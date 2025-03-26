import argparse
import os
import numpy as np
import tables as tb
import matplotlib.pyplot as plt
from scipy.signal import find_peaks_cwt

import invisible_cities.core.fit_functions  as fitf
import invisible_cities.io.channel_param_io as pIO

from invisible_cities.calib.spe_response import poisson_scaled_gaussians, dark_convolution
from invisible_cities.database import load_db as DB
from invisible_cities.calib.calib_functions import seeds_and_bounds, dark_scaler, SensorType
from invisible_cities.cities.components import get_run_number

# Utility function for weighted average and standard deviation
def weighted_av_std(values, weights):
    avg = np.average(values, weights=weights)
    var = np.average((values - avg) ** 2, weights=weights)
    var = weights.sum() * var / (weights.sum() - 1)
    return avg, np.sqrt(var)

# Fit the dark spectrum with a Gaussian to estimate pedestal parameters
def fit_dark_spectrum(bins, dspec):
    gb0 = [(0, -1, 0), (1e99, 1, 10)]  # Bounds: norm, mean, sigma
    av, rms = weighted_av_std(bins[dspec > 100], dspec[dspec > 100])
    sd0 = (dspec.sum(), av, rms)  # Seeds: norm, mean, sigma
    errs = np.sqrt(dspec[dspec > 100])
    errs[errs == 0] = 0.0001  # Avoid zero errors
    gfitRes = fitf.fit(fitf.gauss, bins[dspec > 100], dspec[dspec > 100], sd0, sigma=errs, bounds=gb0)
    return gfitRes.values, gfitRes.errors, gfitRes

# Plot the dark spectrum fit
def plot_dark_fit(bins, dspec, gfitRes, ich, rnum):
    plt.figure()
    # Plot all data initially, but we'll zoom in
    plt.errorbar(bins, dspec, yerr=np.sqrt(dspec), fmt='b.', label='Dark Data')
    fit_x = np.linspace(min(bins), max(bins), 1000)  # Full range for fit curve
    fit_y = fitf.gauss(fit_x, *gfitRes.values)
    plt.plot(fit_x, fit_y, 'r-', label='Gaussian Fit')
    plt.title(f'Dark Spectrum Fit for Channel {ich}')
    plt.xlabel('ADC')
    plt.ylabel('Counts')
    plt.yscale('log')

    # Extract fit parameters (Norm, Mean, Sigma)
    ped_mean = gfitRes.values[1]  # Mean of the pedestal peak
    ped_sigma = gfitRes.values[2]  # Sigma of the pedestal peak

    # Dynamically set x-axis limits
    x_min = ped_mean - 3 * ped_sigma  # 3 sigma covers ~99.7% of the peak
    x_max = ped_mean + 3 * ped_sigma
    plt.xlim(x_min, x_max)

    # Dynamically set y-axis limits within the x-range
    sel = (bins >= x_min) & (bins <= x_max)
    y_min = 0.1  # Minimum for log scale
    y_max = max(dspec[sel]) * 1.5 if max(dspec[sel]) > 0 else 10  # Headroom above max count
    plt.ylim(y_min, y_max)

    # Legend with fit parameters and chi-squared
    param_names = ['Norm', 'Mean', 'Sigma']
    params_str = '\n'.join([f'{name}: {val:.2f} ± {err:.2f}' 
                            for name, val, err in zip(param_names, gfitRes.values, gfitRes.errors)])
    chi2_str = f'Chi2: {gfitRes.chi2:.2f}'
    plt.legend(title=f'Fit Parameters\n{params_str}\n{chi2_str}', loc='upper right')

    # Save plot
    plot_dir = f"sipm_{rnum}/dark"
    if not os.path.isdir(plot_dir):
        os.makedirs(plot_dir)
    plt.savefig(f"{plot_dir}/dark_fit_{ich}.pdf", bbox_inches='tight')
    plt.close()

# Compute seeds and bounds for the fit
def compute_seeds_and_bounds(sensor_type, run_no, ich, func_name, bins, lspec, dspec, ped_vals, ped_errs):
    detector = 'demopp'  # Adjust if necessary
    use_db_gain_seeds = False  # Use database gain seeds
    scaler_func = dark_scaler(dspec[(bins>=-5) & (bins<=5)]) # Note: the selection of bins in dspec is necessary or seeds_and_bounds will fail
    seeds, bounds = seeds_and_bounds(sensor_type, run_no, ich, scaler_func, bins, lspec, ped_vals, detector, ped_errs, func=func_name, use_db_gain_seeds=use_db_gain_seeds)
    return seeds, bounds

# Perform the fit with the specified number of Gaussians
def perform_fit(fit_func, bins, lspec, seeds, bounds):
    errs = np.sqrt(lspec)
    errs[errs == 0] = 1  # Minimum error for empty bins
    rfit = fitf.fit(fit_func, bins, lspec, seeds, sigma=errs, bounds=bounds)
    return rfit

# Plot the spectrum and fit with parameters and chi-squared in the legend
def plot_fit(bins, lspec, rfit, ich, rnum, ped_mean, ped_sigma):
    plt.figure()
    plt.errorbar(bins, lspec, xerr=0.5 * np.diff(bins)[0], yerr=np.sqrt(lspec), fmt='b.', label='Data')
    plt.plot(bins, rfit.fn(bins), 'r-', label='Fit')
    plt.yscale('log')
    plt.title(f'SPE Response Fit to Channel {ich}')
    plt.xlabel('ADC')
    plt.ylabel('Counts')

    # Extract fit parameters from convolution model
    normalization = rfit.values[0]  # Scaling factor
    poisson_mean = rfit.values[1]   # Poisson mean
    gain = rfit.values[2]           # SPE gain
    gain_sigma = rfit.values[3]     # SPE gain sigma

    # Use pedestal mean and sigma from dark spectrum fit
    x_min = ped_mean - 2 * ped_sigma  # Lower limit: 2 sigma below pedestal
    x_max = ped_mean + 5 * gain       # Upper limit: 5 gains above pedestal
    plt.xlim(x_min, x_max)

    # Dynamically set y-axis limits within the x-range
    sel = (bins >= x_min) & (bins <= x_max)
    y_min = 0.1
    y_max = max(lspec[sel]) * 1.5 if max(lspec[sel]) > 0 else 10
    plt.ylim(y_min, y_max)

    # Legend with fit parameters and chi-squared
    param_names = ['Normalization', 'Poisson Mean', 'Gain', 'Gain Sigma']
    params_str = '\n'.join([f'{name}: {val:.3f} ± {err:.3f}' 
                            for name, val, err in zip(param_names, rfit.values, rfit.errors)])
    chi2_str = f'Chi2: {rfit.chi2:.2f}'
    plt.legend(title=f'Fit Parameters\n{params_str}\n{chi2_str}', loc='lower right')

    plot_dir = f"sipm_{rnum}"
    if not os.path.isdir(plot_dir):
        os.mkdir(plot_dir)
    plt.savefig(f"{plot_dir}/fit_{ich}.pdf", bbox_inches='tight')
    plt.close()

# Save all generic parameters
def save_fit_params(param_writer, ich, rfit, ped_vals, ped_errs, bins, b1, b2):
    outDict = {
        "normalization": (rfit.values[0], rfit.errors[0]),     # From convolution fit
        "poisson_mu": (rfit.values[1], rfit.errors[1]),        # From convolution fit
        "pedestal": (ped_vals[1], ped_errs[1]),                # From dark spectrum fit
        "pedestal_sigma": (ped_vals[2], ped_errs[2]),          # From dark spectrum fit
        "gain": (rfit.values[2], rfit.errors[2]),              # From convolution fit
        "gain_sigma": (rfit.values[3], rfit.errors[3]),        # From convolution fit
        "fit_limits": (bins[b1], bins[min(len(bins)-1, b2)]),  # Fit range
        "n_gaussians_chi2": (0, rfit.chi2)                     # No n_gaussians, just chi2
    }
    param_writer(ich, outDict)

# Generate 2D scatter plots for gain and chi-squared
def plot_2d_scatter(pos_x, pos_y, gains, chi2s, rnum):
    plot_dir = f"sipm_{rnum}"
    if not os.path.isdir(plot_dir):
        os.mkdir(plot_dir)

    # Gain scatter plot
    plt.scatter(pos_x, pos_y, c=gains)
    plt.title("Gain Map")
    plt.xlabel("X (mm)")
    plt.ylabel("Y (mm)")
    plt.colorbar(label="Gain (ADC)")
    plt.savefig(f"{plot_dir}/gain_map.pdf")
    plt.close()

    # Chi-squared scatter plot
    plt.scatter(pos_x, pos_y, c=chi2s)
    plt.title("Chi-squared Map")
    plt.xlabel("X (mm)")
    plt.ylabel("Y (mm)")
    plt.colorbar(label="Chi^2")
    plt.savefig(f"{plot_dir}/chi2_map.pdf")
    plt.close()

# Main function to orchestrate the fitting process
def main():

    # Parse command-line arguments
    parser = argparse.ArgumentParser(description='Fit SiPM spectra with n Gaussians.')
    parser.add_argument('file_in', type=str, help='Input file (e.g., sipmCal_R<run#>.h5)')
    parser.add_argument('--rng-low', type=int, default=None, help='Lower channel number range (inclusive)')
    parser.add_argument('--rng-high', type=int, default=None, help='Upper channel number range (inclusive)')
    #parser.add_argument('n_gaussians', type=int, help='Number of Gaussians for the fit (e.g., 2 or 3)')
    args = parser.parse_args()

    file_name = args.file_in
    rng_low = args.rng_low
    rng_high = args.rng_high
    # n_gaussians = args.n_gaussians
    # if n_gaussians < 2:
    #     raise ValueError("Number of Gaussians must be >= 2")

    # Load data
    sipmIn = tb.open_file(file_name, 'r')
    bins = np.array(sipmIn.root.HIST.sipm_spe_bins)
    specsL = np.array(sipmIn.root.HIST.sipm_spe).sum(axis=0)
    specsD = np.array(sipmIn.root.HIST.sipm_dark).sum(axis=0)
    run_no = get_run_number(sipmIn)
    sensor_type = SensorType.SIPM

    # Database file (kept as is)
    db_file = '/home/jrenner/local/jerenner/IC/invisible_cities/database/localdb.DEMOPPDB.sqlite3'
    channs = DB.DataSiPM(db_file, run_no).SensorID.values

    # Set up output file
    func_name = 'conv'
    #func_name = f'ngau_{n_gaussians}'
    out_file = tb.open_file(f'sipmCalParOut_R{run_no}_F{func_name}.h5', 'w')
    param_writer = pIO.channel_param_writer(out_file, sensor_type='sipm', func_name=func_name, param_names=pIO.generic_params)

    # Lists to store results for 2D plots
    gains = []
    chi2s = []

    # Process each SiPM channel
    for ich, (dspec, lspec) in enumerate(zip(specsD, specsL)):

        channel_id = channs[ich]
        # Skip channels outside the specified range if rng_low or rng_high is provided
        if (rng_low is not None and channel_id < rng_low) or (rng_high is not None and channel_id > rng_high):
            #print(f"Skipping sensor {channel_id} (outside range {rng_low}-{rng_high})")
            continue

        print(f"Processing sensor {ich}...")

        # Fit dark spectrum to estimate pedestal
        ped_vals, ped_errs, gfitRes = fit_dark_spectrum(bins[(bins>=-5) & (bins<=5)], dspec[(bins>=-5) & (bins<=5)])
        ped_mean = ped_vals[1]  # Mean from dark spectrum fit
        ped_sigma = ped_vals[2] # Sigma from dark spectrum fit
        plot_dark_fit(bins, dspec, gfitRes, channel_id, run_no)

        # Compute seeds and bounds
        seeds, bounds = compute_seeds_and_bounds(sensor_type, run_no, ich, func_name, bins, lspec, dspec, ped_vals, ped_errs)
        seeds = (seeds[0], 0.1, 18, 4)
        bounds = ((0, 0.099, 10, 2), (bounds[1][0], 0.101, 40, 4))
        print(f"Seeds: {seeds}")
        print(f"Bounds: {bounds}")
        mask = (bins >= ped_mean - 5 * ped_sigma) & (bins <= ped_mean + 4 * seeds[2])
        bins_subset = bins[mask]
        lspec_subset = lspec[mask]
        dspec_subset = dspec[mask]

        # Perform the fit
        # Define the fitting function
        # if('ngau' in func_name):
        #     fit_func = poisson_scaled_gaussians(n_gaussians=n_gaussians)
        # else:
        fit_func = dark_convolution(bins_subset, dspec_subset, 1)
        rfit = perform_fit(fit_func, bins_subset, lspec_subset, seeds, bounds)
        #rfit = perform_fit(fit_func, bins, lspec, seeds, bounds)

        # Plot the results
        plot_fit(bins_subset, lspec_subset, rfit, channel_id, run_no, ped_mean, ped_sigma)

        # Save fit parameters
        save_fit_params(param_writer, channel_id, rfit, ped_vals, ped_errs, bins, 0, len(bins)-1)

        # Store gain and chi-squared for scatter plots
        gains.append(rfit.values[2])  # Gain
        chi2s.append(rfit.chi2)

        print("\n")

    # Generate 2D scatter plots
    if(rng_low is None and rng_high is None):
        pos_x = DB.DataSiPM(db_file, run_no).X.values
        pos_y = DB.DataSiPM(db_file, run_no).Y.values
        plot_2d_scatter(pos_x, pos_y, gains, chi2s, run_no)

    # Clean up
    sipmIn.close()
    out_file.close()

if __name__ == '__main__':
    main()
