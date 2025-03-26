import os
import sys
import numpy as np
import tables as tb
import matplotlib.pyplot as plt

import invisible_cities.core.fit_functions       as fitf
import invisible_cities.calib.spe_response       as speR
import invisible_cities.io.channel_param_io      as pIO

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
    gb0 = [(0, -100, 0), (1e99, 100, 10000)]  # Bounds: norm, mean, sigma
    av, rms = weighted_av_std(bins[dspec > 100], dspec[dspec > 100])
    sd0 = (dspec.sum(), av, rms)  # Seeds: norm, mean, sigma
    errs = np.sqrt(dspec[dspec > 100])
    errs[errs == 0] = 0.0001  # Avoid zero errors
    gfitRes = fitf.fit(fitf.gauss, bins[dspec > 100], dspec[dspec > 100], sd0, sigma=errs, bounds=gb0)
    return gfitRes

# Plot the fit of the dark spectrum
def plot_dark_fit(bins, dspec, gfitRes, ich, rnum):
    plt.figure()
    sel = dspec > 100
    plt.errorbar(bins[sel], dspec[sel], yerr=np.sqrt(dspec[sel]), fmt='b.', label='Dark Data')
    fit_x = np.linspace(bins[sel][0], bins[sel][-1], 1000)
    fit_y = fitf.gauss(fit_x, *gfitRes.values)
    plt.plot(fit_x, fit_y, 'r-', label='Gaussian Fit')
    plt.title(f'Dark Spectrum Fit for Channel {ich}')
    plt.xlabel('ADC')
    plt.ylabel('Counts')

    param_names = ['Norm', 'Mean', 'Sigma']
    params_str = '\n'.join([f'{name}: {val:.2f} ± {err:.2f}' 
                            for name, val, err in zip(param_names, gfitRes.values, gfitRes.errors)])
    chi2_str = f'Chi2: {gfitRes.chi2:.2f}'
    plt.legend(title=f'Fit Parameters\n{params_str}\n{chi2_str}', loc='upper right')

    plot_dir = f"{rnum}_pmt/dark"
    if not os.path.isdir(plot_dir):
        os.makedirs(plot_dir)
    plt.savefig(f"{plot_dir}/dark_fit_ch{ich}.pdf", bbox_inches='tight')
    plt.close()

# Compute seeds and bounds for the fit
def compute_seeds_and_bounds(sensor_type, run_no, ich, bins, lspec, dspec, ped_vals, ped_errs, n_gaussians):
    detector = 'demopp'  # Adjust if necessary
    use_db_gain_seeds = False  # Use database gain seeds; set to False if preferred
    func_name = f'ngau_{n_gaussians}'
    scaler_func = dark_scaler(dspec[bins < 0])  # Precompute scaler with dark spectrum
    seeds, bounds = seeds_and_bounds(sensor_type, run_no, ich, scaler_func, bins, lspec, ped_vals, detector, ped_errs, func=func_name, use_db_gain_seeds=use_db_gain_seeds)
    return seeds, bounds

# Perform the fit with the specified number of Gaussians
def perform_fit(fit_func, bins, lspec, seeds, bounds):
    errs = np.sqrt(lspec)
    errs[errs == 0] = 1  # Minimum error for empty bins
    rfit = fitf.fit(fit_func, bins, lspec, seeds, sigma=errs, bounds=bounds)
    return rfit

# Plot the spectrum and fit, including parameters and chi-squared in the legend
def plot_fit(bins, lspec, rfit, ich, rnum):
    plt.figure()
    # Plot data with error bars
    plt.errorbar(bins, lspec, xerr=0.5 * np.diff(bins)[0], yerr=np.sqrt(lspec), fmt='b.', label='Data')
    plt.plot(bins, rfit.fn(bins), 'r-', label='Fit')
    
    # Set logarithmic y-scale
    plt.yscale('log')
    plt.title(f'SPE Response Fit to Channel {ich}')
    plt.xlabel('ADC')
    plt.ylabel('Counts')

    # Extract fit parameters (adjust indices based on your fit model)
    ped_mean = rfit.values[2]  # 3rd parameter is pedestal mean
    ped_sigma = rfit.values[3] # 4th parameter is pedestal sigma
    gain = rfit.values[4]      # 5th parameter is gain

    # Set x-axis limits dynamically
    x_min = ped_mean - 2 * ped_sigma
    x_max = ped_mean + 5 * gain
    plt.xlim(x_min, x_max)

    # Set y-axis limits
    y_min = 0.1  # Avoid zero for log scale
    y_max = max(lspec) * 1.5  # 50% above max count
    plt.ylim(y_min, y_max)

    # Add legend with fit parameters
    param_names = ['Norm', 'Mu', 'Ped_mean', 'Ped_sigma', 'Gain', 'Gain_sigma']
    params_str = '\n'.join([f'{name}: {val:.2f} ± {err:.2f}' 
                           for name, val, err in zip(param_names, rfit.values, rfit.errors)])
    chi2_str = f'Chi2: {rfit.chi2:.2f}'
    plt.legend(title=f'Fit Parameters\n{params_str}\n{chi2_str}', loc='upper right')

    # Save the plot
    plot_dir = f"{rnum}_pmt"
    if not os.path.isdir(plot_dir):
        os.mkdir(plot_dir)
    plt.savefig(f"{plot_dir}/fit_{ich}.pdf", bbox_inches='tight')
    plt.close()

# Save fit parameters using param_writer
def save_fit_params(param_writer, ich, rfit, n_gaussians, bins, b1, b2):
    outDict = {
        "normalization": (rfit.values[0], rfit.errors[0]),     # normalization
        "poisson_mu": (rfit.values[1], rfit.errors[1]),        # poisson_mu
        "pedestal": (rfit.values[2], rfit.errors[2]),          # pedestal mean
        "pedestal_sigma": (rfit.values[3], rfit.errors[3]),    # pedestal sigma
        "gain": (rfit.values[4], rfit.errors[4]),              # gain
        "gain_sigma": (rfit.values[5], rfit.errors[5]),        # gain sigma
        "fit_limits": (bins[b1], bins[min(len(bins)-1, b2)]),  # fit range
        "n_gaussians_chi2": (n_gaussians, rfit.chi2)           # n_gaussians and chi2
    }
    param_writer(ich, outDict)

# Main function to orchestrate the fitting process
def main():
    # Check command-line arguments
    if len(sys.argv) != 3:
        print("Usage: python pmtCalFit.py <run_number> <n_gaussians>")
        sys.exit(1)

    rnum = sys.argv[1]
    n_gaussians = int(sys.argv[2])
    if n_gaussians < 2:
        raise ValueError("Number of Gaussians must be >= 2")

    # Load data
    file_name = f"pmtCal_R{rnum}.h5"
    h5in = tb.open_file(file_name, 'r')
    bins = np.array(h5in.root.HIST.pmt_dark_bins)
    specsD = np.array(h5in.root.HIST.pmt_dark).sum(axis=0)
    specsL = np.array(h5in.root.HIST.pmt_spe).sum(axis=0)
    run_no = get_run_number(h5in)
    data_pmt_table = h5in.root.Sensors.DataPMT
    sensor_ids = [row['sensorID'] for row in data_pmt_table]  # List of SensorIDs for channels 0 to 255
    sensor_type = SensorType.PMT

    # Set up output file
    func_name = f'ngau_{n_gaussians}'
    pOut = tb.open_file(f'pmtCalParOut_R{rnum}_F{func_name}.h5', 'w')
    param_writer = pIO.channel_param_writer(pOut, sensor_type='pmt', func_name=func_name, param_names=pIO.generic_params)

    # Define the fitting function
    fit_func = speR.poisson_scaled_gaussians(n_gaussians=n_gaussians)

    # Process each PMT channel
    for ich in range(3):  # 3 PMT channels
        print(f"Fitting sensor {sensor_ids[ich]}...")
        
        dspec = specsD[ich]
        lspec = specsL[ich]

        # Fit dark spectrum to estimate pedestal
        gfitRes = fit_dark_spectrum(bins, dspec)
        ped_vals = gfitRes.values  # norm, mean, sigma
        ped_errs = gfitRes.errors
        plot_dark_fit(bins, dspec, gfitRes, sensor_ids[ich], rnum)

        # Compute seeds and bounds
        seeds, bounds = compute_seeds_and_bounds(sensor_type, run_no, sensor_ids[ich], bins, lspec, dspec, ped_vals, ped_errs, n_gaussians)

        # Perform the fit
        rfit = perform_fit(fit_func, bins, lspec, seeds, bounds)

        # Plot the results
        plot_fit(bins, lspec, rfit, sensor_ids[ich], rnum)

        # Save fit parameters
        save_fit_params(param_writer, sensor_ids[ich], rfit, n_gaussians, bins, 0, len(bins)-1)

    pOut.close()

if __name__ == '__main__':
    main()
