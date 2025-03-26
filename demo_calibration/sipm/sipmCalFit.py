import os
import sys
import tables as tb
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.signal import find_peaks
import sqlite3

import invisible_cities.io.channel_param_io as pIO

# Constants
db_path = "/home/jrenner/local/jerenner/IC/invisible_cities/database/localdb.DEMOPPDB.sqlite3"
base_ids = [12000, 13000, 14000, 15000]

# Command-line argument handling
if len(sys.argv) != 2:
    print("Usage: python calib_sipm_spectra.py run_number")
    sys.exit(1)
run_number = int(sys.argv[1])

# Set up input and output paths
input_file = f"sipmSpecs_{run_number}.h5"
output_dir = f"sipm_{run_number}"
os.makedirs(output_dir, exist_ok=True)

# Read input HDF5 file
with tb.open_file(input_file, "r") as h5in:
    spe_table = h5in.root.HIST.sipm_dark[0]  # Dark spectra for all SiPMs
    bin_edges = h5in.root.HIST.sipm_dark_bins[:]  # Bin edges or centers

    # Determine bin centers
    if len(bin_edges) == spe_table.shape[1]:
        bin_centers = bin_edges
    else:
        bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    nsipm = spe_table.shape[0]  # Number of SiPMs

    # Load the DataSiPM dataset for SensorID mapping
    data_sipm_table = h5in.root.Sensors.DataSiPM
    sensor_ids = [row['sensorID'] for row in data_sipm_table]  # List of SensorIDs for channels 0 to 255

# Query the database for SiPM information
conn = sqlite3.connect(db_path)
query = """
    SELECT SensorID, X, Y
    FROM ChannelPosition
    WHERE ? BETWEEN MinRun AND MaxRun
    AND Type = 'Anode'
"""
try:
    sipm_df = pd.read_sql_query(query, conn, params=(run_number,))
    sipm_df = sipm_df.reset_index().rename(columns={"index": "channel"})
except:
    print(f"Failed to query ChannelPosition for run {run_number}. Checking default approach.")
    query = "SELECT SensorID, X, Y FROM ChannelPosition WHERE Type = 'Anode'"
    sipm_df = pd.read_sql_query(query, conn)
    sipm_df = sipm_df.reset_index().rename(columns={"index": "channel"})
conn.close()

# Verify the number of SiPMs matches the HDF5 data
if len(sipm_df) != nsipm:
    print(f"Warning: Database contains {len(sipm_df)} SiPMs, but HDF5 has {nsipm}. Using available data.")

# Define single-Gaussian function
def gaussian(x, amplitude, mean, sigma):
    """Single Gaussian function."""
    return amplitude * np.exp(-((x - mean) ** 2) / (2 * sigma ** 2))

# Function to fit the spectrum
def fit_spectrum(bin_centers, counts):
    """Fit the dark spectrum to two independent single Gaussians with improved peak finding.
    Returns parameters, errors, and chi-squared for each Gaussian separately.
    """
    # Smooth the data to reduce noise (optional)
    counts_smooth = np.convolve(counts, np.ones(3)/3, mode='same')

    # Initial peak detection for noise peak
    peaks, properties = find_peaks(counts_smooth, height=0.5 * np.max(counts), distance=5, prominence=50)
    if len(peaks) < 1:
        return (-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1)  # No significant peaks found

    # Identify noise peak (highest peak)
    noise_peak = peaks[np.argmax(counts_smooth[peaks])]
    noise_mean = bin_centers[noise_peak]

    # Initial fit for noise peak
    p0_noise = [counts_smooth[noise_peak], noise_mean, 2]
    bounds_noise = ([0, -5, 0.5], [np.inf, 5, 5])
    try:
        popt_noise, pcov_noise = curve_fit(gaussian, bin_centers, counts_smooth, p0=p0_noise,
                                           bounds=bounds_noise, maxfev=5000)
        perr_noise = np.sqrt(np.diag(pcov_noise))
    except (RuntimeError, ValueError):
        return (-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1)

    # Define a range for the noise peak (e.g., ±3σ)
    noise_sigma = popt_noise[2]
    noise_range = (bin_centers >= popt_noise[1] - 3 * noise_sigma) & (bin_centers <= popt_noise[1] + 3 * noise_sigma)

    # Calculate chi-squared for noise peak in its range
    fit_noise = gaussian(bin_centers[noise_range], *popt_noise)
    chi2_noise = np.sum(((counts[noise_range] - fit_noise) ** 2) / (counts[noise_range] + 1e-6))

    # Subtract noise fit to isolate SPE peak
    residual = counts - gaussian(bin_centers, *popt_noise)

    # Find SPE peak in residual
    spe_peaks, _ = find_peaks(residual, height=0.2 * np.max(residual), distance=2, prominence=0.02 * np.max(residual), width=1)
    if len(spe_peaks) < 1:
        return (popt_noise[0], popt_noise[1], popt_noise[2], -1, -1, -1, perr_noise[0], perr_noise[1], perr_noise[2], -1, -1, -1, chi2_noise, -1)

    # Select the maximum peak after the noise tail
    valid_spe_peaks = spe_peaks[bin_centers[spe_peaks] > noise_mean + 5]
    if len(valid_spe_peaks) == 0:
        spe_peak = np.argmin(np.abs(bin_centers - (noise_mean + 15)))  # Default to 15 ADC offset
    else:
        spe_peak = valid_spe_peaks[np.argmax(residual[valid_spe_peaks])]
    spe_mean = bin_centers[spe_peak]

    # Construct a cut for fitting after the noise tail
    fit_rng = (bin_centers > (spe_mean - 5)) & (bin_centers < (spe_mean + 5))

    # Fit SPE peak
    p0_spe = [max(1, residual[spe_peak]), spe_mean, 2]
    bounds_spe = ([0, 5, 0.5], [np.inf, 50, 4])
    try:
        popt_spe, pcov_spe = curve_fit(gaussian, bin_centers[fit_rng], residual[fit_rng], p0=p0_spe,
                                       bounds=bounds_spe, maxfev=5000)
        perr_spe = np.sqrt(np.diag(pcov_spe))
    except (RuntimeError, ValueError):
        return (popt_noise[0], popt_noise[1], popt_noise[2], -1, -1, -1, perr_noise[0], perr_noise[1], perr_noise[2], -1, -1, -1, chi2_noise, -1)

    # Calculate chi-squared for SPE peak in its fit range
    fit_spe = gaussian(bin_centers[fit_rng], *popt_spe)
    chi2_spe = np.sum(((residual[fit_rng] - fit_spe) ** 2) / (residual[fit_rng] + 1e-6))

    # Return parameters, errors, and chi-squared for both fits
    return (popt_noise[0], popt_noise[1], popt_noise[2], popt_spe[0], popt_spe[1], popt_spe[2],
            perr_noise[0], perr_noise[1], perr_noise[2], perr_spe[0], perr_spe[1], perr_spe[2], chi2_noise, chi2_spe)

# Function to save fit parameters
def save_fit_params(param_writer, sensor_id, popt, perr, chi2, bin_centers):
    """Save the fit parameters in the correct format using pIO.channel_param_writer.

    Parameters:
    - param_writer: Function to write the parameters (e.g., pIO.channel_param_writer)
    - sensor_id: Identifier for the sensor/channel
    - popt: Array of fit parameters [amp_noise, mean_noise, sigma_noise, amp_spe, mean_spe, sigma_spe]
    - perr: Array of parameter errors [amp_noise_err, mean_noise_err, sigma_noise_err, amp_spe_err, mean_spe_err, sigma_spe_err]
    - chi2: Chi-squared value of the fit
    - bin_centers: Array of bin centers used in the fit
    """
    # Check for failed fit
    if popt[0] == -1 or popt[3] == -1:
        return

    # Extract parameters
    amp_noise, mean_noise, sigma_noise, amp_spe, mean_spe, sigma_spe = popt[:6]
    amp_noise_err, mean_noise_err, sigma_noise_err, amp_spe_err, mean_spe_err, sigma_spe_err = perr[:6]

    # Map to required keys
    normalization = (amp_noise, amp_noise_err)  # Amplitude of noise peak
    pedestal = (mean_noise, mean_noise_err)     # Mean of noise peak
    pedestal_sigma = (sigma_noise, sigma_noise_err)  # Sigma of noise peak

    # Calculate gain as mean_spe - mean_noise
    gain = mean_spe - mean_noise
    gain_error = np.sqrt(mean_spe_err**2 + mean_noise_err**2)

    # Gain sigma is the SPE peak sigma
    gain_sigma = (sigma_spe, sigma_spe_err)

    # Estimate poisson_mu as amp_spe / amp_noise
    if amp_noise != 0:
        poisson_mu = amp_spe / amp_noise
        poisson_mu_error = poisson_mu * np.sqrt((amp_spe_err / amp_spe)**2 + (amp_noise_err / amp_noise)**2)
    else:
        poisson_mu = -1
        poisson_mu_error = -1

    # Fit limits from bin_centers
    fit_limits = (bin_centers[0], bin_centers[-1])

    # Number of Gaussians is 2
    n_gaussians = 2

    # Construct the output dictionary
    outDict = {
        "normalization": normalization,           # Amplitude of noise peak
        "poisson_mu": (poisson_mu, poisson_mu_error),  # Estimated from amplitude ratio
        "pedestal": pedestal,                     # Mean of noise peak
        "pedestal_sigma": pedestal_sigma,         # Sigma of noise peak
        "gain": (gain, gain_error),               # Gain = spe_mean - pedestal
        "gain_sigma": gain_sigma,                 # Sigma of SPE peak
        "fit_limits": fit_limits,                 # Fit range
        "n_gaussians_chi2": (n_gaussians, chi2)   # Number of Gaussians and chi2
    }

    # Write the parameters
    param_writer(sensor_id, outDict)

# Function to plot the fit with parameters in legend
def plot_fit(bin_centers, counts, popt, perr, chi2_noise, chi2_spe, sensor_id, output_dir):
    """Generate and save a plot of the spectrum with fit and parameters in legend."""
    plt.figure(figsize=(8, 6))
    plt.step(bin_centers, counts, where="mid", label="Dark Spectrum")
    if popt[0] != -1 and popt[3] != -1:  # If both fits succeeded
        x_fit = np.linspace(bin_centers[0], 50, 1000)
        y_fit = gaussian(x_fit, *popt[:3]) + gaussian(x_fit, *popt[3:6])
        plt.plot(x_fit, y_fit, 'r-', label="Sum of Fits")
        plt.plot(x_fit, gaussian(x_fit, *popt[:3]), 'g--', label="Noise Fit")
        plt.plot(x_fit, gaussian(x_fit, *popt[3:6]), 'b--', label="SPE Fit")
    elif popt[0] != -1:  # Only noise fit succeeded
        x_fit = np.linspace(bin_centers[0], 50, 1000)
        y_fit = gaussian(x_fit, *popt[:3])
        plt.plot(x_fit, y_fit, 'r-', label="Noise Fit")
    plt.xlim(-10, 50)
    plt.xlabel("Integrated Charge (ADC counts)")
    plt.ylabel("Counts")
    plt.title(f"SiPM {sensor_id} Dark Spectrum Fit")

    # Prepare legend with parameters and chi-squared
    if popt[0] != -1 and popt[3] != -1:
        params_str = (f"Noise: A={popt[0]:.2f}±{perr[0]:.2f}, μ={popt[1]:.2f}±{perr[1]:.2f}\n"
                      f"σ={popt[2]:.2f}±{perr[2]:.2f}, $\\chi^2$ = {chi2_noise:.2f}\n"
                      f"SPE: A={popt[3]:.2f}±{perr[3]:.2f}, μ={popt[4]:.2f}±{perr[4]:.2f}\n"
                      f"σ={popt[5]:.2f}±{perr[5]:.2f}, $\\chi^2$ = {chi2_spe:.2f}\n"
                      f"Gain: {popt[4] - popt[1]:.2f}±{np.sqrt(perr[4]**2 + perr[1]**2):.2f}")
    elif popt[0] != -1:
        params_str = (f"Noise: A={popt[0]:.2f}±{perr[0]:.2f}, μ={popt[1]:.2f}±{perr[1]:.2f}, σ={popt[2]:.2f}±{perr[2]:.2f}\n"
                      "SPE: Fit failed\nChi2: N/A")
    else:
        params_str = "Fit failed"

    plt.legend(title=f"Fit Parameters\n{params_str}", loc='upper right')
    plt.savefig(os.path.join(output_dir, f"fit_{sensor_id}.pdf"))
    plt.close()

# Set up HDF5 output file for parameters
func_name = 'two_gaussians'
out_file = tb.open_file(f'sipmCalParOut_R{run_number}_F{func_name}.h5', 'w')
param_writer = pIO.channel_param_writer(out_file, sensor_type='sipm', func_name=func_name, param_names=pIO.generic_params)

# Process each SiPM
fit_results = []
for channel in range(nsipm):
    sensor_id = sensor_ids[channel]  # Use DataSiPM mapping
    print(f"[sipm {channel}/{nsipm}]: SensorID {sensor_id}")
    spectrum = spe_table[channel]  # Dark spectrum
    popt = fit_spectrum(bin_centers, spectrum)
    if len(popt) == 14:  # Successful fit
        amp_noise, mean_noise, sigma_noise, amp_spe, mean_spe, sigma_spe, \
        amp_noise_err, mean_noise_err, sigma_noise_err, amp_spe_err, mean_spe_err, sigma_spe_err, chi2_noise, chi2_spe = popt
        fit_params = [amp_noise, mean_noise, sigma_noise, amp_spe, mean_spe, sigma_spe]
        fit_errors = [amp_noise_err, mean_noise_err, sigma_noise_err, amp_spe_err, mean_spe_err, sigma_spe_err]
    else:
        print(f"WARNING: invalid fit")
        fit_params = [-1] * 6
        fit_errors = [-1] * 6
        chi2_noise = -1
        chi2_spe = -1

    # Save fit parameters
    save_fit_params(param_writer, sensor_id, fit_params, fit_errors, chi2_noise, bin_centers)

    # Generate fit plot
    plot_fit(bin_centers, spectrum, fit_params, fit_errors, chi2_noise, chi2_spe, sensor_id, output_dir)

    # Store results for CSV
    gain = (mean_spe - mean_noise) if mean_spe != -1 and mean_noise != -1 else np.nan
    fit_results.append({
        "sensor_id": sensor_id,
        "mean_noise": mean_noise if mean_noise != -1 else np.nan,
        "sigma_noise": sigma_noise if sigma_noise != -1 else np.nan,
        "mean_spe": mean_spe if mean_spe != -1 else np.nan,
        "sigma_spe": sigma_spe if sigma_spe != -1 else np.nan,
        "gain": gain
    })

# Save fit results to CSV
fit_df = pd.DataFrame(fit_results)
fit_df.to_csv(os.path.join(output_dir, "fit_results.csv"), index=False)

# Close the HDF5 file
out_file.close()

# Filter out failed fits for summary plots
fit_df_valid = fit_df[fit_df["gain"].notna()]

# Generate DICE board histograms
db_defs = {0: 12000, 1: 13000, 2: 14000, 3: 15000}
for db_num, base_id in db_defs.items():
    db_sensors = fit_df_valid[fit_df_valid["sensor_id"] // 1000 == base_id // 1000]
    if db_sensors.empty:
        continue
    
    plt.figure(figsize=(8, 6))
    x_values = db_sensors["sensor_id"] - base_id
    y_values = db_sensors["gain"]
    plt.bar(x_values, y_values, width=0.8)
    plt.xlabel(f"SensorID - BaseID({base_id})")
    plt.ylabel("Single-Photon Gain (ADC counts)")
    plt.title(f"DICE Board {db_num} Gains")
    plt.xlim(-1, 64)  # Assuming 64 sensors per DB
    plt.savefig(os.path.join(output_dir, f"hist_gains_DB{db_num}.pdf"))
    plt.close()

# Generate spatial gain plot
merged_df = pd.merge(sipm_df, fit_df_valid, left_on="SensorID", right_on="sensor_id")
merged_df = merged_df[["SensorID", "X", "Y", "gain"]]  # Keep relevant columns
plt.figure(figsize=(10, 10))
scatter = plt.scatter(
    merged_df["X"], merged_df["Y"], c=merged_df["gain"],
    s=100, marker="s", cmap="viridis", edgecolors="black"
)
plt.colorbar(scatter, label="Single-Photon Gain (ADC counts)")
plt.xlabel("X Position (mm)")
plt.ylabel("Y Position (mm)")
plt.title(f"SiPM Gains vs. (X, Y) Position - Run {run_number}")
plt.axis("equal")
plt.savefig(os.path.join(output_dir, "gains_xy.pdf"))
plt.close()

print(f"Processing complete. Outputs saved in {output_dir}")