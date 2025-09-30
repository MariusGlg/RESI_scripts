# import h5py
# import matplotlib.pyplot as plt
# import pandas as pd
# import numpy as np
# from scipy.stats import norm
# from configparser import ConfigParser
# from matplotlib.patches import Circle
# import os
#
# # --- Configuration (Kept for context, assuming config.ini is present) ---
# config = ConfigParser()
# file = "config.ini"
# config.read(file)
# path = config["INPUT_FILES"].get("path", "")
# RESI_file_name = config["INPUT_FILES"].get("RESI_file_name", "RESI.hdf5")
# SMLM_cluster_name = config["INPUT_FILES"].get("SMLM_cluster_name", "Cluster.hdf5")
#
#
# # ---------------------------------------------------------------------
#
#
# def load_file(file_path):
#     """load .hdf5_files"""
#     with h5py.File(file_path, "r") as locs_file:
#         key = list(locs_file.keys())[0]
#         locs = locs_file[str(key)][...]
#         data_pd = pd.DataFrame(locs)
#         return data_pd
#
#
# def plot_combined_gaussians_subplots(df_list, column_name, resi_x, resi_sigma):
#     """
#     Plots the histogram in the first subplot and the two scaled Gaussian
#     curves in the second subplot for direct comparison of resolution.
#     """
#
#     locs_cluster = []  # localizations of all clusters (mean set to 0)
#
#     # 1. Process and Recenter Each Cluster
#     for index, (key, value) in enumerate(df_list.items()):  #i in range(len(df_list)):
#         #cluster_data = df_list[i].copy()
#         cluster_data = value
#
#         # Calculate means for recentering
#         x_mean = cluster_data["x"].mean()
#         y_mean = cluster_data["y"].mean()
#
#         # FIX: Use .loc for explicit, unambiguous assignment
#         cluster_data.loc[:, "x"] = cluster_data["x"] - x_mean
#         cluster_data.loc[:, "y"] = cluster_data["y"] - y_mean
#
#         locs_cluster.append(cluster_data)
#
#     # 2. Combine all re-centered localizations
#     final_df = pd.concat(locs_cluster, ignore_index=True)
#     final_df["x"] = final_df["x"] * 130  # rescale to nm
#     final_df["y"] = final_df["y"] * 130  # rescale to nm
#
#     # 3. Fit the combined cluster distribution (uses std_x as the width)
#     # mu_cluster_x will be near 0 because the data was re-centered
#     mu_cluster_x, std_x = norm.fit(final_df[column_name])
#
#     # 4. Prepare for Plotting (Calculate PDFs and Scaling)
#
#     # Determine the x-range
#     x_values = np.linspace(final_df[column_name].min(), final_df[column_name].max(), 200)
#
#     # Calculate unscaled PDF based on the FITTED width (std_x), not the mean sx
#     p_fit_cluster_unscaled = norm.pdf(x_values, mu_cluster_x, std_x)
#
#     # Calculate RESI PDF for comparison
#     p_fit_resi = norm.pdf(x_values, mu_cluster_x, resi_sigma)
#
#     # Determine scaling factor (Target is RESI peak height)
#     cluster_max_height = np.max(p_fit_cluster_unscaled)
#     resi_max_height = np.max(p_fit_resi)
#     scaling_factor = resi_max_height / cluster_max_height
#
#     # Scale the cluster PDF
#     p_fit_resi = p_fit_resi /scaling_factor
#
#     # --------------------------------------------------------------------
#     # 5. Setup Figure and Subplots
#     # --------------------------------------------------------------------
#     fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(14, 6))
#     fig.suptitle("DNA-PAINT vs. RESI resolution", fontsize=14)
#
#     # --------------------------------------------------------------------
#     # Subplot 1: Scaled Histogram + Scaled Cluster Fit
#     # --------------------------------------------------------------------
#     ax1 = axes[1]
#
#     # Plot the histogram, applying the scaling factor via 'weights'
#     # This aligns the bar height with the scaled PDF peak.
#     final_df['x'].hist(
#         ax=ax1,
#         bins=80,
#         density=True,
#         edgecolor='black',
#         alpha=0.2,
#         color='skyblue',
#
#     )
#
#     # Plot the Scaled Cluster Gaussian
#     ax1.plot(x_values, p_fit_cluster_unscaled, 'b-',
#              linewidth=1.5, label=f'Cluster Fit ($\sigma={std_x:.4f}$)')
#     ax1.plot(x_values, p_fit_resi, 'r-',
#              linewidth=1.5, label=f'RESI Fit ($\sigma={resi_sigma:.4f}$)')
#
#     ax1.set_title("Histogram and Scaled Cluster Fit")
#     ax1.set_xlabel(f"{column_name} Position (Re-centered)")
#     ax1.set_ylabel("Relative Density (Scaled)")
#     #ax1.set_xlim(x_values.min(), x_values.max())
#     ax1.set_xlim(left=-60, right=60)
#
#
#     ax1.grid(axis='y', alpha=0.5)
#     ax1.legend()
#
#     # --------------------------------------------------------------------
#     # Subplot 2: Scatter Plot (x vs. y)
#     # --------------------------------------------------------------------
#     ax2 = axes[0]
#     sampled_df = final_df.sample(frac=0.20, random_state=1)  # select x% random from df for plotting
#     ax2.scatter(
#         x=sampled_df["x"],
#         y=sampled_df["y"],
#         s=1,  # Size of the points
#         alpha=0.7,
#         color='darkblue'
#     )
#     circle = Circle(
#         (0, 0),  # (x, y) center coordinates
#         2.35 * resi_sigma,  # Radius (0.05)
#         color='red',  # Fill color
#         alpha=0.6,  # Transparency
#         fill=True  # Make it a solid circle
#     )
#     ax2.add_patch(circle)
#
#     ax2.set_title("Combined DNA-PAINT Cluster (eGFR monomer) (Center of Mass at [0, 0])")
#     ax2.set_xlabel("X Position (nm, Re-centered)")
#     ax2.set_ylabel("Y Position (nm, Re-centered)")
#     ax2.set_xlim(left=-60, right=60)
#
#     # Add center lines
#     ax2.axhline(0, color='r', linestyle=':', alpha=0.5)
#     ax2.axvline(0, color='r', linestyle=':', alpha=0.5)
#     ax2.set_aspect('equal', adjustable='box')  # Added for better visualization of localization data
#
#     # Ensure plots don't overlap
#     plt.tight_layout()
#     plt.show()
#
#
# try:
#     # Check if files exist before attempting to load
#     full_resi_path = os.path.join(path, RESI_file_name)
#     full_cluster_path = os.path.join(path, SMLM_cluster_name)
#
#     if not os.path.exists(full_resi_path) or not os.path.exists(full_cluster_path):
#         raise FileNotFoundError("Check config.ini path and filenames. Files not found.")
#
#     # Load data from HDF5 files
#     resi = load_file(full_resi_path)
#     cluster = load_file(full_cluster_path)
#
#     grouped = cluster.groupby("group")
#     xy_coordinates = {}
#     # extract all values
#     resi_all = []
#     for group_value, group_df in grouped:
#         xy_coordinates[group_value] = group_df[["x", "y", "sx", "sy"]]
#
#     # Extract mean RESI precision
#     #filtered_df = resi[(resi['lpx'] < 0.0090)]
#     resi_sigma = resi["lpx"].mean()
#     resi_sigma = resi_sigma * 130  # rescale to nm
#
#     # Plot the histogram and both Gaussians in subplots
#     plot_combined_gaussians_subplots(xy_coordinates, "x", 0, resi_sigma)
#
# except FileNotFoundError as fnfe:
#     print(f"File Error: {fnfe}. Ensure config.ini is correct and files exist.")
# except Exception as e:
#     print(f"An unexpected error occurred: {e}")

"""
DNA-PAINT Resolution Analysis Script
@author. Marius Glogger
Optical Imaging Competence Centre, FAU Erlangen

This script analyzes the localization precision of DNA-PAINT clusters and corresponding RESI localizations

It performs the following steps:
1. Loads localization data for SMLM clusters and RESI precision (SMLM cluster_center data) data from HDF5 files
   specified in 'config.ini'.
2. Recalculates and re-centers all cluster localizations relative to their center of mass.
3. Rescales all data to nanometers (nm).
4. Fits the combined re-centered localization distribution with a Gaussian curve to
   determine sigma_cluster.
5. Generates a two-panel plot:
    - Left: Scatter plot of re-centered localizations vs. a circle representing RESI precision.
    - Right: Histogram of localizations overlaid with a gaussian of width sigma and
      a gaussian of width sigma_RESI (scaled for comparison).
"""

import h5py
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.stats import norm
from configparser import ConfigParser
from matplotlib.patches import Circle
import os
import sys  # Added for graceful exit on file error


# --- Configuration Management ---
def load_config(file_path="config.ini"):
    """Loads configuration from the specified INI file."""
    config = ConfigParser()
    if not config.read(file_path):
        raise FileNotFoundError(f"Configuration file '{file_path}' not found.")

    path = config["INPUT_FILES"].get("path", "")
    resi_file_name = config["INPUT_FILES"].get("RESI_file_name", "RESI.hdf5")
    smlm_cluster_name = config["INPUT_FILES"].get("SMLM_cluster_name", "Cluster.hdf5")

    return path, resi_file_name, smlm_cluster_name, config


# --- Data Loading Utility ---
def load_hdf5_data(file_path):
    """
    Loads data from the first dataset in an .hdf5 file into a pandas DataFrame.
    """
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"Data file not found at: {file_path}")

    try:
        with h5py.File(file_path, "r") as locs_file:
            # Assumes the first key is the dataset name
            key = next(iter(locs_file.keys()))
            locs = locs_file[str(key)][()]  # Use [()] for cleaner NumPy array extraction
            return pd.DataFrame(locs)
    except StopIteration:
        raise ValueError(f"HDF5 file '{file_path}' is empty or contains no datasets.")
    except Exception as e:
        raise IOError(f"Error reading HDF5 file '{file_path}': {e}")


# --- Plotting Function (Refactored) ---
def plot_combined_gaussians_subplots(cluster_df, resi_sigma_nm, column_name="x", rescale_factor=130):
    """
    Plots the SMLM cluster localization scatter and compares the cluster
    localization distribution to the RESI precision (Gaussian fit).

    Args:
        cluster_df (pd.DataFrame): Combined and grouped SMLM localization data.
        resi_sigma_nm (float): Mean RESI localization precision (sigma) in nm.
        column_name (str): Column ('x' or 'y') to analyze the distribution of.
        rescale_factor (int): Factor to convert raw units to nanometers.
    """

    # 1. Recenter and Rescale All Cluster Localizations

    # Use list comprehension for cleaner and potentially faster processing
    re_centered_dfs = []

    # Process each cluster group
    for group_id, group_df in cluster_df.groupby("group"):
        # Calculate means
        x_mean = group_df["x"].mean()
        y_mean = group_df["y"].mean()

        # Recenter and rescale in one go, creating a new DataFrame slice
        re_centered_df = pd.DataFrame({
            "x": (group_df["x"] - x_mean) * rescale_factor,
            "y": (group_df["y"] - y_mean) * rescale_factor
        })
        re_centered_dfs.append(re_centered_df)

    # Combine all re-centered localizations
    final_df = pd.concat(re_centered_dfs, ignore_index=True)

    # 2. Fit the combined cluster distribution
    # mu_cluster will be near 0
    mu_cluster, std_cluster = norm.fit(final_df[column_name])

    # 3. Prepare for Plotting (Calculate PDFs and Scaling)

    # Determine the x-range for the plot based on the data
    data_min = final_df[column_name].min()
    data_max = final_df[column_name].max()
    x_range = np.linspace(min(-60, data_min), max(60, data_max), 200)  # Use a wider range for smoothness

    # Calculate unscaled PDFs
    p_fit_cluster = norm.pdf(x_range, mu_cluster, std_cluster)
    p_fit_resi_unscaled = norm.pdf(x_range, mu_cluster, resi_sigma_nm)  # Use re-centered mu for alignment

    # Scale the RESI PDF to align with the Cluster PDF's height for direct comparison
    cluster_max_height = np.max(p_fit_cluster)
    resi_max_height = np.max(p_fit_resi_unscaled)
    scaling_factor = cluster_max_height / resi_max_height

    p_fit_resi_scaled = p_fit_resi_unscaled * scaling_factor

    # --------------------------------------------------------------------
    # 4. Setup Figure and Subplots
    # --------------------------------------------------------------------
    fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(14, 6))
    fig.suptitle("SMLM Cluster Localization vs. RESI Precision", fontsize=14)

    # --------------------------------------------------------------------
    # Subplot 1 (Left): Scatter Plot (x vs. y)
    # --------------------------------------------------------------------
    ax1 = axes[0]

    # Sample a fixed number of points (more robust than frac) or keep original logic
    # Using 'frac' as in the original code:
    sampled_df = final_df.sample(frac=0.20, random_state=1)

    ax1.scatter(
        x=sampled_df["x"],
        y=sampled_df["y"],
        s=1,
        alpha=0.7,
        color='darkblue'
    )

    # Add a circle representing the expected RESI localization precision (2.35 * sigma is FWHM/2)
    # FWHM (Full Width at Half Maximum) is ~2.355 * sigma
    fwhm_radius = 2.355 * resi_sigma_nm / 2  # Using the FWHM approximation for a clearer visual boundary

    circle = Circle(
        (0, 0),
        fwhm_radius,  # Use FWHM/2 or FWHM itself depending on desired visualization
        color='red',
        alpha=0.2,  # Use lower alpha for just a boundary hint
        fill=True
    )
    ax1.add_patch(circle)

    # Add a visual reference for the RESI sigma
    ax1.text(fwhm_radius * 1.05, 0, f'$\sigma_{{RESI}}={resi_sigma_nm:.2f}$ nm',
             color='red', fontsize=9, ha='left', va='center')

    ax1.set_title("Combined Localizations (Center of Mass at [0, 0])")
    ax1.set_xlabel("X Position (nm, Re-centered)")
    ax1.set_ylabel("Y Position (nm, Re-centered)")
    ax1.set_xlim(left=-60, right=60)
    ax1.set_ylim(bottom=-60, top=60)
    ax1.axhline(0, color='r', linestyle=':', alpha=0.5)
    ax1.axvline(0, color='r', linestyle=':', alpha=0.5)
    ax1.set_aspect('equal', adjustable='box')

    # --------------------------------------------------------------------
    # Subplot 2 (Right): Scaled Histogram + Gaussian Fits
    # --------------------------------------------------------------------
    ax2 = axes[1]

    # Plot the histogram (already density=True)
    final_df[column_name].hist(
        ax=ax2,
        bins=80,
        density=True,
        edgecolor='black',
        alpha=0.2,
        color='skyblue',
        label='Cluster Data Histogram'
    )

    # Plot the Cluster Gaussian Fit
    ax2.plot(x_range, p_fit_cluster, 'b-',
             linewidth=2, label=f'Cluster Fit ($\sigma={std_cluster:.2f}$ nm)')

    # Plot the Scaled RESI Gaussian Fit
    ax2.plot(x_range, p_fit_resi_scaled, 'r-',
             linewidth=2, label=f'RESI Precision ($\sigma={resi_sigma_nm:.2f}$ nm)')

    ax2.set_title(f"Distribution of $\Delta${column_name} Localizations")
    ax2.set_xlabel(f"{column_name} Position (Re-centered, nm)")
    ax2.set_ylabel("Probability Density (Scaled)")
    ax2.set_xlim(left=-60, right=60)
    ax2.grid(axis='y', alpha=0.5)
    ax2.legend()

    # Final rendering
    plt.tight_layout()

    # Save the figure
    plt.savefig(os.path.join(path, "DNA-PAINT vs RESI resolution.png"))

    plt.show()


# ---------------------------------------------------------------------
# --- Main Execution Block ---
# ---------------------------------------------------------------------
if __name__ == "__main__":
    RESCALE_FACTOR = 130  # The factor to convert pixel units to nanometers

    try:
        # 1. Load Configuration
        path, resi_file_name, smlm_cluster_name, _ = load_config("config.ini")

        full_resi_path = os.path.join(path, resi_file_name)
        full_cluster_path = os.path.join(path, smlm_cluster_name)

        # 2. Load Data
        resi_df = load_hdf5_data(full_resi_path)
        cluster_df = load_hdf5_data(full_cluster_path)

        # 3. Process RESI Data (Localization Precision)
        # Assuming 'lpx' is the localization precision (sigma) in the RESI data
        resi_sigma_mean = resi_df["lpx"].mean()
        resi_sigma_nm = resi_sigma_mean * RESCALE_FACTOR

        # 4. Plot Results
        # The 'cluster_df' contains the grouping column needed for the function's internal processing
        plot_combined_gaussians_subplots(cluster_df, resi_sigma_nm, "x", RESCALE_FACTOR)

    except FileNotFoundError as e:
        print(f"File Error: {e}")
        print("Please ensure 'config.ini' is correctly configured and the HDF5 files exist.")
        sys.exit(1)  # Exit gracefully on error
    except (ValueError, IOError) as e:
        print(f"Data Error: {e}")
        sys.exit(1)
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
        sys.exit(1)