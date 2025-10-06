"""
@author. Marius Glogger
Optical Imaging Competence Centre, FAU Erlangen

Plot nearest-neighbor analysis data (NND) from SMLM RESI experiments (generated with Picasso @Jungmann group)
and compare to complete spatial randomness (CSR) at given cluster densities.
- load .csv files containing NND data
- calculates the cluster density from .yaml file in folder
- performs NND (via KDTree) search on completely spatial randomized dataset that uses density as in
    experimental settings
- plots experimental NND and CSR NND
- Saves plot as .png and CSR NND as .csv in folder defined by foldername

input
folder: path to folder containing NND .csv file
NND_file: name of NND-file (.csv)
"""


import numpy as np
from scipy.spatial import KDTree
import matplotlib.pyplot as plt
import csv
import os
from configparser import ConfigParser
import seaborn as sns
import yaml


# --- Configuration Management ---
def load_config(file_path="config.ini"):
    """Loads configuration from the specified INI file."""
    config = ConfigParser()
    if not config.read(file_path):
        raise FileNotFoundError(f"Configuration file '{file_path}' not found.")

    path = config["INPUT_FILES"].get("path", "")
    resi_file_name = config["INPUT_FILES"].get("RESI_file_name", "RESI.hdf5")
    smlm_cluster_name = config["INPUT_FILES"].get("SMLM_cluster_name", "Cluster.hdf5")
    resi_yaml_name = config["INPUT_FILES"].get("RESI_yaml_name", "RESI.hdf5")
    NND_file_name = config["INPUT_FILES"].get("NND_file_name", "clustered_centers_nn.csv")


    return path, resi_file_name, smlm_cluster_name, config, resi_yaml_name, NND_file_name

path, resi_file_name, smlm_cluster_name, config, resi_yaml_name, NND_file_name = load_config("config.ini")


# Open CSV and plot data
def open_csv(path):
    with open(path, mode ='r') as file:
        csvFile = csv.reader(file)
        data = list(csvFile)
        # extract NND values from list of lists and store as array
        return data

def calc_density(path, target_keys):
    """
    Reads a multi-document YAML file, iterates through all documents,
    and returns a dictionary of all found target keys and their values.
    """
    found_data = {}

    try:
        with open(path, 'r') as file:
            # FIX: Use safe_load_all to correctly handle multiple documents ('---')
            data_documents = list(yaml.safe_load_all(file))

            # Iterate through each document dictionary
            for doc in data_documents:
                if isinstance(doc, dict):
                    # Check each target key against the current document
                    for key in target_keys:
                        if key in doc and key not in found_data:
                            # Store the key/value pair and mark it as found
                            found_data[key] = doc[key]

            return found_data

    except FileNotFoundError:
        print(f"Error: File not found at {os.path.abspath(path)}")
        return None
    except Exception as e:
        print(f"An error occurred while loading or searching the YAML: {e}")
        return None

    except FileNotFoundError:
        print(f"Error: File not found at {os.path.abspath(path)}")
        return None
    except Exception as e:
        print(f"An error occurred while loading or searching the YAML: {e}")
        return None


target_key_names = ["Number of clusters", "Pick Areas (um^2)"]
resi_file = os.path.join(path, resi_yaml_name)
vals = calc_density(resi_file, target_key_names)
cluster = int(vals["Number of clusters"])
area = int(vals["Pick Areas (um^2)"][0])
density = cluster / area
side_length = np.sqrt(area)

NND_file = open_csv(os.path.join(path, NND_file_name))

first_NND = [float(sublist[0]) for sublist in NND_file]
second_NND = [float(sublist[1]) for sublist in NND_file]
third_NND = [float(sublist[2]) for sublist in NND_file]

# Generate random coordinates for n_points  within the calculated area
x = np.random.uniform(0, side_length, cluster)
y = np.random.uniform(0, side_length, cluster)
CSR_coordinates = np.stack((x, y), axis=1)  # stack arrays

tree = KDTree(CSR_coordinates)  # use KDTree for NND search
# Then you find the closest two as the first is the point itself
dists = tree.query(CSR_coordinates, 4)
nn_dist1 = dists[0][:, 1] * 1000  # first NND convert to nm
nn_dist2 = dists[0][:, 2] * 1000  # 2nd NND
nn_dist3 = dists[0][:, 3] * 1000  # 3rd NND

# --- Outlier Removal using IQR ---
def clean_outliers_iqr(data):
    """
    Converts data to a NumPy array (if necessary) and filters out outliers
    using the 1.5 * IQR rule.
    """
    # Ensure data is a NumPy array for proper filtering
    data_array = np.array(data)

    # 1. Calculate Q1, Q3, and IQR
    Q1 = np.percentile(data_array, 25)
    Q3 = np.percentile(data_array, 75)
    IQR = Q3 - Q1

    # 2. Define the bounds
    lower_bound = Q1 - (1.5 * IQR)
    upper_bound = Q3 + (1.5 * IQR)

    # 3. Filter and return the cleaned data
    cleaned_data = data_array[
        (data_array >= lower_bound) & (data_array <= upper_bound)
    ]
    return cleaned_data, lower_bound, upper_bound

# --- APPLY THE FUNCTION TO YOUR DATASETS ---

# This replaces all the repetitive code from your original block.
first_NND_cleaned, lower_bound_1, upper_bound_1 = clean_outliers_iqr(first_NND)
second_NND_cleaned, lower_bound_2, upper_bound_2 = clean_outliers_iqr(second_NND)
third_NND_cleaned, lower_bound_3, upper_bound_3 = clean_outliers_iqr(third_NND)

# Histogram
# Plot simulation
sns.kdeplot(
    data=nn_dist1,
    color="skyblue",
    linewidth=3,
    fill=False,       # Optional: Fills the area under the line
    alpha=1,
    label="sim_NND1_kd")        # Transparency for the fill
sns.kdeplot(
    data=nn_dist2,
    color="navajowhite",
    linewidth=3,
    fill=False,       # Optional: Fills the area under the line
    alpha=1,
    label="sim_NND2_kd")        # Transparency for the fill
sns.kdeplot(
    data=nn_dist3,
    color="mediumaquamarine",
    linewidth=3,
    fill=False,       # Optional: Fills the area under the line
    alpha=1,
    label="sim_NND3_kd")         # Transparency for the fill

plt.title("Experimental vs CSR NND: "'{}'.format(density))
plt.xlabel("nm")
plt.ylabel("rel. frequency")
plt.xlim([0, 300])

# Plot experimental data
plt.hist(first_NND_cleaned, bins='auto', alpha=0.5, label = '1st NND_exp',
         density=True, edgecolor = "black", color = "skyblue")
plt.hist(second_NND_cleaned, bins='auto', alpha=0.5, label = '1st NND_exp',
         density=True, edgecolor = "black", color = "navajowhite")
plt.hist(third_NND_cleaned, bins='auto', alpha=0.5, label = '1st NND_exp',
         density=True, edgecolor = "black", color = "mediumaquamarine")

plt.legend()
plt.savefig(os.path.join(path, "compare_CSR_NND.png"))
plt.show()

# 2. Combine the arrays into a single 2D structure for writing
#    The `np.stack` function combines the arrays as columns
data_to_write = np.stack((nn_dist1, nn_dist2, nn_dist3), axis=1)
# 3. Define the CSV file name and column headers
csv_file = os.path.join(path, "compare_CSR_NND.csv")  # data storage path
headers = ['1st NND', '2nd NND', '3rd NND']  # headers
# 4. Write the data to the CSV file
with open(csv_file, 'w', newline='') as f:
    writer = csv.writer(f)
    # Write the header row
    writer.writerow(headers)
    # Write the data rows
    writer.writerows(data_to_write)
