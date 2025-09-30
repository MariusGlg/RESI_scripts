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
import pandas as pd
from scipy.spatial import KDTree
import matplotlib.pyplot as plt
import csv
import os
from configparser import ConfigParser

#  config.ini file
config = ConfigParser()
file = "config.ini"
config.read(file)
# config parameter
path = config["INPUT_FILES"]["path"]
NND_file_name = config["INPUT_FILES"]["NND_file_name"]
cluster = int(config["PARAMETERS"]["cluster"])
area_ROI = int(config["PARAMETERS"]["area_ROI"])

# calculate cluster density
density = cluster/area_ROI
# Calculate the required side length of the square plot
side_length = np.sqrt(area_ROI)

# Open CSV and plot data
def open_csv(path):
    with open(path, mode ='r') as file:
        csvFile = csv.reader(file)
        data = list(csvFile)
        # extract NND values from list of lists and store as array
        return data

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

# Histogram
plt.hist(nn_dist1, bins = 100, alpha=0.5, label = '1st NND_sim', weights=np.ones_like(nn_dist1) / len(nn_dist1),
         histtype='step', linewidth=2)
plt.hist(nn_dist2, bins = 100, alpha=0.5, label = '2nd NND_sim', weights=np.ones_like(nn_dist2) / len(nn_dist2),
         histtype='step', linewidth=2)
plt.hist(nn_dist3, bins = 100, alpha=0.5, label = '3rd NND_sim', weights=np.ones_like(nn_dist3) / len(nn_dist3),
         histtype='step', linewidth=2)
plt.title("Experimental vs CSR NND: "'{}'.format(density))
plt.xlabel("nm")
plt.ylabel("rel. frequency")
plt.xlim([0, 500])
# Plot experimenta data
plt.hist(first_NND, bins = 100, alpha=0.5, label = '1st NND_exp',
         weights=np.ones_like(first_NND) / len(first_NND), edgecolor = "black", color = "skyblue")
plt.hist(second_NND, bins = 100, alpha=0.5, label = '1st NND_exp',
         weights=np.ones_like(second_NND) / len(second_NND), edgecolor = "black", color = "navajowhite")
plt.hist(third_NND, bins = 100, alpha=0.5, label = '1st NND_exp',
         weights=np.ones_like(third_NND) / len(third_NND), edgecolor = "black", color = "mediumaquamarine")
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
