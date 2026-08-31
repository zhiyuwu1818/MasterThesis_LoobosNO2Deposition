# A3: Data acquisition for leaf area index (LAI) data at Loobos flux tower, this script is a part of data acquisition process.
# In this script, we will get the LAI measurement data for the period of 2023-2025 from a zip downloaded from ICOS website 
# (https://meta.icos-cp.eu/objects/RBbGtSvPuyjgktNpRDUq5fpw), and then extract the data and save as csv file for later use.
# Author: Zhiyu Wu, Date: 11/03/2026

# ------ Step 1: Import necessary libraries ------
import pandas as pd
import glob
import os
import zipfile
import shutil

# ------ Step 2: Define the folder path, unzip the zip file ------
folder_path = 'Data/VeluweInSitu'
file = os.path.join(folder_path, 'ICOSETC_NL-Loo_ARCHIVE_INTERIM_L2.zip')

# Extract the outer zip file
extracted_folder_path = "Data/VeluweInSitu/NL-Loo_ARCHIVE"
os.makedirs(extracted_folder_path, exist_ok=True)

with zipfile.ZipFile(file, "r") as z:
    z.extractall(extracted_folder_path)

print("Outer zip extracted")

# Extract the inner zip file
print(os.listdir(extracted_folder_path))
inner_zip = os.path.join(extracted_folder_path, "ICOSETC_NL-Loo_ARCHIVE_INTERIM_L2.zip")

with zipfile.ZipFile(inner_zip, "r") as z:
    z.extractall(extracted_folder_path)

print("Inner zip extracted")

# ------ Step 3: Read the CSV file, derive information related to LAI, read it into a DataFrame ------
# Select the CSV file and read it into a DataFrame
csv_file = os.path.join(extracted_folder_path, "ICOSETC_NL-Loo_ANCILLARY_INTERIM_L2.csv")
ancillary_data = pd.read_csv(csv_file, na_values=["-9999"])

# Select LAI group in the file
lai_list = []
db = ancillary_data[ancillary_data["VARIABLE_GROUP"] == "GRP_LAI"]

if db.empty:
    print(f"No LAI data in the file {csv_file}")

# Loop through groups (same as GROUP_ID loop in R)
for gid in db["GROUP_ID"].unique():

    g = db[db["GROUP_ID"] == gid]

    # Helper function
    def get_val(var):
        v = g.loc[g["VARIABLE"] == var, "DATAVALUE"]
        return v.iloc[0] if not v.empty else None

    row = {
        "SITE_ID": g["SITE_ID"].iloc[0],
        "GROUP_ID": gid,
        "LAI_TYPE": get_val("LAI_TYPE"),
        "LAI_CANOPY_TYPE": get_val("LAI_CANOPY_TYPE"),
        "LAI_STATISTIC": get_val("LAI_STATISTIC"),
        "LAI": get_val("LAI"),
        "LAI_DATE": get_val("LAI_DATE"),
        "LAI_DATE_START": get_val("LAI_DATE_START"),
        "LAI_DATE_END": get_val("LAI_DATE_END")
    }

    lai_list.append(row)

    # Combine all rows into a DataFrame
lai_df = pd.DataFrame(lai_list)

# Create unified date column like R script
lai_df["LAI_DATEu"] = (
    lai_df["LAI_DATE"]
    .fillna(lai_df["LAI_DATE_START"])
    .astype(str)
    .str[:8]
)

# Print the resulting DataFrame
print(lai_df)

# ------ Step 4: Save the resulting DataFrame to a new CSV file, and delete the extracted folder to save space ------
# Save the resulting DataFrame to a new CSV file
output_csv = os.path.join(folder_path, "Loobos_LAI_23-25.csv")
lai_df.to_csv(output_csv, index=False)
print(f"LAI data saved to {output_csv}")

# Deleteting the extracted folder to save space
shutil.rmtree(extracted_folder_path)
print("Extracted folder deleted")