import numpy as np
import pandas as pd
from netCDF4 import Dataset
from wrf import getvar, to_np, ALL_TIMES
from sklearn.cluster import KMeans
import csv

def extract_initial_cape_shear(ncfile):
    try:
        # CAPE: Convective Available Potential Energy
        cape = to_np(getvar(ncfile, "cape_3d", timeidx=0))
        # Get the max CAPE at each horizontal location (vertical dim is first)
        max_cape_2d = np.nanmax(cape, axis=0)
        max_cape = np.nanmax(max_cape_2d)

        # Wind shear between 1 km and 6 km
        shear = to_np(getvar(ncfile, "wind_shear", timeidx=0, units="m s-1", height=(1000, 6000)))
        shear_magnitude = np.sqrt(shear[0]**2 + shear[1]**2)  # shear[0]=u_diff, shear[1]=v_diff
        max_shear = np.nanmax(shear_magnitude)
        
        return max_cape, max_shear
    except Exception as e:
        print(f"Error extracting CAPE/shear: {e}")
        return np.nan, np.nan
    
def classify_experiment_based_on_max_and_sum(ncfile, variable_name):
    # Extract 4D data (time, lat, lon) using wrf functions
    data = to_np(getvar(ncfile, variable_name, timeidx=ALL_TIMES, method="cat"))

    # Get the maximum value per pixel (lat, lon) over time
    max_per_pixel = np.nanmax(data, axis=0)

    # Calculate the total sum of the 2D field (lat, lon) considering only values > 0
    total_sum = np.nansum(max_per_pixel[max_per_pixel >= 0])

    # Calculate the maximum value in the 2D field (lat, lon)
    max_value = np.nanmax(max_per_pixel)

    # Calculate the ratio of max_value to total_sum
    ratio = max_value / total_sum if total_sum != 0 else np.nan

    # Calculate the mean value for the valid pixels (> 0) in the max_per_pixel field
    mean_value = np.nanmean(max_per_pixel[max_per_pixel >= 0])
    
    return total_sum, max_value, ratio, mean_value


def process_experiment_with_kmeans(csv_file_path, num_clusters=3):
    summary_data = []

    with open(csv_file_path, newline='') as csvfile:
        csv_reader = csv.reader(csvfile)
        next(csv_reader)  # Skip the header row if it exists
        for row_number, row in enumerate(csv_reader):
            conf_num = '{:03d}'.format(row_number)
            input_path = '/home/jorge.gacitua/datosmunin2/EXPERIMENTS_UNWEATHER/DATA/TEST_Multi_4_vars_996/' + conf_num + "/"
            file_in = input_path + 'wrfout_' + conf_num + '.nc'
            
            # Read x1 to x4 values
            try:
                x1, x2, x3, x4 = map(float, row)
            except ValueError:
                print(f"Skipping malformed row {row_number}: {row}")
                continue

            try:
                # Load the WRF file
                ncfile = Dataset(file_in)

                # Calculate reflectivity-based metrics
                total_sum, max_value, ratio, mean_value = classify_experiment_based_on_max_and_sum(ncfile, 'mdbz')

                # Extract max CAPE and shear at the initial time
                max_cape, max_shear = extract_initial_cape_shear(ncfile)

                # Store all values: conf_num, wrf-based metrics, x1–x4
                summary_data.append([
                    conf_num, total_sum, max_value, ratio, mean_value, max_cape, max_shear, x1, x2, x3, x4
                ])
                print(f"Processed experiment {conf_num}")
            except Exception as e:
                print(f"Error processing experiment {conf_num}: {e}")
                continue


    # Create a DataFrame with the summary of all experiments
    columns = ['Experiment', 'Total Sum', 'Max Value', 'Ratio', 'Mean', 'Max CAPE', 'Max Shear', 'x1', 'x2', 'x3', 'x4']
    df_summary = pd.DataFrame(summary_data, columns=columns)

    # Optional: use only meteorological variables, or include x1–x4 too
    features_for_clustering = ['Total Sum', 'Max Value', 'Mean', 'Max CAPE', 'Max Shear']  # or add 'x1'...'x4'

    kmeans = KMeans(n_clusters=num_clusters, random_state=0)
    df_summary['Cluster'] = kmeans.fit_predict(df_summary[features_for_clustering])

    # Save the summary with cluster labels to a CSV file
    df_summary.to_csv('experiment_summary_with_clusters.csv', index=False)
    print("Experiment summary with clusters saved to 'experiment_summary_with_clusters.csv'")


if __name__ == '__main__':
    csv_file_path = 'sampling_batch_996.csv'
    process_experiment_with_kmeans(csv_file_path, num_clusters=4)  # Adjust num_clusters as needed
