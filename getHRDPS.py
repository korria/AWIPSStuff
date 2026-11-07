#!/awips2/python/bin/python
#
#
# getHRDPS.py
#
# Program to download Canadian HRDPS model, clips and prepares for AWIPS ingest.
#
# HISTORY:
# - 2024-05-24 - Initial version - Korri Anderson
#----------------------------------------------------------------
import os
import subprocess
import requests
from requests.exceptions import RequestException
import pygrib
import numpy as np
import xarray as xr
from datetime import datetime, timedelta
import time

# Base URL for HRDPS data

base_url = "https://dd.weather.gc.ca/today/model_hrdps/continental/2.5km/"
variables = [
 "_TMP_AGL-2m",
 #"_DPT_AGL-2m",
 #"_DEPR_AGL-2m",
 #"_RH_AGL-2m", 
 "_HPBL_Sfc",
 "_HGT_Sfc",
 "_HGT_ISBL_0500",
 "_HGT_ISBL_0700",
 "_HGT_ISBL_0750",
 "_HGT_ISBL_0800",
 "_HGT_ISBL_0850",
 "_HGT_ISBL_0875",
 "_HGT_ISBL_0900",
 "_HGT_ISBL_0925",
 #"_HGT_ISBL_0950",
 "_HGT_ISBL_1000",
 "_TMP_ISBL_0500",
 "_TMP_ISBL_0700",
 "_TMP_ISBL_0750",
 "_TMP_ISBL_0800",
 "_TMP_ISBL_0850",
 "_TMP_ISBL_0875",
 "_TMP_ISBL_0900",
 "_TMP_ISBL_0925",
 #"_TMP_ISBL_0950",
 "_TMP_ISBL_1000",
 "_RH_ISBL_0500",
 "_RH_ISBL_0700",
 "_RH_ISBL_0750",
 "_RH_ISBL_0800",
 "_RH_ISBL_0850",
 "_RH_ISBL_0875",
 "_RH_ISBL_0900",
 "_RH_ISBL_0925",
 #"_RH_ISBL_0950",
 "_RH_ISBL_1000",
 #"_HGT_ISBY-1000-500",
 #"_GUST_AGL-10m",
 #"_GUST-Max_AGL-10m",
 #"_VGRD_AGL-10m",
 #"_UGRD_AGL-10m",
 "_WIND_AGL-10m",
 "_WDIR_AGL-10m",
 "_APCP-Accum1h_Sfc",
 "_APCP_Sfc",
 "_TCDC_Sfc",
 "_SNOD_Sfc",
 "_PRMSL_MSL",
 "_PRATE_Sfc",
 "_PTYPE_Sfc",
 #"_WEAFR_Sfc",
 #"_WEAPE_Sfc",
 #"_WEARN_Sfc",
 #"_WEASN_Sfc",
 "-WEonG_DPT_Sfc",
 "-WEonG_GUST_Sfc",
 #"-WEonG_HGTSNLVL_Sfc",
 "-WEonG_CONDAPCPN_Sfc",
 "-WEonG_PCPNTYPE_Sfc",
 "-WEonG_DMNTPCPNTYPE_Sfc",
 "-WEonG_PROBPCPN_Sfc",
 #"-WEonG_PROBLPCPN_Sfc",
 #"-WEonG_PROBFZPCPN_Sfc",
 #"-WEonG_PROBSN_Sfc",
 #"-WEonG_PROBSNSQ_Sfc",
 #"-WEonG_PROBFZRA_Sfc",
 #"-WEonG_PROBRA_Sfc",
 #"-WEonG_PROBTSOCRNC_Sfc",
 #"-WEonG_PROBPL_Sfc",
 #"-WEonG_PROBDZ_Sfc",
 #"-WEonG_PROBFZDZ_Sfc",
 #"_LAND_Sfc",
 
 ]

# Clipping bounds for Idaho, Washington, and Oregon (approximate)
lon_min, lon_max = -125, -110
lat_min, lat_max = 40, 50

# Directory to store downloaded files
download_dir = "/ldad/localapps/getHRDPS/hrdps_data"
os.makedirs(download_dir, exist_ok=True)

# Get the latest date and hour
def get_latest_date_hour():
    current_time = datetime.utcnow()
    yesterday = current_time - timedelta(days=1)
    latest_date = current_time.strftime("%Y%m%d")
    yest_date = yesterday.strftime("%Y%m%d")
    current_hour = current_time.hour
    latest_hour = "00"
    if current_hour >= 22:
       #latest_date = yest_date
       latest_hour = "18"
    elif current_hour >= 16:
       latest_hour = "12"
    elif current_hour >= 10:
       latest_hour = "06"
    elif current_hour >= 4:
       latest_hour = "00"
    elif current_hour < 4:
       latest_hour = "18"
       latest_date = yest_date
    
    return latest_date, latest_hour

latest_date, latest_hour = get_latest_date_hour()

# Function to run a shell command
def run_command(command):
    process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    if process.returncode != 0:
        print(f"Error running command: {command}")
        print(stderr.decode("utf-8"))
    return stdout.decode("utf-8"), stderr.decode("utf-8")

def fetch_with_retries(url, retries=3, backoff_factor=0.3):
    for i in range(retries):
        try:
            response = requests.get(url, timeout=10)
            response.raise_for_status()
            return response
        except (requests.Timeout, requests.ConnectionError) as e:
            if i < retries - 1:
                wait = backoff_factor * (2 ** i)
                print(f"Attempt {i + 1} failed, retrying in {wait} seconds...")
                time.sleep(wait)
            else:
                print("All retry attempts failed.")
                raise
        except RequestException as e:
            print("Request failed: ", e)
            raise

# Function to download, clip, and combine GRIB2 files for each forecast hour
def download_and_combine_files(forecast_hour):
    forecast_str = f"{forecast_hour:03}"
    combined_grib_path = os.path.join(download_dir, f"LDAD-GRIB-CMC-HRDPS_{latest_date}{latest_hour}_P{forecast_str}.grib2")
    print(f"consolidated_hrdps_{latest_date}{latest_hour}_f{forecast_str}.grib2")
    
    with open(combined_grib_path, 'wb') as combined_grib:
        for var in variables:
            file_name = f"{latest_date}T{latest_hour}Z_MSC_HRDPS{var}_RLatLon0.0225_PT{forecast_str}H.grib2"
            file_url = f"{base_url}{latest_hour}/{forecast_hour:03}/{file_name}"
            try:
                response = fetch_with_retries(file_url)
            
                if response.status_code == 200:
                    print(file_name)
                    temp_file_path = os.path.join(download_dir, f"{var}.grib2")
                    with open(temp_file_path, 'wb') as temp_file:
                        temp_file.write(response.content)

                    # Append clipped data to combined file
                    with open(temp_file_path, 'rb') as temp_file:
                        combined_grib.write(temp_file.read())

                    os.remove(temp_file_path)
                else:
                    print(f"Failed to download {file_url}")
            except RequestException as e:
                print("Failed to retrieve URL: ",e)

        # Clip the data using wgrib2
        clipped_grib_path = os.path.join(download_dir, f"LDAD-GRIB-CMC-HRDPS_{latest_date}{latest_hour}_P{forecast_str}.grib2.int")
        clip_command = f"/ldad/localapps/getHRDPS/wgrib2 {combined_grib_path} -set analysis_or_forecast_process_id 50 -GRIB {clipped_grib_path}"
        run_command(clip_command)
            
        # Regrid to polar stereographic projection with grid spacing of 0.025 degrees
        regrid_grib_path = os.path.join(download_dir, f"LDAD-GRIB-BOI-HRDPS_{latest_date}{latest_hour}_P{forecast_str}.grib2")
        regrid_command = (
            f"/ldad/localapps/getHRDPS/wgrib2 {clipped_grib_path} -new_grid_winds earth "
            f"-new_grid nps:-117.0:60.0 -130.0:750:2500 39.0:600:2500 {regrid_grib_path}"
        )
        run_command(regrid_command)
        run_command(f"mv {regrid_grib_path} /data/Incoming/")

        os.remove(clipped_grib_path)
        os.remove(combined_grib_path)


    print(f"Consolidated and clipped GRIB2 file for forecast hour {forecast_hour} saved to {combined_grib_path}")

# Download, clip, and combine GRIB2 files for forecast hours 001 through 048
for hour in range(0, 49):
    download_and_combine_files(hour)

print("Downloaded and consolidated HRDPS data for all forecast hours.")
