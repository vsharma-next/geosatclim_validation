import pandas as pd
import numpy as np
import xarray as xr
from datetime import datetime
import os


def clean_station_name(name):
    """Clean station name by removing spaces and special characters"""
    cleaned = (
        name.replace(" ", "")
        .replace("_", "")
        .replace("/", "")
        .replace("(", "")
        .replace(")", "")
    )
    return cleaned


def get_variable_name(ebcode):
    """Map GEBA codes to standardized variable names"""
    var_mapping = {"GLOBAL": "SIS", "DIFFUS": "SISDIF", "DIRECT": "SISDIR"}
    return var_mapping.get(ebcode, None)


def process_geba_data(metadata_path, data_path, output_path):
    """Single-pass processing of GEBA data with station analysis"""
    metadata = pd.read_csv(metadata_path)
    data = pd.read_csv(data_path)

    # Initialize storage
    station_info = {}
    station_coverage = {}

    # Process metadata and data together
    for _, row in metadata.iterrows():
        station_name = clean_station_name(row["sgname"])
        var_name = get_variable_name(row["ebcode"])

        if var_name is None:
            continue

        # Store station metadata
        if station_name not in station_info:
            station_info[station_name] = {
                "lon": row["sgxlon"],
                "lat": row["sgylat"],
                "elevation": row["sgelevation"],
                "variables": {},
            }

            # Get temporal coverage for this station
            station_data = data[
                data["tskey"].isin(
                    metadata[metadata["sgname"] == row["sgname"]]["tskey"]
                )
            ]
            if not station_data.empty:
                station_coverage[station_name] = {
                    "lon": row["sgxlon"],
                    "lat": row["sgylat"],
                    "elevation": row["sgelevation"],
                    "start_year": station_data["year"].min(),
                    "end_year": station_data["year"].max(),
                }

        station_info[station_name]["variables"][var_name] = row["tskey"]

    # Save station coverage
    coverage_df = pd.DataFrame.from_dict(station_coverage, orient="index")
    coverage_df.index.name = "station"
    coverage_df.to_csv(f"{output_path}/geba_station_coverage.csv")

    # Create netCDF files
    data["timestamp"] = (
        pd.to_datetime(data[["year", "month"]].assign(day=1)).astype("int64") // 1e9
    )
    data_by_tskey = {tskey: group for tskey, group in data.groupby("tskey")}

    # Process each station
    for station_name, station_meta in station_info.items():
        station_data = {}
        time_points = set()

        for var_name, tskey in station_meta["variables"].items():
            if tskey in data_by_tskey:
                var_data = data_by_tskey[tskey]
                station_data[var_name] = var_data
                time_points.update(var_data["timestamp"])

        if not station_data:
            continue

        # Create sorted time array
        time_values = np.array(sorted(time_points))

        # Initialize dataset with coordinates
        ds = xr.Dataset(
            coords={
                "lon": np.array([station_meta["lon"]]),
                "lat": np.array([station_meta["lat"]]),
                "Z": np.array([station_meta["elevation"]]),
                "time": time_values,
            }
        )

        # Add each variable's data
        for var_name, var_data in station_data.items():
            # Create full timeseries with NaN for missing values
            full_data = np.full((len(time_values), 1, 1, 1), np.nan)

            # Create time index mapping
            time_indices = {t: i for i, t in enumerate(time_values)}

            # Fill available data points
            for _, row in var_data.iterrows():
                idx = time_indices[row["timestamp"]]
                full_data[idx, 0, 0, 0] = row["converted_flux_avg"]

            ds[var_name] = (["time", "Z", "lat", "lon"], full_data)
            ds[var_name].attrs["units"] = "W/m²"
            ds[var_name].attrs["_FillValue"] = -32767.0
            ds[var_name].attrs["long_name"] = f"{var_name}; monthly mean"

        # Set coordinate attributes
        ds.lon.attrs["units"] = "degrees_east"
        ds.lon.attrs["long_name"] = "longitude"
        ds.lat.attrs["units"] = "degrees_north"
        ds.lat.attrs["long_name"] = "latitude"
        ds.Z.attrs["units"] = "m"
        ds.Z.attrs["long_name"] = "elevation above sea level"
        ds.time.attrs["units"] = "seconds since 1970-01-01 00:00:00"
        ds.time.attrs["long_name"] = "time"
        ds.time.attrs["calendar"] = "gregorian"

        # Global attributes
        ds.attrs["Title"] = "GEBA Station Data"
        ds.attrs["Institution"] = "GEBA ETH Zurich"
        ds.attrs["Conventions"] = "CF-1.4"
        ds.attrs["Timesteps"] = "Start of compositing period"
        ds.attrs["Coordinates"] = "Pixel center location"

        # Add time coverage attributes
        ds.attrs["time_coverage_start"] = pd.Timestamp(
            time_values[0], unit="s"
        ).isoformat()
        ds.attrs["time_coverage_end"] = pd.Timestamp(
            time_values[-1], unit="s"
        ).isoformat()

        # Save the file
        filename = f"{output_path}/geba/{station_name}.M.nc"
        ds.to_netcdf(filename)

    return station_info, coverage_df


def main():
    metadata_path = "./library/geba_metadata.csv"
    monthly_data_path = "./01_raw_data/geba_monthlydata.zip"
    output_path = "./02_analysis_ready/"

    os.makedirs(f"{output_path}/geba", exist_ok=True)

    station_info, coverage = process_geba_data(
        metadata_path, monthly_data_path, output_path
    )
    print(f"Processed {len(station_info)} stations")
    print("Station coverage summary saved to geba_station_coverage.csv")


if __name__ == "__main__":
    main()
