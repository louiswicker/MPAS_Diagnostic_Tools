#!/usr/bin/env python3
"""
mpas_minmax.py
--------------
Reads an MPAS output NetCDF file and prints the max/min value
and its index location for every variable.

Usage:
    python mpas_minmax.py <mpas_output_file.nc>
"""

import sys
import numpy as np

try:
    from netCDF4 import Dataset
except ImportError:
    sys.exit("ERROR: netCDF4 is not installed. Run: pip install netCDF4")


def format_location(flat_index, shape):
    """Convert a flat index back to a multi-dimensional index tuple."""
    return np.unravel_index(flat_index, shape)


def print_var_stats(nc_file):
    print(f"\nFile: {nc_file}")
    print("=" * 80)

    with Dataset(nc_file, "r") as ds:
        # Print global dimensions for context
        print("\nDimensions:")
        for dim_name, dim in ds.dimensions.items():
            size = "UNLIMITED" if dim.isunlimited() else str(len(dim))
            print(f"  {dim_name}: {size}")

        print("\n" + "=" * 80)
        print(f"{'Variable':<30} {'Dimensions':<35} {'Min':>15} {'Max':>15}")
        print(f"{'':30} {'':35} {'Location':>15} {'Location':>15}")
        print("-" * 80)

        for var_name, var in ds.variables.items():
            try:
                data = var[:]

                # Skip non-numeric or scalar variables
                if not np.issubdtype(data.dtype, np.number):
                    print(f"  {var_name:<28} (skipped — non-numeric, dtype={data.dtype})")
                    continue

                if data.size == 0:
                    print(f"  {var_name:<28} (skipped — empty array)")
                    continue

                # Mask fill/missing values if present
                if hasattr(data, 'mask'):
                    if data.mask.all():
                        print(f"  {var_name:<28} (skipped — all values masked)")
                        continue
                    data_valid = data.compressed()
                    data_flat  = np.ma.filled(data, fill_value=np.nan).ravel()
                else:
                    data_valid = data.ravel()
                    data_flat  = data.ravel()

                min_val = np.nanmin(data_flat)
                max_val = np.nanmax(data_flat)
                min_idx = int(np.nanargmin(data_flat))
                max_idx = int(np.nanargmax(data_flat))

                min_loc = format_location(min_idx, data.shape)
                max_loc = format_location(max_idx, data.shape)

                dim_str = str(var.dimensions)

                print(f"\n  Variable : {var_name}")
                print(f"  Dims     : {dim_str}")
                if hasattr(var, 'units'):
                    print(f"  Units    : {var.units}")
                if hasattr(var, 'long_name'):
                    print(f"  Long name: {var.long_name}")
                print(f"  Shape    : {data.shape}")
                print(f"  Min      : {min_val:.6g}  @ index {min_loc}")
                print(f"  Max      : {max_val:.6g}  @ index {max_loc}")

            except Exception as e:
                print(f"  {var_name:<28} (ERROR: {e})")

    print("\n" + "=" * 80)
    print("Done.")


def main():
    if len(sys.argv) < 2:
        sys.exit(f"Usage: python {sys.argv[0]} <mpas_output_file.nc>")

    nc_file = sys.argv[1]
    print_var_stats(nc_file)


if __name__ == "__main__":
    main()
