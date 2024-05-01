#!/usr/bin/env python

import numpy as np
import os as os
import sys
import argparse
from MPASDomainLib import *
from pathlib import Path

_wps_file = 'namelist.wps'

_xL = 900.e3
_nx = 300
_ny = 300

debug = 11

dir = '/scratch/ywang/MPAS/gnu/mpas_scripts/run_dirs/20240410/dacycles.noise4'
in_grid_file = '/scratch/ywang/MPAS/gnu/mpas_scripts/run_dirs/20240410/init/wofs_mpas_01.init.nc'


output_variables = {'w':'w', 'u':'uReconstructZonal', 'v':'uReconstructMeridional', 
                    'qr': 'qr', 'surface_pressue': 'surface_pressure', 'q2': 'q2', 't2m':'t2m', 'dbz': 'refl10cm'}
                    
#=======================================================================================================
#
#
#  MPAS_LQG:  Lou's Quick reGrider for MPAS history files, and DART output files based on the MPAS grid
#
# 
#=======================================================================================================

def MPAS_lqg(in_grid_file, in_data_file, out_filename, nearest=True):

    ds_data = xr.open_dataset(in_data_file)
    ntimes  = ds_data.Time.shape[0]
    ds_data.close()

    ds_grid = xr.open_dataset(in_grid_file)
    try:
        if ds_grid.on_a_sphere == "YES":
            sphere = True
    except:
        sphere = False
    ds_grid.close()

    print(sphere)


    if debug > 100:
        calc_MPAS_grid_stat( in_grid_file )

    # Most of the work in the main routine is setting up a new horizontal grid

    xg, yg, zg, xC, yC, zC = calc_MPAS_new_grid(in_grid_file, wps_file =_wps_file, 
                                                nx = _nx, ny = _ny, xL_grid = _xL, yL_grid = _xL)

    new_grid = calc_MPAS_quad_grid( in_data_file, xC, yC, xg, yg, out_vars = output_variables, nearest = nearest )

    write_MPAS_quad_netCDF( new_grid, xg, yg, zg, ntimes, out_filename, latlon=sphere )

#=======================================================================================================

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    
    parser.add_argument('-d', dest="in_data_file", type=str,
                        help="Input file from MPAS which has data, if same as grid file, not needed",default="")

    parser.add_argument('-o', dest="outfile", type=str, \
                        help="Filename for interpolated out from MPAS on quad grid", default="")

    parser.add_argument('--nearest', dest="nearest", action='store_false', \
                        help="Flag to turn on 5 pt IDW interpolation", default=True)

    args = parser.parse_args()

    if args.in_data_file == "":
        print(" NO DATA FILE ")
        sys.exit(1)
    else:
        ldir = args.in_data_file[0:5]
        basename = os.path.basename(args.in_data_file)
        in_data_file = os.path.join(dir, args.in_data_file)
        path = Path(ldir)
        path.mkdir(exist_ok=True)
        print("\n Made directory %s" % ldir)

    if args.outfile == "":
        out_filename = ("%s_quad.nc") % os.path.join(ldir, basename[0:-3])
    else:
        out_filename = args.outfile

    nearest = True

    for n in np.arange(18):
    
        print( basename[0:-12] )
        new_data_file = "%s_%2.2i.analysis" % (basename[0:-12],n+1)


        in_data_file = os.path.join(dir, ldir, new_data_file)
        print(" Processing:  %s" % in_data_file)

        MPAS_lqg( in_grid_file, in_data_file, out_filename, nearest=nearest)
    
    print("\n Finished MPAS_LQG process")
