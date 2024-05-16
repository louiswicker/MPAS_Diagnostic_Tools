#!/usr/bin/env python

import numpy as np
import os as os
import sys
import argparse
from MPASDomainLib import *
import yaml
from cbook10 import list2dict, read_yaml

debug = 11

dir = '/scratch/ywang/MPAS/gnu/mpas_scripts/run_dirs/20240410/dacycles.one'


#=======================================================================================================

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    
    parser.add_argument('-i', dest="in_grid_file", type=str,
                        help="Input file from MPAS which has grid information", default="")
                        
    parser.add_argument('-d', dest="in_data_file", type=str,
                        help="Input file from MPAS which has data, if same as grid file, not needed",default="")
                        
    parser.add_argument('-o', dest="outfile", type=str, \
                        help="Filename for interpolated out from MPAS on quad grid", default="")

    parser.add_argument('--interp', dest="interp", action='store_true', \
                        help="Flag to turn on 5 pt IDW interpolation", default=False)

    parser.add_argument('--config', dest="config", type=str, \
                        help="YAML configuration file to read, default is config.yaml", default="config.yaml")

    args = parser.parse_args()

    if args.in_grid_file == "":
        print("\n MPAS_REGRID:  You must specify an MPAS file with grid information!!!")
        parser.print_help()
        sys.exit(1)   
    else:
        in_grid_file = args.in_grid_file

    if args.in_data_file == "":
        in_data_file = in_grid_file
    else:
        in_data_file = args.in_data_file
    
    if args.outfile == "":  
        out_filename = ("%s_quad.nc") % in_data_file[0:-3]
    else:
        out_filename = args.outfile     

    if args.interp == True:
        print("\n MPAS_REGRID:  5-pt IDW interpolation used (1/dis^2)")
        interp = args.interp
    else:
        print("\n MPAS_REGRID:  Nearest neighbor interpolation used ")
        interp = args.interp

    MPAS_lqg( in_grid_file, in_data_file, out_filename, ConfigFile=args.config, interp=interp )
    
    print("\n Finished MPAS_LQG process")
