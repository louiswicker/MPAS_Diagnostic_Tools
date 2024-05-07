#!/usr/bin/env python

import numpy as np
import os as os
import sys
import argparse
from MPASDomainLib import *

_wps_file = 'namelist.wps'

_xL = 900.e3
_nx = 300
_ny = 300

debug = 11

dir = '/scratch/ywang/MPAS/gnu/mpas_scripts/run_dirs/20240410/dacycles.noise4'


varlist = ['w', 'uReconstructZonal', 'uReconstructMeridional', 'refl10cm', 'surface_pressure', 'q2', 't2m',
           'theta', 'rho', 'qv', 'qc', 'qr', 'qi', 'qs', 'qg', 'qh', 
           'volg', 'volh', 'nc', 'nr', 'ni', 'ns', 'ng', 'nh']
                    
#=======================================================================================================
#
def list2dict(list1, list2):

    return dict(zip(list1, list2))

output_variables = list2dict(varlist, varlist)

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

    args = parser.parse_args()

    if args.in_grid_file == "":
        print("\n MPAS_LQG:  You must specify an MPAS file with grid information!!!")
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

    interp = args.interp
    
    MPAS_lqg( in_grid_file, in_data_file, out_filename, interp=interp )
    
    print("\n Finished MPAS_LQG process")
