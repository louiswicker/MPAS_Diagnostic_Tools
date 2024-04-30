#!/usr/bin/env python

import numpy as np
import os as os
import sys
import argparse
from MPASDomainLib import *


#=======================================================================================================
#
#
#  MPAS_GRID_STAT:  Quick stats for MPAS grid shapes
#
# 
#=======================================================================================================

def MPAS_grid_stat(in_grid_file):

    calc_MPAS_grid_stat( in_grid_file )

#=======================================================================================================

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    
    parser.add_argument('-i', dest="in_grid_file", type=str,
                        help="Input file from MPAS which has grid information", default="")
                        
    args = parser.parse_args()

    if args.in_grid_file == "":
        print("\n MPAS_GRID_STAT:  You must specify an MPAS file with grid information!!!")
        parser.print_help()
        sys.exit(1)   
    else:
        in_grid_file = args.in_grid_file

    MPAS_grid_stat( in_grid_file )
    
    print("\n Finished MPAS_GRID_STAT\n")
