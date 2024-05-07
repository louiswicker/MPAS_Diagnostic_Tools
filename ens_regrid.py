#!/usr/bin/env python

import numpy as np
import os as os
import sys
import argparse
from datetime import time, datetime, timedelta
from pathlib import Path

from MPASDomainLib import * 

#=======================================================================================================
#
#
#  MPAS_REGRID_BATCH:  driver script for batch processing of MPAS regridding and can also
#                      create difference files if needed
#
# 
#=======================================================================================================

_members = 1 + np.arange(3)   # generate any list that you want...


src_dir_top = '/scratch/ywang/MPAS/gnu/mpas_scripts/run_dirs'
src_day     = '20240410'

_in_grid_file = '/scratch/ywang/MPAS/gnu/mpas_scripts/run_dirs/20240410/init/wofs_mpas_01.init.nc'

_mpas_file_prefix = "wofs_mpas_%2.2i.restart.2024-04-10_%s.%s.00.nc"

_ncdiff_cmd = "ncdiff -O %s %s %s >/dev/null 2>&1"

#=======================================================================================================
# 
# function to make a complicated path

def parse_mpas_file(dir, member, hour, min):

    if dir == "":

        return _mpas_file_prefix % (member, hour, min)

    else:

        return os.path.join(dir, _mpas_file_prefix % (member, hour, min))

#=======================================================================================================

if __name__ == "__main__":


    parser = argparse.ArgumentParser()
    
    parser.add_argument('-d', dest="in_dir", type=str,
                        help="Input file from MPAS which has data, if same as grid file, not needed",default="")

    parser.add_argument('-t', dest="hhmm", type=str, \
                        help="four digit HHMM to process", default=None)

    parser.add_argument('-r', dest="run_it", action='store_true', default=False, 
                        help="Flag to run script, else, its a test run to check things")

    parser.add_argument('--diff', dest="diff", action='store_true', default=False, 
                        help="Flag to run the ncdiff command on prior and posterior quad files")

    parser.add_argument('--interp', dest="interp", action='store_true', \
                        help="Flag to turn on 5 pt IDW interpolation", default=False)

    parser.add_argument('--config', dest="config", type=str, \
                        help="YAML configuration file to read, default is config.yaml", default="config.yaml")

    args = parser.parse_args()

    if args.in_dir == "":
        print(" NO EXPERIENT DIRECTORY SPECIFIED...Exiting!!! ")
        print(" TOP level source directory is: %s" % src_dir_top)
        sys.exit(1)
    else:
        exp_dir = args.in_dir

        print("\n ----------------------------------")
        print("\n Lets make sure you have the right paths....\n")

        check0 = os.path.isdir(src_dir_top)
        check1 = os.path.isdir(os.path.join(src_dir_top,src_day))
        check2 = os.path.isdir(os.path.join(src_dir_top,src_day,exp_dir))

        if check0:
            print("\n Path to top level src directory is okay %s \n " % src_dir_top)
        else:
            print("\n TOP level source directory is NOT found, exiting... %s" % src_dir_top)
            sys.exit(1)

        if check1:
            print("\n Path to day directory is okay %s \n " % os.path.join(src_dir_top,src_day))
        else:
            print("\n Day directory is NOT found, exiting... %s" % os.path.join(src_dir_top,src_day))
            sys.exit(1)

        if check2:
            print("\n Path to experiment directory is okay %s \n" % os.path.join(src_dir_top,src_day,exp_dir))
        else:
            print("\n Experiment directory is NOT found, exiting... %s" % os.path.join(src_dir_top,src_day,exp_dir))
            sys.exit(1)

        master_path =  os.path.join(src_dir_top,src_day,exp_dir) 

    if args.hhmm == None:
        print("\n MPAS_REGRID_BATCH:  You must specify a 4-digit HHMM!")
        parser.print_help()
        sys.exit(1)
    else:
        hhmm = datetime.strptime("%s %s" %(src_day, args.hhmm), '%Y%m%d %H%M')
        print(" Analysis time: %s \n" % hhmm.strftime("%H%M") )
        print(" ----------------------------------\n")

    if args.run_it:
        print("\n Running script!....")
        run_it = args.run_it
    else:
        print("\n ==============>>>>> HEY!!  -r flag missing, so only doing a test run....")
        run_it = args.run_it
        
    if args.diff:
        print("\n Running NCDIFF....")
        nc_dif = args.diff
    else:
        nc_dif = args.diff

    if args.interp == True:
        print("\n MPAS_REGRID:  5-pt IDW interpolation used (1/dis^2)")
        interp = args.interp
    else:
        print("\n MPAS_REGRID:  Nearest neighbor interpolation used ")
        interp = args.interp

#=======================================================================================================

# Create local prior and posterior directories

    local_prior = os.path.join(src_day, exp_dir, hhmm.strftime("%H%M"), 'prior')

    path = Path(local_prior)
    path.mkdir(parents=True, exist_ok=True)

    print("\n Made local prior directory %s" % local_prior)

    local_post = os.path.join(src_day, exp_dir, hhmm.strftime("%H%M"), 'posterior')

    path = Path(local_post)
    path.mkdir(parents=True, exist_ok=True)

    print("\n Made local posterior directory %s" % local_post)

#=======================================================================================================
#
# Processing section

# Loop over ensemble members

    for n in _members:
    
    # PRIORS FIRST

        # Priors in YW structure are stored in main directory T-15        

        hhmm_m15 = hhmm - timedelta(hours=0, minutes=15)
        hhmm_str = hhmm_m15.strftime("%H%M")

        hh_str   = "%2.2i" % hhmm.hour
        mm_str   = "%2.2i" % hhmm.minute

        prior_filename = parse_mpas_file(os.path.join(master_path, hhmm_str), n, hh_str, mm_str)

        print("\n PRIOR member name:  %s " % prior_filename )

        file_exists =  os.path.isfile(prior_filename)

        if file_exists:

            quad_filename  = os.path.basename( ("%s_quad.nc") % prior_filename[0:-3] )

            prior_quad_filename = os.path.join(local_prior, quad_filename)

            print("\n Now computing %s " % prior_quad_filename)

            if run_it:
                MPAS_lqg( _in_grid_file,  prior_filename, prior_quad_filename, ConfigFile=args.config, interp=interp )

        else:
            print("\n ERROR:  PRIOR member does not exist...Exiting script\n")

    # POSTERIORS

        # Posteriors in YW structure are stored in HHMM/fcst_XX directory...

        hhmm_str = hhmm.strftime("%H%M")
        hh_str   = "%2.2i" % hhmm.hour
        mm_str   = "%2.2i" % hhmm.minute

        fcst_dir = "fcst_%2.2i" % n

        post_filename = parse_mpas_file(os.path.join(master_path, hhmm_str, fcst_dir), n, hh_str, mm_str)

        print("\n POSTERIOR member name:  %s " % post_filename )

        file_exists =  os.path.isfile(post_filename)

        if file_exists:

            quad_filename  = os.path.basename( ("%s_quad.nc") % post_filename[0:-3] )

            post_quad_filename = os.path.join(local_post, quad_filename)

            print("\n Now computing %s " % post_quad_filename)

            if run_it:
                MPAS_lqg( _in_grid_file,  post_filename, post_quad_filename, ConfigFile=args.config, interp=interp )

        else:
            print("\n ERROR:  POSTERIOR member does not exist...Exiting script\n")

        if nc_dif:

            diff_filename = ("anal_incre_%2.2i.nc") % n
            diff_path = os.path.join(src_day, exp_dir, hhmm.strftime("%H%M"), diff_filename)

            print("\n --> NCDIFF is creating: %s\n" % diff_path)

            cmd = _ncdiff_cmd % (post_quad_filename, prior_quad_filename, diff_path)
            print("\n  Running ncdiff command: %s" % cmd)

            os.system(cmd)
            

    print("\n Finished MPAS_REGRID_BATCH process")
