# MPAS Tool: Lou's Quick Gridder
##------------------------------

mpas_lqg.py is a command line tool to quickly interpolate 2D/3D fields associated with the MPAS grid 
onto a quadralateral grid and writes out the fields to a new netCDF4 file.

It is a lightweight interpolator written in Python. A more complete systemn is in 
Larissa Reames's MPASSIT solution.  See https://github.com/LarissaReames-NOAA/MPASSIT

MPASSIT converts an MPAS history file to appear like a WRF history file, so that standard tools can be used (like UPP).


**USAGE:** mpas_lqg.py [-h] [-i IN_GRID_FILE] [-d IN_DATA_FILE] [-o OUTFILE] [--nearest]

options:

* -h, --help       show this help message and exit
* -i IN_GRID_FILE  Input file from MPAS which has grid information
* -d IN_DATA_FILE  Input file from MPAS which has data, if same as grid file, not needed
* -o OUTFILE       Filename for interpolated out from MPAS on quad grid
* --nearest        Flag to turn on 5 pt IDW interpolation

If you are converting an MPAS init file, or a history file, then only the IN_GRID_FILE needs to be specified.

If you are converting a DART diagnostic file, it has no coordinates, so the IN_GRID_FILE provides the coordinate,
but the IN_DATA_FILE (aka the DART diagnostic file) is the file that is converted.

If you dont supply an output file, the suffix "_quad" will be added to the end of the filename.

Default is nearest neighbor interpolation - its really fast, and probably all anyone needs.

Right now, the list of variables that are converted is limited by those at the top of the mpas_lqg.py file.  

output_variables = {'w':'w', 'u':'uReconstructZonal', 'v':'uReconstructMeridional', 'theta': 'theta',
                    'surface_pressue': 'surface_pressure', 'q2': 'q2', 't2m':'t2m'}
                    
Suggestions welcome.
