# MPAS Tool: REGRID
================================

**MPAS_REGRID.py** is a command line tool to quickly interpolate 2D/3D fields associated with the MPAS grid 
onto a quadralateral grid and writes out the fields to a new netCDF4 file.

MPAS_REGRID is a lightweight interpolator written in Python. A more complete converter is in 
Larissa Reames's MPASSIT solution.  See https://github.com/LarissaReames-NOAA/MPASSIT

* **MPASSIT** converts an MPAS history file to appear like a WRF history file, so that standard tools can be used (like UPP).

* **MPAS_REGRID** was originally built to regrid MPAS planer output from idealized runs. The ability to regrid from a real data MPAS run was added to examine output files from DART which have no geographic information.


**USAGE:** mpas_regrid.py [-h] [-i IN_GRID_FILE] [-d IN_DATA_FILE] [-o OUTFILE] [--interp]

options:

* -h, --help       show this help message and exit
* -i IN_GRID_FILE  Input file from MPAS which has grid information
* -d IN_DATA_FILE  Input file from MPAS which has data, if same as grid file, not needed
* -o OUTFILE       Filename for interpolated out from MPAS on quad grid
* --interp         Flag to turn on 5 pt IDW interpolation

If you are converting an MPAS init file, or a history file, then only the IN_GRID_FILE needs to be specified.

If you are converting a DART diagnostic file, it has no coordinates, so the IN_GRID_FILE provides the coordinate,
but the IN_DATA_FILE (aka the DART diagnostic file) is the file that is converted.

If you dont supply an output file, the suffix "_quad" will be added to the end of the filename.

Default is nearest neighbor interpolation - its really fast, and probably all anyone needs.

Right now, the list of variables that are converted is limited by those at the top of the MPASDomain.py file.  

output_variables = {'w':'w', 'u':'uReconstructZonal', 'v':'uReconstructMeridional', 'theta': 'theta',
                    'surface_pressue': 'surface_pressure', 'q2': 'q2', 't2m':'t2m'}

# MPAS Tool: Grid Stat
================================

**MPAS_GRID_STAT.py** is a tools that prints statistics out from your mpas grid 


**USAGE:** mpas_grid_stat.py [-h] [-i IN_GRID_FILE] 

options:

* -h, --help       show this help message and exit
* -i IN_GRID_FILE  Input file from MPAS which has grid information


==================

Python Libs Needed
==================

* argparse
* numpy
* matplotlib
* scipy.spatial
* xarray
* cartopy 
* pyproj
                    
Suggestions welcome.
