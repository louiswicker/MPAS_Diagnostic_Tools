#====================================================================================================
#
#
# This code was initially copied from the GitHub directory: https://github.com/lucas-uw [Xiaodong Chen] from
# Pacific Northest Laboratory. It was written for WRF, and helped inspire me to try and build a similar
# set of tools for MPAS.
#
# Lou Wicker started adapting it to MPAS in April of 2024
#
# Many thanks to lucas-uw!
#
#====================================================================================================

import numpy as np
import cartopy.crs as ccrs
from pyproj import Transformer
import xarray as xr
import re  # regular expression
import os.path
from scipy.spatial import KDTree
from cbook10 import list2dict, read_yaml

debug = 0

_scale = 0.001

_interp = False

default_params = {
 'map_proj'  : 'lambert',
 'ref_lat'   :  35.00,
 'ref_lon'   : -90.00,
 'truelat1'  :  30.0,
 'truelat2'  :  60.0,
 'stand_lon' : -90.0,
                 }

#=======================================================================================================
#
#
#  MPAS_LQG:  Lou's Quick reGrider for MPAS history files, and DART output files based on the MPAS grid
#
#
#=======================================================================================================

def MPAS_lqg( in_grid_file, in_data_file, out_filename, ConfigFile='config.yaml', interp=False ):

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

    if debug > 100:
        calc_MPAS_grid_stat( in_grid_file )

    # Read input configuration file

    config = read_yaml(ConfigFile)

    output_variables = list2dict(config['VarList'], config['VarList'])

    _wps_file, _xL, _nx, _ny = config['QuadGrid']

    # Most of the work in the main routine is setting up a new horizontal grid

    xg, yg, zg, xC, yC, zC = calc_MPAS_new_grid(in_grid_file, wps_file =_wps_file,
                                                nx = _nx, ny = _ny, xL_grid = _xL, yL_grid = _xL)

    new_grid = calc_MPAS_quad_grid( in_data_file, xC, yC, xg, yg, out_vars = output_variables, interp=interp )

    write_MPAS_quad_netCDF( new_grid, xg, yg, zg, ntimes, out_filename, latlon=sphere )

#====================================================================================================
def get_wps_param_value(wps_file, param_name, vartype):

    def read_default():
        try:
            output = default_params[param_name]
            return output
        except KeyError:
            print("\n GET_WPS_PARAM: Invalid parameter name not in defaults, exiting...\n")
            return None
    try:
        with open(wps_file, 'r') as file:
            for line in file.readlines():
                words = re.split('=|\s+|,|\'', line)
                while '' in words:
                    words.remove('')
                if param_name in words:
                    if vartype=='float':
                        output = float(words[1])
                    elif vartype=='int':
                        output = int(words[1])
                    else:
                        output = words[1]

                    return output

                else:
                    if debug > 0:  print("\n GET_WPS_PARAM: No parameter found in file, using defaults...\n")
                    return read_default()

    except FileNotFoundError:
        if debug > 0:  print("\n GET_WPS_PARAM: No file found, using defaults...\n")
        return( read_default() )

#====================================================================================================

def get_transformer(wps_file, latlon2xy = True):

    ref_lat = get_wps_param_value(wps_file, 'ref_lat', 'float')
    ref_lon = get_wps_param_value(wps_file, 'ref_lon', 'float')
    par_lat1 = get_wps_param_value(wps_file, 'truelat1', 'float')
    par_lat2 = get_wps_param_value(wps_file, 'truelat2', 'float')

    proj_daymet = '+proj=lcc +lon_0=%f +lat_1=%f +lat_2=%f' % (ref_lon, par_lat1 , par_lat2)
    
    if latlon2xy:
        return Transformer.from_crs("EPSG:4326", proj_daymet, always_xy=True)
        
    else:
        return Transformer.from_crs(proj_daymet, "EPSG:4326", always_xy=True)

#=======================================================================================================
#
#
#
#=======================================================================================================

def calc_MPAS_grid_stat( filename, ds_in = None ):
    
    if ds_in == None:
        check = os.path.isfile(filename)
        if check:
            ds = xr.open_dataset(filename)
        else:
            print("\n CALC_MPAS_GRID_STAT: mpas file not found, exiting...\n")
            sys.exit(1)
    else:
        ds = ds_in

    try:
        isgridfile = ds.xCell.values
    except KeyError:
        print("\n CALC_MPAS_GRID_STAT: mpas file does not have coordinates, exiting...\n")
        sys.exit(1)   

    latC = np.rad2deg(ds.latCell.values)
    lonC = np.rad2deg(ds.lonCell.values)
    xC = ds.xCell.values
    yC = ds.yCell.values
    zE = ds.zgrid.values[0,:]
    zC = 0.5*(zE[:,1:] + zE[:,:-1])
    areaC = ds.areaCell.values / (3000**2)
    xV = ds.xVertex.values
    yV = ds.yVertex.values
    dis_center = ds.dcEdge.values
    
    print("\n")

    print("MPAS Lats: ",latC.min(), latC.max(), (latC.max()-latC.min()))
    print("MPAS Lons: ",lonC.min()-180., lonC.max()-180., (lonC.max()-lonC.min()))

    print("\n")
    
    print("MPAS xC (KM): ", _scale*xC.min(), _scale*xC.max(), _scale*(xC.max()-xC.min()))
    print("MPAS yC (KM): ", _scale*yC.min(), _scale*yC.max(), _scale*(yC.max()-yC.min()))
    print("MPAS zC (KM): ", _scale*zC.min(), _scale*zC.max())
    
    print("\n")

    print("MPAS xV (KM): ", _scale*xV.min(), _scale*xV.max(), _scale*(xV.max()-xV.min()))
    print("MPAS yV (KM): ", _scale*yV.min(), _scale*yV.max(), _scale*(yV.max()-yV.min()))
    
    print("\n")

    print("MPAS Cell Area: ",areaC.min(), areaC.max())
    print("MPAS cell center distances: ", dis_center.min(), dis_center.max(), dis_center.std())

    print("\n")

    print("MPAS cell distances < 2.9 km ", np.sum(dis_center < 2900.))
    print("MPAS cell distances > 3.1 km ", np.sum(dis_center > 3100.))

    ds.close()

    return
    
#=======================================================================================================
#
#
#
#=======================================================================================================
def calc_MPAS_domain_info( ds, transformer ):

    latC = np.rad2deg(ds.latCell.values)
    lonC = np.rad2deg(ds.lonCell.values)

    lat_min, lat_max = latC.min(), latC.max()
    lon_min, lon_max = lonC.min(), lonC.max()

    lon_corners = [lon_min, lon_min, lon_max, lon_max]
    lat_corners = [lat_min, lat_max, lat_max, lat_min]

    x_corners, y_corners = transformer.transform(lon_corners, lat_corners)
        
    return np.asarray(x_corners, dtype=np.float32), np.asarray(y_corners, dtype=np.float32)

#=======================================================================================================
#
#
#
#=======================================================================================================
def calc_MPAS_new_grid( grid_filename, ds_in = None, wps_file = None, \
                        nx = None, ny = None, dx = 3000., xL_grid = None, yL_grid = None ):

    if ds_in == None:
        
        check = os.path.isfile(grid_filename)

        if check:
            ds = xr.open_dataset(grid_filename)
        else:
            print("\n CALC_MPAS_NEW_GRID: mpas file not found, exiting...\n")
            sys.exit(1)
    else:
        ds = ds_in
        
    try:
        isgridfile = ds.latCell.values
    except KeyError:
        print("\n CALC_MPAS_NEW_GRID: mpas file does not have coordinates, exiting...\n")
        sys.exit(1)   

    if nx == None and ny == None and xL == None and yL == None:
        print("\n CALC_MPAS_NEW_GRID: Not enough domain information...stopping\n")
        sys.exit(1)

    try:
        if ds.on_a_sphere == "YES":
            sphere = True
    except:
        sphere = False

    if sphere:
        
        zE = ds.zgrid.values[0,:]
        zg = 0.5*(zE[1:] + zE[:-1])
        zC = zg
    
        trans        = get_transformer(wps_file)  # gets info for map projection
        x_lcc, y_lcc = calc_MPAS_domain_info(ds, trans)
        
        xL_grid = min(x_lcc[2] - x_lcc[1], x_lcc[3] - x_lcc[0])
        yL_grid = min(y_lcc[0] - y_lcc[1], y_lcc[2] - y_lcc[3])

        if nx < 0:   # create new grid based on _dx and data grid domain
            nx = np.floor( xL_grid / dx ) - 2  # shrink it a bit
            xL = nx * dx
            ny = np.floor( yL_grid / dx ) - 2
            yL = ny * dx
            xg = x_lcc.mean() + dx*(-nx/2 + np.arange(nx))
            yg = y_lcc.mean() - dx*( ny/2 + np.arange(ny))
            
        else:  # divide the data grid domain into _nx & _ny points
            dx = 0.95*xL_grid/float(nx)
            dy = 0.95*yL_grid/float(ny)
            xg = x_lcc.mean() + dx*(-nx/2 + np.arange(nx))
            yg = y_lcc.mean() + dy*(-ny/2 + np.arange(ny))
            nx = nx
            ny = ny

        xC, yC = trans.transform(np.rad2deg(ds.lonCell.values), np.rad2deg(ds.latCell.values))

        if debug > 0:
            print("\nMPAS xC (KM): ", _scale*x_lcc.min(), _scale*x_lcc.max(), _scale*(x_lcc.max()-x_lcc.min()))
            print(" NEW XG (KM): ", _scale*xg.min(), _scale*xg.max(), _scale*(xg.max()-xg.min()))
            print("\nMPAS yC (KM): ", _scale*y_lcc.min(), _scale*y_lcc.max(), _scale*(y_lcc.max()-y_lcc.min()))
            print(" NEW YG (KM): ", _scale*yg.min(), _scale*yg.max(), _scale*(yg.max()-yg.min()))
            print("\nMPAS xCells (KM): ", _scale*xC.min(), _scale*xC.max(), _scale*(xC.max()-xC.min()))
            print("\nMPAS yCells (KM): ", _scale*yC.min(), _scale*yC.max(), _scale*(yC.max()-yC.min()))

    else:   # planner grid from idealized run?
        
        xC = ds.xCell.values
        yC = ds.yCell.values
        zC = ds.zCell.values
        
        zE = ds.zgrid.values[0,:]
        zg = 0.5*(zE[1:] + zE[:-1])
        
        xL_grid = xC.max() - xC.min()
        yL_grid = yC.max() - yC.min()
    
        if nx < 0:   # create new grid based on _dx and data grid domain
            nx = np.floor( xl_grid / dx ) - 2  # shrink it a bit
            xL = nx * dx
            ny = np.floor( yl_grid / dx ) - 2
            yL = ny * dx
            xg = xC.mean() + dx*(-nx/2 + np.arange(nx))
            yg = yC.mean() + dx*(-ny/2 + np.arange(ny))
            
        else:  # divide the data grid domain into _nx & _ny points
            dx = 0.95*xL_grid/float(nx)
            dy = 0.95*yL_grid/float(ny)
            xg = xC.mean() + dx*(-nx/2 + np.arange(nx))
            yg = yC.mean() + dy*(-ny/2 + np.arange(ny))

        if debug > 0:
            print("\nMPAS xC (KM): ", _scale*xC.min(), _scale*xC.max(), _scale*(xC.max()-xC.min()))
            print(" NEW XG (KM): ", _scale*xg.min(), _scale*xg.max(), _scale*(xg.max()-xg.min()))
            print("\nMPAS yC (KM): ", _scale*yC.min(), _scale*yC.max(), _scale*(yC.max()-yC.min()))
            print(" NEW YG (KM): ", _scale*yg.min(), _scale*yg.max(), _scale*(yg.max()-yg.min()))

    ds.close()
    
    return xg, yg, zg, xC, yC, zC
    
#=======================================================================================================
#
#
#
#=======================================================================================================

def calc_MPAS_quad_grid( data_filename, xC, yC, xg, yg, ds_in = None, out_vars = None, interp = False ):

    if ds_in == None:
        check = os.path.isfile(data_filename)
        if check:
            ds = xr.open_dataset(data_filename)
        else:
            print("\n CALC_MPAS_QUAD_GRID: MPAS data file not found, exiting...\n")
            sys.exit(1)
    else:
        ds = ds_in

    if out_vars == None:
       return 

    # create 2D new grid locations

    nx = xg.shape[0]
    ny = yg.shape[0]
    
    xx, yy = np.meshgrid(xg,yg)

    coord_G = list(zip(xx.flatten(), yy.flatten()))

    # create KDTree mapping for MPAS data points

    coord_C = list(zip(xC, yC))

    mtree = KDTree(coord_C)

    if not interp:
        if debug > 0:
            print("\n CALC_MPAS_QUAD_GRID:  Using nearest neighbor interpolation \n")
        dis, index = mtree.query(coord_G, k=1)
    else:
        if debug > 0:
            print("\n CALC_MPAS_QUAD_GRID: Using IDW 5-pt interpolation \n")
        dis, index = mtree.query(coord_G, k=5)
        wght = 1.0 / dis**2
        wsum = 1.0 / np.sum(wght, axis=1)

    ntimes        = ds.Time.shape[0]
    nlevels       = ds.nVertLevels.shape[0]
    nlevels_w     = ds.nVertLevelsP1.shape[0]
    nlevels_soil  = ds.nSoilLevels.shape[0]

    interp_arrays = {}
    
    for key in out_vars:
        
        if key == 'zgrid':  # special case

            ztmp = ds[out_vars[key]].values

            ztmp = 0.5*(ztmp[:,1:] + ztmp[:,:-1])
            fldC  = np.expand_dims(ztmp, axis=0)

        else:

            fldC = ds[out_vars[key]].values

        print(key, fldC.shape)

        if fldC.ndim == 2:
    
            if not interp:
    
                fld_interp = fldC[:,index].reshape(ntimes, ny, nx)
    
            else:
    
                fld_interp = []
    
                for n in np.arange(ntimes):
                    fld_interp.append( np.sum(wght * fldC[n].flatten()[index], axis=1) * wsum )
    
                fld_interp = np.array(fld_interp).reshape(ntimes, ny, nx)
             
            interp_arrays[key] = [2, fld_interp[:,::-1,:], ntimes, ny, nx]
        
        elif fldC.ndim == 3:
    
            fldT = np.moveaxis(fldC, -1, 1)
            mz   = fldT.shape[1]

            match mz:
                        
                case x if x == nlevels_w:
                    fldT = 0.5 * (fldT[:,1:,:] + fldT[:,:-1,:])
                    mz = mz - 1

                case x if x == nlevels:
                    pass      # normal 3D Array

                case x if x == nlevels_soil:
                    print(f"\nNow interpolating soil and albedo, nlevels = {mz}")

                case _:
                    print(f"\nWeird number of vertical levels found, nlevels = {mz}")
                    print(f"\nWeird number of vertical levels found, nlevels = {nlevels_soil}")
            
            if not interp:
    
                fld_interp = fldT[:,:,index].reshape(ntimes,mz,ny,nx)
    
            else:  #IDW
    
                fld_interp = []
                for n in np.arange(ntimes):
                    for k in np.arange(mz):
                        fld_interp.append( np.sum(wght * fldT[n,k].flatten()[index], axis=1) * wsum )
    
                fld_interp = np.array(fld_interp).reshape(ntimes,mz,ny,nx)

            match mz:
                        
                case x if x == nlevels:
                    interp_arrays[key] = [3, fld_interp[:,:,::-1,:], ntimes, mz, ny, nx]
           
                case _:
                    interp_arrays[key] = [22, fld_interp[:,:,::-1,:], ntimes, mz, ny, nx]


        else:
            print("\n CALC_MPAS_QUAD_GRID: %s variable is not yet implemented, dimensions are wrong - DIMS: %i3.3" \
                   % (out_vars[key], len(fldC.shape)))
    
    ds.close()

    return interp_arrays

#=======================================================================================================
#
#
#
#=======================================================================================================

def write_MPAS_quad_netCDF( arrays, xg, yg, zg, ntimes, outfile, latlon=False ):

    # Write to XARRAY file (netCDF4)

    print("\n")
    
    if latlon:
    
        xx, yy     = np.meshgrid(xg,yg)
        trans      = get_transformer( "", latlon2xy=False )  # gets info for map projection
        lonC, latC = trans.transform( xx, yy )
        latC = latC[::-1,:]
    
    for n, key in enumerate(arrays):
            
        if arrays[key][0] == 2:  # 2D spatial data set like precip, t2m, q2m
        
            if latlon:
            
                new = xr.DataArray( arrays[key][1], dims=['nt', 'ny', 'nx'],
                        coords={"time": (["nt"], np.arange(ntimes)),
                                "lon": (["ny","nx"], lonC),
                                "lat": (["ny","nx"], latC) } )

            
            else:
            
                new = xr.DataArray( arrays[key][1], dims=['nt', 'ny', 'nx'],
                                    coords={"time": (["nt"], np.arange(ntimes)),
                                            "x": (["nx"], xg),
                                            "y": (["ny"], yg) } )

        elif arrays[key][0] == 22:  # 2.5D spatial data set like soil

            if latlon:
                new = xr.DataArray( arrays[key][1], dims = ['nt', 'nlayer', 'ny', 'nx'],
                                    coords={"time": (["nt"], np.arange(ntimes)),
                                             "lon": (["ny","nx"], lonC),
                                             "lat": (["ny","nx"], latC),
                                           "layer": (["nlayer"], np.arange(arrays[key][1].shape[1])) } )


            else:
            
                new = xr.DataArray( arrays[key][1], dims = ['nt', 'nz', 'ny', 'nx'],
                                    coords={"time": (["nt"], np.arange(ntimes)),
                                               "x": (["nx"], xg),
                                               "y": (["ny"], yg),
                                           "layer": (["nlayer"], np.arange(arrays[key][1].shape[1])) } )

        else:  # 3D grid

            if latlon:
                new = xr.DataArray( arrays[key][1], dims = ['nt', 'nz', 'ny', 'nx'],
                                    coords={"time": (["nt"], np.arange(ntimes)),
                                             "lon": (["ny","nx"], lonC),
                                             "lat": (["ny","nx"], latC),
                                               "z": (["nz"], zg) } )


            else:
            
                new = xr.DataArray( arrays[key][1], dims = ['nt', 'nz', 'ny', 'nx'],
                                    coords={"time": (["nt"], np.arange(ntimes)),
                                            "x": (["nx"], xg),
                                            "y": (["ny"], yg),
                                            "z": (["nz"], zg) } )
            
        if n == 0:
    
            ds_new = new.to_dataset(name = key)
    
        else:
    
            ds_new[key] = new
            
        print(" Wrote %s" % key)
    
    ds_new.to_netcdf(outfile, mode='w')
    ds_new.close()
    
    print(f'\n Successfully wrote interpolated MPAS data to file: {outfile}')

#====================================================================================================
# WRF CODE
#
#
#
#=======================================================================================================
def calc_corner_point_latlon(center_lat, center_lon, e_we, e_ns, dx, dy, wpsproj, latlonproj, loc):
    
    center_x, center_y = wpsproj.transform_point(center_lon, center_lat, latlonproj)
    if loc=='ll':
        xpt = center_x - dx*e_we/2.0
        ypt = center_y - dy*e_ns/2.0
    elif loc=='lr':
        xpt = center_x + dx*e_we/2.0
        ypt = center_y - dy*e_ns/2.0
    elif loc=='ul':
        xpt = center_x - dx*e_we/2.0
        ypt = center_y + dy*e_ns/2.0
    elif loc=='ur':
        xpt = center_x + dx*e_we/2.0
        ypt = center_y + dy*e_ns/2.0
    corner_lon, corner_lat = latlonproj.transform_point(xpt, ypt, wpsproj)
    
    return corner_lon, corner_lat

def calc_center_point_latlon(corner_lat_parent, corner_lon_parent, dx_parent, dy_parent, e_we, e_ns, dx, dy, i, j, wpsproj, latlonproj):
    corner_x_parent, corner_y_parent = wpsproj.transform_point(corner_lon_parent, corner_lat_parent, latlonproj)
    center_x_child = corner_x_parent + dx_parent*i + dx*e_we/2.0
    center_y_child = corner_y_parent + dy_parent*j + dy*e_ns/2.0
    center_lon_child, center_lat_child = latlonproj.transform_point(center_x_child, center_y_child, wpsproj)
    
    return center_lon_child, center_lat_child



def reproject_corners(corner_lons, corner_lats, wpsproj, latlonproj):
    corner_x = np.zeros((4,1))
    corner_y = np.zeros((4,1))
    corner_x[0], corner_y[0] = wpsproj.transform_point(corner_lons[0], corner_lats[0], latlonproj)
    corner_x[1], corner_y[1] = wpsproj.transform_point(corner_lons[1], corner_lats[1], latlonproj)
    corner_x[2], corner_y[2] = wpsproj.transform_point(corner_lons[2], corner_lats[2], latlonproj)
    corner_x[3], corner_y[3] = wpsproj.transform_point(corner_lons[3], corner_lats[3], latlonproj)

    return corner_x, corner_y

# all these functions below are necessary only when LCC projection is used.
def find_side(ls, side):
    """
    Given a shapely LineString which is assumed to be rectangular, return the
    line corresponding to a given side of the rectangle.
    
    """
    minx, miny, maxx, maxy = ls.bounds
    points = {'left': [(minx, miny), (minx, maxy)],
              'right': [(maxx, miny), (maxx, maxy)],
              'bottom': [(minx, miny), (maxx, miny)],
              'top': [(minx, maxy), (maxx, maxy)],}
    return sgeom.LineString(points[side])


def lambert_xticks(ax, ticks, size):
    """Draw ticks on the bottom x-axis of a Lambert Conformal projection."""
    te = lambda xy: xy[0]
    lc = lambda t, n, b: np.vstack((np.zeros(n) + t, np.linspace(b[2], b[3], n))).T
    xticks, xticklabels = _lambert_ticks(ax, ticks, 'bottom', lc, te)
    ax.xaxis.tick_bottom()
    ax.set_xticks(xticks)
    ax.set_xticklabels([ax.xaxis.get_major_formatter()(xtick) for xtick in xticklabels], size=size)
    

def lambert_yticks_left(ax, ticks, size):
    """Draw ricks on the left y-axis of a Lamber Conformal projection."""
    te = lambda xy: xy[1]
    lc = lambda t, n, b: np.vstack((np.linspace(b[0], b[1], n), np.zeros(n) + t)).T
    yticks, yticklabels = _lambert_ticks(ax, ticks, 'left', lc, te)
    ax.yaxis.tick_left()
    ax.set_yticks(yticks)
    ax.set_yticklabels([ax.yaxis.get_major_formatter()(ytick) for ytick in yticklabels], size=size)

    
def lambert_yticks_right(ax, ticks, size):
    """Draw ricks on the left y-axis of a Lamber Conformal projection."""
    te = lambda xy: xy[1]
    lc = lambda t, n, b: np.vstack((np.linspace(b[0], b[1], n), np.zeros(n) + t)).T
    yticks, yticklabels = _lambert_ticks(ax, ticks, 'right', lc, te)
    ax.yaxis.tick_right()
    ax.set_yticks(yticks)
    ax.set_yticklabels([ax.yaxis.get_major_formatter()(ytick) for ytick in yticklabels], size=size)

def _lambert_ticks(ax, ticks, tick_location, line_constructor, tick_extractor):
    """Get the tick locations and labels for an axis of a Lambert Conformal projection."""
    outline_patch = sgeom.LineString(ax.outline_patch.get_path().vertices.tolist())
    axis = find_side(outline_patch, tick_location)
    n_steps = 30
    extent = ax.get_extent(ccrs.PlateCarree())
    _ticks = []
    for t in ticks:
        xy = line_constructor(t, n_steps, extent)
        proj_xyz = ax.projection.transform_points(ccrs.Geodetic(), xy[:, 0], xy[:, 1])
        xyt = proj_xyz[..., :2]
        ls = sgeom.LineString(xyt.tolist())
        locs = axis.intersection(ls)
        if not locs:
            tick = [None]
        else:
            tick = tick_extractor(locs.xy)
        _ticks.append(tick[0])
    # Remove ticks that aren't visible:    
    ticklabels = copy(ticks)
    while True:
        try:
            index = _ticks.index(None)
        except ValueError:
            break
        _ticks.pop(index)
        ticklabels.pop(index)
    return _ticks, ticklabels

