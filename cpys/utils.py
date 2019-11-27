# THIS LIBRARY CONTAINS SMALL TOOLS HELPFUL IN HANDLING I/0
# Kajetan Chrapkiewicz 2017

## MODULES 
import numpy as np
from math import *
from subprocess import call
import ntpath # FOR path_leave
import itertools
import sys
#from scipy.io import FortranFile
# MY LIBS
from my_signal_processing import *
from constants import *
#from my_plotlib import *
##

## FUNCTIONS

def define_checker(dims, dx, pads, std, ampl):
  if len(dims) == 2:
    nx1, nx3 = dims
    dx1, dx3 = dx # GAPS BETWEEN ANOMALIES' CENTRES
    pads1, pads3 = pads[0 : 2], pads[2 : ] # DISTANCES FROM MODEL MARGINS
    std1, std3 = std
  else:
    eprint('ERROR! 3D NOT YET ADDRESSED')
    quit()
  
  X1, A1 = diracs_comb_1d(nx1, dx1, pads1, ampl)
  X3, A3 = diracs_comb_1d(nx3, dx3, pads3, ampl)
  #l = 50
  #w = 15
  #y = np.zeros((l))
  #y[l / 2 - w - 1 : l / 2 + w] = 1
  ##plt.plot(y)
  ##plt.plot(x, y)
  ##plt.show()
  ##quit()
  #conv_array = np.convolve(A1, y, mode = 'same')
  #plt.plot(conv_array)
  #plt.show()
  #quit()
  A1 = unnormed_filter_gaussian(A1, std1)
  A3 = unnormed_filter_gaussian(A3, std3)
  #plt.plot(A3)
  #plt.show()
  #quit()

  chckr = np.ones((len(X1), len(X3)))
  
  for i in range(len(X1)):
    for j in range(len(X3)):
      chckr[i][j] *= A1[i] * A3[j] / float(ampl) # NORMALIZATION!
  
  X, Z = np.meshgrid(X1, X3)
  chckr = chckr.transpose()
  #plot_model_vtr(nx1, nx3, X, Z, chckr, 'seismic')
  #plt.savefig('chckrbrd.png', format = 'png', dpi = 150, bbox_inches = 'tight')
  #plt.show()
  #quit()
  return X, Z, chckr

# ADD CHECKERBOARD TO THE MODEL OF THE SAME DIMENSIONS
def add_chckr(model, chckr):
  nx3 = len(model)
  nx1 = len(model[0])
  nmodel = np.zeros((nx3, nx1))
  #print len(model), len(chckr)
  for z in np.arange(nx3):
    for x in np.arange(nx1):
      nmodel[z][x] = model[z][x] + chckr[z][x] * model[z][x]
  return nmodel




# FROM NUMBER OF MULTIDIM. LISTS PICK OUT ELEMENTS WITH x AS A FIRST INDEX
def cherrypick_x(x, *lists):
  nlists = []
  for l in lists:
    xy = l[x]
    nlists.append(xy)
  return nlists

# FROM NUMBER OF MULTIDIM. LISTS PICK OUT ELEMENTS WITH (x, y) AS 2 FIRST INDICES 
def cherrypick_xy(x, y, *lists):
  nlists = []
  for l in lists:
    xy = l[x][y]
    nlists.append(xy)
  return nlists

# READ ALL DATA NEEDED TO UPDATE GHOST WAVEFIELD (OUTPUT FILE OF fsprep)
def read_ghost_data(file_name):
  content = read_file(file_name)
  
  # READ HEADER
  x_start, nx, x_step, y_start, ny, y_step = content[0]
  # CAST STRINGS TO CORRECT TYPES
  x_start, nx, y_start, ny = [int(i) for i in [x_start, nx, y_start, ny]]
  x_step, y_step = [float(i) for i in [x_step, y_step]]
  storder = int(content[1][0])
  ##
  
  # READ DATA
  data = content[2: ]
  ##
  
  # DERIVE 'SIZES' OF RESPECTIVE DATA TYPES

  if ny > 1:
    auxs_per_fict = (aux_no) ** 2 # PER FICT. POINT
  else:
    auxs_per_fict = aux_no
  ghost_len = ghosts_no + intersects_no + ficts_no + ficts_no * auxs_per_fict
  ##
  
  all_lvls = []
  line = 0
  for glvl in range(storder):
    ghosts = []
    intersects = []
    ficts = []
    auxs = []    
    for x in range(nx):
      x_ghosts, x_intersects, x_ficts, x_auxs = [], [], [], []
      for y in range(ny):
        ghost_data = data[line : line + ghost_len]
        
        x_ghosts.append(ghost_data[0])
        x_intersects.append(ghost_data[1])
        x_ficts.append(ghost_data[2 : 2 + ficts_no])
        
        all_auxs = ghost_data[2 + ficts_no :]
        xx_auxs = []
        k = 0
        for i in range(ficts_no):
          x_f_auxs = []
          for j in range(auxs_per_fict):
            x_f_auxs.append(all_auxs[k])
            #print 'x_f', len(x_f_auxs)
            k += 1
          #print i + 1, 'final x_f', len(x_f_auxs)
          xx_auxs.append(x_f_auxs)
        #print 'final xx_auxs', len(xx_auxs)
        x_auxs.append(xx_auxs)
        line += ghost_len
        
      ghosts.append(x_ghosts)
      intersects.append(x_intersects)
      ficts.append(x_ficts)
      auxs.append(x_auxs)
    ghosts = [[[float(i) for i in j] for j in k] for k in ghosts]
    intersects = [[[float(i) for i in j] for j in k] for k in intersects]
    ficts = [[[[float(i) for i in j] for j in k] for k in l] for l in ficts]
    auxs = [[[[[float(i) for i in j] for j in k] for k in l] for l in m] for m in auxs]
    
    all_lvls.append([ghosts, intersects, ficts, auxs])
  #ghosts = [[[float(i) for i in j] for j in k] for k in ghosts]
  #intersects = [[[float(i) for i in j] for j in k] for k in intersects]
  #ficts = [[[[float(i) for i in j] for j in k] for k in l] for l in ficts]
  #auxs = [[[[[float(i) for i in j] for j in k] for k in l] for l in m] for m in auxs]
 
  return all_lvls #ghosts, intersects, ficts, auxs

# TURN LIST OF TUPLES INTO 3 LISTS OF THEIR COORDINATES
def split_tuples(tuples):
  X, Y, Z = [], [], []
  for tupl in tuples:
    X.append(tupl[0])
    Y.append(tupl[1])
    Z.append(tupl[2])
  return X, Y, Z

# CUT FS TO LEAVE ONLY AREA WITHIN radius FROM (x, y) <- INDICES (MAKE IT MORE GENERAL)
def cut_fs(X, Y, Z, x, y, radius):
  xmin = x - radius + 1
  if xmin < 0:
    xmin = 0  
  xmax = x + radius 

  
  X = X[xmin : xmax]
  
  if len(Y) > 1:
    ymin = y - radius + 1
    if ymin < 0:
      ymin = 0    
    ymax = y + radius   
    Y = Y[ymin : ymax]
    Z = [i[ymin : ymax] for i in Z[xmin : xmax]]  
  else:
    Z = Z[xmin : xmax]
  return X, Y, Z

# READ FS GRID FILE
def read_fs_grd(file_name):
  content = read_file(file_name)
  x_start, nx, x_step, y_start, ny, y_step = content[0]
  # CAST STRINGS TO CORRECT TYPES
  x_start, nx, y_start, ny = [int(i) for i in [x_start, nx, y_start, ny]]
  x_step, y_step = [float(i) for i in [x_step, y_step]]
  ##
  X = np.arange(x_start, x_start + x_step * nx, x_step)
  if ny > 1:
    Y = np.arange(y_start, y_start + y_step * ny, y_step)
  else:
    Y = [1]
  data = content[1: ]
  data = np.ravel(data) # FLATTEN THE LIST
  data = [float(i) for i in data]
  
  x_slices = get_chunks(data, ny)
  Z = []
  for x_slice in x_slices:
    Z.append(x_slice)
  return X, Y, Z
  
# READ FS GRID FILE
def read_fs_grd_old(file_name):
  content = read_file(file_name)
  x_start, nx, x_step, y_start, ny, y_step = content[0]
  # CAST STRINGS TO CORRECT TYPES
  x_start, nx, y_start, ny = [int(i) for i in [x_start, nx, y_start, ny]]
  x_step, y_step = [float(i) for i in [x_step, y_step]]
  ##
  X = np.arange(x_start, x_start + x_step * nx, x_step)
  if ny > 1:
    Y = np.arange(y_start, y_start + y_step * ny, y_step)
  else:
    Y = [1]
  data = content[1: ]
  print data
  data = np.ravel(data)
  
  data = [float(i) for i in data]
  print 'after'
  print data
  
  if ny > 1:
    y_slices = get_chunks(data, nx)
    Z = []
    for y_slice in y_slices:
      Z.append(y_slice)
  elif ny == 1:
    Z = data
  else:
    eprint('ERROR. STRANGE VALUE OF ny')
  
  return X, Y, Z

# WRITE FULL SegyPrep.key WITH FIXED REGULAR GEOMETRY
def write_segyprep(segyprep_runfile, geom, rec_x, rec_z, sou):
  nx1, nx2, nx3 = geom
  
  recz = int(rec_z[0]) * 10 # INTERIM SOLUTION
  recx0 = int(rec_x[0])  * 10
  recnx = len(rec_x) 
  if len(rec_x) > 1:
    recdx = int(rec_x[1] - rec_x[0])  * 10
  else:
    recdx = 0
  recy0 = 1  
  recny = 1
  recdy = 1  
  
  souz = int(sou[1])  * 10
  soux0 = int(sou[0])  * 10
  sounx = 1
  soudx = int(nx1 / 2)  * 10 #??????
  
  f = open(segyprep_runfile, 'w')
  f.write('io      : fw3d\n')
  f.write('problem : synthetic\n')
  f.write('nx1     : ' + str(nx1) + '\n')
  f.write('nx2     : ' + str(nx2) + '    \n')
  f.write('nx3     : ' + str(nx3) + '  \n')
  f.write('dx      : 10    \n')
  f.write('\n')
  f.write('ttime   : 1000  \n')
  f.write('dtms    : 1.0  \n')
  f.write('\n')
  f.write('geometry : regular\n')
  f.write('fixed array : yes\n')
  f.write('souz   : ' + str(souz) + ' \n')
  f.write('recz   : ' + str(recz) + ' \n')
  f.write('\n')
  f.write('soux0  : ' + str(soux0) + ' \n')
  f.write('souy0  : 1\n')
  f.write('sounx  : ' + str(sounx) + '\n')
  f.write('souny  : 1 \n')
  f.write('soudx  : ' + str(soudx) + '\n')
  f.write('soudy  : 1\n')
  f.write('\n')
  f.write('recx0  : ' + str(recx0) + ' \n')  
  f.write('recy0  : 1\n')
  f.write('recnx  : ' + str(recnx) + '  \n')
  f.write('recny  : 1\n')
  f.write('recdx  : ' + str(recdx) + ' \n')
  f.write('recdy  : 1\n')
  f.close()

# READ SegyPrep.key AND RETURN PARAMS NEEDED IN rotate_fs
def read_segyprep(segyprep_runfile):  
  params = read_runfile(segyprep_runfile)
  #
  nx1 = int(params['nx1'])
  nx2 = int(params['nx2'])
  nx3 = int(params['nx3'])
  
  sou_depth = int(params['souz']) / 10
  soux0 = int(params['soux0']) / 10
  
  rec_depth = int(params['recz']) / 10
  recx0 = int(params['recx0']) / 10
  rec_dx = int(params['recdx']) / 10
  recnx = int(params['recnx'])
  
  if sou_depth != 2 * rec_depth:
    eprint('WARNING. SOURCE SHOULD BE TWICE AS DEEP AS RECEIVERS\n')
  
  if (sou_depth - rec_depth) % rec_dx != 0:
    eprint('WARNING: sou will not be rotated properly!\n')
  
  return nx1, nx2, nx3, sou_depth, soux0, rec_depth, recx0, rec_dx, recnx

# READ Runfile.key AND RETURN RECORDS OF FORMAT: '1-word-key : value' AS A DICTIONARY 
def read_runfile(file_name):
  content = read_file(file_name)
  params = {} # DICTIONARY
  for record in content:
    # GET ONLY RECORDS OF CERTAIN FORMAT:
    if (len(record) == 3) and (record[1] == ':'):
      key = record[0]
      value = record[2]
      params[key] = value
  print 'my_lib.py/read_runfile: ', file_name, ' loaded as a dictionary.'
  return params

# ASSUME btop = etep = 0, take rest_value FOR THE REST
def update_runfile_bounds(file_name, rest_value):
  update_runfile(file_name, ('btop', 0), ('extratop', 0),
                 ('Bbottom', rest_value), ('extrabottom', rest_value), # CAREFUL WITH CAPITALS!
                 ('Bleft', rest_value), ('extraleft', rest_value),
		 ('Bright', rest_value), ('extraright', rest_value),
		 ('Bfront', rest_value), ('extrafront', rest_value),
		 ('Bback', rest_value), ('extraback', rest_value))
		              
# CHANGE RECORDS OF RUNFILE, THE REST REMAINS EXACTLY THE SAME
def update_runfile(file_name, *records): 
  content = read_not_split(file_name)
  f = open(file_name, 'w')
  for line in content:
    nline = line
    line = line.split(None)
    if len(line) == 0:
      continue
    for r in records:
      if line[0] == r[0]:
        nline = r[0] + ' : ' + str(r[1]) + '\n'
        break
    f.write(nline)
  #print records[0]
  #f = open(file_name + 'test', 'w')
  #for key in params:
    #print "key: %s, value: %s" % (key, params[key])
  f.close()
  
# STRAIGHT 2D LINE OF GIVEN SLOPE ANGLE alpha AND FIRST VALUE y0
def straight(x_series, alpha, y0):
  x0 = x_series[0]
  
  A = np.tan(deg2rad(alpha))
  B = y0 - A * x0
  
  return [A * x + B for x in x_series]

# PREPARE 3D GRIDDED FREE SURFACE GAUSSIAN IN X DIRECTION
def prepare_X_sin_fs_grd(file_name, x_start, nx, x_step, y_start, ny, y_step, z_min, amp, lambd):
  X, Y, Z = prepare_flat_fs_grd(file_name, x_start, nx, x_step, y_start, ny, y_step, 1)
  Z = Z.transpose() # TO MODIFY IN X DIRECTION AND KEEP CONSTANT IN Y DIRECTION
  for i in range(len(Y)):
    Z[i] = [(z_min + amp) + amp * np.sin(2 * np.pi * x / lambd) for x in X]
  Z = Z.transpose() # CONVERT BACK TO SHAPE COMPATIBLE WITH OTHER READ/WRITE FUNCTION
  return X, Y, Z

# PREPARE 3D GRIDDED FREE SURFACE GAUSSIAN IN X DIRECTION
def prepare_X_gauss_fs_grd(file_name, x_start, nx, x_step, y_start, ny, y_step, z_min, z_max, mu, sigma):
  X, Y, Z = prepare_flat_fs_grd(file_name, x_start, nx, x_step, y_start, ny, y_step, 1)
  Z = Z.transpose() # TO MODIFY IN X DIRECTION AND KEEP CONSTANT IN Y DIRECTION
  for i in range(len(Y)):
    Z[i] = [z_min + (z_max - z_min) * gauss1d(x, mu, sigma) for x in X]
  Z = Z.transpose() # CONVERT BACK TO SHAPE COMPATIBLE WITH OTHER READ/WRITE FUNCTION
  return X, Y, Z
  
# PREPARE 3D GRIDDED FREE SURFACE TILTED IN X DIRECTION
def prepare_X_tilted_fs_grd(file_name, x_start, nx, x_step, y_start, ny, y_step, alpha, y0):
  X, Y, Z = prepare_flat_fs_grd(file_name, x_start, nx, x_step, y_start, ny, y_step, 1)
  Z = Z.transpose() # TO MODIFY IN X DIRECTION AND KEEP CONSTANT IN Y DIRECTION
  for i in range(len(Y)):
    Z[i] = straight(X, alpha, y0)
  Z = Z.transpose() # CONVERT BACK TO SHAPE COMPATIBLE WITH OTHER READ/WRITE FUNCTION
  return X, Y, Z
  
# PREPARE 3D FLAT (HORIZONTAL) FREE SURFACE GRID FILE
def prepare_flat_fs_grd(file_name, x_start, nx, x_step, y_start, ny, y_step, z_value):
  X = np.arange(x_start, x_step * nx + 1, x_step)
  Y = np.arange(y_start, y_step * ny + 1, y_step) 
  
  Z = np.ones((len(X), len(Y))) * z_value
  return X, Y, Z

# SAVE FILES WITH FS BEFORE ROTATION
def prepare_flat_geometry(fs_depth, rec_depth, rec_dx, n, pad, step):
  x_step, y_step = step
  
  xleft = pad #* rec_dx # MARGINS
  xright = xleft       #
  
  rec_x = np.arange(1, (2 * n + 1) * rec_dx + 1, rec_dx)
  rec_x += xleft
  rec_z = np.zeros(len(rec_x)) + rec_depth + fs_depth 
  rec = zip(rec_x, rec_z)
  #print rec
  
  fs_x = np.arange(1, rec_x[-1] + xright + 1)
  #fs_z = np.zeros(len(fs_x)) + fs_depth
  #fs = zip(fs_x, fs_z)
  
  x_start = 1
  x_end = len(fs_x)
  nx = (x_end - x_start + 1) / x_step
  
  # NOW ONLY 2D
  y_start = 1
  ny = 1  
  ##
  
  X_fs, Y_fs, Z_fs = prepare_flat_fs_grd('dummy', x_start, nx, x_step, y_start, ny, y_step, fs_depth)
  fs = [X_fs, Y_fs, Z_fs]
  
  sou = (rec_x[len(rec) / 2], 2 * rec_depth + fs_depth) # INDICES START FROM 0
  #print sou
  
  return rec, fs, sou, nx, ny

def prepare_layered_model():
  print 'OK'
  #bottoms = np.linspace(0, nx3, nlay)
  ##print vtop, vbot
  #values = np.linspace(vtop, vbot, nlay)
  #bottoms = [int(round(b)) for b in bottoms] # ROUND TO THE NEAREST INTEGER
  ##print bottoms
  ##print values
  ## TMP SOLUTION
  #vels = np.linspace(vtop, vbot, nx3)
  ##print np.linspace 
  ##print len(vels)
  ##print len(np.arange(nx1))
  #
  #for i1 in np.arange(nx1):
  #  for i2 in np.arange(nx2):
  #    lin1 = np.linspace(1500, 4000, 40)
  #    lin2 = np.linspace(4000, 6000, 61)
  #    model[i1, i2, 0:10] = 1500
  #    model[i1, i2, 10:50] = np.linspace(2000, 4000, 40)
  #    model[i1, i2, 50:101] = np.linspace(5000, 6000, 51)
  #  
  #for i3 in np.arange(nx3):
  #  print model[0, 0, i3]
  
  
  #from my_lib import *
  
  ## FUNCTIONS
  #def grad_between2_interfaces(htop, hbot, vtop, vbot, nlays):
    #depths = np.linspace(htop, hbot, nlays)
    #vels = np.linspace(vtop, vbot, nlays)
    
    #print depths
    #print vels
  
  #def continuous(i, depths):
    #if depths[i] == depths[i + 1]
  
  #def read_interfaces(file_name):
    #content = read_file(file_name)
    #depths = [row[0] for row in content]
    #vels = [row[1] for row in content]
    #return depths, vels
  
  ### MAIN
  #file_names = sys.argv[1:] 
  
  #for file_name in file_names:
  
    #grad_between2_interfaces(0, 10, 1.5, 4.5, 10)
    #depths, vels = read_interfaces(file_name)
  print 'OK2'

def make_model_vtr(project_name, nx1, nx2, nx3):
  from scipy.io import FortranFile
  print 'make_model_vtr: nx1, nx2, nx3', nx1, nx2, nx3
  model = np.zeros(shape = (nx1, nx2, nx3), dtype = np.float32)
  
  for i3 in np.arange(nx3):
    for i2 in np.arange(nx2):
      for i1 in np.arange(nx1):
        model[i1, i2, i3] = 3000.0 # m/s
  ###
  
  ## WRITE THE MODEL TO A FILE 
  file_name = project_name + '-TrueVp.vtr'
  f = FortranFile(file_name, 'w')
  
  f.write_record(np.array([1, 2, 0], dtype = np.int32)) # OR 1, 3, 0 WHEN USING NX2 AS WELL 
  f.write_record(np.array([nx3, nx1], dtype = np.int32)) # OR nx3, nx2, nx1 (NOT COMPATIBLE WITH LLUIS)
  
  for i1 in np.arange(nx1):
    for i2 in np.arange(nx2):
      trace = model[i1, i2, :]
      f.write_record(np.array(trace, dtype = np.float32))
  
  f.close()  
  
# SAVE FS IN .vtr FORMAT (2D ONLY)
def save_vtr(file_name, nx, ny, Z):
  from scipy.io import FortranFile
  f = FortranFile(file_name + '.vtr', 'w')
  if ny > 1:
    f.write_record(np.array([1, 2, 0], dtype = np.int32)) 
    f.write_record(np.array([ny, nx], dtype = np.int32)) 
  elif ny == 1:
    f.write_record(np.array([1, 1, 0], dtype = np.int32)) 
    f.write_record(np.array([nx], dtype = np.int32)) 
  else:
    eprint('ERROR. STRANGE VALUE OF ny.')
    sys.exit(1)
  for i in np.arange(nx):
    trace = Z[i, :]
    f.write_record(np.array(trace, dtype = np.float32))  
  
  f.close()  

def read_bin(file_name):
  from scipy.io import FortranFile
  f = FortranFile(file_name, 'r')
  #for i in np.arange(nx):
    #trace = Z[i, :]
    #print i, len(trace)
    #f.write_record(np.array(trace, dtype = np.float32))  
  f.close()
# 2D
def save_bin(file_name, nx, ny, Z):
  f = FortranFile(file_name + '.bin', 'w')
  for i in np.arange(nx):
    trace = Z[i, :]
    print i, len(trace)
    f.write_record(np.array(trace, dtype = np.float32))  
  f.close()

# 2D ONLY AT THE MOMENT
def read_vtr(file_name):
  from scipy.io import FortranFile
  f = FortranFile(file_name, 'r')
  print 'File ', file_name, 'is about to be read.'
  
  header1 = f.read_ints(dtype = np.int32)
  nx3, nx1 = f.read_ints(dtype = np.int32) # OMIT NX2 IN LLUIS CODE?
  print nx3, nx1
  axes_ratio = nx1 / float(nx3)
  model = []
  for x in np.arange(nx1):
    trace = [] # EMPTY JUST IN CASE
    trace = f.read_reals(dtype = np.float32)
    model.append(trace)
  model = np.transpose(model)
  
  f.close()
  print 'File ', file_name, 'is being closed'
  
  return nx1, nx3, model

# SAVE FS WITH HEADER TO .grd FILE (2D/3D)
def save_grd(file_name, x_start, nx, x_step, y_start, ny, y_step, Z):
  f = open(file_name + '.grd', 'w')
  f.write('{:6}'.format(x_start) + '{:6}'.format(nx) + '{:6}'.format(x_step) + 
          '{:6}'.format(y_start) + '{:6}'.format(ny) + '{:6}'.format(y_step) + '\n')
  for i in range(nx):
    for j in range(ny):
      f.write('{:15.8f}'.format(Z[i, j]) + '\n')
  f.close()

# UNNORMALIZED 2D GAUSSIAN FUNCTION (VALUES 0-1)
def gauss2d(x, y, mu_x, mu_y, sig_x, sig_y):
  return np.exp(-((x - mu_x)**2) / (2 * sig_x**2) - ((y - mu_y)**2) / (2 * sig_y**2))

# UNNORMALIZED 1D GAUSSIAN FUNCTION (VALUES 0-1)
def gauss1d(x, mu, sigma):
  return np.exp(-((x - mu)**2) / (2 * sigma**2)) # (1 / np.sqrt(2 * np.pi * sigma**2)) *

# EUCLIDEAN DISTANCE BETWEEN 2 POINTS IN 2D
def dist(point1, point2):
  return np.sqrt( (point1[0] - point2[0])**2 + (point1[1] - point2[1])**2 )

# CONVERT DEGREES TO RADIANS
def deg2rad(angle):
  return angle * np.pi / 180 

# CONVERT RADIANS TO DEGREES
def rad2deg(angle):
  return angle * 180 / np.pi 

# FIND a AND b OF y = ax +b GIVEN TWO POINTS
def find_linear_function(p1, p2):
  x1, y1 = p1 
  x2, y2 = p2
  dy = y2 - y1
  dx = x2 - x1
  a = dy / float(dx)
  b = y1 - x1 * dy / float(dx)
  return a, b

# ROTATE A POINT (px, py) CLOCKWISE BY A GIVEN ANGLE (DEGREES) AROUND A GIVEN ORIGIN (ox, oy)
def rotate2d(origin, point, angle):
  ox, oy = origin
  px, py = point
  angle = deg2rad(angle)
  
  qx = ox + np.cos(angle) * (px - ox) - np.sin(angle) * (py - oy)
  qy = oy + np.sin(angle) * (px - ox) + np.cos(angle) * (py - oy)
  return qx, qy

# SAVE 2D FS TO .pgy FILE 
def save_fs_2d_pgy(file_name, nx1, nx2, nx3, X, Z):
  Z = np.ravel(Z)
  X = [float(i) for i in X]
  Z = [float(i) for i in Z]
  
  
  n = len(X)# nx1 * nx2 # NO. OF POINTS
  sn   = '{:10}'.format(n)
  snx3 = '{:15}'.format(nx3)
  snx2 = '{:15}'.format(nx2)
  snx1 = '{:15}'.format(nx1)
  
  f = open(file_name + '.pgy', 'w')
  f.write(sn + snx3 + snx2 + snx1 + '\n')
  for i in range(len(X)):
    si   = '{:10}'.format(i + 1)
    sz   = '{:15.8f}'.format(Z[i])
    sy   = '{:15.8f}'.format(1.0)
    sx   = '{:15.8f}'.format(X[i])  
    #print si, sz, sy, sx
    f.write(si + sz + sy + sx + '\n')
    i += 1
  f.close()

# PREPARE .pgy FREE SURFACE FILE
def prepare_flat_fs_2d_pgy(file_name, nx1, nx2, nx3):
  n = nx1 * nx2 # NO. OF POINTS
  sn   = '{:10}'.format(n)
  snx3 = '{:15}'.format(nx3)
  snx2 = '{:15}'.format(nx2)
  snx1 = '{:15}'.format(nx1)
  
  f = open(file_name, 'w')
  f.write(sn + snx3 + snx2 + snx1 + '\n')
  i = 1
  for y in range(1, nx2 + 1):
    for x in range(1, nx1 + 1):
      si   = '{:10}'.format(i)
      sz   = '{:15.8f}'.format(1.)
      sy   = '{:15.8f}'.format(y)
      sx   = '{:15.8f}'.format(x)  
      print si, sz, sy, sx
      f.write(si + sz + sy + sx + '\n')
      i += 1
  f.close()

# READ .pgy FILE: HEADER AND 3 COLUMNS OF DATA SEPARATELY
def read_pgy(file_name):
  content = read_file(file_name)
  header = content[0]
  data = content[1: ]
  nx3, nx2, nx1 = [int(i) for i in header[1: ]]
  x = [float(row[3]) for row in data]
  y = [float(row[2]) for row in data]
  z = [float(row[1]) for row in data]
  return nx1, nx2, nx3, x, y, z
  
# READ .pgy FREE SURFACE FILE
def read_fs_2d_pgy(file_name):
  nx1, nx2, nx3, x, y, z = read_pgy(file_name)
  yslices = get_chunks(z, nx1)
  fs = []
  for yslice in yslices:
    fs.append(yslice)
  fs = np.array(fs)
  fs = np.transpose(fs)
  return nx1, nx2, nx3, fs

# DIVIDE LIST l INTO n-ELEMENT PIECES
def get_chunks(l, n):
  for i in range(0, len(l), n):
    yield l[i : i + n]    

# PRINT TO STDERR
def eprint(*args, **kwargs):
  sys.stderr.write(*args, **kwargs)

# EXTRACT A FILE NAME FROM A WHOLE PATH
def path_leave(path):
  head, tail = ntpath.split(path)
  return tail or ntpath.basename(head)

# CHECK IF THE STRING STANDS FOR A FLOAT NUMBER AND RETURN TRUE/FALSE
def is_number(s):
  try:
    float(s)
    return True
  except ValueError:
    return False

# CHECK IF THE STRING STANDS FOR A INT NUMBER AND RETURN TRUE/FALSE
def is_int(s):
  try: 
    int(s)
    return True
  except ValueError:
    return False

# SEARCH THE LIST OF "NUMBERS" TO REPLACE NAN WITH ZEROS AND RETURN CORRECTED LIST
def NaN2zero(content):
  new_content = []
  for c in content:
    if isnan(c):
      c = 0.0
    new_content.append(c)
  return new_content

# SEARCH THE LIST OF "NUMBERS" TO REPLACE NAN WITH ZEROS AND RETURN CORRECTED LIST
def NaN2number(content, number):
  new_content = []
  for c in content:
    if isnan(c):
      c = number
    new_content.append(c)
  return new_content

# SEARCH THE LIST OF NUMBERS FROM THE END AND RETURN THE LENGTH OF ZERO-TAIL
def empty_ending(li):
  count = 0
  for i in reversed(li):
    if i == 0.:
      count += 1
    else:
      break
  return count

# READ A FILE EXACTLY AS IT IS
def read_not_split(file_name):
  f = open(file_name, 'r')
  content = []
  for line in f:
    content.append(line)
  f.close()
  return content

# READ A FILE NON-ZERO LINES, SPLIT THEM WITH None AND RETURN CONTENT 
def read_file(file_name):
  f = open(file_name, 'r')
  print 'my_lib.py/read_file: ', file_name, ' opened.' 
  content = []
  for line in f:
    line = line.split(None)
    if len(line) != 0:
      content.append(line)
  f.close()
  return content

def read_file_1st_char(file_name):
  f = open(file_name, 'r')    
  for line in f:
    line = line.split(None)
    first_char = line[0]
    break
  f.close()
  return first_char

def read_n_columns(file_name, n):
  f = open(file_name, 'r')
  content = []
  for line in f:
    content.append(line.split(None))
  f.close()  
  columns = []
  for i in np.arange(n):
    col = [line[i] for line in content]
    columns.append(col)
  return columns

def read_only_vs_files(fnames):
  depths_list, vs_list = [], []
  for fname in fnames:
    content = read_file(fname)
    depths = [line[0] for line in content]
    vels = [line[1] for line in content]
    # DOUBLE EACH DEPTH TO GET STAIR-CASE LOOK OF PLOT !!!!!!!!!!!
    new_depths = []
    for d in depths:
      new_depths.append(d)
      if float(d) != 0: # DONT DOUBLE FOR DEPTH = O KM
        new_depths.append(d)
    # ASSIGN PROPER VS FOR EACH BOUNDARY
    new_vels = []
    for i in np.arange(len(vels)):
      if i > 0:
        new_vels.append(vels[i - 1])
      new_vels.append(vels[i])
    depths_list.append(new_depths)
    vs_list.append(new_vels)
  return [depths_list, vs_list]

# READ CPS MODEL FILES AND RETURN [DEPTHS_LIST, VS_LIST]
def read_mod(mod_files, d_below_last_boundary): 
  depths_list, vs_list = [], []
  for mod_file in mod_files:
    content = read_file(mod_file)
    data = content[mod_header_length: ]
    
    depths = []
    vs = []
    h = 0.
    i = 0
    
    for line in data:
      depths.append(h)
      vs.append(data[i][2])
      h += float(line[0])
      if i != (len(data) - 1): 
        depths.append(h)
        vs.append(data[i][2])
      else:
        depths.append(h + d_below_last_boundary)
        vs.append(data[i][2])
      i += 1
    depths_list.append(depths)
    vs_list.append(vs)
    
  return [depths_list, vs_list]

# READ CPS RESOLUTION MATRIX
def read_R_file(file_name):
  content = read_file(file_name)
  R, row = [], []
  for i in np.arange(len(content) - 2) + 2: # START FROM FIRST DATA LINE (2nd line) TO SKIP PROBLEMS WITH (i-1) FOR i=0
    if is_int(content[i][0]):  
      row = list(itertools.chain.from_iterable(row)) # MAKE row A FLAT LIST
      row = [float(i) for i in row]
      R.append(row) 
      row = []
    elif not is_int(content[i - 1][0]):
      row.append(content[i])
  R = np.array(R)
  R = np.swapaxes(R, 1, 0)
  R = np.flipud(R)
  R = np.fliplr(R)

  return R # MAKE NUMPY ARRAY FROM LIST OF LISTS

# GET THE CORE OF THE STRING (REMOVE PREFIX AND SUFFIX FROM IT)
def core(string, prefix, suffix): # try variable number of args...
  beginning = filter(None, string.split(suffix))[0]
  return filter(None, beginning.split(prefix))[0] # prefix must by space (' ', not '') if we dont want any

def lengthen_matrix(M, xmin0, xmax0, xmin, xmax):
  print xmin0, xmax0 # INITIAL LIMITS OF EACH MATRIX ROW
  print xmin, xmax # DESIRED LIMITS OF EACH MATRIX ROW (OBJECTIVE OF THIS FUNCTION)
  #print len(M[0]), len(M[-1])
  step = (xmax0 - xmin0) / len(M[0]) # CONSTANT DIFFERENCE BETWEEN SUBSEQUENT ELEMENTS OF THE ROW
  n_add = int(ceil((xmax - xmax0) / step)) # NUMBER OF ELEMENTS TO ADD AT THE END OF EACH ROW
  print 'n_add: ', n_add
  exact_xmax = xmax0 + n_add * step 
  #n_cut = floor((xmin - xmin0) / step) # NUMBER OF ELEMENTS TO SUBTRACT AT THE BEGINNING OF EACH ROW
  end = np.zeros(n_add) # NEW END OF EACH ROW
  new_M = []
  for row in M:
    new_M.append(np.concatenate((row, end))) # GLUE 2 NP.ARRAYS
  new_M = np.array(new_M) # MAKE NP.ARRAY OUT OF LIST
  
  return new_M, exact_xmax

# TURN A LIST OF SUBLISTS INTO A FLAT LIST
def flatten_list(list2d):
  return list(itertools.chain.from_iterable(list2d))

# PRINT UNKNOWN NUMBER OF ARGUMENTS  
def all_args(*arg):
  print "I was called with", len(arg), "arguments:", arg

# iTH COLUMN OF THE LIST
def ith_col(i, li):
  return [l[i] for l in li]
##


#  combined formulas from (Berteussen, 1977) and (Gardner et al., 1974). 
def density(vp):
  rho = 0.77 + 0.32 * vp + 0.09 * (((vp - 6.) ** 4 + 0.1) ** 0.25 - (vp - 6.))
  return rho
