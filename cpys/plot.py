# THIS LIBRARY CONTAINS SMALL TOOLS HELPFUL IN PLOTTING
# Kajetan Chrapkiewicz 2017

## MODULES
from pylab import *
import matplotlib # for log xticks
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.ticker import MaxNLocator
from matplotlib.ticker import FuncFormatter
import matplotlib.colors as mcolors # FOR OWN COLORMAPS
from matplotlib import ticker
from mpl_toolkits.mplot3d import Axes3D
# MY OWN
from utils import *
from constants import *
##

## FUNCTIONS  

# DISPLAY ONLY EVERY nTH TICK FOR BOTH X AND Y AXIS
def display_every_nth_tick(ax, n):
  i = 1
  for tic in ax.xaxis.get_major_ticks():
    if (i % n) != 0:
      tic.label1On = tic.label2On = False
    i += 1
  i = 1
  for tic in ax.yaxis.get_major_ticks():
    if (i % n) != 0:
      tic.label1On = tic.label2On = False
    i += 1


# FOR ALL GHOSTS OF glvl LEVEL, PLOT ALL AUXILIARY POINTS (I, F, A)
def show_all_ghosts(glvl, fs_file, gh_file):
  X_fs, Y_fs, Z_fs = read_fs_grd(fs_file) # DON'T CUT IT, UNLIKE IN show_single_ghost
  all_lvls = read_ghost_data(gh_file)
  ghosts, intersects, ficts, auxs = all_lvls[glvl]
  xmin = margin
  xmax = int(X_fs[-1]) - margin
  if len(Y_fs) > 1:
    ymin = margin
    ymax = int(Y_fs[-1]) - margin
  else:
    ymin = 1
    ymax = 2
  zmin = min(Z_fs)[0] - glvl - 2
  zmax = max(Z_fs)[0] + margin + 1
  
  for x in range(xmin, xmax):
    for y in range(ymin, ymax):
      plot_all_auxs(x, y, X_fs, Y_fs, Z_fs, ghosts, intersects, ficts, auxs)
      plt.ylim(zmax, zmin)
      plt.gca().set_aspect('equal', adjustable='box')
      plt.savefig(str(x).rjust(7, '0') + '-' + str(y) + '.png', format = 'png')
      plt.close()

# PLOT ALL AUXILIARY POINTS (I, F, A) ASSOCIATED WITH A SINGLE GHOST NODE AGAINST A LOCAL FRAGMENT OF FS
def show_single_ghost(x, y, glvl, fs_file, gh_file):
  X_fs, Y_fs, Z_fs = read_fs_grd(fs_file)
  X_fs, Y_fs, Z_fs = cut_fs(X_fs, Y_fs, Z_fs, x, y, node_rad)
  all_lvls = read_ghost_data(gh_file)
  ghosts, intersects, ficts, auxs = all_lvls[glvl]
  plot_all_auxs(x, y, X_fs, Y_fs, Z_fs, ghosts, intersects, ficts, auxs)
  plt.show()

# PLOT ALL AUXILIARY POINTS (I, F, A) ASSOCIATED WITH A SINGLE GHOST NODE
def plot_all_auxs(x, y, X_fs, Y_fs, Z_fs, ghosts, intersects, ficts, auxs):
  if len(Y_fs) > 1: # 3D (HANDLE TICKS)
    fig, ax = plot_fs_3d(X_fs, Y_fs, Z_fs)
    ghost, intersect, ficts, auxs = cherrypick_xy(x, y, ghosts, intersects, ficts, auxs)
    plot_auxs_3d(ax, ghost, intersect, ficts, auxs)
  
  else: # 2D
    Z_fs = np.ravel(Z_fs)
    fig, ax = plot_fs_2d(X_fs, Z_fs)
    plt.xticks(np.arange(min(X_fs), max(X_fs) + 1, 1.0))
    plt.yticks(np.arange(int(min(Z_fs)) - 25 , int(max(Z_fs)) + 25, 1.0)) # CLIPPED BY ylim ANYWAY
    # MAKE XTICK LABELS SPARSER
    i = 1
    for tic in ax.xaxis.get_major_ticks():
      if (i % 5) != 0:
        tic.label1On = tic.label2On = False
      i += 1
    ##
    # MAKE YTICK LABELS SPARSER
    i = 1
    for tic in ax.yaxis.get_major_ticks():
      if (i % 2) != 0:
        tic.label1On = tic.label2On = False
      i += 1
    ##
    
    ghost, intersect, ficts, auxs = cherrypick_x(x, ghosts, intersects, ficts, auxs)
    plot_auxs_2d(ax, ghost[0], intersect[0], ficts[0], auxs[0]) # [0] <= COMPATIBILITY W/ 3D  

# PLOT ALL 3D AUXILIARY POINTS (I, F, A) ASSOCIATED WITH A SINGLE GHOST NODE
def plot_auxs_3d(ax, ghost, intersect, ficts, auxs):
  colors = iter(cm.rainbow(np.linspace(0, 1, len(ficts))))
  ax.scatter(ghost[0], ghost[1], ghost[2], c = 'b')
  ax.scatter(intersect[0], intersect[1], intersect[2], c = 'm')
  
  X, Y, Z = split_tuples(ficts)
  ax.scatter(X, Y, Z, c = 'k', s = psize)
  # PLOT NORMAL
  ax.plot([ghost[0], X[-1]], [ghost[1], Y[-1]], [ghost[2], Z[-1]], c = 'grey', lw = 3, alpha = 0.5)
  ##
  
  for a in auxs:
    X, Y, Z = split_tuples(a)
    ax.scatter(X, Y, Z, c = next(colors))    

# PLOT ALL 2D AUXILIARY POINTS (I, F, A) ASSOCIATED WITH A SINGLE GHOST NODE
def plot_auxs_2d(ax, ghost, intersect, ficts, auxs):
  ax.scatter(ghost[0], ghost[2], c = 'k', s = psize)
  ax.scatter(intersect[0], intersect[2], c = 'w', s = psize)
  # PLOT FICTS
  X, Y, Z = split_tuples(ficts)
  ax.scatter(X, Z, c = 'k', s = psize)
  # PLOT NORMAL
  ax.plot([ghost[0], X[-1]], [ghost[2], Z[-1]], c = 'k')#, lw = nwidth)#, alpha = nopac)
  ##
  # PLOT AUXS
  #colors = iter(cm.rainbow(np.linspace(0, 1, len(ficts)))) # IT SEEMS AUXS ARE IN WRONG ORDER 
  zmax = epsi
  for a in auxs:
    X, Y, Z = split_tuples(a)
    if max(Z) > zmax:
      zmax = max(Z)
    ax.scatter(X, Z, c = 'b', s = psize) #next(colors))  
  plt.ylim(zmax + pad, ghost[2] - pad)
  
def plot_multiple_2d_fs(file_names):
  fig, ax = plt.subplots(1, 1, sharex = True)
  fig.set_size_inches(longer_side, square_side, forward = True) # IN INCHES
  
  plt.xlabel('X [nodes]')
  plt.ylabel('Z [nodes]')
  plt.grid()    
  xmin = big
  xmax = epsi
  zmin = big
  zmax = epsi
  
  
  colors = iter(cm.rainbow(np.linspace(0, 1, len(file_names))))
  for file_name in file_names:
    X, Y, Z = read_fs_grd(file_name)
    Z = np.ravel(Z)
    
    if min(X) < xmin:
      xmin = min(X)
    if max(X) > xmax:
      xmax = max(X)
    if min(Z) < zmin:
      zmin = min(Z)
    if max(Z) > zmax:
      zmax = max(Z)    
    
    ax.plot(X, Z, c = next(colors), lw = 2, label = file_name)
  
  plt.legend(loc = 'upper right', frameon = False, prop = {'size' : 20})
  plt.xlim(xmin - pad, xmax + pad)
  plt.ylim(zmin - pad, zmax + pad)
  ax.set_ylim(ax.get_ylim()[::-1])
  
# PLOT 2D/3D FREE SURFACE
def plot_fs(file_name):
  X, Y, Z = read_fs_grd(file_name)
  
  if len(Y) < 2: 
    Z = np.ravel(Z) # FLATTEN THE LIST
    plot_fs_2d(X, Z) 
  else:
    plot_fs_3d(X, Y, Z)

# MAKE INTERACTIVE 3D SURFACE PLOT GIVEN COORDINATES 1D ARRAYS (X, Y) AND 2D DATA ARRAY (Z)   
def plot_fs_3d(X, Y, Z):
  XX, YY, Z = prepare_3d_mesh(X, Y, Z)
  fig, ax = set_3d_window(X, Y)
  
  ax.plot_surface(XX, YY, Z, color = 'dodgerblue', alpha = 0.1)
  return fig, ax

# PLOT 2D FREE SURFACE PLOT
def plot_fs_2d(X, Z):
  fig, ax = plt.subplots(1, 1, sharex = True)
  fig.set_size_inches(longer_side, square_side, forward = True) # IN INCHES
  
  plt.xlim(min(X) - pad, max(X) + pad)
  #plt.ylim(max(Z) + pad, 0)
  
  plt.xlabel('X [nodes]')
  plt.ylabel('Z [nodes]')
  plt.grid()
  ax.plot(X, Z, c = 'k', lw = 5, alpha = 0.4)
  ax.set_ylim(ax.get_ylim()[::-1])
  #plt.gca().set_aspect('equal', adjustable='box')
  return fig, ax


def plot_shots(x, z):
  plt.scatter(x, z, s = 200, c = 'y', marker = '*')

def plot_receivers(x, z):
  plt.scatter(x, z, s = 200, c = 'aquamarine', marker = 'v')

# NOW ONLY 2D
def plot_model_vtr(nx1, nx3, X, Z, model, colormap, minn, maxx): # 1500 7500
  font = {'size' : 40} #40
  matplotlib.rc('font', **font)
  
  axes_ratio = float(nx1) / nx3
  
  if colormap == 'seismic0':
    seismic = matplotlib.cm.seismic 
    shiftedColorMap(seismic, minn, maxx, name = 'seismic0')
    colormap = 'seismic0'
  
  #vmin = minn, vmax = maxx, cmap = 'seismic0'
  fig, ax = plt.subplots(1, 1, sharex = True)
  fig.set_size_inches(square_side * axes_ratio, square_side, forward = True) # IN INCHES
  cax = pcolormesh(X, Z, model, cmap = colormap, vmin = minn, vmax = maxx) #inferno
  #tickss = np.linspace(minn, maxx, 5, endpoint = True)
  cb = plt.colorbar(cax, orientation = 'vertical')#, ticks = tickss)
  #cb.ax.set_title('m/s')
  #cb.set_ticks(np.linspace(minn, maxx + 1, 5))
  #tick_locator = ticker.MaxNLocator(nbins = 8)
  #cb.locator = tick_locator 
  #cb.update_ticks()
  
  #dx = 0.05 # km
  
  plt.xlim(0, nx1 - 1)
  plt.ylim(0, nx3 - 1)
  #plt.xticks(np.arange(dx * nx1))
  #plt.yticks(np.arange(dx * nx3))
  plt.xlabel('in-line [km]')
  plt.ylabel('depth [km]')
  plt.tight_layout()
  ax.set_ylim(ax.get_ylim()[::-1])
  return fig, ax

# PROBABLY deprecated 
def plot_ghost_data_2d(ghosts, intersects, ficts, auxs):
  ghosts, intersects, ficts, auxs = read_ghost_data(gh_file)
  ghost = ghosts[x][y]
  intersect = intersects[x][y]
  ficts = ficts[x][y]
  auxs = auxs[x][y]
  #print 'intersect: ', intersect
  #print x, y
  #print 'ficts', ficts[x][y] 
  #for i in range(len(ficts)):
    #print 'fict. no ', i + 1, ': ', ficts[i][ :3]
  #print auxs[2]
  plot_1ghost_data_3d(X_fs, Y_fs, Z_fs, ghost, intersect, ficts, auxs)

# PROBABLY deprecated
def plot_fs_ghosts_2d(X_fs, Z_fs, X_gh, Z_gh):
  Z_gh = np.ravel(Z_gh) # BURDEN OF 3D 
  fig, ax = plt.subplots(1, 1, sharex = True)
  fig.set_size_inches(square_side, square_side, forward = True) # IN INCHES
  
  plt.xlim(-pad, max(X_gh) + pad)
  plt.ylim(max(Z_gh) + pad, 0)
  
  plt.xlabel('X [nodes]')
  plt.ylabel('Z [nodes]')
  plt.xticks(np.arange(min(X_fs), max(X_fs) + 1, 1.0))
  plt.grid()
  plt.plot(X_fs, Z_fs, lw = 5, c = 'red', alpha = .5)
  plt.scatter(X_gh, Z_gh) 

# MAKE INTERACTIVE 3D PLOT OF GHOSTS GIVEN COORDINATES 1D ARRAYS (X, Y) AND 2D DATA ARRAY (Z)   
def plot_ghosts_3d(X_fs, Y_fs, Z_fs, X_gh, Y_gh, Z_gh):
  XX_fs, YY_fs = np.meshgrid(X_fs, Y_fs)
  Z_fs = np.array(Z_fs)
  Z_fs = Z_fs.reshape(XX_fs.shape)
  
  XX_gh, YY_gh = np.meshgrid(X_gh, Y_gh)
  Z_gh = np.array(Z_gh)
  Z_gh = Z_gh.reshape(XX_gh.shape)
  
  fig = plt.figure()
  ax = fig.add_subplot(111, projection = '3d')
  
  minn = min(Z_gh.min(axis = 1))
  if minn < 0:
    minn = 0
  maxx = max(Z_gh.max(axis = 1))
  
  plt.xlim(min(X_gh), max(X_gh))
  plt.ylim(min(Y_gh), max(Y_gh))
  ax.set_zlim(minn, maxx)
  
  #ax.set_xlabel("X [nodes]")
  #ax.set_ylabel("Y [nodes]")
  #ax.set_zlabel("Z [nodes]")
  ax.invert_zaxis() # FLIP Z AXIS
  ax.scatter(XX_gh, YY_gh, Z_gh)
  ax.plot_surface(XX_fs, YY_fs, Z_fs, color = 'r', alpha = .1) # 'dodgerblue'

# PREPARE MESH 
def prepare_3d_mesh(X, Y, Z):
  XX, YY = np.meshgrid(X, Y)
  Z = np.array(Z)
  Z = Z.transpose() # CRUCIAL
  Z = Z.reshape(XX.shape)
  return XX, YY, Z

# CREATE AND SET PARAMETERS OF 3D PLOT'S WINDOW 
def set_3d_window(X, Y):
  from mpl_toolkits.mplot3d import Axes3D
  fig = plt.figure()
  ax = fig.add_subplot(111, projection = '3d')
  plt.xlim(min(X), max(X))
  plt.ylim(min(Y), max(Y))
  #ax.set_xlabel("X [nodes]")
  #ax.set_ylabel("Y [nodes]")
  #ax.set_zlabel("Z [nodes]")
  ax.invert_zaxis() # FLIP Z AXIS
  return fig, ax
  
# 2D PROJECTION OF 3D DATA
def plot_2d_projection(X, Y, Z): 
  XX, YY = np.meshgrid(X, Y)
  
  fig, ax = plt.subplots(1, 1, sharex = True)
  fig.set_size_inches(longer_side, square_side, forward = True) # IN INCHES
  cax = pcolormesh(XX, YY, Z, cmap = 'viridis')
  plt.xlim(min(X), max(X))
  plt.ylim(min(Y), max(Y))
  
  plt.xlabel('X')
  plt.ylabel('Y')
  #ax.set_ylim(ax.get_ylim()[::-1]) # FLIP Y-AXIS
  cb = plt.colorbar(cax, orientation = 'horizontal')


def shiftedColorMap(cmap, minn, maxx, start=0, stop=1.0, name='shiftedcmap'):
  print ''
#  (FROM STACK OVERFLOW)
#  Function to offset the "center" of a colormap. Useful for
#  data with a negative min and positive max and you want the
#  middle of the colormap's dynamic range to be at zero
#
#  Input
#  -----
#    cmap : The matplotlib colormap to be altered
#    start : Offset from lowest point in the colormap's range.
#        Defaults to 0.0 (no lower ofset). Should be between
#        0.0 and `midpoint`.
#    midpoint : The new center of the colormap. Defaults to 
#        0.5 (no shift). Should be between 0.0 and 1.0. In
#        general, this should be  1 - vmax/(vmax + abs(vmin))
#        For example if your data range from -15.0 to +5.0 and
#        you want the center of the colormap at 0.0, `midpoint`
#        should be set to  1 - 5/(5 + 15)) or 0.75
#    stop : Offset from highets point in the colormap's range.
#        Defaults to 1.0 (no upper ofset). Should be between
#        `midpoint` and 1.0.
  from mpl_toolkits.axes_grid1 import AxesGrid
  if (maxx == 0) and (minn == 0):
    midpoint = 0.5
  else:
    midpoint = 1 - maxx / (maxx + abs(minn))
  
  cdict = {
      'red': [],
      'green': [],
      'blue': [],
      'alpha': []
  }

  # regular index to compute the colors
  reg_index = np.linspace(start, stop, 257)

  # shifted index to match the data
  shift_index = np.hstack([
      np.linspace(0.0, midpoint, 128, endpoint=False), 
      np.linspace(midpoint, 1.0, 129, endpoint=True)
  ])

  for ri, si in zip(reg_index, shift_index):
      r, g, b, a = cmap(ri)

      cdict['red'].append((si, r, r))
      cdict['green'].append((si, g, g))
      cdict['blue'].append((si, b, b))
      cdict['alpha'].append((si, a, a))

  newcmap = matplotlib.colors.LinearSegmentedColormap(name, cdict)
  plt.register_cmap(cmap=newcmap)

  return newcmap

def forceAspect(ax,aspect=1):
  im = ax.get_images()
  extent =  im[0].get_extent()
  ax.set_aspect(abs((extent[1]-extent[0])/(extent[3]-extent[2]))/aspect)

def my_formatter(x, pos):
  if x.is_integer():
    return str(int(x))
  else:
    return str(x)

def cm2inch(value):
  print ''
  #def cm2inch(*tupl):
  #inch = 2.54
  #if isinstance(tupl[0], tuple):
  #  return tuple(i/inch for i in tupl[0])
  #else:
  #  return tuple(i/inch for i in tupl)
  return value / 2.54

def make_cmap(colors, position=None, bit=False):
  import matplotlib as mpl
#     make_cmap takes a list of tuples which contain RGB values. The RGB
#     values may either be in 8-bit [0 to 255] (in which bit must be set to
#     True when called) or arithmetic [0 to 1] (default). make_cmap returns
#     a cmap with equally spaced colors.
#     Arrange your tuples so that the first color is the lowest value for the
#     colorbar and the last is the highest.
#     position contains values from 0 to 1 to dictate the location of each color.
  
  bit_rgb = np.linspace(0,1,256)
  if position == None:
      position = np.linspace(0,1,len(colors))
  else:
      if len(position) != len(colors):
          sys.exit("position length must be the same as colors")
      elif position[0] != 0 or position[-1] != 1:
          sys.exit("position must start with 0 and end with 1")
  if bit:
      for i in range(len(colors)):
          colors[i] = (bit_rgb[colors[i][0]],
                       bit_rgb[colors[i][1]],
                       bit_rgb[colors[i][2]])
  cdict = {'red':[], 'green':[], 'blue':[]}
  for pos, color in zip(position, colors):
      cdict['red'].append((pos, color[0], color[0]))
      cdict['green'].append((pos, color[1], color[1]))
      cdict['blue'].append((pos, color[2], color[2]))

  cmap = mpl.colors.LinearSegmentedColormap('my_colormap',cdict,256)
  return cmap

def make_colormap(seq):
    """Return a LinearSegmentedColormap
    seq: a sequence of floats and RGB-tuples. The floats should be increasing
    and in the interval (0,1).
    """
    seq = [(None,) * 3, 0.0] + list(seq) + [1.0, (None,) * 3]
    cdict = {'red': [], 'green': [], 'blue': []}
    for i, item in enumerate(seq):
        if isinstance(item, float):
            r1, g1, b1 = seq[i - 1]
            r2, g2, b2 = seq[i + 1]
            cdict['red'].append([item, r1, r2])
            cdict['green'].append([item, g1, g2])
            cdict['blue'].append([item, b1, b2])
    return mcolors.LinearSegmentedColormap('CustomMap', cdict)

def plot_mod_normal(x_list, y_list, ax):
  colors = iter(cm.rainbow(np.linspace(0, 1, len(y_list))))
  for i in np.arange(len(y_list)):
    # SET COLOR AND LABEL FOR LEGEND
    colr = next(colors)
    if i == 0:
      colr = 'k'
      lab = 'mean'
    elif i == 1:
      lab = 'mode'
    elif i == 2:
      lab = 'true'
    elif i == 3:
      lab = 'ref'
    # PLOT
    ax.plot(y_list[i], x_list[i], color = colr, lw = lwidth, label = lab) 

def plot_mod_envelope(x, y1, y2, ax):
  x = [float(s) for s in x]
  y1, y2 = np.array([float(s) for s in y1]), np.array([float(s) for s in y2])

  ax.plot(y1, x, color = 'grey', lw = 1)
  ax.plot(y2, x, color = 'grey', lw = 1)
  ax.fill_betweenx(x, y1, y2, color = 'grey', alpha = 0.1)   

def plot_mod_envelope_color(x, y1, y2, ax, col, alph):
  x = [float(s) for s in x]
  y1, y2 = np.array([float(s) for s in y1]), np.array([float(s) for s in y2])

  ax.plot(y1, x, color = col, lw = 1, ls = 'None')
  ax.plot(y2, x, color = col, lw = 1, ls = 'None')
  ax.fill_betweenx(x, y1, y2, color = col, alpha = alph) 

def plot_stat_post(ax, d_max, files):
  depths_list, vs_list = read_only_vs_files(files)
  # by default first 2 files define borders of the envelope (stdm, stdp in this order), the rest is plotted normally
  x = depths_list[0] # assuming depths are identical in all files!!!!!
  y1, y2 = vs_list[0], vs_list[1]
  x_list = depths_list[2: ]
  y_list = vs_list[2: ]
  
  # To propagate the size change to an existing gui window add forward=True
  
  plot_mod_envelope(x, y1, y2, ax)
  plot_mod_normal(x_list, y_list, ax)
  
  plt.xlim(vmin, vmax)
  plt.ylim(0, d_max)
  ax.set_ylim(ax.get_ylim()[::-1]) # must be below ?lim
  #plt.legend(loc = 'lower left', frameon = False, prop={'size':8}
  
  plt.xlabel('$\mathregular{v_S}$ [km/s]')
  plt.ylabel('depth [km]')
  plt.grid()

def plot_2_envelopes(ax, d_max, files):
  depths_list, vs_list = read_only_vs_files(files)
  # ORDER OF THE FILES stdm1, stdp1, mean1, stdm2, stdp2, mean2
  x1 = depths_list[0] # assuming depths are identical in all files!!!!!
  mean1, stdm1, stdp1  = vs_list[0], vs_list[1],  vs_list[2]
  
  x2 = depths_list[3] # assuming depths are identical in all files!!!!!
  mean2, stdm2, stdp2 = vs_list[3], vs_list[4],  vs_list[5]
  
  plot_mod_envelope_color(x1, stdm1, stdp1, ax, 'red', 0.5)
  ax.plot(mean1, x1, color = 'red', lw = lwidth, alpha = opacity, label = 'option 1')
  
  plot_mod_envelope_color(x2, stdm2, stdp2, ax, 'blue', 0.5)
  ax.plot(mean2, x2, color = 'blue', lw = lwidth, label = 'option 2')
  
  plt.xlim(vmin, vmax)
  plt.ylim(0, d_max)
  ax.set_ylim(ax.get_ylim()[::-1]) 
  plt.legend(loc = 'lower left', frameon = True, prop = {'size' : 15})
  
  plt.xlabel('$\mathregular{v_S}$ [km/s]')
  plt.ylabel('depth [km]')
  plt.grid()

def plot_only_vs(ax, fnames):
  depths_list, vs_list = read_only_vs_files(fnames)
  
  for i in np.arange(len(depths_list)):
    plt.plot(vs_list[i], depths_list[i])
  
  ax.set_ylim(ax.get_ylim()[::-1])

def equalize_series_length(series1, series2):
  l1, l2 = len(series1), len(series2)
  if l1 > l2:
    series1 = series1[ :l2]
  else:
    series2 = series2[ :l1]
  return [series1, series2]
    
def compare_2_series(x, y1, y2, label1, label2, xlabel, ylabel):
  fig, ax = plt.subplots(1, 1, sharex=True)
  fig.set_size_inches(longer_side, shorter_side, forward=True) # in inches
  
  y1, y2 = equalize_series_length(y1, y2) # equalization of lengths!
  y1 = np.array(y1) # list -> array (otherwise 'dimensions are inconsisten')
  y2 = np.array(y2) ## #+ 0.018 <- for comparison of raw kumar and rfr
  
  #samples = np.arange(len(y1)) # converting samples to time domain
  #x = samples / sampling_freq - delay
  #x = samples
  
  ax.plot(x, y1, color = 'r', lw = 2, label = label1)
  ax.plot(x, y2, color = 'b', lw = 2, label = label2)
  ax.fill_between(x, y1, y2, where=y2 >= y1, facecolor='blue', interpolate=True, alpha = .5)
  ax.fill_between(x, y1, y2, where=y2 <= y1, facecolor='red', interpolate=True, alpha = .5)
  plt.xlabel(xlabel)
  plt.ylabel(ylabel)
  plt.legend(loc = 'upper right', frameon = False, prop = {'size' : 15})
  plt.grid()

# RUBBISH

def plot_intersects(x, y, radius, fs_file, gh_file, ghost_data_file):
  xmin = x - radius
  xmax = x + radius - 1
  ymin = y - radius
  ymax = y + radius - 1
  
  X_fs, Y_fs, Z_fs = read_fs_grd(fs_file)
  X_gh, Y_gh, Z_gh = read_fs_grd(gh_file)
  X_in, Y_in, Z_in = read_intersects(ghost_data_file)
  
  XX_fs, YY_fs = np.meshgrid(X_fs, Y_fs)
  Z_fs = np.array(Z_fs)
  Z_fs = Z_fs.reshape(XX_fs.shape)
  
  XX_gh, YY_gh = np.meshgrid(X_gh, Y_gh)
  Z_gh = np.array(Z_gh)
  Z_gh = Z_gh.reshape(XX_gh.shape)  
  
  fig = plt.figure()
  ax = fig.add_subplot(111, projection = '3d')  
  
  plt.xlim(min(X_fs), max(X_fs))
  plt.ylim(min(Y_fs), max(Y_fs))
  
  XX_fs = [i[xmin : xmax] for i in XX_fs[ymin : ymax]]
  YY_fs = [i[xmin : xmax] for i in YY_fs[ymin : ymax]]
  Z_fs = [i[xmin : xmax] for i in  Z_fs[ymin : ymax]]
  
  XX_gh = [i[xmin : xmax] for i in XX_gh[ymin : ymax]]
  YY_gh = [i[xmin : xmax] for i in YY_gh[ymin : ymax]]
  Z_gh = [i[xmin : xmax] for i in  Z_gh[ymin : ymax]]  

  ax.invert_zaxis() # FLIP Z AXIS
  ax.plot_surface(XX_fs, YY_fs, Z_fs, color = 'r', alpha = .5)  
  ax.scatter(X_in, Y_in, Z_in, c = 'r', marker = '.')
  ax.scatter(XX_gh, YY_gh, Z_gh)

