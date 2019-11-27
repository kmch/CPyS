#!/usr/bin/env python

#
## MODULES 
#
import sys
from pylab import *
from utils import *
from plot import *
from constants import *
from scipy.interpolate import *
#
## MAIN
#
matrix_file = sys.argv[1] # tmpsrfi.02
model_file = sys.argv[2] # .DAT FILE !!!
#option = sys.argv[3]

core = filter(None, model_file.split(".dat"))[0]

depths = read_only_vs_files([model_file])[0][0][::2] # [::2] - TAKE EVERY 2nd ITEM
x, y = depths, depths
x = [float(i) for i in x]
y = [float(i) for i in y]
#print y
#x[-1] = 300 # IN FACT IT'S ~298 KM
#y[-1] = 300 # IN FACT IT'S ~298 KM
X, Y = np.meshgrid(x, y)


R = read_R_file(matrix_file)
R = np.clip(R, 0, clipper)#1e6) # VALUES < 0 BECOME 0, VALUES > 1e6 BECOME 1e6

x_max, y_max = [], []

glob_max = 0

for row in np.arange(len(R)):
  if row != 0:
    row_max = max(R[row])
    no_layer = list(R[row]).index(row_max)
    if row_max > glob_max:
      glob_max = row_max
    #print depths[row], depths[no_layer]
    x_max.append(depths[row])
    y_max.append(depths[no_layer])

# PLOT 
fig, ax = plt.subplots(1, 1, sharex = True)
fig.set_size_inches(square_side, square_side, forward = True) # IN INCHES

colors = [(255, 255, 255), (0, 0, 255), (0, 246, 255), (255, 222, 0), (255, 0, 0)]#, (255, 0, 0), (255, 0, 0), (255, 0, 0)] # , (25, 0, 116)
c = mcolors.ColorConverter().to_rgb
c1 = 'white'
c2 = '#00035b' #'#004577'#'#00555a' #'#005f6a' #'petrol' #'#0e87cc' #water lbue#'#047495' #'sea blue''blue' #midnightblue'
c3 = '#41fdfe' #'brightcyan' #'cyan'
c4 = '#9cef43' #(kiwi) #lime
c5 = '#f1da7a' #(sandy) #'yellow'
c6 = '#ffb07c' #'peach' #'orange'
c7 = '#fd8d49' #orangeish'orangered'
c8 = 'red'
c9 = 'maroon'
c10 = 'black'

my_cmap = make_colormap(
    [c(c1), 0.01, c(c1), c(c2), 0.1, c(c2), c(c3), 0.2, c(c3), c(c4), 0.3, c(c4), 
     c(c5), 0.4, c(c5), c(c6), 0.5, c(c6), c(c7), 0.6, c(c7), c(c8), 0.7, c(c8), 
     c(c9), 0.8, c(c9), c(c10), 0.9, c(c10)])

cax = pcolormesh(X, Y, R, cmap = my_cmap, alpha = .5)#'nipy_spectral') #'Greys', 'nipy_spectral'
plt.scatter(y_max, x_max, marker = 'o', facecolors = 'none', edgecolors = 'r', s = 100, label = 'max for a given layer')
plt.plot(np.arange(0, x[-1]), np.arange(0, y[-1]), c = 'k', lw = 2, label = 'main diagonal')

plt.legend(loc = 'lower left', frameon = True, prop={'size' : 12})
plt.axis('square') #!!!!!!!         
ax.tick_params(direction = 'in', length = tick_length, width = tick_width, colors = 'k')
plt.xlim(0, x[-1])
plt.ylim(0, y[-1])
plt.grid()
fig.tight_layout()
ax.set_ylim(ax.get_ylim()[::-1])

#print glob_max
#cax = ax.imshow(R, cmap = 'Greys', interpolation = 'none', aspect = 'auto') #extent=[v_min, v_max, prof, 0], )
cb = plt.colorbar(cax, orientation = 'horizontal')#, ticks=[0, 0.235])# plt.colorbar(im,fraction=0.046, pad=0.04)
tick_locator = ticker.MaxNLocator(nbins = 6)
cb.locator = tick_locator
cb.update_ticks()
plt.xlabel('depth [km]')
plt.ylabel('depth [km]') 

print core

plt.savefig('../resol_' + core + '.png', dpi = 300)
#plt.savefig('../resol_' + core + '.svg', format = 'svg', dpi = 300) # dots per inches
plt.savefig('../resol_' + core + '.ps', format = 'ps', dpi = 300)
plt.show()