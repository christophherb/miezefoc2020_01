import numpy as np
from scipy.interpolate import Rbf
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from matplotlib import cm
from matplotlib import pyplot as plt
import matplotlib.axes._axes as axes
import matplotlib.figure as figure
from ill_datafile import ill_datafile
from generate_figure import generate_figure

fig, ax = generate_figure((10, 4))

dtx_vals = np.array([k for k in range(10, 260, 30)]+[400, 600])+470
sth_vals = np.linspace(-0.2, -1.8, 81)

sth_opt = -1.12
all_amplitudes = np.loadtxt('fit_data/amps')[:-2,:]
all_amplitudesErr = np.loadtxt('fit_data/amp_err')
all_sigmas = np.loadtxt('fit_data/centers')
all_centers = np.loadtxt('fit_data/sigmas')
sth_vals = np.array(sth_vals)-sth_opt
dtx_vals = dtx_vals[:-2]/10
x, y = np.meshgrid(sth_vals, dtx_vals) 
xi = np.linspace(sth_vals[0]+0.01, sth_vals[-1]-0.01, 100)
yi = np.linspace(dtx_vals[0]-15, dtx_vals[-1]+15, 100)
xx, yy = np.meshgrid(xi, yi)

xf = x.flatten()
yf = y.flatten()
zf = all_amplitudes.flatten()
print('xf', xf.shape)
print('xi', xi.shape)
rbf = Rbf(xf, yf, zf, function='linear')
zz = rbf(xx, yy)
zi = griddata((xf, yf), zf, (xi[None, :], yi[:, None]), method='nearest')
print(zz.shape)
#z = x*np.exp(-x**2-y**2)
#ti = np.linspace(-2.0, 2.0, 100)
#scat = ax.scatter(x, y, 100, z, cmap=cm.jet,)
#ax.set_xlabel('x')
#ax.set_ylabel('y')
#xx, yy = np.meshgrid(ti, ti)
#
#
#

im = ax.pcolor(xx, yy, zz)

ax.scatter(x, y, 10, c = 'black')

ax.set_title('RBF interpolation - multiquadrics')



fig.colorbar(im)
plt.show()

