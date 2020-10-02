################################### 3Dplot #################################################
%matplotlib qt
from nicosload import read_nicos_file,read_pad_file
import matplotlib.ticker as plticker
from matplotlib import cm
dtx_vals = [k for k in range(10,260,30)]+[400,600]
sth_vals = np.linspace(-0.2,-1.8,81)

all_amplitudes = np.loadtxt('fit_data/amps')
all_amplitudesErr = np.loadtxt('fit_data/amp_err')
all_sigmas = np.loadtxt('fit_data/centers')
all_centers = np.loadtxt('fit_data/sigmas')

fig = plt.figure()
ax = fig.add_subplot(111,projection = '3d',)
colormap = plt.cm.viridis
#ax.set_prop_cycle(color =[colormap(i) for i in np.linspace(0, 1, len(all_amplitudes))][::-1])  
xs,ys,zs = [],[],[]
for line in range(9):
    xs.append(np.array([dtx_vals[line]+470]*81))
    zs.append(all_amplitudes[line,:]/10**6)
    ys.append(np.linspace(-0.2,-1.8,81))
X,Y = np.meshgrid(xs,ys)
Z = all_amplitudes
surf = ax.plot_surface(X, Y, Z, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)
ax.set_xlabel('dtx (mm)',rotation = 45)
ax.set_ylabel('sth (deg)')
ax.set_zlabel('peak amplitude')

ax.axes.xaxis.labelpad=10
ax.axes.yaxis.labelpad=15
ax.axes.zaxis.labelpad=10
#ax.set_axis_on()
#ax.set_xticklabels([])
ax.axes.xaxis.line.set_color("black")
ax.axes.xaxis.line.set_marker(' ')
ax.axes.yaxis.line.set_color("black")
ax.axes.yaxis.line.set_marker(' ')
ax.axes.zaxis.line.set_color("black")
ax.axes.zaxis.line.set_marker(' ')
loc = plticker.MultipleLocator(base = 0.3)
ax.axes.yaxis.set_major_locator(loc)
ax.dist = 11
plt.tight_layout()
fig.savefig('/home/cherb/Desktop/LRZ Sync+Share/Doktorarbeit/Mieze/Mieze_Foc/images/paper/3D_firstnight.pdf')
import ill_datafile.ill_datafile