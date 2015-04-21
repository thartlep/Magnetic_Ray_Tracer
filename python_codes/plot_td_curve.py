import numpy as np
import matplotlib.cm as cm
import matplotlib.pylab as plt
from matplotlib.ticker import FuncFormatter
from scipy.interpolate import interp1d
import sys



### SET FIGURE SIZE ###############
from pylab import rcParams
rcParams['figure.figsize'] = 8, 4     # 8,6
### END - SET FIGURE SIZE #########

######################################
def gabor_wavelet(x_in_min,frequency):
   return np.cos(2*np.pi*frequency*x_in_min*60.0)*np.exp(-(x_in_min*60.0)**2/(2*(60.*5.)**2))
######################################

filename = 'td_combined_smoothed.txt'
f = open(filename,'r')
lines = f.readlines()
f.close()
distance = []
time = []
for i in range(0,len(lines)):
  splitline = lines[i].split()
  distance.append(float(splitline[0]))
  time.append(float(splitline[1]))

shift_distance_in_Mm = float(sys.argv[1])

distance_range = np.linspace(0.0, 40.0, 200)
time_range = np.linspace(0.0, 30.0, 100)
Y, X = np.meshgrid(time_range,distance_range)
Z = X * 0.0
print Z.shape

frequency = 3e-3
for i in range(0,len(distance_range)):
  d = distance_range[i] + shift_distance_in_Mm
  # get time corresponding to distance d from (distance,time)
  f = interp1d(distance, time, kind='cubic')
  central_time = f(d)
  if i == 0:
     central_time_at_0 = central_time
#  print 'Distance ',d,', time ',central_time
  # compute gabor wavelet 
  Z[i,:] = gabor_wavelet(time_range-central_time+central_time_at_0,frequency)
#  plt.plot(time_range,Z[i,:])
#  plt.show()
 

fig, ax = plt.subplots()
#im = plt.pcolor(X,Y,Z,cmap=cm.gray) 
im = ax.imshow(np.transpose(Z), cmap=cm.gray, vmin=-1, vmax=1, extent=[0., 40., 0., 30.], origin = 'lower',aspect='auto')
im.set_interpolation('bilinear')
plt.ylabel('Time [min]')
plt.xlabel(r'$\Delta$'+' [Mm]')
plt.xlim([0,35])
plt.ylim([0,25])
#plt.xlim([0,40])
#plt.ylim([0,30])
plt.plot([0,40],[0,40/45.3*1000./60.],':',color='white')
plt.savefig('td.eps',bbox_inches='tight')
#plt.show()

