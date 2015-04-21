import numpy as np
import sys
from scipy import interpolate


td = []
for i in range(0,6):
  filename = 'td'+str(i)+'.txt'
  print 'Reading '+filename+' ... '
  f = open(filename,'r')
  lines= f.readlines()
  for line in lines:
    td.append([float(line.split()[0]),float(line.split()[1])])

def getkey(item):
  return item[1]

td_sorted = sorted(td,key=getkey)
t = []
d = []
for item in td_sorted:
 if item[1] > 0.0:
  t.append(item[0])
  d.append(item[1])

f = open('td_combined.txt','w')
for i in range(0,len(t)):
  f.write(str(d[i])+'   '+str(t[i])+'\n')
f.close()



filename = 'td_combined.txt'
f = open(filename,'r')
lines = f.readlines()
f.close()
input_distance = []
input_time = []
for line in lines:
  splitline = line.split()
  input_distance.append(float(splitline[0]))
  input_time.append(float(splitline[1]))

from scipy.interpolate import UnivariateSpline
spl = UnivariateSpline(input_distance, input_time)
spl.set_smoothing_factor(250.)
output_distance = np.arange(0,50.,0.1)
output_time = spl(output_distance)

filename = 'td_combined_smoothed.txt'
f = open(filename,'w')
for i in range(0,len(output_time)):
   f.write(str(output_distance[i])+'   '+str(output_time[i])+'\n')

f.close()

