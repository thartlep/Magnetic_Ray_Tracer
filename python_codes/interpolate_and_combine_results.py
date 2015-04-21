import sys
import numpy as np
import glob

##################################################
def interpolate_raypath(Phi,tau,time_step,absolute_max_time):

 import numpy as np
 from scipy.interpolate import interp1d
 import matplotlib
 matplotlib.use('Agg')
 import pylab as plt

 old_time = Phi[:,6]
 max_time = np.max(old_time)
 if max_time > absolute_max_time:
    max_time = absolute_max_time
 new_time = np.linspace(0.0,max_time,max_time/time_step)

 index_just_past_max = 0
 while Phi[index_just_past_max,6] < max_time and index_just_past_max < len(Phi[:,6])-1:
   index_just_past_max += 1

 old_xvals = Phi[:,0]
 old_yvals = Phi[:,1]
 old_zvals = Phi[:,2]

 old_kxvals = Phi[:,3]
 old_kyvals = Phi[:,4]
 old_kzvals = Phi[:,5]

 fxvals = interp1d(old_time[0:index_just_past_max+1], old_xvals[0:index_just_past_max+1])
 fyvals = interp1d(old_time[0:index_just_past_max+1], old_yvals[0:index_just_past_max+1])
 fzvals = interp1d(old_time[0:index_just_past_max+1], old_zvals[0:index_just_past_max+1])
 fkxvals = interp1d(old_time[0:index_just_past_max+1], old_kxvals[0:index_just_past_max+1])
 fkyvals = interp1d(old_time[0:index_just_past_max+1], old_kyvals[0:index_just_past_max+1])
 fkzvals = interp1d(old_time[0:index_just_past_max+1], old_kzvals[0:index_just_past_max+1])
 ftauvals = interp1d(old_time[0:index_just_past_max+1], tau[0:index_just_past_max+1])

 new_xvals = fxvals(new_time)
 new_yvals = fyvals(new_time)
 new_zvals = fzvals(new_time)
 new_kxvals = fkxvals(new_time)
 new_kyvals = fkyvals(new_time)
 new_kzvals = fkzvals(new_time)
 new_tau = ftauvals(new_time)

 new_Phi = np.zeros([len(new_time),7])
 new_Phi[:,0] = new_xvals
 new_Phi[:,1] = new_yvals
 new_Phi[:,2] = new_zvals
 new_Phi[:,3] = new_kxvals
 new_Phi[:,4] = new_kyvals
 new_Phi[:,5] = new_kzvals
 new_Phi[:,6] = new_time

 return new_Phi, new_tau
###########################################################


# times in seconds here
interpolate_time_step = 60.0 * 0.05
absolute_max_time = 60.0 * 30.0

flip_sign = True

#
files = []
for i in range(1,len(sys.argv)):
  more_files = glob.glob(sys.argv[i].split()[0])
  for another_file in more_files:
    files.append(another_file)

number_of_files = len(files)
results = []
for fileindex in range(0,number_of_files):
  file = files[fileindex]
  print 'Loading ',file,' ... '
  tmp = np.load(file)['arr_0']
  number_of_raypaths = tmp.shape[0]
  for i in range(0,tmp.shape[0]):
    r0 = tmp[i,0]
    k0 = tmp[i,1]
    if not flip_sign:
       Phi = tmp[i,2]
    else:
       Phi = tmp[i,2]*[-1,-1,1,-1,-1,1,1]
    reached = tmp[i,3]
    reached_at_index = tmp[i,4]
    tau = np.array(tmp[i,5])

    if reached:
       reached_time = Phi[reached_at_index,6]
       if reached_time > absolute_max_time:
          reached = False
    Phi, tau = interpolate_raypath(Phi,tau, interpolate_time_step,absolute_max_time)
    if reached:
       reached_at_index = int(reached_time/interpolate_time_step)
       if reached_time > absolute_max_time:
          reached_at_index = len(Phi[:,6])-1
    else:
       reached_at_index = 0

    results.append([r0,k0,Phi,reached,reached_at_index,tau])

## average rapyaths
#results2 = []
#for i in range(0,number_of_raypaths):
#  r0 = results[i][0]
#  k0 = results[i][1]
#  Phi = results[i][2]
#  reached = int(results[i][3])
#  reached_at_index = float(results[i][4])
#  tau = results[i][5]
#
#  for fileindex in range(1,number_of_files):
#     rayindex = i + fileindex*number_of_raypaths
#     r0 += results[rayindex][0]
#     k0 += results[rayindex][1]
#     Phi += results[rayindex][2]
#     reached += results[rayindex][3]
#     reached_at_index += results[rayindex][4]
#     tau += results[rayindex][5]
#
#  norm = 1.0 / number_of_files
#  r0 = r0 * norm
#  k0 = k0 * norm
#  Phi = Phi * norm
#  reached = bool( reached * norm )
#  reached_at_index = int( reached_at_index * norm )
#  tau = tau * norm
#  results2.append([r0,k0,Phi,reached,reached_at_index,tau])
#del results
#results = results2

np.savez('output.interpolated_and_combined.npz',results)


