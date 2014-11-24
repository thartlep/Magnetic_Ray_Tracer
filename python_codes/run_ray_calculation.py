import numpy as np
import matplotlib.pylab as plt
import ray_calculation as ray
import sys
import os, errno, glob

#################################################
def mkdir(dir):
  try:
      os.makedirs(dir)
  except OSError as exc: # Python >2.5
      if exc.errno == errno.EEXIST and os.path.isdir(dir):
          pass
      else: raise
#################################################
def direction(angle_h, angle_xy):
   rad_angle_h = angle_h * 2 * np.pi / 360.
   rad_angle_xy = angle_xy * 2 * np.pi / 360.
   return np.array([ np.cos(rad_angle_h)*np.cos(rad_angle_xy), np.cos(rad_angle_h)*np.sin(rad_angle_xy), np.sin(rad_angle_h) ])
#################################################


dir_simulation_data = sys.argv[1]
case_name = sys.argv[2]
file_index = int(sys.argv[3])

print 'CASE: '+case_name+', DATA: '+dir_simulation_data+'\n'

mkdir('../results/'+case_name)

set_magnetic_field_to_zero = False
set_N_square_to_zero = True
set_omegac_square_to_zero = True

set_to_reach_tau_001 = True       # stop computing when ray reaches tau=0.01
mintime = 15.0                    # minimum travel time to compute
maxtime = 60.0                    # maximum travel time to compute

accuracy = 1e-3

n_S = 10000
omega = 3e-3
additional_argument = (omega,)

depth = 5.0
offset = 7.5

if file_index == 0:
   angle_resolution = 10
   angle_xy_s = range(0,90+angle_resolution,angle_resolution)
   angle_h_s = range(-60,240+angle_resolution,angle_resolution)
if file_index == 1:
   angle_resolution = 10
   angle_xy_s = 5+range(0,80+angle_resolution,angle_resolution)
   angle_h_s = range(-60,240+angle_resolution,angle_resolution)
if file_index == 2:
   angle_xy_s = [0]
#   angle_h_s = [10,-5,-3,-2,-1.5,-1,-0.5,0.,0.5,1,1.5,2,3,5,10,30,60,90,120,150,170,175,177,178,178.5,179,179.5,180,180.5,181,181.5,182,183,185,190,225,270,315]
   angle_h_s = [25,15,5,-3.3,-6.7,-13.3,-16.7,-23.3,-26.7,-32.5,-35,-37.5,-42.5,-45.,-85]
   angle_h_s = angle_h_s + list(180.0 - np.array(angle_h_s))   

# run name
run_name = 'freq_'+str(omega*1e3)+'_mHz__depth_'+str(depth)+'_Mm__offset_'+str(offset)+'_Mm__run'+str(file_index)
print 'Run_name: ',run_name

# initial locations
r0s =           [
                   np.array([   offset,   0.0,   -depth ])*1e8
                ]


# initial directions
k0_directions = []
for angle_h in angle_h_s:
 for angle_xy in angle_xy_s:
   k0_directions.append( direction(angle_h,angle_xy) )

results = ray.propagate_rays(dir_simulation_data,set_magnetic_field_to_zero,set_N_square_to_zero,set_omegac_square_to_zero,r0s,k0_directions,omega,set_to_reach_tau_001,mintime,maxtime,n_S,accuracy)

print 'Saving results ... '
np.savez('results/'+case_name+'/output.raypaths.'+run_name+'.npz',results)

print 'Done.'
