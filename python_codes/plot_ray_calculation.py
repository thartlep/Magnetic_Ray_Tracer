import numpy as np
import matplotlib.pylab as plt
import sys
import ray_calculation as ray

# load raypaths
print 'Loading raypaths ... '
results = []
file_names = sys.argv[1].split()
for file_name in file_names:
  results_tmp = np.load(file_name)['arr_0']
  for i in range(0,results_tmp.shape[0]):
     r0 = np.copy(results_tmp[i,0])
     k0 = np.copy(results_tmp[i,1])
     Phi = np.copy(results_tmp[i,2])
     reached = results_tmp[i,3]
     reached_at_index = results_tmp[i,4]
     tau = results_tmp[i,5]
     results.append([results_tmp[i,0],results_tmp[i,1],results_tmp[i,2],results_tmp[i,3],results_tmp[i,4],results_tmp[i,5]])
#     results.append([r0,k0,Phi,reached,reached_at_index,tau])
     # add mirror version 
     if r0[1] == 0 and abs(k0[1]) > 1e-5:
        k0[1] = -k0[1]
        Phi[:,1] = -Phi[:,1]
        Phi[:,4] = -Phi[:,4]
        results.append([r0,k0,Phi,reached,reached_at_index,tau])

# load model and compute tau = 0.01 iso-surface
#print 'Loading model and computing tau=0.01 iso-surface ... '
#ray.load_model('Rempel_sunspot_data',False,False,False)
#print '.'
#xgrid = np.linspace(-20.0,20.0,1000)
#ygrid = np.array([0.0])
#tau_iso_value = 0.01
#tau_iso_surface = ray.compute_tau_iso_surface(xgrid,ygrid,tau_iso_value)
#print '.'

if len(sys.argv) > 2:
  plot_which = int(sys.argv[2])
else:
  plot_which = 0


##################################################
if plot_which == 0:

 plot_in_3d = False
 if len(sys.argv) > 3:
    if sys.argv[3] == 'True' or sys.argv[3] == 'true':
       plot_in_3d = True

 fig = plt.figure()

 if plot_in_3d:
    from mpl_toolkits.mplot3d import Axes3D
    ax = fig.add_subplot(111, projection='3d')
 else:
    ax = fig.add_subplot(111)

 print 'Plotting raypaths ... '
 for ray_path in results:

    r0 = ray_path[0]
    k0 = ray_path[1]
    Phi = ray_path[2]
    reached = ray_path[3]
    reached_at_index = ray_path[4]

    if np.min(Phi[:,2]) >= -10.0*1e8 and np.max(Phi[:,2]) < 1e8:
      print '   Plotting ray ... '
      print '      Starting point ',r0/1e8
      print '      Starting direction ',k0   
      if reached:
#         x_vals = Phi[0:reached_at_index+1,0]/1e8
#         y_vals = Phi[0:reached_at_index+1,1]/1e8
#         z_vals = Phi[0:reached_at_index+1,2]/1e8 
         x_vals = Phi[:,0]/1e8
         y_vals = Phi[:,1]/1e8
         z_vals = Phi[:,2]/1e8
         if plot_in_3d:
            ax.plot3D([x_vals[reached_at_index]],[y_vals[reached_at_index]],[z_vals[reached_at_index]],'.',color='black')
         else:
            ax.plot([x_vals[reached_at_index]],[z_vals[reached_at_index]],'.',color='black')
         print '      Reached tau = 0.01 after ',Phi[reached_at_index,6]/60.0,' mins at point ',np.array([Phi[reached_at_index,0],Phi[reached_at_index,1],Phi[reached_at_index,2]])/1e8 
      else:
         x_vals = Phi[:,0]/1e8
         y_vals = Phi[:,1]/1e8
         z_vals = Phi[:,2]/1e8
         print '      Did not reach tau = 0.01'
      if plot_in_3d:
         ax.plot3D(x_vals,y_vals,z_vals)
      else:
         ax.plot(x_vals,z_vals)

 if plot_in_3d:
   ax.set_xlabel('x [Mm]')
   ax.set_ylabel('y [Mm]')
   ax.set_zlabel('z [Mm]')
 else:
   ax.set_xlabel('x [Mm]')
   ax.set_ylabel('z [Mm]')

# ax.set_xlim([-25.0,+25.0])
# ax.set_ylim([-25.0,+25.0])
# ax.set_zlim([-10.0,+1.0])

 plt.show()


##################################################
if plot_which == 1:
 print 'Plotting "wave front" ... '

 time_step = 1.0    # mins
 max_time = time_step + 10.0
 for time_index in range(0,int(max_time/time_step),1):

   time_in_min = time_index * time_step
   time_in_sec = time_in_min * 60.0

   x_vals = []
   y_vals = []
   z_vals = []

   for ray_path in results:

      r0 = ray_path[0]
      k0 = ray_path[1]
      Phi = ray_path[2]

      if k0[1] == 0:
         closest_time_index = 0
         smallest_difference = 1e99
         for i in range(0,len(Phi[:,0])):
            if abs(Phi[i,6] - time_in_sec) < smallest_difference:
               closest_time_index = i
               smallest_difference = abs(Phi[i,6] - time_in_sec)
 
         if closest_time_index > 0:     
            x_vals.append(Phi[closest_time_index,0]/1e8)
            y_vals.append(Phi[closest_time_index,1]/1e8)
            z_vals.append(Phi[closest_time_index,2]/1e8)

   if len(x_vals) > 0:
      x_vals.append(x_vals[0])
      y_vals.append(y_vals[0])
      z_vals.append(z_vals[0])
      print '   Plotting "wave front" at t = ',time_in_min,' mins'
      plt.plot(x_vals,z_vals)
      plt.xlabel('x [Mm]')
      plt.ylabel('z [Mm]')
#      plt.title('Time = '+str(time_in_min))
      plt.xlim([-10.0,+15.0])
      plt.ylim([-5.0,+1.0])
 plt.show()
 

##################################################
if plot_which == 2:
   print 'Plotting "wave front" at surface ... '

   x_vals = []
   y_vals = []
   z_vals = []
   t_vals = []

   for ray_path in results:

      r0 = ray_path[0]
      k0 = ray_path[1]
      Phi = ray_path[2]
      reached = ray_path[3]
      reached_at_index = ray_path[4]
      tau = ray_path[5]
 
      if k0[1] == 0:
         if reached:
            xval = Phi[reached_at_index,0]/1e8
            yval = Phi[reached_at_index,1]/1e8
            zval = Phi[reached_at_index,2]/1e8
            tval = Phi[reached_at_index,6]/60.
            x_vals.append(xval)
            y_vals.append(yval)
            z_vals.append(zval)
            t_vals.append(tval)
 
   if len(t_vals) > 0:
      plt.plot(t_vals,x_vals,'.')
      plt.xlabel('t [min]')
      plt.ylabel('x [Mm]')
      plt.xlim([0.0,25.0])
      plt.ylim([-25.0,+25.0])
      plt.show()


##################################################
if plot_which == 3:
   print 'Plotting "wave front" at tau=0.01 surface ... '

   x_vals = []
   y_vals = []
   z_vals = []
   t_vals = []

   for ray_path in results:
#    if ray_path[1][1] == 0.0:
      r0 = ray_path[0]
      k0 = ray_path[1]
      Phi = ray_path[2]
      reached = ray_path[3]
      reached_at_index = ray_path[4]
      tau = ray_path[5]

      if reached:
         xval = Phi[reached_at_index,0]/1e8
         yval = Phi[reached_at_index,1]/1e8
         zval = Phi[reached_at_index,2]/1e8
         tval = Phi[reached_at_index,6]/60.
         x_vals.append(xval)
         y_vals.append(yval)
         z_vals.append(zval)
         t_vals.append(tval)

   time_step = 0.25
   times = np.arange(0.,19.,time_step)
   for time in times:
      out_x_vals = []
      out_y_vals = []
      out_z_vals = []
      for ray_index in range(0,len(x_vals)):
         if abs(t_vals[ray_index] - time) < time_step:
            out_x_vals.append(x_vals[ray_index])
            out_y_vals.append(y_vals[ray_index])
            out_z_vals.append(z_vals[ray_index])

      if len(out_x_vals) > 0:
         plt.plot(out_x_vals,out_y_vals,'.')
         plt.xlabel('x [Mm]')
         plt.ylabel('y [Mm]')
         plt.title('Time = '+str(time)+' min')
         plt.xlim([-20,+20])
         plt.ylim([-20,+20])
         plt.show()

