import numpy as np
import matplotlib.pylab as plt
import sys
import ray_calculation as ray
import numpy as np
from scipy.interpolate import interp1d

interpolate_the_wave_fronts = False
interpolate_the_tau_surface = False

#################################################
def get_angle_from_k(k):
   kx = k[0]
   ky = k[1]
   kz = k[2]
   rad_angle_h = np.arcsin(kz)
   angle_h = rad_angle_h * 360.0 / ( 2 * np.pi )
   if kx >= 0:
      angle_h = angle_h
   else:
      angle_h = -angle_h + 180
   return angle_h
##################################################
def interpolate_raypath(Phi,time_step):

 import numpy as np
 from scipy.interpolate import interp1d
 import matplotlib.pylab as plt

 old_time = Phi[:,6]
 max_time = np.max(old_time)
 new_time = np.linspace(0.0,max_time,max_time/time_step)

 old_xvals = Phi[:,0]
 old_yvals = Phi[:,1]
 old_zvals = Phi[:,2]

 old_kxvals = Phi[:,3]
 old_kyvals = Phi[:,4]
 old_kzvals = Phi[:,5]

 fxvals = interp1d(old_time, old_xvals)
 fyvals = interp1d(old_time, old_yvals)
 fzvals = interp1d(old_time, old_zvals)
 fkxvals = interp1d(old_time, old_kxvals)
 fkyvals = interp1d(old_time, old_kyvals)
 fkzvals = interp1d(old_time, old_kzvals)

 new_xvals = fxvals(new_time)
 new_yvals = fyvals(new_time)
 new_zvals = fzvals(new_time)
 new_kxvals = fkxvals(new_time)
 new_kyvals = fkyvals(new_time)
 new_kzvals = fkzvals(new_time)

# plt.plot(new_time,new_xvals,'o')
# plt.plot(old_time,old_xvals,'-')
# plt.show()

 new_Phi = np.zeros([len(new_time),7])
 new_Phi[:,0] = new_xvals
 new_Phi[:,1] = new_yvals
 new_Phi[:,2] = new_zvals
 new_Phi[:,3] = new_kxvals
 new_Phi[:,4] = new_kyvals
 new_Phi[:,5] = new_kzvals
 new_Phi[:,6] = new_time

 return new_Phi
##################################################
def plot_raypaths(show,plot_in_3d,skip_raypaths,overplot,fig=0,ax=0):

 if not overplot:
    fig = plt.figure()
    if plot_in_3d:
       from mpl_toolkits.mplot3d import Axes3D
       ax = fig.add_subplot(111, projection='3d')
    else:
       ax = fig.add_subplot(111)

 tau_surface_x = []
 tau_surface_y = []
 tau_surface_z = []

 print 'Plotting raypaths ... '
 print 'Number of raypaths: ',len(results)
 for ray_path in results:

    r0 = ray_path[0]
    k0 = ray_path[1]
    Phi = ray_path[2]
    reached = ray_path[3]
    reached_at_index = ray_path[4]

    if plot_in_3d or Phi[1,1] == 0.0:
      print '   Plotting ray ... '
      print '      Starting point ',r0/1e8
      print '      Starting direction ',k0
      if reached:
         x_vals = Phi[0:reached_at_index+1,0]/1e8
         y_vals = Phi[0:reached_at_index+1,1]/1e8
         z_vals = Phi[0:reached_at_index+1,2]/1e8
#         x_vals = Phi[:,0]/1e8
#         y_vals = Phi[:,1]/1e8
#         z_vals = Phi[:,2]/1e8
         tau_surface_x.append(x_vals[reached_at_index])
         tau_surface_y.append(y_vals[reached_at_index])
         tau_surface_z.append(z_vals[reached_at_index])
         if not skip_raypaths:
            if plot_in_3d:
               ax.plot3D([x_vals[reached_at_index]],[y_vals[reached_at_index]],[z_vals[reached_at_index]],'.',color='black')
#            else:
#               ax.plot([x_vals[reached_at_index]],[z_vals[reached_at_index]],'.',color='black')
         print '      Reached tau = 0.01 after ',Phi[reached_at_index,6]/60.0,' mins at point ',np.array([Phi[reached_at_index,0],Phi[reached_at_index,1],Phi[reached_at_index,2]])/1e8
      else:
         x_vals = Phi[:,0]/1e8
         y_vals = Phi[:,1]/1e8
         z_vals = Phi[:,2]/1e8
         print '      Did not reach tau = 0.01'
      if not skip_raypaths:
         if plot_in_3d:
            ax.plot3D(x_vals,y_vals,z_vals)
         else:
            ax.plot(x_vals,z_vals,color='black')

 if plot_in_3d:
   ax.set_xlabel('x [Mm]')
   ax.set_ylabel('y [Mm]')
   ax.set_zlabel('z [Mm]')
 else:
   ax.set_xlabel('x [Mm]')
   ax.set_ylabel('z [Mm]')

 if not plot_in_3d:
    plt.xlim((-5.0,+25.0))
    plt.ylim((-6.0,+1.0))

    # plot tau surface
    if interpolate_the_tau_surface:
       fx = interp1d(np.arange(0,len(tau_surface_x)),tau_surface_x,'slinear')
       fz = interp1d(np.arange(0,len(tau_surface_z)),tau_surface_z,'slinear')
       smooth_tau_surface_x = fx(np.linspace(0,len(tau_surface_x)-1,len(tau_surface_x)*10))
       smooth_tau_surface_z = fz(np.linspace(0,len(tau_surface_z)-1,len(tau_surface_z)*10))
       tau_surface_x = smooth_tau_surface_x
       tau_surface_z = smooth_tau_surface_z
    tau_surface_x, tau_surface_z = (list(x) for x in zip(*sorted(zip(list(tau_surface_x), list(tau_surface_z)))))
    plt.plot(tau_surface_x,tau_surface_z,'-',color='black',linewidth=5)

 if show:
    plt.show()

 if not overplot:
    return fig,ax
##################################################


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
   plot_raypaths(True,False,False,False)

if plot_which == 6:
   plot_raypaths(True,True,False,False)

##################################################
if plot_which == 1 or plot_which == 4 or plot_which == 5:
 print 'Plotting "wave front" ... '

 if plot_which == 4:
    fig,ax = plot_raypaths(False,False,False,False)

 file_counter = 0

 mean_time_step = 0.05        # mins
 time_step = 0.05             # mins

 frequency = 3e-3
 packet_half_time_in_min = 9.0

 wave_packet_times_in_sec = np.arange(-packet_half_time_in_min*60,+packet_half_time_in_min*60,time_step*60)
 gabor_wavelet = np.cos(2*np.pi*frequency*wave_packet_times_in_sec)*np.exp(-np.pi*wave_packet_times_in_sec**2*5e-6)
 gabor_wavelet = gabor_wavelet / np.max(np.abs(gabor_wavelet))

 min_time = -packet_half_time_in_min
 max_time = time_step + 20.0 + packet_half_time_in_min * 2

 for mean_time_index in range(int(min_time/mean_time_step),int(max_time/mean_time_step),1):
   
    if plot_which == 5:
       fig,ax = plot_raypaths(False,False,True,False)
       time_string = "{:5.2f}".format(mean_time_index*mean_time_step)
       plt.text(18.,0.5,'Time = '+time_string+' min')

    for wavelet_index in range(0,len(gabor_wavelet)):

       time_index = wavelet_index - len(gabor_wavelet)/2 + mean_time_index*mean_time_step/time_step

       if time_index >= 0 and time_index < int(max_time/time_step): 
 
          time_in_seconds = time_index*time_step*60.0
          amplitude = gabor_wavelet[wavelet_index]
          if amplitude > 0:
             green = 1.0 - amplitude
             red = 1.0
             blue = 1.0 - amplitude
          else:
             green = 1.0 -abs(amplitude)
             red = 1.0 - abs(amplitude)
             blue = 1.0
          color = [red, green, blue]


          time_in_min = time_index * time_step
          time_in_sec = time_in_min * 60.0

          x_vals = []
          y_vals = []
          z_vals = []
          angle_h_vals = []
          o_vals = []

          for ray_path in results:

             r0 = ray_path[0]
             k0 = ray_path[1]
             Phi = ray_path[2]
             reached = ray_path[3]
             reached_at_index = ray_path[4]
             if reached:
                reached_at_time = Phi[reached_at_index,6]

             if k0[1] == 0:

                closest_time_index = 0
                smallest_difference = 1e99
                for i in range(0,len(Phi[:,0])):
                   if abs(Phi[i,6] - time_in_sec) < smallest_difference:
                      closest_time_index = i
                      smallest_difference = abs(Phi[i,6] - time_in_sec)
 
                if closest_time_index > 0:     
                   if (not reached) or (reached and reached_at_index >= closest_time_index): 
                      x_value = Phi[closest_time_index,0]/1e8
                      y_value = Phi[closest_time_index,1]/1e8
                      z_value = Phi[closest_time_index,2]/1e8
                   else:
                      x_value = float('nan')
                      y_value = float('nan')
                      z_value = float('nan')

                   x_vals.append(x_value)
                   y_vals.append(y_value)
                   z_vals.append(z_value)
                   angle_h = get_angle_from_k(k0)
                   angle_h_vals.append(angle_h)

          if len(x_vals) > 0:
#             print 'BEFORE ORDERING'
#             for i in range(0,len(x_vals)):
#                print angle_h_vals[i],x_vals[i],z_vals[i]
             print '   Plotting "wave front" at t = ',time_in_min,' mins'
             # first, sort them by angle_h
             angle_h_vals,x_vals,z_vals = (list(x) for x in zip(*sorted(zip(list(angle_h_vals), list(x_vals), list(z_vals)))))
#             print 'AFTER ORDERING'
#             for i in range(0,len(x_vals)):
#                print angle_h_vals[i],x_vals[i],z_vals[i]
             # now plot (interpolated or not)
             if interpolate_the_wave_fronts:
                tmp_x_vals = []
                tmp_z_vals = []
                for i in range(0,len(x_vals)):
                   if x_vals[i] < 0 or x_vals[i] >= 0:
                      tmp_x_vals.append(x_vals[i])
                      tmp_z_vals.append(z_vals[i])
                   if not (x_vals[i] < 0 or x_vals[i] >= 0) or i == len(x_vals) - 1:
                      if len(tmp_x_vals) > 1:
                         if len(tmp_x_vals) >= 2:
                            method = 'slinear'
#                         if len(tmp_x_vals) >= 3:
#                            method = 'quadratic'
                         fx = interp1d(np.arange(0,len(tmp_x_vals)),tmp_x_vals,method)
                         fz = interp1d(np.arange(0,len(tmp_z_vals)),tmp_z_vals,method)
                         smooth_x = fx(np.linspace(0,len(tmp_x_vals)-1,len(tmp_x_vals)*10))
                         smooth_z = fz(np.linspace(0,len(tmp_x_vals)-1,len(tmp_x_vals)*10))
                         plt.plot(smooth_x,smooth_z,':',color=color,linewidth=5)
                      tmp_x_vals = []
                      tmp_y_vals = []
             else:
                plt.plot(x_vals,z_vals,'-',color=color,linewidth=5)


    if plot_which == 5:
       plot_raypaths(False,False,True,True,fig,ax) 
    plt.xlabel('x [Mm]')
    plt.ylabel('z [Mm]')
#    plt.title('Time = '+str(time_in_min))
    plt.xlim([0.0,+25.0])
    plt.ylim([-6.0,+1.0])

    if plot_which == 5:
       if file_counter < 10:
          filename = '00'+str(file_counter)+'.png'
       else:
          if file_counter < 100:
             filename = '0'+str(file_counter)+'.png'
          else:
             filename = str(file_counter)+'.png'
       plt.savefig(filename)
       file_counter += 1

 if plot_which == 4:
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
      plt.xlim([0.0,20.0])
      plt.ylim([-20.0,+30.0])
      plt.show()

   for i in range(0,len(t_vals)):
      print t_vals[i],x_vals[i]   


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


