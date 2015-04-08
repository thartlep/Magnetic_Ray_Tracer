import matplotlib
matplotlib.use('Agg')
import pylab as plt

import numpy as np
import sys
import ray_calculation as ray
import numpy as np
from scipy.interpolate import interp1d
import glob

interpolate_the_wave_fronts = True
interpolate_the_tau_surface = True

#################################################
# interpolate
def my_interpolate(x,y,smoothing_factor):

   k = 3
   if len(x) > k:
     from scipy.interpolate import UnivariateSpline
     spl = UnivariateSpline(x, y, k=k)
     spl.set_smoothing_factor(smoothing_factor)
     new_x = np.linspace(np.min(x),np.max(x),200)
     new_y = spl(new_x)
   else:
     new_y = y

   return new_y

#################################################
# 'stolen' from the internets @ http://stackoverflow.com/questions/22988882/how-to-smooth-a-curve-in-python
#
def savitzky_golay(y, window_size, order, deriv=0, rate=1):

    import numpy as np
    from math import factorial

    try:
        window_size = np.abs(np.int(window_size))
        order = np.abs(np.int(order))
    except ValueError, msg:
        raise ValueError("window_size and order have to be of type int")
    if window_size % 2 != 1 or window_size < 1:
        raise TypeError("window_size size must be a positive odd number")
    if window_size < order + 2:
        raise TypeError("window_size is too small for the polynomials order")
    order_range = range(order+1)
    half_window = (window_size -1) // 2
    # precompute coefficients
    b = np.mat([[k**i for i in order_range] for k in range(-half_window, half_window+1)])
    m = np.linalg.pinv(b).A[deriv] * rate**deriv * factorial(deriv)
    # pad the signal at the extremes with
    # values taken from the signal itself
    firstvals = y[0] - np.abs( y[1:half_window+1][::-1] - y[0] )
    lastvals = y[-1] + np.abs(y[-half_window-1:-1][::-1] - y[-1])
    y = np.concatenate((firstvals, y, lastvals))
    return np.convolve( m[::-1], y, mode='valid')
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
 import matplotlib
 matplotlib.use('Agg')
 import pylab as plt

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
    tau_surface_x, tau_surface_z = (list(x) for x in zip(*sorted(zip(list(tau_surface_x), list(tau_surface_z)))))
    if interpolate_the_tau_surface:
      tmp_x = []
      tmp_z = []
      for i in range(0,len(tau_surface_x)):
        if i % 10 == 0:
          tmp_x.append(tau_surface_x[i])
          tmp_z.append(tau_surface_z[i])
      tau_surface_x = tmp_x
      tau_surface_z = tmp_z
      tau_surface_x = savitzky_golay(tau_surface_x, 11, 3)
      tau_surface_z = savitzky_golay(tau_surface_z, 11, 3)
    for add_z in np.arange(1.5,0,-0.1):
      plt.plot(np.array(tau_surface_x),np.array(tau_surface_z)+add_z,'-',color='white',linewidth=10)
    plt.plot(np.array(tau_surface_x),np.array(tau_surface_z)-0.08,'-',color='black',linewidth=8)

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
     tau = list(results_tmp[i,5])
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
 time_step = 2.0              # mins

 frequency = 3e-3
 packet_half_time_in_min = 20.0

 min_time = -packet_half_time_in_min
 max_time = time_step + 20.0 + packet_half_time_in_min * 2

 if len(sys.argv) > 3:
    first_image = int(sys.argv[3])
 else:
    first_image = 0
 if len(sys.argv) > 4:
    last_image = int(sys.argv[4])
 else:
    last_image = 10000000


 for mean_time_index in range(int(min_time/mean_time_step),int(max_time/mean_time_step),1):

    skip_computation = False
    if plot_which == 5:
       if file_counter < first_image or file_counter > last_image:
          skip_computation = True
       if file_counter < 10:
          filename = '00'+str(file_counter)+'.eps'
       else:
          if file_counter < 100:
             filename = '0'+str(file_counter)+'.eps'
          else:
             filename = str(file_counter)+'.eps'
       existing_files = glob.glob(filename)
       if len(existing_files) != 0:
          skip_computation = True

    file_counter += 1
  
    if mean_time_index*mean_time_step != 7.5:
      skip_computation = True

    if skip_computation:

     print ' ... skipping ',filename,' ... '

    else:

     if plot_which == 5:
        fig,ax = plot_raypaths(False,False,True,False)
        time_string = "{:5.2f}".format(mean_time_index*mean_time_step)
        plt.text(15.5,0.5,'Time = '+time_string+' min')

     for case in [0]: #,1]:

      if case == 0:
         line_thickness = 15
         time_step = 1.0 #0.25
      elif case == 1:
         line_thickness = 3
         time_step = 0.05


      print 'Running with line_thickness of ',line_thickness,' and a time_step of ',time_step,' ... '

      wave_packet_times_in_sec = np.arange(-packet_half_time_in_min*60,+packet_half_time_in_min*60,time_step*60)
      gabor_wavelet = np.cos(2*np.pi*frequency*wave_packet_times_in_sec)*np.exp(-wave_packet_times_in_sec**2/(2*300.**2))
      gabor_wavelet = gabor_wavelet / np.max(np.abs(gabor_wavelet))

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
#              print 'BEFORE ORDERING'
#              for i in range(0,len(x_vals)):
#                 print angle_h_vals[i],x_vals[i],z_vals[i]
              print '   Plotting "wave front" at t = ',time_in_min,' mins'
              # first, sort them by angle_h
              angle_h_vals,x_vals,z_vals = (list(x) for x in zip(*sorted(zip(list(angle_h_vals), list(x_vals), list(z_vals)))))
#              print 'AFTER ORDERING'
#              for i in range(0,len(x_vals)):
#                 print angle_h_vals[i],x_vals[i],z_vals[i]




              # remove possible NANs
              xvals = []
              zvals = []
              anglevals = []
              tmp_xvals = []
              tmp_zvals = []
              tmp_anglevals = []

              for i in range(0,len(x_vals)):
                 if x_vals[i] != x_vals[i] or z_vals[i] != z_vals[i]:
                    if len(tmp_zvals) > 0:
                       xvals.append(np.array(tmp_xvals))
                       zvals.append(np.array(tmp_zvals))
                       anglevals.append(np.array(tmp_anglevals))
                       tmp_xvals = []
                       tmp_zvals = []
                       tmp_anglevals = []
                 else:
                    tmp_xvals.append(x_vals[i])
                    tmp_zvals.append(z_vals[i])
                    tmp_anglevals.append(angle_h_vals[i])
              if len(tmp_zvals) > 0:
                 xvals.append(np.array(tmp_xvals))
                 zvals.append(np.array(tmp_zvals))
                 anglevals.append(np.array(tmp_anglevals))

              # now plot (interpolated or not)
              print 'Number of segments: ',len(xvals)
              for i in range(0,len(xvals)):
                if len(xvals[i]) > 10: 


                  new_xvals = []
                  new_zvals = []
                  new_anglevals = []
                  current_xval = 0
                  current_zval = 0
                  current_anglevals = 0
                  current_counter = 0
                  new_val = True
                  for j in range(0,len(anglevals[i])):
#                     print anglevals[i][j],new_anglevals
                     if not new_val:
                        if abs(anglevals[i][j] - (current_angleval/current_counter)) < 2.0:
#                           print 'Less ',current_angleval/current_counter
                           current_xval += xvals[i][j]
                           current_zval += zvals[i][j]
                           current_angleval += anglevals[i][j]
                           current_counter += 1
                        else:
#                           print 'Greater'
                           new_xvals.append(current_xval / current_counter)
                           new_zvals.append(current_zval / current_counter)
                           new_anglevals.append(current_angleval / current_counter)
                           current_xval = 0
                           current_zval = 0
                           current_anglevals = 0
                           current_counter = 0
                           new_val = True
                     if new_val:
#                        print 'New'
                        current_xval = xvals[i][j]
                        current_zval = zvals[i][j]
                        current_angleval = anglevals[i][j]
                        current_counter = 1
                        new_val = False

#                  print new_anglevals
                  plot_xvals = my_interpolate(new_anglevals,new_xvals,1)
#                  plot_xvals = my_interpolate(anglevals[i],xvals[i],2)
#                  plot_xvals = savitzky_golay(xvals[i], 17, 2)
                  plot_zvals = my_interpolate(new_anglevals,new_zvals,1)
#                  plot_zvals = my_interpolate(anglevals[i],zvals[i],2)
#                  plot_zvals = savitzky_golay(zvals[i], 17, 2)
                  plot_xvals = list(plot_xvals)
                  plot_zvals = list(plot_zvals)
                  if 1 == 0 and (plot_xvals[len(plot_xvals)-1] > 10 or plot_xvals[len(plot_xvals)-1] < 5) and plot_zvals[len(plot_zvals)-1] > -1:
                    slope = 0.0
                    xdiff = 0.0
                    n = min(2,len(plot_xvals)-3)
                    if n > 1:
                      counter = 0
                      for j in range(1,n):
                        this_xdiff = plot_xvals[len(plot_xvals)-j]-plot_xvals[len(plot_xvals)-j-1]
                        this_zdiff = plot_zvals[len(plot_zvals)-j]-plot_zvals[len(plot_zvals)-j-1]
                        xdiff += this_xdiff
                        slope += 0.4*this_zdiff/this_xdiff
                        counter += 1.0
                      slope = slope / counter
                      if (plot_xvals[len(plot_xvals)-1] > 10 and slope < 0.) or (plot_xvals[len(plot_xvals)-1] <5 and slope > 0.2):
                        print 'Slope:',slope
                        if plot_xvals[len(plot_xvals)-1] > 10:
                           step = -2.0
                        else:
                           step = 2.0
                        tmp_xvals = plot_xvals[0:len(plot_xvals)-n-1]
                        tmp_zvals = plot_zvals[0:len(plot_zvals)-n-1]
                        tmp_xvals.append(tmp_xvals[len(tmp_xvals)-1]+step)
                        tmp_zvals.append(tmp_zvals[len(tmp_zvals)-1]+slope*step)
                        plot_xvals = tmp_xvals
                        plot_zvals = tmp_zvals
                  plt.plot(plot_xvals,plot_zvals,'-',color=color,linewidth=line_thickness)

     #plt.show()
               
     if plot_which == 5:
        plot_raypaths(False,False,True,True,fig,ax) 
     plt.xlabel('x [Mm]')
     plt.ylabel('z [Mm]')
#     plt.title('Time = '+str(time_in_min))
     plt.xlim([0.0,+22.5])
     plt.ylim([-6.0,+1.0])

     if plot_which == 5:
        plt.savefig(filename)

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


