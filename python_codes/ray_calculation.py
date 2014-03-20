import scipy.integrate as int
import numpy as np
import matplotlib.pylab as plt

#########################################
def load_model(dir):

  from astropy.io import fits as pyfits

  fitsfile = pyfits.open(dir+'/r.fits')
  model_r = fitsfile[0].data
  fitsfile.close()

  fitsfile = pyfits.open(dir+'/z.fits')
  model_z = fitsfile[0].data
  fitsfile.close()

  fitsfile = pyfits.open(dir+'/br.fits')
  model_Br = fitsfile[0].data
  fitsfile.close()

  fitsfile = pyfits.open(dir+'/bz.fits')
  model_Bz = fitsfile[0].data
  fitsfile.close()

  fitsfile = pyfits.open(dir+'/vr.fits')
  model_vr = fitsfile[0].data
  fitsfile.close()

  fitsfile = pyfits.open(dir+'/vz.fits')
  model_vz = fitsfile[0].data
  fitsfile.close()

  fitsfile = pyfits.open(dir+'/T.fits')
  model_T = fitsfile[0].data
  fitsfile.close()

  fitsfile = pyfits.open(dir+'/p.fits')
  model_p = fitsfile[0].data
  fitsfile.close()

  fitsfile = pyfits.open(dir+'/rho.fits')
  model_rho = fitsfile[0].data
  fitsfile.close()

  fitsfile = pyfits.open(dir+'/tau.fits')
  model_tau = fitsfile[0].data
  fitsfile.close()

  fitsfile = pyfits.open(dir+'/gamma1.fits')
  model_gamma1 = fitsfile[0].data
  fitsfile.close()

  return model_r,model_z,model_vr,model_vz,model_Br,model_Bz,model_T,model_p,model_rho,model_tau,model_gamma1

#########################################
def func(Phi,s,omega):

  # input values
  r = np.array([Phi[0],Phi[1],Phi[2]])
  k = np.array([Phi[3],Phi[4],Phi[5]])
  t = Phi[6]
  k_square = k[0]**2 + k[1]**2 + k[2]**2

  # get model values
  B,cs,cs_square,va,va_square,N,N_square,omegac,omegac_square,dcsdxi,dvajdxi,dNdxi,domegacdxi = get_interpolated_model(r)


  # stuff
  delz = np.array([0.,0.,0.,0.,0.0,0.,0.,0.,1.]).reshape((3,3))
  delxy = np.array([1.,0.,0.,0.,1.0,0.,0.,0.,0.]).reshape((3,3))

  # func[0] = d F / d ki
  dFdki = (
            - omega**2*(cs_square+va_square)*2*k
            + cs_square * ( 3*va*k*k + va*(k_square - k*k) + 2*k*(np.dot(va,k)-va*k) )
            + omegac_square*cs_square*2*np.dot(delz,k)
            + cs_square * N_square * 2 * np.dot(delxy,k)
          )

  # -func[1] = d F / d xi
  dFdxi = (
            - omega**2 * k_square * (2*cs*dcsdxi + 2*np.dot(dvajdxi,va))
            + k_square * ( np.dot(k,va) * 2 * cs * dcsdxi + cs_square * np.dot(dvajdxi,k) )
            - (omega**2*2*omegac*domegacdxi -k[2]**2*(omegac_square*2*cs*dcsdxi+cs_square*2*omegac*domegacdxi) )
            + (k[0]**2+k[1]**2)*( N_square*2*cs*dcsdxi + cs_square*2*N*dNdxi ) 
          )

  # -func[2] = d F / d omega
  dFdomega = ( 
               + 4*omega**3
               - 2*omega*(cs_square + va_square)*k_square 
               - omegac_square * 2 * omega
             )

  # return result
  print s,'||',dFdki[0],dFdki[1],dFdki[2],'|',-dFdxi[0],-dFdxi[1],-dFdxi[2],'|',-dFdomega
  return [dFdki[0],dFdki[1],dFdki[2],-dFdxi[0],-dFdxi[1],-dFdxi[2],-dFdomega]

#########################################
def get_interpolated_model(r):

  x = r[0]
  y = r[1]
  z = r[2]

  B = np.array([0.,0.,0.])
  rho = 1.0
  H = 1e99
  cs = 1.0
  g = 0.0


  cs_square= cs**2

  va = B / np.sqrt(4*np.pi*rho)
  va_square = np.dot(va,va)

  omegac_square = cs / ( 2*H )
  omegac = np.sqrt(omegac_square)

  N_square = g/H - g**2/cs_square
  N = np.sqrt(N_square)

  dcsdxi = np.zeros(3)
  dvajdxi = np.zeros([3,3])
  dNdxi = np.zeros(3)
  domegacdxi = np.zeros(3)

  return B,cs,cs_square,va,va_square,N,N_square,omegac,omegac_square,dcsdxi,dvajdxi,dNdxi,domegacdxi 

#########################################
def make_cartesian_model(model_r,model_z,model_vr,model_vz,model_Br,model_Bz,model_T,model_p,model_rho,model_tau,model_gamma1):

  tmp = []
  for i in range(0,len(model_r)):
     tmp.append(-model_r[len(model_r)-1-i])
  for i in range(1,len(model_r)):
     tmp.append(model_r[i])
  cartesian_model_x = np.array(tmp)
  cartesian_model_y = np.array(tmp)
  cartesian_model_z = model_z

  xlen = len(cartesian_model_x)
  ylen = len(cartesian_model_y)
  zlen = len(cartesian_model_z)

  cartesian_model_va_x = np.zeros([xlen,ylen,zlen])
  cartesian_model_va_y = np.zeros([xlen,ylen,zlen])
  cartesian_model_va_z = np.zeros([xlen,ylen,zlen])
  cartesian_model_cs_square = np.zeros([xlen,ylen,zlen])
  cartesian_model_omegac_square = np.zeros([xlen,ylen,zlen])
  cartesian_model_N_square = np.zeros([xlen,ylen,zlen])
  cartesian_model_p = np.zeros([xlen,ylen,zlen])

  print '   Computing va and cs_square ... '
  max_r = np.max(model_r)
  from scipy import interpolate
  for k in range(0,zlen):
     print '      z='+str(k)+' of '+str(zlen)
     interpolator_Br = interpolate.interp1d(model_r[:,0],model_Br[k,:],kind='cubic')
     interpolator_Bz = interpolate.interp1d(model_r[:,0],model_Bz[k,:],kind='cubic')
     interpolator_rho = interpolate.interp1d(model_r[:,0],model_rho[k,:],kind='cubic')
     interpolator_p = interpolate.interp1d(model_r[:,0],model_p[k,:],kind='cubic')
     interpolator_gamma1 = interpolate.interp1d(model_r[:,0],model_gamma1[k,:],kind='cubic')
     for i in range(0,xlen):
        for j in range(0,ylen):
           r_value = np.sqrt(cartesian_model_x[i]**2 + cartesian_model_y[i]**2)
           if r_value > 0:
              if cartesian_model_x[i] != 0:
                 alpha = np.arccos(cartesian_model_x[i]/r_value)
              else:
                 alpha = np.arcsin(cartesian_model_y[i]/r_value)
              if r_value > max_r:
                 r_value = max_r
              Br = interpolator_Br(r_value)
           else:
              alpha = 0.0
              Br = 0.0
           Bz = interpolator_Bz(r_value)
           Bx = Br*np.cos(alpha)
           By = Br*np.sin(alpha)
           rho = interpolator_rho(r_value)
           p = interpolator_p(r_value)
           gamma1 = interpolator_gamma1(r_value)
           factor_va = 1.0/np.sqrt(4*np.pi*rho)
           va_x = Bx * factor_va
           va_y = By * factor_va
           va_z = Bz * factor_va
           cs_square = gamma1 * p / rho
           cartesian_model_va_x[i,j,k] = va_x
           cartesian_model_va_y[i,j,k] = va_y
           cartesian_model_va_z[i,j,k] = va_z
           cartesian_model_cs_square[i,j,k] = cs_square
           cartesian_model_p[i,j,k] = p

  print '   Computing omegac_square and N_square ... '
  for i in range(0,xlen):
     for j in range(0,ylen):
        dpdz = np.gradient(cartesian_model_p[i,j,:],model_z[1]-model_z[0])
        for k in range(0,zlen):        
           g = 2.74e2
           H = -cartesian_model_p[i,j,k] / dpdz[k]
           cs_square = cartesian_model_cs_square[i,j,k]
           omegac_square = np.sqrt(cs_square) / (2*H)
           N_square = g/H - g**2 / cs_square
           cartesian_model_omegac_square[i,j,k] = omegac_square
           cartesian_model_N_square[i,j,k] = N_square

  print '   Computing derivatives ... '
  cartesian_model_domegac_square_dx = np.zeros([xlen,ylen,zlen])
  cartesian_model_dN_square_dx = np.zeros([xlen,ylen,zlen])
  cartesian_model_dcs_square_dx = np.zeros([xlen,ylen,zlen])
  cartesian_model_dva_x_dx = np.zeros([xlen,ylen,zlen])
  cartesian_model_dva_y_dx = np.zeros([xlen,ylen,zlen])
  cartesian_model_dva_z_dx = np.zeros([xlen,ylen,zlen])
  cartesian_model_domegac_square_dy = np.zeros([xlen,ylen,zlen])
  cartesian_model_dN_square_dy = np.zeros([xlen,ylen,zlen])
  cartesian_model_dcs_square_dy = np.zeros([xlen,ylen,zlen])
  cartesian_model_dva_x_dy = np.zeros([xlen,ylen,zlen])
  cartesian_model_dva_y_dy = np.zeros([xlen,ylen,zlen])
  cartesian_model_dva_z_dy = np.zeros([xlen,ylen,zlen])
  cartesian_model_domegac_square_dz = np.zeros([xlen,ylen,zlen])
  cartesian_model_dN_square_dz = np.zeros([xlen,ylen,zlen])
  cartesian_model_dcs_square_dz = np.zeros([xlen,ylen,zlen])
  cartesian_model_dva_x_dz = np.zeros([xlen,ylen,zlen])
  cartesian_model_dva_y_dz = np.zeros([xlen,ylen,zlen])
  cartesian_model_dva_z_dz = np.zeros([xlen,ylen,zlen])
  # compute gradients in x
  xstep = cartesian_model_x[1]-cartesian_model_x[0]
  for j in range(0,ylen):
     for k in range(0,zlen):
        cartesian_model_domegac_square_dx[:,j,k] = np.gradient(cartesian_model_omegac_square[:,j,k],xstep)
        cartesian_model_dN_square_dx[:,j,k] = np.gradient(cartesian_model_N_square[:,j,k],xstep)
        cartesian_model_dcs_square_dx[:,j,k] = np.gradient(cartesian_model_cs_square[:,j,k],xstep)
        cartesian_model_dva_x_dx[:,j,k] = np.gradient(cartesian_model_va_x[:,j,k],xstep)
        cartesian_model_dva_y_dx[:,j,k] = np.gradient(cartesian_model_va_y[:,j,k],xstep)
        cartesian_model_dva_z_dx[:,j,k] = np.gradient(cartesian_model_va_z[:,j,k],xstep)
  # compute gradients in y
  ystep = cartesian_model_y[1]-cartesian_model_y[0]
  for i in range(0,xlen):
     for k in range(0,zlen):
        cartesian_model_domegac_square_dy[i,:,k] = np.gradient(cartesian_model_omegac_square[i,:,k],ystep)
        cartesian_model_dN_square_dy[j,:,k] = np.gradient(cartesian_model_N_square[i,:,k],ystep)
        cartesian_model_dcs_square_dy[j,:,k] = np.gradient(cartesian_model_cs_square[i,:,k],ystep)
        cartesian_model_dva_x_dy[j,:,k] = np.gradient(cartesian_model_va_x[i,:,k],ystep)
        cartesian_model_dva_y_dy[j,:,k] = np.gradient(cartesian_model_va_y[i,:,k],ystep)
        cartesian_model_dva_z_dy[j,:,k] = np.gradient(cartesian_model_va_z[i,:,k],ystep)
  # compute gradients in z
  zstep = cartesian_model_z[1]-cartesian_model_z[0]
  for i in range(0,xlen):
     for j in range(0,ylen):
        cartesian_model_domegac_square_dy[i,j,:] = np.gradient(cartesian_model_omegac_square[i,j,:],ystep)
        cartesian_model_dN_square_dy[i,j,:] = np.gradient(cartesian_model_N_square[i,j,:],ystep)
        cartesian_model_dcs_square_dy[i,j,:] = np.gradient(cartesian_model_cs_square[i,j,:],ystep)
        cartesian_model_dva_x_dz[i,j,:] = np.gradient(cartesian_model_va_x[i,j,:],ystep)
        cartesian_model_dva_y_dz[i,j,:] = np.gradient(cartesian_model_va_y[i,j,:],ystep)
        cartesian_model_dva_z_dz[i,j,:] = np.gradient(cartesian_model_va_z[i,j,:],ystep)

  # return results
  return cartesian_model_x, cartesian_model_y, cartesian_model_z, cartesian_model_cs_square, cartesian_model_va_x, cartesian_model_va_y, cartesian_model_va_z, cartesian_model_omegac_square, cartesian_model_N_square,cartesian_model_domegac_square_dx,cartesian_model_dN_square_dx,cartesian_model_dcs_square_dx,cartesian_model_dva_x_dx,cartesian_model_dva_y_dx,cartesian_model_dva_z_dx,cartesian_model_domegac_square_dy,cartesian_model_dN_square_dy,cartesian_model_dcs_square_dy,cartesian_model_dva_x_dy,cartesian_model_dva_y_dy,cartesian_model_dva_z_dy,cartesian_model_domegac_square_dy,cartesian_model_dN_square_dy,cartesian_model_dcs_square_dy,cartesian_model_dva_x_dz,cartesian_model_dva_y_dz,cartesian_model_dva_z_dz

#########################################

# load model
print 'Loading model in cylindrical coordinates ... '
model_r,model_z,model_vr,model_vz,model_Br,model_Bz,model_T,model_p,model_rho,model_tau,model_gamma1 = load_model('Rempel_sunspot_data')

print 'Make Cartesian model ... '
cartesian_model_x, cartesian_model_y, cartesian_model_z, cartesian_model_cs_square, cartesian_model_va_x, cartesian_model_va_y, cartesian_model_va_z, cartesian_model_omegac_square, cartesian_model_N_square,cartesian_model_domegac_square_dx,cartesian_model_dN_square_dx,cartesian_model_dcs_square_dx,cartesian_model_dva_x_dx,cartesian_model_dva_y_dx,cartesian_model_dva_z_dx,cartesian_model_domegac_square_dy,cartesian_model_dN_square_dy,cartesian_model_dcs_square_dy,cartesian_model_dva_x_dy,cartesian_model_dva_y_dy,cartesian_model_dva_z_dy,cartesian_model_domegac_square_dy,cartesian_model_dN_square_dy,cartesian_model_dcs_square_dy,cartesian_model_dva_x_dz,cartesian_model_dva_y_dz,cartesian_model_dva_z_dz = make_cartesian_model(model_r,model_z,model_vr,model_vz,model_Br,model_Bz,model_T,model_p,model_rho,model_tau,model_gamma1)


# set parameters
min_S = -1e7
n_S = 100


# set intial value
r0 = [    0.0,   0.0,  0.0 ] 
k0 = [   2.2361e-3,   0.0,  -2e-3] 
t0 = 0.0
Phi0 = [ r0[0], r0[1], r0[2], k0[0], k0[1], k0[2], t0 ]
omega = 3e-3

# wave frequency
additional_argument = (omega,)

# range of s
s = np.arange(0,min_S,min_S/n_S)

# call solver
Phi = int.odeint(func,Phi0,s,args=additional_argument)

print '>>>>>>>>>>>>>>>>>>>>>>>'
for i in range(0,n_S):
  print Phi[i,6],Phi[i,0],Phi[i,1],Phi[i,2]

plt.plot(Phi[:,0],Phi[:,2])
plt.xlabel('x')
plt.ylabel('z')
plt.show()


