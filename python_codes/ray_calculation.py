import scipy.integrate as int
import numpy as np
import matplotlib.pylab as plt


#########################################
def cbrt(x):
    from math import pow
    if x >= 0: 
	return pow(x, 1.0/3.0)
    else:
	return -pow(abs(x), 1.0/3.0)

#########################################
def polar(x, y, deg=0):		# radian if deg=0; degree if deg=1
    from math import hypot, atan2, pi
    if deg:
	return hypot(x, y), 180.0 * atan2(y, x) / pi
    else:
	return hypot(x, y), atan2(y, x)

#########################################
def quadratic(a, b, c=None):
    import math, cmath
    if c:		# (ax^2 + bx + c = 0)
	a, b = b / float(a), c / float(a)
    t = a / 2.0
    r = t**2 - b
    if r >= 0:		# real roots
	y1 = math.sqrt(r)
    else:		# complex roots
	y1 = cmath.sqrt(r)
    y2 = -y1
    return y1 - t, y2 - t

#########################################
def cubic(a, b, c, d=None):
    from math import cos
    if d:			# (ax^3 + bx^2 + cx + d = 0)
	a, b, c = b / float(a), c / float(a), d / float(a)
    t = a / 3.0
    p, q = b - 3 * t**2, c - b * t + 2 * t**3
    u, v = quadratic(q, -(p/3.0)**3)
    if type(u) == type(0j):	# complex cubic root
	r, w = polar(u.real, u.imag)
	y1 = 2 * cbrt(r) * cos(w / 3.0)
    else:			# real root
        y1 = cbrt(u) + cbrt(v)
    y2, y3 = quadratic(y1, p + y1**2)
    return y1 - t, y2 - t, y3 - t

#########################################
def load_model(dir):

  from astropy.io import fits as pyfits

  fitsfile = pyfits.open(dir+'/r.fits')
  model_r = fitsfile[0].data * 1e8
  fitsfile.close()

  fitsfile = pyfits.open(dir+'/z.fits')
  model_z = fitsfile[0].data * 1e8
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

  # set magnetic field to zero for now
  model_Br[:,:] = 0.0
  model_Bz[:,:] = 0.0


  g = 2.74e4
  dpdz = np.zeros([len(model_z),len(model_r)])
  for r_index in range(0,len(model_r)):
    dpdz[:,r_index] = (np.gradient(model_p[:,r_index],model_z[1]-model_z[0]))
  model_H = -model_p / dpdz
  model_cs_square = model_gamma1 * model_p / model_rho
  model_va_r = model_Br / np.sqrt(4*np.pi*model_rho)
  model_va_z = model_Bz / np.sqrt(4*np.pi*model_rho)
  model_va_square = model_va_r**2 + model_va_z**2
  model_omegac_square = np.sqrt(model_cs_square) / (2*model_H)
  model_N_square = g/model_H - g**2 / model_cs_square

  # set N_square to zero for now
  model_N_square[:,:] = 0.0
  # set omegac_square to zero for now
  model_omegac_square[:,:] = 0.0

  model_dcssquaredr = np.zeros([len(model_z),len(model_r)])
  for z_index in range(0,len(model_z)):
     model_dcssquaredr[z_index,:] = np.gradient(model_cs_square[z_index,:],model_r[1]-model_r[0])
  model_dcssquaredz = np.zeros([len(model_z),len(model_r)])
  for r_index in range(0,len(model_r)):
     model_dcssquaredz[:,r_index] = np.gradient(model_cs_square[:,r_index],model_z[1]-model_z[0])

  model_dNsquaredr = np.zeros([len(model_z),len(model_r)])
  for z_index in range(0,len(model_z)):
     model_dNsquaredr[z_index,:] = np.gradient(model_N_square[z_index,:],model_r[1]-model_r[0])
  model_dNsquaredz = np.zeros([len(model_z),len(model_r)])
  for r_index in range(0,len(model_r)):
     model_dNsquaredz[:,r_index] = np.gradient(model_N_square[:,r_index],model_z[1]-model_z[0])

  model_domegacsquaredr = np.zeros([len(model_z),len(model_r)])
  for z_index in range(0,len(model_z)):
     model_domegacsquaredr[z_index,:] = np.gradient(model_omegac_square[z_index,:],model_r[1]-model_r[0])
  model_domegacsquaredz = np.zeros([len(model_z),len(model_r)])
  for r_index in range(0,len(model_r)):
     model_domegacsquaredz[:,r_index] = np.gradient(model_omegac_square[:,r_index],model_z[1]-model_z[0])

  model_dvasquaredr = np.zeros([len(model_z),len(model_r)])
  for z_index in range(0,len(model_z)):
     model_dvasquaredr[z_index,:] = np.gradient(model_va_square[z_index,:],model_r[1]-model_r[0])
  model_dvasquaredz = np.zeros([len(model_z),len(model_r)])
  for r_index in range(0,len(model_r)):
     model_dvasquaredz[:,r_index] = np.gradient(model_va_square[:,r_index],model_z[1]-model_z[0])



  return model_r,model_z,model_vr,model_vz,model_Br,model_Bz,model_cs_square,model_va_r,model_va_z,model_omegac_square,model_N_square,model_dcssquaredr,model_dcssquaredz,model_dvasquaredr,model_dvasquaredz,model_dNsquaredr,model_dNsquaredz,model_domegacsquaredr,model_domegacsquaredz

#########################################
def func(Phi,s,omega):

  # input values
  r = np.array([Phi[0],Phi[1],Phi[2]])
  k = np.array([Phi[3],Phi[4],Phi[5]])
  t = Phi[6]
  k_square = k[0]**2 + k[1]**2 + k[2]**2

  # get model values
  B,cs_square,va,va_square,N_square,omegac_square,dcssquaredxi,dvajdxi,dvasquaredxi,dNsquaredxi,domegacsquaredxi = get_model(r)


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
            - omega**2 * k_square * (dcssquaredxi + dvasquaredxi)
            + k_square * ( np.dot(k,va) * dcssquaredxi + cs_square * np.dot(dvajdxi,k) )
            - (omega**2*domegacsquaredxi -k[2]**2*(omegac_square*dcssquaredxi+cs_square*domegacsquaredxi) )
            + (k[0]**2+k[1]**2)*( N_square*dcssquaredxi + cs_square*dNsquaredxi ) 
          )

  # -func[2] = d F / d omega
  dFdomega = ( 
               + 4*omega**3
               - 2*omega*(cs_square + va_square)*k_square 
               - omegac_square * 2 * omega
             )

  # return result
#  print s,'||',dFdki[0],dFdki[1],dFdki[2],'|',-dFdxi[0],-dFdxi[1],-dFdxi[2],'|',-dFdomega
  return [dFdki[0],dFdki[1],dFdki[2],-dFdxi[0],-dFdxi[1],-dFdxi[2],-dFdomega]

#########################################
def get_model(r):

  x_value = r[0]
  y_value = r[1]
  z_value = r[2]

  # compute radius and angle
  r_value = np.sqrt(x_value**2 + y_value**2)
  if r_value > 0:
     cosalpha = x_value / r_value
     sinalpha = y_value / r_value
  else:
     cosalpha = 1.0
     sinalpha = 0.0

  # find nearest r and z in dataset
  best_dr = 1e99
  nearest_r_index = -1
  for r_index in range(0,len(model_r)):
     dr = abs(r_value - model_r[r_index])
     if dr < best_dr:
        best_dr = dr
        nearest_r_index = r_index 
  best_dz = 1e99
  nearest_z_index = -1
  for z_index in range(0,len(model_z)):
     dz = abs(z_value - model_z[z_index])
     if dz < best_dz:
        best_dz = dz
        nearest_z_index = z_index 

  Bz = model_Bz[nearest_z_index,nearest_r_index]
  Br = model_Br[nearest_z_index,nearest_r_index]
  va_r = model_va_r[nearest_z_index,nearest_r_index]
  va_z = model_va_z[nearest_z_index,nearest_r_index]
  cs_square = model_cs_square[nearest_z_index,nearest_r_index]
  omegac_square = model_omegac_square[nearest_z_index,nearest_r_index]
  N_square = model_N_square[nearest_z_index,nearest_r_index]

  Bx = Br * cosalpha
  By = Br * sinalpha
  B = np.array([Bx,By,Bz])

  va_x = abs(va_r * cosalpha)
  va_y = abs(va_r * sinalpha)
  va = np.array([va_x,va_y,va_z])
  va_square = np.dot(va,va)

  dcssquaredr = model_dcssquaredr[nearest_z_index,nearest_r_index]
  dcssquaredz = model_dcssquaredz[nearest_z_index,nearest_r_index]
  dcssquaredx = dcssquaredr * cosalpha
  dcssquaredy = dcssquaredr * sinalpha
  dcssquaredxi = np.array([dcssquaredx,dcssquaredy,dcssquaredz])

  dNsquaredr = model_dNsquaredr[nearest_z_index,nearest_r_index]
  dNsquaredz = model_dNsquaredz[nearest_z_index,nearest_r_index]
  dNsquaredx = dNsquaredr * cosalpha
  dNsquaredy = dNsquaredr * sinalpha
  dNsquaredxi = np.array([dNsquaredx,dNsquaredy,dNsquaredz])

  domegacsquaredr = model_domegacsquaredr[nearest_z_index,nearest_r_index]
  domegacsquaredz = model_domegacsquaredz[nearest_z_index,nearest_r_index]
  domegacsquaredx = domegacsquaredr * cosalpha
  domegacsquaredy = domegacsquaredr * sinalpha
  domegacsquaredxi = np.array([domegacsquaredx,domegacsquaredy,domegacsquaredz])

  dvasquaredr = model_dvasquaredr[nearest_z_index,nearest_r_index]
  dvasquaredz = model_dvasquaredz[nearest_z_index,nearest_r_index]
  dvasquaredx = dvasquaredr * cosalpha
  dvasquaredy = dvasquaredr * sinalpha
  dvasquaredxi = np.array([dvasquaredx,dvasquaredy,dvasquaredz])

  dvajdxi = np.zeros([3,3])

  return B,cs_square,va,va_square,N_square,omegac_square,dcssquaredxi,dvajdxi,dvasquaredxi,dNsquaredxi,domegacsquaredxi 

#########################################

# load model
print 'Loading model in cylindrical coordinates ... '
model_r,model_z,model_vr,model_vz,model_Br,model_Bz,model_cs_square,model_va_r,model_va_z,model_omegac_square,model_N_square,model_dcssquaredr,model_dcssquaredz,model_dvasquaredr,model_dvasquaredz,model_dNsquaredr,model_dNsquaredz,model_domegacsquaredr,model_domegacsquaredz = load_model('Rempel_sunspot_data')


# domain size
print 'Domain size is from ',-np.max(model_r/1e8),' to ',np.max(model_r/1e8),' Mm by ',np.min(model_z/1e8),' to ',np.max(model_z/1e8),' Mm'
print 'Sound speed is in the range of ',np.sqrt(np.min(model_cs_square))/1e5,' to ',np.sqrt(np.max(model_cs_square))/1e5,' km/sec'


# set intial values
#nu = 3e-3
omega = 3e-3
r0 = np.array([   5.0,   0.0,   -3.0 ])*1e8
k0_direction = np.array([  2.0,   0.0,  -1.0])

k0_norm = k0_direction / np.dot(k0_direction,k0_direction)
B,cs_square,va,va_square,N_square,omegac_square,dcssquaredxi,dvajdxi,dvasquaredxi,dNsquaredxi,domegacsquaredxi = get_model(r0)
power0 = omega**4 - omegac_square*omega**2
power2 = -omega**2*(cs_square+va_square) + omegac_square*cs_square*k0_norm[2]**2 + cs_square*N_square*(k0_norm[0]**2+k0_norm[1]**2)
power4 = cs_square*np.dot(va,k0_norm)**2
if power4 == 0:
   k0_mag = np.sqrt( - power0 / power2 )
else:
   k0_mag = np.sqrt(max(quadratic(power4,power2,power0)))
k0 = k0_norm * k0_mag
print 'Initial wave vector is :',k0

t0 = 0.0
Phi0 = [ r0[0], r0[1], r0[2], k0[0], k0[1], k0[2], t0 ]

# wave frequency
additional_argument = (omega,)

# range of s
min_S = -1e5
n_S = 100

desired_maxtime = 30.0 # in mins

print 'Running solver ... '
maxt = 0.0
while maxt < desired_maxtime:
  s = np.arange(0,min_S,min_S/n_S)
  Phi = int.odeint(func,Phi0,s,args=additional_argument)
  maxt = np.max(Phi[:,6])/60.
  min_S = min_S * desired_maxtime / maxt * 1.05
print 'Total travel time:'+str(maxt)+' mins'


#print '>>>>>>>>>>>>>>>>>>>>>>>'
#for i in range(0,n_S):
#  print Phi[i,6],Phi[i,0],Phi[i,1],Phi[i,2]

plt.plot(Phi[:,0]/1e8,Phi[:,2]/1e8)
plt.xlabel('x [Mm]')
plt.ylabel('z [Mm]')
plt.xlim([-np.max(model_r)/1e8,np.max(model_r)/1e8])
plt.ylim([np.min(model_z)/1e8,np.max(model_z)/1e8])
plt.show()


