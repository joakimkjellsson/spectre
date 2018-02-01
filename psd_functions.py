import numpy as np
from netCDF4 import Dataset 
import os,sys,time
from scipy.fftpack import fft, ifft, fftn, ifftn
from scipy.signal import periodogram, hamming, tukey
import scipy.stats as stats

import matplotlib.pyplot as plt
import multiprocessing as mp

def find_time(starttime,endtime,nd):
   startday  = int(starttime[6:8])
   startmon  = int(starttime[4:6])
   startyear = int(starttime[0:4])
   
   endday    = int(endtime[6:8])
   endmon    = int(endtime[4:6])
   endyear   = int(endtime[0:4])
   
   ii = 1
   iday  = startday
   imon  = startmon
   iyear = startyear
   itime = iyear * 10000 + imon * 100 + iday
   
   idaymax = [31,28,31,30,31,30,31,31,30,31,30,31]
   
   idays = np.array([iday])
   imons = np.array([imon])
   iyears = np.array([iyear])
   
   if (int(endtime) - int(starttime) >0 ):
      while( itime < endyear * 10000 + endmon * 100 + endday ):
         
         iday = iday + nd
         if (iday > idaymax[imon-1]):
            iday = iday - idaymax[imon-1]
            imon = imon + 1
         
         if (imon > 12):
            imon  = 1
            iyear = iyear + 1
         
         itime = iyear * 10000 + imon * 100 + iday
         
         idays = np.append(idays,iday)
         imons = np.append(imons,imon)
         iyears = np.append(iyears,iyear)
      
   
   return iyears,imons,idays


def find_kindex(depth,deptht,depthw):
   """
   Find the k index of a given depth 
   """
   
   kdiff = (deptht - depth)
   kt = np.where( np.abs(kdiff) == np.min(np.abs(kdiff)) )[0]
   kdiff = (depthw - depth)
   kw = np.where( np.abs(kdiff) == np.min(np.abs(kdiff)) )[0]
   
   return kt.max(),kw.max()


def calculate_rhines_scale(urms,mean_lat):
   """
   """
   
   b = 2. * 7.2921150 * 1e-5 * np.cos( np.pi/180. * mean_lat ) / (6371.0*1000)
   r = 2 * np.pi * np.sqrt( urms / b )
   return r
   

def set_region(region,tlon,tlat):
   
   info = {}
   info['region'] = region
   
   ## Jet in double-gyre run
   if (region == 'gyre-jet'):
      lon0 = 0
      lon1 = 16
      lat0 = 19
      lat1 = 31
   elif (region == 'gyre-ext'):
      lon0 = 15
      lon1 = 40
      lat0 = 19
      lat1 = 31
   elif (region == 'gyre-global'):
      lon0 = 0
      lon1 = 43
      lat0 = 5
      lat1 = 48
   
   ## Jet in Rob's MITgcm run
   elif (region == 'MITgcm-jet'):
      lon0 = 1.5e5
      lon1 = 1.125e6
      lat0 = 1.5e6
      lat1 = 2.25e6
   elif (region == 'MITgcm-down'):
      lon0 = 200*7500
      lon1 = 350*7500
      lat0 = 1.5e6
      lat1 = 2.25e6
   
   ## Regions on in the real world
   
   ## ACC   
   elif (region == 'acc-full'):
      lon0 = 30.
      lon1 = 70.
      lat0 = -45
      lat1 = -35
   
   ## same region as Scott & Wang JPO,2005
   elif (region == 'acc-scott'):
      lon0 = -128
      lon1 = -112
      lat0 =  -65
      lat1 =  -49
      
   elif (region == 'south-pacific1'):
      lon0 = -160
      lon1 = -100
      lat0 =  -65
      lat1 =  -47
   
   elif (region == 'south-pacific2'):
      lon0 = -175
      lon1 =  -80
      lat0 =  -68
      lat1 =  -47
   
   ## Agulhas retroflection
   elif (region == 'agulhas-retro'):
      lon0 = 25
      lon1 = 68
      lat0 = -45.5
      lat1 = -37
   elif (region == 'agulhas-rings-big'):
      lon0 = -10
      lon1 = 60
      lat0 = -50
      lat1 = 10
      rot = -20
   elif (region == 'agulhas-current-big'):
      lon0 = -5
      lon1 = 67
      lat0 = -50
      lat1 = -10
   elif (region == 'agulhas-retro-2'):
      lon0 = 14
      lon1 = 68
      lat0 = -54
      lat1 = -45.5
            
   ## Gulf Stream
   elif (region[0:10] == 'gulfstream'):
      lat0 = 36.3
      lat1 = 42.
      if (region == 'gulfstream-full'):
         lon0 = -63.5
         lon1 = -34.
      elif (region == 'gulfstream-up'):
         lon0 = -68.5
         lon1 = -40.
      elif (region == 'gulfstream-down'):
         lon0 = -68.5
         lon1 = -50.
      else:
         print ' region ',region,' is not a recognised part of the Gulf Stream!'
         sys.exit()
   
   ## Kuroshio
   elif (region[0:8] == 'kuroshio'):
      lat0 = 30.0
      lat1 = 40.0
      if (region == 'kuroshio-full'):
         lon0 = 141.5
         lon1 = -176.
      ## downstream
      elif (region == 'kuroshio-down'):
         lon0 = 163.
         lon1 = -176 #-165.
      ## upstream
      elif (region == 'kuroshio-up'):
         lon0 = 141.5
         lon1 = 163. 
      elif (region == 'kuroshio-full2'):
         lon0 = 149.
         lon1 = -176.
         lat0 = 28.0
         lat1 = 43.0
      elif (region == 'kuroshio-full3'):
         lon0 = 145.
         lon1 = 180.
         lat0 = 28.0
         lat1 = 43.0
      
   elif (region == 'north-pacific-1'):
      lon0 = 141.5
      lon1 = -176.
      lat0 = 20.0
      lat1 = 32.5
   
   elif (region == 'argentine'):
      lon0 = -52.75
      lon1 = -30.
      lat0 = -53.
      lat1 = -40.
   
   elif (region == 'east-pacific-1'):
      lon0 = -121.5
      lon1 = -82.5
      lat0 = -46.4
      lat1 = -25.
   
   else:
      print ' Region: ',region,' is not a region I recognise! '
      print ' If you have not typed anything stupid (which the esteemed Dr Kjellsson often does) '
      print ' you should add the lon,lats for the region to the function set_region '
      sys.exit()
   
   info['lon0'] = lon0
   info['lon1'] = lon1
   info['lat0'] = lat0
   info['lat1'] = lat1
   
   ##
   ## We have the lon0, lon1, lat0, lat1 of the region
   ## 
   ## Find i,j indices for the region
   ##
   
   #if (grid_type == 'll'):
   if 1:
      ## e.g. if lon1 = 160E and lon0 = 160W
      ## then we convert lon to [0,360] to make lon monotonic
      if (lon1 < lon0):
         tmp_lon1 = lon1+360
         tlon = np.where(tlon<0,tlon+360,tlon)
      else:
         tmp_lon1 = lon1
      
      ## find all indices in [lon0,lon1],[lat0,lat1]
      ind = np.where( (tlon>=lon0) & (tlon<=tmp_lon1) & (tlat>=lat0) & (tlat<=lat1) )
      
      ## if we converted to [0,360] before, we convert back now
      if (lon1 < lon0):
         tlon = np.where(tlon>180,tlon-360,tlon)
         
   ## we should be able to select a subdomain if the grid is x-y rather than lon-lat
   ## I think this is common for idealised MITgcm runs?
   ## In that case, maybe the user can set the x0, x1, y0, y1 variables here?
   #elif (grid_type == 'xy'):
   #   ind = np.where( (tlon >= lon0) & (tlon <=lon1) & (tlat >= lat0) & (tlat <=lat1) )
   
   ## ind[1] has i indices, and ind[0] has j indices   
   i_ind = ind[1]
   j_ind = ind[0]
   ## pick i0,i1,j0,j1 as min/max of i and j arrays
   ## On a wonky grid like e.g. NEMO tripolar grid, 
   ## this may mean that we include places where lon < lon0
   i0 = i_ind.min()
   i1 = i_ind.max()
   j0 = j_ind.min()
   j1 = j_ind.max()
   
   info['i0'] = i0
   info['i1'] = i1
   info['j0'] = j0
   info['j1'] = j1
   
   return info
   

def set_mitgcm_rob_grid(ufile):
   """
   """
   
   print ufile
   ncu = Dataset(ufile,'r')
   
   xx   = ncu.variables['X'][:]
   yy   = ncu.variables['Y'][:]      
   lon  = xx
   lat  = yy
   ulon, ulat = np.meshgrid(lon,lat)
   vlon, vlat = np.meshgrid(lon,lat)
   tlon, tlat = np.meshgrid(lon,lat)
   ncu.close()
   
   return ulon,ulat
   

def read_mitgcm_rob(ufile,tfile,i0,i1,j0,j1,step0=0,step1=1,lw=False):
   """
   """
   
   print ufile
   ncu = Dataset(ufile,'r')
   print tfile
   nct = Dataset(tfile,'r')
   
   xx   = ncu.variables['X'][i0:i1]
   yy   = ncu.variables['Y'][j0:j1]      
   ulon, ulat = np.meshgrid(xx,yy)
   vlon, vlat = np.meshgrid(xx,yy)
   tlon, tlat = np.meshgrid(xx,yy)
   utmp = ncu.variables['UVEL'][step0:step1,:,j0:j1,i0:i1+1]
   uvel_full = 0.5 * (utmp[:,:,:,1:] + utmp[:,:,:,0:-1])
   vtmp = ncu.variables['VVEL'][step0:step1,:,j0:j1+1,i0:i1]
   vvel_full = 0.5 * (vtmp[:,:,1:,:] + vtmp[:,:,0:-1,:])
   etan_full = ncu.variables['ETAN'][step0:step1,:,j0:j1,i0:i1]
   dzt = np.ones(uvel_full[0,:,:,:].shape) * 3000.
   
   for s in range(step0,step1):
      if (s == step0):
         taux_full = nct.variables['taux'][:,j0:j1,i0:i1]
         tauy_full = taux_full * 0.
      else:
         taux_full = np.concatenate( (taux_full, nct.variables['taux'][:,j0:j1,i0:i1]), axis=0 )
         tauy_full = np.concatenate( (tauy_full, taux_full*0), axis=0 )
   
   ## sometimes uvel, vvel are undefined
   ## mask all values larger than 20m/s
   umax = max( np.abs(uvel_full.min()), uvel_full.max() )
   vmax = max( np.abs(vvel_full.min()), vvel_full.max() )
   if (umax > 10 or vmax > 10):
      print ' WARNING!!! uvel > 10 or vmax > 10'
      
   umask = np.where( np.abs(uvel_full) > 20, 1, 0)
   vmask = np.where( np.abs(vvel_full) > 20, 1, 0)
   uvel_full = np.ma.array( uvel_full, mask=umask )
   vvel_full = np.ma.array( vvel_full, mask=vmask )
   etan_full = np.ma.array( etan_full, mask=vmask )
   print taux_full.shape,umask.shape
   taux_full = np.ma.array( taux_full, mask=umask[0,0,:,:] )
   tauy_full = np.ma.array( tauy_full, mask=vmask[0,0,:,:] )
   umax = max( np.abs(uvel_full.min()), uvel_full.max() )
   vmax = max( np.abs(vvel_full.min()), vvel_full.max() )
   if (umax > 10 or vmax > 10):
      print ' uvel > 10 or vmax > 10 or wmax > 10'
      sys.exit()
      
   ncu.close()
   nct.close()
   
   data = {}
   data['ulon'] = ulon
   data['ulat'] = ulat
   data['vlon'] = vlon
   data['vlat'] = vlat
   data['uvel'] = uvel_full
   data['vvel'] = vvel_full
   data['etan'] = etan_full
   data['tlon'] = tlon
   data['tlat'] = tlat
   data['taux'] = taux_full
   data['tauy'] = tauy_full
   data['dzt']  = dzt
   
   return data


def calculate_uv_gradients_xy(uvel,vvel,dxu,dyu,dxv,dyv,dxt,dyt,dxf,dyf):
   """
   Calculates du/dx, du/dy, dv/dx, dv/dy 
   and nabla^2 u, nabla^2 v, 
   and nabla^4 u, nabla^4 v
   at T points if u,v are given at U,V points 
   """
   
   dudx_xy = np.zeros(uvel.shape)
   dvdx_xy = np.zeros(uvel.shape)
   dudy_xy = np.zeros(uvel.shape)
   dvdy_xy = np.zeros(uvel.shape)
   
   lapu = np.zeros(uvel.shape)
   lapv = np.zeros(uvel.shape)
   blpu = np.zeros(uvel.shape)
   blpv = np.zeros(uvel.shape)
   
   for jk in range(0,uvel.shape[0]):
      dudx_xy[jk,1:-1,1:-1] = (uvel[jk,1:-1,1:-1] - uvel[jk,1:-1,0:-2]) / dxu[1:-1,1:-1]
      dudy_xy[jk,1:-1,1:-1] = 0.5 * ( (uvel[jk,2:,1:-1] - uvel[jk,0:-2,1:-1]) / (2*dyu[1:-1,1:-1]) + \
                                      (uvel[jk,2:,0:-2] - uvel[jk,0:-2,0:-2]) / (2*dyu[1:-1,1:-1]) )
      dvdx_xy[jk,1:-1,1:-1] = 0.5 * ( (vvel[jk,1:-1,2:] - vvel[jk,1:-1,0:-2]) / (2*dxv[1:-1,1:-1]) + \
                                      (vvel[jk,0:-2,2:] - vvel[jk,0:-2,0:-2]) / (2*dxv[1:-1,1:-1]) )
      dvdy_xy[jk,1:-1,1:-1] = (vvel[jk,1:-1,1:-1] - vvel[jk,0:-2,1:-1]) / dyv[1:-1,1:-1]
      
      rot  = np.zeros((uvel.shape[1],uvel.shape[2]))
      div  = np.zeros((uvel.shape[1],uvel.shape[2]))
      zuf  = np.zeros((uvel.shape[1],uvel.shape[2]))
      zut  = np.zeros((uvel.shape[1],uvel.shape[2]))
         
      # rot = dv/dx - du/dy at F point
      rot[1:-1,1:-1]  = (vvel[jk,1:-1,2:  ] - vvel[jk,1:-1,1:-1])/dxf[1:-1,1:-1] - (uvel[jk,2:  ,1:-1] - uvel[jk,1:-1,1:-1])/dyf[1:-1,1:-1]
      # div = du/dx + dv/dy at T point
      div[1:-1,1:-1]  = (uvel[jk,1:-1,1:-1] - uvel[jk,1:-1,0:-2])/dxt[1:-1,1:-1] + (vvel[jk,1:-1,1:-1] - vvel[jk,0:-2,1:-1])/dyt[1:-1,1:-1]
      # lap(u) = ddiv/dx - drot/dy at U point
      lapu[jk,1:-1,1:-1] = (div[1:-1,2:  ] - div[1:-1,1:-1])/dxu[1:-1,1:-1] - (rot[1:-1,1:-1] - rot[0:-2,1:-1])/dyu[1:-1,1:-1]
      # lap(v) = drot/dx + ddiv/dy at V point
      lapv[jk,1:-1,1:-1] = (rot[1:-1,1:-1] - rot[1:-1,0:-2])/dxv[1:-1,1:-1] + (div[2:  ,1:-1] - div[1:-1,1:-1])/dyv[1:-1,1:-1]
      # zuf = dlapv/dx - dlapu/dy at F point
      zuf[1:-1,1:-1]  = (lapv[jk,1:-1,2:  ] - lapv[jk,1:-1,1:-1])/dxf[1:-1,1:-1] - (lapu[jk,2:  ,1:-1] - lapu[jk,1:-1,1:-1])/dyf[1:-1,1:-1]
      # zut = dlapu/dx + dlapv/dy at T point
      zut[1:-1,1:-1]  = (lapu[jk,1:-1,1:-1] - lapu[jk,1:-1,0:-2])/dxt[1:-1,1:-1] + (lapv[jk,1:-1,1:-1] - lapv[jk,0:-2,1:-1])/dyt[1:-1,1:-1]
      # blpu = dzut/dx - dzuf/dy at U point
      blpu[jk,1:-1,1:-1] = (zut[1:-1,2:  ] - zut[1:-1,1:-1])/dxu[1:-1,1:-1] - (zuf[1:-1,1:-1] - zuf[0:-2,1:-1])/dyu[1:-1,1:-1]
      # blpv = dzuf/dx + dzut/dy at V point
      blpv[jk,1:-1,1:-1] = (zuf[1:-1,1:-1] - zuf[1:-1,0:-2])/dxv[1:-1,1:-1] + (zut[2:  ,1:-1] - zut[1:-1,1:-1])/dyv[1:-1,1:-1]
   
   data = {}
   data['dudx_xy'] = dudx_xy
   data['dvdx_xy'] = dvdx_xy
   data['dudy_xy'] = dudy_xy
   data['dvdy_xy'] = dvdy_xy
   data['lapu_xy'] = lapu
   data['lapv_xy'] = lapv
   data['blpu_xy'] = blpu
   data['blpv_xy'] = blpv
   
   return data


def calculate_bottom_friction_xy(uvel,vvel,dzt,kbot,Cd=-1e-3,bg_tke=2.5e-3,mode='quadratic'):
   """
   """
   
   ## calculate bottom friction
   ## kbot is the bottom level
   if (0):
      vu = np.zeros(uvel[:,-1,:,:].shape)
      uv = np.zeros(vvel[:,-1,:,:].shape)
      vu[:,1:-1,1:-1] = 0.25 * (vvel[:,kbot,1:-1,1:-1] + vvel[:,kbot,0:-2,1:-1] + vvel[:,kbot,1:-1,2:] + vvel[:,kbot,0:-2,2:])
      uv[:,1:-1,1:-1] = 0.25 * (uvel[:,kbot,1:-1,1:-1] + uvel[:,kbot,1:-1,0:-2] + uvel[:,kbot,2:,1:-1] + uvel[:,kbot,2:,0:-2])
      utrd_bfr2 = -1e-3 * np.sqrt( uvel[:,kbot,:,:]**2 + vu**2 + 2.5 * 1e-3 ) * uvel[:,kbot,:,:] / dzt.sum(axis=1)
      vtrd_bfr2 = -1e-3 * np.sqrt( uv**2 + vvel[:,kbot,:,:]**2 + 2.5 * 1e-3 ) * vvel[:,kbot,:,:] / dzt.sum(axis=1)
   
   Hdep = dzt.sum(axis=0)
   if (mode == 'quadratic'):
      bfru = Cd * np.sqrt( uvel[kbot,:,:]**2 + vvel[kbot,:,:]**2 + bg_tke ) * uvel[kbot,:,:] / Hdep
      bfrv = Cd * np.sqrt( uvel[kbot,:,:]**2 + vvel[kbot,:,:]**2 + bg_tke ) * vvel[kbot,:,:] / Hdep
   
   data = {}
   data['bfru'] = bfru
   data['bfrv'] = bfrv
   
   return data


def calculate_pressure_gradient_xy(ssh,rhd,dxt,dyt,dzt,grav=9.81):
   """
   """
   
   spgu = np.zeros(rhd.shape)
   spgv = np.zeros(rhd.shape)
   hpgu = np.zeros(rhd.shape)
   hpgv = np.zeros(rhd.shape)
   
   for jk in range(0,rhd.shape[0]):
      spgu[jk,1:-1,1:-1] = -grav * (ssh[1:-1,2:] - ssh[1:-1,0:-2])/(2*dxt[1:-1,1:-1])
      spgv[jk,1:-1,1:-1] = -grav * (ssh[2:,1:-1] - ssh[0:-2,1:-1])/(2*dyt[1:-1,1:-1])
      if (jk == 0):
         coeff = -grav * 0.5 * dzt[jk,1:-1,1:-1]
         hpgu[jk,1:-1,1:-1] = coeff * (rhd[jk,1:-1,2:] - rhd[jk,1:-1,0:-2])/(2*dxt[1:-1,1:-1]) 
         hpgv[jk,1:-1,1:-1] = coeff * (rhd[jk,2:,1:-1] - rhd[jk,0:-2,1:-1])/(2*dyt[1:-1,1:-1]) 
      else:
         coeff = -grav * dzt[jk,1:-1,1:-1]
         hpgu[jk,1:-1,1:-1] = hpgu[jk-1,1:-1,1:-1] + coeff * (rhd[jk,1:-1,2:] - rhd[jk,1:-1,0:-2])/(2*dxt[1:-1,1:-1])
         hpgv[jk,1:-1,1:-1] = hpgv[jk-1,1:-1,1:-1] + coeff * (rhd[jk,2:,1:-1] - rhd[jk,0:-2,1:-1])/(2*dyt[1:-1,1:-1])
   
   ## force to zero on boundaries
   spgu[:,0,:] = 0. ; spgu[:,-1,:] = 0. ; spgu[:,:,0] = 0. ; spgu[:,:,-1] = 0.
   spgv[:,0,:] = 0. ; spgv[:,-1,:] = 0. ; spgv[:,:,0] = 0. ; spgv[:,:,-1] = 0.
   hpgu[:,0,:] = 0. ; hpgu[:,-1,:] = 0. ; hpgu[:,:,0] = 0. ; hpgu[:,:,-1] = 0.
   hpgv[:,0,:] = 0. ; hpgv[:,-1,:] = 0. ; hpgv[:,:,0] = 0. ; hpgv[:,:,-1] = 0.
      
   preu = spgu + hpgu
   prev = spgv + hpgv
   
   data = {}
   data['spgu'] = spgu
   data['spgv'] = spgv
   data['hpgu'] = hpgu
   data['hpgv'] = hpgv
   data['preu'] = preu
   data['prev'] = prev
      
   return data
   

def rotate_box(xlist,ylist,alpha):
   """
   """
   
   a = alpha*np.pi/180.
   Lx = xlist[1]-xlist[0]
   Ly = ylist[2]-ylist[1]
   
   x1p = xlist[0]
   y1p = ylist[0]
   
   x2p = xlist[0] + np.cos(a) * Lx
   y2p = ylist[0] + np.sin(a) * Lx
   
   x3p = x2p - np.sin(a) * Ly 
   y3p = y2p + np.cos(a) * Ly
   
   x4p = xlist[0] - np.sin(a) * Lx
   y4p = ylist[0] + np.cos(a) * Lx
   
   return [x1p,x2p,x3p,x4p],[y1p,y2p,y3p,y4p] 
   
   
def calculate_vorticity(uvel,vvel,xx,yy):
   """
   """
   
   rot = np.zeros(uvel.shape)
   
   dvdx = (vvel[1:-1,2:]-vvel[1:-1,0:-2]) / (xx[1:-1,2:]-xx[1:-1,0:-2])
   dudy = (uvel[2:,1:-1]-uvel[0:-2,1:-1]) / (yy[2:,1:-1]-yy[0:-2,1:-1])
   
   rot[1:-1,1:-1] = dvdx[:,:] - dudy[:,:]
   
   return rot
   

def calculate_ke(k2D,l2D,uvel,vvel,laverage=True):
   """
   """
   uhat = fftn(uvel)
   vhat = fftn(vvel)
   z = 0.5 * np.real(uhat * np.conj(uhat) + vhat * np.conj(vhat))
   if (laverage):
      nn = (uvel.shape[1]**2 * uvel.shape[0]**2)
      z = z/float(nn)   
   return z


def calculate_cke(k2D,l2D,uvel,vvel,dz,laverage=True):
   """
   Calculate kinetic energy for each layer
   """
   
   dtot = np.sum(dz,axis=0)   
   for jk in range(0,uvel.shape[0]):
      ## weight is sqrt(dz/H) 
      ## since we multiply u*u, this ends up being
      ## u*u*dz/H
      ## which is what we want to sum
      weight = dz[jk] / dtot
      uhat = fftn(uvel[jk,:,:])
      vhat = fftn(vvel[jk,:,:])
      z = 0.5 * np.real(uhat * np.conj(uhat) + vhat * np.conj(vhat)) * weight
      if (jk == 0):
         cke = z
      else:
         cke += z
   
   if (laverage):
      nn = (uvel.shape[1]**2 * uvel.shape[2]**2)
      cke = cke / float(nn)
         
   return cke


def calculate_ens(k2D,l2D,rot):
   """
   """
   rhat = fftn(rot)
   z = 0.5 * np.real(rhat * np.conj(rhat))
      
   return z

def calculate_spectral_flux(kx,ky,uvel,vvel):
   """
   Calculate spectral flux 
   We assume du/dt = -u du/dx - v du/dy
             dv/dt = -u dv/dx - v dv/dy
   """
   
   uhat = fftn(uvel)
   vhat = fftn(vvel)
   i = np.complex(0,1)
   # du/dx in x,y
   ddx_u = np.real( ifftn(i*kx*uhat) )
   # du/dy in x,y
   ddy_u = np.real( ifftn(i*ky*uhat) )
   # dv/dx in x,y
   ddx_v = np.real( ifftn(i*kx*vhat) )
   # dv/dy in x,y
   ddy_v = np.real( ifftn(i*ky*vhat) )
   
   # adv_u = u * du/dx + v * du/dy
   adv_u = -uvel * ddx_u - vvel * ddy_u
   # adv_v = u * dv/dx + v * dv/dy
   adv_v = -uvel * ddx_v - vvel * ddy_v
   
   # KE trend from advection: 
   # - u * adv_u - v * adv_v
   # in spectral space
   # The minus sign arises as advection 
   # is on the RHS of the momentum eqs. 
   Tkxky = np.real( np.conj(uhat)*fftn(adv_u) + \
                    np.conj(vhat)*fftn(adv_v) )   #[m2/s3]
   
   return Tkxky


def calculate_spectral_flux_baroclinic_barotropic(kx,ky,u_bt,v_bt,u_bc,v_bc):
   """
   Calculate spectral flux for a triad
   i.e. u_bt * u_bc * du_bc/dx + u_bt * v_bc * du_bc/dy + 
        v_bt * u_bc * dv_bc/dx + v_bt * v_bc * dv_bc/dy
   """
   
   print ' WARNING: Gradients for barotropic/baroclinic velocities '
   print '          are calculated in spectral space. '
   print '          You should do the calculations in physical space to agree '
   print '          with a grid-point model '
   
   
   i = np.complex(0,1)
   
   uhat_bt = fftn(u_bt)
   vhat_bt = fftn(v_bt)
   
   nx = u_bc.shape[2]
   ny = u_bc.shape[1]
   nz = u_bc.shape[0]
   
   ddx_u_bt = np.real( ifftn(i*kx*uhat_bt) ) # du_bt/dx
   ddy_u_bt = np.real( ifftn(i*ky*uhat_bt) ) # du_bt/dy
   ddx_v_bt = np.real( ifftn(i*kx*vhat_bt) ) # dv_bt/dx
   ddy_v_bt = np.real( ifftn(i*ky*vhat_bt) ) # dv_bt/dy
   
   for jk in range(0,nz):
      uhat_bc = fftn(u_bc[jk,:,:])
      vhat_bc = fftn(v_bc[jk,:,:])
      
      ddx_u_bc = np.real( ifftn(i*kx*uhat_bc) ) # du_bc/dx
      ddy_u_bc = np.real( ifftn(i*ky*uhat_bc) ) # du_bc/dy
      ddx_v_bc = np.real( ifftn(i*kx*vhat_bc) ) # dv_bc/dx
      ddy_v_bc = np.real( ifftn(i*ky*vhat_bc) ) # dv_bc/dy
      
      if (jk == 0):
         # adv_u = u * du/dx + v * du/dy
         adv_u_bc_bc = u_bc[jk,:,:] * ddx_u_bc + v_bc[jk,:,:] * ddy_u_bc
         adv_u_bt_bt = u_bt[:,:]    * ddx_u_bt + v_bt[:,:]    * ddy_u_bt
         # adv_v = u * dv/dx + v * dv/dy
         adv_v_bc_bc = u_bc[jk,:,:] * ddx_v_bc + v_bc[jk,:,:] * ddy_v_bc
         adv_v_bt_bt = u_bt[:,:]    * ddx_v_bt + v_bt[:,:]    * ddy_v_bt
      else:
         # adv_u = u * du/dx + v * du/dy
         adv_u_bc_bc = adv_u_bc_bc + u_bc[jk,:,:] * ddx_u_bc + v_bc[jk,:,:] * ddy_u_bc
         adv_u_bt_bt = adv_u_bt_bt + u_bt[:,:]    * ddx_u_bt + v_bt[:,:]    * ddy_u_bt
         # adv_v = u * dv/dx + v * dv/dy
         adv_v_bc_bc = adv_v_bc_bc + u_bc[jk,:,:] * ddx_v_bc + v_bc[jk,:,:] * ddy_v_bc
         adv_v_bt_bt = adv_v_bt_bt + u_bt[:,:]    * ddx_v_bt + v_bt[:,:]    * ddy_v_bt
   
   adv_u_bc_bc = adv_u_bc_bc / float(nz)   
   adv_v_bc_bc = adv_v_bc_bc / float(nz)   
   adv_u_bt_bt = adv_u_bt_bt / float(nz)   
   adv_v_bt_bt = adv_v_bt_bt / float(nz)   
   
   # KE trend from advection: 
   # - u * adv_u - v * adv_v
   # in spectral space
   # The minus sign arises as advection 
   # is on the RHS of the momentum eqs. 
   Tk_bt_bc_bc = np.real( -np.conj(uhat_bt)*fftn(adv_u_bc_bc) - \
                           np.conj(vhat_bt)*fftn(adv_v_bc_bc) )   #[m2/s3]
   Tk_bt_bt_bt = np.real( -np.conj(uhat_bt)*fftn(adv_u_bt_bt) - \
                           np.conj(vhat_bt)*fftn(adv_v_bt_bt) )   #[m2/s3]
   Tk_bc_bc_bc = np.real( -np.conj(uhat_bc)*fftn(adv_u_bc_bc) - \
                           np.conj(vhat_bc)*fftn(adv_v_bc_bc) )
   
   nn = (nx**2 * ny**2)                        
   data = {}
   data['Tk_bt_bc_bc'] = Tk_bt_bc_bc / float(nn)
   data['Tk_bt_bt_bt'] = Tk_bt_bt_bt / float(nn)
   data['Tk_bc_bc_bc'] = Tk_bc_bc_bc / float(nn)
   
   return data


def calculate_spectral_flux_baroclinic_barotropic5(kx,ky,u_bt,v_bt,u_bc,v_bc,dzt):
   """
   Calculate spectral flux for a triad
   i.e. u_bt * u_bc * du_bc/dx + u_bt * v_bc * du_bc/dy + 
        v_bt * u_bc * dv_bc/dx + v_bt * v_bc * dv_bc/dy
   
   Used in Kjellsson & Zanna (Fluids, 2017)
   A good reference for the terms is Scott & Arbic (JPO, 2007)
   
   Note: Triads that only include one baroclinic component, 
         e.g. (bt,bt,bc) should be zero. 
         Check that they are! 
   
   Barotropic velocity, (u_bt, v_bt), is 2D
   and baroclinic velocity (u_bc,v_bc) is 3D
   
   Input: 
      kx, ky - 2D arrays with zonal and meridional wavenumbers
      u_bt, v_bt - 2D arrays with barotropic velocity components
      u_bc, v_bc - 3D arrays with baroclinic velocity components
   
   Output:
      data - Dictionary containing 8 different triad interactions
   
   """
   
   i = np.complex(0,1)
   
   uhat_bt = fftn(u_bt)
   vhat_bt = fftn(v_bt)
   
   nx = u_bc.shape[2]
   ny = u_bc.shape[1]
   nz = u_bc.shape[0]
   
   # derivatives for barotropic velocities
   ddx_u_bt = np.real( ifftn(i*kx*uhat_bt) ) # du_bt/dx
   ddy_u_bt = np.real( ifftn(i*ky*uhat_bt) ) # du_bt/dy
   ddx_v_bt = np.real( ifftn(i*kx*vhat_bt) ) # dv_bt/dx
   ddy_v_bt = np.real( ifftn(i*ky*vhat_bt) ) # dv_bt/dy
   
   # Barotropic self-interactions
   # adv_u = u * du/dx + v * du/dy   
   adv_u_bt_bt = u_bt[:,:] * ddx_u_bt + v_bt[:,:] * ddy_u_bt
   # adv_v = u * dv/dx + v * dv/dy
   adv_v_bt_bt = u_bt[:,:] * ddx_v_bt + v_bt[:,:] * ddy_v_bt
   
   adv_u_bc_bc_tot = np.zeros((ny,nx))
   adv_u_bt_bc_tot = np.zeros((ny,nx))
   adv_u_bc_bt_tot = np.zeros((ny,nx))
   adv_v_bc_bc_tot = np.zeros((ny,nx))
   adv_v_bt_bc_tot = np.zeros((ny,nx))
   adv_v_bc_bt_tot = np.zeros((ny,nx))
   
   Tk_bc_bc_bc = np.zeros(uhat_bt.shape)
   Tk_bc_bt_bt = np.zeros(uhat_bt.shape)
   Tk_bc_bc_bt = np.zeros(uhat_bt.shape)
   Tk_bc_bt_bc = np.zeros(uhat_bt.shape)
   
   Hsum = np.sum(dzt,axis=0)
   
   for jk in range(0,nz):
      uhat_bc = fftn(u_bc[jk,:,:])
      vhat_bc = fftn(v_bc[jk,:,:])
      
      ddx_u_bc = np.real( ifftn(i*kx*uhat_bc) ) # du_bc/dx
      ddy_u_bc = np.real( ifftn(i*ky*uhat_bc) ) # du_bc/dy
      ddx_v_bc = np.real( ifftn(i*kx*vhat_bc) ) # dv_bc/dx
      ddy_v_bc = np.real( ifftn(i*ky*vhat_bc) ) # dv_bc/dy
      
      # weight baroclinic advection by layer thickness
      # adv_u = u * du/dx + v * du/dy
      adv_u_bc_bc = (u_bc[jk,:,:] * ddx_u_bc + v_bc[jk,:,:] * ddy_u_bc) * dzt[jk,:,:]/Hsum
      adv_u_bc_bt = (u_bc[jk,:,:] * ddx_u_bt + v_bc[jk,:,:] * ddy_u_bt) * dzt[jk,:,:]/Hsum
      adv_u_bt_bc = (u_bt[:,:]    * ddx_u_bc + v_bt[:,:]    * ddy_u_bc) * dzt[jk,:,:]/Hsum
      # adv_v = u * dv/dx + v * dv/dy
      adv_v_bc_bc = (u_bc[jk,:,:] * ddx_v_bc + v_bc[jk,:,:] * ddy_v_bc) * dzt[jk,:,:]/Hsum
      adv_v_bc_bt = (u_bc[jk,:,:] * ddx_v_bt + v_bc[jk,:,:] * ddy_v_bt) * dzt[jk,:,:]/Hsum
      adv_v_bt_bc = (u_bt[:,:]    * ddx_v_bc + v_bt[:,:]    * ddy_v_bc) * dzt[jk,:,:]/Hsum
      
      # vertical sum of advection terms
      adv_u_bc_bc_tot = adv_u_bc_bc_tot + adv_u_bc_bc
      adv_u_bt_bc_tot = adv_u_bt_bc_tot + adv_u_bt_bc
      adv_u_bc_bt_tot = adv_u_bc_bt_tot + adv_u_bc_bt
      adv_v_bc_bc_tot = adv_v_bc_bc_tot + adv_v_bc_bc
      adv_v_bt_bc_tot = adv_v_bt_bc_tot + adv_v_bt_bc
      adv_v_bc_bt_tot = adv_v_bc_bt_tot + adv_v_bc_bt
      
      # baroclinic budget                                                
      Tk_bc_bc_bc = Tk_bc_bc_bc + np.real( -np.conj(uhat_bc)*fftn(adv_u_bc_bc) - \
                                            np.conj(vhat_bc)*fftn(adv_v_bc_bc) )   
      Tk_bc_bc_bt = Tk_bc_bc_bt + np.real( -np.conj(uhat_bc)*fftn(adv_u_bc_bt) - \
                                            np.conj(vhat_bc)*fftn(adv_v_bc_bt) )   #[m2/s3]
      Tk_bc_bt_bc = Tk_bc_bt_bc + np.real( -np.conj(uhat_bc)*fftn(adv_u_bt_bc) - \
                                            np.conj(vhat_bc)*fftn(adv_v_bt_bc) )   #[m2/s3]
      # weight barotropic advection by layer thickness
      adv_u_bt_bt_new = adv_u_bt_bt * dzt[jk,:,:]/Hsum
      adv_v_bt_bt_new = adv_v_bt_bt * dzt[jk,:,:]/Hsum
      Tk_bc_bt_bt = Tk_bc_bt_bt + np.real( -np.conj(uhat_bc)*fftn(adv_u_bt_bt_new) - \
                                            np.conj(vhat_bc)*fftn(adv_v_bt_bt_new) )   #[m2/s3]
      
   
   # KE trend from advection: 
   # - u * adv_u - v * adv_v
   # The minus sign arises as advection 
   # is on the RHS of the momentum eqs. 
   
   # barotropic budget
   Tk_bt_bc_bc = np.real( -np.conj(uhat_bt)*fftn(adv_u_bc_bc_tot) - \
                           np.conj(vhat_bt)*fftn(adv_v_bc_bc_tot) )   #[m2/s3]
   Tk_bt_bt_bt = np.real( -np.conj(uhat_bt)*fftn(adv_u_bt_bt) - \
                           np.conj(vhat_bt)*fftn(adv_v_bt_bt) )   #[m2/s3]
   Tk_bt_bc_bt = np.real( -np.conj(uhat_bt)*fftn(adv_u_bc_bt_tot) - \
                           np.conj(vhat_bt)*fftn(adv_v_bc_bt_tot) )   #[m2/s3]
   Tk_bt_bt_bc = np.real( -np.conj(uhat_bt)*fftn(adv_u_bt_bc_tot) - \
                           np.conj(vhat_bt)*fftn(adv_v_bt_bc_tot) )   #[m2/s3]
   
   ## The FFT routine in scipy returns
   ## y(j) = (x * exp(-2*pi*sqrt(-1)*j*np.arange(n)/n)).sum()
   ## where j is wavenumber, x is array in gridpoint space, and n is number of grid points
   ## To normalise, we must divide by n
   ## if its a 2D FFT, we must divide by nx*ny
   ## and since we take the square of the transformed variables
   ## we must divide by (nx*ny)^2
   nn2 = (nx**2 * ny**2)                        
   
   data = {}
   data['Tk_bt_bc_bc'] = Tk_bt_bc_bc / float(nn2) #[m2/s3] (energy per second)
   data['Tk_bt_bt_bt'] = Tk_bt_bt_bt / float(nn2)
   data['Tk_bt_bc_bt'] = Tk_bt_bc_bt / float(nn2)
   data['Tk_bt_bt_bc'] = Tk_bt_bt_bc / float(nn2)
   
   data['Tk_bc_bc_bc'] = Tk_bc_bc_bc / float(nn2)
   data['Tk_bc_bc_bt'] = Tk_bc_bc_bt / float(nn2)
   data['Tk_bc_bt_bt'] = Tk_bc_bt_bt / float(nn2)
   data['Tk_bc_bt_bc'] = Tk_bc_bt_bc / float(nn2)
   
   return data


def calculate_spectral_ke_tendency(uvel,vvel,du,dv,win=1.):
   """
   Calculate KE tendency using momentum and momentum tendencies
   """
   
   npoints_sq = uvel.shape[0]**2 * uvel.shape[1]**2
   d_ke = np.real( np.conj( fftn(uvel[:,:]*win) ) * fftn(du[:,:]*win) + \
                   np.conj( fftn(vvel[:,:]*win) ) * fftn(dv[:,:]*win) ) / float(npoints_sq)

   return d_ke
   

def calculate_spectral_ape_flux(kx,ky,uvel,vvel,ape):
   """
   """
   
   i = np.complex(0,1)
   
   nx = uvel.shape[2]
   ny = uvel.shape[1]
   nz = uvel.shape[0]
   
   ## FFT of the square root of APE
   ## We will multiply by conj of APE later, so
   ## we must take square root so that the final 
   ## product is APE
   ahat = fftn(np.sqrt(ape))
   
   for jk in range(0,nz):
      
      uhat = fftn(uvel[jk,:,:])
      vhat = fftn(vvel[jk,:,:])
      
      # dAPE/dx in x,y
      ddx_ape = np.real( ifftn(i*kx*ahat) )
      # dAPE/dy in x,y
      ddy_ape = np.real( ifftn(i*ky*ahat) )
               
      # u * dAPE/dx + v * dAPE/dy
      tmp = uvel[jk,:,:] * ddx_ape + vvel[jk,:,:] * ddy_ape
      adv_tmp = np.real( -np.conj(ahat)*fftn(tmp) )   #[m2/s3]
      if (jk == 0):
         adv_ape = tmp
      else:
         adv_ape += tmp
      
   adv_ape = adv_ape / float(nz)
   nn = (ahat.shape[1]**2 * ahat.shape[0]**2)
   adv_ape = adv_ape / float(nn)
   
   data = {}
   data['adv_ape'] = adv_ape
   
   return data


def calculate_spectral_ape_flux_baroclinic_barotropic(kx,ky,u_bt,v_bt,u_bc,v_bc,ape):
   """
   """
   
   i = np.complex(0,1)
   
   uhat_bt = fftn(u_bt)
   vhat_bt = fftn(v_bt)
   
   nx = u_bt.shape[1]
   ny = u_bt.shape[0]
   nz = u_bc.shape[0]
   
   ## FFT of the square root of APE
   ## We will multiply by conj of APE later, so
   ## we must take square root so that the final 
   ## product is APE
   ahat = fftn(np.sqrt(ape))
   
   # dAPE/dx in x,y
   ddx_ape = np.real( ifftn(i*kx*ahat) )
   # dAPE/dy in x,y
   ddy_ape = np.real( ifftn(i*ky*ahat) )
            
   # u_bt * dAPE/dx + v_bt * dAPE/dy
   adv_bt_ape = u_bt * ddx_ape + v_bt * ddy_ape
   
   for jk in range(0,nz):
      uhat_bc = fftn(u_bc[jk,:,:])
      vhat_bc = fftn(v_bc[jk,:,:])
      # u_bc * dAPE/dx + v_bc * dAPE/dy
      if (jk == 0):
         adv_bc_ape  = u_bc[jk,:,:] * ddx_ape + v_bc[jk,:,:] * ddy_ape
      else:
         adv_bc_ape += u_bc[jk,:,:] * ddx_ape + v_bc[jk,:,:] * ddy_ape
            
   adv_ape_bt_ape = np.real( -np.conj(ahat)*fftn(adv_bt_ape) )   #[m2/s3]
   adv_ape_bc_ape = np.real( -np.conj(ahat)*fftn(adv_bc_ape) )   #[m2/s3]
   nn = (ahat.shape[1]**2 * ahat.shape[0]**2)
   adv_ape_bt_ape = adv_ape_bt_ape / float(nn)
   adv_ape_bc_ape = adv_ape_bc_ape / float(nn)
   
   data = {}
   data['adv_ape_bt_ape'] = adv_ape_bt_ape
   data['adv_ape_bc_ape'] = adv_ape_bc_ape
   
   return data
   

def calculate_spectral_ens_flux(kx,ky,uvel,vvel,rot):
   """
   """
   
   rhat = fftn(rot)
   uhat = fftn(uvel)
   vhat = fftn(vvel)
   i = np.complex(0,1)
   # drot/dx in x,y
   ddx_rot = np.real( ifftn(i*kx*rhat) )
   # drot/dy in x,y
   ddy_rot = np.real( ifftn(i*ky*rhat) )
   
   # adv_rot = u * drot/dx + v * drot/dy
   adv_rot = uvel * ddx_rot + vvel * ddy_rot
   
   # Enstrophy trend from advection: 
   # rot * adv_rot 
   # in spectral space
   Tkxky = np.real( -np.conj(rhat)*fftn(adv_rot) ) 
   
   return Tkxky
   

def calculate_spectral_viscosity(kx,ky,uvel,Ahm,order='4'):
   """
   Calculate spectral flux 
   We assume du/dt = (d4/dx4+d4/dy4) u (order 4) or 
             du/dt = (d2/dx2+d2/dy2) u (order 2) or 
   """
   
   uhat = fftn(uvel)
   i = np.complex(0,1)
   if (order == '4'):
      # d4u/dx4 in x,y
      ddx_u = Ahm * np.real( ifftn(kx**4 * uhat) )
      # d4u/dy4 in x,y
      ddy_u = Ahm * np.real( ifftn(ky**4 * uhat) )
      # 2 d4/dx2dy2 in x,y
      ddxy_u = Ahm * np.real( ifftn(2 * kx**2 * ky**2 * uhat) )
      ##
      Tkxky = np.real( np.conj(uhat) * fftn(ddx_u) + \
                       np.conj(uhat) * fftn(ddxy_u)+ \
                       np.conj(uhat) * fftn(ddy_u)   ) 
                       
   elif (order == '2'):
      # d2u/dx2 in x,y
      ddx_u = -Ahm * np.real( ifftn(kx**2 * uhat) )
      # d2u/dy2 in x,y
      ddy_u = -Ahm * np.real( ifftn(ky**2 * uhat) )
      # KE trend from viscosity: 
      # ddx_u + ddy_u
      # in spectral space
      Tkxky = np.real( np.conj(uhat) * fftn(ddx_u) + \
                       np.conj(uhat) * fftn(ddy_u)   ) 
   
   return Tkxky
   

def calculate_spectral_forcing(kx,ky,uvel,vvel,taux,tauy,rho=1023.,laverage=True):
   """
   Calculate spectral flux from wind forcing
   We assume du/dt = 1/rho*(u_conj * taux + v_conj * tauy)
   """
   
   uhat = fftn(uvel)
   vhat = fftn(vvel)
   txhat = fftn(taux/rho)
   tyhat = fftn(tauy/rho)
   i = np.complex(0,1)
   
   # u_conj * taux
   u_taux = np.conj(uhat) * txhat
   # v_conj * tauy
   v_tauy = np.conj(vhat) * tyhat
   
   # KE trend from wind forcing: 
   Tkxky = np.real( u_taux + v_tauy ) 
   
   nn = uvel.shape[1]**2 * uvel.shape[0]**2
   Tkxky = Tkxky/float(nn)
   
   return Tkxky
   

def integrate_spectrum(psd2D,wvsq,k,dk):
   """
   """
   ## Integrate 2D PSD around lines of constant k
   ## sum(u^2 + v^2) = sum(E)
   ## PSD = d/dk sum(E) [m3/s2]
   ## E = PSD * dk [m2/s2], energy at a given k
   nk = k.shape[0]
   psc = np.zeros((nk))
   for jk in range(0,nk):
      indices = np.where(wvsq >= k[jk]**2)
      psc[jk] = np.sum(psd2D[indices]) 
         
   vpsd   = -(psc[1:] - psc[0:-1]) / dk
   #vpsdk  =  (psc[1:] + psc[0:-1]) * 0.5 
   vpsdk  =  psc[0:-1]
   
   return 0.5*(k[1:]+k[0:-1]),vpsd,vpsdk


def integrate_spectrum2(psd2D,wvsq,k,dk):
   """
   """
   nk = k.shape[0]
   psc_neg = np.zeros((nk))
   psc_pos = np.zeros((nk))
   for jk in range(0,nk):
      indices_neg = np.where( (wvsq >= k[jk]**2) & (psd2D < 0) )
      indices_pos = np.where( (wvsq >= k[jk]**2) & (psd2D > 0) )
      psc_neg[jk] = np.sum(psd2D[indices_neg]) 
      psc_pos[jk] = np.sum(psd2D[indices_pos]) 
         
   vpsd_neg   = -(psc_neg[1:] - psc_neg[0:-1]) / dk
   vpsdk_neg  =  (psc_neg[1:] + psc_neg[0:-1]) * 0.5 
   vpsd_pos   = -(psc_pos[1:] - psc_pos[0:-1]) / dk
   vpsdk_pos  =  (psc_pos[1:] + psc_pos[0:-1]) * 0.5 
   
   return 0.5*(k[1:]+k[0:-1]),vpsd_neg,vpsd_pos,vpsdk_neg,vpsdk_pos


def integrate_spectrum_angle(psd2D,k,l):
   """
   """
   angle = np.linspace(0,np.pi/2.,50)
   wvsq = k**2 + l**2
   nk = angle.shape[0]
   psc = np.zeros((nk))
   for jk in range(0,nk):
      b = np.sqrt(k**2)
      #a = np.arccos(b/np.sqrt(wvsq))
      #a = np.where( c == 0, 1, a )
      indices = np.where( np.arccos(b/np.sqrt(wvsq)) >= angle[jk] )
      psc[jk] = np.sum(psd2D[indices]) 
   
   #box = np.ones(9)/9
   #psc = np.convolve(psc, box, mode='same')
         
   vpsd   = -(psc[1:] - psc[0:-1])
   vpsdk  =  (psc[1:] + psc[0:-1]) * 0.5 
   
   return 0.5*(angle[1:]+angle[0:-1]),vpsd,vpsdk
   

def set_wavenumbers(x0,x1,y0,y1,dx,dy,nx,ny,grid='ll'):
   """
   Set wavenumber arrays given x,y min/max, grid size and resolution
   If grid='ll' (default) it is assumed we have a lon-lat grid
   If grid='xy' we are on a x-y grid
   
   Input:
     x0 - minimum x value (e.g. westernmost longitude)
     x1 - maximum x value (e.g. easternmost longitude)
     y0 - minimum y value (e.g. southern boundary)
     y1 - maximum y value (e.g. northern boundary)
     dx - grid resolution (in degrees if lon-lat or meter if xy)
     dy - meridional resolution
     nx - number of zonal points
     ny - number of meridional points
     
   Output:
     xx, yy - 2D arrays with x,y points on a regular grid
              Can be used if you want to interpolate to this grid
     wn_x, wn_y - 1D arrays of zonal and meridional wavenumbers
     kx, ky - 2D arrays of zonal and meridional wavenumber
     k - 1D array of isotropic wavenumber
     dk - Resolution of k wavenumber array
   
   """
   
   if grid == 'xy':
      ## Set x,y arrays that are monotonic
      xmin = 0
      xmax = x1-x0
      ymin = 0
      ymax = y1-y0
   
      ## Length of domain
      Lx = np.abs(xmax - xmin)
      Ly = np.abs(ymax - ymin)
      
      x = np.arange(x0,x1,nx)
      y = np.arange(y0,y1,ny)
   
   if grid == 'll':
      xmin = x0
      xmax = x1
      ymin = y0
      ymax = y1
      
      # if xmax less than xmin,
      # e.g. xmax = -160 and xmin = 160
      # then add 360 degrees
      if (xmax <= xmin): xmax=xmax+360
      # by checking above, we have 
      # ensured x below is monotonic      
      x = xmin + np.arange(0,nx) * dx
      
      y = ymin + np.arange(0,ny) * dy
      
      ## size of domain in meters
      mean_lat = (y1+y0)/2.
      Lx = 6731000. * (x.max() - x.min()) / 180. * np.pi * np.cos(mean_lat * np.pi/180.)
      Ly = 6731000. * (y.max() - y.min()) / 180. * np.pi
      #print 'Lx,Ly [deg] ',x.max()-x.min(),y.max()-y.min()
      
      # now remove 360 where x > 180 go get 
      # back to original definition
      x = np.where(x>=180,x-360,x)
   
   ## 2D arrays of coordinates   
   xx,yy = np.meshgrid(x,y)
      
   #print 'Lx,Ly [meter]',Lx,Ly
      
   ## Wavenumber vectors
   ## Depend if nx,ny are even or odd numbers
   if (np.mod(nx,2) == 0):
      wn_x = 2.0 * np.pi / Lx * \
             np.concatenate( (np.arange(0,nx/2+1),np.arange(-nx/2+1,0)) )
   else:
      nnx = (nx-1)/2
      wn_x = 2.0 * np.pi / Lx * \
             np.concatenate( (np.arange(0,nnx+1),np.arange(-nnx,0)) )
   if (np.mod(ny,2) == 0):   
      wn_y = 2.0 * np.pi / Ly * \
             np.concatenate( (np.arange(0,ny/2+1),np.arange(-ny/2+1,0)) )       
   else:
      nny = (ny-1)/2
      wn_y = 2.0 * np.pi / Ly * \
             np.concatenate( (np.arange(0,nny+1),np.arange(-nny,0)) )        
   
   ## 2D wavenumber arrays
   kx,ky = np.meshgrid(wn_x,wn_y)
   ## Square total wavenumber
   wvsq  = kx**2 + ky**2
   ## Maximum wavenumber
   wn_max = np.max( np.sqrt(wvsq) )
   
   ## Make an isotropic 1D wavenumber array
   dk = min(2.0*np.pi/Lx,2.0*np.pi/Ly)
   k = np.arange(dk,wn_max+dk,dk)
   nk = k.shape[0]
   
   return xx,yy,wn_x,wn_y,kx,ky,wvsq,k,dk


def smooth_spectrum(data_list,n=20):
   """
   """
   
   data_out = []
   for data in data_list:
      nsmooth = min(data.shape[0],n)
      print np.ones((nsmooth,))/nsmooth
      data = np.convolve(data, np.ones((nsmooth,))/nsmooth, mode='same')
      data_out.append(data)
   
   return data_out
   

def calculate_rho(tem,sal,dep):
   """
   """
   from mod_eos import eos
   eos.tem = tem[:,:,:].T
   eos.sal = sal[:,:,:].T
   eos.dep = dep[:,:,:].T
   tmask = np.where(sal[:,:,:] > 0,1,0)
   eos.tmask = tmask[:,:,:].T
   eos.nn_eos = 0
   eos.eos_init()
   eos.rab_3d()
   alpha = eos.alpha[:,:,:].T
   alpha = np.ma.masked_where( tmask==0, alpha )
   beta  = eos.beta[:,:,:].T
   beta = np.ma.masked_where( tmask==0, beta )
   eos.eos_insitu_pot()
   rho1 = eos.rho[:,:,:].T
   rho1 = np.ma.masked_where( tmask==0, rho1 )
   rho = (rho1+1) * 1026.
   rho = np.ma.masked_where( tmask==0, rho ) 
   eos.clean()
   return rho,rho1


def calculate_means_and_anomalies(uvel,vvel,dzt):      
   """
   """
   
   data = {}
   
   if (1):
      ## uvel_vm - vertical mean
      ## uvel_va - deviations from vertical mean
      uvel_vm = np.mean(uvel,axis=0)
      vvel_vm = np.mean(vvel,axis=0)
      uvel_va = np.zeros(uvel.shape)
      vvel_va = np.zeros(vvel.shape)
      for jk in range(0,uvel.shape[0]):
         uvel_va[jk,:,:] = uvel[jk,:,:] - uvel_vm[:,:]
         vvel_va[jk,:,:] = vvel[jk,:,:] - vvel_vm[:,:]
      
      data['u_vm'] = uvel_vm
      data['v_vm'] = vvel_vm
      data['u_va'] = uvel_va
      data['v_va'] = vvel_va
      
      ## uvel_zm - zonal mean
      uvel_zm = np.mean(uvel,axis=2)
      vvel_zm = np.mean(vvel,axis=2)
      ## uvel_za - deviations from zonal mean
      uvel_za = np.zeros(uvel.shape)
      vvel_za = np.zeros(vvel.shape)
      ## uvel_zm_vm - zonal mean of vertical mean
      uvel_zm_vm = np.mean(uvel_vm,axis=1)
      vvel_zm_vm = np.mean(vvel_vm,axis=1)
      ## uvel_za_vm - zonal deviations of vertical mean
      uvel_za_vm = np.zeros(uvel_vm.shape)
      vvel_za_vm = np.zeros(vvel_vm.shape)
      ## uvel_zm_va - zonal mean of vertical deviations
      uvel_zm_va = np.mean(uvel_va,axis=2)
      vvel_zm_va = np.mean(vvel_va,axis=2)
      ## uvel_za_va - zonal deviations of vertical deviations
      uvel_za_va = np.zeros(uvel.shape)
      vvel_za_va = np.zeros(vvel.shape)
      
      for ji in xrange(0,uvel.shape[2]):
         ## zonal anomalies
         uvel_za[:,:,ji] = uvel[:,:,ji] - uvel_zm[:,:]
         vvel_za[:,:,ji] = vvel[:,:,ji] - vvel_zm[:,:]
         ## zonal anomalies of vertical means (barotropic)
         uvel_za_vm[:,ji] = uvel_vm[:,ji] - uvel_zm_vm[:]
         vvel_za_vm[:,ji] = vvel_vm[:,ji] - vvel_zm_vm[:]
         ## zonal anomalies of vertical anomalies (all baroclinic modes)
         uvel_za_va[:,:,ji] = uvel_va[:,:,ji] - uvel_zm_va[:,:]
         vvel_za_va[:,:,ji] = vvel_va[:,:,ji] - vvel_zm_va[:,:]
      
      data['u_zm'] = uvel_zm
      data['v_zm'] = vvel_zm
      data['u_za'] = uvel_za
      data['v_za'] = vvel_za
      data['u_za_vm'] = uvel_za_vm
      data['v_za_vm'] = vvel_za_vm
      data['u_za_va'] = uvel_za_va
      data['v_za_va'] = vvel_za_va
      
   return data   


def calculate_trend(data_list):
   """
   """
   
   vt   = data_list[0]
   data = data_list[1]

   slope  =  np.ma.zeros((data.shape[1]))
   sign   =  np.ma.zeros(slope.shape)
   
   for i in xrange(0,data.shape[1]):
      ## Masked vector
      tmp = data[:,i]
      
      ## Make vx and vy from valid data
      vx  = vt[:]
      vy  = tmp[:]
      
      k, m, r, p, sig = stats.linregress(vx,vy)
      slope[i] = k
      sign[i]  = p
      
   return slope,sign


def calculate_trend_parallel(data,vt=np.array([])):
   """
   Calculates spatial trends from a masked array of shape (time,x)
   """

   t0 = time.time()

   ## Spawn threads
   ncpu = mp.cpu_count()
   p = mp.Pool(ncpu)
   nump = p._processes

   slope  =  np.ma.zeros((data.shape[1]))
   sign   =  np.ma.zeros(slope.shape)

   if (vt.shape[0] < 1):
      vt = np.arange(0,data.shape[0])

   ##
   ## Split array into ncpu smaller arrays and put
   ## them in a list
   ##
   nx = data.shape[1] ## number of columns in data

   ## Find number of tasks in Pool
   ntask = p._processes

   ## Split array into smaller arrays of size (nt,ny,10)
   ## Each task will make calculations on one array and then go
   ## on the the next
   ##
   ## More sub-arrays gives better load-balance but also more overhead.
   ## For 362 i points, 30-40 sub-arrays seems good.

   s = np.array( np.ceil( np.arange(0,nx,10) ),dtype='int32' ) ##start i for sub arrays
   e = np.array( np.append(s[1:],nx) ,dtype='int32') ## end i for sub arrays

   ## List with all sub-arrays
   l = []
   for j in xrange(0,s.shape[0]):
      l.append([vt,data[:,s[j]:e[j]]])

   ## Let all workers in Pool do calculations asynchronously using routine
   ## calc_spatial_trend_masked
   ##
   ## Results are returned to a list
   results = p.map_async(calculate_trend,l)
   result = results.get()

   ## Put results from sub-arrays into full-size arrays
   for j in xrange(0,s.shape[0]):
      k,m = result[j]
      slope[s[j]:e[j]] = k
      sign[s[j]:e[j]] = m

   p.terminate()
   p.join()

   t1 = time.time()
   print ' calculate_trend_parallel in ',t1-t0,' sec on ',nump,' threads '

   return slope,sign
   

def two_layer(data,dzt,k_interface):
   """
   """
   
   k_upper = np.arange(0,k_interface)
   k_lower = np.arange(k_interface,data.shape[0])
               
   nx = data.shape[2]
   ny = data.shape[1]
   
   data_two = np.zeros((2,ny,nx))
   
   h1 = np.ma.sum(dzt[k_upper,:,:],axis=0)
   h2 = np.ma.sum(dzt[k_lower,:,:],axis=0)
   
   datadz = data[:,:,:]*dzt[:,:,:] 
   data_two[0,:,:] = np.ma.sum( datadz[k_upper,:,:],axis=0) / h1
   data_two[1,:,:] = np.ma.sum( datadz[k_lower,:,:],axis=0) / h2
   
   return data_two,h1,h2
   

def calculate_w(uvel,vvel,dxv,dyu,dzu,dzv,dzt):
   """
   Calculate vertical velocity from continuity
   """
   
   uflux = uvel * dzu * dyu
   vflux = vvel * dzv * dxv
   wflux = np.ma.zeros(uflux.shape)
   for jk in range(nz-1,-1,-1):
      wflux[jk,1:,1:] = wflux[jk+1,1:,1:] -(uflux[jk+1,1:,1:] - uflux[jk+1,1:,0:-1]) - (vflux[jk+1,1:,1:] - vflux[jk+1,0:-1,1:])
      ## add code later to take care of lateral boundaries
      
   return wflux
      
   
         