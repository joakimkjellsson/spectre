import os,sys,time
import numpy as np
from datetime import datetime 
from netCDF4 import Dataset

def read_data(date,run,ddir,mode='data',i0=0,i1=-1,j0=0,j1=-1,levels=np.array([]),readT=False,readW=False,ltend=False):
   """
   Find files with u, v, T, S for the data, and the grid files (lon, lat etc)
   """
   
   data = {} ## load data into a dictionary
   
   if run == 'ORCA1-N406':
      ## ORCA1 simulation by Andrew
      mfile = ddir + '/domain/mask.nc'
      hfile = ddir + '/domain/mesh_hgr.nc'
      zfile = ddir + '/domain/mesh_zgr.nc'
      
      if date.year < 2008:   
         fileprefix = ddir + '/means/%04d/ORCA1-N406_%04d%02d%02dd05' % (date.year,date.year,date.month,date.day)
      elif date.year >= 2008:
         fileprefix = ddir + '/means/%04db/ORCA1-N406_%04d%02d%02dd05' % (date.year,date.year,date.month,date.day)
      
      ufile = fileprefix + 'U.nc'
      vfile = fileprefix + 'V.nc'
      tfile = fileprefix + 'T.nc'
      wfile = fileprefix + 'W.nc'
      
      data['Ahm0'] = 10000. # momentum viscosity
      data['visc_form'] = 'laplacian' # formulation of viscosity
      data['Aht0'] = 1000.  # tracer diffusion
      data['diff_form'] = 'laplacian' # formulation of tracer diffusion
   
      
   if run == 'ORCA025-N401':
      ## ORCA025 simulation by Andrew
      mfile = ddir + '/domain/mask.nc'
      hfile = ddir + '/domain/mesh_hgr.nc'
      zfile = ddir + '/domain/mesh_zgr.nc'
      
      fileprefix = ddir + '/means/%04d/ORCA025-N401_%04d%02d%02dd05' % (date.year,date.year,date.month,date.day)
      
      ufile = fileprefix + 'U.nc'
      vfile = fileprefix + 'V.nc'
      tfile = fileprefix + 'T.nc'
      wfile = fileprefix + 'W.nc'
      
      data['Ahm0'] = -2.2 * 10**11 # momentum viscosity
      data['visc_form'] = 'bilaplacian' # formulation of viscosity
      data['Aht0'] = 300.  # tracer diffusion
      data['diff_form'] = 'laplacian' # formulation of tracer diffusion
      
   if run == 'ORCA0083-N001':
      ## ORCA0083 simulation by Andrew
      mfile = ddir + '/domain/mask.nc'
      hfile = ddir + '/domain/mesh_hgr.nc'
      zfile = ddir + '/domain/mesh_zgr.nc'
      
      fileprefix = ddir + '/means/%04d/ORCA0083-N01_%04d%02d%02dd05' % (date.year,date.year,date.month,date.day)
      
      ufile = fileprefix + 'U.nc'
      vfile = fileprefix + 'V.nc'
      tfile = fileprefix + 'T.nc'
      wfile = fileprefix + 'W.nc'
      
      data['Ahm0'] = -1.25 * 10**10 # momentum viscosity
      data['visc_form'] = 'bilaplacian' # formulation of viscosity
      data['Aht0'] = 125.  # tracer diffusion
      data['diff_form'] = 'laplacian' # formulation of tracer diffusion
   
   if run == 'ORCA0083-N001-CG':
      ## ORCA0083 by Andrew but coarse-grained onto ORCA025 grid
      mfile = ddir + '/mask.nc'
      hfile = ddir + '/mesh_hgr.nc'
      zfile = ddir + '/mesh_zgr.nc'
      
      fileprefix = ddir + '/means/%04d/ORCA0083-N01_%04d%02d%02dd05' % (date.year,date.year,date.month,date.day)
      
      ufile = fileprefix + 'U.nc'
      vfile = fileprefix + 'V.nc'
      tfile = fileprefix + 'T.nc'
      wfile = fileprefix + 'W.nc'
      
      data['Ahm0'] = -1.25 * 10**10 # momentum viscosity
      data['visc_form'] = 'bilaplacian' # formulation of viscosity
      data['Aht0'] = 125.  # tracer diffusion
      data['diff_form'] = 'laplacian' # formulation of tracer diffusion
   
   if run[0:12] == 'GYRE4_square':
      ## My double-gyre runs at 1/4 resolution
      mfile = ddir + '/mask.nc'
      hfile = ddir + '/mesh_hgr.nc'
      zfile = ddir + '/mesh_zgr.nc'
      
      fileprefix = ddir + '/SAVED/'+run+'/'+run+'_5d_%04d0101_%04d1230_grid_' % (date.year,date.year)
      
      ufile = fileprefix + 'U.nc'
      vfile = fileprefix + 'V.nc'
      tfile = fileprefix + 'T.nc'
      wfile = fileprefix + 'W.nc'
      
      data['Ahm0'] = -2.2 * 10**11 # momentum viscosity
      data['visc_form'] = 'bilaplacian' # formulation of viscosity
      data['Aht0'] = 300.  # tracer diffusion
      data['diff_form'] = 'laplacian' # formulation of tracer diffusion
      
   if run[0:13] == 'GYRE12_square':
      ## My double-gyre runs at 1/12 resolution
      mfile = ddir + '/mask.nc'
      hfile = ddir + '/mesh_hgr.nc'
      zfile = ddir + '/mesh_zgr.nc'
      
      fileprefix = ddir + '/SAVED/'+run+'/'+run+'_5d_%04d0101_%04d1230_grid_' % (date.year,date.year)
      
      ufile = fileprefix + 'U.nc'
      vfile = fileprefix + 'V.nc'
      tfile = fileprefix + 'T.nc'
      wfile = fileprefix + 'W.nc'
      
      data['Ahm0'] = -1.25 * 10**10 # momentum viscosity
      data['visc_form'] = 'bilaplacian' # formulation of viscosity
      data['Aht0'] = 125.  # tracer diffusion
      data['diff_form'] = 'laplacian' # formulation of tracer diffusion
   
   if run == 'ORCA05.L46-KJH0004':
      ## Global ORCA05 run by Jan H
      mfile = ddir + '/masks/mask.nc'
      hfile = ddir + '/masks/mesh_hgr.nc'
      zfile = ddir + '/masks/mesh_zgr.nc'
      
      fileprefix = ddir + '/ORCA05.L46-KJH0004_5d_%04d0101_%04d1231_grid_' % (date.year,date.year)
      
      ufile = fileprefix + 'U.nc'
      vfile = fileprefix + 'V.nc'
      tfile = fileprefix + 'T.nc'
      wfile = fileprefix + 'W.nc'
      
      data['Ahm0'] = 0.
      data['visc_form'] = 'bilaplacian'
      data['Aht0'] = 0.
      data['diff_form'] = 'laplacian'
   
   if run == 'INALT10.L46-KJH0017':
      ## Global ORCA05 run with AGRIF nesting (INALT10) by Jan H
      ## This is the data on the ORCA05 grid
      mfile = ddir + '/masks/mask.nc'
      hfile = ddir + '/masks/mesh_hgr.nc'
      zfile = ddir + '/masks/mesh_zgr.nc'
      
      fileprefix = ddir + '/INALT10.L46-KJH0017_5d_%04d0101_%04d1231_grid_' % (date.year,date.year)
      
      ufile = fileprefix + 'U.nc'
      vfile = fileprefix + 'V.nc'
      tfile = fileprefix + 'T.nc'
      wfile = fileprefix + 'W.nc'
      
      data['Ahm0'] = 0.
      data['visc_form'] = 'bilaplacian'
      data['Aht0'] = 0.
      data['diff_form'] = 'laplacian'
   
   if run == 'INALT10.L46-KJH0017-NEST1':
      ## Global ORCA05 run with AGRIF nesting (INALT10) by Jan H
      ## This is the data in the nest
      mfile = ddir + '/masks/1_mask.nc'
      hfile = ddir + '/masks/1_mesh_hgr.nc'
      zfile = ddir + '/masks/1_mesh_zgr.nc'
      
      fileprefix = ddir + '/1_INALT10.L46-KJH0017_5d_%04d0101_%04d1231_grid_' % (date.year,date.year)
      
      ufile = fileprefix + 'U.nc'
      vfile = fileprefix + 'V.nc'
      tfile = fileprefix + 'T.nc'
      wfile = fileprefix + 'W.nc'
      
      data['Ahm0'] = -1.25 * 10**11
      data['visc_form'] = 'bilaplacian'
      data['Aht0'] = 0.
      data['diff_form'] = 'laplacian'
      
      
   if not readT: tfile = None
   if not readW: wfile = None
      
   ##
   ## We use same function to read all NEMO runs
   ## Just make sure filenames are set before
   ## 
   if (run[0:4] == 'ORCA') or (run[0:4] == 'GYRE') or (run[0:5] == 'INALT'):
      if mode == 'grid':
         data_read = read_nemo_grid(mfile,hfile,zfile,i0=i0,i1=i1,j0=j0,j1=j1,levels=levels)
      elif mode == 'data':
         data_read = read_nemo_data(ufile,vfile,tfile=tfile,wfile=wfile,\
                                    i0=i0,i1=i1,j0=j0,j1=j1,levels=levels,ltend=ltend)   
      ## merge settings above with the data we just read
      z = data.update(data_read)
         
   ##
   ## Read AVISO data
   ##
   if run == 'AVISO-sla':
      ## Using anomaly data (mean removed)
      fileprefix = ddir + '/dt_global_allsat_msla_%s_%04d%02d%02d.nc' 
      ufile = fileprefix % ('uv',date.year,date.month,date.day)
      hfile = fileprefix % ('h',date.year,date.month,date.day)
      
   if run == 'AVISO-adt':
      ## Using absolute data (mean included)
      fileprefix = ddir + '/dt_global_allsat_madt_%s_%04d%02d%02d.nc' 
      ufile = fileprefix % ('uv',date.year,date.month,date.day)
      hfile = fileprefix % ('h',date.year,date.month,date.day)
      
   if run[0:5] == 'AVISO':
      if mode == 'grid':
         data = read_aviso_grid()
      elif mode == 'data':
         data = read_aviso_data(ufile,hfile,i0=i0,i1=i1,j0=j0,j1=j1)
      
   ## 
   ## Read MITgcm data (for example)
   ##
   if run == 'MITgcm-test':
      if mode == 'grid':
         data = read_mitgcm_grid()
      elif mode == 'data':
         data = read_mitgcm_data()
      
      
   ## return data
   return data


def read_nemo_grid(mfile,hfile,zfile,i0=0,i1=-1,j0=0,j1=-1,levels=np.array([])):
   """
   Read the NEMO grid 
   
   Here we make a dictionary, data, and populate it with the following
   tlon, tlat - lon, lat at T points
   ulon, ulat - lon, lat at U points
   vlon, vlat - lon, lat at V points
   
   dxu, dyu - dx,dy at U points
   dxv, dyv - dx,dy at V points
   dxt, dyt - dx,dy at T points
   
   nx, ny, nz - number of zonal, meridional and vertical points
   kmt - Number of valid vertical points at each i,j point
   tmask, umask, vmask - Land-sea mask at T, U and V points
   
   dzt - Layer thickness at each i,j,k point
   dzt_1d - Mean layer thickness for each vertical layer
   dept - Depth at each i,j,k T point
   dept_1d - Mean depth at each vertical T point
   depw_1d - Mean depth at each vertical W point
   
   """
   
   print hfile
   data = {}
   nc = Dataset(hfile,'r')
   if i1 == -1:
      i1 = len(nc.dimensions['x'])
   if j1 == -1:
      j1 = len(nc.dimensions['y'])
   if levels.shape[0] == 0:
      levels = np.arange(0,len(nc.dimensions['z']))
   
   data['tlon'] = nc.variables['glamt'][0,j0:j1,i0:i1]
   data['tlat'] = nc.variables['gphit'][0,j0:j1,i0:i1]
   data['ulon'] = nc.variables['glamu'][0,j0:j1,i0:i1]
   data['ulat'] = nc.variables['gphiu'][0,j0:j1,i0:i1]
   data['vlon'] = nc.variables['glamv'][0,j0:j1,i0:i1]
   data['vlat'] = nc.variables['gphiv'][0,j0:j1,i0:i1]
   data['dxu']  = nc.variables['e1u'][0,j0:j1,i0:i1]
   data['dxv']  = nc.variables['e1v'][0,j0:j1,i0:i1]
   data['dxt']  = nc.variables['e1t'][0,j0:j1,i0:i1]
   data['dxf']  = nc.variables['e1f'][0,j0:j1,i0:i1]
   data['dyu']  = nc.variables['e2u'][0,j0:j1,i0:i1]
   data['dyv']  = nc.variables['e2v'][0,j0:j1,i0:i1]
   data['dyt']  = nc.variables['e2t'][0,j0:j1,i0:i1]
   data['dyf']  = nc.variables['e2f'][0,j0:j1,i0:i1]
   nc.close()
   
   print mfile 
   nc = Dataset(mfile,'r')
   data['tmask'] = nc.variables['tmask'][0,levels,j0:j1,i0:i1]
   data['umask'] = nc.variables['umask'][0,levels,j0:j1,i0:i1]
   data['vmask'] = nc.variables['vmask'][0,levels,j0:j1,i0:i1]
   nc.close()
   
   data['nz'] = data['tmask'].shape[0]
   data['ny'] = data['tmask'].shape[1]
   data['nx'] = data['tmask'].shape[2]
   
   print zfile
   nc = Dataset(zfile,'r')
   data['kmt'] = nc.variables['mbathy'][0,j0:j1,i0:i1]
   
   ## if layer thickness is saved as a 3D field
   if ('e3t' in nc.variables.keys()):
      data['dzt'] = nc.variables['e3t'][0,levels,j0:j1,i0:i1]
   elif ('e3t_0' in nc.variables.keys()):
      data['dzt'] = nc.variables['e3t_0'][0,levels,j0:j1,i0:i1]   
         
   ## if not saved, we use the 1D dz to construct 3D field
   elif 'e3t_1d' in nc.variables.keys():
      data['dzt'] = np.zeros((levels.shape[0],kmt.shape[0],kmt.shape[1]))
      for jk in range(0,data['dzt'].shape[0]):
         data['dzt'][jk,:,:] = nc.variables['e3t_1d'][0,levels[jk]]
   
   else:
      print ' Can not find or construct dz at T points '
      sys.exit()
      
   if ('gdept_0' in nc.variables.keys()) and (nc.variables['gdept_0'].ndim == 4):
      data['dept']  = nc.variables['gdept_0'][0,levels,j0:j1,i0:i1]
   elif 'gdept_1d' in nc.variables.keys():
      data['dept'] = np.zeros((levels.shape[0],j1-j0,i1-i0))
      print 
      for jk in range(0,data['dept'].shape[0]):
         data['dept'][jk,:,:] = nc.variables['gdept_1d'][0,levels[jk]]
   else:
      print ' Can not find or construct depth at T points '
      sys.exit()
      
   if ('gdepw_0' in nc.variables.keys()) and (nc.variables['gdepw_0'].ndim == 4):
      data['depw']  = nc.variables['gdepw_0'][0,levels,j0:j1,i0:i1]
   elif 'gdepw_1d':
      data['depw'] = np.zeros((levels.shape[0],j1-j0,i1-i0))
      for jk in range(0,data['depw'].shape[0]):
         data['depw'][jk,:,:] = nc.variables['gdepw_1d'][0,levels[jk]]
   else:
      print ' Can not find or construct dz at W points '
      sys.exit()      
   
   ## make 1D arrays of dz and depth 
   if 'gdept_1d' in nc.variables.keys():
      data['dept_1d']  = nc.variables['gdept_1d'][0,levels]
   elif 'dept' in data.keys():
      dept = np.ma.masked_where(data['tmask'] == 0, data['dept'])
      data['dept_1d'] = np.ma.mean( np.ma.mean(dept,axis=1), axis=1 )
   else: 
      print 'Could not find a 1D depth array'
      sys.exit()
   
   if 'gdepw_1d' in nc.variables.keys():
      data['depw_1d']  = nc.variables['gdepw_1d'][0,levels]
   elif 'depw' in data.keys():
      depw = np.ma.masked_where(data['tmask'] == 0, data['depw'])
      data['depw_1d'] = np.ma.mean( np.ma.mean(depw,axis=1), axis=1 )
   else:
      print 'Could not find a 1D depth array'
      sys.exit()
   
   if 'e3t_1d' in nc.variables.keys():
      data['dzt_1d']  = nc.variables['e3t_1d'][0,levels]
   elif 'dzt' in data.keys():
      dzt = np.ma.masked_where(data['tmask'] == 0, data['dzt'])
      data['dzt_1d'] = np.ma.mean( np.ma.mean(dzt,axis=1), axis=1 )
   else:
      print 'Could not find a 1D dz array'
      sys.exit()
   
   nc.close()
      
   return data


def read_nemo_data(ufile,vfile,tfile=None,wfile=None,ltend=False,i0=0,i1=-1,j0=0,j1=-1,\
                   levels=np.array([]),step0=0,step1=-1):
   """
   Read NEMO model output
   
   Here we make a dictionary, data, and fill it with the following:
   
   uvel, vvel - Zonal and meridional velocities
   wvel - Vertical velocity if wfile exists
   tem, sal - Temperature and salinity (if tfile exists)
   ssh - Sea surface height if t file exists
   taux, tauy - Zonal and meridional wind stress 
   
   tlon, tlat - Lon,lat for T points (and similar for U,V points)
   
   utend_adv, vtend_adv - Online tendencies from momentum advection scheme in NEMO
                          and similar for other tendencies, e.g. visc, hpg, spg
                          (note that utend_zdf is for vertical viscosity and includes
                           wind stress and bottom friction!)
   
   """
   
   print ufile
   ncu = Dataset(ufile,'r')
   ncv = Dataset(vfile,'r')
   if wfile != None:
      ncw = Dataset(wfile,'r')
   if tfile != None:
      nct = Dataset(tfile,'r')
   
   if i1 == -1:
      i1 = len(ncu.dimensions['x'])
   if j1 == -1:
      j1 = len(ncu.dimensions['y'])
   if levels.shape[0] == 0:
      levels = np.arange(0,len(ncu.dimensions['depthu']))
   if step1 == -1:
      step1 = len(ncu.dimensions['time_counter'])
         
   ulon = ncu.variables['nav_lon'][j0:j1,i0:i1]
   ulat = ncu.variables['nav_lat'][j0:j1,i0:i1]
   vlon = ncu.variables['nav_lon'][j0:j1,i0:i1]
   vlat = ncu.variables['nav_lat'][j0:j1,i0:i1]
   
   ## In some NEMO outputs, there is a fillValue = 0
   ## The netCDF library assumes values at or below fillValue are invalid
   ## So we must remove the mask and replace it with the land mask!
   ## We can (safely?) assume that if u=0 and v=0, then it is land
   ## or if abs(u) > 20, or abs(v) > 20 then it is also land. 
   
   ## check if u is named u or vozocrtx etc
   unames = ['u','vozocrtx']
   vnames = ['v','vomecrty']
   tnames = ['tem','votemper']
   snames = ['sal','vosaline']
   enames = ['ssh','sossheig']
   
   ## test if one of the entries in list unames is in the ufile
   match = set(unames).intersection(ncu.variables.keys()); 
   if (len(match) > 0): 
      uname=match.pop() 
   else: 
      print 'could not find zonal velocity in file ',ufile,unames
   
   match = set(vnames).intersection(ncv.variables.keys()); 
   if (len(match) > 0): 
      vname=match.pop() 
   else: 
      print 'could not find meridional velocity in file ',vfile,vnames
   
   if tfile != None:
      match = set(tnames).intersection(nct.variables.keys()); 
      if (len(match) > 0): 
         tname=match.pop() 
      else: 
         print 'could not find temperature in file ',tfile,tnames
      
      match = set(snames).intersection(nct.variables.keys()); 
      if (len(match) > 0): 
         sname=match.pop() 
      else: 
         print 'could not find salinity in file ',tfile,snames
      
      match = set(enames).intersection(nct.variables.keys()); 
      if (len(match) > 0): 
         ename=match.pop() 
      else: 
         print 'could not find ssh in file ',tfile,enames
   
   ## read u,v, and wind stress
   uvel_full = np.array(ncu.variables[uname][step0:step1,levels,j0:j1,i0:i1])
   vvel_full = np.array(ncv.variables[vname][step0:step1,levels,j0:j1,i0:i1])
   taux_full = ncu.variables['sozotaux'][step0:step1,j0:j1,i0:i1]
   tauy_full = ncv.variables['sometauy'][step0:step1,j0:j1,i0:i1]
   
   if wfile != None:
      tlon = ncw.variables['nav_lon'][j0:j1,i0:i1]
      tlat = ncw.variables['nav_lat'][j0:j1,i0:i1]
      tlon = np.array(tlon)
      tlat = np.array(tlat)
      wvel_full = ncw.variables['vovecrtz'][step0:step1,levels,j0:j1,i0:i1]
   
   if tfile != None:
      print ' Read T,S '
      tlon = nct.variables['nav_lon'][j0:j1,i0:i1]
      tlat = nct.variables['nav_lat'][j0:j1,i0:i1]
      tlon = np.array(tlon)
      tlat = np.array(tlat)
      tem_full = nct.variables[tname][step0:step1,levels,j0:j1,i0:i1]
      sal_full = nct.variables[sname][step0:step1,levels,j0:j1,i0:i1]
      ssh_full = nct.variables[ename][step0:step1,j0:j1,i0:i1]
      tmask = np.where(sal_full <= 0,0,1)
   
   ## read GM eddy induced velocities if available
   if ('vozoeivu' in ncu.variables.keys()):
      uvel2 = np.array(ncu.variables['vozoeivu'][step0:step1,levels,j0:j1,i0:i1])
      if isinstance(uvel2,np.ma.MaskedArray):
         uvel2[uvel2.mask] = 0.
      uvel_full = uvel_full + uvel2 
      print ' Added GM u velocities '
   if ('vomeeivv' in ncv.variables.keys()):
      vvel2 = np.array(ncv.variables['vomeeivv'][step0:step1,levels,j0:j1,i0:i1])
      if isinstance(vvel2,np.ma.MaskedArray):
         vvel2[vvel2.mask] = 0.
      vvel_full = vvel_full + vvel2
      print ' Added GM v velocities '
   
   ## sometimes uvel, vvel are undefined
   ## mask all values larger than 20m/s
   umax = max( np.abs(uvel_full.min()), uvel_full.max() )
   vmax = max( np.abs(vvel_full.min()), vvel_full.max() )
   if (umax > 10 or vmax > 10):
      print ' WARNING!!! uvel > 10 or vmax > 10'
      
   umask = np.where( (np.abs(uvel_full) > 20) | (uvel_full == 0.), 1, 0)
   vmask = np.where( (np.abs(vvel_full) > 20) | (vvel_full == 0.), 1, 0)
   uvel_full = np.ma.masked_where( umask == 1, uvel_full )
   vvel_full = np.ma.masked_where( vmask == 1, vvel_full )
   if wfile != None:
      wmask = np.where( np.abs(wvel_full) > 20, 1, 0)
      wvel_full = np.ma.array( wvel_full, mask=wmask )
      wmax = max( np.abs(wvel_full.min()), wvel_full.max() )
   taux_full = np.ma.array( taux_full, mask=umask[:,0,:,:] )
   tauy_full = np.ma.array( tauy_full, mask=vmask[:,0,:,:] )
   umax = max( np.abs(uvel_full.min()), uvel_full.max() )
   vmax = max( np.abs(vvel_full.min()), vvel_full.max() )
   if (umax > 10 or vmax > 10):
      print ' uvel > 10 or vmax > 10 or wmax > 10'
      sys.exit()
      
   #uvel_full = np.ma.masked_where( (uvel_full == 0) & (vvel_full == 0), uvel_full ) 
   #vvel_full = np.ma.masked_where( (uvel_full == 0) & (vvel_full == 0), vvel_full )
   #taux_full = np.ma.masked_where( (uvel_full[:,0,:,:] == 0) & (vvel_full[:,0,:,:] == 0), taux_full )
   #tauy_full = np.ma.masked_where( (uvel_full[:,0,:,:] == 0) & (vvel_full[:,0,:,:] == 0), tauy_full )
         
   data = {}
   data['nt'] = len(ncu.dimensions['time_counter'])
   data['ulon'] = ulon
   data['ulat'] = ulat
   data['vlon'] = vlon
   data['vlat'] = vlat
   data['uvel'] = uvel_full
   data['vvel'] = vvel_full
   if wfile != None:
      data['tlon'] = tlon
      data['tlat'] = tlat
      data['wvel'] = wvel_full
   data['taux'] = taux_full
   data['tauy'] = tauy_full
   if tfile != None:
      data['ssh']  = ssh_full
      data['tlon'] = tlon
      data['tlat'] = tlat
      data['tem'] = tem_full
      data['sal'] = sal_full
   
   if ltend:
      print ' === Read NEMO online tendencies === '
      unc = Dataset(utfile,'r')
      ## advection is rvo (vort x (u,v)) + KE gradient
      ## ldf is horizontal viscosity
      ## zdf is wind + bottom friction + vertical viscosity
      ## cor is planetary vorticity terms
      ## pre is surface pressure gradient + hydrostatic pressure gradient
      data['utend_adv'] = unc.variables['utrd_rvo'][step0:step1,klevels,j0:j1,i0:i1] + \
                          unc.variables['utrd_keg'][step0:step1,klevels,j0:j1,i0:i1]
      data['utend_vsc'] = unc.variables['utrd_ldf'][step0:step1,klevels,j0:j1,i0:i1]
      data['utend_zdf'] = unc.variables['utrd_zdf'][step0:step1,klevels,j0:j1,i0:i1]
      data['utend_cor'] = unc.variables['utrd_pvo'][step0:step1,klevels,j0:j1,i0:i1]
      data['utend_pre'] = unc.variables['utrd_spg'][step0:step1,klevels,j0:j1,i0:i1] +\
                          unc.variables['utrd_hpg'][step0:step1,klevels,j0:j1,i0:i1]
      data['utend_spg'] = unc.variables['utrd_spg'][step0:step1,klevels,j0:j1,i0:i1]
      data['utend_hpg'] = unc.variables['utrd_hpg'][step0:step1,klevels,j0:j1,i0:i1]
      
      data['vtend_adv'] = vnc.variables['vtrd_rvo'][step0:step1,klevels,j0:j1,i0:i1] + \
                          vnc.variables['vtrd_keg'][step0:step1,klevels,j0:j1,i0:i1]
      data['vtend_vsc'] = vnc.variables['vtrd_ldf'][step0:step1,klevels,j0:j1,i0:i1]
      data['vtend_zdf'] = vnc.variables['vtrd_zdf'][step0:step1,klevels,j0:j1,i0:i1]
      data['vtend_cor'] = vnc.variables['vtrd_pvo'][step0:step1,klevels,j0:j1,i0:i1]
      data['vtend_pre'] = vnc.variables['vtrd_spg'][step0:step1,klevels,j0:j1,i0:i1] +\
                          vnc.variables['vtrd_hpg'][step0:step1,klevels,j0:j1,i0:i1]
      data['vtend_spg'] = vnc.variables['vtrd_spg'][step0:step1,klevels,j0:j1,i0:i1]
      data['vtend_hpg'] = vnc.variables['vtrd_hpg'][step0:step1,klevels,j0:j1,i0:i1]
   
   ncu.close()
   ncv.close()
   if wfile != None:
      ncw.close()
   if tfile != None:
      nct.close()
   
   return data


def read_aviso(ufile,hfile,i0=0,i1=-1,j0=0,j1=-1):
   """
   Read AVISO data 
   """
   
   print ufile
   nc = Dataset(ufile,'r')
   lon = nc.variables['lon'][i0:i1]
   lon = np.where( lon>180, lon-360, lon )
   lat = nc.variables['lat'][j0:j1]
   uvel = nc.variables['u'][steps,j0:j1,i0:i1]
   vvel = nc.variables['v'][steps,j0:j1,i0:i1]
   nc.close()
   
   lon,lat = np.meshgrid(lon,lat)
   
   ## sometimes uvel, vvel are undefined
   ## mask all values larger than 20m/s
   umax = max( np.abs(uvel.min()), uvel.max() )
   vmax = max( np.abs(vvel.min()), vvel.max() )
   if (umax > 10 or vmax > 10):
      print ' WARNING!!! uvel > 10 or vmax > 10'
      
   umask = np.where( np.abs(uvel) > 20, 1, 0)
   vmask = np.where( np.abs(vvel) > 20, 1, 0)
   uvel = np.ma.array( uvel, mask=umask )
   vvel = np.ma.array( vvel, mask=vmask )
   umax = max( np.abs(uvel.min()), uvel.max() )
   vmax = max( np.abs(vvel.min()), vvel.max() )
   if (umax > 10 or vmax > 10):
      print ' uvel > 10 or vmax > 10'
      sys.exit()
   
   data['tlon'] = lon
   data['tlat'] = lat
   data['uvel'] = uvel
   data['vvel'] = vvel
   
   return data


def interpolate_alldata(data_list,x2D,y2D,x1D,y1D):
   """
   Given a list of variables
   loop through all variables and interpolate to a regular grid
   given by lon,lat
   
   Input:
      data_list - List where each element is an array
      x2D, y2D  - 2D arrays with x and y points 
      x1D, y1D  - 1D arrays for regular grid
   
   Output:
      data_out  - List where each element is an array
   
   HUGE WARNING: The function interp uses the griddata function in matplotlib!
                 griddata uses nearest-neighbour interpolation. 
                 Thus, it will fill in some missing values, and it is not energy
                 conserving! 
                 You might be adding or removing energy with this scheme!
   
   To do: Write an energy conserving re-gridding routine! 
   
   """
   
   data_out = []
   
   for data in data_list:
      ## get info about the data
      ## supported shapes are (y,x), (t,y,x), (z,y,x), (t,z,y,x)
      
      nd = data.ndim
      if nd == 2:
         ny,nx = data.shape
         nz = 1
         nt = 1
      elif nd == 3:
         nz,ny,nx = data.shape
         nt = 1
      elif nd == 4:
         nt,nz,ny,nx = data.shape
      else:
         print ' shape of data not supported by interpolate_alldata function '
         print ' The data should be 2D, 3D or 4D and the last two dimensions ' 
         print ' should be y and x. '
         print ' shape of data: ',data.shape
         sys.exit()
      
      ## new array with regular gridded data
      data_reg = np.ma.zeros((nt,nz,ny,nx))
      ## new array with irregular data
      data_irr = np.ma.zeros((nt,nz,ny,nx))
      data_irr[:] = data[:].copy()
      
      t0 = time.time()
      for jn in range(0,nt):
         for jk in range(0,nz):
            lst = [x2D, y2D, data_irr[jn,jk,:,:], xx2, yy2]
            data_reg[jn,jk,:,:] = interp(lst)
      t1 = time.time()
      
      if nd == 2:
         data_reg = data_reg[0,0,:,:]
      elif nd == 3:
         data_reg = data_reg[0,:,:,:]
      elif (nd > 4) or (nd < 1):
         print ' not supported shape '
         print data_reg.shape,data.shape
         sys.exit() 
      
      data_out.append(data_reg)
      
   return data_out



   