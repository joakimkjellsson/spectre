##
##
##  Analyse kinetic and potential energy in ocean models
##  as a function of horizontal wavenumber
##
##  Borrows shamlessly from scripts by Rob Scott
##
##  Joakim "Definitely not the real slim shady" Kjellsson, GEOMAR, Jan 2018
##
##  
##  To do:
##    * Fix APE code. It gives odd results sometimes. 
##      Must calculate N2 and buoyancy as in NEMO for it to work! 
##
##

import os,sys,time
import numpy as np

## turns off figures popping up on screen
#import matplotlib
#matplotlib.use('Agg')  

## load neat functions for calendar and time arrays
from datetime import datetime, timedelta
from calendar import monthrange, isleap

## load fft and window function routines
from scipy.fftpack import fft, ifft, fftn, ifftn
from scipy.signal import periodogram, hamming, tukey

## load stats package (for linear regs etc.)
import scipy.stats as stats

## load interpolation routine (not used?)
from scipy.interpolate import griddata

## load plot package
## and added functions (map plotting, etc)
import matplotlib.pyplot as plt
from matplotlib.mlab import stineman_interp
import matplotlib.ticker as mtick
from mpl_toolkits.basemap import Basemap

## load netCDF package
from netCDF4 import Dataset, num2date, date2num

## load my own python functions
from mod_eos import *
from read_data import * 
from read_saved_data import *
from write_data import *
from psd_functions import *  ## function to do spectral calculations
from plot_functions_2 import * ## plots

## load my own cython code
## compile cython code before loading functions
import pyximport
pyximport.install(setup_args={"include_dirs":np.get_include()},
                  reload_support=True)

from psd_functions_cython import *  ## cython code to do some quick loops
         

##
## Define some local functions
##
def make_file(outfile,currentTime,nk,fprec='f4',iprec='i4',cal='noleap'):
   """
   Create an output stream in netCDF
   """
   
   print outfile
   
   nc = Dataset(outfile,'w')
   nc.createDimension('x',grid['nx'])
   nc.createDimension('y',grid['ny'])
   nc.createDimension('z',klevels_keep.shape[0])
   nc.createDimension('k',nk)
   nc.createDimension('t',None)

   ## wavenumber and time variable
   idx = nc.createVariable('tlon',fprec,('x',))
   idy = nc.createVariable('tlat',fprec,('y',))
   idz = nc.createVariable('deptht',fprec,('z',))
   idk = nc.createVariable('k',fprec,('k',))
   idt = nc.createVariable('t',fprec,('t',))
   
   iyear = nc.createVariable('year',fprec,('t',))
   imon  = nc.createVariable('month',fprec,('t',))
   iday  = nc.createVariable('day',fprec,('t',))
   itime = nc.createVariable('time',fprec,('t',))
   itime.units = 'hours since 1850-01-01'
   itime.calendar = cal
   
   ikint = nc.createVariable('interface_depth',fprec,('t',))
   
   ## cross sections
   nc.createVariable('cross_uvel',fprec,('t','z','y'))
   nc.createVariable('cross_rho',fprec,('t','z','y'))
   
   ## KE 
   ipsd_ke = nc.createVariable('psd_ke',fprec,('t','k'))
   ipst_ke = nc.createVariable('psdt_ke',fprec,('t'))
   itt_ke  = nc.createVariable('psd_ke_timescale',fprec,('t','k'))
   
   ## Spectral transfer/flux of KE
   iTk_ke = nc.createVariable('Tk_ke',fprec,('t','k'))
   iPi_ke = nc.createVariable('Pi_ke',fprec,('t','k'))
   iadv = nc.createVariable('adv_xy',fprec,('t',))
   iadv_sp = nc.createVariable('adv_sp',fprec,('t',))
   
   if lbtbc:
      ## KE in bt/bc modes
      ipsd_bke = nc.createVariable('psd_bt_ke',fprec,('t','k'))
      ipsd_beke = nc.createVariable('psd_bt_eke',fprec,('t','k'))
      ipsd_cke = nc.createVariable('psd_bc_ke',fprec,('t','k'))
      ipsd_ceke = nc.createVariable('psd_bc_eke',fprec,('t','k'))
      
      ## Spectral transfer for barotropic KE
      iTk_bt_bc_bc = nc.createVariable('Tk_bt_bc_bc',fprec,('t','k'))
      iTk_bt_bt_bt = nc.createVariable('Tk_bt_bt_bt',fprec,('t','k'))
      iTk_bt_bc_bt = nc.createVariable('Tk_bt_bc_bt',fprec,('t','k'))
      iTk_bt_bt_bc = nc.createVariable('Tk_bt_bt_bc',fprec,('t','k'))
      
      ## Spectral transfer for baroclinic KE
      iTk_bc_bc_bc = nc.createVariable('Tk_bc_bc_bc',fprec,('t','k'))
      iTk_bc_bc_bt = nc.createVariable('Tk_bc_bc_bt',fprec,('t','k'))
      iTk_bc_bt_bt = nc.createVariable('Tk_bc_bt_bt',fprec,('t','k'))
      iTk_bc_bt_bc = nc.createVariable('Tk_bc_bt_bc',fprec,('t','k'))
      
      ## Spectral flux for barotropic KE
      iPi_bt_bc_bc = nc.createVariable('Pi_bt_bc_bc',fprec,('t','k'))
      iPi_bt_bt_bt = nc.createVariable('Pi_bt_bt_bt',fprec,('t','k'))
      iPi_bt_bc_bt = nc.createVariable('Pi_bt_bc_bt',fprec,('t','k'))
      iPi_bt_bt_bc = nc.createVariable('Pi_bt_bt_bc',fprec,('t','k'))
      
      ## Spectral flux for baroclinic KE
      iPi_bc_bc_bc = nc.createVariable('Pi_bc_bc_bc',fprec,('t','k'))
      iPi_bc_bc_bt = nc.createVariable('Pi_bc_bc_bt',fprec,('t','k'))
      iPi_bc_bt_bt = nc.createVariable('Pi_bc_bt_bt',fprec,('t','k'))
      iPi_bc_bt_bc = nc.createVariable('Pi_bc_bt_bc',fprec,('t','k'))
    
   
   ## Planetary vorticity term (Coriolis)
   iTk_cor = nc.createVariable('Tk_cor',fprec,('t','k'))
   iPi_cor = nc.createVariable('Pi_cor',fprec,('t','k'))
   
   ## SSH and surface pressure gradient term
   if lssh:
      iTk_eta = nc.createVariable('Tk_eta',fprec,('t','k'))
      iPi_eta = nc.createVariable('Pi_eta',fprec,('t','k'))
      iTk_spg = nc.createVariable('Tk_spg',fprec,('t','k'))
      iPi_spg = nc.createVariable('Pi_spg',fprec,('t','k'))
   
   ## Hydrostatic pressure gradient term
   if lrho:
      iTk_hpg = nc.createVariable('Tk_hpg',fprec,('t','k'))
      iPi_hpg = nc.createVariable('Pi_hpg',fprec,('t','k'))
   
   ## Total pressure term
   if lssh and lrho:
      iTk_pre = nc.createVariable('Tk_pre',fprec,('t','k'))
      iPi_pre = nc.createVariable('Pi_pre',fprec,('t','k'))
      ipre = nc.createVariable('pre_xy',fprec,('t',))
      ipre_sp = nc.createVariable('pre_sp',fprec,('t',))
   
   ## Horizontal viscosity term
   if lvsc:
      iTk_ke_visc = nc.createVariable('Tk_ke_visc',fprec,('t','k'))
      iPi_ke_visc = nc.createVariable('Pi_ke_visc',fprec,('t','k'))
      ivisc = nc.createVariable('visc_xy',fprec,('t',))
      ivisc_sp = nc.createVariable('visc_sp',fprec,('t',))
      
      if lbtbc:
         iTk_bt_visc = nc.createVariable('Tk_bt_visc',fprec,('t','k'))
         iPi_bt_visc = nc.createVariable('Pi_bt_visc',fprec,('t','k'))
         iTk_bc_visc = nc.createVariable('Tk_bc_visc',fprec,('t','k'))
         iPi_bc_visc = nc.createVariable('Pi_bc_visc',fprec,('t','k'))
         
   ## Wind stress
   if ltau:
      iTk_ke_tau = nc.createVariable('Tk_ke_tau',fprec,('t','k'))
      iPi_ke_tau = nc.createVariable('Pi_ke_tau',fprec,('t','k'))
      itau = nc.createVariable('tau_xy',fprec,('t',))
      itau_sp = nc.createVariable('tau_sp',fprec,('t',))
         
   ## Bottom friction
   if lbfr:
      iTk_bfr = nc.createVariable('Tk_bfr',fprec,('t','k'))
      iPi_bfr = nc.createVariable('Pi_bfr',fprec,('t','k'))
      ibfr = nc.createVariable('bfr_xy',fprec,('t',))
      ibfr_sp = nc.createVariable('bfr_sp',fprec,('t',))
   
   ## Vertical energy fluxes
   if lvertical:
      iTk_ke_vert = nc.createVariable('Tk_ke_vert',fprec,('t','k'))
      iPi_ke_vert = nc.createVariable('Pi_ke_vert',fprec,('t','k'))
      iTk_ke_vvisc = nc.createVariable('Tk_ke_vvisc',fprec,('t','k'))
      iPi_ke_vvisc = nc.createVariable('Pi_ke_vvisc',fprec,('t','k'))
   
   ## Available potential energy terms
   if lape:
      ipsd_ape = nc.createVariable('psd_ape',fprec,('t','k'))
      iTk_ape = nc.createVariable('Tk_ape',fprec,('t','k'))
      iPi_ape = nc.createVariable('Pi_ape',fprec,('t','k'))
      iTk_wb = nc.createVariable('Tk_wb',fprec,('t','k'))
      iPi_wb = nc.createVariable('Pi_wb',fprec,('t','k'))
      iTk_ape_visc = nc.createVariable('Tk_ape_visc',fprec,('t','k'))
      iPi_ape_visc = nc.createVariable('Pi_ape_visc',fprec,('t','k'))
      if lbtbc:
         iTk_bc_ape = nc.createVariable('Tk_bc_ape',fprec,('t','k'))
         iTk_bt_ape = nc.createVariable('Tk_bt_ape',fprec,('t','k'))
         iPi_bc_ape = nc.createVariable('Pi_bc_ape',fprec,('t','k'))
         iPi_bt_ape = nc.createVariable('Pi_bt_ape',fprec,('t','k'))

      
   ## Total domain sums
   ## in wavenumber or physical space
   ip2k = nc.createVariable('pe2ke_xy',fprec,('t',))
   ip2k_sp = nc.createVariable('pe2ke_sp',fprec,('t',))
   
   ## From online tendencies
   if ltend:
      iTk_ke_tend = nc.createVariable('Tk_ke_tend',fprec,('t','k'))
      iPi_ke_tend = nc.createVariable('Pi_ke_tend',fprec,('t','k'))
      iTk_ke_visc_tend = nc.createVariable('Tk_ke_visc_tend',fprec,('t','k'))
      iPi_ke_visc_tend = nc.createVariable('Pi_ke_visc_tend',fprec,('t','k'))
      iTk_spg_tend = nc.createVariable('Tk_spg_tend',fprec,('t','k'))
      iPi_spg_tend = nc.createVariable('Pi_spg_tend',fprec,('t','k'))   
      iTk_hpg_tend = nc.createVariable('Tk_hpg_tend',fprec,('t','k'))
      iPi_hpg_tend = nc.createVariable('Pi_hpg_tend',fprec,('t','k'))
      iTk_pre_tend = nc.createVariable('Tk_pre_tend',fprec,('t','k'))
      iPi_pre_tend = nc.createVariable('Pi_pre_tend',fprec,('t','k'))
      iTk_cor_tend = nc.createVariable('Tk_cor_tend',fprec,('t','k'))
      iPi_cor_tend = nc.createVariable('Pi_cor_tend',fprec,('t','k'))
       
      iadv_tend = nc.createVariable('adv_tend_xy',fprec,('t',))
      iadv_tend_sp = nc.createVariable('adv_tend_sp',fprec,('t',))
      ivisc_tend = nc.createVariable('visc_tend_xy',fprec,('t',))
      ivisc_tend_sp = nc.createVariable('visc_tend_sp',fprec,('t',))
      ipre_tend = nc.createVariable('pre_tend_xy',fprec,('t',))
      ipre_tend_sp = nc.createVariable('pre_tend_sp',fprec,('t',))

      
   ## Enstrophy budget terms   
   ipsd_ens = nc.createVariable('psd_ens',fprec,('t','k'))
   iTk_ens = nc.createVariable('Tk_ens',fprec,('t','k'))
   iPi_ens = nc.createVariable('Pi_ens',fprec,('t','k'))
   
   ## Vertical velocity shear
   ishear   = nc.createVariable('u_shear_rms'  ,fprec,('t',))
   ## Buoyancy anomalies at surface
   if lpe:
      irhoprim = nc.createVariable('rho_prime_rms',fprec,('t',))
   
   ## Return handle for netCDF file
   return nc
   

## ==========================================================================

## 
## Settings for the analysis
##

full_prefix = 'test'  ## unique name for the analysis
grid_type = 'll' ## type of grid (ll - lonlat)
machine   = 'mac'  ## where are you doing the calculations
starttime = datetime(2009,1,5)
endtime   = datetime(2009,12,31)
outputStep = timedelta(days=5)

lcalculate = False  # do calculations and save to files

ltukey  = False  ## use a Tukey window to force zero on boundaries
lrhines = False  ## calculate and plot Rhines scale
lrossby = True    ## read and plot 1st baroclinic Rossby radius
lpsd_freq = False  ## store each step to make frequency spectrum analysis
lbtbc   = False    ## barotropic and baroclinic components
leddy   = False    ## remove mean to get eddy components
lpe     = False  ## calculate potential energy
ltend   = False   ## read full model tendencies
lreadw  = True   ## read vertical velocity (otherwise infer from cont.)
lssh    = False ## calculate ssh and surface pressur gradient
lrho    = False ## calculate buoyancy and hydrostatic pressure gradient
ltau    = True  ## calculate wind stress
lbfr    = True  ## calculate bottom friction
lvsc    = True  ## calculate horizontal viscosity
lvertical = False ## calculate vertical KE flux and viscosity (experimental!)
lape    = False ## APE calculations that dont work!

lsmooth   = False  # smooth data in wavenumber space
nsmooth = 5
lrunningmean = False  # do a running mean in time

linterpolate = False ## interpolate NEMO data to regular grid. May not be necessary

diag_level = 1 ## = 0 no diagnostics plots
               ## = 1 plot KE spectra and map of sfc vorticity each step

ltrend    = False # calculate time trends (not working)
pcrit = 0.1 ## p-limit for significant trend

## if we want potential energy, we need sea surface height and density (if 3D model)
if lpe: 
   lssh=True
   lrho=True

noff = 0

if (machine == 'mac'):
   pdir = '/Users/jkjellsson/Downloads/'
   outdir = '/Users/jkjellsson/Downloads/'
   
elif (machine == 'archer'):
   pdir = '/work/n01/n01/joakim2/frames/'
   outdir = '/work/n01/n01/joakim2/data/psd/'
   
elif (machine == 'jasmin'):
   pdir = './'
   outdir = '/group_workspaces/jasmin2/aopp/joakim/psd/slask/'

pdir = pdir + '/' + full_prefix + '/'
os.system('mkdir -p '+pdir)

## Vertical levels to read in
k0 = 0 ## surface in NEMO
k1 = 75 ## bottom in Andrews ORCA runs

lfull_levels = True ## calculate baroclinic modes from all levels (not working)
ltwo_levels  = False ## average to 2 levels to simplify barotropic/baroclinic calculations

hmin   = 0.   # shallowest level to use 
hmax   = 500. # deepest point to use
hsep   = 500. # depth to split into two levels 

## Diagnostics to plot
lpaper          = False  ## make plots for a paper
lke             = True  ## kinetic energy of top level
lens            = False
lke_freq        = False  ## frequency spectrum of kinetic energy of top level
lbke            = True  ## barotropic kinetic energy
lcke            = True ## baroclinic kinetic energy
lTk_bc_bt       = True  ## spectral transfers of barotropic and baroclinic modes
lTk_wind_visc   = False ## spectral wind forcing and dissipation
lTk_wind_visc_variance = False
lke_budget      = True
lpe_budget      = False
lenergy_budget  = True
lbtbc_budget    = True
lplot_mean_state= False ## plot mean state
lshear          = True
lbtbc_ke_xy     = False
lke_contain     = True ## draw energy containing scale
lbtbc_diff      = False ## difference between resoltuions
plot_surf_ke    = False
##
##
##

## Use viridis colormap if available 
## (only in newer versions of matplotlib)
if 'viridis' in plt.cm.datad.keys():
   cmap = plt.cm.viridis
else:
   cmap = plt.cm.YlGnBu_r

filenameU = []
filenameV = []

name_list = ['INALT10.L46-KJH0017-NEST1']
ddir_list = ['/Users/jkjellsson/data/INALT10.L46/']

regions = ['agulhas-retro'] 

nd = len(name_list) ## number of datasets, not including -5/3 lines etc.
nr = len(regions) ## number of regions to analyse

print ' Regions    : ',regions
print ' Start time : ',starttime
print ' End time   : ',endtime
print ' Data sets  : ',name_list

spectre_list = [] ## list with all analysed data

for jd in range(0,nd):
   
   ## Read global grid
   grid_global = read_data(starttime,name_list[jd],ddir_list[jd],mode='grid')
   
   k0,tmp   = find_kindex(hmin,grid_global['dept_1d'],grid_global['depw_1d'])
   k1,tmp   = find_kindex(hmax,grid_global['dept_1d'],grid_global['depw_1d'])
   ksep,tmp = find_kindex(hsep,grid_global['dept_1d'],grid_global['depw_1d'])
   klevels_keep = np.arange(k0,k1)
   
   region_list = [] ## list with data from all regions for this dataset
   
   for jr in range(0,nr):
      
      ## Size of global grid
      imt = grid_global['tlon'].shape[1]
      jmt = grid_global['tlon'].shape[0]
      
      ## Set lon,lat for current region
      if (regions[jr] != 'global'):
         region_info = set_region(regions[jr],grid_global['tlon'],grid_global['tlat'])
         i0 = region_info['i0']
         i1 = region_info['i1']
         j0 = region_info['j0']
         j1 = region_info['j1']
         lon0 = region_info['lon0']
         lon1 = region_info['lon1']
         lat0 = region_info['lat0']
         lat1 = region_info['lat1']
         print ' lons   : ',region_info['lon0'],region_info['lon1']
         print ' lats   : ',region_info['lat0'],region_info['lat1']
         print ' i0,i1  : ',region_info['i0'],region_info['i1']
         print ' j0,j1  : ',region_info['j0'],region_info['j1']
         print ' Region : ',region_info['region']
         
      else:
         i0 = 0
         i1 = imt
         j0 = 0
         j1 = jmt
      
      ## 
      ## Read the local grid
      ## 
      grid = read_data(starttime,name_list[jd],ddir_list[jd],i0=i0,i1=i1,j0=j0,j1=j1,mode='grid')
      
      ulon = grid['ulon'][:,:]
      ulat = grid['ulat'][:,:]
      vlon = grid['vlon'][:,:]
      vlat = grid['vlat'][:,:]
      tlon = grid['tlon'][:,:]
      tlat = grid['tlat'][:,:]
            
      dxu = grid['dxu'][:,:]
      dxv = grid['dxv'][:,:]
      dxt = grid['dxt'][:,:]
      dxf = grid['dxf'][:,:]
      dyu = grid['dyu'][:,:]
      dyv = grid['dyv'][:,:]
      dyt = grid['dyt'][:,:]
      dyf = grid['dyf'][:,:]
      dzt_full = grid['dzt'][:,:,:]
      dept_full = grid['dept'][:,:,:]
      tmask_full = grid['tmask'][:,:,:]
               
      currentTime = starttime
      
      while currentTime <= endtime:
         
         print ' ====== current time ====='
         print currentTime
         
         ## Read data
         if lrho or lssh:
            readT = True
         else:
            readT = False
         
         if lvertical and not ldiagW:
            readW = True
         else:
            readW = False
         
         ## If you are diagnosing vertical velocity, you must read in all 
         ## vertical levels         
         data = read_data(currentTime,name_list[jd],ddir_list[jd],\
                          i0=i0,i1=i1,j0=j0,j1=j1,\
                          readT=readT, readW=readW, ltend=False)
         
         ## Number of time steps in file
         nt = data['nt']
         
         ## If it is the first time step
         ## we read the global velocity field
         ## and plot all the regions
         if currentTime == starttime and jr == 0:
            data_global = read_data(currentTime,name_list[jd],ddir_list[jd],levels=np.arange(0,2))
            plot_all_regions(grid_global,data_global,[region_info])
         
         
         for jn in range(0,nt):
            
            if currentTime > endtime:
               break
            
            ##
            ##
            uvel = data['uvel'][jn,:,:,:]
            vvel = data['vvel'][jn,:,:,:]
            if lssh:
               etan = data['ssh'][jn,:,:]
            if lrho:
               tem = data['tem'][jn,:,:,:]
               sal = data['sal'][jn,:,:,:]
            
            if lvertical and ldiagW:
               wvel = calculate_w(uvel,vvel,dxv,dyu,\
                      grid['dzt'][:,:,:],grid['dzt'][:,:,:],grid['dzt'][:,:,:])
               
            if ltau:
               ## wind stress 
               taux = data['taux'][jn,:,:]
               tauy = data['tauy'][jn,:,:]
            
            if lpe:
               ## Calculate rho if not AVISO or barotropic MITgcm
               ## NOTE: inside NEMO if statement, can remove if statement below
               rho = np.zeros(tem.shape)
               rhd = np.zeros(tem.shape)
               rhd[:,:,:] = eos_insitu(tem[:,:,:],sal[:,:,:],grid['dept'][:,:,:],grid['tmask'][:,:,:])
               rho[:,:,:] = (rhd[:,:,:]+1)*1026.
               
               fig = plt.figure()
               ax1 = fig.add_subplot(111)
               cf1 = ax1.contourf(tlon,tlat,rho[0,:,:])
               plt.colorbar(cf1,ax=ax1)
            
            ## Load online tendencies from the model 
            if ltend:
               utend_adv = data['utend_adv'][jn,:,:,:] 
               utend_vsc = data['utend_vsc'][jn,:,:,:]
               utend_zdf = data['utend_zdf'][jn,:,:,:]
               utend_cor = data['utend_cor'][jn,:,:,:]
               utend_pre = data['utend_pre'][jn,:,:,:]
               utend_spg = data['utend_spg'][jn,:,:,:]
               utend_hpg = data['utend_hpg'][jn,:,:,:]
               
               vtend_adv = data['vtend_adv'][jn,:,:,:] 
               vtend_vsc = data['vtend_vsc'][jn,:,:,:]
               vtend_zdf = data['vtend_zdf'][jn,:,:,:] 
               vtend_cor = data['vtend_cor'][jn,:,:,:] 
               vtend_pre = data['vtend_pre'][jn,:,:,:] 
               vtend_spg = data['vtend_spg'][jn,:,:,:] 
               vtend_hpg = data['vtend_hpg'][jn,:,:,:] 
               
               if (leddy):
                  print ' === Havent written code to use eddy trends. Stopping === '
                  sys.exit()
               
            if leddy:
               ## Remove time mean to get eddy field
               for nn in range(0,uvel.shape[0]):
                  uvel[:,:,:] = uvel[:,:,:] - u_mean[:,:,:]
                  vvel[:,:,:] = vvel[:,:,:] - v_mean[:,:,:]
                  if lssh:
                     etan[:,:] = etan_full[:,:] - h_mean[:,:]
                  if ltau:
                     taux[:,:] = taux_full[:,:] - taux_mean[:,:]
                     tauy[:,:] = tauy_full[:,:] - tauy_mean[:,:]
                  if lvertical:
                     wvel[:,:,:] = wvel[:,:,:] - w_mean_full[:,:,:]
            
            ## Plot the time averaged fields
            #if (leddy and lplot_mean_state and name[jd] != 'AVISO' and name[jd][0:10] != 'MITgcm-Rob'):
            #   plot_mean_fields()
            
            ## Calculate first baroclinic Rossby radius 
            #if lrossby:
            #   Lross2D = calculate_rossby_radius()
            #   Lross = Lross2D.mean()
            
            ## if region crosses date line, then
            ## dlon will reach big numbers
            ## so we add or remove 360 degree when we
            ## find those points
            dlon = tlon[:,1:]-tlon[:,0:-1]
            dlon = np.where(dlon >  180, dlon-360, dlon)
            dlon = np.where(dlon < -180, dlon+360, dlon)
            dlon = dlon.mean()
            dlat = tlat[1:,:]-tlat[0:-1,:]
            dlat = dlat.mean()
            
            ## Here we set the wavenumber arrays
            ## The code has only been used for grid='ll', i.e. lon-lat grid
            ## but it should be possible and in fact more accurate to use grid='xy'
            xx,yy,wn_x,wn_y,kx,ky,wvsq,k,dk = set_wavenumbers(tlon.min(),tlon.max(),tlat.min(),tlat.max(),\
                                                         dlon,dlat,tlon.shape[1],tlon.shape[0],grid='ll')
            
            ## Make sure the fields are free from masks
            uvel[uvel.mask] = 0.
            vvel[vvel.mask] = 0.
            if lssh:
               etan[etan.mask] = 0.
            if ltau:
               taux[taux.mask] = 0.
               tauy[tauy.mask] = 0. 
            if lvertical:
               wvel[wvel.mask] = 0.
            if lrho:
               tem[tem.mask] = 0.
               sal[sal.mask] = 0.
            
            ##
            ## === Interpolate all data to T points if needed and calculate some gradients ===
            ##

            ## calculate first order, second order, and fourth order gradients of velocity
            grad_data = calculate_uv_gradients_xy(uvel,vvel,dxu,dyu,\
                                                            dxv,dyv,\
                                                            dxt,dyt,\
                                                            dxf,dyf)
            ## interpolate u,v and online tendencies to T points
            uvel[:,1:-1,1:-1] = 0.5 * (uvel[:,1:-1,1:-1] + uvel[:,1:-1,0:-2])
            vvel[:,1:-1,1:-1] = 0.5 * (vvel[:,1:-1,1:-1] + vvel[:,0:-2,1:-1])
            if ltend:
               utend_adv[:,1:-1,1:-1] = 0.5 * (utend_adv[:,1:-1,1:-1] + utend_adv[:,1:-1,0:-2])
               utend_cor[:,1:-1,1:-1] = 0.5 * (utend_cor[:,1:-1,1:-1] + utend_cor[:,1:-1,0:-2])
               utend_pre[:,1:-1,1:-1] = 0.5 * (utend_pre[:,1:-1,1:-1] + utend_pre[:,1:-1,0:-2])
               utend_spg[:,1:-1,1:-1] = 0.5 * (utend_spg[:,1:-1,1:-1] + utend_spg[:,1:-1,0:-2])
               utend_hpg[:,1:-1,1:-1] = 0.5 * (utend_hpg[:,1:-1,1:-1] + utend_hpg[:,1:-1,0:-2])
               utend_vsc[:,1:-1,1:-1] = 0.5 * (utend_vsc[:,1:-1,1:-1] + utend_vsc[:,1:-1,0:-2])
               utend_zdf[:,1:-1,1:-1] = 0.5 * (utend_zdf[:,1:-1,1:-1] + utend_zdf[:,1:-1,0:-2])

               vtend_adv[:,1:-1,1:-1] = 0.5 * (vtend_adv[:,1:-1,1:-1] + vtend_adv[:,0:-2,1:-1])
               vtend_cor[:,1:-1,1:-1] = 0.5 * (vtend_cor[:,1:-1,1:-1] + vtend_cor[:,0:-2,1:-1])
               vtend_pre[:,1:-1,1:-1] = 0.5 * (vtend_pre[:,1:-1,1:-1] + vtend_pre[:,0:-2,1:-1])
               vtend_spg[:,1:-1,1:-1] = 0.5 * (vtend_spg[:,1:-1,1:-1] + vtend_spg[:,0:-2,1:-1])
               vtend_hpg[:,1:-1,1:-1] = 0.5 * (vtend_hpg[:,1:-1,1:-1] + vtend_hpg[:,0:-2,1:-1])
               vtend_vsc[:,1:-1,1:-1] = 0.5 * (vtend_vsc[:,1:-1,1:-1] + vtend_vsc[:,0:-2,1:-1])
               vtend_zdf[:,1:-1,1:-1] = 0.5 * (vtend_zdf[:,1:-1,1:-1] + vtend_zdf[:,0:-2,1:-1])

            ## total depth we are using (used for the wind forcing calculations)
            Hdep = dzt_full[klevels_keep,:,:].sum(axis=0)
            
            ## calculate advection
            u_adv = -uvel * grad_data['dudx_xy'] - vvel * grad_data['dudy_xy']
            v_adv = -uvel * grad_data['dvdx_xy'] - vvel * grad_data['dvdy_xy']
            u_adv = u_adv[klevels_keep,:,:]
            v_adv = v_adv[klevels_keep,:,:]

            ## calculate Coriolis terms
            ## if we are on lon,lat grid => easy
            ## if we are on x-y, I just assume f=1e-4 constant
            omega = 7.2921159e-5
            if grid_type == 'll': ft = 2 * omega * np.sin(np.pi/180. * tlat)
            else: ft = 1e-4
            u_cor =  ft * vvel[:,:,:]
            v_cor = -ft * uvel[:,:,:]
            u_cor = u_cor[klevels_keep,:,:]
            v_cor = v_cor[klevels_keep,:,:]

            if ltau:
               ## calculate wind forcing
               rho0 = 1026.
               u_tau = taux[:,:] / (rho0 * Hdep)
               v_tau = tauy[:,:] / (rho0 * Hdep)
            
            ## calculate bottom friction
            ## default: quadratic
            ## find bottom velocities first
            kmt = np.array(grid['kmt'],dtype='int32')-1
            utmp = np.array(uvel,dtype='float32')
            vtmp = np.array(vvel,dtype='float32')
            dztmp = np.array(dzt_full,dtype='float32')
            u_bot = find_bottom_velocity(utmp,kmt)
            v_bot = find_bottom_velocity(vtmp,kmt)
            dz_bot = find_bottom_velocity(dztmp,kmt)
            dztmp = np.ones((1,u_bot.shape[0],u_bot.shape[1]))
            bfr_data = calculate_bottom_friction_xy(u_bot[np.newaxis,:,:],v_bot[np.newaxis,:,:],dztmp,0,Cd=-1e-3,bg_tke=2.5e-3)
            u_bfr = bfr_data['bfru'] 
            v_bfr = bfr_data['bfrv'] 

            ## calculate hor. viscosity
            if data['visc_form'] == 'laplacian':
               u_vsc = data['Ahm0'] * grad_data['lapu_xy']
               v_vsc = data['Ahm0'] * grad_data['lapv_xy']
            elif data['visc_form'] == 'bilaplacian':
               u_vsc = data['Ahm0'] * grad_data['blpu_xy']
               v_vsc = data['Ahm0'] * grad_data['blpv_xy']
            
            u_vsc = u_vsc[klevels_keep,:,:]
            v_vsc = v_vsc[klevels_keep,:,:]
            
            if lssh and lrho:
               ## calculate pressure gradient
               if leddy:
                  pre_data = calculate_pressure_gradient_xy(etan,rhd-rd_mean_full,\
                                                            dxt,dyt,dzt_full,grav=9.81)
               else:
                  pre_data = calculate_pressure_gradient_xy(etan,rhd,dxt,dyt,dzt_full,grav=9.81)
               
               u_spg = pre_data['spgu'] ## u trend for surface pressure gradient
               v_spg = pre_data['spgv'] ## v trend for surface pressure gradient
               u_hpg = pre_data['hpgu'] ## u trend for hydrostatic pressure gradient
               v_hpg = pre_data['hpgv'] ## v trend for hydrostatic pressure gradient
               u_pre = pre_data['preu'] ## u trend for total pressure gradient
               v_pre = pre_data['prev'] ## v trend for total pressure gradient
   
               u_spg = u_spg[klevels_keep,:,:]
               v_spg = v_spg[klevels_keep,:,:]
               u_hpg = u_hpg[klevels_keep,:,:]
               v_hpg = v_hpg[klevels_keep,:,:]
               u_pre = u_pre[klevels_keep,:,:]
               v_pre = v_pre[klevels_keep,:,:]

            uvel     = uvel[klevels_keep,:,:]
            vvel     = vvel[klevels_keep,:,:]
            if lvertical:
               wvel     = wvel[klevels_keep,:,:]
            if lpe:
               tem      = tem[klevels_keep,:,:]
               sal      = sal[klevels_keep,:,:]
               rho      = rho[klevels_keep,:,:]
               rhd      = rhd[klevels_keep,:,:]
            if leddy:
               t_mean   = t_mean_full[klevels_keep,:,:]
               s_mean   = s_mean_full[klevels_keep,:,:]
               r_mean   = r_mean_full[klevels_keep,:,:]
               rd_mean  = rd_mean_full[klevels_keep,:,:]
               u_mean   = u_mean_full[klevels_keep,:,:]
               v_mean   = v_mean_full[klevels_keep,:,:]
               w_mean   = w_mean_full[klevels_keep,:,:]
            
            ##
            ##  If we are using data with several vertical levels
            ##  i.e. not a barotropic model
            ##
            if grid['nz'] > 1:
               
               ## Temporary depth variable 
               ## only depth of the layers we are studying
               deptht_tmp = dept_full[klevels_keep,:,:].copy()
               
               ## same for level thickness
               dzt = grid['dzt'][klevels_keep,:,:].copy()
               dz  = grid['dept_1d'][klevels_keep].copy()
               
               ## Save square of density anomalies at west most edge
               if lpe:
                  tmp1 = (rho[np.newaxis,:,:,0] - r_mean[np.newaxis,:,:,0])**2
                  if currentTime == starttime:
                     cross_rho = tmp1
                  else:
                     cross_rho = np.concatenate( (cross_rho,tmp1),axis=0 )
                     
               tmp2 = uvel[np.newaxis,:,:,0]
               if currentTime == starttime:
                  cross_uvel = tmp2
               else:
                  cross_uvel = np.concatenate( (cross_uvel,tmp2), axis=0 )
               
               if lpe:
                  ## Calculate buoyancy as negative density anomalies,
                  ## i.e. anomalouly light fluid is positive buoyancy anomaly
                  rho0 = 1023.
                  bprim = -9.81/rho0 * (rho - r_mean) 

               if lfull_levels:
                  if lvertical:
                     ## interpolate vertical velocity to T point
                     wtmp = np.ma.zeros(wvel.shape)
                     wtmp[0:-1,:,:] = 0.5 * (wvel[0:-1,:,:]+wvel[1:,:,:])
                     wtmp[-1,:,:] = np.ma.masked
                     wvel = wtmp.copy()
                     
                     ##
                     ## Calculate vertical velocity gradients on full levels
                     ## (I'm not too sure that this is the best way...)
                     ##
                     dudz = np.ma.zeros(uvel.shape)
                     dvdz = np.ma.zeros(vvel.shape)
                     dudz[0 ,:,:] = np.ma.masked_where( dudz[0 ,:,:] == 0, dudz[0 ,:,:] ) ## mask
                     dvdz[0 ,:,:] = np.ma.masked_where( dvdz[0 ,:,:] == 0, dvdz[0 ,:,:] ) ## top
                     dudz[-1,:,:] = np.ma.masked_where( dudz[-1,:,:] == 0, dudz[-1,:,:] ) ## and bottom
                     dvdz[-1,:,:] = np.ma.masked_where( dvdz[-1,:,:] == 0, dvdz[-1,:,:] ) ## levels
                     dudz[1:-1,:,:] = (uvel[0:-2,:,:] - uvel[2:,:,:])/(deptht_tmp[0:-2,:,:]-deptht_tmp[2:,:,:])
                     dvdz[1:-1,:,:] = (vvel[0:-2,:,:] - vvel[2:,:,:])/(deptht_tmp[0:-2,:,:]-deptht_tmp[2:,:,:])

                  ## RMS of shear (diff between top and bottom level) and 
                  ## density (top level buoyancy) anomalies
                  ushear  = np.sqrt( np.mean( (uvel[0,:,:] - uvel[-1,:,:])**2 + (vvel[0,:,:] - vvel[-1,:,:])**2 ) )
                  if lpe:
                     rhoprim = np.sqrt( np.mean( (rho[0,:,:] - r_mean[0,:,:])**2 ) )
                  
                  ## Calculate APE
                  ## Exactly how to calculate APE is debatable
                  ## Here I adopt the rather simple version of b2/N2
                  ## as in Capet et al. 2008.
                  if lpe:
                     N1 = np.zeros(uvel.shape)
                     N2 = np.zeros(uvel.shape)
                     ape = np.zeros((uvel.shape[1],uvel.shape[2]))
                     Hdep = np.sum(dzt,axis=0)
                     ## calculate N2
                     ## WE SHOULD REALLY DO THIS USING THE TEOS-10 ALGORITHMS!
                     for jk in range(1,uvel.shape[0]-1):
                        N2[jk,:,:] = 9.81/rho0 * (r_mean[jk-1,:,:] - r_mean[jk+1,:,:]) / (deptht_tmp[jk-1,:,:] - deptht_tmp[jk+1,:,:])
                     
                     N2_mean = N2.mean(axis=0) * np.ones(N2.shape)
                     ## ensure N2 not negative
                     ## if it is negative, replace with time mean value
                     N1 = np.ma.where( N2 <= 0., np.sqrt(N2_mean), N1 )
                     N2 = N2_mean
                     N1 = np.sqrt(N2_mean)
                     ## vertical integral 
                     for jk in range(0,uvel.shape[0]-1):
                        ape  += 0.5 * bprim[jk,:,:]**2 / N2[jk,:,:] * dzt[jk,:,:] #/ Hdep
                        print ' === N1 min,max, jk ',N1[jk,:,:].min(),N1[jk,:,:].max(),jk
                        
                     #fig = plt.figure()
                     #ax1 = fig.add_subplot(111)
                     #cf1 = ax1.contourf(np.sqrt(N2_mean[0,:,:]))
                     #plt.colorbar(cf1)

                     #plt.show()
                     #sys.exit()


               if (ltwo_levels):
                  ##
                  ## Convert to two-layer model
                  ##
                  
                  ## the upper and lower levels
                  k_upper = np.arange(0,ksep)
                  k_lower = np.arange(ksep,grid['dept_1d'][klevels_keep].shape[0])
                  
                  z_two,h1,h2 = two_layer(deptht_tmp,dzt,k_interface)
                  z1 = z_two[0,:,:].mean()
                  z2 = z_two[1,:,:].mean()

                  #print ' Average all levels to two levels at z1, z2 ',z1,z2
                  #print ' mean depth of interface = ',deptht_tmp[k_interface,:,:].mean()

                  u_two = np.zeros(uvel.shape)
                  u_two = u_two[0:2,:,:]
                  v_two = np.zeros(u_two.shape)
                  w_two = np.zeros(u_two.shape)
                  if (1):
                     t_two = np.zeros(u_two.shape)
                     s_two = np.zeros(u_two.shape)
                     r_two = np.zeros(u_two.shape)

                  ## Interpolate u,v to two layers separated at level = k_interface
                  u_two,h1,h2 = two_layer(uvel[:,:,:],dzt,k_interface)
                  v_two,h1,h2 = two_layer(vvel[:,:,:],dzt,k_interface)
                  if (jn == 0 and leddy):
                     u_mean_two,h1,h2 = two_layer(u_mean[:,:,:],dzt,k_interface)
                     v_mean_two,h1,h2 = two_layer(v_mean[:,:,:],dzt,k_interface)
                     if (grid_type == 'xy'):
                        u_mean_xy = u_mean_two.copy()
                        v_mean_xy = v_mean_two.copy()

                  ## Interpolate w to mid-layers and then average for each layer as done for u,v
                  wtmp = np.ma.zeros(wvel.shape)
                  wtmp[0:-1,:,:] = 0.5 * (wvel[0:-1,:,:]+wvel[1:,:,:])
                  wtmp[-1,:,:] = np.ma.masked
                  w_two,h1,h2 = two_layer(wtmp[:,:,:],dzt,k_interface)
                  ## Only keep vertical velocity at interface
                  #w_two[1,:,:] = wvel[k_interface,:,:]
                  if (1):
                     t_two,h1,h2 = two_layer(tem[:,:,:] ,dzt,k_interface)
                     s_two,h1,h2 = two_layer(sal[:,:,:] ,dzt,k_interface)
                     ##
                     ## Averaging rho onto two levels
                     ## We could also calculate rho from two-level averages of T,S
                     ## which would seem less accurate.
                     ## Perhaps test both?
                     ##
                     r_two,h1,h2 = two_layer(rho[:,:,:] ,dzt,k_interface)
                     r_two_mean,h1,h2 = two_layer(r_mean[:,:,:], dzt,k_interface)

                     ## Calculate APE locally (Molemaker & McWilliams, J. Fluid Mech., 2010)
                     N2 = np.zeros(uvel.shape)
                     #bprim = np.zeros(uvel.shape)
                     ape = np.zeros(uvel.shape)
                     rho0 = 1023.
                     N2[1:,:,:] = 9.81/rho0 * (r_mean[0:-1,:,:] - r_mean[1:,:,:]) / (deptht_tmp[0:-1,:,:] - deptht_tmp[1:,:,:])
                     N2 = np.ma.masked_where(N2 <= 0,N2)
                     #bprim[:,:,:] = 9.81/rho0 * (rho[:,:,:] - r_mean[:,:,:])
                     b2 = bprim**2
                     ape[1:,:,:] = b2[1:,:,:] / N2[1:,:,:]
                     print ' APE, b2, N2 min,max ',ape.min(),ape.max(),b2.min(),b2.max(),N2.min(),N2.max()
                     ape_two,h1,h2 = two_layer(ape[:,:,:] ,dzt,k_interface)
                     N2_two,h1,h2 = two_layer(N2[:,:,:] ,dzt,k_interface)

                     if (0):
                        fig = plt.figure()
                        ax1 = fig.add_subplot(221)
                        ax1.set_title('N2')
                        cf1 = ax1.contourf(N2_two[0,:,:])
                        plt.colorbar(cf1,ax=ax1)

                        ax2 = fig.add_subplot(222)
                        ax2.set_title('bprim')
                        cf2 = ax2.contourf(9.81/rho0 * (r_two[0,:,:] - r_two_mean[0,:,:]))
                        plt.colorbar(cf2,ax=ax2)

                        ax3 = fig.add_subplot(223)
                        ax3.set_title('bprim^2/2N2')
                        cf3 = ax3.contourf(0.5*(9.81/rho0 * (r_two[0,:,:] - r_two_mean[0,:,:]))**2/N2_two[0,:,:])
                        plt.colorbar(cf3,ax=ax3)

                        ax4 = fig.add_subplot(224)
                        ax4.set_title('1/2 u^2 + v^2')
                        cf4 = ax4.contourf(0.5*(u_two[0,:,:]**2 + v_two[0,:,:]))
                        plt.colorbar(cf4,ax=ax4)

                  ## Interpolate tendencies onto two levels
                  if (ltend):
                     utend_adv,h1,h2 = two_layer(utend_adv[:,:,:] ,dzt,k_interface)
                     utend_cor,h1,h2 = two_layer(utend_cor[:,:,:] ,dzt,k_interface)
                     utend_spg,h1,h2 = two_layer(utend_spg[:,:,:] ,dzt,k_interface)
                     utend_hpg,h1,h2 = two_layer(utend_hpg[:,:,:] ,dzt,k_interface)
                     utend_pre,h1,h2 = two_layer(utend_pre[:,:,:] ,dzt,k_interface)
                     utend_vsc,h1,h2 = two_layer(utend_vsc[:,:,:] ,dzt,k_interface)
                     utend_zdf,h1,h2 = two_layer(utend_zdf[:,:,:] ,dzt,k_interface)

                     vtend_adv,h1,h2 = two_layer(vtend_adv[:,:,:] ,dzt,k_interface)
                     vtend_cor,h1,h2 = two_layer(vtend_cor[:,:,:] ,dzt,k_interface)
                     vtend_spg,h1,h2 = two_layer(vtend_spg[:,:,:] ,dzt,k_interface)
                     vtend_hpg,h1,h2 = two_layer(vtend_hpg[:,:,:] ,dzt,k_interface)
                     vtend_pre,h1,h2 = two_layer(vtend_pre[:,:,:] ,dzt,k_interface)
                     vtend_vsc,h1,h2 = two_layer(vtend_vsc[:,:,:] ,dzt,k_interface)
                     vtend_zdf,h1,h2 = two_layer(vtend_zdf[:,:,:] ,dzt,k_interface)

                  u_adv,h1,h2 = two_layer(u_adv[:,:,:] ,dzt,k_interface)
                  u_cor,h1,h2 = two_layer(u_cor[:,:,:] ,dzt,k_interface)
                  u_spg,h1,h2 = two_layer(u_spg[:,:,:] ,dzt,k_interface)
                  u_hpg,h1,h2 = two_layer(u_hpg[:,:,:] ,dzt,k_interface)
                  u_pre,h1,h2 = two_layer(u_pre[:,:,:] ,dzt,k_interface)
                  u_vsc,h1,h2 = two_layer(u_vsc[:,:,:] ,dzt,k_interface)

                  v_adv,h1,h2 = two_layer(v_adv[:,:,:] ,dzt,k_interface)
                  v_cor,h1,h2 = two_layer(v_cor[:,:,:] ,dzt,k_interface)
                  v_spg,h1,h2 = two_layer(v_spg[:,:,:] ,dzt,k_interface)
                  v_hpg,h1,h2 = two_layer(v_hpg[:,:,:] ,dzt,k_interface)
                  v_pre,h1,h2 = two_layer(v_pre[:,:,:] ,dzt,k_interface)
                  v_vsc,h1,h2 = two_layer(v_vsc[:,:,:] ,dzt,k_interface)

                  ## Calculate vertical velocity gradients on full levels,
                  ## then average du/dz and dv/dz on two layers to get diffusivity on two layers
                  dudz = np.ma.zeros(uvel.shape)
                  dvdz = np.ma.zeros(vvel.shape)
                  dudz[0 ,:,:] = np.ma.masked_where( dudz[0 ,:,:] == 0, dudz[0 ,:,:] ) ## mask
                  dvdz[0 ,:,:] = np.ma.masked_where( dvdz[0 ,:,:] == 0, dvdz[0 ,:,:] ) ## top
                  dudz[-1,:,:] = np.ma.masked_where( dudz[-1,:,:] == 0, dudz[-1,:,:] ) ## and bottom
                  dvdz[-1,:,:] = np.ma.masked_where( dvdz[-1,:,:] == 0, dvdz[-1,:,:] ) ## levels
                  dudz[1:-1,:,:] = (uvel[0:-2,:,:] - uvel[2:,:,:])/(deptht_tmp[0:-2,:,:]-deptht_tmp[2:,:,:])
                  dvdz[1:-1,:,:] = (vvel[0:-2,:,:] - vvel[2:,:,:])/(deptht_tmp[0:-2,:,:]-deptht_tmp[2:,:,:])
                  dudz_two,h1,h2 = two_layer(dudz[:,:,:], dzt,k_interface)
                  dvdz_two,h1,h2 = two_layer(dvdz[:,:,:], dzt,k_interface)

                  ## RMS of shear and density anomalies
                  ushear  = np.sqrt( np.mean( (u_two[0,:,:] - u_two[1,:,:])**2 + (v_two[0,:,:] - v_two[1,:,:])**2 ) )
                  if (lpe):
                     rhoprim = np.sqrt( np.mean( (r_two[0,:,:] - r_two_mean[0,:,:])**2 ) )

                  uvel = u_two.copy()
                  vvel = v_two.copy()
                  wvel = w_two.copy()
                  u_bot = u_two[-1,:,:]
                  v_bot = v_two[-1,:,:]
                  dz   = np.array( [np.sum(dz[:k_interface]), np.sum(dz[k_interface:])] )
                  dzt  = np.concatenate( (h1[np.newaxis,:,:],h2[np.newaxis,:,:]), axis=0 )
                  deptht_tmp = z_two.copy()

                  ## Calculate APE
                  if (lpe):
                     rho0 = 1023.
                     N2 = np.zeros(uvel.shape)
                     bprim = np.zeros(uvel.shape)
                     N1 = np.zeros(uvel.shape)
                     ape = np.zeros((uvel.shape[1],uvel.shape[2]))
                     for jk in range(0,uvel.shape[0]):
                        N2[jk,:,:] = 9.81/rho0 * (r_two_mean[0,:,:] - r_two_mean[1,:,:]) / (z_two[0,:,:] - z_two[1,:,:])
                        bprim[jk,:,:] = 9.81/rho0 * (r_two[jk,:,:] - r_two_mean[jk,:,:])

                     N2_mean = N2.mean(axis=0) * np.ones(N2.shape)
                     #N1 = np.ma.where( N2 <= 0., np.sqrt(N2_mean), N1 )
                     N2 = N2_mean
                     N1 = np.sqrt(N2_mean)
                     for jk in range(0,uvel.shape[0]):
                        #N1[jk,:,:] = np.sqrt(N2[jk,:,:])
                        ape  += 0.5 * bprim[jk,:,:]**2 #/ N1[jk,:,:]**2 * dzt[jk,:,:] / np.sum(dzt,axis=0)
                        print ' === APE min,max, jk ',ape[:,:].min(),ape[:,:].max(),N1[jk,:,:].min(),N1[jk,:,:].max()

                        if (0):
                           fig = plt.figure()
                           ax1 = fig.add_subplot(221)
                           ax2 = fig.add_subplot(222)
                           ax3 = fig.add_subplot(223)
                           ax4 = fig.add_subplot(224)
                           ax1.set_title('rho_prim')
                           ax2.set_title('N1')
                           ax3.set_title('dz')
                           ax4.set_title('w')
                           cf1 = ax1.contourf(bprim[jk,:,:])
                           cf2 = ax2.contourf(N1[jk,:,:])
                           cf3 = ax3.contourf(z_two[0,:,:]-z_two[1,:,:])
                           cf4 = ax4.contourf(wvel[jk,:,:])
                           plt.colorbar(cf1,ax=ax1)
                           plt.colorbar(cf2,ax=ax2)
                           plt.colorbar(cf3,ax=ax3)
                           plt.colorbar(cf4,ax=ax4)

                           fig = plt.figure()
                           ax1 = fig.add_subplot(121)
                           ax2 = fig.add_subplot(122)
                           ke = 0.5*(uvel[jk,:,:]**2 + vvel[jk,:,:]**2)
                           title1 = 'APE, '+regions[jr]+' min=%5f, max=%5f' % (ape.min(),ape.max())
                           title2 = ' KE, '+regions[jr]+' min=%5f, max=%5f' % (ke.min(),ke.max())
                           ax1.set_title(title1,fontsize=10)
                           ax2.set_title(title2,fontsize=10)
                           cf1 = ax1.contourf(ape)
                           cf2 = ax2.contourf(ke)
                           plt.colorbar(cf1,ax=ax1)
                           plt.colorbar(cf2,ax=ax2)



               if (lrhines):
                  Lrhines = calculate_rhines_scale(np.sqrt(np.mean(uvel**2+vvel**2)),tlat.mean())

            ## If grid is lon-lat,
            ## we interpolate to regular lon-lat grid
            ## and then use the wavenumbers in x-y from above
            ## In Kjellsson & Zanna (Fluids, 2017) it was not really necessary
            if linterpolate:
               data_irr = []
               data_irr.append(uvel)
               data_irr.append(vvel)
               
               data_irr.append(dzt)
               data_irr.append(deptht_tmp)
               
               if 1:
                  data_irr.append(u_adv)
                  data_irr.append(v_adv)
                  data_irr.append(u_cor)
                  data_irr.append(v_cor)
               
               if lvertical:
                  data_irr.append(wvel)
                  data_irr.append(dudz)
                  data_irr.append(dvdz)
                  
               if lssh:
                  data_irr.append(u_spg)
                  data_irr.append(v_spg)
                  
               if lrho:
                  data_irr.append(u_hpg)
                  data_irr.append(v_hpg)
               
               if ltau:
                  data_irr.append(u_tau)
                  data_irr.append(v_tau)
               
               if lbfr:
                  data_irr.append(u_bfr)
                  data_irr.append(v_bfr)
               
               if lvsc:
                  data_irr.append(u_vsc)
                  data_irr.append(v_vsc)
               
               if ltend:
                  data_irr.append(utend_adv)
                  data_irr.append(utend_cor)
                  data_irr.append(utend_pre)
                  data_irr.append(utend_vsc)
                  data_irr.append(utend_zdf)
                  
                  data_irr.append(vtend_adv)
                  data_irr.append(vtend_cor)
                  data_irr.append(vtend_pre)
                  data_irr.append(vtend_vsc)
                  data_irr.append(vtend_zdf)
               
               if lape:
                  data_irr.append(bprim)
                  data_irr.append(N1)
                  data_irr.append(ape)
                  data_irr.append(ape2)
               
               data_reg = interpolate_alldata(data_irr,tlon[:,:],tlat[:,:],xx,yy)
            
            ## =====================================================
            ##
            ##  Spectral fluxes and KE in the frequency domain
            ##  c.f. Arbic et al. 2012 (https://doi.org/10.1175/JPO-D-11-0151.1)
            ##
            ##  This requires storing data from ALL time steps 
            ##  which is very memory expensive
            ##  I think the best way to go in this case is to 
            ##  split the full domain into smaller ones 
            ##
            ##  So, I would recommend making lots of small "regions"
            ##  and then aggreating all data using this option
            ##
            ## =====================================================
            
            ## Save uvel, vvel, taux, tauy
            ## This is necessary for frequency spectrum calculations
            ## But saving all time steps in memory is expensive!
            if (lpsd_freq):
               if (jn == 0):
                  uvel_store = uvel[np.newaxis,:,:,:]
                  vvel_store = vvel[np.newaxis,:,:,:]
                  if (name[jd] != 'AVISO'):
                     taux_store = taux[np.newaxis,:,:]
                     tauy_store = tauy[np.newaxis,:,:]
                     if (name[jd] != 'MITgcm-Rob'):
                        if (1):
                           wb_store  = bprim[np.newaxis,0,:,:] * wvel[np.newaxis,0,:,:]
                        if (lpe):
                           ape_store = ape[np.newaxis,:,:]
                           #rhop_store= (rho[np.newaxis,:,:,0] - r_mean[np.newaxis,:,:,0])**2

               else:
                  uvel_store = np.concatenate( (uvel[np.newaxis,:,:,:],uvel_store), axis=0 )
                  vvel_store = np.concatenate( (vvel[np.newaxis,:,:,:],vvel_store), axis=0 )

                  if (name[jd] != 'AVISO'):
                     taux_store = np.concatenate( (taux[np.newaxis,:,:],taux_store), axis=0 )
                     tauy_store = np.concatenate( (tauy[np.newaxis,:,:],tauy_store), axis=0 )
                     if (name[jd] != 'MITgcm-Rob'):
                        if (1):
                           wb_store  = np.concatenate( (bprim[np.newaxis,0,:,:]*wvel[np.newaxis,0,:,:], wb_store), axis=0 )
                        if (lpe):
                           ape_store = np.concatenate( (ape[np.newaxis,:,:],ape_store), axis=0 )
                           #rhop_store= np.concatenate( ((rho[np.newaxis,:,:,0]-r_mean[np.newaxis,:,:,0])**2,rho_store), axis=0 )


            ## Calculate relative vorticity
            ## needed to calculate enstrophy fluxes
            vort = np.zeros((uvel.shape[0],uvel.shape[1],uvel.shape[2]))
            for jk in range(0,uvel.shape[0]):
               vort[jk,:,:] = calculate_vorticity(uvel[jk,:,:],vvel[jk,:,:],xx,yy)
            
            if diag_level >= 1:
               fig = plt.figure()
               ax = fig.add_subplot(1,1,1)
               ax.set_title('uppermost vorticity')
               cf = ax.contourf(vort[0,:,:],levels=np.arange(-0.2,0.22,0.02),cmap=cmap,extend='both')
               plt.colorbar(cf,ax=ax)
               fig.savefig(pdir+'/surf_vorticity.png',format='png')
            
            ## Default: No window function
            window_x = np.ones(tlon.shape[1])
            window_y = np.ones(tlon.shape[0])

            if ltukey:
               ## Use a Tukey window
               ## alpha between 0 and 1 is the fraction of the domain that is tapered
               alpha = 0.2
               window_x = tukey(nx,alpha=alpha)
               window_y = tukey(ny,alpha=alpha)
            
            win = np.zeros(tlon.shape)
            for jj in range(0,tlon.shape[0]):
               win[jj,:] = window_x[:]**2 * window_y[jj]**2
            
            if diag_level >= 1:
               fig = plt.figure()
               ax1 = fig.add_subplot(111)
               ax1.set_title('Min, Max of window function: %3f,%3f' %(win.min(),win.max()))
               cf1 = ax1.contourf(win)
               plt.colorbar(cf1,ax=ax1)
               fig.savefig('window_function.png',format='png')
            
            ##
            ## Barotropic and baroclinic calculations
            ## Only used if we have more than one vertical level
            ## (otherwise everything is barotropic!)
            ##

            if lbtbc and grid['nz'] > 1:
               if lfull_levels:
                  ## calculate vertical mean as barotropic component
                  ## and the anomalies from the vertical mean as baroclinic
                  ## (see Chemke & Kaspi, GRL 2016)
                  means_anomalies = calculate_means_and_anomalies(uvel,vvel,dzt)
                  um = means_anomalies['u_vm'][:,:]
                  vm = means_anomalies['v_vm'][:,:]
                  ut = means_anomalies['u_va'][:,:]
                  vt = means_anomalies['v_va'][:,:]
                  
                  ## Barotropic KE (includes zonal mean)
                  bke_2D = calculate_ke(kx,ky,means_anomalies['u_vm']*win,means_anomalies['v_vm']*win)
                  ## Barotropic EKE (deviations from zonal mean)
                  beke_2D = calculate_ke(kx,ky,means_anomalies['u_za_vm']*win,means_anomalies['v_za_vm']*win)
                  ## Baroclinic KE and EKE
                  cke_2D = np.zeros((beke_2D.shape))
                  ceke_2D = np.zeros((beke_2D.shape))
                  ## sum vertically
                  for jk in range(0,uvel.shape[0]):
                     cke_tmp = calculate_ke(kx,ky,means_anomalies['u_va'][jk,:,:]*win,\
                                                  means_anomalies['v_va'][jk,:,:]*win) / float(uvel.shape[0])
                     ceke_tmp = calculate_ke(kx,ky,means_anomalies['u_za_va'][jk,:,:]*win,\
                                                   means_anomalies['v_za_va'][jk,:,:]*win) / float(uvel.shape[0])
                     cke_2D  = cke_2D  + cke_tmp
                     ceke_2D = ceke_2D + ceke_tmp
                  
                  ## Calculate spectral fluxes of triad interactions
                  ## Change in barotropic EKE from baroclinic interactions etc.
                  ## NOTE: Here we calculate gradients as du/dx = i*k*uhat
                  ## which is not the same way as NEMO does it
                  ## We should really calculate the gradients in physical space first!
                  
                  ## Warning: This function assumes layers are equally thick! 
                  ## Not the case for most ocean models
                  ## In Kjellsson & Zanna (Fluids, 2017) we used 
                  ## calculate_spectral_flux_baroclinic_barotropic5 
                  ## which does not make this assumpion!
                  Tk_data = calculate_spectral_flux_baroclinic_barotropic(kx,ky,\
                            means_anomalies['u_za_vm']*win,means_anomalies['v_za_vm']*win,\
                            means_anomalies['u_za_va']*win,means_anomalies['v_za_va']*win)

                  adv_bt_bc_bc = Tk_data['Tk_bt_bc_bc']
                  adv_bt_bt_bt = Tk_data['Tk_bt_bt_bt']
                  adv_bc_bc_bc = Tk_data['Tk_bc_bc_bc']
                  
                  if lape:
                     # APE advection
                     Pk_data = calculate_spectral_ape_flux(kx,ky,uvel*win,vvel*win,ape)
                     adv_bt_ape = Pk_data['adv_ape']
                     adv_bc_ape = Pk_data['adv_ape']


               elif (ltwo_levels):

                  ## Barotropic velocity
                  um = np.sum(uvel[:,:,:]*dzt[:,:,:],axis=0) / np.sum(dzt,axis=0)
                  vm = np.sum(vvel[:,:,:]*dzt[:,:,:],axis=0) / np.sum(dzt,axis=0)
                  ## Baroclinic velocity as difference from barotropic for each layer
                  ut = np.zeros(uvel.shape)
                  vt = np.zeros(vvel.shape)
                  for jk in range(0,uvel.shape[0]):
                     ut[jk,:,:] = (uvel[jk,:,:] - um[:,:])
                     vt[jk,:,:] = (vvel[jk,:,:] - vm[:,:])

                  ## Plot mean barotropic and baroclinic KE
                  if (jn == 0 and leddy):
                     um_mean = np.sum(u_mean_xy[:,:,:]*dzt[:,:,:],axis=0) / np.sum(dzt,axis=0)
                     vm_mean = np.sum(v_mean_xy[:,:,:]*dzt[:,:,:],axis=0) / np.sum(dzt,axis=0)
                     ut_mean = np.zeros(uvel.shape)
                     vt_mean = np.zeros(vvel.shape)
                     for jk in range(0,uvel.shape[0]):
                        ut_mean[jk,:,:] = (u_mean_xy[jk,:,:] - um_mean[:,:])
                        vt_mean[jk,:,:] = (v_mean_xy[jk,:,:] - vm_mean[:,:])
                     
                     if ldiag_lev >= 1:
                        fig = plt.figure()
                        ax1 = fig.add_subplot(121)
                        ax2 = fig.add_subplot(122)
                        zplot1 = 0.5*(um_mean**2 + vm_mean**2)
                        zplot2 = 0.5*(ut_mean[0,:,:]**2 + vt_mean[0,:,:]**2)
                        cf1 = ax1.contourf( zplot1 )
                        cf2 = ax2.contourf( zplot2 )
                        plt.colorbar(cf1,ax=ax1)
                        plt.colorbar(cf2,ax=ax2)
                        fig.savefig('bt_bc_KE_xy_'+name[jd]+'_'+regions[jr]+'.pdf',format='pdf')

                  if (lbtbc_ke_xy):
                     if (jn == 0):
                        btke_store = 0.5 * (um[np.newaxis,  :,:]**2 + vm[np.newaxis,  :,:]**2)
                        bcke_store = 0.5 * (ut[np.newaxis,0,:,:]**2 + vt[np.newaxis,0,:,:]**2)
                     else:
                        btke_store = np.concatenate( (0.5*(um[np.newaxis,:,:]**2+vm[np.newaxis,:,:]**2),btke_store),axis=0 )
                        bcke_store = np.concatenate( (0.5*(ut[np.newaxis,0,:,:]**2+vt[np.newaxis,0,:,:]**2),bcke_store),axis=0 )


                  if(0):
                     fig = plt.figure()
                     ax1 = fig.add_subplot(221)
                     ax2 = fig.add_subplot(222)
                     ax3 = fig.add_subplot(223)
                     ax4 = fig.add_subplot(224)
                     ax1.set_title('u_bc_1')
                     cf1 = ax1.contourf(ut[0,:,:])
                     plt.colorbar(cf1,ax=ax1)
                     ax2.set_title('u_bc_2')
                     cf2 = ax2.contourf(ut[1,:,:])
                     plt.colorbar(cf2,ax=ax2)
                     ax3.set_title('u_bc_mean')
                     cf3 = ax3.contourf(np.sum(ut[:,:,:]*dzt[:,:,:],axis=0)/np.sum(dzt,axis=0))
                     plt.colorbar(cf3,ax=ax3)
                     ax4.set_title('u_bt')
                     cf4 = ax4.contourf(um[:,:])
                     plt.colorbar(cf4,ax=ax4)
                     fig.savefig('u_bt_bc_'+name[jd]+'.png',format='png')

                     i = np.complex(0,1)
                     uhat_bt = fftn(um)
                     vhat_bt = fftn(vm)
                     # derivatives for barotropic velocities
                     ddx_u_bt = np.real( ifftn(i*kx*uhat_bt) ) # du_bt/dx
                     ddy_u_bt = np.real( ifftn(i*ky*uhat_bt) ) # du_bt/dy
                     ddx_v_bt = np.real( ifftn(i*kx*vhat_bt) ) # dv_bt/dx
                     ddy_v_bt = np.real( ifftn(i*ky*vhat_bt) ) # dv_bt/dy

                     # adv_u = u * du/dx + v * du/dy
                     adv_u_bt_bt = um[:,:] * ddx_u_bt + vm[:,:] * ddy_u_bt
                     # adv_v = u * dv/dx + v * dv/dy
                     adv_v_bt_bt = um[:,:] * ddx_v_bt + vm[:,:] * ddy_v_bt

                     fig,ax = plt.subplots(2,2)
                     Tk_bc_bt_bt1 = np.zeros(uhat_bt.shape)
                     Tk_bc_bt_bt2 = np.zeros(uhat_bt.shape)
                     for kkk in range(0,2):
                        uhat_bc = fftn(ut[kkk,:,:])
                        vhat_bc = fftn(vt[kkk,:,:])
                        dzhat   = fftn(dzt[kkk,:,:])
                        adv_u_bt_bt_new = adv_u_bt_bt * dzt[kkk,:,:]
                        adv_v_bt_bt_new = adv_v_bt_bt * dzt[kkk,:,:]
                        zplot1 = -np.conj(uhat_bc) * fftn(adv_u_bt_bt_new) - np.conj(vhat_bc) * fftn(adv_v_bt_bt_new)
                        zplot2 = -(np.conj(uhat_bc) * fftn(adv_u_bt_bt) - np.conj(vhat_bc) * fftn(adv_v_bt_bt)) * dzhat
                        cf = ax[0,kkk].contourf(np.log10( np.abs(zplot1) ))
                        plt.colorbar(cf,ax=ax[0,kkk])
                        ax[0,kkk].set_title('transfer jk=%2d' % (kkk,))
                        Tk_bc_bt_bt1 = Tk_bc_bt_bt1 + np.real( zplot1 )
                        Tk_bc_bt_bt2 = Tk_bc_bt_bt2 + np.real( zplot2 )

                     cf = ax[1,0].contourf(np.log10( np.abs(Tk_bc_bt_bt1) ))
                     plt.colorbar(cf,ax=ax[1,0])
                     ax[1,0].set_title('transfer total 1')
                     cf = ax[1,1].contourf(np.log10( np.abs(Tk_bc_bt_bt2) ))
                     plt.colorbar(cf,ax=ax[1,1])
                     ax[1,1].set_title('transfer total 2')
                     fig.savefig('transfer_'+name[jd]+'.png',format='png')

                  ## Barotropic KE
                  ## bke is barotropic component
                  ## beke is barotropic but with zonal mean removed
                  bke_2D  = calculate_ke(kx,ky,um*win,vm*win)
                  beke_2D = calculate_ke(kx,ky,um*win,vm*win)

                  ## Baroclinic KE
                  print ' this is the thickness vector to calculate_cke '
                  print ' is it ok? ',dz
                  sys.exit()
                  cke_2D  = calculate_cke(kx,ky,ut*win,vt*win,dz)
                  ceke_2D = calculate_cke(kx,ky,ut*win,vt*win,dz)
                  cke_2D = cke_2D
                  ceke_2D = ceke_2D
                  
                  ## triad interactions
                  Tk_data = calculate_spectral_flux_baroclinic_barotropic5(kx,ky,um*win,vm*win,ut*win,vt*win,dzt)
                  
                  adv_bt_bc_bc = Tk_data['Tk_bt_bc_bc']  ## KE bt -> KE bc, barotropisation
                  adv_bt_bc_bt = Tk_data['Tk_bt_bc_bt']  ##
                  adv_bt_bt_bc = Tk_data['Tk_bt_bt_bc']  ##
                  adv_bt_bt_bt = Tk_data['Tk_bt_bt_bt']  ## barotropic self-interaction

                  adv_bc_bc_bc = Tk_data['Tk_bc_bc_bc']  ## baroclinic self-interaction
                  adv_bc_bc_bt = Tk_data['Tk_bc_bc_bt']  ## catalytic flux
                  adv_bc_bt_bt = Tk_data['Tk_bc_bt_bt']  ##
                  adv_bc_bt_bc = Tk_data['Tk_bc_bt_bc']  ##
                  

                  ## APE advection
                  ## Not working very well!
                  Pk_data = calculate_spectral_ape_flux_baroclinic_barotropic(kx,ky,um*win,vm*win,ut*win,vt*win,ape)
                  adv_bt_ape = Pk_data['adv_ape_bt_ape'] ## Adv. of APE by barotropic mode
                  adv_bc_ape = Pk_data['adv_ape_bc_ape'] ## Adv. of APE by baroclinic mode

               else:
                  print ' To calculate barotropic and baroclinic components '
                  print ' you must either specify lfull_levels, or ltwo_levels '
                  sys.exit()
               
               
               ## calculate first order, second order, and fourth order gradients
               ## for barotropic and baroclinic velocities
               grad_data_bt = calculate_uv_gradients_xy(um[np.newaxis,:,:],vm[np.newaxis,:,:],\
                                                                  dxu,dyu,\
                                                                  dxv,dyv,\
                                                                  dxt,dyt,\
                                                                  dxf,dyf)
               grad_data_bc = calculate_uv_gradients_xy(ut,vt,\
                                                                  dxu,dyu,\
                                                                  dxv,dyv,\
                                                                  dxt,dyt,\
                                                                  dxf,dyf)                                                   
               ## Viscosity tendencies in barotropic/baroclinic modes
               if data['visc_form'] == 'bilaplacian':
                  um_visc = data['Ahm0'] * grad_data_bt['blpu_xy'][0,:,:]
                  vm_visc = data['Ahm0'] * grad_data_bt['blpv_xy'][0,:,:]
                  ut_visc = data['Ahm0'] * grad_data_bc['blpu_xy'][:,:,:]
                  vt_visc = data['Ahm0'] * grad_data_bc['blpv_xy'][:,:,:]
                  
               elif data['visc_form'] == 'laplacian':
                  um_visc = data['Ahm0'] * grad_data_bt['lapu_xy'][0,:,:]
                  vm_visc = data['Ahm0'] * grad_data_bt['lapv_xy'][0,:,:]
                  ut_visc = data['Ahm0'] * grad_data_bc['lapu_xy'][0,:,:]
                  vt_visc = data['Ahm0'] * grad_data_bc['lapv_xy'][0,:,:]
               else:
                  print ' the visc_form = ',data['visc_form']
                  print ' is not recognised so no horizontal viscosity for '
                  print ' barotropic/baroclinic KE can not be calculated '
                  sys.exit()
               

            if lvertical:
               ## Vertical momentum flux between top and second layer
               ## This is very experimental! 
               ## Results may be crap. 
               if grid['nz'] > 1:
                  if (ltend):
                     ddz_uhat = fftn(utend_vert[0,kk,:,:])
                     ddz_vhat = fftn(vtend_vert[0,kk,:,:])
                  else:
                     ## w * du/dz
                     ddz_u = -wvel[0,:,:] * (uvel[0,:,:]-uvel[1,:,:])/(deptht_tmp[0,:,:]-deptht_tmp[1,:,:])
                     ddz_v = -wvel[0,:,:] * (vvel[0,:,:]-vvel[1,:,:])/(deptht_tmp[0,:,:]-deptht_tmp[1,:,:])
                     ddz_uhat = fftn(ddz_u)
                     ddz_vhat = fftn(ddz_v)
                  ## u * w * du/dz
                  ke_vert = np.conj(fftn(uvel[0,:,:])) * ddz_uhat + \
                            np.conj(fftn(vvel[0,:,:])) * ddz_vhat

                  nn = ddz_uhat.shape[0]**2 * ddz_vhat.shape[1]**2
                  ke_vert = ke_vert / float(nn)

               ## Vertical viscosity
               ## Here I have assumed vertical viscosity is
               ## d/dz (A_z du/dz)
               ## But A_z is spatially non-uniform, so it must have been stored
               ## You could maybe make something up? Some profile from climatology? 
               if (grid['nz'] > 1) and lvvisc:
                  if (ltend):
                     u_visc = utend_vvsc[0,kk,:,:]
                     v_visc = vtend_vvsc[0,kk,:,:]
                  else:
                     u_visc = (Avm_s*dudz[0,:,:] - Avm_b*dudz[1,:,:]) / (deptht_tmp[0,:,:] - deptht_tmp[1,:,:])
                     v_visc = (Avm_s*dvdz[0,:,:] - Avm_b*dvdz[1,:,:]) / (deptht_tmp[0,:,:] - deptht_tmp[1,:,:])
                  ke_vvisc = np.conj( fftn(uvel[0,:,:]) ) * fftn(u_visc) + \
                             np.conj( fftn(vvel[0,:,:]) ) * fftn(v_visc)
                  nn = u_visc.shape[0]**2 * u_visc.shape[1]**2
                  ke_vvisc = ke_vvisc / float(nn)
            else:
               ke_vert  = np.zeros((grid['ny'],grid['nx']))
               ke_vvisc = np.zeros((grid['ny'],grid['nx']))

            ## Advection of APE
            ## This is something I am very unsure of
            ## Requires a lot more thinking and testing
            if lape and lrho and grid['nz'] > 1:
               uhat = fftn(uvel)
               vhat = fftn(vvel)
               if (0):
                  ## Calculate APE in physical space
                  ## integrate over all levels
                  ape  = (bprim[0,:,:]/N1[0,:,:])**2 * dzt[0,:,:] / np.sum(dzt,axis=0)
                  for jk in range(1,bprim.shape[0]):
                     ape += (bprim[jk,:,:]/N1[jk,:,:])**2 * dzt[jk,:,:] / np.sum(dzt,axis=0)
                  ## Take square root and Fourier transform
                  ## Then APE in Fourier space is np.conj(ahat) * ahat
                  ## I think this is wrong
                  ## Best way is probably to have b' and N1 at each level
                  ## then ahat = FFT(b'/N1) at each level
                  ## Then calculate np.conj(ahat) * ahat and integrate vertically
                  ape1 = np.sqrt(ape)
                  ahat = fftn(ape1*win)
                  i = np.complex(0,1)
                  # dAPE/dx in x,y
                  ddx_ape = np.real( ifftn(i*kx*ahat) )
                  # dAPE/dy in x,y
                  ddy_ape = np.real( ifftn(i*ky*ahat) )

                  # adv_ape = u * dAPE/dx + v * dAPE/dy
                  adv_ape = (uvel[0,:,:] * ddx_ape + vvel[0,:,:] * ddy_ape) * dzt[0,:,:]/dzt.sum(axis=0) + \
                            (uvel[1,:,:] * ddx_ape + vvel[1,:,:] * ddy_ape) * dzt[1,:,:]/dzt.sum(axis=0)

                  # KE trend from advection = - ape * adv_ape
                  # in spectral space
                  # The minus sign arises as advection
                  # is on the RHS of the momentum eqs.
                  ape_adv = np.real( -np.conj(ahat)*fftn(adv_ape) )   #[m2/s3]

               ## following Molemaker & McWilliams (J Fluid Mech, 2010)
               ## This is the best method I've found so far
               ## But the code below is still very experimental
               if (1):
                  for kk in range(0,bprim.shape[0]):
                     #atmp = np.sign(bprim[kk,:,:]) * np.sqrt(2. * ape)
                     atmp = bprim[kk,:,:]
                     
                     if diag_level >= 1:
                        fig = plt.figure()
                        ax1 = fig.add_subplot(221)
                        ax1.set_title('bprim')
                        cf1 = ax1.contourf(atmp)
                        plt.colorbar(cf1,ax=ax1)
                        fig.savefig(pdir+'/bprim.png',format='png')

                     ahat = fftn(atmp)
                     i = np.complex(0,1)
                     # dAPE/dx in x,y
                     ddx_ape = np.real( ifftn(i*kx*ahat) )
                     # dAPE/dy in x,y
                     ddy_ape = np.real( ifftn(i*ky*ahat) )

                     # adv_ape = u * dAPE/dx + v * dAPE/dy
                     adv_tmp = (uvel[kk,:,:] * ddx_ape + vvel[kk,:,:] * ddy_ape) #* dzt[0,:,:]/dzt.sum(axis=0)

                     # KE trend from advection = - ape * adv_ape
                     # in spectral space
                     # The minus sign arises as advection
                     # is on the RHS of the momentum eqs.
                     ape_adv = np.real( -np.conj(ahat)*fftn(adv_tmp) )   #[m2/s3]

                     uhat = fftn(uvel[kk,:,:]*win)
                     vhat = fftn(vvel[kk,:,:]*win)
                     i = np.complex(0,1)
                     # du/dx in x,y
                     ddx_u = np.real( ifftn(i*kx*uhat) )
                     ddx_v = np.real( ifftn(i*kx*vhat) )
                     # du/dy in x,y
                     ddy_u = np.real( ifftn(i*ky*uhat) )
                     ddy_v = np.real( ifftn(i*ky*vhat) )

                     # adv_ape = u * du/dx + v * du/dy
                     adv_u = (uvel[kk,:,:] * ddx_u + vvel[kk,:,:] * ddy_u) #* dzt[0,:,:]/dzt.sum(axis=0)

               
            ##
            ## Here we calculate all the KE fluxes in spectral space
            ## and vertically integrate if needed
            ##
            ## The above code has given us 2D arrays of tendency or KE 
            ## at each wavenumber. 
            ## We now assume isotropy, and make cumulative sums. 
            ## I.e. psd_ke(K) = sum KE(k,l) where k**2+l**2 > K**2
            ## This gives us psd_ke(K) as integral energy in units [m2/s2]
            ## and spectral fluxes Pi_ke(K) in units [m2/s3]
            ## The spectral density KE or transfer at wavenumber K is then
            ## d/dk psd_ke in units [m3/s2]
            ## and spectral transfer is Tk_ke = d/dk Pi_ke in units [m3/s3]
            ## 
            
            ## Mean KE in physical space
            ke_xy = 0.5 * np.mean(uvel[0,:,:]**2 + vvel[0,:,:]**2)
            
            ## Enstrophy and enstrophy flux
            ## Probably needs more work.
            ens_2D = calculate_ens(kx,ky,vort[0,:,:]*win)
            ens_adv = calculate_spectral_ens_flux(kx,ky,uvel[0,:,:]*win,vvel[0,:,:]*win,vort[0,:,:]*win)

            ## Power spectral density of KE
            psd_2D = 0.
            if ltwo_levels:
               nz = 1
            else:
               nz = uvel.shape[0]
            for kk in range(0,uvel.shape[0]):
               psd_tmp = calculate_ke(kx,ky,uvel[kk,:,:]*win,vvel[kk,:,:]*win)
               psd_2D += psd_tmp * dz[kk]/dz.sum()
            
            ## Spectral flux (i.e. momentum advection in spectral space)
            ## Done in function calculate_spectral_ke_tendency
            ##
            ## We do: fft(u) * fft(du/dx) 
            ## Not the same as fft(u) * i * k * fft(u) as in Scott & Wang 2005
            
            ke_adv = 0.
            ke_adv_tend = 0
            Hdep = np.sum(dzt,axis=0)
            if ltwo_levels: 
               nz = 1
            else:
               nz = uvel.shape[0]
               
            for kk in range(0,nz):
               adv_tmp = calculate_spectral_ke_tendency(uvel[kk,:,:], vvel[kk,:,:],\
                                                        u_adv[kk,:,:],v_adv[kk,:,:],win=win)
               ke_adv += adv_tmp * dz[kk]/dz.sum()
               if ltend:
                  adv_tmp = calculate_spectral_ke_tendency(uvel[kk,:,:], vvel[kk,:,:],\
                                                           utend_adv[kk,:,:],vtend_adv[kk,:,:],win=win)
                  ## add to vertical sum and multiply by layer depth over total depth                                         
                  ke_adv_tend += adv_tmp * dz[kk]/dz.sum()
            
            ##
            ## Forcing terms in spectral space
            ## 
            ## For wind forcing, we do: fft(u) * fft(tau)
            ## For viscosity, fft(u) * fft(A [d2/dx2+d2/dy2]u)
            ##
            ## Wind forcing u,v tendencies have been multiplied by rho/H
            ## where H is the depth we are considering. 
            ## 
            
            if ltau:
               ke_tau = calculate_spectral_ke_tendency(uvel[0,:,:], vvel[0,:,:],\
                                                       u_tau[:,:],v_tau[:,:],win=win)
               if (ltwo_levels):
                  ke_tau = ke_tau * 1 

            ## Viscosity for barotropic/baroclinic modes
            if lvsc:
               ke_visc = 0.
               ke_visc_tend = 0.
               if lbtbc:
                  ke_bt_visc = 0.
                  ## Horizontal viscosity for barotropic mode
                  ke_bt_visc += calculate_spectral_ke_tendency(um[:,:], vm[:,:],\
                                                               um_visc[:,:],vm_visc[:,:],win=win)
                  ## Horizontal viscosity for baroclinic modes
                  ke_bc_visc = 0.
                  for kk in range(0,ut.shape[0]):
                     ke_bc_visc += calculate_spectral_ke_tendency(ut[kk,:,:], vt[kk,:,:],\
                                                        ut_visc[kk,:,:],vt_visc[kk,:,:],win=win)
               
               ## Viscosity for full velocity fields   
               if ltwo_levels:
                  nz = 1
               else:
                  nz = uvel.shape[0]
               
               for kk in range(0,uvel.shape[0]):
                  visc_tmp = calculate_spectral_ke_tendency(uvel[kk,:,:], vvel[kk,:,:],\
                                                        u_vsc[kk,:,:],v_vsc[kk,:,:],win=win) 
                  ke_visc += visc_tmp * dz[kk]/dz.sum()
                  if ltend:
                     visc_tmp = calculate_spectral_ke_tendency(uvel[kk,:,:], vvel[kk,:,:],\
                                                        utend_vsc[kk,:,:],vtend_vsc[kk,:,:],win=win) 
                     ke_visc += visc_tmp * dz[kk]/dz.sum()
            
            ## Bottom drag
            if lbfr:
               if 1:
                  uhat = fftn(u_bot[:,:]*win)
                  vhat = fftn(v_bot[:,:]*win)
                  ## u_bfr and v_bfr are Cd |u| u / H, so ke_bfric is a vertical mean quantity
                  ke_bfric = calculate_spectral_ke_tendency(u_bot[:,:], v_bot[:,:],\
                                                            u_bfr[:,:], v_bfr[:,:],win=win) 
                  if lbtbc:
                     ## barotropic flux from bottom friction
                     ## [u] * [u_bfr], where [] is 1/H int()dz
                     dep = np.sum(dzt,axis=0)
                     ## the tendency in u is u_bfr/dz_bot
                     ## the vertical mean then becomes u_bfr/dz_bot * dz_bot / dep
                     ## so we skip dz_bot alltogether
                     ufr = u_bfr[:,:] / dep[:,:]
                     vfr = v_bfr[:,:] / dep[:,:]
                     ke_bfric = calculate_spectral_ke_tendency(um[:,:], vm[:,:],\
                                                               ufr[:,:],vfr[:,:],win=win)
            else:
               ke_bfric = 0.

            ## Planetary vorticity
            if 1:
               ke_cor = 0
               ke_cor_tend = 0
               if ltwo_levels:
                  nz = 1
               else:
                  nz = uvel.shape[0]
               if 1:
                  for kk in range(0,uvel.shape[0]):
                     coeff = dz[kk]/dz.sum()
                     cor_tmp = calculate_spectral_ke_tendency(uvel[kk,:,:], vvel[kk,:,:],\
                                                        u_cor[kk,:,:],v_cor[kk,:,:],win=win) 
                     ke_cor += cor_tmp * coeff
                     if (ltend):
                        cor_tmp = calculate_spectral_ke_tendency(uvel[kk,:,:], vvel[kk,:,:],\
                                                        utend_cor[kk,:,:],vtend_cor[kk,:,:],win=win) 
                        ke_cor_tend += cor_tmp * coeff

            ## Sea surface height (surface pressure gradient)
            if lssh:
               ke_spg = 0
               ke_spg_tend = 0
               if ltwo_levels:
                  nz = 1
               else:
                  nz = uvel.shape[0]
               for kk in range(0,uvel.shape[0]):
                  spg_tmp = calculate_spectral_ke_tendency(uvel[kk,:,:], vvel[kk,:,:],\
                                                           u_spg[kk,:,:],v_spg[kk,:,:],win=win) 
                  ke_spg += spg_tmp * coeff
                  if (ltend):
                     spg_tmp = calculate_spectral_ke_tendency(uvel[kk,:,:], vvel[kk,:,:],\
                                                           utend_spg[kk,:,:],vtend_spg[kk,:,:],win=win) 
                     ke_spg_tend += spg_tmp * coeff
                     
            ## Hydrostatic pressure gradient (only for  baroclinic models)
            if lrho:
               ke_hpg = 0
               ke_hpg_tend = 0
               if ltwo_levels:
                  nz = 1
               else:
                  nz = uvel.shape[0]
               for kk in range(0,uvel.shape[0]):
                  hpg_tmp = calculate_spectral_ke_tendency(uvel[kk,:,:], vvel[kk,:,:],\
                                                           u_hpg[kk,:,:],v_hpg[kk,:,:],win=win) 
                  ke_hpg += hpg_tmp * coeff
                  if ltend:
                     hpg_tmp = calculate_spectral_ke_tendency(uvel[kk,:,:], vvel[kk,:,:],\
                                                              utend_hpg[kk,:,:],vtend_hpg[kk,:,:],win=win) 
                     ke_hpg_tend += hpg_tmp * coeff
               
            ## Total pressure gradient (i.e. PE->KE conversion)
            if lrho and lssh:
               ke_pre = 0
               ke_pre_tend = 0
               if ltwo_levels:
                  nz = 1
               else:
                  nz = uvel.shape[0]
               for kk in range(0,uvel.shape[0]):
                  pre_tmp = calculate_spectral_ke_tendency(uvel[kk,:,:], vvel[kk,:,:],\
                                                           u_pre[kk,:,:],v_pre[kk,:,:],win=win) 
                  ke_pre += pre_tmp * coeff
                  if (ltend):
                     pre_tmp = calculate_spectral_ke_tendency(uvel[kk,:,:], vvel[kk,:,:],\
                                                           utend_pre[kk,:,:],vtend_pre[kk,:,:],win=win) 
                     ke_pre_tend += pre_tmp * coeff
               
               
            ## APE
            if lape:
               if (ltwo_levels):
                  for jk in range(0,2):
                     #apetmp = bprim[jk,:,:]/N1[jk,:,:] * np.sqrt(dzt[jk,:,:]/np.sum(dzt,axis=0)) * win
                     apetmp = bprim[jk,:,:]
                     ahat = fftn(apetmp)
                     z = 0.5 * np.real(ahat * np.conj(ahat))
                     nn = (ape.shape[1]**2 * ape.shape[0]**2)
                     if (jk==0):
                        pe2D  = z/float(nn)
                     else:
                        pe2D += z/float(nn)

               if (lfull_levels):
                  for jk in range(0,bprim.shape[0]-2):
                     tmp = 0.5 * (bprim[jk,:,:]**2)/(N1[jk,:,:]**2) #* dzt[jk,:,:]/np.ma.sum(dzt,axis=0)
                     apetmp = np.sign(bprim[jk,:,:]) * np.ma.sqrt( 2*tmp )
                     ahat = fftn(apetmp)

                     if (0):
                        fig = plt.figure()
                        ax1 = fig.add_subplot(221)
                        ax1.set_title('ape')
                        cf1 = ax1.contourf(apetmp)
                        plt.colorbar(cf1,ax=ax1)

                        ax2 = fig.add_subplot(222)
                        ax2.set_title('dz/H')
                        cf2 = ax2.contourf(dzt[jk,:,:]/np.ma.sum(dzt,axis=0))
                        plt.colorbar(cf2,ax=ax2)

                        ax3 = fig.add_subplot(223)
                        ax3.set_title('bprim^2')
                        cf3 = ax3.contourf(bprim[jk,:,:]**2)
                        plt.colorbar(cf3,ax=ax3)

                        ax4 = fig.add_subplot(224)
                        ax4.set_title('N^2')
                        cf4 = ax4.contourf(N1[jk,:,:]**2)
                        plt.colorbar(cf4,ax=ax4)

                     if (jk==0):
                        z = np.conj(ahat) * ahat
                     else:
                        z += np.conj(ahat) * ahat

                  nn = (ape.shape[1]**2 * ape.shape[0]**2)
                  pe2D  = z/float(nn)


            ## Conversion APE to PE
            if lape:
               if (ltwo_levels):
                  for kk in range(0,wvel.shape[0]):
                     what = fftn(wvel[kk,:,:]*win)
                     bhat = fftn(bprim[kk,:,:]*win)
                     ## bouyancy is defined as b = -g*(rho-rho_mean)/rho0
                     ## conversion APE -> KE is then w * b
                     z = np.real(np.conj(what) * bhat)
                     nn = (nx**2 * ny**2)
                     pe2ke =  z/float(nn) * dz[kk]/dz.sum()
                     ke2pe = -z/float(nn) * dz[kk]/dz.sum()
                     print ' calculated w*b at level '

               else:
                  pe2ke = 0.
                  ke2pe = 0.
                  for kk in range(0,wvel.shape[0]):
                     if (ltend):
                        z = np.conj(fftn(uvel[kk,:,:])) * fftn(utend_spg[kk,:,:]) +\
                            np.conj(fftn(vvel[kk,:,:])) * fftn(vtend_spg[kk,:,:])
                     else:
                        what = fftn(wvel[kk,:,:]*win)
                        bhat = fftn(bprim[kk,:,:]*win)
                        ## bouyancy is defined as b = -g*(rho-rho_mean)/rho0
                        ## conversion APE -> KE is then w * b
                        ## so here I ignore both minus signs
                        z = np.real(np.conj(what) * bhat)
                     nn = (nx**2 * ny**2)
                     pe2ke +=  z/float(nn) * dz[kk]/dz.sum()
                     ke2pe += -z/float(nn) * dz[kk]/dz.sum()
                     print ' calculated grad(p) at level ',kk

                  #fig,ax = plt.subplots(2,2)
                  #kk = np.array([[2,7],[15,25]])
                  #for j in range(0,2):
                  #   for i in range(0,2):
                  #      ax[i,j].set_title('w(%2d,:,:)' % (kk[i,j],))
                  #      cf = ax[i,j].contourf(wvel[kk[i,j],:,:])
                  #      plt.colorbar(cf,ax=ax[i,j])
                  #
                  #plt.show()
                  #sys.exit()


               #elif (lfull_levels):
               #   for jk in range(0,wvel.shape[0]):
               #      what = fftn(wvel[jk,:,:]*win)
               #      ## bouyancy is defined as b = -g*(rho-rho_mean)/rho0
               #      ## conversion APE -> KE is then -w * b
               #      ## so here I ignore both minus signs
               #      #btmp = np.sum(bprim * dzt,axis=0) / np.sum(dzt,axis=0)
               #      btmp = bprim[jk,:,:] #* dzt[0,:,:]/dzt.sum(axis=0)
               #      bhat = fftn(btmp*win)
               #      if (jk == 0):
               #         z = np.real(np.conj(what) * bhat) * dzt[jk,:,:]
               #      else:
               #         z += np.real(np.conj(what) * bhat) * dzt[jk,:,:]
               #
               #      #what = fftn(wvel[1,:,:]*win)
               #      #btmp = bprim[1,:,:] * dzt[1,:,:]/dzt.sum(axis=0)
               #      #bhat = fftn(btmp*win)
               #      #z += np.real(np.conj(what) * bhat)
               #      nn = (what.shape[1]**2 * what.shape[0]**2)
               #
               #   pe2ke =  z/(float(nn) * dzt.sum(axis=0))
               #   ke2pe = -z/(float(nn) * dzt.sum(axis=0))


            ## Diffusion of APE
            if lape:
               #ape  = (bprim[0,:,:]/N1[0,:,:])**2 * dzt[0,:,:] / np.sum(dzt,axis=0)
               ape  = bprim[0,:,:]
               for jk in range(1,bprim.shape[0]):
                  #ape += (bprim[jk,:,:]/N1[jk,:,:])**2 * dzt[jk,:,:] / np.sum(dzt,axis=0)
                  ape += bprim[jk,:,:]
               #ape = np.sqrt(ape)
               ahat = fftn(ape * win)
               # d2APE/dx2 in x,y
               ddx_APE = np.real( ifftn(-Aht0 * kx**2 * ahat) )
               # d2APE/dx2 in x,y
               ddy_APE = np.real( ifftn(-Aht0 * ky**2 * ahat) )
               # APE trend from diffusion:
               # ddx_u + ddy_u
               # in spectral space
               Tkxky = np.real( np.conj(ahat) * fftn(ddx_APE) + \
                                np.conj(ahat) * fftn(ddy_APE)   )

               nn = (ape.shape[1]**2 * ape.shape[0]**2)
               ape_visc = Tkxky / float(nn)


            ## Budget in physical space
            if (1):
               adv_xy = np.ma.zeros((uvel.shape[0]))
               tau_xy = np.ma.zeros((uvel.shape[0]))
               bfr_xy = np.ma.zeros((uvel.shape[0]))
               visc_xy = np.ma.zeros((uvel.shape[0]))
               pre_xy = np.ma.zeros((uvel.shape[0]))
               cor_xy = np.ma.zeros((uvel.shape[0]))
               p2k_xy = np.ma.zeros((uvel.shape[0]))

               adv_tend_xy = np.ma.zeros((uvel.shape[0]))
               visc_tend_xy = np.ma.zeros((uvel.shape[0]))
               pre_tend_xy = np.ma.zeros((uvel.shape[0]))
               cor_tend_xy = np.ma.zeros((uvel.shape[0]))

               for kk in range(0,uvel.shape[0]):
                  ## adv
                  if (ltend):
                     tmp = uvel[kk,:,:] * utend_adv[kk,:,:] + vvel[kk,:,:] * vtend_adv[kk,:,:]
                     adv_tend_xy[kk] = np.mean(tmp) * dz[kk]/dz.sum()

                  tmp = uvel[kk,:,:] * u_adv[kk,:,:] + vvel[kk,:,:] * v_adv[kk,:,:]
                  adv_xy[kk] = np.mean(tmp) * dz[kk]/dz.sum()

                  ## cor
                  if (ltend):
                     tmp = uvel[kk,:,:] * utend_cor[kk,:,:] + vvel[kk,:,:] * vtend_cor[kk,:,:]
                     cor_tend_xy[kk] = np.mean(tmp) * dz[kk]/dz.sum()

                  tmp = uvel[kk,:,:] * u_cor[kk,:,:] + vvel[kk,:,:] * v_cor[kk,:,:]
                  cor_xy[kk] = np.mean(tmp) * dz[kk]/dz.sum()

                  if ltau:
                     ## wind
                     tmp = uvel[0,:,:] * u_tau[:,:] + vvel[0,:,:] * v_tau[:,:]
                     tau_xy[kk] = np.mean(tmp) * dz[kk]/dz.sum()

                     ## bfr
                     tmp = u_bot[:,:] * u_bfr[:,:] + v_bot[:,:] * v_bfr[:,:]
                     bfr_xy[kk] = np.mean(tmp) * dz[kk]/dz.sum()
                     
                  ## visc
                  if (ltend):
                     tmp = uvel[kk,:,:] * utend_vsc[kk,:,:] + vvel[kk,:,:] * vtend_vsc[kk,:,:]
                     visc_tend_xy[kk] = np.mean(tmp) * dz[kk]/dz.sum()

                  tmp = uvel[kk,:,:] * u_vsc[kk,:,:] + vvel[kk,:,:] * v_vsc[kk,:,:]
                  visc_xy[kk] = np.mean(tmp) * dz[kk]/dz.sum()

                  ## pre
                  if lpe:
                     if (ltend):
                        tmp = uvel[kk,:,:] * utend_pre[kk,:,:] + vvel[kk,:,:] * vtend_pre[kk,:,:]
                        pre_tend_xy[kk] = np.mean(tmp) * dz[kk]/dz.sum()
                     
                     tmp = uvel[kk,:,:] * u_pre[kk,:,:] + vvel[kk,:,:] * v_pre[kk,:,:]
                     pre_xy[kk] = np.mean(tmp) * dz[kk]/dz.sum()

               adv_xy = np.ma.sum(adv_xy,axis=0)
               if ltau:
                  tau_xy = np.ma.sum(tau_xy,axis=0)
                  bfr_xy = np.ma.sum(bfr_xy,axis=0)
               if lpe:
                  pre_xy = np.ma.sum(pre_xy,axis=0)
               cor_xy = np.ma.sum(cor_xy,axis=0)
               visc_xy = np.ma.sum(visc_xy,axis=0)
               
               adv_tend_xy = np.ma.sum(adv_tend_xy,axis=0)
               pre_tend_xy = np.ma.sum(pre_tend_xy,axis=0)
               cor_tend_xy = np.ma.sum(cor_tend_xy,axis=0)
               visc_tend_xy = np.ma.sum(visc_tend_xy,axis=0)
               

            ##
            ## Integrate 2D PSD around lines of constant k
            ##
            ## For spectral KE flux, we here define PI(k)
            ## which is the transfer of energy from all
            ## scales of wavenumber < k (large scales) to
            ## those of wavenumber > k (small scales)
            ## PI = \int_k^inf T(k) dk or
            ## PI = -\int_0^k T(k) dk
            ## but since we don't have the largest scales (smallest k)
            ## we need to integrate from infinity to k
            ## where T(k) is e.g. -u_conj*u*du/dx
            ## and dE/dt = u_conj*du/dt
            ## It follows that dE/dt = -dPI/dk
            ##

            #
            # psd_ke [m3/s2]
            # Tk_ke  [m3/s3]
            # Tk_tau [m3/s3]
            # Tk_visc[m3/s3]
            #
            # Spectral energy densiy budget
            # d/dt psd_ke = Tk_ke + Tk_tau + Tk_visc
            #
            # Pi_ke  [m2/s3]
            #
            # Spectral energy budget
            # d/dt KE = Pi_ke + Pi_tau + Pi_visc
            #
            nk = k.shape[0]
            
            vk,psd_ke,ke_1D   = integrate_spectrum(psd_2D,wvsq,k,dk)
            va,psa_ke,kea_1D  = integrate_spectrum_angle(psd_2D,kx,ky)
            vk,psd_ens,ens_1D = integrate_spectrum(ens_2D,wvsq,k,dk)
            vk,Tk_ke,Pi_ke    = integrate_spectrum(ke_adv,wvsq,k,dk)

            if (ltend):
               vk,Tk_ke_tend,Pi_ke_tend = integrate_spectrum(ke_adv_tend,wvsq,k,dk)
            vk,Tk_ens,Pi_ens  = integrate_spectrum(ens_adv,wvsq,k,dk)

            if 'Ahm0' in data.keys():
               vk,Tk_visc,Pi_visc = integrate_spectrum(ke_visc,wvsq,k,dk)
               if (ltend):
                  vk,Tk_visc_tend,Pi_visc_tend = integrate_spectrum(ke_visc,wvsq,k,dk)
            if ltau:
               vk,Tk_bfr,Pi_bfr   = integrate_spectrum(ke_bfric,wvsq,k,dk)
               vk,Tk_tau,Pi_tau   = integrate_spectrum(ke_tau,wvsq,k,dk)
            if lpe:
               vk,Tk_spg,Pi_spg   = integrate_spectrum(ke_spg,wvsq,k,dk)
               vk,Tk_pre,Pi_pre   = integrate_spectrum(ke_pre,wvsq,k,dk)
               if (ltend):
                  vk,Tk_spg_tend,Pi_spg_tend = integrate_spectrum(ke_spg_tend,wvsq,k,dk)
                  vk,Tk_pre_tend,Pi_pre_tend = integrate_spectrum(ke_pre_tend,wvsq,k,dk)
            
            vk,Tk_cor,Pi_cor   = integrate_spectrum(ke_cor,wvsq,k,dk)
            if (ltend):
               vk,Tk_cor_tend,Pi_cor_tend = integrate_spectrum(ke_cor_tend,wvsq,k,dk)
                  
            if lrho:
               vk,Tk_hpg,Pi_hpg   = integrate_spectrum(ke_hpg,wvsq,k,dk)
               if (ltend):
                  vk,Tk_hpg_tend,Pi_hpg_tend = integrate_spectrum(ke_hpg_tend,wvsq,k,dk)

            if lvertical:
               vk,Tk_vert,Pi_vert = integrate_spectrum(ke_vert,wvsq,k,dk)
               vk,Tk_vvisc,Pi_vvisc = integrate_spectrum(ke_vvisc,wvsq,k,dk)
               print ' Tk_vvisc',Tk_vvisc

            if lape:
               vk,Tk_wb,Pi_wb = integrate_spectrum(pe2ke,wvsq,k,dk)
               if (lpe and name[jd][0:10] != 'MITgcm-Rob'):
                  vk,psd_ape,ape_1D = integrate_spectrum(pe2D,wvsq,k,dk)
                  vk,Tk_ape,Pi_ape = integrate_spectrum(ape_adv,wvsq,k,dk)
                  vk,Tk_ape_visc,Pi_ape_visc = integrate_spectrum(ape_visc,wvsq,k,dk)

                  if (lbtbc):
                     vk,Tk_bc_ape,Pi_bc_ape = integrate_spectrum(adv_bc_ape,wvsq,k,dk)
                     vk,Tk_bt_ape,Pi_bt_ape = integrate_spectrum(adv_bt_ape,wvsq,k,dk)

            ## split Tk_ke into negative fluxes and positive fluxes
            vk,Tk_ke_neg,Tk_ke_pos,tmp1,tmp2 = integrate_spectrum2(ke_adv,wvsq,k,dk)

            if lbtbc and grid['nz'] > 1:
               vk,psd_bke ,bke  = integrate_spectrum(bke_2D ,wvsq,k,dk)
               vk,psd_beke,beke = integrate_spectrum(beke_2D,wvsq,k,dk)
               vk,psd_cke ,cke  = integrate_spectrum(cke_2D ,wvsq,k,dk)
               vk,psd_ceke,ceke = integrate_spectrum(ceke_2D,wvsq,k,dk)

               vk,Tk_bt_bc_bc,Pi_bt_bc_bc = integrate_spectrum(adv_bt_bc_bc,wvsq,k,dk)
               vk,Tk_bt_bt_bt,Pi_bt_bt_bt = integrate_spectrum(adv_bt_bt_bt,wvsq,k,dk)
               vk,Tk_bt_bc_bt,Pi_bt_bc_bt = integrate_spectrum(adv_bt_bc_bt,wvsq,k,dk)
               vk,Tk_bt_bt_bc,Pi_bt_bt_bc = integrate_spectrum(adv_bt_bt_bc,wvsq,k,dk)

               vk,Tk_bc_bc_bc,Pi_bc_bc_bc = integrate_spectrum(adv_bc_bc_bc,wvsq,k,dk)
               vk,Tk_bc_bc_bt,Pi_bc_bc_bt = integrate_spectrum(adv_bc_bc_bt,wvsq,k,dk)
               vk,Tk_bc_bt_bt,Pi_bc_bt_bt = integrate_spectrum(adv_bc_bt_bt,wvsq,k,dk)
               vk,Tk_bc_bt_bc,Pi_bc_bt_bc = integrate_spectrum(adv_bc_bt_bc,wvsq,k,dk)

               vk,Tk_bt_visc,Pi_bt_visc = integrate_spectrum(ke_bt_visc,wvsq,k,dk)
               vk,Tk_bc_visc,Pi_bc_visc = integrate_spectrum(ke_bc_visc,wvsq,k,dk)
            
            ens = np.zeros((nk))
            ens2 = np.zeros((nk))
            ## Calculate PI(k)
            for jk in range(0,nk):
               indices = np.where(wvsq >= k[jk]**2)
               ens[jk] = np.sum(ens_2D[indices])
               ens2[jk] = np.sum(ens_adv[indices])

            ens1D  = -(ens[1:] - ens[0:-1]) / dk
            ens_adv1D = -(ens2[1:]-ens2[0:-1]) / dk
            ens_adv1D2 = 0.5*(ens2[1:]+ens2[0:-1])

            ## Calculate time scale,
            tt_ke = np.zeros((nk-1))
            for jk in range(0,nk-1):
               tt_ke[jk] = 1./(np.sqrt(psd_ke[jk]*dk)*vk[jk]) / 86400.
            
            ## Open output stream
            outfile = outdir + '/psd_'+full_prefix+'_'+name_list[jd]+'_'+regions[jr]+'_%04d%02d%02d.nc' % (currentTime.year,currentTime.month,currentTime.day)
            nc = make_file(outfile,currentTime,vk.shape[0])
            
            ## Write isotropic wavenumber vector
            nc.variables['k'][:] = vk[:]
            ## write grid into 
            nc.variables['tlon'][:] = grid['tlon'][grid['ny']/2,:]
            nc.variables['tlat'][:] = grid['tlat'][:,grid['nx']/2]
            nc.variables['deptht'][:] = grid['dept_1d'][klevels_keep]
            
            if ltwo_levels:
               ## write interface depth
               nc.variables['interface_depth'] = deptht[k_interface,:,:].mean()
            if leddy:
               nc.eddy = 'Velocities with subtracted mean '
            else:
               nc.eddy = 'Absolute velocities '
            
            ## make a time array
            ## in fractions of years   
            date = currentTime
            tt = date.timetuple()
            d0 = datetime(currentTime.year  ,1,1)
            d1 = datetime(currentTime.year+1,1,1)
            ndays = d1-d0
            #idt[0] = currentTime.year + (tt.tm_yday-1)/float(ndays.days)
            nc.variables['time'][0] = date2num(	currentTime, nc.variables['time'].units, \
                                                             calendar=nc.variables['time'].calendar)
            if diag_level >= 1:
               # control print the time variable
               print ' Wrote time value ',nc.variables['time'][0],' which corresponds to '
               print ' actual time: ',num2date(nc.variables['time'][0],nc.variables['time'].units, \
                                            calendar=nc.variables['time'].calendar)
            
            ## year, mon, day
            nc.variables['year'][0]  = currentTime.year
            nc.variables['month'][0] = currentTime.month
            nc.variables['day'][0]   = currentTime.day

            ## write power spectrum of KE 
            ## and spectral transfers and fluxes
            nc.variables['psd_ke'][0,:] = psd_ke[:]
            nc.variables['psdt_ke'][0]   = ke_xy
            nc.variables['psd_ke_timescale'][0,:]  = tt_ke[:]
            nc.variables['Tk_ke'][0,:] = Tk_ke[:]
            nc.variables['Pi_ke'][0,:] = Pi_ke[:]
            nc.variables['adv_xy'][0] = adv_xy
            nc.variables['adv_sp'][0] = Pi_ke[0]
            if ltend:
               nc.variables['Tk_ke_tend'][0,:] = Tk_ke_tend[:]
               nc.variables['Pi_ke_tend'][0,:] = Pi_ke_tend[:]
               nc.variables['adv_tend_xy'][0] = adv_tend_xy
               nc.variables['adv_tend_sp'][0] = Pi_ke_tend[0]
            
            ## write enstrophy and enstrophy transfers   
            nc.variables['psd_ens'][0,:] = psd_ens[:]
            nc.variables['Tk_ens'][0,:] = Tk_ens[:]
            nc.variables['Pi_ens'][0,:] = Pi_ens[:]
            
            if lvsc:
               ## spectral transfer and flux from viscosity 
               nc.variables['Tk_ke_visc'][0,:] = Tk_visc[:]
               nc.variables['Pi_ke_visc'][0,:] = Pi_visc[:]
               nc.variables['visc_xy'][0] = visc_xy
               nc.variables['visc_sp'][0] = Pi_visc[0]
               
               if ltend:
                  nc.variables['Tk_ke_visc_tend'][0,:] = Tk_visc_tend[:]
                  nc.variables['Pi_ke_visc_tend'][0,:] = Pi_visc_tend[:]
                  nc.variables['visc_tend_xy'][0] = visc_tend_xy
                  nc.variables['visc_tend_sp'][0] = Pi_visc_tend[0]
                  
            if ltau:
               ## spectral transfer from wind forcing
               nc.variables['Tk_ke_tau'][0,:] = Tk_tau[:]
               nc.variables['Pi_ke_tau'][0,:] = Pi_tau[:]
               nc.variables['tau_xy'][0] = tau_xy
               nc.variables['tau_sp'][0] = Pi_tau[0]
            
            if lbfr:
               ## spectral transfer from bottom friction   
               nc.variables['Tk_bfr'][0,:] = Tk_bfr[:]
               nc.variables['Pi_bfr'][0,:] = Pi_bfr[:]
               nc.variables['bfr_xy'][0] = bfr_xy
               nc.variables['bfr_sp'][0] = Pi_bfr[0]

            if lssh:
               ## surface pressure gradient
               nc.variables['Tk_spg'][0,:] = Tk_spg
               nc.variables['Pi_spg'][0,:] = Pi_spg
               if ltend:
                  nc.variables['Tk_spg_tend'][0,:] = Tk_spg_tend
                  nc.variables['Pi_spg_tend'][0,:] = Pi_spg_tend
            
            if lrho:
               ## hydrostatic pressure gradient
               nc.variables['Tk_hpg'][0,:] = Tk_hpg
               nc.variables['Pi_hpg'][0,:] = Pi_hpg
               if ltend:
                  nc.variables['Tk_hpg_tend'][0,:] = Tk_hpg_tend
                  nc.variables['Pi_hpg_tend'][0,:] = Pi_hpg_tend
                     
            if lrho and lssh:
               ## total pressure gradient
               nc.variables['Tk_pre'][0,:] = Tk_pre
               nc.variables['Pi_pre'][0,:] = Pi_pre
               nc.variables['pre_xy'][0] = pre_xy
               nc.variables['pre_sp'][0] = Pi_pre[0]
               if ltend:
                  nc.variables['Tk_pre_tend'][0,:] = Tk_pre_tend
                  nc.variables['Pi_pre_tend'][0,:] = Pi_pre_tend
                  nc.variables['pre_tend_xy'][0] = pre_tend_xy
                  nc.variables['pre_tend_sp'][0] = Pi_pre_tend[0]
                  
            if 1:
               ## planetary vorticity
               nc.variables['Tk_cor'][0,:] = Tk_cor
               nc.variables['Pi_cor'][0,:] = Pi_cor
               if ltend:
                  nc.variables['Tk_cor_tend'][0,:] = Tk_cor_tend
                  nc.variables['Pi_cor_tend'][0,:] = Pi_cor_tend
                  
            if lape:
               nc.variables['p2k_xy'][0] = p2k_xy
               nc.variables['p2k_sp'][0] = Pi_wb[0]
               
            if lvertical:
               nc.variables['Tk_ke_vert'][0,:] = Tk_vert[:]
               nc.variables['Pi_ke_vert'][0,:] = Pi_vert[:]
               nc.variables['Tk_ke_vvisc'][0,:] = Tk_vvisc[:]
               nc.variables['Pi_ke_vvisc'][0,:] = Pi_vvisc[:]

            if lape:
               nc.variables['Tk_wb'][0,:] = Tk_wb[:]
               nc.variables['Pi_wb'][0,:] = Pi_wb[:]
               if 1:
                  ## power spectrum, spectral transfer from advection and diffusion
                  ## for APE. Also conversion APE-KE
                  nc.variables['psd_ape'][0,:] = psd_ape[:]
                  nc.variables['Tk_ape'][0,:] = Tk_ape[:]
                  nc.variables['Pi_ape'][0,:] = Pi_ape[:]
                  nc.variables['Tk_ape_visc'][0,:] = Tk_ape_visc[:]
                  nc.variables['Pi_ape_visc'][0,:] = Pi_ape_visc[:]
                  nc.variables['shear'][0] = ushear
                  nc.variables['rhoprim'][0] = rhoprim

                  if lbtbc:
                     ## spectral transfer and flux of KE in
                     ## barotropic and baroclinic modes
                     nc.variables['Tk_bc_ape'][0,:] = Tk_bc_ape[:]
                     nc.variables['Tk_bt_ape'][0,:] = Tk_bt_ape[:]
                     nc.variables['Pi_bc_ape'][0,:] = Pi_bc_ape[:]
                     nc.variables['Pi_bt_ape'][0,:] = Pi_bt_ape[:]

            if lbtbc and grid['nz'] >1:
               ## power spectrum of barotropic and baroclinic KE
               ## Also transfers and fluxes for triad interactions
               nc.variables['psd_bke'][0,:] = psd_bke[:]
               nc.variables['psd_beke'][0,:] = psd_beke[:]
               nc.variables['psd_cke'][0,:] = psd_cke[:]
               nc.variables['psd_ceke'][0,:] = psd_ceke[:]
               
               nc.variables['Tk_bt_bc_bc'][0,:] = Tk_bt_bc_bc[:]
               nc.variables['Tk_bt_bt_bt'][0,:] = Tk_bt_bt_bt[:]
               nc.variables['Tk_bt_bc_bt'][0,:] = Tk_bt_bc_bt[:]
               nc.variables['Tk_bt_bt_bc'][0,:] = Tk_bt_bt_bc[:]
               nc.variables['Tk_bc_bc_bc'][0,:] = Tk_bc_bc_bc[:]
               nc.variables['Tk_bc_bc_bt'][0,:] = Tk_bc_bc_bt[:]
               nc.variables['Tk_bc_bt_bt'][0,:] = Tk_bc_bt_bt[:]
               nc.variables['Tk_bc_bt_bc'][0,:] = Tk_bc_bt_bc[:]

               nc.variables['Pi_bt_bc_bc'][0,:] = Pi_bt_bc_bc[:]
               nc.variables['Pi_bt_bt_bt'][0,:] = Pi_bt_bt_bt[:]
               nc.variables['Pi_bt_bc_bt'][0,:] = Pi_bt_bc_bt[:]
               nc.variables['Pi_bt_bt_bc'][0,:] = Pi_bt_bt_bc[:]
               nc.variables['Pi_bc_bc_bc'][0,:] = Pi_bc_bc_bc[:]
               nc.variables['Pi_bc_bc_bt'][0,:] = Pi_bc_bc_bt[:]
               nc.variables['Pi_bc_bt_bt'][0,:] = Pi_bc_bt_bt[:]
               nc.variables['Pi_bc_bt_bc'][0,:] = Pi_bc_bt_bc[:]

               nc.variables['Tk_bt_visc'][0,:] = Tk_bt_visc[:]
               nc.variables['Pi_bt_visc'][0,:] = Pi_bt_visc[:]
               nc.variables['Tk_bc_visc'][0,:] = Tk_bc_visc[:]
               nc.variables['Pi_bc_visc'][0,:] = Pi_bc_visc[:]
            
            if 1:
               nc.variables['cross_uvel'][0,:,:] = cross_uvel[0,:,:]
            if lrho:
               nc.variables['cross_rho'][0,:,:] = cross_rho[0,:,:]

            ## Store means
            if currentTime == starttime:
               print ' Set averages to zero '
               
               data_mean = {}
               
               #data_mean['vk'] = vk[:]
               #data_mean['tlat'] = grid['tlat'][:,:]
               #data_mean['deptht'] = grid['dept_1d'][klevels_keep]
               data_mean['numsteps'] = 0
               
            data_mean['numsteps'] += 1
            
            varnames = nc.variables.keys()
               
            ## Concatenate all data along the time dimension
            ## (only used for plotting KE below)
            data_mean = store_nc_data(nc,data_mean)
            
            ## close file
            nc.close()
            
            
            ## Calculate total energy or flux in spectral space
            adv_sp = Pi_ke[0]
            if 'Ahm0' in data.keys():
               visc_sp = Pi_visc[0]
            if ltau:
               tau_sp = Pi_tau[0]
            if lssh and lrho:
               pre_sp = Pi_pre[0]
            if lape:
               p2k_sp = Pi_wb[0]
            
            ##
            ## Add time step to currentTime
            ##
            currentTime += outputStep
            
         ## end loop over steps in file
      ## end time loop
      
      ## Plot KE and fluxes if needed
      if diag_level >= 1:
         fig = plt.figure()
         ax1 = fig.add_subplot(111)
         ax1.set_title('KE ('+name_list[jd]+' '+regions[jr]+')')
         ax1.loglog(vk,data_mean['psd_ke'].mean(axis=0))
         
         fig = plt.figure()
         ax1 = fig.add_subplot(111)
         ax1.set_title('Spectral transfer ('+name_list[jd]+' '+regions[jr]+')')
         ax1.semilogx(vk,data_mean['Tk_ke'].mean(axis=0))
         
         fig = plt.figure()
         ax1 = fig.add_subplot(111)
         ax1.set_title('Spectral flux ('+name_list[jd]+' '+regions[jr]+')')
         ax1.semilogx(vk,data_mean['Pi_ke'].mean(axis=0))


print ' ========================================================== '
print ' '
print ' End of analysis for runs: '
for name in name_list:
   print ' '+name
print ' '
print ' Time: ',time.time()
print ' '
print ' Data has been stored in: '+outdir
print ' Calculated KE and spectral fluxes ' 
if ltau:
   print ' Calculated wind forcing '
if lbfr:
   print ' Calculated bottom friction '
if lvsc:
   print ' Calculate horizontal viscosity '
if lssh:
   print ' Calculated surface pressure gradient '
if lrho:
   print ' Calculated hydrostatic pressure gradient '
if lpe:
   print ' Calculated potential energy (EXPERIMENTAL!) '   
if lvertical:
   print ' Calculated vertical energy fluxes and viscosity (EXPERIMENTAL!)'
if lbtbc:
   print ' Calculated barotropic/baroclinic KE and spectral fluxes '
print ' '   
print ' ========================================================== '

plt.show()
sys.exit()


