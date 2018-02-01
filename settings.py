from datetime import datetime, timedelta

def setup_analysis():
   """
   Settings for spectral energy analysis
   """
   
   setup = {}
   
   ## Prefix for the analysis
   setup['prefix'] = 'test'
   
   ## Grid type, start time, end time, output time step
   setup['grid_type']  = 'll' ## type of grid (ll - lonlat, xy - cartesian)
   setup['starttime']  = datetime(2009,1,5)
   setup['endtime']    = datetime(2009,12,31)
   setup['outputStep'] = timedelta(days=5)
   
   ## Lists with simulation names, directories, and regions to study
   ## also where to store plots (pdir) and results (outdir)
   name_list = ['ORCA05.L46-KJH0004','INALT10.L46-KJH0017','INALT10.L46-KJH0017-NEST1']
   ddir_list = ['/Users/jkjellsson/data/ORCA05.L46/',\
                '/Users/jkjellsson/data/INALT10.L46/',\
                '/Users/jkjellsson/data/INALT10.L46/']
   regions   = ['agulhas-retro'] 
   pdir      = './'
   outdir    = './'
   
   setup['names']   = name_list
   setup['dirs']    = ddir_list
   setup['regions'] = regions
   setup['pdir']    = pdir
   setup['outdir']  = outdir
   
   ## What to calculate 
   setup['lcalculate'] = False  # do calculations and save to files
   setup['lrhines']    = False  ## calculate and plot Rhines scale
   setup['lrossby']    = True    ## read and plot 1st baroclinic Rossby radius
   setup['lpsd_freq']  = False  ## store each step to make frequency spectrum analysis 
   setup['lbtbc']      = False    ## barotropic and baroclinic components
   setup['lpe']        = False  ## calculate potential energy
   setup['lape']       = False ## APE calculations that dont work!
   setup['ltend']      = False   ## read full model tendencies
   setup['lreadw']     = True   ## read vertical velocity (otherwise infer from cont.)
   setup['lssh']       = False ## calculate ssh and surface pressur gradient
   setup['lrho']       = False ## calculate buoyancy and hydrostatic pressure gradient
   setup['ltau']       = True  ## calculate wind stress
   setup['lbfr']       = True  ## calculate bottom friction
   setup['lvsc']       = True  ## calculate horizontal viscosity
   setup['lvertical']  = False ## calculate vertical KE flux and viscosity (experimental!)
   
   ## How to calculate
   setup['ltukey']       = False  ## use a Tukey window to force zero on boundaries
   setup['leddy']        = False  ## remove mean to get eddy components
   setup['linterpolate'] = False  ## interpolate NEMO data to regular grid. May not be necessary
   setup['lrotate']      = False  ## rotate grid (experimental!)
   setup['ltrend']       = False # calculate time trends (not working)
   
   ## Diagnostics (=0 no diagnostics, =1 print some control stuff, some plots, =2 plot lots of stuff)
   setup['diag_level'] = 1 ## = 0 no diagnostics plots
                           ## = 1 plot KE spectra and map of sfc vorticity each step
   
   ## Depths to read in 
   setup['hmin'] = 0.   # shallowest level to use 
   setup['hmax'] = 500. # deepest point to use
   setup['hsep'] = 500. # depth to split into two levels 
   
   return setup