import numpy as np
from netCDF4 import Dataset

def store_nc_data(nc,data_dict):
   """
   """
   
   for varname in nc.variables.keys():
      vardata = nc.variables[varname]
      
      ## if the data already has been stored once, and it has time as a coordinate
      ## then we concatenate along the time coordinate
      if (varname in data_dict.keys()) and ('t' in vardata.dimensions):
         data_dict[varname] = np.concatenate( (data_dict[varname],vardata[:]), axis=0 )
      else:
         data_dict[varname] = vardata[:]
         
   return data_dict


def read_stored_data(run_list,region_list,prefix,starttime,endtime,outputStep,ddir):
   """
   """
   
   ## number of data sets
   nd = len(run_list)
   ## number of regions
   nr = len(region_list)
   
   data_list = []
   for jd in range(0,nd):
      reg_list = []
      for jr in range(0,nr):
         currentTime = starttime
         while currentTime <= endtime:
            psdfile = ddir + '/psd_'+prefix+'_'+run_list[jd]+'_'+region_list[jr]+'_%04d%02d%02d.nc' % (currentTime.year,currentTime.month,currentTime.day)
            
            print psdfile
            nc = Dataset(psdfile,'r')
            
            ## Store means
            if currentTime == starttime:
               print ' Set averages to zero '
               
               data_mean = {}
               
               #data_mean['vk'] = vk[:]
               #data_mean['tlat'] = grid['tlat'][:,:]
               #data_mean['deptht'] = grid['dept_1d'][klevels_keep]
               data_mean['numsteps'] = 0
               data_mean['region'] = region_list[jr]
               data_mean['run'] = run_list[jd]
               
            data_mean['numsteps'] += 1
               
            ## Add all data above
            data_mean = store_nc_data(nc,data_mean)
            
            ## close file
            nc.close()
            
            currentTime += outputStep
            
         reg_list.append(data_mean)
      
      data_list.append(reg_list)
            
   return data_list
