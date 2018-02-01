import os,sys
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from netCDF4 import Dataset

## control fonts for plots
from matplotlib import rc
rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)


def plot_mean_fields(grid,data,lcrossy=True,lcrossx=True,lbarotropic=True):
   """
   Plot some mean fields
   """
   
   ## Plot mean fields
   figname = pdir+'/'+name[jd]+'_'+regions[jr]
   umean_max = np.zeros((nr))
   vmean_max = np.zeros((nr))
   wmean_max = np.zeros((nr))
   tmean_max = np.zeros((nr))
   smean_max = np.zeros((nr))
   rmean_max = np.zeros((nr))
   
   u = np.ma.masked_where( (u==0) & (v==0), u )
   v = np.ma.masked_where( (u==0) & (v==0), v )
   w = np.ma.masked_where( (u==0) & (v==0), w )
   t = np.ma.masked_where( s==0, t )
   s = np.ma.masked_where( s==0, s )
   rho = np.ma.masked_where( s==0, rho )
   
   
   ## Plot vertical mean u
   zplot = np.ma.mean(u,axis=0)
   fig = plt.figure()
   ax1 = fig.add_subplot(111)
   ax1.set_xlabel('i index')
   ax1.set_ylabel('j index')
   ax1.set_title('Vertical mean zonal velocity')
   if (umax == 0): umax=zplot.max()
   levels = np.linspace(-umax/4.,umax/4.,21)
   cf1 = ax1.contourf(zplot,levels=levels,extend='both')
   plt.colorbar(cf1,ax=ax1)
   
   fig.savefig(name+'_vertical_mean_u.pdf',format='pdf')
   
   ## Plot area mean vertical structure
   vu = np.ma.array([])
   vv = np.ma.array([])
   vw = np.ma.array([])
   vt = np.ma.array([])
   vs = np.ma.array([])
   vr = np.ma.array([])
   vzt = np.ma.array([])
   vzw = np.ma.array([])
   for jk in range(0,u.shape[0]):
      vu = np.ma.append(vu,np.ma.mean(u[jk,:,:]))
      vv = np.ma.append(vv,np.ma.mean(v[jk,:,:]))
      vw = np.ma.append(vw,np.ma.mean(w[jk,:,:]))
      vt = np.ma.append(vt,np.ma.mean(t[jk,:,:]))
      vs = np.ma.append(vs,np.ma.mean(s[jk,:,:]))
      vr = np.ma.append(vr,np.ma.mean(rho[jk,:,:]))
      vzt = np.ma.append(vzt,np.ma.mean(tdep[jk,:,:]))
      vzw = np.ma.append(vzw,np.ma.mean(wdep[jk,:,:]))
   
   zticks = np.array([5,50,100,500,1000,2000,5000])
   zlabel = []
   for zt in zticks:
      zlabel.append('%4d' % zt)
         
   fig,ax = plt.subplots(2,3)
   ax[0,0].set_title('u')
   ax[0,0].semilogy(vu,vzt,basey=2)
   ax[0,1].set_title('v')
   ax[0,1].semilogy(vv,vzt,basey=2)
   ax[0,2].set_title('w')
   ax[0,2].semilogy(vw,vzw,basey=2)
   ax[1,0].set_title('T')
   ax[1,0].semilogy(vt,vzt,basey=2)
   ax[1,1].set_title('S')
   ax[1,1].semilogy(vs,vzt,basey=2)
   ax[1,2].set_title('rho')
   ax[1,2].semilogy(vr,vzt,basey=2)
   
   for jj in range(0,2):
      for ji in range(0,3):
         xlims = ax[jj,ji].get_xlim()
         xticks = np.array([xlims[0],0.5*(xlims[0]+xlims[1]),xlims[1]])
         ax[jj,ji].set_xticks(xticks)
         ax[jj,ji].invert_yaxis()
         ax[jj,ji].set_yticks(zticks)
         ax[jj,ji].set_yticklabels(zlabel)
   
   fig.savefig(name+'_uvwtsr_profiles.pdf',format='pdf')   
   
   
   ## Plot zonal mean u
   i1 = int(np.floor(u.shape[2]/2))
   zplot1 = u[:,:,:i1].mean(axis=2)
   zplot2 = rho[:,:,:i1].mean(axis=2)
   fig = plt.figure()
   ax1 = fig.add_subplot(111)
   ax1.set_title('Zonal mean zonal velocity')
   #cf1 = ax1.contourf(zplot1,levels=np.linspace(-umax/3.,umax/3.,21),extend='both')
   #ct1 = ax1.contour(zplot2,colors='k')
   cf1 = ax1.contourf(ulat[:,0],vzt,zplot1,levels=np.linspace(-umax/3.,umax/3.,21),extend='both')
   ct1 = ax1.contour( ulat[:,0],vzt,zplot2,colors='k')
   plt.colorbar(cf1,ax=ax1)
   plt.clabel(ct1, fontsize=9, inline=1)
   #zticks = np.arange(0,vzt.shape[0],5)
   #yticks = np.arange(0,ulat[:,0].shape[0],np.floor(ulat[:,0].shape[0]/5.))
   #zlabel = []
   #ylabel = []
   #for k in zticks:
   #   zlabel.append( '%4d' % (vzt[k],) )
   #for j in yticks:
   #   ylabel.append( '%6.2f' % (ulat[j,0],) )
   #ax1.set_yticks(zticks)
   #ax1.set_yticklabels(zlabel)
   ax1.set_ylim([0,2000])
   ax1.invert_yaxis()
   #ax1.set_xticks(yticks)
   #ax1.set_xticklabels(ylabel)
   fig.savefig(name+'_zonal_mean_u.pdf',format='pdf')
   
   
   ## Plot zonal mean u shear
   #zplot = np.ma.mean( -(u[0:-1,:,:]-u[1:,:,:])/(tdep[0:-1,:,:]-tdep[1:,:,:]),axis=2 )
   #levels = np.linspace(-0.006,0.006,21)
   #zz = 0.5*(vzt[1:]+vzt[0:-1])
   #fig = plt.figure()
   #ax1 = fig.add_subplot(111)
   #ax1.set_title('Zonal mean zonal velocity shear')
   #cf1 = ax1.contourf(ulat[:,0],zz,zplot,levels=levels,extend='both')
   #ax1.invert_yaxis()
   #ax1.set_yscale('log', nonposy='clip', basey=2)
   #ax1.set_yticks(zticks)
   #ax1.set_yticklabels(zlabel)
   #plt.colorbar(cf1,ax=ax1)
   #fig.savefig(name+'_zonal_mean_shear_u.pdf',format='pdf')
   
   ## Plot meridional cross section
   i1 = int(np.floor(u.shape[2]/2))
   zplot1 = u[:,:,i1]
   zplot2 = rho[:,:,i1]
   fig = plt.figure()
   ax1 = fig.add_subplot(111)
   cf1 = ax1.contourf(ulat[:,0],vzt,zplot1,levels=np.linspace(-umax,umax,21),extend='both')
   plt.colorbar(cf1,ax=ax1)
   ct1 = ax1.contour( ulat[:,0],vzt,zplot2,colors='k')
   plt.clabel(ct1, fontsize=9, inline=1)
   #zticks = np.arange(0,vzt.shape[0],5)
   #yticks = np.arange(0,ulat[:,0].shape[0],np.floor(ulat[:,0].shape[0]/5))
   #zlabel = []
   #ylabel = []
   #for k in zticks:
   #   zlabel.append( '%4d' % (vzt[k],) )
   #for j in yticks:
   #   ylabel.append( '%6.2f' % (ulat[j,0],) )
   #ax1.set_yticks(zticks)
   #ax1.set_yticklabels(zlabel)
   ax1.set_ylim([0,2000])
   ax1.invert_yaxis()
   #ax1.set_xticks(yticks)
   #ax1.set_xticklabels(ylabel)
   fig.savefig(name+'_shear_and_rho.pdf',format='pdf')
   

def plot_all_regions(grid,data,regions):
   """
   Plot a global map and draw the regions
   """
   
   ## use the new viridis colormap, if available
   if 'viridis' in plt.cm.datad.keys():
      cmap = plt.cm.viridis
   else:
      cmap = plt.cm.YlGnBu_r
   
   if 1:
      ## Plot map of KE and mark region of interest
      uvel = data['uvel'][0,0,:,:]
      vvel = data['vvel'][0,0,:,:]
      ke   = 0.5 * (uvel**2 + vvel**2)
      
      print ' Plot global map of surface KE '
      fig_ke_glo = plt.figure()
      ax_ke_glo  = fig_ke_glo.add_axes([0.1,0.2,0.6,0.6])
      cax_ke_glo = fig_ke_glo.add_axes([0.74,0.3,0.02,0.4])
      ax_ke_glo.set_title(r'Surface kinetic energy',fontsize=10)
      bmap = Basemap(projection='eck4',lon_0=0.,resolution='c',ax=ax_ke_glo)
      levs = np.array([0,0.001,0.002,0.005,0.01,0.02,0.05,0.1,0.2])
      levels_norm = np.linspace(0,1,levs.shape[0])
      ke_colors = cmap(levels_norm)
      
      x,y = bmap(grid['tlon'],grid['tlat'])
      x = x.flatten()
      y = y.flatten()
      z = ke.flatten()
      cf  = bmap.contourf(x,y,z,tri=True,extend='max',colors=ke_colors,levels=levs)
      bmap.fillcontinents(color='Gainsboro',lake_color=ke_colors[0])
      cb = fig_ke_glo.colorbar(cf,cax=cax_ke_glo,spacing='uniform')
      cb.set_label(r'Kinetic energy $[\mathrm{m}^2$ $\mathrm{s}^{-2}]$',fontsize=8)
      cb.ax.tick_params(labelsize=8)
      
      ## draw regions
      for region in regions:
         i0 = region['i0']
         i1 = region['i1']
         j0 = region['j0']
         j1 = region['j1']
         print ' draw region ',region['region'],i0,i1,j0,j1
         lam = grid['tlon'][j0,i0:i1] 
         phi = grid['tlat'][j0,i0:i1]
         x,y = bmap(lam,phi)
         bmap.plot(x,y,'-r')
         
         lam = grid['tlon'][j0:j1,i1]
         phi = grid['tlat'][j0:j1,i1]
         x,y = bmap(lam,phi)
         bmap.plot(x,y,'-r')
         
         lam = grid['tlon'][j1,i0:i1]
         phi = grid['tlat'][j1,i0:i1]
         x,y = bmap(lam,phi)
         bmap.plot(x,y,'-r')
         
         lam = grid['tlon'][j0:j1,i0]
         phi = grid['tlat'][j0:j1,i0]
         x,y = bmap(lam,phi)
         bmap.plot(x,y,'-r')
      
      ## save fig   
      fig_ke_glo.savefig('map_all_regions.png',format='png',dpi=600)
      
      if 0:
         if 1:
            if 1:
               ## Plot local Kuroshio map
               fig_ke_ks = plt.figure()
               ax_ke_ks  = fig_ke_ks.add_axes([0.1,0.2,0.6,0.6])
               cax_ke_ks = fig_ke_ks.add_axes([0.74,0.3,0.02,0.4])
               ax_ke_ks.set_title(r'Surface kinetic energy in '+name[jd][:-5],fontsize=10)
               #bmap = Basemap(llcrnrlon=135.,llcrnrlat=28,urcrnrlon=180,urcrnrlat=43,\
               #                  resolution='l',projection='merc',ax=ax_ke_ks)
               #m.drawcoastlines()
               bmap = Basemap(llcrnrlon=135.,llcrnrlat=28.,urcrnrlon=170.,urcrnrlat=43.,\
                           rsphere=(6378137.00,6356752.3142),\
                           resolution='i',projection='merc',\
                           lat_0=35.5,lon_0=150,lat_ts=20.,ax=ax_ke_ks)

               bmap.drawparallels(np.arange(10,90,5),labels=[0,1,0,1])
               bmap.drawmeridians(np.arange(-180,180,10),labels=[0,1,0,1])

               levs = np.array([0,0.001,0.002,0.005,0.01,0.02,0.05,0.1,0.2])
               levels_norm = np.linspace(0,1,levs.shape[0])
               ke_colors = cmap(levels_norm)
               ke = 0.5 * ( uplot**2 + vplot**2   )
               x,y = bmap(lon_glo,lat_glo)
               ke[ke.mask] = 0.
               x = x.flatten()
               y = y.flatten()
               z = ke.flatten()
               cf  = bmap.contourf(x,y,z,tri=True,extend='max',colors=ke_colors,levels=levs)
               bmap.fillcontinents(color='Gainsboro',lake_color=ke_colors[0])
               cb = fig_ke_ks.colorbar(cf,cax=cax_ke_ks,spacing='uniform')
               cb.set_label(r'Kinetic energy $[\mathrm{m}^2$ $\mathrm{s}^{-2}]$',fontsize=8)
               cb.ax.tick_params(labelsize=8)


def plot_stored_data(data_list,format='png'):
   """
   Plot the data stored in the data_list
   data_list is assumed to be generated by read_stored_data
   so that it has one element for each simulation analysed, 
   and each of the simulations has one element for each region analysed. 
   Within each of those elements, data is stored for all time steps
   
   So to get the density of KE from the first simulation and second region do
   psd = data_list[0][1]['psd_ke']
   which is then an array of size (t,k) where t is time and k isotropic wavenumber. 
   Storing all time steps means we can use time filtering etc. 
   """
   
   nd = len(data_list)
   nr = len(data_list[0])
   
   ## Loop over regions
   for jr in range(0,nr):
      
      ##
      ## Plot KE spectra
      ##
      
      fig = plt.figure()
      ax1 = fig.add_subplot(111)
      ax1.set_title('Spectral density of KE in '+data_list[0][jr]['region'])
      
      for jd in range(0,nd):
         vk  = data_list[jd][jr]['k']
         psd = data_list[jd][jr]['psd_ke'].mean(axis=0)
         print psd.shape
         ax1.loglog(vk,psd,label=data_list[jd][jr]['run'])
      
      ax1.legend(loc=3)
      fig.savefig('psd_KE.'+format,format=format)
      
      
      ##
      ## Plot spectral flux
      ##
      
      fig = plt.figure()
      ax1 = fig.add_subplot(111)
      ax1.set_title('Spectral flux of KE in '+data_list[0][jr]['region'])
      
      for jd in range(0,nd):
         vk = data_list[jd][jr]['k']
         pi = data_list[jd][jr]['Pi_ke'].mean(axis=0)
         ax1.semilogx(vk,pi,label=data_list[jd][jr]['run'])
      
      ax1.legend(loc=4)
      fig.savefig('pi_KE.'+format,format=format)
      
      
      ##
      ## Plot wind forcing
      ##
      
      fig = plt.figure()
      ax1 = fig.add_subplot(111)
      ax1.set_title('Wind forcing in '+data_list[0][jr]['region'])
      
      for jd in range(0,nd):
         vk  = data_list[jd][jr]['k']
         print data_list[jd][jr].keys()
         if 'Pi_ke_tau' in data_list[jd][jr].keys():
            tau = data_list[jd][jr]['Pi_ke_tau'].mean(axis=0)
            ax1.semilogx(vk,tau,label=data_list[jd][jr]['run'])
      
      ax1.legend(loc=4)
      fig.savefig('pi_tau_ke.'+format,format=format)
      
      
      ##
      ## Plot horizontal viscosity
      ##
      
      fig = plt.figure()
      ax1 = fig.add_subplot(111)
      ax1.set_title('Horizontal viscosity in '+data_list[0][jr]['region'])
      
      for jd in range(0,nd):
         vk  = data_list[jd][jr]['k']
         print data_list[jd][jr].keys()
         if 'Pi_ke_visc' in data_list[jd][jr].keys():
            visc = data_list[jd][jr]['Pi_ke_visc'].mean(axis=0)
            ax1.semilogx(vk,visc,label=data_list[jd][jr]['run'])
      
      ax1.legend(loc=4)
      fig.savefig('pi_visc_ke.'+format,format=format)
      

def plot_stored_data_old(data_list):
   """
   """
   
   for jd in range(0,nd):
      for jr in range(0,nr):      
         data = data_list[jd][jr] 
         
         lat = data['tlat'][:]
         dep = data['deptht'][:]
         
         fig = plt.figure()
         ax1 = fig.add_subplot(111)
         rho  = np.ma.mean(data['cross_rho'],axis=0)
         uvel = np.ma.mean(data['cross_uvel'],axis=0)
         
         cf1 = ax1.contourf(lat, dep, uvel, levels=np.linspace(0,uvel.max(),20), extend='both' )
         ct1 = ax1.contour( lat, dep, rho, colors='k',levels=np.linspace(-0.3,0.3,6) )
         plt.clabel(ct1, fontsize=9, inline=1)
         plt.colorbar(cf1,ax=ax1)
         ax1.set_ylim([0,dep.max()])
         ax1.invert_yaxis()
         
         fig.savefig(pdir+'/cross_rho_'+data['name']+'_'+data['region']+'_'+starttime+'-'+endtime+'.pdf',format='pdf')
         
         
         if 1:
            zplot1 = np.mean(btke_store[:,:,:],axis=0)
            zplot2 = np.mean(bcke_store[:,:,:],axis=0)
            
            fig = plt.figure()
            ax1 = fig.add_subplot(121)
            ax1.set_title('Barotropic KE')
            ax2 = fig.add_subplot(122)
            ax2.set_title('Baroclinic KE (top layer)')
            levs = np.array([0,0.001,0.002,0.005,0.01,0.02,0.05,0.1,0.2])/4.
            levels_norm = np.linspace(0,1,levs.shape[0])
            ke_colors = cmap(levels_norm)
            vx = np.arange(0,xx2.shape[1])
            vy = np.arange(0,yy2.shape[0])
            xticks = np.arange(0,xx2.shape[1],np.floor(xx2.shape[1]/3))
            yticks = np.arange(0,yy2.shape[0],np.floor(yy2.shape[0]/3))
            xlabel = []
            ylabel = []
            for i in xticks:
               xlabel.append( '%6.2f' % (xx2[0,i],) )
            for j in yticks:
               ylabel.append( '%6.2f' % (yy2[j,0],) )
            cf1 = ax1.contourf(vx,vy,zplot1,extend='max',colors=ke_colors,levels=levs)
            cf2 = ax2.contourf(vx,vy,zplot2,extend='max',colors=ke_colors,levels=levs)
            ax1.set_xticks(xticks)
            ax1.set_xticklabels(xlabel)
            ax1.set_yticks(yticks)
            ax1.set_yticklabels(ylabel)
            ax2.set_xticks(xticks)
            ax2.set_xticklabels(xlabel)
            ax2.set_yticks(yticks)
            ax2.set_yticklabels(ylabel)
            cb2 = fig.colorbar(cf2,ax=ax2,spacing='uniform')
            fig.savefig('bt_bc_KE_xy_'+name[jd]+'_'+regions[jr]+'_'+starttime+'-'+endtime+'.pdf',format='pdf')
      
         ##
         ## Running mean
         ##
         if lrunningmean:
            run_data = np.zeros(data.shape)
            
            nmean = 12
            intmean = 73
            intn = run_mean_store.shape[2]/intmean
            int_mean_store = run_mean_store[:,:,0:intn,:].copy()
            
            varnames = data.keys()
            for jv in range(0,len(varnames)):
               tmp = data[varnames[jv]][:]
               print ' calc running mean ',varnames[jv],tmp.shape
               nsmooth = min(tmp.shape[0],nmean)
               
               run_tmp = convolve1d(tmp, np.ones((nsmooth,))/nsmooth, axis=0, mode='same')
               #tmp = np.convolve(tmp, np.ones((nsmooth,))/nsmooth, mode='same')
               run_data[:] = run_tmp[:]
               
         ##
         ## Calculate frequency spectrum
         ##

         if lpsd_freq:
            Lt = nt * 5. * 86400.
            if (np.mod(nt,2) == 0):
               arr = np.concatenate( (np.arange(0,nt/2+1),np.arange(-nt/2+1,0)) )
            else:
               nnt = (nt-1)/2
               arr = np.concatenate( (np.arange(0,nnt+1),np.arange(-nnt,0)) )
            wn_t = 2.0 * np.pi / Lt * arr
            wn_max = np.max( np.sqrt(wn_t**2) )
            dk = 2.0*np.pi/Lt
            k_t = np.arange(dk,wn_max+dk,dk)

            #use top layer
            uhat = fft(uvel_store[:,0,:,:],axis=0)
            vhat = fft(vvel_store[:,0,:,:],axis=0)
            #print uhat.shape
            ke_t = np.real(np.conj(uhat) * uhat) + np.real(np.conj(vhat) * vhat)

            nk_t = k_t.shape[0]
            psc = np.zeros((nk_t))
            ke1D_t = np.zeros((nk_t))
            for jk in range(0,nk_t):
               steps = np.where(wn_t**2 >= k_t[jk]**2)[0]
               psc[jk] = np.sum(ke_t[steps,:,:])
               steps = np.where(wn_t**2 == k_t[jk]**2)[0]
               ke1D_t[jk] = np.sum( ke_t[steps,:,:] ) / (nx*ny)
               #print steps

            psc = psc / (nx*ny)
            #print ke_t.shape,psc.shape
            psd_t = -(psc[1:] - psc[0:-1]) / dk
            psd_t[psd_t == 0] = 0.
            print 'psd_t min, max ',psd_t.min(),psd_t.max()
            print 'psd_t <= 0',psd_t[psd_t<=0]
            #ke1D_t = np.mean( np.mean(ke_t,axis=1), axis=1 )
            print 'ke_t <= 0',ke1D_t[ke1D_t<=0]
            for jk in range(0,nk_t-1):
               if (psc[jk] == psc[jk+1]):
                  print jk,psc[jk],ke1D_t[jk],ke1D_t[jk+1],psd_t[jk]

            ke_t  =  (psc[1:] + psc[0:-1]) * 0.5
            ke1D_t  =  (ke1D_t[1:] + ke1D_t[0:-1]) * 0.5
            vk_t = 0.5 * (k_t[1:] + k_t[0:-1])
            print ke1D_t.shape
            print wn_t.shape

            ##
            ## Also, plot surface plot of KE and std
            ##

            fig = plt.figure()
            ax1 = fig.add_subplot(2,1,1)
            ax2 = fig.add_subplot(2,1,2)
            ax1.set_title('Mean KE')
            ax2.set_title('Std KE')
            cf1 = ax1.pcolormesh(np.mean(uvel_store[:,0,:,:]**2 + vvel_store[:,0,:,:]**2,axis=0),vmin=-0.1,vmax=0.1)
            plt.colorbar(cf1,ax=ax1)
            cf2 = ax2.pcolormesh(np.std( uvel_store[:,0,:,:]**2 + vvel_store[:,0,:,:]**2,axis=0),vmin=0,vmax=0.1)
            plt.colorbar(cf2,ax=ax2)

            if (lpe and name[jd] != 'AVISO' and name[jd][0:10] != 'MITgcm-Rob'):
               fig = plt.figure()
               ax1 = fig.add_subplot(111)
               ax1.set_title('Conversion APE->KE',fontsize=12)
               zplot = np.mean(wb_store,axis=0)
               zmin = zplot.min()
               zmax = max(abs(zmin),zplot.max())
               levels = np.linspace(-zmax*0.8,zmax*0.8,11)
               cf1 = ax1.contourf(xx2,yy2,zplot,extend='both',cmap=plt.cm.RdBu_r)
               ax1.set_xlim([xx2.min(),xx2.max()])
               ax1.set_ylim([yy2.min(),yy2.max()])
               plt.colorbar(cf1)
               fig.savefig('wb_mean_'+regions[jr]+'_'+name[jd]+'.pdf',format='pdf')



def set_psd_fig(x='k'):
   fig = plt.figure()
   ax1  = fig.add_subplot(111)
   ax11 = ax1.twiny()

   if (x == 'k'):
      ax1.set_xlabel(r'Wavenumber $K$ $[\mathrm{m}^{-1}]$')
      ax11.set_xlabel(r'Wavelength $L = 2\pi/K$ $[\mathrm{km}]$')
   elif (x == 'o'):
      ax1.set_xlabel(r'Frequency $\omega$ $[\mathrm{days}^{-1}]$')
      ax11.set_xlabel(r'Time scale $[\mathrm{days}]$')

   return fig,ax1,ax11


def set_flux_fig():
   fig_1 = plt.figure(figsize=(8,10))
   ax_1  = fig_1.add_axes([0.2,0.7,0.4,0.18])
   ax_2  = fig_1.add_axes([0.2,0.5,0.4,0.18])
   ax_22 = fig_1.add_axes([0.2,0.3,0.4,0.18])
   ax_3  = fig_1.add_axes([0.2,0.1,0.4,0.18])
   ax_11 = ax_1.twiny()

   ax_1.spines['right'].set_visible(False)
   ax_1.spines['bottom'].set_visible(False)
   ax_1.yaxis.set_ticks_position('left')
   ax_1.tick_params(axis='x',which='both',bottom='off',top='off',labelbottom='off')

   ax_11.spines['right'].set_visible(False)
   ax_11.spines['bottom'].set_visible(False)
   ax_11.xaxis.set_ticks_position('top')
   ax_11.tick_params(axis='x',bottom='off',labelbottom='off')

   ax_2.spines['right'].set_visible(False)
   ax_2.spines['top'].set_visible(False)
   ax_2.spines['bottom'].set_visible(False)
   ax_2.yaxis.set_ticks_position('left')
   ax_2.tick_params(axis='x',which='both',bottom='off',top='off',labelbottom='off')
   #ax_22.spines['top'].set_visible(False)
   ax_22.spines['right'].set_visible(False)
   ax_22.yaxis.set_ticks_position('left')
   ax_22.spines['top'].set_visible(False)
   ax_22.spines['bottom'].set_visible(False)
   ax_22.tick_params(axis='x',which='both',bottom='off',top='off',labelbottom='off')
   #ax_22.xaxis.set_ticks_position('bottom')

   ax_3.spines['top'].set_visible(False)
   ax_3.spines['right'].set_visible(False)
   ax_3.yaxis.set_ticks_position('left')
   ax_3.xaxis.set_ticks_position('bottom')

   ax_11.set_xlabel(r'Wavelength $2\pi/K$ $[\mathrm{km}]$')
   ax_3.set_xlabel(r'Wavenumber $K=\sqrt{k^2 + l^2}$ $[\mathrm{m}^{-1}]$')

   return fig_1,ax_1,ax_11,ax_2,ax_22,ax_3


def k2l(k):
    """
    Returns length given wavenumber k.
    """
    return (2.0 * np.pi / k)/1000. ## [km]

def k2t(k):
    """
    Returns time scale given wavenumber k.
    """
    return 2.0 * np.pi / k / 86400.


def convert_k_to_len(ax_f):
    """
    Update second axis according with first axis.
    """
    y1, y2 = ax1.get_xlim()
    ax11.set_xlim(k2l(y1), k2l(y2))
    ax11.figure.canvas.draw()

def convert_k_to_len2(ax_f2):
    """
    Update second axis according with first axis.
    """
    y1, y2 = ax1.get_xlim()
    ax11.set_xlim(k2l(y1), k2l(y2))
    ax11.figure.canvas.draw()

    y1, y2 = ax2.get_xlim()
    ax22.set_xlim(k2l(y1), k2l(y2))
    ax22.figure.canvas.draw()

def convert_k_to_len4(ax_f4):
    """
    Update second axis according with first axis.
    """
    y1, y2 = ax11.get_xlim()
    ax_11.set_xlim(k2l(y1), k2l(y2))
    ax_11.figure.canvas.draw()

    y1, y2 = ax12.get_xlim()
    ax_12.set_xlim(k2l(y1), k2l(y2))
    ax_12.figure.canvas.draw()

    y1, y2 = ax21.get_xlim()
    ax_21.set_xlim(k2l(y1), k2l(y2))
    ax_21.figure.canvas.draw()

    y1, y2 = ax22.get_xlim()
    ax_22.set_xlim(k2l(y1), k2l(y2))
    ax_22.figure.canvas.draw()

def convert_k_to_time(ax_f):
    y1, y2 = ax1.get_xlim()
    ax11.set_xlim(k2t(y1), k2t(y2))
    ax11.figure.canvas.draw()



def draw_lines(ax1,vk,psd,m0=False,m53=False,m3=False,m23=False):
   """
   """

   ##
   ## Plot some example slopes
   ##
   if (1):
      if (m0):
         ## draw zero line
         ax1.semilogx(vk,vk*0,'-k')

      if (m53):
         ## draw -5/3 slope
         psdmax_m53 = vk.max()
         psdmin_m53 = psdmax_m53/50.
         kmax_m53 = vk[psd==psdmax_m53] * 1.5
         C_m53 = psdmax_m53/kmax_m53**(-5./3.)
         kmin_m53 = (psdmin_m53/C_m53)**(-3./5.)
         ax1.loglog([kmin_m53,kmax_m53],[psdmin_m53,psdmax_m53],\
                     '--k',linewidth=2.,label='k^(-5/3)')

      if (m3):
         psdmax_m3 = psd.max()
         psdmin_m3 = psdmax_m3/50.
         kmax_m3 = vk[psd==psdmax_m3] * 2.
         C_m3 = psdmax_m3/kmax_m3**(-3.)
         kmin_m3 = (psdmin_m3/C_m3)**(-1./3.)
         ax1.loglog([kmin_m3,kmax_m3],[psdmin_m3,psdmax_m3],\
                    ':k',linewidth=2.,label='k^(-3)')

               
               
               