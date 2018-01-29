##
##
##  Equation of state calculations taken from NEMO
##
##  Joakim Kjellsson, AOPP, 2016
##
##

import numpy as np
   
def eos_insitu(tem,sal,dep,tmask):
   """
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE eos_insitu  ***
      !!
      !! ** Purpose :   Compute the in situ density (ratio rho/rau0) from
      !!       potential temperature and salinity using an equation of state
      !!       defined through the namelist parameter nn_eos.
      !!
      !! ** Method  :   prd(t,s,z) = ( rho(t,s,z) - rau0 ) / rau0
      !!         with   prd    in situ density anomaly      no units
      !!                t      TEOS10: CT or EOS80: PT      Celsius
      !!                s      TEOS10: SA or EOS80: SP      TEOS10: g/kg or EOS80: psu
      !!                z      depth                        meters
      !!                rho    in situ density              kg/m^3
      !!                rau0   reference density            kg/m^3
      !!
      !!     nn_eos = -1 : polynomial TEOS-10 equation of state is used for rho(t,s,z).
      !!         Check value: rho = 1028.21993233072 kg/m^3 for z=3000 dbar, ct=3 Celcius, sa=35.5 g/kg
      !!
      !!     nn_eos =  0 : polynomial EOS-80 equation of state is used for rho(t,s,z).
      !!         Check value: rho = 1028.35011066567 kg/m^3 for z=3000 dbar, pt=3 Celcius, sp=35.5 psu
      !!
      !!     nn_eos =  1 : simplified equation of state
      !!              prd(t,s,z) = ( -a0*(1+lambda/2*(T-T0)+mu*z+nu*(S-S0))*(T-T0) + b0*(S-S0) ) / rau0
      !!              linear case function of T only: rn_alpha<>0, other coefficients = 0
      !!              linear eos function of T and S: rn_alpha and rn_beta<>0, other coefficients=0
      !!              Vallis like equation: use default values of coefficients
      !!
      !! ** Action  :   compute prd , the in situ density (no units)
      !!
      !! References :   Roquet et al, Ocean Modelling, in preparation (2014)
      !!                Vallis, Atmospheric and Oceanic Fluid Dynamics, 2006
      !!                TEOS-10 Manual, 2010
      !!----------------------------------------------------------------------
   """
   
   jpi = tem.shape[2]
   jpj = tem.shape[1]
   jpk = tem.shape[0]
   jpkm1 = jpk-1
   
   eos = eos_vars()
   
   rhd = np.zeros((jpk,jpj,jpi))
   
   if( eos.nn_eos == -1 or eos.nn_eos == 0 ):
      ##          !==  polynomial TEOS-10 / EOS-80 ==!
      zh  = dep[:,:,:] * eos.r1_Z0                                  # depth
      zt  = tem[:,:,:] * eos.r1_T0                                  # temperature
      zs  = np.sqrt( np.abs( sal[:,:,:] + eos.rdeltaS ) * eos.r1_S0 )   # square root salinity
      ztm = tmask[:,:,:]                                        # tmask
      
      
      zn3 = eos.EOS013*zt \
          + eos.EOS103*zs+eos.EOS003
      
      zn2 = (eos.EOS022*zt   \
          + eos.EOS112*zs+eos.EOS012)*zt   \
          + (eos.EOS202*zs+eos.EOS102)*zs+eos.EOS002
      
      zn1 = (((eos.EOS041*zt   \
          + eos.EOS131*zs+eos.EOS031)*zt   \
          + (eos.EOS221*zs+eos.EOS121)*zs+eos.EOS021)*zt   \
          + ((eos.EOS311*zs+eos.EOS211)*zs+eos.EOS111)*zs+eos.EOS011)*zt   \
          + (((eos.EOS401*zs+eos.EOS301)*zs+eos.EOS201)*zs+eos.EOS101)*zs+eos.EOS001
      
      zn0 = (((((eos.EOS060*zt   \
          + eos.EOS150*zs+eos.EOS050)*zt   \
          + (eos.EOS240*zs+eos.EOS140)*zs+eos.EOS040)*zt   \
          + ((eos.EOS330*zs+eos.EOS230)*zs+eos.EOS130)*zs+eos.EOS030)*zt   \
          + (((eos.EOS420*zs+eos.EOS320)*zs+eos.EOS220)*zs+eos.EOS120)*zs+eos.EOS020)*zt   \
          + ((((eos.EOS510*zs+eos.EOS410)*zs+eos.EOS310)*zs+eos.EOS210)*zs+eos.EOS110)*zs+eos.EOS010)*zt   \
          + (((((eos.EOS600*zs+eos.EOS500)*zs+eos.EOS400)*zs+eos.EOS300)*zs+eos.EOS200)*zs+eos.EOS100)*zs+eos.EOS000
      
      zn  = ( ( zn3 * zh + zn2 ) * zh + zn1 ) * zh + zn0
      
      rhd[:,:,:] = (  zn * eos.r1_rau0 - 1.   ) * ztm  # ! density anomaly (masked)
      
   if ( eos.nn_eos == 1 ):               # !==  simplified EOS  ==!
      zt  = tem[:,:,:] - 10. 
      zs  = sal[:,:,:] - 35. 
      zh  = dep[:,:,:]
      ztm = tmask[:,:,:]
         
      zn =  - eos.rn_a0 * ( 1.  + 0.5 * eos.rn_lambda1*zt + eos.rn_mu1*zh ) * zt   \
         + eos.rn_b0 * ( 1.  - 0.5 *eos.rn_lambda2*zs - eos.rn_mu2*zh ) * zs   \
         - eos.rn_nu * zt * zs
                                            
      rhd[:,:,:] = zn * eos.r1_rau0 * ztm              #  ! density anomaly (masked)
   
   return rhd   
   
   
def eos_insitu_pot(tem,sal,dep,tmask):
   """
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE eos_insitu_pot  ***
      !!
      !! ** Purpose :   Compute the in situ density (ratio rho/rau0) and the
      !!      potential volumic mass (Kg/m3) from potential temperature and
      !!      salinity fields using an equation of state defined through the
      !!     namelist parameter nn_eos.
      !!
      !! ** Action  : - prd  , the in situ density (no units)
      !!              - prhop, the potential volumic mass (Kg/m3)
      !!
      !!----------------------------------------------------------------------
      INTEGER  ::   ji, jj, jk, jsmp             ! dummy loop indices
      INTEGER  ::   jdof
      REAL(wp) ::   zt , zh , zstemp, zs , ztm   ! local scalars
      REAL(wp) ::   zn , zn0, zn1, zn2, zn3      !   -      -
      REAL(wp), DIMENSION(:), ALLOCATABLE :: zn0_sto, zn_sto, zsign    ! local vectors
      !!----------------------------------------------------------------------
   """
   
   jpi = tem.shape[2]
   jpj = tem.shape[1]
   jpk = tem.shape[0]
   jpkm1 = jpk
   
   rho  = np.zeros((jpk,jpj,jpi))
   rhop = np.zeros((jpk,jpj,jpi))
   
   eos = eos_vars()
   
   if (eos.nn_eos == -1 or eos.nn_eos == 0): #               !==  polynomial TEOS-10 / EOS-80 ==!
      zh  = dep[:,:,:] * eos.r1_Z0                       #           ! depth
      zt  = tem[:,:,:] * eos.r1_T0                       #    ! temperature
      zs  = np.sqrt( np.abs( sal[:,:,:] + eos.rdeltaS ) * eos.r1_S0 ) #  ! square root salinity
      ztm = tmask[:,:,:]                                    #     ! tmask
      
      zn3 = eos.EOS013*zt   \
          + eos.EOS103*zs+eos.EOS003
      
      zn2 = (eos.EOS022*zt   \
          + eos.EOS112*zs+eos.EOS012)*zt   \
          + (eos.EOS202*zs+eos.EOS102)*zs+eos.EOS002
       
      zn1 = (((eos.EOS041*zt   \
          + eos.EOS131*zs+eos.EOS031)*zt   \
          + (eos.EOS221*zs+eos.EOS121)*zs+eos.EOS021)*zt   \
          + ((eos.EOS311*zs+eos.EOS211)*zs+eos.EOS111)*zs+eos.EOS011)*zt   \
          + (((eos.EOS401*zs+eos.EOS301)*zs+eos.EOS201)*zs+eos.EOS101)*zs+eos.EOS001
        
      zn0 = (((((eos.EOS060*zt   \
          + eos.EOS150*zs+eos.EOS050)*zt   \
          + (eos.EOS240*zs+eos.EOS140)*zs+eos.EOS040)*zt   \
          + ((eos.EOS330*zs+eos.EOS230)*zs+eos.EOS130)*zs+eos.EOS030)*zt   \
          + (((eos.EOS420*zs+eos.EOS320)*zs+eos.EOS220)*zs+eos.EOS120)*zs+eos.EOS020)*zt   \
          + ((((eos.EOS510*zs+eos.EOS410)*zs+eos.EOS310)*zs+eos.EOS210)*zs+eos.EOS110)*zs+eos.EOS010)*zt   \
          + (((((eos.EOS600*zs+eos.EOS500)*zs+eos.EOS400)*zs+eos.EOS300)*zs+eos.EOS200)*zs+eos.EOS100)*zs+eos.EOS000
      
      zn  = ( ( zn3 * zh + zn2 ) * zh + zn1 ) * zh + zn0
      
      rhop[:,:,:] = zn0 * ztm                          # ! potential density referenced at the surface
      
      rho[:,:,:] = (  zn * eos.r1_rau0 - 1.   ) * ztm    #  ! density anomaly (masked)
         
   if (nn_eos == 1):                # !==  simplified EOS  ==!
      zt  = tem[:,:,:] - 10. 
      zs  = sal[:,:,:] - 35. 
      zh  = dep[:,:,:]
      ztm = tmask[:,:,:]
      
      #                                                     ! potential density referenced at the surface
      zn =  - eos.rn_a0 * ( 1.  + 0.5 *eos.rn_lambda1*zt ) * zt   \
         + eos.rn_b0 * ( 1.  - 0.5 *eos.rn_lambda2*zs ) * zs   \
         - eos.rn_nu * zt * zs
      rhop[:,:,:] = ( eos.rau0 + zn ) * ztm
      #                                                     ! density anomaly (masked)
      zn = zn - ( eos.rn_a0 * eos.rn_mu1 * zt + eos.rn_b0 * eos.rn_mu2 * zs ) * zh
      rho[:,:,:] = zn * eos.r1_rau0 * ztm
   
   return rhop
         
   
def rab_3d(tem,sal,dep,tmask):
   """
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE rab_3d  ***
      !!
      !! ** Purpose :   Calculates thermal/haline expansion ratio at T-points
      !!
      !! ** Method  :   calculates alpha / beta at T-points
      !!
      !! ** Action  : - pab     : thermal/haline expansion ratio at T-points
      !!----------------------------------------------------------------------
      INTEGER  ::   ji, jj, jk                ! dummy loop indices
      REAL(wp) ::   zt , zh , zs , ztm        ! local scalars
      REAL(wp) ::   zn , zn0, zn1, zn2, zn3   !   -      -
      !!----------------------------------------------------------------------
   """
   
   jpi = tem.shape[2]
   jpj = tem.shape[1]
   jpk = tem.shape[0]
   jpkm1 = jpk-1
   
   alpha = np.zeros((jpk,jpj,jpi))
   beta  = np.zeros((jpk,jpj,jpi))
   
   eos = eos_vars()
   
   if (eos.nn_eos == -1 or eos.nn_eos == 0): #                !==  polynomial TEOS-10 / EOS-80 ==!
      zh  = dep[:,:,:] * eos.r1_Z0                               # ! depth
      zt  = tem[:,:,:] * eos.r1_T0                           # ! temperature
      zs  = np.sqrt( np.abs( sal[:,:,:] + eos.rdeltaS ) * eos.r1_S0 ) #  ! square root salinity
      ztm = tmask[:,:,:]                                #         ! tmask
      
      # alpha
      zn3 = eos.ALP003
      zn2 = eos.ALP012*zt + eos.ALP102*zs+eos.ALP002
      zn1 = ((eos.ALP031*zt   \
          + eos.ALP121*zs+eos.ALP021)*zt   \
          + (eos.ALP211*zs+eos.ALP111)*zs+eos.ALP011)*zt   \
          + ((eos.ALP301*zs+eos.ALP201)*zs+eos.ALP101)*zs+eos.ALP001
      
      zn0 = ((((eos.ALP050*zt   \
          + eos.ALP140*zs+eos.ALP040)*zt   \
          + (eos.ALP230*zs+eos.ALP130)*zs+eos.ALP030)*zt   \
          + ((eos.ALP320*zs+eos.ALP220)*zs+eos.ALP120)*zs+eos.ALP020)*zt   \
          + (((eos.ALP410*zs+eos.ALP310)*zs+eos.ALP210)*zs+eos.ALP110)*zs+eos.ALP010)*zt   \
          + ((((eos.ALP500*zs+eos.ALP400)*zs+eos.ALP300)*zs+eos.ALP200)*zs+eos.ALP100)*zs+eos.ALP000
        
      zn  = ( ( zn3 * zh + zn2 ) * zh + zn1 ) * zh + zn0
      
      alpha[:,:,:] = zn * eos.r1_rau0 * ztm
      
      # beta
      zn3 = eos.BET003
      zn2 = eos.BET012*zt + eos.BET102*zs+eos.BET002
      zn1 = ((eos.BET031*zt   \
          + eos.BET121*zs+eos.BET021)*zt   \
          + (eos.BET211*zs+eos.BET111)*zs+eos.BET011)*zt   \
          + ((eos.BET301*zs+eos.BET201)*zs+eos.BET101)*zs+eos.BET001
      
      zn0 = ((((eos.BET050*zt   \
          + eos.BET140*zs+eos.BET040)*zt   \
          + (eos.BET230*zs+eos.BET130)*zs+eos.BET030)*zt   \
          + ((eos.BET320*zs+eos.BET220)*zs+eos.BET120)*zs+eos.BET020)*zt   \
          + (((eos.BET410*zs+eos.BET310)*zs+eos.BET210)*zs+eos.BET110)*zs+eos.BET010)*zt   \
          + ((((eos.BET500*zs+eos.BET400)*zs+eos.BET300)*zs+eos.BET200)*zs+eos.BET100)*zs+eos.BET000
      
      zn  = ( ( zn3 * zh + zn2 ) * zh + zn1 ) * zh + zn0
      
      beta[:,:,:] = zn / zs * eos.r1_rau0 * ztm
      
   if (eos.nn_eos == 1): #                  !==  simplified EOS  ==!
      zt  = tem[:,:,:] - 10.   # ! pot. temperature anomaly (t-T0)
      zs  = sal[:,:,:] - 35.   # ! abs. salinity anomaly (s-S0)
      zh  = dep[:,:,:]            #     ! depth in meters at t-point
      ztm = tmask[:,:,:]          #        ! land/sea bottom mask = surf. mask
      
      zn  = eos.rn_a0 * ( 1.  + eos.rn_lambda1*zt + eos.rn_mu1*zh ) + eos.rn_nu*zs
      alpha[:,:,:] = zn * eos.r1_rau0 * ztm  # ! alpha
         
      zn  = eos.rn_b0 * ( 1.  - eos.rn_lambda2*zs - eos.rn_mu2*zh ) - eos.rn_nu*zt
      beta[:,:,:] = zn * eos.r1_rau0 * ztm  # ! beta
   
   return alpha,beta
         
         
def bn2(tem,sal,alpha,beta,dep,depw,dzw,tmask):
   """
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE bn2  ***
      !!
      !! ** Purpose :   Compute the local Brunt-Vaisala frequency at the 
      !!                time-step of the input arguments
      !!
      !! ** Method  :   pn2 = grav * (alpha dk[T] + beta dk[S] ) / e3w
      !!      where alpha and beta are given in pab, and computed on T-points.
      !!      N.B. N^2 is set one for all to zero at jk=1 in istate module.
      !!
      !! ** Action  :   pn2 : square of the brunt-vaisala frequency at w-point 
      !!
      !!----------------------------------------------------------------------
   """
   
   jpi = tem.shape[2]
   jpj = tem.shape[1]
   jpk = tem.shape[0]
   jpkm1 = jpk-1
   
   N2 = np.zeros((jpk,jpj,jpi))
   
   zrw =   ( depw[1:-1,:,:] - dep[1:-1,:,:] ) / ( dep[0:-2,:,:] - dep[1:-1,:,:] ) 
   zaw = alpha[1:-1,:,:] * (1. - zrw) + alpha[0:-2,:,:] * zrw 
   zbw = beta[1:-1,:,:]  * (1. - zrw) + beta[0:-2,:,:]  * zrw
   
   # interior points only (2=< jk =< jpkm1 )
   # surface and bottom value set to zero one for all in istate.F90
   N2[1:-1,:,:] = 9.81 * (   zaw * ( tem[0:-2,:,:] - tem[1:-1,:,:] )     \
                           - zbw * ( sal[0:-2,:,:] - sal[1:-1,:,:] )  )  \
                         / dzw[1:-1,:,:] * tmask[1:-1,:,:]
   
   return N2
   
   
def eos_pen(tem,sal,dep,tmask):
   """
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE eos_pen  ***
      !!
      !! ** Purpose :   Calculates nonlinear anomalies of alpha_PE, beta_PE and PE at T-points
      !!
      !! ** Method  :   PE is defined analytically as the vertical 
      !!                   primitive of EOS times -g integrated between 0 and z>0.
      !!                pen is the nonlinear bsq-PE anomaly: pen = ( PE - rau0 gz ) / rau0 gz - rd
      !!                                                      = 1/z * /int_0^z rd dz - rd 
      !!                                where rd is the density anomaly (see eos_rhd function)
      !!                ab_pe are partial derivatives of PE anomaly with respect to T and S:
      !!                    ab_pe(1) = - 1/(rau0 gz) * dPE/dT + drd/dT = - d(pen)/dT
      !!                    ab_pe(2) =   1/(rau0 gz) * dPE/dS + drd/dS =   d(pen)/dS
      !!
      !! ** Action  : - pen         : PE anomaly given at T-points
      !!            : - pab_pe  : given at T-points
      !!                    pab_pe(:,:,:,jp_tem) is alpha_pe
      !!                    pab_pe(:,:,:,jp_sal) is beta_pe
      !!----------------------------------------------------------------------
      INTEGER  ::   ji, jj, jk                ! dummy loop indices
      REAL(wp) ::   zt , zh , zs , ztm        ! local scalars
      REAL(wp) ::   zn , zn0, zn1, zn2        !   -      -
      !!----------------------------------------------------------------------
   """
   
   jpi = tem.shape[2]
   jpj = tem.shape[1]
   jpk = tem.shape[0]
   jpkm1 = jpk-1
   
   eos = eos_vars()
   
   PE_an = np.zeros((jpk,jpj,jpi))
   dPE_dT = np.zeros((jpk,jpj,jpi))
   dPE_dS = np.zeros((jpk,jpj,jpi))
   
   if (eos.nn_eos == -1 or eos.nn_eos == 0): #                !==  polynomial TEOS-10 / EOS-80 ==!
      zh  = dep[:,:,:] * eos.r1_Z0                        #        ! depth
      zt  = tem[:,:,:] * eos.r1_T0                        #   ! temperature
      zs  = np.sqrt( np.abs( sal[:,:,:] + eos.rdeltaS ) * eos.r1_S0 ) #  ! square root salinity
      ztm = tmask[:,:,:]                        #                 ! tmask
      
      zn2 = (eos.PEN012)*zt   \
          + eos.PEN102*zs+eos.PEN002
      
      zn1 = ((eos.PEN021)*zt   \
          + eos.PEN111*zs+eos.PEN011)*zt   \
          + (eos.PEN201*zs+eos.PEN101)*zs+eos.PEN001
      
      zn0 = ((((eos.PEN040)*zt   \
          + eos.PEN130*zs+eos.PEN030)*zt   \
          + (eos.PEN220*zs+eos.PEN120)*zs+eos.PEN020)*zt   \
          + ((eos.PEN310*zs+eos.PEN210)*zs+eos.PEN110)*zs+eos.PEN010)*zt   \
          + (((eos.PEN400*zs+eos.PEN300)*zs+eos.PEN200)*zs+eos.PEN100)*zs+eos.PEN000
         
      zn  = ( zn2 * zh + zn1 ) * zh + zn0
      
      PE_an[:,:,:]  = zn * zh * eos.r1_rau0 * ztm
      
      #
      # alphaPE non-linear anomaly
      zn2 = eos.APE002
      
      zn1 = (eos.APE011)*zt   \
          + eos.APE101*zs+eos.APE001
      
      zn0 = (((eos.APE030)*zt   \
          + eos.APE120*zs+eos.APE020)*zt   \
          + (eos.APE210*zs+eos.APE110)*zs+eos.APE010)*zt   \
          + ((eos.APE300*zs+eos.APE200)*zs+eos.APE100)*zs+eos.APE000
      
      zn  = ( zn2 * zh + zn1 ) * zh + zn0
                                    
      dPE_dT[:,:,:] = zn * zh * eos.r1_rau0 * ztm
      
      #
      # betaPE non-linear anomaly
      zn2 = eos.BPE002
      
      zn1 = (eos.BPE011)*zt   \
          + eos.BPE101*zs+eos.BPE001
      
      zn0 = (((eos.BPE030)*zt   \
          + eos.BPE120*zs+eos.BPE020)*zt   \
          + (eos.BPE210*zs+eos.BPE110)*zs+eos.BPE010)*zt   \
          + ((eos.BPE300*zs+eos.BPE200)*zs+eos.BPE100)*zs+eos.BPE000
      
      zn  = ( zn2 * zh + zn1 ) * zh + zn0
      
      dPE_dS[:,:,:] = zn / zs * zh * eos.r1_rau0 * ztm
      
   if (nn_eos == 1): #                !==  Vallis (2006) simplified EOS  ==!
      zt  = tem[:,:,:] - 10.  # ! temperature anomaly (t-T0)
      zs  = sal[:,:,:] - 35.  # ! abs. salinity anomaly (s-S0)
      zh  = dep[:,:,:]           #    ! depth in meters  at t-point
      ztm = tmask[:,:,:]         #       ! tmask
      zn  = 0.5 * zh * eos.r1_rau0 * ztm
      #                                    ! Potential Energy
      PE_an[:,:,:] = ( eos.rn_a0 * eos.rn_mu1 * zt + eos.rn_b0 * eos.rn_mu2 * zs ) * zn
      #                                    ! alphaPE
      dPE_dT[:,:,:] = - eos.rn_a0 * eos.rn_mu1 * zn
      dPE_dS[:,:,:] =   eos.rn_b0 * eos.rn_mu2 * zn
   
   return PE_an,dPE_dT,dPE_dS
         
         
class eos_vars:
   
   #def __init__(self,nn_eos):
   #   self.nn_eos = nn_eos
   
   if (1):   
      rau0        = 1026.               #  !: volumic mass of reference     [kg/m3]
      rcp         = 3991.86795711963    #  !: heat capacity     [J/K]
      nn_eos      = -1
      
      if (nn_eos == -1): #                       !==  polynomial TEOS-10  ==!
         #
         rdeltaS = 32. 
         r1_S0  = 0.875 /35.16504 
         r1_T0  = 1. /40. 
         r1_Z0  = 1.e-4 
         #
         EOS000 = 8.0189615746e+02 
         EOS100 = 8.6672408165e+02 
         EOS200 = -1.7864682637e+03 
         EOS300 = 2.0375295546e+03 
         EOS400 = -1.2849161071e+03 
         EOS500 = 4.3227585684e+02 
         EOS600 = -6.0579916612e+01 
         EOS010 = 2.6010145068e+01 
         EOS110 = -6.5281885265e+01 
         EOS210 = 8.1770425108e+01 
         EOS310 = -5.6888046321e+01 
         EOS410 = 1.7681814114e+01 
         EOS510 = -1.9193502195 
         EOS020 = -3.7074170417e+01 
         EOS120 = 6.1548258127e+01 
         EOS220 = -6.0362551501e+01 
         EOS320 = 2.9130021253e+01 
         EOS420 = -5.4723692739 
         EOS030 = 2.1661789529e+01 
         EOS130 = -3.3449108469e+01 
         EOS230 = 1.9717078466e+01 
         EOS330 = -3.1742946532 
         EOS040 = -8.3627885467 
         EOS140 = 1.1311538584e+01 
         EOS240 = -5.3563304045 
         EOS050 = 5.4048723791e-01 
         EOS150 = 4.8169980163e-01 
         EOS060 = -1.9083568888e-01 
         EOS001 = 1.9681925209e+01 
         EOS101 = -4.2549998214e+01 
         EOS201 = 5.0774768218e+01 
         EOS301 = -3.0938076334e+01 
         EOS401 = 6.6051753097 
         EOS011 = -1.3336301113e+01 
         EOS111 = -4.4870114575 
         EOS211 = 5.0042598061 
         EOS311 = -6.5399043664e-01 
         EOS021 = 6.7080479603 
         EOS121 = 3.5063081279 
         EOS221 = -1.8795372996 
         EOS031 = -2.4649669534 
         EOS131 = -5.5077101279e-01 
         EOS041 = 5.5927935970e-01 
         EOS002 = 2.0660924175 
         EOS102 = -4.9527603989 
         EOS202 = 2.5019633244 
         EOS012 = 2.0564311499 
         EOS112 = -2.1311365518e-01 
         EOS022 = -1.2419983026 
         EOS003 = -2.3342758797e-02 
         EOS103 = -1.8507636718e-02 
         EOS013 = 3.7969820455e-01 
         #
         ALP000 = -6.5025362670e-01 
         ALP100 = 1.6320471316 
         ALP200 = -2.0442606277 
         ALP300 = 1.4222011580 
         ALP400 = -4.4204535284e-01 
         ALP500 = 4.7983755487e-02 
         ALP010 = 1.8537085209 
         ALP110 = -3.0774129064 
         ALP210 = 3.0181275751 
         ALP310 = -1.4565010626 
         ALP410 = 2.7361846370e-01 
         ALP020 = -1.6246342147 
         ALP120 = 2.5086831352 
         ALP220 = -1.4787808849 
         ALP320 = 2.3807209899e-01 
         ALP030 = 8.3627885467e-01 
         ALP130 = -1.1311538584 
         ALP230 = 5.3563304045e-01 
         ALP040 = -6.7560904739e-02 
         ALP140 = -6.0212475204e-02 
         ALP050 = 2.8625353333e-02 
         ALP001 = 3.3340752782e-01 
         ALP101 = 1.1217528644e-01 
         ALP201 = -1.2510649515e-01 
         ALP301 = 1.6349760916e-02 
         ALP011 = -3.3540239802e-01 
         ALP111 = -1.7531540640e-01 
         ALP211 = 9.3976864981e-02 
         ALP021 = 1.8487252150e-01 
         ALP121 = 4.1307825959e-02 
         ALP031 = -5.5927935970e-02 
         ALP002 = -5.1410778748e-02 
         ALP102 = 5.3278413794e-03 
         ALP012 = 6.2099915132e-02 
         ALP003 = -9.4924551138e-03 
         #
         BET000 = 1.0783203594e+01 
         BET100 = -4.4452095908e+01 
         BET200 = 7.6048755820e+01 
         BET300 = -6.3944280668e+01 
         BET400 = 2.6890441098e+01 
         BET500 = -4.5221697773 
         BET010 = -8.1219372432e-01 
         BET110 = 2.0346663041 
         BET210 = -2.1232895170 
         BET310 = 8.7994140485e-01 
         BET410 = -1.1939638360e-01 
         BET020 = 7.6574242289e-01 
         BET120 = -1.5019813020 
         BET220 = 1.0872489522 
         BET320 = -2.7233429080e-01 
         BET030 = -4.1615152308e-01 
         BET130 = 4.9061350869e-01 
         BET230 = -1.1847737788e-01 
         BET040 = 1.4073062708e-01 
         BET140 = -1.3327978879e-01 
         BET050 = 5.9929880134e-03 
         BET001 = -5.2937873009e-01 
         BET101 = 1.2634116779 
         BET201 = -1.1547328025 
         BET301 = 3.2870876279e-01 
         BET011 = -5.5824407214e-02 
         BET111 = 1.2451933313e-01 
         BET211 = -2.4409539932e-02 
         BET021 = 4.3623149752e-02 
         BET121 = -4.6767901790e-02 
         BET031 = -6.8523260060e-03 
         BET002 = -6.1618945251e-02 
         BET102 = 6.2255521644e-02 
         BET012 = -2.6514181169e-03 
         BET003 = -2.3025968587e-04 
         #
         PEN000 = -9.8409626043 
         PEN100 = 2.1274999107e+01 
         PEN200 = -2.5387384109e+01 
         PEN300 = 1.5469038167e+01 
         PEN400 = -3.3025876549 
         PEN010 = 6.6681505563 
         PEN110 = 2.2435057288 
         PEN210 = -2.5021299030 
         PEN310 = 3.2699521832e-01 
         PEN020 = -3.3540239802 
         PEN120 = -1.7531540640 
         PEN220 = 9.3976864981e-01 
         PEN030 = 1.2324834767 
         PEN130 = 2.7538550639e-01 
         PEN040 = -2.7963967985e-01 
         PEN001 = -1.3773949450 
         PEN101 = 3.3018402659 
         PEN201 = -1.6679755496 
         PEN011 = -1.3709540999 
         PEN111 = 1.4207577012e-01 
         PEN021 = 8.2799886843e-01 
         PEN002 = 1.7507069098e-02 
         PEN102 = 1.3880727538e-02 
         PEN012 = -2.8477365341e-01 
         #
         APE000 = -1.6670376391e-01 
         APE100 = -5.6087643219e-02 
         APE200 = 6.2553247576e-02 
         APE300 = -8.1748804580e-03 
         APE010 = 1.6770119901e-01 
         APE110 = 8.7657703198e-02 
         APE210 = -4.6988432490e-02 
         APE020 = -9.2436260751e-02 
         APE120 = -2.0653912979e-02 
         APE030 = 2.7963967985e-02 
         APE001 = 3.4273852498e-02 
         APE101 = -3.5518942529e-03 
         APE011 = -4.1399943421e-02 
         APE002 = 7.1193413354e-03 
         #
         BPE000 = 2.6468936504e-01 
         BPE100 = -6.3170583896e-01 
         BPE200 = 5.7736640125e-01 
         BPE300 = -1.6435438140e-01 
         BPE010 = 2.7912203607e-02 
         BPE110 = -6.2259666565e-02 
         BPE210 = 1.2204769966e-02 
         BPE020 = -2.1811574876e-02 
         BPE120 = 2.3383950895e-02 
         BPE030 = 3.4261630030e-03 
         BPE001 = 4.1079296834e-02 
         BPE101 = -4.1503681096e-02 
         BPE011 = 1.7676120780e-03 
         BPE002 = 1.7269476440e-04 
         #
      if (nn_eos == 0): #                        !==  polynomial EOS-80 formulation  ==!
         #
         rdeltaS = 20. 
         r1_S0  = 1. /40. 
         r1_T0  = 1. /40. 
         r1_Z0  = 1.e-4 
         #
         EOS000 = 9.5356891948e+02 
         EOS100 = 1.7136499189e+02 
         EOS200 = -3.7501039454e+02 
         EOS300 = 5.1856810420e+02 
         EOS400 = -3.7264470465e+02 
         EOS500 = 1.4302533998e+02 
         EOS600 = -2.2856621162e+01 
         EOS010 = 1.0087518651e+01 
         EOS110 = -1.3647741861e+01 
         EOS210 = 8.8478359933 
         EOS310 = -7.2329388377 
         EOS410 = 1.4774410611 
         EOS510 = 2.0036720553e-01 
         EOS020 = -2.5579830599e+01 
         EOS120 = 2.4043512327e+01 
         EOS220 = -1.6807503990e+01 
         EOS320 = 8.3811577084 
         EOS420 = -1.9771060192 
         EOS030 = 1.6846451198e+01 
         EOS130 = -2.1482926901e+01 
         EOS230 = 1.0108954054e+01 
         EOS330 = -6.2675951440e-01 
         EOS040 = -8.0812310102 
         EOS140 = 1.0102374985e+01 
         EOS240 = -4.8340368631 
         EOS050 = 1.2079167803 
         EOS150 = 1.1515380987e-01 
         EOS060 = -2.4520288837e-01 
         EOS001 = 1.0748601068e+01 
         EOS101 = -1.7817043500e+01 
         EOS201 = 2.2181366768e+01 
         EOS301 = -1.6750916338e+01 
         EOS401 = 4.1202230403 
         EOS011 = -1.5852644587e+01 
         EOS111 = -7.6639383522e-01 
         EOS211 = 4.1144627302 
         EOS311 = -6.6955877448e-01 
         EOS021 = 9.9994861860 
         EOS121 = -1.9467067787e-01 
         EOS221 = -1.2177554330 
         EOS031 = -3.4866102017 
         EOS131 = 2.2229155620e-01 
         EOS041 = 5.9503008642e-01 
         EOS002 = 1.0375676547 
         EOS102 = -3.4249470629 
         EOS202 = 2.0542026429 
         EOS012 = 2.1836324814 
         EOS112 = -3.4453674320e-01 
         EOS022 = -1.2548163097 
         EOS003 = 1.8729078427e-02 
         EOS103 = -5.7238495240e-02 
         EOS013 = 3.8306136687e-01 
         #
         ALP000 = -2.5218796628e-01 
         ALP100 = 3.4119354654e-01 
         ALP200 = -2.2119589983e-01 
         ALP300 = 1.8082347094e-01 
         ALP400 = -3.6936026529e-02 
         ALP500 = -5.0091801383e-03 
         ALP010 = 1.2789915300 
         ALP110 = -1.2021756164 
         ALP210 = 8.4037519952e-01 
         ALP310 = -4.1905788542e-01 
         ALP410 = 9.8855300959e-02 
         ALP020 = -1.2634838399 
         ALP120 = 1.6112195176 
         ALP220 = -7.5817155402e-01 
         ALP320 = 4.7006963580e-02 
         ALP030 = 8.0812310102e-01 
         ALP130 = -1.0102374985 
         ALP230 = 4.8340368631e-01 
         ALP040 = -1.5098959754e-01 
         ALP140 = -1.4394226233e-02 
         ALP050 = 3.6780433255e-02 
         ALP001 = 3.9631611467e-01 
         ALP101 = 1.9159845880e-02 
         ALP201 = -1.0286156825e-01 
         ALP301 = 1.6738969362e-02 
         ALP011 = -4.9997430930e-01 
         ALP111 = 9.7335338937e-03 
         ALP211 = 6.0887771651e-02 
         ALP021 = 2.6149576513e-01 
         ALP121 = -1.6671866715e-02 
         ALP031 = -5.9503008642e-02 
         ALP002 = -5.4590812035e-02 
         ALP102 = 8.6134185799e-03 
         ALP012 = 6.2740815484e-02 
         ALP003 = -9.5765341718e-03 
         #
         BET000 = 2.1420623987 
         BET100 = -9.3752598635 
         BET200 = 1.9446303907e+01 
         BET300 = -1.8632235232e+01 
         BET400 = 8.9390837485 
         BET500 = -1.7142465871 
         BET010 = -1.7059677327e-01 
         BET110 = 2.2119589983e-01 
         BET210 = -2.7123520642e-01 
         BET310 = 7.3872053057e-02 
         BET410 = 1.2522950346e-02 
         BET020 = 3.0054390409e-01 
         BET120 = -4.2018759976e-01 
         BET220 = 3.1429341406e-01 
         BET320 = -9.8855300959e-02 
         BET030 = -2.6853658626e-01 
         BET130 = 2.5272385134e-01 
         BET230 = -2.3503481790e-02 
         BET040 = 1.2627968731e-01 
         BET140 = -1.2085092158e-01 
         BET050 = 1.4394226233e-03 
         BET001 = -2.2271304375e-01 
         BET101 = 5.5453416919e-01 
         BET201 = -6.2815936268e-01 
         BET301 = 2.0601115202e-01 
         BET011 = -9.5799229402e-03 
         BET111 = 1.0286156825e-01 
         BET211 = -2.5108454043e-02 
         BET021 = -2.4333834734e-03 
         BET121 = -3.0443885826e-02 
         BET031 = 2.7786444526e-03 
         BET002 = -4.2811838287e-02 
         BET102 = 5.1355066072e-02 
         BET012 = -4.3067092900e-03 
         BET003 = -7.1548119050e-04 
         #
         PEN000 = -5.3743005340 
         PEN100 = 8.9085217499 
         PEN200 = -1.1090683384e+01 
         PEN300 = 8.3754581690 
         PEN400 = -2.0601115202 
         PEN010 = 7.9263222935 
         PEN110 = 3.8319691761e-01 
         PEN210 = -2.0572313651 
         PEN310 = 3.3477938724e-01 
         PEN020 = -4.9997430930 
         PEN120 = 9.7335338937e-02 
         PEN220 = 6.0887771651e-01 
         PEN030 = 1.7433051009 
         PEN130 = -1.1114577810e-01 
         PEN040 = -2.9751504321e-01 
         PEN001 = -6.9171176978e-01 
         PEN101 = 2.2832980419 
         PEN201 = -1.3694684286 
         PEN011 = -1.4557549876 
         PEN111 = 2.2969116213e-01 
         PEN021 = 8.3654420645e-01 
         PEN002 = -1.4046808820e-02 
         PEN102 = 4.2928871430e-02 
         PEN012 = -2.8729602515e-01 
         #
         APE000 = -1.9815805734e-01 
         APE100 = -9.5799229402e-03 
         APE200 = 5.1430784127e-02 
         APE300 = -8.3694846809e-03 
         APE010 = 2.4998715465e-01 
         APE110 = -4.8667669469e-03 
         APE210 = -3.0443885826e-02 
         APE020 = -1.3074788257e-01 
         APE120 = 8.3359333577e-03 
         APE030 = 2.9751504321e-02 
         APE001 = 3.6393874690e-02 
         APE101 = -5.7422790533e-03 
         APE011 = -4.1827210323e-02 
         APE002 = 7.1824006288e-03 
         #
         BPE000 = 1.1135652187e-01 
         BPE100 = -2.7726708459e-01 
         BPE200 = 3.1407968134e-01 
         BPE300 = -1.0300557601e-01 
         BPE010 = 4.7899614701e-03 
         BPE110 = -5.1430784127e-02 
         BPE210 = 1.2554227021e-02 
         BPE020 = 1.2166917367e-03 
         BPE120 = 1.5221942913e-02 
         BPE030 = -1.3893222263e-03 
         BPE001 = 2.8541225524e-02 
         BPE101 = -3.4236710714e-02 
         BPE011 = 2.8711395266e-03 
         BPE002 = 5.3661089288e-04 
         #
      
      rau0_rcp    = rau0 * rcp 
      r1_rau0     = 1.  / rau0
      r1_rcp      = 1.  / rcp
      r1_rau0_rcp = 1.  / rau0_rcp 
      
      
      