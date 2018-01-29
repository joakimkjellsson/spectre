MODULE eos
!!
!!
!!  Equation of state calculations taken from NEMO
!!
!!  Joakim Kjellsson, AOPP, 2016
!!
!!
   
   IMPLICIT NONE
   
   INTEGER, PARAMETER :: wp = SELECTED_REAL_KIND(15,307)
   
   REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) ::  tem, sal, dep, depw, &
   &                                           alpha, beta, rho, rhop, N2, &
   &                                           PE_an, dPE_dT, dPE_dS, &
   &                                           tmask, dzw
   REAL(wp), ALLOCATABLE, DIMENSION(:,:,:,:) ::  tem_ran, sal_ran
   INTEGER                                 ::  nn_eos     ! = 0/1/2 type of eq. of state and Brunt-Vaisala frequ.
   INTEGER                                 ::  nn_sto_eos 
   INTEGER                                 ::  jpi, fs_2, fs_jpim1, jpj, jpjm1, jpk, jpkm1
   LOGICAL                                 ::  ln_useCT   ! determine if eos_pt_from_ct is used to compute sst_m
   LOGICAL                                 ::  ln_sto_eos ! stochastic EOS
   
   REAL(wp) ::   rn_a0      = 1.6550e-1_wp     ! thermal expansion coeff. 
   REAL(wp) ::   rn_b0      = 7.6554e-1_wp     ! saline  expansion coeff. 
   REAL(wp) ::   rn_lambda1 = 5.9520e-2_wp     ! cabbeling coeff. in T^2        
   REAL(wp) ::   rn_lambda2 = 5.4914e-4_wp     ! cabbeling coeff. in S^2        
   REAL(wp) ::   rn_mu1     = 1.4970e-4_wp     ! thermobaric coeff. in T  
   REAL(wp) ::   rn_mu2     = 1.1090e-5_wp     ! thermobaric coeff. in S  
   REAL(wp) ::   rn_nu      = 2.4341e-3_wp     ! cabbeling coeff. in theta*salt  
   
   ! TEOS10/EOS80 parameters
   REAL(wp) ::   rau0, rcp, r1_S0, r1_T0, r1_Z0, rdeltaS, r1_rau0, r1_rau0_rcp, rau0_rcp, r1_rcp
   
   ! EOS parameters
   REAL(wp) ::   EOS000 , EOS100 , EOS200 , EOS300 , EOS400 , EOS500 , EOS600
   REAL(wp) ::   EOS010 , EOS110 , EOS210 , EOS310 , EOS410 , EOS510
   REAL(wp) ::   EOS020 , EOS120 , EOS220 , EOS320 , EOS420
   REAL(wp) ::   EOS030 , EOS130 , EOS230 , EOS330
   REAL(wp) ::   EOS040 , EOS140 , EOS240
   REAL(wp) ::   EOS050 , EOS150
   REAL(wp) ::   EOS060
   REAL(wp) ::   EOS001 , EOS101 , EOS201 , EOS301 , EOS401
   REAL(wp) ::   EOS011 , EOS111 , EOS211 , EOS311
   REAL(wp) ::   EOS021 , EOS121 , EOS221
   REAL(wp) ::   EOS031 , EOS131
   REAL(wp) ::   EOS041
   REAL(wp) ::   EOS002 , EOS102 , EOS202
   REAL(wp) ::   EOS012 , EOS112
   REAL(wp) ::   EOS022
   REAL(wp) ::   EOS003 , EOS103
   REAL(wp) ::   EOS013 
   
   ! ALPHA parameters
   REAL(wp) ::   ALP000 , ALP100 , ALP200 , ALP300 , ALP400 , ALP500
   REAL(wp) ::   ALP010 , ALP110 , ALP210 , ALP310 , ALP410
   REAL(wp) ::   ALP020 , ALP120 , ALP220 , ALP320
   REAL(wp) ::   ALP030 , ALP130 , ALP230
   REAL(wp) ::   ALP040 , ALP140
   REAL(wp) ::   ALP050
   REAL(wp) ::   ALP001 , ALP101 , ALP201 , ALP301
   REAL(wp) ::   ALP011 , ALP111 , ALP211
   REAL(wp) ::   ALP021 , ALP121
   REAL(wp) ::   ALP031
   REAL(wp) ::   ALP002 , ALP102
   REAL(wp) ::   ALP012
   REAL(wp) ::   ALP003
   
   ! BETA parameters
   REAL(wp) ::   BET000 , BET100 , BET200 , BET300 , BET400 , BET500
   REAL(wp) ::   BET010 , BET110 , BET210 , BET310 , BET410
   REAL(wp) ::   BET020 , BET120 , BET220 , BET320
   REAL(wp) ::   BET030 , BET130 , BET230
   REAL(wp) ::   BET040 , BET140
   REAL(wp) ::   BET050
   REAL(wp) ::   BET001 , BET101 , BET201 , BET301
   REAL(wp) ::   BET011 , BET111 , BET211
   REAL(wp) ::   BET021 , BET121
   REAL(wp) ::   BET031
   REAL(wp) ::   BET002 , BET102
   REAL(wp) ::   BET012
   REAL(wp) ::   BET003

   ! PEN parameters
   REAL(wp) ::   PEN000 , PEN100 , PEN200 , PEN300 , PEN400
   REAL(wp) ::   PEN010 , PEN110 , PEN210 , PEN310
   REAL(wp) ::   PEN020 , PEN120 , PEN220
   REAL(wp) ::   PEN030 , PEN130
   REAL(wp) ::   PEN040
   REAL(wp) ::   PEN001 , PEN101 , PEN201
   REAL(wp) ::   PEN011 , PEN111
   REAL(wp) ::   PEN021
   REAL(wp) ::   PEN002 , PEN102
   REAL(wp) ::   PEN012
   
   ! ALPHA_PEN parameters
   REAL(wp) ::   APE000 , APE100 , APE200 , APE300
   REAL(wp) ::   APE010 , APE110 , APE210
   REAL(wp) ::   APE020 , APE120
   REAL(wp) ::   APE030
   REAL(wp) ::   APE001 , APE101
   REAL(wp) ::   APE011
   REAL(wp) ::   APE002

   ! BETA_PEN parameters
   REAL(wp) ::   BPE000 , BPE100 , BPE200 , BPE300
   REAL(wp) ::   BPE010 , BPE110 , BPE210
   REAL(wp) ::   BPE020 , BPE120
   REAL(wp) ::   BPE030
   REAL(wp) ::   BPE001 , BPE101
   REAL(wp) ::   BPE011
   REAL(wp) ::   BPE002
   
   CONTAINS

   SUBROUTINE eos_insitu
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
      INTEGER  ::   ji, jj, jk                ! dummy loop indices
      REAL(wp) ::   zt , zh , zs , ztm        ! local scalars
      REAL(wp) ::   zn , zn0, zn1, zn2, zn3   !   -      -
      !!----------------------------------------------------------------------
      !
      jpi = size(tem,1)
      jpj = size(tem,2)
      jpk = size(tem,3)
      jpkm1 = jpk-1
      !
      IF (.NOT. ALLOCATED(rho)) THEN 
         ALLOCATE( rho(jpi,jpj,jpk) )
      END IF
      !
      SELECT CASE( nn_eos )
      !
      CASE( -1, 0 )                !==  polynomial TEOS-10 / EOS-80 ==!
         !
         DO jk = 1, jpkm1
            DO jj = 1, jpj
               DO ji = 1, jpi
                  !
                  zh  = dep (ji,jj,jk) * r1_Z0                                  ! depth
                  zt  = tem (ji,jj,jk) * r1_T0                           ! temperature
                  zs  = SQRT( ABS( sal(ji,jj,jk) + rdeltaS ) * r1_S0 )   ! square root salinity
                  ztm = tmask(ji,jj,jk)                                         ! tmask
                  !
                  zn3 = EOS013*zt   &
                     &   + EOS103*zs+EOS003
                     !
                  zn2 = (EOS022*zt   &
                     &   + EOS112*zs+EOS012)*zt   &
                     &   + (EOS202*zs+EOS102)*zs+EOS002
                     !
                  zn1 = (((EOS041*zt   &
                     &   + EOS131*zs+EOS031)*zt   &
                     &   + (EOS221*zs+EOS121)*zs+EOS021)*zt   &
                     &   + ((EOS311*zs+EOS211)*zs+EOS111)*zs+EOS011)*zt   &
                     &   + (((EOS401*zs+EOS301)*zs+EOS201)*zs+EOS101)*zs+EOS001
                     !
                  zn0 = (((((EOS060*zt   &
                     &   + EOS150*zs+EOS050)*zt   &
                     &   + (EOS240*zs+EOS140)*zs+EOS040)*zt   &
                     &   + ((EOS330*zs+EOS230)*zs+EOS130)*zs+EOS030)*zt   &
                     &   + (((EOS420*zs+EOS320)*zs+EOS220)*zs+EOS120)*zs+EOS020)*zt   &
                     &   + ((((EOS510*zs+EOS410)*zs+EOS310)*zs+EOS210)*zs+EOS110)*zs+EOS010)*zt   &
                     &   + (((((EOS600*zs+EOS500)*zs+EOS400)*zs+EOS300)*zs+EOS200)*zs+EOS100)*zs+EOS000
                     !
                  zn  = ( ( zn3 * zh + zn2 ) * zh + zn1 ) * zh + zn0
                  !
                  rho(ji,jj,jk) = (  zn * r1_rau0 - 1._wp  ) * ztm  ! density anomaly (masked)
                  !
               END DO
            END DO
         END DO
         !
      CASE( 1 )                !==  simplified EOS  ==!
         !
         DO jk = 1, jpkm1
            DO jj = 1, jpj
               DO ji = 1, jpi
                  zt  = tem (ji,jj,jk) - 10._wp
                  zs  = sal (ji,jj,jk) - 35._wp
                  zh  = dep (ji,jj,jk)
                  ztm = tmask(ji,jj,jk)
                  !
                  zn =  - rn_a0 * ( 1._wp + 0.5_wp*rn_lambda1*zt + rn_mu1*zh ) * zt   &
                     &  + rn_b0 * ( 1._wp - 0.5_wp*rn_lambda2*zs - rn_mu2*zh ) * zs   &
                     &  - rn_nu * zt * zs
                     !                                 
                  rho(ji,jj,jk) = zn * r1_rau0 * ztm                ! density anomaly (masked)
               END DO
            END DO
         END DO
         !
      END SELECT
      !
   END SUBROUTINE eos_insitu
   
   
   SUBROUTINE eos_insitu_pot
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
      !
      jpi = size(tem,1)
      jpj = size(tem,2)
      jpk = size(tem,3)
      jpkm1 = jpk
      !
      IF (.NOT. ALLOCATED(rho)) THEN 
         ALLOCATE( rho(jpi,jpj,jpk) )
      END IF
      IF (.NOT. ALLOCATED(rhop)) THEN 
         ALLOCATE( rhop(jpi,jpj,jpk) )
      END IF
      !
      SELECT CASE ( nn_eos )
      !
      CASE( -1, 0 )                !==  polynomial TEOS-10 / EOS-80 ==!
         !
         ! Stochastic equation of state
         IF ( ln_sto_eos ) THEN
            ALLOCATE(zn0_sto(1:2*nn_sto_eos))
            ALLOCATE(zn_sto(1:2*nn_sto_eos))
            ALLOCATE(zsign(1:2*nn_sto_eos))
            DO jsmp = 1, 2*nn_sto_eos, 2
              zsign(jsmp)   = 1._wp
              zsign(jsmp+1) = -1._wp
            END DO
            !
            DO jk = 1, jpkm1
               DO jj = 1, jpj
                  DO ji = 1, jpi
                     !
                     ! compute density (2*nn_sto_eos) times:
                     ! (1) for t+dt, s+ds (with the random TS fluctutation computed in sto_pts)
                     ! (2) for t-dt, s-ds (with the opposite fluctuation)
                     DO jsmp = 1, nn_sto_eos*2
                        jdof   = (jsmp + 1) / 2
                        zh     = dep (ji,jj,jk) * r1_Z0                                  ! depth
                        zt     = (tem (ji,jj,jk) + tem_ran(ji,jj,jk,jdof) * zsign(jsmp)) * r1_T0    ! temperature
                        zstemp = sal  (ji,jj,jk) + sal_ran(ji,jj,jk,jdof) * zsign(jsmp)
                        zs     = SQRT( ABS( zstemp + rdeltaS ) * r1_S0 )   ! square root salinity
                        ztm    = tmask(ji,jj,jk)                                         ! tmask
                        !
                        zn3 = EOS013*zt   &
                           &   + EOS103*zs+EOS003
                           !
                        zn2 = (EOS022*zt   &
                           &   + EOS112*zs+EOS012)*zt   &
                           &   + (EOS202*zs+EOS102)*zs+EOS002
                           !
                        zn1 = (((EOS041*zt   &
                           &   + EOS131*zs+EOS031)*zt   &
                           &   + (EOS221*zs+EOS121)*zs+EOS021)*zt   &
                           &   + ((EOS311*zs+EOS211)*zs+EOS111)*zs+EOS011)*zt   &
                           &   + (((EOS401*zs+EOS301)*zs+EOS201)*zs+EOS101)*zs+EOS001
                           !
                        zn0_sto(jsmp) = (((((EOS060*zt   &
                           &   + EOS150*zs+EOS050)*zt   &
                           &   + (EOS240*zs+EOS140)*zs+EOS040)*zt   &
                           &   + ((EOS330*zs+EOS230)*zs+EOS130)*zs+EOS030)*zt   &
                           &   + (((EOS420*zs+EOS320)*zs+EOS220)*zs+EOS120)*zs+EOS020)*zt   &
                           &   + ((((EOS510*zs+EOS410)*zs+EOS310)*zs+EOS210)*zs+EOS110)*zs+EOS010)*zt   &
                           &   + (((((EOS600*zs+EOS500)*zs+EOS400)*zs+EOS300)*zs+EOS200)*zs+EOS100)*zs+EOS000
                           !
                        zn_sto(jsmp)  = ( ( zn3 * zh + zn2 ) * zh + zn1 ) * zh + zn0_sto(jsmp)
                     END DO
                     !
                     ! compute stochastic density as the mean of the (2*nn_sto_eos) densities
                     rhop(ji,jj,jk) = 0._wp ; rho(ji,jj,jk) = 0._wp
                     DO jsmp = 1, nn_sto_eos*2
                        rhop(ji,jj,jk) = rhop(ji,jj,jk) + zn0_sto(jsmp)  ! potential density referenced at the surface
                        !
                        rho(ji,jj,jk) = rho(ji,jj,jk) + (  zn_sto(jsmp) * r1_rau0 - 1._wp  )   ! density anomaly (masked)
                     END DO
                     rhop (ji,jj,jk) = 0.5_wp * rhop (ji,jj,jk) * ztm / nn_sto_eos
                     rho  (ji,jj,jk) = 0.5_wp * rho  (ji,jj,jk) * ztm / nn_sto_eos
                  END DO
               END DO
            END DO
            DEALLOCATE(zn0_sto,zn_sto,zsign)
            
         ! Non-stochastic equation of state
         ELSE
            DO jk = 1, jpkm1
               DO jj = 1, jpj
                  DO ji = 1, jpi
                     !
                     zh  = dep (ji,jj,jk) * r1_Z0                                  ! depth
                     zt  = tem (ji,jj,jk) * r1_T0                           ! temperature
                     zs  = SQRT( ABS( sal(ji,jj,jk) + rdeltaS ) * r1_S0 )   ! square root salinity
                     ztm = tmask(ji,jj,jk)                                         ! tmask
                     !
                     zn3 = EOS013*zt   &
                        &   + EOS103*zs+EOS003
                        !
                     zn2 = (EOS022*zt   &
                        &   + EOS112*zs+EOS012)*zt   &
                        &   + (EOS202*zs+EOS102)*zs+EOS002
                        !
                     zn1 = (((EOS041*zt   &
                        &   + EOS131*zs+EOS031)*zt   &
                        &   + (EOS221*zs+EOS121)*zs+EOS021)*zt   &
                        &   + ((EOS311*zs+EOS211)*zs+EOS111)*zs+EOS011)*zt   &
                        &   + (((EOS401*zs+EOS301)*zs+EOS201)*zs+EOS101)*zs+EOS001
                        !
                     zn0 = (((((EOS060*zt   &
                        &   + EOS150*zs+EOS050)*zt   &
                        &   + (EOS240*zs+EOS140)*zs+EOS040)*zt   &
                        &   + ((EOS330*zs+EOS230)*zs+EOS130)*zs+EOS030)*zt   &
                        &   + (((EOS420*zs+EOS320)*zs+EOS220)*zs+EOS120)*zs+EOS020)*zt   &
                        &   + ((((EOS510*zs+EOS410)*zs+EOS310)*zs+EOS210)*zs+EOS110)*zs+EOS010)*zt   &
                        &   + (((((EOS600*zs+EOS500)*zs+EOS400)*zs+EOS300)*zs+EOS200)*zs+EOS100)*zs+EOS000
                        !
                     zn  = ( ( zn3 * zh + zn2 ) * zh + zn1 ) * zh + zn0
                     !
                     rhop(ji,jj,jk) = zn0 * ztm                           ! potential density referenced at the surface
                     !
                     rho(ji,jj,jk) = (  zn * r1_rau0 - 1._wp  ) * ztm      ! density anomaly (masked)
                  END DO
               END DO
            END DO
         ENDIF
         
      CASE( 1 )                !==  simplified EOS  ==!
         !
         DO jk = 1, jpkm1
            DO jj = 1, jpj
               DO ji = 1, jpi
                  zt  = tem  (ji,jj,jk) - 10._wp
                  zs  = sal  (ji,jj,jk) - 35._wp
                  zh  = dep  (ji,jj,jk)
                  ztm = tmask(ji,jj,jk)
                  !                                                     ! potential density referenced at the surface
                  zn =  - rn_a0 * ( 1._wp + 0.5_wp*rn_lambda1*zt ) * zt   &
                     &  + rn_b0 * ( 1._wp - 0.5_wp*rn_lambda2*zs ) * zs   &
                     &  - rn_nu * zt * zs
                  rhop(ji,jj,jk) = ( rau0 + zn ) * ztm
                  !                                                     ! density anomaly (masked)
                  zn = zn - ( rn_a0 * rn_mu1 * zt + rn_b0 * rn_mu2 * zs ) * zh
                  rho(ji,jj,jk) = zn * r1_rau0 * ztm
                  !
               END DO
            END DO
         END DO
         !
      END SELECT
      !
   END SUBROUTINE eos_insitu_pot
   
   
   SUBROUTINE rab_3d
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
      !
      jpi = size(tem,1)
      jpj = size(tem,2)
      jpk = size(tem,3)
      jpkm1 = jpk-1
      !
      IF (.NOT. ALLOCATED(alpha)) THEN 
         ALLOCATE( alpha(jpi,jpj,jpk) )
      END IF
      IF (.NOT. ALLOCATED(beta)) THEN 
         ALLOCATE( beta(jpi,jpj,jpk) )
      END IF
      !
      SELECT CASE ( nn_eos )
      !
      CASE( -1, 0 )                !==  polynomial TEOS-10 / EOS-80 ==!
         !
         DO jk = 1, jpkm1
            DO jj = 1, jpj
               DO ji = 1, jpi
                  !
                  zh  = dep (ji,jj,jk) * r1_Z0                                ! depth
                  zt  = tem (ji,jj,jk) * r1_T0                           ! temperature
                  zs  = SQRT( ABS( sal(ji,jj,jk) + rdeltaS ) * r1_S0 )   ! square root salinity
                  ztm = tmask(ji,jj,jk)                                         ! tmask
                  !
                  ! alpha
                  zn3 = ALP003
                  !
                  zn2 = ALP012*zt + ALP102*zs+ALP002
                  !
                  zn1 = ((ALP031*zt   &
                     &   + ALP121*zs+ALP021)*zt   &
                     &   + (ALP211*zs+ALP111)*zs+ALP011)*zt   &
                     &   + ((ALP301*zs+ALP201)*zs+ALP101)*zs+ALP001
                     !
                  zn0 = ((((ALP050*zt   &
                     &   + ALP140*zs+ALP040)*zt   &
                     &   + (ALP230*zs+ALP130)*zs+ALP030)*zt   &
                     &   + ((ALP320*zs+ALP220)*zs+ALP120)*zs+ALP020)*zt   &
                     &   + (((ALP410*zs+ALP310)*zs+ALP210)*zs+ALP110)*zs+ALP010)*zt   &
                     &   + ((((ALP500*zs+ALP400)*zs+ALP300)*zs+ALP200)*zs+ALP100)*zs+ALP000
                     !
                  zn  = ( ( zn3 * zh + zn2 ) * zh + zn1 ) * zh + zn0
                  !
                  alpha (ji,jj,jk) = zn * r1_rau0 * ztm
                  if (ji == 250 .and. jj == 50 .and. jk == 1) then
                     print*,zn3,zn2,zn1,zn0,zt,zs,zh,alpha(ji,jj,jk)
                  end if   
                  !
                  ! beta
                  zn3 = BET003
                  !
                  zn2 = BET012*zt + BET102*zs+BET002
                  !
                  zn1 = ((BET031*zt   &
                     &   + BET121*zs+BET021)*zt   &
                     &   + (BET211*zs+BET111)*zs+BET011)*zt   &
                     &   + ((BET301*zs+BET201)*zs+BET101)*zs+BET001
                     !
                  zn0 = ((((BET050*zt   &
                     &   + BET140*zs+BET040)*zt   &
                     &   + (BET230*zs+BET130)*zs+BET030)*zt   &
                     &   + ((BET320*zs+BET220)*zs+BET120)*zs+BET020)*zt   &
                     &   + (((BET410*zs+BET310)*zs+BET210)*zs+BET110)*zs+BET010)*zt   &
                     &   + ((((BET500*zs+BET400)*zs+BET300)*zs+BET200)*zs+BET100)*zs+BET000
                     !
                  zn  = ( ( zn3 * zh + zn2 ) * zh + zn1 ) * zh + zn0
                  !
                  beta (ji,jj,jk) = zn / zs * r1_rau0 * ztm
                  !
               END DO
            END DO
         END DO
         !
      CASE( 1 )                  !==  simplified EOS  ==!
         !
         DO jk = 1, jpkm1
            DO jj = 1, jpj
               DO ji = 1, jpi
                  zt  = tem (ji,jj,jk) - 10._wp   ! pot. temperature anomaly (t-T0)
                  zs  = sal (ji,jj,jk) - 35._wp   ! abs. salinity anomaly (s-S0)
                  zh  = dep (ji,jj,jk)                 ! depth in meters at t-point
                  ztm = tmask(ji,jj,jk)                  ! land/sea bottom mask = surf. mask
                  !
                  zn  = rn_a0 * ( 1._wp + rn_lambda1*zt + rn_mu1*zh ) + rn_nu*zs
                  alpha (ji,jj,jk) = zn * r1_rau0 * ztm   ! alpha
                  !
                  zn  = rn_b0 * ( 1._wp - rn_lambda2*zs - rn_mu2*zh ) - rn_nu*zt
                  beta (ji,jj,jk) = zn * r1_rau0 * ztm   ! beta
                  !
               END DO
            END DO
         END DO
         !
      CASE DEFAULT
         PRINT*, '          bad flag value for nn_eos = ', nn_eos
         STOP
         !
      END SELECT
      !
   END SUBROUTINE rab_3d
   
   
   SUBROUTINE bn2
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
      INTEGER  ::   ji, jj, jk      ! dummy loop indices
      REAL(wp) ::   zaw, zbw, zrw   ! local scalars
      !!----------------------------------------------------------------------
      !
      jpi = size(tem,1)
      jpj = size(tem,2)
      jpk = size(tem,3)
      jpkm1 = jpk-1
      !
      IF (.NOT. ALLOCATED(N2)) THEN 
         ALLOCATE( N2(jpi,jpj,jpk) )
      END IF
      !
      DO jk = 2, jpkm1           ! interior points only (2=< jk =< jpkm1 )
         DO jj = 1, jpj          ! surface and bottom value set to zero one for all in istate.F90
            DO ji = 1, jpi
               zrw =   ( depw(ji,jj,jk  ) - dep(ji,jj,jk) )   &
                  &  / ( dep (ji,jj,jk-1) - dep(ji,jj,jk) ) 
                  !
               zaw = alpha(ji,jj,jk) * (1. - zrw) + alpha(ji,jj,jk-1) * zrw 
               zbw = beta (ji,jj,jk) * (1. - zrw) + beta (ji,jj,jk-1) * zrw
               !
               N2(ji,jj,jk) = 9.81 * (   zaw * ( tem(ji,jj,jk-1) - tem(ji,jj,jk) )     &
                  &                    - zbw * ( sal(ji,jj,jk-1) - sal(ji,jj,jk) )  )  &
                  &            / dzw(ji,jj,jk) * tmask(ji,jj,jk)
            END DO
         END DO
      END DO
      !
   END SUBROUTINE bn2

   
   SUBROUTINE eos_pen
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
      !
      jpi = size(tem,1)
      jpj = size(tem,2)
      jpk = size(tem,3)
      jpkm1 = jpk-1
      !
      IF (.NOT. ALLOCATED(PE_an)) THEN 
         ALLOCATE( PE_an(jpi,jpj,jpk) )
      END IF
      IF (.NOT. ALLOCATED(dPE_dT)) THEN 
         ALLOCATE( dPE_dT(jpi,jpj,jpk) )
      END IF
      IF (.NOT. ALLOCATED(dPE_dS)) THEN 
         ALLOCATE( dPE_dS(jpi,jpj,jpk) )
      END IF
      !
      SELECT CASE ( nn_eos )
      !
      CASE( -1, 0 )                !==  polynomial TEOS-10 / EOS-80 ==!
         !
         DO jk = 1, jpkm1
            DO jj = 1, jpj
               DO ji = 1, jpi
                  !
                  zh  = dep (ji,jj,jk) * r1_Z0                                ! depth
                  zt  = tem (ji,jj,jk) * r1_T0                           ! temperature
                  zs  = SQRT( ABS( sal(ji,jj,jk) + rdeltaS ) * r1_S0 )   ! square root salinity
                  ztm = tmask(ji,jj,jk)                                         ! tmask
                  !
                  ! potential energy non-linear anomaly
                  zn2 = (PEN012)*zt   &
                     &   + PEN102*zs+PEN002
                     !
                  zn1 = ((PEN021)*zt   &
                     &   + PEN111*zs+PEN011)*zt   &
                     &   + (PEN201*zs+PEN101)*zs+PEN001
                     !
                  zn0 = ((((PEN040)*zt   &
                     &   + PEN130*zs+PEN030)*zt   &
                     &   + (PEN220*zs+PEN120)*zs+PEN020)*zt   &
                     &   + ((PEN310*zs+PEN210)*zs+PEN110)*zs+PEN010)*zt   &
                     &   + (((PEN400*zs+PEN300)*zs+PEN200)*zs+PEN100)*zs+PEN000
                     !
                  zn  = ( zn2 * zh + zn1 ) * zh + zn0
                  !
                  PE_an(ji,jj,jk)  = zn * zh * r1_rau0 * ztm
                  !
                  ! alphaPE non-linear anomaly
                  zn2 = APE002
                  !
                  zn1 = (APE011)*zt   &
                     &   + APE101*zs+APE001
                     !
                  zn0 = (((APE030)*zt   &
                     &   + APE120*zs+APE020)*zt   &
                     &   + (APE210*zs+APE110)*zs+APE010)*zt   &
                     &   + ((APE300*zs+APE200)*zs+APE100)*zs+APE000
                     !
                  zn  = ( zn2 * zh + zn1 ) * zh + zn0
                  !                              
                  dPE_dT(ji,jj,jk) = zn * zh * r1_rau0 * ztm
                  !
                  ! betaPE non-linear anomaly
                  zn2 = BPE002
                  !
                  zn1 = (BPE011)*zt   &
                     &   + BPE101*zs+BPE001
                     !
                  zn0 = (((BPE030)*zt   &
                     &   + BPE120*zs+BPE020)*zt   &
                     &   + (BPE210*zs+BPE110)*zs+BPE010)*zt   &
                     &   + ((BPE300*zs+BPE200)*zs+BPE100)*zs+BPE000
                     !
                  zn  = ( zn2 * zh + zn1 ) * zh + zn0
                  !                              
                  dPE_dS(ji,jj,jk) = zn / zs * zh * r1_rau0 * ztm
                  !
               END DO
            END DO
         END DO
         !
      CASE( 1 )                !==  Vallis (2006) simplified EOS  ==!
         !
         DO jk = 1, jpkm1
            DO jj = 1, jpj
               DO ji = 1, jpi
                  zt  = tem(ji,jj,jk) - 10._wp  ! temperature anomaly (t-T0)
                  zs  = sal(ji,jj,jk) - 35._wp  ! abs. salinity anomaly (s-S0)
                  zh  = dep(ji,jj,jk)               ! depth in meters  at t-point
                  ztm = tmask(ji,jj,jk)                ! tmask
                  zn  = 0.5_wp * zh * r1_rau0 * ztm
                  !                                    ! Potential Energy
                  PE_an(ji,jj,jk) = ( rn_a0 * rn_mu1 * zt + rn_b0 * rn_mu2 * zs ) * zn
                  !                                    ! alphaPE
                  dPE_dT(ji,jj,jk) = - rn_a0 * rn_mu1 * zn
                  dPE_dS(ji,jj,jk) =   rn_b0 * rn_mu2 * zn
                  !
               END DO
            END DO
         END DO
         !
      CASE DEFAULT
         PRINT*, '          bad flag value for nn_eos = ', nn_eos
         STOP
         !
      END SELECT
      !
   END SUBROUTINE eos_pen

   SUBROUTINE eos_init
      
      rau0        = 1026._wp                 !: volumic mass of reference     [kg/m3]
      rcp         = 3991.86795711963_wp      !: heat capacity     [J/K]
      
      SELECT CASE( nn_eos )         ! check option
      !
      CASE( -1 )                       !==  polynomial TEOS-10  ==!
         !
         rdeltaS = 32._wp
         r1_S0  = 0.875_wp/35.16504_wp
         r1_T0  = 1._wp/40._wp
         r1_Z0  = 1.e-4_wp
         !
         EOS000 = 8.0189615746e+02_wp
         EOS100 = 8.6672408165e+02_wp
         EOS200 = -1.7864682637e+03_wp
         EOS300 = 2.0375295546e+03_wp
         EOS400 = -1.2849161071e+03_wp
         EOS500 = 4.3227585684e+02_wp
         EOS600 = -6.0579916612e+01_wp
         EOS010 = 2.6010145068e+01_wp
         EOS110 = -6.5281885265e+01_wp
         EOS210 = 8.1770425108e+01_wp
         EOS310 = -5.6888046321e+01_wp
         EOS410 = 1.7681814114e+01_wp
         EOS510 = -1.9193502195_wp
         EOS020 = -3.7074170417e+01_wp
         EOS120 = 6.1548258127e+01_wp
         EOS220 = -6.0362551501e+01_wp
         EOS320 = 2.9130021253e+01_wp
         EOS420 = -5.4723692739_wp
         EOS030 = 2.1661789529e+01_wp
         EOS130 = -3.3449108469e+01_wp
         EOS230 = 1.9717078466e+01_wp
         EOS330 = -3.1742946532_wp
         EOS040 = -8.3627885467_wp
         EOS140 = 1.1311538584e+01_wp
         EOS240 = -5.3563304045_wp
         EOS050 = 5.4048723791e-01_wp
         EOS150 = 4.8169980163e-01_wp
         EOS060 = -1.9083568888e-01_wp
         EOS001 = 1.9681925209e+01_wp
         EOS101 = -4.2549998214e+01_wp
         EOS201 = 5.0774768218e+01_wp
         EOS301 = -3.0938076334e+01_wp
         EOS401 = 6.6051753097_wp
         EOS011 = -1.3336301113e+01_wp
         EOS111 = -4.4870114575_wp
         EOS211 = 5.0042598061_wp
         EOS311 = -6.5399043664e-01_wp
         EOS021 = 6.7080479603_wp
         EOS121 = 3.5063081279_wp
         EOS221 = -1.8795372996_wp
         EOS031 = -2.4649669534_wp
         EOS131 = -5.5077101279e-01_wp
         EOS041 = 5.5927935970e-01_wp
         EOS002 = 2.0660924175_wp
         EOS102 = -4.9527603989_wp
         EOS202 = 2.5019633244_wp
         EOS012 = 2.0564311499_wp
         EOS112 = -2.1311365518e-01_wp
         EOS022 = -1.2419983026_wp
         EOS003 = -2.3342758797e-02_wp
         EOS103 = -1.8507636718e-02_wp
         EOS013 = 3.7969820455e-01_wp
         !
         ALP000 = -6.5025362670e-01_wp
         ALP100 = 1.6320471316_wp
         ALP200 = -2.0442606277_wp
         ALP300 = 1.4222011580_wp
         ALP400 = -4.4204535284e-01_wp
         ALP500 = 4.7983755487e-02_wp
         ALP010 = 1.8537085209_wp
         ALP110 = -3.0774129064_wp
         ALP210 = 3.0181275751_wp
         ALP310 = -1.4565010626_wp
         ALP410 = 2.7361846370e-01_wp
         ALP020 = -1.6246342147_wp
         ALP120 = 2.5086831352_wp
         ALP220 = -1.4787808849_wp
         ALP320 = 2.3807209899e-01_wp
         ALP030 = 8.3627885467e-01_wp
         ALP130 = -1.1311538584_wp
         ALP230 = 5.3563304045e-01_wp
         ALP040 = -6.7560904739e-02_wp
         ALP140 = -6.0212475204e-02_wp
         ALP050 = 2.8625353333e-02_wp
         ALP001 = 3.3340752782e-01_wp
         ALP101 = 1.1217528644e-01_wp
         ALP201 = -1.2510649515e-01_wp
         ALP301 = 1.6349760916e-02_wp
         ALP011 = -3.3540239802e-01_wp
         ALP111 = -1.7531540640e-01_wp
         ALP211 = 9.3976864981e-02_wp
         ALP021 = 1.8487252150e-01_wp
         ALP121 = 4.1307825959e-02_wp
         ALP031 = -5.5927935970e-02_wp
         ALP002 = -5.1410778748e-02_wp
         ALP102 = 5.3278413794e-03_wp
         ALP012 = 6.2099915132e-02_wp
         ALP003 = -9.4924551138e-03_wp
         !
         BET000 = 1.0783203594e+01_wp
         BET100 = -4.4452095908e+01_wp
         BET200 = 7.6048755820e+01_wp
         BET300 = -6.3944280668e+01_wp
         BET400 = 2.6890441098e+01_wp
         BET500 = -4.5221697773_wp
         BET010 = -8.1219372432e-01_wp
         BET110 = 2.0346663041_wp
         BET210 = -2.1232895170_wp
         BET310 = 8.7994140485e-01_wp
         BET410 = -1.1939638360e-01_wp
         BET020 = 7.6574242289e-01_wp
         BET120 = -1.5019813020_wp
         BET220 = 1.0872489522_wp
         BET320 = -2.7233429080e-01_wp
         BET030 = -4.1615152308e-01_wp
         BET130 = 4.9061350869e-01_wp
         BET230 = -1.1847737788e-01_wp
         BET040 = 1.4073062708e-01_wp
         BET140 = -1.3327978879e-01_wp
         BET050 = 5.9929880134e-03_wp
         BET001 = -5.2937873009e-01_wp
         BET101 = 1.2634116779_wp
         BET201 = -1.1547328025_wp
         BET301 = 3.2870876279e-01_wp
         BET011 = -5.5824407214e-02_wp
         BET111 = 1.2451933313e-01_wp
         BET211 = -2.4409539932e-02_wp
         BET021 = 4.3623149752e-02_wp
         BET121 = -4.6767901790e-02_wp
         BET031 = -6.8523260060e-03_wp
         BET002 = -6.1618945251e-02_wp
         BET102 = 6.2255521644e-02_wp
         BET012 = -2.6514181169e-03_wp
         BET003 = -2.3025968587e-04_wp
         !
         PEN000 = -9.8409626043_wp
         PEN100 = 2.1274999107e+01_wp
         PEN200 = -2.5387384109e+01_wp
         PEN300 = 1.5469038167e+01_wp
         PEN400 = -3.3025876549_wp
         PEN010 = 6.6681505563_wp
         PEN110 = 2.2435057288_wp
         PEN210 = -2.5021299030_wp
         PEN310 = 3.2699521832e-01_wp
         PEN020 = -3.3540239802_wp
         PEN120 = -1.7531540640_wp
         PEN220 = 9.3976864981e-01_wp
         PEN030 = 1.2324834767_wp
         PEN130 = 2.7538550639e-01_wp
         PEN040 = -2.7963967985e-01_wp
         PEN001 = -1.3773949450_wp
         PEN101 = 3.3018402659_wp
         PEN201 = -1.6679755496_wp
         PEN011 = -1.3709540999_wp
         PEN111 = 1.4207577012e-01_wp
         PEN021 = 8.2799886843e-01_wp
         PEN002 = 1.7507069098e-02_wp
         PEN102 = 1.3880727538e-02_wp
         PEN012 = -2.8477365341e-01_wp
         !
         APE000 = -1.6670376391e-01_wp
         APE100 = -5.6087643219e-02_wp
         APE200 = 6.2553247576e-02_wp
         APE300 = -8.1748804580e-03_wp
         APE010 = 1.6770119901e-01_wp
         APE110 = 8.7657703198e-02_wp
         APE210 = -4.6988432490e-02_wp
         APE020 = -9.2436260751e-02_wp
         APE120 = -2.0653912979e-02_wp
         APE030 = 2.7963967985e-02_wp
         APE001 = 3.4273852498e-02_wp
         APE101 = -3.5518942529e-03_wp
         APE011 = -4.1399943421e-02_wp
         APE002 = 7.1193413354e-03_wp
         !
         BPE000 = 2.6468936504e-01_wp
         BPE100 = -6.3170583896e-01_wp
         BPE200 = 5.7736640125e-01_wp
         BPE300 = -1.6435438140e-01_wp
         BPE010 = 2.7912203607e-02_wp
         BPE110 = -6.2259666565e-02_wp
         BPE210 = 1.2204769966e-02_wp
         BPE020 = -2.1811574876e-02_wp
         BPE120 = 2.3383950895e-02_wp
         BPE030 = 3.4261630030e-03_wp
         BPE001 = 4.1079296834e-02_wp
         BPE101 = -4.1503681096e-02_wp
         BPE011 = 1.7676120780e-03_wp
         BPE002 = 1.7269476440e-04_wp
         !
      CASE( 0 )                        !==  polynomial EOS-80 formulation  ==!
         !
         rdeltaS = 20._wp
         r1_S0  = 1._wp/40._wp
         r1_T0  = 1._wp/40._wp
         r1_Z0  = 1.e-4_wp
         !
         EOS000 = 9.5356891948e+02_wp
         EOS100 = 1.7136499189e+02_wp
         EOS200 = -3.7501039454e+02_wp
         EOS300 = 5.1856810420e+02_wp
         EOS400 = -3.7264470465e+02_wp
         EOS500 = 1.4302533998e+02_wp
         EOS600 = -2.2856621162e+01_wp
         EOS010 = 1.0087518651e+01_wp
         EOS110 = -1.3647741861e+01_wp
         EOS210 = 8.8478359933_wp
         EOS310 = -7.2329388377_wp
         EOS410 = 1.4774410611_wp
         EOS510 = 2.0036720553e-01_wp
         EOS020 = -2.5579830599e+01_wp
         EOS120 = 2.4043512327e+01_wp
         EOS220 = -1.6807503990e+01_wp
         EOS320 = 8.3811577084_wp
         EOS420 = -1.9771060192_wp
         EOS030 = 1.6846451198e+01_wp
         EOS130 = -2.1482926901e+01_wp
         EOS230 = 1.0108954054e+01_wp
         EOS330 = -6.2675951440e-01_wp
         EOS040 = -8.0812310102_wp
         EOS140 = 1.0102374985e+01_wp
         EOS240 = -4.8340368631_wp
         EOS050 = 1.2079167803_wp
         EOS150 = 1.1515380987e-01_wp
         EOS060 = -2.4520288837e-01_wp
         EOS001 = 1.0748601068e+01_wp
         EOS101 = -1.7817043500e+01_wp
         EOS201 = 2.2181366768e+01_wp
         EOS301 = -1.6750916338e+01_wp
         EOS401 = 4.1202230403_wp
         EOS011 = -1.5852644587e+01_wp
         EOS111 = -7.6639383522e-01_wp
         EOS211 = 4.1144627302_wp
         EOS311 = -6.6955877448e-01_wp
         EOS021 = 9.9994861860_wp
         EOS121 = -1.9467067787e-01_wp
         EOS221 = -1.2177554330_wp
         EOS031 = -3.4866102017_wp
         EOS131 = 2.2229155620e-01_wp
         EOS041 = 5.9503008642e-01_wp
         EOS002 = 1.0375676547_wp
         EOS102 = -3.4249470629_wp
         EOS202 = 2.0542026429_wp
         EOS012 = 2.1836324814_wp
         EOS112 = -3.4453674320e-01_wp
         EOS022 = -1.2548163097_wp
         EOS003 = 1.8729078427e-02_wp
         EOS103 = -5.7238495240e-02_wp
         EOS013 = 3.8306136687e-01_wp
         !
         ALP000 = -2.5218796628e-01_wp
         ALP100 = 3.4119354654e-01_wp
         ALP200 = -2.2119589983e-01_wp
         ALP300 = 1.8082347094e-01_wp
         ALP400 = -3.6936026529e-02_wp
         ALP500 = -5.0091801383e-03_wp
         ALP010 = 1.2789915300_wp
         ALP110 = -1.2021756164_wp
         ALP210 = 8.4037519952e-01_wp
         ALP310 = -4.1905788542e-01_wp
         ALP410 = 9.8855300959e-02_wp
         ALP020 = -1.2634838399_wp
         ALP120 = 1.6112195176_wp
         ALP220 = -7.5817155402e-01_wp
         ALP320 = 4.7006963580e-02_wp
         ALP030 = 8.0812310102e-01_wp
         ALP130 = -1.0102374985_wp
         ALP230 = 4.8340368631e-01_wp
         ALP040 = -1.5098959754e-01_wp
         ALP140 = -1.4394226233e-02_wp
         ALP050 = 3.6780433255e-02_wp
         ALP001 = 3.9631611467e-01_wp
         ALP101 = 1.9159845880e-02_wp
         ALP201 = -1.0286156825e-01_wp
         ALP301 = 1.6738969362e-02_wp
         ALP011 = -4.9997430930e-01_wp
         ALP111 = 9.7335338937e-03_wp
         ALP211 = 6.0887771651e-02_wp
         ALP021 = 2.6149576513e-01_wp
         ALP121 = -1.6671866715e-02_wp
         ALP031 = -5.9503008642e-02_wp
         ALP002 = -5.4590812035e-02_wp
         ALP102 = 8.6134185799e-03_wp
         ALP012 = 6.2740815484e-02_wp
         ALP003 = -9.5765341718e-03_wp
         !
         BET000 = 2.1420623987_wp
         BET100 = -9.3752598635_wp
         BET200 = 1.9446303907e+01_wp
         BET300 = -1.8632235232e+01_wp
         BET400 = 8.9390837485_wp
         BET500 = -1.7142465871_wp
         BET010 = -1.7059677327e-01_wp
         BET110 = 2.2119589983e-01_wp
         BET210 = -2.7123520642e-01_wp
         BET310 = 7.3872053057e-02_wp
         BET410 = 1.2522950346e-02_wp
         BET020 = 3.0054390409e-01_wp
         BET120 = -4.2018759976e-01_wp
         BET220 = 3.1429341406e-01_wp
         BET320 = -9.8855300959e-02_wp
         BET030 = -2.6853658626e-01_wp
         BET130 = 2.5272385134e-01_wp
         BET230 = -2.3503481790e-02_wp
         BET040 = 1.2627968731e-01_wp
         BET140 = -1.2085092158e-01_wp
         BET050 = 1.4394226233e-03_wp
         BET001 = -2.2271304375e-01_wp
         BET101 = 5.5453416919e-01_wp
         BET201 = -6.2815936268e-01_wp
         BET301 = 2.0601115202e-01_wp
         BET011 = -9.5799229402e-03_wp
         BET111 = 1.0286156825e-01_wp
         BET211 = -2.5108454043e-02_wp
         BET021 = -2.4333834734e-03_wp
         BET121 = -3.0443885826e-02_wp
         BET031 = 2.7786444526e-03_wp
         BET002 = -4.2811838287e-02_wp
         BET102 = 5.1355066072e-02_wp
         BET012 = -4.3067092900e-03_wp
         BET003 = -7.1548119050e-04_wp
         !
         PEN000 = -5.3743005340_wp
         PEN100 = 8.9085217499_wp
         PEN200 = -1.1090683384e+01_wp
         PEN300 = 8.3754581690_wp
         PEN400 = -2.0601115202_wp
         PEN010 = 7.9263222935_wp
         PEN110 = 3.8319691761e-01_wp
         PEN210 = -2.0572313651_wp
         PEN310 = 3.3477938724e-01_wp
         PEN020 = -4.9997430930_wp
         PEN120 = 9.7335338937e-02_wp
         PEN220 = 6.0887771651e-01_wp
         PEN030 = 1.7433051009_wp
         PEN130 = -1.1114577810e-01_wp
         PEN040 = -2.9751504321e-01_wp
         PEN001 = -6.9171176978e-01_wp
         PEN101 = 2.2832980419_wp
         PEN201 = -1.3694684286_wp
         PEN011 = -1.4557549876_wp
         PEN111 = 2.2969116213e-01_wp
         PEN021 = 8.3654420645e-01_wp
         PEN002 = -1.4046808820e-02_wp
         PEN102 = 4.2928871430e-02_wp
         PEN012 = -2.8729602515e-01_wp
         !
         APE000 = -1.9815805734e-01_wp
         APE100 = -9.5799229402e-03_wp
         APE200 = 5.1430784127e-02_wp
         APE300 = -8.3694846809e-03_wp
         APE010 = 2.4998715465e-01_wp
         APE110 = -4.8667669469e-03_wp
         APE210 = -3.0443885826e-02_wp
         APE020 = -1.3074788257e-01_wp
         APE120 = 8.3359333577e-03_wp
         APE030 = 2.9751504321e-02_wp
         APE001 = 3.6393874690e-02_wp
         APE101 = -5.7422790533e-03_wp
         APE011 = -4.1827210323e-02_wp
         APE002 = 7.1824006288e-03_wp
         !
         BPE000 = 1.1135652187e-01_wp
         BPE100 = -2.7726708459e-01_wp
         BPE200 = 3.1407968134e-01_wp
         BPE300 = -1.0300557601e-01_wp
         BPE010 = 4.7899614701e-03_wp
         BPE110 = -5.1430784127e-02_wp
         BPE210 = 1.2554227021e-02_wp
         BPE020 = 1.2166917367e-03_wp
         BPE120 = 1.5221942913e-02_wp
         BPE030 = -1.3893222263e-03_wp
         BPE001 = 2.8541225524e-02_wp
         BPE101 = -3.4236710714e-02_wp
         BPE011 = 2.8711395266e-03_wp
         BPE002 = 5.3661089288e-04_wp
         !
      CASE (1)
         !
      CASE DEFAULT                     !==  ERROR in nn_eos  ==!
         PRINT*, '          bad flag value for nn_eos = ', nn_eos
         STOP
         !
      END SELECT
      !
      rau0_rcp    = rau0 * rcp 
      r1_rau0     = 1._wp / rau0
      r1_rcp      = 1._wp / rcp
      r1_rau0_rcp = 1._wp / rau0_rcp 
      !
   END SUBROUTINE eos_init
   
   SUBROUTINE clean
      
      if (allocated(tem)) deallocate(tem)
      if (allocated(sal)) deallocate(sal)
      if (allocated(dep)) deallocate(dep)
      if (allocated(depw)) deallocate(depw)
      if (allocated(alpha)) deallocate(alpha)
      if (allocated(beta)) deallocate(beta)
      if (allocated(rho)) deallocate(rho)
      if (allocated(rhop)) deallocate(rhop)
      if (allocated(N2)) deallocate(N2)
      if (allocated(PE_an)) deallocate(PE_an)
      if (allocated(dPE_dT)) deallocate(dPE_dS)
      if (allocated(tmask)) deallocate(tmask)
      if (allocated(dzw)) deallocate(dzw)
      
   
   END SUBROUTINE
   
END MODULE