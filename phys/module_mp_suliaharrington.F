      MODULE MODULE_MP_SULIAHARRINGTON
      USE module_wrf_error
!     USE module_utility, ONLY: WRFU_Clock, WRFU_Alarm  ! GT
!     USE module_domain, ONLY : HISTORY_ALARM, Is_alarm_tstep  ! GT

!     USE module_state_description
      USE module_dm
      USE module_comm_dm

      USE ieee_arithmetic
      IMPLICIT NONE

!     INCLUDE 'mpif.h'

      REAL, PARAMETER :: PI = 3.1415926535897932384626434
      REAL, PARAMETER :: SQRTPI = 0.9189385332046727417803297
      LOGICAL, PARAMETER :: PIRE_CHEM = .FALSE.

      PUBLIC  ::  MP_SULIAHARRINGTON
!     PUBLIC  ::  SULIAHARRINGTON_INIT
      PRIVATE  ::  POLYSVP1
      PRIVATE :: GAMMA
      PRIVATE :: PI, SQRTPI
      PRIVATE :: GAMMLN
!     PRIVATE :: EVOLVE
      PRIVATE :: FINDGTP
      PRIVATE :: LENCONV
      PRIVATE :: FLHS
      PRIVATE :: CAPACITANCE_GAMMA
      PRIVATE :: SULIAHARRINGTON_MICRO

      REAL, PRIVATE ::  rhoi    !BULK DENSITY OF CLOUD ICE
      REAL, PRIVATE ::  nu      !DISTRIBUTION SHAPE FACTOR
      REAL, PRIVATE ::  nuc     !ICE NUCLEATION CONCENTRAION (#/L)
      REAL, PRIVATE ::  rd      !GAS CONSTANT OF DRY AIR
      REAL, PRIVATE ::  cp      !SPECIFIC HEAT FOR DRY AIR (CONST P)
      REAL, PRIVATE ::  ao      !INITIAL CHARACTERISTIC A-AXIS LENGTH
      REAL, PRIVATE ::  li0     !INITIAL SEMI-MAJOR AXIS LENGTH
      REAL, PRIVATE ::  mi0     !INITIAL PARTICLE MASS
      REAL, PRIVATE ::  gammnu  !GAMMA DIST WITH SHAPE, NU
      REAL, PRIVATE ::  qsmall  !THRES MIN MIXING RATIO VALUE
      REAL, PRIVATE ::  rv      !GAS CONSTANT OF WATER VAPOR
      REAL, PRIVATE ::  bif     !FALL SPEED EXPONENTIAL COEFF
      REAL, PRIVATE ::  aif     !FALL SPEED LEADING COEFF
      REAL, PRIVATE ::  g       !ACCELERATION OF GRAVITY
      REAL, PRIVATE ::  cpw
      REAL, PRIVATE :: RIN      !ICE NUCLEAITON RADIUS,CONTACT
!                                 FREEZING, MEYERS

!     parameters for snow microphysics
      REAL, PRIVATE ::  DCS     !THRESHOLD AGGREGATION RADIUS
      REAL, PRIVATE ::  BS      !SNOW FALLSPEED PARAMETER
      REAL, PRIVATE ::  RHOSN   !SNOW/AGGREGATE DENSITY
      REAL, PRIVATE ::  EII     !ICE AGGREGATION EFFICIENCY
      REAL, PRIVATE ::  DS      !SNOW SIZE DIST PARAM
      REAL, PRIVATE ::  CS
      REAL, PRIVATE ::  f1s
      REAL, PRIVATE ::  f2s
      REAL, PRIVATE ::  lammins, lammaxs

!     parameters for rain microphysics
      REAL, PRIVATE :: f1r      ! VENTILATION PARAMETER FOR RAIN
      REAL, PRIVATE :: f2r      ! " "
      REAL, PRIVATE :: ar       ! FALL SPEED PARAMETER FOR RAIN
      REAL, PRIVATE :: as       ! FALL SPEED PARAMETER FOR SNOW
      REAL, PRIVATE :: rhosu    ! APPROXIMATE AIR DENSITY NEAR 850MB
      REAL, PRIVATE :: rhow     ! DENSITY OF LIQUID WATER
      REAL, PRIVATE :: br       ! FALLSPEED PARAMETER FOR RAIN
      REAL, PRIVATE :: bc       ! FALLSPEED PARAMETER FOR CLOUD WATER
      REAL, PRIVATE :: aimm     ! FREEZING PARAMETER FOR RAIN
      REAL, PRIVATE :: bimm     ! FREEZING PARAMETER FOR RAIN
      REAL, PRIVATE :: lamminr, lammaxr

!     parameters for graupel microphysics
      REAL, PRIVATE :: AG       !FALL SPEED PARAMETER FOR GRAUPEL
      REAL, PRIVATE :: BG       !FALL SPEED PARAMETER FOR GRAUPEL
      REAL, PRIVATE :: RHOG     !GRAUPEL DENSITY
      REAL, PRIVATE :: MG0      !MASS OF EMBRYO GRAUPEL
      REAL, PRIVATE :: CG       !SIZE DISTRIBUTION PARAMETER FOR GRAUPEL
      REAL, PRIVATE :: DG       !SIZE DISTRIBUTION PARAMETER FOR GRAUPEL
      REAL, PRIVATE :: ECI      !COLLECTION EFFICIENCY,ICE-DROPLET
      REAL, PRIVATE :: ECR      !COLLECTION EFFICIENCI B/N DROPLETS/RAIN & SNOW/RAIN
      REAL, PRIVATE :: LAMMING
      REAL, PRIVATE :: LAMMAXG
      REAL, PRIVATE :: MMULT

      REAL, PRIVATE :: cons1,cons2,cons3,cons4,cons5,cons6,cons7,cons8,cons9    
      REAL, PRIVATE :: cons10,cons11,cons13,cons14,cons15,cons17,cons18,cons19
      REAL, PRIVATE :: cons20,cons23,cons24,cons25,cons26,cons29
      REAL, PRIVATE :: cons31,cons32,cons34,cons35,cons36,cons39   
      REAL, PRIVATE :: cons40,cons41
      REAL, PRIVATE :: FUDGE

  
CONTAINS
!*************************************************************************************
!     THIS SUBROUTINE INITIALIZES ALL PHYSICAL CONSTANTS AMND PARAMETERS 
!     NEEDED BY THE MICROPHYSICS SCHEME.
!     NEEDS TO BE CALLED AT FIRST TIMESTEP,
!     PRIOR TO CALL TO MAIN MICROPHYSICS INTERFACE
!-------------------------------------------------------------------------------------
      SUBROUTINE SULIAHARRINGTON_INIT
!************************************************************************************

      IMPLICIT NONE

      rhoi = 920.
!     rhoi = 500.
      nu = 4.
      nuc = 1.
      rd = 287.15
      cp = 1005.
      cpw = 4187.
      ao = .1e-6
      li0 = 1.e-6
      mi0 = 4./3.*pi*rhoi*(li0)**3
      gammnu = exp(gammln(nu))
      qsmall = 1.e-14
      rv = 461.5
      RIN = 0.1E-6

      bif = 0.99
      aif = 0.81*(1000.)**bif

!     add constants for snow microphysics !kjs 02/10/2015
      DCS = 125.e-6
      BS = 0.41
      DS = 3.
      rhosn = 100.
      CS = rhosn*pi/6.
      EII = 0.1
      lammins = 1./2000.e-6
      lammaxs = 1./10.e-6
      f1s = 0.86
      f2s = 0.28

!     add constants for graupel microphysics !kjs 02/10/2015
      RHOG = 400.
      AG = 19.3
      BG = 0.37
      CG = RHOG*pi/6.
      DG = 3.
      ECI = 0.7
      ECR = 1.0
      MG0 = 1.6E-10
      LAMMING = 1./2000.E-6
      LAMMAXG = 1./20.E-6
      MMULT = 4./3.*PI*RHOI*(5.e-6)**3

!     add constants for rain microphysics
      ar = 841.99667
      as = 11.72
      rhosu = 85000./(287.15*273.15)
      rhow = 997.
      f1r = 0.78
      f2r = 0.308
      f1s = 0.86
      f2s = 0.28
      g = 9.806
      br = 0.8
      bc = 2.0
      bimm = 100.
      aimm = 0.66
      lamminr = 1./2800.e-6
      lammaxr = 1./20.e-6

      cons1 = gamma(1.+DS)*CS
      cons2 = gamma(1.+DG)*CG
      cons3 = gamma(4.+BS)/6.
      cons4 = gamma(4.+ br)/6.
      cons5 = gamma(1.+BS)
      cons6 = gamma(1.+BR)
      cons7 = gamma(4.+BG)/6.
      cons8 = gamma(1.+BG)
      cons9 = gamma(5./2.+ br/2.)
      cons10 = gamma(5./2.+BS/2.)
      cons11 = gamma(5./2.+BG/2.)
      cons13 = gamma(BS+3.)*pi/4.*ECI
      cons14 = gamma(BG+3.)*pi/4.*ECI
      cons15 = -1108.*EII*pi**((1.-BS)/3.)*rhosn**((-2.-BS)/3.)/(4.*720.)
      cons17 = 4.*2.*3.*rhosu*pi*ECI*ECI*gamma(2.*BS+2.)/(8.*(rhog-rhosn))
      cons18 = rhosn*rhosn
      cons19 = rhow*rhow
      cons20 = 20.*pi*pi*rhow*bimm
      cons23 = pi/4.*EII*gamma(BS+3.) 
      cons24 = pi/4.*ECR*gamma(BR+3.)
      cons25 = pi*pi/24.*rhow*ECR*gamma(BR+6.)
      cons26 = pi/6.*rhow
      cons29 = 4./3.*pi*rhow*(25.e-6)**3
      cons31 = pi*pi*ECR*rhosn
      cons32 = pi/2.*ecr
      cons34 = 5./2.+ br/2.
      cons35 = 5./2.+BS/2.
      cons36 = 5./2.+BG/2.
      cons39 = pi*pi/36.*rhow*bimm
      cons40 = pi/6.*bimm
      cons41 = pi*pi*ecr*rhow

      FUDGE = 0.9999

      END SUBROUTINE SULIAHARRINGTON_INIT

!********************************************************************************
!     THIS SUBROUTINE IS MAIN INTERFACE WITH THE TWO-MOMENT MICROPHYSICS SCHEME
!     THIS INTERFACE TAKES IN 3D VARIABLES FROM DRIVER MODEL, CONVERTS TO 2D FOR
!     CALL TO THE MAIN MICROPHYSICS SUBROUTINE (SUBROUTINE SULIAHARRINGTON_MICRO) 
!     WHICH OPERATES ON A 2D (X-Z) GRID.
!     2D VARIABLES FROM THE MAIN MICROPHYSICS SUBROUTINE ARE THEN REASSIGNED BACK 
!     TO 3D FOR OUTPUT BACK TO DRIVER MODEL USING THIS INTERFACE.
!     MICROPHYSICS TENDENCIES ARE ADDED TO VARIABLES HERE BEFORE BEING PASSED BACK 
!     TO DRIVER MODEL.

!     THIS CODE WAS WRITTEN BY JERRY HARRINGTON (PSU,JYH10@PSU.EDU), HUGH MORRISON (NCAR), 
!     AND KARA SULIA (PSU), AND IMPLEMENTED BY KARA SULIA (KJS5066@GMAIL.COM).
!-------------------------------------------------------------------------------
      SUBROUTINE MP_SULIAHARRINGTON(ID,ITIMESTEP,&
      TH, QV, QC, QR, QI, QS, QG, NC, NR, NI, NS,&
      NG, AI, CI, RHO, PII, P, DT, DZ, HT, W,    &
      ICEDEP, ICESUB, RAINEVAP, SNOWEVAP, &
      SNOWMELT, SNOWDEP, SNOWSUB, SNOWACCR,      &
      CLOUDCOND, CLOUDEVAP, ICEMELT, ICENUC,     &
      RAINFRZ, CLOUDFRZ, SEDR, SEDI, SEDS,    &
      SEDG, NAGGOUT, NRAGGOUT, NSAGGOUT, PRAOUT, &
      CONTACT,PRCOUT, DEP, PRCIOUT,&
      PSACWSOUT, QMULTSOUT,QMULTROUT,PIACROUT,&
      PRACIOUT, PIACRSOUT, PRACISOUT, PRACGOUT,&
      PSACWGOUT,PGSACWOUT,PGRACSOUT, PRDGOUT,&
      EPRDGOUT, EVPMGOUT, PGMLTOUT, QMULTGOUT,QMULTRGOUT,&
      IMMERSION, PRAIOUT, CFRZR               &
      ,RAINACC,RAINACCT,SNOWACC,SNOWACCT &
      ,GRAUPELACC,GRAUPELACCT,FRZNACC,FRZNACCT &
      ,PHI,RHOICE,RELH      &
      ,IDS,IDE, JDS,JDE, KDS,KDE              & ! domain dims
      ,IMS,IME, JMS,JME, KMS,KME              & ! memory dims
      ,ITS,ITE, JTS,JTE, KTS,KTE              & ! tile   dims
      )
!*******************************************************************************

!     QV - water vapor mixing ratio (kg/kg)
!     QC - cloud water mixing ratio (kg/kg)
!     QI - cloud ice mixing ratio (kg/kg)
!     NI - cloud ice number concentration (1/kg)
!     AI - a-axis length mixing ratio (m/kg)
!     CI - c-axis length mixing ratio (m/kg)
!     NOTE: HT NOT USED BY THIS SCHEME AND NOT NEEDED TO BE PASSED INTO SCHEME
!     P - AIR PRESSURE (PA)
!     W - VERTICAL AIR VELOCITY (M/S)
!     TH - POTENTIAL TEMPERATURE (K)
!     PII - exner function - used to convert potential temp to temp
!     DZ - difference in height over interface (m)
!     DT - model time step (sec)
!     ITIMESTEP - time step counter
!     RAINACC - accumulated grid-scale precipitation (mm)
!     RAINACCT - one time step grid scale precipitation (mm/time step)
!     SNOWACC - accumulated grid-scale snow plus cloud ice (mm)
!     SNOWACCT - one time step grid scale snow plus cloud ice (mm/time step)
!     GRAUPELACC - accumulated grid-scale graupel (mm)
!     GRAUPELACCT - one time step grid scale graupel (mm/time step)
!     FRZNACC - accumulated grid-scale snow plus cloud ice plus graupel (mm)
!     FRZNACCT - one time step grid scale snow plus cloud ice plus graupel (mm/time step)

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: ids, ide, jds, jde, kds, kde, &
      ims, ime, jms, jme, kms, kme, &
      its, ite, jts, jte, kts, kte, &
      ITIMESTEP,ID

!     Temporary changed from INOUT to IN

      REAL, DIMENSION(ims:ime, kms:kme, jms:jme), INTENT(INOUT):: &
      QV, QC, QR, QS, QG, QI, NC, NR, NS, NG, NI, AI, CI, TH
!     , effcs, effis

      REAL, DIMENSION(ims:ime, kms:kme, jms:jme), INTENT(IN):: &
      PII, P, DZ, RHO, W

      REAL, DIMENSION(ims:ime, kms:kme, jms:jme), INTENT(INOUT) ::   &
      ICEDEP, ICESUB, RAINEVAP, SNOWEVAP, SNOWMELT,         &
      SNOWDEP, SNOWSUB, SNOWACCR, &
      ICEMELT, ICENUC, RAINFRZ, CLOUDFRZ, PHI, RHOICE, RELH,&
      SEDI, SEDS, SEDR, SEDG, NAGGOUT, NRAGGOUT, NSAGGOUT, PRAOUT, &
      CLOUDCOND, CLOUDEVAP, CONTACT,PRCOUT, DEP, PRCIOUT,&
      PSACWSOUT, QMULTSOUT,QMULTROUT,PIACROUT,&
      PRACIOUT, PIACRSOUT, PRACISOUT, PRACGOUT,&
      PSACWGOUT,PGSACWOUT,PGRACSOUT, PRDGOUT,&
      EPRDGOUT, EVPMGOUT, PGMLTOUT, QMULTGOUT,QMULTRGOUT, IMMERSION, &
      PRAIOUT, CFRZR

      REAL, DIMENSION(ims:ime, jms:jme), INTENT(IN) :: HT
      REAL, DIMENSION(ims:ime, jms:jme), INTENT(INOUT) ::   RAINACC, &
      RAINACCT,SNOWACC, SNOWACCT, GRAUPELACC, GRAUPELACCT, FRZNACC, &
      FRZNACCT

      REAL, INTENT(IN):: DT

      REAL, DIMENSION(its:ite, kts:kte, jts:jte) ::T,THL

      REAL, DIMENSION(its:ite, kts:kte) :: QV2D, QC2D, QR2D, QS2D, &
      QG2D, QI2D, NC2D, NR2D, NS2D, NG2D, NI2D, AI2D, CI2D, T2D,   &
      RHO2D, P2D, ICEDEP2D, ICESUB2D, RAINEVAP2D, SNOWEVAP2D, SNOWMELT2D,&
      SNOWDEP2D, SNOWSUB2D, SNOWACCR2D, CLOUDCOND2D, CLOUDEVAP2D,  &
      ICEMELT2D, ICENUC2D, RAINFRZ2D, CLOUDFRZ2D, PHI2D, RHOICE2D, &
      RELH2D, NAGGOUT2D, NRAGGOUT2D, NSAGGOUT2D, PRAOUT2D, RSED, GSED,&
      ISED, SSED, CONTACT2D, PRCOUT2D, DEP2D, PRCIOUT2D, &
      PSACWSOUT2D,QMULTSOUT2D,QMULTROUT2D,PIACROUT2D, PRACIOUT2D, &
      PIACRSOUT2D, PRACISOUT2D, PRACGOUT2D, PSACWGOUT2D,PGSACWOUT2D,&
      PGRACSOUT2D, PRDGOUT2D, EPRDGOUT2D, EVPMGOUT2D, PGMLTOUT2D,&
      QMULTGOUT2D, QMULTRGOUT2D, IMMERSION2D, PRAIOUT2D, CFRZR2D

      REAL, DIMENSION(its:ite, jts:jte) :: LWP,IWP,PATHVAR,PATHVAR2
      REAL, DIMENSION(kts:kte) ::  avgthl,TOT,TOT2,DZQ

      INTEGER i, j, k, ihr
      REAL, DIMENSION(its:ite) :: PRECPRT, SNOWRT, SNOWPRT, GRPLPRT

      INTEGER, PARAMETER :: nin = 15, nccn = 3, nhrs = 8, & !nhrs=8 for 3-hourly
!                            start_day = 1, end_day=9, size_opt = 8
                            start_day = 11, end_day = 19, size_opt = 8
      REAL IIN(ide,jde,kde,nhrs),CCN(ide,kde,nhrs,8),&
           IIN_SUM(ide,jde,kde,nhrs), IIN_SUMJ(ide,kde,nhrs) !NN(ide-1,jde-1,kde-1,18)


      REAL time_new, ts
      INTEGER, SAVE :: time_save = 0, hr_save = 0
      REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:,:) :: CCN1, CCN2, CCN3
      REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:,:) :: IIN_SUM1,IIN_SUM2,&
      IIN_SUM3
      REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: IIN_SUMJ1,IIN_SUMJ2,&
      IIN_SUMJ3

      ihr = 1
      IF(PIRE_CHEM)THEN
         time_new = real(start_day) + real(itimestep*dt)/86400.
         IF(int(time_new) .gt. int(time_save))THEN
            CALL APM(itimestep,dt,id,ids,ide,jds,jde,kds,kde,nhrs,&
                 IIN,IIN_SUM,CCN,start_day,size_opt,IIN_SUMJ)
            IF(id == 1)THEN
               IF(.not.allocated(CCN1))THEN
                  allocate(CCN1(ide,kde,nhrs,8))
                  allocate(IIN_SUM1(ide,jde,kde,nhrs))
                  allocate(IIN_SUMJ1(ide,kde,nhrs))
               END IF
               CCN1 = CCN
               IIN_SUM1 = IIN_SUM
               IIN_SUMJ1 = IIN_SUMJ
            ELSE IF(id == 2)THEN
               IF(.not.allocated(CCN2))THEN
                  allocate(CCN2(ide,kde,nhrs,8))
                  allocate(IIN_SUM2(ide,jde,kde,nhrs))
                  allocate(IIN_SUMJ2(ide,kde,nhrs))
               END IF
               CCN2 = CCN
               IIN_SUM2 = IIN_SUM
               IIN_SUMJ2 = IIN_SUMJ
            ELSE IF(id == 3)THEN
               IF(.not.allocated(CCN3))THEN
                  allocate(CCN3(ide,kde,nhrs,8))
                  allocate(IIN_SUM3(ide,jde,kde,nhrs))
                  allocate(IIN_SUMJ3(ide,kde,nhrs))
               END IF
               CCN3 = CCN
               IIN_SUM3 = IIN_SUM
               IIN_SUMJ3 = IIN_SUMJ
               time_save = int(time_new)
            END IF
         END IF
         IF(id == 1)THEN
            CCN = CCN1
            IIN_SUM = IIN_SUM1
            IIN_SUMJ = IIN_SUMJ1
         ELSE IF(id == 2)THEN
            CCN = CCN2
            IIN_SUM = IIN_SUM2
            IIN_SUMJ = IIN_SUMJ2
         ELSE IF(id == 3)THEN
            CCN = CCN3
            IIN_SUM = IIN_SUM3
            IIN_SUMJ = IIN_SUMJ3
         END IF
         
         if(nhrs .eq. 8) then
            ts = (time_new - int(time_new))*24.
            ihr = int(ts/3)+1 
         else
            ihr = 1
         end if
      END IF


!     LOGICAL, EXTERNAL :: wrf_dm_on_monitor

      DO i=its,ite
         DO j=jts,jte
            DO k=kts,kte

               T(i,k,j) = TH(i,k,j) * PII(i,k,j)

            END DO
         END DO
      END DO
!     CONVERT 3D VARIABLES TO 2D FOR MICROPHYSICS SCHEME

      DO j=jts,jte              !(east-west)
         DO i=its,ite           !(north-south)
            DO k=kts,kte        !(vertical)

               QV2D(i,k) = QV(i,k,j)
               QC2D(i,k) = QC(i,k,j)
               QR2D(i,k) = QR(i,k,j)
               QI2D(i,k) = QI(i,k,j)
               QS2D(i,k) = QS(i,k,j)
               QG2D(i,k) = QG(i,k,j)

               NC2D(i,k) = NC(i,k,j)
               NR2D(i,k) = NR(i,k,j)
               NI2D(i,k) = NI(i,k,j)
               NS2D(i,k) = NS(i,k,j)
               NG2D(i,k) = NG(i,k,j)
               AI2D(i,k) = AI(i,k,j)
               CI2D(i,k) = CI(i,k,j)

               T2D(i,k) = T(i,k,j)

               ICEDEP2D(i,k) = ICEDEP(i,k,j)
               ICESUB2D(i,k) = ICESUB(i,k,j)
               RAINEVAP2D(i,k) = RAINEVAP(i,k,j)
               SNOWEVAP2D(i,k) = SNOWEVAP(i,k,j)
               SNOWMELT2D(i,k) = SNOWMELT(i,k,j)
               SNOWDEP2D(i,k) = SNOWDEP(i,k,j)
               SNOWSUB2D(i,k) = SNOWSUB(i,k,j)
               SNOWACCR2D(i,k) = SNOWACCR(i,k,j)
               CLOUDCOND2D(i,k) = CLOUDCOND(i,k,j)
               CLOUDEVAP2D(i,k) = CLOUDEVAP(i,k,j)
               ICEMELT2D(i,k) = ICEMELT(i,k,j)
               ICENUC2D(i,k) = ICENUC(i,k,j)
               RAINFRZ2D(i,k) = RAINFRZ(i,k,j)
               CLOUDFRZ2D(i,k) = CLOUDFRZ(i,k,j)
!not sure if the sed and process blocks need to be here
               RSED(i,k) = SEDR(i,k,j)
               GSED(i,k) = SEDG(i,k,j)
               SSED(i,k) = SEDS(i,k,j)
               ISED(i,k) = SEDI(i,k,j)

               NAGGOUT2D(i,k) = NAGGOUT(i,k,j) 
               NRAGGOUT2D(i,k) = NRAGGOUT(i,k,j)
               NSAGGOUT2D(i,k) = NSAGGOUT(i,k,j)
               PRAOUT2D(i,k) = PRAOUT(i,k,j)
               CONTACT2D(i,k) = CONTACT(i,k,j)
               PRCOUT2D(i,k) = PRCOUT(i,k,j)
               DEP2D(i,k) = DEP(i,k,j)
               PRCIOUT2D(i,k) = PRCIOUT(i,k,j) 
               PSACWSOUT2D(i,k) = PSACWSOUT(i,k,j)
               PIACROUT2D(i,k) = PIACROUT(i,k,j)
               PRACIOUT2D(i,k) = PRACIOUT(i,k,j)
               PIACRSOUT2D(i,k) = PIACRSOUT(i,k,j)
               PRACISOUT2D(i,k) = PRACISOUT(i,k,j)
               PRACGOUT2D(i,k) = PRACGOUT(i,k,j)
               PSACWGOUT2D(i,k) = PSACWGOUT(i,k,j)
               PGSACWOUT2D(i,k) = PGSACWOUT(i,k,j)
               PGRACSOUT2D(i,k) = PGRACSOUT(i,k,j)
               PRDGOUT2D(i,k) = PRDGOUT(i,k,j)
               EPRDGOUT2D(i,k) = EPRDGOUT(i,k,j)
               EVPMGOUT2D(i,k) = EVPMGOUT(i,k,j)
               PGMLTOUT2D(i,k) = PGMLTOUT(i,k,j)
               QMULTGOUT2D(i,k) = QMULTGOUT(i,k,j)
               QMULTRGOUT2D(i,k) = QMULTRGOUT(i,k,j)
               IMMERSION2D(i,k) = IMMERSION(i,k,j)
               PRAIOUT2D(i,k) = PRAIOUT(i,k,j)
               CFRZR2D(i,k) = CFRZR(i,k,j)



               P2D(i,k) = P(i,k,j)
               RHO2D(i,k) = RHO(i,k,j)
               DZQ(k) = DZ(i,k,j)

            END DO
         END DO

         !print*,'sulia_harrington'
         CALL SULIAHARRINGTON_MICRO(T2D, DT, QV2D, QC2D, QR2D, QS2D,&
         QG2D, QI2D, NC2D, NR2D, NS2D, NG2D, NI2D, AI2D, CI2D, P2D, &
         RHO2D, DZQ, IMS, IME, JMS, JME, KMS, KME, ITS, ITE, JTS,   &
         JTE, KTS, KTE, IDS, IDE, JDS, JDE, KDS, KDE, &
         ITIMESTEP, ICEDEP2D, ICESUB2D, RAINEVAP2D, SNOWEVAP2D,      &
         SNOWMELT2D, SNOWDEP2D, SNOWSUB2D, SNOWACCR2D, CLOUDCOND2D,  &
         CLOUDEVAP2D, ICEMELT2D, ICENUC2D, RAINFRZ2D, CLOUDFRZ2D,    &
         RSED, ISED, SSED, GSED, NAGGOUT2D, NRAGGOUT2D, NSAGGOUT2D, &
         PRAOUT2D, PHI2D, RHOICE2D, RELH2D, IIN, CCN, IIN_SUM,IIN_SUMJ,&
         PRECPRT, SNOWRT, SNOWPRT, GRPLPRT, IHR, NHRS,&
         CONTACT2D,PRCOUT2D, DEP2D, PRCIOUT2D, PSACWSOUT2D, &
         QMULTSOUT2D,QMULTROUT2D,PIACROUT2D, PRACIOUT2D, PIACRSOUT2D,&
         PRACISOUT2D, PRACGOUT2D, PSACWGOUT2D,PGSACWOUT2D,&
         PGRACSOUT2D, PRDGOUT2D, EPRDGOUT2D, EVPMGOUT2D, PGMLTOUT2D,&
         QMULTGOUT2D, QMULTRGOUT2D, IMMERSION2D, PRAIOUT2D, CFRZR2D)

!     TRANSFER 2D ARRAYS BACK TO 3D FOR WRF

         DO i=its,ite
            DO k=kts,kte

               QV(i,k,j) = QV2D(i,k)
               QC(i,k,j) = QC2D(i,k)
               QI(i,k,j) = QI2D(i,k)
               QR(i,k,j) = QR2D(i,k)
               QS(i,k,j) = QS2D(i,k)
               QG(i,k,j) = QG2D(i,k)

               NC(i,k,j) = NC2D(i,k)
               NR(i,k,j) = NR2D(i,k)
               NI(i,k,j) = NI2D(i,k)
               NS(i,k,j) = NS2D(i,k)
               NG(i,k,j) = NG2D(i,k)
               AI(i,k,j) = AI2D(i,k)
               CI(i,k,j) = CI2D(i,k)
               T(i,k,j) = T2D(i,k)
               TH(i,k,j) = T(i,k,j)/PII(i,k,j)

               ICEDEP(i,k,j) = ICEDEP2D(i,k)
               ICESUB(i,k,j) = ICESUB2D(i,k)
               RAINEVAP(i,k,j) = RAINEVAP2D(i,k)
               SNOWEVAP(i,k,j) = SNOWEVAP2D(i,k)
               SNOWMELT(i,k,j) = SNOWMELT2D(i,k)
               SNOWDEP(i,k,j) = SNOWDEP2D(i,k)
               SNOWSUB(i,k,j) = SNOWSUB2D(i,k)
               SNOWACCR(i,k,j) = SNOWACCR2D(i,k)
               CLOUDCOND(i,k,j) = CLOUDCOND2D(i,k)
               CLOUDEVAP(i,k,j) = CLOUDEVAP2D(i,k)
               ICEMELT(i,k,j) = ICEMELT2D(i,k)
               ICENUC(i,k,j) = ICENUC2D(i,k)
               RAINFRZ(i,k,j) = RAINFRZ2D(i,k)
               CLOUDFRZ(i,k,j) = CLOUDFRZ2D(i,k)

               SEDR(i,k,j) = RSED(i,k)
               SEDG(i,k,j) = GSED(i,k)
               SEDS(i,k,j) = SSED(i,k)
               SEDI(i,k,j) = ISED(i,k)

               NAGGOUT(i,k,j) = NAGGOUT2D(i,k)
               NRAGGOUT(i,k,j) = NRAGGOUT2D(i,k)
               NSAGGOUT(i,k,j) = NSAGGOUT2D(i,k)
               PRAOUT(i,k,j) = PRAOUT2D(i,k)
               CONTACT(i,k,j) = CONTACT2D(i,k)
               PRCOUT(i,k,j) = PRCOUT2D(i,k)
               DEP(i,k,j) = DEP2D(i,k)
               PRCIOUT(i,k,j) = PRCIOUT2D(i,k)
               PSACWSOUT(i,k,j) = PSACWSOUT2D(i,k)
               PIACROUT(i,k,j) = PIACROUT2D(i,k)
               PRACIOUT(i,k,j) = PRACIOUT2D(i,k)
               PIACRSOUT(i,k,j) = PIACRSOUT2D(i,k)
               PRACISOUT(i,k,j) = PRACISOUT2D(i,k)
               PRACGOUT(i,k,j) = PRACGOUT2D(i,k)
               PSACWGOUT(i,k,j) = PSACWGOUT2D(i,k)
               PGSACWOUT(i,k,j) = PGSACWOUT2D(i,k)
               PGRACSOUT(i,k,j) = PGRACSOUT2D(i,k)
               PRDGOUT(i,k,j) = PRDGOUT2D(i,k)
               EPRDGOUT(i,k,j) = EPRDGOUT2D(i,k)
               EVPMGOUT(i,k,j) = EVPMGOUT2D(i,k)
               PGMLTOUT(i,k,j) = PGMLTOUT2D(i,k)
               QMULTGOUT(i,k,j) = QMULTGOUT2D(i,k)
               QMULTRGOUT(i,k,j) = QMULTRGOUT2D(i,k)
               IMMERSION(i,k,j) = IMMERSION2D(i,k)
               PRAIOUT(i,k,j) = PRAIOUT2D(i,k)
               CFRZR(i,k,j) = CFRZR2D(i,k)

               RHOICE(i,k,j) = RHOICE2D(i,k)
               RELH(i,k,j) = RELH2D(i,k)
               PHI(i,k,j) = PHI2D(i,k)
               THL(i,k,j)=TH(i,k,j)*exp((-QC(i,k,j)*2.46E6)/&
               (cp*T(i,k,j)))


            END DO              !end k loop

         RAINACC(i,j) = RAINACC(i,j)+PRECPRT(i)
         if(RAINACC(i,j).gt.1000.)print*,real(itimestep)*dt,i,j,RAINACC(i,j),PRECPRT(i)

         RAINACCT(i,j) = PRECPRT(i)
         SNOWACC(i,j) = SNOWACC(i,j)+SNOWPRT(i)
         SNOWACCT(i,j) = SNOWPRT(i)
         GRAUPELACC(i,j) = GRAUPELACC(i,j)+GRPLPRT(i)
         GRAUPELACCT(i,j) = GRPLPRT(i)
         FRZNACC(i,j) = FRZNACC(i,j) + SNOWRT(i)
         FRZNACCT(i,j) = SNOWRT(i)
         END DO
      END DO

      END SUBROUTINE MP_SULIAHARRINGTON
      
!*****************************************************************************
      SUBROUTINE SULIAHARRINGTON_MICRO(t, dt, qv, qc, qr, qs, qg, qi, nc,  &
      nr, ns, ng, ni, ai, ci, p, rho, dzq, ims, ime, jms, jme, kms, kme,   &
      its, ite, jts, jte, kts, kte, ids, ide, jds, jde, kds, kde,      &
      itimestep,icedep, icesub, rainevap,&
      snowevap, snowmelt, snowdep, snowsub, snowaccr, cloudcond,       &
      cloudevap, icemelt, icenuc, rainfrz, cloudfrz, rsed, ised, ssed, &
      gsed, naggout, nraggout, nsaggout,   &
      praout, phi, rhoice, relh, iin, ccn, iin_sum, iin_sumj, &
      precprt, snowrt, snowprt, grplprt,ihr,nhrs, contact,&
      prcout, dep, prciout, psacwsout, qmultsout, qmultrout, piacrout,&
      praciout, piacrsout, pracisout, pracgout, psacwgout, pgsacwout,&
      pgracsout, prdgout, eprdgout, evpmgout, pgmltout, qmultgout,&
      qmultrgout, immersion, praiout, cfrzr)
!****************************************************************************

      IMPLICIT NONE

      INTEGER, INTENT(IN) ::  &
      ims,ime,jms,jme,kms,kme,&
      its,ite,jts,jte,kts,kte,&
      ids,ide,jds,jde,kds,kde,&
      itimestep, nhrs

      REAL, DIMENSION(its:ite,kts:kte) :: &
      qv, qc, qr, qs, qg, qi, nc, nr, ns, ng, ni, ai, ci, t
 
      REAL, DIMENSION(its:ite,kts:kte) :: &
      icedep, icesub, rainevap, snowevap, snowmelt, snowdep, snowsub,&
      snowaccr, cloudcond, cloudevap, icemelt, icenuc, rainfrz, &
      cloudfrz, dep, contact

      REAL, DIMENSION(kts:kte) :: dzq
      REAL, DIMENSION(its:ite,kts:kte) :: ised, ssed, gsed, rsed

      REAL, DIMENSION(its:ite,kts:kte) :: p, rho, rhoice,phi,relh
      REAL :: dt

      INTEGER i, k, iflag, nstep, n, iaspect, ISDACflag, homofreeze,&
      masssizeflag, ipart, resupply, sphrflag, SEDON, RAINON,       &
      DSTRCHECKS, ICE_CALCS, EVOLVE_ON, snow_on, ice_start_time,    &
      processes, LTRUE, CNUM, redden, nucleation, icontactflag,     &
      ccnflag, graupel, demottflag, ihr

!     process rates: calculates changes due to nucleation,deposition,cond/evap
      REAL mnuccd               !change in qi freezing aerosol
      REAL nnuccd               !change in ni freezing aerosol
      REAL prd, ard, crd        !change in qi, ai, and ci, deposition of ice
      REAL pcc                  !condendation/evaporation cloud droplets
      REAL prc                  !autoconversion of droplets
      REAL pra                  !accretion of droplets by rain
      REAL pre                  !evaporation of rain
      REAL nprci, prci          !change in ni, qi autoconversion cloud ice by snow
      REAL psmlt                !change in q melting snow to rain
      REAL evpms                !change in q melting snow evaporating (snow -->melt-->evap)
      REAL pracs                !change in q rain-snow collection
      REAL pracs1               !change in q rain-snow collection
      REAL npracs1              !change in n rain-snow collection
      REAL prai                 !change in q accretion of cloud ice by snow
      REAL npra                 !change in n due to droplet accretion by rain
      REAL nprai                !change in n accretion of cloud ice
      REAL nprc                 !change in nc autoconversion droplets
      REAL nprc1                !change in nr autoconversion droplets
      REAL nsmlts               !change in n melting snow
      REAL nsmltr               !change in n melting snow to rain
      REAL prds, eprds          !change in q deposition/sublimation snow
      REAL nsubr, nsubs         !loss of nr,ns during evap,sub
      REAL agg, nagg, nsagg, nragg
      REAL nnew, prd1

      REAL pracg                !change in q collection rain by graupel
      REAL pracg1               !change in q collection rain by graupel
      REAL psacwg               !change in q collection droplets by graupel
      REAL pgsacw               !conversion q to graupel due to collection droplets by snow
      REAL pgracs               !conversion q to graupel due to collection rain by snow
      REAL prdg                 !dep of graupel
      REAL eprdg                !sub of graupel
      REAL evpmg                !change q melting of graupel and evap
      REAL pgmlt                !change q melting of graupel
      REAL npracg               !change n collection rain by graupel
      REAL npracg1              !change n collection rain by graupel
      REAL npsacwg              !change n collection droplets by graupel
      REAL nscng                !change n conversion to graupel due to collection droplets by snow
      REAL ngracs               !change n conversion to graupel due to collection rain by snow
      REAL ngmltg               !change n melting graupel
      REAL ngmltr               !change n melting graupel to rain
      REAL nsubg                !change n sub/dep graupel
      REAL psacr                !convertion due to coll of snow by rain
      REAL nmultg               !ice mult due to acc droplets by graupel
      REAL nmultrg              !ice mult due to acc rain by graupel
      REAL qmultg               !change q due to ice mult droplets/graup
      REAL qmultrg              !change q due to ice mult rain/graupel
      REAL qmultr               !change q due to ice rain/snow
      REAL qmults               !change q due to ice mult droplets/snow
      REAL nmults               !change n due to ice mult droplets/snow
      REAL nmultr               !change n due to ice rain/snow
      REAL piacr                !change qr, ice-rain collection
      REAL piacrs               !change qr, ice rain collision, added to snow
      REAL niacr                !change n, ice-rain collection
      REAL praci                !change qi, ice-rain collection
      REAL niacrs               !change n, ice rain collision, added to snow
      REAL pracis               !change qi, ice rain collision, added to snow
      REAL psacws               !change q droplet accretion by snow
      REAL npsacws              !change n droplet accretion by snow
!     ice characteristics
      REAL rhobar
      REAL qitend, nitend
      REAL fmult

      REAL dum, dum1, dum2, temp, theta, celsius, press, rhoa, cpm, qsdum
      REAL evs, evi, qvs, qvi, qvv, xxls, xxlv, xxlf, dqsidt, abi, epss
      REAL epsg
      REAL ani, cni, rni, anf, cnf, rnf, phii, phif
      REAL deltastr, alphstr, alphv, betam
      REAL iwci, iwcf, vi, nidum, voltmp
      REAL weight, igr1, igr2, igr, losstot
      REAL vtbarb, vtbarbm, vtbarblen, vtdumm, vtdumn
      REAL alpham, alpha_v, beta_v, alpha_a, beta_a, lmean, cap
      REAL qt_adv, qt_sed, qi_frac, qiold, niold
      REAL*8 sui, sup, qvqvs, qvqvsi

      REAL lammin, lammax, lamc, pgam

      REAL, DIMENSION(its:ite,kts:kte) :: vtrmi1, vtrni1, vtrli1, &
      effi1, vtrmc

!     sedimenation terms
      REAL, DIMENSION(kts:kte) ::                                      &
      fc, fr, fs, fg, fi, fnc, fnr, fns, fng,fni, fci, fai, dumc, dumr,&
      dums,dumg, dumi, dumnc, dumnr, dumns,dumng, dumni, dumai, dumci, &
      falloutc, falloutr, fallouts, falloutg, fallouti, falloutnc,     &
      falloutnr, falloutns, falloutng, falloutni, falloutai, falloutci,&
      qcsten, qrsten, qssten, qisten, ncsten, nrsten, nssten, nisten,  &
      aisten, cisten, qgsten, ngsten, arn, acn, asn, agn,              &
      qiloss2,qcloss,qiloss,qrloss,qcloss2,qrloss2
      
      REAL, DIMENSION(its:ite,kts:kte) :: naggout, nraggout, &
      nsaggout, praout, prciout, prcout, psacwsout, qmultsout,&
      qmultrout, piacrout, praciout, piacrsout, pracisout,&
      pracgout, psacwgout, pgsacwout, pgracsout, prdgout, eprdgout,&
      evpmgout, pgmltout, qmultgout, qmultrgout, immersion, praiout, &
      cfrzr

      REAL umc, unc, ums, umr, umg, uns, unr, ung
      
      REAL, DIMENSION(kts:kte) :: thetal,dthl

      REAL falltndc, falltndr, falltnds,falltndg,  falltndi, falltndnc, &
      falltndnr, falltndns, falltndng, falltndni, falltndai, falltndci,  &
      rgvm

      REAL, DIMENSION(its:ite) :: precprt, snowrt, snowprt, grplprt

      REAL ncloud, cdist1

!     parameters for rain microphysics

      REAL n0rr, dv, lamr, mu, kap, sc, ab, epsr, dqsdt, ratio, mnuccr
      REAL nnuccr, lams, n0s, lamg, n0g

!     parameters for cloud water processes

      REAL mnuccc, nnuccc, mnucci, nnucci

      REAL alphap,alphamsp,betamp

      REAL gammaFindSumTrip(60)
      DATA gammaFindSumTrip/  0.910547, 0.81807, 0.6874, 0.60127, 1.59767&
      ,2.32423,  2.08818,  1.61921, 1.15865, 0.863071, 0.617586, 0.453917&
      ,0.351975, 0.28794, 0.269298, 0.28794, 0.333623, 0.418883,  0.56992&
      ,0.796458, 1.14325,  1.64103, 1.90138,  1.82653, 1.61921,   1.47436&
      ,1.32463,  1.25556,  1.22239,   1.206,  1.11522,  1.10751,  1.10738&
      ,1.11484,  1.12234,  1.12221, 1.14529,  1.16884,  1.20104,  1.22573&
      ,1.25094,  1.27666,  1.31183,  1.3388,  1.35704,  1.37553,  1.38479&
      ,1.39411,  1.40349,  1.41294, 1.42245,  1.43202,  1.44166,  1.45137&
      ,1.46114,  1.47097,  1.48087, 1.50105,  1.50087,  1.51098/

      INTEGER, PARAMETER :: nin = 15, nccn = 8
      INTEGER size_opt
      REAL :: IIN(ide,jde,kde,nhrs),CCN(ide,kde,nhrs,nccn),&
                          IIN_SUM(ide,jde,kde,nhrs),IIN_SUMJ(ide,kde,nhrs)
      
      REAL rdry(nin),ccnsup(nccn), denom, nuc1, n1,&
           n2, nuc_in
      DATA rdry/1.500E-02, 3.182E-02, 6.364E-02, 1.255E-01,2.291E-01&
      ,3.873E-01, 6.225E-01, 9.843E-01, 1.531E+00, 2.312E+00&
      ,3.269E+00, 5.214E+00, 9.708E+00, 1.632E+01, 2.541E+01/ !microns
      DATA ccnsup/0.001,0.002,0.003,0.004,0.005,0.006,0.007,0.008/

      resupply=0
      masssizeflag=0
      ipart=1                   !0 for RAMS dendrites and needles (Walko et al., 1995)
!     1 for RAMS plates and columns (Walko et al., 1995)
!     2 for Mitchell's mass relations (hex plates, columns)
!     3 for Mitchell's mass relations (stellars, columns)
!     4 for Wood's (2007) mass-size relations
      iaspect    = 0            !set constant aspect ratio (sensitivity study)
      ISDACflag  = 0            !test for ISDAC intercomparison
      sphrflag   = 0            !all ice assumed spheres
      redden     = 0            !reduces density of spheres
      homofreeze = 1            !homogeneous freezing
      snow_on    = 1            !ice --> snow & aggregation
      SEDON      = 1            !sedimentation
      EVOLVE_ON  = 1            !depositional growth
      RAINON     = 2            !0==all rain off; 1==basic processes on; 2==all rain processes
      ICE_CALCS  = 1            !all ice calculations
      ice_start_time = 60.*60.*0.0 !time to begin ice nucleation & homogeneous freezing !LG change to 4
      LTRUE = 0
      graupel = 1               !1 on
      processes = 1
      nucleation = 1 !0 simple, 1 meyers, 2 demott
      demottflag = 1 !0 DeMott 2010, 1 DeMott 2015
      icontactflag = 0 !0 contact freezing off, 1 on
      ccnflag = 0 !0 CCN from aerosol data used, 1 CCN from aerosol data NOT used

      IF(ISDACflag .EQ. 1) homofreeze=0
      IF(ISDACflag .eq. 1) rhoi = 88.4
      IF(sphrflag  .eq. 1) rhoi = 920.
      If(redden .eq. 1) rhoi = 500.
  
      alphap = 0.0
      alphamsp = (4./3.)*pi*rhoi
      betamp = 3.0
      
      qt_adv=0.0
      qt_sed=0.0

      DO i = its,ite
         DO k = kts,kte
            IF(qi(i,k).gt.qsmall.and.ni(i,k).gt.qsmall)THEN
               ni(i,k) = max(ni(i,k),qsmall)
               ai(i,k) = max(ai(i,k),qsmall)
               ci(i,k) = max(ci(i,k),qsmall)

               ani = ((ai(i,k)**2)/(ci(i,k)*nu*ni(i,k)))**(1./3.)
               cni = ((ci(i,k)**2)/(ai(i,k)*nu*ni(i,k)))**(1./3.)

            ELSE 
               qi(i,k) = 0.0
               ni(i,k) = 0.0
               ani = 0.0
               cni = 0.0
            END IF
            rni = (ani**2*cni)**(1./3.)
            ai(i,k) = nu*ni(i,k)*ani
            ci(i,k) = nu*ni(i,k)*cni

            IF(qc(i,k).lt.qsmall)THEN
               qc(i,k) = 0.
               nc(i,k) = 0.
            END IF
            IF(qr(i,k).lt.qsmall)THEN
               qr(i,k) = 0.
               nr(i,k) = 0.
            END IF
            IF(qs(i,k).lt.qsmall)THEN
               qs(i,k) = 0.
               ns(i,k) = 0.
            END IF
            IF(qg(i,k).lt.qsmall)THEN
               qg(i,k) = 0.
               ng(i,k) = 0.
            END IF
!     initialize process rates
            mnuccd = 0.;mnuccr = 0.;mnuccc = 0.;mnucci = 0.
            nnuccd = 0.;nnuccr = 0.;nnuccc = 0.;nnucci = 0.

            pcc = 0.0; prc = 0.0; nprc = 0.0; nprc1 = 0.0
            pra = 0.0; npra = 0.0; nragg = 0.0
            pre = 0.0; pracs = 0.0; pracg = 0.0
            npracg = 0.0; evpms = 0.0; psmlt = 0.0
            evpmg = 0.0; pgmlt = 0.0
            nsubr = 0.0; nsmlts = 0.0; nsmltr = 0.0
            ngmltg = 0.0; ngmltr = 0.0;
            nsagg = 0.0; psacws = 0.0; npsacws = 0.0
            psacwg = 0.0; npsacwg = 0.0
            pracs1 = 0.0; npracs1 = 0.0; psacr = 0.0
            pracg1 = 0.0; npracg1 = 0.0
            qmultr = 0.0; nmultr = 0.0; qmults = 0.0
            nmults = 0.0; qmultg = 0.0; nmultg = 0.0
            qmultrg = 0.0; nmultrg = 0.0
            pgsacw = 0.0; nscng = 0.0
            pgracs = 0.0; ngracs = 0.0
            niacr = 0.0; piacr = 0.0; praci = 0.0
            niacrs = 0.0; piacrs = 0.0; pracis = 0.0
            prd = 0.0; ard = 0.0; crd = 0.0
            prds = 0.0; eprds = 0.0
            prdg = 0.0; eprdg = 0.0
            nsubs = 0.0; nsubg = 0.0
            qitend = 0.0; nitend = 0.0       
            nprci = 0.; nprai = 0.; prai = 0.; prci = 0.
            agg = 0.; nagg = 0.

            vtbarbm = 0.;vtbarb = 0.;vtbarblen = 0.;

            qrsten(k)  = 0.0
            qcsten(k)  = 0.0
            qssten(k) = 0.0
            qgsten(k) = 0.0
            qisten(k) = 0.0
            ncsten(k) = 0.0
            nrsten(k) = 0.0
            nssten(k) = 0.0
            ngsten(k) = 0.0
            nisten(k) = 0.0
            aisten(k) = 0.0
            cisten(k) = 0.0
            qiloss(k) = 0.0
            qrloss(k) = 0.0
            qcloss(k) = 0.0
            qiloss2(k)=0.0
            qcloss2(k)=0.0
            qrloss2(k)=0.0
            phi(i,k)=1.0

            ised(i,k) = 0.0; rsed(i,k) = 0.0; gsed(i,k) = 0.0
            ssed(i,k) = 0.0  

            naggout(i,k) = 0.0; nraggout(i,k) = 0.0
            nsaggout(i,k) = 0.0; praout(i,k) = 0.0

            rainevap(i,k) = 0.0; snowevap(i,k) = 0.0
            snowmelt(i,k) = 0.0; snowdep(i,k) = 0.0; snowsub(i,k) = 0.0
            snowaccr(i,k) = 0.0;

            prcout(i,k) = 0.0; psacwsout(i,k) = 0.0; qmultsout(i,k) = 0.0
            qmultrout(i,k) = 0.0; piacrout(i,k) = 0.0; praciout(i,k) = 0.0
            piacrsout(i,k) = 0.0; pracisout(i,k) = 0.0; pracgout(i,k) = 0.0
            psacwgout(i,k) = 0.0; pgsacwout(i,k) = 0.0; pgracsout(i,k) = 0.0
            prdgout(i,k) = 0.0; eprdgout(i,k) = 0.0; evpmgout(i,k) = 0.0
            pgmltout(i,k) = 0.0; qmultgout(i,k) = 0.0; qmultrgout(i,k) = 0.0
            icenuc(i,k) = 0.0; contact(i,k) = 0.0; dep(i,k) = 0.0 
            immersion(i,k) = 0.0; praiout(i,k) = 0.0; cfrzr(i,k) = 0.0

            thetal(k)=0.0
            dthl(k)=0.0

!     initialize ice characteristics in case first time ice nucleation

            deltastr = 1.
            rhobar = 920.
            IF(ISDACflag .eq. 1) rhobar = 88.4
            If(redden .eq. 1) rhoi = 500.
            temp = t(i,k)
            celsius = temp-273.15
            press = p(i,k)
            rhoa = rho(i,k)
            theta = temp*(100000./press)**(rd/cp)

!     set parameters that vary in space and time
!     add fix for low pressure
            evs = min(FUDGE*press,polysvp1(temp,0)) !Pa
            evi = min(FUDGE*press,polysvp1(temp,1)) !pa
            qvs = 0.622*evs/(press-evs)
            qvi = 0.622*evi/(press-evi)
            qvqvs = qv(i,k)/qvs
            qvqvsi = qv(i,k)/qvi
            qvv = qvs/qvi - 1. 
            sup = qv(i,k)/qvs-1.
            relh(i,k) = qv(i,k)/qvs*100.
            sui = qv(i,k)/qvi-1.
            xxls = 3.15e6-2370.*temp+0.3337e6
            xxlv = 3.1484E6-2370.*temp
            xxlf = xxls-xxlv
            dqsdt = xxlv*qvs/(rv*temp**2)
            dqsidt = xxls*qvi/(rv*temp**2)
            cpm = cp*(1.+0.887*qv(i,k))
            abi = 1.+dqsidt*xxls/cpm
!     diffusivity of water vapor
            dv = 8.794E-5*temp**1.81/press
!     viscocity of air
            mu = (1.496E-6*temp**1.5/(temp+120.))
            kap = 1.414E3*mu

!     schmidt
            sc = mu/(rho(i,k)*dv)

!     parameter for rain evaporation
            ab = 1.+dqsdt*xxlv/cpm

!     "a" parameter for rain fallspeed
            arn(k) = (rhosu/rho(i,k))**0.54*ar
!     "a" parameter for cloud droplet fallspeed, based on Stokes Law
            acn(k) = g*rhow/(18.*mu)
!     "a" parameter for snow fall speed
            asn(k) = (rhosu/rho(i,k))**0.54*as
!     "a" parameter for graupel fall speed
            agn(k) = (rhosu/rho(i,k))**0.54*ag
!     set droplet concentration (units of kg-1)

            IF(PIRE_CHEM)THEN
               IF(sup .gt. ccnsup(1) .and. sup .lt. ccnsup(nccn))THEN
                  
                  weight = (real(int(sup*1000.))+1.0) - (sup*1000.)
                  n1 = int(sup*1000.)
                  n2 = n1+1
                  nc(i,k) = weight*ccn(i,k,ihr,n1) + (1.-weight)*ccn(i,k,ihr,n2)
                  
               ELSE IF(sup .ge. ccnsup(nccn))THEN
                  nc(i,k) = ccn(i,k,ihr,nccn)
               ELSE
                  nc(i,k) = 0.0
               END IF
            ELSE
               nc(i,k)=250.
            END IF
            nc(i,k) = nc(i,k)*1.e6/rho(i,k)
           
            IF(qc(i,k).lt.qsmall.and.qi(i,k).lt.qsmall.and.&
               qs(i,k).lt.qsmall.and.qr(i,k).lt.qsmall.and.&
               qg(i,k).lt.qsmall)THEN
               IF((temp.lt.273.15.and.qvqvsi.lt.FUDGE).or.&
                  (temp.ge.273.15.and.qvqvs.lt.FUDGE))THEN
                  !print*,'there is no condensate'
                  !print*,'qc',qc(i,k),'qi',qi(i,k),'qs',qs(i,k),'qr',&
                  !qr(i,k),temp,qvqvs,qvqvsi
                  GOTO 200
               END IF
            END IF

!     RAIN PROCESSES ----------------------------------------------------------------------
            
!     calculate rain size distribution parameters
            nr(i,k) = max(0.,nr(i,k))
            nc(i,k) = max(0.,nc(i,k))
            ns(i,k) = max(0.,ns(i,k))
            ng(i,k) = max(0.,ng(i,k))
            
            IF(qr(i,k).ge.qsmall)THEN
               lamr = (pi*rhow*nr(i,k)/qr(i,k))**(1./3.)
               n0rr = nr(i,k)*lamr
               
               IF(lamr .lt. lamminr)THEN
                  lamr = lamminr
                  n0rr = lamr**4*qr(i,k)/(pi*rhow)
                  nr(i,k) = n0rr/lamr
               ELSE IF(lamr.gt.lammaxr)THEN
                  lamr = lammaxr
                  n0rr = lamr**4*qr(i,k)/(pi*rhow)
                  nr(i,k) = n0rr/lamr
               END IF
            END IF

!     CLOUD DROPLETS ---------------------------------------------------------------------
!     Martin et al. (1994) formula for pgam
            pgam = 0.0
            lamc = 0.0
            CDIST1 = 0.0
            IF(qc(i,k).gt.qsmall)THEN
               dum = press/(287.15*temp)
               pgam = 0.0005714*(nc(i,k)/1.e6*dum)+0.2714
               pgam = 1./(pgam**2)-1.
               pgam = max(pgam,2.)
               pgam = min(pgam,10.)
               
               lamc = (cons26*nc(i,k)*GAMMA(pgam+4.)/&
               (qc(i,k)*GAMMA(pgam+1.)))**(1./3.)
               
               lammin = (pgam+1.)/60.e-6
               lammax = (pgam+1.)/1.e-6
               
               IF(lamc.lt.lammin)THEN
                  lamc = lammin
                  nc(i,k) = exp(3.*log(lamc)+log(qc(i,k))+&
                  log(GAMMA(pgam+1.))-log(GAMMA(pgam+4.)))/cons26
               ELSE IF(lamc.gt.lammax)THEN
                  lamc = lammax
                  nc(i,k) = exp(3.*log(lamc)+log(qc(i,k))+&
                  log(GAMMA(pgam+1.))-log(GAMMA(pgam+4.)))/cons26
               END IF
               CDIST1 = nc(i,k)/gamma(pgam+1.)
            END IF              !qc >=qsmall
!     SNOW---------------------------------------------------------------------------------
            IF(qs(i,k).ge.qsmall)THEN
               lams = (cons1*ns(i,k)/qs(i,k))**(1./DS)
               n0s = ns(i,k)*lams
                  
               IF(lams.lt.lammins)THEN
                  lams = lammins
                  n0s = lams**4*qs(i,k)/cons1
                  ns(i,k) = n0s/lams
               ELSE IF(lams.gt.lammaxs)THEN
                  lams = lammaxs
                  n0s = lams**4*qs(i,k)/cons1
                  ns(i,k) = n0s/lams
               END IF
               qsdum = 2.*pi*n0s*rho(i,k)*(f1s/(lams*lams)+&
               f2s*SQRT(asn(k)*rho(i,k)/mu)*sc**(1./3.)*cons10/&
               (lams**cons35))
            END IF           !qs >= qsmall

            if(qg(i,k).ge.qsmall)then
                lamg = (cons2*ng(i,k)/qg(i,k))**(1./DG)
                n0g = ng(i,k)*lamg
                if(lamg.lt.lamming)then
                   lamg = lamming
                   n0g = lamg**4*qg(i,k)/cons2
                   ng(i,k) = n0g/lamg
                else if(lamg.gt.lammaxg)then
                   lamg = lammaxg
                   n0g = lamg**4*qg(i,k)/cons2
                   ng(i,k) = n0g/lamg
                end if
            end if
            IF(processes .eq. 1)THEN
   
            IF(temp.ge.273.15)THEN


!       melt all very small snow and graupel

               IF(RAINON.ge.1)THEN
               IF(qs(i,k).lt.1.e-6)THEN
                  qr(i,k) = qr(i,k)+qs(i,k)
                  nr(i,k) = nr(i,k)+ns(i,k)
                  temp=temp-qs(i,k)*xxlf/cpm
                  qs(i,k) = 0.0
                  ns(i,k) = 0.0
               ENDIF

               IF(qg(i,k).lt.1.e-6)THEN
                  qr(i,k) = qr(i,k)+qg(i,k)
                  nr(i,k) = nr(i,k)+ng(i,k)
                  temp=temp-qg(i,k)*xxlf/cpm
                  qg(i,k) = 0.0
                  ng(i,k) = 0.0
               ENDIF
               ENDIF
!       recalculate rain distribution parameters if rain now exists
               IF(qr(i,k).ge.qsmall)THEN
                  lamr = (pi*rhow*nr(i,k)/qr(i,k))**(1./3.)
                  n0rr = nr(i,k)*lamr

                  IF(lamr .lt. lamminr)THEN
                     lamr = lamminr
                     n0rr = lamr**4*qr(i,k)/(pi*rhow)
                     nr(i,k) = n0rr/lamr
                  ELSE IF(lamr.gt.lammaxr)THEN
                     lamr = lammaxr
                     n0rr = lamr**4*qr(i,k)/(pi*rhow)
                     nr(i,k) = n0rr/lamr
                  END IF
               END IF

!     autoconversion of cloud liquid water to rain
               IF(RAINON.ge.1)THEN
               IF(qc(i,k).ge.1.E-6)THEN
!     formula from khairoutdinov and kogan 2000, mwr
                  prc=1350.*qc(i,k)**2.47*(nc(i,k)/1.e6*&
                  rho(i,k))**(-1.79)
              
                  nprc1 = prc/cons29
                  nprc = prc/qc(i,k)*nc(i,k)
                  
                  nprc = min(nprc,nc(i,k)/dt)
                  nprc1 = min(nprc1,nprc)
               END IF
               

!     accretion of cloud water by rain
!     continuous collection eq with grav. coll. kernel, droplet fall spd neglected

               IF(qr(i,k).ge.1.E-8 .and. qc(i,k).ge.1.E-8)THEN
                  dum=(qc(i,k)*qr(i,k))
                  pra = 67.*dum**1.15
                  npra = pra/qc(i,k)*nc(i,k)
               END IF
               praout(i,k) = pra
!     self-collection of rain drops
               IF(qr(i,k).ge.1.e-8)THEN
                  dum1 = 300.e-6
                  IF(1./lamr.lt.dum1)THEN
                     dum = 1.
                  ELSE IF(1./lamr.ge.dum1)THEN
                     dum = 2.-exp(2300.*(1./lamr-dum1))
                  END IF
                  nragg = -5.78*dum*nr(i,k)*qr(i,k)*rho(i,k)
               END IF

!     evaporation of rain (RUTLEDGE AND HOBBS 1983)
               
               IF(qr(i,k).ge.qsmall)THEN
                  epsr = 2.*pi*n0rr*rho(i,k)*dv*(f1r/(lamr*lamr)+f2r*&
                  SQRT(arn(k)*rho(i,k)/mu)*sc**(1./3.)*cons9/&
                  (lamr**cons34))
               ELSE
                  epsr = 0.
               END IF
!     no condensation onto rain, only evap allowed
               
               IF (qv(i,k).lt.qvs) THEN
                  pre = min((epsr*(qv(i,k)-qvs)/ab),0.)
               ELSE
                  pre = 0.
               END IF        
               END IF !RAINON>=1
            !     collection of snow by rain above freezing
!     formula from Ikawa and Saito (1991)

               IF(RAINON.eq.2)THEN
               IF(qr(i,k).ge.1.e-8.and.qs(i,k).ge.1.e-8)THEN
                  ums = asn(k)*cons3/(lams**bs)
                  umr = arn(k)*cons4/lamr**br
                  uns = asn(k)*cons5/lams**bs
                  unr = arn(k)*cons6/lamr**br
!     set realistic limits on fallspeeds
                     
                  dum = (rhosu/rho(i,k))**0.54
                  ums = min(ums,1.2*dum)
                  uns = min(uns,1.2*dum)
                  umr = min(umr,9.1*dum)
                  unr = min(unr,9.1*dum)
                     
                  pracs = cons41*(SQRT((1.2*umr-0.95*ums)**2 + 0.08*&
                  ums*umr)*rho(i,k)*n0rr*n0s/lamr**3* &
                  (5./(lamr**3*lams)+2./(lamr*lamr*lams*lams)+0.5/&
                  (lamr*lams**3)))
               END IF        !qr & qs >=1.e-8
!     collection of graupel by rain above freezing
               IF(graupel .eq. 1)then
                  IF(qr(i,k).ge.1.e-8.and.qg(i,k).ge.1.e-8)then
                     umg = agn(k)*cons7/(lamg**BG)
                     umr = arn(k)*cons4/(lamr**BR)
                     ung = agn(k)*cons8/(lamg**BG)
                     unr = arn(k)*cons6/(lamr**BR)
                     
                     dum = (rhosu/rho(i,k))**0.54
                     umg = min(umg,20.*dum)
                     ung = min(ung,20.*dum)
                     umr = min(umr,9.1*dum)
                     unr = min(unr,9.1*dum)

!       pracg is mixing ratio of rain per sec collected by graupel
                     pracg = cons41*(((1.2*umr-0.95*umg)**2+0.08*umg*&
                        umr)**0.5*rho(i,k)*n0rr*n0g/lamr**3*(5./&
                        (lamr**3*lamg)+2./(lamr**2*lamg**2)+0.5/&
                        (lamr*lamg**3)))
!       assume 1 mm drops are shed, get number shed per sec
                     npracg = cons32*rho(i,k)*(1.7*(unr-ung)**2+0.3*&
                        unr*ung)**0.5*n0rr*n0g*(1./(lamr**3*lamg)+&
                        1./(lamr**2*lamg**2)+1./(lamr*lamg**3))

                     dum = pracg/5.2E-7
                     npracg = npracg-dum

                  END IF
               END IF
!     melting of snow
!     snow may persist above freezing, from Rutledge and Hobbs (1984)
!     if supersat, snow melts to form rain
               IF(qs(i,k).ge.1.e-8)THEN
                  dum = -cpw/xxlf*(temp-273.15)*pracs
                  psmlt = qsdum*kap*(273.15-temp)/xxlf+dum
                  IF(qvqvs.lt.1.)THEN
                     epss = qsdum*dv
                     evpms = (qv(i,k)-qvs)*epss/ab
                     evpms = max(evpms,psmlt)
                     psmlt = psmlt-evpms
                  END IF
               END IF        !qs>=1.e-8
               pracs = 0.

!       melting of graupel
!       graupel may persist above freezing, from Rutledge and Hobbs
!       (1984)
!       if supersat, graupel melts to form rain
               if(graupel.eq.1)then
                  if(qg(i,k).ge.1.e-8)then
                     dum = -cpw/xxlf*(temp-273.15)*pracg
                     pgmlt = 2.*pi*n0g*kap*(273.15-temp)/xxlf*(f1s/&
                        (lamg*lamg)+f2s*(agn(k)*rho(i,k)/mu)**0.5*&
                        sc**(1./3.)*cons11/(lamg**cons36))+dum   !melting of graupel to form rain
                     if(qvqvs.lt.1)then !if subsat, graupel sublimates
                        epsg = 2.*pi*n0g*rho(i,k)*dv*(f1s/(lamg*lamg)+&
                          f2s*(agn(k)*rho(i,k)/mu)**0.5*sc**(1./3.)*&
                          cons11/(lamg**cons36))
                        evpmg = (qv(i,k)-qvs)*epsg/ab
                        evpmg = max(evpmg,pgmlt) !graupel sublimating
                        pgmlt = pgmlt-evpmg !remainder of graupel melts to rain
                     end if
                  end if
                  pracg = 0.
               end if
               END IF!RAINON==2
!     check cloud water conservation, update qc and qr with process rates
               
               dum = (prc+pra)*dt
               IF(dum.gt.qc(i,k).and.qc(i,k).ge.qsmall)THEN
                  ratio = qc(i,k)/dum
                  prc = prc*ratio
                  pra = pra*ratio
               END IF
            
!     conservation of qr, NOTE: pre is a negative number

               dum = (-pracs-pracg-pre-pra-prc+psmlt+pgmlt)*dt
               IF(dum.gt.qr(i,k).and.qr(i,k).ge.qsmall) THEN
                  ratio = (qr(i,k)/dt+pra+prc+pracs+pracg-psmlt-pgmlt)/(-pre)
                  pre = pre*ratio
               END IF

!     conservation of qs
               dum = (pracs-evpms-psmlt)*dt !melting, evap, & accretion of snow
               IF(dum.gt.qs(i,k).and.qs(i,k).ge.qsmall)THEN
                  ratio = qs(i,k)/dum
                  pracs = pracs*ratio
                  psmlt = psmlt*ratio
                  evpms = evpms*ratio
               END IF

!     conservation of qg
               dum = (-pgmlt-evpmg-pracg)*dt
               if (dum.gt.qg(i,k).and.qg(i,k).ge.qsmall)then
                  ratio = qg(i,k)/dum
                  pgmlt = pgmlt*ratio
                  evpmg = evpmg*ratio
                  pracg = pracg*ratio
               end if

               IF(pre.lt.0.0)THEN
                  dum = pre*dt/qr(i,k)
                  dum = max(-1.,dum)
                  nsubr = dum*nr(i,k)/dt
               END IF
               IF(evpms+psmlt.lt.0.0)THEN
                  dum = (evpms+psmlt)*dt/qs(i,k)
                  dum = max(-1.,dum)
                  nsmlts = dum*ns(i,k)/dt
               END IF
               IF(psmlt.lt.0.0)THEN
                  dum = (psmlt)*dt/qs(i,k)
                  dum = max(-1.,dum)
                  nsmltr = dum*ns(i,k)/dt
               END IF
               IF(evpmg+pgmlt.lt.0.0)THEN
                  dum = (evpmg+pgmlt)*dt/qg(i,k)
                  dum = max(-1.,dum)
                  ngmltg = dum*ng(i,k)/dt
               END IF
               IF(pgmlt.lt.0.0)THEN
                  dum = (pgmlt)*dt/qg(i,k)
                  dum = max(-1.,dum)
                  ngmltr = dum*ng(i,k)/dt
               END IF


            ELSE !temp
!     autoconversion of cloud liquid water to rain
               IF(RAINON.ge.1)THEN
               IF(qc(i,k).ge.1.E-6)THEN
!     formula from khairoutdinov and kogan 2000, mwr
                  prc=1350.*qc(i,k)**2.47*(nc(i,k)/1.e6*&
                  rho(i,k))**(-1.79)

                  nprc1 = prc/cons29
                  nprc = prc/(qc(i,k)/nc(i,k))
             
                  nprc = min(nprc,nc(i,k)/dt)
                  nprc1 = min(nprc1,nprc)
               END IF
               END IF !RAINON>=1
!     aggregation of qs
               IF(qs(i,k).ge.1.e-8)THEN
                  nsagg = cons15*asn(k)*(qs(i,k)*rho(i,k))**&
                  ((2.+BS)/3.)*(ns(i,k)*rho(i,k))**((4.-BS)/3.)&
                  /rho(i,k)
               END IF

!     accretion of cloud droples onto snow/graupel
!     use continuous growth equations with
!     simple gravitational collection kernel ignoring snow
               if(graupel .eq. 1)then
                  if(qs(i,k) .ge.1.e-8.and.qc(i,k).ge.qsmall)then
                     psacws = cons13*asn(k)*qc(i,k)*rho(i,k)*n0s/&
                         lams**(BS+3.)
                     
                     npsacws = cons13*asn(k)*nc(i,k)*rho(i,k)*n0s/&
                         lams**(BS+3.)
                  end if
!     collection of cloud water by graupel
                  if(qg(i,k).ge.1.e-8.and.qc(i,k).ge.qsmall)then
                     psacwg = cons14*agn(k)*qc(i,k)*rho(i,k)*n0g/&
                         lamg**(BG+3.)
                  
                     npsacwg = cons14*agn(k)*nc(i,k)*rho(i,k)*n0g/&
                         lamg**(BG+3.)
                  end if
               end if
               
               IF(RAINON.eq.2)THEN
!     accretion of rain water by snow
               IF(qr(i,k).ge.1.e-8.and.qs(i,k).ge.1.e-8)THEN
                  ums = asn(k)*cons3/(lams**bs)
                  umr = arn(k)*cons4/lamr**br
                  uns = asn(k)*cons5/lams**bs
                  unr = arn(k)*cons6/lamr**br
!     set realistic limits on fallspeeds
                     
                  dum = (rhosu/rho(i,k))**0.54
                  ums = min(ums,1.2*dum)
                  uns = min(uns,1.2*dum)
                  umr = min(umr,9.1*dum)
                  unr = min(unr,9.1*dum)
                     
                  pracs1 = cons41*(SQRT((1.2*umr-0.95*ums)**2 + 0.08*&
                  ums*umr)*rho(i,k)*n0rr*n0s/lamr**3* &
                  (5./(lamr**3*lams)+2./(lamr*lamr*lams*lams)+0.5/&
                  (lamr*lams**3)))

                  npracs1 = cons32*rho(i,k)*(1.7*(unr-uns)**2 + 0.3*&
                  unr*uns)**0.5*n0rr*n0s*(1./(lamr**3*lams) + 1./&
                  (lamr**2*lams**2)+1./(lamr*lams**3))

                  pracs1 = min(pracs1,qr(i,k)/dt)

!     collection of snow by rain -- needed for graupel conversion calcs
!     only calculate if snow and rain mixing ratios exceed 0.1 g/kg

                  if(qs(i,k).ge.0.1e-3.and.qr(i,k).ge.0.1e-3)then
                     psacr = cons31*(((1.2*umr-0.95*ums)**2+0.08*ums*&
                        umr)**0.5*rho(i,k)*n0rr*n0s/lams**3*(5./(lams**3*&
                        lamr)+2./(lams**2*lamr**2)+0.5/(lams*lamr**3)))
                  end if

               END IF

!     collection of rainwater by graupel, Ikawa and Saito 1990
               IF(graupel.eq.1)then
                  IF(qr(i,k).ge.1.e-8.and.qg(i,k).ge.1.e-8)then
                     umg = agn(k)*cons7/(lamg**BG)
                     umr = arn(k)*cons4/(lamr**BR)
                     ung = agn(k)*cons8/(lamg**BG)
                     unr = arn(k)*cons6/(lamr**BR)
                    
                     dum = (rhosu/rho(i,k))**0.54
                     umg = min(umg,20.*dum)
                     ung = min(ung,20.*dum)
                     umr = min(umr,9.1*dum)
                     unr = min(unr,9.1*dum)
                     
!     pracg is mixing ratio of rain per sec collected by graupel
                     pracg1 = cons41*(((1.2*umr-0.95*umg)**2+0.08*umg*&
                        umr)**0.5*rho(i,k)*n0rr*n0g/lamr**3*(5./&
                        (lamr**3*lamg)+2./(lamr**2*lamg**2)+0.5/&
                        (lamr*lamg**3)))
!     assume 1 mm drops are shed, get number shed per sec
                     npracg1 = cons32*rho(i,k)*(1.7*(unr-ung)**2+0.3*&
                        unr*ung)**0.5*n0rr*n0g*(1./(lamr**3*lamg)+&
                        1./(lamr**2*lamg**2)+1./(lamr*lamg**3))
                     
                     pracg1 = min(pracg1,qr(i,k)/dt)
                     
                  END IF
               END IF !graupel
               END IF!RAINON==2

!     rime-splintering of snow
!     hallet-mossop 1974
!     number of splinters formed is based on mass of rimed water
               IF (graupel .eq. 1)THEN
                  if(qs(i,k).ge.0.1e-3)then
                     if(qc(i,k).ge.0.5e-3.or.qr(i,k).ge.0.1e-3)then
                        if(psacws.gt.0.0.or.pracs1.gt.0.0)then
                           if(temp.le.270.16.and.temp.gt.268.16)then
                              fmult = (270.16-temp)/2.
                           else if(temp.ge.265.16.and.temp.le.268.16)then
                              fmult = (temp-265.16)/3.
                           else
                              fmult = 0.0
                           end if
                        end if
!     splintering from droplets accreted onto snow
                        if(psacws.gt.0.0)then
                           nmults = 35.e4*psacws*fmult*1000.
                           qmults = nmults*mmult
!     constrain so that transfer of mass from snow to ice cannot be more mass
!     than was rimed onto snow
                           qmults = min(qmults,psacws)
                           psacws = psacws-qmults
                        end if
                   
                        if(pracs1.gt.0.)then
                           nmultr = 35.e4*pracs1*fmult*1000.
                           qmultr = nmultr*mmult
                           
                           qmultr = min(qmultr,pracs1)
                           pracs1 = pracs1-qmultr
                        end if
                     end if
                  end if
              
!     rime-splintering of graupel
                  if(qg(i,k).ge.0.1e-3)then
                     if(qc(i,k).ge.0.5e-3.or.qr(i,k).ge.0.1e-3)then
                        if(psacwg.gt.0.0.or.pracg1.gt.0.0)then
                           if(temp.le.270.16.and.temp.gt.268.16)then
                              fmult = (270.16-temp)/2.
                           else if(temp.ge.265.16.and.temp.le.268.16)then
                              fmult = (temp-265.16)/3.
                           else
                              fmult = 0.
                           end if
!     splintering from droplets accreted onto graupel
                           if(psacwg.gt.0.0)then
                              nmultg = 35.e4*psacwg*fmult*1000.
                              qmultg = nmultg*mmult
!     constrain so that transfer of mass from snow to ice cannot be more mass
!     than was rimed onto snow
                              qmultg = min(qmultg,psacwg)
                              psacwg = psacwg-qmultg
                           end if
                           
                           if(pracg1.gt.0.)then
                              nmultrg = 35.e4*pracg1*fmult*1000.
                              qmultrg = nmultrg*mmult
                              
                              qmultrg = min(qmultrg,pracg1)
                              pracg1 = pracg1-qmultrg
                           end if
                        end if
                     end if
                  end if

!     conversion of rimed cloud water onto snow to graupel
                  if(psacws.gt.0.0)then
                     if(qs(i,k).ge.0.1e-3.and.qc(i,k).ge.0.5e-3)then
!     portion of riming converted to graupel
                        pgsacw = min(psacws,cons17*dt*n0s*qc(i,k)*&
                           qc(i,k)*asn(k)*asn(k)/(rho(i,k)*lams**&
                           (2.*BS+2.)))
                        
                        dum = max(rhosn/(rhog-rhosn)*pgsacw,0.)
                        nscng = dum/mg0*rho(i,k)
                        nscng = min(nscng,ns(i,k)/dt)
!     portion of riming left for snow
                        psacws = psacws - pgsacw 
                     end if
                  end if
                  
!     conversion of rimed rain water onto snow converted to graupel
                  if(pracs1.gt.0.0)then
                     if(qs(i,k).ge.0.1e-3.and.qr(i,k).ge.0.1e-3)then
                        dum = cons18*(4./lams)**3*(4./lams)**3 &
                           /(cons18*(4./lams)**3*(4./lams)**3+&
                           cons19*(4./lamr)**3*(4./lamr)**3)
                        dum = min(dum,1.)
                        dum = max(dum,0.)
                        pgracs = (1.-dum)*pracs1
                        ngracs = (1.-dum)*npracs1
                        ngracs = min(ngracs,nr(i,k)/dt)
                        ngracs = min(ngracs,ns(i,k)/dt)
                        
                        pracs1 = pracs1 - pgracs
                        npracs1 = npracs1 - ngracs
                        
                        psacr = psacr*(1.-dum)
                     end if
                  end if


               END IF           !graupel

!     accretion of cloud water by rain
!     continuous collection eq with grav. coll. kernel, droplet fall spd neglected
               IF(RAINON.ge.1)THEN
               IF(qr(i,k).ge.1.E-8 .and. qc(i,k).ge.1.E-8)THEN
                  dum=(qc(i,k)*qr(i,k))
                  pra = 67.*dum**1.15
                  npra = pra/qc(i,k)*nc(i,k)
               END IF
               
!     self-collection of rain drops
               IF(qr(i,k).ge.1.e-8)THEN
                  dum1 = 300.e-6
                  IF(1./lamr.lt.dum1)THEN
                     dum = 1.
                  ELSE IF(1./lamr.ge.dum1)THEN
                     dum = 2.-exp(2300.*(1./lamr-dum1))
                  END IF
                  nragg = -5.78*dum*nr(i,k)*qr(i,k)*rho(i,k)
               END IF
!     evaporation of rain (RUTLEDGE AND HOBBS 1983)
               
               IF(qr(i,k).ge.qsmall)THEN
                  epsr = 2.*pi*n0rr*rho(i,k)*dv*(f1r/(lamr*lamr)+f2r*&
                  SQRT(arn(k)*rho(i,k)/mu)*sc**(1./3.)*cons9/&
                  (lamr**cons34))
               ELSE
                  epsr = 0.
               END IF
               
!     no condensation onto rain, only evap allowed
               
               IF (qv(i,k).lt.qvs) THEN
                  pre = min((epsr*(qv(i,k)-qvs)/ab),0.)
               ELSE
                  pre = 0.
               END IF    
               END IF!RAINON>=1
!     collision of rain and ice to produce snow or graupel
               IF(RAINON.eq.2)THEN
               IF(graupel.eq.1)then
                  if(qr(i,k).ge.1.e-8.and.qi(i,k).ge.1.e-8)then
                     if(qr(i,k).ge.0.1e-3)then
                        niacr = cons24*ni(i,k)*n0rr*arn(k)/lamr**&
                           (BR+3.)*rho(i,k)
                        piacr = cons25*ni(i,k)*n0rr*arn(k)/lamr**&
                           (BR+3.)/lamr**3*rho(i,k)
                        praci = cons24*qi(i,k)*n0rr*arn(k)/lamr**&
                           (BR+3.)*rho(i,k)
                        niacr = min(niacr,nr(i,k)/dt)
                        niacr = min(niacr,ni(i,k)/dt)
                     else
                        niacrs = cons24*ni(i,k)*n0rr*arn(k)/lamr**&
                           (BR+3.)*rho(i,k)
                        piacrs = cons25*ni(i,k)*n0rr*arn(k)/lamr**&
                           (BR+3.)/lamr**3*rho(i,k)
                        pracis = cons24*qi(i,k)*n0rr*arn(k)/lamr**&
                           (BR+3.)*rho(i,k)
                        niacrs = min(niacrs,nr(i,k)/dt)
                        niacrs = min(niacrs,ni(i,k)/dt)
                     end if
                  end if
               END IF
    
!     deposition of qs
               IF(qs(i,k).ge.qsmall)THEN
                  epss = 2.*pi*n0s*rho(i,k)*dv*(f1s/(lams*lams)+&
                     f2s*(asn(k)*rho(i,k)/mu)**0.5*sc**(1./3.)*cons10/&
                     (lams**cons35))
               ELSE
                  epss = 0.0
               END IF      
               prds = epss*(qv(i,k)-qvi)/abi
               END IF!RAINON==2
!     deposition of graupel
               IF(graupel.eq.1)then
                  if(qg(i,k).ge.qsmall)then
                     epsg = 2.*pi*n0g*rho(i,k)*dv*(f1s/(lamg*lamg)+&
                        f2s*(agn(k)*rho(i,k)/mu)**0.5*sc**(1./3.)*cons11/&
                        (lamg**cons36))
                  else
                     epsg = 0.0
                  end if
                  prdg = epsg*(qv(i,k)-qvi)/abi
               END IF!graupel
               
               dum = (qv(i,k)-qvi)/dt
               IF((dum.gt.0..and.prds+prdg.gt.dum*FUDGE).or.&
                  (dum.lt.0..and.prds+prdg.lt.dum*FUDGE))THEN
                  prds = FUDGE*prds*dum/(prds+prdg)
                  prdg = FUDGE*prdg*dum/(prds+prdg)
               END IF
               IF(prds.lt.0.)THEN
                  eprds = prds
                  prds = 0.
               END IF
               IF(prdg.lt.0.0)then
                  eprdg = prdg
                  prdg = 0.0
               END IF
               
!     check cloud water conservation, update qc and qr with process rates
               dum = (prc+pra+psacws+qmults+psacwg+pgsacw+qmultg)*dt
               IF(dum.gt.qc(i,k).and.qc(i,k).ge.qsmall)THEN
                  ratio = qc(i,k)/dum
                  prc = prc*ratio
                  pra = pra*ratio
                  psacws = psacws*ratio
                  qmults = qmults*ratio
                  qmultg = qmultg*ratio
                  psacwg = psacwg*ratio
                  pgsacw = pgsacw*ratio
               END IF
               
!     conservation of qr, NOTE: pre is a negative number
               
               dum = ((pracs1-pre)+(qmultr+qmultrg-prc)+(-pra)+&
                  piacr+piacrs+pgracs+pracg1)*dt
               IF(dum.gt.qr(i,k).and.qr(i,k).ge.qsmall) THEN
                  ratio = (qr(i,k)/dt+pra+prc)/(-pre+qmultr+qmultrg+&
                     pracs1+piacr+piacrs+pgracs+pracg1)
                  pre = pre*ratio   
                  pracs1 = pracs1*ratio
                  qmultr = qmultr*ratio
                  qmultrg = qmultrg*ratio
                  piacr = piacr*ratio
                  piacrs = piacrs*ratio
                  pgracs = pgracs*ratio
                  pracg1 = pracg1*ratio
               END IF

               IF(pre.lt.0.)THEN
                  dum = pre*dt/qr(i,k)
                  dum = max(-1.,dum)
                  nsubr = dum*nr(i,k)/dt
               END IF
               
!     conservation of qs
               dum = (-prds-psacws-pracs1-eprds+psacr-piacrs-pracis)*dt 
               IF(dum.gt.qs(i,k).and.qs(i,k).ge.qsmall)THEN
                  ratio = (qs(i,k)/dt+prds+psacws+pracs1+piacrs+pracis)&
                     /(-eprds+psacr)
                  eprds = eprds*ratio
                  psacr = psacr*ratio
               END IF

               IF(eprds.lt.0.)THEN
                  dum = eprds*dt/qs(i,k)
                  dum = max(-1.,dum)
                  nsubs = dum*ns(i,k)/dt
               END IF

!     conservation of qg
               dum = (-psacwg-pracg1-pgsacw-pgracs-prdg-eprdg-&
                  piacr-praci-psacr)*dt
               if(dum.gt.qg(i,k).and.qg(i,k).ge.qsmall)then
                  ratio = (qg(i,k)/dt+psacwg+pracg1+pgsacw+pgracs+prdg+&
                     psacr+piacr+praci)/(-eprdg)

                  eprdg = eprdg*ratio
               end if


!     sublimate, melt, or evaporate number concentration

               if(eprdg.lt.0.0)then
                  dum = eprdg*dt/qg(i,k)
                  dum = max(-1.,dum)
                  nsubg = dum*ng(i,k)/dt
               end if

            END IF              !temp<273.15
            
            END IF

!     IF(ICE_CALCS .eq. 1)THEN
            IF(qi(i,k).gt.qsmall.and.ni(i,k).gt.qsmall)THEN
               
!     set minimum values for ni,ai,ci, otherwise
               
               ni(i,k) = max(ni(i,k),qsmall)
               ai(i,k) = max(ai(i,k),qsmall)
               ci(i,k) = max(ci(i,k),qsmall)
               
!     get characteristic ci,ai (i.e. cni,ani) from ci,ai
               
               ani = ai(i,k)/(nu*ni(i,k))
               cni = ci(i,k)/(nu*ni(i,k))
               
!     get deltastr from cni and ani
!     deltastr = 1 for ice particles pre-diagnosed as spheres
               
               IF((log(ani)-log(ao)).gt. 0.01 &
                  .and.(log(cni)-log(ao)).gt.0.001)THEN
                  deltastr = (log(cni)-log(ao))/(log(ani)-log(ao))
               ELSE
                  deltastr = 1.
               ENDIF
               IF(ISDACflag .eq. 1)deltastr = 1.0
               IF(iaspect .eq. 1) deltastr = 0.8
               IF(masssizeflag .eq. 1) deltastr=1.0
        IF(sphrflag .eq. 1) deltastr = 1.0  
!     make sure deltastr within reasonable limits
!     IF(deltastr.lt.0.7)THEN
               
!     deltastr = 0.7
!     cni=ao**(1.-deltastr)*ani**deltastr
!     ci(i,k)=nu*ni(i,k)*cni
               
!     ELSE IF (deltastr.gt.1.5)THEN
               
!     deltastr=1.5
!     cni=ao**(1.-deltastr)*ani**deltastr
!     ci(i,k)=nu*ni(i,k)*cni
               
!     END IF
               
               if(deltastr.lt.0.55) then
                  voltmp=(4./3.)*pi*ao**(1.-deltastr)*ani**(2.+deltastr)* &
                  (exp(gammln(NU+deltastr+2.)))/gammnu
                  deltastr=0.55
                  ani=((3.*voltmp*gammnu)/ &
                  (4.*pi*ao**(1.-deltastr)*(exp(gammln(NU+deltastr+2.)))))** &
                  (1./(2.+deltastr))
                  ai(i,k)=max((NU*ni(i,k)*ani),1.e-20)
                  ani=ai(i,k)/(NU*ni(i,k))
               else if (deltastr.gt.1.5) then
                  voltmp=(4./3.)*pi*ao**(1.-deltastr)*ani**(2.+deltastr)* &
                  (exp(gammln(NU+deltastr+2.)))/gammnu
                  deltastr=1.5
                  ani=((3.*voltmp*gammnu)/ &
                  (4.*pi*ao**(1.-deltastr)*(exp(gammln(NU+deltastr+2.)))))** &
                  (1./(2.+deltastr))
                  ai(i,k)=max((NU*ni(i,k)*ani),1.e-20)
                  ani=ai(i,k)/(NU*ni(i,k))
                  cni=ao**(1.-deltastr)*ani**deltastr
                  ci(i,k)=max((NU*ni(i,k)*cni),1.e-20)
                  cni=ci(i,k)/(nu*ni(i,k))
               endif
               
               betam = 2.+deltastr
               alphstr = ao**(1.-deltastr)
               alphv = (4./3.)*pi*alphstr
               
!     get avg ice density 
               rhobar = qi(i,k)*gammnu/(ni(i,k)*alphv* &
                        ani**betam*exp(gammln(nu+betam)))
               IF(ISDACflag .eq. 1) rhobar = 88.4
               IF(masssizeflag .eq. 1.or.sphrflag.eq.1) rhobar = 920.
               If(redden .eq. 1) rhoi = 500.

!     rhobar=500.
!     check bounds for ice density
               IF(rhobar.gt.920.)THEN 
                  
                  rhobar=920.
                  ani=((qi(i,k)*gammnu)/(rhobar*ni(i,k)*alphv*&
                  exp(gammln(nu+betam))))**(1./betam)
                  cni=ao**(1.-deltastr)*ani**deltastr
                  ci(i,k)=nu*ni(i,k)*cni
                  ai(i,k)=nu*ni(i,k)*ani
                  
               ELSE IF(rhobar.lt.50.)THEN

                  rhobar=50.
                  ani=((qi(i,k)*gammnu)/(rhobar*ni(i,k)*alphv*&
                  exp(gammln(nu+betam))))**(1./betam)
                  cni=ao**(1.-deltastr)*ani**deltastr
                  ci(i,k)=nu*ni(i,k)*cni
                  ai(i,k)=nu*ni(i,k)*ani
                  
               END IF
!     get rni (characteristic equivalent volume ice radius)


               rni = (qi(i,k)*3./(ni(i,k)*rhobar*4.*pi*&
               (exp(gammln(nu+deltastr+2.))/gammnu)))**(1./3.)

               IF(masssizeflag .eq. 1)THEN
                  rni=ani
                  IF(cni.gt.ani)rni=cni
               ENDIF
               
!     make sure rni is within reasonable bounds,
               
               IF(rni.lt.2.e-6)THEN
                  
                  rni=2.e-6
                  ni(i,k)=3.*qi(i,k)*gammnu/(4.*pi*rhobar*rni**3.*&
                  (exp(gammln(nu+deltastr+2.))))
                  ani=((qi(i,k)*gammnu)/(rhobar*ni(i,k)*alphv* &
                  exp(gammln(nu+betam))))**(1./betam) 
!                  IF(masssizeflag .eq. 1)THEN
!                     IF(ani.ge.cni)ani=rni
!                     IF(cni.gt.ani)cni=rni
!                     betam=betamp
!                     ni(i,k)=(qi(i,k)*exp(gammln(nu+betam)))/&
!                     (rni**betam*alphamsp*gammnu)
!                     IF(ani.ge.cni)cni=ao**(1.-deltastr)*ani**deltastr
!                     IF(cni.gt.ani)ani=ao**(deltastr-1.)*cni**&
!                     (1./deltastr)
!                  ELSE
                     cni=ao**(1.-deltastr)*ani**deltastr
!                  ENDIF
                  ci(i,k)=nu*ni(i,k)*cni
                  ai(i,k)=nu*ni(i,k)*ani
                  
               ELSE IF(rni.gt.2.e-3)THEN
                  
                  rni=2.e-3
                  ni(i,k)=3.*qi(i,k)*gammnu/(4.*pi*rhobar*rni**3.* &
                  (exp(gammln(nu+deltastr+2.))))
                  ani=((qi(i,k)*gammnu)/(rhobar*ni(i,k)*alphv* &
                  exp(gammln(nu+betam))))**(1./betam)
!                  IF(masssizeflag .eq. 1)THEN
!                     IF(ani.ge.cni)ani=rni
!                     IF(cni.gt.ani)cni=rni
!                     betam=betamp
!                     ni(i,k)=(qi(i,k)*exp(gammln(nu+betam)))/&
!                     (rni**betam*alphamsp*gammnu)
!                     IF(ani.ge.cni)cni=ao**(1.-deltastr)*ani**deltastr
!                     IF(cni.gt.ani)ani=ao**(deltastr-1.)*ani**&
!                     (1./deltastr)
!                  ELSE
                     cni=ao**(1.-deltastr)*ani**deltastr
!                  END IF
                  ci(i,k)=nu*ni(i,k)*cni
                  ai(i,k)=nu*ni(i,k)*ani
               END IF  
               
!     get iwc to calculate iwc tendency
!     by substracting final and initial values
               vi = 4./3.*pi*rni**3.*exp(gammln(nu+deltastr+2.))/gammnu 
               iwci = ni(i,k)*rhobar*vi*rho(i,k)
!     IF(masssizeflag.eq.1)&
!     iwci=ni(i,k)*rho(i,k)*alpham*rni**betam*&
!     (gammnu/exp(gammln(nu+deltastr+2.)))**(3./betam) 
               
!     calculate the particle inherent growth ratio (IGR) based on temp
!     igr=1 for ice particles pre-diagnosed as spheres
               IF(celsius .le. -1. .and. celsius .ge. -59.)THEN 
                  weight = (ABS(real(int(celsius))) + 1.0) - ABS(celsius)
                  igr1 = gammaFindSumTrip(int(celsius)*(-1))
                  igr2 = gammaFindSumTrip((int(celsius)*(-1))+1)
                  igr = weight*igr1 + (1.-weight)*igr2    
               ELSE IF(celsius .lt. -60.)THEN
                  igr = 1.5
               ELSE
                  igr = 1.0
               END IF

               IF(iaspect .eq. 1) igr = .27
               IF(ISDACflag .eq. 1) igr = 1.0
               IF(sphrflag .eq. 1) igr=1.0
               If(redden .eq. 1) rhoi = 500.
!     calculate number concentration from number mixing ratio

               IF(EVOLVE_ON .eq. 1) THEN
                  
                  nidum = ni(i,k)*rho(i,k)
                  CALL EVOLVE(ani,nidum,sui,sup,qvv,temp,press,igr,dt,iwcf,&
                  cnf,iwci,phii,phif,cni,rni,rnf,anf,deltastr,mu,&
                  rhobar,vtbarb,vtbarbm,alphstr,vtbarblen,rhoa,i,k,iaspect&
                  ,ISDACflag,masssizeflag,ipart,sphrflag,redden,itimestep) 
                  
                  betam = deltastr + 2.0
                  alphstr = ao**(1.-deltastr)
                  alphv = 4./3.*pi*alphstr
                  
                  phi(i,k)=phif
                  
!     get deposition/sublimation process rates for qi, ai, and ci
!     deposition rate for ice [=] kg/kg/s (mixing ratio rate)
                  
                  prd=(iwcf-iwci)/rho(i,k)/dt              
                  ard=(anf-ani)*nu*ni(i,k)/dt
                  crd=(cnf-cni)*nu*ni(i,k)/dt 
               END IF           !EVOLVE_ON
            END IF              !ice is present
!     END IF                 !ICE_CALCS
!     get ice nucleation
            IF(REAL(itimestep)*dt.gt.ice_start_time)THEN !.and.& 
               IF(PIRE_CHEM)THEN

                  nuc1 = iin_sumj(i,k,ihr)
                  !nuc1 = sum(iin_sum(i,:,k,ihr))/real(jde) !#/kg !changed iin to iin_sum for sum of bins from 4-15
                  !iin_sumj = sum(iin_sum(i,:,k,ihr),dim=2) !summing iin_sum through j dimension
                  !nuc_in = iin_sumj(i,k,ihr)
               ELSE

                  nuc1 = nuc !#/L
                  !nuc_in = nuc
               END IF

                IF(nucleation .eq. 0) then
                   CALL SIMPLENUC(temp,sui,ISDACflag,rho(i,k),ni(i,k),&
                     dt,mnuccd,nnuccd,nuc1,rdry,nin)
                ELSE IF (nucleation .eq. 1) then
                   CALL MEYERS(temp,sui,rho(i,k),ni(i,k),dt,mnuccd,pgam,&
                     lamc,CDIST1,mu,nnuccc,mnuccc,nnuccd,press,icontactflag,nuc1,&
                     rdry,nin)
                ELSE IF(nucleation .eq. 2) then
                   CALL DEMOTT(temp,sui,nnuccd,mnuccd,rho(i,k),nuc1,dt,&
                        rdry,nin,demottflag,nnucci,mnucci)                
                END IF
             END IF
 
!     make sure doesn't push into subsat or supersat
            
            iflag=0
            prd1=prd
            IF(prd.lt.0.0.and.sui.ge.0.0)THEN
               IF(prd.lt.-FUDGE*sui*qvi/abi/dt)THEN
                   prd = -FUDGE*sui*qvi/abi/dt
                   iflag = 1
               END IF
            END IF 
            IF(prd.lt.0.0.and.sui.lt.0.0)THEN
               IF(prd.lt.FUDGE*sui*qvi/abi/dt)THEN
                 prd = FUDGE*sui*qvi/abi/dt
                 iflag = 1
               END IF
            END IF
            IF((prd.gt.0..and.prd.gt.FUDGE*sui*qvi/abi/dt))THEN
               prd=FUDGE*sui*qvi/abi/dt
               iflag=1
            END IF

!     conservation checks
            prd=max(prd,-qi(i,k)/dt)
!     add growth tendencies to get updated variables
            qi(i,k)=qi(i,k)+(prd)*dt
            ai(i,k)=ai(i,k)+(ard)*dt
            ci(i,k)=ci(i,k)+(crd)*dt

            IF(prd .gt. 0.0)THEN
               icedep(i,k) = prd
            ELSE
               icesub(i,k) = prd
            ENDIF

            IF(qi(i,k).ge.qsmall.and.ni(i,k).gt.qsmall)THEN

!     if sublimation, reduce crystal number,
!     NOTE: this should not impact rhobar, since
!     rhobar contains terms with qi/ni and ratio
!     of qi/ni is assumed constant during loss of ni
               
               IF(prd.lt.0.)THEN
                  ni(i,k)=ni(i,k)+prd*ni(i,k)/qi(i,k)*dt
               END IF

!     set minimum ni to avoid taking root of a negative number
               ni(i,k)=max(ni(i,k),qsmall)

!     if iflag=1, then recalculate ai and ci assuming same deltastr
!     and same rhobar
!     this is needed for consistency between qi, ai, ci, etc.
!     since growth rate prd must be scaled back
               
               IF(iflag.eq.1)THEN
                  alphstr=ao**(1.-deltastr)
                  alphv=4./3.*pi*alphstr
                  betam=2.+deltastr   
                  
                  ani=((qi(i,k)*gammnu)/(rhobar*ni(i,k)*alphv*&
                  exp(gammln(nu+betam))))**(1./betam)
!                  IF(masssizeflag .eq. 1)THEN
!                     betam=betamp
!                     ani=((qi(i,k)*gammnu)/(ni(i,k)*alphamsp*&
!                     exp(gammln(nu+betam))))**(1./betam)
!                  END IF
                  cni=ao**(1.-deltastr)*ani**deltastr
                  ci(i,k)=nu*ni(i,k)*cni
                  ai(i,k)=nu*ni(i,k)*ani
               END IF           ! iflag = 1

               rni = (qi(i,k)*3./(ni(i,k)*rhobar*4.*pi*&
               (exp(gammln(nu+deltastr+2.))/gammnu)))**(1./3.)
!     make sure rni is within reasonable bounds,
               
               IF(rni.lt.2.e-6)THEN
                  
                  rni=2.e-6
                  ni(i,k)=3.*qi(i,k)*gammnu/(4.*pi*rhobar*rni**3*&
                  (exp(gammln(nu+deltastr+2.))))
                  ani=((qi(i,k)*gammnu)/(rhobar*ni(i,k)*alphv* &
                  exp(gammln(nu+betam))))**(1./betam) 
                  ci(i,k)=nu*ni(i,k)*cni
                  ai(i,k)=nu*ni(i,k)*ani
                  
               ELSE IF(rni.gt.2.e-3)THEN
                  
                  rni=2.e-3
                  ni(i,k)=3.*qi(i,k)*gammnu/(4.*pi*rhobar*rni**3* &
                  (exp(gammln(nu+deltastr+2.))))
                  ani=((qi(i,k)*gammnu)/(rhobar*ni(i,k)*alphv* &
                  exp(gammln(nu+betam))))**(1./betam)
                  ci(i,k)=nu*ni(i,k)*cni
                  ai(i,k)=nu*ni(i,k)*ani
               END IF  
            END IF              ! q > qsmall
            
            betam=2.+deltastr
            alphstr=ao**(1.-deltastr)
            alphv=4./3.*pi*alphstr     

!     calculate simplified aggregation for snow category 
!     kjs 02/2015
            
!     First, calulcate the aggregate-available ice,
!     which is the amount of ice that would autoconvert to snow
!     in traditional schemes. This is for r_i > 125 um.

            IF(qi(i,k).ge.1.e-8.and.qv(i,k)/qvi.gt.1..and.rni.gt.0.)THEN
               nprci = 4./DCS/rhobar*(qv(i,k)-qvi)*rho(i,k)*ni(i,k)&
               /rni*exp(-DCS/rni)*dv/abi ! #/kg/s
               prci = 4./3.*pi*rhobar*DCS**3*nprci !kg/kg/s
!     Second, if the aggregate-available ice mass mixing ratio is
!     > 1.e-8 kg/kg, then calculate total aggregate number and mass
!     and add it to the snow category.  The aggeagates are assumed 
!     to be spheres and have a density of rhosn=100 kg/m3.
               agg = 0.
               nagg = 0.
               IF(prci*dt.gt.1.e-8)THEN
                  dum2 = cons15*asn(k)
                  nagg = dum2*(rho(i,k)*prci*dt)**((2.+BS)/3.)*&
                  (nprci*dt*rho(i,k))**((4.-BS)/3.)/rho(i,k) !#/kg/s
                  agg = 4./3.*pi*rhosn*DCS**3*nagg !kg/kg/s
               END IF
!     Third, update qs, ns, qi, ni, ai, and ci due to aggregation. 
!     Assume rhobar remains the same.
!     Note: agg & nagg are negative numbers
               
               IF(nagg .ne. 0.0)THEN 
                  qs(i,k) = qs(i,k) - agg*dt
                  ns(i,k) = ns(i,k) - nagg*dt
                  
                  qi(i,k) = qi(i,k) + agg*dt
                  ni(i,k) = ni(i,k) + nagg*dt
               END IF
              
               !adding collection rates for output
               naggout(i,k) = nagg !#kg^-1 s^-1
               nraggout(i,k) = nragg
               nsaggout(i,k) = nsagg
               prciout(i,k) = prci
            END IF
!.......................................................................
! ACCRETION OF CLOUD ICE BY SNOW
! FOR THIS CALCULATION, IT IS ASSUMED THAT THE VS >> VI
! AND DS >> DI FOR CONTINUOUS COLLECTION

            IF(qs(i,k).ge.qsmall.and.ni(i,k).ge.1.e-8.and.qi(i,k).ge.qsmall)THEN

               lams = (cons1*ns(i,k)/qs(i,k))**(1./DS)
               n0s = ns(i,k)*lams

               IF(lams.lt.lammins)THEN
                  lams = lammins
                  n0s = lams**4*qs(i,k)/cons1
                  ns(i,k) = n0s/lams
               ELSE IF(lams.gt.lammaxs)THEN
                  lams = lammaxs
                  n0s = lams**4*qs(i,k)/cons1
                  ns(i,k) = n0s/lams
               END IF

               prai = cons23*asn(k)*qi(i,k)*rho(i,k)*n0s/&
                      lams**(bs+3.)
               nprai = cons23*asn(k)*ni(i,k)*rho(i,k)*n0s/lams**(bs+3.)
               nprai = min(nprai,ni(i,k)/dt)

               qs(i,k) = qs(i,k) + prai*dt
               ns(i,k) = ns(i,k) + nprai*dt

               qi(i,k) = qi(i,k) - prai*dt
               ni(i,k) = ni(i,k) - nprai*dt

               praiout(i,k) = prai
            END IF
 
            IF(qi(i,k).gt.qsmall.and.ni(i,k).gt.qsmall)THEN
               ani = ((qi(i,k)*gammnu)/(rhobar*ni(i,k)*alphv*&
               exp(gammln(nu+betam))))**(1./betam)
               cni=ao**(1.-deltastr)*ani**deltastr
               ci(i,k)=nu*ni(i,k)*cni
               ai(i,k)=nu*ni(i,k)*ani  
            END IF
            ns(i,k)=max(ns(i,k),qsmall)
            IF(qs(i,k).lt.qsmall.or.ns(i,k).lt.qsmall)THEN
               qs(i,k) = 0.0
               ns(i,k) = 0.0
            END IF
!..................................................................
!     now update qi, ni, ai, and ci due to ice nucleation
!     for simplicity, assume that rhobar and deltastr are 
!     constant during nucleation (probably not the best assumption,
!     should look into future modifications....)

!     add tendencies to temp, water vapor
!     if qi < qsmall, then zero out ni, ai, ci


!     make sure not to evaporate more liquid than available
            IF(qc(i,k).ge.qsmall)THEN
               IF((mnuccc+mnucci)*dt .lt. qc(i,k))then
                  qi(i,k)=qi(i,k)+(mnuccc+mnucci)*dt
                  ni(i,k)=ni(i,k)+(nnuccc+nnucci)*dt
               ELSE
                  qi(i,k)=qi(i,k)+qc(i,k)
                  ni(i,k)=ni(i,k)+nc(i,k)
               END IF
            END IF

            qi(i,k)=qi(i,k)+(mnuccd)*dt
            ni(i,k)=ni(i,k)+(nnuccd)*dt

!     set minimum ni to avoid division by zero
            ni(i,k)=max(ni(i,k),qsmall)

            IF(qi(i,k).ge.qsmall.and.ni(i,k).ge.qsmall)THEN

               qi_frac = min(max(qiold/qi(i,k),0.0),1.0)
               rhobar = rhoi*(1.-qi_frac)+rhobar*(qi_frac)

               ani=((qi(i,k)*gammnu)/(rhobar*ni(i,k)*alphv*&
                   exp(gammln(nu+betam))))**(1./betam)
               cni=ao**(1.-deltastr)*ani**deltastr
               ci(i,k)=nu*ni(i,k)*cni
               ai(i,k)=nu*ni(i,k)*ani

               ai(i,k)=max(ai(i,k),qsmall)
               ci(i,k)=max(ci(i,k),qsmall)

!     final check on limits for deltastr, rhobar, and rn
               
!     get cni and ani from ci and ai

               ani=ai(i,k)/(nu*ni(i,k))
               cni=ci(i,k)/(nu*ni(i,k))

!     get deltastr from ci and ai

               IF((log(ani)-log(ao)).gt. 0.01 .and.&
                  (log(cni)-log(ao)).gt.0.001)THEN
                  deltastr = (log(cni)-log(ao))/(log(ani)-log(ao))

               ELSE
                  deltastr = 1.
               ENDIF

               IF(ISDACflag .eq. 1 .or.masssizeflag .eq. 1 .or. &
               sphrflag .eq. 1) deltastr = 1.0
               IF(iaspect .eq. 1) deltastr = 0.8

!     make sure deltastr is in reasonable limits
!     if adjustment is needed, keep ai the same, adjust ci
               
!              IF(deltastr.lt.0.7)THEN
!                  deltastr=0.7           
!                  cni=ao**(1.-deltastr)*ani**deltastr
!                  ci(i,k)=nu*ni(i,k)*cni

!               ELSE IF(deltastr.gt.1.3)THEN
                  
!                  deltastr=1.3
!                  cni=ao**(1.-deltastr)*ani**deltastr
!                  ci(i,k)=nu*ni(i,k)*cni
                  
!               END IF

              if(deltastr.lt.0.55) then
                 voltmp=(4./3.)*pi*ao**(1.-deltastr)*ani**(2.+deltastr)* &
                 (exp(gammln(NU+deltastr+2.)))/gammnu
                 deltastr=0.55
                 ani=((3.*voltmp*gammnu)/ &
                      (4.*pi*ao**(1.-deltastr)*(exp(gammln(NU+deltastr+2.)))))** &
                      (1./(2.+deltastr))
                 ai(i,k)=max((NU*ni(i,k)*ani),1.e-20)
                 ani=ai(i,k)/(NU*ni(i,k))
              else if (deltastr.gt.1.5) then
                 voltmp=(4./3.)*pi*ao**(1.-deltastr)*ani**(2.+deltastr)* &
                       (exp(gammln(NU+deltastr+2.)))/gammnu
                 deltastr=1.5
                 ani=((3.*voltmp*gammnu)/ &
                     (4.*pi*ao**(1.-deltastr)*(exp(gammln(NU+deltastr+2.)))))** &
                     (1./(2.+deltastr))
                 ai(i,k)=max((NU*ni(i,k)*ani),1.e-20)
                 ani=ai(i,k)/(NU*ni(i,k))
                 cni=ao**(1.-deltastr)*ani**deltastr
                 ci(i,k)=max((NU*ni(i,k)*cni),1.e-20)
                 cni=ci(i,k)/(nu*ni(i,k))
              endif
               betam=2.+deltastr
               alphstr=ao**(1.-deltastr)
               alphv=4./3.*pi*alphstr              
               
!     get avg ice density
               rhobar = qi(i,k)*gammnu/(ni(i,k)*alphv* &
               ani**betam*exp(gammln(nu+betam)))

               IF(ISDACflag .eq. 1) rhobar = 88.4
               IF(masssizeflag .eq. 1) rhobar = 920.
               IF(sphrflag .eq. 1) rhobar = 920.
               If(redden .eq. 1) rhoi = 500.
!     check to make sure ice density in bounds
!     if necessary adjust ani, then we also need to recalculate cni,
!     assuming the same deltastr

               IF(rhobar.gt.920.)THEN

                  rhobar=920.
                  ani=((qi(i,k)*gammnu)/(rhobar*ni(i,k)*alphv*&
                  exp(gammln(nu+betam))))**(1./betam)
                  cni=ao**(1.-deltastr)*ani**deltastr
                  ci(i,k)=nu*ni(i,k)*cni
                  ai(i,k)=nu*ni(i,k)*ani

               ELSE IF(rhobar.lt.50.)THEN

                  rhobar=50.
                  ani=((qi(i,k)*gammnu)/(rhobar*ni(i,k)*alphv*&
                  exp(gammln(nu+betam))))**(1./betam)
                  cni=ao**(1.-deltastr)*ani**deltastr
                  ci(i,k)=nu*ni(i,k)*cni
                  ai(i,k)=nu*ni(i,k)*ani
                  
               END IF
!     get rni

            END IF ! qi > qsmall

            IF(qi(i,k).lt.qsmall)THEN
               qi(i,k)=0.
               ai(i,k)=0.
               ni(i,k)=0.
               ci(i,k)=0.
               rhobar=920.
            ELSE
               qv(i,k)=qv(i,k)-(prd+mnuccd)*dt
               temp=temp+((prd+mnuccd)*xxls/cpm*dt)
               icenuc(i,k) = icenuc(i,k) + mnuccd
               dep(i,k) = mnuccd
               cfrzr(i,k) = mnuccr
               qiold = qi(i,k)
               IF(qc(i,k).ge.qsmall)THEN
                  IF((mnuccc+mnucci)*dt .lt. qc(i,k))then
                     qc(i,k)=qc(i,k)-(mnuccc+mnucci)*dt
                     nc(i,k)=nc(i,k)-(nnuccc+nnucci)*dt
                     temp=temp+((mnuccc+mnucci)*xxlf/cpm*dt)
                     icenuc(i,k) = mnuccc+mnucci 
                     contact(i,k) = mnuccc
                     immersion(i,k) = mnucci
                  ELSE
                     temp=temp+qc(i,k)*xxlf/cpm
                     icenuc(i,k)=qc(i,k)/dt
                     contact(i,k)=qc(i,k)/dt*mnuccc/(mnuccc+mnucci)
                     immersion(i,k)=qc(i,k)/dt*mnucci/(mnuccc+mnucci)
                     qc(i,k)=0.0
                     nc(i,k)=0.0
                  END IF
               END IF
            END IF

            IF(qs(i,k).lt.qsmall)THEN
               qs(i,k)=0.
               ns(i,k)=0.
            END IF

!     get fallspeed parameters
!     for now, use formulation from LH74, side planes
!     express v-D relationship in terms of 'a' length for simplicity
            
!     note: no air density correction factor applied yet!
            
            IF(qi(i,k).gt.qsmall.and.ni(i,k).gt.qsmall) THEN
               
!     updated ani, cni are calculated right before final checks above
               ni(i,k) = max(ni(i,k),qsmall)
               ani=ai(i,k)/(nu*ni(i,k))
               cni=ci(i,k)/(nu*ni(i,k))
               dum=exp(gammln(nu+deltastr+2.+bif))/&
               exp(gammln(nu+deltastr+2.))
               vtdumm=aif*dum*ani**bif
               dum=exp(gammln(nu+bif))/gammnu
               vtdumn=aif*dum*ani**bif
                   
!     limit fallspeed to 5 m/s
               
               vtrmi1(i,k)=max(min(vtbarbm,5.),0.)
               vtrni1(i,k)=max(min(vtbarb,5.),0.)
               vtrli1(i,k)=max(min(vtbarblen,5.),0.)
            ELSE
               
               vtrli1(i,k)=0.
               vtrmi1(i,k)=0.
               vtrni1(i,k)=0.
            END IF
!     get ice effective radius
!     for now just set to 10 micron

            effi1(i,k)=10.
            
!     get phi
            
            IF(qi(i,k).ge.qsmall.and.ani.gt.qsmall)THEN
               phii = cni/ani*exp(gammln(nu+deltastr-1.))/gammnu
            ELSE
               phii=1.
               ani=0.
               cni=0.
               rhobar=920.
            END IF
            IF(iaspect .eq. 1) phii = 0.27
            IF(ISDACflag .eq. 1) rhobar=88.4 
            IF(masssizeflag .eq. 1) rhobar=920.
            IF(sphrflag .eq. 1) rhobar =920.
            If(redden .eq. 1) rhoi = 500.
            !rhobar=500.

!defining output variables as 2D
            prcout(i,k) = prc
            psacwsout(i,k) = psacws
            qmultsout(i,k) = qmults
            qmultrout(i,k) = qmultr
            piacrout(i,k) = piacr
            praciout(i,k) = praci
            piacrsout(i,k) = piacrs
            pracisout(i,k) = pracis
            pracgout(i,k) = pracg
            psacwgout(i,k) = psacwg
            pgsacwout(i,k) = pgsacw
            pgracsout(i,k) = pgracs
            prdgout(i,k) = prdg
            eprdgout(i,k) = eprdg
            evpmgout(i,k) = evpmg
            pgmltout(i,k) = pgmlt
            qmultgout(i,k) = qmultg
            qmultrgout(i,k) = qmultrg

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     add tendencies to cloud, rain, and snow

            qv(i,k) = qv(i,k) + &
                      (-pre-evpms-evpmg)*dt + &
                      (-pcc-prds-eprds-prdg-eprdg)*dt
            temp = temp + &
                   (pre*xxlv + (evpms+evpmg)*xxls + &
                   (psmlt+pgmlt-pracs-pracg)*xxlf)/cpm*dt +&
                   ((prds+eprds+prdg+eprdg)*xxls + pcc*xxlv + &
                   (psacws+qmults+qmultg+qmultr+qmultrg+pracs1+&
                   +psacwg+pracg1+pgsacw+pgracs+piacr+piacrs)*&
                   xxlf)/cpm*dt

            qc(i,k) = qc(i,k) + &
                      (-pra-prc)*dt +&
                      (pcc-psacws-qmults-qmultg-psacwg-pgsacw)*dt
            
            nc(i,k) = nc(i,k) + &
                      (-npra-nprc)*dt + &
                      (-npsacws-npsacwg)*dt

            qr(i,k) = qr(i,k) + &
                      (pre+prc+pra-psmlt-pgmlt+pracs+pracg)*dt + &
                      (-pracs1-qmultr-qmultrg-piacr-piacrs-&
                      pracg1-pgracs)*dt
            nr(i,k) = nr(i,k) + &
                      (nprc1+nragg-npracg-nsmltr-ngmltr)*dt + &
                      (-npracs1-niacr-niacrs-npracg1-ngracs+nsubr)*dt
            qs(i,k) = qs(i,k) + &
                      (psmlt+evpms-pracs)*dt +&
                      (psacws+prds+pracs1+eprds-psacr+piacrs+&
                      pracis)*dt
            ns(i,k) = ns(i,k) + &
                      nsmlts*dt + &
                      (nsagg-nscng-ngracs+niacrs+nsubs)*dt

            qg(i,k) = qg(i,k) + &
                      (pgmlt+evpmg-pracg)*dt + &
                      (pracg1+psacwg+pgsacw+pgracs+prdg+&
                      eprdg+piacr+psacr)*dt
            ng(i,k) = ng(i,k) + &
                      (ngmltg)*dt + &
                      (nscng+ngracs+niacr+nsubg)*dt

            qitend = (qmults+qmultg+qmultr+qmultrg-praci-pracis)
            nitend = (nmults+nmultg+nmultr+nmultrg-niacr-niacrs)
            qi(i,k) = qi(i,k) + qitend*dt
            ni(i,k) = ni(i,k) + nitend*dt


            rainevap(i,k) = -1.*pre
            snowevap(i,k) = -1.*evpms
            snowmelt(i,k) = -1.*psmlt
            snowdep(i,k) = prds
            snowsub(i,k) = eprds
            snowaccr(i,k) = pracs
!     add liquid, do saturation adjustment to get updated qc
!     get updated thermodynamics (updated after ice microphysics)
            
!     get temp from 0theta (in K)
           
            dum=min(0.99*press,polysvp1(temp,0))
            qvs=0.622*dum/(press-dum) 
            dum=qv(i,k)-qvs
            pcc=dum/(1.+xxlv**2*qvs/(cpm*rv*temp**2))/dt
            IF(pcc*dt+qc(i,k).lt.0.)THEN
               pcc=-qc(i,k)/dt
            END IF
            !condrate(i,k)=pcc
!     update variables due to condensation/evaporation

            temp=temp+pcc*xxlv/cpm*dt
            qv(i,k)=qv(i,k)-pcc*dt
            qc(i,k)=qc(i,k)+pcc*dt
            IF(pcc .gt. 0.0)THEN
               cloudcond(i,k) = pcc
            ELSE
               cloudevap(i,k) = pcc
            END IF
!     for now, assume droplet fallspeed is zero           
            vtrmc(i,k)=0.

            !temp3(i,k)=temp
            t(i,k) = temp

!     add melting, melt all ice within one time step above freezing
            IF(t(i,k).ge.273.15)THEN
               IF(qi(i,k).ge.qsmall)THEN
                  qr(i,k)=qr(i,k)+qi(i,k)
                  nr(i,k)=nr(i,k)+ni(i,k)
                  t(i,k)=t(i,k)-qi(i,k)*xxlf/cpm
                  icemelt(i,k) = qi(i,k)/dt
                  qi(i,k)=0.
                  ni(i,k)=0.
                  ai(i,k)=0.
                  ci(i,k)=0.
               END IF
!               IF(qs(i,k).ge.qsmall)THEN
!                  qr(i,k)=qr(i,k)+qs(i,k)
!                  nr(i,k)=nr(i,k)+ns(i,k)
!                  t(i,k)=t(i,k)-qs(i,k)*xxlf/cpm
!                  qs(i,k)=0.
!                  ns(i,k)=0.
!               END IF
            END IF
!     homogeneous freezing, freeze andd cloud and rain water within one time-step below -40
            IF(homofreeze .eq. 1.and.t(i,k).le.233.15.and.&
               real(itimestep)*dt.gt.ice_start_time)THEN
               IF(nr(i,k).ge.qsmall.and.qr(i,k).ge.qsmall)THEN

                  qi(i,k)=qi(i,k)+qr(i,k)
                  ni(i,k)=ni(i,k)+nr(i,k)
                  t(i,k)=t(i,k)+qr(i,k)*(xxls-xxlv)/cpm
!                  nnew = qr(i,k)/((4./3.)*pi*rhoi*(200.e-6)**3)
                  ai(i,k) = ai(i,k) + nu*nr(i,k)*(qr(i,k)*3./(4.*pi*&
                  rhoi*nr(i,k))*gammnu/exp(gammln(nu+3)))**(1./3.)
                  ci(i,k) = ci(i,k) + nu*nr(i,k)*(qr(i,k)*3./(4.*pi*&
                  rhoi*nr(i,k))*gammnu/exp(gammln(nu+3)))**(1./3.)
                  rainfrz(i,k) = qr(i,k)/dt 
                  qr(i,k)=0.
                  nr(i,k)=0.

               END IF
               IF(nc(i,k).ge.qsmall.and.qc(i,k).ge.qsmall)THEN

                  qi(i,k)=qi(i,k)+qc(i,k)
                  ni(i,k)=ni(i,k)+nc(i,k)
                  t(i,k)=t(i,k)+qc(i,k)*(xxls-xxlv)/cpm
                  !nnew = qc(i,k)/((4./3.)*pi*rhoi*(20.e-6)**3)
                  ai(i,k) = ai(i,k) + nu*nc(i,k)*(qc(i,k)*3./(4.*pi*&
                  rhoi*nc(i,k))*gammnu/exp(gammln(nu+3)))**(1./3.)
                  ci(i,k) = ci(i,k) + nu*nc(i,k)*(qc(i,k)*3./(4.*pi*&
                  rhoi*nc(i,k))*gammnu/exp(gammln(nu+3)))**(1./3.)
                  !ni(i,k) = ni(i,k) + nnew
                  cloudfrz(i,k) = qc(i,k)/dt
                  qc(i,k)=0.
                  nc(i,k)=0.

               END IF
            END IF
            LTRUE = 1

 200        CONTINUE

            IF(qi(i,k) .gt. qsmall)THEN
               ani = ai(i,k)/(ni(i,k)*nu)
               cni = ci(i,k)/(ni(i,k)*nu)

               ai(i,k) = ani**2*cni*nu*ni(i,k)
               ci(i,k) = cni**2*ani*nu*ni(i,k)
            END IF
            rhoice(i,k)=rhobar

         END DO !end k

         precprt(i) = 0.0
         snowrt(i) = 0.0
         snowprt(i) = 0.0
         grplprt(i) = 0.0


         IF(LTRUE.eq.0)GOTO 400
!     SEDIMENTATION --------------------------------------------------------------
        IF(SEDON .eq. 1)THEN
         nstep = 1
         
         DO k = kte,kts,-1
           
            IF(qi(i,k).ge.qsmall)THEN
               fi(k) = vtrmi1(i,k)
               fni(k) = vtrni1(i,k)
               fai(k) = vtrmi1(i,k)
               fci(k) = vtrmi1(i,k)
            ELSE
               fi(k) = 0.0
               fni(k) = 0.0
               fai(k) = 0.0
               fci(k) = 0.0
            END IF
            !sedm_i(i,k) = fi(k)

!     calculate rain fallspeed
            nr(i,k) = max(nr(i,k),0.0)
            IF(qr(i,k).ge.qsmall)THEN!.and.nr(i,k).gt.0.0)THEN
               lamr = (pi*rhow*nr(i,k)/qr(i,k))**(1./3.)
               lamr = max(lamr,lamminr)
               lamr = min(lamr,lammaxr)
               fr(k) = arn(k)*cons4/lamr**br
               fnr(k) = arn(k)*cons6/lamr**br
            ELSE
               fr(k) = 0.
               fnr(k) = 0.
            END IF
            fr(k) = min(fr(k),9.1*(rhosu/rho(i,k))**0.54)
            fnr(k) = min(fnr(k),9.1*(rhosu/rho(i,k))**0.54)

!     hm add cloud water sedimentation

            fc(k) = 0.0
            fnc(k) = 0.0
            nc(i,k) = max(nc(i,k),0.)
            IF(qc(i,k).ge.qsmall)THEN!.and.nc(i,k).gt.0.0)THEN
               DUM = p(i,k)/(287.15*t(i,k))
               PGAM = 0.0005714*(nc(i,k)/1.E6*DUM)+0.2741
               PGAM = 1./(PGAM**2)-1
               PGAM = MAX(PGAM,2.)
               PGAM = MIN(PGAM,10.)

               lamc = (CONS26*nc(i,k)*GAMMA(PGAM+4.)/&
               (qc(i,k)*GAMMA(PGAM+1.)))**(1./3.)
               LAMMIN = (PGAM+1.)/60.E-6
               LAMMAX = (PGAM+1.)/1.E-6
               lamc=MAX(lamc,LAMMIN)
               lamc=MIN(lamc,LAMMAX)

               fc(k) = acn(k)*GAMMA(4.+BC+PGAM)/(lamc**BC*GAMMA(PGAM+4.))
               fnc(k) = acn(k)*gamma(1.+BC+pgam)/(lamc**BC*gamma(pgam+1.))

             ELSE
               fc(k) = 0.
               fnc(k) = 0.
             END IF
            fc(k) = min(fc(k),9.1*(rhosu/rho(i,k))**0.54)
            fnc(k) = min(fnc(k),9.1*(rhosu/rho(i,k))**0.54)

!     calculate snow/aggregate sedimentation
            IF(qs(i,k).ge.qsmall)THEN
               dum = (cons1*ns(i,k)/qs(i,k))**(1./DS)
               fs(k) = asn(k)*cons3/dum**BS
               fns(k) = asn(k)*cons5/dum**BS
            ELSE
               fs(k) = 0.
               fns(k) = 0.
            END IF
            fs(k) = min(fs(k),1.2*(rhosu/rho(i,k))**0.54)
            fns(k) = min(fns(k),1.2*(rhosu/rho(i,k))**0.54)

            IF (qg(i,k).ge.qsmall) THEN
               dum = (cons2*ng(i,k)/qg(i,k))**(1./DG)
               dum=MAX(dum,LAMMING)
               dum=MIN(dum,LAMMAXG)
               fg(k) = agn(k)*cons7/(dum**BG)
               fng(k) = agn(k)*cons8/dum**BG
            ELSE
               fg(k) = 0.
               fng(k) = 0.
            END IF
            fg(k) = min(fg(k),20.*(rhosu/rho(i,k))**0.54)
            fng(k) = min(fng(k),20.*(rhosu/rho(i,k))**0.54)

!modify fallspeed below level of precip
            IF(k.le.kte-1)THEN
               IF(fi(k).lt.1.e-10) fi(k) = fi(k+1)
               IF(fni(k).lt.1.e-10) fni(k) = fni(k+1)
               IF(fai(k).lt.1.e-10) fai(k) = fai(k+1)
               IF(fci(k).lt.1.e-10) fci(k) = fci(k+1)
               IF(fs(k).lt.1.e-10) fs(k) = fs(k+1)
               IF(fns(k).lt.1.e-10) fns(k) = fns(k+1)
               IF(fr(k).lt.1.e-10) fr(k) = fr(k+1)
               IF(fnr(k).lt.1.e-10) fnr(k) = fnr(k+1)
               IF(fc(k).lt.1.e-10) fc(k) = fc(k+1)
               IF(fnc(k).lt.1.e-10) fnc(k) = fnc(k+1)
               IF(fg(k).lt.1.e-10) fg(k) = fg(k+1)
               IF(fng(k).lt.1.e-10) fng(k) = fng(k+1)
            END IF
!     calculate number of split time steps

            rgvm = max(fi(k),fni(k),fai(k),fci(k),fr(k),fnr(k),fc(k),&
               fnc(k),fs(k),fns(k),fg(k),fng(k))
            nstep = max(int(rgvm*dt/dzq(k)+1.),nstep)
            
!     multiply variables by rho for right units
            dumi(k) = qi(i,k)*rho(i,k) !kgm^-3
            dumni(k) = ni(i,k)*rho(i,k) !#m^-3
            dumai(k) = ai(i,k)*rho(i,k) !unitless (ai is volume now)
            dumci(k) = ci(i,k)*rho(i,k) !unitless (ci is volume now)
            dumr(k) = qr(i,k)*rho(i,k) !kgm^-3
            dumc(k) = qc(i,k)*rho(i,k) !kgm^-3
            dumnr(k) = nr(i,k)*rho(i,k) !kgm^-3
            dumnc(k) = nc(i,k)*rho(i,k) !kgm^-3
            dums(k) = qs(i,k)*rho(i,k)
            dumns(k) = ns(i,k)*rho(i,k)
            dumg(k) = qg(i,k)*rho(i,k)
            dumng(k) = ng(i,k)*rho(i,k)
         END DO

         DO n = 1, nstep
            DO k = kts,kte

               fallouti(k) = fi(k)*dumi(k) !kgm^-2 s^-1
               falloutni(k) = fni(k)*dumni(k) !#m^-2 s^-1
               falloutai(k) = fai(k)*dumai(k) !m s^-1
               falloutci(k) = fci(k)*dumci(k) !m s^-1
               falloutr(k) = fr(k)*dumr(k) !kgm^-2 s^-1
               falloutc(k) = fc(k)*dumc(k) !kgm^-2 s^-1
               falloutnr(k) = fnr(k)*dumnr(k) !kgm^-2 s^-1
               falloutnc(k) = fnc(k)*dumnc(k) !kgm^-2 s^-1
               fallouts(k) = fs(k)*dums(k)
               falloutns(k) = fns(k)*dumns(k)
               falloutg(k) = fg(k)*dumg(k)
               falloutng(k) = fng(k)*dumng(k)

            END DO
            
!     top of model
                      
            k = kte

            falltndi = fallouti(k)/dzq(k) !kgm^-3 s^-1
            falltndni = falloutni(k)/dzq(k) !#m^-3 s^-1
            falltndai = falloutai(k)/dzq(k) !s^-1
            falltndci = falloutci(k)/dzq(k) !s^-1
            falltndr = falloutr(k)/dzq(k) !kgm^-3 s^-1
            falltndc = falloutc(k)/dzq(k) !kgm^-3 s^-1
            falltndnr = falloutnr(k)/dzq(k) !kgm^-3 s^-1
            falltndnc = falloutnc(k)/dzq(k) !kgm^-3 s^-1
            falltnds = fallouts(k)/dzq(k)
            falltndns = falloutns(k)/dzq(k)
            falltndg = falloutg(k)/dzq(k)
            falltndng = falloutng(k)/dzq(k)
            
!     sedimentation tendencies should be renewed every time step
!     so be sure to initialize to zero at beginning of code (kjs)
           
            qisten(k) = qisten(k) - falltndi/nstep/rho(i,k) !s^-1
            nisten(k) = nisten(k) - falltndni/nstep/rho(i,k) !#kg^-1 s^-1
            aisten(k) = aisten(k) - falltndai/nstep/rho(i,k) !m^3 kg^-1 s^-1
            cisten(k) = cisten(k) - falltndci/nstep/rho(i,k) !m^3 kg^-1 s^-1
            qrsten(k) = qrsten(k) - falltndr/nstep/rho(i,k) !s^-1
            qcsten(k) = qcsten(k) - falltndc/nstep/rho(i,k) !s^-1
            nrsten(k) = nrsten(k) - falltndnr/nstep/rho(i,k) !s^-1
            ncsten(k) = ncsten(k) - falltndnc/nstep/rho(i,k) !s^-1
            qssten(k) = qssten(k) - falltnds/nstep/rho(i,k)
            nssten(k) = nssten(k) - falltndns/nstep/rho(i,k)
            qgsten(k) = qgsten(k) - falltndg/nstep/rho(i,k)
            ngsten(k) = ngsten(k) - falltndng/nstep/rho(i,k)
            
            dumi(k) = dumi(k) - falltndi*dt/nstep             
            dumni(k) = dumni(k) - falltndni*dt/nstep
            dumai(k) = dumai(k) - falltndai*dt/nstep
            dumci(k) = dumci(k) - falltndci*dt/nstep
            dumr(k) = dumr(k) - falltndr*dt/nstep
            dumc(k) = dumc(k) - falltndc*dt/nstep
            dumnr(k) = dumnr(k) - falltndnr*dt/nstep
            dumnc(k) = dumnc(k) - falltndnc*dt/nstep
            dums(k) = dums(k) - falltnds*dt/nstep
            dumns(k) = dumns(k) - falltndns*dt/nstep
            dumg(k) = dumg(k) - falltndg*dt/nstep
            dumng(k) = dumng(k) - falltndng*dt/nstep            
            
            DO k = kte-1,kts,-1
               
               falltndi = (fallouti(k+1) - fallouti(k))/dzq(k)
               falltndni = (falloutni(k+1) - falloutni(k))/dzq(k)
               falltndai = (falloutai(k+1) - falloutai(k))/dzq(k)
               falltndci = (falloutci(k+1) - falloutci(k))/dzq(k)
               falltndr = (falloutr(k+1) - falloutr(k))/dzq(k)
               falltndc = (falloutc(k+1) - falloutc(k))/dzq(k)
               falltndnr = (falloutnr(k+1) - falloutnr(k))/dzq(k)
               falltndnc = (falloutnc(k+1) - falloutnc(k))/dzq(k)
               falltnds = (fallouts(k+1) - fallouts(k))/dzq(k)
               falltndns = (falloutns(k+1) - falloutns(k))/dzq(k)
               falltndg = (falloutg(k+1) - falloutg(k))/dzq(k)
               falltndng = (falloutng(k+1) - falloutng(k))/dzq(k)
               

               qisten(k) = qisten(k) + falltndi/nstep/rho(i,k)
               nisten(k) = nisten(k) + falltndni/nstep/rho(i,k)
               aisten(k) = aisten(k) + falltndai/nstep/rho(i,k)
               cisten(k) = cisten(k) + falltndci/nstep/rho(i,k)
               qrsten(k) = qrsten(k) + falltndr/nstep/rho(i,k)
               qcsten(k) = qcsten(k) + falltndc/nstep/rho(i,k)
               nrsten(k) = nrsten(k) + falltndnr/nstep/rho(i,k)
               ncsten(k) = ncsten(k) + falltndnc/nstep/rho(i,k)
               qssten(k) = qssten(k) + falltnds/nstep/rho(i,k)
               nssten(k) = nssten(k) + falltndns/nstep/rho(i,k)
               qgsten(k) = qgsten(k) + falltndg/nstep/rho(i,k)
               ngsten(k) = ngsten(k) + falltndng/nstep/rho(i,k)

               dumi(k) = dumi(k) + falltndi*dt/nstep
               dumni(k) = dumni(k) + falltndni*dt/nstep
               dumai(k) = dumai(k) + falltndai*dt/nstep
               dumci(k) = dumci(k) + falltndci*dt/nstep
               dumr(k) = dumr(k) + falltndr*dt/nstep
               dumc(k) = dumc(k) + falltndc*dt/nstep
               dumnr(k) = dumnr(k) + falltndnr*dt/nstep
               dumnc(k) = dumnc(k) + falltndnc*dt/nstep
               dums(k) = dums(k) + falltnds*dt/nstep
               dumns(k) = dumns(k) + falltndns*dt/nstep
               dumg(k) = dumg(k) + falltndg*dt/nstep
               dumng(k) = dumng(k) + falltndng*dt/nstep
  

               qiloss(k)=qiloss(k)+fallouti(k)
               qrloss(k)=qrloss(k)+falloutr(k)+falloutc(k)


               ised(i,k) = ised(i,k) + fallouti(k)/nstep
               ssed(i,k) = ssed(i,k) + fallouts(k)/nstep
               gsed(i,k) = gsed(i,k) + falloutg(k)/nstep
               rsed(i,k) = rsed(i,k) + falloutr(k)/nstep
            END DO
!     get precipitation and snowfall accumulation during time step
!     NOTE: factor of 1000 converts m to mm, but division by density cancels this
            
            precprt(i) = precprt(i) + (fallouti(kts)+falloutr(kts)+&
            falloutc(kts)+fallouts(kts)+falloutg(kts))*dt/nstep !kgm^-2

            snowrt(i) = snowrt(i) + (fallouti(kts)+fallouts(kts)+falloutg(kts))*dt/nstep
            snowprt(i) = snowprt(i) + (fallouti(kts)+fallouts(kts))*dt/nstep
            grplprt(i) = grplprt(i) + (falloutg(kts))*dt/nstep

         END DO
                           !end nstep loop
         END IF  !SEDON       
!     add on sedimenation tendencies for mixing ratio to rest of tendencies
         DO k=kts,kte
!     add new tendencies to mixing ratios


            qi(i,k) = qi(i,k) + qisten(k)*dt !kg kg^-1
            ni(i,k) = ni(i,k) + nisten(k)*dt !#kg^-1
            ai(i,k) = ai(i,k) + aisten(k)*dt !m kg^-1
            ci(i,k) = ci(i,k) + cisten(k)*dt !m kg^-1
            qr(i,k) = qr(i,k) + qrsten(k)*dt !kg kg^-1
            qc(i,k) = qc(i,k) + qcsten(k)*dt !kg kg^-1
            nr(i,k) = nr(i,k) + nrsten(k)*dt !kg kg^-1
            nc(i,k) = nc(i,k) + ncsten(k)*dt !kg kg^-1
            qs(i,k) = qs(i,k) + qssten(k)*dt
            ns(i,k) = ns(i,k) + nssten(k)*dt
            qg(i,k) = qg(i,k) + qgsten(k)*dt
            ng(i,k) = ng(i,k) + ngsten(k)*dt

            evs = min(FUDGE*p(i,k),polysvp1(t(i,k),0)) !Pa
            evi = min(FUDGE*p(i,k),polysvp1(t(i,k),1)) !pa
            qvs = 0.622*evs/(p(i,k)-evs)
            qvi = 0.622*evi/(p(i,k)-evi)
            qvqvs = qv(i,k)/qvs
            qvqvsi = qv(i,k)/qvi

            IF(qvqvs.lt.0.9)THEN
               IF(qr(i,k).lt.1.e-8)THEN
                  qv(i,k) = qv(i,k) + qr(i,k)
                  t(i,k) = t(i,k) - qr(i,k)*xxlv/cpm
                  qr(i,k) = 0.0
               END IF
               IF(qc(i,k).lt.1.e-8)THEN
                  qv(i,k) = qv(i,k) + qc(i,k)
                  t(i,k) = t(i,k) - qc(i,k)*xxlv/cpm
                  qc(i,k) = 0.0
               END IF
            END IF

            IF(qvqvsi.lt.0.9)THEN
               IF(qi(i,k).lt.1.e-8)THEN
                  qv(i,k) = qv(i,k) + qi(i,k)
                  t(i,k) = t(i,k) - qi(i,k)*xxls/cpm
                  qi(i,k) = 0.0
               END IF
               IF(qs(i,k).lt.1.e-8)THEN
                  qv(i,k) = qv(i,k) + qs(i,k)
                  t(i,k) = t(i,k) - qs(i,k)*xxls/cpm
                  qs(i,k) = 0.0
               END IF
               IF(qg(i,k).lt.1.e-8)THEN
                  qv(i,k) = qv(i,k) + qg(i,k)
                  t(i,k) = t(i,k) - qg(i,k)*xxls/cpm
                  qg(i,k) = 0.0
               END IF
            END IF
            IF(qi(i,k) .lt. qsmall .or. ni(i,k) .lt. qsmall)THEN
               qi(i,k) = 0.
               ni(i,k) = 0.
               ai(i,k) = 0.
               ci(i,k) = 0.
               qiloss(k) = 0.
            END IF
            IF(qr(i,k) .lt. qsmall .or. nr(i,k) .lt. qsmall)THEN
               qr(i,k) = 0.
               nr(i,k) = 0.
               qrloss(k) = 0.
            END IF
            IF(qc(i,k) .lt. qsmall .or. nc(i,k) .lt. qsmall)THEN
               qc(i,k) = 0.
               nc(i,k) = 0.
               qcloss(k) = 0.
            END IF
            IF(qs(i,k).lt.qsmall .or. ns(i,k) .lt. qsmall)THEN
               qs(i,k) = 0.0
               ns(i,k) = 0.0
            END IF
            IF(qg(i,k).lt.qsmall .or. ng(i,k) .lt. qsmall)THEN
               qg(i,k) = 0.0
               ng(i,k) = 0.0
            END IF

            ni(i,k) = max(0.,ni(i,k))
            ns(i,k) = max(0.,ns(i,k))
            ng(i,k) = max(0.,ng(i,k))
            nc(i,k) = max(0.,nc(i,k))
            nr(i,k) = max(0.,nr(i,k))

!     recalculate cloud, rain, and snow distributions
            IF(qr(i,k).ge.qsmall)THEN
               lamr = (pi*rhow*nr(i,k)/qr(i,k))**(1./3.)
               
               IF(lamr.lt.lamminr)THEN
                  lamr = lamminr
                  n0rr = lamr**4*qr(i,k)/(pi*rhow)
                  nr(i,k) = n0rr/lamr
               ELSE IF(lamr.gt.lammaxr)THEN
                  lamr = lammaxr
                  n0rr = lamr**4*qr(i,k)/(pi*rhow)
                  nr(i,k) = n0rr/lamr
               END IF
            END IF
            IF(qc(i,k) .ge. qsmall)THEN
               dum = p(i,k)/(287.15*t(i,k))
               pgam = 0.0005714*(nc(i,k)/1.e6*dum)+0.2714
               pgam = 1./(pgam**2)-1.
               pgam = max(pgam,2.)
               pgam = min(pgam,10.)

               lamc = (cons26*nc(i,k)*gamma(pgam+4.)/&
               (qc(i,k)*gamma(pgam+1.)))**(1./3.)

               lammin = (pgam+1.)/60.e-6
               lammax = (pgam+1.)/1.e-6

               IF(lamc.lt.lammin)THEN
                  lamc = lammin
                  nc(i,k) = exp(3.*log(lamc)+log(qc(i,k))+&
                  log(gamma(pgam+1.))-log(gamma(pgam+4.)))/cons26
               ELSE IF(lamc.gt.lammax)THEN
                  lamc = lammax
                  nc(i,k) = exp(3.*log(lamc)+log(qc(i,k))+&
                  log(gamma(pgam+1.))-log(gamma(pgam+4.)))/cons26
               END IF
            END IF

            IF(qs(i,k).ge.qsmall)THEN
               lams = (cons1*ns(i,k)/qs(i,k))**(1./DS)

               IF(lams.lt.lammins)THEN
                  lams = lammins
                  n0s = lams**4*qs(i,k)/cons1
                  ns(i,k) = n0s/lams
               ELSE IF(lams.gt.lammaxs)THEN
                  lams = lammaxs
                  n0s = lams**4*qs(i,k)/cons1
                  ns(i,k) = n0s/lams
               END IF
            END IF

            if(qg(i,k).ge.qsmall)then
                lamg = (cons2*ng(i,k)/qg(i,k))**(1./DG)
                n0g = ng(i,k)*lamg
                if(lamg.lt.lamming)then
                   lamg = lamming
                   n0g = lamg**4*qg(i,k)/cons2
                   ng(i,k) = n0g/lamg
                else if(lamg.gt.lammaxg)then
                   lamg = lammaxg
                   n0g = lamg**4*qg(i,k)/cons2
                   ng(i,k) = n0g/lamg
                end if
            end if

            ni(i,k) = max(0.,ni(i,k))
            ns(i,k) = max(0.,ns(i,k))
            ng(i,k) = max(0.,ng(i,k))
            nc(i,k) = max(0.,nc(i,k))
            nr(i,k) = max(0.,nr(i,k))
         END DO
 400     CONTINUE
      END DO     
               !end i loop
      RETURN

      END SUBROUTINE SULIAHARRINGTON_MICRO

!THIS SUBROUTINE:
!  READS IN BINARY DATA FROM THE ADVANCED PARTICLE MICROPHYSICS SCHEME (FANGQUN YU) 
!    FOR A GIVEN DOMAIN AND DAY 
!  LOADS ICE NUCLEI NUMBER FOR A GIVEN SIZE (15 SIZES: 1.500E-02, 3.182E-02, 
!    6.364E-02, 1.255E-01, 2.291E-01, 3.873E-01, 6.225E-01, 9.843E-01, 1.531E+00, 
!    2.312E+00, 3.269E+00, 5.214E+00, 9.708E+00, 1.632E+01, 2.541E+01 MICRONS)
!    N_DIST 1-15
!  LOADS CCN NUMBER FOR A GIVEN SUPERSATURATION (3 SUPERSATURATIONS: 0.8, 0.4, AND 0.2%)
!    N_DIST 16-18
      SUBROUTINE APM(itimestep,dt,id,ids,ide,jds,jde,kds,kde,nhrs,&
                     IIN,IIN_SUM,CCNOUT,start_day,size_opt,IIN_SUMJ)
      IMPLICIT NONE

      INTEGER, PARAMETER :: nccn = 3
      INTEGER, INTENT(IN) :: itimestep, id, start_day,nhrs,&
                             ids, ide, jds, jde, kds, kde, size_opt
      REAL, INTENT(IN) :: dt
      INTEGER i, j, k, n, h, ff1
      REAL elapsed_time, NN(ide-1,jde-1,kde-1,18,nhrs)
      REAL, INTENT(OUT) :: IIN(ide,jde,kde,nhrs), CCNOUT(ide,kde,nhrs,8),&
                           IIN_SUM(ide,jde,kde,nhrs),IIN_SUMJ(ide,kde,nhrs)

      REAL CCN(ide,jde,kde,nhrs,nccn)

      CHARACTER(len=1) :: dd,ff
      CHARACTER(len=2) :: fff
      CHARACTER(len=68) :: filename1
      CHARACTER(len=88) :: filename2
      CHARACTER(len=62) :: path

      elapsed_time = real(start_day) + real(itimestep*dt)/86400.
      ff1 = int(elapsed_time)

      path ='/network/asrc/scratch/asrcall/PIRE_APM/&
             wrfdustout-30L-3hr-NAM/' !LG 30 layer NAM data
!      path ='/network/asrc/scratch/asrcall/PIRE_APM/&
!             wrfdustout-40L-3hr-NAM/' !LG 40 layer NAM data
!      path ='/network/asrc/scratch/asrcall/PIRE_APM/&
!             wrfdustout-30L-3hr/' !LG 30 layer GFS data

      write(dd,'(i1)') id
      if(ff1.lt.10)then
         write(ff,'(i1)') ff1
         filename2 = path//'WRFOUTd0'//dd//'-2014-01-0'//ff//'.bin' !3-hourly average
         print*,'grid id: ',dd,'; day: ',ff,'; file: ',filename2
      else
         write(fff,'(i2)') ff1
         filename2 = path//'WRFOUTd0'//dd//'-2013-12-'//fff//'.bin'
         print*,'grid id: ',dd,'; day: ',fff,'; file: ',filename2
      end if

      open(100,file=filename2,access='direct',form='unformatted',recl=ide*jde*kde*18*nhrs)
      read(100,rec=1) NN

      IIN(:,:,:,:) = 0.0
      IIN_SUM(:,:,:,:) = 0.0
      IIN_SUMJ(:,:,:) = 0.0
      CCN(:,:,:,:,:) = 0.0

      DO i = ids,ide-1
         DO k = kds,kde-1
            DO h = 1, nhrs
               DO j = jds,jde-1 
               
                  IIN(i,j,k,h) = NN(i,j,k,size_opt,h)
                  IIN_SUM(i,j,k,h) = SUM(NN(i,j,k,4:15,h)) 
                  CCN(i,j,k,h,1) = NN(i,j,k,18,h)
                  CCN(i,j,k,h,2) = NN(i,j,k,17,h)
                  CCN(i,j,k,h,3) = NN(i,j,k,16,h)
               END DO!j
            

               IIN_SUMJ(i,k,h) = SUM(IIN_SUM(i,:,k,h))/real(jde)
               CCNOUT(i,k,h,1) = 0.0
               CCNOUT(i,k,h,2) = sum(CCN(i,:,k,h,1))/real(jde)
               CCNOUT(i,k,h,4) = sum(CCN(i,:,k,h,2))/real(jde)
               CCNOUT(i,k,h,8) = sum(CCN(i,:,k,h,3))/real(jde)
               CCNOUT(i,k,h,3) = (CCNOUT(i,k,h,2)+CCNOUT(i,k,h,4))/2.0
               CCNOUT(i,k,h,6) = (CCNOUT(i,k,h,4)+CCNOUT(i,k,h,8))/2.0
               CCNOUT(i,k,h,5) = (CCNOUT(i,k,h,4)+CCNOUT(i,k,h,6))/2.0
               CCNOUT(i,k,h,7) = (CCNOUT(i,k,h,6)+CCNOUT(i,k,h,8))/2.0
            END DO!h
         END DO!k
      END DO!i


      close(100)
      END SUBROUTINE APM

!---------------------------------------------------------------------
      SUBROUTINE SIMPLENUC(temp,sui,ISDACflag,rhodum,nidum,dt,mnuccd,&
        nnuccd,nuc1,rdry,nin)
!---------------------------------------------------------------------

      IMPLICIT NONE

      REAL, INTENT(IN) :: temp, dt, nuc1, rhodum, nidum, rdry(nin)
      REAL, INTENT(OUT) :: mnuccd, nnuccd
      REAL :: dum, mi, wght, rdry1(nin), denom
      INTEGER, INTENT(IN) :: ISDACflag, nin
      INTEGER :: k
      REAL*8 :: sui
            
      IF(temp.lt.268.15 .and.sui.ge.0.05)THEN
      !qc(i,k).gt.1.e-7)THEN

        IF (PIRE_CHEM) THEN
           dum = nuc1 !#/kg

           denom = 0.
           DO k = nin, 4, -1
              denom = denom + real(k)
           END DO

           mi = 0.
           DO k = 4,nin
              wght = real(k)/114.
              mi = mi + 4./3.*pi*rhoi*((rdry1(k)*1.e-6)**3)*wght
           END DO
        ELSE           
           dum = nuc1*1000./rhodum !1/L*1000/rho ~ 1000/kg
           mi = mi0
        END IF

        IF(ISDACflag .eq. 0)THEN                  
!           dum = nuc*1000./rhodum !1/L*1000/rho ~ 1000/kg
           !IF(nidum.lt.dum)THEN
           IF(nidum .lt. dum) THEN !dum used to be nuc1
              nnuccd=(dum-nidum)/dt
           END IF

         ELSE IF(ISDACflag .eq. 1) THEN
!           dum = nuc*1000./rhodum !1/L*1000/rho ~ 1000/kg
           mi = 4./3.*pi*88.4*(5.e-6)**3
           nnuccd = max(0.0,(dum-nidum)/dt)
        END IF

        mnuccd = nnuccd*mi
      END IF 

      END SUBROUTINE SIMPLENUC

!----------------------------------------------------------------------
      SUBROUTINE MEYERS(temp,sui,rhodum,nidum,dt,mnuccd,pgam,&
        lamc,cdist,mu,nnuccc,mnuccc,nnuccd,press,icontactflag,nuc1,rdry,&
        nin)

!----------------------------------------------------------------------
! NEED TO FIGURE OUT HOW TO USE NUC1 HERE LG        
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: nin
      REAL :: temp, dum, nnuccd, dt, mnuccd, nuc1, mi
      REAL*8 :: sui
      REAL :: rhodum, nidum, mnuccc, DAP, nnuccc, cdist, PGAM, lamc,&
       Nic,mu, press, nucdum, rdry(nin),denom, wght
      INTEGER :: icontactflag, k

      IF(temp.lt.268.15 .and.sui.ge.0.05 .and.temp.gt.238.15)THEN
      !qc(i,k).gt.1.e-7)THEN

!     depositional freezing                  
         dum=exp(-0.639+0.1296*100.*sui)*1000./rhodum !number mixing
!         ratio

         denom = 0.
         DO k = nin, 4, -1
            denom = denom + real(k)
         END DO

         mi = 0.
         DO k = 4,nin
            wght = real(k)/denom
            mi = mi + 4./3.*pi*rhoi*((rdry(k)*1.e-6)**3)*wght
         END DO


         IF (PIRE_CHEM) THEN
            nucdum = nuc1 !keep as #/kg
         ELSE
            nucdum = nuc1*(1000./rhodum) !converting #/L to #/kg for min comparison below
         END IF

         dum = min(dum,nucdum) !limits ice mixing ratio to aerosol amount

         nnuccd=dum/dt
     
         mnuccd = (nnuccd*mi)

!     contact freezing
        if (icontactflag .eq. 1) then
        
           dum = 7.37*temp/(288.*10.*press)/100.

           Nic = exp(-2.8 + 0.262*(273.15-temp))

           DAP = 4.*pi*1.38e-23/(6.*pi*RIN)*temp*(1.+dum/RIN)/mu
           
           mnuccc = pi*pi/3.*rhow*DAP*Nic*&
              exp(log(cdist)+log(gamma(PGAM+5.))-4.*log(lamc))

           nnuccc = 2.*pi*DAP*Nic*cdist*gamma(pgam+2.)/lamc
        
        end if
      END IF

      END SUBROUTINE MEYERS


!---------------------------------------------------------------------
      SUBROUTINE DEMOTT(temp,sui,nnuccd,mnuccd,rhodum,nuc1,dt,&
        rdry,nin, demottflag,nnucci,mnucci)
!DeMott et al. (2010)
!---------------------------------------------------------------------
      IMPLICIT NONE

      INTEGER :: k
      INTEGER, INTENT(IN) :: demottflag, nin
      REAL :: denom, dum, wght, mi, cf, alpha, beta, lamda, delta, rdry1(nin)
      REAL, INTENT(IN) :: temp, rhodum, nuc1, dt, rdry(nin)
      REAL, INTENT(OUT) :: nnuccd, mnuccd, nnucci, mnucci
      REAL*8 :: sui

      cf = 3. !mineral dust
      alpha = 0.
      beta = 1.25
      lamda = 0.46
      delta = -11.6

      !rdry = rdry*(1.E-6) !converting um to m

      IF(temp.lt.268.15 .and.sui.ge.0.05.and.temp.gt.238.15)THEN

        IF (PIRE_CHEM) THEN
           dum = nuc1*rhodum*(1./(100.**3)) !converting nuc (#/kg) to !na(#/cm^3) 
        ELSE
           dum = nuc1*0.001 !converting nuc (#/L) to na (#/cm^3)
        END IF

        denom = 0.
        DO k = nin, 4, -1
           denom = denom + real(k)
        END DO

        mi = 0.
        DO k = 4,nin
           wght = real(k)/denom
           !mi = mi + 4./3.*pi*rhoi*((25.*1.e-6)**3)*wght
           mi = mi + 4./3.*pi*rhoi*((rdry(k)*1.e-6)**3)*wght
        END DO

        IF (demottflag .eq. 0) then !DeMott 2010
           nnuccd = (0.0000594*((273.16-temp)**(10./3.))*(dum**&
              (0.0264*(273.16-temp)+0.0033))*1000./rhodum)/dt!1/L*1000/rho ~1000/kg
           mnuccd = nnuccd*mi

!     DeMott 2015 !parameterization can be treated as immersion or condensation freezing
        ELSE IF (demottflag .eq. 1) then 
           nnucci = (cf*(dum**(alpha*(273.16-temp)+beta))*exp(lamda*&
              (273.16-temp)+delta)*1000./rhodum)/dt !1/L*1000/rho ~1000/kg
           mnucci = nnucci*mi 

        END IF

      END IF

      END SUBROUTINE DEMOTT

!_____________________________________________________________________
      REAL FUNCTION GAMMLN(XX)
!---------------------------------------------------------------------
      INTEGER J
      REAL*8 COF(6),STP,HALF,ONE,FPF,X,TMP,SER
      REAL XX
      DATA COF,STP/76.18009173D0,-86.50532033D0,24.01409822D0,&
          -1.231739516D0,.120858003D-2,-.536382D-5,2.50662827465D0/
      DATA HALF,ONE,FPF/0.5D0,1.0D0,5.5D0/
      X=XX-ONE
      TMP=X+FPF
      TMP=(X+HALF)*LOG(TMP)-TMP
      SER=ONE
      DO J=1,6
         X=X+ONE
         SER=SER+COF(J)/X
      ENDDO
      GAMMLN=TMP+LOG(STP*SER)
      RETURN
      END FUNCTION GAMMLN

!___________________________________________________________________
      REAL FUNCTION POLYSVP1 (T,type)
!--------------------------------------------------------------------

!  Compute saturation vapor pressure by using 
! function from Goff and Gatch (1946)

!  Polysvp returned in units of pa.
!  T is input in units of K.
!  type refers to saturation with respect to liquid (0) or ice (1)

      implicit none

      real dum
      real T,t0
      integer type

      t0=273.16

! ice


! ice

      if (type.eq.1) then
         
!     Goff Gatch equation (good down to -100 C)

            polysvp1 = 10.**(-9.09718*(273.16/t-1.)-3.56654*&
            alog10(273.16/t)+0.876793*(1.-t/273.16)+&
            alog10(6.1071))*100.

        
            
         
      end if
      
      
! liquid
      if (type.eq.0) then

!     Goff Gatch equation, uncertain below -70 C

         polysvp1 = 10.**(-7.90298*(373.16/t-1.)+&
          5.02808*alog10(373.16/t)- &
          1.3816e-7*(10**(11.344*(1.-t/373.16))-1.)+&
          8.1328e-3*(10**(-3.49149*(373.16/t-1.))-1.)+ &
          alog10(1013.246))*100. 

    
      end if
      
      end function polysvp1
      

      subroutine getmasssize(celsiushab,lmean,alpham,betam,&
      alpha_a,beta_a,alpha_v,beta_v,qi,ipart)


      implicit none
      real celsiushab,lmean,alpham,betam,&
      alpha_a,beta_a,alpha_v,beta_v,qi
      integer ipart

      if(lmean*1.e6.lt.125.0)then
         if(ipart.eq.0)then
            alpham = 1.23e-3    ! needles
            betam = 1.8
         else if(ipart.eq.1)then
            alpham = 110.8      ! columns
            betam = 2.91
         endif
         beta_a = 2.0
         alpha_a = 0.684*100.**(beta_a-2.)
      else
         if(ipart.eq.0)then
            alpham = 1.23e-3    ! needles
            betam = 1.8
         else if(ipart.eq.1)then
            alpham =  2.739e-3  ! columns
            betam = 1.74
         endif
         beta_a = 1.5
         alpha_a = 0.0696*100.**(beta_a-2.)
      endif

      if(ipart.eq.2.or.ipart.eq.3)then
         if(lmean*1.e6.le.100.)then
            betam = 2.91
            alpham = 0.167*100.**betam/1000.
            beta_a = 2.0
            alpha_a = 0.684*100.**(beta_a-2.)
         else if(lmean*1.e6.gt.100..and.lmean*1.e6.le.300.)then
            betam = 1.91
            alpham = 0.00166*100.**betam/1000.
            beta_a = 1.5
            alpha_a = 0.0696*100.**(beta_a-2.)
         else if(lmean*1.e6.gt.300.)then
            betam = 1.74
            alpham = 0.000907*100.**betam/1000.
            beta_a = 1.414
            alpha_a = 0.0512*100.**(beta_a-2.)
         endif
      endif

      if(celsiushab.lt.-10.0.and.celsiushab.gt.-20.0)then
         if(ipart.eq.0)then
            alpham = .377e-2   !dendrites
            betam = 2.
            if(lmean*1.e6.le.90.)then
               beta_a = 1.85
               alpha_a = 0.24*100.**(beta_a-2.)
            else
               beta_a = 1.63
               alpha_a = 0.11*100.**(beta_a-2.)
            endif
! Alex's data
!     alpha_m = 2.423e-1 !stellar dendrites
!     beta_m = 2.53
!     alpha_m = 2.327e-2  ! classic dendrites
!     beta_m = 2.29
         else if(ipart.eq.1)then
            alpham = .8854     ! hex plates
            betam =  2.5
            if(lmean*1.e6.le.90.)then
               beta_a = 1.85
               alpha_a = 0.24*100.**(beta_a-2.)
            else
               beta_a = 2.0
               alpha_a = 0.65*100.**(beta_a-2.)
            endif
         endif

         if(ipart.eq.2)then
            betam = 2.45
            alpham = 0.00739*100.**betam/1000. ! hex plates
            if(lmean*1.e6.le.90.)then
               beta_a = 1.85
               alpha_a = 0.24*100.**(beta_a-2.)
            else
               beta_a = 1.63
               alpha_a = 0.11*100.**(beta_a-2.)
            endif
         endif

         if(ipart.eq.3)then
            if(lmean*1.e6.le.90.)then ! dendrites
               betam = 2.42
               alpham = 0.00583*100.**betam/1000.
               beta_a = 1.85
               alpha_a = 0.24*100.**(beta_a-2.)
            else
               betam = 1.67
               alpham = 0.00027*100.**betam/1000.
               beta_a = 1.63
               alpha_a = 0.11*100.**(beta_a-2.)
            endif
         endif
      endif


      if(ipart.eq.4)then
      if(celsiushab.lt.0.0.and.celsiushab.ge.-3)then
         qi = 0.179*2.          !cap for a column
         betam = 3.0            ! columns
         alpham = 0.0450*1000.**betam/(1000.*1000.)
         beta_v = 0.650
         alpha_v = 0.825*1000.**beta_v
         if(lmean*1.e6.le.100.)then
            beta_a = 2.0
            alpha_a = 0.684*100.**(beta_a-2.)
         else if(lmean*1.e6.gt.100..and.lmean*1.e6.le.300.)then
            beta_a = 1.5
            alpha_a = 0.0696*100.**(beta_a-2.)
         else if(lmean*1.e6.gt.300.)then
            beta_a = 1.414
            alpha_a = 0.0512*100.**(beta_a-2.)
         endif
      else if(celsiushab.lt.-3.0.and.celsiushab.ge.-8.)then
         qi = .1803*2.
         betam = 2.01           ! needles
         alpham = 0.0092*1000.**betam/(1000.*1000.)
         beta_v = 0.590
         alpha_v = 0.444*1000.**beta_v
         if(lmean*1.e6.le.100.)then
            beta_a = 2.0
            alpha_a = 0.684*100.**(beta_a-2.)
         else if(lmean*1.e6.gt.100..and.lmean*1.e6.le.300.)then
            beta_a = 1.5
            alpha_a = 0.0696*100.**(beta_a-2.)
         else if(lmean*1.e6.gt.300.)then
            beta_a = 1.414
            alpha_a = 0.0512*100.**(beta_a-2.)
         endif
      else if(celsiushab.lt.-8.0.and.celsiushab.ge.-11.0)then
         qi = 0.429*2.
         betam = 1.9            ! sector, broad-branched
         alpham = 0.037*1000.**betam/(1000.*1000.)
         beta_v = 0.410
         alpha_v = 0.690*1000.**beta_v
         if(lmean*1.e6.le.90.)then
            beta_a = 1.85
            alpha_a = 0.24*100.**(beta_a-2.)
         else
            beta_a = 2.0
            alpha_a = 0.65*100.**(beta_a-2.)
         endif
      else if(celsiushab.lt.-11.0.and.celsiushab.ge.-17.0)then
         qi = .3183*2.
         betam = 2.19           ! dendrites
         alpham = 0.0141*1000.**betam/(1000.*1000.)
         beta_v = 0.217
         alpha_v = 0.376*1000.**beta_v
         if(lmean*1.e6.le.90.)then
            beta_a = 1.85
            alpha_a = 0.24*100.**(beta_a-2.)
         else
            beta_a = 1.63
            alpha_a = 0.11*100.**(beta_a-2.)
         endif
      else if(celsiushab.lt.-17.0.and.celsiushab.ge.-21.0)then
         qi = .3183*2.
         betam = 2.19           ! radiating assemblage of dendrites
         alpham = 0.0141*1000.**betam/(1000.*1000.)
         beta_v = 0.217
         alpha_v = 0.376*1000.**beta_v
         if(lmean*1.e6.le.90.)then
            beta_a = 1.85
            alpha_a = 0.24*100.**(beta_a-2.)
         else
            beta_a = 1.63
            alpha_a = 0.11*100.**(beta_a-2.)
         endif
      else if(celsiushab.lt.-21.0)then
         qi = 0.429*2.
         betam = 1.9            ! cold type
         alpham = 0.037*1000.**betam/(1000.*1000.)
         beta_v = 0.410
         alpha_v = 0.690*1000.**beta_v
         if(lmean*1.e6.le.90.)then
            beta_a = 1.85
            alpha_a = 0.24*100.**(beta_a-2.)
         else
            beta_a = 2.0
            alpha_a = 0.65*100.**(beta_a-2.)
         endif
      endif
      endif
      
      return
      END SUBROUTINE GETMASSSIZE

                                                            


!************************************************************************
! SUBROUTINE EVOLVE
!   Purpose: This routine solves the equations for the growth of 
!            non-spherical ice following the Chen and Lamb (1994)
!            method as per the document of Harrington and Sulia (2010)
!
!   Variables: (NOTE: All lengths are characteristic lengths!!!)
!           ani = initial a-axis length at time t
!           nu = shape of the distribution
!           ni = number density of ice crystals
!           sui = ice supersaturation
!           temp = temperaure
!           press = pressure
!           igr = inherent growth ratio
!           deltt = time-step 
!           iwc = final ice water content at time t+deltt
!           cf = final c-axis length at time t+deltt
!           iwci = initial ice water content at time t
!           phii = initial aspect ratio at time t
!           phif = final aspect ratio at time t+deltt
!           cni = initial c-axis length at time t
!           rni = initial equivalent volume radius at time t
!           rnf = final equivalent volume radius at time t+deltt
!           anf = final a-axis length at time t+deltt
!           deltastr = historical factor relating the a and c axes
!
!------------------------------------------------------------------------

      SUBROUTINE EVOLVE(ani,ni,sui,sup,qvv,temp,press,igr,deltt,iwc,cf,&
          iwci,phii,phif,cni,rni,rnf,anf,deltastr,mu,rhoavg,vtbarb,&
          vtbarbm,alphstr,vtbarblen,rhoa,i,k,iaspect,ISDACflag,&
          masssizeflag,ipart,sphrflag,redden,itimestep)
!***********************************************************************

      IMPLICIT NONE

      INTEGER ipart,ivent,i,k,iaspect,ISDACflag,redden,masssizeflag
      INTEGER sphrflag,itimestep,MICRO_ON
      REAL ani,ni,temp,press,qvv
      REAL igr,gi,capgam,dmdt,deltt,qe
      REAL dlndt,mu
      REAL l,lmean,fs,gammnu1,gammnubet,anf,betam
      REAL iwc,cf,alphstr,deltastr,iwci,phii,phif
      REAL cni,gammnu2delt,gammnu1delt,r,rf,vf
      REAL gamrats,alphan,rni,rnf,alphanr,rhoavg,rhodep,rhosub
      REAL drho,igr_den,vi,betavol,Vmin,videp,fallcheck
      REAL fv,dvs,kts,grav,rhoa,am,bm,bv1,bv2
      REAL bt1,bt2,bx,gv,gt,fh,ntherm,xvent,nre,nsch,npr,xm
      REAL bl,ba,aa,al,etaa,celsius,wghtv1
      REAL igrvent,wghtv2,wghtv3,vtbarb,vtbarbm,vtbarblen,vt
      REAL cv, bv, bm1
      REAL alpham,alpha_a,beta_a,alpha_v,beta_v,cap
      REAL alphams, alpha_as, alphap, alphamsp, betamp
      REAL*8 sui,sup,xn

      MICRO_ON = 1
      bv = 0.5
      cv = 12.*2.**bv
      bm1 = 3.0
      ivent=1
      IF(ISDACflag .EQ. 1) ivent=0
      grav = 9.81
      celsius = temp-273.15
      fv = 1.0
      fh = 1.0
      gi = findgtp(temp,press,fv,fh)
      dvs = 0.211*(temp/273.15)**1.94 * &
          (1013.25*100./press)*1.0/100.0**2.0 ! vapor diffusivity
      kts = 2.3823e-2 +7.1177e-5*(temp-273.15) ! thermal conductivity
!      if(MICRO_ON.eq.1)THEN
      l=ani
!      IF(masssizeflag.EQ.1)THEN
!         IF(cni.GT.ani)l=cni         
!      END IF
!      IF(MICRO_ON.eq.1)THEN
      drho = ((polysvp1(temp,0)-polysvp1(temp,1))/(Rv*temp))*1000.0
!      IF(MICRO_ON.eq.1)THEN
!      IGR_DEN = IGR

      betam = 2. + deltastr
      gammnu=exp(gammln(nu))
      gammnu1=exp(gammln(nu+1.0))
      gammnubet=exp(gammln(nu+betam))

      capgam = capacitance_gamma(l,deltastr,nu,alphstr)
      IF(ISDACflag .EQ. 1) capgam=capgam*2./pi
      lmean = l*nu

      cap = 0.
      alpham = 0.
      alphams = 0.
      alpha_as = 0.
      alphap = 0.
      alphamsp = 0.
      betamp = 0. 
      alpha_a = 0.
      beta_a = 0.
      alpha_v = 0.
      beta_v = 0.

!      IF(masssizeflag .EQ. 1)THEN

!         IF(ipart .eq. 0)THEN
!            capgam = 2.*.1803*l*gammnu1/gammnu ! cap for a needle in RAMS
            
!         ELSE IF(ipart .ge. 1)THEN
!            capgam = 2.*0.179*l*gammnu1/gammnu !cap for a column in RAMS
!         END IF
         
!         IF(celsius .lt. -10.0 .and. celsius .gt. -20.0)THEN
            
!            IF(ipart .eq. 0 .or. ipart .eq. 3)THEN
!               capgam = 2.*.3183*l*gammnu1/gammnu !cap for a dendrite
               
!            ELSE IF(ipart.eq.1.or.ipart.eq.2)THEN
!               capgam = 2.*0.429*l*gammnu1/gammnu ! cap for a plate in RAMS
!            END IF
!            
!         END IF  

!         cap=capgam/(l*gammnu1/gammnu)

!         CALL getmasssize(celsius,l,alpham,&
!         betam,alpha_a,beta_a,alpha_v,beta_v,cap,ipart)

!         alphams=alpham*2.**betam
!         alpha_as=alpha_a*2.**beta_a
!         IF(alphap.eq.0.0)THEN
!            alphap=alphams
!            alphamsp=alphams
!            betamp=betam
!         END IF
!     deltastr=betam-2.
!     alphstr=ao**(1.-deltastr)
!         capgam=cap   
!         gammnubet=exp(gammln(nu+betam))
!         iwci=ni*alphams*l**betam*gammnubet/gammnu
         
!      END IF

!      dmdt = 4.*pi*gi*sui*ni*capgam

        !     get alpham and betam from deltastr
!      dlndt = dmdt/(ni*lenconv(ani,deltastr,nu,alpham,&
!      betam,alphstr,rhoavg))


      gammnubet = exp(gammln(nu+betam))
      gammnu2delt = exp(gammln(nu+2.+deltastr))
      gammnu1delt = exp(gammln(nu+deltastr-1.))
      gamrats = exp(gammln(nu+betam/3.))/gammnu
      
      phii = cni/ani*gammnu1delt/gammnu ! diagnose initial aspect ratio
      IF(iaspect .eq. 1) phii = 0.27
      
      r = alphstr**(1./3.) * (l)**(betam/3.)
!      fs = capgam/(r*gamrats)
      fs = capgam/(rni*gamrats)
      alphan = cni*ani**(-igr)
      alphanr = ani/rni**(3./(2.+igr))

      if(phii.lt.1.0)then
         bl = 1.
         al = 2.
         aa = pi
         ba = 2.
         qe = (rhoavg/920.)**(2./3.)
      else if(phii .gt. 1.0)then
         al = 2.*alphstr**(1./3.)
         bl = (deltastr + 2.)/3.
         aa = pi*alphstr
         ba = deltastr + 1.
         qe = 1.0
      else if(phii .eq. 1.0)then
         bl = 1.
         al = 2.
         aa = pi
         ba = 2.
         qe = 1.0
      endif
      
      etaa = (1.718 + 0.0049*celsius - 1.2e-5*celsius**2)*1.e-4 &
      *100./1000.               ! dynamic viscisoty
      nsch = etaa/(rhoa*dvs)    ! Schmidt number
      npr = etaa/(rhoa*kts)     ! Prandlt number

!      xn =  2./rhoavg*(rhoavg-rhoa)*grav*rhoa/etaa**2 & 
!     *(4./3.*pi*rhoavg)*alphstr*al**2/(aa*qe) * qe**(3./4.)
      xn =  2./rhoavg*(rhoavg-rhoa)*grav*rhoa/etaa**2 & 
      *(4./3.*pi*rhoavg)*10.**(alog10(alphstr)+2.* &
      alog10(al)-alog10(aa))/qe * qe**(3./4.)
      
      bx = deltastr+2.+2.*bl-ba
      
!      IF(masssizeflag.eq.1)THEN
!         xn = (8.*alphams*grav*rhoa)/(alpha_as*etaa**2)
!         bx = betam + 2.0 - beta_a
!      END IF
      wghtv1 = exp(gammln(nu+deltastr+2.+2.*bl-ba))/gammnu
      xm = xn*ani**bx * wghtv1  ! average Best number
      
      if(xm.le.10.)then
         am = 0.04394
         bm = 0.970
      else if(xm.gt.10.0.and.xm.le.585.)then
         am = 0.06049
         bm = 0.831
      else if(xm.gt.585.and.xm.le.1.56e5)then
         am = 0.2072
         bm = 0.638
      else if(xm.gt.1.56e5.and.xm.lt.1.e8)then
         am = 1.0865
         bm = 0.499
      else
         am = 1.0865
         bm = 0.499

      endif
      
      Nre = am*xm**bm
      
      wghtv2 = exp(gammln(nu+bx*bm-1.))/gammnu
      vtbarb = etaa/rhoa*0.5 * am*(xn)**bm*ani**(bx*bm-1.) * wghtv2
      wghtv3 = exp(gammln(nu+bx*bm-1.+2.+deltastr))&
      /exp(gammln(nu+2.+deltastr))
      vtbarbm = etaa/rhoa*0.5 * am*(xn)**bm*ani**(bx*bm-1.) * wghtv3

 
! add length
      if(phii.lt.1.0)then
         wghtv3 = exp(gammln(nu+bx*bm-1.+1)) &
         /exp(gammln(nu+1.))
         vtbarblen = etaa/rhoa*0.5 * am*(xn)**bm*ani**(bx*bm-1.) &
         * wghtv3
      else if(phii.gt.1.0)then
         wghtv3 = exp(gammln(nu+bx*bm-1.+deltastr))&
         /exp(gammln(nu+deltastr))
         vtbarblen = etaa/rhoa*0.5 * am*(xn)**bm*ani**(bx*bm-1.) &
         * wghtv3
      else if (phii.eq.1.0)then
         wghtv3 = exp(gammln(nu+bx*bm-1.+1))&
         /exp(gammln(nu+1.))
         vtbarblen = etaa/rhoa*0.5 * am*(xn)**bm*ani**(bx*bm-1.) &
         * wghtv3
      endif
      
      IF(ISDACflag .eq. 1)THEN
         vtbarb = cv*ani**bv*exp(gammln(nu+bv))/gammnu
         vtbarbm = cv*ani**bv*exp(gammln(nu+bv+1.))/exp(gammln(nu+1.))
         vtbarblen = cv*ani**bv*exp(gammln(nu+bv+bm1))/&
         exp(gammln(nu+bm1))
      END IF
      
!      IF(masssizeflag .EQ. 1)THEN
!         vtbarb = alpha_v*2.**beta_v*l**beta_v*&
!         exp(gammln(nu+beta_v))/gammnu
!         vtbarbm = alpha_v*2.**beta_v*l**beta_v*&
!         exp(gammln(nu+beta_v+betam))/exp(gammln(nu+betam))
!         vtbarblen = alpha_v*2.**beta_v*l**beta_v*&
!         exp(gammln(nu+beta_v+1.))/exp(gammln(nu+1.))
!      END IF

         
      xvent = nsch**(1./3.)*SQRT(nre)
      ntherm = SQRT(nre)*npr**(1./3.)

      if(xvent.le.1.0)then
         bv1 = 1.0
         bv2 = 0.14
         gv = 2.
      else
         bv1 = 0.86
         bv2 = 0.28
         gv = 1.
      endif

      if(ntherm.lt.1.4)then
         bt1 = 1.0
         bt2 = 0.108
         gt = 2.0
      else
         bt1 = 0.78
         bt2 = 0.308
         gt = 1.0
      endif

      if(ivent.eq.1)then
         fv = bv1 + bv2*xvent**gv
!         fva = bv1 + bv2*xvent**gv*(ani/rni)**(gv/2.)
!         fva = bv1 + bv2*xvent**gv*(alphstr**(-1./3.)*
!     +        ln**(1.-(deltastr+2.)/3.) * 
!     +        exp(gammln(nu + 1. -(deltastr+2.)/3.))/gammnu)**(gv/2.)
!         fvc = bv1 + bv2*xvent**gv*(cni/rni)**(gv/2.)
!         fvc = bv1 + bv2*xvent**gv*(alphstr**(2./3.)*
!     +        ln**(deltastr-(deltastr+2.)/3.) * 
!     +        exp(gammln(nu + deltastr -(deltastr+2.)/3.))/gammnu)
!     +        **(gv/2.)
         fh = bt1 + bt2*ntherm**gt

         igrvent = igr
      else
         fv = 1.0
         fh = 1.0
         igrvent = igr
      endif

      gi = findgtp(temp,press,fv,fh)
      

!      rhodep = 920.*exp((-3.0*max(drho-0.05,0.0))/IGR) ! get density of ice deposited
 
!      rhodep = 920.                 
                                             ! during growth
      vt = mu/rhoa*0.5*Nre/ani
      fallcheck = sqrt(dvs*pi*2*cni/(vt*nu))
      rhosub = rhoi*((1.-sui/qvv)+igr*(sui/qvv))
      IF(igr .lt. 1.0)THEN
         IF(sup .ge. 0.0)THEN
            IF(ani .gt. fallcheck)THEN
                rhodep = rhoi*igr
            ELSE
               rhodep = rhoi
            END IF
         ELSE
            IF(ani .gt. fallcheck)THEN
               rhodep = rhosub
            ELSE
               rhodep = rhoi
            END IF
         END IF
      ELSE IF(igr.gt.1.0)THEN 
         IF(sup .ge. 0.0)THEN
            rhodep = rhoi/igr
         ELSE
            rhodep = rhosub
         END IF
      ELSE
         rhodep = 920.
      END IF
      if(sui.lt.0.0)then        ! OR get density of ice removed 
         rhodep = rhoavg        ! during sublimation using polynomial
         videp = 4./3.*pi*rni**3 & 
         *exp(gammln(nu+deltastr+2.))/gammnu
         Vmin = 4./3.*pi*ao**3
         if(Vmin.lt.videp)then
            betavol = alog(rhoi/rhoavg)*1./(alog(Vmin/videp))
            rhodep = rhoavg*(1.+betavol)
         else
            rhodep = rhoavg 
         endif
      endif


      IF(ISDACflag.eq.1) rhodep = 88.4
!      IF(masssizeflag .eq. 1.or.sphrflag.eq.1) rhodep=920.
      If(redden .eq. 1) rhodep = 500.
      !rhodep=500.
!     rf = ((rni*gamrats)**2 + 2.*gi*sui*fs/rhodep &
!           *gammnu/gammnubet*gamrats**3*deltt)**(0.5) ! rmean after deltat

      rf = (max(((rni*gamrats)**2 + 2.*gi*sui*fs/rhodep &
           *gammnu/gammnubet*gamrats**3*deltt),1.e-20))**(0.5) ! rmean after deltat
      rnf = rf/gamrats           ! convert to rn (characteristic r)
      anf = alphanr*rnf**(3./(2.+igr)) ! characteristic a-axis after deltat
      vi = 4./3.*pi*rni**3 * exp(gammln(nu+deltastr+2.))/gammnu ! mean initial volume
      vf = 4./3.*pi*rnf**3 * exp(gammln(nu+deltastr+2.))/gammnu ! mean volume after deltat
      rhoavg = rhoavg*(vi/vf) + rhodep*(1.-vi/vf) ! new average ice density
      rhoavg = max(min(rhoavg,920.),50.)

      IF(ISDACflag .eq.1) rhoavg = 88.4
!      IF(masssizeflag.eq.1.or.sphrflag.eq.1)rhoavg=920.
      !rhoavg=500.
      iwc = ni*rhoavg*vf        ! IWC after deltat
!      IF(masssizeflag.eq.1) iwc=ni*alpham*rnf**betam*&
!      (gammnu/exp(gammln(nu+deltastr+2.)))**(3./betam)

      if(igr.ne.1.0.and.(log(anf)-log(ao)).gt.0.01 .and.&
      (3.*log(rnf)-2.*log(anf)-log(ao)).gt.0.001)THEN
         deltastr = (3.*log(rnf)-2.*log(anf)-log(ao)) &
         /(log(anf)-log(ao))    ! diagnose new deltastar
         IF(masssizeflag.eq.1)deltastr=1.0
      else
         deltastr = 1.
      endif
      IF(iaspect .eq. 1) deltastr = 0.8
      deltastr = min(max(deltastr,0.55),1.5)
  
      if(deltastr.ge.1.0)then   ! if columns, get c from c-a relation
         cf = (ao**(1.-deltastr))*anf**deltastr
      endif

      phif = phii*(rnf**3/r**3)**((igr-1.)/(igr+2.)) ! if plates, get c from aspect ratio after deltat
      IF(iaspect .eq.1) phif = 0.27

!      rnf = rf
      if(deltastr.lt.1.0)then
         cf = phif*anf*gammnu/exp(gammln(nu+deltastr-1.))
      endif

      IF(masssizeflag.eq.1)THEN
         anf = (ani**(betam-1.0) + 4.*pi*gi*sui*cap*(betam-1.)/&
         (alphams*betam)*exp(gammln(nu+1.))/&
         exp(gammln(nu+betam)) * deltt)**(1./(betam-1.))
         rnf = anf
!     iwc = ni*alphams*lf**betam * gammnubet/gammnu
         iwc = ni*alphams*anf**betam * gammnubet/gammnu
!     dlndt = dmdt/(ni*lenconv(ln,deltastr,nu,alpham,betam,alphstr))
      END IF
!      END IF !MICRO_ON

 
      return
      END SUBROUTINE EVOLVE




!------------------------------------------------------------------
!     FUNCTION FINDGTP
!     
!     This function calculates Gi(T,P) for ice growth based on the
!     function in Pruppacher and Klett (1997)
!     
!------------------------------------------------------------------
      
      real function findgtp(temp,p0,fv,fh)
      
      implicit none
      real temp,p0,xdvstar,xk,eeq,esi_new,rv,L,gtp1,gtp2,dvstar
      real xkstar,dv,fv,fh
!      real polysvp1
      
      L = FLHS(temp)            !2.833e6               ! latent heat of sublimation
      dv = 0.211*(temp/273.15)**1.94 *(101300.25/p0)*1.0/100.0**2.0 ! vapor diffusivity
      xdvstar = dv*fv              ! vapor diffusivity in cm^2/s
!     xk = xkstar              ! thermal conductivity in W/m/K
      xk = (2.3823e-2 +7.1177e-5*(temp-273.15))*fh ! thermal conductivity
      eeq =  polysvp1(temp,1)           ! equilibrium vp
      rv = 461.5                ! H2O gas constant
!      WRITE(*,*) 'T', temp, 'P', p0, 'Ls', L, 'Kt', xk, 'Dv', dv, 'ei', 
!     +     eeq
      if(xdvstar.eq.0.)then
         findgtp=0.0
      else
         gtp1 = rv*temp/(xdvstar*eeq)
         gtp2 = (L**2.0/(xk*rv*temp**2.0)) - (L/(xk*temp))
         findgtp = 1.0/(gtp1+gtp2)
      endif
      return
      end function findgtp

!------------------------------------------------------------------------------

      real function capacitance_gamma(axis,deltastr,nu,alpha_s)
      
      implicit none
      real l,phi,a,c,axis,deltastr
      real a1,a2,b1,capone,c1,c2,d1,d2,nu
      real li,phicoli,phiplai,bcol,bpla
      real acol,apla,b2,alpha_s,beta_s
      real gammnu,gammad1,gammad2!,gammln

      a1 = 0.58
      a2 = 0.95
      b1 = 0.75

      l = axis

      li = ao/nu
      phicoli = 1.0
      phiplai = 1.0
      bcol = 1./deltastr
      bpla = deltastr
      acol = 1./phicoli * li**(1.-bcol)
      apla = phiplai * li**(1.-bpla)

!    Oblate Spheroid
      if (deltastr.le.1.) then
       a1 = 0.6369427        !0.58/0.9105999
       b1 = 0.0
       a2 = 0.57*a1   ! fit to spheroid
       b2 = 0.95   ! fit to spheroid
       alpha_s = apla
       beta_s = bpla-1.0

!     Prolate Spheroid
      elseif (deltastr.gt.1.) then
       a1 = 0.5714285   ! 0.58/1.015000
       b1 = -1.0
       a2 = 0.75*a1      ! fit to spheroid
       b2 = 0.82-1.0  ! fit to spheroid 
       alpha_s = apla
       beta_s = bpla-1.0

      endif

 
      if (deltastr.eq.1..and.nu.eq.1.) then
       gammnu = exp(gammln(nu))
       gammad1 = exp(gammln(nu+1.0))
       capacitance_gamma = l*gammad1/gammnu
       
      else if(deltastr.le.1.)then
       c1 = a1*alpha_s**b1
       c2 = a2*alpha_s**b2
       d1 = beta_s*b1+1.
       d2 = beta_s*b2+1.

       gammnu = exp(gammln(nu))
       gammad1 = exp(gammln(nu+d1))
       gammad2 = exp(gammln(nu+d2))
       capacitance_gamma = c1*l**d1 * gammad1/gammnu &
          + c2*l**d2 * gammad2/gammnu


      else if(deltastr.gt.1.) then
       c1 = a1*alpha_s**(b1+1.)
       c2 = a2*alpha_s**(b2+1.)
       d1 = beta_s*b1+deltastr
       d2 = beta_s*b2+deltastr

       gammnu = exp(gammln(nu))
       gammad1 = exp(gammln(nu+d1))
       gammad2 = exp(gammln(nu+d2))
       capacitance_gamma = c1*l**d1 * gammad1/gammnu &
          + c2*l**d2 * gammad2/gammnu

      endif

      return
      end function capacitance_gamma

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real function lenconv(axis,deltastr,nu,alpha_m,beta_m,alphstr,&
      rhoi)

      implicit none
      real l,phi,a,c,axis
      real a1,a2,b1,capone,c1,c2,d1,d2
      real deltastr,li,phicoli,phiplai,bcol,bpla
      real acol,apla,b2!,pi
      real gammnu,gammabeta,alpha_m,beta_m
      real alphstr,nu,rhoi

! NOTE: USE SEMI-MAJOR AXIS NOT MAJOR AXIS LENGTH!!

!      pi = 3.14159

      
      l = axis

      li = ao/nu
      phicoli = 1.0
      phiplai = 1.0
      bcol = 1./deltastr
      bpla = deltastr
      acol = 1./phicoli * li**(1.-bcol)
      apla = phiplai * li**(1.-bpla)

      if(deltastr.ne.1.) then
       beta_m = 2.+bpla
       alpha_m = 4./3.*pi*apla*rhoi
       alphstr = apla

      elseif (deltastr.eq.1.)then
       l = axis
       beta_m = 3.0
       alpha_m = 4./3.*pi*rhoi
       alphstr = 1.0
      endif

      gammnu = exp(gammln(nu))
      gammabeta = exp(gammln(nu+beta_m))
      lenconv = alpha_m*beta_m*l**(beta_m-1.)*gammabeta/gammnu

      return
      end function lenconv


!**********************************************************************
!     This function calculates the enthalpy of sublimation for       
!     water as per Bohren and Albrecht (1999, pg 197, eq 5.64). 
!     Arguement is in kelvin and the returned result is in mks 
!     Coded: JYH Feb 2005    
!**********************************************************************
      REAL FUNCTION FLHS (tmp)

      REAL tmp,cpv,ci,t0,xlvs0,delcv,delt,xlvs
        
      cpv = 1850.               ! specific heat of vapor const p
      ci = 2106.                ! specific heat of ice
      t0 = 273.15
      xlvs0 = 2.834e6           ! enthalpy of sublimation at 0 C 
      delcv = cpv - ci
      delt = t0 - tmp
      xlvs = (xlvs0 + delcv*delt) ! B&A equation 5.64
      FLHS = xlvs
  
      RETURN
      END FUNCTION FLHS


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      REAL FUNCTION GAMMA(X)

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!CD    DOUBLE PRECISION FUNCTION DGAMMA(X)
!C----------------------------------------------------------------------
!C
!C THIS ROUTINE CALCULATES THE GAMMA FUNCTION FOR A REAL ARGUMENT X.
!C   COMPUTATION IS BASED ON AN ALGORITHM OUTLINED IN REFERENCE 1.
!C   THE PROGRAM USES RATIONAL FUNCTIONS THAT APPROXIMATE THE GAMMA
!C   FUNCTION TO AT LEAST 20 SIGNIFICANT DECIMAL DIGITS.  COEFFICIENTS
!C   FOR THE APPROXIMATION OVER THE INTERVAL (1,2) ARE UNPULBISHED.
!C   THOSE FOR THE APPROXIMATION FOR X .GE. 12 ARE FROM REFERENCE 2.
!C   THE ACCURACY ACHIEVED DEPENDS ON THE ARITHMETIC SYSTEM, THE
!C   COMPILER, THE INTRINSIC FUNCTIONS, AND PROPER SELECTION OF THE
!C   MACHINE-DEPENDENT CONSTANTS.
!C
!C
!C*******************************************************************
!C*******************************************************************
!C
!C EXPLANATION OF MACHINE-DEPENDENT CONSTANTS
!C
!C BETA   - RADIX FOR THE FLOATING-POINT REPRESENTATION
!C MAXEXP - THE SMALLEST POSITIVE POWER OF BETA THAT OVERFLOWS
!C XBIG   - THE LARGEST ARGUMENT FOR WHICH GAMMA(X) IS REPRESENTABLE
!C          IN THE MACHINE, I.E., THE SOLUTION TO THE EQUATION
!C                  GAMMA(XBIG) = BETA**MAXEXP
!C XINF   - THE LARGEST MACHINE REPRESENTABLE FLOATING-POINT NUMBER;
!C          APPROXIMATELY BETA**MAXEXP
!C EPS    - THE SMALLEST POSITIVE FLOATING-POINT NUMBER SUCH THAT
!C          1.0+EPS .GT. 1.0
!C XMININ - THE SMALLEST POSITIVE FLOATING-POINT NUMBER SUCH THAT
!C          1/XMININ IS MACHINE REPRESENTABLE
!C
!C     APPROXIMATE VALUES FOR SOME IMPORTANT MACHINES ARE:
!C
!C                            BETA       MAXEXP        XBIG
!C
!C CRAY-1         (S.P.)        2         8191        966.961
!C CYBER 180/855
!C   UNDER NOS    (S.P.)        2         1070        177.803
!C IEEE (IBM/XT,
!C   SUN, ETC.)   (S.P.)        2          128        35.040
!C IEEE (IBM/XT,
!C   SUN, ETC.)   (D.P.)        2         1024        171.624
!C IBM 3033       (D.P.)       16           63        57.574
!C VAX D-FORMAT   (D.P.)        2          127        34.844
!C VAX G-FORMAT   (D.P.)        2         1023        171.489
!C
!C                            XINF         EPS        XMININ
!C
!C CRAY-1         (S.P.)   5.45E+2465   7.11E-15    1.84E-2466
!C CYBER 180/855
!C   UNDER NOS    (S.P.)   1.26E+322    3.55E-15    3.14E-294
!C IEEE (IBM/XT,
!C   SUN, ETC.)   (S.P.)   3.40E+38     1.19E-7     1.18E-38
!C IEEE (IBM/XT,
!C   SUN, ETC.)   (D.P.)   1.79D+308    2.22D-16    2.23D-308
!C IBM 3033       (D.P.)   7.23D+75     2.22D-16    1.39D-76
!C VAX D-FORMAT   (D.P.)   1.70D+38     1.39D-17    5.88D-39
!C VAX G-FORMAT   (D.P.)   8.98D+307    1.11D-16    1.12D-308
!C
!C*******************************************************************
!C*******************************************************************
!C
!C ERROR RETURNS
!C
!C  THE PROGRAM RETURNS THE VALUE XINF FOR SINGULARITIES OR
!C     WHEN OVERFLOW WOULD OCCUR.  THE COMPUTATION IS BELIEVED
!C     TO BE FREE OF UNDERFLOW AND OVERFLOW.
!C
!C
!C  INTRINSIC FUNCTIONS REQUIRED ARE:
!C
!C     INT, DBLE, EXP, LOG, REAL, SIN
!C
!C
!C REFERENCES:  AN OVERVIEW OF SOFTWARE DEVELOPMENT FOR SPECIAL
!C              FUNCTIONS   W. J. CODY, LECTURE NOTES IN MATHEMATICS,
!C              506, NUMERICAL ANALYSIS DUNDEE, 1975, G. A. WATSON
!C              (ED.), SPRINGER VERLAG, BERLIN, 1976.
!C
!C              COMPUTER APPROXIMATIONS, HART, ET. AL., WILEY AND
!C              SONS, NEW YORK, 1968.
!C
!C  LATEST MODIFICATION: OCTOBER 12, 1989
!C
!C  AUTHORS: W. J. CODY AND L. STOLTZ
!C           APPLIED MATHEMATICS DIVISION
!C           ARGONNE NATIONAL LABORATORY
!C           ARGONNE, IL 60439
!C
!C----------------------------------------------------------------------
      implicit none
      integer :: I,N
      logical :: l_parity
      real ::                                                       &
          CONV,EPS,FACT,HALF,ONE,res,sum,TWELVE,                    &
          TWO,X,XBIG,XDEN,XINF,XMININ,XNUM,Y,Y1,YSQ,Z,ZERO
      real, dimension(7) :: C
      real, dimension(8) :: P
      real, dimension(8) :: Q
      real, parameter    :: constant1 = 0.9189385332046727417803297

!----------------------------------------------------------------------
!  MATHEMATICAL CONSTANTS
!----------------------------------------------------------------------
      data ONE,HALF,TWELVE,TWO,ZERO/1.0E0,0.5E0,12.0E0,2.0E0,0.0E0/
!----------------------------------------------------------------------
!  MACHINE DEPENDENT PARAMETERS
!----------------------------------------------------------------------
      data XBIG,XMININ,EPS/35.040E0,1.18E-38,1.19E-7/,XINF/3.4E38/
!----------------------------------------------------------------------
!  NUMERATOR AND DENOMINATOR COEFFICIENTS FOR RATIONAL MINIMAX
!     APPROXIMATION OVER (1,2).
!----------------------------------------------------------------------
      data P/-1.71618513886549492533811E+0,2.47656508055759199108314E+1,&
             -3.79804256470945635097577E+2,6.29331155312818442661052E+2,&
             8.66966202790413211295064E+2,-3.14512729688483675254357E+4,&
             -3.61444134186911729807069E+4,6.64561438202405440627855E+4/
      data Q/-3.08402300119738975254353E+1,3.15350626979604161529144E+2,&
             -1.01515636749021914166146E+3,-3.10777167157231109440444E+3,&
              2.25381184209801510330112E+4,4.75584627752788110767815E+3,&
            -1.34659959864969306392456E+5,-1.15132259675553483497211E+5/
!----------------------------------------------------------------------
!  COEFFICIENTS FOR MINIMAX APPROXIMATION OVER (12, INF).
!----------------------------------------------------------------------
      data C/-1.910444077728E-03,8.4171387781295E-04,&
           -5.952379913043012E-04,7.93650793500350248E-04,&
           -2.777777777777681622553E-03,8.333333333333333331554247E-02,&
            5.7083835261E-03/
!----------------------------------------------------------------------
!  STATEMENT FUNCTIONS FOR CONVERSION BETWEEN INTEGER AND FLOAT
!----------------------------------------------------------------------
      CONV(I) = REAL(I)
      l_parity=.FALSE.
      FACT=ONE
      N=0
      Y=X
      if (Y.LE.ZERO) then
!----------------------------------------------------------------------
!  ARGUMENT IS NEGATIVE
!----------------------------------------------------------------------
        Y=-X
        Y1=AINT(Y)
        res=Y-Y1
        if (res.NE.ZERO) then
          if(Y1.NE.AINT(Y1*HALF)*TWO)l_parity=.TRUE.
          FACT=-PI/SIN(PI*res)
          Y=Y+ONE
        else
          res=XINF
          goto 900
        endif
      endif
!----------------------------------------------------------------------
!  ARGUMENT IS POSITIVE
!----------------------------------------------------------------------
      if (Y.LT.EPS) then
!----------------------------------------------------------------------
!  ARGUMENT .LT. EPS
!----------------------------------------------------------------------
        if (Y.GE.XMININ) then
          res=ONE/Y
        else
          res=XINF
          goto 900
        endif
      elseif (Y.LT.TWELVE) then
        Y1=Y
        if (Y.LT.ONE) then
!----------------------------------------------------------------------
!  0.0 .LT. ARGUMENT .LT. 1.0
!----------------------------------------------------------------------
          Z=Y
          Y=Y+ONE
        else
!----------------------------------------------------------------------
!  1.0 .LT. ARGUMENT .LT. 12.0, REDUCE ARGUMENT IF NECESSARY
!----------------------------------------------------------------------
          N=INT(Y)-1
          Y=Y-CONV(N)
          Z=Y-ONE
        endif
!----------------------------------------------------------------------
!  EVALUATE APPROXIMATION FOR 1.0 .LT. ARGUMENT .LT. 2.0
!----------------------------------------------------------------------
        XNUM=ZERO
        XDEN=ONE
        do I=1,8
          XNUM=(XNUM+P(I))*Z
          XDEN=XDEN*Z+Q(I)
        enddo
        res=XNUM/XDEN+ONE
        if (Y1.LT.Y) then
!----------------------------------------------------------------------
!  ADJUST RESULT FOR CASE  0.0 .LT. ARGUMENT .LT. 1.0
!----------------------------------------------------------------------
          res=res/Y1
        elseif (Y1.GT.Y) then
!----------------------------------------------------------------------
!  ADJUST RESULT FOR CASE  2.0 .LT. ARGUMENT .LT. 12.0
!----------------------------------------------------------------------
          do I=1,N
            res=res*Y
            Y=Y+ONE
          enddo
        endif
      else
!----------------------------------------------------------------------
!  EVALUATE FOR ARGUMENT .GE. 12.0,
!----------------------------------------------------------------------
        if (Y.LE.XBIG) then
          YSQ=Y*Y
          sum=C(7)
          do I=1,6
            sum=sum/YSQ+C(I)
          enddo
          sum=sum/Y-Y+constant1
          sum=sum+(Y-HALF)*log(Y)
          res=exp(sum)
        else
          res=XINF
          goto 900
        endif
      endif
!----------------------------------------------------------------------
!  FINAL ADJUSTMENTS AND RETURN
!----------------------------------------------------------------------
      if (l_parity)res=-res
      if (FACT.NE.ONE)res=FACT/res
  900 gamma=res
      return

! ---------- LAST LINE OF GAMMA ----------
      END FUNCTION GAMMA


!      LOGICAL FUNCTION  wrf_dm_on_monitor()
!      IMPLICIT NONE
!#ifndef STUBMPI
!      INCLUDE 'mpif.h'
!      INTEGER tsk, ierr, mpi_comm_local
!      CALL wrf_get_dm_communicator( mpi_comm_local )
!      CALL mpi_comm_rank ( mpi_comm_local, tsk , ierr )
!      wrf_dm_on_monitor = tsk .EQ. 0
!#else
!      wrf_dm_on_monitor = .TRUE.
!#endif
!      RETURN
!      END FUNCTION wrf_dm_on_monitor



      END MODULE MODULE_MP_SULIAHARRINGTON

!
