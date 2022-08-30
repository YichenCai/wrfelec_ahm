      MODULE MODULE_MP_AHM
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

      PUBLIC  ::  MP_AHM
!     PUBLIC  ::  SULIAHARRINGTON_INIT
      PRIVATE  ::  POLYSVP1
      PRIVATE :: GAMMA
      PRIVATE :: PI, SQRTPI
      PRIVATE :: GAMMLN
!     PRIVATE :: EVOLVE
      PRIVATE :: FINDGTP
      PRIVATE :: LENCONV
      PRIVATE :: FLHS
      PRIVATE :: GET_X
      PRIVATE :: CAPACITANCE_GAMMA
      PRIVATE :: AHM_MICRO
      PRIVATE :: DELBK, DELABK, gamma_dp

      REAL, PRIVATE ::  rhoi    !BULK DENSITY OF CLOUD ICE
      REAL, PRIVATE ::  nu      !DISTRIBUTION SHAPE FACTOR
      REAL, PRIVATE ::  nuc     !ICE NUCLEATION CONCENTRAION (#/L)
      REAL, PRIVATE ::  rd      !GAS CONSTANT OF DRY AIR
      REAL, PRIVATE ::  cp      !SPECIFIC HEAT FOR DRY AIR (CONST P)
      REAL, PRIVATE ::  ao      !INITIAL CHARACTERISTIC A-AXIS LENGTH
      REAL, PRIVATE ::  co      !INITIAL CHARACTERISTIC C-AXIS LENGTH
      REAL, PRIVATE ::  li0     !INITIAL SEMI-MAJOR AXIS LENGTH
      REAL, PRIVATE ::  mi0     !INITIAL PARTICLE MASS
      !REAL, PRIVATE ::  gammnu  !GAMMA DIST WITH SHAPE, NU
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
      REAL, PRIVATE :: ars       ! FALL SPEED PARAMETER FOR SNOW
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

      REAL, PRIVATE :: coll_ni(5), coll_an(6), coll_cn(6), coll_nu(8), coll_rho(9)
      REAL, PRIVATE :: coll(5,6,6,8,9), ncoll(5,6,6,8,9)
      REAL, PRIVATE :: coll_phi(20), coll_r(27)
      REAL, PRIVATE :: acoll(20,27), ancoll(20,27), cncoll(20,27), ddcoll(20,27)
      INTEGER, PRIVATE :: ii, jj, kk, ll, mm, nn, oo, iii, jjj, kkk, lll, mmm, nnn, ooo

      !---------------------------------------------------------------------------
      !Explicit lightning public vars (add with caution)--------------------------
      integer :: ipelec = 0
      integer :: isaund = 0
      logical :: idonic = .false.

      ! For "an" indices
      integer :: na = 9 ! to be changed later
      integer :: ltmp
      integer, parameter :: lqmx = 30
      integer, parameter :: lt = 1
      integer, parameter :: lv = 2
      integer, parameter :: lc = 3
      integer, parameter :: lr = 4
      integer, parameter :: li = 5
      !   integer, private :: lis = 0
      integer, private :: ls = 6
      integer, private :: lh = 7
      !   integer, private :: lhl = 0
      integer, private :: lhab = 7
      integer, private :: lg = 7
      integer, private  :: lccn = 0 ! 0 or 9, other indices adjusted accordingly
      integer, private :: lccna = 0
      integer, private :: lcina = 0
      integer, private :: lcin = 0
      integer, private :: lnc = 0
      integer, private :: lnr = 0
      integer, private :: lni = 0
      !   integer, private :: lnis = 0
      integer, private :: lns = 0
      integer, private :: lnh = 0
      !   integer, private :: lnhl = 0
      ! Space charge
      integer :: lscw = 0
      integer :: lscr = 0
      integer :: lsci = 0
      !   integer :: lscis = 0
      integer :: lscs = 0
      integer :: lsch = 0
      !   integer :: lschl = 0
      !   integer :: lscwi = 0
      integer :: lscpi = 0
      integer :: lscni = 0
      integer :: lscpli = 0
      integer :: lscnli = 0
      integer :: lne = 0 ! last varible for transforming

      integer lsc(lc:lqmx)
      integer ln(lc:lqmx)
      integer denscale(lc:lqmx) ! flag for density scaling (mixing ratio conversion)
      real qxmin(lc:lqmx)
      ! End "an" indices

      integer :: nonigrd = 0

      ! Takahashi lookup table
      integer, parameter :: nlwc=30, ntem=31
      real takalu(0:ntem,0:nlwc)

      ! fundamental unit of charge
      real ec,ec_i
      parameter (ec = 1.602e-19) ! Coulomb / proton or electron
      parameter (ec_i = 1.0/ec) ! #proton or electron / Coulomb
      !---------------------------------------------------------------------------
      !End explicit lightning ----------------------------------------------------
CONTAINS
!*************************************************************************************
!     THIS SUBROUTINE INITIALIZES ALL PHYSICAL CONSTANTS AMND PARAMETERS 
!     NEEDED BY THE MICROPHYSICS SCHEME.
!     NEEDS TO BE CALLED AT FIRST TIMESTEP,
!     PRIOR TO CALL TO MAIN MICROPHYSICS INTERFACE
!-------------------------------------------------------------------------------------
      SUBROUTINE AHM_INIT(nssl_params,idonictmp)
!************************************************************************************

      IMPLICIT NONE
      !---------------------------------------------------------------------------
      !Explicit lightning --------------------------------------------------------
      real, intent(in), dimension(20) :: nssl_params ! for lightning
      logical, optional, intent(in) :: idonictmp
      integer :: igam
      double precision :: arg
      integer :: i,j

      ipelec   = Nint(nssl_params(11))
      isaund   = Nint(nssl_params(12))
      IF ( present(idonictmp) ) THEN
          idonic = idonictmp
      ENDIF

      IF ( ipelec > 0 .or. idonic ) THEN
          print*, 'ipelec, isaund, idonic=',ipelec, isaund, idonic
      ENDIF

      if (ipelec>0) then
          lhab = lh ! 7
          lccn = lhab+1 ! 8
          lnc = lhab+2 ! 9
          lnr = lhab+3 ! 10
          lni = lhab+4 ! 11
          lns = lhab+5 ! 12
          lnh = lhab+6 ! 13
          ltmp = lnh
          lscw  = ltmp+1 ! 14
          lscr  = ltmp+2 ! 15
          lsci  = ltmp+3 ! 16
          lscs  = ltmp+4 ! 17
          lsch  = ltmp+5 ! 18
          ltmp = lsch
          lscpi = ltmp+1 ! 19
          lscni = ltmp+2 ! 20
          ltmp = lscni

          na = ltmp ! 20

          ln(:) = 0
          ln(lc) = lnc
          ln(lr) = lnr
          ln(li) = lni
          ln(ls) = lns
          ln(lh) = lnh

          lsc(:) = 0
          lsc(lc) = lscw
          lsc(lr) = lscr
          lsc(li) = lsci
          lsc(ls) = lscs
          lsc(lh) = lsch

          ! Define qxmin from nssl
          qxmin(lc) = 1.0e-13
          qxmin(li) = 1.0e-13
          qxmin(ls) = 1.0e-13
          qxmin(lr) = 1.0e-12
          qxmin(lh) = 1.0e-12

          IF ( nonigrd .eq. -1 .and. ipelec .gt. 0 ) THEN
              open(unit=90,file='takahashi.txt',form='formatted')
              DO i=nlwc,0,-1
                  read(90,*) (takalu(j,i), j=0,31)
              END DO
              close(90)
          END IF
      ENDIF
      !---------------------------------------------------------------------------
      !End explicit lightning --------------------------------------------------------

      rhoi = 920.
!     rhoi = 500.
      nu = 4.
      nuc = 1000.
      rd = 287.15
      cp = 1005.
      cpw = 4187.
      ao = 1.e-6
      co = 1.e-6
      li0 = 1.e-6
      mi0 = 4./3.*pi*rhoi*(li0)**3
      !gammnu = exp(gammln(nu))
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
      ars = 11.72
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

      iii = 5
      jjj = 6
      kkk = 6
      lll = 8
      mmm = 9
      nnn = 20 
      ooo = 27

 
      OPEN(1,FILE="COLLV5_E0.05.bin",form='unformatted')!!Lookup table for aggregation mass and number
      READ(1) (coll_ni(ii),ii=1,iii) !ni = 1, 10, 100, 1000, 10000 L-1
      READ(1) (coll_an(jj),jj=1,jjj) !an = 1, 10, 100, 1000, 10000 um
      READ(1) (coll_cn(kk),kk=1,kkk) !cn = 1, 10, 100, 1000, 10000 um
      READ(1) (coll_nu(ll),ll=1,lll) !nu = 1, 2, 3, 4, 5, 6, 7, 8
      READ(1) (coll_rho(mm),mm=1,mmm)!rho = 100, 200, 300, 400, 500, 600, 700, 800, 900 kg/m3
      READ(1) (((((coll(ii,jj,kk,ll,mm),ii=1,iii),jj=1,jjj),kk=1,kkk),ll=1,lll),mm=1,mmm)
      !READ(1) (((((ncoll(ii,jj,kk,ll,mm),ii=1,iii),jj=1,jjj),kk=1,kkk),ll=1,lll),mm=1,mmm)
      CLOSE(1)
      
      OPEN(1,FILE="COLLV6_E0.05.bin",form='unformatted')!!Lookup table for aggregation mass and number
      READ(1) (coll_ni(ii),ii=1,iii) !ni = 1, 10, 100, 1000, 10000 L-1
      READ(1) (coll_an(jj),jj=1,jjj) !an = 1, 10, 100, 1000, 10000 um
      READ(1) (coll_cn(kk),kk=1,kkk) !cn = 1, 10, 100, 1000, 10000 um
      READ(1) (coll_nu(ll),ll=1,lll) !nu = 1, 2, 3, 4, 5, 6, 7, 8
      READ(1) (coll_rho(mm),mm=1,mmm)!rho = 100, 200, 300, 400, 500, 600, 700, 800, 900 kg/m3
      !READ(1) (((((coll(ii,jj,kk,ll,mm),ii=1,iii),jj=1,jjj),kk=1,kkk),ll=1,lll),mm=1,mmm)
      READ(1) (((((ncoll(ii,jj,kk,ll,mm),ii=1,iii),jj=1,jjj),kk=1,kkk),ll=1,lll),mm=1,mmm)
      CLOSE(1)
      
      print*,'INIT', coll(1,1,1,1,1),ncoll(1,1,1,1,1)

      OPEN(1,FILE="IPAS_newphi_flat_newform.bin",form='unformatted')!!Lookup table for aggregation a, an, and cn
      READ(1) (coll_phi(nn),nn=1,nnn) !phi = 0.01 --> 100.0 logarithmically spaced
      READ(1) (coll_r(oo),oo=1,ooo)   !r = 1 -> 10, 20 -> 100, 200 -> 1000 microns
      READ(1) ((acoll(nn,oo),nn=1,nnn),oo=1,ooo)  !a_avg
      READ(1) ((ancoll(nn,oo),nn=1,nnn),oo=1,ooo) !an
      READ(1) ((cncoll(nn,oo),nn=1,nnn),oo=1,ooo) !cn
      READ(1) ((ddcoll(nn,oo),nn=1,nnn),oo=1,ooo) !delta density
      CLOSE(1)

      print*,'INIT2', coll_phi, ddcoll(nnn,ooo)

      END SUBROUTINE AHM_INIT

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
      SUBROUTINE MP_AHM(ID,ITIMESTEP,&
      TH, QV, QC, QR, QI, QS, QG, NC, NR, NI, NS,&
      NG, AI, CI, AS, CS, RHO, PII, P, DT, DZ, HT, W,    &
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
      ,QAGG, QNAGG, CPLX, PHIS, RHOS, CPLXS    &
      ,PHI,RHOICE,RELH      &
      ,IDS,IDE, JDS,JDE, KDS,KDE              & ! domain dims
      ,IMS,IME, JMS,JME, KMS,KME              & ! memory dims
      ,ITS,ITE, JTS,JTE, KTS,KTE              & ! tile   dims
!!                ------------------------------ explicite lightning ------------------------------
!!                ---------------------------------------------------------------------------------
              ,SCR                            & !
              ,SCW                            & !
              ,SCI                            & !
              ,SCS                            & !
              ,SCH                            & !
      ,SCTOT                          & !
      ,csaci03d, chaci03d, chacs03d   & !
      ,NONINDUC                       & ! A
      ,rscghis_2d                     & ! D
      ,INDUC                          & ! D
      ,ELEC                           & ! E
      ,SCION                          & ! D
      ,SCIONA                         & !
      ,dx                             & !
      ,dy                             & !
      ,ipelectmp                      & !
      ,umc,umr,umi,ums,umg            & !
!!                ------------------------------ end lightning ------------------------------------
!!                ---------------------------------------------------------------------------------
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

#if ( defined(DM_PARALLEL)  &&   ! defined(STUBMPI) )
#define MPI
              USE module_dm, ONLY : local_communicator, mytask
              ! keep a spacing line here to keep Apple cpp from adding a space in front of the endif
#endif
      IMPLICIT NONE
#if ( defined(DM_PARALLEL)  &&   ! defined(STUBMPI) ) || defined(MPI)
              INCLUDE 'mpif.h'
#else
      integer :: mytask = 0
#endif

      INTEGER, INTENT(IN) :: ids, ide, jds, jde, kds, kde, &
      ims, ime, jms, jme, kms, kme, &
      its, ite, jts, jte, kts, kte, &
      ITIMESTEP,ID

      REAL, DIMENSION(ims:ime, kms:kme, jms:jme), INTENT(INOUT):: &
      QV, QC, QR, QS, QG, QI, NC, NR, NS, NG, NI, AI, CI, AS, CS, TH

      REAL, DIMENSION(ims:ime, kms:kme, jms:jme), INTENT(IN):: &
      PII, P, DZ, RHO, W

      REAL, DIMENSION(ims:ime, kms:kme, jms:jme), INTENT(OUT) ::      &
      ICEDEP, ICESUB, RAINEVAP, SNOWEVAP, SNOWMELT, SNOWDEP, SNOWSUB, &
      SNOWACCR, ICEMELT, ICENUC, RAINFRZ, CLOUDFRZ, PHI, RHOICE, RELH,&
      SEDI, SEDS, SEDR, SEDG, NAGGOUT, NRAGGOUT, NSAGGOUT, PRAOUT,    &
      CLOUDCOND, CLOUDEVAP, CONTACT, PRCOUT, DEP, PRCIOUT, PSACWSOUT, &
      QMULTSOUT, QMULTROUT, PIACROUT, PRACIOUT, PIACRSOUT, PRACISOUT, &
      PRACGOUT, PSACWGOUT, PGSACWOUT, PGRACSOUT, PRDGOUT, EPRDGOUT,   &
      EVPMGOUT, PGMLTOUT, QMULTGOUT,QMULTRGOUT, IMMERSION, PRAIOUT,   &
      CFRZR, QAGG, QNAGG, CPLX, CPLXS, RHOS, PHIS

      REAL, DIMENSION(ims:ime, jms:jme), INTENT(IN) :: HT
      REAL, DIMENSION(ims:ime, jms:jme), INTENT(INOUT) ::   RAINACC, &
      RAINACCT,SNOWACC, SNOWACCT, GRAUPELACC, GRAUPELACCT, FRZNACC, &
      FRZNACCT

      REAL, INTENT(IN):: DT

      REAL, DIMENSION(its:ite, kts:kte, jts:jte) ::T,THL

      REAL, DIMENSION(its:ite, kts:kte) :: QV2D, QC2D, QR2D, QS2D,   &
      QG2D, QI2D, NC2D, NR2D, NS2D, NG2D, NI2D, AI2D, CI2D, AS2D,    &
      CS2D, T2D, RHO2D, P2D, ICEDEP2D, ICESUB2D, RAINEVAP2D,         &
      SNOWEVAP2D, SNOWMELT2D, SNOWDEP2D, SNOWSUB2D, SNOWACCR2D,      &
      CLOUDCOND2D, CLOUDEVAP2D, ICEMELT2D, ICENUC2D, RAINFRZ2D,      &
      CLOUDFRZ2D, PHI2D, RHOICE2D, CPLX2D, PHIS2D, RHOS2D, CPLXS2D,  &
      RELH2D, NAGGOUT2D, NRAGGOUT2D, NSAGGOUT2D, PRAOUT2D, RSED,     &
      GSED, ISED, SSED, CONTACT2D, PRCOUT2D, DEP2D, PRCIOUT2D,       &
      PSACWSOUT2D,QMULTSOUT2D,QMULTROUT2D,PIACROUT2D, PRACIOUT2D,    &
      PIACRSOUT2D, PRACISOUT2D, PRACGOUT2D, PSACWGOUT2D, PGSACWOUT2D,&
      PGRACSOUT2D, PRDGOUT2D, EPRDGOUT2D, EVPMGOUT2D, PGMLTOUT2D,    &
      QMULTGOUT2D, QMULTRGOUT2D, IMMERSION2D, PRAIOUT2D, CFRZR2D,    &
      QAGG2D, QNAGG2D

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
      
      REAL :: start, finish

!                ------------------------------ explicite lightning ------------------------------
!                ---------------------------------------------------------------------------------
      real, dimension(ims:ime, kms:kme, jms:jme),  optional, intent(inout)::  &
              scr,scw,sci,scs,sch,sctot,scion,sciona &
              ,csaci03d, chaci03d, chacs03d &
              ,umc,umr,umi,ums,umg
      real, dimension(ims:ime, kms:kme, jms:jme),  optional, intent(inout)::  &
              noninduc,induc,elec
      real, dimension(ims:ime , jms:jme),  optional, intent(inout)::  &
              rscghis_2d
      real, optional, intent(in) :: dx,dy
      integer, optional, intent(in) :: ipelectmp

      real, dimension(its:ite, 1, kts:kte, na) :: an

      double precision :: chgneg1,chgpos1,chgneg2,chgpos2,chgneg3,chgpos3
      double precision :: chgiona1,chgiona2,chgiona3,sctot3,scwtot,scrtot,scitot,scstot,schtot,schltot
      double precision :: scwmax,scrmax,scimax,scsmax,schmax,schlmax
      double precision :: scwmin,scrmin,scimin,scsmin,schmin,schlmin
      real, dimension(its:ite, kts:kte) :: scsacw,scsaci,schacw,schaci,schacs
      real, dimension(its:ite) :: nonind_sum
      real, dimension(its:ite, 1, kts:kte) :: dz2d
      real, dimension(its:ite,kts:kte) :: elec2d
      real, dimension(its:ite,kts:kte) :: csaci02d, chaci02d, chacs02d
      real umx(its:ite,kts:kte,lc:lhab) ! fall speed

      real :: tmp_elec, dv
      integer :: nx,ny,nz


#ifdef MPI
#if defined(MPI)
              integer, parameter :: ntot = 50
              double precision  mpitotindp(ntot), mpitotoutdp(ntot)
              INTEGER :: mpi_error_code = 1
#endif
#endif
      !---------------------- elec vars description --------------------------------
      ! SCX -  total net space charge on hydrometeor species x (C/kg)
      ! NONINDUC - total non inductive charging rates (C/m3*s)
      ! rscghis_2d - maximum non inductive charging rate in the column (C/m3*s)
      ! INDUC - total inductive/polarization charging rates in C/m3*s (Ziegler et al. 1991,MWR)
      ! ELEC - electric field vector (V/m)
      ! csaci03d = collision rate for ice-snow, nssl eq (#/kg/s)
      ! chaci03d = collision rate for ice-graupel (#/kg/s)
      ! chacs03d = collision rate for snow-graupel (#/kg/s)
      ! csaci03d_morr = collision rate for ice-snow, morr eq  (#/kg/s)
      ! umx = hydrometeor fall speed
      ! qxsten = sedimentation rates
      !---------------------- end elec vars description -----------------------------

      IF ( present( ipelectmp ) ) THEN
          ipelec = ipelectmp
      ELSE
          ipelec = 0
      ENDIF
      IF ( present(scw) .and. ipelec > 0 ) THEN
!          scr = itimestep * 1e-18 ! testing

          chgiona1 = 0.0d0
          chgiona2 = 0.0d0
          chgiona3 = 0.0d0
          chgneg1 = 0.0
          chgpos1 = 0.0
          chgneg2 = 0.0
          chgpos2 = 0.0
          chgneg3 = 0.0
          chgpos3 = 0.0
          sctot3  = 0.0d0
          scwtot  = 0.0d0
          scrtot  = 0.0d0
          scitot  = 0.0d0
          scstot  = 0.0d0
          schtot  = 0.0d0
          schltot  = 0.0d0

          scwmax  = 0.0d0
          scrmax  = 0.0d0
          scimax  = 0.0d0
          scsmax  = 0.0d0
          schmax  = 0.0d0
          schlmax  = 0.0d0

          scwmin  = 0.0d0
          scrmin  = 0.0d0
          scimin  = 0.0d0
          scsmin  = 0.0d0
          schmin  = 0.0d0
          schlmin  = 0.0d0

          scsacw = 0.0
          scsaci = 0.0
          schacw = 0.0
          schaci = 0.0
          schacs = 0.0
          nonind_sum = 0.0

          csaci03d = 0.0
          chaci03d = 0.0
          chacs03d = 0.0
          csaci02d = 0.0
          chaci02d = 0.0
          chacs02d = 0.0

          umc = 0.0
          umr = 0.0
          umi = 0.0
          ums = 0.0
          umg = 0.0
      END IF
!                ------------------------------ end lightning ------------------------------------
!                ---------------------------------------------------------------------------------

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
               AS2D(i,k) = AS(i,k,j)
               CS2D(i,k) = CS(i,k,j)

               T2D(i,k) = T(i,k,j)

!not sure if the sed and process blocks need to be here
               RSED(i,k) = SEDR(i,k,j)
               GSED(i,k) = SEDG(i,k,j)
               SSED(i,k) = SEDS(i,k,j)
               ISED(i,k) = SEDI(i,k,j)

               P2D(i,k) = P(i,k,j)
               RHO2D(i,k) = RHO(i,k,j)
               DZQ(k) = DZ(i,k,j)

            END DO
         END DO

         !     ------------------------------ explicite lightning ------------------------------
         !     ---------------------------------------------------------------------------------

         IF ( present(scw) .and. ipelec > 0 ) THEN
             an(:,1,:,:) = 0.0
             an(its:ite,1,kts:kte,lt)   = th(its:ite,kts:kte,j)
             an(its:ite,1,kts:kte,lscw)   = scw(its:ite,kts:kte,j) * RHO2D(its:ite,kts:kte)
             an(its:ite,1,kts:kte,lscr)   = scr(its:ite,kts:kte,j) * RHO2D(its:ite,kts:kte)
             an(its:ite,1,kts:kte,lsci)   = sci(its:ite,kts:kte,j) * RHO2D(its:ite,kts:kte)
             an(its:ite,1,kts:kte,lscs)   = scs(its:ite,kts:kte,j) * RHO2D(its:ite,kts:kte)
             an(its:ite,1,kts:kte,lsch)   = sch(its:ite,kts:kte,j) * RHO2D(its:ite,kts:kte)
             !        IF ( lhl > 1 ) an(:,1,:,lschl)   = schl(:,:,j)
             an(its:ite,1,kts:kte,lscpi) = Max(0.0, sciona(its:ite,kts:kte,j)) * ec_i * RHO2D(its:ite,kts:kte) ! Convert sciona[C/kg] -> an[#/m3]
             an(its:ite,1,kts:kte,lscni) = -Min(0.0,sciona(its:ite,kts:kte,j)) * ec_i * RHO2D(its:ite,kts:kte)

             DO k = kts,kte
                 DO i = its,ite
                     ! Fake ions for now -- will need to change this if real ions ever get used in WRF
                     ! convert net ion density to positive or negative ions, so that underlying code is
                     ! still ion-friendly
                     tmp_elec =  (an(i,1,k,lscw)+ &
                             an(i,1,k,lscr)+ &
                             an(i,1,k,lsci)+ &
                             an(i,1,k,lscs)+ &
                             an(i,1,k,lsch)+ &
                             ! an(i,1,k,lschl)+ &
                             sciona(i,k,j)) *dx*dy*dz(i,k,j)  ! [kg/m3]*[C/kg]*[m3] = [C]
                     chgneg1 = chgneg1 + Min(0.0, tmp_elec)
                     chgpos1 = chgpos1 + Max(0.0, tmp_elec)
                     chgiona1 = chgiona1 + (an(i,1,k,lscpi)-an(i,1,k,lscni)) * ec * RHO2D(i,k) * dx*dy*dz(i,k,j) ! an(lscions) are in #/m3, ec in C/#. !!! rho_igl looks like a mistake here
                 END DO
             END DO
             if (ipelec>1) then
                 elec2d(its:ite,kts:kte)   =  elec(its:ite,kts:kte,j) ! Ez here =elecsave
                 dz2d (its:ite,1,kts:kte) = dz(its:ite,kts:kte,j)
             end if
         END IF
         !     ------------------------------ end lightning ------------------------------------
         !     ---------------------------------------------------------------------------------
         call CPU_TIME(start)
         !print*,'sulia_harrington'
         CALL AHM_MICRO(T2D, DT, QV2D, QC2D, QR2D, QS2D, QG2D, QI2D,  &
         NC2D, NR2D, NS2D, NG2D, NI2D, AI2D, CI2D, AS2D, CS2D, P2D,   &
         RHO2D, DZQ, IMS, IME, JMS, JME, KMS, KME, ITS, ITE, JTS, JTE,& 
         KTS, KTE, IDS, IDE, JDS, JDE, KDS, KDE, ITIMESTEP, ICEDEP2D, &
         ICESUB2D, RAINEVAP2D, SNOWEVAP2D, SNOWMELT2D, SNOWDEP2D,     &
         SNOWSUB2D, SNOWACCR2D, CLOUDCOND2D, CLOUDEVAP2D, ICEMELT2D,  &
         ICENUC2D, RAINFRZ2D, CLOUDFRZ2D, RSED, ISED, SSED, GSED,     &
         NAGGOUT2D, NRAGGOUT2D, NSAGGOUT2D, PRAOUT2D, PHI2D, RHOICE2D,&
         CPLX2D, PHIS2D, RHOS2D, CPLXS2D, RELH2D, IIN, CCN, IIN_SUM,  &
         IIN_SUMJ, PRECPRT, SNOWRT, SNOWPRT, GRPLPRT, IHR, NHRS,      &
         CONTACT2D,PRCOUT2D, DEP2D, PRCIOUT2D, PSACWSOUT2D,           &
         QMULTSOUT2D,QMULTROUT2D,PIACROUT2D, PRACIOUT2D, PIACRSOUT2D, &
         PRACISOUT2D, PRACGOUT2D, PSACWGOUT2D,PGSACWOUT2D,PGRACSOUT2D,&
         PRDGOUT2D, EPRDGOUT2D, EVPMGOUT2D, PGMLTOUT2D, QMULTGOUT2D,  &
         QMULTRGOUT2D, IMMERSION2D, PRAIOUT2D, CFRZR2D, QAGG2D, QNAGG2D,&
!     ------------------------------   lightning   ------------------------------------
         an, elec2d, ipelec, j, &
         scsacw,scsaci,schacw,schaci,schacs &
         ,CSACI02d,CHACI02d,CHACS02d &
         ,umx &
!     ------------------------------ end lightning ------------------------------------
         )
!     ------------------------------   LIGHTNING --------------------------------------
          IF (IPELEC>0) THEN
             nx = ite-its+1
             ny = 1         ! set up as 2D slabs
             nz = kte-kts+1

             ! ------------------ compute ind/nonind charge rates ------------------
             IF ( present( noninduc ) .and. ipelec >= 1 ) THEN
                 DO k = kts,kte
                     DO i = its,ite
                         induc(i,k,j)  = schacw(i,k) ! schacw
                         noninduc(i,k,j) = schaci(i,k) + schacs(i,k) ! schaci + schacs
                         IF ( .not. ( noninduc(i,k,j) > -1.e-3 .and. noninduc(i,k,j) < 1.e-3 ) ) THEN
                             write(0,*) 'problem with noninduc! i,j,k,noninduc = ',i,j,k,noninduc(i,k,j)
                             write (0,*) 'stopping'
                             print*, 'stopping'
                             STOP
                         ENDIF
                         IF ( .not. ( induc(i,k,j) > -1.e-3 .and. induc(i,k,j) < 1.e-3 ) ) THEN
                             write(0,*) 'problem with induc! i,j,k,induc = ',i,j,k,induc(i,k,j)
                             write (0,*) 'stopping'
                             print*, 'stopping'
                             STOP
                         ENDIF
                         IF ( present( rscghis_2d ) ) THEN
                             nonind_sum(i) = nonind_sum(i) + Abs( noninduc(i,k,j) )*dz(i,k,j)
                         ENDIF
                     ENDDO
                 ENDDO
                 IF ( present( rscghis_2d ) ) THEN
                     DO i = its,ite
                         rscghis_2d(i,j) = Max( rscghis_2d(i,j), nonind_sum(i) )
                     ENDDO
                 ENDIF
                 DO k = kts,kte
                     DO i = its,ite
                         induc(i,k,j) = induc(i,k,j)/rho(i,k,j) ! convert to C/kg/s
                         noninduc(i,k,j) = noninduc(i,k,j)/rho(i,k,j) ! convert to C/kg/s
                     END DO
                 END DO
             ENDIF

             IF ( present( sciona ) .and. lnh > 1) THEN
                 CALL SCIONATT    &
                         &  (nx,ny,nz,na,j &
                         &  ,0,0 &
                         &  ,dt,dz2d &
                         &  ,an &
                         &  )
             ENDIF

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             !!!  Collect into chg...2 (an[m-3])
             DO k = kts,kte
                 DO i = its,ite
                     tmp_elec =   (an(i,1,k,lscw)+ &
                             an(i,1,k,lscr)+ &
                             an(i,1,k,lsci)+ &
                             an(i,1,k,lscs)+ &
                             an(i,1,k,lsch)+ &
                             (an(i,1,k,lscpi)-an(i,1,k,lscni))*ec) *dx*dy*dz(i,k,j)
                     chgneg2 = chgneg2 + Min(0.0, tmp_elec)
                     chgpos2 = chgpos2 + Max(0.0, tmp_elec)
                     chgiona2 = chgiona2 + (an(i,1,k,lscpi)-an(i,1,k,lscni)) * ec * RHO2D(i,k) * dx*dy*dz(i,k,j)
                 ENDDO
             ENDDO
             !!!  End chg...2
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             scw(its:ite,kts:kte,j) = an(its:ite,1,kts:kte,lscw) / RHO2D(its:ite,kts:kte)
             scr(its:ite,kts:kte,j) = an(its:ite,1,kts:kte,lscr) / RHO2D(its:ite,kts:kte)
             sci(its:ite,kts:kte,j) = an(its:ite,1,kts:kte,lsci) / RHO2D(its:ite,kts:kte)
             scs(its:ite,kts:kte,j) = an(its:ite,1,kts:kte,lscs) / RHO2D(its:ite,kts:kte)
             sch(its:ite,kts:kte,j) = an(its:ite,1,kts:kte,lsch) / RHO2D(its:ite,kts:kte)
             sctot(its:ite,kts:kte,j)  = scw(its:ite,kts:kte,j) &
                                        +scr(its:ite,kts:kte,j) &
                                        +sci(its:ite,kts:kte,j) &
                                        +scs(its:ite,kts:kte,j) &
                                        +sch(its:ite,kts:kte,j) &
                                        +sciona(its:ite,kts:kte,j) ! to get C/m3 for RHS calc later on
             csaci03d(its:ite,kts:kte,j) = csaci02d(its:ite,kts:kte)/RHO2D(its:ite,kts:kte)
             chaci03d(its:ite,kts:kte,j) = chaci02d(its:ite,kts:kte)/RHO2D(its:ite,kts:kte)
             chacs03d(its:ite,kts:kte,j) = chacs02d(its:ite,kts:kte)/RHO2D(its:ite,kts:kte)
             umc(its:ite,kts:kte,j) = umx(its:ite,kts:kte,lc)
             umr(its:ite,kts:kte,j) = umx(its:ite,kts:kte,lr)
             umi(its:ite,kts:kte,j) = umx(its:ite,kts:kte,li)
             ums(its:ite,kts:kte,j) = umx(its:ite,kts:kte,ls)
             umg(its:ite,kts:kte,j) = umx(its:ite,kts:kte,lg)

             DO i = its,ite
                 IF ( i.le.(ids+2) .or. i.ge.(ide-2)  ) THEN
                     IF ( j.le.(jds+2) .or. j.ge.(jde-2) ) THEN
                         scw(i,kts:kte,j) = 0.
                         scr(i,kts:kte,j) = 0.
                         sci(i,kts:kte,j) = 0.
                         scs(i,kts:kte,j) = 0.
                         sch(i,kts:kte,j) = 0.
                         sctot(i,kts:kte,j) = 0.
                         csaci03d(i,kts:kte,j) = 0.
                         chaci03d(i,kts:kte,j) = 0.
                         chacs03d(i,kts:kte,j) = 0.
                     ENDIF
                 endif
             enddo

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             !!!  Collect into chg...3 (an[kg-1])
             DO k = kts,kte
                 DO i = its,ite
                     tmp_elec =  rho(i,k,j)*(an(i,1,k,lscw)+ &
                             an(i,1,k,lscr)+ &
                             an(i,1,k,lsci)+ &
                             an(i,1,k,lscs)+ &
                             an(i,1,k,lsch)+ &
                             sciona(i,k,j)) *dx*dy*dz(i,k,j)
                     chgneg3 = chgneg3 + Min(0.0, tmp_elec)
                     chgpos3 = chgpos3 + Max(0.0, tmp_elec)
                     chgiona3 = chgiona3 + sciona(i,k,j)*rho(i,k,j)*dx*dy*dz(i,k,j)
                 END DO
             end do
             !!!  End chg...3
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             ! ------------------ For mpi calculation ------------------------
             dv = rho(i,k,j)*dx*dy*dz(i,k,j)
             IF ( present(elec) ) THEN
                 sctot3 = sctot3 + sctot(i,k,j)*dx*dy*dz(i,k,j)
                 scwtot  = scwtot + scw(i,k,j)*dv
                 scrtot  = scrtot + scr(i,k,j)*dv
                 scitot  = scitot + sci(i,k,j)*dv
                 scstot  = scstot + scs(i,k,j)*dv
                 schtot  = schtot + sch(i,k,j)*dv

                 scwmax = Max( scwmax, scw(i,k,j) )
                 scrmax = Max( scrmax, scr(i,k,j) )
                 scimax = Max( scimax, sci(i,k,j) )
                 scsmax = Max( scsmax, scs(i,k,j) )
                 schmax = Max( schmax, sch(i,k,j) )

                 scwmin = Min( scwmin, scw(i,k,j) )
                 scrmin = Min( scrmin, scr(i,k,j) )
                 scimin = Min( scimin, sci(i,k,j) )
                 scsmin = Min( scsmin, scs(i,k,j) )
                 schmin = Min( schmin, sch(i,k,j) )
             ENDIF

             IF ( present( scw ) ) THEN
#ifdef DM_PARALLEL
                   mpitotindp(1)  = chgpos1
                   mpitotindp(2)  = chgneg1
                   mpitotindp(3)  = chgpos2
                   mpitotindp(4)  = chgneg2
                   mpitotindp(5)  = chgpos3
                   mpitotindp(6)  = chgneg3
                   mpitotindp(7)  = chgiona1
                   mpitotindp(8)  = chgiona2
                   mpitotindp(9)  = chgiona3
                   mpitotindp(10)  = sctot3
                   mpitotindp(11)  = scwtot
                   mpitotindp(12)  = scrtot
                   mpitotindp(13)  = scitot
                   mpitotindp(14)  = scstot
                   mpitotindp(15)  = schtot
                   mpitotindp(16)  = schltot
!                   CALL MPI_Reduce(mpitotindp, mpitotoutdp, 16, MPI_DOUBLE_PRECISION, MPI_SUM, 0, local_communicator, mpi_error_code)
                   IF ( mytask == 0 ) THEN
                       chgpos1 = mpitotoutdp(1)
                       chgneg1 = mpitotoutdp(2)
                       chgpos2 = mpitotoutdp(3)
                       chgneg2 = mpitotoutdp(4)
                       chgpos3 = mpitotoutdp(5)
                       chgneg3 = mpitotoutdp(6)
                       chgiona1 = mpitotoutdp(7)
                       chgiona2 = mpitotoutdp(8)
                       chgiona3 = mpitotoutdp(9)
                       sctot3   = mpitotoutdp(10)
                       scwtot  = mpitotoutdp(11)
                       scrtot  = mpitotoutdp(12)
                       scitot  = mpitotoutdp(13)
                       scstot  = mpitotoutdp(14)
                       schtot  = mpitotoutdp(15)
                       schltot = mpitotoutdp(16)
                   ENDIF
                   mpitotindp(1)  = scwmax
                   mpitotindp(2)  = scrmax
                   mpitotindp(3)  = scimax
                   mpitotindp(4)  = scsmax
                   mpitotindp(5)  = schmax
                   mpitotindp(6)  = schlmax
                   mpitotindp(7)  = -scwmin
                   mpitotindp(8)  = -scrmin
                   mpitotindp(9)  = -scimin
                   mpitotindp(10) = -scsmin
                   mpitotindp(11) = -schmin
                   mpitotindp(12) = -schlmin
!                   CALL MPI_Reduce(mpitotindp, mpitotoutdp, 16, MPI_DOUBLE_PRECISION, MPI_MAX, 0, local_communicator, mpi_error_code)
                   IF ( mytask == 0 ) THEN
                       scwmax  = mpitotoutdp(1)
                       scrmax  = mpitotoutdp(2)
                       scimax  = mpitotoutdp(3)
                       scsmax  = mpitotoutdp(4)
                       schmax  = mpitotoutdp(5)
                       schlmax = mpitotoutdp(6)
                       scwmin  = -mpitotoutdp(7)
                       scrmin  = -mpitotoutdp(8)
                       scimin  = -mpitotoutdp(9)
                       scsmin  = -mpitotoutdp(10)
                       schmin  = -mpitotoutdp(11)
                       schlmin = -mpitotoutdp(12)
                   ENDIF
#endif
              ENDIF ! present( scw )
         ENDIF ! IPELEC>0
!     ------------------------------ end lightning ------------------------------------

!     TRANSFER 2D ARRAYS BACK TO 3D FOR WRF
         CALL CPU_TIME(finish)
         !print '("Time = ",f6.3," seconds.")',finish-start

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
               AS(i,k,j) = AS2D(i,k)
               CS(i,k,j) = CS2D(i,k)
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
               QAGG(i,k,j) = QAGG2D(i,k)
               QNAGG(i,k,j) = QNAGG2D(i,k)
               CPLX(i,k,j) = CPLX2D(i,k)
               CPLXS(i,k,j) = CPLXS2D(i,k)
               RHOS(i,k,j) = RHOS2D(i,k)
               PHIS(i,k,j) = PHIS2D(i,k)

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
               QAGG(i,k,j) = QAGG2D(i,k)
               QNAGG(i,k,j) = QNAGG2D(i,k)

               RHOICE(i,k,j) = RHOICE2D(i,k)
               RELH(i,k,j) = RELH2D(i,k)
               PHI(i,k,j) = PHI2D(i,k)
               CPLX(i,k,j) = CPLX2D(i,k)

               RHOS(i,k,j) = RHOS2D(i,k)
               PHIS(i,k,j) = PHIS2D(i,k)
               CPLXS(i,k,j) = CPLXS2D(i,k)

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

      END SUBROUTINE MP_AHM
      
!*****************************************************************************
      SUBROUTINE AHM_MICRO(t, dt, qv, qc, qr, qs, qg, qi, nc,      &
      nr, ns, ng, ni, ai, ci, as, cs, p, rho, dzq, ims, ime, jms,  &
      jme, kms, kme, its, ite, jts, jte, kts, kte, ids, ide, jds,  &
      jde, kds, kde, itimestep,icedep, icesub, rainevap, snowevap, &
      snowmelt, snowdep, snowsub, snowaccr, cloudcond, cloudevap,  &
      icemelt, icenuc, rainfrz, cloudfrz, rsed, ised, ssed, gsed,  &
      naggout, nraggout, nsaggout, praout, phi, rhoice, cplx, phis,&
      rhos, cplxs, relh, iin, ccn, iin_sum, iin_sumj, precprt,     &
      snowrt, snowprt, grplprt,ihr,nhrs, contact, prcout, dep,     &
      prciout, psacwsout, qmultsout, qmultrout, piacrout, praciout,&
      piacrsout, pracisout, pracgout, psacwgout, pgsacwout,        &
      pgracsout, prdgout, eprdgout, evpmgout, pgmltout, qmultgout, &
      qmultrgout, immersion, praiout, cfrzr, qagg, qnagg, &
      an, elec2d, ipelec, j,&
      scsacw,scsaci,schacw,schaci,schacs &
      ,CSACI0,CHACI0,CHACS0 &
      ,umx &
              )
!****************************************************************************

      IMPLICIT NONE

      INTEGER, INTENT(IN) ::  &
      ims,ime,jms,jme,kms,kme,&
      its,ite,jts,jte,kts,kte,&
      ids,ide,jds,jde,kds,kde,&
      itimestep, nhrs

      REAL, SAVE :: nus = 3.

      REAL, DIMENSION(its:ite,kts:kte), INTENT(INOUT) :: &
      qv, qc, qr, qs, qg, qi, nc, nr, ns, ng, ni, ai, ci, as, cs, t
 
      REAL, DIMENSION(its:ite,kts:kte), INTENT(OUT) :: &
      icedep, icesub, rainevap, snowevap, snowmelt, snowdep, snowsub,&
      snowaccr, cloudcond, cloudevap, icemelt, icenuc, rainfrz,      &
      cloudfrz, dep, contact, ised, ssed, gsed, rsed, rhoice, cplx,  &
      phis, rhos, cplxs, phi, relh

      REAL, DIMENSION(its:ite,kts:kte), INTENT(OUT) :: &
      naggout, nraggout, nsaggout, praout, prciout, prcout, psacwsout,&
      qmultsout, qmultrout, piacrout, praciout, piacrsout, pracisout, &
      pracgout, psacwgout, pgsacwout, pgracsout, prdgout, eprdgout,   &
      evpmgout, pgmltout, qmultgout, qmultrgout, immersion, praiout,  &
      cfrzr, qagg, qnagg

      REAL, DIMENSION(kts:kte), INTENT(IN) :: dzq
      REAL, DIMENSION(its:ite,kts:kte), INTENT(IN) :: p, rho 
      REAL :: dt

      INTEGER i, k, iflag, nstep, n, iaspect, homofreeze,&
      masssizeflag, ipart, resupply, sphrflag, SEDON, RAINON,       &
      DSTRCHECKS, ICE_CALCS, EVOLVE_ON, snowflag, ice_start_time,  &
      processes, LTRUE, CNUM, redden, nucleation, icontactflag,     &
      ccnflag, graupel, demottflag, ihr, meltflag, effflag

!     process rates: calculates changes due to nucleation,deposition,cond/evap
      REAL mnuccd               !change in qi freezing aerosol
      REAL nnuccd               !change in ni freezing aerosol
      REAL prd, ard, crd        !change in qi, ai, and ci, deposition of ice
      REAL qmlt, nmlt           !change in qi and ni due to melting
      REAL amlt, cmlt           !change in ai and ci due to melting
      REAL qsmlt, nsmlt         !change in qs and ns due to melting
      REAL asmlt, csmlt         !change in as and cs due to melting
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
      REAL agg, nagg, aagg, anagg, cnagg, nsagg, nragg, ddagg !change in q and n aggregation
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
      REAL epsg, kt, fh, fv
      REAL ani, cni, rni, anf, cnf, rnf, phii, phif
      REAL deltastr, alphstr, alphv, betam, gammnu
      REAL iwci, iwcf, vi, nidum, voltmp
      REAL weight, igr1, igr2, igr, losstot
      REAL vtbarb, vtbarbm, vtbarblen, vtdumm, vtdumn
      REAL alpham, alpha_v, beta_v, alpha_a, beta_a, lmean, cap
      REAL qt_adv, qt_sed, qi_frac, qiold, niold, dsold
      REAL*8 sui, sup, qvqvs, qvqvsi
      REAL tmpmelt, gi

      !snow characteristics (aggregate microphysics)
      REAL ans, cns, rns, rhobars, deltastrs, vs, swci, swcf, phisf, rnsf, &
           vtbarbs, vtbarbms, vtbarblens, ards, ansf, crds, cnsf, phiis, &
           qs_new, ans_new, cns_new, ns_new, as_new, dd_new, frac

      REAL lammin, lammax, lamc, pgam

      REAL, DIMENSION(its:ite,kts:kte) :: vtrmi1, vtrni1, vtrli1, &
      effi1, vtrmc, vtrms, vtrns, vtrls

!     sedimenation terms
      REAL, DIMENSION(kts:kte) ::                                      &
      fc, fr, fs, fg, fi, fnc, fnr, fns, fng,fni, fci, fai, dumc, dumr,&
      dums,dumg, dumi, dumnc, dumnr, dumns,dumng, dumni, dumai, dumci, &
      falloutc, falloutr, fallouts, falloutg, fallouti, falloutnc,     &
      falloutnr, falloutns, falloutng, falloutni, falloutai, falloutci,&
      qcsten, qrsten, qssten, qisten, ncsten, nrsten, nssten, nisten,  &
      aisten, cisten, qgsten, ngsten, arn, acn, asn, agn,              &
      qiloss2,qcloss,qiloss,qrloss,qcloss2,qrloss2,                     &
      fas,fcs,dumas,dumcs,falloutas,falloutcs,assten,cssten

      REAL umc, unc, ums, umr, umg, uns, unr, ung
      
      REAL, DIMENSION(kts:kte) :: thetal,dthl

      REAL falltndc, falltndr, falltnds,falltndg,  falltndi, falltndnc, &
      falltndnr, falltndns, falltndng, falltndni, falltndai, falltndci,  &
      rgvm, falltndas, falltndcs

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
      
      INTEGER eii1,eii2
      REAL EII_C(8), iieff
      DATA EII_C/0.09,0.21,0.6,0.1,0.08,0.02/ !! ice-ice aggregation efficiency for 
                                                !! -5,-10,-15,-20,-25,-30
                                                !! from Connolly et al. 2012
      
      REAL start, finish

      !     ------------------------------   lightning   ------------------------------------
      ! passed in
      real, dimension(its:ite, 1, kts:kte, na), INTENT(INOUT) :: an
      real, dimension(its:ite,kts:kte), INTENT(INOUT) :: elec2d
!      real, dimension(its:ite, 1, kts:kte), INTENT(INOUT) :: t2_2d,t3_2d,t5_2d,t6_2d,t8_2d
      INTEGER, INTENT(IN)  :: ipelec
      INTEGER, INTENT(IN)  :: j

      ! constants
      logical :: largeion = .false.   ! true: use large ion category
      REAL DIAG,DIAS,DIAI ! # WEIGHTED MEAN DIAMETER
      REAL VT ! RELATIVE FALL SPEED
      REAL :: EHS_COLLSN = 0.5, EHI_COLLSN = 1.0, ESI_COLLSN = 1.0 ! MAX COLLISION EFFICIENCY
      REAL ESICLSN,EHICLSN,EHSCLSN ! COLLISION EFFICIENCY
      REAL BB_G,BB_I,BB_S
      REAL XNU_G,XNU_I,XNU_S
      REAL XMU_G,XMU_I,XMU_S
      REAL DA0LH,DA0LI,DA0LS,DAB0LH_IH,DAB0LH_SH
      REAL RVT,AA2      ! RATIO OF COLLECTION KERNELS (ZRNIC ET AL, 1993)
      REAL XV_I,XV_S ! VOLUMN OF ICE,SNOW. NEEDED FOR CSACI0
      PARAMETER ( RVT = 0.104, AA2 = 5.78E3 )
      REAL, PARAMETER :: CINU = 0.0
      REAL :: DTPINV, DTP
      INTEGER, PARAMETER :: ITEST=1
      REAL CLIONPMX,CLIONNMX
      PARAMETER (CLIONPMX=1.E9,CLIONNMX=1.E9) ! TAKAHASHI 84
      real mixedphasefac
      real qconkq,qconm,qconn,qsign,ftelwc,elecfac
      real cfce15,gf8,gf4i,gf3p5,gf1a,gf1p5,gf4p5
      logical :: csi_flag = .True.
      real ni_elec
      real, parameter :: tfr = 273.15, thnuc = 235.15
      real    :: delqnia = -50.0e-14
      real    :: delqxia =  50.0e-14
      real    :: delqnsa = -50.0e-14
      real    :: delqxsa =  50.0e-14
      real    :: scwppmx = 20.0e-12
      real    :: scippmx = 20.0e-12
      integer, parameter :: icediaflag=4 ! ice diameter used for charging: =1 ani; =2 cni; =3 max axis; =4 rni
      integer, parameter :: snowdiaflag=4
      real   :: ecollmx = 0.5 ! Maximum collision efficiency for graup/hail with ice; used only for charging rates
      real :: scxacymax = 3000.e-12 ! Maximum charging rate magnitude (C s-1 m-3)


      ! tmps..
      REAL TMP, tmppre
      REAL :: tmp1 = 0.
      REAL CCWTMP,CCITMP ! ,CIPTMP,CIRTMP
      DOUBLE PRECISION FRAC_elec
      INTEGER IL, IC
      REAL SCFAC
      REAL CPQC,CPCI,CPQC0,CPCI0
      REAL QCWTMP
      REAL QCITMP
      REAL IL5(its:ite,kts:kte)
      real delqnxa(lc:lqmx)
      real delqxxa(lc:lqmx)
      integer idelq
!      real :: scworig=0.0
      integer iltest
      integer :: iremoveqwfrz = 1    ! Whether to remove (=1) or not (=0) the newly-frozen cloud droplets (ibfc=1) from the CWC used for charge separation
      real :: rarfac  = 1.0
      real lamselec,n0selec,anselec,cnselec,rhoselec,deltastrselec,rnselec
      real lamgelec,n0gelec
      real rhoielec,anielec,cnielec,rhobarelec,deltastrelec,rnielec
      real fspeed

      REAL, INTENT(INOUT), DIMENSION(its:ite,kts:kte) ::   CSACI0 ! SNOW-ICE COLLISION RATE (#COLLI./M3/S)
      REAL, INTENT(INOUT), DIMENSION(its:ite,kts:kte) ::   CHACI0 ! GRAUP-ICE COLLISION RATE (#COLLI./M3/S)
      REAL, INTENT(INOUT), DIMENSION(its:ite,kts:kte) ::   CHACS0 ! GRAUP-SNOW COLLISION RATE (#COLLI./M3/S)
      ! csaci0_morr is simply 9 times NPRAI for now, since eii=0.1
      REAL, DIMENSION(its:ite,kts:kte) ::   csaci0_morr ! SNOW-ICE COLLISION RATE, morr (#COLLI./kg/S)
      REAL, DIMENSION(kts:kte, lscw:lsch) :: dumsc,tmpsc ! for use in sedimentation
      REAL :: SCX(its:ite,kts:kte,LC:LHAB) ! space charge per m3
      REAL, DIMENSION(its:ite,kts:kte) :: CIONP,CIONN,CLIONP,CLIONN,CIONPMXD,CIONNMXD
      REAL temcg(its:ite,kts:kte) ! T in Celcius
      REAL dezcomp(its:ite,kts:kte)
      integer iexy(lc:lqmx,lc:lqmx)
      real exy(ls:lhab,lc:ls) ! sticking efficiency
      real cxacy(its:ite,kts:kte,ls:lhab,lc:ls) ! collision rates
      real xdia(its:ite,kts:kte,lc:lhab) ! diameter; =1/lamx
      real umx(its:ite,kts:kte,lc:lhab) ! fall speed
      real rarx(its:ite,kts:kte,ls:lhab)

      ! scxacy series
      real, intent(out) :: scsacw(its:ite,kts:kte) ! inductive
      real, intent(out) :: scsaci(its:ite,kts:kte) ! noninductive
      !
      real, intent(out) :: schacw(its:ite,kts:kte) ! inductive
      real, intent(out) :: schaci(its:ite,kts:kte), schacs(its:ite,kts:kte) ! noninductive
      !
      real :: scxacy(its:ite,kts:kte,ls:lhab,lc:ls)

      ! charge per kg hdmt
      real :: fscrw
      real :: fsccw
      real :: fscci
      real :: fscsw
      real :: fschw

      ! production of space charge term  + / -
      real ::  psctot
      real ::  psccwi, psccwd ! cloud water
      real ::  psccii, psccid ! ice crystals
      real ::  pscrwi, pscrwd ! rain water
      real ::  pscswi, pscswd ! snow
      real ::  pschwi, pschwd ! graupel

      ! production of space charge from mass transfer
      real ::  psccwmi, psccwmd
      real ::  psccimi, psccimd
      real ::  pscrwmi, pscrwmd
      real ::  pscswmi, pscswmd
      real ::  pschwmi, pschwmd

      ! ions: assume mass = 0
      real ::  pscpii, pscpid !positive ions
      real ::  pscnii, pscnid !negative ions
      real ::  pscplii, pscplid !positive (large) ions
      real ::  pscnlii, pscnlid !negative (large) ions

      ! sedimentation rates
      real ::  pscc_sed(kts:kte)
      real ::  pscr_sed(kts:kte)
      real ::  psci_sed(kts:kte)
      real ::  pscs_sed(kts:kte)
      real ::  pscg_sed(kts:kte)
      real, dimension(kts:kte) :: faloutscc,faloutscr,faloutsci,faloutscs,faloutscg
      real :: faltndscc,faltndscr,faltndsci,faltndscs,faltndscg

      ! temp mp rates for elec
      real :: praci_e,pracis_e
      real :: pra_e,prc_e,psacws_e,qmults_e,qmultg_e,psacwg_e,pgsacw_e
      real :: psmlt_e,psacr_e
      real :: pgmlt_e
      real :: pracs1_e,qmultr_e,qmultrg_e,piacr_e,piacrs_e,pracg1_e,pgracs_e,pracs_e,pracg_e

      IF (IPELEC>0) THEN
          ! INITIALIZE VALUES NEEDED FOR EXPLICIT LIGHTNING

          DTP = DT
          DTPINV = 1.D0/DTP

          BB_G = 1./3.
          BB_I = 0.3429
          BB_S = 1./3.

          XNU_G = (0.-2.)/3.
          XNU_I = 0.
          XNU_S = -0.8

          XMU_G = 1./3.
          XMU_I = 1.
          XMU_S = 1.

          DA0LH = DELBK(BB_G, XNU_G, XMU_G, 0)
          DA0LI = DELBK(BB_I, XNU_I, XMU_I, 0)
          DA0LS = DELBK(BB_S, XNU_S, XMU_S, 0)
          DAB0LH_IH = DELABK(BB_I, BB_G, XNU_I, XNU_G, XMU_I, XMU_G, 0)
          DAB0LH_SH = DELABK(BB_S, BB_G, XNU_S, XNU_G, XMU_S, XMU_G, 0)
          gf4p5 = 11.63172839656745 ! gamma(4.0+0.5)
          gf3p5 = 3.32335097 ! gamma(3.5)
          gf1p5 = 0.8862269255  ! gamma(1.5)
!          print*, 'DA0LH,DA0LI,DA0LS',DA0LH,DA0LI,DA0LS
!          print*, 'DAB0LH_IH,DAB0LH_SH',DAB0LH_IH,DAB0LH_SH

          DO i = its,ite
              DO k = kts,kte
                  dezcomp(i,k) = Max(-1.00e5, Min(elec2d(i,k),1.00e5) ) ! Ez is passed (elecsave)
              END DO
          END DO

          mixedphasefac = 1.0

          iexy(:,:)=0 ! sets to zero the ones I might have forgotten
          !     snow
          iexy(ls,li) = 1
          iexy(ls,lc) = 1 ; iexy(ls,lr) = 0 ;
          !     graupel
          iexy(lh,ls)  = 1 ; iexy(lh,li) = 1 ;
          iexy(lh,lc) = 1 ; iexy(lh,lr)  = 0 ;

          delqnxa(lc:lqmx) = 0.0 ! delqnia
          delqxxa(lc:lqmx) = 0.0 ! delqxia
          delqnxa(li) = delqnia
          delqxxa(li) = delqxia
          delqnxa(ls) = delqnsa
          delqxxa(ls) = delqxsa

          elecfac = 1.0

          ! Initialization

          !        LAMC(kts:kte) = 0.0 ! safety precaution for use in induc charging
          scxacy = 0.0
          cxacy  = 0.0
          xdia   = 0.0
          scx    = 0.0
          exy    = 0.0
          umx    = 0.0

          schacw = 0.0
          scsacw = 0.0
          schaci = 0.0
          scsaci = 0.0
          schacs = 0.0

          fscrw = 0.0
          fsccw = 0.0
          fscci = 0.0
          fscsw = 0.0
          fschw = 0.0

          psccwi = 0.0 ! cloud water
          psccwd = 0.0
          psccwmi = 0.0
          psccwmd = 0.0
          psccii = 0.0 ! cloud ice
          psccid = 0.0
          psccimi = 0.0
          psccimd = 0.0
          pscrwi = 0.0 ! rain water
          pscrwd = 0.0
          pscrwmi = 0.0
          pscrwmd = 0.0
          pscswi = 0.0 ! snow water
          pscswd = 0.0
          pscswmi = 0.0
          pscswmd = 0.0
          pschwi = 0.0 ! graupel
          pschwd = 0.0
          pschwmi = 0.0
          pschwmd = 0.0

          pscpii = 0.0
          pscpid = 0.0
          pscnii = 0.0
          pscnid = 0.0
          pscplii = 0.0
          pscplid = 0.0
          pscnlii = 0.0
          pscnlid = 0.0

          pscc_sed = 0.0
          pscr_sed = 0.0
          psci_sed = 0.0
          pscs_sed = 0.0
          pscg_sed = 0.0

          CSACI0 = 0.0
          CHACI0 = 0.0
          CHACS0 = 0.0
          csaci0_morr = 0.0

          faloutscc = 0.0
          faloutscr = 0.0
          faloutsci = 0.0
          faloutscs = 0.0
          faloutscg = 0.0

          faltndscc = 0.0
          faltndscr = 0.0
          faltndsci = 0.0
          faltndscs = 0.0
          faltndscg = 0.0

          CIONP = 0.0
          CIONN = 0.0
          CLIONP = 0.0
          CLIONN = 0.0
          CIONPMXD = 0.0
          CIONNMXD = 0.0

          DUMSC = 0.0

          vtrms = 0.0
          vtrmi1 = 0.0
      END IF !(IPELEC>0)
      !     ------------------------------ end lightning ------------------------------------
      resupply=0
      masssizeflag=0
      ipart=1                   !0 for RAMS dendrites and needles (Walko et al., 1995)
!     1 for RAMS plates and columns (Walko et al., 1995)
!     2 for Mitchell's mass relations (hex plates, columns)
!     3 for Mitchell's mass relations (stellars, columns)
!     4 for Wood's (2007) mass-size relations
      iaspect    = 0            !set constant aspect ratio (sensitivity study)
      sphrflag   = 0            !all ice assumed spheres
      redden     = 0            !reduces density of spheres
      homofreeze = 1            !homogeneous freezing
      snowflag  = 1             !ice --> snow & aggregation
      meltflag = 1              !1 Jensen_Ishmael melting scheme, o/w instantaneous
      effflag = 1               !1 for ice-ice agg efficiency from connolly et al 2012
      SEDON      = 1            !sedimentation
      EVOLVE_ON  = 1            !depositional growth
      RAINON     = 2            !0==all rain off; 1==basic processes on; 2==all rain processes
      ICE_CALCS  = 1            !all ice calculations
      ice_start_time = 60.*60.*0.0 !time to begin ice nucleation & homogeneous freezing !LG change to 4
      LTRUE = 0
      graupel = 1               !1 on
      processes = 1
      nucleation = 0 !0 simple, 1 meyers, 2 demott
      demottflag = 1 !0 DeMott 2010, 1 DeMott 2015
      icontactflag = 0 !0 contact freezing off, 1 on
      ccnflag = 0 !0 CCN from aerosol data used, 1 CCN from aerosol data NOT used

      IF(sphrflag  .eq. 1) rhoi = 920.
      If(redden .eq. 1) rhoi = 500.
  
      gammnu = exp(gammln(nu))
      alphap = 0.0
      alphamsp = (4./3.)*pi*rhoi
      betamp = 3.0
      
      qt_adv=0.0
      qt_sed=0.0

      DO i = its,ite
         DO k = kts,kte

            !ICE
            IF(qi(i,k).gt.qsmall.and.ni(i,k).gt.qsmall)THEN
               ni(i,k) = max(ni(i,k),qsmall)
               ai(i,k) = max(ai(i,k),qsmall)
               ci(i,k) = max(ci(i,k),qsmall)
               if(isnan(ai(i,k)))print*,'ai1',real(itimestep)*dt,ai(i,k),ani,nu,ni(i,k)
               if(ci(i,k).gt.qsmall)ani = ((ai(i,k)**2)/(ci(i,k)*nu*ni(i,k)))**(1./3.)
               if(ai(i,k).gt.qsmall)cni = ((ci(i,k)**2)/(ai(i,k)*nu*ni(i,k)))**(1./3.)
               if(isnan(ai(i,k)))print*,'ai2',real(itimestep)*dt,ai(i,k),ani,nu,ni(i,k)
            ELSE 
               qi(i,k) = 0.0
               ni(i,k) = 0.0
               ai(i,k) = 0.0
               ci(i,k) = 0.0
               ani = 0.0
               cni = 0.0
               if (ipelec>0) then
               an(i,1,k,lsc(li)) = 0.0
               end if
            END IF
            rni = (ani**2*cni)**(1./3.)
            ai(i,k) = nu*ni(i,k)*ani
            ci(i,k) = nu*ni(i,k)*cni

            !CLOUD WATER
            IF(qc(i,k).lt.qsmall)THEN
               qc(i,k) = 0.
               nc(i,k) = 0.
               if (ipelec>0) then
               an(i,1,k,lsc(lc)) = 0.0
               end if
            END IF

            !RAIN WATER
            IF(qr(i,k).lt.qsmall)THEN
               qr(i,k) = 0.
               nr(i,k) = 0.
               if (ipelec>0) then
               an(i,1,k,lsc(lr)) = 0.0
               end if
            END IF

            !SNOW (AGGREGATES)
            IF(qs(i,k).gt.qsmall.and.ns(i,k).gt.qsmall)THEN  
               IF(snowflag .eq. 2)THEN               
                  ns(i,k) = max(ns(i,k),qsmall)
                  
                  if(cs(i,k).gt.qsmall)ans = ((as(i,k)**2)/(cs(i,k)*nus*ns(i,k)))**(1./3.)
                  if(as(i,k).gt.qsmall)cns = ((cs(i,k)**2)/(as(i,k)*nus*ns(i,k)))**(1./3.)
                  
                  as(i,k) = nus*ns(i,k)*ans
                  cs(i,k) = nus*ns(i,k)*cns
                  rns = (ans**2*cns)**(1./3.)
                  
                  !if(as(i,k).gt.0.0)print*,'as0',real(itimestep)*dt,i,k,as(i,k),cs(i,k),((cs(i,k)**2)/(as(i,k)*nus*ns(i,k)))**(1./3.)
               END IF
            ELSE
               rhos(i,k) = 920.
               qs(i,k) = 0.
               ns(i,k) = 0.
               as(i,k) = 0.
               cs(i,k) = 0.
               ans = 0.0
               cns = 0.0
               rns = 0.0
               if (ipelec>0) then
               an(i,1,k,lsc(ls)) = 0.0
               end if
            END IF
            

            !GRAUPEL
            IF(qg(i,k).lt.qsmall)THEN
               qg(i,k) = 0.
               ng(i,k) = 0.
               if (ipelec>0) then
               an(i,1,k,lsc(lg)) = 0.0
               end if
            END IF
!     initialize process rates
            mnuccd = 0.;mnuccr = 0.;mnuccc = 0.;mnucci = 0.
            nnuccd = 0.;nnuccr = 0.;nnuccc = 0.;nnucci = 0.

            prd = 0.0; ard = 0.0; crd = 0.0; pcc = 0.0; prc = 0.0
            pra = 0.0; prds = 0.0; eprds = 0.0; nsubr = 0.0; nsubs = 0.0
            pre = 0.0; npra = 0.0; nprci = 0.; nprai = 0.; nprc = 0.0
            nprc1 = 0.0; prci = 0.0; agg = 0.0; nagg = 0.0; aagg = 0.0
            anagg = 0.0; cnagg = 0.0; ddagg = 0.0; nsagg = 0.0
            nragg = 0.0; psmlt = 0.0; pracs = 0.0; evpms = 0.0
            nsmlts = 0.0; nsmltr = 0.0; 
            
            qmlt = 0.0; amlt = 0.0; cmlt = 0.0; nmlt = 0.0
            qsmlt = 0.0; asmlt = 0.0; csmlt = 0.0; nsmlt = 0.0

            pracg = 0.0; npracg = 0.0; evpmg = 0.0; pgmlt = 0.0
            ngmltg = 0.0; ngmltr = 0.0; psacws = 0.0; npsacws = 0.0
            psacwg = 0.0; npsacwg = 0.0; pracs1 = 0.0; npracs1 = 0.0
            psacr = 0.0; pracg1 = 0.0; npracg1 = 0.0; qmultr = 0.0
            nmultr = 0.0; qmults = 0.0; nmults = 0.0; qmultg = 0.0
            nmultg = 0.0; qmultrg = 0.0; nmultrg = 0.0; pgsacw = 0.0
            nscng = 0.0; pgracs = 0.0; ngracs = 0.0; niacr = 0.0
            piacr = 0.0; praci = 0.0; niacrs = 0.0; piacrs = 0.0
            pracis = 0.0; prdg = 0.0; eprdg = 0.0; nsubg = 0.0
            qitend = 0.0; nitend = 0.0; prai = 0.0

            vtbarbm = 0.;vtbarb = 0.;vtbarblen = 0.
            vtbarbms = 0.;vtbarbs = 0.;vtbarblens = 0.

            qrsten(k) = 0.0; qcsten(k)  = 0.0; qssten(k) = 0.0
            qgsten(k) = 0.0; qisten(k) = 0.0; ncsten(k) = 0.0
            nrsten(k) = 0.0; nssten(k) = 0.0; ngsten(k) = 0.0
            nisten(k) = 0.0; aisten(k) = 0.0; cisten(k) = 0.0
            assten(k) = 0.0; cssten(k) = 0.0; qiloss(k) = 0.0
            qrloss(k) = 0.0; qcloss(k) = 0.0; qiloss2(k) = 0.0
            qcloss2(k) = 0.0; qrloss2(k) = 0.0; phi(i,k) = 1.0
            phis(i,k) = 1.0; cplx(i,k) = 0.0; cplxs(i,k) = 0.0


            icedep(i,k) = 0.0; icesub(i,k) = 0.0;
            cloudcond(i,k) = 0.0; cloudevap(i,k) = 0.0; icemelt(i,k) = 0.0;
            rainfrz(i,k) = 0.0; cloudfrz(i,k) = 0.0;

            prciout(i,k) = 0.0; qagg(i,k) = 0.0; qnagg(i,k) = 0.0

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
            deltastrs = 1.
            rhobar = 920.
            rhobars = 920.
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
            kt = 2.3823e-2 +7.1177e-5*(temp-273.15)
            fh = 1.0
            fv = 1.0
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
            asn(k) = (rhosu/rho(i,k))**0.54*ars
!     "a" parameter for graupel fall speed
            agn(k) = (rhosu/rho(i,k))**0.54*ag
!     set droplet concentration (units of kg-1)

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

            IF(sphrflag .eq. 1) igr=1.0

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
            IF(snowflag .eq. 1)THEN
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
            END IF              !snowflag
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
                        if (ipelec>0) then
                        an(i,1,k,lsc(lr)) = an(i,1,k,lsc(lr)) + an(i,1,k,lsc(ls))
                        an(i,1,k,lsc(ls)) = 0.0
                        end if
                     ENDIF
                     
                     IF(qg(i,k).lt.1.e-6)THEN
                        qr(i,k) = qr(i,k)+qg(i,k)
                        nr(i,k) = nr(i,k)+ng(i,k)
                        temp=temp-qg(i,k)*xxlf/cpm
                        qg(i,k) = 0.0
                        ng(i,k) = 0.0
                        if (ipelec>0) then
                        an(i,1,k,lsc(lr)) = an(i,1,k,lsc(lr)) + an(i,1,k,lsc(lg))
                        an(i,1,k,lsc(lg)) = 0.0
                        end if
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
                  END IF        !RAINON>=1
!     collection of snow by rain above freezing
!     formula from Ikawa and Saito (1991)
                  
                  IF(RAINON.eq.2)THEN
                     IF(snowflag.eq.1)THEN
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
                        END IF  !qr & qs >=1.e-8
                     END IF     !snowflag
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
                           
!     pracg is mixing ratio of rain per sec collected by graupel
                           pracg = cons41*(((1.2*umr-0.95*umg)**2+0.08*umg*&
                           umr)**0.5*rho(i,k)*n0rr*n0g/lamr**3*(5./&
                           (lamr**3*lamg)+2./(lamr**2*lamg**2)+0.5/&
                           (lamr*lamg**3)))
!     assume 1 mm drops are shed, get number shed per sec
                           npracg = cons32*rho(i,k)*(1.7*(unr-ung)**2+0.3*&
                           unr*ung)**0.5*n0rr*n0g*(1./(lamr**3*lamg)+1.0&
                           /(lamr**2*lamg**2)+1./(lamr*lamg**3))
                           
                           dum = pracg/5.2E-7
                           npracg = npracg-dum
                           
                        END IF
                     END IF
!     melting of snow
!     snow may persist above freezing, from Rutledge and Hobbs (1984)
!     if supersat, snow melts to form rain
                     IF(snowflag .eq. 1)THEN
                        IF(qs(i,k).ge.1.e-8)THEN
                           dum = -cpw/xxlf*(temp-273.15)*pracs
                           psmlt = qsdum*kap*(273.15-temp)/xxlf+dum
                           IF(qvqvs.lt.1.)THEN
                              epss = qsdum*dv
                              evpms = (qv(i,k)-qvs)*epss/ab
                              evpms = max(evpms,psmlt)
                              psmlt = psmlt-evpms
                           END IF
                        END IF  !qs>=1.e-8
                        pracs = 0.
                     END IF     !snowflag
                     
!     melting of graupel
!     graupel may persist above freezing, from Rutledge and Hobbs
!     (1984)
!     if supersat, graupel melts to form rain
                     if(graupel.eq.1)then
                        if(qg(i,k).ge.1.e-8)then
                           dum = -cpw/xxlf*(temp-273.15)*pracg
                           pgmlt = 2.*pi*n0g*kap*(273.15-temp)/xxlf*(f1s/&
                           (lamg*lamg)+f2s*(agn(k)*rho(i,k)/mu)**0.5*&
                           sc**(1./3.)*cons11/(lamg**cons36))+dum !melting of graupel to form rain
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
                  END IF        !RAINON==2
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
                  IF(snowflag.eq.1)THEN
                     dum = (pracs-evpms-psmlt)*dt !melting, evap, & accretion of snow
                     IF(dum.gt.qs(i,k).and.qs(i,k).ge.qsmall)THEN
                        ratio = qs(i,k)/dum
                        pracs = pracs*ratio
                        psmlt = psmlt*ratio
                        evpms = evpms*ratio
                     END IF
                  END IF        !snowflag
                  
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
                  
                  
               ELSE             !temp separation
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
                  END IF        !RAINON>=1
!     aggregation of qs
                  IF(snowflag.eq.1)THEN
                     IF(qs(i,k).ge.1.e-8)THEN
                        nsagg = cons15*asn(k)*(qs(i,k)*rho(i,k))**&
                        ((2.+BS)/3.)*(ns(i,k)*rho(i,k))**((4.-BS)/3.)&
                        /rho(i,k)
                     END IF
                  END IF        !snowflag
                  
!     accretion of cloud droples onto snow/graupel
!     use continuous growth equations with
!     simple gravitational collection kernel ignoring snow
                  if(graupel .eq. 1)then
                     IF(snowflag.eq.1)THEN
                        if(qs(i,k) .ge.1.e-8.and.qc(i,k).ge.qsmall)then
                           psacws = cons13*asn(k)*qc(i,k)*rho(i,k)*n0s/&
                           lams**(BS+3.)
                           
                           npsacws = cons13*asn(k)*nc(i,k)*rho(i,k)*n0s/&
                           lams**(BS+3.)
                        end if
                     END IF     !snowflag
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
                     IF(snowflag.eq.1)THEN
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
                     END IF     !snowflag
                     
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
                           unr*ung)**0.5*n0rr*n0g*(1./(lamr**3*lamg)+1.0/&
                           (lamr**2*lamg**2)+1./(lamr*lamg**3))
                           
                           pracg1 = min(pracg1,qr(i,k)/dt)
                           
                        END IF
                     END IF     !graupel
                  END IF        !RAINON==2
                  
!     rime-splintering of snow
!     hallet-mossop 1974
!     number of splinters formed is based on mass of rimed water
                  IF (graupel .eq. 1)THEN
                     IF(snowflag.eq.1)THEN
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
                     END IF     !snowflag
                     
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
                     IF(snowflag.eq.1)THEN
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
                     END IF     !snowflag
                  END IF        !graupel
                  
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
                  END IF        !RAINON>=1
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
                     IF(snowflag.eq.1)THEN
                        IF(qs(i,k).ge.qsmall)THEN
                           epss = 2.*pi*n0s*rho(i,k)*dv*(f1s/(lams*lams)+&
                           f2s*(asn(k)*rho(i,k)/mu)**0.5*sc**(1./3.)*cons10/&
                           (lams**cons35))
                        ELSE
                           epss = 0.0
                        END IF      
                        prds = epss*(qv(i,k)-qvi)/abi
                     END IF     !snowflag
                  END IF        !RAINON==2
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
                  END IF        !graupel
                  
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
                  IF(snowflag.eq.1)THEN
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
                  END IF        !snowflag
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
                  
               END IF           !temp<273.15
            END IF              !processes

!     IF(ICE_CALCS .eq. 1)THEN
            IF(qi(i,k).gt.qsmall.and.ni(i,k).gt.qsmall)THEN
               
!     set minimum values for ni,ai,ci, otherwise
               
               ni(i,k) = max(ni(i,k),qsmall)
               ai(i,k) = max(ai(i,k),qsmall)
               ci(i,k) = max(ci(i,k),qsmall)
               
!     get characteristic ci,ai (i.e. cni,ani) from ci,ai
               
               ani = ai(i,k)/(nu*ni(i,k))
               cni = ci(i,k)/(nu*ni(i,k))
               

               CALL ICE_CHECKS(1, nu,ni(i,k),qi(i,k),ani,cni,rni,deltastr,rhobar,&
               iaspect,sphrflag,redden)
               
               
               ci(i,k)=nu*ni(i,k)*cni
               ai(i,k)=nu*ni(i,k)*ani
!     get iwc to calculate iwc tendency
!     by substracting final and initial values
               vi = 4./3.*pi*rni**3.*exp(gammln(nu+deltastr+2.))/gammnu 
               iwci = ni(i,k)*rhobar*vi*rho(i,k)

!     calculate number concentration from number mixing ratio

               IF(EVOLVE_ON .eq. 1) THEN
                  
                  nidum = ni(i,k)*rho(i,k)
                  CALL EVOLVE(nu,ani,nidum,sui,sup,qvv,temp,press,igr,dt,iwcf,&
                  cnf,iwci,phii,phif,cni,rni,rnf,anf,deltastr,mu,&
                  rhobar,vtbarb,vtbarbm,vtbarblen,rhoa,i,k,iaspect&
                  ,ipart,sphrflag,redden,itimestep) 
                  

                  betam = deltastr + 2.0
                  alphstr = co/ao**(deltastr)
                  alphv = 4./3.*pi*alphstr
                  
                  phi(i,k)=phif
                  
!     get deposition/sublimation process rates for qi, ai, and ci
!     deposition rate for ice [=] kg/kg/s (mixing ratio rate)
                  
                  prd=(iwcf-iwci)/rho(i,k)/dt              
                  ard=(anf-ani)*nu*ni(i,k)/dt
                  crd=(cnf-cni)*nu*ni(i,k)/dt 
               END IF           !EVOLVE_ON
            
            !melting of ice follows that in jensen_ishmael
                IF(meltflag.eq.1)THEN
                    IF(temp.ge.273.15)THEN
                !.. Melting: ice shapes become more spherical

                       qmlt = 2.*pi*(kt*fh*(273.15-temp)+rho(i,k)*xxlv*dv*fv*(qvs-qv(i,k)))/   & 
                            xxlf*(ni(i,k)*nu*max(ani,cni)) 
                            !-(cpw/xxlf*(temp-273.15)*(rimetotal(cc)/rhoair(k) + dQImltri(cc)))

                       qmlt = min(qmlt,0.)
                       qmlt = max(qmlt,(-qi(i,k)/dt))

!.. Don't let small very number mixing ratios (large sizes) cause ice precipitation
                       if(qmlt.lt.0.) then
                          if(ai(i,k).lt.1.e-12.or.ci(i,k).lt.1.e-12) then
                             qmlt = -qi(i,k)/dt
                          endif
                       endif

                       nmlt = max((-ni(i,k)/dt),(ni(i,k)*qmlt/qi(i,k)))

                       gi=MIN((MAX((NINT(((nu+2.+deltastr)*100000.)-355000+1)),1)),505001)
                       tmpmelt= 4./3.*pi*alphstr*exp(gammln(gi))/gammnu

                       amlt = (1./tmpmelt)*(1./rhobar)*qmlt + ai(i,k)*nmlt/ni(i,k)
                       cmlt = amlt*cni/ani*(1.+(2.*deltastr))/(2.+deltastr)

                       icemelt(i,k) = qmlt*-1.
                    END IF !temp
                END IF !melt
               
            END IF              !ice is present
!     END IF                 !ICE_CALCS


            IF(snowflag .eq.2)THEN
               IF(qs(i,k).gt.qsmall.and.ns(i,k).gt.qsmall)THEN
                  if(isnan(ns(i,k)).or..not.ieee_is_finite(ns(i,k)))&
                  print*,'ns2',real(itimestep)*dt,i,k,ns(i,k)
                  ns(i,k) = max(ns(i,k),qsmall)
                  
                  ans = as(i,k)/(nus*ns(i,k))
                  cns = cs(i,k)/(nus*ns(i,k))

                  if(isnan(ans).or..not.ieee_is_finite(ans))&
                  print*,'ans',real(itimestep)*dt,i,k,ans,as(i,k),nus,ns(i,k)
                  

                  !print*,'snow1',itimestep, rhos(i,k),qs(i,k),ns(i,k),ans,cns,nus,as(i,k),cs(i,k)
                  !check that deltastr, rhobar, and rni are within reasonable bounds
                  IF(isnan(rns))print*,"BEFORECHECKS",rns
                  CALL ICE_CHECKS(2, nus,ns(i,k),qs(i,k),ans,cns,rns,deltastrs,rhos(i,k),&
                       iaspect,sphrflag,redden)
                  IF(isnan(rns))print*,"AFTERCHECKS",rns
                    
                  !CALL DSTR_CHECK(2,nus,ns(i,k),ans,cns,deltastrs,iaspect,sphrflag)
                      !CALL RHO_CHECK(flag,n,deltastr,qi,ni,ani,cni,sphrflag,redden,rhobar)
                  !CALL R_CHECK(nus,qs(i,k),ns(i,k),cns,ans,rns,rhos(i,k),deltastrs)
                  !print*,'snow2',itimestep, rhos(i,k),deltastrs
                  cs(i,k)=nus*ns(i,k)*cns
                  as(i,k)=nus*ns(i,k)*ans
                  if(isnan(ns(i,k)).or..not.ieee_is_finite(ns(i,k)))&
                  print*,'ns3',real(itimestep)*dt,i,k,ns(i,k)
                  !     get iwc to calculate iwc tendency
                  !     by substracting final and initial values
                  vs = 4./3.*pi*rns**3.*exp(gammln(nus+deltastrs+2.))/exp(gammln(nus)) 
                  swci = ns(i,k)*rhos(i,k)*vs*rho(i,k)
                  if(isnan(swci))print*,'swci',real(itimestep)*dt,i,k,swci,&
                  ns(i,k),rhos(i,k),vs,rho(i,k),deltastrs,rns,nus
                  IF(iaspect .eq. 1) igr = .27
                  IF(sphrflag .eq. 1) igr=1.0
                  If(redden .eq. 1) rhobars = 100.
                  !     calculate number concentration from number mixing ratio
                  
                  IF(EVOLVE_ON .eq. 1) THEN
                     
                     nidum = ns(i,k)*rho(i,k)
                     
                     IF(isnan(rns))print*,"BEFOREEVOLVE",rns
                     CALL EVOLVE(nus,ans,nidum,sui,sup,qvv,temp,press,igr,dt,swcf,&
                          cnsf,swci,phiis,phisf,cns,rns,rnsf,ansf,deltastrs,mu,&
                          rhos(i,k),vtbarbs,vtbarbms,vtbarblens,rhoa,i,k,iaspect&
                          ,ipart,sphrflag,redden,itimestep) 
                     IF(isnan(rns))print*,"AFTEREVOLVE",rns
                           
                     phis(i,k)=phisf
                     
                     !     get deposition/sublimation process rates for qi, ai, and ci
                     !     deposition rate for ice [=] kg/kg/s (mixing ratio rate)
                     
                     prds=(swcf-swci)/rho(i,k)/dt              
                     ards=(ansf-ans)*nus*ns(i,k)/dt
                     crds=(cnsf-cns)*nus*ns(i,k)/dt
                     !if(isnan(prds).or.prds.lt.0.0)print*,'prds0',real(itimestep)*dt,i,k,&
                     !prds,swcf,swci,rho(i,k),dt
                  END IF           !EVOLVE_ON
                  
                  !melting of snow follows that of melting ice jensen_ishmael
                  IF(meltflag.eq.1)THEN
                    IF(temp.ge.273.15)THEN
                !.. Melting: ice shapes become more spherical

                       qsmlt = 2.*pi*(kt*fh*(273.15-temp)+rho(i,k)*xxlv*dv*fv*(qvs-qv(i,k)))/   & 
                            xxlf*(ns(i,k)*nu*max(ans,cns)) 
                            !-(cpw/xxlf*(temp-273.15)*(rimetotal(cc)/rhoair(k) + dQImltri(cc)))

                       qsmlt = min(qsmlt,0.)
                       qsmlt = max(qsmlt,(-qs(i,k)/dt))

!.. Don't let small very number mixing ratios (large sizes) cause ice precipitation
                       if(qsmlt.lt.0.) then
                          if(as(i,k).lt.1.e-12.or.cs(i,k).lt.1.e-12) then
                             qsmlt = -qs(i,k)/dt
                          endif
                       endif

                       nsmlt = max((-ns(i,k)/dt),(ns(i,k)*qsmlt/qs(i,k)))

                       alphstr = co/ao**deltastrs
                       gi=MIN((MAX((NINT(((nus+2.+deltastrs)*100000.)-355000+1)),1)),505001)
                       tmpmelt= 4./3.*pi*alphstr*exp(gammln(gi))/exp(gammln(nus))

                       asmlt = (1./tmpmelt)*(1./rhos(i,k))*qsmlt + as(i,k)*nsmlt/ns(i,k)
                       csmlt = asmlt*cns/ans*(1.+(2.*deltastrs))/(2.+deltastrs)

                    END IF!temp
                END IF!melt
                  
    
               END IF              !snow is present
            END IF!snowflag




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
                   CALL SIMPLENUC(temp,sui,rho(i,k),ni(i,k),&
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
                  alphstr=co/ao**(deltastr)
                  alphv=4./3.*pi*alphstr
                  betam=2.+deltastr   
                  
                  ani=((qi(i,k)*gammnu)/(rhobar*ni(i,k)*alphv*&
                  exp(gammln(nu+betam))))**(1./betam)

                  cni=co*(ani/ao)**deltastr
                  ci(i,k)=nu*ni(i,k)*cni
                  ai(i,k)=nu*ni(i,k)*ani
               END IF           ! iflag = 1

               CALL R_CHECK(1,nu,qi(i,k),ni(i,k),cni,ani,rni,rhobar,deltastr)
               
               ci(i,k)=nu*ni(i,k)*cni
               ai(i,k)=nu*ni(i,k)*ani
              
            END IF              ! q > qsmall
                

            if(snowflag.eq.2)then
               !prds_test = prds
               prds=max(prds,-qs(i,k)/dt)
               
               qs(i,k)=qs(i,k)+(prds)*dt
               as(i,k)=as(i,k)+(ards)*dt
               cs(i,k)=cs(i,k)+(crds)*dt
               !if(as(i,k).gt.0.0)print*,'as1',real(itimestep)*dt,i,k,as(i,k),cs(i,k),((cs(i,k)**2)/(as(i,k)*nus*ns(i,k)))**(1./3.)
               IF(prds .gt. 0.0)THEN
                  snowdep(i,k) = prds
               ELSE
                  snowsub(i,k) = prds
               ENDIF

               !     if sublimation, reduce snow number,
               !     NOTE: this should not impact rhobars, since
               !     rhobars contains terms with qs/ns and ratio
               !     of qs/ns is assumed constant during loss of ns
               
               IF(qs(i,k).ge.qsmall.and.ns(i,k).gt.qsmall)THEN
               
                  IF(prds.lt.0.)THEN
                     ns(i,k)=ns(i,k)+prds*ns(i,k)/qs(i,k)*dt
                  END IF

                  !     set minimum ni to avoid taking root of a negative number
                  ns(i,k)=max(ns(i,k),qsmall)
                  if(isnan(ns(i,k)).or..not.ieee_is_finite(ns(i,k)))print*,'ns4',&
               real(itimestep)*dt,i,k,ns(i,k)
                  !     if iflag=1, then recalculate as and cs assuming same deltastrs
                  !     and same rhobars
                  !     this is needed for consistency between qs, as, cs, etc.
                  !     since growth rate prd must be scaled back
               
                  IF(iflag.eq.1)THEN
                     alphstr=co/ao**(deltastrs)
                     alphv=4./3.*pi*alphstr
                     betam=2.+deltastrs   
                     
                     ans=((qs(i,k)*exp(gammln(nus+betam)))/(rhos(i,k)*ns(i,k)*alphv*&
                          exp(gammln(nus+betam))))**(1./betam)
                     
                     cns=co*(ans/ao)**deltastrs
                  END IF           ! iflag = 1
                  
                  IF(isnan(rns))print*,"BEFORE-R-CHECK",rns
                  CALL R_CHECK(2,nus,qs(i,k),ns(i,k),cns,ans,rns,rhos(i,k),deltastrs)
                  IF(isnan(rns))print*,"AFTER-R-CHECK",rns
               
                  cs(i,k)=nus*ns(i,k)*cns
                  as(i,k)=nus*ns(i,k)*ans
                  !if(as(i,k).gt.0.0)print*,'as2',real(itimestep)*dt,i,k,as(i,k),cs(i,k),((cs(i,k)**2)/(as(i,k)*nus*ns(i,k)))**(1./3.)
               END IF              ! q > qsmall
            end if !snowflag

!     calculate simplified aggregation for snow category 
!     kjs 02/2015
            
!     First, calulcate the aggregate-available ice,
!     which is the amount of ice that would autoconvert to snow
!     in traditional schemes. This is for r_i > 125 um.
            betam=2.+deltastr
            alphstr=co/ao**(deltastr)
            alphv=4./3.*pi*alphstr 
            agg = 0.
            nagg = 0.
            IF(snowflag.eq.1)THEN
               IF(qi(i,k).ge.1.e-8.and.qv(i,k)/qvi.gt.1..and.rni.gt.0.)THEN
                  nprci = 4./DCS/rhobar*(qv(i,k)-qvi)*rho(i,k)*ni(i,k)&
                  /rni*exp(-DCS/rni)*dv/abi ! #/kg/s
                  prci = 4./3.*pi*rhobar*DCS**3*nprci !kg/kg/s
!     Second, if the aggregate-available ice mass mixing ratio is
!     > 1.e-8 kg/kg, then calculate total aggregate number and mass
!     and add it to the snow category.  The aggeagates are assumed 
!     to be spheres and have a density of rhosn=100 kg/m3.
                 
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
                      if (ipelec>0) then
                          tmp = an(i,1,k,lsc(li))/qi(i,k)
                          an(i,1,k,lsc(ls)) = an(i,1,k,lsc(ls)) - tmp*agg*dt
                          an(i,1,k,lsc(li)) = an(i,1,k,lsc(li)) + tmp*agg*dt
                     end if
                     qs(i,k) = qs(i,k) - agg*dt
                     ns(i,k) = ns(i,k) - nagg*dt

                     qi(i,k) = qi(i,k) + agg*dt
                     ni(i,k) = ni(i,k) + nagg*dt
                  END IF
                  
!adding collection rates for output
                  qagg(i,k) = agg
                  qnagg(i,k) = nagg
                  naggout(i,k) = nagg !#kg^-1 s^-1
                  nsaggout(i,k) = nsagg
                  prciout(i,k) = prci
               END IF
!.......................................................................
! ACCRETION OF CLOUD ICE BY SNOW
! FOR THIS CALCULATION, IT IS ASSUMED THAT THE VS >> VI
!     AND DS >> DI FOR CONTINUOUS COLLECTION
               
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

                  if (ipelec>0) then
                      tmp = an(i,1,k,lsc(li))/qi(i,k)
                      an(i,1,k,lsc(ls)) = an(i,1,k,lsc(ls)) + tmp*prai*dt
                      an(i,1,k,lsc(li)) = an(i,1,k,lsc(li)) - tmp*prai*dt
                  end if
                  qs(i,k) = qs(i,k) + prai*dt
                  ns(i,k) = ns(i,k) + nprai*dt

                  qi(i,k) = qi(i,k) - prai*dt
                  ni(i,k) = ni(i,k) - nprai*dt

                  praiout(i,k) = prai
               END IF
            ELSE IF(snowflag.eq.2)THEN
               if(ni(i,k).gt.qsmall)then
                  
                  nidum = ni(i,k)*rho(i,k)  !#/kg --> #/m3
                  CALL COLL_LOOKUP(nidum,ani,cni,rhobar,agg,nagg,aagg,anagg,cnagg,ddagg)
                  !CALL CPU_TIME(start)
                  !CALL COLL_COMPUTE(nidum,ani,cni,rhobar,nu,ns(i,k)*rho(i,k),ans,cns,rhobars,nus)
                  !CALL CPU_TIME(finish)
                  !print '("Time = ",f6.3," seconds.")',start,finish,finish-start
                  
                  if (effflag.eq.1)then
                      !!the aggregation efficiency used in the lookup table is 0.05 (5%), so need to divide out this
                      !!value and raplace with temp dependant EII from Connolly et al 2012
                      if(celsius .ge. -5.0)then 
                          iieff = EII_C(1)
                      else if(celsius .le. -30.0)then
                          iieff = EII_C(6)
                      else !celsius between -5 and -30
                          eii1 = abs(int(celsius/-5.0))
                          eii2 = abs(int(celsius/-5.0)) + 1
                          weight = celsius/((eii1+eii2)*-5.0)
                          iieff = EII_C(eii1)*(1.-weight) + EII_C(eii2)*weight
                      end if

                      !print*,'eff',celsius,iieff,agg
                      agg = agg/0.05*iieff
                      nagg = nagg/0.05*iieff
                      !print*,'eff2',agg
                  end if


                  !print*,"AFTER LOOKUP", nidum,ani,cni,rhobar,agg,nagg,ns(i,k)
                  qs_new = agg*dt/rho(i,k)  !kg/m3/s --> kg/kg
                  ns_new = nagg*dt/rho(i,k) !#/m3/s --> #/kg
                  !if(qs_new.gt.qsmall)print*,'qs',qs_new,agg
                  !if(ns_new.gt.qsmall)print*,'ns',ns_new,nagg
                  
                  if(qs(i,k)+qs_new.gt.qsmall.and.ns(i,k)+ns_new.gt.qsmall)then
                     ans_new = anagg*1.e-6
                     cns_new = cnagg*1.e-6
                     as_new = aagg*1.e-6*ns_new
                     dd_new = min(max((1.+ddagg)*rhobar,0.),920.)
                     !dd in IPAS is rho_newagg - rho_monomer (rho_newagg should be smaller)
                     !so is the percentage increase or decrease of new agg relative to original density,
                     !which in this case is the density of the monomer ice crystal
                     !so to get agg density, multiply the ice density by (1+%change),
                     !%change will most likely be negative
                     !print*,'density',dd_new

                     frac = 0.
                     if(ns(i,k)+ns_new.gt.qsmall) frac = ns(i,k)/(ns(i,k)+ns_new)
                     if(qs(i,k)+qs_new.gt.qsmall) frac = qs(i,k)/(qs(i,k)+qs_new)

                     if (ipelec>0.and.qi(i,k).gt.qsmall) then
                         tmp = an(i,1,k,lsc(li))/qi(i,k) ! charge per mass of ice
                         an(i,1,k,lsc(ls)) = an(i,1,k,lsc(ls)) + tmp*agg*dt/rho(i,k)
                         an(i,1,k,lsc(li)) = an(i,1,k,lsc(li)) - tmp*agg*dt/rho(i,k)
                     end if
                     qs(i,k) = qs(i,k) + qs_new
                     ns(i,k) = ns(i,k) + ns_new
                     
                     !if(ns(i,k).gt.10)print*,'ns>10',ns(i,k),ns_new,nagg
                     !if(as(i,k).gt.10.0)print*,'as before', as(i,k),as_new, nus
                     
                     as(i,k) = as(i,k)*frac + as_new*(1.-frac)
                     ans = ans*frac +  ans_new*(1.-frac)
                     cns = cns*frac +  cns_new*(1.-frac)
                     
                     !if(as(i,k).gt.0.0)print*,'as3',real(itimestep)*dt,i,k,as(i,k),cs(i,k),&
                     !((cs(i,k)**2)/(as(i,k)*nus*ns(i,k)))**(1./3.)

                     rhos(i,k) = rhos(i,k)*frac + dd_new*(1.-frac)
                     
                     if(ans.gt.0.)cs(i,k) = as(i,k)*cns/ans
                
                     if(ns(i,k)*ans.gt.0.)nus = max(as(i,k)/(ns(i,k)*ans),1.)
                     !if(ns(i,k)*ans.gt.0.)nus = min(max(as(i,k)/(ns(i,k)*ans),1.),8.)
                     IF(nus .gt. 10)print*,'nus',nus,as(i,k),ns(i,k),ans, as_new, frac,&
                     qs_new,qs(i,k)

                     IF(isnan(rns))print*,"BEFORECHECKAGAIN",rns
                     CALL DSTR_CHECK(2,nus,ns(i,k),ans,cns,deltastrs,iaspect,sphrflag)
                     CALL R_CHECK(2,nus,qs(i,k),ns(i,k),cns,ans,rns,rhos(i,k),deltastrs)
                     IF(isnan(rns))print*,"AFTERCHECKAGAIN",rns

                     qi(i,k) = qi(i,k) - agg*dt/rho(i,k)
                     ni(i,k) = ni(i,k) - nagg*dt/rho(i,k)

                     qagg(i,k) = agg/rho(i,k)
                     qnagg(i,k) = nagg/rho(i,k)
  
                  end if

               end if
               
            END IF!snowflag

            IF(qi(i,k).gt.qsmall.and.ni(i,k).gt.qsmall)THEN
               ani = ((qi(i,k)*gammnu)/(rhobar*ni(i,k)*alphv*&
               exp(gammln(nu+betam))))**(1./betam)
               cni=co*(ani/ao)**deltastr
               ci(i,k)=nu*ni(i,k)*cni
               ai(i,k)=nu*ni(i,k)*ani  
            END IF
            ns(i,k)=max(ns(i,k),qsmall)
            IF(qs(i,k).lt.qsmall.or.ns(i,k).lt.qsmall)THEN
               qs(i,k) = 0.0
               ns(i,k) = 0.0
               IF(snowflag.eq.2)then
                  as(i,k) = 0.0
                  cs(i,k) = 0.0
               END IF
               if (ipelec>0) then
                   an(i,1,k,lsc(ls)) = 0.0
               end if
            END IF
!..................................................................
!     now update qi, ni, ai, and ci due to ice nucleation
!     for simplicity, assume that rhobar and deltastr are 
!     constant during nucleation (probably not the best assumption,
!     should look into future modifications....)


!     make sure not to evaporate more liquid than available
            IF(qc(i,k).ge.qsmall)THEN
               IF((mnuccc+mnucci)*dt .lt. qc(i,k))then
                  qi(i,k)=qi(i,k)+(mnuccc+mnucci)*dt
                  ni(i,k)=ni(i,k)+(nnuccc+nnucci)*dt
                  if (ipelec>0) then
                      tmp = an(i,1,k,lsc(lc))/qc(i,k) ! charge per mass of water
                      an(i,1,k,lsc(li)) = an(i,1,k,lsc(li)) + tmp*(mnuccc+mnucci)*dt
!                      an(i,1,k,lsc(lc)) = an(i,1,k,lsc(lc)) - tmp*(mnuccc+mnucci)*dt ! done later
                  end if
               ELSE
                  qi(i,k)=qi(i,k)+qc(i,k)
                  ni(i,k)=ni(i,k)+nc(i,k)
                  if (ipelec>0) then
                      an(i,1,k,lsc(li)) = an(i,1,k,lsc(li)) + an(i,1,k,lsc(lc))
                  end if
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
               cni=co*(ani/ao)**deltastr
               ci(i,k)=nu*ni(i,k)*cni
               ai(i,k)=nu*ni(i,k)*ani

               ai(i,k)=max(ai(i,k),qsmall)
               ci(i,k)=max(ci(i,k),qsmall)

               ani=ai(i,k)/(nu*ni(i,k))
               cni=ci(i,k)/(nu*ni(i,k))

!     final check on limits for deltastr, rhobar
               CALL DSTR_CHECK(1,nu,ni(i,k),ani,cni,deltastr,iaspect,sphrflag)
               CALL RHO_CHECK(1,nu,deltastr,qi(i,k),ni(i,k),ani,cni,sphrflag,redden,rhobar)
               ci(i,k)=nu*ni(i,k)*cni
               ai(i,k)=nu*ni(i,k)*ani

            END IF ! qi > qsmall

!     add tendencies to temp, water vapor
!     if qi < qsmall, then zero out ni, ai, ci
            IF(qi(i,k).lt.qsmall)THEN
               qi(i,k)=0.
               ai(i,k)=0.
               ni(i,k)=0.
               ci(i,k)=0.
               rhobar=920.
               if (ipelec>0) then
                   an(i,1,k,lsc(li)) = 0.0
               end if
            ELSE
               qv(i,k)=qv(i,k)-(prd+mnuccd)*dt
               temp=temp+((prd+mnuccd)*xxls/cpm*dt)
               icenuc(i,k) = icenuc(i,k) + mnuccd
               dep(i,k) = mnuccd
               cfrzr(i,k) = mnuccr
               qiold = qi(i,k)
               IF(qc(i,k).ge.qsmall)THEN
                  IF((mnuccc+mnucci)*dt .lt. qc(i,k))then
                     if (ipelec>0) then
                         tmp = an(i,1,k,lsc(lc))/qc(i,k) ! charge per mass of water
                         an(i,1,k,lsc(lc)) = an(i,1,k,lsc(lc)) - tmp*(mnuccc+mnucci)*dt
                     end if
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
                     if (ipelec>0) then
                         an(i,1,k,lsc(lc)) = 0.0
                     end if
                  END IF
               END IF
            END IF
            
            
            !update to ice and snow after melting
            IF(meltflag .eq. 1)THEN
                IF(temp.ge.273.15)THEN

                    IF(qi(i,k).gt.qsmall)THEN
                       rhobar = rhobar*((qi(i,k)+(qmlt*dt))/qi(i,k)) + &
                            RHOW*((-qmlt*dt)/qi(i,k))
                       rhobar = min(rhobar,RHOI)
                    ELSE
                       rhobar = RHOI
                    END IF

                    if (ipelec>0.and.qi(i,k).ge.qsmall) then
                        tmp = an(i,1,k,lsc(li))/qi(i,k) ! charge per mass of ice
                        an(i,1,k,lsc(li)) = an(i,1,k,lsc(li)) + tmp*qmlt*dt
                        an(i,1,k,lsc(lr)) = an(i,1,k,lsc(lr)) - tmp*qmlt*dt
                    end if
                    qi(i,k) = qi(i,k) + qmlt*dt
                    ni(i,k) = ni(i,k) + nmlt*dt
                    ai(i,k) = ai(i,k) + amlt*dt
                    ci(i,k) = ci(i,k) + cmlt*dt
                    
                    qr(i,k)=qr(i,k)-qmlt*dt
                    nr(i,k)=nr(i,k)-nmlt*dt

                    !.. Small masses of leftover ice can become rain
                    if(qi(i,k).le.QSMALL) then
                       qr(i,k) = qr(i,k) + qi(i,k)
                       qi(i,k)=0.
                       if (ipelec>0) then
                           an(i,1,k,lsc(lr)) = an(i,1,k,lsc(lr)) + an(i,1,k,lsc(li))
                           an(i,1,k,lsc(li)) = 0.0
                       end if
                    endif

                    IF(qi(i,k).gt.QSMALL)THEN

                        ni(i,k)=max(ni(i,k),QSMALL)
                        ai(i,k)=max(ai(i,k),QSMALL)
                        ci(i,k)=max(ci(i,k),QSMALL)

                        ani=ai(i,k)/(nu*ni(i,k))
                        cni=ci(i,k)/(nu*ni(i,k))

                        dsold = deltastr
                        deltastr = (log(cni)-log(ao))/(log(ani)-log(ao))

                        if(dsold.le.1.) then
                           if(deltastr.lt.dsold) then
                              deltastr=dsold
                           endif
                           if(deltastr.gt.1.) then
                              deltastr=dsold
                           endif
                        else
                           if(deltastr.gt.dsold) then
                              deltastr=dsold
                           endif
                           if(deltastr.lt.1.) then
                              deltastr=dsold
                           endif
                        endif
                        cni=ao**(1.-deltastr)*ani**deltastr

                        CALL ICE_CHECKS(1, nu,ni(i,k),qi(i,k),ani,cni,rni,deltastr,rhobar,&
                        iaspect,sphrflag,redden)

                    END IF !qi > qsmall

                    IF(snowflag.eq.2)THEN
                        IF(qs(i,k).gt.qsmall)THEN
                           rhos(i,k) = rhos(i,k)*((qs(i,k)+(qsmlt*dt))/qs(i,k)) + &
                                RHOW*((-qsmlt*dt)/qs(i,k))
                           rhos(i,k) = min(rhos(i,k),RHOI)
                        ELSE
                           rhos(i,k) = RHOI
                        END IF

                        if (ipelec>0.and.qs(i,k).ge.qsmall) then
                            tmp = an(i,1,k,lsc(ls))/qs(i,k) ! charge per mass of snow
                            an(i,1,k,lsc(ls)) = an(i,1,k,lsc(ls)) + tmp*qsmlt*dt
                            an(i,1,k,lsc(lr)) = an(i,1,k,lsc(lr)) - tmp*qsmlt*dt
                        end if
                        qs(i,k) = qs(i,k) + qsmlt*dt
                        ns(i,k) = ns(i,k) + nsmlt*dt
                        as(i,k) = as(i,k) + asmlt*dt
                        cs(i,k) = cs(i,k) + csmlt*dt
                        
                        qr(i,k)=qr(i,k)-qsmlt*dt
                        nr(i,k)=nr(i,k)-nsmlt*dt


                        !.. Small masses of leftover snow can become rain
                        if(qs(i,k).le.QSMALL) then
                           qr(i,k) = qr(i,k) + qs(i,k)
                           qs(i,k)=0.
                           if (ipelec>0) then
                               an(i,1,k,lsc(lr)) = an(i,1,k,lsc(lr)) + an(i,1,k,lsc(ls))
                               an(i,1,k,lsc(ls)) = 0.0
                           end if
                        endif

                        IF(qs(i,k).gt.QSMALL)THEN

                            ns(i,k)=max(ns(i,k),QSMALL)
                            as(i,k)=max(as(i,k),QSMALL)
                            cs(i,k)=max(cs(i,k),QSMALL)

                            ans=as(i,k)/(nus*ns(i,k))
                            cns=cs(i,k)/(nus*ns(i,k))

                            dsold = deltastrs
                            deltastrs = (log(cns)-log(ao))/(log(ans)-log(ao))

                            if(dsold.le.1.) then
                               if(deltastrs.lt.dsold) then
                                  deltastrs=dsold
                               endif
                               if(deltastrs.gt.1.) then
                                  deltastrs=dsold
                               endif
                            else
                               if(deltastrs.gt.dsold) then
                                  deltastrs=dsold
                               endif
                               if(deltastrs.lt.1.) then
                                  deltastrs=dsold
                               endif
                            endif
                            cns=ao**(1.-deltastrs)*ans**deltastrs

                            CALL ICE_CHECKS(1, nus,ns(i,k),qs(i,k),ans,cns,rns,deltastrs,rhos(i,k),&
                            iaspect,sphrflag,redden)

                        END IF !qs > qsmall
                    END IF !snowflag
                END IF !temp>273
            END IF

            IF(qs(i,k).lt.qsmall)THEN
               qs(i,k)=0.
               ns(i,k)=0.
               IF(snowflag.eq.2)then
                  as(i,k) = 0.0
                  cs(i,k) = 0.0
               	END IF
               if (ipelec>0) then
                   an(i,1,k,lsc(ls)) = 0.0
               end if
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
            


            ns(i,k) = max(ns(i,k),qsmall)
            IF(qs(i,k).gt.qsmall.and.ns(i,k).gt.qsmall) THEN

!     limit fallspeed to 5 m/s
               vtrms(i,k)=max(min(vtbarbms,5.),0.)
               vtrns(i,k)=max(min(vtbarbs,5.),0.)
               vtrls(i,k)=max(min(vtbarblens,5.),0.)

            ELSE
               
               vtrls(i,k)=0.
               vtrms(i,k)=0.
               vtrns(i,k)=0.
               
            END IF



!     get phi
            
            IF(qi(i,k).ge.qsmall.and.ani.gt.qsmall)THEN
               phii = cni/ani*exp(gammln(nu+deltastr-1.))/gammnu
            ELSE
               phii=1.
               ani=0.
               cni=0.
               rhobar=920.
            END IF
            IF(snowflag.eq.2)THEN
               IF(qs(i,k).ge.qsmall.and.ans.gt.qsmall)THEN
                  phis(i,k) = cns/ans*exp(gammln(nus+deltastrs-1.))/exp(gammln(nus))
               ELSE
                  phis(i,k)=1.
                  ans=0.
                  cns=0.
                  rhos(i,k)=920.
               END IF
            END IF
            IF(iaspect .eq. 1) phii = 0.27
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
            IF(snowflag.eq.1)THEN
               qs(i,k) = qs(i,k) + &
                      (psmlt+evpms-pracs)*dt +&
                      (psacws+prds+pracs1+eprds-psacr+piacrs+&
                      pracis)*dt
               ns(i,k) = ns(i,k) + &
                      nsmlts*dt + &
                      (nsagg-nscng-ngracs+niacrs+nsubs)*dt
            END IF
            qg(i,k) = qg(i,k) + &
                      (pgmlt+evpmg-pracg)*dt + &
                      (pracg1+psacwg+pgsacw+pgracs+prdg+&
                      eprdg+piacr+psacr &
                              +praci & ! yc added 6/22/2022
                              )*dt
            ng(i,k) = ng(i,k) + &
                      (ngmltg)*dt + &
                      (nscng+ngracs+niacr+nsubg)*dt

            qitend = (qmults+qmultg+qmultr+qmultrg-praci-pracis)
            nitend = (nmults+nmultg+nmultr+nmultrg-niacr-niacrs)
            qi(i,k) = qi(i,k) + qitend*dt
            ni(i,k) = ni(i,k) + nitend*dt

            ! ⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡
            ! yc added 6/30/2022. for elec
            IF (IPELEC>0) THEN
                ! update an "q and n" part (C/m3 for charge)
                an(i,1,k,lc) = qc(i,k) ! kg/kg
                an(i,1,k,lr) = qr(i,k)
                an(i,1,k,li) = qi(i,k)
                an(i,1,k,ls) = qs(i,k)
                an(i,1,k,lh) = qg(i,k)

                an(i,1,k,lnc) = nc(i,k)
                an(i,1,k,lnr) = nr(i,k)
                an(i,1,k,lni) = ni(i,k)
                an(i,1,k,lns) = ns(i,k)
                an(i,1,k,lnh) = ng(i,k)

                DO il = lc,lhab
                    scx(i,k,il) = an(i,1,k,lsc(il))
                    IF ( Abs(scx(i,k,il)) > 1000.e-9 ) THEN
                    write(0,*) 'Problem0a with scx il = ',il,lsc(il),scx(i,k,il),an(i,1,k,il),an(i,1,k,ln(il))
                    write(0,*) 'at it, ix,jy,kz = ',itimestep ,i, j, k
                    ENDIF
                ENDDO
                cionp(i,k) = an(i,1,k,lscpi) ! #/m3
                cionn(i,k) = an(i,1,k,lscni)
                if ( largeion ) then
                    clionp(i,k) = an(i,1,k,lscpli)
                    clionn(i,k) = an(i,1,k,lscnli)
                end if


                ! calculate size & fall speed
                dum = (rhosu/rho(i,k))**0.54
                ! SNOW---------------------------------------------------------------------------------
                IF(qs(i,k).ge.qsmall.and.ns(i,k).gt.qsmall)THEN
                    if (snowflag.eq.1) then
                        lamselec = (cons1*ns(i,k)/qs(i,k))**(1./DS)
                        n0selec = ns(i,k)*lamselec
                        ! only updating lams that is needed, in order not to change ns before it gets updated later
                        IF(lamselec.lt.lammins)THEN
                            lamselec = lammins
                            n0selec = lamselec**4*qs(i,k)/cons1
                            an(i,1,k,lns) = n0selec/lamselec
                        ELSE IF(lamselec.gt.lammaxs)THEN
                            lamselec = lammaxs
                            n0selec = lamselec**4*qs(i,k)/cons1
                            an(i,1,k,lns) = n0selec/lamselec
                        END IF
                        xdia(i,k,ls) = 1.0/lamselec

                        fspeed = asn(k)*cons3/((cons1*ns(i,k)/qs(i,k))**(1./DS))**BS
                        umx(i,k,ls) = min(fspeed,1.2*dum)
                    else if (snowflag.eq.2) then
                        anselec = as(i,k)/(nus*ns(i,k))
                        cnselec = cs(i,k)/(nus*ns(i,k))
                        if (snowdiaflag.eq.1) then
                            xdia(i,k,ls) = 2.*anselec
                        else if (snowdiaflag.eq.2) then
                            xdia(i,k,ls) = 2.*cnselec
                        else if (snowdiaflag.eq.3) then
                            xdia(i,k,ls) = 2.*max(anselec,cnselec)
                        else if (snowdiaflag.eq.4) then ! eq volumne radius
                            rhoselec = rhos(i,k)
                            deltastrselec = deltastrs
                            CALL ICE_CHECKS(1,nus,an(i,1,k,lns),an(i,1,k,ls),anselec,cnselec,rnselec,deltastrselec,rhoselec,&
                                    iaspect,sphrflag,redden)
                            xdia(i,k,ls) = 2.*rnselec
                        end if

                        umx(i,k,ls) = vtrms(i,k)
                    end if
                END IF

                ! GRAUPEL---------------------------------------------------------------------------------
                if(qg(i,k).ge.qsmall .and. ng(i,k).ge.qsmall)then
                    ! only updating lamg that is needed, in order not to change ng before it gets updated later
                    lamgelec = (cons2*ng(i,k)/qg(i,k))**(1./DG)
                    n0gelec = ng(i,k)*lamgelec
                    if(lamgelec.lt.lamming)then
                        lamgelec = lamming
                        n0gelec = lamgelec**4*qg(i,k)/cons2
                        an(i,1,k,lnh) = n0gelec/lamgelec
                    else if(lamgelec.gt.lammaxg)then
                        lamgelec = lammaxg
                        n0gelec = lamgelec**4*qg(i,k)/cons2
                        an(i,1,k,lnh) = n0gelec/lamgelec
                    end if
                    xdia(i,k,lh) = 1.0/lamgelec

                    fspeed = agn(k)*cons7/(lamgelec**BG)
                    umx(i,k,lh) = min(fspeed,20.*dum)
                end if

                ! CLOUD DROPLET---------------------------------------------------------------------------------
                IF(qc(i,k).ge.qsmall .and. nc(i,k).ge.qsmall)THEN ! here using number weighted
                    pgam = 0.0005714*(nc(i,k)/1.e6* (press/(287.15*temp)) )+0.2714
                    pgam = 1./(pgam**2)-1.
                    pgam = min(max(pgam,2.),10.)
                    lamc = (cons26*nc(i,k)*GAMMA(pgam+4.)/&
                                (qc(i,k)*GAMMA(pgam+1.)))**(1./3.)
                    lammin = (pgam+1.)/60.e-6
                    lammax = (pgam+1.)/1.e-6
                    IF(lamc.lt.lammin)THEN
                        lamc = lammin
                        an(i,1,k,lnc) = exp(3.*log(lamc)+log(qc(i,k))+&
                                log(GAMMA(pgam+1.))-log(GAMMA(pgam+4.)))/cons26
                    ELSE IF(lamc.gt.lammax)THEN
                        lamc = lammax
                        an(i,1,k,lnc) = exp(3.*log(lamc)+log(qc(i,k))+&
                                log(GAMMA(pgam+1.))-log(GAMMA(pgam+4.)))/cons26
                    END IF

                    xdia(i,k,lc) = (6/(PI*RHOW)*qc(i,k)/nc(i,k))**(1./3.)
                    xdia(i,k,lc) = min(max(xdia(i,k,lc),1.E-6),60.E-6)

                    fspeed = acn(k)*GAMMA(4.+BC+PGAM)/(lamc**BC*GAMMA(PGAM+4.))
                    umx(i,k,lc) = min(fspeed,9.1*dum)
                END IF

                ! CLOUD ICE---------------------------------------------------------------------------------
                IF(qi(i,k).gt.qsmall.and.ni(i,k).gt.qsmall) THEN
                    if (sphrflag.eq.1) then
                        rhoielec = 920.
                        If(redden .eq. 1) rhoielec = 500.
                        xdia(i,k,lc) = (6/(PI*rhoielec)*qi(i,k)/ni(i,k))**(1./3.)
                        xdia(i,k,lc) = min(max(xdia(i,k,lc),1.E-6),(2.*125.E-6+100.E-6)) ! adapted from morr
                    else
                        anielec=ai(i,k)/(nu*ni(i,k))
                        cnielec=ci(i,k)/(nu*ni(i,k))
                        if (icediaflag.eq.1) then
                            xdia(i,k,li) = 2.*anielec
                        else if (icediaflag.eq.2) then
                            xdia(i,k,li) = 2.*cnielec
                        else if (icediaflag.eq.3) then
                            xdia(i,k,li) = 2.* max(anielec,cnielec)
                        else if (icediaflag.eq.4) then
                            rhobarelec = rhobar
                            deltastrelec = deltastr
                            CALL ICE_CHECKS(1,nu,an(i,1,k,lni),an(i,1,k,li),anielec,cnielec,rnielec,deltastrelec,rhoielec,&
                                    iaspect,sphrflag,redden)
                            xdia(i,k,li) = 2.*rnielec
                        end if
                    end if

                    umx(i,k,li) = vtrmi1(i,k)
                END IF

                IF(qr(i,k).gt.qsmall.and.nr(i,k).gt.qsmall) THEN
                    fspeed = arn(k)*cons4/(lamr**BR)
                    umx(i,k,lr) = min(fspeed,9.1*dum)
                END IF

                an(i,1,k,lnc) = an(i,1,k,lnc) * rho(i,k) ! #/m3 consistent w/ nssl
                an(i,1,k,lnr) = an(i,1,k,lnr) * rho(i,k)
                an(i,1,k,lni) = an(i,1,k,lni) * rho(i,k)
                an(i,1,k,lns) = an(i,1,k,lns) * rho(i,k)
                an(i,1,k,lnh) = an(i,1,k,lnh) * rho(i,k)

                !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
                ! Collision rate of snow - ice (CSACI0: #/m3/s)
                IF (an(i,1,k,ls).gt.qsmall.and.an(i,1,k,lns).gt.qsmall .and. &
                    an(i,1,k,li).ge.qsmall.and.an(i,1,k,lni).gt.qsmall) THEN
                    if (snowflag.eq.1) then
                        csaci0(i,k) = pi/4.*gamma(bs+3.)*asn(k)*an(i,1,k,lni)*n0selec/lamselec**(bs+3.)*rho(i,k) ! same eq as for nprai except no eii
                    else if (snowflag.eq.2) then
                        csaci0(i,k) = pi/4.*xdia(i,k,ls)*umx(i,k,ls)*an(i,1,k,lni)*an(i,1,k,lns)
                    end if
                END IF

                !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
                ! Collision rate of graupel-ice (CHACI0)
                IF (an(i,1,k,lh).gt.qsmall.and.an(i,1,k,lnh).gt.qsmall .and. &
                    an(i,1,k,li).gt.qsmall.and.an(i,1,k,lni).gt.qsmall) THEN
                    vt = sqrt((umx(i,k,lh)-umx(i,k,li))**2 + 0.04*umx(i,k,lh)*umx(i,k,li))
                    EHICLSN = EHI_COLLSN ! Graupel-ice collision efficiency, currently used as constant 1.
                    CHACI0(i,K) = 0.25 * PI * EHICLSN * &
                            (an(i,1,k,lnh)) * (an(i,1,k,lni)) * VT * &
                            (DA0LH * xdia(i,k,lh)**2 + &
                                    DAB0LH_IH * xdia(i,k,lh) * xdia(i,k,li) + &
                                    DA0LI * xdia(i,k,li)**2)
                END IF

                !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
                ! Collision rate of graupel-snow (CHACS0)
                IF (an(i,1,k,lh).gt.qsmall.and.an(i,1,k,lnh).gt.qsmall .and. &
                    an(i,1,k,ls).ge.qsmall.and.an(i,1,k,lns).gt.qsmall) THEN
                    vt = sqrt((umx(i,k,lh)-umx(i,k,ls))**2 + 0.04*umx(i,k,lh)*umx(i,k,ls))

                    ! Graupel-snow collision efficiency.
                    ! Defined as a linear function of snow size.
                    ! Zero at Ds<=40um, 0.5 at Ds>=150um
                    EHSCLSN = EHS_COLLSN
                    IF ( xdia(i,k,ls) < 40.E-6 ) THEN
                        EHSCLSN = 0.0
                    ELSEIF ( xdia(i,k,ls) < 150.E-6 ) THEN
                        EHSCLSN =  EHS_COLLSN*(xdia(i,k,ls) - 40.E-6)/(150.E-6 - 40.E-6)
                    ELSE
                        EHSCLSN = EHS_COLLSN
                    ENDIF

                    CHACS0(i,K) = 0.25 * PI * EHSCLSN * &
                            (an(i,1,k,lnh)) * (an(i,1,k,lns)) * VT * &
                            (DA0LH * xdia(i,k,lh)**2 + &
                                    DAB0LH_SH * xdia(i,k,lh) * xdia(i,k,ls) + &
                                    DA0LS * xdia(i,k,ls)**2)
                END IF

                cxacy(i,k,ls,li) = csaci0(i,k)
                cxacy(i,k,lh,li) = chaci0(i,k)
                cxacy(i,k,lh,ls) = chacs0(i,k)

                !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
                ! Collisional charging rates
                ! Sticking eff. Collision eff assumed 1
                exy(ls,lc) = eci ! 0.7
                if (snowflag.eq.2) then
                exy(ls,li) = 0.0 ! since no ice-snow collection considered
                else if (snowflag.eq.1) then
                    exy(ls,li) = EII ! same as collection eff since assumed 100% collision
                end if
                exy(lh,lc) = eci ! 0.7
                exy(lh,li) = 0.0 ! since no graup-ice collection considered
                exy(lh,ls) = 0.0 ! since no graup-snow collection considered

                temcg(i,k) = t(i,k)-273.15

                IF ( ipelec .ge. 2 .or. idonic ) THEN
                    IF ( j.gt.(jds+2) .and. j.lt.(jde-2) ) THEN
                        DO il=ls,lhab ! Collector: loop through all ice habits (s,h) Collector
                            DO ic=li,ls ! Collectee: loop through all mixed phase particles (i,s)
                                ! il=lh
                                ! ic=li
                                IF ( iexy(il,ic) .eq. 1  ) THEN ! (si,hi,hs)
                                    IF (an(i,1,k,il) .gt. qxmin(il)  .and. an(i,1,k,ic) .gt. qxmin(ic) &
                                            .and. an(i,1,k,ln(il)).gt.qsmall .and. an(i,1,k,ln(ic)).gt.qsmall) THEN  ! make sure there are particles
                                        IF ( t(i,k) .lt. tfr .and. t(i,k) .gt. thnuc ) THEN ! 0.0=tfreezing and -38C= assumed T for homogeneous ice nucleation
                                            IF ( iremoveqwfrz == 1 ) THEN
                                                qcwtmp = Max(0.0, an(i,1,k,lc) - dtp*MNUCCC)
                                            ELSE
                                                qcwtmp = an(i,1,k,lc)
                                            ENDIF

                                            ! Currently set nonigrd = 0 (saundx), default isaund = 12
                                            IF (nonigrd .eq. -1 ) THEN
                                                vt = umx(i,k,il) - umx(i,k,ic)
                                                CALL takax(isaund,ntem,nlwc,takalu,temcg(i,k),qcwtmp,vt,   &
                                                        &     xdia(i,k,ic),rho(i,k),ftelwc,exy(il,lc),an(i,1,k,lr),rarfac)
                                                scxacy(i,k,il,ic) = ftelwc
                                                if ( scxacy(i,k,il,ic) .gt. 0.0 ) then
                                                    scxacy(i,k,il,ic) = min(delqxxa(ic),scxacy(i,k,il,ic))
                                                end if
                                                if ( scxacy(i,k,il,ic) .lt. 0.0 ) then
                                                    scxacy(i,k,il,ic) = max(delqnxa(ic),scxacy(i,k,il,ic))
                                                end if
                                            ELSEIF ( nonigrd .eq. 1 .or. nonigrd .eq. 0 ) THEN
                                                IF (nonigrd .eq. 0) THEN ! currently used
                                                    qconkq = 0.0
                                                    qconm = 1.0
                                                    qconn = 1.0
                                                    ftelwc = 0.0
                                                    vt = umx(i,k,il) - umx(i,k,lc)
                                                    call saundx(isaund,temcg(i,k),qcwtmp,exy(il,lc),vt,xdia(i,k,ic),   &
                                                            &   rho(i,k),qsign,ftelwc,qconkq,qconm,qconn,idelq,rarfac)
                                                ELSE ! nonigrd=1
                                                    vt = umx(i,k,il) - umx(i,k,lc)
                                                    ! rarx only done for graup (as in nssl), doesn't matter since nonigrd is not 1 for now
                                                    rarx(i,k,lh) = (PSACWG*rho(i,k))*1.0e3*rho(i,k)/((pi/2.0)*xdia(i,k,lh)*an(i,1,k,lnh)) ! ! PSACWG: equate to qhacw cw collected by gw, cw->gw
                                                    call saundy(isaund,temcg(i,k),qcwtmp,rarx(i,k,il),vt,xdia(i,k,ic),   &
                                                            &   rho(i,k),qsign,ftelwc,qconkq,qconm,qconn,idelq)
                                                ENDIF

                                                cfce15 = 1.0e-15
                                                scxacy(i,k,il,ic) = qconkq*(xdia(i,k,ic)**qconm)   &
                                                        &  *((abs(umx(i,k,il)-umx(i,k,ic)))**qconn)   &
                                                        &  *ftelwc*cfce15
                                                if ( scxacy(i,k,il,ic) .gt. 0.0 ) then ! {
                                                    scxacy(i,k,il,ic) = min(delqxxa(ic),scxacy(i,k,il,ic))
                                                end if
                                                if ( scxacy(i,k,il,ic) .lt. 0.0 ) then
                                                    scxacy(i,k,il,ic) = max(delqnxa(ic),scxacy(i,k,il,ic))
                                                end if
                                            end if

                                            ! Formula 7 in M05 to get total space charge per collisions within a grid box
                                            tmp = scxacy(i,k,il,ic)
                                            scxacy(i,k,il,ic) =   &
                                                    &  cxacy(i,k,il,ic)*ecollmx*(1.0-exy(il,ic))*scxacy(i,k,il,ic)
                                            IF ( Abs( scxacy(i,k,il,ic) ) .gt. scxacymax ) THEN ! scxacymax=3000e-12
                                                scxacy(i,k,il,ic) = Min(scxacymax, Max( -scxacymax, scxacy(i,k,il,ic) ) )
                                            ENDIF

                                        END IF ! temperture allows mixed phase
                                    END IF ! an mass/number > qsmall
                                END IF ! iexy(il,ic) .eq. 1
                            END DO ! ic
                        END DO ! il
                        !
                        !  End of ice-ice loop
                        !
                    ELSE
                        scsacw(i,k) = 0.0
                        scsaci(i,k) = 0.0

                        schacw(i,k) = 0.0
                        schaci(i,k) = 0.0
                        schacs(i,k) = 0.0
                        scxacy(i,k,:,:)=0.0
                    END IF ! boundary check on y
                END IF ! ipelec>=2

                IF ( i.gt.(ids+2) .and. i.lt.(ide-2)  ) THEN
                    scsacw(i,k) = elecfac*scxacy(i,k,ls,lc) ! induc
                    scsaci(i,k) = elecfac*scxacy(i,k,ls,li) ! non-induc

                    schacw(i,k) = elecfac*scxacy(i,k,lh,lc) ! induc
                    schaci(i,k) = elecfac*scxacy(i,k,lh,li) ! non-induc
                    schacs(i,k) = elecfac*scxacy(i,k,lh,ls) ! non-induc
                ELSE
                    scsacw(i,k) = 0.0
                    scsaci(i,k) = 0.0

                    schacw(i,k) = 0.0
                    schaci(i,k) = 0.0
                    schacs(i,k) = 0.0
                    scxacy(i,k,:,:)=0.0
                END IF ! x boundary check

                !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
                ! Charge per hydrometeor content
                if ( ipelec .ge. 1) then
                    if ( an(i,1,k,lr) .gt. qxmin(lr) ) then
                        fscrw = scx(i,k,lr)/an(i,1,k,lr)
                    end if
                    if ( an(i,1,k,lc) .gt. qxmin(lc) ) then
                        fsccw = scx(i,k,lc)/an(i,1,k,lc)
                    end if
                    if ( an(i,1,k,li) .gt. qxmin(li) ) then
                        fscci = scx(i,k,li)/an(i,1,k,li)
                    end if
                    if ( an(i,1,k,ls) .gt. qxmin(ls) ) then
                        fscsw = scx(i,k,ls)/an(i,1,k,ls)
                    end if
                    if ( an(i,1,k,lh) .gt. qxmin(lh) ) then
                        fschw = scx(i,k,lh)/an(i,1,k,lh)
                    end if
                end if
                ! il5 =0, warm; =1, frozen
                il5 = 0
                if ( t(i,k) .lt. 273.15 ) then
                    il5 = 1
                end if

                !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
                ! Charge production rates (pscxxi,pscxxd,pscxxmi,pscxxmd)
                !!
                !!
                !!  Charge production from mass transfer
                !!
                !!
                !-------------- First get sink rates to ensure conservation
                !-------------- i.e., loss of charge due to mass loss shouldn't lead to sign reversal
                !-------------- i.e., pscxxmd shouldn't be more than scx/dt
                !
                !  cloud ice
                !
                if ( an(i,1,k,li) .gt. qxmin(li) ) then
!                    psccimi = fsccw*(qmults+qmultg) +&
!                            fscrw*(qmultr+qmultrg)
                    praci_e = praci
                    pracis_e = pracis
                    psccimd = fscci*(-praci-pracis)
                    if (-praci.gt.0.0 .or.&
                            -pracis.gt.0.0) then
                        write(0,*) 'sink terms pos for ice'
                        write(0,*) 'temper = ',T(i,k)
                        write(0,*) 'i,j,k,t=',i,j,k,itimestep
                        write(0,*) '-praci', -praci
                        write(0,*) '-pracis', -pracis
                    end if
                    if (psccimd*an(i,1,k,lsc(li)).gt.0.) then ! pscxxmd should lead to decrease of scx, so pscxxmd*scx should <=0
                        write(0,*) 'ice psccimd same sign as an(lsci)'
                        write(0,*) 'temper = ',T(i,k)
                        write(0,*) 'i,j,k,t=',i,j,k,itimestep
                        write(0,*) 'psccimd, an(i,1,k,lsc(li))',psccimd, an(i,1,k,lsc(li))
                        write(0,*) 'fscci,praci_e,pracis_e'
                        write(0,*) fscci,praci_e,pracis_e
                        write(0,*) '(-praci_e-pracis_e)'
                        write(0,*) (-praci_e-pracis_e)
                        write(0,*) 'qx,cx = ',an(i,1,k,li),an(i,1,k,ln(li))
                        write(0,*) 'stopping'
                        print*, 'stopping'
                        STOP
                    end if
                    if ((abs(psccimd)*dt) .gt. abs(an(i,1,k,lsc(li)))) then
                        tmp = abs(an(i,1,k,lsc(li)) / (psccimd*dt))
                        write(0,*) 'conserving ice sci, tmp= ',tmp
                        write(0,*) 'temper = ',T(i,k)
                        write(0,*) 'i,j,k,t=',i,j,k,itimestep
                        write(0,*) 'psccimd, dt',psccimd, dt
                        write(0,*) 'psccimd*dt, an(i,1,k,lsc(li))', psccimd*dt, an(i,1,k,lsc(li))
                        praci_e = tmp * praci_e
                        pracis_e = tmp * pracis_e
                        psccimd = tmp * psccimd
                        write(0,*) 'psccimd after conserving',psccimd, psccimd*dt
                    end if
                end if

                !
                !  cloud water
                !
                if ( an(i,1,k,lc) .gt. qxmin(lc) ) then
!                    psccwmi = 0.0
                    pra_e = pra
                    prc_e = prc
                    psacws_e = psacws
                    qmults_e = qmults
                    qmultg_e = qmultg
                    psacwg_e = psacwg
                    pgsacw_e = pgsacw
                    psccwmd = fsccw*(-pra-prc-psacws-qmults-qmultg-psacwg-pgsacw)
                    if (-pra.gt.0.0 .or.&
                        -prc.gt.0.0 .or.&
                        -psacws.gt.0.0 .or.&
                        -qmults.gt.0.0 .or.&
                        -qmultg.gt.0.0 .or.&
                        -psacwg.gt.0.0 .or.&
                        -pgsacw.gt.0.0) then
                        write(0,*) 'sink terms pos for cloud'
                        write(0,*) 'temper = ',T(i,k)
                        write(0,*) 'i,j,k,t=',i,j,k,itimestep
                        write(0,*) '-pra',-pra
                        write(0,*) '-prc',-prc
                        write(0,*) '-psacws',-psacws
                        write(0,*) '-qmults',-qmults
                        write(0,*) '-qmultg',-qmultg
                        write(0,*) '-psacwg',-psacwg
                        write(0,*) '-pgsacw',-pgsacw
                    end if
                    if (psccwmd*an(i,1,k,lsc(lc)).gt.0.) then ! pscxxmd should lead to decrease of scx, so pscxxmd*scx should be neg
                        write(0,*) 'cloud psccwmd same sign as an(lscw)'
                        write(0,*) 'temper = ',T(i,k)
                        write(0,*) 'i,j,k,t=',i,j,k,itimestep
                        write(0,*) 'psccwmd, an(i,1,k,lsc(lc))',psccwmd, an(i,1,k,lsc(lc))
                        write(0,*) 'fsccw,pra,prc,psacws,qmults,qmultg,psacwg,pgsacw'
                        write(0,*) fsccw,pra,prc,psacws,qmults,qmultg,psacwg,pgsacw
                        write(0,*) '(-pra_e-prc_e-psacws_e-qmults_e-qmultg_e-psacwg_e-pgsacw_e)'
                        write(0,*) (-pra_e-prc_e-psacws_e-qmults_e-qmultg_e-psacwg_e-pgsacw_e)
                        write(0,*) 'qx,cx = ',an(i,1,k,lc),an(i,1,k,ln(lc))
                        write(0,*) 'stopping'
                        print*, 'stopping'
                        STOP
                    end if
                    if ((abs(psccwmd)*dt) .gt. abs(an(i,1,k,lsc(lc)))) then
                        tmp = abs(an(i,1,k,lsc(lc)) / (psccwmd*dt))
                        write(0,*) 'conserving cloud scw, tmp= ',tmp
                        write(0,*) 'temper = ',T(i,k)
                        write(0,*) 'i,j,k,t=',i,j,k,itimestep
                        write(0,*) 'psccwmd, dt',psccwmd, dt
                        write(0,*) 'psccwmd*dt, an(i,1,k,lsc(lc))', psccwmd*dt, an(i,1,k,lsc(lc))
                        pra_e = tmp * pra_e
                        prc_e = tmp * prc_e
                        psacws_e = tmp * psacws_e
                        qmults_e = tmp * qmults_e
                        qmultg_e = tmp * qmultg_e
                        psacwg_e = tmp * psacwg_e
                        pgsacw_e = tmp * pgsacw_e
                        psccwmd = tmp * psccwmd
                        write(0,*) 'psccwmd after conserving',psccwmd, psccwmd*dt
                    end if
                end if

                !
                !  snow
                !
                if ( an(i,1,k,ls) .gt. qxmin(ls) ) then
!                    pscswmi = fsccw*(psacws) +&
!                            fscrw*(pracs1+piacrs+max(0.,-pracs)) +&
!                            fscci*(pracis)
                    psmlt_e = psmlt
                    psacr_e = psacr
                    pscswmd = fscsw*(min(0.,psmlt)-psacr)
                    if (min(0.,psmlt).gt.0.0 .or.&
                            -psacr.gt.0.0) then
                        write(0,*) 'sink terms pos for snow'
                        write(0,*) 'temper = ',T(i,k)
                        write(0,*) 'i,j,k,t=',i,j,k,itimestep
                        write(0,*) 'min(0.,psmlt)', min(0.,psmlt)
                        write(0,*) '-psacr', -psacr
                    end if
                    if (pscswmd*an(i,1,k,lsc(ls)).gt.0.) then ! pscxxmd should lead to decrease of scx, so pscxxmd*scx should be neg
                        write(0,*) 'snow pscswmd same sign as an(lscs)'
                        write(0,*) 'temper = ',T(i,k)
                        write(0,*) 'i,j,k,t=',i,j,k,itimestep
                        write(0,*) 'pscswmd, an(i,1,k,lsc(ls))',pscswmd, an(i,1,k,lsc(ls))
                        write(0,*) 'fscsw,psmlt,psacr'
                        write(0,*) fscsw,psmlt,psacr
                        write(0,*) '(min(0.,psmlt_e)-psacr_e)'
                        write(0,*) (min(0.,psmlt_e)-psacr_e)
                        write(0,*) 'qx,cx = ',an(i,1,k,ls),an(i,1,k,ln(ls))
                        write(0,*) 'stopping'
                        print*, 'stopping'
                        STOP
                    end if
                    if ((abs(pscswmd)*dt) .gt. abs(an(i,1,k,lsc(ls)))) then
                        tmp = abs(an(i,1,k,lsc(ls)) / (pscswmd*dt))
                        write(0,*) 'conserving snow scs, tmp= ',tmp
                        write(0,*) 'temper = ',T(i,k)
                        write(0,*) 'i,j,k,t=',i,j,k,itimestep
                        write(0,*) 'pscswmd, dt',pscswmd, dt
                        write(0,*) 'pscswmd*dt, an(i,1,k,lsc(ls))', pscswmd*dt, an(i,1,k,lsc(ls))
                        psmlt_e = tmp * psmlt_e
                        psacr_e = tmp * psacr_e
                        pscswmd = tmp * pscswmd
                        write(0,*) 'pscswmd after conserving',pscswmd, pscswmd*dt
                    end if
                end if

                !
                !  graupel
                !
                if ( an(i,1,k,lh) .gt. qxmin(lh) ) then
!                    pschwmi = fscrw*(max(0.,-pracg)+pgracs+piacr+pracg1) +&
!                            fsccw*(psacwg+pgsacw) +&
!                            fscsw*(psacr) +&
!                            fscci*(praci)
                    pgmlt_e = pgmlt
                    pschwmd = fschw*(min(0.,pgmlt_e))
                    if ((min(0.,pgmlt_e)).gt.0.0) then
                        write(0,*) 'sink terms pos for graupel'
                        write(0,*) 'temper = ',T(i,k)
                        write(0,*) 'i,j,k,t=',i,j,k,itimestep
                        write(0,*) '(min(0.,pgmlt_e))', (min(0.,pgmlt_e))
                    end if
                    if (pschwmd*an(i,1,k,lsc(lh)).gt.0.) then ! pscxxmd should lead to decrease of scx, so pscxxmd*scx should be neg
                        write(0,*) 'graupel pschwmd same sign as an(lsch)'
                        write(0,*) 'temper = ',T(i,k)
                        write(0,*) 'i,j,k,t=',i,j,k,itimestep
                        write(0,*) 'pschwmd, an(i,1,k,lsc(lh))',pschwmd, an(i,1,k,lsc(lh))
                        write(0,*) 'fschw,pgmlt'
                        write(0,*) fschw,pgmlt
                        write(0,*) '(min(0.,pgmlt_e))'
                        write(0,*) (min(0.,pgmlt_e))
                        write(0,*) 'qx,cx = ',an(i,1,k,lh),an(i,1,k,ln(lh))
                        write(0,*) 'stopping'
                        print*, 'stopping'
                        STOP
                    end if
                    if ((abs(pschwmd)*dt) .gt. abs(an(i,1,k,lsc(lh)))) then
                        tmp = abs(an(i,1,k,lsc(lh)) / (pschwmd*dt))
                        write(0,*) 'conserving graupel sch, tmp= ',tmp
                        write(0,*) 'temper = ',T(i,k)
                        write(0,*) 'i,j,k,t=',i,j,k,itimestep
                        write(0,*) 'pschwmd, dt',pschwmd, dt
                        write(0,*) 'pschwmd*dt, an(i,1,k,lsc(lh))', pschwmd*dt, an(i,1,k,lsc(lh))
                        pgmlt_e = tmp * pgmlt_e
                        pschwmd = tmp * pschwmd
                        write(0,*) 'pschwmd after conserving',pschwmd, pschwmd*dt
                    end if
                end if

                !
                !  rain
                !
                if ( an(i,1,k,lr) .gt. qxmin(lr) ) then
!                    pscrwmi = fsccw*(prc+pra) +&
!                            fscsw*(max(0.,-psmlt)) +&
!                            fschw*(max(0.,-pgmlt))
                    pracs1_e = pracs1
                    qmultr_e = qmultr
                    qmultrg_e = qmultrg
                    piacr_e = piacr
                    piacrs_e = piacrs
                    pracg1_e = pracg1
                    pgracs_e = pgracs
                    pracs_e = pracs
                    pracg_e = pracg
                    pscrwmd = fscrw*(-pracs1-qmultr-qmultrg-piacr-piacrs &
                            -pracg1-pgracs &
                                     +min(0.,pracs)+min(0.,pracg))
                    if (-pracs1.gt.0.0 .or.&
                            -qmultr.gt.0.0 .or.&
                            -qmultrg.gt.0.0 .or.&
                            -piacr.gt.0.0 .or.&
                            -piacrs.gt.0.0 .or.&
                            -pracg1.gt.0.0 .or.&
                            -pgracs.gt.0.0 .or.&
                            min(0.,pracs).ne.0.0 .or.&
                            min(0.,pracg).ne.0.0) then
                        write(0,*) 'sink terms pos for rain'
                        write(0,*) 'temper = ',T(i,k)
                        write(0,*) 'i,j,k,t=',i,j,k,itimestep
                        write(0,*) '-pracs1', -pracs1
                        write(0,*) '-qmultr', -qmultr
                        write(0,*) '-qmultrg', -qmultrg
                        write(0,*) '-piacr', -piacr
                        write(0,*) '-piacrs', -piacrs
                        write(0,*) '-pracg1', -pracg1
                        write(0,*) '-pgracs', -pgracs
                        write(0,*) 'min(0.,pracs)', min(0.,pracs)
                        write(0,*) 'min(0.,pracg)', min(0.,pracg)
                    end if
                    if (pscrwmd*an(i,1,k,lsc(lr)).gt.0.) then ! pscxxmd should lead to decrease of scx, so pscxxmd*scx should be neg
                        write(0,*) 'rain pscrwmd same sign as an(lscr)'
                        write(0,*) 'temper = ',T(i,k)
                        write(0,*) 'i,j,k,t=',i,j,k,itimestep
                        write(0,*) 'pscrwmd, an(i,1,k,lsc(lr))',pscrwmd, an(i,1,k,lsc(lr))
                        write(0,*) 'fschw,pracs1,qmultr,qmultrg,piacr,piacrs,pracg1,pgracs,pracs,pracg'
                        write(0,*) fschw,pracs1,qmultr,qmultrg,piacr,piacrs,pracg1,pgracs,pracs,pracg
                        write(0,*) '(-pracs1-qmultr-qmultrg-piacr-piacrs -pracg1-pgracs +min(0.,pracs)+min(0.,pracg))'
                        write(0,*) (-pracs1-qmultr-qmultrg-piacr-piacrs -pracg1-pgracs +min(0.,pracs)+min(0.,pracg))
                        write(0,*) 'qx,cx = ',an(i,1,k,lr),an(i,1,k,ln(lr))
                        write(0,*) 'stopping'
                        print*, 'stopping'
                        STOP
                    end if
                    if ((abs(pscrwmd)*dt) .gt. abs(an(i,1,k,lsc(lr)))) then
                        tmp = abs(an(i,1,k,lsc(lr)) / (pscrwmd*dt))
                        write(0,*) 'conserving rain scr, tmp= ',tmp
                        write(0,*) 'temper = ',T(i,k)
                        write(0,*) 'i,j,k,t=',i,j,k,itimestep
                        write(0,*) 'pscrwmd, dt',pscrwmd, dt
                        write(0,*) 'pscrwmd*dt, an(i,1,k,lsc(lr))', pscrwmd*dt, an(i,1,k,lsc(lr))
                        pracs1_e = tmp * pracs1_e
                        qmultr_e = tmp * qmultr_e
                        qmultrg_e = tmp * qmultrg_e
                        piacr_e = tmp * piacr_e
                        piacrs_e = tmp * piacrs_e
                        pracg1_e = tmp * pracg1_e
                        pgracs_e = tmp * pgracs_e
                        pracs_e = tmp * pracs_e
                        pracg_e = tmp * pracg_e
                        pscrwmd = tmp * pscrwmd
                        write(0,*) 'pscrwmd after conserving',pscrwmd, pscrwmd*dt
                    end if
                end if

                !-------------- Now get source rates from conserved sink rates
                !
                !  cloud ice
                !
                if ( an(i,1,k,li) .gt. qxmin(li) ) then
                    psccimi = fsccw*(qmults_e+qmultg_e) +&
                            fscrw*(qmultr_e+qmultrg_e)
                end if

                !
                !  cloud water
                !
                if ( an(i,1,k,lc) .gt. qxmin(lc) ) then
                    psccwmi = 0.0
                end if

                !
                !  snow
                !
                if ( an(i,1,k,ls) .gt. qxmin(ls) ) then
                    pscswmi = fsccw*(psacws_e) +&
                            fscrw*(pracs1_e+piacrs_e+max(0.,-pracs_e)) +&
                            fscci*(pracis_e)
                end if

                !
                !  graupel
                !
                if ( an(i,1,k,lh) .gt. qxmin(lh) ) then
                    pschwmi = fscrw*(max(0.,-pracg_e)+pgracs_e+piacr_e+pracg1_e) +&
                            fsccw*(psacwg_e+pgsacw_e) +&
                            fscsw*(psacr_e) +&
                            fscci*(praci_e)
                end if

                !
                !  rain
                !
                if ( an(i,1,k,lr) .gt. qxmin(lr) ) then
                    pscrwmi = fsccw*(prc_e+pra_e) +&
                            fscsw*(max(0.,-psmlt_e)) +&
                            fschw*(max(0.,-pgmlt_e))
                end if

                !!
                !!
                !! Collisional Charging transfers (inductive/noninductive)
                !!
                !!
                !
                !  cloud water
                !
                psccwi =    max(-schacw(i,k),0.0)   &
                        & + max(-scsacw(i,k),0.0)
                psccwd =    min(-schacw(i,k),0.0)   &
                        & + min(-scsacw(i,k),0.0)

                qcwtmp = Max(0.0, an(i,1,k,lc))
                IF ( qcwtmp .gt. qxmin(lc) ) then
                    ccwtmp = Max(0.001, an(i,1,k,lnc))
                    cpqc = (scx(i,k,lc) + dtp*(psccwi+psccwd + psccwmi+psccwmd)) &
                            /ccwtmp
                    IF ( Abs(cpqc) .gt. scwppmx .and.   &
                            & (psccwi+psccwd) .ne. 0.0 .and.   &
                            & Sign(1.0, scx(i,k,lc)) .eq.   &
                                    & Sign(1.0, (psccwi+psccwd))) THEN
                        cpqc0 = Sign(scwppmx,cpqc)
                        scfac = (cpqc0*ccwtmp - scx(i,k,lc) -   &
                                &       dtp*(psccwmi+psccwmd))/   &
                                &              ( dtp*(psccwi+psccwd) )
                        IF ( scfac .gt. 1.0 ) THEN
                            write(0,*) 'OUCH: scwfac = ',scfac
                            write(0,*) 'OUCH: scwfac = ',scfac
                        ENDIF
                        scfac = Max(scfac,0.0)
                        scfac = Min(scfac,1.0)
                    ELSE
                        scfac = 1.0
                    ENDIF
                ELSE
                    scfac = 0.0
                ENDIF ! ( qcwtmp .gt. qxmin(lc) )

                IF ( scfac .lt. 1.0 .and. an(i,1,k,lc) .gt. qxmin(lc) ) THEN
                    schacw(i,k) = scfac*schacw(i,k)
                    scsacw(i,k) = scfac*scsacw(i,k)
                    psccwi =  max(-schacw(i,k),0.0)   &
                          & + max(-scsacw(i,k),0.0)
                    psccwd =  min(-schacw(i,k),0.0)   &
                          & + min(-scsacw(i,k),0.0)
                ENDIF

                !
                !  cloud ice
                !
                psccii =    max(-schaci(i,k),0.0)   &
                        & + max(-scsaci(i,k),0.0)
                psccid =    min(-schaci(i,k),0.0)   &
                        & + min(-scsaci(i,k),0.0)

                qcitmp = Max(0.0, an(i,1,k,li))
                IF ( qcitmp .gt. qxmin(li) ) then
                    ccitmp = Max(0.001, an(i,1,k,lni) )
                    cpci = (scx(i,k,li) + dtp*(psccii+psccid + psccimi+psccimd)) &
                            /ccitmp
                    IF ( Abs(cpci) .gt. scippmx .and.   &
                            & (psccii+psccid ) .ne. 0.0 .and.   &
                            & Sign(1.0,scx(i,k,li)) .eq.   &
                                    & Sign(1.0,(psccii+psccid))) THEN
                        cpci0 = Sign(scippmx,cpci)
                        scfac = (cpci0*ccitmp - scx(i,k,li) -    &
                                &       dtp*(psccimi+psccimd))/   &
                                &              ( dtp*(psccii+psccid) )
                        IF ( scfac .gt. 1.0 ) THEN
                            write(0,*) 'OUCH: sccifac = ',scfac,   &
                                    & dtp*(psccimi+psccimd),cpci,scx(i,k,li),ccitmp,   &
                                    & i,j,k,t(i,k)
                        ENDIF
                        scfac = Max(scfac,0.0)
                        scfac = Min(scfac,1.0)
                    ELSE
                        scfac = 1.0
                    ENDIF
                ELSE
                    scfac = 0.0
                ENDIF ! ( qcitmp .gt. qxmin(li) )

                IF ( scfac .lt. 1.0 ) THEN
                    schaci(i,k) = scfac*schaci(i,k)
                    scsaci(i,k) = scfac*scsaci(i,k)
                    psccii =    max(-schaci(i,k),0.0)   &
                            & + max(-scsaci(i,k),0.0)
                    psccid =    min(-schaci(i,k),0.0)   &
                            & + min(-scsaci(i,k),0.0)
                ENDIF

                !
                !  rain
                !
!                pscrwi = max(-scsacr(i,k),0.0)
!                pscrwd = min(-scsacr(i,k),0.0)

                !
                !  snow
                !
                pscswi =    max(scsaci(i,k),0.0)   &
                        & + max(scsacw(i,k),0.0)   &
                        & + max(-schacs(i,k),0.0)
                pscswd =    min(scsaci(i,k),0.0)   &
                        & + min(scsacw(i,k),0.0)   &
                        & + min(-schacs(i,k),0.0)

                !
                !  graupel
                !
                pschwi =    max(schaci(i,k),0.0)   &
                        & + max(schacs(i,k),0.0)   &
                        & + max(schacw(i,k),0.0)
                pschwd =    min(schaci(i,k),0.0)   &
                        & + min(schacs(i,k),0.0)   &
                        & + min(schacw(i,k),0.0)

                do iltest=lc,lh
                    if (isnan(scx(i,k,iltest))) then
                        write (0,*) 'scx=nan, iltest=',iltest
                        write (0,*) 'i,j,k,t=',i,j,k,itimestep
                        write (0,*) 'scx(k,iltest)=',scx(i,k,iltest)
                        write (0,*) 'an(i,1,k,lsc(il))=',an(i,1,k,lsc(il))
                    end if
                end do

                ! update scx array
                scx(i,k,lc) = scx(i,k,lc) +   &
                        &   dt*(psccwi + psccwd + psccwmi + psccwmd)
                scx(i,k,li) = scx(i,k,li) +   &
                        &   dt*(psccii + psccid + psccimi + psccimd)
                scx(i,k,lr) = scx(i,k,lr) +   &
                        &   dt*(pscrwi + pscrwd + pscrwmi + pscrwmd)
                scx(i,k,ls) = scx(i,k,ls) +   &
                        &   dt*(pscswi + pscswd + pscswmi + pscswmd)
                scx(i,k,lh) = scx(i,k,lh) +   &
                        &   dt*(pschwi + pschwd + pschwmi + pschwmd)

                cionp(i,k) = cionp(i,k) +  &
                        &   dt*ec_i*(pscpii + pscpid) ! #/m3 * C / C/m3. ! each ion has electron charge
                cionn(i,k) = cionn(i,k) -  &
                        &   dt*ec_i*(pscnii + pscnid)
                if ( largeion ) then
                    clionp(i,k) = clionp(i,k) +  &
                            &   dt*ec_i*(pscplii + pscplid)
                    clionn(i,k) = clionn(i,k) -  &
                            &   dt*ec_i*(pscnlii + pscnlid)
                endif

                do iltest=lc,lh
                    if (isnan(scx(i,k,iltest))) then
                        write (0,*) 'scx=nan, iltest=',iltest
                        write (0,*) 'i,j,k,t=',i,j,k,itimestep
                        write (0,*) 'scx(k,iltest)=',scx(i,k,iltest)
                        write (0,*) 'psccwi psccwd psccwmi psccwmd '
                        write (0,*) psccwi, psccwd, psccwmi, psccwmd
                        write (0,*) 'psccii psccid psccimi psccimd '
                        write (0,*) psccii, psccid, psccimi, psccimd
                        write (0,*) 'pscrwi pscrwd pscrwmi pscrwmd '
                        write (0,*) pscrwi, pscrwd, pscrwmi, pscrwmd
                        write (0,*) 'pscswi pscswd pscswmi pscswmd '
                        write (0,*) pscswi, pscswd, pscswmi, pscswmd
                        write (0,*) 'pschwi pschwd pschwmi pschwmd '
                        write (0,*) pschwi, pschwd, pschwmi, pschwmd
                        write (0,*) 'stopping'
                        print*, 'stopping'
                        STOP
                    end if
                end do

                ! put scx arrays back into an array
                an(i,1,k,lscpi) = cionp(i,k)
                an(i,1,k,lscni) = cionn(i,k)

                DO il = lc,lhab
                    tmp = an(i,1,k,lsc(il))
                    an(i,1,k,lsc(il)) = scx(i,k,il)
                    IF (.not.( Abs (scx(i,k,il) ) < 1000.e-9  )) THEN
                        write(0,*) 'Problem0a with scx il = ',il,lsc(il),scx(i,k,il)
                        write(0,*) 'Problem2 with scx il = ',tmp,tmp-scx(i,k,il),(tmp-scx(i,k,il))*dtpinv
                        write(0,*) 'temper = ',T(i,k)
                        write(0,*) 'i,j,k,t=',i,j,k,itimestep
                        write(0,*) 'qx,cx = ',an(i,1,k,il),an(i,1,k,ln(il))
                        write(0,*)  psccwi, psccwd
                        write(0,*)  psccii, psccid
                        write(0,*)  pscrwi, pscrwd
                        write(0,*)  pscswi, pscswd
                        write(0,*)  pschwi, pschwd

                        write(0,*)  pscpii ,pscpid
                        write(0,*)  pscnii ,pscnid
                        write(0,*)
                        write(0,*)  psccwmi, psccwmd
                        write(0,*)  psccimi, psccimd
                        write(0,*)  pscrwmi, pscrwmd
                        write(0,*)  pscswmi, pscswmd
                        write(0,*)  pschwmi, pschwmd

                        write(0,*) 'fscrw',fscrw,max(0.,-pracg),pgracs,piacr,pracg1
                        write(0,*) 'fsccw',fsccw,psacwg,pgsacw
                        write(0,*) 'fscsw',fscsw,psacr
                        write(0,*) 'fscci',fscci,praci
                        write(0,*) fscrw*(max(0.,-pracg)+pgracs+piacr+pracg1)
                        write(0,*) fsccw*(psacwg+pgsacw)
                        write(0,*) fscsw*(psacr), fscci*(praci)
                        write(0,*)
                    ENDIF
                ENDDO

            END IF ! IPELEC
            ! end charging
            ! ⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡


            rainevap(i,k) = -1.*pre
            nraggout(i,k) = nragg
            IF(snowflag.eq.1)THEN
               snowevap(i,k) = -1.*evpms
               snowmelt(i,k) = -1.*psmlt
               snowdep(i,k) = prds
               snowsub(i,k) = eprds
               snowaccr(i,k) = pracs
            ELSE IF(snowflag.eq.2)THEN
               snowmelt(i,k) = qsmlt*-1.
            END IF
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
            !relh(i,k) = qv(i,k)/qvs*100.
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

            IF(meltflag.ne.1)THEN
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
                      if (ipelec>0) then
                          an(i,1,k,lsc(lr)) = an(i,1,k,lsc(lr)) + an(i,1,k,lsc(li))
                          an(i,1,k,lsc(li)) = 0.0
                      end if
                   END IF
                   IF(snowflag.eq.2)THEN
                      IF(qs(i,k).ge.qsmall)THEN
                         qr(i,k)=qr(i,k)+qs(i,k)
                         nr(i,k)=nr(i,k)+ns(i,k)
                         t(i,k)=t(i,k)-qs(i,k)*xxlf/cpm
                         snowmelt(i,k) = qs(i,k)/dt
                         qs(i,k)=0.
                         ns(i,k)=0.
                         as(i,k)=0.
                         cs(i,k)=0.
                         if (ipelec>0) then
                             an(i,1,k,lsc(lr)) = an(i,1,k,lsc(lr)) + an(i,1,k,lsc(ls))
                             an(i,1,k,lsc(ls)) = 0.0
                         end if
                      END IF
                   END IF
                END IF
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
                  if (ipelec>0) then
                      an(i,1,k,lsc(li)) = an(i,1,k,lsc(li)) + an(i,1,k,lsc(lr))
                      an(i,1,k,lsc(lr)) = 0.0
                  end if

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
                  if (ipelec>0) then
                      an(i,1,k,lsc(li)) = an(i,1,k,lsc(li)) + an(i,1,k,lsc(lc))
                      an(i,1,k,lsc(lc)) = 0.0
                  end if

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

            IF(qs(i,k) .gt. qsmall)THEN
               ans = as(i,k)/(ns(i,k)*nus)
               cns = cs(i,k)/(ns(i,k)*nus)
               as(i,k) = ans**2*cns*nus*ns(i,k)
               cs(i,k) = cns**2*ans*nus*ns(i,k)
            END IF

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
            IF(snowflag.eq.1)THEN
               IF(qs(i,k).ge.qsmall)THEN
                  dum = (cons1*ns(i,k)/qs(i,k))**(1./DS)
                  fs(k) = asn(k)*cons3/dum**BS
                  fns(k) = asn(k)*cons5/dum**BS
                  fs(k) = min(fs(k),1.2*(rhosu/rho(i,k))**0.54)
                  fns(k) = min(fns(k),1.2*(rhosu/rho(i,k))**0.54)
               ELSE
                  fs(k) = 0.
                  fns(k) = 0.
               END IF
            ELSE IF(snowflag.eq.2)THEN
                IF(qs(i,k).ge.qsmall)THEN
                   fs(k) = vtrms(i,k)
                   fns(k) = vtrns(i,k)
                   fas(k) = vtrms(i,k)
                   fcs(k) = vtrms(i,k)
                ELSE
                   fs(k) = 0.0
                   fns(k) = 0.0
                   fas(k) = 0.0
                   fcs(k) = 0.0
                END IF
            END IF !snowflag

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
               IF(fas(k).lt.1.e-10) fas(k) = fas(k+1)
               IF(fcs(k).lt.1.e-10) fcs(k) = fcs(k+1)
               IF(fr(k).lt.1.e-10) fr(k) = fr(k+1)
               IF(fnr(k).lt.1.e-10) fnr(k) = fnr(k+1)
               IF(fc(k).lt.1.e-10) fc(k) = fc(k+1)
               IF(fnc(k).lt.1.e-10) fnc(k) = fnc(k+1)
               IF(fg(k).lt.1.e-10) fg(k) = fg(k+1)
               IF(fng(k).lt.1.e-10) fng(k) = fng(k+1)
            END IF
!     calculate number of split time steps

            rgvm = max(fi(k),fni(k),fai(k),fci(k),fr(k),fnr(k),fc(k),&
               fnc(k),fs(k),fns(k),fas(k),fcs(k),fg(k),fng(k))
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
            dumas(k) = as(i,k)*rho(i,k) !unitless (as is volume now)
            dumcs(k) = cs(i,k)*rho(i,k) !unitless (cs is volume now)
            dumg(k) = qg(i,k)*rho(i,k)
            dumng(k) = ng(i,k)*rho(i,k)
            !-----------------------------------------------------------------------
            !Explicit lightning ----------------------------------------------------
            if (ipelec>0) then
                dumsc(k,lsc(lc)) = an(i,1,k,lsc(lc))
                dumsc(k,lsc(lr)) = an(i,1,k,lsc(lr))
                dumsc(k,lsc(li)) = an(i,1,k,lsc(li))
                dumsc(k,lsc(ls)) = an(i,1,k,lsc(ls))
                dumsc(k,lsc(lg)) = an(i,1,k,lsc(lg))
            end if
            !End explicit lightning ------------------------------------------------
            !-----------------------------------------------------------------------
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
               falloutas(k) = fas(k)*dumas(k) !m s^-1
               falloutcs(k) = fcs(k)*dumcs(k) !m s^-1
               falloutg(k) = fg(k)*dumg(k)
               falloutng(k) = fng(k)*dumng(k)
               !-----------------------------------------------------------------------
               !Explicit lightning ----------------------------------------------------
               if (ipelec>0) then
                   faloutscc(k) = fc(k)*dumsc(k,lsc(lc))
                   faloutscr(k) = fr(k)*dumsc(k,lsc(lr))
                   faloutsci(k) = fi(k)*dumsc(k,lsc(li))
                   faloutscs(k) = fs(k)*dumsc(k,lsc(ls))
                   faloutscg(k) = fg(k)*dumsc(k,lsc(lh))
               end if
               !End explicit lightning ------------------------------------------------
               !-----------------------------------------------------------------------

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
            falltndas = falloutas(k)/dzq(k) !s^-1
            falltndcs = falloutcs(k)/dzq(k) !s^-1
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
            assten(k) = assten(k) - falltndas/nstep/rho(i,k) !m^3 kg^-1 s^-1
            cssten(k) = cssten(k) - falltndcs/nstep/rho(i,k) !m^3 kg^-1 s^-1
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
            dumas(k) = dumas(k) - falltndas*dt/nstep
            dumcs(k) = dumcs(k) - falltndcs*dt/nstep
            dumg(k) = dumg(k) - falltndg*dt/nstep
            dumng(k) = dumng(k) - falltndng*dt/nstep            
            
            !-----------------------------------------------------------------------
            !Explicit lightning ----------------------------------------------------
            if (ipelec>0) then
                faltndscc = faloutscc(k)/dzq(k)
                faltndscr = faloutscr(k)/dzq(k)
                faltndsci = faloutsci(k)/dzq(k)
                faltndscs = faloutscs(k)/dzq(k)
                faltndscg = faloutscg(k)/dzq(k)

                pscc_sed(k) = pscc_sed(k)-faltndscc/nstep
                pscr_sed(k) = pscr_sed(k)-faltndscr/nstep
                psci_sed(k) = psci_sed(k)-faltndsci/nstep
                pscs_sed(k) = pscs_sed(k)-faltndscs/nstep
                pscg_sed(k) = pscg_sed(k)-faltndscg/nstep

                dumsc(k,lsc(lc)) = dumsc(k,lsc(lc)) - faltndscc*dt/nstep
                dumsc(k,lsc(lr)) = dumsc(k,lsc(lr)) - faltndscr*dt/nstep
                dumsc(k,lsc(li)) = dumsc(k,lsc(li)) - faltndsci*dt/nstep
                dumsc(k,lsc(ls)) = dumsc(k,lsc(ls)) - faltndscs*dt/nstep
                dumsc(k,lsc(lg)) = dumsc(k,lsc(lg)) - faltndscg*dt/nstep
            end if
            !End explicit lightning ------------------------------------------------
            !-----------------------------------------------------------------------

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
               falltndas = (falloutas(k+1) - falloutas(k))/dzq(k)
               falltndcs = (falloutcs(k+1) - falloutcs(k))/dzq(k)
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
               assten(k) = assten(k) + falltndas/nstep/rho(i,k)
               cssten(k) = cssten(k) + falltndcs/nstep/rho(i,k)
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
               dumas(k) = dumas(k) + falltndas*dt/nstep
               dumcs(k) = dumcs(k) + falltndcs*dt/nstep
               dumg(k) = dumg(k) + falltndg*dt/nstep
               dumng(k) = dumng(k) + falltndng*dt/nstep
  

               qiloss(k)=qiloss(k)+fallouti(k)
               qrloss(k)=qrloss(k)+falloutr(k)+falloutc(k)

               !-----------------------------------------------------------------------
               !Explicit lightning ----------------------------------------------------
               if (ipelec>0) then
                   FALTNDSCC = (FALOUTSCC(K+1)-FALOUTSCC(K))/DZQ(K)
                   FALTNDSCR = (FALOUTSCR(K+1)-FALOUTSCR(K))/DZQ(K)
                   FALTNDSCI = (FALOUTSCI(K+1)-FALOUTSCI(K))/DZQ(K)
                   FALTNDSCS = (FALOUTSCS(K+1)-FALOUTSCS(K))/DZQ(K)
                   FALTNDSCG = (FALOUTSCG(K+1)-FALOUTSCG(K))/DZQ(K)

                   PSCC_SED(K) = PSCC_SED(K)+FALTNDSCC/NSTEP
                   PSCR_SED(K) = PSCR_SED(K)+FALTNDSCR/NSTEP
                   PSCI_SED(K) = PSCI_SED(K)+FALTNDSCI/NSTEP
                   PSCS_SED(K) = PSCS_SED(K)+FALTNDSCS/NSTEP
                   PSCG_SED(K) = PSCG_SED(K)+FALTNDSCG/NSTEP

                   DUMSC(K,lsc(LC)) = DUMSC(K,lsc(LC)) + FALTNDSCC*DT/NSTEP
                   DUMSC(K,lsc(LR)) = DUMSC(K,lsc(LR)) + FALTNDSCR*DT/NSTEP
                   DUMSC(K,lsc(LI)) = DUMSC(K,lsc(LI)) + FALTNDSCI*DT/NSTEP
                   DUMSC(K,lsc(LS)) = DUMSC(K,lsc(LS)) + FALTNDSCS*DT/NSTEP
                   DUMSC(K,lsc(LG)) = DUMSC(K,lsc(LG)) + FALTNDSCG*DT/NSTEP
               end if
               !End explicit lightning ------------------------------------------------
               !-----------------------------------------------------------------------

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
            as(i,k) = as(i,k) + assten(k)*dt !m kg^-1
            cs(i,k) = cs(i,k) + cssten(k)*dt !m kg^-1
            qg(i,k) = qg(i,k) + qgsten(k)*dt
            ng(i,k) = ng(i,k) + ngsten(k)*dt
            !-----------------------------------------------------------------------
            !Explicit lightning ----------------------------------------------------
            if (ipelec>0) then
                tmpsc(k,lsc(lc)) = an(i,1,k,lsc(lc))
                tmpsc(k,lsc(lr)) = an(i,1,k,lsc(lr))
                tmpsc(k,lsc(li)) = an(i,1,k,lsc(li))
                tmpsc(k,lsc(ls)) = an(i,1,k,lsc(ls))
                tmpsc(k,lsc(lg)) = an(i,1,k,lsc(lh))
                an(i,1,k,lsc(lc)) = an(i,1,k,lsc(lc)) + pscc_sed(k)*dt
                an(i,1,k,lsc(lr)) = an(i,1,k,lsc(lr)) + pscr_sed(k)*dt
                an(i,1,k,lsc(li)) = an(i,1,k,lsc(li)) + psci_sed(k)*dt
                an(i,1,k,lsc(ls)) = an(i,1,k,lsc(ls)) + pscs_sed(k)*dt
                an(i,1,k,lsc(lh)) = an(i,1,k,lsc(lh)) + pscg_sed(k)*dt
                do il = lc,lhab
                    if (.not.( abs (an(i,1,k,lsc(il)) ) < 1000.e-9  )) then
                        write(0,*) 'problem0a with scx after sedi.'
                        write(0,*) 'il, lsc(il), = ',il,lsc(il),an(i,1,k,lsc(il))
                        write(0,*) 'new, old = ',an(i,1,k,lsc(il)), tmpsc(k,lsc(il))
                        write(0,*) 'pscc_sed(k), pscc_sed(k)*dt',pscc_sed(k), pscc_sed(k)*dt
                        write(0,*) 'temper = ',t(i,k)
                        write(0,*) 'i,j,k,t=',i,j,k,itimestep
                        write(0,*) 'qx,cx = ',an(i,1,k,il),an(i,1,k,ln(il))
                        write(0,*)  pscc_sed(k)
                        write(0,*)  psci_sed(k)
                        write(0,*)  pscr_sed(k)
                        write(0,*)  pscs_sed(k)
                        write(0,*)  pscg_sed(k)
                        if (abs(an(i,1,k,lsc(il))).gt.100000.e-9) then
                            write (0,*) 'Exceeding 1e-4'
                            write (0,*) 'stopping'
                            print*, 'stopping'
                            stop
                        end if
                    end if ! >1e-6
                enddo ! il
            end if ! ipelec            !End explicit lightning ------------------------------------------------
            !-----------------------------------------------------------------------


            ! yc added 2/24/2022
!            if(qi(i,k).lt.0.0) write(0,*) 'qi neg',qi(i,k),i,k
!            if(ni(i,k).lt.0.0) write(0,*) 'ni neg',ni(i,k),i,k
!            if(ai(i,k).lt.0.0) write(0,*) 'ai neg',ai(i,k),i,k
!            if(ci(i,k).lt.0.0) write(0,*) 'ci neg',ci(i,k),i,k
!            if(qr(i,k).lt.0.0) write(0,*) 'qr neg',qr(i,k),i,k
!            if(qc(i,k).lt.0.0) write(0,*) 'qc neg',qc(i,k),i,k
!            if(nr(i,k).lt.0.0) write(0,*) 'nr neg',nr(i,k),i,k
!            if(nc(i,k).lt.0.0) write(0,*) 'nc neg',nc(i,k),i,k
!            if(qs(i,k).lt.0.0) write(0,*) 'qs neg',qs(i,k),i,k
!            if(ns(i,k).lt.0.0) write(0,*) 'ns neg',ns(i,k),i,k
!            if(as(i,k).lt.0.0) write(0,*) 'as neg',as(i,k),i,k
!            if(cs(i,k).lt.0.0) write(0,*) 'cs neg',cs(i,k),i,k
!            if(qg(i,k).lt.0.0) write(0,*) 'qg neg',qg(i,k),i,k
!            if(ng(i,k).lt.0.0) write(0,*) 'ng neg',ng(i,k),i,k
            ! get rid of artificial neg q/n, prevent adding to qv
!            qi(i,k) = max(qi(i,k),0.0)
!            ni(i,k) = max(ni(i,k),0.0)
!            ai(i,k) = max(ai(i,k),0.0)
!            ci(i,k) = max(ci(i,k),0.0)
!            qr(i,k) = max(qr(i,k),0.0)
!            qc(i,k) = max(qc(i,k),0.0)
!            nr(i,k) = max(nr(i,k),0.0)
!            nc(i,k) = max(nc(i,k),0.0)
!            qs(i,k) = max(qs(i,k),0.0)
!            ns(i,k) = max(ns(i,k),0.0)
!            as(i,k) = max(as(i,k),0.0)
!            cs(i,k) = max(cs(i,k),0.0)
!            qg(i,k) = max(qg(i,k),0.0)
!            ng(i,k) = max(ng(i,k),0.0)
            !if(as(i,k).gt.0.0)print*,'as4',real(itimestep)*dt,i,k,as(i,k),cs(i,k),((cs(i,k)**2)/(as(i,k)*nus*ns(i,k)))**(1./3.)

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
                  if (ipelec>0) then
                      IF ( an(i,1,k,lscr) .gt. 0.0 ) THEN
                          an(i,1,k,lscpi) = an(i,1,k,lscpi) + an(i,1,k,lscr)*ec_i
                      ELSE
                          an(i,1,k,lscni) = an(i,1,k,lscni) - an(i,1,k,lscr)*ec_i
                      ENDIF
                      an(i,1,k,lscr) = 0.0
                  end if
               END IF
               IF(qc(i,k).lt.1.e-8)THEN
                  qv(i,k) = qv(i,k) + qc(i,k)
                  t(i,k) = t(i,k) - qc(i,k)*xxlv/cpm
                  qc(i,k) = 0.0
                  if (ipelec>0) then
                      IF ( an(i,1,k,lscw) .gt. 0.0 ) THEN
                          an(i,1,k,lscpi) = an(i,1,k,lscpi) + an(i,1,k,lscw)*ec_i
                      ELSE
                          an(i,1,k,lscni) = an(i,1,k,lscni) - an(i,1,k,lscw)*ec_i
                      ENDIF
                      an(i,1,k,lscw) = 0.0
                  end if
               END IF
            END IF

            IF(qvqvsi.lt.0.9)THEN
               IF(qi(i,k).lt.1.e-8)THEN
                  qv(i,k) = qv(i,k) + qi(i,k)
                  t(i,k) = t(i,k) - qi(i,k)*xxls/cpm
                  qi(i,k) = 0.0
                  if (ipelec>0) then
                      IF ( an(i,1,k,lsc(li)) .gt. 0.0 ) THEN
                          an(i,1,k,lscpi) = an(i,1,k,lscpi) + an(i,1,k,lsc(li))*ec_i
                      ELSE
                          an(i,1,k,lscni) = an(i,1,k,lscni) - an(i,1,k,lsc(li))*ec_i
                      ENDIF
                      an(i,1,k,lsc(li)) = 0.0
                  end if
               END IF
               IF(qs(i,k).lt.1.e-8)THEN
                  qv(i,k) = qv(i,k) + qs(i,k)
                  t(i,k) = t(i,k) - qs(i,k)*xxls/cpm
                  qs(i,k) = 0.0
                  if (ipelec>0) then
                      IF ( an(i,1,k,lsc(ls)) .gt. 0.0 ) THEN
                          an(i,1,k,lscpi) = an(i,1,k,lscpi) + an(i,1,k,lsc(ls))*ec_i
                      ELSE
                          an(i,1,k,lscni) = an(i,1,k,lscni) - an(i,1,k,lsc(ls))*ec_i
                      ENDIF
                      an(i,1,k,lsc(ls)) = 0.0
                  end if
               END IF
               IF(qg(i,k).lt.1.e-8)THEN
                  qv(i,k) = qv(i,k) + qg(i,k)
                  t(i,k) = t(i,k) - qg(i,k)*xxls/cpm
                  qg(i,k) = 0.0
                  if (ipelec>0) then
                      IF ( an(i,1,k,lsc(lh)) .gt. 0.0 ) THEN
                          an(i,1,k,lscpi) = an(i,1,k,lscpi) + an(i,1,k,lsc(lh))*ec_i
                      ELSE
                          an(i,1,k,lscni) = an(i,1,k,lscni) - an(i,1,k,lsc(lh))*ec_i
                      ENDIF
                      an(i,1,k,lsc(lh)) = 0.0
                  end if
               END IF
            END IF
            IF(qi(i,k) .lt. qsmall .or. ni(i,k) .lt. qsmall)THEN
               qi(i,k) = 0.
               ni(i,k) = 0.
               ai(i,k) = 0.
               ci(i,k) = 0.
               qiloss(k) = 0.
               if (ipelec>0) then
                   an(i,1,k,lsc(li)) = 0.0
               end if
            END IF
            IF(qr(i,k) .lt. qsmall .or. nr(i,k) .lt. qsmall)THEN
               qr(i,k) = 0.
               nr(i,k) = 0.
               qrloss(k) = 0.
               if (ipelec>0) then
                   an(i,1,k,lsc(lr)) = 0.0
               end if
            END IF
            IF(qc(i,k) .lt. qsmall .or. nc(i,k) .lt. qsmall)THEN
               qc(i,k) = 0.
               nc(i,k) = 0.
               qcloss(k) = 0.
               if (ipelec>0) then
                   an(i,1,k,lsc(lc)) = 0.0
               end if
            END IF
            IF(qs(i,k).lt.qsmall .or. ns(i,k) .lt. qsmall)THEN
               qs(i,k) = 0.0
               ns(i,k) = 0.0
               as(i,k) = 0.0
               cs(i,k) = 0.0
               if (ipelec>0) then
                   an(i,1,k,lsc(ls)) = 0.0
               end if
            END IF
            IF(qg(i,k).lt.qsmall .or. ng(i,k) .lt. qsmall)THEN
               qg(i,k) = 0.0
               ng(i,k) = 0.0
               if (ipelec>0) then
                   an(i,1,k,lsc(lh)) = 0.0
               end if
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

            IF(snowflag.eq.1)THEN
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
            END IF!snowflag

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

      END SUBROUTINE AHM_MICRO

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
      SUBROUTINE SIMPLENUC(temp,sui,rhodum,nidum,dt,mnuccd,&
        nnuccd,nuc1,rdry,nin)
!---------------------------------------------------------------------

      IMPLICIT NONE

      REAL, INTENT(IN) :: temp, dt, nuc1, rhodum, nidum, rdry(nin)
      REAL, INTENT(OUT) :: mnuccd, nnuccd
      REAL :: dum, mi, wght, rdry1(nin), denom
      INTEGER, INTENT(IN) ::  nin
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
        
!     dum = nuc*1000./rhodum !1/L*1000/rho ~ 1000/kg
!IF(nidum.lt.dum)THEN
        IF(nidum .lt. dum) THEN !dum used to be nuc1
           nnuccd=(dum-nidum)/dt
        END IF


        mnuccd = nnuccd*mi
      END IF 

      END SUBROUTINE SIMPLENUC
!----------------------------------------------------------------------
      SUBROUTINE MEYERS(temp,sui,rhodum,nidum,dt,mnuccd,pgam,&
        lamc,cdist,mu,nnuccc,mnuccc,nnuccd,press,icontactflag,nuc1,&
        rdry,nin)
!---------------------------------------------------------------------
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


    !********************************************************
    !THIS SUBROUTINE COMPUTES DELTASTR, ICE DENSITY, AND ICE
    !EQUIVALENT VOL RADIUS AND SUBSEQUENTLY PERFORMS LIMIT
    !CHECKS.
    !********************************************************
    SUBROUTINE ICE_CHECKS(flag,n,ni,qi,ani,cni,rni,deltastr,rhobar,iaspect,&
         sphrflag,redden)
                 
      IMPLICIT NONE
      INTEGER iaspect, sphrflag, redden, flag
      REAL qi, ni, ani, cni, rni, deltastr, rhobar, n
      
      IF(flag.eq.2 .and. isnan(rhobar))print*,'BEFORE1', rhobar
      IF(flag.eq.2 .and. isnan(rni))print*,'BEFORE2', rni
      
      CALL DSTR_CHECK(flag,n,ni,ani,cni,deltastr,iaspect,sphrflag)
      IF(flag.eq.2 .and. isnan(rni))print*,'AFTERDSTR', rni
      CALL RHO_CHECK(flag,n,deltastr,qi,ni,ani,cni,sphrflag,redden,rhobar)
      IF(flag.eq.2 .and. isnan(rni))print*,'AFTERRHO1', rni
      IF(flag.eq.2 .and. isnan(rhobar))print*,'AFTERRHO2', rhobar
      CALL R_CHECK(flag,n,qi,ni,cni,ani,rni,rhobar,deltastr)
      
      IF(flag.eq.2 .and. isnan(rhobar))print*,'AFTER1', rhobar
      IF(flag.eq.2 .and. isnan(rni))print*,'AFTER2', rni
                 

    END SUBROUTINE ICE_CHECKS

    SUBROUTINE DSTR_CHECK(flag,n,ni,ani,cni,deltastr,iaspect,sphrflag)
      
      IMPLICIT NONE
      INTEGER :: iaspect,sphrflag,flag
      REAL ::  ani, cni, deltastr, voltmp, ai, ci, ni, gn, n
      
      IF(flag.eq.2)THEN
      IF(isnan(n))print*,'n1D',flag,n,ni,ani,cni,deltastr,iaspect,sphrflag
      IF(isnan(deltastr))print*,'deltastr1D',flag,n,ni,ani,cni,deltastr,iaspect,sphrflag
      IF(isnan(ni))print*,'ni1D',flag,n,ni,ani,cni,deltastr,iaspect,sphrflag
      IF(isnan(ani))print*,'ani1D',flag,n,ni,ani,cni,deltastr,iaspect,sphrflag
      IF(isnan(cni))print*,'cni1D',flag,n,ni,ani,cni,deltastr,iaspect,sphrflag
      END IF

      !     get deltastr from cni and ani
      !     deltastr = 1 for ice particles pre-diagnosed as spheres

      gn = exp(gammln(n))
      
      IF((log(ani)-log(ao)).gt. 0.01 &
           .and.(log(cni)-log(co)).gt.0.001)THEN
         deltastr = (log(cni)-log(co))/(log(ani)-log(ao))
      ELSE
         deltastr = 1.
      ENDIF
      
      IF(iaspect .eq. 1) deltastr = 0.8
      IF(sphrflag .eq. 1) deltastr = 1.0
      
               
      if(deltastr.lt.0.55) then
         voltmp=(4./3.)*pi*co/ao**(deltastr)*ani**(2.+deltastr)* &
              (exp(gammln(n+deltastr+2.)))/gn
         deltastr=0.55
         ani=((3.*voltmp*gn)/ &
              (4.*pi*co/ao**(deltastr)*(exp(gammln(n+deltastr+2.)))))** &
              (1./(2.+deltastr))
         ai=max((n*ni*ani),1.e-20)
         ani=ai/(n*ni)
      else if (deltastr.gt.1.5) then
         voltmp=(4./3.)*pi*co/ao**(deltastr)*ani**(2.+deltastr)* &
              (exp(gammln(n+deltastr+2.)))/gn
         deltastr=1.5
         ani=((3.*voltmp*gn)/ &
              (4.*pi*co/ao**(deltastr)*(exp(gammln(n+deltastr+2.)))))** &
              (1./(2.+deltastr))
         ai=max((n*ni*ani),1.e-20)
         ani=ai/(n*ni)
         cni=co*(ani/ao)**deltastr
         ci=max((n*ni*cni),1.e-20)
         cni=ci/(n*ni)
      endif
      
      IF(flag.eq.2)THEN
      IF(isnan(n))print*,'n2D',flag,n,ni,ani,cni,deltastr,iaspect,sphrflag
      IF(isnan(deltastr))print*,'deltastr2D',flag,n,ni,ani,cni,deltastr,iaspect,sphrflag
      IF(isnan(ni))print*,'ni2D',flag,n,ni,ani,cni,deltastr,iaspect,sphrflag
      IF(isnan(ani))print*,'ani2D',flag,n,ni,ani,cni,deltastr,iaspect,sphrflag
      IF(isnan(cni))print*,'cni2D',flag,n,ni,ani,cni,deltastr,iaspect,sphrflag
      END IF
      
    END SUBROUTINE DSTR_CHECK

    SUBROUTINE RHO_CHECK(flag,n,deltastr,qi,ni,ani,cni,sphrflag,redden,rhobar)
      IMPLICIT NONE
      INTEGER sphrflag, redden, flag
      REAL deltastr, qi, ni, ani, cni, n, gn
      REAL alphstr, alphv, betam
      REAL rhobar
      
      IF(flag.eq.2)THEN
      IF(isnan(n))print*,'n1RHO',flag,n,deltastr,qi,ni,ani,cni,sphrflag,redden,rhobar
      IF(isnan(deltastr))print*,'deltastr1RHO',flag,n,deltastr,qi,ni,ani,cni,sphrflag,redden,rhobar
      IF(isnan(qi))print*,'qi1RHO',flag,n,deltastr,qi,ni,ani,cni,sphrflag,redden,rhobar
      IF(isnan(ni))print*,'ni1RHO',flag,n,deltastr,qi,ni,ani,cni,sphrflag,redden,rhobar
      IF(isnan(ani))print*,'ani1RHO',flag,n,deltastr,qi,ni,ani,cni,sphrflag,redden,rhobar
      IF(isnan(cni))print*,'cni1RHO',flag,n,deltastr,qi,ni,ani,cni,sphrflag,redden,rhobar
      IF(isnan(rhobar))print*,'rhobar1RHO',flag,n,deltastr,qi,ni,ani,cni,sphrflag,redden,rhobar
      END IF

      gn = exp(gammln(n))
      
      betam = 2.+deltastr
      alphstr = co/ao**(deltastr)
      alphv = (4./3.)*pi*alphstr
      
      !     get avg ice density 
               
      rhobar = qi*gn/(ni*alphv* &
           ani**betam*exp(gammln(n+betam)))
               
      IF(sphrflag.eq.1) rhobar = 920.
      If(redden .eq. 1) rhoi = 500.

      IF(rhobar.gt.920.)THEN 
                      
         rhobar=920.
         ani=((qi*gn)/(rhobar*ni*alphv*&
              exp(gammln(n+betam))))**(1./betam)
         cni=co*(ani/ao)**deltastr
                  
      ELSE IF(rhobar.lt.50.)THEN
                      
         rhobar=50.
         ani=((qi*gn)/(rhobar*ni*alphv*&
              exp(gammln(n+betam))))**(1./betam)
         cni=co*(ani/ao)**deltastr
                      
      END IF
      
      IF(flag.eq.2)THEN
      IF(isnan(n))print*,'n2RHO',flag,n,deltastr,qi,ni,ani,cni,sphrflag,redden,rhobar
      IF(isnan(deltastr))print*,'deltastr2RHO',flag,n,deltastr,qi,ni,ani,cni,sphrflag,redden,rhobar
      IF(isnan(qi))print*,'qi2RHO',flag,n,deltastr,qi,ni,ani,cni,sphrflag,redden,rhobar
      IF(isnan(ni))print*,'ni2RHO',flag,n,deltastr,qi,ni,ani,cni,sphrflag,redden,rhobar
      IF(isnan(ani))print*,'ani2RHO',flag,n,deltastr,qi,ni,ani,cni,sphrflag,redden,rhobar
      IF(isnan(cni))print*,'cni2RHO',flag,n,deltastr,qi,ni,ani,cni,sphrflag,redden,rhobar
      IF(isnan(rhobar))print*,'rhobar2RHO',flag,n,deltastr,qi,ni,ani,cni,sphrflag,redden,rhobar
      END IF

      
    END SUBROUTINE RHO_CHECK

    !     get rni (characteristic equivalent volume ice radius)
    SUBROUTINE R_CHECK(flag,n,qi,ni,cni,ani,rni,rhobar,deltastr)
      IMPLICIT NONE
      INTEGER flag
      REAL qi, ni, cni, ani, rhobar, deltastr, n, gn
      REAL rni
      REAL alphstr, alphv, betam
      
      IF(flag.eq.2)THEN
      IF(isnan(n))print*,'n1R',flag,n,qi,ni,cni,ani,rni,rhobar,deltastr
      IF(isnan(deltastr))print*,'deltastr1R',flag,n,qi,ni,cni,ani,rni,rhobar,deltastr
      IF(isnan(qi))print*,'qi1R',flag,n,qi,ni,cni,ani,rni,rhobar,deltastr
      IF(isnan(ni))print*,'ni1R',flag,n,qi,ni,cni,ani,rni,rhobar,deltastr
      IF(isnan(ani))print*,'ani1R',flag,n,qi,ni,cni,ani,rni,rhobar,deltastr
      IF(isnan(cni))print*,'cni1R',flag,n,qi,ni,cni,ani,rni,rhobar,deltastr
      IF(isnan(rni))print*,'rni1R',flag,n,qi,ni,cni,ani,rni,rhobar,deltastr
      IF(isnan(rhobar))print*,'rhobar1R',flag,n,qi,ni,cni,ani,rni,rhobar,deltastr
      END IF

      gn = exp(gammln(n))
      
      betam = 2.+deltastr
      alphstr = co/ao**(deltastr)
      alphv = (4./3.)*pi*alphstr

      rni = (qi*3./(ni*rhobar*4.*pi*&
           (exp(gammln(n+deltastr+2.))/gn)))**(1./3.)
      
      
      !     make sure rni is within reasonable bounds,
      

      IF(rni.lt.2.e-6)THEN
         
         rni=2.e-6
         ni=3.*qi*gn/(4.*pi*rhobar*rni**3.*&
              (exp(gammln(n+deltastr+2.))))
         ani=((qi*gn)/(rhobar*ni*alphv* &
              exp(gammln(n+betam))))**(1./betam) 
         cni=co*(ani/ao)**deltastr
         
      ELSE IF(rni.gt.2.e-3)THEN
         
         rni=2.e-3
         ni=3.*qi*gn/(4.*pi*rhobar*rni**3.* &
              (exp(gammln(n+deltastr+2.))))
         ani=((qi*gn)/(rhobar*ni*alphv* &
              exp(gammln(n+betam))))**(1./betam)
         cni=co*(ani/ao)**deltastr
         
      END IF
      
     IF(flag.eq.2)THEN
      IF(isnan(n))print*,'n2R',flag,n,qi,ni,cni,ani,rni,rhobar,deltastr
      IF(isnan(deltastr))print*,'deltastr2R',flag,n,qi,ni,cni,ani,rni,rhobar,deltastr
      IF(isnan(qi))print*,'qi2R',flag,n,qi,ni,cni,ani,rni,rhobar,deltastr
      IF(isnan(ni))print*,'ni2R',flag,n,qi,ni,cni,ani,rni,rhobar,deltastr
      IF(isnan(ani))print*,'ani2R',flag,n,qi,ni,cni,ani,rni,rhobar,deltastr
      IF(isnan(cni))print*,'cni2R',flag,n,qi,ni,cni,ani,rni,rhobar,deltastr
      IF(isnan(rni))print*,'rni2R',flag,n,qi,ni,cni,ani,rni,rhobar,deltastr,gn, &
      (exp(gammln(n+deltastr+2.))/gn)
      IF(isnan(rhobar))print*,'rhobar2R',flag,n,qi,ni,cni,ani,rni,rhobar,deltastr
      END IF

    END SUBROUTINE R_CHECK


!This subroutine is used to calculate the collection of two different hydrometeor species,
!e.g., the collection of ice by snow. This is leiu of a 10D lookup table.
    SUBROUTINE COLL_COMPUTE(ni,ai,ci,rhoi,nui,ns,as,cs,rhos,nus)
    IMPLICIT NONE
    
    INTEGER :: n,o,p,q
    INTEGER, PARAMETER :: nn = 74
    REAL, PARAMETER :: E = 0.05
    REAL :: r(92)
    DATA r/0.5e-6,1.e-6,2.e-6,3.e-6,4.e-6,5.e-6,6.e-6,7.e-6,8.e-6,9.e-6,&
                  1.e-5,2.e-5,3.e-5,4.e-5,5.e-5,6.e-5,7.e-5,8.e-5,9.e-5,&
                  1.e-4,2.e-4,3.e-4,4.e-4,5.e-4,6.e-4,7.e-4,8.e-4,9.e-4,&
                  1.e-3,2.e-3,3.e-3,4.e-3,5.e-3,6.e-3,7.e-3,8.e-3,9.e-3,&
                  1.e-2,2.e-2,3.e-2,4.e-2,5.e-2,6.e-2,7.e-2,8.e-2,9.e-2,&
                  0.1,  0.2,  0.3,  0.4,  0.5,  0.6,  0.7,  0.8,  0.9,&
                  1.,   2.,   3.,   4.,   5.,   6.,   7.,   8.,   9.,&
                  10.,  20.,  30.,  40.,  50.,  60.,  70.,  80.,  90.,&
                  100., 200., 300., 400., 500., 600., 700., 800., 900.,&
                  1000.,2000.,3000.,4000.,5000.,6000.,7000.,8000.,9000.,10000./
    REAL, INTENT(IN) :: ni,ai,ci,rhoi,nui,ns,as,cs,rhos,nus
    REAL :: a1(nn),a2(nn),c1(nn),c2(nn),da1(nn),dr(nn),nrate(nn,nn,nn,nn),mrate(nn,nn,nn,nn)
    REAL :: maxx,maxy,vol1,vol2,n1,n2,v1,v2,maxx2,maxy2,coll,ncoll,area
    
    dr(1) = 0.5e-6
    DO n=2,nn,1
        dr(n) = (r(n)-r(n-1))
    END DO
    !print*,'COMPUTE'
    
    
    DO n=1,nn,1
        a1(n) = r(n)
        DO o=1,nn,1
            c1(o) = r(o)
            
            maxx = max(a1(n),c1(o))
            vol1 = a1(n)**2*c1(o)
            n1 = ni/exp(gammln(nui))**2 * (a1(n)/ai)**(nui-1)*(exp(-a1(n)/ai)/ai) *&
                                          (c1(o)/ci)**(nui-1)*(exp(-c1(o)/ci)/ci)
            maxx2 = maxx*2.e6                              
            if(maxx2 .lt. 41)then 
                v1 = 0.0028*(maxx2)**2.
            else if(maxx2 .ge. 41 .and. maxx2 .lt. 839)then 
                v1 = 0.0791*(maxx2)**1.101
            else 
                v1 = 62.29*(maxx2)**0.1098
            end if
            v1 = v1/100.
            
            
            DO p=1,nn,1
              a2(p) = r(p)
              DO q=1,nn,1
                  c2(q) = r(q)
                  
                  maxy = max(a2(p),c2(q))
                  area = (maxx+maxy)**2
                  vol2 = a2(p)**2*c2(q)
                  n2 = ns/exp(gammln(nus))**2 * (a2(p)/as)**(nus-1)*(exp(-a2(p)/as)/as) *&
                                                (c2(q)/cs)**(nus-1)*(exp(-c2(q)/cs)/cs)
                  maxy2 = maxy*2.e6                              
                  if(maxy2 .lt. 41)then 
                      v2 = 0.0028*(maxy2)**2.
                  else if(maxy2 .ge. 41 .and. maxy2 .lt. 839)then 
                      v2 = 0.0791*(maxy2)**1.101
                  else 
                      v2 = 62.29*(maxx2)**0.1098
                  end if
                  v2 = v2/100.
                  
                  nrate(q,p,o,n) = pi*area*abs(v1-v2)*n1*n2*E*dr(n)*dr(o)*dr(p)*dr(q)
                  mrate(q,p,o,n) = rhos*(vol1+vol2)*nrate(q,p,o,n)
              END DO!c2     
            END DO!a2
        END DO!c1
    END DO!a1
    ncoll = sum(sum(sum(sum(nrate,dim=4),dim=3),dim=2),dim=1)
    coll = sum(sum(sum(sum(mrate,dim=4),dim=3),dim=2),dim=1)
    
    
    
    END SUBROUTINE COLL_COMPUTE





    SUBROUTINE COLL_LOOKUP(ni,an,cn,rho,agg,nagg,aagg,anagg,cnagg,ddagg)

      IMPLICIT NONE
      INTEGER :: i, j, k, l, m, n, o, flag
      INTEGER :: si, sj, sk, sl, sm, sn, so
      REAL :: phi, rn
      REAL :: wt_ni, wt_an, wt_cn, wt_nu, wt_rho, wt_phi, wt_rn
      REAL, INTENT(IN) :: ni, an, cn, rho
      REAL, INTENT(OUT) :: agg, nagg, aagg, anagg, cnagg, ddagg


      phi = cn/an
      rn = (an**2*cn)**(1./3.)*1.e6 !!equivalent vol spherical radius in microns

      si = size(coll_ni)
      sj = size(coll_an)
      sk = size(coll_cn)
      sl = size(coll_nu)
      sm = size(coll_rho)
      sn = size(coll_phi)
      so = size(coll_r)
     
      !!ITERATE_LOOP RETURNS THE LOOKUP TABLE LOCATION AND THE ASSOCIATED WEIGHTING
      !!FOR THE MODELED VALUE
      CALL ITERATE_LOOKUP(1,coll_ni,  ni,  si,  wt_ni,  i)
      CALL ITERATE_LOOKUP(2,coll_an,  an,  sj,  wt_an,  j)
      CALL ITERATE_LOOKUP(3,coll_cn,  cn,  sk,  wt_cn,  k)
      CALL ITERATE_LOOKUP(4,coll_nu,  nu,  sl,  wt_nu,  l)
      CALL ITERATE_LOOKUP(5,coll_rho, rho, sm,  wt_rho, m)
      CALL ITERATE_LOOKUP(6,coll_phi, phi, sn,  wt_phi, n)
      CALL ITERATE_LOOKUP(7,coll_r,   rn,  so,  wt_rn,  o)
      
      !print*,"IN LOOKUP",ni,an,cn,rho,nu
      !print*,"LOOKUP VALS", i,coll_ni(i), j,coll_an(j),j,coll_cn(k),l,coll_nu(l),m,coll_rho(m)
      !print*,"WGHT",wt_ni,wt_an,wt_cn,wt_nu,wt_rho
      !print*,'COLL VALS', coll(i,j,k,l,m), coll(i,j-1,k-1,l,m), ncoll(i,j,k,l,m), ncoll(i,j-1,k-1,l,m),ncoll(1,2,2,4,9),ncoll(1,1,1,1,1)
      !print*,"IPAS VALS", acoll(1,1),ancoll(1,1),cncoll(1,1),ddcoll(1,1)

      CALL WGHTED_LOOKUP5D(i,j,k,l,m,si,sj,sk,sl,sm,wt_ni,wt_an,wt_cn,wt_nu,wt_rho,coll,agg)
      CALL WGHTED_LOOKUP5D(i,j,k,l,m,si,sj,sk,sl,sm,wt_ni,wt_an,wt_cn,wt_nu,wt_rho,ncoll,nagg)

      CALL WGHTED_LOOKUP2D(n,o,sn,so,wt_phi,wt_rn,acoll,aagg)
      CALL WGHTED_LOOKUP2D(n,o,sn,so,wt_phi,wt_rn,ancoll,anagg)
      CALL WGHTED_LOOKUP2D(n,o,sn,so,wt_phi,wt_rn,cncoll,cnagg)
      CALL WGHTED_LOOKUP2D(n,o,sn,so,wt_phi,wt_rn,ddcoll,ddagg)
      
    END SUBROUTINE COLL_LOOKUP


    !ITERATE_LOOKUP: THIS SUBROUTINE DETERMINES THE WEIGHTING BETWEEN TWO INDICES OF A
    !PARTICULAR DIMENSION BASED ON THE LOOKUP TABLE VALUES AND THE MODELED VALUES
    SUBROUTINE ITERATE_LOOKUP(f,lookup_value,modeled_value,lookup_size,wght,i)
      
      IMPLICIT NONE
      INTEGER :: n, value_found_flag, f
      INTEGER, INTENT(IN) :: lookup_size
      INTEGER, INTENT(OUT) :: i 
      REAL, INTENT(IN) :: lookup_value(lookup_size), modeled_value
      REAL, INTENT(OUT) :: wght
      
      if(f.eq.1)then
         !print*,'lookup values: ',x,' actual value: ',y,' array size: ',s
      end if

      value_found_flag = 0
      i = lookup_size  
      DO n = 1, lookup_size   !increment lookup table until lookup value >= modeled value
         if(lookup_value(n) .ge. modeled_value .and. value_found_flag .eq. 0)then
            i = n
            value_found_flag = 1
         end if
      END DO

      !if(f.eq.1)print*,'upper index: ',i
      
      if(i .gt. 1)then                  !when the modeled value exsits w/in bounds of lookup
         wght = (lookup_value(i)-modeled_value)/(lookup_value(i)-lookup_value(i-1))
      else if(i.eq.1)then               !when the modeled value <= the smallest lookup value
         wght = (lookup_value(i)-modeled_value)/(lookup_value(i)-0.0)
      else if(i.eq.lookup_size.and.value_found_flag.eq.0)then !when the modeled value >= the largest lookup value
         wght = 0.0
      end if
      
    END SUBROUTINE ITERATE_LOOKUP


    !!WGHTED_LOOKUP2D: THIS SUBROUTINE TAKES THE WEIGHTING OF A PARTICULR DIMENSION AND 
    !!AND DETERMINES THE LOOKUP TABLE-WEIGHTED VALUE BASED ON THE CHOSEN LOOKUP TABLE
    !INDEX. THIS IS COMPLETED IN 2 DIMENSIONS FOR A 2D ARRAY.
    SUBROUTINE WGHTED_LOOKUP2D(i,j,si,sj,wghti,wghtj,x,y)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: i,j,si,sj
      REAL, INTENT(IN) :: wghti,wghtj,x(si,sj)
      REAL, INTENT(OUT) :: y
      REAL tmpi1,tmpi2,tmpj,x1,x2,x3,x4

      !i1 = wght i at j1
      !i2 = wght i at j2

      x1 = x(i-1,j-1)
      x2 = x(i,  j-1)
      x3 = x(i-1,j  )
      x4 = x(i,  j  )
      if(i.eq.1)then
         x1 = 0
         x3 = 0
      else if(j.eq.1)then
         x1 = 0
         x2 = 0
      end if
         
      tmpi1 = x1*wghti + x2*(1.-wghti) !i wght
      tmpi2 = x3*wghti + x4*(1.-wghti)     !i wght
         
      !i-wghted j
      tmpj = tmpi1*wghtj + tmpi2*(1.-wghtj)                                  !j wght

      y = tmpj

    END SUBROUTINE WGHTED_LOOKUP2D

    !!WGHTED_LOOKUP5D: THIS SUBROUTINE TAKES THE WEIGHTING OF A PARTICULR DIMENSION AND 
    !!AND DETERMINES THE LOOKUP TABLE-WEIGHTED VALUE BASED ON THE CHOSEN LOOKUP TABLE
    !INDEX. THIS IS COMPLETED IN 5 DIMENSIONS FOR A 5D ARRAY.
    SUBROUTINE WGHTED_LOOKUP5D(i,j,k,l,m,si,sj,sk,sl,sm,wghti,wghtj,wghtk,wghtl,wghtm,x,y)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: i,j,k,l,m,si,sj,sk,sl,sm
      REAL, INTENT(IN) :: wghti,wghtj,wghtk,wghtl,wghtm,x(si,sj,sk,sl,sm)
      REAL, INTENT(OUT) :: y
      REAL tmpi1,tmpi2,tmpj1,tmpj2,tmpk1,tmpk2,tmpl1,tmpl2,tmpm
      REAL x1,x2,x3,x4

      !i1 = wght i at j1,k1,l1,m1
      !i2 = wght i at j2,k1,l1,m1

      x1 = get_x(i-1,j-1,k-1,l-1,m-1,si,sj,sk,sl,sm,x)
      x2 = get_x(i,  j-1,k-1,l-1,m-1,si,sj,sk,sl,sm,x)
      x3 = get_x(i-1,j,  k-1,l-1,m-1,si,sj,sk,sl,sm,x)
      x4 = get_x(i,  j,  k-1,l-1,m-1,si,sj,sk,sl,sm,x)

      tmpi1 = x1*wghti + x2*(1.-wghti) !i wght
      tmpi2 = x3*wghti + x4*(1.-wghti)     !i wght
         
      !j1 = i-wghted j at k1,l1,m1
      tmpj1 = tmpi1*wghtj + tmpi2*(1.-wghtj)                                  !j wght

      !----------------------------------------------------------------------------------------

      !i1 = wght i at j1,k2,l1,m1
      !i2 = wght i at j2,k2,l1,m1

      x1 = get_x(i-1,j-1,k,  l-1,m-1,si,sj,sk,sl,sm,x)
      x2 = get_x(i,  j-1,k,  l-1,m-1,si,sj,sk,sl,sm,x)
      x3 = get_x(i-1,j,  k,  l-1,m-1,si,sj,sk,sl,sm,x)
      x4 = get_x(i,  j,  k,  l-1,m-1,si,sj,sk,sl,sm,x)

      tmpi1 = x1*wghti + x2*(1.-wghti)!i wght
      tmpi2 = x3*wghti + x4*(1.-wghti)    !i wght

      !j2 = i-wghted j at k2,l1,m1
      tmpj2 = tmpi1*wghtj + tmpi2*(1.-wghtj)                                  !j wght

      !----------------------------------------------------------------------------------------

      !k1 = (i,j)-wghted k at l1,m1
      tmpk1 = tmpj1*wghtk + tmpj2*(1.-wghtk)                                  !k wght

      !========================================================================================      

      !i1 = wght i at j1,k1,l2,m1
      !i2 = wght i at j2,k1,l2,m1

      x1 = get_x(i-1,j-1,k-1,l,  m-1,si,sj,sk,sl,sm,x)
      x2 = get_x(i,  j-1,k-1,l,  m-1,si,sj,sk,sl,sm,x)
      x3 = get_x(i-1,j,  k-1,l,  m-1,si,sj,sk,sl,sm,x)
      x4 = get_x(i,  j,  k-1,l,  m-1,si,sj,sk,sl,sm,x)

      tmpi1 = x1*wghti + x2*(1.-wghti)     !i wght
      tmpi2 = x3*wghti + x4*(1.-wghti)         !i wght

      !j1 = i-wghted j at k1,l2,m1
      tmpj1 = tmpi1*wghtj + tmpi2*(1.-wghtj)                                 !j wght

      !----------------------------------------------------------------------------------------

      !i1 = wght i at j1,k2,l2,m1
      !i2 = wght i at j2,k2,l2,m1

      x1 = get_x(i-1,j-1,k,  l,  m-1,si,sj,sk,sl,sm,x)
      x2 = get_x(i,  j-1,k,  l,  m-1,si,sj,sk,sl,sm,x)
      x3 = get_x(i-1,j,  k,  l,  m-1,si,sj,sk,sl,sm,x)
      x4 = get_x(i,  j,  k,  l,  m-1,si,sj,sk,sl,sm,x)

      tmpi1 = x1*wghti + x2*(1.-wghti)          !i wght
      tmpi2 = x3*wghti + x4*(1.-wghti)              !i wght

      !j2 = i-wghted j at k2,l2,m1
      tmpj2 = tmpi1*wghtj + tmpi2*(1.-wghtj)                                  !j wght

      !----------------------------------------------------------------------------------------

      !k2 = (i,j)-wghted k at l2,m1
      tmpk2 = tmpj1*wghtk + tmpj2*(1.-wghtk)                                  !k wght

      !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      !l1 = (i,j,k)-wghted l at m1
      tmpl1 = tmpk1*wghtl + tmpk2*(1.-wghtl)                                  !l wght

      !||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

      !i1 = wght i at j1,k1,l1,m2
      !i2 = wght i at j2,k1,l1,m2

      x1 = get_x(i-1,j-1,k-1,l-1,m,  si,sj,sk,sl,sm,x)
      x2 = get_x(i,  j-1,k-1,l-1,m,  si,sj,sk,sl,sm,x)
      x3 = get_x(i-1,j,  k-1,l-1,m,  si,sj,sk,sl,sm,x)
      x4 = get_x(i,  j,  k-1,l-1,m,  si,sj,sk,sl,sm,x)

      tmpi1 = x1*wghti + x2*(1.-wghti)      !i wght
      tmpi2 = x3*wghti + x4*(1.-wghti)          !i wght

      !j1 = i-wghted j at k1,l1,m2
      tmpj1 = tmpi1*wghtj + tmpi2*(1.-wghtj)                                  !j wght

      !----------------------------------------------------------------------------------------

      !i1 = wght i at j1,k2,l1,m2
      !i2 = wght i at j2,k2,l1,m2

      x1 = get_x(i-1,j-1,k,  l-1,m,  si,sj,sk,sl,sm,x)
      x2 = get_x(i,  j-1,k,  l-1,m,  si,sj,sk,sl,sm,x)
      x3 = get_x(i-1,j,  k,  l-1,m,  si,sj,sk,sl,sm,x)
      x4 = get_x(i,  j,  k,  l-1,m,  si,sj,sk,sl,sm,x)

      tmpi1 = x1*wghti + x2*(1.-wghti)          !i wght
      tmpi2 = x3*wghti + x4*(1.-wghti)              !i wght

      !j2 = i-wghted j at k2,l1,m2
      tmpj2 = tmpi1*wghtj + tmpi2*(1.-wghtj)                                  !j wght

      !----------------------------------------------------------------------------------------

      !k1 = (i,j)-wghted k at l1,m2
      tmpk1 = tmpj1*wghtk + tmpj2*(1.-wghtk)                                  !k wght

      !========================================================================================      

      !i1 = wght i at j1,k1,l2,m2
      !i2 = wght i at j2,k1,l2,m2

      x1 = get_x(i-1,j-1,k-1,l,  m,  si,sj,sk,sl,sm,x)
      x2 = get_x(i,  j-1,k-1,l,  m,  si,sj,sk,sl,sm,x)
      x3 = get_x(i-1,j,  k-1,l,  m,  si,sj,sk,sl,sm,x)
      x4 = get_x(i,  j,  k-1,l,  m,  si,sj,sk,sl,sm,x)

      tmpi1 = x1*wghti + x2*(1.-wghti)          !i wght
      tmpi2 = x3*wghti + x4*(1.-wghti)              !i wght

      !j1 = i-wghted j at k1,l2,m2
      tmpj1 = tmpi1*wghtj + tmpi2*(1.-wghtj)                                  !j wght

      !----------------------------------------------------------------------------------------

      !i1 = wght i at j1,k2,l2,m2
      !i2 = wght i at j2,k2,l2,m2

      x1 = get_x(i-1,j-1,k,  l,  m,  si,sj,sk,sl,sm,x)
      x2 = get_x(i,  j-1,k,  l,  m,  si,sj,sk,sl,sm,x)
      x3 = get_x(i-1,j,  k,  l,  m,  si,sj,sk,sl,sm,x)
      x4 = get_x(i,  j,  k,  l,  m,  si,sj,sk,sl,sm,x)

      tmpi1 = x1*wghti + x2*(1.-wghti)              !i wght
      tmpi2 = x3*wghti + x4*(1.-wghti)                  !i wght

      !j2 = i-wghted j at k2,l2,m2
      tmpj2 = tmpi1*wghtj + tmpi2*(1.-wghtj)                                  !j wght

      !----------------------------------------------------------------------------------------

      !k2 = (i,j)-wghted k at l2,m2
      tmpk2 = tmpj1*wghtk + tmpj2*(1.-wghtk)                                  !k wght

      !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      !l2 = (i,j,k)-wghted l at m2
      tmpl2 = tmpk1*wghtl + tmpk2*(1.-wghtl)                                  !l wght

      !m = (i,j,k,l)-wghted m
      tmpm = tmpl1*wghtm + tmpl2*(1.-wghtm)                                   !m wght

      y = tmpm

    END SUBROUTINE WGHTED_LOOKUP5D


      real function get_x(i,j,k,l,m,si,sj,sk,sl,sm,x)
      implicit none
      integer i,j,k,l,m,si,sj,sk,sl,sm
      real x(si,sj,sk,sl,sm)
      
      if(i.eq.0.or.j.eq.0.or.k.eq.0.or.l.eq.0.or.m.eq.0)then
         get_x = 0
      else 
         get_x = x(i,j,k,l,m)
      end if
      return 
      end function get_x

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

      SUBROUTINE EVOLVE(nn,ani,ni,sui,sup,qvv,temp,press,igr,deltt,iwc,cf,&
          iwci,phii,phif,cni,rni,rnf,anf,deltastr,mu,rhoavg,vtbarb,&
          vtbarbm,vtbarblen,rhoa,i,k,iaspect,&
          ipart,sphrflag,redden,itimestep)
!***********************************************************************

      IMPLICIT NONE

      INTEGER ipart,ivent,i,k,iaspect,redden
      INTEGER sphrflag,itimestep,MICRO_ON
      REAL ani,ni,temp,press,qvv,nn
      REAL igr,gi,capgam,dmdt,deltt,qe
      REAL dlndt,mu
      REAL l,lmean,fs,gammnu,gammnu1,gammnubet,anf,betam
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
      grav = 9.81
      celsius = temp-273.15
      fv = 1.0
      fh = 1.0
      gi = findgtp(temp,press,fv,fh)
      dvs = 0.211*(temp/273.15)**1.94 * &
          (1013.25*100./press)*1.0/100.0**2.0 ! vapor diffusivity
      kts = 2.3823e-2 +7.1177e-5*(temp-273.15) ! thermal conductivity

      l=ani
      drho = ((polysvp1(temp,0)-polysvp1(temp,1))/(Rv*temp))*1000.0

      alphstr=co/ao**(deltastr)
      betam = 2. + deltastr
      gammnu=exp(gammln(nn))
      gammnu1=exp(gammln(nn+1.0))
      gammnubet=exp(gammln(nn+betam))

      capgam = capacitance_gamma(l,deltastr,nn,alphstr)
      lmean = l*nn

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

      gammnubet = exp(gammln(nn+betam))
      gammnu2delt = exp(gammln(nn+2.+deltastr))
      gammnu1delt = exp(gammln(nn+deltastr-1.))
      gamrats = exp(gammln(nn+betam/3.))/gammnu
      
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
      
      wghtv1 = exp(gammln(nn+deltastr+2.+2.*bl-ba))/gammnu
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
      
      wghtv2 = exp(gammln(nn+bx*bm-1.))/gammnu
      vtbarb = etaa/rhoa*0.5 * am*(xn)**bm*ani**(bx*bm-1.) * wghtv2
      wghtv3 = exp(gammln(nn+bx*bm-1.+2.+deltastr))&
      /exp(gammln(nn+2.+deltastr))
      vtbarbm = etaa/rhoa*0.5 * am*(xn)**bm*ani**(bx*bm-1.) * wghtv3

 
! add length
      if(phii.lt.1.0)then
         wghtv3 = exp(gammln(nn+bx*bm-1.+1)) &
         /exp(gammln(nn+1.))
         vtbarblen = etaa/rhoa*0.5 * am*(xn)**bm*ani**(bx*bm-1.) &
         * wghtv3
      else if(phii.gt.1.0)then
         wghtv3 = exp(gammln(nn+bx*bm-1.+deltastr))&
         /exp(gammln(nn+deltastr))
         vtbarblen = etaa/rhoa*0.5 * am*(xn)**bm*ani**(bx*bm-1.) &
         * wghtv3
      else if (phii.eq.1.0)then
         wghtv3 = exp(gammln(nn+bx*bm-1.+1))&
         /exp(gammln(nn+1.))
         vtbarblen = etaa/rhoa*0.5 * am*(xn)**bm*ani**(bx*bm-1.) &
         * wghtv3
      endif

         
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
      fallcheck = sqrt(dvs*pi*2*cni/(vt*nn))
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
         *exp(gammln(nn+deltastr+2.))/gammnu
         Vmin = 4./3.*pi*ao**3
         if(Vmin.lt.videp)then
            betavol = alog(rhoi/rhoavg)*1./(alog(Vmin/videp))
            rhodep = rhoavg*(1.+betavol)
         else
            rhodep = rhoavg 
         endif
      endif


      If(redden .eq. 1) rhodep = 500.
      !rhodep=500.
!     rf = ((rni*gamrats)**2 + 2.*gi*sui*fs/rhodep &
!           *gammnu/gammnubet*gamrats**3*deltt)**(0.5) ! rmean after deltat

      rf = (max(((rni*gamrats)**2 + 2.*gi*sui*fs/rhodep &
           *gammnu/gammnubet*gamrats**3*deltt),1.e-20))**(0.5) ! rmean after deltat
      rnf = rf/gamrats           ! convert to rn (characteristic r)
      anf = alphanr*rnf**(3./(2.+igr)) ! characteristic a-axis after deltat
      vi = 4./3.*pi*rni**3 * exp(gammln(nn+deltastr+2.))/gammnu ! mean initial volume
      vf = 4./3.*pi*rnf**3 * exp(gammln(nn+deltastr+2.))/gammnu ! mean volume after deltat
      rhoavg = rhoavg*(vi/vf) + rhodep*(1.-vi/vf) ! new average ice density
      rhoavg = max(min(rhoavg,920.),50.)

      !rhoavg=500.
      iwc = ni*rhoavg*vf        ! IWC after deltat

      if(igr.ne.1.0.and.(log(anf)-log(ao)).gt.0.01 .and.&
      (3.*log(rnf)-2.*log(anf)-log(ao)).gt.0.001)THEN
         deltastr = (3.*log(rnf)-2.*log(anf)-log(ao)) &
         /(log(anf)-log(ao))    ! diagnose new deltastar
      else
         deltastr = 1.
      endif
      IF(iaspect .eq. 1) deltastr = 0.8
      deltastr = min(max(deltastr,0.55),1.5)
  
      if(deltastr.ge.1.0)then   ! if columns, get c from c-a relation
         cf = (co/ao**(deltastr))*anf**deltastr
      endif

      phif = phii*(rnf**3/r**3)**((igr-1.)/(igr+2.)) ! if plates, get c from aspect ratio after deltat
      IF(iaspect .eq.1) phif = 0.27

!      rnf = rf
      if(deltastr.lt.1.0)then
         cf = phif*anf*gammnu/exp(gammln(nn+deltastr-1.))
      endif


 
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

          !    -----------------------------------------------------------------------!
          !    ------------------------  Explicit Lightning  -------------------------!
          Function delbk(bb,nu,mu,k)
              !
              !  Purpose: Caluculates collection coefficients following Siefert (2006)
              !
              !  delbk is equation (90) (b collecting b -- self-collection)
              !  mass-diameter relationship: D = a*x**(b), where x = particle mass
              !  general distribution: n(x) = A*x**(nu)*Exp(-lam*x**(mu))
              !  where
              !      A = mu*N/(Gamma((nu+1)/mu)) *lam**((nu+1)/mu)
              !
              !      lam = ( Gamma((nu+1)/mu)/Gamma((nu+2)/mu) * xbar )**(-mu)
              !
              !     where  xbar = L/N  (mass content)/(number concentration) = q*rhoa/N
              !

              implicit none
              real delbk
              real nu, mu, bb
        integer k, igam
        integer ngm0,ngm1,ngm2
        parameter (ngm0=3001,ngm1=500,ngm2=500)
        double precision, parameter :: dgam = 0.01, dgami = 100.
        double precision :: arg
        double precision gmoi(0:ngm0) ! ,gmod(0:ngm1,0:ngm2),gmdi(0:ngm1,0:ngm2)

              real tmp, del
              real x1, x2, x3, x4
              integer i

        gmoi(0) = 1.d32
        do igam = 1,ngm0
            arg = dgam*igam
            gmoi(igam) = gamma_dp(arg)
        end do

              tmp = ((1.0 + nu)/mu)
              i = Int(dgami*(tmp))
              del = tmp - dgam*i
              x1 = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami

              tmp = ((2.0 + nu)/mu)
              i = Int(dgami*(tmp))
              del = tmp - dgam*i
              x2 = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami

              tmp = ((1.0 + 2.0*bb + k + nu)/mu)
              i = Int(dgami*(tmp))
              del = tmp - dgam*i
              x3 = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami

              !      delbk =  &
              !     &  ((Gamma_sp((1.0 + nu)/mu)/Gamma_sp((2.0 + nu)/mu))**(2.0*bb + k)* &
              !     &    Gamma_sp((1.0 + 2.0*bb + k + nu)/mu))/Gamma_sp((1.0 + nu)/mu)

              delbk =  &
                      &  ((x1/x2)**(2.0*bb + k)* &
                      &    x3)/x1

              RETURN
          END  Function delbk
          !
          !    !+---+-----------------------------------------------------------------+
          !
          Function delabk(ba,bb,nua,nub,mua,mub,k)

              implicit none
              real delabk
              real nua, mua, ba
        integer k, igam
              real nub, mub, bb

              integer i
              real tmp,del

              real g1pnua, g2pnua, g1pbapnua, g1pbbpk, g1pnub, g2pnub

        integer ngm0,ngm1,ngm2
        parameter (ngm0=3001,ngm1=500,ngm2=500)
        double precision, parameter :: dgam = 0.01, dgami = 100.
        double precision :: arg
        double precision gmoi(0:ngm0) ! ,gmod(0:ngm1,0:ngm2),gmdi(0:ngm1,0:ngm2)

        gmoi(0) = 1.d32
        do igam = 1,ngm0
            arg = dgam*igam
            gmoi(igam) = gamma_dp(arg)
        end do

              tmp = (1. + nua)/mua
              i = Int(dgami*(tmp))
              del = tmp - dgam*i
              IF ( i+1 > ngm0 ) THEN
                  write(0,*) 'delabk: i+1 > ngm0!!!!',i,ngm0,nua,mua,tmp
            write (0,*) 'stopping'
            print*, 'stopping'
                  STOP
              ENDIF
              g1pnua = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami
              !        write(91,*) 'delabk: g1pnua,gamma = ',g1pnua,Gamma_sp((1. + nua)/mua)

              tmp = ((2. + nua)/mua)
              i = Int(dgami*(tmp))
              del = tmp - dgam*i
              g2pnua = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami

              tmp = ((1. + ba + nua)/mua)
              i = Int(dgami*(tmp))
              del = tmp - dgam*i
              g1pbapnua = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami

              tmp = ((1. + nub)/mub)
              i = Int(dgami*(tmp))
              del = tmp - dgam*i
              g1pnub = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami

              tmp = ((2 + nub)/mub)
              i = Int(dgami*(tmp))
              del = tmp - dgam*i
              g2pnub = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami

              tmp = ((1. + bb + k + nub)/mub)
              i = Int(dgami*(tmp))
              del = tmp - dgam*i
              g1pbbpk = gmoi(i) + (gmoi(i+1) - gmoi(i))*del*dgami

              delabk =  &
                      &  (2.*(g1pnua/g2pnua)**ba*     &
                      &    g1pbapnua*                                               &
                      &    (g1pnub/g2pnub)**(bb + k)*                                &
                      &    g1pbbpk)/                                                &
                      &  (g1pnua*g1pnub)

              RETURN
          END Function delabk
!
!    !+---+-----------------------------------------------------------------+
!
    DOUBLE PRECISION FUNCTION GAMMA_DP(xx)

        implicit none
        double precision xx
        integer j

        ! Double precision ser,stp,tmp,x,y,cof(6)

        real*8 ser,stp,tmp,x,y,cof(6)
        SAVE cof,stp
        DATA cof,stp/76.18009172947146d+0,  &
                &            -86.50532032941677d0,   &
                &             24.01409824083091d0,   &
                &             -1.231739572450155d0,  &
                &              0.1208650973866179d-2,&
                &             -0.5395239384953d-5,   &
                &              2.5066282746310005d0/

        x = xx
        y = x
        tmp = x + 5.5d0
        tmp = (x + 0.5d0)*Log(tmp) - tmp
        ser = 1.000000000190015d0
        DO j=1,6
            y = y + 1.0d0
            ser = ser + cof(j)/y
        END DO
        gamma_dp = Exp(tmp + log(stp*ser/x))

        RETURN
    END function gamma_dp

    ! ####################################################################
    !                SUBROUTINE TAKAX
    ! ####################################################################
    subroutine takax(isaund0,nt,nlc,lookup,temcg,qcw,vt,awdia,     &
            &                 rho0,ftlwc,exw,qrw,rarfac)

        ! ####################################################################
!
        !  Purpose: For calling different versions of Takahashi charging
        !           routine
!
!
        ! ####################################################################

        implicit none

        integer isaund0,isaund,nt,nlc
        real lookup(0:nt,0:nlc)
        real qcw, qrw, lwc
        real rho0
        real ftlwc
        real temcg
        real vt
        real awdia,awdia1
        real exw
        real rarfac

        lwc = qcw
        isaund = Abs(isaund0)
        IF ( isaund0 .eq. -1 ) lwc = (qcw + qrw)
        IF ( isaund0 .eq. -2 ) lwc = (qcw + qrw)*rarfac

        !  BEGIN EXECUTABLE CODE
        IF ( isaund .eq. 1 ) THEN ! taka size/vel. depend.
            call taka(nt,nlc,lookup,temcg,lwc,vt,awdia,rho0,ftlwc)
        ELSEIF ( isaund .eq. 2 ) THEN ! saunders size dependence
            awdia1 = awdia !  3.67*awdia
            call taka2(nt,nlc,lookup,temcg,lwc,vt,awdia1,rho0,ftlwc)
        ELSEIF ( isaund .eq. 3 ) THEN ! no size or velocity dependence
            call taka3(nt,nlc,lookup,temcg,lwc,vt,awdia,rho0,ftlwc)
        ELSEIF ( isaund .eq. 4 ) THEN ! taka size/vel. depend. w/ D = 3.67/lambda
            call taka4(nt,nlc,lookup,temcg,lwc,vt,awdia,rho0,ftlwc)
        ELSEIF ( isaund .eq. 5 ) THEN ! taka size/vel. depend. and using RAR
            call takarar(nt,nlc,lookup,temcg,lwc,vt,awdia,rho0,ftlwc,exw)
        ELSE
            call taka(nt,nlc,lookup,temcg,lwc,vt,awdia,rho0,ftlwc)
        ENDIF

    END subroutine takax
!
    ! ####################################################################
    !                SUBROUTINE TAKA2
    ! ####################################################################
    subroutine TAKA2(nt,nlc,lookup,temcg,qcw,vt,awdia,rho0,ftlwc)

        ! ####################################################################
!
        !  Purpose:
        !           Noninductive charging per Takahashi 1978 lab results.
        !           Using lookup table from Wojcik (1994).  For temperatures
        !           lower than -30 C, the values at -30 C are used.
!
        !  7/31/99 Now limit low-temp charging rates by a parabolic function
        !          with a value of zero at -30 C and goes to zero at -43 C.
        !
        !  8/10/99 Fixed bug in the interpolation: was accidentally
        !          _extrapolating_.   Oops.
        !
        !  10/26/99  Test using a Saunders et al. size dependence on the crystal size
        !
        ! ####################################################################

        implicit none

        integer nt,nlc
        real lookup(0:nt,0:nlc)
        real qlwc
        real qcw
        real rho0
        real ftlwc
        real temt,temcg
        integer item,ilwc
        real fact,facl
        real a1,a2,a3
        real alf
        real vt
        real awdia
        real d0
        parameter(d0 = 100.0e-6)
        real v0
        parameter(v0 = 8.0)
        real fac
        real qconkq, qconm, qconn

        !  BEGIN EXECUTABLE CODE

        IF ( temcg .gt. -30 ) THEN
            fac = 1.0
        ELSEIF ( temcg .gt. -40.0 ) THEN
            fac = 1.0 - ((temcg+30.0)/(40.0 - 30.0))**2
        ELSE
            fac = 0.0
        END IF


        temt = Max (temcg, -30.0)
        item = Int(-temt) + 1
        qlwc = qcw*(1.e3)*rho0

        qlwc = Min ( qlwc, 30.0)
        ilwc = 0
        ftlwc = 0.0

        IF (qlwc .ge. 10.0) THEN
            ilwc = Int( qlwc/10.0 + 1.0e-5) + 1 + 27
            facl = (qlwc - lookup(0,ilwc) )*0.1
        ELSEIF (qlwc .ge. 1.0) THEN
            ilwc = Int( qlwc + 1.0e-5) + 1 + 18
            facl = (qlwc - lookup(0,ilwc) )*1.0
        ELSEIF (qlwc .ge. 0.1) THEN
            ilwc = Int( qlwc*10.0 + 1.0e-5 ) + 1 + 9
            facl = (qlwc - lookup(0,ilwc) )*10.0
        ELSEIF (qlwc .ge. 0.01) THEN
            ilwc = Int ( qlwc*100.0 + 1.0e-5) + 1
            facl = (qlwc - lookup(0,ilwc) )*100.0
        ELSE
            ilwc = 0
            ftlwc = 0.0
        END IF

        IF (ilwc .gt. 0) THEN
            fact = -temt - Float(item) + 1.0
            !        a1 = (1.0-fact)*lookup(item,ilwc) + fact*lookup(item+1,ilwc)
            !        a2 = (1.0-fact)*lookup(item,ilwc+1) + fact*lookup(item+1,ilwc+1)
            a1 = (1.0-fact)*lookup(item,ilwc) +      &
                    &          fact*lookup(Min(nt,item+1),ilwc)
            a2 = (1.0-fact)*lookup(item,ilwc+1) +      &
                    &          fact*lookup(Min(nt,item+1),ilwc+1)
            ftlwc = fac*(a1 + facl*(a2-a1))
            IF ( ftlwc .gt. 0.0 ) THEN

                if ( awdia*1.e6.lt.155. ) then
                    qconkq = 4.9e13
                    qconm  = 3.76
                    qconn  = 2.5
                end if
                if ( awdia*1.e6.ge.155. .and. awdia*1.e6.le.452. ) then
                    qconkq = 4.0e6
                    qconm  = 1.9
                    qconn  = 2.5
                end if
                if ( awdia*1.e6.gt.452. ) then
                    qconkq = 52.8
                    qconm  = 0.44
                    qconn  = 2.5
                end if

            ELSEIF ( ftlwc .lt. 0.0 ) THEN

                if ( awdia*1.e6.lt.253. ) then
                    qconkq = 5.24e8
                    qconm  = 2.54
                    qconn  = 2.8
                end if
                if ( awdia*1.e6.gt.253. ) then
                    qconkq = 24.0
                    qconm  = 0.50
                    qconn  = 2.8
                end if

            ENDIF

            !        alf = Min(5.0*(awdia/d0)**2*Abs(vt)/v0,10.0)
            !        alf = qconkq*(3.67*awdia)**qconm*(0.333*Abs(vt))**qconn ! 11/9/01 erm
            !        alf = qconkq*(awdia)**qconm*(Abs(vt))**qconn        ! 11/9/01 erm
            ! Need to figure out what to do with qconkq, if anything...
            alf = (awdia/d0)**qconm*(Abs(vt)/v0)**qconn        ! 12/9/04 erm
            ftlwc = 1.0e-15*alf*ftlwc
        END IF

        RETURN
    END subroutine TAKA2

    ! ####################################################################
    !                SUBROUTINE TAKARAR
    ! ####################################################################
    subroutine TAKARAR(nt,nlc,lookup,temcg,qcw,vt,awdia,rho0,ftlwc,     &
            &                   exw)

        ! ####################################################################
        !
        !  Purpose:
        !           Noninductive charging per Takahashi 1978 lab results.
        !           Using lookup table from Wojcik (1994).  For temperatures
        !           lower than -30 C, the values at -30 C are used.
        !
        !  7/31/99 Now limit low-temp charging rates by a parabolic function
        !          with a value of zero at -30 C and goes to zero at -43 C.
!
        !  8/10/99 Fixed bug in the interpolation: was accidentally
        !          _extrapolating_.   Oops.
!
        ! ####################################################################

        implicit none

        integer nt,nlc
        real lookup(0:nt,0:nlc)
        real qlwc
        real qcw
        real rho0
        real ftlwc
        real exw
        real temt,temcg
        integer item,ilwc
        real fact,facl
        real a1,a2,a3
        real alf
        real vt
        real awdia
        real d0
        parameter(d0 = 100.0e-6)
        real v0
        parameter(v0 = 8.0)
        real fac
        real rar

        !  BEGIN EXECUTABLE CODE

        IF ( temcg .gt. -30 ) THEN
            fac = 1.0
        ELSEIF ( temcg .gt. -40.0 ) THEN
            fac = 1.0 - ((temcg+30.0)/(40.0 - 30.0))**2
        ELSE
            fac = 0.0
        END IF
!
        ! here, qlwc is scaled by vt/9.0, where 9.0m/s is the riming rod speed in Takahashi (1978)
!
        qlwc = exw*qcw*1.0e3*rho0*vt/9.0

        temt = Max (temcg, -30.0)
        item = Int(-temt) + 1
!        qlwc = qcw*(1.e3)*rho0

        qlwc = Min ( qlwc, 30.0)
        ilwc = 0
        ftlwc = 0.0

        IF (qlwc .ge. 10.0) THEN
            ilwc = Int( qlwc/10.0 + 1.0e-5) + 1 + 27
            facl = (qlwc - lookup(0,ilwc) )*0.1
        ELSEIF (qlwc .ge. 1.0) THEN
            ilwc = Int( qlwc + 1.0e-5) + 1 + 18
            facl = (qlwc - lookup(0,ilwc) )*1.0
        ELSEIF (qlwc .ge. 0.1) THEN
            ilwc = Int( qlwc*10.0 + 1.0e-5 ) + 1 + 9
            facl = (qlwc - lookup(0,ilwc) )*10.0
        ELSEIF (qlwc .ge. 0.01) THEN
            ilwc = Int ( qlwc*100.0 + 1.0e-5) + 1
            facl = (qlwc - lookup(0,ilwc) )*100.0
        ELSE
            ilwc = 0
            ftlwc = 0.0
        END IF

        IF ( ilwc .eq. nlc ) THEN
            write(0,*) 'Warning: ilwc = nlc! setting to nlc-1'
            ilwc = nlc - 1
        ENDIF
        IF ( item .gt. nt ) THEN
            write(0,*) 'Warning: item > nt! setting to nt'
            item = nt
        ENDIF
        IF ( ilwc .lt. 0 ) THEN
            write(0,*) 'Warning: ilwc < 0! setting to 0'
            ilwc = 0
        ENDIF
        IF ( item .lt. 0 ) THEN
            write(0,*) 'Warning: item < 0! setting to 0'
            item = 0
        ENDIF

        IF (ilwc .gt. 0) THEN
            fact = -temt - Float(item) + 1.0
            a1 = (1.0-fact)*lookup(item,ilwc) +      &
                    &          fact*lookup(Min(nt,item+1),ilwc)
            a2 = (1.0-fact)*lookup(item,ilwc+1) +      &
                    &          fact*lookup(Min(nt,item+1),ilwc+1)
            ftlwc = fac*(a1 + facl*(a2-a1))

            alf = Min(5.0*(awdia/d0)**2*Abs(vt)/v0,10.0)
            ftlwc = 1.0e-15*alf*ftlwc
        END IF

        RETURN
    END subroutine takarar
    ! ####################################################################
    !                SUBROUTINE TAKA
    ! ####################################################################
    subroutine TAKA(nt,nlc,lookup,temcg,qcw,vt,awdia,rho0,ftlwc)

        ! ####################################################################
!
        !  Purpose:
        !           Noninductive charging per Takahashi 1978 lab results.
        !           Using lookup table from Wojcik (1994).  For temperatures
        !           lower than -30 C, the values at -30 C are used.
        !
        !  7/31/99 Now limit low-temp charging rates by a parabolic function
        !          with a value of zero at -30 C and goes to zero at -43 C.
        !
        !  8/10/99 Fixed bug in the interpolation: was accidentally
        !          _extrapolating_.   Oops.
        !
        ! ####################################################################

        implicit none

        integer nt,nlc
        real lookup(0:nt,0:nlc)
        real qlwc
        real qcw
        real rho0
        real ftlwc
        real temt,temcg
        integer item,ilwc
        real fact,facl
        real a1,a2,a3
        real alf
        real vt
        real awdia
        real d0
        parameter(d0 = 100.0e-6)
        real v0
        parameter(v0 = 8.0)
        real fac

        !  BEGIN EXECUTABLE CODE

        IF ( temcg .gt. -30 ) THEN
            fac = 1.0
        ELSEIF ( temcg .gt. -40.0 ) THEN
            fac = 1.0 - ((temcg+30.0)/(40.0 - 30.0))**2
        ELSE
            fac = 0.0
        END IF


        temt = Max (temcg, -30.0)
        item = Int(-temt) + 1
        qlwc = qcw*(1.e3)*rho0

        qlwc = Min ( qlwc, 30.0)
        ilwc = 0
        ftlwc = 0.0

        IF (qlwc .ge. 10.0) THEN
            ilwc = Int( qlwc/10.0 + 1.0e-5) + 1 + 27
            facl = (qlwc - lookup(0,ilwc) )*0.1
        ELSEIF (qlwc .ge. 1.0) THEN
            ilwc = Int( qlwc + 1.0e-5) + 1 + 18
            facl = (qlwc - lookup(0,ilwc) )*1.0
        ELSEIF (qlwc .ge. 0.1) THEN
            ilwc = Int( qlwc*10.0 + 1.0e-5 ) + 1 + 9
            facl = (qlwc - lookup(0,ilwc) )*10.0
        ELSEIF (qlwc .ge. 0.01) THEN
            ilwc = Int ( qlwc*100.0 + 1.0e-5) + 1
            facl = (qlwc - lookup(0,ilwc) )*100.0
        ELSE
            ilwc = 0
            ftlwc = 0.0
        END IF

        IF ( ilwc .eq. nlc ) THEN
            write(0,*) 'Warning: ilwc = nlc! setting to nlc-1'
            ilwc = nlc - 1
        ENDIF
        IF ( item .gt. nt ) THEN
            write(0,*) 'Warning: item > nt! setting to nt'
            item = nt
        ENDIF
        IF ( ilwc .lt. 0 ) THEN
            write(0,*) 'Warning: ilwc < 0! setting to 0'
            ilwc = 0
        ENDIF
        IF ( item .lt. 0 ) THEN
            write(0,*) 'Warning: item < 0! setting to 0'
            item = 0
        ENDIF

        IF (ilwc .gt. 0) THEN
            fact = -temt - Float(item) + 1.0
            a1 = (1.0-fact)*lookup(item,ilwc) +      &
                    &          fact*lookup(Min(nt,item+1),ilwc)
            a2 = (1.0-fact)*lookup(item,ilwc+1) +      &
                    &          fact*lookup(Min(nt,item+1),ilwc+1)
            ftlwc = fac*(a1 + facl*(a2-a1))

            alf = Min(5.0*(awdia/d0)**2*Abs(vt)/v0,10.0)
            ftlwc = 1.0e-15*alf*ftlwc
        END IF

        RETURN
    END subroutine TAKA

    ! ####################################################################
    !                SUBROUTINE TAKA3
    ! ####################################################################
    subroutine TAKA3(nt,nlc,lookup,temcg,qcw,vt,awdia,rho0,ftlwc)

        ! ####################################################################
!
        !  Purpose:
        !           Noninductive charging per Takahashi 1978 lab results.
        !           Using lookup table from Wojcik (1994).  For temperatures
        !           lower than -30 C, the values at -30 C are used.
!
        !  7/31/99 Now limit low-temp charging rates by a parabolic function
        !          with a value of zero at -30 C and goes to zero at -43 C.
!
        !  8/10/99 Fixed bug in the interpolation: was accidentally
        !          _extrapolating_.   Oops.
!
        ! ####################################################################

        implicit none

        integer nt,nlc
        real lookup(0:nt,0:nlc)
        real qlwc
        real qcw
        real rho0
        real ftlwc
        real temt,temcg
        integer item,ilwc
        real fact,facl
        real a1,a2,a3
        real alf
        real vt
        real awdia
        real d0
        parameter(d0 = 100.0e-6)
        real v0
        parameter(v0 = 8.0)
        real fac

        !  BEGIN EXECUTABLE CODE

        IF ( temcg .gt. -30 ) THEN
            fac = 1.0
        ELSEIF ( temcg .gt. -40.0 ) THEN
            fac = 1.0 - ((temcg+30.0)/(40.0 - 30.0))**2
        ELSE
            fac = 0.0
        END IF


        temt = Max (temcg, -30.0)
        item = Int(-temt) + 1
        qlwc = qcw*(1.e3)*rho0

        qlwc = Min ( qlwc, 30.0)
        ilwc = 0
        ftlwc = 0.0

        IF (qlwc .ge. 10.0) THEN
            ilwc = Int( qlwc/10.0 + 1.0e-5) + 1 + 27
            facl = (qlwc - lookup(0,ilwc) )*0.1
        ELSEIF (qlwc .ge. 1.0) THEN
            ilwc = Int( qlwc + 1.0e-5) + 1 + 18
            facl = (qlwc - lookup(0,ilwc) )*1.0
        ELSEIF (qlwc .ge. 0.1) THEN
            ilwc = Int( qlwc*10.0 + 1.0e-5 ) + 1 + 9
            facl = (qlwc - lookup(0,ilwc) )*10.0
        ELSEIF (qlwc .ge. 0.01) THEN
            ilwc = Int ( qlwc*100.0 + 1.0e-5) + 1
            facl = (qlwc - lookup(0,ilwc) )*100.0
        ELSE
            ilwc = 0
            ftlwc = 0.0
        END IF

        IF (ilwc .gt. 0) THEN
            fact = -temt - Float(item) + 1.0
            !        a1 = (1.0-fact)*lookup(item,ilwc) + fact*lookup(item+1,ilwc)
            !        a2 = (1.0-fact)*lookup(item,ilwc+1) + fact*lookup(item+1,ilwc+1)
            a1 = (1.0-fact)*lookup(item,ilwc) +      &
                    &          fact*lookup(Min(nt,item+1),ilwc)
            a2 = (1.0-fact)*lookup(item,ilwc+1) +      &
                    &          fact*lookup(Min(nt,item+1),ilwc+1)
            ftlwc = fac*(a1 + facl*(a2-a1))

!            alf = Min(5.0*(awdia/d0)**2*Abs(vt)/v0,10.0)
            ftlwc = 1.0e-15*ftlwc
        END IF

        RETURN
    END subroutine TAKA3

    ! ####################################################################
    !                SUBROUTINE TAKA4
    ! ####################################################################
    subroutine TAKA4(nt,nlc,lookup,temcg,qcw,vt,awdia,rho0,ftlwc)

        ! ####################################################################
        !
        !  Purpose:
        !           Noninductive charging per Takahashi 1978 lab results.
        !           Using lookup table from Wojcik (1994).  For temperatures
        !           lower than -30 C, the values at -30 C are used.
        !
        !  7/31/99 Now limit low-temp charging rates by a parabolic function
        !          with a value of zero at -30 C and goes to zero at -43 C.
        !
        !  8/10/99 Fixed bug in the interpolation: was accidentally
        !          _extrapolating_.   Oops.
        !
        ! ####################################################################

        implicit none

        integer nt,nlc
        real lookup(0:nt,0:nlc)
        real qlwc
        real qcw
        real rho0
        real ftlwc
        real temt,temcg
        integer item,ilwc
        real fact,facl
        real a1,a2,a3
        real alf
        real vt
        real awdia
        real d0
        parameter(d0 = 100.0e-6)
        real v0
        parameter(v0 = 8.0)
        real fac

        !  BEGIN EXECUTABLE CODE

        IF ( temcg .gt. -30 ) THEN
            fac = 1.0
        ELSEIF ( temcg .gt. -40.0 ) THEN
            fac = 1.0 - ((temcg+30.0)/(40.0 - 30.0))**2
        ELSE
            fac = 0.0
        END IF


        temt = Max (temcg, -30.0)
        item = Int(-temt) + 1
        qlwc = qcw*(1.e3)*rho0

        qlwc = Min ( qlwc, 30.0)
        ilwc = 0
        ftlwc = 0.0

        IF (qlwc .ge. 10.0) THEN
            ilwc = Int( qlwc/10.0 + 1.0e-5) + 1 + 27
            facl = (qlwc - lookup(0,ilwc) )*0.1
        ELSEIF (qlwc .ge. 1.0) THEN
            ilwc = Int( qlwc + 1.0e-5) + 1 + 18
            facl = (qlwc - lookup(0,ilwc) )*1.0
        ELSEIF (qlwc .ge. 0.1) THEN
            ilwc = Int( qlwc*10.0 + 1.0e-5 ) + 1 + 9
            facl = (qlwc - lookup(0,ilwc) )*10.0
        ELSEIF (qlwc .ge. 0.01) THEN
            ilwc = Int ( qlwc*100.0 + 1.0e-5) + 1
            facl = (qlwc - lookup(0,ilwc) )*100.0
        ELSE
            ilwc = 0
            ftlwc = 0.0
        END IF

        IF (ilwc .gt. 0) THEN
            fact = -temt - Float(item) + 1.0
            !        a1 = (1.0-fact)*lookup(item,ilwc) + fact*lookup(item+1,ilwc)
            !        a2 = (1.0-fact)*lookup(item,ilwc+1) + fact*lookup(item+1,ilwc+1)
            a1 = (1.0-fact)*lookup(item,ilwc) +      &
                    &          fact*lookup(Min(nt,item+1),ilwc)
            a2 = (1.0-fact)*lookup(item,ilwc+1) +      &
                    &          fact*lookup(Min(nt,item+1),ilwc+1)
            ftlwc = fac*(a1 + facl*(a2-a1))

            alf = Min(5.0*(3.67*awdia/d0)**2*Abs(vt)/v0,10.0)
            ftlwc = 1.0e-15*alf*ftlwc
        END IF

        RETURN
    END subroutine TAKA4
!
    ! ####################################################################
    ! ####################################################################
    !                SUBROUTINE SAUNDX
    ! ####################################################################
    subroutine saundx(isaund,temcg,qcw,exw,vt,awdia,rho0,qsign,ftrar,     &
            &                 qconkq,qconm,qconn,idelq,rarfac)

        ! ####################################################################
        !
        ! 05.04.2004  Fixed another error in the S91 code (in region S8 of
        !             Helsdon et al.)
        !             Added isaund = -2 option to replace anom. pos. zone with normal
        !              negative charging
        !             Added isaund = -3 option to use same equations as Wojcik 1994
!
        !  11/9/2001  Fixed a discrepancy between the crystal diameter used
        !             to calculate charge separation and the diameter used
        !             to find the size range.  Previously used awdia for finding
        !             the size range, but then used fac*awdia for calculating the
        !             charge separation.  (OK if fac = 1, but not for fac=3.67)
        !
        !
        !
        !  Purpose:
        !   Calls the appropriate version of the Saunders et al. scheme
        !
        !  -5 : Saunders etal 1991 (following Helsdon et al. 2001, but use normal charging instead of 'anomalous' zones)
        !  -4 : Saunders etal 1991 (following Helsdon et al. 2001, but remove positive 'anomalous' zone)
        !  -3 : Saunders etal 1991 (following Helsdon et al. 2001)
        !  -2 : Saunders etal 1991 (remove positive 'anomalous' zone)
        !  -1 : Saunders etal 1991 (unmodified 'anomalous' regions)
        !   0 : Saunders 1991 (modified as in Wojcik 1994)
        !   1 : RR scheme (with extra factor 3.67**qconm)
        !   2 : RR scheme ( no extra factor )
        !   3 : Saunders and Peck Scheme (with extra factor 3.67**qconm)
        !   4 : Saunders and Peck Scheme ( no extra factor )
        !   5 : Saunders and Peck Scheme (with extra factor 3.67**qconm, 0.5*rar)
        !   6 : Saunders and Peck Scheme (with extra factor 3.67**qconm, 0.75*rar)
        !   7 : Saunders and Peck Scheme ( no extra factor, 0.5*rar)
        !   8 : RR scheme ( no extra factor, , 0.5*rar )
        !   9 : Saunders and Peck Scheme ( no extra factor, cutoff at -32.47 as orig eq. from sp98 )
        !  10 : Brooks et al. RARcrit for T > -15 using saund2 (otherwise same as isaund=2) (set rarfac to negative in saund2)
        !  11 : Brooks et al. RARcrit for T > -15 using saund6 (otherwise same as isaund=4)
        !  12 : Brooks et al. RARcrit for T > -15 using saund6 (otherwise same as isaund=9)
        !  13 : Brooks et al. RARcrit for T > -15 using saund8
        ! ####################################################################

        implicit none
        integer isaund  ! scheme choice
        real temcg ! temperature
        real qcw ! cloud water mixing ratio
        real exw ! cloud water collection efficiency
        real vt ! terminal speed difference between x and cw
        real awdia ! crystal diameter
        real rho0  ! air density
        real qsign ! sign of charge acquired by rimer (not used by calling prog)
        real ftrar ! charge factor based on temp and RAR
        real qconkq ! factor kq ( or 'B')
        real qconm  ! exponent on crystal diameter ('a')
        real qconn  ! exponent on speed ('b')
        integer idelq ! charge sign, not used by calling program
        real awdia1
        real rarfac

        ! ####################################################################

        idelq = 0
        awdia1 = awdia

        IF ( isaund .eq. 1 ) THEN
            awdia1 = 3.67*awdia
            call saund2(temcg,qcw,exw,vt,awdia1,rho0,qsign,ftrar,     &
                    &                 qconkq,qconm,qconn,idelq,3.67,1.0)
        ELSEIF ( isaund .eq. 2 ) THEN
            call saund2(temcg,qcw,exw,vt,awdia,rho0,qsign,ftrar,     &
                    &                 qconkq,qconm,qconn,idelq,0.0,rarfac)
        ELSEIF ( isaund .eq. 8 ) THEN
            call saund2(temcg,qcw,exw,vt,awdia,rho0,qsign,ftrar,     &
                    &                 qconkq,qconm,qconn,idelq,0.0,0.5)
        ELSEIF ( isaund .eq. 3 ) THEN
            awdia1 = 3.67*awdia
            call saund6(temcg,qcw,exw,vt,awdia1,rho0,qsign,ftrar,     &
                    &                 qconkq,qconm,qconn,idelq,3.67,1.0)
        ELSEIF ( isaund .eq. 4 ) THEN
            call saund6(temcg,qcw,exw,vt,awdia,rho0,qsign,ftrar,     &
                    &                 qconkq,qconm,qconn,idelq,0.0,rarfac)
        ELSEIF ( isaund .eq. 5 ) THEN
            awdia1 = 3.67*awdia
            call saund6(temcg,qcw,exw,vt,awdia1,rho0,qsign,ftrar,     &
                    &                 qconkq,qconm,qconn,idelq,3.67,0.5)
        ELSEIF ( isaund .eq. 6 ) THEN
            awdia1 = 3.67*awdia
            call saund6(temcg,qcw,exw,vt,awdia1,rho0,qsign,ftrar,     &
                    &                 qconkq,qconm,qconn,idelq,3.67,0.75)
        ELSEIF ( isaund .eq. 7 ) THEN
            call saund6(temcg,qcw,exw,vt,awdia,rho0,qsign,ftrar,     &
                    &                 qconkq,qconm,qconn,idelq,0.0,0.5)
        ELSEIF ( isaund .eq. 9 ) THEN
            idelq = 1
            call saund6(temcg,qcw,exw,vt,awdia,rho0,qsign,ftrar,     &
                    &                 qconkq,qconm,qconn,idelq,0.0,rarfac)
        ELSEIF ( isaund .eq. 10 ) THEN
            call saund2(temcg,qcw,exw,vt,awdia,rho0,qsign,ftrar,     &
                    &                 qconkq,qconm,qconn,idelq,0.0,-rarfac)
        ELSEIF ( isaund .eq. 11 ) THEN
            call saund6(temcg,qcw,exw,vt,awdia,rho0,qsign,ftrar,     &
                    &                 qconkq,qconm,qconn,idelq,0.0,-rarfac)
        ELSEIF ( isaund .eq. 12 ) THEN ! default
            idelq = 1
            call saund6(temcg,qcw,exw,vt,awdia,rho0,qsign,ftrar,     &
                    &                 qconkq,qconm,qconn,idelq,0.0,-rarfac)
        ELSEIF ( isaund .eq. 13 ) THEN
            call saund8(temcg,qcw,exw,vt,awdia,rho0,qsign,ftrar,     &
                    &                 qconkq,qconm,qconn,idelq,0.0,-rarfac)
        ELSEIF ( isaund .eq. 0 ) THEN
            call saund(temcg,qcw,exw,awdia,rho0,qsign,ftrar,     &
                    &                 qconkq,qconm,qconn,idelq,0)
        ELSEIF ( isaund .eq. -1 ) THEN
            call saund(temcg,qcw,exw,awdia,rho0,qsign,ftrar,     &
                    &                 qconkq,qconm,qconn,idelq,1)
        ELSEIF ( isaund .eq. -2 ) THEN
            call saund(temcg,qcw,exw,awdia,rho0,qsign,ftrar,     &
                    &                 qconkq,qconm,qconn,idelq,2)
        ELSEIF ( isaund .eq. -3 ) THEN
            call saund(temcg,qcw,exw,awdia,rho0,qsign,ftrar,     &
                    &                 qconkq,qconm,qconn,idelq,3)
        ELSEIF ( isaund .eq. -4 ) THEN
            call saund(temcg,qcw,exw,awdia,rho0,qsign,ftrar,     &
                    &                 qconkq,qconm,qconn,idelq,4)
        ELSEIF ( isaund .eq. -5 ) THEN
            call saund(temcg,qcw,exw,awdia,rho0,qsign,ftrar,     &
                    &                 qconkq,qconm,qconn,idelq,5)
        ELSE
            write(6,*) 'STOP! unsupported value of isaund = ',isaund
            write (0,*) 'stopping'
            print*, 'stopping'
            STOP
        ENDIF
        RETURN
    END subroutine saundx

    ! ####################################################################
    !                SUBROUTINE SAUNDY
    ! ####################################################################
    subroutine saundy(isaund,temcg,qcw,rar,vt,awdia,rho0,qsign,ftrar,     &
            &                 qconkq,qconm,qconn,idelq)

        ! ####################################################################
        !
        !  4.1.2002   altered version of saundx to take rar as input instead
        !             of exw.
        !
        !  11/9/2001  Fixed a discrepancy between the crystal diameter used
        !             to calculate charge separation and the diameter used
        !             to find the size range.  Previously used awdia for finding
        !             the size range, but then used fac*awdia for calculating the
        !             charge separation.  (OK if fac = 1, but not for fac=3.67)
        !
        !
        !
        !  Purpose:
        !   Calls the appropriate version of the Saunders et al. scheme
        !
        !   1 : RAR scheme (with extra factor 3.67**qconm)
        !   2 : RAR scheme ( no extra factor )
        !   3 : Saunders and Peck Scheme (with extra factor 3.67**qconm)
        !   4 : Saunders and Peck Scheme ( no extra factor )
        !   5 : Saunders and Peck Scheme (with extra factor 3.67**qconm, 0.5*rar)
        !   6 : Saunders and Peck Scheme (with extra factor 3.67**qconm, 0.75*rar)
        !
        !
        ! ####################################################################

        implicit none
        integer isaund  ! scheme choice
        real temcg ! temperature
        real qcw ! cloud water mixing ratio
        real exw ! cloud water collection efficiency
        real rar
        real vt ! terminal speed difference between x and cw
        real awdia ! crystal diameter
        real rho0  ! air density
        real qsign ! sign of charge acquired by rimer (not used by calling prog)
        real ftrar ! charge factor based on temp and RAR
        real qconkq ! factor kq ( or 'B')
        real qconm  ! exponent on crystal diameter ('a')
        real qconn  ! exponent on speed ('b')
        integer idelq ! charge sign, not used by calling program
        real awdia1

        ! ####################################################################

        awdia1 = awdia

        IF ( isaund .eq. 4 ) THEN
            call saund7(temcg,qcw,rar,vt,awdia,rho0,qsign,ftrar,     &
                    &                 qconkq,qconm,qconn,idelq,0.0,1.0)
        ENDIF
        RETURN
    END subroutine saundy
    ! ####################################################################
    !                SUBROUTINE SAUND7
    ! ####################################################################
    subroutine saund7(temcg,qcw,rar,vt,awdia,rho0,qsign,ftrar,     &
            &                 qconkq,qconm,qconn,idelq,qfac,rarfac)

        ! ####################################################################
        !
        ! 4.1.2002  New version that uses the model-calculated RAR (derived
        !           from the collection rate qxacw)
        !
        !cPurpose:
        !   Saunders charging scheme based on rime accretion rate (RAR)
        !   as in Brooks et al. (1997) and Saunders and Peck (1998)
        !
        ! 7/27/2000 Try with original Saunders and Peck curve (1998) plus the
        !           drop-off starting at -23.7
!
!
        ! 7/24/99  new adjustment to RARcritical curve: make a cubic drop-off
        !          from -23.7 to -40 (and zero charging below -40)
!
        ! 12/8/99 temporary test of putting in a factor of (3.67)**qconm to
        !         approximate using the mass-weighted mean crystal diameter
!
        ! 12/17/99 for positive charging, just use ftrar = 6.74*( rar - rarc )
        !          for all temperatures (took out alternate eq. for T<-34)
        !
        ! ####################################################################

        implicit none
        real temcg ! temperature
        real qcw ! cloud water mixing ratio
        !      real exw ! cloud water collection efficiency
        real vt ! terminal speed difference between x and cw
        real awdia ! crystal diameter
        real rho0  ! air density
        real qsign ! sign of charge acquired by rimer (not used by calling prog)
        real ftrar ! charge factor based on temp and RAR
        real qconkq ! factor kq ( or 'B')
        real qconm  ! exponent on crystal diameter ('a')
        real qconn  ! exponent on speed ('b')
        integer idelq ! charge sign, not used by calling program
        real qfac
        real rarfac  ! factor to reduce rar for purposes of charge calculation

        real rar ! rime accretion rate
        real rarc ! critical RAR
        real t,tc
        parameter (tc = -23.7)
        real delr,delri
        parameter (delr=0.5, delri=1.0/delr)
        real rar0  ! lower limit of RAR where charging goes to zero
        parameter (rar0 = 0.1)
        real q0
        parameter (q0 = 6.48)
        real tema,tmin
        parameter (tema = -7.0 , tmin = -43.)

        real fac

        ! ####################################################################
        ! Begin Executable code
        ! ####################################################################

        IF ( temcg .gt. -30 ) THEN
            fac = 1.0
        ELSEIF ( temcg .gt. -40.0 ) THEN
            fac = 1.0 - ((temcg+30.0)/(40.0 - 30.0))**2
        ELSE
            fac = 0.0
        END IF

        ftrar = 0.0

        !      rar = exw*qcw*1.0e3*rho0*vt*rarfac
        t = temcg

        IF (t .gt. tc ) THEN
            rarc = 1.0 + t*(7.9262e-2 + t*(4.4847e-2 +     &
                    &  t*(7.4754e-3 + t*(5.4686e-4 +     &
                            &  t*(1.6737e-5 + t*1.7613e-7)))))
            !      ELSE ! IF ( t .le. -23.0) THEN
            !       rarc = 3.27
            !      ELSE ! IF ( t .le. -30.0) THEN
            !       rarc = 1.795
        ELSE ! IF ( t .le. tc) THEN
            rarc = 3.39608*( 1.0 - Abs( ( (t - tc)/(tc + 40) )**3) )
        END IF



        ! check for RAR below threshold:
        IF ( rar .le. rar0 .or. temcg .le. tmin) THEN
            qsign = 0.0
            ftrar = 0.0
            qconkq = 0.0
            qconn = 1.0
            qconm = 1.0
            GOTO 999
        END IF

        ! for now, use charging values at -8.0 for warmer temps
        ! and values at -23.0 for colder temps (Saunders and Peck, 1998)
        IF (temcg .gt. -8.0 ) THEN
            t = -8.0
        ELSEIF (temcg .lt. -23.0) THEN
            t = -23.0
        END IF

        qsign = -1.0
        IF (rar .gt. rarc) qsign = 1.0

        IF (qsign .gt. 0.5) THEN

            !        IF ( Abs(rar - rarc) .lt. delr) THEN
            !         ftrar = delri*Abs(rar - rarc)*
            !     :           (6.74*(rarc + delr) + 1.36*t + 10.5 )
!        ELSE
            !         ftrar = 6.74*rar + 1.36*t + 10.5
!        END IF
            !        ftrar = Max(0.0, ftrar)

            !        IF ( Abs(rar - rarc) .lt. delr .and. temcg .gt. tc) THEN
            !         ftrar = Max( 0.0, delri*Abs(rar - rarc)*
            !     :           (6.74*(rarc + delr) + 1.36*t + 10.5 ) )
!        ELSE
            !         ftrar = Max( 0.0, 6.74*rar + 1.36*t + 10.5 )
!        ENDIF

            !        IF ( temcg .gt. -34.0 ) THEN
            ftrar = 6.74*( rar - rarc )
!        ELSE
            !          ftrar = Max ( 0.0, 6.74*( rar -
            !     :     4*( 1.0 - Abs( ( (-34.0 + 25.0 )/(-25.0 - tmin) )**2) ) ) )
!        END IF

            if ( awdia*1.e6.lt.155. ) then
                qconkq = 4.9e13
                qconm  = 3.76
                qconn  = 2.5
            end if
            if ( awdia*1.e6.ge.155. .and. awdia*1.e6.le.452. ) then
                qconkq = 4.0e6
                qconm  = 1.9
                qconn  = 2.5
            end if
            if ( awdia*1.e6.gt.452. ) then
                qconkq = 52.8
                qconm  = 0.44
                qconn  = 2.5
            end if

            IF ( qfac .gt. 1.0 ) THEN
                ftrar = fac*ftrar*qfac**qconm
            ELSE
                ftrar = fac*ftrar
            ENDIF

        ELSEIF (qsign .lt. -0.5) THEN

            !        IF ( Abs(rar - rarc) .lt. delr) THEN
            !         ftrar = delri*Abs(rar - rarc)*
            !     :           (3.02 - 10.59*(rarc-delr) + 2.95*(rarc-delr)**2 )
            !        ELSEIF (rar .lt. 0.4) THEN
            !         ftrar = ((rar-0.1)/0.3)*(3.02 - 10.59*(0.4) + 2.95*(0.4)**2)
!        ELSE
            !         ftrar = 3.02 - 10.59*rar + 2.95*rar**2
!        END IF

            ftrar = Min(1.0, 0.6*Abs(rarc - rar0) )*(6.5)*(-1.0 +      &
                    &    4.0/(rarc - rar0)**2 * (rar - (rarc + rar0)/2.0 )**2  )
            ftrar = Min( 0.0, ftrar )

            if ( awdia*1.e6.lt.253. ) then
                qconkq = 5.24e8
                qconm  = 2.54
                qconn  = 2.8
            end if
            if ( awdia*1.e6.gt.253. ) then
                qconkq = 24.0
                qconm  = 0.50
                qconn  = 2.8
            end if

            IF ( qfac .gt. 1.0 ) THEN
                ftrar = fac*ftrar*qfac**qconm
            ELSE
                ftrar = fac*ftrar
            ENDIF

        END IF


        999  CONTINUE
        RETURN
    END subroutine saund7

    ! ####################################################################
    !                SUBROUTINE SAUND6
    ! ####################################################################
    subroutine saund6(temcg,qcw,exw,vt,awdia,rho0,qsign,ftrar,     &
            &                 qconkq,qconm,qconn,idelq,qfac,rarfac)

        ! ####################################################################
!
        !  Purpose:
        !   Saunders charging scheme based on rime accretion rate (RAR)
        !   as in Brooks et al. (1997) and Saunders and Peck (1998)
        !
        ! 10/12/2003:  Add line 'rarc = Max(rarc, 0.0)' to prevent negative
        !              values of rarc for -43 < T < -40
        !
        !
        ! 7/27/2000 Try with original Saunders and Peck curve (1998) plus the
        !           drop-off starting at -23.7
        !
        !
        ! 7/24/99  new adjustment to RARcritical curve: make a cubic drop-off
        !          from -23.7 to -40 (and zero charging below -40)
        !         (correction, 10/12/2003: charging was set zero below -43, not -40)
        !
        ! 12/8/99 temporary test of putting in a factor of (3.67)**qconm to
        !         approximate using the mass-weighted mean crystal diameter
        !
        ! 12/17/99 for positive charging, just use ftrar = 6.74*( rar - rarc )
        !          for all temperatures (took out alternate eq. for T<-34)
        !
        ! ####################################################################

        implicit none
        real temcg ! temperature
        real qcw ! cloud water mixing ratio
        real exw ! cloud water collection efficiency
        real vt ! terminal speed difference between x and cw
        real awdia ! crystal diameter
        real rho0  ! air density
        real qsign ! sign of charge acquired by rimer (not used by calling prog)
        real ftrar ! charge factor based on temp and RAR
        real qconkq ! factor kq ( or 'B')
        real qconm  ! exponent on crystal diameter ('a')
        real qconn  ! exponent on speed ('b')
        integer idelq ! charge sign, not used by calling program
        real qfac
        real rarfac  ! factor to reduce rar for purposes of charge calculation

        real rar ! rime accretion rate
        real rarc ! critical RAR
        real t,tc
        parameter (tc = -23.7)
        real delr,delri
        parameter (delr=0.5, delri=1.0/delr)
        real rar0  ! lower limit of RAR where charging goes to zero
        parameter (rar0 = 0.1)
        real q0
        parameter (q0 = 6.48)
        real tema,tmin
        parameter (tema = -7.0 , tmin = -37.0)
        integer ibs

        real fac

        ! ####################################################################
        ! Begin Executable code
        ! ####################################################################

        IF ( temcg .gt. -30. ) THEN
            fac = 1.0
!        ELSEIF ( temcg .gt. -40.0 ) THEN
!            fac = 1.0 - ((temcg+30.0)/(40.0 - 30.0))**2
        ELSEIF ( temcg .gt. tmin ) THEN
            fac = 1.0 - ((temcg+30.0)/(-tmin - 30.0))**2
        ELSE
            fac = 0.0
        END IF

        IF ( rarfac .gt. 0.0 ) THEN
            ibs = 0
        ELSE
            ibs = 1
        ENDIF

        ftrar = 0.0

        rar = exw*qcw*1.0e3*rho0*vt*Abs(rarfac)
        t = temcg

        IF ( ibs .eq. 1 .and. t .gt. -15.0 ) THEN
            rarc = Max( 0.0, Min( 3.29, -1.47 - 0.2*t ) )
        ELSE
            IF (t .gt. tc ) THEN
                rarc = 1.0 + t*(7.9262e-2 + t*(4.4847e-2 +     &
                        &    t*(7.4754e-3 + t*(5.4686e-4 +     &
                                &    t*(1.6737e-5 + t*1.7613e-7)))))
                !      ELSE ! IF ( t .le. -23.0) THEN
                !       rarc = 3.27
                !      ELSE ! IF ( t .le. -30.0) THEN
                !       rarc = 1.795
            ELSE ! IF ( t .le. tc) THEN
                rarc = 3.39608*( 1.0 - Abs( ( (t - tc)/(tc - tmin) )**3) )
            ENDIF
        ENDIF

        rarc = Max(rarc, rar0)  ! ERM 10/12/2003; changed 0 to rar0 5/8/2004

!
        ! New option (2/9/2005) to use original SP98 equation and cut off charging where
        !  it hits zero
!
        IF ( idelq .eq. 1 ) THEN

            IF ( t .le. tc) THEN
                rarc = Max(0.0, 1.0 + t*(7.9262e-2 + t*(4.4847e-2 +     &
                        &    t*(7.4754e-3 + t*(5.4686e-4 +     &
                                &    t*(1.6737e-5 + t*1.7613e-7))))) )
            END IF

            IF ( t .lt. -33. ) rarc = rar0

            IF ( temcg .gt. -25.0 ) THEN
                fac = 1.0
            ELSEIF ( temcg .gt. -32.47 ) THEN
                fac = 1.0 + ((temcg+25.0)/(32.47 - 25.00))
            ELSE
                fac = 0.0
                qsign = 0.0
                ftrar = 0.0
                qconkq = 0.0
                qconn = 1.0
                qconm = 1.0
                GOTO 999
            END IF

        ENDIF


        ! check for RAR below threshold:
        IF ( rar .le. rar0 .or. temcg .le. tmin) THEN
            qsign = 0.0
            ftrar = 0.0
            qconkq = 0.0
            qconn = 1.0
            qconm = 1.0
            GOTO 999
        END IF

        ! for now, use charging values at -8.0 for warmer temps
        ! and values at -23.0 for colder temps (Saunders and Peck, 1998)
        IF (temcg .gt. -8.0 ) THEN
            t = -8.0
        ELSEIF (temcg .lt. -23.0) THEN
            t = -23.0
        END IF

        qsign = -1.0
        IF (rar .gt. rarc) qsign = 1.0

        IF (qsign .gt. 0.5) THEN

            ftrar = 6.74*( rar - rarc )

            if ( awdia*1.e6.lt.155. ) then
                qconkq = 4.9e13
                qconm  = 3.76
                qconn  = 2.5
            end if
            if ( awdia*1.e6.ge.155. .and. awdia*1.e6.le.452. ) then
                qconkq = 4.0e6
                qconm  = 1.9
                qconn  = 2.5
            end if
            if ( awdia*1.e6.gt.452. ) then
                qconkq = 52.8
                qconm  = 0.44
                qconn  = 2.5
            end if

            IF ( qfac .gt. 1.0 ) THEN
                ftrar = fac*ftrar*qfac**qconm
            ELSE
                ftrar = fac*ftrar
            ENDIF

        ELSEIF (qsign .lt. -0.5) THEN


            ftrar = Min(1.0, 0.6*Abs(rarc - rar0) )*(6.5)*(-1.0 +      &
                    &    4.0/(rarc - rar0)**2 * (rar - (rarc + rar0)/2.0 )**2  )
            ftrar = Min( 0.0, ftrar )

            if ( awdia*1.e6.lt.253. ) then
                qconkq = 5.24e8
                qconm  = 2.54
                qconn  = 2.8
            end if
            if ( awdia*1.e6.gt.253. ) then
                qconkq = 24.0
                qconm  = 0.50
                qconn  = 2.8
            end if

            IF ( qfac .gt. 1.0 ) THEN
                ftrar = fac*ftrar*qfac**qconm
            ELSE
                ftrar = fac*ftrar
            ENDIF

        END IF


        999  CONTINUE
        RETURN
    END subroutine saund6


    ! ####################################################################
    !                SUBROUTINE SAUND8
    ! ####################################################################
    subroutine saund8(temcg,qcw,exw,vt,awdia,rho0,qsign,ftrar,     &
            &                 qconkq,qconm,qconn,idelq,qfac,rarfac)

        ! ####################################################################
!
        !  Purpose:
        !   Saunders charging scheme based on rime accretion rate (RAR)
        !   as in Brooks et al. (1997) and Saunders and Peck (1998)
!
        ! 10/12/2003:  Add line 'rarc = Max(rarc, 0.0)' to prevent negative
        !              values of rarc for -43 < T < -40
!
!
        ! 7/27/2000 Try with original Saunders and Peck curve (1998) plus the
        !           drop-off starting at -23.7
!
!
        ! 7/24/99  new adjustment to RARcritical curve: make a cubic drop-off
        !          from -23.7 to -40 (and zero charging below -40)
        !         (correction, 10/12/2003: charging was set zero below -43, not -40)
!
        ! 12/8/99 temporary test of putting in a factor of (3.67)**qconm to
        !         approximate using the mass-weighted mean crystal diameter
!
        ! 12/17/99 for positive charging, just use ftrar = 6.74*( rar - rarc )
        !          for all temperatures (took out alternate eq. for T<-34)
!
        ! ####################################################################

        implicit none
        real temcg ! temperature
        real qcw ! cloud water mixing ratio
        real exw ! cloud water collection efficiency
        real vt ! terminal speed difference between x and cw
        real awdia ! crystal diameter
        real rho0  ! air density
        real qsign ! sign of charge acquired by rimer (not used by calling prog)
        real ftrar ! charge factor based on temp and RAR
        real qconkq ! factor kq ( or 'B')
        real qconm  ! exponent on crystal diameter ('a')
        real qconn  ! exponent on speed ('b')
        integer idelq ! charge sign, not used by calling program
        real qfac
        real rarfac  ! factor to reduce rar for purposes of charge calculation

        real rar ! rime accretion rate
        real rarc ! critical RAR
        real t,tc
!        parameter (tc = -23.7)
        parameter (tc = -24.5)
        real delr,delri
        parameter (delr=0.5, delri=1.0/delr)
        real rar0  ! lower limit of RAR where charging goes to zero
        parameter (rar0 = 0.1)
        real q0
        parameter (q0 = 6.48)
        real tema,tmin,tcc
        parameter (tema = -7.0 , tmin = -37.0)
        real :: tmincutoff
        integer ibs
        real, parameter :: pi = 3.141592654

        real fac

        ! ####################################################################
        ! Begin Executable code
        ! ####################################################################

        tcc = -30.
        tmincutoff = -37.
        IF ( temcg .gt. tcc ) THEN
            fac = 1.0
            !      ELSEIF ( temcg .gt. -40.0 ) THEN
            !        fac = 1.0 - ((temcg+30.0)/(40.0 - 30.0))**2
        ELSEIF ( temcg .gt. tmincutoff ) THEN
            !        fac = 1.0 - ((temcg+30.0)/(-tmin - 30.0))**2
            ! Cosine funtion roll-off
            fac = 0.5 *(1. + Cos(pi *(((temcg - tcc)/(-tmincutoff + tcc)))))
        ELSE
            fac = 0.0
        END IF

        IF ( rarfac .gt. 0.0 ) THEN
            ibs = 0
        ELSE
            ibs = 1
        ENDIF

        ftrar = 0.0

        rar = exw*qcw*1.0e3*rho0*vt*Abs(rarfac)
        t = temcg

        IF ( ibs .eq. 1 .and. t .gt. -15.0 ) THEN
            rarc = Max( 0.0, Min( 3.29, -1.47 - 0.2*t ) )
        ELSE
            IF (t .gt. tc ) THEN
                rarc = 1.0 + t*(7.9262e-2 + t*(4.4847e-2 +     &
                        &    t*(7.4754e-3 + t*(5.4686e-4 +     &
                                &    t*(1.6737e-5 + t*1.7613e-7)))))
                !      ELSE ! IF ( t .le. -23.0) THEN
                !       rarc = 3.27
                !      ELSE ! IF ( t .le. -30.0) THEN
                !       rarc = 1.795
            ELSEIF ( t > tmin) THEN
                !         rarc = 3.39608*( 1.0 - Abs( ( (t - tc)/(tc - tmin) )**3) )
                rarc = 3.42533* 0.5 *(1. + Cos(pi *(((t - tc)/(-tmin + tc))))) ! 3.42533 is the value on the curve at T = -24.5284
            ELSE
                rarc = 0
            ENDIF
        ENDIF

        rarc = Max(rarc, rar0)  ! ERM 10/12/2003; changed 0 to rar0 5/8/2004

        !
        ! New option (2/9/2005) to use original SP98 equation and cut off charging where
        !  it hits zero
        !
        IF ( idelq .eq. 1 ) THEN

            IF ( t .le. tc) THEN
                rarc = Max(0.0, 1.0 + t*(7.9262e-2 + t*(4.4847e-2 +     &
                        &    t*(7.4754e-3 + t*(5.4686e-4 +     &
                                &    t*(1.6737e-5 + t*1.7613e-7))))) )
            END IF

            IF ( t .lt. -33. ) rarc = rar0

            IF ( temcg .gt. -25.0 ) THEN
                fac = 1.0
            ELSEIF ( temcg .gt. -32.47 ) THEN
                fac = 1.0 + ((temcg+25.0)/(32.47 - 25.00))
            ELSE
                fac = 0.0
                qsign = 0.0
                ftrar = 0.0
                qconkq = 0.0
                qconn = 1.0
                qconm = 1.0
                GOTO 999
            END IF

        ENDIF


        ! check for RAR below threshold:
        IF ( rar .le. rar0 .or. temcg .le. tmin) THEN
            qsign = 0.0
            ftrar = 0.0
            qconkq = 0.0
            qconn = 1.0
            qconm = 1.0
            GOTO 999
        END IF

        ! for now, use charging values at -8.0 for warmer temps
        ! and values at -23.0 for colder temps (Saunders and Peck, 1998)
        IF (temcg .gt. -8.0 ) THEN
            t = -8.0
        ELSEIF (temcg .lt. -23.0) THEN
            t = -23.0
        END IF

        qsign = -1.0
        IF (rar .gt. rarc) qsign = 1.0

        IF (qsign .gt. 0.5) THEN

            ftrar = 6.74*( rar - rarc )

            if ( awdia*1.e6.lt.155. ) then
                qconkq = 4.9e13
                qconm  = 3.76
                qconn  = 2.5
            end if
            if ( awdia*1.e6.ge.155. .and. awdia*1.e6.le.452. ) then
                qconkq = 4.0e6
                qconm  = 1.9
                qconn  = 2.5
            end if
            if ( awdia*1.e6.gt.452. ) then
                qconkq = 52.8
                qconm  = 0.44
                qconn  = 2.5
            end if

            IF ( qfac .gt. 1.0 ) THEN
                ftrar = fac*ftrar*qfac**qconm
            ELSE
                ftrar = fac*ftrar
            ENDIF

        ELSEIF (qsign .lt. -0.5) THEN


            ftrar = Min(1.0, 0.6*Abs(rarc - rar0) )*(6.5)*(-1.0 +      &
                    &    4.0/(rarc - rar0)**2 * (rar - (rarc + rar0)/2.0 )**2  )
            ftrar = Min( 0.0, ftrar )

            if ( awdia*1.e6.lt.253. ) then
                qconkq = 5.24e8
                qconm  = 2.54
                qconn  = 2.8
            end if
            if ( awdia*1.e6.gt.253. ) then
                qconkq = 24.0
                qconm  = 0.50
                qconn  = 2.8
            end if

            IF ( qfac .gt. 1.0 ) THEN
                ftrar = fac*ftrar*qfac**qconm
            ELSE
                ftrar = fac*ftrar
            ENDIF

        END IF


        999  CONTINUE
        RETURN
    END subroutine saund8

          ! ####################################################################
          !                SUBROUTINE SAUND2
          ! ####################################################################
          subroutine saund2(temcg,qcw,exw,vt,awdia,rho0,qsign,ftrar,     &
                  &                 qconkq,qconm,qconn,idelq,qfac,rarfac)

              ! ####################################################################
              !
              !  Purpose:
              !   Saunders charging scheme based on rime accretion rate (RAR)
              !   as in Brooks et al. (1997) and Saunders and Peck (1998)
              !
              ! 8.27.2002 added rarfac input
              !
              ! 7/24/99  new adjustment to RARcritical curve: make a cubic drop-off
              !          from -23.7 to -40 (and zero charging below -40)
              !  version for small droplet curve (Saunders, ICAE 1999)
              ! 8/15/99 use higher RARcrit at low temp and keep small droplet
              !
              ! 12/8/99 temporary test of putting in a factor of (3.67)**qconm to
              !         approximate using the mass-weighted mean crystal diameter
              !
              ! 12/17/99 for positive charging, just use ftrar = 6.74*( rar - rarc )
              !          for all temperatures (took out alternate eq. for T<-34)
              !
              ! ####################################################################

              implicit none
              real temcg ! temperature
              real qcw ! cloud water mixing ratio
              real exw ! cloud water collection efficiency
              real vt ! terminal speed difference between x and cw
              real awdia ! crystal diameter
              real rho0  ! air density
              real qsign ! sign of charge acquired by rimer (not used by calling prog)
              real ftrar ! charge factor based on temp and RAR
              real qconkq ! factor kq ( or 'B')
              real qconm  ! exponent on crystal diameter ('a')
              real qconn  ! exponent on speed ('b')
              integer idelq ! charge sign, not used by calling program
              real qfac

              real rar ! rime accretion rate
              real rarc ! critical RAR
              real t,tc
              parameter (tc = -23.7)
              real delr,delri
              parameter (delr=0.5, delri=1.0/delr)
              real rar0  ! lower limit of RAR where charging goes to zero
              parameter (rar0 = 0.1)
              real q0
              parameter (q0 = 6.48)
              real tema,tmin
              parameter (tema = -7.0 , tmin = -43.)
              integer ibs

              real fac,rarfac

              ! ####################################################################
              ! Begin Executable code
              ! ####################################################################

              IF ( temcg .gt. -30 ) THEN
                  fac = 1.0
              ELSEIF ( temcg .gt. -40.0 ) THEN
                  fac = 1.0 - ((temcg+30.0)/(40.0 - 30.0))**2
              ELSE
                  fac = 0.0
              END IF

              IF ( rarfac .gt. 0.0 ) THEN
                  ibs = 0
              ELSE
                  ibs = 1
              ENDIF

              ftrar = 0.0

              rar = exw*qcw*1.0e3*rho0*vt*Abs(rarfac)
              t = temcg
              IF ( ibs .eq. 1 .and. t .gt. -15.0 ) THEN
                  rarc = Max( rar0, Min( 3.29, -1.47 - 0.2*t ) )
              ELSE
                  IF (t .gt. -7.0 ) THEN
                      rarc = 1.0 + t*(7.9262e-2 + t*(4.4847e-2 +     &
                              &    t*(7.4754e-3 + t*(5.4686e-4 +     &
                                      &    t*(1.6737e-5 + t*1.7613e-7)))))
                  ELSEIF ( t .gt. -16.0 ) THEN
                      rarc =  8.0*Abs((tema-t)/(10.0+tema))*     &
                              &      Exp(- Abs((tema-t)/(10.0+tema))) +     &
                              &     1.0 + t*(7.9262e-2 + t*(4.4847e-2 +     &
                              &    t*(7.4754e-3 + t*(5.4686e-4 +     &
                                      &    t*(1.6737e-5 + t*1.7613e-7)))))
                  ELSEIF ( t .gt. -21.7 ) THEN
                      rarc = 4*( 1.0 - Abs( ( (t + 25.0 )/(-25.0 - tmin) )**2) )
                  ELSE
                      rarc = 4*( 1.0 - Abs( ( (-21.7 + 25.0 )/(-25.0 - tmin) )**2) )
                  ENDIF
              ENDIF ! ibs .eq. 1

              ! check for RAR below threshold:
              IF ( rar .le. rar0 .or. temcg .le. tmin) THEN
                  qsign = 0.0
                  ftrar = 0.0
                  qconkq = 0.0
                  qconn = 1.0
                  qconm = 1.0
                  GOTO 999
              END IF

              ! for now, use charging values at -8.0 for warmer temps
              ! and values at -23.0 for colder temps (Saunder and Peck, 1998)
              IF (temcg .gt. -8.0 ) THEN
                  t = -8.0
              ELSEIF (temcg .lt. -23.0) THEN
                  t = -23.0
              END IF

              qsign = -1.0
              IF (rar .gt. rarc) qsign = 1.0

              IF (qsign .gt. 0.5) THEN

                  !        IF ( Abs(rar - rarc) .lt. delr) THEN
                  !         ftrar = delri*Abs(rar - rarc)*
                  !     :           (6.74*(rarc + delr) + 1.36*t + 10.5 )
                  !        ELSE
                  !         ftrar = 6.74*rar + 1.36*t + 10.5
                  !        END IF
                  !        ftrar = Max(0.0, ftrar)

                  !        IF ( Abs(rar - rarc) .lt. delr .and. temcg .gt. tc) THEN
                  !         ftrar = Max( 0.0, delri*Abs(rar - rarc)*
                  !     :           (6.74*(rarc + delr) + 1.36*t + 10.5 ) )
                  !        ELSE
                  !         ftrar = Max( 0.0, 6.74*rar + 1.36*t + 10.5 )
                  !        END IF

                  !        IF ( temcg .gt. -34.0 ) THEN
                  ftrar = 6.74*( rar - rarc )
                  !        ELSE
                  !          ftrar = Max ( 0.0, 6.74*( rar -
                  !     :     4*( 1.0 - Abs( ( (-34.0 + 25.0 )/(-25.0 - tmin) )**2) ) ) )
                  !        END IF

                  if ( awdia*1.e6.lt.155. ) then
                      qconkq = 4.9e13
                      qconm  = 3.76
                      qconn  = 2.5
                  end if
                  if ( awdia*1.e6.ge.155. .and. awdia*1.e6.le.452. ) then
                      qconkq = 4.0e6
                      qconm  = 1.9
                      qconn  = 2.5
                  end if
                  if ( awdia*1.e6.gt.452. ) then
                      qconkq = 52.8
                      qconm  = 0.44
                      qconn  = 2.5
                  end if

                  IF ( qfac .gt. 1.0 ) THEN
                      ftrar = fac*ftrar*qfac**qconm
                  ELSE
                      ftrar = fac*ftrar
                  ENDIF

              ELSEIF (qsign .lt. -0.5) THEN

                  !        IF ( Abs(rar - rarc) .lt. delr) THEN
                  !         ftrar = delri*Abs(rar - rarc)*
                  !     :           (3.02 - 10.59*(rarc-delr) + 2.95*(rarc-delr)**2 )
                  !        ELSEIF (rar .lt. 0.4) THEN
                  !         ftrar = ((rar-0.1)/0.3)*(3.02 - 10.59*(0.4) + 2.95*(0.4)**2)
                  !        ELSE
                  !         ftrar = 3.02 - 10.59*rar + 2.95*rar**2
                  !        END IF

                  ftrar = Min(1.0, 0.6*Abs(rarc - rar0) )*(6.5)*(-1.0 +      &
                          &    4.0/(rarc - rar0)**2 * (rar - (rarc + rar0)/2.0 )**2  )
                  ftrar = Min( 0.0, ftrar )

                  if ( awdia*1.e6.lt.253. ) then
                      qconkq = 5.24e8
                      qconm  = 2.54
                      qconn  = 2.8
                  end if
                  if ( awdia*1.e6.gt.253. ) then
                      qconkq = 24.0
                      qconm  = 0.50
                      qconn  = 2.8
                  end if

                  IF ( qfac .gt. 1.0 ) THEN
                      ftrar = fac*ftrar*qfac**qconm
                  ELSE
                      ftrar = fac*ftrar
                  ENDIF

              END IF


              999  CONTINUE
              RETURN
          END subroutine saund2


          ! ####################################################################
          !                SUBROUTINE SAUND5
          ! ####################################################################
          subroutine saund5(temcg,qcw,exw,vt,awdia,rho0,qsign,ftrar,     &
                  &                 qconkq,qconm,qconn,idelq)

              ! ####################################################################
              !
              !  Purpose:
              !   Saunders charging scheme based on rime accretion rate (RAR)
              !   as in Brooks et al. (1997) and Saunders and Peck (1998)
              !
              ! 7/24/99  new adjustment to RARcritical curve: make a cubic drop-off
              !          from -23.7 to -40 (and zero charging below -40)
              !  version for small droplet curve (Saunders, ICAE 1999)
              ! ####################################################################

              implicit none
              real temcg ! temperature
              real qcw ! cloud water mixing ratio
              real exw ! cloud water collection efficiency
              real vt ! terminal speed difference between x and cw
              real awdia ! crystal diameter
              real rho0  ! air density
              real qsign ! sign of charge acquired by rimer (not used by calling prog)
              real ftrar ! charge factor based on temp and RAR
              real qconkq ! factor kq ( or 'B')
              real qconm  ! exponent on crystal diameter ('a')
              real qconn  ! exponent on speed ('b')
              integer idelq ! charge sign, not used by calling program

              real rar ! rime accretion rate
              real rarc ! critical RAR
              real t,tc
              parameter (tc = -23.7)
              real delr,delri
              parameter (delr=0.5, delri=1.0/delr)
              real rar0  ! lower limit of RAR where charging goes to zero
              parameter (rar0 = 0.1)
              real q0
              parameter (q0 = 6.48)
              real tema,tmin
              parameter (tema = -7.0 , tmin = -43.)

              ! ####################################################################
              ! Begin Executable code
              ! ####################################################################

              ftrar = 0.0

              rar = exw*qcw*1.0e3*rho0*vt
              t = temcg
              IF (t .gt. -7.0 ) THEN
                  rarc = 1.0 + t*(7.9262e-2 + t*(4.4847e-2 +     &
                          &  t*(7.4754e-3 + t*(5.4686e-4 +     &
                                  &  t*(1.6737e-5 + t*1.7613e-7)))))
              ELSEIF ( t .gt. -16.0 ) THEN
                  rarc =  8.0*Abs((tema-t)/(10.0+tema))*     &
                          &      Exp(- Abs((tema-t)/(10.0+tema))) +     &
                          &     1.0 + t*(7.9262e-2 + t*(4.4847e-2 +     &
                          &  t*(7.4754e-3 + t*(5.4686e-4 +     &
                                  &  t*(1.6737e-5 + t*1.7613e-7)))))
              ELSE
                  rarc = 4*( 1.0 - Abs( ( (t + 25.0 )/(-25.0 - tmin) )**2) )
              END IF

              ! check for RAR below threshold:
              IF ( rar .le. rar0 .or. temcg .le. tmin) THEN
                  qsign = 0.0
                  ftrar = 0.0
                  qconkq = 0.0
                  qconn = 1.0
                  qconm = 1.0
                  GOTO 999
              END IF

              ! for now, use charging values at -8.0 for warmer temps
              ! and values at -23.0 for colder temps (Saunder and Peck, 1998)
              IF (temcg .gt. -8.0 ) THEN
                  t = -8.0
              ELSEIF (temcg .lt. -23.0) THEN
                  t = -23.0
              END IF

              qsign = -1.0
              IF (rar .gt. rarc) qsign = 1.0

              IF (qsign .gt. 0.5) THEN

                  !        IF ( Abs(rar - rarc) .lt. delr) THEN
                  !         ftrar = delri*Abs(rar - rarc)*
                  !     :           (6.74*(rarc + delr) + 1.36*t + 10.5 )
                  !        ELSE
                  !         ftrar = 6.74*rar + 1.36*t + 10.5
                  !        END IF
                  !        ftrar = Max(0.0, ftrar)

                  !        IF ( Abs(rar - rarc) .lt. delr .and. temcg .gt. tc) THEN
                  !         ftrar = Max( 0.0, delri*Abs(rar - rarc)*
                  !     :           (6.74*(rarc + delr) + 1.36*t + 10.5 ) )
                  !        ELSE
                  !         ftrar = Max( 0.0, 6.74*rar + 1.36*t + 10.5 )
                  !        END IF

                  IF ( temcg .gt. -34.0 ) THEN
                      ftrar = 6.74*( rar - rarc )
                  ELSE
                      ftrar = Max ( 0.0, 6.74*( rar -      &
                              &     4*( 1.0 - Abs( ( (-34.0 + 25.0 )/(-25.0 - tmin) )**2) ) ) )
                  END IF

                  if ( awdia*1.e6.lt.155. ) then
                      qconkq = 4.9e13
                      qconm  = 3.76
                      qconn  = 2.5
                  end if
                  if ( awdia*1.e6.ge.155. .and. awdia*1.e6.le.452. ) then
                      qconkq = 4.0e6
                      qconm  = 1.9
                      qconn  = 2.5
                  end if
                  if ( awdia*1.e6.gt.452. ) then
                      qconkq = 52.8
                      qconm  = 0.44
                      qconn  = 2.5
                  end if

              ELSEIF (qsign .lt. -0.5) THEN

                  !        IF ( Abs(rar - rarc) .lt. delr) THEN
                  !         ftrar = delri*Abs(rar - rarc)*
                  !     :           (3.02 - 10.59*(rarc-delr) + 2.95*(rarc-delr)**2 )
                  !        ELSEIF (rar .lt. 0.4) THEN
                  !         ftrar = ((rar-0.1)/0.3)*(3.02 - 10.59*(0.4) + 2.95*(0.4)**2)
                  !        ELSE
                  !         ftrar = 3.02 - 10.59*rar + 2.95*rar**2
                  !        END IF

                  ftrar = Min(1.0, 0.6*Abs(rarc - rar0) )*(6.5)*(-1.0 +      &
                          &    4.0/(rarc - rar0)**2 * (rar - (rarc + rar0)/2.0 )**2  )
                  ftrar = Min( 0.0, ftrar )

                  if ( awdia*1.e6.lt.253. ) then
                      qconkq = 5.24e8
                      qconm  = 2.54
                      qconn  = 2.8
                  end if
                  if ( awdia*1.e6.gt.253. ) then
                      qconkq = 24.0
                      qconm  = 0.50
                      qconn  = 2.8
                  end if
              END IF


              999  CONTINUE
              RETURN
          END subroutine saund5

          ! ####################################################################
          !                SUBROUTINE SAUND4
          ! ####################################################################
          subroutine saund4(temcg,qcw,exw,vt,awdia,rho0,qsign,ftrar,     &
                  &                 qconkq,qconm,qconn,idelq)

              ! ####################################################################
              !
              !  Purpose:
              !   Saunders charging scheme based on rime accretion rate (RAR)
              !   as in Brooks et al. (1997) and Saunders and Peck (1998)
              !
              ! 10/12/2003:  Add line 'rarc = Max(rarc, 0.0)' to prevent negative
              !              rarc for -43 < T < -40
              !
              !
              ! 7/24/99  new adjustment to RARcritical curve: make a cubic drop-off
              !          from -23.7 to -40 (and zero charging below -40)
              !
              ! ####################################################################

              implicit none
              real temcg ! temperature
              real qcw ! cloud water mixing ratio
              real exw ! cloud water collection efficiency
              real vt ! terminal speed difference between x and cw
              real awdia ! crystal diameter
              real rho0  ! air density
              real qsign ! sign of charge acquired by rimer (not used by calling prog)
              real ftrar ! charge factor based on temp and RAR
              real qconkq ! factor kq ( or 'B')
              real qconm  ! exponent on crystal diameter ('a')
              real qconn  ! exponent on speed ('b')
              integer idelq ! charge sign, not used by calling program

              real rar ! rime accretion rate
              real rarc ! critical RAR
              real t,tc
              parameter (tc = -23.7)
              real delr,delri
              parameter (delr=0.5, delri=1.0/delr)
              real rar0  ! lower limit of RAR where charging goes to zero
              parameter (rar0 = 0.1)
              real q0
              parameter (q0 = 6.48)

              ! ####################################################################
              ! Begin Executable code
              ! ####################################################################

              ftrar = 0.0

              rar = exw*qcw*1.0e3*rho0*vt
              t = temcg
              !      IF (t .gt. -23.0 ) THEN
              !      IF (t .gt. -30.0 ) THEN
              IF (t .gt. tc ) THEN
                  rarc = 1.0 + t*(7.9262e-2 + t*(4.4847e-2 +     &
                          &  t*(7.4754e-3 + t*(5.4686e-4 +     &
                                  &  t*(1.6737e-5 + t*1.7613e-7)))))
                  !      ELSE ! IF ( t .le. -23.0) THEN
                  !       rarc = 3.27
                  !      ELSE ! IF ( t .le. -30.0) THEN
                  !       rarc = 1.795
              ELSE ! IF ( t .le. tc) THEN
                  rarc = 3.39608*( 1.0 - Abs( ( (t - tc)/(tc + 40) )**3) )
              END IF


              rarc = Max(rarc, 0.0)


              ! check for RAR below threshold:
              IF ( rar .le. rar0 .or. temcg .le. -40.0) THEN
                  qsign = 0.0
                  ftrar = 0.0
                  qconkq = 0.0
                  qconn = 1.0
                  qconm = 1.0
                  GOTO 999
              END IF

              ! for now, use charging values at -8.0 for warmer temps
              ! and values at -23.0 for colder temps (Saunder and Peck, 1998)
              IF (temcg .gt. -8.0 ) THEN
                  t = -8.0
              ELSEIF (temcg .lt. -23.0) THEN
                  t = -23.0
              END IF

              qsign = -1.0
              IF (rar .gt. rarc) qsign = 1.0

              IF (qsign .gt. 0.5) THEN

                  !        IF ( Abs(rar - rarc) .lt. delr) THEN
                  !         ftrar = delri*Abs(rar - rarc)*
                  !     :           (6.74*(rarc + delr) + 1.36*t + 10.5 )
                  !        ELSE
                  !         ftrar = 6.74*rar + 1.36*t + 10.5
                  !        END IF
                  !        ftrar = Max(0.0, ftrar)

                  IF ( Abs(rar - rarc) .lt. delr .and. temcg .gt. tc) THEN
                      ftrar = Max( 0.0, delri*Abs(rar - rarc)*     &
                              &           (6.74*(rarc + delr) + 1.36*t + 10.5 ) )
                  ELSE
                      ftrar = Max( 0.0, 6.74*rar + 1.36*t + 10.5 )
                  END IF

                  if ( awdia*1.e6.lt.155. ) then
                      qconkq = 4.9e13
                      qconm  = 3.76
                      qconn  = 2.5
                  end if
                  if ( awdia*1.e6.ge.155. .and. awdia*1.e6.le.452. ) then
                      qconkq = 4.0e6
                      qconm  = 1.9
                      qconn  = 2.5
                  end if
                  if ( awdia*1.e6.gt.452. ) then
                      qconkq = 52.8
                      qconm  = 0.44
                      qconn  = 2.5
                  end if

              ELSEIF (qsign .lt. -0.5) THEN

                  !        IF ( Abs(rar - rarc) .lt. delr) THEN
                  !         ftrar = delri*Abs(rar - rarc)*
                  !     :           (3.02 - 10.59*(rarc-delr) + 2.95*(rarc-delr)**2 )
                  !        ELSEIF (rar .lt. 0.4) THEN
                  !         ftrar = ((rar-0.1)/0.3)*(3.02 - 10.59*(0.4) + 2.95*(0.4)**2)
                  !        ELSE
                  !         ftrar = 3.02 - 10.59*rar + 2.95*rar**2
                  !        END IF

                  ftrar = Min(1.0, 0.6*Abs(rarc - rar0) )*(7.0)*(-1.0 +      &
                          &    4.0/(rarc - rar0)**2 * (rar - (rarc + rar0)/2.0 )**2  )
                  ftrar = Min( 0.0, ftrar )

                  if ( awdia*1.e6.lt.253. ) then
                      qconkq = 5.24e8
                      qconm  = 2.54
                      qconn  = 2.8
                  end if
                  if ( awdia*1.e6.gt.253. ) then
                      qconkq = 24.0
                      qconm  = 0.50
                      qconn  = 2.8
                  end if
              END IF


              999  CONTINUE
              RETURN
          END subroutine saund4

          ! ####################################################################
          !                SUBROUTINE SAUND3
          ! ####################################################################
          subroutine saund3(temcg,qcw,exw,vt,awdia,rho0,qsign,ftrar,     &
                  &                 qconkq,qconm,qconn,idelq)

              ! ####################################################################
              !
              !  Purpose:
              !   Saunders charging scheme based on rime accretion rate (RAR)
              !   as in Brooks et al. (1997) and Saunders and Peck (1998)
              !
              !  8/2/99  Put in quadratic drop-off in transferred charge
              !          starting at -30 (fac=1.0) down to -43 (fac=0.0)
              !
              ! ####################################################################

              implicit none
              real temcg ! temperature
              real qcw ! cloud water mixing ratio
              real exw ! cloud water collection efficiency
              real vt ! terminal speed difference between x and cw
              real awdia ! crystal diameter
              real rho0  ! air density
              real qsign ! sign of charge acquired by rimer (not used by calling prog)
              real ftrar ! charge factor based on temp and RAR
              real qconkq ! factor kq ( or 'B')
              real qconm  ! exponent on crystal diameter ('a')
              real qconn  ! exponent on speed ('b')
              integer idelq ! charge sign, not used by calling program

              real rar ! rime accretion rate
              real rarc ! critical RAR
              real t
              real delr,delri
              parameter (delr=0.5, delri=1.0/delr)
              real rar0  ! lower limit of RAR where charging goes to zero
              parameter (rar0 = 0.1)
              real q0
              parameter (q0 = 6.48)
              real fac

              ! ####################################################################
              ! Begin Executable code
              ! ####################################################################

              ftrar = 0.0

              IF ( temcg .gt. -30 ) THEN
                  fac = 1.0
              ELSEIF ( temcg .gt. -40.0 ) THEN
                  fac = 1.0 - ((temcg+30.0)/(40.0 - 30.0))**2
              ELSE
                  fac = 0.0
              END IF

              rar = exw*qcw*1.0e3*rho0*vt
              t = temcg
              IF (t .gt. -23.0 ) THEN
                  !      IF (t .gt. -30.0 ) THEN
                  rarc = 1.0 + t*(7.9262e-2 + t*(4.4847e-2 +     &
                          &  t*(7.4754e-3 + t*(5.4686e-4 +     &
                                  &  t*(1.6737e-5 + t*1.7613e-7)))))
              ELSE ! IF ( t .le. -23.0) THEN
                  rarc = 3.27
                  !      ELSE ! IF ( t .le. -30.0) THEN
                  !       rarc = 1.795
              END IF


              ! check for RAR below threshold:
              IF ( rar .le. rar0 ) THEN
                  qsign = 0.0
                  ftrar = 0.0
                  qconkq = 0.0
                  qconn = 1.0
                  qconm = 1.0
                  GOTO 999
              END IF

              ! for now, use charging values at -8.0 for warmer temps
              ! and values at -23.0 for colder temps (Saunder and Peck, 1998)
              IF (temcg .gt. -8.0 ) THEN
                  t = -8.0
              ELSEIF (temcg .lt. -23.0) THEN
                  t = -23.0
              END IF

              qsign = -1.0
              IF (rar .gt. rarc) qsign = 1.0

              IF (qsign .gt. 0.5) THEN

                  !        IF ( Abs(rar - rarc) .lt. delr) THEN
                  !         ftrar = delri*Abs(rar - rarc)*
                  !     :           (6.74*(rarc + delr) + 1.36*t + 10.5 )
                  !        ELSE
                  !         ftrar = 6.74*rar + 1.36*t + 10.5
                  !        END IF
                  !        ftrar = Max(0.0, ftrar)

                  ftrar = fac*Max( 0.0, 6.74*rar + 1.36*t + 10.5 )

                  if ( awdia*1.e6.lt.155. ) then
                      qconkq = 4.9e13
                      qconm  = 3.76
                      qconn  = 2.5
                  end if
                  if ( awdia*1.e6.ge.155. .and. awdia*1.e6.le.452. ) then
                      qconkq = 4.0e6
                      qconm  = 1.9
                      qconn  = 2.5
                  end if
                  if ( awdia*1.e6.gt.452. ) then
                      qconkq = 52.8
                      qconm  = 0.44
                      qconn  = 2.5
                  end if

              ELSEIF (qsign .lt. -0.5) THEN

                  !        IF ( Abs(rar - rarc) .lt. delr) THEN
                  !         ftrar = delri*Abs(rar - rarc)*
                  !     :           (3.02 - 10.59*(rarc-delr) + 2.95*(rarc-delr)**2 )
                  !        ELSEIF (rar .lt. 0.4) THEN
                  !         ftrar = ((rar-0.1)/0.3)*(3.02 - 10.59*(0.4) + 2.95*(0.4)**2)
                  !        ELSE
                  !         ftrar = 3.02 - 10.59*rar + 2.95*rar**2
                  !        END IF

                  ftrar = fac*q0*(-1.0 +      &
                          &    4.0/(rarc - rar0)**2 * (rar - (rarc + rar0)/2.0 )**2  )

                  if ( awdia*1.e6.lt.253. ) then
                      qconkq = 5.24e8
                      qconm  = 2.54
                      qconn  = 2.8
                  end if
                  if ( awdia*1.e6.gt.253. ) then
                      qconkq = 24.0
                      qconm  = 0.50
                      qconn  = 2.8
                  end if
              END IF


              999  CONTINUE
              RETURN
          END subroutine saund3

          !
          ! 6/21/99
          !         Also put in a linear drop-off for -7.4 < T < 0
          !       try setting elwc = 0.5*lwc
          !
          !  11.11.03  ERM fixed some bugs that allowed charging at EW below
          !            threshold for T < -20.  Also corrected charging at low EW
          !            at -16 > T > -20.
          !
          ! ####################################################################
          !                SUBROUTINE SAUND
          ! ####################################################################
          subroutine saund(temcg,qcw,exw,awdia,rho0,qsign,ftelwc,     &
                  &                 qconkq,qconm,qconn,idelq,ianom)

              implicit none

              integer ianom  ! =0 to reduce 'anomalous' regions
              ! =1 to keep 'anonmalous' regions unchanged
              ! =2 to reduce neg. anom. zone and
              !     remove the pos. anom. zone (use normal neg. instead)
              ! =3 to run same as Helsdon et al. 2001
              ! =4 to run Helsdon et al. 2001 without pos. anom. zone
              ! =5 to run Helsdon et al. 2001 without pos. or neg. anom. zones (use normal instead)
              ! =(anything else) to remove the anomalous zones altogether
              real anom
              integer ieq, iseq

              real temcg ! temperature
              real qcw ! cloud water mixing ratio
              real exw ! cloud water collection efficiency
              !      real vt ! terminal speed difference between x and cw
              real awdia ! crystal diameter
              real rho0  ! air density
              real qsign ! sign of charge acquired by rimer (not used by calling prog)
              real ftrar ! charge factor based on temp and RAR
              real qconkq ! factor kq ( or 'B')
              real qconm  ! exponent on crystal diameter ('a')
              real qconn  ! exponent on speed ('b')
              integer idelq ! charge sign, not used by calling program
              real ftelwc, fac

              real elwc, trevsau, qlwc, pi, cew
              !
              ! ####################################################
              !

              IF ( temcg .gt. -30 ) THEN
                  fac = 1.0
              ELSEIF ( temcg .gt. -40.0 ) THEN
                  fac = 1.0 - ((temcg+30.0)/(40.0 - 30.0))**2
              ELSE
                  fac = 0.0
              END IF

              ieq = 0
              pi    = 4.*atan(1.0)
              !
              qlwc = qcw*(1.e3)*rho0
              elwc = qlwc*Min(1.0,exw)
              !      elwc = qlwc*0.5
              !
              !
              !   Saunders et al. 1991 non-inductive parameterizations
              !
              !
              idelq = 0
              ftelwc = 0.0
              qsign = 0.0
              qconkq = 0.0
              qconm = 1.0
              qconn = 1.0

              IF ( ( temcg .lt. -20.0 .and. elwc .le. 0.061 ) .or.     &
                      &     (  elwc .le. 0.026 ) ) THEN
                  RETURN
              ENDIF
              !
              !   case I  (positive anomalous zone)
              !
              IF ( ( temcg.lt.(-20.0) .or.      &
                      &    ( temcg.lt.(-18.0) .and. ianom .ne. 3 ) )     &
                      &      .and. elwc.lt.(0.1596) .and. ianom .ne. 2      &
                      &       .and. ianom .ne. 4 .and. ianom .ne. 5 ) then
                  !
                  ieq = 1
                  idelq = 1
                  !
                  !  sign
                  !
                  qsign = 1.0
                  !
                  !  f(t,lwc)
                  !
                  ftelwc = 0.0

                  if (elwc.lt.(0.12).and.elwc.gt.(0.061)) then
                      ftelwc = ( (2041.76)*elwc-(128.70) )
                  end if
                  if (elwc.lt.(0.1596).and.elwc.ge.(0.12)) then
                      ftelwc = ( (-2900.22)*elwc+(462.91) )
                  end if

                  IF ( temcg .gt. -20.0 ) THEN
                      ftelwc = ftelwc*( -temcg - 18.0 )*0.5
                  ENDIF
                  !
                  ! f(awdia)
                  !

                  IF ( ianom .eq. 0 .or. ianom .eq. 3 ) THEN
                      anom = 0.1
                  ELSEIF ( ianom .eq. 1 ) THEN
                      anom = 1.0
                  ELSE
                      anom = 0.0
                  ENDIF

                  if ( awdia*1.e6.lt.155. ) then
                      qconkq = anom*4.9e13
                      qconm  = 3.76
                      qconn  = 2.5
                  end if
                  if ( awdia*1.e6.ge.155. .and.      &
                          &     awdia*1.e6.le.452. ) then
                      qconkq = anom*4.9e6
                      qconm  = 1.9
                      qconn  = 2.5
                  end if
                  if ( awdia*1.e6.gt.452. ) then
                      qconkq = anom*52.8
                      qconm  = 0.44
                      qconn  = 2.5
                  end if
                  !
                  !  end case I
                  !
              end if
              !
              !
              !   case II  (zero charging at low EW)
              !
              if ( ( (temcg.ge.(-20.0) .or. (ianom.eq.2) ) .and.      &
                      &      temcg.lt.(-16.0)) .and.      &
                      &     elwc.lt.(0.06) ) then
                  !
                  ieq = 2
                  idelq = 1
                  !
                  !  sign
                  !
                  qsign = 0.0
                  !
                  !  f(t,lwc)
                  !
                  ftelwc = 0.0
                  !
                  ! f(awdia)
                  !
                  qconkq = 0.0
                  qconm  = 1.0
                  qconn  = 1.0
                  !
                  !  end case II
                  !
              end if
              !
              !   case III  (negative anomalous zone)
              !
              if ( ( temcg.ge.(-16.0) .or.      &
                      &      ( temcg.ge.(-18.0) .and. ianom .ne. 3 .and. ianom .ne. 4) )     &
                      &        .and. elwc.lt.(0.22) .and. ianom .ne. 5 ) then
                  !
                  ieq = 3
                  idelq = 1
                  !
                  !  sign
                  !
                  qsign = -1.
                  !
                  !  f(t,lwc)
                  !
                  !  If elwc.lt.0.03:
                  ftelwc = 0.0                                  ! ERM bug fix 11.11.2003

                  if (elwc.lt.(0.14).and.elwc.gt.(0.026)) then
                      ftelwc = ( (-314.40)*elwc+(7.92) ) !iseq=2
                  end if
                  if (elwc.lt.(0.22).and.elwc.ge.(0.14)) then
                      ftelwc = ( (419.4)*elwc-(92.64) ) !iseq=3
                  end if

                  ftelwc = Min ( 0.0, ftelwc )

                  IF ( temcg .lt. -16.0 ) THEN
                      ftelwc = ftelwc*( temcg + 18.0 )*0.5
                  ENDIF

                  IF ( temcg .gt. -7.38 ) ftelwc = ftelwc*Abs(temcg/7.38)
                  !
                  ! f(awdia)
                  !

                  IF ( ianom .eq. 0 .or. (ianom .ge. 2 .and. ianom .le. 4) ) THEN
                      anom = 0.2
                  ELSEIF ( ianom .eq. 1 ) THEN
                      anom = 1.0
                  ELSE
                      anom = 0.0
                  ENDIF

                  if ( awdia*1.e6.lt.253. ) then
                      qconkq = anom*5.24e8
                      qconm  = 2.54
                      qconn  = 2.8
                  end if
                  if ( awdia*1.e6.gt.253. ) then
                      qconkq = anom*24.0
                      qconm  = 0.50
                      qconn  = 2.8
                  end if
                  !
                  !  end case III
                  !
              end if
              !
              !  ASSIGN NORMAL CHARGING ZONES
              !
              trevsau = -15.06*elwc-7.38

              iseq = 7 ! set default for [something]
              IF ( temcg .gt. -10.69 .and. elwc .gt. 0.22 ) THEN
                  cew = 0.22
                  iseq = 1
              ELSEIF (temcg .lt. -24. .and. elwc .gt. 1.1 ) THEN
                  cew = 1.1
                  iseq = 8
              ELSEIF (temcg .lt. -24.0 .and. elwc .le. 1.1 ) THEN
                  cew = 1.1
                  iseq = 7
              ELSEIF ( temcg .le. -10.69 .and. temcg .ge. -24.0 ) THEN
                  cew = -0.49 - (6.64e-2)*temcg
                  IF ( elwc .gt. cew ) THEN
                      iseq = 4
                  ELSE
                      iseq = 7
                  ENDIF
              ENDIF

              IF ( temcg .le. -7.38 .and. temcg .ge. -24.0 .and. ianom .eq. 5 ) THEN
                  cew = -0.49 - (6.64e-2)*temcg
                  IF ( elwc .gt. cew ) THEN
                      iseq = 4
                  ELSE
                      iseq = 7
                  ENDIF
              ENDIF


              !
              !   case IV (normal positive zone, above the critical line)
              !

              IF ( iseq .eq. 1 .or. iseq .eq. 8 .or. iseq .eq. 4 ) THEN
                  !
                  ieq = 4
                  idelq = 1
                  !
                  !  sign
                  !
                  qsign = 1.0
                  !
                  !  f(t,lwc)
                  !
                  IF ( iseq .eq. 4 ) THEN
                      ftelwc = (20.22)*elwc+(1.36)*temcg+(10.05)
                  ELSEIF ( iseq .eq. 1 .and. ianom /= 5 ) THEN
                      ftelwc = 20.22*elwc - 4.4484 ! value of T set to -10.69
                      !cccc      ENDIF
                  ELSEIF ( iseq .eq. 1 .and. ianom == 5 ) THEN
                      ftelwc = 20.22*(elwc - cew)
                  ELSEIF ( iseq .eq. 8 ) THEN
                      IF ( ianom .eq. 3 .or. ianom .eq. 4 ) THEN
                          ftelwc = Max(0.0, (20.22)*elwc - 22.24 )
                      ELSE  ! create region of zero charging
                          ftelwc = Max(0.0, (20.22)*elwc+(1.36)*temcg+(10.05))
                      ENDIF
                  ENDIF
                  IF ( temcg .gt. -7.38 )     &
                          &   ftelwc = Abs(temcg/7.38)*     &
                          &            ((20.22)*elwc + (1.36)*(-7.38) + 10.05)

                  !      IF ( ftelwc .lt. 0 ) ftelwc = 0.0
                  ftelwc = Max ( 0.0 , ftelwc )
                  !
                  ! f(awdia)
                  !
                  if ( awdia*1.e6.lt.155. ) then
                      qconkq = 4.9e13
                      qconm  = 3.76
                      qconn  = 2.5
                  end if
                  if ( awdia*1.e6.ge.155. .and.      &
                          &     awdia*1.e6.le.452. ) then
                      qconkq = 4.9e6
                      qconm  = 1.9
                      qconn  = 2.5
                  end if
                  if ( awdia*1.e6.gt.452. ) then
                      qconkq = 52.8
                      qconm  = 0.44
                      qconn  = 2.5
                  end if
                  !
                  !  end case IV
                  !
              END IF
              !
              !
              !   case V (normal negative zone, below the critical line)
              !
              !      if ( temcg.le.(trevsau) .and.
              !     >     elwc.lt.(1.1) .and.
              !     >     idelq .eq. 0 ) then
              IF ( idelq .eq. 0 .and. iseq .eq. 7 ) THEN
                  !
                  ieq = 7
                  idelq = 1
                  !
                  !  sign
                  !
                  qsign = -1.
                  !
                  !  f(t,lwc)
                  !
                  ftelwc = Min( 0.0, (3.02)-(31.76)*elwc+(26.53)*elwc**2 )
                  !
                  ! f(awdia)
                  !
                  if ( awdia*1.e6.lt.253. ) then
                      qconkq = 5.24e8
                      qconm  = 2.54
                      qconn  = 2.8
                  end if
                  if ( awdia*1.e6.gt.253. ) then
                      qconkq = 24.0
                      qconm  = 0.50
                      qconn  = 2.8
                  end if
                  !
                  !  end case V
                  !
              end if

              ftelwc = fac*ftelwc
              return
          end subroutine saund
          ! ####################################################################
          !
          !    ! #####################################################################
          !    ! #####################################################################
          !    !
          !    ! Subroutine for reattachment of floating ion charge
          !    !
          SUBROUTINE SCIONATT    &
                  &  (nx,ny,nz,na,jyslab &
                  &  ,nor,norz &
                  &  ,dtp,dz3d &
                  &  ,an &
                  &  )

              implicit none

              integer :: nx,ny,nz,na
              integer :: ng
              integer :: nor,norz, jyslab ! ,nht,ngt,igsr
              real    :: dtp  ! time step


              !
              ! external temporary arrays
              !

              real an(-nor+1:nx+nor,-nor+1:ny+nor,-norz+1:nz+norz,na)

              real dz3d(-nor+1:nx+nor,-nor+1:ny+nor,-norz+1:nz+norz)

              ! local

              integer ix,jy,kz,i,n

              double precision frachl
              double precision frach
              double precision fracr
              double precision fracw
              double precision fraci
              double precision fracs
              double precision fractot,frac
              double precision dv, deninv, tmp, tmpion

              ! density defined in nssl
              !        real    :: rhos = 100.
              !        real    :: rhog = 500. ! not correct, but for now
              !        real    :: rhoh = 900.
              !        real    :: rhor = 1000.
              !        real    :: rhoi = 900.
              !        real    :: rhow = 1000.
              ! density defined in morr
              real    :: rhos = 100.
              real    :: rhog = 400.
              real    :: rhor = 997.
              real    :: rhoi = 500. ! the most different
              real    :: rhow = 997.

              logical, parameter :: dampicecharge = .false.    ! whether to apply charge damping to cloud ice and snow (WRF ONLY)
              real, parameter    :: dampchargethresh = 0.09e-9 ! charge density magnitude at which damping stops
              real, parameter    :: dampchargetimeconst = 3600. ! time constant for charge damping

              IF ( lnh < 1 ) THEN ! not double-moment
                  RETURN
              ENDIF
              !        IF (IHAIL.EQ.1)  RHOG = 900.

              DO kz = 1,nz
                  DO jy = 1,1
                      DO ix = 1,nx
                          IF ( an(ix,jy,kz,lscpi) - an(ix,jy,kz,lscni) /= 0.0 ) THEN
                              fracs = 0.0
                              fraci = 0.0
                              frach = 0.0
                              fracw = 0.0
                              fracr = 0.0
                              frachl = 0.0
                              if (an(ix,jy,kz,lc) > qsmall) fracw=((an(ix,jy,kz,lc)/rhow)**(2./3.) ) * Max(an(ix,jy,kz,lnc), 0.0)**(1./3.)
                              if (an(ix,jy,kz,lr) > qsmall) fracr=((an(ix,jy,kz,lr)/rhor)**(2./3.) ) * Max(an(ix,jy,kz,lnr), 0.0)**(1./3.)
                              if (an(ix,jy,kz,li) > qsmall) fraci=((an(ix,jy,kz,li)/rhoi)**(2./3.) ) * Max(an(ix,jy,kz,lni), 0.0)**(1./3.)
                              if (an(ix,jy,kz,ls) > qsmall) fracs=((an(ix,jy,kz,ls)/rhos)**(2./3.) ) * Max(an(ix,jy,kz,lns), 0.0)**(1./3.)
                              if (an(ix,jy,kz,lh) > qsmall) frach=((an(ix,jy,kz,lh)/rhog)**(2./3.) ) * Max(an(ix,jy,kz,lnh), 0.0)**(1./3.)

                              fractot = fracw + fracr + fraci + fracs + frach + frachl

                              IF ( fractot > 0.0 ) THEN
                                  tmpion = (an(ix,jy,kz,lscpi) - an(ix,jy,kz,lscni))*ec ! C/m3

                                  IF (fracw > 0) then
                                      an(ix,jy,kz,lscw) = an(ix,jy,kz,lscw) + tmpion*(fracw/fractot)
                                  ENDIF

                                  IF (fracr > 0) then
                                      an(ix,jy,kz,lscr) = an(ix,jy,kz,lscr) + tmpion*(fracr/fractot)
                                  ENDIF

                                  IF (fraci > 0) then
                                      an(ix,jy,kz,lsci) = an(ix,jy,kz,lsci) + tmpion*(fraci/fractot)
                                  ENDIF

                                  IF (fracs > 0) then ! avoid division by zero IF fractot = 0 because qx=0
                                      an(ix,jy,kz,lscs) = an(ix,jy,kz,lscs) + tmpion*(fracs/fractot)
                                  ENDIF

                                  IF (frach > 0) then
                                      an(ix,jy,kz,lsch) = an(ix,jy,kz,lsch) + tmpion*(frach/fractot)
                                  ENDIF

                                  !                            IF (lnhl > 1 .and. frachl > 0 ) then
                                  !                                an(ix,jy,kz,lschl) = an(ix,jy,kz,lschl) + tmpion*(frachl/fractot)
                                  !                            ENDIF


                                  an(ix,jy,kz,lscpi) = 0.0
                                  an(ix,jy,kz,lscni) = 0.0

                              ENDIF ! fractot

                          ENDIF ! non-zero ion charge

                      ENDDO
                  ENDDO
              ENDDO

              IF ( dampicecharge ) THEN
                  frac = Exp(-dtp/dampchargetimeconst)
                  DO kz = 1,nz
                      DO jy = 1,1
                          DO ix = 1,nx
                              IF ( Abs( an(ix,jy,kz,lsci) ) > dampchargethresh ) THEN
                                  an(ix,jy,kz,lsci) = frac*an(ix,jy,kz,lsci)
                              ENDIF

                              IF ( Abs( an(ix,jy,kz,lscs) ) > dampchargethresh ) THEN
                                  an(ix,jy,kz,lscs) = frac*an(ix,jy,kz,lscs)
                              ENDIF

                          ENDDO
                      ENDDO
                  ENDDO
              ENDIF

              RETURN
          END SUBROUTINE SCIONATT

          ! ##############################################################################

          !    -----------------------------------------------------------------------!
          !    ----------------------  End Explicit Lightning  -----------------------!
      END MODULE MODULE_MP_AHM

      !
