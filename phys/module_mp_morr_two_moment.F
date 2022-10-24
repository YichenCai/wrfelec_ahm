!WRF:MODEL_LAYER:PHYSICS
!

! THIS MODULE CONTAINS THE TWO-MOMENT MICROPHYSICS CODE DESCRIBED BY
!     MORRISON ET AL. (2009, MWR)

! CHANGES FOR V3.2, RELATIVE TO MOST RECENT (BUG-FIX) CODE FOR V3.1

! 1) ADDED ACCELERATED MELTING OF GRAUPEL/SNOW DUE TO COLLISION WITH RAIN, FOLLOWING LIN ET AL. (1983)
! 2) INCREASED MINIMUM LAMBDA FOR RAIN, AND ADDED RAIN DROP BREAKUP FOLLOWING MODIFIED VERSION
!     OF VERLINDE AND COTTON (1993)
! 3) CHANGE MINIMUM ALLOWED MIXING RATIOS IN DRY CONDITIONS (RH < 90%), THIS IMPROVES RADAR REFLECTIIVITY
!     IN LOW REFLECTIVITY REGIONS
! 4) BUG FIX TO MAXIMUM ALLOWED PARTICLE FALLSPEEDS AS A FUNCTION OF AIR DENSITY
! 5) BUG FIX TO CALCULATION OF LIQUID WATER SATURATION VAPOR PRESSURE (CHANGE IS VERY MINOR)
! 6) INCLUDE WRF CONSTANTS PER SUGGESTION OF JIMY

! bug fix, 5/12/10
! 7) bug fix for saturation vapor pressure in low pressure, to avoid division by zero
! 8) include 'EP2' WRF constant for saturation mixing ratio calculation, instead of hardwire constant

! CHANGES FOR V3.3
! 1) MODIFICATION FOR COUPLING WITH WRF-CHEM (PREDICTED DROPLET NUMBER CONCENTRATION) AS AN OPTION
! 2) MODIFY FALLSPEED BELOW THE LOWEST LEVEL OF PRECIPITATION, WHICH PREVENTS
!      POTENTIAL FOR SPURIOUS ACCUMULATION OF PRECIPITATION DURING SUB-STEPPING FOR SEDIMENTATION
! 3) BUG FIX TO LATENT HEAT RELEASE DUE TO COLLISIONS OF CLOUD ICE WITH RAIN
! 4) CLEAN UP OF COMMENTS IN THE CODE

! additional minor bug fixes and small changes, 5/30/2011
! minor revisions by A. Ackerman April 2011:
! 1) replaced kinematic with dynamic viscosity
! 2) replaced scaling by air density for cloud droplet sedimentation
!    with viscosity-dependent Stokes expression
! 3) use Ikawa and Saito (1991) air-density scaling for cloud ice
! 4) corrected typo in 2nd digit of ventilation constant F2R

! additional fixes:
! 5) TEMPERATURE FOR ACCELERATED MELTING DUE TO COLLIIONS OF SNOW AND GRAUPEL
!    WITH RAIN SHOULD USE CELSIUS, NOT KELVIN (BUG REPORTED BY K. VAN WEVERBERG)
! 6) NPRACS IS NOT SUBTRACTED FROM SNOW NUMBER CONCENTRATION, SINCE
!    DECREASE IN SNOW NUMBER IS ALREADY ACCOUNTED FOR BY NSMLTS
! 7) fix for switch for running w/o graupel/hail (cloud ice and snow only)

! hm bug fix 3/16/12

! 1) very minor change to limits on autoconversion source of rain number when cloud water is depleted

! WRFV3.5
! hm/A. Ackerman bug fix 11/08/12

! 1) for accelerated melting from collisions, should use rain mass collected by snow, not snow mass
!    collected by rain
! 2) minor changes to some comments
! 3) reduction of maximum-allowed ice concentration from 10 cm-3 to 0.3
!    cm-3. This was done to address the problem of excessive and persistent
!    anvil cirrus produced by the scheme.

! CHANGES FOR WRFV3.5.1
! 1) added output for snow+cloud ice and graupel time step and accumulated
!    surface precipitation
! 2) bug fix to option w/o graupel/hail (IGRAUP = 1), include PRACI, PGSACW,
!    and PGRACS as sources for snow instead of graupel/hail, bug reported by
!    Hailong Wang (PNNL)
! 3) very minor fix to immersion freezing rate formulation (negligible impact)
! 4) clarifications to code comments
! 5) minor change to shedding of rain, remove limit so that the number of
!    collected drops can smaller than number of shed drops
! 6) change of specific heat of liquid water from 4218 to 4187 J/kg/K

! CHANGES FOR WRFV3.6.1
! 1) minor bug fix to melting of snow and graupel, an extra factor of air density (RHO) was removed
!    from the calculation of PSMLT and PGMLT
! 2) redundant initialization of PSMLT (non answer-changing)

! CHANGES FOR WRFV3.8.1
! 1) changes and cleanup of code comments
! 2) correction to universal gas constant (very small change)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! THIS SCHEME IS A BULK DOUBLE-MOMENT SCHEME THAT PREDICTS MIXING
! RATIOS AND NUMBER CONCENTRATIONS OF FIVE HYDROMETEOR SPECIES:
! CLOUD DROPLETS, CLOUD (SMALL) ICE, RAIN, SNOW, AND GRAUPEL/HAIL.

MODULE MODULE_MP_MORR_TWO_MOMENT
    USE     module_wrf_error
    USE module_utility, ONLY: WRFU_Clock, WRFU_Alarm  ! GT
    USE module_domain, ONLY : HISTORY_ALARM, Is_alarm_tstep  ! GT
    USE module_mp_radar

    ! USE WRF PHYSICS CONSTANTS
    use module_model_constants, ONLY: CP, G, R => r_d, RV => r_v, EP_2

    !  USE module_state_description

    IMPLICIT NONE

    REAL, PARAMETER :: PI = 3.1415926535897932384626434
    REAL, PARAMETER :: xxx = 0.9189385332046727417803297

    PUBLIC  ::  MP_MORR_TWO_MOMENT
    PUBLIC  ::  POLYSVP

    PRIVATE :: GAMMA, DERF1
    PRIVATE :: PI, xxx
    PRIVATE :: MORR_TWO_MOMENT_MICRO
    PRIVATE :: DELBK, DELABK, GAMMA_DP

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! SWITCHES FOR MICROPHYSICS SCHEME
    ! IACT = 1, USE POWER-LAW CCN SPECTRA, NCCN = CS^K
    ! IACT = 2, USE LOGNORMAL AEROSOL SIZE DIST TO DERIVE CCN SPECTRA
    ! IACT = 3, ACTIVATION CALCULATED IN MODULE_MIXACTIVATE

    INTEGER, PRIVATE ::  IACT

    ! INUM = 0, PREDICT DROPLET CONCENTRATION
    ! INUM = 1, ASSUME CONSTANT DROPLET CONCENTRATION
    ! !!!NOTE: PREDICTED DROPLET CONCENTRATION NOT AVAILABLE IN THIS VERSION
    ! CONTACT HUGH MORRISON (morrison@ucar.edu) FOR FURTHER INFORMATION

    INTEGER, PRIVATE ::  INUM

    ! FOR INUM = 1, SET CONSTANT DROPLET CONCENTRATION (CM-3)
    REAL, PRIVATE ::      NDCNST

    ! SWITCH FOR LIQUID-ONLY RUN
    ! ILIQ = 0, INCLUDE ICE
    ! ILIQ = 1, LIQUID ONLY, NO ICE

    INTEGER, PRIVATE ::  ILIQ

    ! SWITCH FOR ICE NUCLEATION
    ! INUC = 0, USE FORMULA FROM RASMUSSEN ET AL. 2002 (MID-LATITUDE)
    !      = 1, USE MPACE OBSERVATIONS

    INTEGER, PRIVATE ::  INUC

    ! IBASE = 1, NEGLECT DROPLET ACTIVATION AT LATERAL CLOUD EDGES DUE TO
    !             UNRESOLVED ENTRAINMENT AND MIXING, ACTIVATE
    !             AT CLOUD BASE OR IN REGION WITH LITTLE CLOUD WATER USING
    !             NON-EQULIBRIUM SUPERSATURATION,
    !             IN CLOUD INTERIOR ACTIVATE USING EQUILIBRIUM SUPERSATURATION
    ! IBASE = 2, ASSUME DROPLET ACTIVATION AT LATERAL CLOUD EDGES DUE TO
    !             UNRESOLVED ENTRAINMENT AND MIXING DOMINATES,
    !             ACTIVATE DROPLETS EVERYWHERE IN THE CLOUD USING NON-EQUILIBRIUM
    !             SUPERSATURATION, BASED ON THE
    !             LOCAL SUB-GRID AND/OR GRID-SCALE VERTICAL VELOCITY
    !             AT THE GRID POINT

    ! NOTE: ONLY USED FOR PREDICTED DROPLET CONCENTRATION (INUM = 0) IN NON-WRF-CHEM VERSION OF CODE

    INTEGER, PRIVATE ::  IBASE

    ! INCLUDE SUB-GRID VERTICAL VELOCITY IN DROPLET ACTIVATION
    ! ISUB = 0, INCLUDE SUB-GRID W (RECOMMENDED FOR LOWER RESOLUTION)
    ! ISUB = 1, EXCLUDE SUB-GRID W, ONLY USE GRID-SCALE W

    ! NOTE: ONLY USED FOR PREDICTED DROPLET CONCENTRATION (INUM = 0) IN NON-WRF-CHEM VERSION OF CODE

    INTEGER, PRIVATE ::  ISUB

    ! SWITCH FOR GRAUPEL/NO GRAUPEL
    ! IGRAUP = 0, INCLUDE GRAUPEL
    ! IGRAUP = 1, NO GRAUPEL

    INTEGER, PRIVATE ::  IGRAUP

    ! HM ADDED NEW OPTION FOR HAIL
    ! SWITCH FOR HAIL/GRAUPEL
    ! IHAIL = 0, DENSE PRECIPITATING ICE IS GRAUPEL
    ! IHAIL = 1, DENSE PRECIPITATING GICE IS HAIL

    INTEGER, PRIVATE ::  IHAIL

    ! CLOUD MICROPHYSICS CONSTANTS

    REAL, PRIVATE ::      AI,AC,AS,AR,AG ! 'A' PARAMETER IN FALLSPEED-DIAM RELATIONSHIP
    REAL, PRIVATE ::      BI,BC,BS,BR,BG ! 'B' PARAMETER IN FALLSPEED-DIAM RELATIONSHIP
    !     REAL, PRIVATE ::      R           ! GAS CONSTANT FOR AIR
    !     REAL, PRIVATE ::      RV          ! GAS CONSTANT FOR WATER VAPOR
    !     REAL, PRIVATE ::      CP          ! SPECIFIC HEAT AT CONSTANT PRESSURE FOR DRY AIR
    REAL, PRIVATE ::      RHOSU       ! STANDARD AIR DENSITY AT 850 MB
    REAL, PRIVATE ::      RHOW        ! DENSITY OF LIQUID WATER
    REAL, PRIVATE ::      RHOI        ! BULK DENSITY OF CLOUD ICE
    REAL, PRIVATE ::      RHOSN       ! BULK DENSITY OF SNOW
    REAL, PRIVATE ::      RHOG        ! BULK DENSITY OF GRAUPEL
    REAL, PRIVATE ::      AIMM        ! PARAMETER IN BIGG IMMERSION FREEZING
    REAL, PRIVATE ::      BIMM        ! PARAMETER IN BIGG IMMERSION FREEZING
    REAL, PRIVATE ::      ECR         ! COLLECTION EFFICIENCY BETWEEN DROPLETS/RAIN AND SNOW/RAIN
    REAL, PRIVATE ::      DCS         ! THRESHOLD SIZE FOR CLOUD ICE AUTOCONVERSION
    REAL, PRIVATE ::      MI0         ! INITIAL SIZE OF NUCLEATED CRYSTAL
    REAL, PRIVATE ::      MG0         ! MASS OF EMBRYO GRAUPEL
    REAL, PRIVATE ::      F1S         ! VENTILATION PARAMETER FOR SNOW
    REAL, PRIVATE ::      F2S         ! VENTILATION PARAMETER FOR SNOW
    REAL, PRIVATE ::      F1R         ! VENTILATION PARAMETER FOR RAIN
    REAL, PRIVATE ::      F2R         ! VENTILATION PARAMETER FOR RAIN
    !     REAL, PRIVATE ::      G           ! GRAVITATIONAL ACCELERATION
    REAL, PRIVATE ::      QSMALL      ! SMALLEST ALLOWED HYDROMETEOR MIXING RATIO
    REAL, PRIVATE ::      CI,DI,CS,DS,CG,DG ! SIZE DISTRIBUTION PARAMETERS FOR CLOUD ICE, SNOW, GRAUPEL
    REAL, PRIVATE ::      EII         ! COLLECTION EFFICIENCY, ICE-ICE COLLISIONS
    REAL, PRIVATE ::      ECI         ! COLLECTION EFFICIENCY, ICE-DROPLET COLLISIONS
    REAL, PRIVATE ::      RIN     ! RADIUS OF CONTACT NUCLEI (M)
    ! hm, add for V3.2
    REAL, PRIVATE ::      CPW     ! SPECIFIC HEAT OF LIQUID WATER

    ! CCN SPECTRA FOR IACT = 1

    REAL, PRIVATE ::      C1     ! 'C' IN NCCN = CS^K (CM-3)
    REAL, PRIVATE ::      K1     ! 'K' IN NCCN = CS^K

    ! AEROSOL PARAMETERS FOR IACT = 2

    REAL, PRIVATE ::      MW      ! MOLECULAR WEIGHT WATER (KG/MOL)
    REAL, PRIVATE ::      OSM     ! OSMOTIC COEFFICIENT
    REAL, PRIVATE ::      VI      ! NUMBER OF ION DISSOCIATED IN SOLUTION
    REAL, PRIVATE ::      EPSM    ! AEROSOL SOLUBLE FRACTION
    REAL, PRIVATE ::      RHOA    ! AEROSOL BULK DENSITY (KG/M3)
    REAL, PRIVATE ::      MAP     ! MOLECULAR WEIGHT AEROSOL (KG/MOL)
    REAL, PRIVATE ::      MA      ! MOLECULAR WEIGHT OF 'AIR' (KG/MOL)
    REAL, PRIVATE ::      RR      ! UNIVERSAL GAS CONSTANT
    REAL, PRIVATE ::      BACT    ! ACTIVATION PARAMETER
    REAL, PRIVATE ::      RM1     ! GEOMETRIC MEAN RADIUS, MODE 1 (M)
    REAL, PRIVATE ::      RM2     ! GEOMETRIC MEAN RADIUS, MODE 2 (M)
    REAL, PRIVATE ::      NANEW1  ! TOTAL AEROSOL CONCENTRATION, MODE 1 (M^-3)
    REAL, PRIVATE ::      NANEW2  ! TOTAL AEROSOL CONCENTRATION, MODE 2 (M^-3)
    REAL, PRIVATE ::      SIG1    ! STANDARD DEVIATION OF AEROSOL S.D., MODE 1
    REAL, PRIVATE ::      SIG2    ! STANDARD DEVIATION OF AEROSOL S.D., MODE 2
    REAL, PRIVATE ::      F11     ! CORRECTION FACTOR FOR ACTIVATION, MODE 1
    REAL, PRIVATE ::      F12     ! CORRECTION FACTOR FOR ACTIVATION, MODE 1
    REAL, PRIVATE ::      F21     ! CORRECTION FACTOR FOR ACTIVATION, MODE 2
    REAL, PRIVATE ::      F22     ! CORRECTION FACTOR FOR ACTIVATION, MODE 2
    REAL, PRIVATE ::      MMULT   ! MASS OF SPLINTERED ICE PARTICLE
    REAL, PRIVATE ::      LAMMAXI,LAMMINI,LAMMAXR,LAMMINR,LAMMAXS,LAMMINS,LAMMAXG,LAMMING

    ! CONSTANTS TO IMPROVE EFFICIENCY

    REAL, PRIVATE :: CONS1,CONS2,CONS3,CONS4,CONS5,CONS6,CONS7,CONS8,CONS9,CONS10
    REAL, PRIVATE :: CONS11,CONS12,CONS13,CONS14,CONS15,CONS16,CONS17,CONS18,CONS19,CONS20
    REAL, PRIVATE :: CONS21,CONS22,CONS23,CONS24,CONS25,CONS26,CONS27,CONS28,CONS29,CONS30
    REAL, PRIVATE :: CONS31,CONS32,CONS33,CONS34,CONS35,CONS36,CONS37,CONS38,CONS39,CONS40
    REAL, PRIVATE :: CONS41

    !-----------------------------------------------------------------------
    !Explicit lightning shared vars-----------------------------------------

    integer :: ipelec = 0
    integer :: isaund = 0
    logical :: idonic = .false.
    logical, private :: cleardiag = .false.

    !  electrical permitivity of air C**2 / (N m**2)
    real, private :: eperao
    parameter (eperao  = 8.8592e-12 )

    ! fundamental unit of charge
    real ec,ec_i
    parameter (ec = 1.602e-19) ! Coulomb / proton or electron
    parameter (ec_i = 1.0/ec) ! #proton or electron / Coulomb

    real    :: scwppmx = 20.0e-12
    real    :: scippmx = 20.0e-12

    integer :: nic_noliq = 0  ! switch for NIC without riming (MST, 2006)
    real    :: qc_noliq = 0.1 ! charge per collision (in fC) for non-riming NIC.
    logical :: largeion = .false.   ! true: use large ion category
    integer :: nonigrd = 0
    integer :: ftauopt = 1
    real    :: rgard1  = 1.0
    real    :: rarfac  = 1.0
    real    :: trever = -15.0
    real    :: cidiamin = 0.0  ! minimum crystal size for charging
    integer, private :: iremoveqwfrz = 1    ! Whether to remove (=1) or not (=0) the newly-frozen cloud droplets (ibfc=1) from the CWC used for charge separation

    real    :: erbnd  = 0.01
    real    :: fdgt   = 1.0
    real    :: costhe = 0.35
    real    :: delqnsa = -50.0e-14
    real    :: delqxsa =  50.0e-14
    real    :: delqnsb = -2.0e-14
    real    :: delqxsb =  2.0e-14
    real    :: delqnia = -50.0e-14
    real    :: delqxia =  50.0e-14
    real, private    :: delqnw = -1.0e-10!-1.0e-12 !
    real, private    :: delqxw =  1.0e-10! 1.0e-12 !
    real    :: scxacymax = 3000.e-12 ! Maximum charging rate magnitude (C s-1 m-3)
    real   , private :: ecollmx = 0.5 ! Maximum collision efficiency for graup/hail with ice; used only for charging rates
    real    :: tindmn = 233, tindmx = 298.0  ! min and max temperatures where inductive charging is allowed
    real gf4p5, gf4ds, gf4br


    ! Takahashi lookup table
    integer, parameter :: nlwc=30, ntem=31
    real takalu(0:ntem,0:nlwc)

    !  gamma function lookup table
    integer ngm0,ngm1,ngm2
    parameter (ngm0=3001,ngm1=500,ngm2=500)
    double precision, parameter :: dgam = 0.01, dgami = 100.
    double precision gmoi(0:ngm0) ! ,gmod(0:ngm1,0:ngm2),gmdi(0:ngm1,0:ngm2)

    ! Maximum collision efficiency
    real :: ehs_collsn = 0.5, ehi_collsn = 1.0, esi_collsn = 1.0

    ! "an" indices
    integer :: na = 9
    integer :: ltmp

    integer, parameter :: icespheres = 0 ! turn ice spheres (frozen droplets) on (1) or off (0). NOT COMPLETE IN WRF/ARPS/CM1 CODE!
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
    !   integer :: lschab = 0

    !   integer :: lscb = 0
    !   integer :: lsce = 0
    !   integer :: lsceq = 0

    integer :: lne = 0 ! last varible for transforming

    integer lsc(lc:lqmx)
    integer ln(lc:lqmx)
    integer denscale(lc:lqmx) ! flag for density scaling (mixing ratio conversion)
    ! End "an" indices

    real qxmin(lc:lqmx)

    !---------------------------------------------------------------------------
    !End Explicit lightning shared vars-----------------------------------------

CONTAINS

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE MORR_TWO_MOMENT_INIT(hail_opt,nssl_params,idonictmp) ! RAS
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! THIS SUBROUTINE INITIALIZES ALL PHYSICAL CONSTANTS AMND PARAMETERS
        ! NEEDED BY THE MICROPHYSICS SCHEME.
        ! NEEDS TO BE CALLED AT FIRST TIME STEP, PRIOR TO CALL TO MAIN MICROPHYSICS INTERFACE
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        IMPLICIT NONE

        INTEGER, INTENT(IN):: hail_opt ! RAS

        integer n,i
        !---------------------------------------------------------------------------
        !Explicit lightning --------------------------------------------------------
        real, intent(in), dimension(20) :: nssl_params ! for lightning
        logical, optional, intent(in) :: idonictmp
        integer :: igam
        double precision :: arg
        integer :: j
        !---------------------------------------------------------------------------
        !End explicit lightning --------------------------------------------------------

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        ! THE FOLLOWING PARAMETERS ARE USER-DEFINED SWITCHES AND NEED TO BE
        ! SET PRIOR TO CODE COMPILATION

        ! INUM IS AUTOMATICALLY SET TO 0 FOR WRF-CHEM BELOW,
        ! ALLOWING PREDICTION OF DROPLET CONCENTRATION
        ! THUS, THIS PARAMETER SHOULD NOT BE CHANGED HERE
        ! AND SHOULD BE LEFT TO 1

        !---------------------------------------------------------------------------
        !Explicit lightning --------------------------------------------------------

        ipelec   = Nint(nssl_params(11))
        isaund   = Nint(nssl_params(12))
        IF ( present(idonictmp) ) THEN
            idonic = idonictmp
        ENDIF

        ! Test print ipelec
        !      print *, "ipelec, isaund: ", ipelec, isaund

        !      Index calculation (corresponds to nssl's "ipconc == 5")
        lhab = lh ! 7
        lccn = lhab+1 ! 8
        lnc = lhab+2 ! 9
        lnr = lhab+3 ! 10
        lni = lhab+4 ! 11
        lns = lhab+5 ! 12
        lnh = lhab+6 ! 13
        ltmp = lnh

        IF ( ipelec > 0 .or. idonic ) THEN
            lscw  = ltmp+1 ! 14
            lscr  = ltmp+2 ! 15
            lsci  = ltmp+3 ! 16
            lscs  = ltmp+4 ! 17
            lsch  = ltmp+5 ! 18
            ltmp = lsch

            lscpi = ltmp+1 ! 19
            lscni = ltmp+2 ! 20
            ltmp = lscni
        ENDIF
        na = ltmp ! 20

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

        ! fill in the complete gamma function lookup table
        gmoi(0) = 1.d32
        do igam = 1,ngm0
            arg = dgam*igam
            gmoi(igam) = gamma_dp(arg)
        end do

        ! Define qxmin from nssl (later changed in micro)
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

        !---------------------------------------------------------------------------
        !End Explicit lightning ----------------------------------------------------

        INUM = 1

        ! SET CONSTANT DROPLET CONCENTRATION (UNITS OF CM-3)
        ! IF NO COUPLING WITH WRF-CHEM

        NDCNST = 250.

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! NOTE, THE FOLLOWING OPTIONS RELATED TO DROPLET ACTIVATION
        ! (IACT, IBASE, ISUB) ARE NOT AVAILABLE IN CURRENT VERSION
        ! FOR WRF-CHEM, DROPLET ACTIVATION IS PERFORMED
        ! IN 'MIX_ACTIVATE', NOT IN MICROPHYSICS SCHEME


        ! IACT = 1, USE POWER-LAW CCN SPECTRA, NCCN = CS^K
        ! IACT = 2, USE LOGNORMAL AEROSOL SIZE DIST TO DERIVE CCN SPECTRA

        IACT = 2

        ! IBASE = 1, NEGLECT DROPLET ACTIVATION AT LATERAL CLOUD EDGES DUE TO
        !             UNRESOLVED ENTRAINMENT AND MIXING, ACTIVATE
        !             AT CLOUD BASE OR IN REGION WITH LITTLE CLOUD WATER USING
        !             NON-EQULIBRIUM SUPERSATURATION ASSUMING NO INITIAL CLOUD WATER,
        !             IN CLOUD INTERIOR ACTIVATE USING EQUILIBRIUM SUPERSATURATION
        ! IBASE = 2, ASSUME DROPLET ACTIVATION AT LATERAL CLOUD EDGES DUE TO
        !             UNRESOLVED ENTRAINMENT AND MIXING DOMINATES,
        !             ACTIVATE DROPLETS EVERYWHERE IN THE CLOUD USING NON-EQUILIBRIUM
        !             SUPERSATURATION ASSUMING NO INITIAL CLOUD WATER, BASED ON THE
        !             LOCAL SUB-GRID AND/OR GRID-SCALE VERTICAL VELOCITY
        !             AT THE GRID POINT

        ! NOTE: ONLY USED FOR PREDICTED DROPLET CONCENTRATION (INUM = 0)

        IBASE = 2

        ! INCLUDE SUB-GRID VERTICAL VELOCITY (standard deviation of w) IN DROPLET ACTIVATION
        ! ISUB = 0, INCLUDE SUB-GRID W (RECOMMENDED FOR LOWER RESOLUTION)
        ! currently, sub-grid w is constant of 0.5 m/s (not coupled with PBL/turbulence scheme)
        ! ISUB = 1, EXCLUDE SUB-GRID W, ONLY USE GRID-SCALE W

        ! NOTE: ONLY USED FOR PREDICTED DROPLET CONCENTRATION (INUM = 0)

        ISUB = 0
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


        ! SWITCH FOR LIQUID-ONLY RUN
        ! ILIQ = 0, INCLUDE ICE
        ! ILIQ = 1, LIQUID ONLY, NO ICE

        ILIQ = 0

        ! SWITCH FOR ICE NUCLEATION
        ! INUC = 0, USE FORMULA FROM RASMUSSEN ET AL. 2002 (MID-LATITUDE)
        !      = 1, USE MPACE OBSERVATIONS (ARCTIC ONLY)

        INUC = 0

        ! SWITCH FOR GRAUPEL/HAIL NO GRAUPEL/HAIL
        ! IGRAUP = 0, INCLUDE GRAUPEL/HAIL
        ! IGRAUP = 1, NO GRAUPEL/HAIL

        IGRAUP = 0

        ! HM ADDED 11/7/07
        ! SWITCH FOR HAIL/GRAUPEL
        ! IHAIL = 0, DENSE PRECIPITATING ICE IS GRAUPEL
        ! IHAIL = 1, DENSE PRECIPITATING ICE IS HAIL
        ! NOTE ---> RECOMMEND IHAIL = 1 FOR CONTINENTAL DEEP CONVECTION

        !IHAIL = 0 !changed to namelist option (hail_opt) by RAS
        ! Check if namelist option is feasible, otherwise default to graupel - RAS
        IF (hail_opt .eq. 1) THEN
            IHAIL = 1
        ELSE
            IHAIL = 0
        ENDIF

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! SET PHYSICAL CONSTANTS

        ! FALLSPEED PARAMETERS (V=AD^B)
        AI = 700.
        AC = 3.E7
        AS = 11.72
        AR = 841.99667
        BI = 1.
        BC = 2.
        BS = 0.41
        BR = 0.8
        IF (IHAIL.EQ.0) THEN
            AG = 19.3
            BG = 0.37
        ELSE ! (MATSUN AND HUGGINS 1980)
            AG = 114.5
            BG = 0.5
        END IF

        ! CONSTANTS AND PARAMETERS
        !         R = 287.15
        !         RV = 461.5
        !         CP = 1005.
        RHOSU = 85000./(287.15*273.15)
        RHOW = 997.
        RHOI = 500.
        RHOSN = 100.
        IF (IHAIL.EQ.0) THEN
            RHOG = 400.
            !            ! Test increase rhog
            !            RHOG = 800.
        ELSE
            RHOG = 900.
        END IF
        AIMM = 0.66
        BIMM = 100.
        ECR = 1.
        DCS = 125.E-6
        MI0 = 4./3.*PI*RHOI*(10.E-6)**3
        MG0 = 1.6E-10
        F1S = 0.86
        F2S = 0.28
        F1R = 0.78
        !         F2R = 0.32
        ! fix 053011
        F2R = 0.308
        !         G = 9.806
        QSMALL = 1.E-14
        EII = 0.1
        ECI = 0.7
        ! HM, ADD FOR V3.2
        ! hm, 7/23/13
        !         CPW = 4218.
        CPW = 4187.

        ! SIZE DISTRIBUTION PARAMETERS

        CI = RHOI*PI/6.
        DI = 3.
        CS = RHOSN*PI/6.
        DS = 3.
        CG = RHOG*PI/6.
        DG = 3.

        ! RADIUS OF CONTACT NUCLEI
        RIN = 0.1E-6

        MMULT = 4./3.*PI*RHOI*(5.E-6)**3

        ! SIZE LIMITS FOR LAMBDA

        LAMMAXI = 1./1.E-6
        LAMMINI = 1./(2.*DCS+100.E-6)
        LAMMAXR = 1./20.E-6
        !         LAMMINR = 1./500.E-6
        LAMMINR = 1./2800.E-6

        LAMMAXS = 1./10.E-6
        LAMMINS = 1./2000.E-6
        LAMMAXG = 1./20.E-6
        LAMMING = 1./2000.E-6

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! note: these parameters only used by the non-wrf-chem version of the
        !       scheme with predicted droplet number

        ! CCN SPECTRA FOR IACT = 1

        ! MARITIME
        ! MODIFIED FROM RASMUSSEN ET AL. 2002
        ! NCCN = C*S^K, NCCN IS IN CM-3, S IS SUPERSATURATION RATIO IN %

        K1 = 0.4
        C1 = 120.

        ! CONTINENTAL

        !              K1 = 0.5
        !              C1 = 1000.

        ! AEROSOL ACTIVATION PARAMETERS FOR IACT = 2
        ! PARAMETERS CURRENTLY SET FOR AMMONIUM SULFATE

        MW = 0.018
        OSM = 1.
        VI = 3.
        EPSM = 0.7
        RHOA = 1777.
        MAP = 0.132
        MA = 0.0284
        ! hm fix 6/23/16
        !         RR = 8.3187
        RR = 8.3145
        BACT = VI*OSM*EPSM*MW*RHOA/(MAP*RHOW)

        ! AEROSOL SIZE DISTRIBUTION PARAMETERS CURRENTLY SET FOR MPACE
        ! (see morrison et al. 2007, JGR)
        ! MODE 1

        RM1 = 0.052E-6
        SIG1 = 2.04
        NANEW1 = 72.2E6
        F11 = 0.5*EXP(2.5*(LOG(SIG1))**2)
        F21 = 1.+0.25*LOG(SIG1)

        ! MODE 2

        RM2 = 1.3E-6
        SIG2 = 2.5
        NANEW2 = 1.8E6
        F12 = 0.5*EXP(2.5*(LOG(SIG2))**2)
        F22 = 1.+0.25*LOG(SIG2)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        ! CONSTANTS FOR EFFICIENCY

        CONS1=GAMMA(1.+DS)*CS
        CONS2=GAMMA(1.+DG)*CG
        CONS3=GAMMA(4.+BS)/6.
        CONS4=GAMMA(4.+BR)/6.
        CONS5=GAMMA(1.+BS)
        CONS6=GAMMA(1.+BR)
        CONS7=GAMMA(4.+BG)/6.
        CONS8=GAMMA(1.+BG)
        CONS9=GAMMA(5./2.+BR/2.)
        CONS10=GAMMA(5./2.+BS/2.)
        CONS11=GAMMA(5./2.+BG/2.)
        CONS12=GAMMA(1.+DI)*CI
        CONS13=GAMMA(BS+3.)*PI/4.*ECI
        CONS14=GAMMA(BG+3.)*PI/4.*ECI
        CONS15=-1108.*EII*PI**((1.-BS)/3.)*RHOSN**((-2.-BS)/3.)/(4.*720.)
        CONS16=GAMMA(BI+3.)*PI/4.*ECI
        CONS17=4.*2.*3.*RHOSU*PI*ECI*ECI*GAMMA(2.*BS+2.)/(8.*(RHOG-RHOSN))
        CONS18=RHOSN*RHOSN
        CONS19=RHOW*RHOW
        CONS20=20.*PI*PI*RHOW*BIMM
        CONS21=4./(DCS*RHOI)
        CONS22=PI*RHOI*DCS**3/6.
        CONS23=PI/4.*EII*GAMMA(BS+3.)
        CONS24=PI/4.*ECR*GAMMA(BR+3.)
        CONS25=PI*PI/24.*RHOW*ECR*GAMMA(BR+6.)
        CONS26=PI/6.*RHOW
        CONS27=GAMMA(1.+BI)
        CONS28=GAMMA(4.+BI)/6.
        CONS29=4./3.*PI*RHOW*(25.E-6)**3
        CONS30=4./3.*PI*RHOW
        CONS31=PI*PI*ECR*RHOSN
        CONS32=PI/2.*ECR
        CONS33=PI*PI*ECR*RHOG
        CONS34=5./2.+BR/2.
        CONS35=5./2.+BS/2.
        CONS36=5./2.+BG/2.
        CONS37=4.*PI*1.38E-23/(6.*PI*RIN)
        CONS38=PI*PI/3.*RHOW
        CONS39=PI*PI/36.*RHOW*BIMM
        CONS40=PI/6.*BIMM
        CONS41=PI*PI*ECR*RHOW

        !+---+-----------------------------------------------------------------+
        !..Set these variables needed for computing radar reflectivity.  These
        !.. get used within radar_init to create other variables used in the
        !.. radar module.

        xam_r = PI*RHOW/6.
        xbm_r = 3.
        xmu_r = 0.
        xam_s = CS
        xbm_s = DS
        xmu_s = 0.
        xam_g = CG
        xbm_g = DG
        xmu_g = 0.

        call radar_init
        !+---+-----------------------------------------------------------------+


    END SUBROUTINE MORR_TWO_MOMENT_INIT

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! THIS SUBROUTINE IS MAIN INTERFACE WITH THE TWO-MOMENT MICROPHYSICS SCHEME
    ! THIS INTERFACE TAKES IN 3D VARIABLES FROM DRIVER MODEL, CONVERTS TO 1D FOR
    ! CALL TO THE MAIN MICROPHYSICS SUBROUTINE (SUBROUTINE MORR_TWO_MOMENT_MICRO)
    ! WHICH OPERATES ON 1D VERTICAL COLUMNS.
    ! 1D VARIABLES FROM THE MAIN MICROPHYSICS SUBROUTINE ARE THEN REASSIGNED BACK TO 3D FOR OUTPUT
    ! BACK TO DRIVER MODEL USING THIS INTERFACE.
    ! MICROPHYSICS TENDENCIES ARE ADDED TO VARIABLES HERE BEFORE BEING PASSED BACK TO DRIVER MODEL.

    ! THIS CODE WAS WRITTEN BY HUGH MORRISON (NCAR) AND SLAVA TATARSKII (GEORGIA TECH).

    ! FOR QUESTIONS, CONTACT: HUGH MORRISON, E-MAIL: MORRISON@UCAR.EDU, PHONE:303-497-8916

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    SUBROUTINE MP_MORR_TWO_MOMENT(ITIMESTEP,                       &
            TH, QV, QC, QR, QI, QS, QG, NI, NS, NR, NG, &
            RHO, PII, P, DT_IN, DZ, HT, W,          &
            RAINNC, RAINNCV, SR,                    &
            SNOWNC,SNOWNCV,GRAUPELNC,GRAUPELNCV,    & ! hm added 7/13/13
            refl_10cm, diagflag, do_radar_ref,      & ! GT added for reflectivity calcs
            qrcuten, qscuten, qicuten               & ! hm added
            ,F_QNDROP, qndrop                        & ! hm added, wrf-chem
            ,IDS,IDE, JDS,JDE, KDS,KDE               & ! domain dims
            ,IMS,IME, JMS,JME, KMS,KME               & ! memory dims
            ,ITS,ITE, JTS,JTE, KTS,KTE               & ! tile   dims            )
            !jdf		   ,C2PREC3D,CSED3D,ISED3D,SSED3D,GSED3D,RSED3D & ! HM ADD, WRF-CHEM
            ,rainprod, evapprod                      &
            ,QLSINK,PRECR,PRECI,PRECS,PRECG &        ! HM ADD, WRF-CHEM
            !-----------------------------------------------------------------------
            !Explicit lightning ----------------------------------------------------Add Arguments
            !                     ,cn                       & ! modified from "CCW,"
            !                     ,CCW                      & ! modified from "CCW,"
            ,SCR                       & !
            ,SCW                       & !
            ,SCI                       & ! P
            ,SCS                       & ! U
            ,SCH                       & ! R
            !                     ,SCHL                      & ! E
            ,SCTOT                     & !
            ,NONINDUC                  & ! A
            ,rscghis_2d                & ! D
            ,INDUC                     & ! D
            ,ELEC                      & ! E
            ,SCION                     & ! D
            ,SCIONA                    & !
            ,dx                        & !
            ,dy                        & !
            ,ipelectmp                 & !
            ,csaci03d                  & !
            ,chaci03d                  & !
            ,chacs03d                  & !
            ,csaci03d_morr             & !
            ,umc,umr,umi,ums,umg       & !
            ,QGSTEN3D,QRSTEN3D,QISTEN3D,QNISTEN3D,QCSTEN3D &
            ,PRD3d &
            ,PRE3d &
            ,PRDS3d &
            ,MNUCCC3d &
            ,PRA3d &
            ,PRC3d &
            ,PCC3d &
            ,MNUCCD3d &
            ,MNUCCR3d &
            ,PRAI3d &
            ,PRCI3d &
            ,PSACWS3d &
            ,PSACWI3d &
            ,QMULTS3d &
            ,QMULTR3d &
            ,PRACS3d &
            ,PCCN3d &
            ,PSMLT3d &
            ,EVPMS3d &
            ,PIACR3d &
            ,PRACI3d &
            ,PIACRS3d &
            ,PRACIS3d &
            ,EPRD3d &
            ,EPRDS3d &
            ,PRACG3d &
            ,PSACWG3d &
            ,PGSACW3d &
            ,PGRACS3d &
            ,PRDG3d &
            ,EPRDG3d &
            ,EVPMG3d &
            ,PGMLT3d &
            ,PSACR3d &
            ,QMULTG3d &
            ,QMULTRG3d &
            !-----------------------------------------------------------------------
            !End Explicit lightning ------------------------------------------------
            )

        ! QV - water vapor mixing ratio (kg/kg)
        ! QC - cloud water mixing ratio (kg/kg)
        ! QR - rain water mixing ratio (kg/kg)
        ! QI - cloud ice mixing ratio (kg/kg)
        ! QS - snow mixing ratio (kg/kg)
        ! QG - graupel mixing ratio (KG/KG)
        ! NI - cloud ice number concentration (1/kg)
        ! NS - Snow Number concentration (1/kg)
        ! NR - Rain Number concentration (1/kg)
        ! NG - Graupel number concentration (1/kg)
        ! NOTE: RHO AND HT NOT USED BY THIS SCHEME AND DO NOT NEED TO BE PASSED INTO SCHEME!!!!
        ! P - AIR PRESSURE (PA)
        ! W - VERTICAL AIR VELOCITY (M/S)
        ! TH - POTENTIAL TEMPERATURE (K)
        ! PII - exner function - used to convert potential temp to temp
        ! DZ - difference in height over interface (m)
        ! DT_IN - model time step (sec)
        ! ITIMESTEP - time step counter
        ! RAINNC - accumulated grid-scale precipitation (mm)
        ! RAINNCV - one time step grid scale precipitation (mm/time step)
        ! SNOWNC - accumulated grid-scale snow plus cloud ice (mm)
        ! SNOWNCV - one time step grid scale snow plus cloud ice (mm/time step)
        ! GRAUPELNC - accumulated grid-scale graupel (mm)
        ! GRAUPELNCV - one time step grid scale graupel (mm/time step)
        ! SR - one time step mass ratio of snow to total precip
        ! qrcuten, rain tendency from parameterized cumulus convection
        ! qscuten, snow tendency from parameterized cumulus convection
        ! qicuten, cloud ice tendency from parameterized cumulus convection

        ! variables below currently not in use, not coupled to PBL or radiation codes
        ! TKE - turbulence kinetic energy (m^2 s-2), NEEDED FOR DROPLET ACTIVATION (SEE CODE BELOW)
        ! NCTEND - droplet concentration tendency from pbl (kg-1 s-1)
        ! NCTEND - CLOUD ICE concentration tendency from pbl (kg-1 s-1)
        ! KZH - heat eddy diffusion coefficient from YSU scheme (M^2 S-1), NEEDED FOR DROPLET ACTIVATION (SEE CODE BELOW)
        ! EFFCS - CLOUD DROPLET EFFECTIVE RADIUS OUTPUT TO RADIATION CODE (micron)
        ! EFFIS - CLOUD DROPLET EFFECTIVE RADIUS OUTPUT TO RADIATION CODE (micron)
        ! HM, ADDED FOR WRF-CHEM COUPLING
        ! QLSINK - TENDENCY OF CLOUD WATER TO RAIN, SNOW, GRAUPEL (KG/KG/S)
        ! CSED,ISED,SSED,GSED,RSED - SEDIMENTATION FLUXES (KG/M^2/S) FOR CLOUD WATER, ICE, SNOW, GRAUPEL, RAIN
        ! PRECI,PRECS,PRECG,PRECR - SEDIMENTATION FLUXES (KG/M^2/S) FOR ICE, SNOW, GRAUPEL, RAIN

        ! rainprod - total tendency of conversion of cloud water/ice and graupel to rain (kg kg-1 s-1)
        ! evapprod - tendency of evaporation of rain (kg kg-1 s-1)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        ! reflectivity currently not included!!!!
        ! REFL_10CM - CALCULATED RADAR REFLECTIVITY AT 10 CM (DBZ)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        ! EFFC - DROPLET EFFECTIVE RADIUS (MICRON)
        ! EFFR - RAIN EFFECTIVE RADIUS (MICRON)
        ! EFFS - SNOW EFFECTIVE RADIUS (MICRON)
        ! EFFI - CLOUD ICE EFFECTIVE RADIUS (MICRON)

        ! ADDITIONAL OUTPUT FROM MICRO - SEDIMENTATION TENDENCIES, NEEDED FOR LIQUID-ICE STATIC ENERGY

        ! QGSTEN - GRAUPEL SEDIMENTATION TEND (KG/KG/S)
        ! QRSTEN - RAIN SEDIMENTATION TEND (KG/KG/S)
        ! QISTEN - CLOUD ICE SEDIMENTATION TEND (KG/KG/S)
        ! QNISTEN - SNOW SEDIMENTATION TEND (KG/KG/S)
        ! QCSTEN - CLOUD WATER SEDIMENTATION TEND (KG/KG/S)

        ! WVAR - STANDARD DEVIATION OF SUB-GRID VERTICAL VELOCITY (M/S)

        !-----------------------------------------------------------------------
        !Explicit lightning ----------------------------------------------------
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
        !-----------------------------------------------------------------------
        !End Explicit lightning ------------------------------------------------



#if ( defined(DM_PARALLEL)  &&   ! defined(STUBMPI) )
#define MPI
        USE module_dm, ONLY : &
                local_communicator, mytask
        ! keep a spacing line here to keep Apple cpp from adding a space in front of the endif
#endif

        IMPLICIT NONE

#if ( defined(DM_PARALLEL)  &&   ! defined(STUBMPI) ) || defined(MPI)
        INCLUDE 'mpif.h'
#else
      integer :: mytask = 0
#endif

        INTEGER,      INTENT(IN   )    ::   ids, ide, jds, jde, kds, kde , &
                ims, ime, jms, jme, kms, kme , &
                its, ite, jts, jte, kts, kte
        ! Temporary changed from INOUT to IN

        REAL, DIMENSION(ims:ime, kms:kme, jms:jme), INTENT(INOUT):: &
                qv, qc, qr, qi, qs, qg, ni, ns, nr, TH, NG
        !jdf                      qndrop ! hm added, wrf-chem
        REAL, DIMENSION(ims:ime, kms:kme, jms:jme), optional,INTENT(INOUT):: qndrop
        !jdf  REAL, DIMENSION(ims:ime, kms:kme, jms:jme),INTENT(INOUT):: CSED3D, &
        REAL, DIMENSION(ims:ime, kms:kme, jms:jme), optional,INTENT(INOUT):: QLSINK, &
                rainprod, evapprod, &
                PRECI,PRECS,PRECG,PRECR ! HM, WRF-CHEM
        !, effcs, effis

        REAL, DIMENSION(ims:ime, kms:kme, jms:jme), INTENT(IN):: &
                pii, p, dz, rho, w !, tke, nctend, nitend,kzh
        REAL, INTENT(IN):: dt_in
        INTEGER, INTENT(IN):: ITIMESTEP

        REAL, DIMENSION(ims:ime, jms:jme), INTENT(INOUT):: &
                RAINNC, RAINNCV, SR, &
                ! hm added 7/13/13
                SNOWNC,SNOWNCV,GRAUPELNC,GRAUPELNCV

        REAL, DIMENSION(ims:ime, kms:kme, jms:jme), INTENT(INOUT)::       &  ! GT
                refl_10cm

        REAL , DIMENSION( ims:ime , jms:jme ) , INTENT(IN) ::       ht

        ! LOCAL VARIABLES

        REAL, DIMENSION(its:ite, kts:kte, jts:jte)::                     &
                effi, effs, effr, EFFG

        REAL, DIMENSION(its:ite, kts:kte, jts:jte)::                     &
                T, WVAR, EFFC

        REAL, DIMENSION(kts:kte) ::                                                                &
                QC_TEND1D, QI_TEND1D, QNI_TEND1D, QR_TEND1D,                      &
                NI_TEND1D, NS_TEND1D, NR_TEND1D,                                  &
                QC1D, QI1D, QR1D,NI1D, NS1D, NR1D, QS1D,                          &
                T_TEND1D,QV_TEND1D, T1D, QV1D, P1D, W1D, WVAR1D,         &
                EFFC1D, EFFI1D, EFFS1D, EFFR1D,DZ1D,   &
                ! HM ADD GRAUPEL
                QG_TEND1D, NG_TEND1D, QG1D, NG1D, EFFG1D, &

                ! ADD SEDIMENTATION TENDENCIES (UNITS OF KG/KG/S)
                QGSTEN,QRSTEN, QISTEN, QNISTEN, QCSTEN, &
                ! ADD CUMULUS TENDENCIES
                QRCU1D, QSCU1D, QICU1D

        ! add cumulus tendencies

        REAL, DIMENSION(ims:ime, kms:kme, jms:jme), INTENT(IN):: &
                qrcuten, qscuten, qicuten

        LOGICAL, INTENT(IN), OPTIONAL ::                F_QNDROP  ! wrf-chem
        LOGICAL :: flag_qndrop  ! wrf-chem
        integer :: iinum ! wrf-chem

        ! wrf-chem
        REAL, DIMENSION(kts:kte) :: nc1d, nc_tend1d,C2PREC,CSED,ISED,SSED,GSED,RSED
        REAL, DIMENSION(kts:kte) :: rainprod1d, evapprod1d
        ! HM add reflectivity
        REAL, DIMENSION(kts:kte) :: dBZ

        REAL PRECPRT1D, SNOWRT1D, SNOWPRT1D, GRPLPRT1D ! hm added 7/13/13

        INTEGER I,K,J

        REAL DT

        LOGICAL, OPTIONAL, INTENT(IN) :: diagflag
        INTEGER, OPTIONAL, INTENT(IN) :: do_radar_ref

        !-----------------------------------------------------------------------
        !Explicit lightning ----------------------------------------------------Define vars
        real, dimension(ims:ime, kms:kme, jms:jme),  optional, intent(inout)::  &
                scr,scw,sci,scs,sch,sciona,sctot,induc,noninduc
        real, dimension(ims:ime, kms:kme, jms:jme,2),optional, intent(inout) :: scion
        real, dimension(ims:ime, jms:jme), optional, intent(inout):: rscghis_2d
        real, dimension(ims:ime, kms:kme, jms:jme),  optional, intent(inout):: csaci03d, chaci03d, chacs03d, csaci03d_morr ! to be passed out [#/kg/s]
        real, dimension(ims:ime, kms:kme, jms:jme),  optional, intent(inout):: umc, umr, umi, ums, umg ! fall speed
        real, dimension(ims:ime, kms:kme, jms:jme),  optional, intent(inout):: QGSTEN3D,QRSTEN3D,QISTEN3D,QNISTEN3D,QCSTEN3D ! sed rates
        real, dimension(ims:ime, kms:kme, jms:jme),  optional, intent(in) :: elec ! elecsave = Ez
        real, optional, intent(in) :: dx,dy
        integer, optional, intent(in) :: ipelectmp

        real, dimension(ims:ime, kms:kme, jms:jme),  optional, intent(inout):: PRD3d
        real, dimension(ims:ime, kms:kme, jms:jme),  optional, intent(inout):: PRE3d
        real, dimension(ims:ime, kms:kme, jms:jme),  optional, intent(inout):: PRDS3d
        real, dimension(ims:ime, kms:kme, jms:jme),  optional, intent(inout):: MNUCCC3d
        real, dimension(ims:ime, kms:kme, jms:jme),  optional, intent(inout):: PRA3d
        real, dimension(ims:ime, kms:kme, jms:jme),  optional, intent(inout):: PRC3d
        real, dimension(ims:ime, kms:kme, jms:jme),  optional, intent(inout):: PCC3d
        real, dimension(ims:ime, kms:kme, jms:jme),  optional, intent(inout):: MNUCCD3d
        real, dimension(ims:ime, kms:kme, jms:jme),  optional, intent(inout):: MNUCCR3d
        real, dimension(ims:ime, kms:kme, jms:jme),  optional, intent(inout):: PRAI3d
        real, dimension(ims:ime, kms:kme, jms:jme),  optional, intent(inout):: PRCI3d
        real, dimension(ims:ime, kms:kme, jms:jme),  optional, intent(inout):: PSACWS3d
        real, dimension(ims:ime, kms:kme, jms:jme),  optional, intent(inout):: PSACWI3d
        real, dimension(ims:ime, kms:kme, jms:jme),  optional, intent(inout):: QMULTS3d
        real, dimension(ims:ime, kms:kme, jms:jme),  optional, intent(inout):: QMULTR3d
        real, dimension(ims:ime, kms:kme, jms:jme),  optional, intent(inout):: PRACS3d
        real, dimension(ims:ime, kms:kme, jms:jme),  optional, intent(inout):: PCCN3d
        real, dimension(ims:ime, kms:kme, jms:jme),  optional, intent(inout):: PSMLT3d
        real, dimension(ims:ime, kms:kme, jms:jme),  optional, intent(inout):: EVPMS3d
        real, dimension(ims:ime, kms:kme, jms:jme),  optional, intent(inout):: PIACR3d
        real, dimension(ims:ime, kms:kme, jms:jme),  optional, intent(inout):: PRACI3d
        real, dimension(ims:ime, kms:kme, jms:jme),  optional, intent(inout):: PIACRS3d
        real, dimension(ims:ime, kms:kme, jms:jme),  optional, intent(inout):: PRACIS3d
        real, dimension(ims:ime, kms:kme, jms:jme),  optional, intent(inout):: EPRD3d
        real, dimension(ims:ime, kms:kme, jms:jme),  optional, intent(inout):: EPRDS3d
        real, dimension(ims:ime, kms:kme, jms:jme),  optional, intent(inout):: PRACG3d
        real, dimension(ims:ime, kms:kme, jms:jme),  optional, intent(inout):: PSACWG3d
        real, dimension(ims:ime, kms:kme, jms:jme),  optional, intent(inout):: PGSACW3d
        real, dimension(ims:ime, kms:kme, jms:jme),  optional, intent(inout):: PGRACS3d
        real, dimension(ims:ime, kms:kme, jms:jme),  optional, intent(inout):: PRDG3d
        real, dimension(ims:ime, kms:kme, jms:jme),  optional, intent(inout):: EPRDG3d
        real, dimension(ims:ime, kms:kme, jms:jme),  optional, intent(inout):: EVPMG3d
        real, dimension(ims:ime, kms:kme, jms:jme),  optional, intent(inout):: PGMLT3d
        real, dimension(ims:ime, kms:kme, jms:jme),  optional, intent(inout):: PSACR3d
        real, dimension(ims:ime, kms:kme, jms:jme),  optional, intent(inout):: QMULTG3d
        real, dimension(ims:ime, kms:kme, jms:jme),  optional, intent(inout):: QMULTRG3d

        !        Local vars.
        real, dimension(kts:kte) :: elec1d ! ez = elecsave bar
        real, dimension(its:ite, 1, kts:kte) :: sciona2d,dz2d,rho2d

        REAL, DIMENSION(its:ite, kts:kte, jts:jte) :: rho_igl ! air density calculated using ideal gas law
        real, dimension(its:ite, 1, kts:kte) :: t0,t1,t2,t3,t4,t5,t6,t7,t8,t9
        real, dimension(kts:kte) :: t2_1d,t3_1d,t5_1d,t6_1d,t8_1d
        !        REAL, DIMENSION(its:ite, kts:kte, jts:jte, 2, lc:lhab) :: fx ! =1: mass weigted, =2: number weighted
        !                fr, fnr, fi, fni, fs, fns, fc, fnc, fg, fng   ! terminal fall velocity
        REAL, DIMENSION(kts:kte) ::                                                                &
                !                fr1d, fnr1d, fi1d, fni1d, fs1d, fns1d, fc1d, fnc1d, fg1d, fng1d   & ! terminal fall velocity
                umc1d, umr1d, umi1d, ums1d, umg1d ! 1d fall speed, pass inout micro
        double precision :: chgneg1,chgpos1,chgneg2,chgpos2,chgneg3,chgpos3
        double precision :: chgiona1,chgiona2,chgiona3,sctot3,scwtot,scrtot,scitot,scstot,schtot,schltot
        double precision :: scwmax,scrmax,scimax,scsmax,schmax,schlmax
        double precision :: scwmin,scrmin,scimin,scsmin,schmin,schlmin
        logical :: f_cnatmp

        real, dimension(its:ite, 1, kts:kte, na) :: an
        real, dimension(kts:kte, na) :: an1d
        real, dimension(kts:kte) :: to_sed
        real, dimension(kts:kte) :: csaci0, chaci0, chacs0, csaci0_morr ! local 1d array [#/m3/s]

        real, dimension(kts:kte) :: PRD
        real, dimension(kts:kte) :: PRE
        real, dimension(kts:kte) :: PRDS
        real, dimension(kts:kte) :: MNUCCC
        real, dimension(kts:kte) :: PRA
        real, dimension(kts:kte) :: PRC
        real, dimension(kts:kte) :: PCC
        real, dimension(kts:kte) :: MNUCCD
        real, dimension(kts:kte) :: MNUCCR
        real, dimension(kts:kte) :: PRAI
        real, dimension(kts:kte) :: PRCI
        real, dimension(kts:kte) :: PSACWS
        real, dimension(kts:kte) :: PSACWI
        real, dimension(kts:kte) :: QMULTS
        real, dimension(kts:kte) :: QMULTR
        real, dimension(kts:kte) :: PRACS
        real, dimension(kts:kte) :: PCCN
        real, dimension(kts:kte) :: PSMLT
        real, dimension(kts:kte) :: EVPMS
        real, dimension(kts:kte) :: PIACR
        real, dimension(kts:kte) :: PRACI
        real, dimension(kts:kte) :: PIACRS
        real, dimension(kts:kte) :: PRACIS
        real, dimension(kts:kte) :: EPRD
        real, dimension(kts:kte) :: EPRDS
        real, dimension(kts:kte) :: PRACG
        real, dimension(kts:kte) :: PSACWG
        real, dimension(kts:kte) :: PGSACW
        real, dimension(kts:kte) :: PGRACS
        real, dimension(kts:kte) :: PRDG
        real, dimension(kts:kte) :: EPRDG
        real, dimension(kts:kte) :: EVPMG
        real, dimension(kts:kte) :: PGMLT
        real, dimension(kts:kte) :: PSACR
        real, dimension(kts:kte) :: QMULTG
        real, dimension(kts:kte) :: QMULTRG

        integer :: lnb, il
        integer :: nx,ny,nz
        integer :: nstep, ltrue
        integer, dimension(its:ite, jts:jte) :: nstep2d, ltrue2d
        integer, parameter :: nor = 0
        real :: tmp, dv
        integer ix,jy,kz

#ifdef MPI

#if defined(MPI)
        integer, parameter :: ntot = 50
        double precision  mpitotindp(ntot), mpitotoutdp(ntot)
        INTEGER :: mpi_error_code = 1
#endif
#endif
        !-----------------------------------------------------------------------
        !End Explicit lightning ------------------------------------------------


        ! below for wrf-chem
        flag_qndrop = .false.
        IF ( PRESENT ( f_qndrop ) ) flag_qndrop = f_qndrop
        !!!!!!!!!!!!!!!!!!!!!!

        ! Initialize tendencies (all set to 0) and transfer
        ! array to local variables
        DT = DT_IN

        DO I=ITS,ITE
            DO J=JTS,JTE
                DO K=KTS,KTE
                    T(I,K,J)        = TH(i,k,j)*PII(i,k,j)

                    ! NOTE: WVAR NOT CURRENTLY USED IN CODE !!!!!!!!!!
                    ! currently assign wvar to 0.5 m/s (not coupled with PBL scheme)

                    WVAR(I,K,J)     = 0.5
                    rho_igl(i,k,j) = P(i,k,j)/(R*T(i,k,j))

                    ! currently mixing of number concentrations also is neglected (not coupled with PBL schemes)

                END DO
            END DO
        END DO

        !-----------------------------------------------------------------------
        !Explicit lightning ----------------------------------------------------Init vars
        IF ( present( ipelectmp ) ) THEN
            ipelec = ipelectmp
        ELSE
            ipelec = 0
        ENDIF

        if (ipelec>0) then
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

        t2_1d(:) = 0.0
        t3_1d(:) = 0.0
        t5_1d(:) = 0.0
        t6_1d(:) = 0.0
        t8_1d(:) = 0.0

        PRD3d = 0.0
        PRE3d = 0.0
        PRDS3d = 0.0
        MNUCCC3d = 0.0
        PRA3d = 0.0
        PRC3d = 0.0
        PCC3d = 0.0
        MNUCCD3d = 0.0
        MNUCCR3d = 0.0
        PRAI3d = 0.0
        PRCI3d = 0.0
        PSACWS3d = 0.0
        PSACWI3d = 0.0
        QMULTS3d = 0.0
        QMULTR3d = 0.0
        PRACS3d = 0.0
        PCCN3d = 0.0
        PSMLT3d = 0.0
        EVPMS3d = 0.0
        PIACR3d = 0.0
        PRACI3d = 0.0
        PIACRS3d = 0.0
        PRACIS3d = 0.0
        EPRD3d = 0.0
        EPRDS3d = 0.0
        PRACG3d = 0.0
        PSACWG3d = 0.0
        PGSACW3d = 0.0
        PGRACS3d = 0.0
        PRDG3d = 0.0
        EPRDG3d = 0.0
        EVPMG3d = 0.0
        PGMLT3d = 0.0
        PSACR3d = 0.0
        QMULTG3d = 0.0
        QMULTRG3d = 0.0


            print*, "-----------------------------------------------------------"
            print*, "at it = ", itimestep, "dt = ", dt_in, " pre print"
            print*, "ims:ime, kms:kme, jms:jme"
            print*, ims,ime, kms,kme, jms,jme
            print*, "its:ite, kts:kte, jts:jte"
            print*, its,ite, kts,kte, jts,jte

            ! print*, "max of scw = ", maxval(scw), "at ", maxloc(scw)
            ! print*, "max of scr = ", maxval(scr), "at ", maxloc(scr)
            ! print*, "max of sci = ", maxval(sci), "at ", maxloc(sci)
            ! print*, "max of scs = ", maxval(scs), "at ", maxloc(scs)
            ! print*, "max of sch = ", maxval(sch), "at ", maxloc(sch)

        endif
        !-----------------------------------------------------------------------
        !End Explicit lightning ------------------------------------------------

        do j=jts,jte      ! j loop (north-south)
            !-----------------------------------------------------------------------
            !Explicit lightning ---------------------------------------------------- init for the slab

            IF ( present(scw) .and. ipelec > 0 ) THEN
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                !!!!        Construct "an" array BEFORED UPDATED BY MORR_2MOM_MICRO
                ! only for space charge, others done in MORR_2MOM_MICRO

                an(:,1,:,:) = 0.0
                an(its:ite,1,kts:kte,lt)   = th(its:ite,kts:kte,j)
                !
                !            !    mixing ratio (kg/kg)
                !            an(its:ite,1,kts:kte,lv)   = qv(its:ite,kts:kte,j)
                !            an(its:ite,1,kts:kte,lc)   = qc(its:ite,kts:kte,j)
                !            an(its:ite,1,kts:kte,lr)   = qr(its:ite,kts:kte,j)
                !            an(its:ite,1,kts:kte,li)   = qi(its:ite,kts:kte,j)
                !            an(its:ite,1,kts:kte,ls)   = qs(its:ite,kts:kte,j)
                !            an(its:ite,1,kts:kte,lh)   = qg(its:ite,kts:kte,j)
                !
                !            !    number concentration (#/kg)
                !            an(its:ite,1,kts:kte,lnc)  = qndrop(its:ite,kts:kte,j)
                !            an(its:ite,1,kts:kte,lnr)  = nr(its:ite,kts:kte,j)
                !            an(its:ite,1,kts:kte,lni)  = ni(its:ite,kts:kte,j)
                !            an(its:ite,1,kts:kte,lns)  = ns(its:ite,kts:kte,j)
                !            an(its:ite,1,kts:kte,lnh)  = ng(its:ite,kts:kte,j)
                !            !        IF ( lhl > 1 ) an(:,1,:,lnhl) = nhl(:,:,j)

                !!!!  need to add ccn
                !                an(:,1,:,lccn) = 0.0

                ! for use in nssl_sc
                dz2d (its:ite,1,kts:kte) = dz (its:ite,kts:kte,j)
                rho2d(its:ite,1,kts:kte) = rho_igl(its:ite,kts:kte,j)

                !     space charge (C/kg)
                an(its:ite,1,kts:kte,lscw)   = scw(its:ite,kts:kte,j)
                an(its:ite,1,kts:kte,lscr)   = scr(its:ite,kts:kte,j)
                an(its:ite,1,kts:kte,lsci)   = sci(its:ite,kts:kte,j)
                an(its:ite,1,kts:kte,lscs)   = scs(its:ite,kts:kte,j)
                an(its:ite,1,kts:kte,lsch)   = sch(its:ite,kts:kte,j)
                !        IF ( lhl > 1 ) an(:,1,:,lschl)   = schl(:,:,j)
                an(its:ite,1,kts:kte,lscpi) = Max(0.0, sciona(its:ite,kts:kte,j)) * ec_i * rho_igl(its:ite,kts:kte,j) ! Convert sciona[C/kg] -> an[#/m3]
                an(its:ite,1,kts:kte,lscni) = -Min(0.0,sciona(its:ite,kts:kte,j)) * ec_i * rho_igl(its:ite,kts:kte,j)
                sciona2d(its:ite,1,kts:kte) = sciona(its:ite,kts:kte,j) * rho_igl(its:ite,kts:kte,j)

                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                !!!  Collect into chg...1 (an[kg-1])
                DO k = kts,kte
                    DO i = its,ite
                        ! Fake ions for now -- will need to change this if real ions ever get used in WRF
                        ! convert net ion density to positive or negative ions, so that underlying code is
                        ! still ion-friendly
                        tmp =  rho_igl(i,k,j) * (an(i,1,k,lscw)+ &
                                an(i,1,k,lscr)+ &
                                an(i,1,k,lsci)+ &
                                an(i,1,k,lscs)+ &
                                an(i,1,k,lsch)+ &
                                ! an(i,1,k,lschl)+ &
                                sciona(i,k,j)) *dx*dy*dz(i,k,j)  ! [kg/m3]*[C/kg]*[m3] = [C]
                        chgneg1 = chgneg1 + Min(0.0, tmp)
                        chgpos1 = chgpos1 + Max(0.0, tmp)
                        chgiona1 = chgiona1 + (an(i,1,k,lscpi)-an(i,1,k,lscni)) * ec * rho_igl(i,k,j) * dx*dy*dz(i,k,j) ! rho_igl looks like a mistake here
                    END DO
                END DO
                !!!  End chg...1
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                !!!!       "denscale" list: =1, need to convert from kg-1 -> m-3
                !!!!       Here includes lsc(c,r,i,s,h) only. # concen are converted in micro.
                denscale(:) = 0
                denscale(lccn:lsch) = 1
                lnb = lhab + 1 ! lccn
                DO il = lscw,na
                    IF ( denscale(il) == 1 ) THEN ! true for lsc(w,r,i,s,h), no ion
                        an(its:ite,1,kts:kte,il) = an(its:ite,1,kts:kte,il) * rho_igl(its:ite,kts:kte,j)
                    ENDIF
                ENDDO ! il
                !!!!      End "denscale"
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            END IF
            !---------------------------------------------------------------------------
            !End Explicit lightning ----------------------------------------------------

            do i=its,ite      ! i loop (east-west)
                !
                ! Transfer 3D arrays into 1D for microphysical calculations
                !

                ! hm , initialize 1d tendency arrays to zero

                do k=kts,kte   ! k loop (vertical)

                    QC_TEND1D(k)  = 0.
                    QI_TEND1D(k)  = 0.
                    QNI_TEND1D(k) = 0.
                    QR_TEND1D(k)  = 0.
                    NI_TEND1D(k)  = 0.
                    NS_TEND1D(k)  = 0.
                    NR_TEND1D(k)  = 0.
                    T_TEND1D(k)   = 0.
                    QV_TEND1D(k)  = 0.
                    nc_tend1d(k) = 0. ! wrf-chem

                    QC1D(k)       = QC(i,k,j)
                    QI1D(k)       = QI(i,k,j)
                    QS1D(k)       = QS(i,k,j)
                    QR1D(k)       = QR(i,k,j)

                    NI1D(k)       = NI(i,k,j)

                    NS1D(k)       = NS(i,k,j)
                    NR1D(k)       = NR(i,k,j)
                    ! HM ADD GRAUPEL
                    QG1D(K)       = QG(I,K,j)
                    NG1D(K)       = NG(I,K,j)
                    QG_TEND1D(K)  = 0.
                    NG_TEND1D(K)  = 0.

                    T1D(k)        = T(i,k,j)
                    QV1D(k)       = QV(i,k,j)
                    P1D(k)        = P(i,k,j)
                    DZ1D(k)       = DZ(i,k,j)
                    W1D(k)        = W(i,k,j)
                    WVAR1D(k)     = WVAR(i,k,j)
                    ! add cumulus tendencies, already decoupled
                    qrcu1d(k)     = qrcuten(i,k,j)
                    qscu1d(k)     = qscuten(i,k,j)
                    qicu1d(k)     = qicuten(i,k,j)
                end do  !jdf added this
                ! below for wrf-chem
                IF (flag_qndrop .AND. PRESENT( qndrop )) THEN
                    iact = 3
                    DO k = kts, kte
                        nc1d(k)=qndrop(i,k,j)
                        iinum=0
                    ENDDO
                ELSE
                    DO k = kts, kte
                        nc1d(k)=0. ! temporary placeholder, set to constant in microphysics subroutine
                        iinum=1
                    ENDDO
                ENDIF

                !jdf  end do

                !-----------------------------------------------------------------------
                !Explicit lightning ----------------------------------------------------
                if (ipelec>0) then
                an1d(kts:kte,1:na) = an(i, 1, kts:kte, 1:na)
                if (ipelec>1) then
                    elec1d(kts:kte)   =  elec(i,kts:kte,j) ! Ez here =elecsave
                end if
                CHACI0(:) = 0.0
                CSACI0(:) = 0.0
                CHACS0(:) = 0.0
                csaci0_morr(:) = 0.0
                umc1d(:) = 0.0
                umr1d(:) = 0.0
                umi1d(:) = 0.0
                ums1d(:) = 0.0
                umg1d(:) = 0.0
                end if
                !End Explicit lightning ------------------------------------------------
                !-----------------------------------------------------------------------
                call MORR_TWO_MOMENT_MICRO(QC_TEND1D, QI_TEND1D, QNI_TEND1D, QR_TEND1D,            &
                NI_TEND1D, NS_TEND1D, NR_TEND1D,                                                  &
                QC1D, QI1D, QS1D, QR1D,NI1D, NS1D, NR1D,                                          &
                T_TEND1D,QV_TEND1D, T1D, QV1D, P1D, DZ1D, W1D, WVAR1D,                   &
                PRECPRT1D,SNOWRT1D,                                                               &
                SNOWPRT1D,GRPLPRT1D,                 & ! hm added 7/13/13
                EFFC1D,EFFI1D,EFFS1D,EFFR1D,DT,                                                   &
                IMS,IME, JMS,JME, KMS,KME,                   &
                ITS,ITE, JTS,JTE, KTS,KTE,                   & ! HM ADD GRAUPEL
                QG_TEND1D,NG_TEND1D,QG1D,NG1D,EFFG1D, &
                qrcu1d, qscu1d, qicu1d, &
                ! ADD SEDIMENTATION TENDENCIES
                QGSTEN,QRSTEN,QISTEN,QNISTEN,QCSTEN, &
                nc1d, nc_tend1d, iinum, C2PREC,CSED,ISED,SSED,GSED,RSED & !wrf-chem
                !-----------------------------------------------------------------------
                !Explicit lightning ----------------------------------------------------
                !                , fr1d, fnr1d, fi1d, fni1d, fs1d, fns1d, fc1d, fnc1d, fg1d, fng1d  & ! terminal fall velocity
                , umc1d, umr1d, umi1d, ums1d, umg1d & ! fall speed
                , nstep, ipelec, ltrue &
                , I, J, ITIMESTEP &
                , IDS,IDE, JDS,JDE, KDS,KDE & ! domain dims
                , CSACI0, CHACI0, CHACS0, csaci0_morr &
                , an1d, elec1d &
                , t2_1d,t3_1d,t5_1d,t6_1d,t8_1d &
                ,PRD &
                ,PRE &
                ,PRDS &
                ,MNUCCC &
                ,PRA &
                ,PRC &
                ,PCC &
                ,MNUCCD &
                ,MNUCCR &
                ,PRAI &
                ,PRCI &
                ,PSACWS &
                ,PSACWI &
                ,QMULTS &
                ,QMULTR &
                ,PRACS &
                ,PCCN &
                ,PSMLT &
                ,EVPMS &
                ,PIACR &
                ,PRACI &
                ,PIACRS &
                ,PRACIS &
                ,EPRD &
                ,EPRDS &
                ,PRACG &
                ,PSACWG &
                ,PGSACW &
                ,PGRACS &
                ,PRDG &
                ,EPRDG &
                ,EVPMG &
                ,PGMLT &
                ,PSACR &
                ,QMULTG &
                ,QMULTRG &
                !-----------------------------------------------------------------------
                !End Explicit lightning ------------------------------------------------
#if (WRF_CHEM == 1)
                ,rainprod1d, evapprod1d & !wrf-chem
#endif
                )

                !
                ! Transfer 1D arrays back into 3D arrays
                !

                !-----------------------------------------------------------------------
                !Explicit lightning ----------------------------------------------------
                if (present(scw) .and. ipelec>0) then
                    umc(i,kts:kte,j) = umc1d(kts:kte)
                    umr(i,kts:kte,j) = umr1d(kts:kte)
                    umi(i,kts:kte,j) = umi1d(kts:kte)
                    ums(i,kts:kte,j) = ums1d(kts:kte)
                    umg(i,kts:kte,j) = umg1d(kts:kte)

                    csaci03d(i,kts:kte,j) = csaci0(kts:kte)/rho_igl(i,kts:kte,j)
                    !                    csaci03d(i,kts:kte,j) = csaci0(kts:kte)
                    chaci03d(i,kts:kte,j) = chaci0(kts:kte)/rho_igl(i,kts:kte,j)
                    chacs03d(i,kts:kte,j) = chacs0(kts:kte)/rho_igl(i,kts:kte,j)
                    csaci03d_morr(i,kts:kte,j) = csaci0_morr(kts:kte)

                    an(i, 1, kts:kte, 1:na) = an1d(kts:kte,1:na)

                    if (ltrue) then
                        nstep2d(i,j) = nstep
                    else
                        nstep2d(i,j) = 0
                    end if
                    ltrue2d (i,j) = ltrue

                    t2(i,1, kts:kte) = t2_1d(kts:kte)
                    t3(i,1, kts:kte) = t3_1d(kts:kte)
                    t5(i,1, kts:kte) = t5_1d(kts:kte)
                    t6(i,1, kts:kte) = t6_1d(kts:kte)
                    t8(i,1, kts:kte) = t8_1d(kts:kte)

                    QGSTEN3D(i,kts:kte,j) = QGSTEN(kts:kte)
                    QRSTEN3D(i,kts:kte,j) = QRSTEN(kts:kte)
                    QISTEN3D(i,kts:kte,j) = QISTEN(kts:kte)
                    QNISTEN3D(i,kts:kte,j) = QNISTEN(kts:kte)
                    QCSTEN3D(i,kts:kte,j) = QCSTEN(kts:kte)

                    PRD3d(i,kts:kte,j) = PRD(kts:kte)
                    PRE3d(i,kts:kte,j) = PRE(kts:kte)
                    PRDS3d(i,kts:kte,j) = PRDS(kts:kte)
                    MNUCCC3d(i,kts:kte,j) = MNUCCC(kts:kte)
                    PRA3d(i,kts:kte,j) = PRA(kts:kte)
                    PRC3d(i,kts:kte,j) = PRC(kts:kte)
                    PCC3d(i,kts:kte,j) = PCC(kts:kte)
                    MNUCCD3d(i,kts:kte,j) = MNUCCD(kts:kte)
                    MNUCCR3d(i,kts:kte,j) = MNUCCR(kts:kte)
                    PRAI3d(i,kts:kte,j) = PRAI(kts:kte)
                    PRCI3d(i,kts:kte,j) = PRCI(kts:kte)
                    PSACWS3d(i,kts:kte,j) = PSACWS(kts:kte)
                    PSACWI3d(i,kts:kte,j) = PSACWI(kts:kte)
                    QMULTS3d(i,kts:kte,j) = QMULTS(kts:kte)
                    QMULTR3d(i,kts:kte,j) = QMULTR(kts:kte)
                    PRACS3d(i,kts:kte,j) = PRACS(kts:kte)
                    PCCN3d(i,kts:kte,j) = PCCN(kts:kte)
                    PSMLT3d(i,kts:kte,j) = PSMLT(kts:kte)
                    EVPMS3d(i,kts:kte,j) = EVPMS(kts:kte)
                    PIACR3d(i,kts:kte,j) = PIACR(kts:kte)
                    PRACI3d(i,kts:kte,j) = PRACI(kts:kte)
                    PIACRS3d(i,kts:kte,j) = PIACRS(kts:kte)
                    PRACIS3d(i,kts:kte,j) = PRACIS(kts:kte)
                    EPRD3d(i,kts:kte,j) = EPRD(kts:kte)
                    EPRDS3d(i,kts:kte,j) = EPRDS(kts:kte)
                    PRACG3d(i,kts:kte,j) = PRACG(kts:kte)
                    PSACWG3d(i,kts:kte,j) = PSACWG(kts:kte)
                    PGSACW3d(i,kts:kte,j) = PGSACW(kts:kte)
                    PGRACS3d(i,kts:kte,j) = PGRACS(kts:kte)
                    PRDG3d(i,kts:kte,j) = PRDG(kts:kte)
                    EPRDG3d(i,kts:kte,j) = EPRDG(kts:kte)
                    EVPMG3d(i,kts:kte,j) = EVPMG(kts:kte)
                    PGMLT3d(i,kts:kte,j) = PGMLT(kts:kte)
                    PSACR3d(i,kts:kte,j) = PSACR(kts:kte)
                    QMULTG3d(i,kts:kte,j) = QMULTG(kts:kte)
                    QMULTRG3d(i,kts:kte,j) = QMULTRG(kts:kte)

                end if
                !-----------------------------------------------------------------------
                !End Explicit lightning ------------------------------------------------

                do k=kts,kte

                    ! hm, add tendencies to update global variables
                    ! HM, TENDENCIES FOR Q AND N NOW ADDED IN M2005MICRO, SO WE
                    ! ONLY NEED TO TRANSFER 1D VARIABLES BACK TO 3D

                    QC(i,k,j)        = QC1D(k)
                    QI(i,k,j)        = QI1D(k)
                    QS(i,k,j)        = QS1D(k)
                    QR(i,k,j)        = QR1D(k)
                    NI(i,k,j)        = NI1D(k)
                    NS(i,k,j)        = NS1D(k)
                    NR(i,k,j)        = NR1D(k)
                    QG(I,K,j)        = QG1D(K)
                    NG(I,K,j)        = NG1D(K)

                    T(i,k,j)         = T1D(k)
                    TH(I,K,J)        = T(i,k,j)/PII(i,k,j) ! CONVERT TEMP BACK TO POTENTIAL TEMP
                    QV(i,k,j)        = QV1D(k)

                    EFFC(i,k,j)      = EFFC1D(k)
                    EFFI(i,k,j)      = EFFI1D(k)
                    EFFS(i,k,j)      = EFFS1D(k)
                    EFFR(i,k,j)      = EFFR1D(k)
                    EFFG(I,K,j)      = EFFG1D(K)

                    ! wrf-chem
                    IF (flag_qndrop .AND. PRESENT( qndrop )) THEN
                        qndrop(i,k,j) = nc1d(k)
                        !jdf         CSED3D(I,K,J) = CSED(K)
                    END IF
                    IF ( PRESENT( QLSINK ) ) THEN
                        if(qc(i,k,j)>1.e-10) then
                            QLSINK(I,K,J)  = C2PREC(K)/QC(I,K,J)
                        else
                            QLSINK(I,K,J)  = 0.0
                        endif
                    END IF
                    IF ( PRESENT( PRECR ) ) PRECR(I,K,J) = RSED(K)
                    IF ( PRESENT( PRECI ) ) PRECI(I,K,J) = ISED(K)
                    IF ( PRESENT( PRECS ) ) PRECS(I,K,J) = SSED(K)
                    IF ( PRESENT( PRECG ) ) PRECG(I,K,J) = GSED(K)
                    ! EFFECTIVE RADIUS FOR RADIATION CODE (currently not coupled)
                    ! HM, ADD LIMIT TO PREVENT BLOWING UP OPTICAL PROPERTIES, 8/18/07
                    !          EFFCS(I,K,J)     = MIN(EFFC(I,K,J),50.)
                    !          EFFCS(I,K,J)     = MAX(EFFCS(I,K,J),1.)
                    !          EFFIS(I,K,J)     = MIN(EFFI(I,K,J),130.)
                    !          EFFIS(I,K,J)     = MAX(EFFIS(I,K,J),13.)

#if ( WRF_CHEM == 1)
                    IF ( PRESENT( rainprod ) ) rainprod(i,k,j) = rainprod1d(k)
                    IF ( PRESENT( evapprod ) ) evapprod(i,k,j) = evapprod1d(k)
#endif

                end do

                ! hm modified so that m2005 precip variables correctly match wrf precip variables
                RAINNC(i,j) = RAINNC(I,J)+PRECPRT1D
                RAINNCV(i,j) = PRECPRT1D
                ! hm, added 7/13/13
                SNOWNC(i,j) = SNOWNC(I,J)+SNOWPRT1D
                SNOWNCV(i,j) = SNOWPRT1D
                GRAUPELNC(i,j) = GRAUPELNC(I,J)+GRPLPRT1D
                GRAUPELNCV(i,j) = GRPLPRT1D
                SR(i,j) = SNOWRT1D/(PRECPRT1D+1.E-12)

                !+---+-----------------------------------------------------------------+
                IF ( PRESENT (diagflag) ) THEN
                    if (diagflag .and. do_radar_ref == 1) then
                        call refl10cm_hm (qv1d, qr1d, nr1d, qs1d, ns1d, qg1d, ng1d,   &
                                t1d, p1d, dBZ, kts, kte, i, j)
                        do k = kts, kte
                            refl_10cm(i,k,j) = MAX(-35., dBZ(k))
                        enddo
                    endif
                ENDIF
                !+---+-----------------------------------------------------------------+

            end do

            !-----------------------------------------------------------------------
            !Explicit lightning ---------------------------------------------------- remaining processes after micro (2D)

            nx = ite-its+1
            ny = 1         ! set up as 2D slabs
            nz = kte-kts+1

            IF ( cleardiag ) THEN
                IF ( present( rscghis_2d ) ) THEN
                    rscghis_2d(its:ite,j) = 0.0
                ENDIF
            ENDIF


            IF ( present( noninduc ) .and. ipelec >= 1 ) THEN
                DO k = kts,kte
                    DO i = its,ite
                        induc(i,k,j)  = t5(i,1,k) ! schacw
                        noninduc(i,k,j) = (t2(i,1,k) + t3(i,1,k) + t6(i,1,k) + t8(i,1,k)) ! schaci + schacs
                        IF ( .not. ( noninduc(i,k,j) > -1.e-3 .and. noninduc(i,k,j) < 1.e-3 ) ) THEN
                            write(0,*) 'problem with noninduc! i,j,k,noninduc = ',i,j,k,noninduc(i,k,j)
                            STOP
                        ENDIF
                        IF ( .not. ( induc(i,k,j) > -1.e-3 .and. induc(i,k,j) < 1.e-3 ) ) THEN
                            write(0,*) 'problem with induc! i,j,k,induc = ',i,j,k,induc(i,k,j),t5(i,1,k)
                            !                            STOP
                        ENDIF
                        IF ( present( rscghis_2d ) ) THEN
                            t9(i,1,1) = t9(i,1,1) + Abs( noninduc(i,k,j) )*dz(i,k,j)
                        ENDIF
                    ENDDO
                ENDDO
                IF ( present( rscghis_2d ) ) THEN
                    DO i = its,ite
                        rscghis_2d(i,j) = Max( rscghis_2d(i,j), t9(i,1,1) )
                    ENDDO
                ENDIF
            ENDIF

            IF ( present( sciona ) .and. lnh > 1 .and. ipelec > 0) THEN

                CALL SCIONATT    &
                        &  (nx,ny,nz,na,j &
                        &  ,nor,nor &
                        &  ,dt,dz2d &
                        &  ,an &
                        &  )
            ENDIF

            IF ( present(scw) .and. ipelec > 0 .and. lscw > 1 ) THEN
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                !!!  Collect into chg...2 (an[m-3])
                DO k = kts,kte
                    DO i = its,ite
                        !                IF ( lhl > 1 ) THEN
                        !                    tmp =   (an(i,1,k,lscw)+an(i,1,k,lscr)+an(i,1,k,lsci)+ &
                        !                            &             an(i,1,k,lscs)+an(i,1,k,lsch)+an(i,1,k,lschl)+ &
                        !                            &             (an(i,1,k,lscpi)-an(i,1,k,lscni))*ec )*dx*dy*dz(i,k,j)
                        !                ELSE
                        tmp =   (an(i,1,k,lscw)+ &
                                an(i,1,k,lscr)+ &
                                an(i,1,k,lsci)+ &
                                an(i,1,k,lscs)+ &
                                an(i,1,k,lsch)+ &
                                (an(i,1,k,lscpi)-an(i,1,k,lscni))*ec) *dx*dy*dz(i,k,j)
                        !                ENDIF
                        chgneg2 = chgneg2 + Min(0.0, tmp)
                        chgpos2 = chgpos2 + Max(0.0, tmp)
                        chgiona2 = chgiona2 + (an(i,1,k,lscpi)-an(i,1,k,lscni)) * ec * rho_igl(i,k,j) * dx*dy*dz(i,k,j)
                    ENDDO
                ENDDO
                !!!  End chg...2
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                !!!!       "denscale" list: =1, need to convert from m-3 back to kg-1
                !!!!       Here includes ln(c,r,i,s,h), lsc(c,r,i,s,h)
                DO il = lnb,na
                    IF ( denscale(il) == 1 ) THEN
                        an(its:ite,1,kts:kte,il) = an(its:ite,1,kts:kte,il) / rho_igl(its:ite,kts:kte,j)
                    ENDIF
                ENDDO ! il
                !!!!      End "denscale" list
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                sciona(its:ite,kts:kte,j) = (an(ix,1,kz,lscpi) - an(ix,1,kz,lscni))&
                        * ec / rho_igl(its:ite,kts:kte,j)
                scw(its:ite,kts:kte,j) = an(its:ite,1,kts:kte,lscw)
                scr(its:ite,kts:kte,j) = an(its:ite,1,kts:kte,lscr)
                sci(its:ite,kts:kte,j) = an(its:ite,1,kts:kte,lsci)
                scs(its:ite,kts:kte,j) = an(its:ite,1,kts:kte,lscs)
                sch(its:ite,kts:kte,j) = an(its:ite,1,kts:kte,lsch)
                sctot(its:ite,kts:kte,j)  = rho_igl(its:ite,kts:kte,j)*(scw(its:ite,kts:kte,j) &
                        +scr(its:ite,kts:kte,j) &
                        +sci(its:ite,kts:kte,j) &
                        +scs(its:ite,kts:kte,j) &
                        +sch(its:ite,kts:kte,j) &
                        +sciona(its:ite,kts:kte,j) ) ! to get C/m3 for RHS calc later on

                DO i = its,ite
                    IF ( i.le.(ids+2) .or. i.ge.(ide-2)  ) THEN
                        IF ( j.le.(jds+2) .or. j.ge.(jde-2) ) THEN
                            scw(i,kts:kte,j) = 0.
                            scr(i,kts:kte,j) = 0.
                            sci(i,kts:kte,j) = 0.
                            scs(i,kts:kte,j) = 0.
                            sch(i,kts:kte,j) = 0.
                            sctot(i,kts:kte,j) = 0.
                        END IF
                    end if
                end do


                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                !!!  Collect into chg...3 (an[kg-1])
                DO k = kts,kte
                    DO i = its,ite
                        tmp =  rho_igl(i,k,j)*(an(i,1,k,lscw)+ &
                                an(i,1,k,lscr)+ &
                                an(i,1,k,lsci)+ &
                                an(i,1,k,lscs)+ &
                                an(i,1,k,lsch)+ &
                                sciona(i,k,j)) *dx*dy*dz(i,k,j)
                        chgneg3 = chgneg3 + Min(0.0, tmp)
                        chgpos3 = chgpos3 + Max(0.0, tmp)
                        chgiona3 = chgiona3 + sciona(i,k,j)*rho_igl(i,k,j)*dx*dy*dz(i,k,j)
                    END DO
                end do
                !!!  End chg...3
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                dv = rho_igl(i,k,j)*dx*dy*dz(i,k,j)

                IF ( present(elec) ) THEN
                    sctot3 = sctot3 + sctot(i,k,j)*dx*dy*dz(i,k,j)
                    scwtot  = scwtot + scw(i,k,j)*dv
                    scrtot  = scrtot + scr(i,k,j)*dv
                    scitot  = scitot + sci(i,k,j)*dv
                    scstot  = scstot + scs(i,k,j)*dv
                    schtot  = schtot + sch(i,k,j)*dv
                    !                    schltot = schltot + schl(i,k,j)*dv

                    scwmax = Max( scwmax, scw(i,k,j) )
                    scrmax = Max( scrmax, scr(i,k,j) )
                    scimax = Max( scimax, sci(i,k,j) )
                    scsmax = Max( scsmax, scs(i,k,j) )
                    schmax = Max( schmax, sch(i,k,j) )
                    !                    IF ( lhl > 1 ) schlmax= Max(schlmax, schl(i,k,j) )

                    scwmin = Min( scwmin, scw(i,k,j) )
                    scrmin = Min( scrmin, scr(i,k,j) )
                    scimin = Min( scimin, sci(i,k,j) )
                    scsmin = Min( scsmin, scs(i,k,j) )
                    schmin = Min( schmin, sch(i,k,j) )
                    !                    IF ( lhl > 1 ) schlmin= Min(schlmin, schl(i,k,j) )

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
                   IF ( present(scw) .and. ipelec > 0 .and. mytask == 0 ) THEN !{
#else
                    IF ( present(scw) .and. ipelec > 0 ) THEN !{
#endif
                    ENDIF
                ENDIF
            ENDIF ! ipelec
            !-----------------------------------------------------------------------
            !End Explicit lightning ------------------------------------------------

        end do !j

    END SUBROUTINE MP_MORR_TWO_MOMENT

    SUBROUTINE SC_SED1D(sc1d, fv, dzq, nstep, dt                     &
            ,ITS,ITE, JTS,JTE, KTS,KTE&
            ,i,j,itimestep)
        !        Purpose: Sedimentation of space charge

        implicit none

        real, INTENT(IN)  :: dt
        integer, INTENT(IN)  :: nstep
        integer, INTENT(IN)  :: ITS,ITE, JTS,JTE, KTS,KTE
        integer, intent(in):: i,j,itimestep
        real, INTENT(INOUT) :: sc1d (kts:kte) ! space charge (C/m3)
        real, INTENT(IN) :: fv (kts:kte)  ! terminal fall velocity (m/s)
        real, INTENT(IN) :: dzq (kts:kte) ! difference in height across level (m)

        ! local vars
        integer :: N, k
        real :: scfalout(kts:kte) ! fallout flux (C/m2/s)
        real :: scfaltnd(kts:kte) ! net fallout tendency (C/m3/s)
        real :: scsten(kts:kte) ! sedimentation tendency (C/m3/s)
        real :: sc1dori(kts:kte)

        scfalout(kts:kte) = 0
        scfaltnd(kts:kte) = 0
        scsten(kts:kte)   = 0

        DO N = 1,NSTEP
            scfalout(kts:kte) = fv(kts:kte) * sc1d(kts:kte)

            k = kte
            scfaltnd(k) = scfalout(k)/dzq(k)
            scsten(k)   = -scfaltnd(k)/NSTEP

            do k = kte-1, kts, -1
                scfaltnd(k) = (scfalout(k+1)-scfalout(k))/dzq(k)
                scsten(k)   = scfaltnd(k)/NSTEP
            end do
        END DO

        sc1dori(kts:kte) = sc1d(kts:kte)
        sc1d(kts:kte) = sc1d(kts:kte) + dt*scsten(kts:kte)
        if (maxval(abs(sc1d)).ge.1.) then
            write(0,*) 'sed exceeded at i,j,it=', i,j,itimestep
            write(0,*) 'sc1d=',sc1d
            write(0,*) 'sc1dori=',sc1dori
            write(0,*) 'fv=',fv
            write(0,*) 'scfalout=',scfalout
            write(0,*) 'dzq=',dzq
            write(0,*) 'scfaltnd=',scfaltnd
            write(0,*) 'nstep=',nstep
            write(0,*) 'scsten=',scsten
            write(0,*) 'dt=',dt
            !            stop
        end if
    END SUBROUTINE SC_SED1D

    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    SUBROUTINE MORR_TWO_MOMENT_MICRO(QC3DTEN,QI3DTEN,QNI3DTEN,QR3DTEN,         &
            NI3DTEN,NS3DTEN,NR3DTEN,QC3D,QI3D,QNI3D,QR3D,NI3D,NS3D,NR3D,              &
            T3DTEN,QV3DTEN,T3D,QV3D,PRES,DZQ,W3D,WVAR,PRECRT,SNOWRT,            &
            SNOWPRT,GRPLPRT,                & ! hm added 7/13/13
            EFFC,EFFI,EFFS,EFFR,DT,                                                   &
            IMS,IME, JMS,JME, KMS,KME,           &
            ITS,ITE, JTS,JTE, KTS,KTE,           & ! ADD GRAUPEL
            QG3DTEN,NG3DTEN,QG3D,NG3D,EFFG,qrcu1d,qscu1d, qicu1d,    &
            QGSTEN,QRSTEN,QISTEN,QNISTEN,QCSTEN, &
            nc3d,nc3dten,iinum, & ! wrf-chem
            c2prec,CSED,ISED,SSED,GSED,RSED  &  ! hm added, wrf-chem
            !            , fr, fnr, fi, fni, fs, fns, fc, fnc, fg, fng  & ! terminal fall velocity
            , umc1d, umr1d, umi1d, ums1d, umg1d & ! fall speed
            , nstep, ipelec, ltrue &
            , I, J, ITIMESTEP &
            , IDS,IDE, JDS,JDE, KDS,KDE & ! domain dims
            , CSACI0, CHACI0, CHACS0, csaci0_morr &
            , an1d, elec1d &
            , t2,t3,t5,t6,t8 &
            ,PRD &
            ,PRE &
            ,PRDS &
            ,MNUCCC &
            ,PRA &
            ,PRC &
            ,PCC &
            ,MNUCCD &
            ,MNUCCR &
            ,PRAI &
            ,PRCI &
            ,PSACWS &
            ,PSACWI &
            ,QMULTS &
            ,QMULTR &
            ,PRACS &
            ,PCCN &
            ,PSMLT &
            ,EVPMS &
            ,PIACR &
            ,PRACI &
            ,PIACRS &
            ,PRACIS &
            ,EPRD &
            ,EPRDS &
            ,PRACG &
            ,PSACWG &
            ,PGSACW &
            ,PGRACS &
            ,PRDG &
            ,EPRDG &
            ,EVPMG &
            ,PGMLT &
            ,PSACR &
            ,QMULTG &
            ,QMULTRG &
#if (WRF_CHEM == 1)
            ,rainprod, evapprod &
#endif
            )

    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    ! THIS PROGRAM IS THE MAIN TWO-MOMENT MICROPHYSICS SUBROUTINE DESCRIBED BY
    ! MORRISON ET AL. 2005 JAS AND MORRISON ET AL. 2009 MWR

    ! THIS SCHEME IS A BULK DOUBLE-MOMENT SCHEME THAT PREDICTS MIXING
    ! RATIOS AND NUMBER CONCENTRATIONS OF FIVE HYDROMETEOR SPECIES:
    ! CLOUD DROPLETS, CLOUD (SMALL) ICE, RAIN, SNOW, AND GRAUPEL/HAIL.

    ! CODE STRUCTURE: MAIN SUBROUTINE IS 'MORR_TWO_MOMENT'. ALSO INCLUDED IN THIS FILE IS
    ! 'FUNCTION POLYSVP', 'FUNCTION DERF1', AND
    ! 'FUNCTION GAMMA'.

    ! NOTE: THIS SUBROUTINE USES 1D ARRAY IN VERTICAL (COLUMN), EVEN THOUGH VARIABLES ARE CALLED '3D'......

    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

    ! DECLARATIONS

    IMPLICIT NONE

    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    ! THESE VARIABLES BELOW MUST BE LINKED WITH THE MAIN MODEL.
    ! DEFINE ARRAY SIZES

    ! INPUT NUMBER OF GRID CELLS

    ! INPUT/OUTPUT PARAMETERS                                 ! DESCRIPTION (UNITS)
    INTEGER, INTENT( IN)  :: IMS,IME, JMS,JME, KMS,KME,          &
            ITS,ITE, JTS,JTE, KTS,KTE

    REAL, DIMENSION(KTS:KTE) ::  QC3DTEN            ! CLOUD WATER MIXING RATIO TENDENCY (KG/KG/S)
    REAL, DIMENSION(KTS:KTE) ::  QI3DTEN            ! CLOUD ICE MIXING RATIO TENDENCY (KG/KG/S)
    REAL, DIMENSION(KTS:KTE) ::  QNI3DTEN           ! SNOW MIXING RATIO TENDENCY (KG/KG/S)
    REAL, DIMENSION(KTS:KTE) ::  QR3DTEN            ! RAIN MIXING RATIO TENDENCY (KG/KG/S)
    REAL, DIMENSION(KTS:KTE) ::  NI3DTEN            ! CLOUD ICE NUMBER CONCENTRATION (1/KG/S)
    REAL, DIMENSION(KTS:KTE) ::  NS3DTEN            ! SNOW NUMBER CONCENTRATION (1/KG/S)
    REAL, DIMENSION(KTS:KTE) ::  NR3DTEN            ! RAIN NUMBER CONCENTRATION (1/KG/S)
    REAL, DIMENSION(KTS:KTE) ::  QC3D               ! CLOUD WATER MIXING RATIO (KG/KG)
    REAL, DIMENSION(KTS:KTE) ::  QI3D               ! CLOUD ICE MIXING RATIO (KG/KG)
    REAL, DIMENSION(KTS:KTE) ::  QNI3D              ! SNOW MIXING RATIO (KG/KG)
    REAL, DIMENSION(KTS:KTE) ::  QR3D               ! RAIN MIXING RATIO (KG/KG)
    REAL, DIMENSION(KTS:KTE) ::  NI3D               ! CLOUD ICE NUMBER CONCENTRATION (1/KG)
    REAL, DIMENSION(KTS:KTE) ::  NS3D               ! SNOW NUMBER CONCENTRATION (1/KG)
    REAL, DIMENSION(KTS:KTE) ::  NR3D               ! RAIN NUMBER CONCENTRATION (1/KG)
    REAL, DIMENSION(KTS:KTE) ::  T3DTEN             ! TEMPERATURE TENDENCY (K/S)
    REAL, DIMENSION(KTS:KTE) ::  QV3DTEN            ! WATER VAPOR MIXING RATIO TENDENCY (KG/KG/S)
    REAL, DIMENSION(KTS:KTE) ::  T3D                ! TEMPERATURE (K)
    REAL, DIMENSION(KTS:KTE) ::  QV3D               ! WATER VAPOR MIXING RATIO (KG/KG)
    REAL, DIMENSION(KTS:KTE) ::  PRES               ! ATMOSPHERIC PRESSURE (PA)
    REAL, DIMENSION(KTS:KTE) ::  DZQ                ! DIFFERENCE IN HEIGHT ACROSS LEVEL (m)
    REAL, DIMENSION(KTS:KTE) ::  W3D                ! GRID-SCALE VERTICAL VELOCITY (M/S)
    REAL, DIMENSION(KTS:KTE) ::  WVAR               ! SUB-GRID VERTICAL VELOCITY (M/S)
    ! below for wrf-chem
    REAL, DIMENSION(KTS:KTE) ::  nc3d
    REAL, DIMENSION(KTS:KTE) ::  nc3dten
    integer, intent(in) :: iinum

    ! HM ADDED GRAUPEL VARIABLES
    REAL, DIMENSION(KTS:KTE) ::  QG3DTEN            ! GRAUPEL MIX RATIO TENDENCY (KG/KG/S)
    REAL, DIMENSION(KTS:KTE) ::  NG3DTEN            ! GRAUPEL NUMB CONC TENDENCY (1/KG/S)
    REAL, DIMENSION(KTS:KTE) ::  QG3D            ! GRAUPEL MIX RATIO (KG/KG)
    REAL, DIMENSION(KTS:KTE) ::  NG3D            ! GRAUPEL NUMBER CONC (1/KG)

    ! HM, ADD 1/16/07, SEDIMENTATION TENDENCIES FOR MIXING RATIO

    REAL, DIMENSION(KTS:KTE) ::  QGSTEN            ! GRAUPEL SED TEND (KG/KG/S)
    REAL, DIMENSION(KTS:KTE) ::  QRSTEN            ! RAIN SED TEND (KG/KG/S)
    REAL, DIMENSION(KTS:KTE) ::  QISTEN            ! CLOUD ICE SED TEND (KG/KG/S)
    REAL, DIMENSION(KTS:KTE) ::  QNISTEN           ! SNOW SED TEND (KG/KG/S)
    REAL, DIMENSION(KTS:KTE) ::  QCSTEN            ! CLOUD WAT SED TEND (KG/KG/S)

    ! hm add cumulus tendencies for precip
    REAL, DIMENSION(KTS:KTE) ::   qrcu1d
    REAL, DIMENSION(KTS:KTE) ::   qscu1d
    REAL, DIMENSION(KTS:KTE) ::   qicu1d

    ! OUTPUT VARIABLES

    REAL PRECRT                ! TOTAL PRECIP PER TIME STEP (mm)
    REAL SNOWRT                ! SNOW PER TIME STEP (mm)
    ! hm added 7/13/13
    REAL SNOWPRT      ! TOTAL CLOUD ICE PLUS SNOW PER TIME STEP (mm)
    REAL GRPLPRT	  ! TOTAL GRAUPEL PER TIME STEP (mm)

    REAL, DIMENSION(KTS:KTE) ::   EFFC            ! DROPLET EFFECTIVE RADIUS (MICRON)
    REAL, DIMENSION(KTS:KTE) ::   EFFI            ! CLOUD ICE EFFECTIVE RADIUS (MICRON)
    REAL, DIMENSION(KTS:KTE) ::   EFFS            ! SNOW EFFECTIVE RADIUS (MICRON)
    REAL, DIMENSION(KTS:KTE) ::   EFFR            ! RAIN EFFECTIVE RADIUS (MICRON)
    REAL, DIMENSION(KTS:KTE) ::   EFFG            ! GRAUPEL EFFECTIVE RADIUS (MICRON)

    ! MODEL INPUT PARAMETERS (FORMERLY IN COMMON BLOCKS)

    REAL DT         ! MODEL TIME STEP (SEC)

    !.....................................................................................................
    ! LOCAL VARIABLES: ALL PARAMETERS BELOW ARE LOCAL TO SCHEME AND DON'T NEED TO COMMUNICATE WITH THE
    ! REST OF THE MODEL.

    ! SIZE PARAMETER VARIABLES

    REAL, DIMENSION(KTS:KTE) :: LAMC          ! SLOPE PARAMETER FOR DROPLETS (M-1)
    REAL, DIMENSION(KTS:KTE) :: LAMI          ! SLOPE PARAMETER FOR CLOUD ICE (M-1)
    REAL, DIMENSION(KTS:KTE) :: LAMS          ! SLOPE PARAMETER FOR SNOW (M-1)
    REAL, DIMENSION(KTS:KTE) :: LAMR          ! SLOPE PARAMETER FOR RAIN (M-1)
    REAL, DIMENSION(KTS:KTE) :: LAMG          ! SLOPE PARAMETER FOR GRAUPEL (M-1)
    REAL, DIMENSION(KTS:KTE) :: CDIST1        ! PSD PARAMETER FOR DROPLETS
    REAL, DIMENSION(KTS:KTE) :: N0I           ! INTERCEPT PARAMETER FOR CLOUD ICE (KG-1 M-1)
    REAL, DIMENSION(KTS:KTE) :: N0S           ! INTERCEPT PARAMETER FOR SNOW (KG-1 M-1)
    REAL, DIMENSION(KTS:KTE) :: N0RR          ! INTERCEPT PARAMETER FOR RAIN (KG-1 M-1)
    REAL, DIMENSION(KTS:KTE) :: N0G           ! INTERCEPT PARAMETER FOR GRAUPEL (KG-1 M-1)
    REAL, DIMENSION(KTS:KTE) :: PGAM          ! SPECTRAL SHAPE PARAMETER FOR DROPLETS

    ! MICROPHYSICAL PROCESSES

    REAL, DIMENSION(KTS:KTE) ::  NSUBC     ! LOSS OF NC DURING EVAP
    REAL, DIMENSION(KTS:KTE) ::  NSUBI     ! LOSS OF NI DURING SUB.
    REAL, DIMENSION(KTS:KTE) ::  NSUBS     ! LOSS OF NS DURING SUB.
    REAL, DIMENSION(KTS:KTE) ::  NSUBR     ! LOSS OF NR DURING EVAP
    REAL, DIMENSION(KTS:KTE) ::  PRD       ! DEP CLOUD ICE
    REAL, DIMENSION(KTS:KTE) ::  PRE       ! EVAP OF RAIN
    REAL, DIMENSION(KTS:KTE) ::  PRDS      ! DEP SNOW
    REAL, DIMENSION(KTS:KTE) ::  NNUCCC    ! CHANGE N DUE TO CONTACT FREEZ DROPLETS
    REAL, DIMENSION(KTS:KTE) ::  MNUCCC    ! CHANGE Q DUE TO CONTACT FREEZ DROPLETS
    REAL, DIMENSION(KTS:KTE) ::  PRA       ! ACCRETION DROPLETS BY RAIN
    REAL, DIMENSION(KTS:KTE) ::  PRC       ! AUTOCONVERSION DROPLETS
    REAL, DIMENSION(KTS:KTE) ::  PCC       ! COND/EVAP DROPLETS
    REAL, DIMENSION(KTS:KTE) ::  NNUCCD    ! CHANGE N FREEZING AEROSOL (PRIM ICE NUCLEATION)
    REAL, DIMENSION(KTS:KTE) ::  MNUCCD    ! CHANGE Q FREEZING AEROSOL (PRIM ICE NUCLEATION)
    REAL, DIMENSION(KTS:KTE) ::  MNUCCR    ! CHANGE Q DUE TO CONTACT FREEZ RAIN
    REAL, DIMENSION(KTS:KTE) ::  NNUCCR    ! CHANGE N DUE TO CONTACT FREEZ RAIN
    REAL, DIMENSION(KTS:KTE) ::  NPRA      ! CHANGE IN N DUE TO DROPLET ACC BY RAIN
    REAL, DIMENSION(KTS:KTE) ::  NRAGG     ! SELF-COLLECTION/BREAKUP OF RAIN
    REAL, DIMENSION(KTS:KTE) ::  NSAGG     ! SELF-COLLECTION OF SNOW
    REAL, DIMENSION(KTS:KTE) ::  NPRC      ! CHANGE NC AUTOCONVERSION DROPLETS
    REAL, DIMENSION(KTS:KTE) ::  NPRC1      ! CHANGE NR AUTOCONVERSION DROPLETS
    REAL, DIMENSION(KTS:KTE) ::  PRAI      ! CHANGE Q ACCRETION CLOUD ICE BY SNOW
    REAL, DIMENSION(KTS:KTE) ::  PRCI      ! CHANGE Q AUTOCONVERSIN CLOUD ICE TO SNOW
    REAL, DIMENSION(KTS:KTE) ::  PSACWS    ! CHANGE Q DROPLET ACCRETION BY SNOW
    REAL, DIMENSION(KTS:KTE) ::  NPSACWS   ! CHANGE N DROPLET ACCRETION BY SNOW
    REAL, DIMENSION(KTS:KTE) ::  PSACWI    ! CHANGE Q DROPLET ACCRETION BY CLOUD ICE
    REAL, DIMENSION(KTS:KTE) ::  NPSACWI   ! CHANGE N DROPLET ACCRETION BY CLOUD ICE
    REAL, DIMENSION(KTS:KTE) ::  NPRCI     ! CHANGE N AUTOCONVERSION CLOUD ICE BY SNOW
    REAL, DIMENSION(KTS:KTE) ::  NPRAI     ! CHANGE N ACCRETION CLOUD ICE
    REAL, DIMENSION(KTS:KTE) ::  NMULTS    ! ICE MULT DUE TO RIMING DROPLETS BY SNOW
    REAL, DIMENSION(KTS:KTE) ::  NMULTR    ! ICE MULT DUE TO RIMING RAIN BY SNOW
    REAL, DIMENSION(KTS:KTE) ::  QMULTS    ! CHANGE Q DUE TO ICE MULT DROPLETS/SNOW
    REAL, DIMENSION(KTS:KTE) ::  QMULTR    ! CHANGE Q DUE TO ICE RAIN/SNOW
    REAL, DIMENSION(KTS:KTE) ::  PRACS     ! CHANGE Q RAIN-SNOW COLLECTION
    REAL, DIMENSION(KTS:KTE) ::  NPRACS    ! CHANGE N RAIN-SNOW COLLECTION
    REAL, DIMENSION(KTS:KTE) ::  PCCN      ! CHANGE Q DROPLET ACTIVATION
    REAL, DIMENSION(KTS:KTE) ::  PSMLT     ! CHANGE Q MELTING SNOW TO RAIN
    REAL, DIMENSION(KTS:KTE) ::  EVPMS     ! CHNAGE Q MELTING SNOW EVAPORATING
    REAL, DIMENSION(KTS:KTE) ::  NSMLTS    ! CHANGE N MELTING SNOW
    REAL, DIMENSION(KTS:KTE) ::  NSMLTR    ! CHANGE N MELTING SNOW TO RAIN
    ! HM ADDED 12/13/06
    REAL, DIMENSION(KTS:KTE) ::  PIACR     ! CHANGE QR, ICE-RAIN COLLECTION
    REAL, DIMENSION(KTS:KTE) ::  NIACR     ! CHANGE N, ICE-RAIN COLLECTION
    REAL, DIMENSION(KTS:KTE) ::  PRACI     ! CHANGE QI, ICE-RAIN COLLECTION
    REAL, DIMENSION(KTS:KTE) ::  PIACRS     ! CHANGE QR, ICE RAIN COLLISION, ADDED TO SNOW
    REAL, DIMENSION(KTS:KTE) ::  NIACRS     ! CHANGE N, ICE RAIN COLLISION, ADDED TO SNOW
    REAL, DIMENSION(KTS:KTE) ::  PRACIS     ! CHANGE QI, ICE RAIN COLLISION, ADDED TO SNOW
    REAL, DIMENSION(KTS:KTE) ::  EPRD      ! SUBLIMATION CLOUD ICE
    REAL, DIMENSION(KTS:KTE) ::  EPRDS     ! SUBLIMATION SNOW
    ! HM ADDED GRAUPEL PROCESSES
    REAL, DIMENSION(KTS:KTE) ::  PRACG    ! CHANGE IN Q COLLECTION RAIN BY GRAUPEL
    REAL, DIMENSION(KTS:KTE) ::  PSACWG    ! CHANGE IN Q COLLECTION DROPLETS BY GRAUPEL
    REAL, DIMENSION(KTS:KTE) ::  PGSACW    ! CONVERSION Q TO GRAUPEL DUE TO COLLECTION DROPLETS BY SNOW
    REAL, DIMENSION(KTS:KTE) ::  PGRACS    ! CONVERSION Q TO GRAUPEL DUE TO COLLECTION RAIN BY SNOW
    REAL, DIMENSION(KTS:KTE) ::  PRDG    ! DEP OF GRAUPEL
    REAL, DIMENSION(KTS:KTE) ::  EPRDG    ! SUB OF GRAUPEL
    REAL, DIMENSION(KTS:KTE) ::  EVPMG    ! CHANGE Q MELTING OF GRAUPEL AND EVAPORATION
    REAL, DIMENSION(KTS:KTE) ::  PGMLT    ! CHANGE Q MELTING OF GRAUPEL
    REAL, DIMENSION(KTS:KTE) ::  NPRACG    ! CHANGE N COLLECTION RAIN BY GRAUPEL
    REAL, DIMENSION(KTS:KTE) ::  NPSACWG    ! CHANGE N COLLECTION DROPLETS BY GRAUPEL
    REAL, DIMENSION(KTS:KTE) ::  NSCNG    ! CHANGE N CONVERSION TO GRAUPEL DUE TO COLLECTION DROPLETS BY SNOW
    REAL, DIMENSION(KTS:KTE) ::  NGRACS    ! CHANGE N CONVERSION TO GRAUPEL DUE TO COLLECTION RAIN BY SNOW
    REAL, DIMENSION(KTS:KTE) ::  NGMLTG    ! CHANGE N MELTING GRAUPEL
    REAL, DIMENSION(KTS:KTE) ::  NGMLTR    ! CHANGE N MELTING GRAUPEL TO RAIN
    REAL, DIMENSION(KTS:KTE) ::  NSUBG    ! CHANGE N SUB/DEP OF GRAUPEL
    REAL, DIMENSION(KTS:KTE) ::  PSACR    ! CONVERSION DUE TO COLL OF SNOW BY RAIN
    REAL, DIMENSION(KTS:KTE) ::  NMULTG    ! ICE MULT DUE TO ACC DROPLETS BY GRAUPEL
    REAL, DIMENSION(KTS:KTE) ::  NMULTRG    ! ICE MULT DUE TO ACC RAIN BY GRAUPEL
    REAL, DIMENSION(KTS:KTE) ::  QMULTG    ! CHANGE Q DUE TO ICE MULT DROPLETS/GRAUPEL
    REAL, DIMENSION(KTS:KTE) ::  QMULTRG    ! CHANGE Q DUE TO ICE MULT RAIN/GRAUPEL

    ! TIME-VARYING ATMOSPHERIC PARAMETERS

    REAL, DIMENSION(KTS:KTE) ::   KAP   ! THERMAL CONDUCTIVITY OF AIR
    REAL, DIMENSION(KTS:KTE) ::   EVS   ! SATURATION VAPOR PRESSURE
    REAL, DIMENSION(KTS:KTE) ::   EIS   ! ICE SATURATION VAPOR PRESSURE
    REAL, DIMENSION(KTS:KTE) ::   QVS   ! SATURATION MIXING RATIO
    REAL, DIMENSION(KTS:KTE) ::   QVI   ! ICE SATURATION MIXING RATIO
    REAL, DIMENSION(KTS:KTE) ::   QVQVS ! SAUTRATION RATIO
    REAL, DIMENSION(KTS:KTE) ::   QVQVSI! ICE SATURAION RATIO
    REAL, DIMENSION(KTS:KTE) ::   DV    ! DIFFUSIVITY OF WATER VAPOR IN AIR
    REAL, DIMENSION(KTS:KTE) ::   XXLS  ! LATENT HEAT OF SUBLIMATION
    REAL, DIMENSION(KTS:KTE) ::   XXLV  ! LATENT HEAT OF VAPORIZATION
    REAL, DIMENSION(KTS:KTE) ::   CPM   ! SPECIFIC HEAT AT CONST PRESSURE FOR MOIST AIR
    REAL, DIMENSION(KTS:KTE) ::   MU    ! VISCOCITY OF AIR
    REAL, DIMENSION(KTS:KTE) ::   SC    ! SCHMIDT NUMBER
    REAL, DIMENSION(KTS:KTE) ::   XLF   ! LATENT HEAT OF FREEZING
    REAL, DIMENSION(KTS:KTE) ::   RHO   ! AIR DENSITY
    REAL, DIMENSION(KTS:KTE) ::   AB    ! CORRECTION TO CONDENSATION RATE DUE TO LATENT HEATING
    REAL, DIMENSION(KTS:KTE) ::   ABI    ! CORRECTION TO DEPOSITION RATE DUE TO LATENT HEATING

    ! TIME-VARYING MICROPHYSICS PARAMETERS

    REAL, DIMENSION(KTS:KTE) ::   DAP    ! DIFFUSIVITY OF AEROSOL
    REAL    NACNT                    ! NUMBER OF CONTACT IN
    REAL    FMULT                    ! TEMP.-DEP. PARAMETER FOR RIME-SPLINTERING
    REAL    COFFI                    ! ICE AUTOCONVERSION PARAMETER

    ! FALL SPEED WORKING VARIABLES (DEFINED IN CODE)

    REAL, DIMENSION(KTS:KTE) ::    DUMI,DUMR,DUMFNI,DUMG,DUMFNG
    REAL UNI, UMI,UMR
    REAL, DIMENSION(KTS:KTE) ::    FR, FI, FNI,FG,FNG
    REAL RGVM
    REAL, DIMENSION(KTS:KTE) ::   FALOUTR,FALOUTI,FALOUTNI
    REAL FALTNDR,FALTNDI,FALTNDNI,RHO2
    REAL, DIMENSION(KTS:KTE) ::   DUMQS,DUMFNS
    REAL UMS,UNS
    REAL, DIMENSION(KTS:KTE) ::   FS,FNS, FALOUTS,FALOUTNS,FALOUTG,FALOUTNG
    REAL FALTNDS,FALTNDNS,UNR,FALTNDG,FALTNDNG
    REAL, DIMENSION(KTS:KTE) ::    DUMC,DUMFNC
    REAL UNC,UMC,UNG,UMG
    REAL, DIMENSION(KTS:KTE) ::   FC,FALOUTC,FALOUTNC
    REAL FALTNDC,FALTNDNC
    REAL, DIMENSION(KTS:KTE) ::   FNC,DUMFNR,FALOUTNR
    REAL FALTNDNR
    REAL, DIMENSION(KTS:KTE) ::   FNR

    ! FALL-SPEED PARAMETER 'A' WITH AIR DENSITY CORRECTION

    REAL, DIMENSION(KTS:KTE) ::    AIN,ARN,ASN,ACN,AGN

    ! EXTERNAL FUNCTION CALL RETURN VARIABLES

    !      REAL GAMMA,      ! EULER GAMMA FUNCTION
    !      REAL POLYSVP,    ! SAT. PRESSURE FUNCTION
    !      REAL DERF1        ! ERROR FUNCTION

    ! DUMMY VARIABLES

    REAL DUM,DUM1,DUM2,DUMT,DUMQV,DUMQSS,DUMQSI,DUMS

    ! PROGNOSTIC SUPERSATURATION

    REAL DQSDT    ! CHANGE OF SAT. MIX. RAT. WITH TEMPERATURE
    REAL DQSIDT   ! CHANGE IN ICE SAT. MIXING RAT. WITH T
    REAL EPSI     ! 1/PHASE REL. TIME (SEE M2005), ICE
    REAL EPSS     ! 1/PHASE REL. TIME (SEE M2005), SNOW
    REAL EPSR     ! 1/PHASE REL. TIME (SEE M2005), RAIN
    REAL EPSG     ! 1/PHASE REL. TIME (SEE M2005), GRAUPEL

    ! NEW DROPLET ACTIVATION VARIABLES
    REAL TAUC     ! PHASE REL. TIME (SEE M2005), DROPLETS
    REAL TAUR     ! PHASE REL. TIME (SEE M2005), RAIN
    REAL TAUI     ! PHASE REL. TIME (SEE M2005), CLOUD ICE
    REAL TAUS     ! PHASE REL. TIME (SEE M2005), SNOW
    REAL TAUG     ! PHASE REL. TIME (SEE M2005), GRAUPEL
    REAL DUMACT,DUM3

    ! COUNTING/INDEX VARIABLES

    INTEGER K,NSTEP,N ! ,I

    ! LTRUE IS ONLY USED TO SPEED UP THE CODE !!
    ! LTRUE, SWITCH = 0, NO HYDROMETEORS IN COLUMN,
    !               = 1, HYDROMETEORS IN COLUMN

    INTEGER LTRUE

    ! DROPLET ACTIVATION/FREEZING AEROSOL


    REAL    CT      ! DROPLET ACTIVATION PARAMETER
    REAL    TEMP1   ! DUMMY TEMPERATURE
    REAL    SAT1    ! DUMMY SATURATION
    REAL    SIGVL   ! SURFACE TENSION LIQ/VAPOR
    REAL    KEL     ! KELVIN PARAMETER
    REAL    KC2     ! TOTAL ICE NUCLEATION RATE

    REAL CRY,KRY   ! AEROSOL ACTIVATION PARAMETERS

    ! MORE WORKING/DUMMY VARIABLES

    REAL DUMQI,DUMNI,DC0,DS0,DG0
    REAL DUMQC,DUMQR,RATIO,SUM_DEP,FUDGEF

    ! EFFECTIVE VERTICAL VELOCITY  (M/S)
    REAL WEF

    ! WORKING PARAMETERS FOR ICE NUCLEATION

    REAL ANUC,BNUC

    ! WORKING PARAMETERS FOR AEROSOL ACTIVATION

    REAL AACT,GAMM,GG,PSI,ETA1,ETA2,SM1,SM2,SMAX,UU1,UU2,ALPHA

    ! DUMMY SIZE DISTRIBUTION PARAMETERS

    REAL DLAMS,DLAMR,DLAMI,DLAMC,DLAMG,LAMMAX,LAMMIN

    INTEGER IDROP

    ! FOR WRF-CHEM
    REAL, DIMENSION(KTS:KTE)::C2PREC,CSED,ISED,SSED,GSED,RSED
#if (WRF_CHEM == 1)
    REAL, DIMENSION(KTS:KTE), INTENT(INOUT) :: rainprod, evapprod
#endif
    REAL, DIMENSION(KTS:KTE)                :: tqimelt ! melting of cloud ice (tendency)

    ! comment lines for wrf-chem since these are intent(in) in that case
    !       REAL, DIMENSION(KTS:KTE) ::  NC3DTEN            ! CLOUD DROPLET NUMBER CONCENTRATION (1/KG/S)
    !       REAL, DIMENSION(KTS:KTE) ::  NC3D               ! CLOUD DROPLET NUMBER CONCENTRATION (1/KG)

    !-----------------------------------------------------------------------
    ! yc added 9/22/2020 for lightning
    INTEGER, INTENT(IN)  :: ipelec, i, j, itimestep

    ! consts
    REAL DIAG,DIAS,DIAI ! # WEIGHTED MEAN DIAMETER
    REAL VT ! RELATIVE FALL SPEED
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
    real cfce15,gf8,gf4i,gf3p5,gf1a,gf1p5
    logical :: csi_morr = .True.
    real ni
    real, parameter :: tfr = 273.15, thnuc = 235.15
    integer, intent(in):: IDS,IDE, JDS,JDE, KDS,KDE ! domain boundaries

    ! tmps..
    REAL TMP, tmppre
    REAL :: tmp1 = 0., dumtmp
    REAL :: tmpscr(kts:kte), dumrtmp(kts:kte)
    REAL CCWTMP,CCITMP ! ,CIPTMP,CIRTMP
    DOUBLE PRECISION FRAC
    INTEGER IL, IC
    REAL SCFAC
    REAL CPQC,CPCI
    REAL CPQC0,CPCI0
    REAL QCWTMP(KTS:KTE)
    REAL QCITMP
    REAL IL5(KTS:KTE)
    real delqnxa(lc:lqmx)
    real delqxxa(lc:lqmx)
    integer idelq
    real :: scworig=0.0
    integer iltest

    REAL, INTENT(INOUT), DIMENSION(KTS:KTE) ::   CSACI0 ! SNOW-ICE COLLISION RATE (#COLLI./M3/S)
    REAL, INTENT(INOUT), DIMENSION(KTS:KTE) ::   CHACI0 ! GRAUP-ICE COLLISION RATE (#COLLI./M3/S)
    REAL, INTENT(INOUT), DIMENSION(KTS:KTE) ::   CHACS0 ! GRAUP-SNOW COLLISION RATE (#COLLI./M3/S)
    ! csaci0_morr is simply 9 times NPRAI for now, since eii=0.1
    REAL, INTENT(INOUT), DIMENSION(KTS:KTE) ::   csaci0_morr ! SNOW-ICE COLLISION RATE, morr (#COLLI./kg/S)
    REAL, INTENT(INOUT), DIMENSION(KTS:KTE) ::  umc1d, umr1d, umi1d, ums1d, umg1d ! 1d fallspeed to pass out

    REAL, DIMENSION(KTS:KTE, NA) :: AN1D ! CONTAINS Q, QN, SC per volumn m3
    REAL, DIMENSION(KTS:KTE, lscw:lsch) :: DUMSC ! for use in sedimentation
    REAL :: SCX(KTS:KTE,LC:LHAB) ! space charge per m3
    REAL :: CIONP(KTS:KTE),CIONN(KTS:KTE),CLIONP(KTS:KTE),CLIONN(KTS:KTE)
    REAL CIONPMXD(KTS:KTE),CIONNMXD(KTS:KTE)
    REAL CLIONPMXD(KTS:KTE),CLIONNMXD(KTS:KTE)
    REAL temcg(KTS:KTE) ! T in Celcius
    REAL elec1d(KTS:KTE)
    REAL dezcomp(KTS:KTE)
    REAL rhoinv(KTS:KTE)
    integer iexy(lc:lqmx,lc:lqmx)
    real exy(KTS:KTE,ls:lhab,lc:ls) ! sticking efficiency
    real cxacy(KTS:KTE,ls:lhab,lc:ls) ! collision rates
    real xdia(KTS:KTE,lc:lhab) ! diameter; =1/lamx
    real umx(KTS:KTE,lc:lhab) ! fall speed
    real rarx(KTS:KTE,ls:lhab)
    real, dimension(kts:kte) :: t2,t3,t5,t6,t8

    ! scxacy series
    real :: scsacw(kts:kte)
    real :: scsacr(kts:kte)
    real :: scsaci(kts:kte)
    !
    real :: schacw(kts:kte), schacr(kts:kte)
    real :: schaci(kts:kte), schacs(kts:kte)
    !
    real :: schmlr(kts:kte), scsmlr(kts:kte)
    real :: scsshr(kts:kte), schshr(kts:kte)
    real :: scsdep(kts:kte)
    !
    real :: schlacw(kts:kte), schlacr(kts:kte)
    real :: schlaci(kts:kte), schlacs(kts:kte)
    !
    real :: scxacy(kts:kte,ls:lhab,lc:ls)

    ! charge per kg hdmt
    real :: fscrw(kts:kte)
    real :: fsccw(kts:kte)
    real :: fscci(kts:kte)
    real :: fscsw(kts:kte)
    real :: fschw(kts:kte)

    ! production of space charge term  + / -
    real ::  psctot(KTS:KTE)
    real ::  psccwi(KTS:KTE), psccwd(KTS:KTE) ! cloud water
    real ::  psccii(KTS:KTE), psccid(KTS:KTE) ! ice crystals
    real ::  pscisi(KTS:KTE), pscisd(KTS:KTE) ! ice spheres
    real ::  pscrwi(KTS:KTE), pscrwd(KTS:KTE) ! rain water
    real ::  pscswi(KTS:KTE), pscswd(KTS:KTE) ! snow
    real ::  pschwi(KTS:KTE), pschwd(KTS:KTE) ! graupel
    real ::  pschli(KTS:KTE), pschld(KTS:KTE) ! hail

    ! production of space charge from mass transfer
    real ::  psccwmi(KTS:KTE), psccwmd(KTS:KTE)
    real ::  psccimi(KTS:KTE), psccimd(KTS:KTE)
    real ::  pscismi(KTS:KTE), pscismd(KTS:KTE)
    real ::  pscrwmi(KTS:KTE), pscrwmd(KTS:KTE)
    real ::  pscswmi(KTS:KTE), pscswmd(KTS:KTE)
    real ::  pschwmi(KTS:KTE), pschwmd(KTS:KTE)
    real ::  pschlmi(KTS:KTE), pschlmd(KTS:KTE)

    ! ions: assume mass = 0
    real ::  pscpii(KTS:KTE), pscpid(KTS:KTE) !positive ions
    real ::  pscnii(KTS:KTE), pscnid(KTS:KTE) !negative ions
    real ::  pscplii(KTS:KTE), pscplid(KTS:KTE) !positive (large) ions
    real ::  pscnlii(KTS:KTE), pscnlid(KTS:KTE) !negative (large) ions

    ! sedimentation rates
    real ::  pscc_sed(kts:kte)
    real ::  pscr_sed(kts:kte)
    real ::  psci_sed(kts:kte)
    real ::  pscs_sed(kts:kte)
    real ::  pscg_sed(kts:kte)
    real, dimension(kts:kte) :: faloutscc,faloutscr,faloutsci,faloutscs,faloutscg
    real :: faltndscc,faltndscr,faltndsci,faltndscs,faltndscg

    !-----------------------------------------------------------------------
    ! end lightning

    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

    ! SET LTRUE INITIALLY TO 0

    LTRUE = 0

    IF (IPELEC>0) THEN
        ! VALUES NEEDED FOR EXPLICIT LIGHTNING

        DTP = DT
        DTPINV = 1.D0/DTP
        RHOINV(KTS:KTE) = (R*T3D(K))/PRES(K)

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

        DO K = KTS,KTE
        dezcomp(k) = Max(-1.00e5, Min(elec1d(k),1.00e5) ) ! Ez is passed (elecsave)
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
        umc1d  = 0.0
        umr1d  = 0.0
        umi1d  = 0.0
        ums1d  = 0.0
        umg1d  = 0.0

        schacw(kts:kte) = 0.0
        scsacw(kts:kte) = 0.0
        schaci(kts:kte) = 0.0
        scsaci(kts:kte) = 0.0
        scsacr(kts:kte) = 0.0
        schacr(kts:kte) = 0.0
        schacs(kts:kte) = 0.0
        schmlr(kts:kte) = 0.0
        scsmlr(kts:kte) = 0.0
        scsdep(kts:kte) = 0.0
        scsshr(kts:kte) = 0.0
        schshr(kts:kte) = 0.0
        schlacw(kts:kte) = 0.0
        schlacr(kts:kte) = 0.0
        schlaci(kts:kte) = 0.0
        schlacs(kts:kte) = 0.0

        fscrw(kts:kte) = 0.0
        fsccw(kts:kte) = 0.0
        fscci(kts:kte) = 0.0
        fscsw(kts:kte) = 0.0
        fschw(kts:kte) = 0.0

        psccwi(kts:kte) = 0.0 ! cloud water
        psccwd(kts:kte) = 0.0
        psccwmi(kts:kte) = 0.0
        psccwmd(kts:kte) = 0.0
        psccii(kts:kte) = 0.0 ! cloud ice
        psccid(kts:kte) = 0.0
        psccimi(kts:kte) = 0.0
        psccimd(kts:kte) = 0.0
        pscrwi(kts:kte) = 0.0 ! rain water
        pscrwd(kts:kte) = 0.0
        pscrwmi(kts:kte) = 0.0
        pscrwmd(kts:kte) = 0.0
        pscswi(kts:kte) = 0.0 ! snow water
        pscswd(kts:kte) = 0.0
        pscswmi(kts:kte) = 0.0
        pscswmd(kts:kte) = 0.0
        pschwi(kts:kte) = 0.0 ! graupel
        pschwd(kts:kte) = 0.0
        pschwmi(kts:kte) = 0.0
        pschwmd(kts:kte) = 0.0

        pscpii(KTS:KTE) = 0.0
        pscpid(KTS:KTE) = 0.0
        pscnii(KTS:KTE) = 0.0
        pscnid(KTS:KTE) = 0.0
        pscplii(KTS:KTE) = 0.0
        pscplid(KTS:KTE) = 0.0
        pscnlii(KTS:KTE) = 0.0
        pscnlid(KTS:KTE) = 0.0

        pscc_sed(kts:kte) = 0.0
        pscr_sed(kts:kte) = 0.0
        psci_sed(kts:kte) = 0.0
        pscs_sed(kts:kte) = 0.0
        pscg_sed(kts:kte) = 0.0

        CSACI0(kts:kte) = 0.0
        CHACI0(kts:kte) = 0.0
        CHACS0(kts:kte) = 0.0
        csaci0_morr(kts:kte) = 0.0

        t2(kts:kte) = 0.0
        t3(kts:kte) = 0.0
        t5(kts:kte) = 0.0
        t6(kts:kte) = 0.0
        t8(kts:kte) = 0.0

        faloutscc(kts:kte) = 0.0
        faloutscr(kts:kte) = 0.0
        faloutsci(kts:kte) = 0.0
        faloutscs(kts:kte) = 0.0
        faloutscg(kts:kte) = 0.0

        faltndscc = 0.0
        faltndscr = 0.0
        faltndsci = 0.0
        faltndscs = 0.0
        faltndscg = 0.0

        PRD(kts:kte) = 0.0
        PRE(kts:kte) = 0.0
        PRDS(kts:kte) = 0.0
        MNUCCC(kts:kte) = 0.0
        PRA(kts:kte) = 0.0
        PRC(kts:kte) = 0.0
        PCC(kts:kte) = 0.0
        MNUCCD(kts:kte) = 0.0
        MNUCCR(kts:kte) = 0.0
        PRAI(kts:kte) = 0.0
        PRCI(kts:kte) = 0.0
        PSACWS(kts:kte) = 0.0
        PSACWI(kts:kte) = 0.0
        QMULTS(kts:kte) = 0.0
        QMULTR(kts:kte) = 0.0
        PRACS(kts:kte) = 0.0
        PCCN(kts:kte) = 0.0
        PSMLT(kts:kte) = 0.0
        EVPMS(kts:kte) = 0.0
        PIACR(kts:kte) = 0.0
        PRACI(kts:kte) = 0.0
        PIACRS(kts:kte) = 0.0
        PRACIS(kts:kte) = 0.0
        EPRD(kts:kte) = 0.0
        EPRDS(kts:kte) = 0.0
        PRACG(kts:kte) = 0.0
        PSACWG(kts:kte) = 0.0
        PGSACW(kts:kte) = 0.0
        PGRACS(kts:kte) = 0.0
        PRDG(kts:kte) = 0.0
        EPRDG(kts:kte) = 0.0
        EVPMG(kts:kte) = 0.0
        PGMLT(kts:kte) = 0.0
        PSACR(kts:kte) = 0.0
        QMULTG(kts:kte) = 0.0
        QMULTRG(kts:kte) = 0.0

        DUMSC(:,:) = 0.0

    END IF !(IPELEC>0)

    ! ATMOSPHERIC PARAMETERS THAT VARY IN TIME AND HEIGHT
    DO K = KTS,KTE

        ! NC3DTEN LOCAL ARRAY INITIALIZED
        NC3DTEN(K) = 0.
        ! INITIALIZE VARIABLES FOR WRF-CHEM OUTPUT TO ZERO

        C2PREC(K)=0.
        CSED(K)=0.
        ISED(K)=0.
        SSED(K)=0.
        GSED(K)=0.
        RSED(K)=0.

#if (WRF_CHEM == 1)
        rainprod(K) = 0.
        evapprod(K) = 0.
        tqimelt(K)  = 0.
        PRC(K)      = 0.
        PRA(K)      = 0.
#endif

        ! LATENT HEAT OF VAPORATION

        XXLV(K) = 3.1484E6-2370.*T3D(K)

        ! LATENT HEAT OF SUBLIMATION

        XXLS(K) = 3.15E6-2370.*T3D(K)+0.3337E6

        CPM(K) = CP*(1.+0.887*QV3D(K))

        ! SATURATION VAPOR PRESSURE AND MIXING RATIO

        ! hm, add fix for low pressure, 5/12/10
        EVS(K) = min(0.99*pres(k),POLYSVP(T3D(K),0))   ! PA
        EIS(K) = min(0.99*pres(k),POLYSVP(T3D(K),1))   ! PA

        ! MAKE SURE ICE SATURATION DOESN'T EXCEED WATER SAT. NEAR FREEZING

        IF (EIS(K).GT.EVS(K)) EIS(K) = EVS(K)

        QVS(K) = EP_2*EVS(K)/(PRES(K)-EVS(K))
        QVI(K) = EP_2*EIS(K)/(PRES(K)-EIS(K))

        QVQVS(K) = QV3D(K)/QVS(K)
        QVQVSI(K) = QV3D(K)/QVI(K)

        ! AIR DENSITY

        RHO(K) = PRES(K)/(R*T3D(K))

        ! ADD NUMBER CONCENTRATION DUE TO CUMULUS TENDENCY
        ! ASSUME N0 ASSOCIATED WITH CUMULUS PARAM RAIN IS 10^7 M^-4
        ! ASSUME N0 ASSOCIATED WITH CUMULUS PARAM SNOW IS 2 X 10^7 M^-4
        ! FOR DETRAINED CLOUD ICE, ASSUME MEAN VOLUME DIAM OF 80 MICRON

        IF (QRCU1D(K).GE.1.E-10) THEN
            DUM=1.8e5*(QRCU1D(K)*DT/(PI*RHOW*RHO(K)**3))**0.25
            NR3D(K)=NR3D(K)+DUM
        END IF
        IF (QSCU1D(K).GE.1.E-10) THEN
            DUM=3.e5*(QSCU1D(K)*DT/(CONS1*RHO(K)**3))**(1./(DS+1.))
            NS3D(K)=NS3D(K)+DUM
        END IF
        IF (QICU1D(K).GE.1.E-10) THEN
            DUM=QICU1D(K)*DT/(CI*(80.E-6)**DI)
            NI3D(K)=NI3D(K)+DUM
        END IF

        ! AT SUBSATURATION, REMOVE SMALL AMOUNTS OF CLOUD/PRECIP WATER
        ! hm modify 7/0/09 change limit to 1.e-8

        IF (QVQVS(K).LT.0.9) THEN
            IF (QR3D(K).LT.1.E-8) THEN
                QV3D(K)=QV3D(K)+QR3D(K)
                T3D(K)=T3D(K)-QR3D(K)*XXLV(K)/CPM(K)
                QR3D(K)=0.
                IF ( an1d(k,lscr) .gt. 0.0 ) THEN
                    an1d(k,lscpi) = an1d(k,lscpi) + an1d(k,lscr)*ec_i
                ELSE
                    an1d(k,lscni) = an1d(k,lscni) - an1d(k,lscr)*ec_i
                ENDIF
                an1d(k,lscr) = 0.0
            END IF
            IF (QC3D(K).LT.1.E-8) THEN
                QV3D(K)=QV3D(K)+QC3D(K)
                T3D(K)=T3D(K)-QC3D(K)*XXLV(K)/CPM(K)
                QC3D(K)=0.
                IF ( an1d(k,lscw) .gt. 0.0 ) THEN
                    an1d(k,lscpi) = an1d(k,lscpi) + an1d(k,lscw)*ec_i
                ELSE
                    an1d(k,lscni) = an1d(k,lscni) - an1d(k,lscw)*ec_i
                ENDIF
                an1d(k,lscw) = 0.0
            END IF
        END IF

        IF (QVQVSI(K).LT.0.9) THEN
            IF (QI3D(K).LT.1.E-8) THEN
                QV3D(K)=QV3D(K)+QI3D(K)
                T3D(K)=T3D(K)-QI3D(K)*XXLS(K)/CPM(K)
                QI3D(K)=0.
                IF ( an1d(k,lsci) .gt. 0.0 ) THEN
                    an1d(k,lscpi) = an1d(k,lscpi) + an1d(k,lsci)*ec_i
                ELSE
                    an1d(k,lscni) = an1d(k,lscni) - an1d(k,lsci)*ec_i
                ENDIF
                an1d(k,lsci) = 0.0
            END IF
            IF (QNI3D(K).LT.1.E-8) THEN
                QV3D(K)=QV3D(K)+QNI3D(K)
                T3D(K)=T3D(K)-QNI3D(K)*XXLS(K)/CPM(K)
                QNI3D(K)=0.
                IF ( an1d(k,lscs) .gt. 0.0 ) THEN
                    an1d(k,lscpi) = an1d(k,lscpi) + an1d(k,lscs)*ec_i
                ELSE
                    an1d(k,lscni) = an1d(k,lscni) - an1d(k,lscs)*ec_i
                ENDIF
                an1d(k,lscs) = 0.0
            END IF
            IF (QG3D(K).LT.1.E-8) THEN
                QV3D(K)=QV3D(K)+QG3D(K)
                T3D(K)=T3D(K)-QG3D(K)*XXLS(K)/CPM(K)
                QG3D(K)=0.
                IF ( an1d(k,lsch) .gt. 0.0 ) THEN
                    an1d(k,lscpi) = an1d(k,lscpi) + an1d(k,lsch)*ec_i
                ELSE
                    an1d(k,lscni) = an1d(k,lscni) - an1d(k,lsch)*ec_i
                ENDIF
                an1d(k,lsch) = 0.0
            END IF
        END IF

        ! HEAT OF FUSION

        XLF(K) = XXLS(K)-XXLV(K)

        !..................................................................
        ! IF MIXING RATIO < QSMALL SET MIXING RATIO AND NUMBER CONC TO ZERO

        IF (QC3D(K).LT.QSMALL) THEN
            QC3D(K) = 0.
            NC3D(K) = 0.
            EFFC(K) = 0.
            an1d(k,lscw) = 0.0
        END IF
        IF (QR3D(K).LT.QSMALL) THEN
            QR3D(K) = 0.
            NR3D(K) = 0.
            EFFR(K) = 0.
            an1d(k,lscr) = 0.0
        END IF
        IF (QI3D(K).LT.QSMALL) THEN
            QI3D(K) = 0.
            NI3D(K) = 0.
            EFFI(K) = 0.
            an1d(k,lsci) = 0.0
        END IF
        IF (QNI3D(K).LT.QSMALL) THEN
            QNI3D(K) = 0.
            NS3D(K) = 0.
            EFFS(K) = 0.
            an1d(k,lscs) = 0.0
        END IF
        IF (QG3D(K).LT.QSMALL) THEN
            QG3D(K) = 0.
            NG3D(K) = 0.
            EFFG(K) = 0.
            an1d(k,lsch) = 0.0
        END IF

        ! INITIALIZE SEDIMENTATION TENDENCIES FOR MIXING RATIO

        QRSTEN(K) = 0.
        QISTEN(K) = 0.
        QNISTEN(K) = 0.
        QCSTEN(K) = 0.
        QGSTEN(K) = 0.

        !..................................................................
        ! MICROPHYSICS PARAMETERS VARYING IN TIME/HEIGHT

        ! fix 053011
        MU(K) = 1.496E-6*T3D(K)**1.5/(T3D(K)+120.)

        ! FALL SPEED WITH DENSITY CORRECTION (HEYMSFIELD AND BENSSEMER 2006)

        DUM = (RHOSU/RHO(K))**0.54

        ! fix 053011
        !            AIN(K) = DUM*AI
        ! AA revision 4/1/11: Ikawa and Saito 1991 air-density correction
        AIN(K) = (RHOSU/RHO(K))**0.35*AI
        ARN(K) = DUM*AR
        ASN(K) = DUM*AS
        !            ACN(K) = DUM*AC
        ! AA revision 4/1/11: temperature-dependent Stokes fall speed
        ACN(K) = G*RHOW/(18.*MU(K))
        ! HM ADD GRAUPEL 8/28/06
        AGN(K) = DUM*AG

        ! hm 4/7/09 bug fix, initialize lami to prevent later division by zero
        LAMI(K)=0.

        !..................................
        ! IF THERE IS NO CLOUD/PRECIP WATER, AND IF SUBSATURATED, THEN SKIP MICROPHYSICS
        ! FOR THIS LEVEL

        IF (QC3D(K).LT.QSMALL.AND.QI3D(K).LT.QSMALL.AND.QNI3D(K).LT.QSMALL &
                .AND.QR3D(K).LT.QSMALL.AND.QG3D(K).LT.QSMALL) THEN
            IF (T3D(K).LT.273.15.AND.QVQVSI(K).LT.0.999) GOTO 200
            IF (T3D(K).GE.273.15.AND.QVQVS(K).LT.0.999) GOTO 200
        END IF

        ! THERMAL CONDUCTIVITY FOR AIR

        ! fix 053011
        KAP(K) = 1.414E3*MU(K)

        ! DIFFUSIVITY OF WATER VAPOR

        DV(K) = 8.794E-5*T3D(K)**1.81/PRES(K)

        ! SCHMIT NUMBER

        ! fix 053011
        SC(K) = MU(K)/(RHO(K)*DV(K))

        ! PSYCHOMETIC CORRECTIONS

        ! RATE OF CHANGE SAT. MIX. RATIO WITH TEMPERATURE

        DUM = (RV*T3D(K)**2)

        DQSDT = XXLV(K)*QVS(K)/DUM
        DQSIDT =  XXLS(K)*QVI(K)/DUM

        ABI(K) = 1.+DQSIDT*XXLS(K)/CPM(K)
        AB(K) = 1.+DQSDT*XXLV(K)/CPM(K)

        !
        !.....................................................................
        !.....................................................................
        ! CASE FOR TEMPERATURE ABOVE FREEZING

        IF (T3D(K).GE.273.15) THEN

            !......................................................................
            !HM ADD, ALLOW FOR CONSTANT DROPLET NUMBER
            ! INUM = 0, PREDICT DROPLET NUMBER
            ! INUM = 1, SET CONSTANT DROPLET NUMBER

            IF (iinum.EQ.1) THEN
                ! CONVERT NDCNST FROM CM-3 TO KG-1
                NC3D(K)=NDCNST*1.E6/RHO(K)
            END IF

            ! GET SIZE DISTRIBUTION PARAMETERS

            ! MELT VERY SMALL SNOW AND GRAUPEL MIXING RATIOS, ADD TO RAIN
            IF (QNI3D(K).LT.1.E-6) THEN
                QR3D(K)=QR3D(K)+QNI3D(K)
                NR3D(K)=NR3D(K)+NS3D(K)
                T3D(K)=T3D(K)-QNI3D(K)*XLF(K)/CPM(K)
                QNI3D(K) = 0.
                NS3D(K) = 0.
                an1d(k,lsc(lr)) = an1d(k,lsc(lr))+an1d(k,lsc(ls))
                an1d(k,lsc(ls)) = 0.
            END IF
            IF (QG3D(K).LT.1.E-6) THEN
                QR3D(K)=QR3D(K)+QG3D(K)
                NR3D(K)=NR3D(K)+NG3D(K)
                T3D(K)=T3D(K)-QG3D(K)*XLF(K)/CPM(K)
                QG3D(K) = 0.
                NG3D(K) = 0.
                an1d(k,lsc(lr)) = an1d(k,lsc(lr))+an1d(k,lsc(lh))
                an1d(k,lsc(lh)) = 0.
            END IF

            IF (QC3D(K).LT.QSMALL.AND.QNI3D(K).LT.1.E-8.AND.QR3D(K).LT.QSMALL.AND.QG3D(K).LT.1.E-8) GOTO 300

            ! MAKE SURE NUMBER CONCENTRATIONS AREN'T NEGATIVE

            NS3D(K) = MAX(0.,NS3D(K))
            NC3D(K) = MAX(0.,NC3D(K))
            NR3D(K) = MAX(0.,NR3D(K))
            NG3D(K) = MAX(0.,NG3D(K))

            !......................................................................
            ! RAIN

            IF (QR3D(K).GE.QSMALL) THEN
                LAMR(K) = (PI*RHOW*NR3D(K)/QR3D(K))**(1./3.)
                N0RR(K) = NR3D(K)*LAMR(K)

                ! CHECK FOR SLOPE

                ! ADJUST VARS

                IF (LAMR(K).LT.LAMMINR) THEN

                    LAMR(K) = LAMMINR

                    N0RR(K) = LAMR(K)**4*QR3D(K)/(PI*RHOW)

                    NR3D(K) = N0RR(K)/LAMR(K)
                ELSE IF (LAMR(K).GT.LAMMAXR) THEN
                    LAMR(K) = LAMMAXR
                    N0RR(K) = LAMR(K)**4*QR3D(K)/(PI*RHOW)

                    NR3D(K) = N0RR(K)/LAMR(K)
                END IF
            END IF

            !......................................................................
            ! CLOUD DROPLETS

            ! MARTIN ET AL. (1994) FORMULA FOR PGAM

            IF (QC3D(K).GE.QSMALL) THEN

                DUM = PRES(K)/(287.15*T3D(K))
                PGAM(K)=0.0005714*(NC3D(K)/1.E6*DUM)+0.2714
                PGAM(K)=1./(PGAM(K)**2)-1.
                PGAM(K)=MAX(PGAM(K),2.)
                PGAM(K)=MIN(PGAM(K),10.)

                ! CALCULATE LAMC

                LAMC(K) = (CONS26*NC3D(K)*GAMMA(PGAM(K)+4.)/   &
                        (QC3D(K)*GAMMA(PGAM(K)+1.)))**(1./3.)

                ! LAMMIN, 60 MICRON DIAMETER
                ! LAMMAX, 1 MICRON

                LAMMIN = (PGAM(K)+1.)/60.E-6
                LAMMAX = (PGAM(K)+1.)/1.E-6

                IF (LAMC(K).LT.LAMMIN) THEN
                    LAMC(K) = LAMMIN

                    NC3D(K) = EXP(3.*LOG(LAMC(K))+LOG(QC3D(K))+              &
                            LOG(GAMMA(PGAM(K)+1.))-LOG(GAMMA(PGAM(K)+4.)))/CONS26
                ELSE IF (LAMC(K).GT.LAMMAX) THEN
                    LAMC(K) = LAMMAX

                    NC3D(K) = EXP(3.*LOG(LAMC(K))+LOG(QC3D(K))+              &
                            LOG(GAMMA(PGAM(K)+1.))-LOG(GAMMA(PGAM(K)+4.)))/CONS26

                END IF
!                xdia(k,lc) = 1.0/LAMC(K)
!                xdia(k,lc) = (PGAM(K)+1.)/LAMC(K)
                xdia(k,lc) = (6/(PI*RHOW)*QC3D(k)/NC3D(k))**(1./3.)
                xdia(k,lc) = min(max(xdia(k,lc),1.E-6),60.E-6)

            END IF

            !......................................................................
            ! SNOW

            IF (QNI3D(K).GE.QSMALL) THEN
                LAMS(K) = (CONS1*NS3D(K)/QNI3D(K))**(1./DS)
                N0S(K) = NS3D(K)*LAMS(K)

                ! CHECK FOR SLOPE

                ! ADJUST VARS

                IF (LAMS(K).LT.LAMMINS) THEN
                    LAMS(K) = LAMMINS
                    N0S(K) = LAMS(K)**4*QNI3D(K)/CONS1

                    NS3D(K) = N0S(K)/LAMS(K)

                ELSE IF (LAMS(K).GT.LAMMAXS) THEN

                    LAMS(K) = LAMMAXS
                    N0S(K) = LAMS(K)**4*QNI3D(K)/CONS1

                    NS3D(K) = N0S(K)/LAMS(K)
                END IF
                xdia(k,ls) = 1.0/LAMS(K)
            END IF

            !......................................................................
            ! GRAUPEL

            IF (QG3D(K).GE.QSMALL) THEN
                LAMG(K) = (CONS2*NG3D(K)/QG3D(K))**(1./DG)
                N0G(K) = NG3D(K)*LAMG(K)

                ! ADJUST VARS

                IF (LAMG(K).LT.LAMMING) THEN
                    LAMG(K) = LAMMING
                    N0G(K) = LAMG(K)**4*QG3D(K)/CONS2

                    NG3D(K) = N0G(K)/LAMG(K)

                ELSE IF (LAMG(K).GT.LAMMAXG) THEN

                    LAMG(K) = LAMMAXG
                    N0G(K) = LAMG(K)**4*QG3D(K)/CONS2

                    NG3D(K) = N0G(K)/LAMG(K)
                END IF
                xdia(k,lh) = 1.0/LAMG(K)
            END IF

            !.....................................................................
            ! ZERO OUT PROCESS RATES

            PRC(K) = 0.
            NPRC(K) = 0.
            NPRC1(K) = 0.
            PRA(K) = 0.
            NPRA(K) = 0.
            NRAGG(K) = 0.
            NSMLTS(K) = 0.
            NSMLTR(K) = 0.
            EVPMS(K) = 0.
            PCC(K) = 0.
            PRE(K) = 0.
            NSUBC(K) = 0.
            NSUBR(K) = 0.
            PRACG(K) = 0.
            NPRACG(K) = 0.
            PSMLT(K) = 0.
            PGMLT(K) = 0.
            EVPMG(K) = 0.
            PRACS(K) = 0.
            NPRACS(K) = 0.
            NGMLTG(K) = 0.
            NGMLTR(K) = 0.
            ! yc: add cold rates
            MNUCCC(K) = 0.
            NNUCCC(K) = 0.
            PRC(K) = 0.
            NPRC(K) = 0.
            NPRC1(K) = 0.
            NSAGG(K) = 0.
            PSACWS(K) = 0.
            NPSACWS(K) = 0.
            PSACWI(K) = 0.
            NPSACWI(K) = 0.
            PRACS(K) = 0.
            NPRACS(K) = 0.
            NMULTS(K) = 0.
            QMULTS(K) = 0.
            NMULTR(K) = 0.
            QMULTR(K) = 0.
            NMULTG(K) = 0.
            QMULTG(K) = 0.
            NMULTRG(K) = 0.
            QMULTRG(K) = 0.
            MNUCCR(K) = 0.
            NNUCCR(K) = 0.
            PRA(K) = 0.
            NPRA(K) = 0.
            NRAGG(K) = 0.
            PRCI(K) = 0.
            NPRCI(K) = 0.
            PRAI(K) = 0.
            NPRAI(K) = 0.
            NNUCCD(K) = 0.
            MNUCCD(K) = 0.
            PCC(K) = 0.
            PRE(K) = 0.
            PRD(K) = 0.
            PRDS(K) = 0.
            EPRD(K) = 0.
            EPRDS(K) = 0.
            NSUBC(K) = 0.
            NSUBI(K) = 0.
            NSUBS(K) = 0.
            NSUBR(K) = 0.
            PIACR(K) = 0.
            NIACR(K) = 0.
            PRACI(K) = 0.
            PIACRS(K) = 0.
            NIACRS(K) = 0.
            PRACIS(K) = 0.
            ! HM: ADD GRAUPEL PROCESSES
            PRACG(K) = 0.
            PSACR(K) = 0.
            PSACWG(K) = 0.
            PGSACW(K) = 0.
            PGRACS(K) = 0.
            PRDG(K) = 0.
            EPRDG(K) = 0.
            NPRACG(K) = 0.
            NPSACWG(K) = 0.
            NSCNG(K) = 0.
            NGRACS(K) = 0.
            NSUBG(K) = 0.

            !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            ! CALCULATION OF MICROPHYSICAL PROCESS RATES, T > 273.15 K

            !.................................................................
            !.......................................................................
            ! AUTOCONVERSION OF CLOUD LIQUID WATER TO RAIN
            ! FORMULA FROM BEHENG (1994)
            ! USING NUMERICAL SIMULATION OF STOCHASTIC COLLECTION EQUATION
            ! AND INITIAL CLOUD DROPLET SIZE DISTRIBUTION SPECIFIED
            ! AS A GAMMA DISTRIBUTION

            ! USE MINIMUM VALUE OF 1.E-6 TO PREVENT FLOATING POINT ERROR

            IF (QC3D(K).GE.1.E-6) THEN

                ! HM ADD 12/13/06, REPLACE WITH NEWER FORMULA
                ! FROM KHAIROUTDINOV AND KOGAN 2000, MWR

                PRC(K)=1350.*QC3D(K)**2.47*  &
                        (NC3D(K)/1.e6*RHO(K))**(-1.79)

                ! note: nprc1 is change in Nr,
                ! nprc is change in Nc

                NPRC1(K) = PRC(K)/CONS29
                NPRC(K) = PRC(K)/(QC3D(k)/NC3D(K))

                ! hm bug fix 3/20/12
                NPRC(K) = MIN(NPRC(K),NC3D(K)/DT)
                NPRC1(K) = MIN(NPRC1(K),NPRC(K))

            END IF

            !.......................................................................
            ! HM ADD 12/13/06, COLLECTION OF SNOW BY RAIN ABOVE FREEZING
            ! FORMULA FROM IKAWA AND SAITO (1991)

            IF (QR3D(K).GE.1.E-8.AND.QNI3D(K).GE.1.E-8) THEN

                UMS = ASN(K)*CONS3/(LAMS(K)**BS)
                UMR = ARN(K)*CONS4/(LAMR(K)**BR)
                UNS = ASN(K)*CONS5/LAMS(K)**BS
                UNR = ARN(K)*CONS6/LAMR(K)**BR

                ! SET REASLISTIC LIMITS ON FALLSPEEDS

                ! bug fix, 10/08/09
                dum=(rhosu/rho(k))**0.54
                UMS=MIN(UMS,1.2*dum)
                UNS=MIN(UNS,1.2*dum)
                UMR=MIN(UMR,9.1*dum)
                UNR=MIN(UNR,9.1*dum)

                ! hm fix, 2/12/13
                ! for above freezing conditions to get accelerated melting of snow,
                ! we need collection of rain by snow (following Lin et al. 1983)
                !            PRACS(K) = CONS31*(((1.2*UMR-0.95*UMS)**2+              &
                !                  0.08*UMS*UMR)**0.5*RHO(K)*                     &
                !                 N0RR(K)*N0S(K)/LAMS(K)**3*                    &
                !                  (5./(LAMS(K)**3*LAMR(K))+                    &
                !                  2./(LAMS(K)**2*LAMR(K)**2)+                  &
                !                  0.5/(LAMS(K)*LAMR(K)**3)))

                PRACS(K) = CONS41*(((1.2*UMR-0.95*UMS)**2+                   &
                        0.08*UMS*UMR)**0.5*RHO(K)*                      &
                        N0RR(K)*N0S(K)/LAMR(K)**3*                              &
                        (5./(LAMR(K)**3*LAMS(K))+                    &
                                2./(LAMR(K)**2*LAMS(K)**2)+                  &
                                0.5/(LAMR(k)*LAMS(k)**3)))

                ! fix 053011, npracs no longer subtracted from snow
                !            NPRACS(K) = CONS32*RHO(K)*(1.7*(UNR-UNS)**2+            &
                !                0.3*UNR*UNS)**0.5*N0RR(K)*N0S(K)*              &
                !                (1./(LAMR(K)**3*LAMS(K))+                      &
                !                 1./(LAMR(K)**2*LAMS(K)**2)+                   &
                !                 1./(LAMR(K)*LAMS(K)**3))

            END IF

            ! ADD COLLECTION OF GRAUPEL BY RAIN ABOVE FREEZING
            ! ASSUME ALL RAIN COLLECTION BY GRAUPEL ABOVE FREEZING IS SHED
            ! ASSUME SHED DROPS ARE 1 MM IN SIZE

            IF (QR3D(K).GE.1.E-8.AND.QG3D(K).GE.1.E-8) THEN

                UMG = AGN(K)*CONS7/(LAMG(K)**BG)
                UMR = ARN(K)*CONS4/(LAMR(K)**BR)
                UNG = AGN(K)*CONS8/LAMG(K)**BG
                UNR = ARN(K)*CONS6/LAMR(K)**BR

                ! SET REASLISTIC LIMITS ON FALLSPEEDS
                ! bug fix, 10/08/09
                dum=(rhosu/rho(k))**0.54
                UMG=MIN(UMG,20.*dum)
                UNG=MIN(UNG,20.*dum)
                UMR=MIN(UMR,9.1*dum)
                UNR=MIN(UNR,9.1*dum)

                ! PRACG IS MIXING RATIO OF RAIN PER SEC COLLECTED BY GRAUPEL/HAIL
                PRACG(K) = CONS41*(((1.2*UMR-0.95*UMG)**2+                   &
                        0.08*UMG*UMR)**0.5*RHO(K)*                      &
                        N0RR(K)*N0G(K)/LAMR(K)**3*                              &
                        (5./(LAMR(K)**3*LAMG(K))+                    &
                                2./(LAMR(K)**2*LAMG(K)**2)+				   &
                                0.5/(LAMR(k)*LAMG(k)**3)))

                ! ASSUME 1 MM DROPS ARE SHED, GET NUMBER SHED PER SEC

                DUM = PRACG(K)/5.2E-7

                NPRACG(K) = CONS32*RHO(K)*(1.7*(UNR-UNG)**2+            &
                        0.3*UNR*UNG)**0.5*N0RR(K)*N0G(K)*              &
                        (1./(LAMR(K)**3*LAMG(K))+                      &
                                1./(LAMR(K)**2*LAMG(K)**2)+                   &
                                1./(LAMR(K)*LAMG(K)**3))

                ! hm 7/15/13, remove limit so that the number of collected drops can smaller than
                ! number of shed drops
                !            NPRACG(K)=MAX(NPRACG(K)-DUM,0.)
                NPRACG(K)=NPRACG(K)-DUM

            END IF

            !.......................................................................
            ! ACCRETION OF CLOUD LIQUID WATER BY RAIN
            ! CONTINUOUS COLLECTION EQUATION WITH
            ! GRAVITATIONAL COLLECTION KERNEL, DROPLET FALL SPEED NEGLECTED

            IF (QR3D(K).GE.1.E-8 .AND. QC3D(K).GE.1.E-8) THEN

                ! 12/13/06 HM ADD, REPLACE WITH NEWER FORMULA FROM
                ! KHAIROUTDINOV AND KOGAN 2000, MWR

                DUM=(QC3D(K)*QR3D(K))
                PRA(K) = 67.*(DUM)**1.15
                NPRA(K) = PRA(K)/(QC3D(K)/NC3D(K))

            END IF
            !.......................................................................
            ! SELF-COLLECTION OF RAIN DROPS
            ! FROM BEHENG(1994)
            ! FROM NUMERICAL SIMULATION OF THE STOCHASTIC COLLECTION EQUATION
            ! AS DESCRINED ABOVE FOR AUTOCONVERSION

            IF (QR3D(K).GE.1.E-8) THEN
                ! include breakup add 10/09/09
                dum1=300.e-6
                if (1./lamr(k).lt.dum1) then
                    dum=1.
                else if (1./lamr(k).ge.dum1) then
                    dum=2.-exp(2300.*(1./lamr(k)-dum1))
                end if
                !            NRAGG(K) = -8.*NR3D(K)*QR3D(K)*RHO(K)
                NRAGG(K) = -5.78*dum*NR3D(K)*QR3D(K)*RHO(K)
            END IF

            !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            ! CALCULATE EVAP OF RAIN (RUTLEDGE AND HOBBS 1983)

            IF (QR3D(K).GE.QSMALL) THEN
                EPSR = 2.*PI*N0RR(K)*RHO(K)*DV(K)*                           &
                        (F1R/(LAMR(K)*LAMR(K))+                       &
                                F2R*(ARN(K)*RHO(K)/MU(K))**0.5*                      &
                                        SC(K)**(1./3.)*CONS9/                   &
                                        (LAMR(K)**CONS34))
            ELSE
                EPSR = 0.
            END IF

            ! NO CONDENSATION ONTO RAIN, ONLY EVAP ALLOWED

            IF (QV3D(K).LT.QVS(K)) THEN
                PRE(K) = EPSR*(QV3D(K)-QVS(K))/AB(K)
                PRE(K) = MIN(PRE(K),0.)
            ELSE
                PRE(K) = 0.
            END IF

            !.......................................................................
            ! MELTING OF SNOW

            ! SNOW MAY PERSITS ABOVE FREEZING, FORMULA FROM RUTLEDGE AND HOBBS, 1984
            ! IF WATER SUPERSATURATION, SNOW MELTS TO FORM RAIN

            IF (QNI3D(K).GE.1.E-8) THEN

                ! fix 053011
                ! HM, MODIFY FOR V3.2, ADD ACCELERATED MELTING DUE TO COLLISION WITH RAIN
                !             DUM = -CPW/XLF(K)*T3D(K)*PRACS(K)
                DUM = -CPW/XLF(K)*(T3D(K)-273.15)*PRACS(K)

                ! hm fix 1/20/15
                !             PSMLT(K)=2.*PI*N0S(K)*KAP(K)*(273.15-T3D(K))/       &
                !                    XLF(K)*RHO(K)*(F1S/(LAMS(K)*LAMS(K))+        &
                !                    F2S*(ASN(K)*RHO(K)/MU(K))**0.5*                      &
                !                    SC(K)**(1./3.)*CONS10/                   &
                !                   (LAMS(K)**CONS35))+DUM
                PSMLT(K)=2.*PI*N0S(K)*KAP(K)*(273.15-T3D(K))/       &
                        XLF(K)*(F1S/(LAMS(K)*LAMS(K))+        &
                        F2S*(ASN(K)*RHO(K)/MU(K))**0.5*                      &
                                SC(K)**(1./3.)*CONS10/                   &
                                (LAMS(K)**CONS35))+DUM

                ! IN WATER SUBSATURATION, SNOW MELTS AND EVAPORATES

                IF (QVQVS(K).LT.1.) THEN
                    EPSS = 2.*PI*N0S(K)*RHO(K)*DV(K)*                            &
                            (F1S/(LAMS(K)*LAMS(K))+                       &
                                    F2S*(ASN(K)*RHO(K)/MU(K))**0.5*                      &
                                            SC(K)**(1./3.)*CONS10/                   &
                                            (LAMS(K)**CONS35))
                    ! hm fix 8/4/08
                    EVPMS(K) = (QV3D(K)-QVS(K))*EPSS/AB(K)
                    EVPMS(K) = MAX(EVPMS(K),PSMLT(K))
                    PSMLT(K) = PSMLT(K)-EVPMS(K)
                END IF
            END IF

            !.......................................................................
            ! MELTING OF GRAUPEL

            ! GRAUPEL MAY PERSITS ABOVE FREEZING, FORMULA FROM RUTLEDGE AND HOBBS, 1984
            ! IF WATER SUPERSATURATION, GRAUPEL MELTS TO FORM RAIN

            IF (QG3D(K).GE.1.E-8) THEN

                ! fix 053011
                ! HM, MODIFY FOR V3.2, ADD ACCELERATED MELTING DUE TO COLLISION WITH RAIN
                !             DUM = -CPW/XLF(K)*T3D(K)*PRACG(K)
                DUM = -CPW/XLF(K)*(T3D(K)-273.15)*PRACG(K)

                ! hm fix 1/20/15
                !             PGMLT(K)=2.*PI*N0G(K)*KAP(K)*(273.15-T3D(K))/ 		 &
                !                    XLF(K)*RHO(K)*(F1S/(LAMG(K)*LAMG(K))+                &
                !                    F2S*(AGN(K)*RHO(K)/MU(K))**0.5*                      &
                !                    SC(K)**(1./3.)*CONS11/                   &
                !                   (LAMG(K)**CONS36))+DUM
                PGMLT(K)=2.*PI*N0G(K)*KAP(K)*(273.15-T3D(K))/ 		 &
                        XLF(K)*(F1S/(LAMG(K)*LAMG(K))+                &
                        F2S*(AGN(K)*RHO(K)/MU(K))**0.5*                      &
                                SC(K)**(1./3.)*CONS11/                   &
                                (LAMG(K)**CONS36))+DUM

                ! IN WATER SUBSATURATION, GRAUPEL MELTS AND EVAPORATES

                IF (QVQVS(K).LT.1.) THEN
                    EPSG = 2.*PI*N0G(K)*RHO(K)*DV(K)*                                &
                            (F1S/(LAMG(K)*LAMG(K))+                               &
                                    F2S*(AGN(K)*RHO(K)/MU(K))**0.5*                      &
                                            SC(K)**(1./3.)*CONS11/                   &
                                            (LAMG(K)**CONS36))
                    ! hm fix 8/4/08
                    EVPMG(K) = (QV3D(K)-QVS(K))*EPSG/AB(K)
                    EVPMG(K) = MAX(EVPMG(K),PGMLT(K))
                    PGMLT(K) = PGMLT(K)-EVPMG(K)
                END IF
            END IF

            ! HM, V3.2
            ! RESET PRACG AND PRACS TO ZERO, THIS IS DONE BECAUSE THERE IS NO
            ! TRANSFER OF MASS FROM SNOW AND GRAUPEL TO RAIN DIRECTLY FROM COLLECTION
            ! ABOVE FREEZING, IT IS ONLY USED FOR ENHANCEMENT OF MELTING AND SHEDDING

            PRACG(K) = 0.
            PRACS(K) = 0.

            !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

            ! FOR CLOUD ICE, ONLY PROCESSES OPERATING AT T > 273.15 IS
            ! MELTING, WHICH IS ALREADY CONSERVED DURING PROCESS
            ! CALCULATION

            ! CONSERVATION OF QC

            DUM = (PRC(K)+PRA(K))*DT

            IF (DUM.GT.QC3D(K).AND.QC3D(K).GE.QSMALL) THEN

                RATIO = QC3D(K)/DUM

                PRC(K) = PRC(K)*RATIO
                PRA(K) = PRA(K)*RATIO

            END IF

            ! CONSERVATION OF SNOW

            DUM = (-PSMLT(K)-EVPMS(K)+PRACS(K))*DT

            IF (DUM.GT.QNI3D(K).AND.QNI3D(K).GE.QSMALL) THEN

                ! NO SOURCE TERMS FOR SNOW AT T > FREEZING
                RATIO = QNI3D(K)/DUM

                PSMLT(K) = PSMLT(K)*RATIO
                EVPMS(K) = EVPMS(K)*RATIO
                PRACS(K) = PRACS(K)*RATIO

            END IF

            ! CONSERVATION OF GRAUPEL

            DUM = (-PGMLT(K)-EVPMG(K)+PRACG(K))*DT

            IF (DUM.GT.QG3D(K).AND.QG3D(K).GE.QSMALL) THEN

                ! NO SOURCE TERM FOR GRAUPEL ABOVE FREEZING
                RATIO = QG3D(K)/DUM

                PGMLT(K) = PGMLT(K)*RATIO
                EVPMG(K) = EVPMG(K)*RATIO
                PRACG(K) = PRACG(K)*RATIO

            END IF

            ! CONSERVATION OF QR
            ! HM 12/13/06, ADDED CONSERVATION OF RAIN SINCE PRE IS NEGATIVE

            DUM = (-PRACS(K)-PRACG(K)-PRE(K)-PRA(K)-PRC(K)+PSMLT(K)+PGMLT(K))*DT

            IF (DUM.GT.QR3D(K).AND.QR3D(K).GE.QSMALL) THEN

                RATIO = (QR3D(K)/DT+PRACS(K)+PRACG(K)+PRA(K)+PRC(K)-PSMLT(K)-PGMLT(K))/ &
                        (-PRE(K))
                PRE(K) = PRE(K)*RATIO

            END IF

            !....................................

            QV3DTEN(K) = QV3DTEN(K)+(-PRE(K)-EVPMS(K)-EVPMG(K))

            T3DTEN(K) = T3DTEN(K)+(PRE(K)*XXLV(K)+(EVPMS(K)+EVPMG(K))*XXLS(K)+&
                    (PSMLT(K)+PGMLT(K)-PRACS(K)-PRACG(K))*XLF(K))/CPM(K)

            QC3DTEN(K) = QC3DTEN(K)+(-PRA(K)-PRC(K))
            QR3DTEN(K) = QR3DTEN(K)+(PRE(K)+PRA(K)+PRC(K)-PSMLT(K)-PGMLT(K)+PRACS(K)+PRACG(K))
            QNI3DTEN(K) = QNI3DTEN(K)+(PSMLT(K)+EVPMS(K)-PRACS(K))
            QG3DTEN(K) = QG3DTEN(K)+(PGMLT(K)+EVPMG(K)-PRACG(K))
            ! fix 053011
            !      NS3DTEN(K) = NS3DTEN(K)-NPRACS(K)
            ! HM, bug fix 5/12/08, npracg is subtracted from nr not ng
            !      NG3DTEN(K) = NG3DTEN(K)
            NC3DTEN(K) = NC3DTEN(K)+ (-NPRA(K)-NPRC(K))
            NR3DTEN(K) = NR3DTEN(K)+ (NPRC1(K)+NRAGG(K)-NPRACG(K))

            ! HM ADD, WRF-CHEM, ADD TENDENCIES FOR C2PREC

            C2PREC(K) = PRA(K)+PRC(K)
            IF (PRE(K).LT.0.) THEN
                DUM = PRE(K)*DT/QR3D(K)
                DUM = MAX(-1.,DUM)
                NSUBR(K) = DUM*NR3D(K)/DT
            END IF

            IF (EVPMS(K)+PSMLT(K).LT.0.) THEN
                DUM = (EVPMS(K)+PSMLT(K))*DT/QNI3D(K)
                DUM = MAX(-1.,DUM)
                NSMLTS(K) = DUM*NS3D(K)/DT
            END IF
            IF (PSMLT(K).LT.0.) THEN
                DUM = PSMLT(K)*DT/QNI3D(K)
                DUM = MAX(-1.0,DUM)
                NSMLTR(K) = DUM*NS3D(K)/DT
            END IF
            IF (EVPMG(K)+PGMLT(K).LT.0.) THEN
                DUM = (EVPMG(K)+PGMLT(K))*DT/QG3D(K)
                DUM = MAX(-1.,DUM)
                NGMLTG(K) = DUM*NG3D(K)/DT
            END IF
            IF (PGMLT(K).LT.0.) THEN
                DUM = PGMLT(K)*DT/QG3D(K)
                DUM = MAX(-1.0,DUM)
                NGMLTR(K) = DUM*NG3D(K)/DT
            END IF

            NS3DTEN(K) = NS3DTEN(K)+(NSMLTS(K))
            NG3DTEN(K) = NG3DTEN(K)+(NGMLTG(K))
            NR3DTEN(K) = NR3DTEN(K)+(NSUBR(K)-NSMLTR(K)-NGMLTR(K))

            300  CONTINUE

            !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            ! NOW CALCULATE SATURATION ADJUSTMENT TO CONDENSE EXTRA VAPOR ABOVE
            ! WATER SATURATION

            DUMT = T3D(K)+DT*T3DTEN(K)
            DUMQV = QV3D(K)+DT*QV3DTEN(K)
            ! hm, add fix for low pressure, 5/12/10
            dum=min(0.99*pres(k),POLYSVP(DUMT,0))
            DUMQSS = EP_2*dum/(PRES(K)-dum)
            DUMQC = QC3D(K)+DT*QC3DTEN(K)
            DUMQC = MAX(DUMQC,0.)

            ! SATURATION ADJUSTMENT FOR LIQUID

            DUMS = DUMQV-DUMQSS
            PCC(K) = DUMS/(1.+XXLV(K)**2*DUMQSS/(CPM(K)*RV*DUMT**2))/DT
            IF (PCC(K)*DT+DUMQC.LT.0.) THEN
                PCC(K) = -DUMQC/DT
            END IF

            QV3DTEN(K) = QV3DTEN(K)-PCC(K)
            T3DTEN(K) = T3DTEN(K)+PCC(K)*XXLV(K)/CPM(K)
            QC3DTEN(K) = QC3DTEN(K)+PCC(K)

#if (WRF_CHEM == 1)
            evapprod(k) = - PRE(K) - EVPMS(K) - EVPMG(K)
            rainprod(k) = PRA(K) + PRC(K) + tqimelt(K)
#endif

            !.......................................................................
            ! ACTIVATION OF CLOUD DROPLETS
            ! ACTIVATION OF DROPLET CURRENTLY NOT CALCULATED
            ! DROPLET CONCENTRATION IS SPECIFIED !!!!!

            !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            ! SUBLIMATE, MELT, OR EVAPORATE NUMBER CONCENTRATION
            ! THIS FORMULATION ASSUMES 1:1 RATIO BETWEEN MASS LOSS AND
            ! LOSS OF NUMBER CONCENTRATION

            !     IF (PCC(K).LT.0.) THEN
            !        DUM = PCC(K)*DT/QC3D(K)
            !           DUM = MAX(-1.,DUM)
            !        NSUBC(K) = DUM*NC3D(K)/DT
            !     END IF

            ! UPDATE TENDENCIES

            !        NC3DTEN(K) = NC3DTEN(K)+NSUBC(K)

            !.....................................................................
            !.....................................................................
        ELSE  ! TEMPERATURE < 273.15

            !......................................................................
            !HM ADD, ALLOW FOR CONSTANT DROPLET NUMBER
            ! INUM = 0, PREDICT DROPLET NUMBER
            ! INUM = 1, SET CONSTANT DROPLET NUMBER

            IF (iinum.EQ.1) THEN
                ! CONVERT NDCNST FROM CM-3 TO KG-1
                NC3D(K)=NDCNST*1.E6/RHO(K)
            END IF

            ! CALCULATE SIZE DISTRIBUTION PARAMETERS
            ! MAKE SURE NUMBER CONCENTRATIONS AREN'T NEGATIVE

            NI3D(K) = MAX(0.,NI3D(K))
            NS3D(K) = MAX(0.,NS3D(K))
            NC3D(K) = MAX(0.,NC3D(K))
            NR3D(K) = MAX(0.,NR3D(K))
            NG3D(K) = MAX(0.,NG3D(K))

            !......................................................................
            ! CLOUD ICE

            IF (QI3D(K).GE.QSMALL) THEN
                LAMI(K) = (CONS12*                 &
                        NI3D(K)/QI3D(K))**(1./DI)
                N0I(K) = NI3D(K)*LAMI(K)

                ! CHECK FOR SLOPE

                ! ADJUST VARS

                IF (LAMI(K).LT.LAMMINI) THEN

                    LAMI(K) = LAMMINI

                    N0I(K) = LAMI(K)**4*QI3D(K)/CONS12

                    NI3D(K) = N0I(K)/LAMI(K)
                ELSE IF (LAMI(K).GT.LAMMAXI) THEN
                    LAMI(K) = LAMMAXI
                    N0I(K) = LAMI(K)**4*QI3D(K)/CONS12

                    NI3D(K) = N0I(K)/LAMI(K)
                END IF
                xdia(k,li) = 1.0/LAMI(K)
            END IF

            !......................................................................
            ! RAIN

            IF (QR3D(K).GE.QSMALL) THEN
                LAMR(K) = (PI*RHOW*NR3D(K)/QR3D(K))**(1./3.)
                N0RR(K) = NR3D(K)*LAMR(K)

                ! CHECK FOR SLOPE

                ! ADJUST VARS

                IF (LAMR(K).LT.LAMMINR) THEN

                    LAMR(K) = LAMMINR

                    N0RR(K) = LAMR(K)**4*QR3D(K)/(PI*RHOW)

                    NR3D(K) = N0RR(K)/LAMR(K)
                ELSE IF (LAMR(K).GT.LAMMAXR) THEN
                    LAMR(K) = LAMMAXR
                    N0RR(K) = LAMR(K)**4*QR3D(K)/(PI*RHOW)

                    NR3D(K) = N0RR(K)/LAMR(K)
                END IF
            END IF

            !......................................................................
            ! CLOUD DROPLETS

            ! MARTIN ET AL. (1994) FORMULA FOR PGAM

            IF (QC3D(K).GE.QSMALL) THEN

                DUM = PRES(K)/(287.15*T3D(K))
                PGAM(K)=0.0005714*(NC3D(K)/1.E6*DUM)+0.2714
                PGAM(K)=1./(PGAM(K)**2)-1.
                PGAM(K)=MAX(PGAM(K),2.)
                PGAM(K)=MIN(PGAM(K),10.)

                ! CALCULATE LAMC

                LAMC(K) = (CONS26*NC3D(K)*GAMMA(PGAM(K)+4.)/   &
                        (QC3D(K)*GAMMA(PGAM(K)+1.)))**(1./3.)

                ! LAMMIN, 60 MICRON DIAMETER
                ! LAMMAX, 1 MICRON

                LAMMIN = (PGAM(K)+1.)/60.E-6
                LAMMAX = (PGAM(K)+1.)/1.E-6

                IF (LAMC(K).LT.LAMMIN) THEN
                    LAMC(K) = LAMMIN

                    NC3D(K) = EXP(3.*LOG(LAMC(K))+LOG(QC3D(K))+              &
                            LOG(GAMMA(PGAM(K)+1.))-LOG(GAMMA(PGAM(K)+4.)))/CONS26
                ELSE IF (LAMC(K).GT.LAMMAX) THEN
                    LAMC(K) = LAMMAX
                    NC3D(K) = EXP(3.*LOG(LAMC(K))+LOG(QC3D(K))+              &
                            LOG(GAMMA(PGAM(K)+1.))-LOG(GAMMA(PGAM(K)+4.)))/CONS26

                END IF
!                xdia(k,lc) = 1.0/LAMC(K)
!                xdia(k,lc) = (PGAM(K)+1.)/LAMC(K)
                xdia(k,lc) = (6/(PI*RHOW)*QC3D(k)/NC3D(k))**(1./3.)
                xdia(k,lc) = min(max(xdia(k,lc),1.E-6),60.E-6)

                ! TO CALCULATE DROPLET FREEZING

                CDIST1(K) = NC3D(K)/GAMMA(PGAM(K)+1.)

            END IF

            !......................................................................
            ! SNOW

            IF (QNI3D(K).GE.QSMALL) THEN
                LAMS(K) = (CONS1*NS3D(K)/QNI3D(K))**(1./DS)
                N0S(K) = NS3D(K)*LAMS(K)

                ! CHECK FOR SLOPE

                ! ADJUST VARS

                IF (LAMS(K).LT.LAMMINS) THEN
                    LAMS(K) = LAMMINS
                    N0S(K) = LAMS(K)**4*QNI3D(K)/CONS1

                    NS3D(K) = N0S(K)/LAMS(K)

                ELSE IF (LAMS(K).GT.LAMMAXS) THEN

                    LAMS(K) = LAMMAXS
                    N0S(K) = LAMS(K)**4*QNI3D(K)/CONS1

                    NS3D(K) = N0S(K)/LAMS(K)
                END IF
                xdia(k,ls) = 1.0/LAMS(K)
            END IF

            !......................................................................
            ! GRAUPEL

            IF (QG3D(K).GE.QSMALL) THEN
                LAMG(K) = (CONS2*NG3D(K)/QG3D(K))**(1./DG)
                N0G(K) = NG3D(K)*LAMG(K)

                ! CHECK FOR SLOPE

                ! ADJUST VARS

                IF (LAMG(K).LT.LAMMING) THEN
                    LAMG(K) = LAMMING
                    N0G(K) = LAMG(K)**4*QG3D(K)/CONS2

                    NG3D(K) = N0G(K)/LAMG(K)

                ELSE IF (LAMG(K).GT.LAMMAXG) THEN

                    LAMG(K) = LAMMAXG
                    N0G(K) = LAMG(K)**4*QG3D(K)/CONS2

                    NG3D(K) = N0G(K)/LAMG(K)
                END IF
                xdia(k,lh) = 1.0/LAMG(K)
            END IF

            !.....................................................................
            ! ZERO OUT PROCESS RATES

            MNUCCC(K) = 0.
            NNUCCC(K) = 0.
            PRC(K) = 0.
            NPRC(K) = 0.
            NPRC1(K) = 0.
            NSAGG(K) = 0.
            PSACWS(K) = 0.
            NPSACWS(K) = 0.
            PSACWI(K) = 0.
            NPSACWI(K) = 0.
            PRACS(K) = 0.
            NPRACS(K) = 0.
            NMULTS(K) = 0.
            QMULTS(K) = 0.
            NMULTR(K) = 0.
            QMULTR(K) = 0.
            NMULTG(K) = 0.
            QMULTG(K) = 0.
            NMULTRG(K) = 0.
            QMULTRG(K) = 0.
            MNUCCR(K) = 0.
            NNUCCR(K) = 0.
            PRA(K) = 0.
            NPRA(K) = 0.
            NRAGG(K) = 0.
            PRCI(K) = 0.
            NPRCI(K) = 0.
            PRAI(K) = 0.
            NPRAI(K) = 0.
            NNUCCD(K) = 0.
            MNUCCD(K) = 0.
            PCC(K) = 0.
            PRE(K) = 0.
            PRD(K) = 0.
            PRDS(K) = 0.
            EPRD(K) = 0.
            EPRDS(K) = 0.
            NSUBC(K) = 0.
            NSUBI(K) = 0.
            NSUBS(K) = 0.
            NSUBR(K) = 0.
            PIACR(K) = 0.
            NIACR(K) = 0.
            PRACI(K) = 0.
            PIACRS(K) = 0.
            NIACRS(K) = 0.
            PRACIS(K) = 0.
            ! HM: ADD GRAUPEL PROCESSES
            PRACG(K) = 0.
            PSACR(K) = 0.
            PSACWG(K) = 0.
            PGSACW(K) = 0.
            PGRACS(K) = 0.
            PRDG(K) = 0.
            EPRDG(K) = 0.
            NPRACG(K) = 0.
            NPSACWG(K) = 0.
            NSCNG(K) = 0.
            NGRACS(K) = 0.
            NSUBG(K) = 0.
            ! yc: add warm rates
            PRC(K) = 0.
            NPRC(K) = 0.
            NPRC1(K) = 0.
            PRA(K) = 0.
            NPRA(K) = 0.
            NRAGG(K) = 0.
            NSMLTS(K) = 0.
            NSMLTR(K) = 0.
            EVPMS(K) = 0.
            PCC(K) = 0.
            PRE(K) = 0.
            NSUBC(K) = 0.
            NSUBR(K) = 0.
            PRACG(K) = 0.
            NPRACG(K) = 0.
            PSMLT(K) = 0.
            PGMLT(K) = 0.
            EVPMG(K) = 0.
            PRACS(K) = 0.
            NPRACS(K) = 0.
            NGMLTG(K) = 0.
            NGMLTR(K) = 0.

            !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            ! CALCULATION OF MICROPHYSICAL PROCESS RATES
            ! ACCRETION/AUTOCONVERSION/FREEZING/MELTING/COAG.
            !.......................................................................
            ! FREEZING OF CLOUD DROPLETS
            ! ONLY ALLOWED BELOW -4 C
            IF (QC3D(K).GE.QSMALL .AND. T3D(K).LT.269.15) THEN

                ! NUMBER OF CONTACT NUCLEI (M^-3) FROM MEYERS ET AL., 1992
                ! FACTOR OF 1000 IS TO CONVERT FROM L^-1 TO M^-3

                ! MEYERS CURVE

                NACNT = EXP(-2.80+0.262*(273.15-T3D(K)))*1000.

                ! COOPER CURVE
                !        NACNT =  5.*EXP(0.304*(273.15-T3D(K)))

                ! FLECTHER
                !     NACNT = 0.01*EXP(0.6*(273.15-T3D(K)))

                ! CONTACT FREEZING

                ! MEAN FREE PATH

                DUM = 7.37*T3D(K)/(288.*10.*PRES(K))/100.

                ! EFFECTIVE DIFFUSIVITY OF CONTACT NUCLEI
                ! BASED ON BROWNIAN DIFFUSION

                DAP(K) = CONS37*T3D(K)*(1.+DUM/RIN)/MU(K)

                MNUCCC(K) = CONS38*DAP(K)*NACNT*EXP(LOG(CDIST1(K))+   &
                        LOG(GAMMA(PGAM(K)+5.))-4.*LOG(LAMC(K)))
                NNUCCC(K) = 2.*PI*DAP(K)*NACNT*CDIST1(K)*           &
                        GAMMA(PGAM(K)+2.)/                         &
                        LAMC(K)

                ! IMMERSION FREEZING (BIGG 1953)

                !           MNUCCC(K) = MNUCCC(K)+CONS39*                   &
                !                  EXP(LOG(CDIST1(K))+LOG(GAMMA(7.+PGAM(K)))-6.*LOG(LAMC(K)))*             &
                !                   EXP(AIMM*(273.15-T3D(K)))

                !           NNUCCC(K) = NNUCCC(K)+                                  &
                !            CONS40*EXP(LOG(CDIST1(K))+LOG(GAMMA(PGAM(K)+4.))-3.*LOG(LAMC(K)))              &
                !                *EXP(AIMM*(273.15-T3D(K)))

                ! hm 7/15/13 fix for consistency w/ original formula
                MNUCCC(K) = MNUCCC(K)+CONS39*                   &
                        EXP(LOG(CDIST1(K))+LOG(GAMMA(7.+PGAM(K)))-6.*LOG(LAMC(K)))*             &
                        (EXP(AIMM*(273.15-T3D(K)))-1.)

                NNUCCC(K) = NNUCCC(K)+                                  &
                        CONS40*EXP(LOG(CDIST1(K))+LOG(GAMMA(PGAM(K)+4.))-3.*LOG(LAMC(K)))              &
                                *(EXP(AIMM*(273.15-T3D(K)))-1.)

                ! PUT IN A CATCH HERE TO PREVENT DIVERGENCE BETWEEN NUMBER CONC. AND
                ! MIXING RATIO, SINCE STRICT CONSERVATION NOT CHECKED FOR NUMBER CONC

                NNUCCC(K) = MIN(NNUCCC(K),NC3D(K)/DT)

            END IF

            !.................................................................
            !.......................................................................
            ! AUTOCONVERSION OF CLOUD LIQUID WATER TO RAIN
            ! FORMULA FROM BEHENG (1994)
            ! USING NUMERICAL SIMULATION OF STOCHASTIC COLLECTION EQUATION
            ! AND INITIAL CLOUD DROPLET SIZE DISTRIBUTION SPECIFIED
            ! AS A GAMMA DISTRIBUTION

            ! USE MINIMUM VALUE OF 1.E-6 TO PREVENT FLOATING POINT ERROR

            IF (QC3D(K).GE.1.E-6) THEN

                ! HM ADD 12/13/06, REPLACE WITH NEWER FORMULA
                ! FROM KHAIROUTDINOV AND KOGAN 2000, MWR

                PRC(K)=1350.*QC3D(K)**2.47*  &
                        (NC3D(K)/1.e6*RHO(K))**(-1.79)

                ! note: nprc1 is change in Nr,
                ! nprc is change in Nc

                NPRC1(K) = PRC(K)/CONS29
                NPRC(K) = PRC(K)/(QC3D(K)/NC3D(K))

                ! hm bug fix 3/20/12
                NPRC(K) = MIN(NPRC(K),NC3D(K)/DT)
                NPRC1(K) = MIN(NPRC1(K),NPRC(K))

            END IF

            !.......................................................................
            ! SELF-COLLECTION OF DROPLET NOT INCLUDED IN KK2000 SCHEME

            ! SNOW AGGREGATION FROM PASSARELLI, 1978, USED BY REISNER, 1998
            ! THIS IS HARD-WIRED FOR BS = 0.4 FOR NOW

            IF (QNI3D(K).GE.1.E-8) THEN
                NSAGG(K) = CONS15*ASN(K)*RHO(K)**            &
                        ((2.+BS)/3.)*QNI3D(K)**((2.+BS)/3.)*                  &
                        (NS3D(K)*RHO(K))**((4.-BS)/3.)/                       &
                        (RHO(K))
            END IF

            !.......................................................................
            ! ACCRETION OF CLOUD DROPLETS ONTO SNOW/GRAUPEL
            ! HERE USE CONTINUOUS COLLECTION EQUATION WITH
            ! SIMPLE GRAVITATIONAL COLLECTION KERNEL IGNORING

            ! SNOW

            IF (QNI3D(K).GE.1.E-8 .AND. QC3D(K).GE.QSMALL) THEN

                PSACWS(K) = CONS13*ASN(K)*QC3D(K)*RHO(K)*               &
                        N0S(K)/                        &
                        LAMS(K)**(BS+3.)
                NPSACWS(K) = CONS13*ASN(K)*NC3D(K)*RHO(K)*              &
                        N0S(K)/                        &
                        LAMS(K)**(BS+3.)

            END IF

            !............................................................................
            ! COLLECTION OF CLOUD WATER BY GRAUPEL

            IF (QG3D(K).GE.1.E-8 .AND. QC3D(K).GE.QSMALL) THEN

                PSACWG(K) = CONS14*AGN(K)*QC3D(K)*RHO(K)*               &
                        N0G(K)/                        &
                        LAMG(K)**(BG+3.)
                NPSACWG(K) = CONS14*AGN(K)*NC3D(K)*RHO(K)*              &
                        N0G(K)/                        &
                        LAMG(K)**(BG+3.)
            END IF

            !.......................................................................
            ! HM, ADD 12/13/06
            ! CLOUD ICE COLLECTING DROPLETS, ASSUME THAT CLOUD ICE MEAN DIAM > 100 MICRON
            ! BEFORE RIMING CAN OCCUR
            ! ASSUME THAT RIME COLLECTED ON CLOUD ICE DOES NOT LEAD
            ! TO HALLET-MOSSOP SPLINTERING

            IF (QI3D(K).GE.1.E-8 .AND. QC3D(K).GE.QSMALL) THEN

                ! PUT IN SIZE DEPENDENT COLLECTION EFFICIENCY BASED ON STOKES LAW
                ! FROM THOMPSON ET AL. 2004, MWR

                IF (1./LAMI(K).GE.100.E-6) THEN

                    PSACWI(K) = CONS16*AIN(K)*QC3D(K)*RHO(K)*               &
                            N0I(K)/                        &
                            LAMI(K)**(BI+3.)
                    NPSACWI(K) = CONS16*AIN(K)*NC3D(K)*RHO(K)*              &
                            N0I(K)/                        &
                            LAMI(K)**(BI+3.)
                END IF
            END IF

            !.......................................................................
            ! ACCRETION OF RAIN WATER BY SNOW
            ! FORMULA FROM IKAWA AND SAITO, 1991, USED BY REISNER ET AL, 1998

            IF (QR3D(K).GE.1.E-8.AND.QNI3D(K).GE.1.E-8) THEN

                UMS = ASN(K)*CONS3/(LAMS(K)**BS)
                UMR = ARN(K)*CONS4/(LAMR(K)**BR)
                UNS = ASN(K)*CONS5/LAMS(K)**BS
                UNR = ARN(K)*CONS6/LAMR(K)**BR

                ! SET REASLISTIC LIMITS ON FALLSPEEDS

                ! bug fix, 10/08/09
                dum=(rhosu/rho(k))**0.54
                UMS=MIN(UMS,1.2*dum)
                UNS=MIN(UNS,1.2*dum)
                UMR=MIN(UMR,9.1*dum)
                UNR=MIN(UNR,9.1*dum)

                PRACS(K) = CONS41*(((1.2*UMR-0.95*UMS)**2+                   &
                        0.08*UMS*UMR)**0.5*RHO(K)*                      &
                        N0RR(K)*N0S(K)/LAMR(K)**3*                              &
                        (5./(LAMR(K)**3*LAMS(K))+                    &
                                2./(LAMR(K)**2*LAMS(K)**2)+                  &
                                0.5/(LAMR(k)*LAMS(k)**3)))

                NPRACS(K) = CONS32*RHO(K)*(1.7*(UNR-UNS)**2+            &
                        0.3*UNR*UNS)**0.5*N0RR(K)*N0S(K)*              &
                        (1./(LAMR(K)**3*LAMS(K))+                      &
                                1./(LAMR(K)**2*LAMS(K)**2)+                   &
                                1./(LAMR(K)*LAMS(K)**3))

                ! MAKE SURE PRACS DOESN'T EXCEED TOTAL RAIN MIXING RATIO
                ! AS THIS MAY OTHERWISE RESULT IN TOO MUCH TRANSFER OF WATER DURING
                ! RIME-SPLINTERING

                PRACS(K) = MIN(PRACS(K),QR3D(K)/DT)

                ! COLLECTION OF SNOW BY RAIN - NEEDED FOR GRAUPEL CONVERSION CALCULATIONS
                ! ONLY CALCULATE IF SNOW AND RAIN MIXING RATIOS EXCEED 0.1 G/KG

                ! HM MODIFY FOR WRFV3.1
                !            IF (IHAIL.EQ.0) THEN
                IF (QNI3D(K).GE.0.1E-3.AND.QR3D(K).GE.0.1E-3) THEN
                    PSACR(K) = CONS31*(((1.2*UMR-0.95*UMS)**2+              &
                            0.08*UMS*UMR)**0.5*RHO(K)*                     &
                            N0RR(K)*N0S(K)/LAMS(K)**3*                               &
                            (5./(LAMS(K)**3*LAMR(K))+                    &
                                    2./(LAMS(K)**2*LAMR(K)**2)+                  &
                                    0.5/(LAMS(K)*LAMR(K)**3)))
                END IF
                !            END IF

            END IF

            !.......................................................................

            ! COLLECTION OF RAINWATER BY GRAUPEL, FROM IKAWA AND SAITO 1990,
            ! USED BY REISNER ET AL 1998
            IF (QR3D(K).GE.1.E-8.AND.QG3D(K).GE.1.E-8) THEN

                UMG = AGN(K)*CONS7/(LAMG(K)**BG)
                UMR = ARN(K)*CONS4/(LAMR(K)**BR)
                UNG = AGN(K)*CONS8/LAMG(K)**BG
                UNR = ARN(K)*CONS6/LAMR(K)**BR

                ! SET REASLISTIC LIMITS ON FALLSPEEDS
                ! bug fix, 10/08/09
                dum=(rhosu/rho(k))**0.54
                UMG=MIN(UMG,20.*dum)
                UNG=MIN(UNG,20.*dum)
                UMR=MIN(UMR,9.1*dum)
                UNR=MIN(UNR,9.1*dum)

                PRACG(K) = CONS41*(((1.2*UMR-0.95*UMG)**2+                   &
                        0.08*UMG*UMR)**0.5*RHO(K)*                      &
                        N0RR(K)*N0G(K)/LAMR(K)**3*                              &
                        (5./(LAMR(K)**3*LAMG(K))+                    &
                                2./(LAMR(K)**2*LAMG(K)**2)+				   &
                                0.5/(LAMR(k)*LAMG(k)**3)))

                NPRACG(K) = CONS32*RHO(K)*(1.7*(UNR-UNG)**2+            &
                        0.3*UNR*UNG)**0.5*N0RR(K)*N0G(K)*              &
                        (1./(LAMR(K)**3*LAMG(K))+                      &
                                1./(LAMR(K)**2*LAMG(K)**2)+                   &
                                1./(LAMR(K)*LAMG(K)**3))

                ! MAKE SURE PRACG DOESN'T EXCEED TOTAL RAIN MIXING RATIO
                ! AS THIS MAY OTHERWISE RESULT IN TOO MUCH TRANSFER OF WATER DURING
                ! RIME-SPLINTERING

                PRACG(K) = MIN(PRACG(K),QR3D(K)/DT)

            END IF

            !.......................................................................
            ! RIME-SPLINTERING - SNOW
            ! HALLET-MOSSOP (1974)
            ! NUMBER OF SPLINTERS FORMED IS BASED ON MASS OF RIMED WATER

            ! DUM1 = MASS OF INDIVIDUAL SPLINTERS

            ! HM ADD THRESHOLD SNOW AND DROPLET MIXING RATIO FOR RIME-SPLINTERING
            ! TO LIMIT RIME-SPLINTERING IN STRATIFORM CLOUDS
            ! THESE THRESHOLDS CORRESPOND WITH GRAUPEL THRESHOLDS IN RH 1984

            !v1.4
            IF (QNI3D(K).GE.0.1E-3) THEN
                IF (QC3D(K).GE.0.5E-3.OR.QR3D(K).GE.0.1E-3) THEN
                    IF (PSACWS(K).GT.0..OR.PRACS(K).GT.0.) THEN
                        IF (T3D(K).LT.270.16 .AND. T3D(K).GT.265.16) THEN

                            IF (T3D(K).GT.270.16) THEN
                                FMULT = 0.
                            ELSE IF (T3D(K).LE.270.16.AND.T3D(K).GT.268.16)  THEN
                                FMULT = (270.16-T3D(K))/2.
                            ELSE IF (T3D(K).GE.265.16.AND.T3D(K).LE.268.16)   THEN
                                FMULT = (T3D(K)-265.16)/3.
                            ELSE IF (T3D(K).LT.265.16) THEN
                                FMULT = 0.
                            END IF

                            ! 1000 IS TO CONVERT FROM KG TO G

                            ! SPLINTERING FROM DROPLETS ACCRETED ONTO SNOW

                            IF (PSACWS(K).GT.0.) THEN
                                NMULTS(K) = 35.E4*PSACWS(K)*FMULT*1000.
                                QMULTS(K) = NMULTS(K)*MMULT

                                ! CONSTRAIN SO THAT TRANSFER OF MASS FROM SNOW TO ICE CANNOT BE MORE MASS
                                ! THAN WAS RIMED ONTO SNOW

                                QMULTS(K) = MIN(QMULTS(K),PSACWS(K))
                                PSACWS(K) = PSACWS(K)-QMULTS(K)

                            END IF

                            ! RIMING AND SPLINTERING FROM ACCRETED RAINDROPS

                            IF (PRACS(K).GT.0.) THEN
                                NMULTR(K) = 35.E4*PRACS(K)*FMULT*1000.
                                QMULTR(K) = NMULTR(K)*MMULT

                                ! CONSTRAIN SO THAT TRANSFER OF MASS FROM SNOW TO ICE CANNOT BE MORE MASS
                                ! THAN WAS RIMED ONTO SNOW

                                QMULTR(K) = MIN(QMULTR(K),PRACS(K))

                                PRACS(K) = PRACS(K)-QMULTR(K)

                            END IF

                        END IF
                    END IF
                END IF
            END IF

            !.......................................................................
            ! RIME-SPLINTERING - GRAUPEL
            ! HALLET-MOSSOP (1974)
            ! NUMBER OF SPLINTERS FORMED IS BASED ON MASS OF RIMED WATER

            ! DUM1 = MASS OF INDIVIDUAL SPLINTERS

            ! HM ADD THRESHOLD SNOW MIXING RATIO FOR RIME-SPLINTERING
            ! TO LIMIT RIME-SPLINTERING IN STRATIFORM CLOUDS

            !         IF (IHAIL.EQ.0) THEN
            ! v1.4
            IF (QG3D(K).GE.0.1E-3) THEN
                IF (QC3D(K).GE.0.5E-3.OR.QR3D(K).GE.0.1E-3) THEN
                    IF (PSACWG(K).GT.0..OR.PRACG(K).GT.0.) THEN
                        IF (T3D(K).LT.270.16 .AND. T3D(K).GT.265.16) THEN

                            IF (T3D(K).GT.270.16) THEN
                                FMULT = 0.
                            ELSE IF (T3D(K).LE.270.16.AND.T3D(K).GT.268.16)  THEN
                                FMULT = (270.16-T3D(K))/2.
                            ELSE IF (T3D(K).GE.265.16.AND.T3D(K).LE.268.16)   THEN
                                FMULT = (T3D(K)-265.16)/3.
                            ELSE IF (T3D(K).LT.265.16) THEN
                                FMULT = 0.
                            END IF

                            ! 1000 IS TO CONVERT FROM KG TO G

                            ! SPLINTERING FROM DROPLETS ACCRETED ONTO GRAUPEL

                            IF (PSACWG(K).GT.0.) THEN
                                NMULTG(K) = 35.E4*PSACWG(K)*FMULT*1000.
                                QMULTG(K) = NMULTG(K)*MMULT

                                ! CONSTRAIN SO THAT TRANSFER OF MASS FROM GRAUPEL TO ICE CANNOT BE MORE MASS
                                ! THAN WAS RIMED ONTO GRAUPEL

                                QMULTG(K) = MIN(QMULTG(K),PSACWG(K))
                                PSACWG(K) = PSACWG(K)-QMULTG(K)

                            END IF

                            ! RIMING AND SPLINTERING FROM ACCRETED RAINDROPS

                            IF (PRACG(K).GT.0.) THEN
                                NMULTRG(K) = 35.E4*PRACG(K)*FMULT*1000.
                                QMULTRG(K) = NMULTRG(K)*MMULT

                                ! CONSTRAIN SO THAT TRANSFER OF MASS FROM GRAUPEL TO ICE CANNOT BE MORE MASS
                                ! THAN WAS RIMED ONTO GRAUPEL

                                QMULTRG(K) = MIN(QMULTRG(K),PRACG(K))
                                PRACG(K) = PRACG(K)-QMULTRG(K)

                            END IF
                        END IF
                    END IF
                END IF
            END IF
            !         END IF

            !........................................................................
            ! CONVERSION OF RIMED CLOUD WATER ONTO SNOW TO GRAUPEL/HAIL

            !           IF (IHAIL.EQ.0) THEN
            IF (PSACWS(K).GT.0.) THEN
                ! ONLY ALLOW CONVERSION IF QNI > 0.1 AND QC > 0.5 G/KG FOLLOWING RUTLEDGE AND HOBBS (1984)
                IF (QNI3D(K).GE.0.1E-3.AND.QC3D(K).GE.0.5E-3) THEN

                    ! PORTION OF RIMING CONVERTED TO GRAUPEL (REISNER ET AL. 1998, ORIGINALLY IS1991)
                    PGSACW(K) = MIN(PSACWS(K),CONS17*DT*N0S(K)*QC3D(K)*QC3D(K)* &
                            ASN(K)*ASN(K)/ &
                            (RHO(K)*LAMS(K)**(2.*BS+2.)))

                    ! MIX RAT CONVERTED INTO GRAUPEL AS EMBRYO (REISNER ET AL. 1998, ORIG M1990)
                    DUM = MAX(RHOSN/(RHOG-RHOSN)*PGSACW(K),0.)

                    ! NUMBER CONCENTRAITON OF EMBRYO GRAUPEL FROM RIMING OF SNOW
                    NSCNG(K) = DUM/MG0*RHO(K)
                    ! LIMIT MAX NUMBER CONVERTED TO SNOW NUMBER
                    NSCNG(K) = MIN(NSCNG(K),NS3D(K)/DT)

                    ! PORTION OF RIMING LEFT FOR SNOW
                    PSACWS(K) = PSACWS(K) - PGSACW(K)
                END IF
            END IF

            ! CONVERSION OF RIMED RAINWATER ONTO SNOW CONVERTED TO GRAUPEL

            IF (PRACS(K).GT.0.) THEN
                ! ONLY ALLOW CONVERSION IF QNI > 0.1 AND QR > 0.1 G/KG FOLLOWING RUTLEDGE AND HOBBS (1984)
                IF (QNI3D(K).GE.0.1E-3.AND.QR3D(K).GE.0.1E-3) THEN
                    ! PORTION OF COLLECTED RAINWATER CONVERTED TO GRAUPEL (REISNER ET AL. 1998)
                    DUM = CONS18*(4./LAMS(K))**3*(4./LAMS(K))**3 &
                            /(CONS18*(4./LAMS(K))**3*(4./LAMS(K))**3+ &
                                    CONS19*(4./LAMR(K))**3*(4./LAMR(K))**3)
                    DUM=MIN(DUM,1.)
                    DUM=MAX(DUM,0.)
                    PGRACS(K) = (1.-DUM)*PRACS(K)
                    NGRACS(K) = (1.-DUM)*NPRACS(K)
                    ! LIMIT MAX NUMBER CONVERTED TO MIN OF EITHER RAIN OR SNOW NUMBER CONCENTRATION
                    NGRACS(K) = MIN(NGRACS(K),NR3D(K)/DT)
                    NGRACS(K) = MIN(NGRACS(K),NS3D(K)/DT)

                    ! AMOUNT LEFT FOR SNOW PRODUCTION
                    PRACS(K) = PRACS(K) - PGRACS(K)
                    NPRACS(K) = NPRACS(K) - NGRACS(K)
                    ! CONVERSION TO GRAUPEL DUE TO COLLECTION OF SNOW BY RAIN
                    PSACR(K)=PSACR(K)*(1.-DUM)
                END IF
            END IF
            !           END IF

            !.......................................................................
            ! FREEZING OF RAIN DROPS
            ! FREEZING ALLOWED BELOW -4 C

            IF (T3D(K).LT.269.15.AND.QR3D(K).GE.QSMALL) THEN

                ! IMMERSION FREEZING (BIGG 1953)
                !            MNUCCR(K) = CONS20*NR3D(K)*EXP(AIMM*(273.15-T3D(K)))/LAMR(K)**3 &
                !                 /LAMR(K)**3

                !            NNUCCR(K) = PI*NR3D(K)*BIMM*EXP(AIMM*(273.15-T3D(K)))/LAMR(K)**3

                ! hm fix 7/15/13 for consistency w/ original formula
                MNUCCR(K) = CONS20*NR3D(K)*(EXP(AIMM*(273.15-T3D(K)))-1.)/LAMR(K)**3 &
                        /LAMR(K)**3

                NNUCCR(K) = PI*NR3D(K)*BIMM*(EXP(AIMM*(273.15-T3D(K)))-1.)/LAMR(K)**3

                ! PREVENT DIVERGENCE BETWEEN MIXING RATIO AND NUMBER CONC
                NNUCCR(K) = MIN(NNUCCR(K),NR3D(K)/DT)

            END IF

            !.......................................................................
            ! ACCRETION OF CLOUD LIQUID WATER BY RAIN
            ! CONTINUOUS COLLECTION EQUATION WITH
            ! GRAVITATIONAL COLLECTION KERNEL, DROPLET FALL SPEED NEGLECTED

            IF (QR3D(K).GE.1.E-8 .AND. QC3D(K).GE.1.E-8) THEN

                ! 12/13/06 HM ADD, REPLACE WITH NEWER FORMULA FROM
                ! KHAIROUTDINOV AND KOGAN 2000, MWR

                DUM=(QC3D(K)*QR3D(K))
                PRA(K) = 67.*(DUM)**1.15
                NPRA(K) = PRA(K)/(QC3D(K)/NC3D(K))

            END IF
            !.......................................................................
            ! SELF-COLLECTION OF RAIN DROPS
            ! FROM BEHENG(1994)
            ! FROM NUMERICAL SIMULATION OF THE STOCHASTIC COLLECTION EQUATION
            ! AS DESCRINED ABOVE FOR AUTOCONVERSION

            IF (QR3D(K).GE.1.E-8) THEN
                ! include breakup add 10/09/09
                dum1=300.e-6
                if (1./lamr(k).lt.dum1) then
                    dum=1.
                else if (1./lamr(k).ge.dum1) then
                    dum=2.-exp(2300.*(1./lamr(k)-dum1))
                end if
                !            NRAGG(K) = -8.*NR3D(K)*QR3D(K)*RHO(K)
                NRAGG(K) = -5.78*dum*NR3D(K)*QR3D(K)*RHO(K)
            END IF

            !.......................................................................
            ! AUTOCONVERSION OF CLOUD ICE TO SNOW
            ! FOLLOWING HARRINGTON ET AL. (1995) WITH MODIFICATION
            ! HERE IT IS ASSUMED THAT AUTOCONVERSION CAN ONLY OCCUR WHEN THE
            ! ICE IS GROWING, I.E. IN CONDITIONS OF ICE SUPERSATURATION

            IF (QI3D(K).GE.1.E-8 .AND.QVQVSI(K).GE.1.) THEN

                !           COFFI = 2./LAMI(K)
                !           IF (COFFI.GE.DCS) THEN
                NPRCI(K) = CONS21*(QV3D(K)-QVI(K))*RHO(K)                         &
                        *N0I(K)*EXP(-LAMI(K)*DCS)*DV(K)/ABI(K)
                PRCI(K) = CONS22*NPRCI(K)
                NPRCI(K) = MIN(NPRCI(K),NI3D(K)/DT)

                !           END IF
            END IF

            !.......................................................................
            ! ACCRETION OF CLOUD ICE BY SNOW
            ! FOR THIS CALCULATION, IT IS ASSUMED THAT THE VS >> VI
            ! AND DS >> DI FOR CONTINUOUS COLLECTION

            IF (QNI3D(K).GE.1.E-8 .AND. QI3D(K).GE.QSMALL) THEN
                PRAI(K) = CONS23*ASN(K)*QI3D(K)*RHO(K)*N0S(K)/     &
                        LAMS(K)**(BS+3.)
                NPRAI(K) = CONS23*ASN(K)*NI3D(K)*                                       &
                        RHO(K)*N0S(K)/                                 &
                        LAMS(K)**(BS+3.)
                NPRAI(K)=MIN(NPRAI(K),NI3D(K)/DT)
                if (ipelec>0) then
                    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
                    ! Collision rate of snow-ice, morr (CSACI0_MORR)
                    CSACI0_MORR(K) = NPRAI(K)/EII
                    !                CSACI0_MORR(K) = PI/4.* GAMMA(BS+3.)*ASN(K)/LAMS(K)**(BS+2.) ! Collision kernel ***CONS23=PI/4.*EII*GAMMA(BS+3.)
                    !                CSACI0_MORR(K) = NI3D(K)*RHO(K)*N0S(K)/LAMS(K) ! Everything else
                    !                print*, 'CSACI0_MORR=',CSACI0_MORR(k), 'i,j,k=', i,j,k
                    !                print*, 'ni3d=',ni3d(k),'rho=',rho(k),'n0s=',n0s(k),'lams=',lams(k)
                end if
            END IF

            !.......................................................................
            ! HM, ADD 12/13/06, COLLISION OF RAIN AND ICE TO PRODUCE SNOW OR GRAUPEL
            ! FOLLOWS REISNER ET AL. 1998
            ! ASSUMED FALLSPEED AND SIZE OF ICE CRYSTAL << THAN FOR RAIN

            IF (QR3D(K).GE.1.E-8.AND.QI3D(K).GE.1.E-8.AND.T3D(K).LE.273.15) THEN

                ! ALLOW GRAUPEL FORMATION FROM RAIN-ICE COLLISIONS ONLY IF RAIN MIXING RATIO > 0.1 G/KG,
                ! OTHERWISE ADD TO SNOW

                IF (QR3D(K).GE.0.1E-3) THEN
                    NIACR(K)=CONS24*NI3D(K)*N0RR(K)*ARN(K) &
                            /LAMR(K)**(BR+3.)*RHO(K)
                    PIACR(K)=CONS25*NI3D(K)*N0RR(K)*ARN(K) &
                            /LAMR(K)**(BR+3.)/LAMR(K)**3*RHO(K)
                    PRACI(K)=CONS24*QI3D(K)*N0RR(K)*ARN(K)/ &
                            LAMR(K)**(BR+3.)*RHO(K)
                    NIACR(K)=MIN(NIACR(K),NR3D(K)/DT)
                    NIACR(K)=MIN(NIACR(K),NI3D(K)/DT)
                ELSE
                    NIACRS(K)=CONS24*NI3D(K)*N0RR(K)*ARN(K) &
                            /LAMR(K)**(BR+3.)*RHO(K)
                    PIACRS(K)=CONS25*NI3D(K)*N0RR(K)*ARN(K) &
                            /LAMR(K)**(BR+3.)/LAMR(K)**3*RHO(K)
                    PRACIS(K)=CONS24*QI3D(K)*N0RR(K)*ARN(K)/ &
                            LAMR(K)**(BR+3.)*RHO(K)
                    NIACRS(K)=MIN(NIACRS(K),NR3D(K)/DT)
                    NIACRS(K)=MIN(NIACRS(K),NI3D(K)/DT)
                END IF
            END IF

            !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            ! NUCLEATION OF CLOUD ICE FROM HOMOGENEOUS AND HETEROGENEOUS FREEZING ON AEROSOL

            IF (INUC.EQ.0) THEN

                ! add threshold according to Greg Thomspon

                if ((QVQVS(K).GE.0.999.and.T3D(K).le.265.15).or. &
                        QVQVSI(K).ge.1.08) then

                    ! hm, modify dec. 5, 2006, replace with cooper curve
                    kc2 = 0.005*exp(0.304*(273.15-T3D(K)))*1000. ! convert from L-1 to m-3
                    ! limit to 500 L-1
                    kc2 = min(kc2,500.e3)
                    kc2=MAX(kc2/rho(k),0.)  ! convert to kg-1

                    IF (KC2.GT.NI3D(K)+NS3D(K)+NG3D(K)) THEN
                        NNUCCD(K) = (KC2-NI3D(K)-NS3D(K)-NG3D(K))/DT
                        MNUCCD(K) = NNUCCD(K)*MI0
                    END IF

                END IF

            ELSE IF (INUC.EQ.1) THEN

                IF (T3D(K).LT.273.15.AND.QVQVSI(K).GT.1.) THEN

                    KC2 = 0.16*1000./RHO(K)  ! CONVERT FROM L-1 TO KG-1
                    IF (KC2.GT.NI3D(K)+NS3D(K)+NG3D(K)) THEN
                        NNUCCD(K) = (KC2-NI3D(K)-NS3D(K)-NG3D(K))/DT
                        MNUCCD(K) = NNUCCD(K)*MI0
                    END IF
                END IF

            END IF

            !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

            101      CONTINUE

            !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            ! CALCULATE EVAP/SUB/DEP TERMS FOR QI,QNI,QR

            ! NO VENTILATION FOR CLOUD ICE

            IF (QI3D(K).GE.QSMALL) THEN

                EPSI = 2.*PI*N0I(K)*RHO(K)*DV(K)/(LAMI(K)*LAMI(K))

            ELSE
                EPSI = 0.
            END IF

            IF (QNI3D(K).GE.QSMALL) THEN
                EPSS = 2.*PI*N0S(K)*RHO(K)*DV(K)*                            &
                        (F1S/(LAMS(K)*LAMS(K))+                       &
                                F2S*(ASN(K)*RHO(K)/MU(K))**0.5*                      &
                                        SC(K)**(1./3.)*CONS10/                   &
                                        (LAMS(K)**CONS35))
            ELSE
                EPSS = 0.
            END IF

            IF (QG3D(K).GE.QSMALL) THEN
                EPSG = 2.*PI*N0G(K)*RHO(K)*DV(K)*                                &
                        (F1S/(LAMG(K)*LAMG(K))+                               &
                                F2S*(AGN(K)*RHO(K)/MU(K))**0.5*                      &
                                        SC(K)**(1./3.)*CONS11/                   &
                                        (LAMG(K)**CONS36))


            ELSE
                EPSG = 0.
            END IF

            IF (QR3D(K).GE.QSMALL) THEN
                EPSR = 2.*PI*N0RR(K)*RHO(K)*DV(K)*                           &
                        (F1R/(LAMR(K)*LAMR(K))+                       &
                                F2R*(ARN(K)*RHO(K)/MU(K))**0.5*                      &
                                        SC(K)**(1./3.)*CONS9/                   &
                                        (LAMR(K)**CONS34))
            ELSE
                EPSR = 0.
            END IF

            ! ONLY INCLUDE REGION OF ICE SIZE DIST < DCS
            ! DUM IS FRACTION OF D*N(D) < DCS

            ! LOGIC BELOW FOLLOWS THAT OF HARRINGTON ET AL. 1995 (JAS)
            IF (QI3D(K).GE.QSMALL) THEN
                DUM=(1.-EXP(-LAMI(K)*DCS)*(1.+LAMI(K)*DCS))
                PRD(K) = EPSI*(QV3D(K)-QVI(K))/ABI(K)*DUM
            ELSE
                DUM=0.
            END IF
            ! ADD DEPOSITION IN TAIL OF ICE SIZE DIST TO SNOW IF SNOW IS PRESENT
            IF (QNI3D(K).GE.QSMALL) THEN
                PRDS(K) = EPSS*(QV3D(K)-QVI(K))/ABI(K)+ &
                        EPSI*(QV3D(K)-QVI(K))/ABI(K)*(1.-DUM)
                ! OTHERWISE ADD TO CLOUD ICE
            ELSE
                PRD(K) = PRD(K)+EPSI*(QV3D(K)-QVI(K))/ABI(K)*(1.-DUM)
            END IF
            ! VAPOR DPEOSITION ON GRAUPEL
            PRDG(K) = EPSG*(QV3D(K)-QVI(K))/ABI(K)

            ! NO CONDENSATION ONTO RAIN, ONLY EVAP

            IF (QV3D(K).LT.QVS(K)) THEN
                PRE(K) = EPSR*(QV3D(K)-QVS(K))/AB(K)
                PRE(K) = MIN(PRE(K),0.)
            ELSE
                PRE(K) = 0.
            END IF

            ! MAKE SURE NOT PUSHED INTO ICE SUPERSAT/SUBSAT
            ! FORMULA FROM REISNER 2 SCHEME

            DUM = (QV3D(K)-QVI(K))/DT

            FUDGEF = 0.9999
            SUM_DEP = PRD(K)+PRDS(K)+MNUCCD(K)+PRDG(K)

            IF( (DUM.GT.0. .AND. SUM_DEP.GT.DUM*FUDGEF) .OR.                      &
                    (DUM.LT.0. .AND. SUM_DEP.LT.DUM*FUDGEF) ) THEN
                MNUCCD(K) = FUDGEF*MNUCCD(K)*DUM/SUM_DEP
                PRD(K) = FUDGEF*PRD(K)*DUM/SUM_DEP
                PRDS(K) = FUDGEF*PRDS(K)*DUM/SUM_DEP
                PRDG(K) = FUDGEF*PRDG(K)*DUM/SUM_DEP
            ENDIF

            ! IF CLOUD ICE/SNOW/GRAUPEL VAP DEPOSITION IS NEG, THEN ASSIGN TO SUBLIMATION PROCESSES

            IF (PRD(K).LT.0.) THEN
                EPRD(K)=PRD(K)
                PRD(K)=0.
            END IF
            IF (PRDS(K).LT.0.) THEN
                EPRDS(K)=PRDS(K)
                PRDS(K)=0.
            END IF
            IF (PRDG(K).LT.0.) THEN
                EPRDG(K)=PRDG(K)
                PRDG(K)=0.
            END IF

            !.......................................................................
            !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

            ! CONSERVATION OF WATER
            ! THIS IS ADOPTED LOOSELY FROM MM5 RESINER CODE. HOWEVER, HERE WE
            ! ONLY ADJUST PROCESSES THAT ARE NEGATIVE, RATHER THAN ALL PROCESSES.

            ! IF MIXING RATIOS LESS THAN QSMALL, THEN NO DEPLETION OF WATER
            ! THROUGH MICROPHYSICAL PROCESSES, SKIP CONSERVATION

            ! NOTE: CONSERVATION CHECK NOT APPLIED TO NUMBER CONCENTRATION SPECIES. ADDITIONAL CATCH
            ! BELOW WILL PREVENT NEGATIVE NUMBER CONCENTRATION
            ! FOR EACH MICROPHYSICAL PROCESS WHICH PROVIDES A SOURCE FOR NUMBER, THERE IS A CHECK
            ! TO MAKE SURE THAT CAN'T EXCEED TOTAL NUMBER OF DEPLETED SPECIES WITH THE TIME
            ! STEP

            !****SENSITIVITY - NO ICE

            IF (ILIQ.EQ.1) THEN
                MNUCCC(K)=0.
                NNUCCC(K)=0.
                MNUCCR(K)=0.
                NNUCCR(K)=0.
                MNUCCD(K)=0.
                NNUCCD(K)=0.
            END IF

            ! ****SENSITIVITY - NO GRAUPEL
            IF (IGRAUP.EQ.1) THEN
                PRACG(K) = 0.
                PSACR(K) = 0.
                PSACWG(K) = 0.
                PRDG(K) = 0.
                EPRDG(K) = 0.
                EVPMG(K) = 0.
                PGMLT(K) = 0.
                NPRACG(K) = 0.
                NPSACWG(K) = 0.
                NSCNG(K) = 0.
                NGRACS(K) = 0.
                NSUBG(K) = 0.
                NGMLTG(K) = 0.
                NGMLTR(K) = 0.
                ! fix 053011
                PIACRS(K)=PIACRS(K)+PIACR(K)
                PIACR(K) = 0.
                ! fix 070713
                PRACIS(K)=PRACIS(K)+PRACI(K)
                PRACI(K) = 0.
                PSACWS(K)=PSACWS(K)+PGSACW(K)
                PGSACW(K) = 0.
                PRACS(K)=PRACS(K)+PGRACS(K)
                PGRACS(K) = 0.
            END IF

            ! CONSERVATION OF QC

            DUM = (PRC(K)+PRA(K)+MNUCCC(K)+PSACWS(K)+PSACWI(K)+QMULTS(K)+PSACWG(K)+PGSACW(K)+QMULTG(K))*DT

            IF (DUM.GT.QC3D(K).AND.QC3D(K).GE.QSMALL) THEN
                RATIO = QC3D(K)/DUM

                PRC(K) = PRC(K)*RATIO
                PRA(K) = PRA(K)*RATIO
                MNUCCC(K) = MNUCCC(K)*RATIO
                PSACWS(K) = PSACWS(K)*RATIO
                PSACWI(K) = PSACWI(K)*RATIO
                QMULTS(K) = QMULTS(K)*RATIO
                QMULTG(K) = QMULTG(K)*RATIO
                PSACWG(K) = PSACWG(K)*RATIO
                PGSACW(K) = PGSACW(K)*RATIO
            END IF

            ! CONSERVATION OF QI

            DUM = (-PRD(K)-MNUCCC(K)+PRCI(K)+PRAI(K)-QMULTS(K)-QMULTG(K)-QMULTR(K)-QMULTRG(K) &
                    -MNUCCD(K)+PRACI(K)+PRACIS(K)-EPRD(K)-PSACWI(K))*DT

            IF (DUM.GT.QI3D(K).AND.QI3D(K).GE.QSMALL) THEN

                RATIO = (QI3D(K)/DT+PRD(K)+MNUCCC(K)+QMULTS(K)+QMULTG(K)+QMULTR(K)+QMULTRG(K)+ &
                        MNUCCD(K)+PSACWI(K))/ &
                        (PRCI(K)+PRAI(K)+PRACI(K)+PRACIS(K)-EPRD(K))

                PRCI(K) = PRCI(K)*RATIO
                PRAI(K) = PRAI(K)*RATIO
                PRACI(K) = PRACI(K)*RATIO
                PRACIS(K) = PRACIS(K)*RATIO
                EPRD(K) = EPRD(K)*RATIO

            END IF

            ! CONSERVATION OF QR

            DUM=((PRACS(K)-PRE(K))+(QMULTR(K)+QMULTRG(K)-PRC(K))+(MNUCCR(K)-PRA(K))+ &
                    PIACR(K)+PIACRS(K)+PGRACS(K)+PRACG(K))*DT

            IF (DUM.GT.QR3D(K).AND.QR3D(K).GE.QSMALL) THEN

                RATIO = (QR3D(K)/DT+PRC(K)+PRA(K))/ &
                        (-PRE(K)+QMULTR(K)+QMULTRG(K)+PRACS(K)+MNUCCR(K)+PIACR(K)+PIACRS(K)+PGRACS(K)+PRACG(K))

                PRE(K) = PRE(K)*RATIO
                PRACS(K) = PRACS(K)*RATIO
                QMULTR(K) = QMULTR(K)*RATIO
                QMULTRG(K) = QMULTRG(K)*RATIO
                MNUCCR(K) = MNUCCR(K)*RATIO
                PIACR(K) = PIACR(K)*RATIO
                PIACRS(K) = PIACRS(K)*RATIO
                PGRACS(K) = PGRACS(K)*RATIO
                PRACG(K) = PRACG(K)*RATIO

            END IF

            ! CONSERVATION OF QNI
            ! CONSERVATION FOR GRAUPEL SCHEME

            IF (IGRAUP.EQ.0) THEN

                DUM = (-PRDS(K)-PSACWS(K)-PRAI(K)-PRCI(K)-PRACS(K)-EPRDS(K)+PSACR(K)-PIACRS(K)-PRACIS(K))*DT

                IF (DUM.GT.QNI3D(K).AND.QNI3D(K).GE.QSMALL) THEN

                    RATIO = (QNI3D(K)/DT+PRDS(K)+PSACWS(K)+PRAI(K)+PRCI(K)+PRACS(K)+PIACRS(K)+PRACIS(K))/(-EPRDS(K)+PSACR(K))

                    EPRDS(K) = EPRDS(K)*RATIO
                    PSACR(K) = PSACR(K)*RATIO

                END IF

                ! FOR NO GRAUPEL, NEED TO INCLUDE FREEZING OF RAIN FOR SNOW
            ELSE IF (IGRAUP.EQ.1) THEN

                DUM = (-PRDS(K)-PSACWS(K)-PRAI(K)-PRCI(K)-PRACS(K)-EPRDS(K)+PSACR(K)-PIACRS(K)-PRACIS(K)-MNUCCR(K))*DT

                IF (DUM.GT.QNI3D(K).AND.QNI3D(K).GE.QSMALL) THEN

                    RATIO = (QNI3D(K)/DT+PRDS(K)+PSACWS(K)+PRAI(K)+PRCI(K)+PRACS(K)+PIACRS(K)+PRACIS(K)+MNUCCR(K))/(-EPRDS(K)+PSACR(K))

                    EPRDS(K) = EPRDS(K)*RATIO
                    PSACR(K) = PSACR(K)*RATIO

                END IF

            END IF

            ! CONSERVATION OF QG

            DUM = (-PSACWG(K)-PRACG(K)-PGSACW(K)-PGRACS(K)-PRDG(K)-MNUCCR(K)-EPRDG(K)-PIACR(K)-PRACI(K)-PSACR(K))*DT

            IF (DUM.GT.QG3D(K).AND.QG3D(K).GE.QSMALL) THEN

                RATIO = (QG3D(K)/DT+PSACWG(K)+PRACG(K)+PGSACW(K)+PGRACS(K)+PRDG(K)+MNUCCR(K)+PSACR(K)+&
                        PIACR(K)+PRACI(K))/(-EPRDG(K))

                EPRDG(K) = EPRDG(K)*RATIO

            END IF

            ! TENDENCIES

            QV3DTEN(K) = QV3DTEN(K)+(-PRE(K)-PRD(K)-PRDS(K)-MNUCCD(K)-EPRD(K)-EPRDS(K)-PRDG(K)-EPRDG(K))

            ! BUG FIX HM, 3/1/11, INCLUDE PIACR AND PIACRS
            T3DTEN(K) = T3DTEN(K)+(PRE(K)                                 &
                    *XXLV(K)+(PRD(K)+PRDS(K)+                            &
                    MNUCCD(K)+EPRD(K)+EPRDS(K)+PRDG(K)+EPRDG(K))*XXLS(K)+         &
                    (PSACWS(K)+PSACWI(K)+MNUCCC(K)+MNUCCR(K)+                      &
                            QMULTS(K)+QMULTG(K)+QMULTR(K)+QMULTRG(K)+PRACS(K) &
                            +PSACWG(K)+PRACG(K)+PGSACW(K)+PGRACS(K)+PIACR(K)+PIACRS(K))*XLF(K))/CPM(K)

            QC3DTEN(K) = QC3DTEN(K)+                                      &
                    (-PRA(K)-PRC(K)-MNUCCC(K)+PCC(K)-                  &
                            PSACWS(K)-PSACWI(K)-QMULTS(K)-QMULTG(K)-PSACWG(K)-PGSACW(K))
            QI3DTEN(K) = QI3DTEN(K)+                                      &
                    (PRD(K)+EPRD(K)+PSACWI(K)+MNUCCC(K)-PRCI(K)-                                 &
                            PRAI(K)+QMULTS(K)+QMULTG(K)+QMULTR(K)+QMULTRG(K)+MNUCCD(K)-PRACI(K)-PRACIS(K))
            QR3DTEN(K) = QR3DTEN(K)+                                      &
                    (PRE(K)+PRA(K)+PRC(K)-PRACS(K)-MNUCCR(K)-QMULTR(K)-QMULTRG(K) &
                            -PIACR(K)-PIACRS(K)-PRACG(K)-PGRACS(K))

            IF (IGRAUP.EQ.0) THEN

                QNI3DTEN(K) = QNI3DTEN(K)+                                    &
                        (PRAI(K)+PSACWS(K)+PRDS(K)+PRACS(K)+PRCI(K)+EPRDS(K)-PSACR(K)+PIACRS(K)+PRACIS(K))
                NS3DTEN(K) = NS3DTEN(K)+(NSAGG(K)+NPRCI(K)-NSCNG(K)-NGRACS(K)+NIACRS(K))
                QG3DTEN(K) = QG3DTEN(K)+(PRACG(K)+PSACWG(K)+PGSACW(K)+PGRACS(K)+ &
                        PRDG(K)+EPRDG(K)+MNUCCR(K)+PIACR(K)+PRACI(K)+PSACR(K))
                NG3DTEN(K) = NG3DTEN(K)+(NSCNG(K)+NGRACS(K)+NNUCCR(K)+NIACR(K))

                ! FOR NO GRAUPEL, NEED TO INCLUDE FREEZING OF RAIN FOR SNOW
            ELSE IF (IGRAUP.EQ.1) THEN

                QNI3DTEN(K) = QNI3DTEN(K)+                                    &
                        (PRAI(K)+PSACWS(K)+PRDS(K)+PRACS(K)+PRCI(K)+EPRDS(K)-PSACR(K)+PIACRS(K)+PRACIS(K)+MNUCCR(K))
                NS3DTEN(K) = NS3DTEN(K)+(NSAGG(K)+NPRCI(K)-NSCNG(K)-NGRACS(K)+NIACRS(K)+NNUCCR(K))

            END IF

            NC3DTEN(K) = NC3DTEN(K)+(-NNUCCC(K)-NPSACWS(K)                &
                    -NPRA(K)-NPRC(K)-NPSACWI(K)-NPSACWG(K))

            NI3DTEN(K) = NI3DTEN(K)+                                      &
                    (NNUCCC(K)-NPRCI(K)-NPRAI(K)+NMULTS(K)+NMULTG(K)+NMULTR(K)+NMULTRG(K)+ &
                            NNUCCD(K)-NIACR(K)-NIACRS(K))

            NR3DTEN(K) = NR3DTEN(K)+(NPRC1(K)-NPRACS(K)-NNUCCR(K)      &
                    +NRAGG(K)-NIACR(K)-NIACRS(K)-NPRACG(K)-NGRACS(K))

            ! HM ADD, WRF-CHEM, ADD TENDENCIES FOR C2PREC

            C2PREC(K) = PRA(K)+PRC(K)+PSACWS(K)+QMULTS(K)+QMULTG(K)+PSACWG(K)+ &
                    PGSACW(K)+MNUCCC(K)+PSACWI(K)
            !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            ! NOW CALCULATE SATURATION ADJUSTMENT TO CONDENSE EXTRA VAPOR ABOVE
            ! WATER SATURATION

            DUMT = T3D(K)+DT*T3DTEN(K)
            DUMQV = QV3D(K)+DT*QV3DTEN(K)
            ! hm, add fix for low pressure, 5/12/10
            dum=min(0.99*pres(k),POLYSVP(DUMT,0))
            DUMQSS = EP_2*dum/(PRES(K)-dum)
            DUMQC = QC3D(K)+DT*QC3DTEN(K)
            DUMQC = MAX(DUMQC,0.)

            ! SATURATION ADJUSTMENT FOR LIQUID

            DUMS = DUMQV-DUMQSS
            PCC(K) = DUMS/(1.+XXLV(K)**2*DUMQSS/(CPM(K)*RV*DUMT**2))/DT
            IF (PCC(K)*DT+DUMQC.LT.0.) THEN
                PCC(K) = -DUMQC/DT
            END IF

            QV3DTEN(K) = QV3DTEN(K)-PCC(K)
            T3DTEN(K) = T3DTEN(K)+PCC(K)*XXLV(K)/CPM(K)
            QC3DTEN(K) = QC3DTEN(K)+PCC(K)

            !.......................................................................
            ! ACTIVATION OF CLOUD DROPLETS
            ! ACTIVATION OF DROPLET CURRENTLY NOT CALCULATED
            ! DROPLET CONCENTRATION IS SPECIFIED !!!!!

            !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            ! SUBLIMATE, MELT, OR EVAPORATE NUMBER CONCENTRATION
            ! THIS FORMULATION ASSUMES 1:1 RATIO BETWEEN MASS LOSS AND
            ! LOSS OF NUMBER CONCENTRATION

            !     IF (PCC(K).LT.0.) THEN
            !        DUM = PCC(K)*DT/QC3D(K)
            !           DUM = MAX(-1.,DUM)
            !        NSUBC(K) = DUM*NC3D(K)/DT
            !     END IF

            IF (EPRD(K).LT.0.) THEN
                DUM = EPRD(K)*DT/QI3D(K)
                DUM = MAX(-1.,DUM)
                NSUBI(K) = DUM*NI3D(K)/DT
            END IF
            IF (EPRDS(K).LT.0.) THEN
                DUM = EPRDS(K)*DT/QNI3D(K)
                DUM = MAX(-1.,DUM)
                NSUBS(K) = DUM*NS3D(K)/DT
            END IF
            IF (PRE(K).LT.0.) THEN
                DUM = PRE(K)*DT/QR3D(K)
                DUM = MAX(-1.,DUM)
                NSUBR(K) = DUM*NR3D(K)/DT
            END IF
            IF (EPRDG(K).LT.0.) THEN
                DUM = EPRDG(K)*DT/QG3D(K)
                DUM = MAX(-1.,DUM)
                NSUBG(K) = DUM*NG3D(K)/DT
            END IF

            !        nsubr(k)=0.
            !        nsubs(k)=0.
            !        nsubg(k)=0.

            ! UPDATE TENDENCIES

            !        NC3DTEN(K) = NC3DTEN(K)+NSUBC(K)
            NI3DTEN(K) = NI3DTEN(K)+NSUBI(K)
            NS3DTEN(K) = NS3DTEN(K)+NSUBS(K)
            NG3DTEN(K) = NG3DTEN(K)+NSUBG(K)
            NR3DTEN(K) = NR3DTEN(K)+NSUBR(K)

#if (WRF_CHEM == 1)
            evapprod(k) = - PRE(K) - EPRDS(K) - EPRDG(K)
            rainprod(k) = PRA(K) + PRC(K) + PSACWS(K) + PSACWG(K) + PGSACW(K) &
                    + PRAI(K) + PRCI(K) + PRACI(K) + PRACIS(K)  &
                    + PRDS(K) + PRDG(K)
#endif

        END IF !!!!!! TEMPERATURE

        ! SWITCH LTRUE TO 1, SINCE HYDROMETEORS ARE PRESENT
        LTRUE = 1

        200     CONTINUE

    END DO

    ! INITIALIZE PRECIP AND SNOW RATES
    PRECRT = 0.
    SNOWRT = 0.
    ! hm added 7/13/13
    SNOWPRT = 0.
    GRPLPRT = 0.

    ! IF THERE ARE NO HYDROMETEORS, THEN SKIP TO END OF SUBROUTINE

    IF (LTRUE.EQ.0) GOTO 400

    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    !.......................................................................
    ! CALCULATE SEDIMENATION
    ! THE NUMERICS HERE FOLLOW FROM REISNER ET AL. (1998)
    ! FALLOUT TERMS ARE CALCULATED ON SPLIT TIME STEPS TO ENSURE NUMERICAL
    ! STABILITY, I.E. COURANT# < 1

    !.......................................................................

    NSTEP = 1

    DO K = KTE,KTS,-1

        DUMI(K) = QI3D(K)+QI3DTEN(K)*DT
        DUMQS(K) = QNI3D(K)+QNI3DTEN(K)*DT
        DUMR(K) = QR3D(K)+QR3DTEN(K)*DT
        dumrtmp(k) = dumr(k)
        DUMFNI(K) = NI3D(K)+NI3DTEN(K)*DT
        DUMFNS(K) = NS3D(K)+NS3DTEN(K)*DT
        DUMFNR(K) = NR3D(K)+NR3DTEN(K)*DT
        DUMC(K) = QC3D(K)+QC3DTEN(K)*DT
        DUMFNC(K) = NC3D(K)+NC3DTEN(K)*DT
        DUMG(K) = QG3D(K)+QG3DTEN(K)*DT
        DUMFNG(K) = NG3D(K)+NG3DTEN(K)*DT

        ! SWITCH FOR CONSTANT DROPLET NUMBER
        IF (iinum.EQ.1) THEN
            DUMFNC(K) = NC3D(K)
        END IF

        ! GET DUMMY LAMDA FOR SEDIMENTATION CALCULATIONS

        ! MAKE SURE NUMBER CONCENTRATIONS ARE POSITIVE
        DUMFNI(K) = MAX(0.,DUMFNI(K))
        DUMFNS(K) = MAX(0.,DUMFNS(K))
        DUMFNC(K) = MAX(0.,DUMFNC(K))
        DUMFNR(K) = MAX(0.,DUMFNR(K))
        DUMFNG(K) = MAX(0.,DUMFNG(K))

        !......................................................................
        ! CLOUD ICE

        IF (DUMI(K).GE.QSMALL) THEN
            DLAMI = (CONS12*DUMFNI(K)/DUMI(K))**(1./DI)
            DLAMI=MAX(DLAMI,LAMMINI)
            DLAMI=MIN(DLAMI,LAMMAXI)
        END IF
        !......................................................................
        ! RAIN

        IF (DUMR(K).GE.QSMALL) THEN
            DLAMR = (PI*RHOW*DUMFNR(K)/DUMR(K))**(1./3.)
            DLAMR=MAX(DLAMR,LAMMINR)
            DLAMR=MIN(DLAMR,LAMMAXR)
        END IF
        !......................................................................
        ! CLOUD DROPLETS

        IF (DUMC(K).GE.QSMALL) THEN
            DUM = PRES(K)/(287.15*T3D(K))
            PGAM(K)=0.0005714*(NC3D(K)/1.E6*DUM)+0.2714
            PGAM(K)=1./(PGAM(K)**2)-1.
            PGAM(K)=MAX(PGAM(K),2.)
            PGAM(K)=MIN(PGAM(K),10.)

            DLAMC = (CONS26*DUMFNC(K)*GAMMA(PGAM(K)+4.)/(DUMC(K)*GAMMA(PGAM(K)+1.)))**(1./3.)
            LAMMIN = (PGAM(K)+1.)/60.E-6
            LAMMAX = (PGAM(K)+1.)/1.E-6
            DLAMC=MAX(DLAMC,LAMMIN)
            DLAMC=MIN(DLAMC,LAMMAX)
        END IF
        !......................................................................
        ! SNOW

        IF (DUMQS(K).GE.QSMALL) THEN
            DLAMS = (CONS1*DUMFNS(K)/ DUMQS(K))**(1./DS)
            DLAMS=MAX(DLAMS,LAMMINS)
            DLAMS=MIN(DLAMS,LAMMAXS)
        END IF
        !......................................................................
        ! GRAUPEL

        IF (DUMG(K).GE.QSMALL) THEN
            DLAMG = (CONS2*DUMFNG(K)/ DUMG(K))**(1./DG)
            DLAMG=MAX(DLAMG,LAMMING)
            DLAMG=MIN(DLAMG,LAMMAXG)
        END IF

        !......................................................................
        ! CALCULATE NUMBER-WEIGHTED AND MASS-WEIGHTED TERMINAL FALL SPEEDS

        ! CLOUD WATER

        IF (DUMC(K).GE.QSMALL) THEN
            UNC =  ACN(K)*GAMMA(1.+BC+PGAM(K))/ (DLAMC**BC*GAMMA(PGAM(K)+1.))
            UMC = ACN(K)*GAMMA(4.+BC+PGAM(K))/  (DLAMC**BC*GAMMA(PGAM(K)+4.))
        ELSE
            UMC = 0.
            UNC = 0.
        END IF

        IF (DUMI(K).GE.QSMALL) THEN
            UNI =  AIN(K)*CONS27/DLAMI**BI
            UMI = AIN(K)*CONS28/(DLAMI**BI)
        ELSE
            UMI = 0.
            UNI = 0.
        END IF

        IF (DUMR(K).GE.QSMALL) THEN
            UNR = ARN(K)*CONS6/DLAMR**BR
            UMR = ARN(K)*CONS4/(DLAMR**BR)
        ELSE
            UMR = 0.
            UNR = 0.
        END IF

        IF (DUMQS(K).GE.QSMALL) THEN
            UMS = ASN(K)*CONS3/(DLAMS**BS)
            UNS = ASN(K)*CONS5/DLAMS**BS
        ELSE
            UMS = 0.
            UNS = 0.
        END IF

        IF (DUMG(K).GE.QSMALL) THEN
            UMG = AGN(K)*CONS7/(DLAMG**BG)
            UNG = AGN(K)*CONS8/DLAMG**BG
        ELSE
            UMG = 0.
            UNG = 0.
        END IF

        ! SET REALISTIC LIMITS ON FALLSPEED

        ! bug fix, 10/08/09
        dum=(rhosu/rho(k))**0.54
        UMS=MIN(UMS,1.2*dum)
        UNS=MIN(UNS,1.2*dum)
        ! fix 053011
        ! fix for correction by AA 4/6/11
        UMI=MIN(UMI,1.2*(rhosu/rho(k))**0.35)
        UNI=MIN(UNI,1.2*(rhosu/rho(k))**0.35)
        UMR=MIN(UMR,9.1*dum)
        UNR=MIN(UNR,9.1*dum)
        UMC=MIN(UMC,9.1*dum)
        UNC=MIN(UNC,9.1*dum)
        UMG=MIN(UMG,20.*dum)
        UNG=MIN(UNG,20.*dum)

        FR(K) = UMR
        FI(K) = UMI
        FNI(K) = UNI
        FS(K) = UMS
        FNS(K) = UNS
        FNR(K) = UNR
        FC(K) = UMC
        FNC(K) = UNC
        FG(K) = UMG
        FNG(K) = UNG

        ! Assign to umx
        if (ipelec>0) then
            umx(k,lc) = umc
            umx(k,lr) = umr
            umx(k,li) = umi
            umx(k,ls) = ums
            !            umx(k,lh) = umg*8.
            umx(k,lh) = umg
            umc1d(k)  = umc
            umr1d(k)  = umr
            umi1d(k)  = umi
            ums1d(k)  = ums
            umg1d(k)  = umg
        end if

        ! V3.3 MODIFY FALLSPEED BELOW LEVEL OF PRECIP

        IF (K.LE.KTE-1) THEN
            IF (FR(K).LT.1.E-10) THEN
                FR(K)=FR(K+1)
            END IF
            IF (FI(K).LT.1.E-10) THEN
                FI(K)=FI(K+1)
            END IF
            IF (FNI(K).LT.1.E-10) THEN
                FNI(K)=FNI(K+1)
            END IF
            IF (FS(K).LT.1.E-10) THEN
                FS(K)=FS(K+1)
            END IF
            IF (FNS(K).LT.1.E-10) THEN
                FNS(K)=FNS(K+1)
            END IF
            IF (FNR(K).LT.1.E-10) THEN
                FNR(K)=FNR(K+1)
            END IF
            IF (FC(K).LT.1.E-10) THEN
                FC(K)=FC(K+1)
            END IF
            IF (FNC(K).LT.1.E-10) THEN
                FNC(K)=FNC(K+1)
            END IF
            IF (FG(K).LT.1.E-10) THEN
                FG(K)=FG(K+1)
            END IF
            IF (FNG(K).LT.1.E-10) THEN
                FNG(K)=FNG(K+1)
            END IF
        END IF ! K LE KTE-1

        ! CALCULATE NUMBER OF SPLIT TIME STEPS

        RGVM = MAX(FR(K),FI(K),FS(K),FC(K),FNI(K),FNR(K),FNS(K),FNC(K),FG(K),FNG(K))
        ! VVT CHANGED IFIX -> INT (GENERIC FUNCTION)
        NSTEP = MAX(INT(RGVM*DT/DZQ(K)+1.),NSTEP)

        ! MULTIPLY VARIABLES BY RHO
        DUMR(k) = DUMR(k)*RHO(K)
        DUMI(k) = DUMI(k)*RHO(K)
        DUMFNI(k) = DUMFNI(K)*RHO(K)
        DUMQS(k) = DUMQS(K)*RHO(K)
        DUMFNS(k) = DUMFNS(K)*RHO(K)
        DUMFNR(k) = DUMFNR(K)*RHO(K)
        DUMC(k) = DUMC(K)*RHO(K)
        DUMFNC(k) = DUMFNC(K)*RHO(K)
        DUMG(k) = DUMG(K)*RHO(K)
        DUMFNG(k) = DUMFNG(K)*RHO(K)

        !-----------------------------------------------------------------------
        !Explicit lightning ----------------------------------------------------
        if (ipelec>0) then
            DUMSC(K,lsc(LC)) = an1d(K,lsc(LC))
            DUMSC(K,lsc(LR)) = an1d(K,lsc(LR))
            DUMSC(K,lsc(LI)) = an1d(K,lsc(LI))
            DUMSC(K,lsc(LS)) = an1d(K,lsc(LS))
            DUMSC(K,lsc(LG)) = an1d(K,lsc(LG))
        end if
        !End explicit lightning ------------------------------------------------
        !-----------------------------------------------------------------------
    END DO

    DO N = 1,NSTEP

        DO K = KTS,KTE
            FALOUTR(K) = FR(K)*DUMR(K)
            FALOUTI(K) = FI(K)*DUMI(K)
            FALOUTNI(K) = FNI(K)*DUMFNI(K)
            FALOUTS(K) = FS(K)*DUMQS(K)
            FALOUTNS(K) = FNS(K)*DUMFNS(K)
            FALOUTNR(K) = FNR(K)*DUMFNR(K)
            FALOUTC(K) = FC(K)*DUMC(K)
            FALOUTNC(K) = FNC(K)*DUMFNC(K)
            FALOUTG(K) = FG(K)*DUMG(K)
            FALOUTNG(K) = FNG(K)*DUMFNG(K)
            !-----------------------------------------------------------------------
            !Explicit lightning ----------------------------------------------------
            if (ipelec>0) then
                FALOUTSCC(K) = FC(K)*DUMSC(K,lsc(LC))
                FALOUTSCR(K) = FR(K)*DUMSC(K,lsc(LR))
                FALOUTSCI(K) = FI(K)*DUMSC(K,lsc(LI))
                FALOUTSCS(K) = FS(K)*DUMSC(K,lsc(LS))
                FALOUTSCG(K) = FG(K)*DUMSC(K,lsc(LH))
                tmpscr(k) = DUMSC(K,lsc(LR))
            end if
            !End explicit lightning ------------------------------------------------
            !-----------------------------------------------------------------------
        END DO

        ! TOP OF MODEL

        K = KTE
        FALTNDR = FALOUTR(K)/DZQ(k)
        FALTNDI = FALOUTI(K)/DZQ(k)
        FALTNDNI = FALOUTNI(K)/DZQ(k)
        FALTNDS = FALOUTS(K)/DZQ(k)
        FALTNDNS = FALOUTNS(K)/DZQ(k)
        FALTNDNR = FALOUTNR(K)/DZQ(k)
        FALTNDC = FALOUTC(K)/DZQ(k)
        FALTNDNC = FALOUTNC(K)/DZQ(k)
        FALTNDG = FALOUTG(K)/DZQ(k)
        FALTNDNG = FALOUTNG(K)/DZQ(k)
        ! ADD FALLOUT TERMS TO EULERIAN TENDENCIES

        QRSTEN(K) = QRSTEN(K)-FALTNDR/NSTEP/RHO(k)
        QISTEN(K) = QISTEN(K)-FALTNDI/NSTEP/RHO(k)
        NI3DTEN(K) = NI3DTEN(K)-FALTNDNI/NSTEP/RHO(k)
        QNISTEN(K) = QNISTEN(K)-FALTNDS/NSTEP/RHO(k)
        NS3DTEN(K) = NS3DTEN(K)-FALTNDNS/NSTEP/RHO(k)
        NR3DTEN(K) = NR3DTEN(K)-FALTNDNR/NSTEP/RHO(k)
        QCSTEN(K) = QCSTEN(K)-FALTNDC/NSTEP/RHO(k)
        NC3DTEN(K) = NC3DTEN(K)-FALTNDNC/NSTEP/RHO(k)
        QGSTEN(K) = QGSTEN(K)-FALTNDG/NSTEP/RHO(k)
        NG3DTEN(K) = NG3DTEN(K)-FALTNDNG/NSTEP/RHO(k)

        DUMR(K) = DUMR(K)-FALTNDR*DT/NSTEP
        DUMI(K) = DUMI(K)-FALTNDI*DT/NSTEP
        DUMFNI(K) = DUMFNI(K)-FALTNDNI*DT/NSTEP
        DUMQS(K) = DUMQS(K)-FALTNDS*DT/NSTEP
        DUMFNS(K) = DUMFNS(K)-FALTNDNS*DT/NSTEP
        DUMFNR(K) = DUMFNR(K)-FALTNDNR*DT/NSTEP
        DUMC(K) = DUMC(K)-FALTNDC*DT/NSTEP
        DUMFNC(K) = DUMFNC(K)-FALTNDNC*DT/NSTEP
        DUMG(K) = DUMG(K)-FALTNDG*DT/NSTEP
        DUMFNG(K) = DUMFNG(K)-FALTNDNG*DT/NSTEP

        !-----------------------------------------------------------------------
        !Explicit lightning ----------------------------------------------------
        if (ipelec>0) then
            FALTNDSCC = FALOUTSCC(K)/DZQ(K)
            FALTNDSCR = FALOUTSCR(K)/DZQ(K)
            FALTNDSCI = FALOUTSCI(K)/DZQ(K)
            FALTNDSCS = FALOUTSCS(K)/DZQ(K)
            FALTNDSCG = FALOUTSCG(K)/DZQ(K)

            PSCC_SED(K) = PSCC_SED(K)-FALTNDSCC/NSTEP
            PSCR_SED(K) = PSCR_SED(K)-FALTNDSCR/NSTEP
            PSCI_SED(K) = PSCI_SED(K)-FALTNDSCI/NSTEP
            PSCS_SED(K) = PSCS_SED(K)-FALTNDSCS/NSTEP
            PSCG_SED(K) = PSCG_SED(K)-FALTNDSCG/NSTEP

            DUMSC(K,lsc(LC)) = DUMSC(K,lsc(LC)) - FALTNDSCC*DT/NSTEP
            DUMSC(K,lsc(LR)) = DUMSC(K,lsc(LR)) - FALTNDSCR*DT/NSTEP
            DUMSC(K,lsc(LI)) = DUMSC(K,lsc(LI)) - FALTNDSCI*DT/NSTEP
            DUMSC(K,lsc(LS)) = DUMSC(K,lsc(LS)) - FALTNDSCS*DT/NSTEP
            DUMSC(K,lsc(LG)) = DUMSC(K,lsc(LG)) - FALTNDSCG*DT/NSTEP
        end if
        !End explicit lightning ------------------------------------------------
        !-----------------------------------------------------------------------

        DO K = KTE-1,KTS,-1
            FALTNDR = (FALOUTR(K+1)-FALOUTR(K))/DZQ(K)
            FALTNDI = (FALOUTI(K+1)-FALOUTI(K))/DZQ(K)
            FALTNDNI = (FALOUTNI(K+1)-FALOUTNI(K))/DZQ(K)
            FALTNDS = (FALOUTS(K+1)-FALOUTS(K))/DZQ(K)
            FALTNDNS = (FALOUTNS(K+1)-FALOUTNS(K))/DZQ(K)
            FALTNDNR = (FALOUTNR(K+1)-FALOUTNR(K))/DZQ(K)
            FALTNDC = (FALOUTC(K+1)-FALOUTC(K))/DZQ(K)
            FALTNDNC = (FALOUTNC(K+1)-FALOUTNC(K))/DZQ(K)
            FALTNDG = (FALOUTG(K+1)-FALOUTG(K))/DZQ(K)
            FALTNDNG = (FALOUTNG(K+1)-FALOUTNG(K))/DZQ(K)

            ! ADD FALLOUT TERMS TO EULERIAN TENDENCIES

            QRSTEN(K) = QRSTEN(K)+FALTNDR/NSTEP/RHO(k)
            QISTEN(K) = QISTEN(K)+FALTNDI/NSTEP/RHO(k)
            NI3DTEN(K) = NI3DTEN(K)+FALTNDNI/NSTEP/RHO(k)
            QNISTEN(K) = QNISTEN(K)+FALTNDS/NSTEP/RHO(k)
            NS3DTEN(K) = NS3DTEN(K)+FALTNDNS/NSTEP/RHO(k)
            NR3DTEN(K) = NR3DTEN(K)+FALTNDNR/NSTEP/RHO(k)
            QCSTEN(K) = QCSTEN(K)+FALTNDC/NSTEP/RHO(k)
            NC3DTEN(K) = NC3DTEN(K)+FALTNDNC/NSTEP/RHO(k)
            QGSTEN(K) = QGSTEN(K)+FALTNDG/NSTEP/RHO(k)
            NG3DTEN(K) = NG3DTEN(K)+FALTNDNG/NSTEP/RHO(k)

            DUMR(K) = DUMR(K)+FALTNDR*DT/NSTEP
            DUMI(K) = DUMI(K)+FALTNDI*DT/NSTEP
            DUMFNI(K) = DUMFNI(K)+FALTNDNI*DT/NSTEP
            DUMQS(K) = DUMQS(K)+FALTNDS*DT/NSTEP
            DUMFNS(K) = DUMFNS(K)+FALTNDNS*DT/NSTEP
            DUMFNR(K) = DUMFNR(K)+FALTNDNR*DT/NSTEP
            DUMC(K) = DUMC(K)+FALTNDC*DT/NSTEP
            DUMFNC(K) = DUMFNC(K)+FALTNDNC*DT/NSTEP
            DUMG(K) = DUMG(K)+FALTNDG*DT/NSTEP
            DUMFNG(K) = DUMFNG(K)+FALTNDNG*DT/NSTEP


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

            ! FOR WRF-CHEM, NEED PRECIP RATES (UNITS OF KG/M^2/S)
            CSED(K)=CSED(K)+FALOUTC(K)/NSTEP
            ISED(K)=ISED(K)+FALOUTI(K)/NSTEP
            SSED(K)=SSED(K)+FALOUTS(K)/NSTEP
            GSED(K)=GSED(K)+FALOUTG(K)/NSTEP
            RSED(K)=RSED(K)+FALOUTR(K)/NSTEP
        END DO

        ! GET PRECIPITATION AND SNOWFALL ACCUMULATION DURING THE TIME STEP
        ! FACTOR OF 1000 CONVERTS FROM M TO MM, BUT DIVISION BY DENSITY
        ! OF LIQUID WATER CANCELS THIS FACTOR OF 1000

        PRECRT = PRECRT+(FALOUTR(KTS)+FALOUTC(KTS)+FALOUTS(KTS)+FALOUTI(KTS)+FALOUTG(KTS))  &
                *DT/NSTEP
        SNOWRT = SNOWRT+(FALOUTS(KTS)+FALOUTI(KTS)+FALOUTG(KTS))*DT/NSTEP
        ! hm added 7/13/13
        SNOWPRT = SNOWPRT+(FALOUTI(KTS)+FALOUTS(KTS))*DT/NSTEP
        GRPLPRT = GRPLPRT+(FALOUTG(KTS))*DT/NSTEP

    END DO

    DO K=KTS,KTE

        ! ADD ON SEDIMENTATION TENDENCIES FOR MIXING RATIO TO REST OF TENDENCIES

        QR3DTEN(K)=QR3DTEN(K)+QRSTEN(K)
        QI3DTEN(K)=QI3DTEN(K)+QISTEN(K)
        QC3DTEN(K)=QC3DTEN(K)+QCSTEN(K)
        QG3DTEN(K)=QG3DTEN(K)+QGSTEN(K)
        QNI3DTEN(K)=QNI3DTEN(K)+QNISTEN(K)

        ! PUT ALL CLOUD ICE IN SNOW CATEGORY IF MEAN DIAMETER EXCEEDS 2 * dcs

        !hm 4/7/09 bug fix
        !        IF (QI3D(K).GE.QSMALL.AND.T3D(K).LT.273.15) THEN
        IF (QI3D(K).GE.QSMALL.AND.T3D(K).LT.273.15.AND.LAMI(K).GE.1.E-10) THEN
            IF (1./LAMI(K).GE.2.*DCS) THEN
                QNI3DTEN(K) = QNI3DTEN(K)+QI3D(K)/DT+ QI3DTEN(K)
                NS3DTEN(K) = NS3DTEN(K)+NI3D(K)/DT+   NI3DTEN(K)
                QI3DTEN(K) = -QI3D(K)/DT
                NI3DTEN(K) = -NI3D(K)/DT
            END IF
        END IF

        ! ⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡
        ! yc added 9/21/2020. for elec
        IF (IPELEC>0) THEN

            if (csi_morr) then
                !                qxmin(ls) = 1.E-8
                !                qxmin(lh) = 1.E-8
                !                qxmin(lc) = QSMALL
                !!                qxmin(lr) = 1.E-8
                !                qxmin(li) = QSMALL
            end if

            !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            IF (QI3D(K).GE.QXMIN(LI) .AND. QNI3D(K).GE.QXMIN(LS) .AND. not(csi_morr)) THEN
                !            IF (QNI3D(K).GE.1.E-8 .AND. QI3D(K).GE.QSMALL) THEN
                IF (NI3D(K).GT.QSMALL .AND. NS3D(K).GT.QSMALL) THEN
                    XV_I = QI3D(K)/NI3D(K)/RHOI
                    XV_S = QNI3D(K)/NS3D(K)/RHOSN

                    ! Graupel-ice collision efficiency
                    ! Currently used as constant 1.
                    ESICLSN = ESI_COLLSN

                    CSACI0(K) = ESICLSN * RVT * AA2 * &
                            (NS3D(K)*RHO(K)) * (NI3D(K)*RHO(K)) * &
                            ((CINU + 2.)*XV_I/(CINU + 1.) + XV_S)
                END IF
            END IF

            !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            ! Collision rate of graupel-ice (CHACI0)
            IF (QI3D(K).GE.QXMIN(LI) .AND. QG3D(K).GE.QXMIN(LH)) THEN
                IF (NG3D(K).GT.QSMALL .AND. NI3D(K).GT.QSMALL) THEN
                    vt = sqrt((umx(k,lh)-umx(k,li))**2 + 0.04*umx(k,lh)*umx(k,li))

                    DIAG = xdia(k,lh)
                    DIAI = xdia(k,li)

                    ! Graupel-ice collision efficiency
                    ! Currently used as constant 1.
                    EHICLSN = EHI_COLLSN

                    CHACI0(K) = 0.25 * PI * EHICLSN * &
                            (NG3D(K)*RHO(K)) * (NI3D(K)*RHO(K)) * VT * &
                            (DA0LH*DIAG**2 + &
                                    DAB0LH_IH*DIAG*DIAI + &
                                    DA0LI*DIAI**2)
                END IF
            END IF

            !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            ! Collision rate of graupel-snow (CHACS0)
            IF (QNI3D(K).GE.QXMIN(LS) .AND. QG3D(K).GE.QXMIN(LH)) THEN
                IF (NG3D(K).GT.QSMALL .AND. NS3D(K).GT.QSMALL) THEN
                    vt = sqrt((umx(k,lh)-umx(k,ls))**2 + 0.04*umx(k,lh)*umx(k,ls))

                    DIAG = xdia(k,lh)
                    DIAS = xdia(k,ls)

                    ! Graupel-snow collision efficiency.
                    ! Defined as a linear function of snow size.
                    ! Zero at Ds<=40um, 0.5 at Ds>=150um
                    EHSCLSN = EHS_COLLSN
                    IF ( DIAS < 40.E-6 ) THEN
                        EHSCLSN = 0.0
                    ELSEIF ( DIAS < 150.E-6 ) THEN
                        EHSCLSN =  EHS_COLLSN*(DIAS - 40.E-6)/(150.E-6 - 40.E-6)
                    ELSE
                        EHSCLSN = EHS_COLLSN
                    ENDIF

                    CHACS0(K) = 0.25 * PI * EHSCLSN * &
                            (NG3D(K)*RHO(K)) * (NS3D(K)*RHO(K)) * VT * &
                            (DA0LH*DIAG**2 + &
                                    DAB0LH_SH*DIAG*DIAS + &
                                    DA0LS*DIAS**2)
                END IF
            END IF


            !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            ! update an1d "q and n" part (C/m3 for charge)
            an1d(k,lc) = qc3d(k) ! kg/kg
            an1d(k,lr) = qr3d(k)
            an1d(k,li) = qi3d(k)
            an1d(k,ls) = qni3d(k)
            an1d(k,lh) = qg3d(k)

            an1d(k,lnc) = nc3d(k)*RHO(K) ! #/m3 consistent w/ nssl
            an1d(k,lnr) = nr3d(k)*RHO(K)
            an1d(k,lni) = ni3d(k)*RHO(K)
            an1d(k,lns) = ns3d(k)*RHO(K)
            an1d(k,lnh) = ng3d(k)*RHO(K)

            !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            ! Load charge into scx array
            DO il = lc,lhab
                scx(k,il) = an1d(k,lsc(il)) ! C/m3
                IF ( Abs(scx(k,il)) > 1000.e-9 ) THEN
                    !                    write(0,*) 'Problem0a with scx il = ',il,lsc(il),scx(k,il),an1d(k,il),an1d(k,ln(il))
                    !                    write(0,*) 'at it, ix,jy,kz = ',itimestep ,i, j, k
                    !                    STOP
                ENDIF
            ENDDO
            cionp(k) = an1d(k,lscpi) ! #/m3
            cionn(k) = an1d(k,lscni)
            if ( largeion ) then
                clionp(k) = an1d(k,lscpli)
                clionn(k) = an1d(k,lscnli)
            end if

            !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            ! Threshold for ion charge
            frac = 0.1d0
            cionpmxd(k) = frac*cionp(k)*dtpinv
            cionnmxd(k) = frac*cionn(k)*dtpinv
            if ( largeion ) then
                clionpmxd(k) = frac*clionp(k)*dtpinv
                clionnmxd(k) = frac*clionn(k)*dtpinv
            end if

            !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            ! Collisional charging rates
            ! Sticking eff. Collision eff assumed 1
            exy(k,ls,lc) = eci ! 0.7
            exy(k,ls,li) = eii ! 0.1
            exy(k,lh,lc) = eci ! 0.7
            exy(k,lh,li) = 0.0
            exy(k,lh,ls) = 0.0

            ni = 1.0
            cxacy(k,ls,lc) = 0.0 ! not used in nssl anyway
            cxacy(k,lh,lc) = 0.0 ! not used in nssl anyway
            if (csi_morr) then
                cxacy(k,ls,li) = ni*csaci0_morr(k)*rho(k) ! if use morr_eq, convert to m-3
            else
                cxacy(k,ls,li) = ni*csaci0(k) ! if use nssl_eq, no need to convert
            end if
            cxacy(k,lh,li) = ni*chaci0(k)
            cxacy(k,lh,ls) = chacs0(k)

            temcg(k) = t3d(k)-273.15

            IF ( ipelec .ge. 2 .or. idonic ) THEN
                IF ( j.gt.(jds+2) .and. j.lt.(jde-2) ) THEN
                    DO il=ls,lhab ! loop through all ice habits (s,h)
                        DO ic=li,ls ! loop through all mixed phase particles (i,s)
                            !                    il=lh
                            !                    ic=li
                            IF ( iexy(il,ic) .eq. 1  ) THEN ! (si,hi,hs)
                                IF (an1d(k,il) .gt. qxmin(il)  .and. an1d(k,ic) .gt. qxmin(ic) &
                                        .and. an1d(k,ln(il)).gt.qsmall .and. an1d(k,ln(ic)).gt.qsmall) THEN  ! make sure there are particles
                                    IF ( t3d(k) .lt. tfr .and. t3d(k) .gt. thnuc ) THEN ! 0.0=tfreezing and -38C= assumed T for homogeneous ice nucleation
                                        IF ( iremoveqwfrz == 1 ) THEN
                                            qcwtmp(k) = Max(0.0, an1d(k,lc) - dtp*MNUCCC(k) )
                                        ELSE
                                            qcwtmp(k) = an1d(k,lc)
                                        ENDIF

                                        ! Currently set nonigrd = 0 (saundx)
                                        IF (nonigrd .eq. -1 ) THEN
                                            vt = umx(k,il) - umx(k,ic)
                                            CALL takax(isaund,ntem,nlwc,takalu,temcg(k),qcwtmp(k),vt,   &
                                                    &     xdia(k,ic),rho(k),ftelwc,exy(k,il,lc),an1d(k,lr),rarfac)
                                            scxacy(k,il,ic) = ftelwc
                                            if ( scxacy(k,il,ic) .gt. 0.0 ) then
                                                scxacy(k,il,ic) = min(delqxxa(ic),scxacy(k,il,ic))
                                            end if
                                            if ( scxacy(k,il,ic) .lt. 0.0 ) then
                                                scxacy(k,il,ic) = max(delqnxa(ic),scxacy(k,il,ic))
                                            end if
                                        ELSEIF ( nonigrd .eq. 1 .or. nonigrd .eq. 0 ) THEN
                                            IF (nonigrd .eq. 0) THEN ! currently used
                                                qconkq = 0.0
                                                qconm = 1.0
                                                qconn = 1.0
                                                ftelwc = 0.0
                                                vt = umx(k,il) - umx(k,lc)
                                                call saundx(isaund,temcg(k),qcwtmp(k),exy(k,il,lc),vt,xdia(k,ic),   &
                                                        &   rho(k),qsign,ftelwc,qconkq,qconm,qconn,idelq,rarfac)
                                            ELSE ! nonigrd=1
                                                vt = umx(k,il) - umx(k,lc)
                                                ! rarx only done for graup (as in nssl), doesn't matter since nonigrd is not 1 for now
                                                rarx(k,lh) = (PSACWG(k)*rho(k))*1.0e3*rho(k)/((pi/2.0)*xdia(k,lh)*an1d(k,lnh)) ! ! PSACWG: equate to qhacw cw collected by gw, cw->gw
                                                call saundy(isaund,temcg(k),qcwtmp(k),rarx(k,il),vt,xdia(k,ic),   &
                                                        &   rho(k),qsign,ftelwc,qconkq,qconm,qconn,idelq)
                                            ENDIF

                                            cfce15 = 1.0e-15
                                            scxacy(k,il,ic) = qconkq*(xdia(k,ic)**qconm)   &
                                                    &  *((abs(umx(k,il)-umx(k,ic)))**qconn)   &
                                                    &  *ftelwc*cfce15
                                            if ( scxacy(k,il,ic) .gt. 0.0 ) then ! {
                                                scxacy(k,il,ic) = min(delqxxa(ic),scxacy(k,il,ic))
                                            end if
                                            if ( scxacy(k,il,ic) .lt. 0.0 ) then
                                                scxacy(k,il,ic) = max(delqnxa(ic),scxacy(k,il,ic))
                                            end if
                                        end if

                                        ! Formula 7 in M05 to get total space charge per collisions within a grid box
                                        tmp = scxacy(k,il,ic)
                                        scxacy(k,il,ic) =   &
                                                &  cxacy(k,il,ic)*ecollmx*(1.0-exy(k,il,ic))*scxacy(k,il,ic)
                                        IF ( Abs( scxacy(k,il,ic) ) .gt. scxacymax ) THEN ! scxacymax=3000e-12
                                            scxacy(k,il,ic) = Min(scxacymax, Max( -scxacymax, scxacy(k,il,ic) ) )
                                        ENDIF

                                    END IF ! temperture allows mixed phase

                                    ! charging with no liquid water (Mitzeva).
                                    ! not used in nssl

                                END IF ! an1d
                            END IF ! iexy(il,ic) .eq. 1
                        END DO ! ic
                    END DO ! il
                    !
                    !  End of ice-ice loop
                    !

                    ! Loop for graupel-droplet inductive charging
                    ! Mansell et al 2005, adapted from Ziegler et al 1991
                    IF ( erbnd .gt. 0.0 .and. fdgt .gt. 0.0 .and. ipelec .ge. 3 ) THEN ! erbnd=0.01, fdgt=1
                        DO il=ls,lhab
                            !                        il=lh
                            IF ( iexy(il,lc) .eq. 1 ) THEN
                                qcwtmp(k) = Max(0.0, an1d(k,lc) - dtp*MNUCCC(k) )
                                if ( an1d(k,il) .gt. qxmin(il) .and. qcwtmp(k) .gt. qxmin(lc)    &
                                        & .and. xdia(k,lc).gt.0.0 &
                                        & .and. exy(k,il,lc) .gt. 0.0 .and. an1d(k,ln(il)) > 1.e-6 &
                                        & .and. an1d(k,lnc)>qsmall) THEN
                                    IF ( t3d(k) .gt. tindmn .and. t3d(k) .lt. tindmx ) then ! -40 to 25 C
                                        scxacy(k,il,lc) =  &
                                                &  0.125*pi**3*exy(k,il,lc)*erbnd*fdgt*an1d(k,ln(il))   &
                                                &  *an1d(k,lnc)*umx(k,il)*6.0/gf4p5*(xdia(k,lc)**2)   &
                                                &  *( pi*eperao*costhe*dezcomp(k)   &
                                                        &      *gf3p5*(xdia(k,il)**2)   &
                                                        &     - (1.0/3.0)*gf1p5*an1d(k,lsc(il))/an1d(k,ln(il)) ) ! must use from previous time step

                                        if (abs(scxacy(k,il,lc)).gt.delqxw) then
                                            write(0,*) 'i,j,k,t=', i,j,k,itimestep
                                            write(0,*) 'il=',il,' scxacy(k,il,lc)=',scxacy(k,il,lc)
                                            write(0,*) '0.125*pi**3=',0.125*pi**3
                                            write(0,*) 'exy(k,il,lc)=',exy(k,il,lc)
                                            write(0,*) 'erbnd*fdgt=',erbnd*fdgt
                                            write(0,*) '0.125*pi**3*exy(k,il,lc)*erbnd*fdgt=',0.125*pi**3*exy(k,il,lc)*erbnd*fdgt
                                            write(0,*) 'an1d(k,ln(il)), an1d(k,lnc)=',an1d(k,ln(il)), an1d(k,lnc)
                                            write(0,*) 'umx(k,il), xdia(k,lc), xdia(k,il)', umx(k,il), xdia(k,lc), xdia(k,il)
                                            write(0,*) '6.0/gf4p5=', 6.0/gf4p5
                                            write(0,*) 'pi*eperao*costhe=',pi*eperao*costhe
                                            write(0,*) 'dezcomp(k)=',dezcomp(k)
                                            write(0,*) 'gf3p5*(xdia(k,il)**2=',gf3p5*xdia(k,il)**2
                                            write(0,*) '(1.0/3.0)*gf1p5=',(1.0/3.0)*gf1p5
                                            write(0,*) 'an1d(k,lsc(il)), an1d(k,ln(il))',an1d(k,lsc(il)), an1d(k,ln(il))
                                            write(0,*) 'an1d(k,lsc(il))/an1d(k,ln(il))=',an1d(k,lsc(il))/an1d(k,ln(il))
                                            !!                                            STOP
                                        end if
                                        if ( scxacy(k,il,lc) .gt. 0.0 ) then
                                            scxacy(k,il,lc) = min(delqxw,scxacy(k,il,lc))
                                        elseif ( scxacy(k,il,lc) .lt. 0.0 ) then
                                            scxacy(k,il,lc) = max(delqnw,scxacy(k,il,lc))
                                        end if

                                    ELSE
                                        scxacy(k,il,lc) = 0.0
                                    end if  ! temg
                                ENDIF ! an1d
                            ENDIF ! iexy
                        ENDDO  ! il
                    ENDIF ! erbnd

                ELSE
                    scsacw(k) = 0.0
                    scsaci(k) = 0.0
                    scsacr(k) = 0.0

                    schacw(k) = 0.0
                    schaci(k) = 0.0
                    schacs(k) = 0.0
                    schacr(k) = 0.0
                    scxacy(k,:,:)=0.0
                END IF ! boundary check on y
            END IF ! ipelec>=2

            IF ( i.gt.(ids+2) .and. i.lt.(ide-2)  ) THEN
                scsacw(k) = elecfac*scxacy(k,ls,lc) ! induc
                scsaci(k) = elecfac*scxacy(k,ls,li) ! non-induc
                scsacr(k) = elecfac*scxacy(k,ls,lr) ! =0

                schacw(k) = elecfac*scxacy(k,lh,lc) ! induc
                schaci(k) = elecfac*scxacy(k,lh,li) ! non-induc
                schacs(k) = elecfac*scxacy(k,lh,ls) ! non-induc
                schacr(k) = elecfac*scxacy(k,lh,lr) ! =0


            ELSE

                scsacw(k) = 0.0
                scsaci(k) = 0.0
                scsacr(k) = 0.0

                schacw(k) = 0.0
                schaci(k) = 0.0
                schacs(k) = 0.0
                schacr(k) = 0.0
                scxacy(k,:,:)=0.0
            END IF ! x boundary check


            !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            ! Charge per hydrometeor content
            if ( ipelec .ge. 1 .and. itest .eq. 1 ) then
                if ( an1d(k,lr) .gt. qxmin(lr) ) then
                    fscrw(k) = scx(k,lr)/an1d(k,lr)
                    tmp1=scx(k,lr)
                end if
                if ( an1d(k,lc) .gt. qxmin(lc) ) then
                    fsccw(k) = scx(k,lc)/an1d(k,lc)
                end if
                if ( an1d(k,li) .gt. qxmin(li) ) then
                    fscci(k) = scx(k,li)/an1d(k,li)
                end if
                if ( an1d(k,ls) .gt. qxmin(ls) ) then
                    fscsw(k) = scx(k,ls)/an1d(k,ls)
                end if
                if ( an1d(k,lh) .gt. qxmin(lh) ) then
                    fschw(k) = scx(k,lh)/an1d(k,lh)
                end if
            end if

            !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            ! il5 =0, warm; =1, frozen
            il5(k) = 0
            if ( t3d(k) .lt. 273.15 ) then
                il5(k) = 1
            end if


            !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            ! Charge production rates (pscxxi,pscxxd,pscxxmi,pscxxmd)
            !
            !
            !  Charge production from mass transfer
            !
            !
            if ( itest .eq. 1 ) then
                !
                !  Cloud ice
                !
                psccimi(k) =    il5(k)*(fsccw(k)*PSACWI(k) & ! PSACWI: equate to qiacw: collection of cw by ci, cw->ci
                        !                        &  +fsccw(k)*qwfrzc(k)   & ! qwfrzc: homogeneous frz in nssl, not used since this is done later instantaneously
                        &  +fsccw(k)*mnuccc(k) & ! mnuccc: equate to qwctfzc, here includs contact and immersion frz, cw->ci
                        &  +fsccw(k)*(qmults(k)+qmultg(k)) & ! h-m ice multip. on sw & gw, cw->ci
                        &  +fscrw(k)*(qmultr(k)+qmultrg(k)) & ! h-m ice multip. on sw & gw, rw->ci
                        )
                !                        &  +fsccw(k)*qicichr(k))   & ! qicichr: Hobbs-Rangno ice enhancement, not used in nssl either
                !                        &  +il5(k)*fscrw(k)*(1. - ifrzs)*qrfrzs(k)   & ! not used since ifrzs=1 in nssl
                !                        &  +fschw(k)*qhmul1(k)   & ! qhmul1: in nssl, hw->ci, hallet-mossop ice multip. from graup
                !                        &  +fschl(k)*qhlmul1(k)  & ! hail not considered
                !                        &  +fscsw(k)*qsmul(k) & ! qsmul: in nssl, sw->ci, hallet-mossop ice multip. from snow

                psccimd(k) =    fscci(k)*(&
                        &   il5(k)*(-PRCI(k)   & ! prci: equate to qscni, autoconversion, ci->sw
                                &           -PRAI(k)   & ! prai: equate to qsaci, ci collected by sw, ci->sw
                                &           -praci(k)  & ! praci: equate to qraci, ci collected by rw to produce gw, ci->gw
                                &           -pracis(k) & ! PRACIS: equate to qracis, ci collected by rw to produce sw, ci->sw
                                &           +min(0.,eprd(k)) & ! EPRD: neg, sublimation of ice, v->ci
                                &          ) &
                        !                        &  -qraci(k)   & ! qraci: ci->rw ci collected by rain
                        !                        &  -qhaci(k)   & ! qraci: ci->gw ci collected by graup
                        !                        &  -qhlaci(k) & ! collection of ci by hail, no hail, not used
                        )
                ! conservation of ice charge
                if (abs(psccimd(k)).gt.0.0) then
                    if (sign(1.0,psccimd(k)).eq.sign(1.0,scx(k,li))) then
                        write(0,*) 'psccimd scx same sign', psccimd(k), scx(k,li)
                        write(0,*) 'i,j,k,t=', i,j,k,itimestep
                        write(0,*) 'fscci(k),scx(k,li),il5(k)',fscci(k),scx(k,li),il5(k)
                        write(0,*) 'stopping for psccimd scxi same sign'
                        stop
                    end if
                    if (abs(psccimd(k)*dt).gt.abs(scx(k,li))) then
                        write(0,*) 'conserving pscci', psccimd(k), dt, psccimd(k)*dt, scx(k,li)
                        write(0,*) 'fscci(k)',fscci(k)
                        write(0,*) 'i,j,k,t=', i,j,k,itimestep
                        write(0,*) 'il5(k)',il5(k)
                        write(0,*) 'an1d(k,li)',an1d(k,li)
                        dumtmp = il5(k)*(-PRCI(k)   & ! prci: equate to qscni, autoconversion, ci->sw
                                &           -PRAI(k)   & ! prai: equate to qsaci, ci collected by sw, ci->sw
                                &           -praci(k)  & ! praci: equate to qraci, ci collected by rw to produce gw, ci->gw
                                &           -pracis(k) & ! PRACIS: equate to qracis, ci collected by rw to produce sw, ci->sw
                                &           +min(0.,eprd(k)) & ! EPRD: neg, sublimation of ice, v->ci
                                &          )
                        write(0,*) 'dumtmp',dumtmp
                        if (an1d(k,li)+dumtmp*dt.lt.0.0) then
                            write(0,*) 'neg qi in psc conserve',an1d(k,li)+dumtmp*dt
                            write(0,*) 'an1d(k,li),dumtmp*dt',an1d(k,li),dumtmp*dt
                        else if (an1d(k,li)+dumtmp*dt.gt.0.0) then
                            write(0,*) 'pos qi in pscci conserve',an1d(k,li)+dumtmp*dt
                            write(0,*) 'an1d(k,li),dumtmp*dt',an1d(k,li),dumtmp*dt
                        else
                            write(0,*) 'zero qi+dumtmp*dt',an1d(k,li)+dumtmp*dt
                        end if
                        psccimd(k) = -scx(k,li)*dtpinv
                        write(0,*) 'conserved pscci', psccimd(k), dt, psccimd(k)*dt, scx(k,li)
                    end if
                end if

                !
                !  Cloud water
                !
                psccwmi(k) = 0.0
                psccwmd(k) =    fsccw(k)*(&
                        &  -PRA(k) & ! PRA: equate to qracw, cw collected by rw, cw->rw
                                &  -PRC(k) & ! PRC: equate to qrcnw, cw autoconvert to rw, cw->rw
                                &  +il5(k)*(-PSACWI(k)   & ! PSACWI: equate to qiacw: collection of cw by ci
                                        &  -mnuccc(k)    & ! mnuccc: equate to qwctfzc, here includs contact and immersion frz
                                        &  -PGSACW(k)  & ! PGSACW: cw collected by sw to produce gw, cw->gw
                                        &  -qmults(k)-qmultg(k) & ! h-m ice multip. on sw & gw, cw->ci
                                        &  -PSACWS(k)  & ! PSACWS: equate to qsacw, cw collected by sw, after deducted by PGSACW and QMULTS, cw->sw
                                        &  -PSACWG(k)  & ! PSACWG: equate to qhacw, cw collected by gw, after deducted by QMULTG, cw->gw
                                        !                                 &  -qwfrz(k)    & ! homo freezing: done later
                                        ) &
                                & +min(0.,pcc(k)) & ! evap of cw, v->cw
                        !                        &  -qiihr(k)  & ! qiihr: Hobbs-Rangno ice enhancement, not used in nssl either
                        !                        &  -qhlacw(k) & ! hail not used
                        )
                ! conservation of droplet charge
                if (abs(psccwmd(k)).gt.0.0 ) then
                    if (sign(1.0,psccwmd(k)).eq.sign(1.0,scx(k,lc))) then
                        write(0,*) 'psccwmd scx same sign', psccwmd(k), scx(k,lc)
                        write(0,*) 'i,j,k,t=', i,j,k,itimestep
                        write(0,*) 'fsccw(k),scx(k,lc),il5(k)',fsccw(k),scx(k,lc),il5(k)
                        write(0,*) 'stopping for psccwmd scxc same sign'
                        stop
                    end if
                    if (abs(psccwmd(k)*dt).gt.abs(scx(k,lc))) then
                        write(0,*) 'conserving psccw', psccwmd(k), dt, psccwmd(k)*dt, scx(k,lc)
                        write(0,*) 'fsccw(k)',fsccw(k)
                        write(0,*) 'i,j,k,t=', i,j,k,itimestep
                        dumtmp = -PRA(k) & ! PRA: equate to qracw, cw collected by rw, cw->rw
                                &  -PRC(k) & ! PRC: equate to qrcnw, cw autoconvert to rw, cw->rw
                                &  +il5(k)*(-PSACWI(k)   & ! PSACWI: equate to qiacw: collection of cw by ci
                                        &  -mnuccc(k)    & ! mnuccc: equate to qwctfzc, here includs contact and immersion frz
                                        &  -PGSACW(k)  & ! PGSACW: cw collected by sw to produce gw, cw->gw
                                        &  -qmults(k)-qmultg(k) & ! h-m ice multip. on sw & gw, cw->ci
                                        &  -PSACWS(k)  & ! PSACWS: equate to qsacw, cw collected by sw, after deducted by PGSACW and QMULTS, cw->sw
                                        &  -PSACWG(k)  & ! PSACWG: equate to qhacw, cw collected by gw, after deducted by QMULTG, cw->gw
                                        !                                 &  -qwfrz(k)    & ! homo freezing: done later
                                        ) &
                                & +min(0.,pcc(k))
                        write(0,*) 'dumtmp',dumtmp
                        if (an1d(k,lc)+dumtmp*dt.lt.0.0) then
                            write(0,*) 'neg qc in psc conserve',an1d(k,lc)+dumtmp*dt
                            write(0,*) 'an1d(k,lc),dumtmp*dt',an1d(k,lc),dumtmp*dt
                            !                            write(0,*) 'MNUCCR(k),pre(k)',MNUCCR(k),pre(k)
                        else if (an1d(k,lc)+dumtmp*dt.gt.0.0) then
                            write(0,*) 'pos qc in psccw conserve',an1d(k,lc)+dumtmp*dt
                            write(0,*) 'an1d(k,lc),dumtmp*dt',an1d(k,lc),dumtmp*dt
                        else
                            write(0,*) 'zero qc+dumtmp*dt',an1d(k,lc)+dumtmp*dt
                        end if
                        psccwmd(k) = -scx(k,lc)*dtpinv
                        write(0,*) 'conserved psccw', psccwmd(k), dt, psccwmd(k)*dt, scx(k,lc)
                    end if
                end if

                !
                !  Snow
                !
                pscswmi(k) =    il5(k)*(fscci(k)*PRCI(k) & ! PRCI: equate to qscni, autoconversion, ci->sw
                        &  +fscci(k)*PRAI(k)   & ! PRAI: equate to qsaci, ci collected by sw, ci->sw
                        !                                &  +fscrw(k)*ifrzs*qrfrzs(k)   & ! rw frozen to form sw, not used in morr since all rw frozen to gw
                        &  +fscci(k)*PRACIS(k) & ! PRACIS: equate to qracis, ci collected by rw to produce sw, ci->sw
                        &  +fsccw(k)*PSACWS(k) & ! PSACWS: equate to qsacw, cw collected by sw, after deducted by PGSACW and QMULTS, cw->sw
                        &  +fscrw(k)*PIACRS(k) & ! PIACRS: equate to qracis, rw collected by ci to produce sw, rw->sw
                        &  +fscrw(k)*PRACS(k)  & ! PRACS: equate to qsacr, rw collected by sw, rw->sw
                        )
                pscswmd(k) =    fscsw(k)*(&
                        &  -il5(k)*PSACR(k) & ! PSACR: sw collected by rw to produce gw, sw->gw
                                !                        &  - qhacs(k) & ! qhacs: sw collected by gw, not used in morr
                                !                        &  - qhlacs(k)  & ! hail not included
                                &  + (1-il5(k))*psmlt(k)*mixedphasefac & ! psmlt: equate to qsmlr, negative, sw melting, rw->sw
                                &  + (1-il5(k))*min(0.,evpms(k)) & ! evpms: evaporation of melted snow, negative, v-sw
                                &  + il5(k)*min(0.,eprds(k)) & ! eprds: sublimation snow, negative, v-sw
                        !                        &  - qsmul(k)  & ! h-m multiplication, but qs not involved in morr
                        !                        &  + min(0.0,qscev(k)) & ! =0 in nssl anyway
                        )
                ! conservation of snow charge
                if (abs(pscswmd(k)).gt.0.0 ) then
                    if (sign(1.0,pscswmd(k)).eq.sign(1.0,scx(k,ls))) then
                        write(0,*) 'pscswmd scx same sign', pscswmd(k), scx(k,ls)
                        write(0,*) 'i,j,k,t=', i,j,k,itimestep
                        write(0,*) 'fscsw(k),scx(k,ls),il5(k)',fscsw(k),scx(k,ls),il5(k)
                        write(0,*) 'stopping for pscswmd scxs same sign'
                        stop
                    end if
                    if (abs(pscswmd(k)*dt).gt.abs(scx(k,ls))) then
                        write(0,*) 'conserving pscsw', pscswmd(k), dt, pscswmd(k)*dt, scx(k,ls)
                        write(0,*) 'fscsw(k)',fscsw(k)
                        write(0,*) 'i,j,k,t=', i,j,k,itimestep
                        dumtmp = -il5(k)*PSACR(k) & ! PSACR: sw collected by rw to produce gw, sw->gw
                                &  + (1-il5(k))*psmlt(k)*mixedphasefac & ! psmlt: equate to qsmlr, negative, sw melting, rw->sw
                                &  + (1-il5(k))*min(0.,evpms(k)) & ! evpms: evaporation of melted snow, negative, v-sw
                                &  + il5(k)*min(0.,eprds(k))
                        write(0,*) 'dumtmp',dumtmp
                        if (an1d(k,ls)+dumtmp*dt.lt.0.0) then
                            write(0,*) 'neg qs in psc conserve',an1d(k,ls)+dumtmp*dt
                            write(0,*) 'an1d(k,ls),dumtmp*dt',an1d(k,ls),dumtmp*dt
                        else if (an1d(k,ls)+dumtmp*dt.gt.0.0) then
                            write(0,*) 'pos qs in pscsw conserve',an1d(k,ls)+dumtmp*dt
                            write(0,*) 'an1d(k,ls),dumtmp*dt',an1d(k,ls),dumtmp*dt
                        else
                            write(0,*) 'zero qs+dumtmp*dt',an1d(k,ls)+dumtmp*dt
                        end if
                        pscswmd(k) = -scx(k,ls)*dtpinv
                        write(0,*) 'conserved pscsw', pscswmd(k), dt, pscswmd(k)*dt, scx(k,ls)
                    end if
                end if

                !                fscsw2(k) = 0.0
                !                tmp =  qx(k,ls) + dtp*(pqswi(k) + pqswd(k) - qsshr(k)) ! should this add back the shed rain mass? probably!
                !                IF ( tmp > qxmin(ls) .and. qsshr(k) < 0.0 ) THEN
                !                    chgtmp = scx(k,ls) + dtp*(pscswmi(k) - pscswmd(k))
                !                    fscsw2(k) = chgtmp/tmp
                !                    pscswmd(k) =  pscswmd(k) + fscsw2(k)*qsshr(k)
                !                ENDIF

                !
                !  Graupel
                !
                pschwmi(k) =    il5(k)*(fscci(k)*PRACI(k) & ! PRACI: equate to qracif, ci collide w/ rw to produce gw, ci->gw
                        &  +fscrw(k)*PIACR(k) & ! PIACR: equate to qiacrf, rw collide w/ ci to produce gw, rw->gw
                        &  +fscrw(k)*MNUCCR(k) & ! MNUCCR: equate to qrfrz, rw frozen into gw, rw->gw
                        &  +fsccw(k)*PGSACW(k) & ! PGSACW: cw collected by sw to produce gw, cw->gw
                        &  +fsccw(k)*PSACWG(k) & ! PSACWG: equate to qhacw, cw collected by gw, cw->gw
                        &  +fscrw(k)*PGRACS(k) & ! PGRACS: rw collected by sw to produce gw, rw->gw
                        &  +fscsw(k)*PSACR(k) & ! PSACR: sw collected by rw to produce gw, sw->gw
                        &  +fscrw(k)*PRACG(k) & ! PRACG: equate to qhacr, rw collected by gw, rw->gw
                        )
                !                        &  +fscsw(k)*qhacs(k) & ! sw collected by gw, not included in morr
                !                        &  +fscci(k)*qhaci(k) & ! sw collected by ci, not included in morr
                pschwmd(k) =    fschw(k)*(&
                        !                        &  - qhlcnh(k)   & ! hw converted to hl, hail not used
                        &  +(1-il5(k))*pgmlt(k)*mixedphasefac   & ! pgmlt: equate to qhmlr, negative, gw melting, rw->gw
                                &  +(1-il5(k))*min(0.,evpmg(k))   & ! evpmg: negative, evap of melted graup, v->gw
                                &  +il5(k)*min(0.,eprdg(k))   & ! eprdg: negative, sublimation of graup, v->gw
                        !                        &  -qhmul1(k) ! h-m multiplication, but qg not involved in morr
                        )
                ! conservation of graupel charge
                if (abs(pschwmd(k)).gt.0.0 ) then
                    if (sign(1.0,pschwmd(k)).eq.sign(1.0,scx(k,lh))) then
                        write(0,*) 'pschwmd scx same sign', pschwmd(k), scx(k,lh)
                        write(0,*) 'i,j,k,t=', i,j,k,itimestep
                        write(0,*) 'fschw(k),scx(k,lh),il5(k)',fschw(k),scx(k,lh),il5(k)
                        write(0,*) 'stopping for pschwmd scxh same sign'
                        stop
                    end if
                    if (abs(pschwmd(k)*dt).gt.abs(scx(k,lh))) then
                        write(0,*) 'conserving pschw', pschwmd(k), dt, pschwmd(k)*dt, scx(k,lh)
                        write(0,*) 'fschw(k)',fschw(k)
                        write(0,*) 'i,j,k,t=', i,j,k,itimestep
                        dumtmp = (1-il5(k))*pgmlt(k)*mixedphasefac   & ! pgmlt: equate to qhmlr, negative, gw melting, rw->gw
                                &  +(1-il5(k))*min(0.,evpmg(k))   & ! evpmg: negative, evap of melted graup, v->gw
                                &  +il5(k)*min(0.,eprdg(k))
                        write(0,*) 'dumtmp',dumtmp
                        if (an1d(k,lh)+dumtmp*dt.lt.0.0) then
                            write(0,*) 'neg qh in psc conserve',an1d(k,lh)+dumtmp*dt
                            write(0,*) 'an1d(k,lh),dumtmp*dt',an1d(k,lh),dumtmp*dt
                            !                            write(0,*) 'MNUCCR(k),pre(k)',MNUCCR(k),pre(k)
                        else if (an1d(k,lh)+dumtmp*dt.gt.0.0) then
                            write(0,*) 'pos qh in pschw conserve',an1d(k,lh)+dumtmp*dt
                            write(0,*) 'an1d(k,lh),dumtmp*dt',an1d(k,lh),dumtmp*dt
                        else
                            write(0,*) 'zero qh+dumtmp*dt',an1d(k,lh)+dumtmp*dt
                        end if
                        pschwmd(k) = -scx(k,lh)*dtpinv
                        write(0,*) 'conserved pschw', pschwmd(k), dt, pschwmd(k)*dt, scx(k,lh)
                    end if
                end if

                !                fschw2(k) = 0.0
                !                tmp =  qx(k,lh) + dtp*(pqhwi(k) + pqhwd(k) - qhshr(k)) ! should this add back the shed rain mass? probably!
                !                IF ( tmp > qxmin(lh) .and. qhshr(k) < 0.0 ) THEN
                !                    chgtmp = scx(k,lh) + dtp*(pschwmi(k) - pschwmd(k))
                !                    fschw2(k) = chgtmp/tmp
                !                    pschwmd(k) =  pschwmd(k) + fschw2(k)*qhshr(k)
                !                ENDIF

                !
                !  Rain
                !
                pscrwmi(k) =    &
                        &  fsccw(k)*PRA(k) & ! PRA: equate to qracw, cw collected by rw, cw->rw
                        &  +fsccw(k)*PRC(k) & ! PRC: equate to qrcnw, cw autoconvert to rw, cw->rw
                        &  +(1-il5(k))*mixedphasefac*(-fscsw(k)*psmlt(k)  & ! psmlt: equate to qsmlr, negative, sw melting, rw->sw
                                &  -fschw(k)*pgmlt(k)  & ! pgmlt: equate to qhmlr, negative, gw melting, rw->gw
                                !                                                   &  -fschl(k)*qhlmlr(k) & ! hail not used
                                )
                pscrwmd(k) =     fscrw(k)*(&
                        &  il5(k)*(-PIACR(k) & ! PIACR: equate to qiacr, rw collected by ci to produce gw, rw->gw
                                &  -MNUCCR(k) & ! MNUCCR: equate to qrfrz, rw frozen into gw, rw->gw
                                &  -qmultr(k)-qmultrg(k) & ! h-m ice multip. on sw & gw, rw->ci
                                &  -PIACRS(k) & ! PIACRS: equate to qracis, rw collected by ci to produce sw, rw->sw
                                &  -PGRACS(k) & ! PGRACS: rw collected by sw to produce gw, rw->gw
                                &  -PRACS(k) & ! PRACS: equate to qsacr, rw collected by sw, rw->sw, =0 for warm
                                &  -PRACG(k) & ! PRACG: equate to qhacr, rw collected by gw, rw->gw, =0 for warm
                                ) &
                                !                        &  -qhlacr(k)  & ! hail not used
                                &  +min(0.0,pre(k)) & ! pre: equate to qrcev, evap of rw, negative only, rw->qv
                        )
                ! conservation of rain charge
                if (abs(pscrwmd(k)).gt.0.0 ) then
                    if (sign(1.0,pscrwmd(k)).eq.sign(1.0,scx(k,lr))) then
                        write(0,*) 'pscrwmd scx same sign', pscrwmd(k), scx(k,lr)
                        write(0,*) 'i,j,k,t=', i,j,k,itimestep
                        write(0,*) 'fscrw(k),scx(k,lr),il5(k)',fscrw(k),scx(k,lr),il5(k)
                        write(0,*) 'stopping for pscrwmd scxr same sign'
                        stop
                    end if
                    if (abs(pscrwmd(k)*dt).gt.abs(scx(k,lr))) then
                        write(0,*) 'conserving pscrw', pscrwmd(k), dt, pscrwmd(k)*dt, scx(k,lr)
                        write(0,*) 'fscrw(k)',fscrw(k)
                        write(0,*) 'i,j,k,t=', i,j,k,itimestep
                        dumtmp =  il5(k)*(-PIACR(k) & ! PIACR: equate to qiacr, rw collected by ci to produce gw, rw->gw
                                &  -MNUCCR(k) & ! MNUCCR: equate to qrfrz, rw frozen into gw, rw->gw
                                &  -qmultr(k)-qmultrg(k) & ! h-m ice multip. on sw & gw, rw->ci
                                &  -PIACRS(k) & ! PIACRS: equate to qracis, rw collected by ci to produce sw, rw->sw
                                &  -PGRACS(k) & ! PGRACS: rw collected by sw to produce gw, rw->gw
                                ) &
                                &  -PRACS(k) & ! PRACS: equate to qsacr, rw collected by sw, rw->sw
                                &  -PRACG(k) & ! PRACG: equate to qhacr, rw collected by gw, rw->gw
                                !                        &  -qhlacr(k)  & ! hail not used
                                &  +min(0.0,pre(k))  ! pre: equate to qrcev, evap of rw, negative only, rw->qv
                        write(0,*) 'dumtmp',dumtmp
                        if (an1d(k,lr)+dumtmp*dt.lt.0.0) then
                            write(0,*) 'neg qr in psc conserve',an1d(k,lr)+dumtmp*dt
                            write(0,*) 'an1d(k,lr),dumtmp*dt',an1d(k,lr),dumtmp*dt
                            write(0,*) 'MNUCCR(k),pre(k)',MNUCCR(k),pre(k)
                        else if (an1d(k,lr)+dumtmp*dt.gt.0.0) then
                            write(0,*) 'pos qr in pscrw conserve',an1d(k,lr)+dumtmp*dt
                            write(0,*) 'an1d(k,lr),dumtmp*dt',an1d(k,lr),dumtmp*dt
                            write(0,*) 'MNUCCR(k),pre(k)',MNUCCR(k),pre(k)
                        else
                            write(0,*) 'zero qr+dumtmp*dt',an1d(k,lr)+dumtmp*dt
                        end if

                        pscrwmd(k) = -scx(k,lr)*dtpinv
                        write(0,*) 'conserved pscrw', pscrwmd(k), dt, pscrwmd(k)*dt, scx(k,lr)
                    end if
                end if

            end if

            !
            !
            ! Collisional Charging transfers (inductive/noninductive)
            !
            !
            if ( itest .eq. 1 ) then
                !
                !  cloud water
                !
                psccwi(k) =  psccwi(k)   &
                        & + max(-schacw(k),0.0)   &
                        !                        & + max(-schlacw(k),0.0)   &
                        & + max(-scsacw(k),0.0)
                psccwd(k) =  psccwd(k)   &
                        & + min(-schacw(k),0.0)   &
                        !                        & + min(-schlacw(k),0.0)   &
                        & + min(-scsacw(k),0.0)

                qcwtmp(k) = Max(0.0, qc3d(k) + dtp*qc3dten(k))
                IF ( qcwtmp(k) .gt. qxmin(lc) ) then
                    ccwtmp = Max(0.001, (nc3d(k) + dtp*nc3dten(k)) *rho(k) )
                    cpqc = (scx(k,lc) + dtp*(psccwi(k)+psccwd(k) + psccwmi(k)+psccwmd(k))) &
                            /ccwtmp
                    IF ( Abs(cpqc) .gt. scwppmx .and.   &
                            & (psccwi(k)+psccwd(k)) .ne. 0.0 .and.   &
                            & Sign(1.0, scx(k,lc)) .eq.   &
                                    & Sign(1.0, (psccwi(k)+psccwd(k)))) THEN
                        cpqc0 = Sign(scwppmx,cpqc)
                        scfac = (cpqc0*ccwtmp - scx(k,lc) -   &
                                &       dtp*(psccwmi(k)+psccwmd(k)))/   &
                                &              ( dtp*(psccwi(k)+psccwd(k)) )
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

                IF ( scfac .lt. 1.0 .and. QC3D(k) .gt. qxmin(lc) ) THEN
                    schacw(k) = scfac*schacw(k)
                    !                    schlacw(k) = scfac*schlacw(k)
                    scsacw(k) = scfac*scsacw(k)
                    psccwi(k) = 0.0
                    psccwd(k) = 0.0
                    psccwi(k) =  psccwi(k)   &
                            & + max(-schacw(k),0.0)   &
                            !                            & + max(-schlacw(k),0.0)   &
                            & + max(-scsacw(k),0.0)
                    psccwd(k) =  psccwd(k)   &
                            & + min(-schacw(k),0.0)   &
                            !                            & + min(-schlacw(k),0.0)   &
                            & + min(-scsacw(k),0.0)
                ENDIF

                !
                !  cloud ice
                !
                psccii(k) =  psccii(k)   &
                        & + max(-schaci(k),0.0)   &
                        !                    & + max(-schlaci(k),0.0)   &
                        & + max(-scsaci(k),0.0)
                psccid(k) =  psccid(k)   &
                        & + min(-schaci(k),0.0)   &
                        !                    & + min(-schlaci(k),0.0)   &
                        & + min(-scsaci(k),0.0)

                qcitmp = Max(0.0, qi3d(k) + dtp*QI3DTEN(k))
                IF ( qcitmp .gt. qxmin(li) ) then
                    ccitmp = Max(0.001, (NI3D(k) + dtp*NI3DTEN(k)) *rho(k) )
                    cpci = (scx(k,li) + dtp*(psccii(k)+psccid(k) + psccimi(k)+psccimd(k))) &
                            /ccitmp
                    IF ( Abs(cpci) .gt. scippmx .and.   &
                            & (psccii(k)+psccid(k) ) .ne. 0.0 .and.   &
                            & Sign(1.0,scx(k,li)) .eq.   &
                                    & Sign(1.0,(psccii(k)+psccid(k)))) THEN
                        cpci0 = Sign(scippmx,cpci)
                        scfac = (cpci0*ccitmp - scx(k,li) -    &
                                &       dtp*(psccimi(k)+psccimd(k)))/   &
                                &              ( dtp*(psccii(k)+psccid(k)) )
                        IF ( scfac .gt. 1.0 ) THEN
                            write(0,*) 'OUCH: sccifac = ',scfac,   &
                                    & dtp*(psccimi(k)+psccimd(k)),cpci,scx(k,li),ccitmp
                            write(0,*) 'OUCH: sccifac = ',scfac,   &
                                    & dtp*(psccimi(k)+psccimd(k)),cpci,scx(k,li),ccitmp,   &
                                    & k,t3d(k)
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
                    schaci(k) = scfac*schaci(k)
                    !                    schlaci(k) = scfac*schlaci(k)
                    scsaci(k) = scfac*scsaci(k)
                    psccii(k) = 0.0
                    psccid(k) = 0.0
                    psccii(k) =  psccii(k)   &
                            & + max(-schaci(k),0.0)   &
                            !                            & + max(-schlaci(k),0.0)   &
                            & + max(-scsaci(k),0.0)
                    psccid(k) =  psccid(k)   &
                            & + min(-schaci(k),0.0)   &
                            !                            & + min(-schlaci(k),0.0)   &
                            & + min(-scsaci(k),0.0)
                ENDIF

                !
                !  rain
                !
                pscrwi(k) =  pscrwi(k)  &
                        !                        & + max(-schacr(k),0.0)    &
                        & + max(-scsacr(k),0.0)
                !                        & + max(scsshr(k),0.0)     & !Canosa & List 93
                !                        & + max(schshr(k),0.0)       !Canosa & List 93

                pscrwd(k) =  pscrwd(k)  &
                        !                        & + min(-schacr(k),0.0)    &
                        & + min(-scsacr(k),0.0)
                !                        & + min(scsshr(k),0.0)     & !Canosa & List 93
                !                        & + min(schshr(k),0.0)       !Canosa & List 93

                !
                !  snow
                !
                pscswi(k) =  pscswi(k)   &
                        & + max(scsaci(k),0.0)   &
                        & + max(scsacw(k),0.0)   &
                        & + max(scsacr(k),0.0)   &
                        & + max(-schacs(k),0.0)
                !                        & + max(-schlacs(k),0.0)  &
                !                        & + max(scsmlr(k),0.0)    &  !Drake
                !                        & + max(scsdep(k),0.0)    &  !Dong&Hallett
                !                        & - min(scsshr(k),0.0)       !Canosa & List 93
                pscswd(k) =  pscswd(k)   &
                        & + min(scsaci(k),0.0)   &
                        & + min(scsacw(k),0.0)   &
                        & + min(scsacr(k),0.0)   &
                        & + min(-schacs(k),0.0)
                !                        & + min(-schlacs(k),0.0)  &
                !     > + min(-scfacs(k),0.0)
                !                        & + min(scsmlr(k),0.0)    &   !Drake, but always 0.? (since scsmlr >= 0)
                !                        & + min(scsdep(k),0.0)    &   !Dong&Hallett
                !                        & - max(scsshr(k),0.0)        !Canosa & List 93

                !
                !  graupel
                !
                pschwi(k) =  pschwi(k)   &
                        & + max(schaci(k),0.0)   &
                        & + max(schacs(k),0.0)   &
                        & + max(schacw(k),0.0)
                !                        & + max(schacr(k),0.0)   &
                !                        & + max(schmlr(k),0.0)   & !Drake
                !                        & - min(schshr(k),0.0)     !Canosa & List 93
                pschwd(k) =  pschwd(k)   &
                        & + min(schaci(k),0.0)   &
                        & + min(schacs(k),0.0)   &
                        & + min(schacw(k),0.0)
                !                        & + min(schacr(k),0.0)   &
                !                        & + min(schmlr(k),0.0)   & !Drake, but always 0.? (since schmlr >= 0)
                !                        & - max(schshr(k),0.0)     !Canosa & List 93

                !                !  ions not changed as in micro suggested by Ted
                !                !
                !                pscnii(k) =  pscnii(k)  &
                !                        &               - min(fscrw(k),0.0)*min(qrcev(k),0.0)  &
                !                        &               - min(fscsw(k),0.0)*min(qscev(k),0.0)  &
                !                        &               + min(-schmlr(k),0.0)    &
                !                        &               + min(-scsmlr(k),0.0)
                !                pscpii(k) =  pscpii(k)  &
                !                        &               - max(fscrw(k),0.0)*min(qrcev(k),0.0)  &
                !                        &               - max(fscsw(k),0.0)*min(qscev(k),0.0)  &
                !                        &               + max(-schmlr(k),0.0)    &
                !                        &               + max(-scsmlr(k),0.0)
                !                pscnid(k) =  pscnid(k)  &
                !                        &               - max(-scsdep(k),0.0)
                !                pscpid(k) =  pscpid(k)  &
                !                        &               - min(-scsdep(k),0.0)

            end if ! itest




            !            !
            !            ! EVERYTHING ELSE
            !            !

            scworig = scx(k,lc)
            ! update scx array (complete)
            scx(k,lc) = scx(k,lc) +   &
                    &   dt*(psccwi(k)+psccwd(k)+psccwmi(k)+psccwmd(k)+pscc_sed(k))
            scx(k,li) = scx(k,li) +   &
                    &   dt*(psccii(k)+psccid(k)+psccimi(k)+psccimd(k)+psci_sed(k))
            scx(k,lr) = scx(k,lr) +   &
                    &   dt*(pscrwi(k)+pscrwd(k)+pscrwmi(k)+pscrwmd(k)+pscr_sed(k))
            scx(k,ls) = scx(k,ls) +   &
                    &   dt*(pscswi(k)+pscswd(k)+pscswmi(k)+pscswmd(k)+pscs_sed(k))
            scx(k,lh) = scx(k,lh) +   &
                    &   dt*(pschwi(k)+pschwd(k)+pschwmi(k)+pschwmd(k)+pscg_sed(k))


            cionp(k) = cionp(k) +  &
                    &   dt*ec_i*(pscpii(k)+pscpid(k)) ! #/m3 * C / C/m3. ! each ion has electron charge
            cionn(k) = cionn(k) -  &
                    &   dt*ec_i*(pscnii(k)+pscnid(k))
            if ( largeion ) then
                clionp(k) = clionp(k) +  &
                        &   dt*ec_i*(pscplii(k)+pscplid(k))
                clionn(k) = clionn(k) -  &
                        &   dt*ec_i*(pscnlii(k)+pscnlid(k))
            endif

            ! put scx arrays back into an array
            an1d(k,lscpi) = cionp(k)
            an1d(k,lscni) = cionn(k)


            DO il = lc,lhab
                tmp = an1d(k,lsc(il))
                an1d(k,lsc(il)) = scx(k,il)
                IF (.not.( Abs (scx(k,il) ) < 1000000.e-9  )) THEN
                    write(0,*) 'Problem0a with scx il = ',il,lsc(il),scx(k,il),an1d(k,il),an1d(k,ln(il))
                    write(0,*) 'Problem2 with scx il = ',il,scx(k,il),tmp,k,tmp-scx(k,il),(tmp-scx(k,il))*dtpinv
                    write(0,*) 'temper = ',T3D(k)
                    write(0,*) 'i,j,k,t=',i,j,k,itimestep
                    write(0,*) 'qx,cx = ',an1d(k,il),an1d(k,ln(il))
                    write(0,*)  psccwi(k), psccwd(k)
                    write(0,*)  psccii(k), psccid(k)
                    write(0,*)  pscrwi(k), pscrwd(k)
                    write(0,*)  pscswi(k), pscswd(k)
                    write(0,*)  pschwi(k), pschwd(k)
                    !                                        write(0,*)  pschli(k), pschld(k)
                    write(0,*)  pscpii(k) ,pscpid(k)
                    write(0,*)  pscnii(k) ,pscnid(k)
                    write(0,*)
                    write(0,*)  psccwmi(k), psccwmd(k)
                    write(0,*)  psccimi(k), psccimd(k)
                    write(0,*)  pscrwmi(k), pscrwmd(k)
                    write(0,*)  pscswmi(k), pscswmd(k)
                    write(0,*)  pschwmi(k), pschwmd(k)
                    !                                        write(0,*)  pschlmi(k), pschlmd(k)
                    write(0,*)
                    write(0,*)  pscc_sed(k)
                    write(0,*)  psci_sed(k)
                    write(0,*)  pscr_sed(k)
                    write(0,*)  pscs_sed(k)
                    write(0,*)  pscg_sed(k)

                ENDIF
            ENDDO
            do iltest=lc,lh
                if (isnan(scx(k,iltest))) then
                    write (0,*) 'scx=nan, iltest=',iltest
                    write (0,*) 'i,j,k,t=',i,j,k,itimestep
                    write (0,*) 'scx(k,iltest)=',scx(k,iltest)
                    write (0,*) 'psccwi(k) psccwd(k) psccwmi(k) psccwmd(k) pscc_sed(k)'
                    write (0,*) psccwi(k), psccwd(k), psccwmi(k), psccwmd(k), pscc_sed(k)
                    write (0,*) 'psccii(k) psccid(k) psccimi(k) psccimd(k) psci_sed(k)'
                    write (0,*) psccii(k), psccid(k), psccimi(k), psccimd(k), psci_sed(k)
                    write (0,*) 'pscrwi(k) pscrwd(k) pscrwmi(k) pscrwmd(k) pscr_sed(k)'
                    write (0,*) pscrwi(k), pscrwd(k), pscrwmi(k), pscrwmd(k), pscr_sed(k)
                    write (0,*) 'pscswi(k) pscswd(k) pscswmi(k) pscswmd(k) pscs_sed(k)'
                    write (0,*) pscswi(k), pscswd(k), pscswmi(k), pscswmd(k), pscs_sed(k)
                    write (0,*) 'pschwi(k) pschwd(k) pschwmi(k) pschwmd(k) pscg_sed(k)'
                    write (0,*) pschwi(k), pschwd(k), pschwmi(k), pschwmd(k), pscg_sed(k)
                    write (0,*) 'stopping'
                    print*, 'stopping'
                    STOP
                end if
            end do

            ! PUT ALL CLOUD ICE IN SNOW CATEGORY IF MEAN DIAMETER EXCEEDS 2 * dcs
            ! To be consistent with same thing done earlier fro q and n
            IF (QI3D(K).GE.QSMALL.AND.T3D(K).LT.273.15.AND.LAMI(K).GE.1.E-10) THEN
                IF (1./LAMI(K).GE.2.*DCS) THEN
                    an1d(k,lsc(ls)) = an1d(k,lsc(ls)) + an1d(k,lsc(li))
                    an1d(k,lsc(li)) = 0.0
                END IF
            END IF

            !
            ! Load the save arrays
            !
            !            t1(k) = scsaci(k)
            t2(k) = schacs(k)
            t3(k) = schaci(k)
            !            t4(k) = scsacw(k)
            t5(k) = schacw(k)
            t6(k) = schlacs(k) ! =0
            t8(k) = schlaci(k) ! =0


        else
            scx(k,lsc(lc):lsc(lhab)) = 0.0
            cionp(k) = 0.0
            cionn(k) = 0.0
            CSACI0(k) = 0.0
            CHACI0(k) = 0.0
            CHACS0(k) = 0.0
            an1d(k,lsc(lc):lsc(lhab)) = 0.0
            umc1d(k) = 0.0
            umr1d(k) = 0.0
            umi1d(k) = 0.0
            ums1d(k) = 0.0
            umg1d(k) = 0.0

        END IF ! ipelec

        ! end charging
        ! ⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡⚡


        ! hm add tendencies here, then call sizeparameter
        ! to ensure consisitency between mixing ratio and number concentration

        QC3D(k)        = QC3D(k)+QC3DTEN(k)*DT
        QI3D(k)        = QI3D(k)+QI3DTEN(k)*DT
        QNI3D(k)       = QNI3D(k)+QNI3DTEN(k)*DT
        QR3D(k)        = QR3D(k)+QR3DTEN(k)*DT
        NC3D(k)        = NC3D(k)+NC3DTEN(k)*DT
        NI3D(k)        = NI3D(k)+NI3DTEN(k)*DT
        NS3D(k)        = NS3D(k)+NS3DTEN(k)*DT
        NR3D(k)        = NR3D(k)+NR3DTEN(k)*DT


        IF (IGRAUP.EQ.0) THEN
            QG3D(k)        = QG3D(k)+QG3DTEN(k)*DT
            NG3D(k)        = NG3D(k)+NG3DTEN(k)*DT
        END IF

        ! ADD TEMPERATURE AND WATER VAPOR TENDENCIES FROM MICROPHYSICS
        T3D(K)         = T3D(K)+T3DTEN(k)*DT
        QV3D(K)        = QV3D(K)+QV3DTEN(k)*DT

        ! yc added 2/24/2022
        ! get rid of artificial neg q/n, prevent adding to qv
        !        qc3d(k) = max(qc3d(k),0.0)
        !        qi3d(k) = max(qi3d(k),0.0)
        !        qni3d(k) = max(qni3d(k),0.0)
        !        qr3d(k) = max(qr3d(k),0.0)
        !        nc3d(k) = max(nc3d(k),0.0)
        !        ni3d(k) = max(ni3d(k),0.0)
        !        ns3d(k) = max(ns3d(k),0.0)
        !        nr3d(k) = max(nr3d(k),0.0)
        !        qg3d(k) = max(qg3d(k),0.0)
        !        ng3d(k) = max(ng3d(k),0.0)
        !        qv3d(k) = max(qv3d(k),0.0)
        ! SATURATION VAPOR PRESSURE AND MIXING RATIO

        ! hm, add fix for low pressure, 5/12/10
        EVS(K) = min(0.99*pres(k),POLYSVP(T3D(K),0))   ! PA
        EIS(K) = min(0.99*pres(k),POLYSVP(T3D(K),1))   ! PA

        ! MAKE SURE ICE SATURATION DOESN'T EXCEED WATER SAT. NEAR FREEZING

        IF (EIS(K).GT.EVS(K)) EIS(K) = EVS(K)

        QVS(K) = EP_2*EVS(K)/(PRES(K)-EVS(K))
        QVI(K) = EP_2*EIS(K)/(PRES(K)-EIS(K))

        QVQVS(K) = QV3D(K)/QVS(K)
        QVQVSI(K) = QV3D(K)/QVI(K)

        ! AT SUBSATURATION, REMOVE SMALL AMOUNTS OF CLOUD/PRECIP WATER
        ! hm 7/9/09 change limit to 1.e-8

        IF (QVQVS(K).LT.0.9) THEN
            IF (QR3D(K).LT.1.E-8) THEN
                if(ipelec>0 .and.qr3d(k)<0.0) then
                    write(0,*) '-------------------------------------------'
                    write(0,*) 'saturation adjust w/ neg qr: '
                    write(0,*) 'i,j,k,t', i,j,k,itimestep
                    write(0,*) 'qr3d(k)',qr3d(k)
                    write(0,*) 'an1d(k,lscr)', an1d(k,lscr)
                    !                    write(0,*)
                end if
                QV3D(K)=QV3D(K)+QR3D(K)
                T3D(K)=T3D(K)-QR3D(K)*XXLV(K)/CPM(K)
                QR3D(K)=0.
                IF ( an1d(k,lscr) .gt. 0.0 ) THEN
                    an1d(k,lscpi) = an1d(k,lscpi) + an1d(k,lscr)*ec_i
                ELSE
                    an1d(k,lscni) = an1d(k,lscni) - an1d(k,lscr)*ec_i
                ENDIF
                an1d(k,lscr) = 0.0
            END IF
            IF (QC3D(K).LT.1.E-8) THEN
                QV3D(K)=QV3D(K)+QC3D(K)
                T3D(K)=T3D(K)-QC3D(K)*XXLV(K)/CPM(K)
                QC3D(K)=0.
                IF ( an1d(k,lscw) .gt. 0.0 ) THEN
                    an1d(k,lscpi) = an1d(k,lscpi) + an1d(k,lscw)*ec_i
                ELSE
                    an1d(k,lscni) = an1d(k,lscni) - an1d(k,lscw)*ec_i
                ENDIF
                an1d(k,lscw) = 0.0
            END IF
        END IF

        IF (QVQVSI(K).LT.0.9) THEN
            IF (QI3D(K).LT.1.E-8) THEN
                QV3D(K)=QV3D(K)+QI3D(K)
                T3D(K)=T3D(K)-QI3D(K)*XXLS(K)/CPM(K)
                QI3D(K)=0.
                IF ( an1d(k,lsci) .gt. 0.0 ) THEN
                    an1d(k,lscpi) = an1d(k,lscpi) + an1d(k,lsci)*ec_i
                ELSE
                    an1d(k,lscni) = an1d(k,lscni) - an1d(k,lsci)*ec_i
                ENDIF
                an1d(k,lsci) = 0.0
            END IF
            IF (QNI3D(K).LT.1.E-8) THEN
                QV3D(K)=QV3D(K)+QNI3D(K)
                T3D(K)=T3D(K)-QNI3D(K)*XXLS(K)/CPM(K)
                QNI3D(K)=0.
                IF ( an1d(k,lscs) .gt. 0.0 ) THEN
                    an1d(k,lscpi) = an1d(k,lscpi) + an1d(k,lscs)*ec_i
                ELSE
                    an1d(k,lscni) = an1d(k,lscni) - an1d(k,lscs)*ec_i
                ENDIF
                an1d(k,lscs) = 0.0
            END IF
            IF (QG3D(K).LT.1.E-8) THEN
                QV3D(K)=QV3D(K)+QG3D(K)
                T3D(K)=T3D(K)-QG3D(K)*XXLS(K)/CPM(K)
                QG3D(K)=0.
                IF ( an1d(k,lsch) .gt. 0.0 ) THEN
                    an1d(k,lscpi) = an1d(k,lscpi) + an1d(k,lsch)*ec_i
                ELSE
                    an1d(k,lscni) = an1d(k,lscni) - an1d(k,lsch)*ec_i
                ENDIF
                an1d(k,lsch) = 0.0
            END IF
        END IF

        !..................................................................
        ! IF MIXING RATIO < QSMALL SET MIXING RATIO AND NUMBER CONC TO ZERO

        IF (QC3D(K).LT.QSMALL) THEN
            QC3D(K) = 0.
            NC3D(K) = 0.
            EFFC(K) = 0.
            an1d(k,lscw) = 0.0
        END IF
        IF (QR3D(K).LT.QSMALL) THEN
            QR3D(K) = 0.
            NR3D(K) = 0.
            EFFR(K) = 0.
            an1d(k,lscr) = 0.0
        END IF
        IF (QI3D(K).LT.QSMALL) THEN
            QI3D(K) = 0.
            NI3D(K) = 0.
            EFFI(K) = 0.
            an1d(k,lsci) = 0.0
        END IF
        IF (QNI3D(K).LT.QSMALL) THEN
            QNI3D(K) = 0.
            NS3D(K) = 0.
            EFFS(K) = 0.
            an1d(k,lscs) = 0.0
        END IF
        IF (QG3D(K).LT.QSMALL) THEN
            QG3D(K) = 0.
            NG3D(K) = 0.
            EFFG(K) = 0.
            an1d(k,lsch) = 0.0
        END IF

        !..................................
        ! IF THERE IS NO CLOUD/PRECIP WATER, THEN SKIP CALCULATIONS

        IF (QC3D(K).LT.QSMALL.AND.QI3D(K).LT.QSMALL.AND.QNI3D(K).LT.QSMALL &
                .AND.QR3D(K).LT.QSMALL.AND.QG3D(K).LT.QSMALL) GOTO 500

        !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        ! CALCULATE INSTANTANEOUS PROCESSES

        ! ADD MELTING OF CLOUD ICE TO FORM RAIN

        IF (QI3D(K).GE.QSMALL.AND.T3D(K).GE.273.15) THEN
            QR3D(K) = QR3D(K)+QI3D(K)
            an1d(k,lscr)=an1d(k,lscr)+an1d(k,lsci)
            T3D(K) = T3D(K)-QI3D(K)*XLF(K)/CPM(K)
#if (WRF_CHEM == 1)
            tqimelt(K)=QI3D(K)/DT
#endif
            QI3D(K) = 0.
            an1d(k,lsci) = 0.
            NR3D(K) = NR3D(K)+NI3D(K)
            NI3D(K) = 0.
        END IF

        ! ****SENSITIVITY - NO ICE
        IF (ILIQ.EQ.1) GOTO 778

        ! HOMOGENEOUS FREEZING OF CLOUD WATER

        IF (T3D(K).LE.233.15.AND.QC3D(K).GE.QSMALL) THEN
            QI3D(K)=QI3D(K)+QC3D(K)
            an1d(k,lsci)=an1d(k,lsci)+an1d(k,lscw)
            T3D(K)=T3D(K)+QC3D(K)*XLF(K)/CPM(K)
            QC3D(K)=0.
            an1d(k,lscw)=0.
            NI3D(K)=NI3D(K)+NC3D(K)
            NC3D(K)=0.
        END IF

        ! HOMOGENEOUS FREEZING OF RAIN

        IF (IGRAUP.EQ.0) THEN

            IF (T3D(K).LE.233.15.AND.QR3D(K).GE.QSMALL) THEN
                QG3D(K) = QG3D(K)+QR3D(K)
                an1d(k,lsch)=an1d(k,lsch)+an1d(k,lscr)
                T3D(K) = T3D(K)+QR3D(K)*XLF(K)/CPM(K)
                QR3D(K) = 0.
                an1d(k,lscr) = 0.
                NG3D(K) = NG3D(K)+ NR3D(K)
                NR3D(K) = 0.
            END IF

        ELSE IF (IGRAUP.EQ.1) THEN

            IF (T3D(K).LE.233.15.AND.QR3D(K).GE.QSMALL) THEN
                QNI3D(K) = QNI3D(K)+QR3D(K)
                an1d(k,lscs)=an1d(k,lscs)+an1d(k,lscr)
                T3D(K) = T3D(K)+QR3D(K)*XLF(K)/CPM(K)
                QR3D(K) = 0.
                an1d(k,lscr) = 0.
                NS3D(K) = NS3D(K)+NR3D(K)
                NR3D(K) = 0.
            END IF

        END IF

        778    CONTINUE

        ! MAKE SURE NUMBER CONCENTRATIONS AREN'T NEGATIVE

        NI3D(K) = MAX(0.,NI3D(K))
        NS3D(K) = MAX(0.,NS3D(K))
        NC3D(K) = MAX(0.,NC3D(K))
        NR3D(K) = MAX(0.,NR3D(K))
        NG3D(K) = MAX(0.,NG3D(K))

        !......................................................................
        ! CLOUD ICE

        IF (QI3D(K).GE.QSMALL) THEN
            LAMI(K) = (CONS12*                 &
                    NI3D(K)/QI3D(K))**(1./DI)

            ! CHECK FOR SLOPE

            ! ADJUST VARS

            IF (LAMI(K).LT.LAMMINI) THEN

                LAMI(K) = LAMMINI

                N0I(K) = LAMI(K)**4*QI3D(K)/CONS12

                NI3D(K) = N0I(K)/LAMI(K)
            ELSE IF (LAMI(K).GT.LAMMAXI) THEN
                LAMI(K) = LAMMAXI
                N0I(K) = LAMI(K)**4*QI3D(K)/CONS12

                NI3D(K) = N0I(K)/LAMI(K)
            END IF
            xdia(k,li) = 1.0/LAMI(K)
        END IF

        !......................................................................
        ! RAIN

        IF (QR3D(K).GE.QSMALL) THEN
            LAMR(K) = (PI*RHOW*NR3D(K)/QR3D(K))**(1./3.)

            ! CHECK FOR SLOPE

            ! ADJUST VARS

            IF (LAMR(K).LT.LAMMINR) THEN

                LAMR(K) = LAMMINR

                N0RR(K) = LAMR(K)**4*QR3D(K)/(PI*RHOW)

                NR3D(K) = N0RR(K)/LAMR(K)
            ELSE IF (LAMR(K).GT.LAMMAXR) THEN
                LAMR(K) = LAMMAXR
                N0RR(K) = LAMR(K)**4*QR3D(K)/(PI*RHOW)

                NR3D(K) = N0RR(K)/LAMR(K)
            END IF

        END IF

        !......................................................................
        ! CLOUD DROPLETS

        ! MARTIN ET AL. (1994) FORMULA FOR PGAM

        IF (QC3D(K).GE.QSMALL) THEN

            DUM = PRES(K)/(287.15*T3D(K))
            PGAM(K)=0.0005714*(NC3D(K)/1.E6*DUM)+0.2714
            PGAM(K)=1./(PGAM(K)**2)-1.
            PGAM(K)=MAX(PGAM(K),2.)
            PGAM(K)=MIN(PGAM(K),10.)

            ! CALCULATE LAMC

            LAMC(K) = (CONS26*NC3D(K)*GAMMA(PGAM(K)+4.)/   &
                    (QC3D(K)*GAMMA(PGAM(K)+1.)))**(1./3.)

            ! LAMMIN, 60 MICRON DIAMETER
            ! LAMMAX, 1 MICRON

            LAMMIN = (PGAM(K)+1.)/60.E-6
            LAMMAX = (PGAM(K)+1.)/1.E-6

            IF (LAMC(K).LT.LAMMIN) THEN
                LAMC(K) = LAMMIN
                NC3D(K) = EXP(3.*LOG(LAMC(K))+LOG(QC3D(K))+              &
                        LOG(GAMMA(PGAM(K)+1.))-LOG(GAMMA(PGAM(K)+4.)))/CONS26

            ELSE IF (LAMC(K).GT.LAMMAX) THEN
                LAMC(K) = LAMMAX
                NC3D(K) = EXP(3.*LOG(LAMC(K))+LOG(QC3D(K))+              &
                        LOG(GAMMA(PGAM(K)+1.))-LOG(GAMMA(PGAM(K)+4.)))/CONS26

            END IF
            !                xdia(k,lc) = 1.0/LAMC(K)
!            xdia(k,lc) = (PGAM(K)+1.)/LAMC(K)
            xdia(k,lc) = (6/(PI*RHOW)*QC3D(k)/NC3D(k))**(1./3.)
            xdia(k,lc) = min(max(xdia(k,lc),1.E-6),60.E-6)

        END IF

        !......................................................................
        ! SNOW

        IF (QNI3D(K).GE.QSMALL) THEN
            LAMS(K) = (CONS1*NS3D(K)/QNI3D(K))**(1./DS)

            ! CHECK FOR SLOPE

            ! ADJUST VARS

            IF (LAMS(K).LT.LAMMINS) THEN
                LAMS(K) = LAMMINS
                N0S(K) = LAMS(K)**4*QNI3D(K)/CONS1

                NS3D(K) = N0S(K)/LAMS(K)

            ELSE IF (LAMS(K).GT.LAMMAXS) THEN

                LAMS(K) = LAMMAXS
                N0S(K) = LAMS(K)**4*QNI3D(K)/CONS1
                NS3D(K) = N0S(K)/LAMS(K)
            END IF
            xdia(k,ls) = 1.0/LAMS(K)

        END IF

        !......................................................................
        ! GRAUPEL

        IF (QG3D(K).GE.QSMALL) THEN
            LAMG(K) = (CONS2*NG3D(K)/QG3D(K))**(1./DG)

            ! CHECK FOR SLOPE

            ! ADJUST VARS

            IF (LAMG(K).LT.LAMMING) THEN
                LAMG(K) = LAMMING
                N0G(K) = LAMG(K)**4*QG3D(K)/CONS2

                NG3D(K) = N0G(K)/LAMG(K)

            ELSE IF (LAMG(K).GT.LAMMAXG) THEN

                LAMG(K) = LAMMAXG
                N0G(K) = LAMG(K)**4*QG3D(K)/CONS2

                NG3D(K) = N0G(K)/LAMG(K)
            END IF
            xdia(k,lh) = 1.0/LAMG(K)

        END IF

        500  CONTINUE

        ! CALCULATE EFFECTIVE RADIUS

        IF (QI3D(K).GE.QSMALL) THEN
            EFFI(K) = 3./LAMI(K)/2.*1.E6
        ELSE
            EFFI(K) = 25.
        END IF

        IF (QNI3D(K).GE.QSMALL) THEN
            EFFS(K) = 3./LAMS(K)/2.*1.E6
        ELSE
            EFFS(K) = 25.
        END IF

        IF (QR3D(K).GE.QSMALL) THEN
            EFFR(K) = 3./LAMR(K)/2.*1.E6
        ELSE
            EFFR(K) = 25.
        END IF

        IF (QC3D(K).GE.QSMALL) THEN
            EFFC(K) = GAMMA(PGAM(K)+4.)/                        &
                    GAMMA(PGAM(K)+3.)/LAMC(K)/2.*1.E6
        ELSE
            EFFC(K) = 25.
        END IF

        IF (QG3D(K).GE.QSMALL) THEN
            EFFG(K) = 3./LAMG(K)/2.*1.E6
        ELSE
            EFFG(K) = 25.
        END IF

        ! HM ADD 1/10/06, ADD UPPER BOUND ON ICE NUMBER, THIS IS NEEDED
        ! TO PREVENT VERY LARGE ICE NUMBER DUE TO HOMOGENEOUS FREEZING
        ! OF DROPLETS, ESPECIALLY WHEN INUM = 1, SET MAX AT 10 CM-3
        !          NI3D(K) = MIN(NI3D(K),10.E6/RHO(K))
        ! HM, 12/28/12, LOWER MAXIMUM ICE CONCENTRATION TO ADDRESS PROBLEM
        ! OF EXCESSIVE AND PERSISTENT ANVIL
        ! NOTE: THIS MAY CHANGE/REDUCE SENSITIVITY TO AEROSOL/CCN CONCENTRATION
        NI3D(K) = MIN(NI3D(K),0.3E6/RHO(K))

        ! ADD BOUND ON DROPLET NUMBER - CANNOT EXCEED AEROSOL CONCENTRATION
        IF (iinum.EQ.0.AND.IACT.EQ.2) THEN
            NC3D(K) = MIN(NC3D(K),(NANEW1+NANEW2)/RHO(K))
        END IF
        ! SWITCH FOR CONSTANT DROPLET NUMBER
        IF (iinum.EQ.1) THEN
            ! CHANGE NDCNST FROM CM-3 TO KG-1
            NC3D(K) = NDCNST*1.E6/RHO(K)
        END IF

    END DO !!! K LOOP

    400         CONTINUE

    IF (IPELEC>0) THEN
        !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        ! update an1d "q and n" part for use in SCIONA later
        an1d(kts:kte,lc) = qc3d(kts:kte)
        an1d(kts:kte,lr) = qr3d(kts:kte)
        an1d(kts:kte,li) = qi3d(kts:kte)
        an1d(kts:kte,ls) = qni3d(kts:kte)
        an1d(kts:kte,lh) = qg3d(kts:kte)

        an1d(kts:kte,lnc) = nc3d(kts:kte)*RHO(kts:kte)
        an1d(kts:kte,lnr) = nr3d(kts:kte)*RHO(kts:kte)
        an1d(kts:kte,lni) = ni3d(kts:kte)*RHO(kts:kte)
        an1d(kts:kte,lns) = ns3d(kts:kte)*RHO(kts:kte)
        an1d(kts:kte,lnh) = ng3d(kts:kte)*RHO(kts:kte)
    END IF
    ! ALL DONE !!!!!!!!!!!
    RETURN
    END SUBROUTINE MORR_TWO_MOMENT_MICRO

    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

    REAL FUNCTION POLYSVP (T,TYPE)

        !-------------------------------------------

        !  COMPUTE SATURATION VAPOR PRESSURE

        !  POLYSVP RETURNED IN UNITS OF PA.
        !  T IS INPUT IN UNITS OF K.
        !  TYPE REFERS TO SATURATION WITH RESPECT TO LIQUID (0) OR ICE (1)

        ! REPLACE GOFF-GRATCH WITH FASTER FORMULATION FROM FLATAU ET AL. 1992, TABLE 4 (RIGHT-HAND COLUMN)

        IMPLICIT NONE

        REAL DUM
        REAL T
        INTEGER TYPE
        ! ice
        real a0i,a1i,a2i,a3i,a4i,a5i,a6i,a7i,a8i
        data a0i,a1i,a2i,a3i,a4i,a5i,a6i,a7i,a8i /&
                6.11147274, 0.503160820, 0.188439774e-1, &
                0.420895665e-3, 0.615021634e-5,0.602588177e-7, &
                0.385852041e-9, 0.146898966e-11, 0.252751365e-14/

        ! liquid
        real a0,a1,a2,a3,a4,a5,a6,a7,a8

        ! V1.7
        data a0,a1,a2,a3,a4,a5,a6,a7,a8 /&
                6.11239921, 0.443987641, 0.142986287e-1, &
                0.264847430e-3, 0.302950461e-5, 0.206739458e-7, &
                0.640689451e-10,-0.952447341e-13,-0.976195544e-15/
        real dt

        ! ICE

        IF (TYPE.EQ.1) THEN

            !         POLYSVP = 10.**(-9.09718*(273.16/T-1.)-3.56654*                &
            !          LOG10(273.16/T)+0.876793*(1.-T/273.16)+						&
            !          LOG10(6.1071))*100.


            dt = max(-80.,t-273.16)
            polysvp = a0i + dt*(a1i+dt*(a2i+dt*(a3i+dt*(a4i+dt*(a5i+dt*(a6i+dt*(a7i+a8i*dt)))))))
            polysvp = polysvp*100.

        END IF

        ! LIQUID

        IF (TYPE.EQ.0) THEN

            dt = max(-80.,t-273.16)
            polysvp = a0 + dt*(a1+dt*(a2+dt*(a3+dt*(a4+dt*(a5+dt*(a6+dt*(a7+a8*dt)))))))
            polysvp = polysvp*100.

            !         POLYSVP = 10.**(-7.90298*(373.16/T-1.)+                        &
            !             5.02808*LOG10(373.16/T)-									&
            !             1.3816E-7*(10**(11.344*(1.-T/373.16))-1.)+				&
            !             8.1328E-3*(10**(-3.49149*(373.16/T-1.))-1.)+				&
            !             LOG10(1013.246))*100.

        END IF


    END FUNCTION POLYSVP

    !------------------------------------------------------------------------------

    REAL FUNCTION GAMMA(X)
        !----------------------------------------------------------------------
        !
        ! THIS ROUTINE CALCULATES THE GAMMA FUNCTION FOR A REAL ARGUMENT X.
        !   COMPUTATION IS BASED ON AN ALGORITHM OUTLINED IN REFERENCE 1.
        !   THE PROGRAM USES RATIONAL FUNCTIONS THAT APPROXIMATE THE GAMMA
        !   FUNCTION TO AT LEAST 20 SIGNIFICANT DECIMAL DIGITS.  COEFFICIENTS
        !   FOR THE APPROXIMATION OVER THE INTERVAL (1,2) ARE UNPUBLISHED.
        !   THOSE FOR THE APPROXIMATION FOR X .GE. 12 ARE FROM REFERENCE 2.
        !   THE ACCURACY ACHIEVED DEPENDS ON THE ARITHMETIC SYSTEM, THE
        !   COMPILER, THE INTRINSIC FUNCTIONS, AND PROPER SELECTION OF THE
        !   MACHINE-DEPENDENT CONSTANTS.
        !
        !
        !*******************************************************************
        !*******************************************************************
        !
        ! EXPLANATION OF MACHINE-DEPENDENT CONSTANTS
        !
        ! BETA   - RADIX FOR THE FLOATING-POINT REPRESENTATION
        ! MAXEXP - THE SMALLEST POSITIVE POWER OF BETA THAT OVERFLOWS
        ! XBIG   - THE LARGEST ARGUMENT FOR WHICH GAMMA(X) IS REPRESENTABLE
        !          IN THE MACHINE, I.E., THE SOLUTION TO THE EQUATION
        !                  GAMMA(XBIG) = BETA**MAXEXP
        ! XINF   - THE LARGEST MACHINE REPRESENTABLE FLOATING-POINT NUMBER;
        !          APPROXIMATELY BETA**MAXEXP
        ! EPS    - THE SMALLEST POSITIVE FLOATING-POINT NUMBER SUCH THAT
        !          1.0+EPS .GT. 1.0
        ! XMININ - THE SMALLEST POSITIVE FLOATING-POINT NUMBER SUCH THAT
        !          1/XMININ IS MACHINE REPRESENTABLE
        !
        !     APPROXIMATE VALUES FOR SOME IMPORTANT MACHINES ARE:
        !
        !                            BETA       MAXEXP        XBIG
        !
        ! CRAY-1         (S.P.)        2         8191        966.961
        ! CYBER 180/855
        !   UNDER NOS    (S.P.)        2         1070        177.803
        ! IEEE (IBM/XT,
        !   SUN, ETC.)   (S.P.)        2          128        35.040
        ! IEEE (IBM/XT,
        !   SUN, ETC.)   (D.P.)        2         1024        171.624
        ! IBM 3033       (D.P.)       16           63        57.574
        ! VAX D-FORMAT   (D.P.)        2          127        34.844
        ! VAX G-FORMAT   (D.P.)        2         1023        171.489
        !
        !                            XINF         EPS        XMININ
        !
        ! CRAY-1         (S.P.)   5.45E+2465   7.11E-15    1.84E-2466
        ! CYBER 180/855
        !   UNDER NOS    (S.P.)   1.26E+322    3.55E-15    3.14E-294
        ! IEEE (IBM/XT,
        !   SUN, ETC.)   (S.P.)   3.40E+38     1.19E-7     1.18E-38
        ! IEEE (IBM/XT,
        !   SUN, ETC.)   (D.P.)   1.79D+308    2.22D-16    2.23D-308
        ! IBM 3033       (D.P.)   7.23D+75     2.22D-16    1.39D-76
        ! VAX D-FORMAT   (D.P.)   1.70D+38     1.39D-17    5.88D-39
        ! VAX G-FORMAT   (D.P.)   8.98D+307    1.11D-16    1.12D-308
        !
        !*******************************************************************
        !*******************************************************************
        !
        ! ERROR RETURNS
        !
        !  THE PROGRAM RETURNS THE VALUE XINF FOR SINGULARITIES OR
        !     WHEN OVERFLOW WOULD OCCUR.  THE COMPUTATION IS BELIEVED
        !     TO BE FREE OF UNDERFLOW AND OVERFLOW.
        !
        !
        !  INTRINSIC FUNCTIONS REQUIRED ARE:
        !
        !     INT, DBLE, EXP, LOG, REAL, SIN
        !
        !
        ! REFERENCES:  AN OVERVIEW OF SOFTWARE DEVELOPMENT FOR SPECIAL
        !              FUNCTIONS   W. J. CODY, LECTURE NOTES IN MATHEMATICS,
        !              506, NUMERICAL ANALYSIS DUNDEE, 1975, G. A. WATSON
        !              (ED.), SPRINGER VERLAG, BERLIN, 1976.
        !
        !              COMPUTER APPROXIMATIONS, HART, ET. AL., WILEY AND
        !              SONS, NEW YORK, 1968.
        !
        !  LATEST MODIFICATION: OCTOBER 12, 1989
        !
        !  AUTHORS: W. J. CODY AND L. STOLTZ
        !           APPLIED MATHEMATICS DIVISION
        !           ARGONNE NATIONAL LABORATORY
        !           ARGONNE, IL 60439
        !
        !----------------------------------------------------------------------
        implicit none
        INTEGER I,N
        LOGICAL PARITY
        REAL                                                          &
                CONV,EPS,FACT,HALF,ONE,RES,SUM,TWELVE,                    &
                TWO,X,XBIG,XDEN,XINF,XMININ,XNUM,Y,Y1,YSQ,Z,ZERO
        REAL, DIMENSION(7) :: C
        REAL, DIMENSION(8) :: P
        REAL, DIMENSION(8) :: Q
        !----------------------------------------------------------------------
        !  MATHEMATICAL CONSTANTS
        !----------------------------------------------------------------------
        DATA ONE,HALF,TWELVE,TWO,ZERO/1.0E0,0.5E0,12.0E0,2.0E0,0.0E0/


        !----------------------------------------------------------------------
        !  MACHINE DEPENDENT PARAMETERS
        !----------------------------------------------------------------------
        DATA XBIG,XMININ,EPS/35.040E0,1.18E-38,1.19E-7/,XINF/3.4E38/
        !----------------------------------------------------------------------
        !  NUMERATOR AND DENOMINATOR COEFFICIENTS FOR RATIONAL MINIMAX
        !     APPROXIMATION OVER (1,2).
        !----------------------------------------------------------------------
        DATA P/-1.71618513886549492533811E+0,2.47656508055759199108314E+1,  &
                -3.79804256470945635097577E+2,6.29331155312818442661052E+2,  &
                8.66966202790413211295064E+2,-3.14512729688483675254357E+4,  &
                -3.61444134186911729807069E+4,6.64561438202405440627855E+4/
        DATA Q/-3.08402300119738975254353E+1,3.15350626979604161529144E+2,  &
                -1.01515636749021914166146E+3,-3.10777167157231109440444E+3, &
                2.25381184209801510330112E+4,4.75584627752788110767815E+3,  &
                -1.34659959864969306392456E+5,-1.15132259675553483497211E+5/
        !----------------------------------------------------------------------
        !  COEFFICIENTS FOR MINIMAX APPROXIMATION OVER (12, INF).
        !----------------------------------------------------------------------
        DATA C/-1.910444077728E-03,8.4171387781295E-04,                      &
                -5.952379913043012E-04,7.93650793500350248E-04,				   &
                -2.777777777777681622553E-03,8.333333333333333331554247E-02,	   &
                5.7083835261E-03/
        !----------------------------------------------------------------------
        !  STATEMENT FUNCTIONS FOR CONVERSION BETWEEN INTEGER AND FLOAT
        !----------------------------------------------------------------------
        CONV(I) = REAL(I)
        PARITY=.FALSE.
        FACT=ONE
        N=0
        Y=X
        IF(Y.LE.ZERO)THEN
            !----------------------------------------------------------------------
            !  ARGUMENT IS NEGATIVE
            !----------------------------------------------------------------------
            Y=-X
            Y1=AINT(Y)
            RES=Y-Y1
            IF(RES.NE.ZERO)THEN
                IF(Y1.NE.AINT(Y1*HALF)*TWO)PARITY=.TRUE.
                FACT=-PI/SIN(PI*RES)
                Y=Y+ONE
            ELSE
                RES=XINF
                GOTO 900
            ENDIF
        ENDIF
        !----------------------------------------------------------------------
        !  ARGUMENT IS POSITIVE
        !----------------------------------------------------------------------
        IF(Y.LT.EPS)THEN
            !----------------------------------------------------------------------
            !  ARGUMENT .LT. EPS
            !----------------------------------------------------------------------
            IF(Y.GE.XMININ)THEN
                RES=ONE/Y
            ELSE
                RES=XINF
                GOTO 900
            ENDIF
        ELSEIF(Y.LT.TWELVE)THEN
            Y1=Y
            IF(Y.LT.ONE)THEN
                !----------------------------------------------------------------------
                !  0.0 .LT. ARGUMENT .LT. 1.0
                !----------------------------------------------------------------------
                Z=Y
                Y=Y+ONE
            ELSE
                !----------------------------------------------------------------------
                !  1.0 .LT. ARGUMENT .LT. 12.0, REDUCE ARGUMENT IF NECESSARY
                !----------------------------------------------------------------------
                N=INT(Y)-1
                Y=Y-CONV(N)
                Z=Y-ONE
            ENDIF
            !----------------------------------------------------------------------
            !  EVALUATE APPROXIMATION FOR 1.0 .LT. ARGUMENT .LT. 2.0
            !----------------------------------------------------------------------
            XNUM=ZERO
            XDEN=ONE
            DO I=1,8
                XNUM=(XNUM+P(I))*Z
                XDEN=XDEN*Z+Q(I)
            END DO
            RES=XNUM/XDEN+ONE
            IF(Y1.LT.Y)THEN
                !----------------------------------------------------------------------
                !  ADJUST RESULT FOR CASE  0.0 .LT. ARGUMENT .LT. 1.0
                !----------------------------------------------------------------------
                RES=RES/Y1
            ELSEIF(Y1.GT.Y)THEN
                !----------------------------------------------------------------------
                !  ADJUST RESULT FOR CASE  2.0 .LT. ARGUMENT .LT. 12.0
                !----------------------------------------------------------------------
                DO I=1,N
                    RES=RES*Y
                    Y=Y+ONE
                END DO
            ENDIF
        ELSE
            !----------------------------------------------------------------------
            !  EVALUATE FOR ARGUMENT .GE. 12.0,
            !----------------------------------------------------------------------
            IF(Y.LE.XBIG)THEN
                YSQ=Y*Y
                SUM=C(7)
                DO I=1,6
                    SUM=SUM/YSQ+C(I)
                END DO
                SUM=SUM/Y-Y+xxx
                SUM=SUM+(Y-HALF)*LOG(Y)
                RES=EXP(SUM)
            ELSE
                RES=XINF
                GOTO 900
            ENDIF
        ENDIF
        !----------------------------------------------------------------------
        !  FINAL ADJUSTMENTS AND RETURN
        !----------------------------------------------------------------------
        IF(PARITY)RES=-RES
        IF(FACT.NE.ONE)RES=FACT/RES
        900 GAMMA=RES
        RETURN
        ! ---------- LAST LINE OF GAMMA ----------
    END FUNCTION GAMMA


    REAL FUNCTION DERF1(X)
        IMPLICIT NONE
        REAL X
        REAL, DIMENSION(0 : 64) :: A, B
        REAL W,T,Y
        INTEGER K,I
        DATA A/                                                 &
                0.00000000005958930743E0, -0.00000000113739022964E0, &
                0.00000001466005199839E0, -0.00000016350354461960E0, &
                0.00000164610044809620E0, -0.00001492559551950604E0, &
                0.00012055331122299265E0, -0.00085483269811296660E0, &
                0.00522397762482322257E0, -0.02686617064507733420E0, &
                0.11283791670954881569E0, -0.37612638903183748117E0, &
                1.12837916709551257377E0,	                          &
                0.00000000002372510631E0, -0.00000000045493253732E0, &
                0.00000000590362766598E0, -0.00000006642090827576E0, &
                0.00000067595634268133E0, -0.00000621188515924000E0, &
                0.00005103883009709690E0, -0.00037015410692956173E0, &
                0.00233307631218880978E0, -0.01254988477182192210E0, &
                0.05657061146827041994E0, -0.21379664776456006580E0, &
                0.84270079294971486929E0,							  &
                0.00000000000949905026E0, -0.00000000018310229805E0, &
                0.00000000239463074000E0, -0.00000002721444369609E0, &
                0.00000028045522331686E0, -0.00000261830022482897E0, &
                0.00002195455056768781E0, -0.00016358986921372656E0, &
                0.00107052153564110318E0, -0.00608284718113590151E0, &
                0.02986978465246258244E0, -0.13055593046562267625E0, &
                0.67493323603965504676E0, 							  &
                0.00000000000382722073E0, -0.00000000007421598602E0, &
                0.00000000097930574080E0, -0.00000001126008898854E0, &
                0.00000011775134830784E0, -0.00000111992758382650E0, &
                0.00000962023443095201E0, -0.00007404402135070773E0, &
                0.00050689993654144881E0, -0.00307553051439272889E0, &
                0.01668977892553165586E0, -0.08548534594781312114E0, &
                0.56909076642393639985E0,							  &
                0.00000000000155296588E0, -0.00000000003032205868E0, &
                0.00000000040424830707E0, -0.00000000471135111493E0, &
                0.00000005011915876293E0, -0.00000048722516178974E0, &
                0.00000430683284629395E0, -0.00003445026145385764E0, &
                0.00024879276133931664E0, -0.00162940941748079288E0, &
                0.00988786373932350462E0, -0.05962426839442303805E0, &
                0.49766113250947636708E0 /
        DATA (B(I), I = 0, 12) /                                  &
                -0.00000000029734388465E0,  0.00000000269776334046E0, 	&
                -0.00000000640788827665E0, -0.00000001667820132100E0,  &
                -0.00000021854388148686E0,  0.00000266246030457984E0, 	&
                0.00001612722157047886E0, -0.00025616361025506629E0, 	&
                0.00015380842432375365E0,  0.00815533022524927908E0, 	&
                -0.01402283663896319337E0, -0.19746892495383021487E0,  &
                0.71511720328842845913E0 /
        DATA (B(I), I = 13, 25) /                                 &
                -0.00000000001951073787E0, -0.00000000032302692214E0,  &
                0.00000000522461866919E0,  0.00000000342940918551E0, 	&
                -0.00000035772874310272E0,  0.00000019999935792654E0, 	&
                0.00002687044575042908E0, -0.00011843240273775776E0, 	&
                -0.00080991728956032271E0,  0.00661062970502241174E0, 	&
                0.00909530922354827295E0, -0.20160072778491013140E0, 	&
                0.51169696718727644908E0 /
        DATA (B(I), I = 26, 38) /                                 &
                0.00000000003147682272E0, -0.00000000048465972408E0,   &
                0.00000000063675740242E0,  0.00000003377623323271E0, 	&
                -0.00000015451139637086E0, -0.00000203340624738438E0, 	&
                0.00001947204525295057E0,  0.00002854147231653228E0, 	&
                -0.00101565063152200272E0,  0.00271187003520095655E0, 	&
                0.02328095035422810727E0, -0.16725021123116877197E0, 	&
                0.32490054966649436974E0 /
        DATA (B(I), I = 39, 51) /                                 &
                0.00000000002319363370E0, -0.00000000006303206648E0,   &
                -0.00000000264888267434E0,  0.00000002050708040581E0, 	&
                0.00000011371857327578E0, -0.00000211211337219663E0, 	&
                0.00000368797328322935E0,  0.00009823686253424796E0, 	&
                -0.00065860243990455368E0, -0.00075285814895230877E0, 	&
                0.02585434424202960464E0, -0.11637092784486193258E0, 	&
                0.18267336775296612024E0 /
        DATA (B(I), I = 52, 64) /                                 &
                -0.00000000000367789363E0,  0.00000000020876046746E0, 	&
                -0.00000000193319027226E0, -0.00000000435953392472E0, 	&
                0.00000018006992266137E0, -0.00000078441223763969E0, 	&
                -0.00000675407647949153E0,  0.00008428418334440096E0, 	&
                -0.00017604388937031815E0, -0.00239729611435071610E0, 	&
                0.02064129023876022970E0, -0.06905562880005864105E0,   &
                0.09084526782065478489E0 /
        W = ABS(X)
        IF (W .LT. 2.2D0) THEN
            T = W * W
            K = INT(T)
            T = T - K
            K = K * 13
            Y = ((((((((((((A(K) * T + A(K + 1)) * T +              &
                    A(K + 2)) * T + A(K + 3)) * T + A(K + 4)) * T +     &
                    A(K + 5)) * T + A(K + 6)) * T + A(K + 7)) * T +     &
                    A(K + 8)) * T + A(K + 9)) * T + A(K + 10)) * T + 	  &
                    A(K + 11)) * T + A(K + 12)) * W
        ELSE IF (W .LT. 6.9D0) THEN
            K = INT(W)
            T = W - K
            K = 13 * (K - 2)
            Y = (((((((((((B(K) * T + B(K + 1)) * T +               &
                    B(K + 2)) * T + B(K + 3)) * T + B(K + 4)) * T + 	  &
                    B(K + 5)) * T + B(K + 6)) * T + B(K + 7)) * T + 	  &
                    B(K + 8)) * T + B(K + 9)) * T + B(K + 10)) * T + 	  &
                    B(K + 11)) * T + B(K + 12)
            Y = Y * Y
            Y = Y * Y
            Y = Y * Y
            Y = 1 - Y * Y
        ELSE
            Y = 1
        END IF
        IF (X .LT. 0) Y = -Y
        DERF1 = Y
    END FUNCTION DERF1

    !+---+-----------------------------------------------------------------+

    subroutine refl10cm_hm (qv1d, qr1d, nr1d, qs1d, ns1d, qg1d, ng1d, &
            t1d, p1d, dBZ, kts, kte, ii, jj)

        IMPLICIT NONE

        !..Sub arguments
        INTEGER, INTENT(IN):: kts, kte, ii, jj
        REAL, DIMENSION(kts:kte), INTENT(IN)::                            &
                qv1d, qr1d, nr1d, qs1d, ns1d, qg1d, ng1d, t1d, p1d
        REAL, DIMENSION(kts:kte), INTENT(INOUT):: dBZ

        !..Local variables
        REAL, DIMENSION(kts:kte):: temp, pres, qv, rho
        REAL, DIMENSION(kts:kte):: rr, nr, rs, ns, rg, ng

        DOUBLE PRECISION, DIMENSION(kts:kte):: ilamr, ilamg, ilams
        DOUBLE PRECISION, DIMENSION(kts:kte):: N0_r, N0_g, N0_s
        DOUBLE PRECISION:: lamr, lamg, lams
        LOGICAL, DIMENSION(kts:kte):: L_qr, L_qs, L_qg

        REAL, DIMENSION(kts:kte):: ze_rain, ze_snow, ze_graupel
        DOUBLE PRECISION:: fmelt_s, fmelt_g
        DOUBLE PRECISION:: cback, x, eta, f_d

        INTEGER:: i, k, k_0, kbot, n
        LOGICAL:: melti

        !+---+

        do k = kts, kte
            dBZ(k) = -35.0
        enddo

        !+---+-----------------------------------------------------------------+
        !..Put column of data into local arrays.
        !+---+-----------------------------------------------------------------+
        do k = kts, kte
            temp(k) = t1d(k)
            qv(k) = MAX(1.E-10, qv1d(k))
            pres(k) = p1d(k)
            rho(k) = 0.622*pres(k)/(R*temp(k)*(qv(k)+0.622))

            if (qr1d(k) .gt. 1.E-9) then
                rr(k) = qr1d(k)*rho(k)
                nr(k) = nr1d(k)*rho(k)
                lamr = (xam_r*xcrg(3)*xorg2*nr(k)/rr(k))**xobmr
                ilamr(k) = 1./lamr
                N0_r(k) = nr(k)*xorg2*lamr**xcre(2)
                L_qr(k) = .true.
            else
                rr(k) = 1.E-12
                nr(k) = 1.E-12
                L_qr(k) = .false.
            endif

            if (qs1d(k) .gt. 1.E-9) then
                rs(k) = qs1d(k)*rho(k)
                ns(k) = ns1d(k)*rho(k)
                lams = (xam_s*xcsg(3)*xosg2*ns(k)/rs(k))**xobms
                ilams(k) = 1./lams
                N0_s(k) = ns(k)*xosg2*lams**xcse(2)
                L_qs(k) = .true.
            else
                rs(k) = 1.E-12
                ns(k) = 1.E-12
                L_qs(k) = .false.
            endif

            if (qg1d(k) .gt. 1.E-9) then
                rg(k) = qg1d(k)*rho(k)
                ng(k) = ng1d(k)*rho(k)
                lamg = (xam_g*xcgg(3)*xogg2*ng(k)/rg(k))**xobmg
                ilamg(k) = 1./lamg
                N0_g(k) = ng(k)*xogg2*lamg**xcge(2)
                L_qg(k) = .true.
            else
                rg(k) = 1.E-12
                ng(k) = 1.E-12
                L_qg(k) = .false.
            endif
        enddo

        !+---+-----------------------------------------------------------------+
        !..Locate K-level of start of melting (k_0 is level above).
        !+---+-----------------------------------------------------------------+
        melti = .false.
        k_0 = kts
        do k = kte-1, kts, -1
            if ( (temp(k).gt.273.15) .and. L_qr(k)                         &
                    .and. (L_qs(k+1).or.L_qg(k+1)) ) then
                k_0 = MAX(k+1, k_0)
                melti=.true.
                goto 195
            endif
        enddo
        195  continue

        !+---+-----------------------------------------------------------------+
        !..Assume Rayleigh approximation at 10 cm wavelength. Rain (all temps)
        !.. and non-water-coated snow and graupel when below freezing are
        !.. simple. Integrations of m(D)*m(D)*N(D)*dD.
        !+---+-----------------------------------------------------------------+

        do k = kts, kte
            ze_rain(k) = 1.e-22
            ze_snow(k) = 1.e-22
            ze_graupel(k) = 1.e-22
            if (L_qr(k)) ze_rain(k) = N0_r(k)*xcrg(4)*ilamr(k)**xcre(4)
            if (L_qs(k)) ze_snow(k) = (0.176/0.93) * (6.0/PI)*(6.0/PI)     &
                    * (xam_s/900.0)*(xam_s/900.0)          &
                    * N0_s(k)*xcsg(4)*ilams(k)**xcse(4)
            if (L_qg(k)) ze_graupel(k) = (0.176/0.93) * (6.0/PI)*(6.0/PI)  &
                    * (xam_g/900.0)*(xam_g/900.0)       &
                    * N0_g(k)*xcgg(4)*ilamg(k)**xcge(4)
        enddo

        !+---+-----------------------------------------------------------------+
        !..Special case of melting ice (snow/graupel) particles.  Assume the
        !.. ice is surrounded by the liquid water.  Fraction of meltwater is
        !.. extremely simple based on amount found above the melting level.
        !.. Uses code from Uli Blahak (rayleigh_soak_wetgraupel and supporting
        !.. routines).
        !+---+-----------------------------------------------------------------+

        if (melti .and. k_0.ge.kts+1) then
            do k = k_0-1, kts, -1

                !..Reflectivity contributed by melting snow
                if (L_qs(k) .and. L_qs(k_0) ) then
                    fmelt_s = MAX(0.005d0, MIN(1.0d0-rs(k)/rs(k_0), 0.99d0))
                    eta = 0.d0
                    lams = 1./ilams(k)
                    do n = 1, nrbins
                        x = xam_s * xxDs(n)**xbm_s
                        call rayleigh_soak_wetgraupel (x,DBLE(xocms),DBLE(xobms), &
                                fmelt_s, melt_outside_s, m_w_0, m_i_0, lamda_radar, &
                                CBACK, mixingrulestring_s, matrixstring_s,          &
                                inclusionstring_s, hoststring_s,                    &
                                hostmatrixstring_s, hostinclusionstring_s)
                        f_d = N0_s(k)*xxDs(n)**xmu_s * DEXP(-lams*xxDs(n))
                        eta = eta + f_d * CBACK * simpson(n) * xdts(n)
                    enddo
                    ze_snow(k) = SNGL(lamda4 / (pi5 * K_w) * eta)
                endif


                !..Reflectivity contributed by melting graupel

                if (L_qg(k) .and. L_qg(k_0) ) then
                    fmelt_g = MAX(0.005d0, MIN(1.0d0-rg(k)/rg(k_0), 0.99d0))
                    eta = 0.d0
                    lamg = 1./ilamg(k)
                    do n = 1, nrbins
                        x = xam_g * xxDg(n)**xbm_g
                        call rayleigh_soak_wetgraupel (x,DBLE(xocmg),DBLE(xobmg), &
                                fmelt_g, melt_outside_g, m_w_0, m_i_0, lamda_radar, &
                                CBACK, mixingrulestring_g, matrixstring_g,          &
                                inclusionstring_g, hoststring_g,                    &
                                hostmatrixstring_g, hostinclusionstring_g)
                        f_d = N0_g(k)*xxDg(n)**xmu_g * DEXP(-lamg*xxDg(n))
                        eta = eta + f_d * CBACK * simpson(n) * xdtg(n)
                    enddo
                    ze_graupel(k) = SNGL(lamda4 / (pi5 * K_w) * eta)
                endif

            enddo
        endif

        do k = kte, kts, -1
            dBZ(k) = 10.*log10((ze_rain(k)+ze_snow(k)+ze_graupel(k))*1.d18)
        enddo


    end subroutine refl10cm_hm

    !+---+-----------------------------------------------------------------+

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
        integer k

        real tmp, del
        real x1, x2, x3, x4
        integer i

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

    !+---+-----------------------------------------------------------------+

    Function delabk(ba,bb,nua,nub,mua,mub,k)

        implicit none
        real delabk
        real nua, mua, ba
        integer k
        real nub, mub, bb

        integer i
        real tmp,del

        real g1pnua, g2pnua, g1pbapnua, g1pbbpk, g1pnub, g2pnub

        tmp = (1. + nua)/mua
        i = Int(dgami*(tmp))
        del = tmp - dgam*i
        IF ( i+1 > ngm0 ) THEN
            write(0,*) 'delabk: i+1 > ngm0!!!!',i,ngm0,nua,mua,tmp
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

    !+---+-----------------------------------------------------------------+

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
        !      qlwc = qcw*(1.e3)*rho0

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

            !        alf = Min(5.0*(awdia/d0)**2*Abs(vt)/v0,10.0)
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
            !      ELSEIF ( temcg .gt. -40.0 ) THEN
            !        fac = 1.0 - ((temcg+30.0)/(40.0 - 30.0))**2
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
        !      parameter (tc = -23.7)
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

    ! #####################################################################
    ! #####################################################################
    !
    ! Subroutine for reattachment of floating ion charge
    !
    SUBROUTINE SCIONATT    &
            &  (nx,ny,nz,na,jyslab &
            &  ,nor,norz &
            &  ,dtp,dz3d &
            &  ,an &
            !            &  ,dn,p2 &
            !            &  ,pn,w &
            !            &  ,tmp3d &
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
        IF (IHAIL.EQ.1)  RHOG = 900.

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
                        if (an(ix,jy,kz,lc) > qxmin(lc)) fracw=((an(ix,jy,kz,lc)/rhow)**(2./3.) ) * Max(an(ix,jy,kz,lnc), 0.0)**(1./3.)
                        if (an(ix,jy,kz,lr) > qxmin(lr)) fracr=((an(ix,jy,kz,lr)/rhor)**(2./3.) ) * Max(an(ix,jy,kz,lnr), 0.0)**(1./3.)
                        if (an(ix,jy,kz,li) > qxmin(li)) fraci=((an(ix,jy,kz,li)/rhoi)**(2./3.) ) * Max(an(ix,jy,kz,lni), 0.0)**(1./3.)
                        if (an(ix,jy,kz,ls) > qxmin(ls)) fracs=((an(ix,jy,kz,ls)/rhos)**(2./3.) ) * Max(an(ix,jy,kz,lns), 0.0)**(1./3.)
                        if (an(ix,jy,kz,lh) > qxmin(lh)) frach=((an(ix,jy,kz,lh)/rhog)**(2./3.) ) * Max(an(ix,jy,kz,lnh), 0.0)**(1./3.)
                        !                        IF ( lhl > 1 ) THEN
                        !                            if (an(ix,jy,kz,lhl) > qxmin(lhl)) frachl=((an(ix,jy,kz,lhl)/rhoh)**(2./3.) ) * Max(an(ix,jy,kz,lnhl), 0.0)**(1./3.)
                        !                        ENDIF


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
END MODULE module_mp_morr_two_moment
!+---+-----------------------------------------------------------------+
!+---+-----------------------------------------------------------------+
