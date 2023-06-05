!*==IDENT.spg  processed by SPAG 6.72Dc at 22:10 on 27 May 2023
!
!-----------------------------------------------------------------------
!
! ****** Calculate the polarization brightness and brightness
! ****** from the electron density (produced by the MAS code).
!
!-----------------------------------------------------------------------
!
! ****** Updates and bug fixes:
!
!        05/19/94, ZM, Version 1.00:
!
!         - Original version of program.
!
!        01/24/95, ZM, Version 1.01:
!
!         - Generalized to handle 3D fields.
!         - Added the capability to enhance the density perturbation.
!         - Added the capability to apply a "vignetting function"
!           to pB.
!
!        02/06/95, ZM, Version 1.02:
!
!         - Added the ability to rotate the output pB image
!           in the theta direction.
!
!        04/28/95, ZM, Version 1.03:
!
!         - Added the ability to multiply pB by a radial vignetting
!           function.
!
!        03/03/97, ZM, Version 1.04:
!
!         - Added the ability to specify the solar B0 and P angles.
!
!        06/15/2004, ZM, Version 1.05:
!
!         - Updated the program to be able to deal with
!           nonuniform meshes in phi.
!
!        05/20/2009, ZM, Version 1.06:
!
!         - Allowed file names to be longer.
!
!        07/06/2010, ZM, Version 1.07:
!
!         - Corrected a bug in routine RDHDF that caused the program
!           to crash if the file being read from is not present.
!           This is a well-know bug in the HDF library that has been
!           corrected in our modern tools libraries.
!
!        08/30/2010, ZM, Version 1.08:
!
!         - Updated the program to use ZM's tools libraries,
!           and to use FORTRAN90.  This was a much-awaited
!           upgrade.  The new program uses a simpler line-of-sight
!           integration scheme, modeled after the one used in
!           GETEIT, with an adaptive integration step size
!           based on the local density mesh.
!         - The computation of the brightness B was added as
!           an option.
!         - Note that the usage line for this version is not
!           compatible with previous versions.
!         - Clarified the normalization of the output pB and B.
!           The output pB and B are in units of I0, the central
!           disk brightness, 2.49e10 [erg/cm^2/s/sr].
!
!        10/25/2018, ZM, Version 1.09:
!
!         - Generalized the program to allow the use of non-parallel
!           lines of sight to compute pB and B.  This allows
!           computations in regions that are large compared to the
!           distance from the Sun to the observer (such as for the
!           STEREO HI imagers and for WISPR on Parker Solar Probe).
!           The previous version assumed that the lines of sight
!           were parallel (i.e., corresponding to an observer at
!           infinity).  This new version has the capability to work
!           as before (i.e., an infinitely distant observer), which
!           is the default behavior.  To use the new feature (i.e.,
!           to calculate along non-parallel lines of sight), specify
!           the distance of the observer from the Sun in AU
!           [astronomical units] using the option -r <dist>.  The
!           previous results can be obtained by omitting this option.
!           When using this new option, for domains that are small
!           compared to the Sun-observer distance, the results
!           ought to be close to those calculated previously.  When
!           using the new option, the image domain [X0,X1] x [Y0,Y1]
!           is specified using elongation and altitude angles
!           [in degrees].
!         - Added shortcuts for the image view for the WISPR cameras
!           on PSP.  Use the flags -wispr1 and -wispr2 to specify
!           the image coordinates for these cameras.
!         - Removed the ability to specify the limits of the
!           integration along the line of sight (-z0 and -z1).
!           These had dubious merit.  The integration is now
!           always performed throughout the whole simulation
!           domain.
!         - This version has OpenMP commands to parallelize the
!           main loop of the computation.  These were adapted from
!           RL's version 1.08b from 04/12/2015.
!
!        08/20/2019, CD, Version 1.10:
!
!         - Updated the integrators to start and end only within the
!           bounds of the density mesh. This improves accuracy.
!         - Added the option for cubic spline interpolation, activated
!           with the -cubic flag.
!         - Cubic interpolation is used for rho and vignette function.
!         - These updates help if one applies special processing
!           to the images, which might otherwise bring out grid and/or
!           interpolation artifacts buried in the image.
!         - getpb now requires the zm_spline library as a dependency.
!
!        02/27/2020, CD, Version 1.11:
!
!         - Bugfix: Integrator update in v1.10 didn't account for LOSs
!           that pass outside the density mesh. Should be fine now.
!
!        11/16/2020, CD, Version 1.12:
!
!         - Added -he_frac flag to get ne from MAS rho if he_frac != 0.
!
!        04/01/2021, CD, Version 1.13:
!
!         - MAJOR UPDATE: Add option to compute emissivity weighted
!           integrals of scalar and vector fields.
!           Set the -help flag for a description of the new options.
!         - Add -mu flag to specify a custom limb-darkening constant.
!         - Rework how files are read and written to be more generic.
!         - SIDE EFFECT: If -scalar or -vr, vt, vp files are read in,
!           then the integration limits are adjusted to be the minimum
!           radius of every file read in (including rho). This can
!           change the integration limits from the half to main mesh,
!           leading to small floating point differences in pb & b
!           compared to previous versions (but that is OK!).
!
!-----------------------------------------------------------------------
!
!#######################################################################
!
! ****** Note that this tool uses modules from ZM's tools library.
!
!#######################################################################
      MODULE IDENT
      IMPLICIT NONE
!*--IDENT147
!
      CHARACTER(*) , PARAMETER :: CNAME = 'GETPB'
      CHARACTER(*) , PARAMETER :: CVERS = '1.13'
      CHARACTER(*) , PARAMETER :: CDATE = '04/22/2021'
!
      END
!*==TYPES.spg  processed by SPAG 6.72Dc at 22:10 on 27 May 2023
!#######################################################################
      MODULE TYPES
!
!-----------------------------------------------------------------------
! ****** Definition of data structures.
!-----------------------------------------------------------------------
!
!
      USE NUMBER_TYPES
      USE SDS_DEF
      USE INVINT_DEF
      USE SPLINE_DEF
      IMPLICIT NONE
!*--TYPES168
!
! ****** Inverse interpolation table structure definitions.
!
      STRUCTURE ::VTAB
         RECORD /ITAB  /  , DIMENSION(3)::C
      ENDSTRUCTURE
!
! ****** Container for 2D/3D MAS style fields and interpolation tables.
!
      STRUCTURE ::MAS_FIELD
         RECORD /SDS   / ::SDS
         RECORD /VTAB  / ::INVTAB
         RECORD /SPL2D / ::SPL2
         RECORD /SPL3D / ::SPL3
      ENDSTRUCTURE
!
      END
!*==CONSTANTS.spg  processed by SPAG 6.72Dc at 22:10 on 27 May 2023
!#######################################################################
      MODULE CONSTANTS
!
!
      USE NUMBER_TYPES
      IMPLICIT NONE
!*--CONSTANTS193
!
      REAL(r_typ) , PARAMETER :: PI = 3.14159265358979323846_R_TYP
      REAL(r_typ) , PARAMETER :: TWOPI = PI*2._R_TYP
      REAL(r_typ) , PARAMETER :: HALFPI = PI*.5_R_TYP
      REAL(r_typ) , PARAMETER :: DEGTORAD = PI/180._R_TYP
!
! ****** One AU in solar radii.
!
      REAL(r_typ) , PARAMETER :: AU_IN_RS = 214.939469_R_TYP
!
! ****** Offset to ensure integration is inside domain.
!
      REAL(r_typ) , PARAMETER :: INTEGRATION_OFFSET = 1.E-10_R_TYP
!
      END
!*==PARAMS.spg  processed by SPAG 6.72Dc at 22:10 on 27 May 2023
!#######################################################################
      MODULE PARAMS
!
!-----------------------------------------------------------------------
! ****** Parameters.
!-----------------------------------------------------------------------
!
!
      USE NUMBER_TYPES
      IMPLICIT NONE
!*--PARAMS220
!
      CHARACTER(512) :: rho_file
!
      CHARACTER(512) :: vr_file , vt_file , vp_file
      CHARACTER(512) :: scalar_file
!
      CHARACTER(512) :: vlos_file_out , vx_file_out , vy_file_out
      CHARACTER(512) :: scalar_file_out
      CHARACTER(512) :: angle_file_out
!
      CHARACTER(512) :: pb_file
      CHARACTER(512) :: b_file
!
      CHARACTER(512) :: vf_file
!
      LOGICAL :: verbose
!
      LOGICAL :: compute_pb
      LOGICAL :: compute_b
!
! ****** Options for emissivity weighted integrals.
!
      LOGICAL :: do_integral_weighting = .FALSE.
      LOGICAL :: compute_vector_integration = .FALSE.
      LOGICAL :: compute_scalar_integration = .FALSE.
      LOGICAL :: compute_angle_integration = .FALSE.
!
! ****** Weight the scalar/vector integrals by B instead (default pB).
!
      LOGICAL :: weight_integral_by_b = .FALSE.
!
! ****** Flag indicating that the density file is from the
! ****** old MAS code.
!
      LOGICAL :: oldmas
!
      REAL(r_typ) :: crlong
      REAL(r_typ) :: b0angle
      REAL(r_typ) :: pangle
      REAL(r_typ) :: rocc
      REAL(r_typ) :: power
      REAL(r_typ) :: disk_value
!
      LOGICAL :: power_set
!
! ****** Distance of the observer from the Sun.
!
      LOGICAL :: r_obs_set
      REAL(r_typ) :: r_obs
      REAL(r_typ) :: r_obs_rs
!
! ****** Step size multiplier.
!
      REAL(r_typ) :: dsmult
!
! ****** Shorthand flags to set the WISPR camera domains.
!
      LOGICAL :: wispr1
      LOGICAL :: wispr2
!
! ****** Helium Fraction Parameters (like in MAS).
! ****** HE_RHO is defined by rho(norm)=HE_RHO*n_e(norm),
!
      REAL(r_typ) :: he_frac = 0.0
      REAL(r_typ) :: he_rho
!
! ****** Option to use cubic interpolation
!
      LOGICAL :: cubic
!
! ****** Limb darkening coefficient at 520 nm. Default is from
! ****** Altschuler & Perry, Solar Physics, 23, 410 (1972).
!
      REAL(r_typ) :: mu_ld = 0.63_R_TYP
!
! ****** Flag for the precision of the output files.
!
      LOGICAL :: hdf32
!
      END
!*==IMAGE_REGION.spg  processed by SPAG 6.72Dc at 22:10 on 27 May 2023
!#######################################################################
      MODULE IMAGE_REGION
!
!
      USE NUMBER_TYPES
      IMPLICIT NONE
!*--IMAGE_REGION308
!
! ****** Image dimensions.
!
      INTEGER :: nx
      INTEGER :: ny
!
! ****** Image region limits.
!
      REAL(r_typ) :: x0
      REAL(r_typ) :: x1
      REAL(r_typ) :: y0
      REAL(r_typ) :: y1
!
      LOGICAL :: image_limits_set
!
! ****** Number of OpenMP iterations to do in each thread.
!
      INTEGER , PARAMETER :: ITERATIONS_PER_THREAD = 500
!
      END
!*==GEOMETRY.spg  processed by SPAG 6.72Dc at 22:10 on 27 May 2023
!#######################################################################
      MODULE GEOMETRY
!
!
      USE NUMBER_TYPES
      IMPLICIT NONE
!*--GEOMETRY336
!
! ****** Sin and cos of the P and B0 angles.
!
      REAL(r_typ) :: sin_p , cos_p
      REAL(r_typ) :: sin_b0 , cos_b0
!
! ****** Central meridian longitude.
!
      REAL(r_typ) :: cml
!
! ****** Flag to indicate an axisymmetric density (i.e., 2D).
!
      LOGICAL :: twodee
!
! ****** Maximum radius for all of the input meshes.
!
      REAL(r_typ) :: rmax = -1.0
!
      END
!*==MAS_FIELDS.spg  processed by SPAG 6.72Dc at 22:10 on 27 May 2023
!#######################################################################
      MODULE MAS_FIELDS
!
!-----------------------------------------------------------------------
! ****** Storage for the density.
!-----------------------------------------------------------------------
!
!
      USE NUMBER_TYPES
      USE TYPES
      USE SPLINE_DEF
      IMPLICIT NONE
!*--MAS_FIELDS369
!
! ****** Density field container.
!
      RECORD /MAS_FIELD/ ::rho
!
! ****** Containers for arbitrary vector field, v, in spherical coords.
!
      RECORD /MAS_FIELD/ ::vr
      RECORD /MAS_FIELD/ ::vt
      RECORD /MAS_FIELD/ ::vp
!
! ****** Containers for arbitrary scalar field.
!
      RECORD /MAS_FIELD/ ::sfield
!
! ****** Normalization factor for electron density.
!
! ****** Multiply the normalized density read in [in MAS code units]
! ****** by FN_N to get electron density in [/cm^3].
!
      REAL(r_typ) , PARAMETER :: FN_N = 1.E8_R_TYP
!
      END
!*==IMAGE_FIELDS.spg  processed by SPAG 6.72Dc at 22:10 on 27 May 2023
!#######################################################################
      MODULE IMAGE_FIELDS
!
!-----------------------------------------------------------------------
! ****** Storage for the computed pB and B, and their associated
! ****** arrays.
!-----------------------------------------------------------------------
!
!
      USE NUMBER_TYPES
      USE TYPES
      IMPLICIT NONE
!*--IMAGE_FIELDS406
!
! ****** Image axes for the case when the observer distance
! ****** has been specified (elongation, altitude).
!
      REAL(r_typ) , DIMENSION(:) , POINTER :: elong , alt
!
! ****** Image axes for the case when the observer distance
! ****** is infinite (POS x and y).
!
      REAL(r_typ) , DIMENSION(:) , POINTER :: x_pos , y_pos
!
! ****** Polarization brightness and brightness.
!
      REAL(r_typ) , DIMENSION(:,:) , ALLOCATABLE :: pb
      REAL(r_typ) , DIMENSION(:,:) , ALLOCATABLE :: b
!
! ****** Emissivity weighted averages of the vector field components
! ****** in the frame of the image plane (los is TOWARDS observer).
!
      REAL(r_typ) , DIMENSION(:,:) , ALLOCATABLE :: vlos_avg
      REAL(r_typ) , DIMENSION(:,:) , ALLOCATABLE :: vx_avg
      REAL(r_typ) , DIMENSION(:,:) , ALLOCATABLE :: vy_avg
!
! ****** Emissivity weighted averages of the scalar field.
!
      REAL(r_typ) , DIMENSION(:,:) , ALLOCATABLE :: sf_avg
!
! ****** Emissivity weighted average elongation angle (0 when r=rmin).
!
      REAL(r_typ) , DIMENSION(:,:) , ALLOCATABLE :: angle_avg
!
! ****** Closest distance to the Sun.  This is used for vignetting
! ****** and multiplication by a radial power.
!
      REAL(r_typ) , DIMENSION(:,:) , ALLOCATABLE :: r_min
!
      END
!*==VIGNETTING_FUNCTION.spg  processed by SPAG 6.72Dc at 22:10 on 27 May 2023
!#######################################################################
      MODULE VIGNETTING_FUNCTION
!
!-----------------------------------------------------------------------
! ****** Vignetting function definition.
!-----------------------------------------------------------------------
!
!
      USE NUMBER_TYPES
      USE SPLINE_DEF
      IMPLICIT NONE
!*--VIGNETTING_FUNCTION456
!
! ****** Vignetting function.
!
      INTEGER :: nvf
      REAL(r_typ) , DIMENSION(:) , ALLOCATABLE :: r_vf
      REAL(r_typ) , DIMENSION(:) , ALLOCATABLE :: vf
!
! ****** Spline structure for cubic interpolation.
!
      RECORD /SPL1D / ::vf_spl1
!
      END
!*==GET_IMAGE_OBS_DISTANCE_SPECIFIED.spg  processed by SPAG 6.72Dc at 22:10 on 27 May 2023
!#######################################################################
      SUBROUTINE GET_IMAGE_OBS_DISTANCE_SPECIFIED
!
!-----------------------------------------------------------------------
!
! ****** Compute pB and/or B for the case when the observer distance
! ****** has been specified.
!
! ****** This uses non-parallel lines of sight for the computation,
! ****** and should be used when the image domain is of significant
! ****** size compared to the Sun-observer distance.
!
! ****** In this case, the image coordinates [X0,X1] x [Y0,Y1] are
! ****** specified in terms of elongation and altitude angles with
! ****** respect to the ecliptic plane.
!
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!
      USE IDENT
      USE NUMBER_TYPES
      USE CONSTANTS
      USE PARAMS
      USE IMAGE_REGION
      USE GEOMETRY
      USE IMAGE_FIELDS
      IMPLICIT NONE
!*--GET_IMAGE_OBS_DISTANCE_SPECIFIED499
!
!-----------------------------------------------------------------------
!
      REAL(R_Typ) , PARAMETER :: HALF = .5_R_TYP
      REAL(R_Typ) , PARAMETER :: ONE = 1._R_TYP
!
!-----------------------------------------------------------------------
!
      INTEGER :: i , j , nstop
      REAL(R_Typ) :: dx , dy
      REAL(R_Typ) :: s0 , s1
      REAL(R_Typ) :: dmax
      REAL(R_Typ) :: e_rad , a_rad , cos_a
      REAL(R_Typ) :: s , t , ds , dsm
      REAL(R_Typ) :: t_r_min
      REAL(R_Typ) :: ds_local
      REAL(R_Typ) :: rv , tv , pv
      REAL(R_Typ) :: pb_local , b_local
      REAL(R_Typ) , DIMENSION(3) :: v_obs , v_ref , v_los , v_r_min
!
!-----------------------------------------------------------------------
!
      REAL(R_Typ) :: xo , yo , zo , xr , yr , zr , xv , yv , zv
      REAL(R_Typ) :: ds_fac , x_los
!
! ****** Unit vectors of observer frame in heliographic coordinates.
!
      REAL(R_Typ) , DIMENSION(3) :: ex , ey , ez
!
!-----------------------------------------------------------------------
!
! ****** Note that the (x,y,z) coordinates are defined such
! ****** that the x-axis points from the Sun to the observer,
! ****** and the z-axis points vertically upward.
!
! ****** Construct the coordinates for the axes of the image.
! ****** These are specified in terms of [elongation,altitude]
! ****** angles, in degrees.
!
! ****** The elongation is the angle between the Sun-observer line
! ****** and the LOS, measured positive towards the West limb
! ****** (with the Sun being at elongation=0).
!
! ****** The altitude is the angle expressing the vertical
! ****** inclination of the LOS with respect to the ecliptic
! ****** plane (with positive angles specifying a LOS above the
! ****** ecliptic plane, and negative angles specifying a LOS
! ****** below the ecliptic plane.
!
      ALLOCATE (ELOng(NX))
      ALLOCATE (ALT(NY))
!
      IF ( NX.EQ.1 ) THEN
         ELOng(1) = X0
      ELSE
         dx = (X1-X0)/(NX-1)
         DO i = 1 , NX
            ELOng(i) = X0 + (i-1)*dx
         ENDDO
      ENDIF
!
      IF ( NY.EQ.1 ) THEN
         ALT(1) = Y0
      ELSE
         dy = (Y1-Y0)/(NY-1)
         DO j = 1 , NY
            ALT(j) = Y0 + (j-1)*dy
         ENDDO
      ENDIF
!
! ****** Cartesian vector from the Sun to the observer.
!
      v_obs(1) = R_Obs_rs
      v_obs(2) = 0.
      v_obs(3) = 0.
!
! ****** Calculate pB and/or B.
!
!$omp parallel do shared(pb,b,r_min)
!$omp& shared(sf_avg,angle_avg,vlos_avg,vx_avg,vy_avg)
!$omp& shared(v_obs)
!$omp& private(ex,ey,ez)
!$omp& private(i,j,ds_fac,x_los,xo,yo,zo,xr,yr,zr,xv,yv,zv)
!$omp& private(e_rad,a_rad,cos_a,v_ref,t_r_min,v_r_min)
!$omp& private(t,s,dsm,ds,ds_local,nstop,v_los,rv,tv,pv,s0,s1,dmax)
!$omp& private(pb_local,b_local)
!$omp& collapse(2)
!$omp& schedule(dynamic,iterations_per_thread)
      DO j = 1 , NY
         DO i = 1 , NX
!
! ****** Cartesian vector from the Sun to the reference point.
!
! ****** This is a point on a sphere of radius R_OBS centered
! ****** at the observer, where R_OBS is the distance of the
! ****** observer from the Sun.
!
            e_rad = DEGTORAD*ELOng(i)
            a_rad = DEGTORAD*ALT(j)
            cos_a = COS(a_rad)
            v_ref(1) = R_Obs_rs*(ONE+COS(PI-e_rad)*cos_a)
            v_ref(2) = R_Obs_rs*SIN(PI-e_rad)*cos_a
            v_ref(3) = R_Obs_rs*SIN(a_rad)
!
! ****** Find the point along this LOS that is closest to the
! ****** Sun.  This is the point on the LOS at which it
! ****** intersects the Thomson sphere.  The distance from this
! ****** point to the Sun is needed in the calculation of pB
! ****** and B.
!
            t_r_min = -(v_obs(1)*(v_ref(1)-v_obs(1))+v_obs(2)           &
                    & *(v_ref(2)-v_obs(2))+v_obs(3)*(v_ref(3)-v_obs(3)))&
                    & /R_Obs_rs**2
!
            v_r_min(1) = t_r_min*v_ref(1) + (ONE-t_r_min)*v_obs(1)
            v_r_min(2) = t_r_min*v_ref(2) + (ONE-t_r_min)*v_obs(2)
            v_r_min(3) = t_r_min*v_ref(3) + (ONE-t_r_min)*v_obs(3)
!
            R_Min(i,j) = SQRT(v_r_min(1)**2+v_r_min(2)**2+v_r_min(3)**2)
!
! ****** Check if this LOS intersects the occulting disk.
! ****** If so, skip the calculation along this LOS.
!
            IF ( R_Min(i,j).LE.MAX(ROCc,ONE) ) THEN
               IF ( COMpute_pb ) PB(i,j) = DISk_value
               IF ( COMpute_b ) B(i,j) = DISk_value
               GOTO 50
            ENDIF
!
! ****** The distance s is along the LOS, starting at the observer
! ****** (s = 0), and increasing toward the reference point
! ****** (s = R_OBS).  The parameter t is proportional to s,
! ****** with t = 0 specifying the observer, and t = 1 specifying
! ****** the reference point.
!
            IF ( COMpute_pb ) PB(i,j) = 0.
            IF ( COMpute_b ) B(i,j) = 0.
!
! ****** Check if the LOS doesn't intersect the grid, leave it zeroed.
!
            IF ( R_Min(i,j).GE.RMAx ) GOTO 50
!
! ****** Compute unit vectors for LOS if doing vector/angle projection.
! ****** Here the projection changes for each pixel (non-parallel).
!
            IF ( COMpute_vector_integration .OR.                        &
               & COMpute_angle_integration ) THEN
               CALL TRANSFORM_POSITION(v_obs(1),v_obs(2),v_obs(3),rv,tv,&
                                     & pv)
               CALL RTP2XYZ(rv,tv,pv,xo,yo,zo)
               CALL TRANSFORM_POSITION(v_r_min(1),v_r_min(2),v_r_min(3),&
                                     & rv,tv,pv)
               CALL RTP2XYZ(rv,tv,pv,xr,yr,zr)
               xv = xo - xr
               yv = yo - yr
               zv = zo - zr
               CALL GET_UNIT_VECTORS(xv,yv,zv,ex,ey,ez)
            ENDIF
!
! ****** Integrate along the LOS, starting at first intersection
! ****** of the LOS with the density mesh.
!
            dmax = SQRT(RMAx**2-R_Min(i,j)**2)
            s0 = t_r_min*R_Obs_rs - dmax + INTEGRATION_OFFSET
            s1 = t_r_min*R_Obs_rs + dmax - INTEGRATION_OFFSET
!
            s = s0
            dsm = 0.
            nstop = 0
            DO
!
! ****** Set the position along the LOS.
!
               t = s/R_Obs_rs
!
               v_los(1) = t*v_ref(1) + (ONE-t)*v_obs(1)
               v_los(2) = t*v_ref(2) + (ONE-t)*v_obs(2)
               v_los(3) = t*v_ref(3) + (ONE-t)*v_obs(3)
!
               CALL TRANSFORM_POSITION(v_los(1),v_los(2),v_los(3),rv,tv,&
                                     & pv)
!
               CALL GET_KERNELS(R_Min(i,j),rv,tv,pv,COMpute_pb,         &
                              & COMpute_b,pb_local,b_local,ds_local)
!
               ds = ds_local*DSMult
               s = s + ds
               IF ( s.GE.s1 ) THEN
                  ds = ds - (s-s1)
                  s = s1
                  nstop = nstop + 1
               ENDIF
               ds_fac = HALF*(dsm+ds)
!
               IF ( COMpute_pb ) PB(i,j) = PB(i,j) + ds_fac*pb_local
!
               IF ( COMpute_b ) B(i,j) = B(i,j) + ds_fac*b_local
!
! ****** Compute emissivity weighted integrals if needed.
! ****** For non-plane-parallel angle need to dot position with x-hat.
!
               IF ( DO_integral_weighting ) THEN
                  IF ( COMpute_angle_integration ) THEN
                     CALL RTP2XYZ(rv,tv,pv,xv,yv,zv)
                     x_los = ex(1)*xv + ex(2)*yv + ex(3)*zv
                  ELSE
                     x_los = 0.
                  ENDIF
!
                  CALL ADD_TO_INTEGRALS(i,j,pb_local,b_local,ds_fac,ex, &
                                      & ey,ez,rv,tv,pv,x_los)
               ENDIF
!
               dsm = ds
               IF ( nstop.GE.2 ) GOTO 50
            ENDDO
 50      ENDDO
      ENDDO
!$omp end parallel do
!
      END
!*==GET_IMAGE_OBS_AT_INFINITY.spg  processed by SPAG 6.72Dc at 22:10 on 27 May 2023
!#######################################################################
      SUBROUTINE GET_IMAGE_OBS_AT_INFINITY
!
!-----------------------------------------------------------------------
!
! ****** Compute pB and/or B for the case when the observer is at
! ****** infinity.
!
! ****** This uses parallel lines of sight for the computation.
!
! ****** In this case, the image coordinates [X0,X1] x [Y0,Y1] are
! ****** specified in terms of horizontal and vertical coordinates
! ****** in the plane of the sky.
!
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!
      USE IDENT
      USE NUMBER_TYPES
      USE CONSTANTS
      USE PARAMS
      USE IMAGE_REGION
      USE GEOMETRY
      USE IMAGE_FIELDS
      USE MAS_FIELDS
      IMPLICIT NONE
!*--GET_IMAGE_OBS_AT_INFINITY750
!
!-----------------------------------------------------------------------
!
      REAL(R_Typ) , PARAMETER :: ZERO = 0._R_TYP
      REAL(R_Typ) , PARAMETER :: HALF = .5_R_TYP
      REAL(R_Typ) , PARAMETER :: ONE = 1._R_TYP
!
!-----------------------------------------------------------------------
!
      INTEGER :: i , j , nstop
      REAL(R_Typ) :: dx , dy
      REAL(R_Typ) :: s0 , s1
      REAL(R_Typ) :: dmax
      REAL(R_Typ) :: s , ds , dsm
      REAL(R_Typ) :: ds_local
      REAL(R_Typ) :: rv , tv , pv
      REAL(R_Typ) :: pb_local , b_local
      REAL(R_Typ) , DIMENSION(3) :: v_los
!
!-----------------------------------------------------------------------
!
      REAL(R_Typ) :: xv , yv , zv
      REAL(R_Typ) :: ds_fac
!
! ****** Unit vectors of observer frame in heliographic coordinates.
!
      REAL(R_Typ) , DIMENSION(3) :: ex , ey , ez
!
!-----------------------------------------------------------------------
!
! ****** Note that the (x,y,z) coordinates are defined such
! ****** that the x-axis points from the Sun to the observer,
! ****** and the z-axis points vertically upward.
!
! ****** Construct the coordinates for the axes of the image.
! ****** These are specified in terms of horizontal and vertical
! ****** plane-of-sky positions, in solar radii.
!
      ALLOCATE (X_Pos(NX))
      ALLOCATE (Y_Pos(NY))
!
      IF ( NX.EQ.1 ) THEN
         X_Pos(1) = X0
      ELSE
         dx = (X1-X0)/(NX-1)
         DO i = 1 , NX
            X_Pos(i) = X0 + (i-1)*dx
         ENDDO
      ENDIF
!
      IF ( NY.EQ.1 ) THEN
         Y_Pos(1) = Y0
      ELSE
         dy = (Y1-Y0)/(NY-1)
         DO j = 1 , NY
            Y_Pos(j) = Y0 + (j-1)*dy
         ENDDO
      ENDIF
!
! ****** Compute unit vectors for LOS if doing vector projection.
! ****** For the plane-parallel assumption, they never change.
!
      IF ( COMpute_vector_integration ) THEN
         CALL TRANSFORM_POSITION(ONE,ZERO,ZERO,rv,tv,pv)
         CALL RTP2XYZ(rv,tv,pv,xv,yv,zv)
         CALL GET_UNIT_VECTORS(xv,yv,zv,ex,ey,ez)
      ENDIF
!
! ****** Calculate pB and/or B.
!
!$omp parallel do shared(pb,b,r_min)
!$omp& shared(sf_avg,angle_avg,vlos_avg,vx_avg,vy_avg)
!$omp& shared(ex,ey,ez)
!$omp& private(i,j,ds_fac)
!$omp& private(s,dsm,ds,ds_local,nstop,v_los,rv,tv,pv,s0,s1,dmax)
!$omp& private(pb_local,b_local)
!$omp& collapse(2)
!$omp& schedule(dynamic,iterations_per_thread)
      DO j = 1 , NY
         DO i = 1 , NX
!
! ****** Calculate the plane-of-sky radius.
!
            R_Min(i,j) = SQRT(X_Pos(i)**2+Y_Pos(j)**2)
!
! ****** Check if this LOS intersects the occulting disk.
! ****** If so, skip the calculation along this LOS.
!
            IF ( R_Min(i,j).LE.MAX(ROCc,ONE) ) THEN
               IF ( COMpute_pb ) PB(i,j) = DISk_value
               IF ( COMpute_b ) B(i,j) = DISk_value
               GOTO 50
            ENDIF
!
! ****** The distance s is along the LOS, referenced with respect
! ****** to the solar limb (s = 0), and increasing away from the
! ****** observer (so s = -x == v_los(1) in the observer frame).
!
            IF ( COMpute_pb ) PB(i,j) = 0.
            IF ( COMpute_b ) B(i,j) = 0.
!
! ****** Check if the LOS doesn't intersect the grid, leave it zeroed.
!
            IF ( R_Min(i,j).GE.RMAx ) GOTO 50
!
! ****** Integrate along the LOS, starting at first intersection
! ****** of the LOS with the density mesh.
!
            dmax = SQRT(RMAx**2-R_Min(i,j)**2)
            s0 = -dmax + INTEGRATION_OFFSET
            s1 = dmax - INTEGRATION_OFFSET
!
            s = s0
            dsm = 0.
            nstop = 0
            DO
!
! ****** Set the position along the LOS.
!
               v_los(1) = -s
               v_los(2) = X_Pos(i)
               v_los(3) = Y_Pos(j)
!
               CALL TRANSFORM_POSITION(v_los(1),v_los(2),v_los(3),rv,tv,&
                                     & pv)
!
               CALL GET_KERNELS(R_Min(i,j),rv,tv,pv,COMpute_pb,         &
                              & COMpute_b,pb_local,b_local,ds_local)
!
               ds = ds_local*DSMult
               s = s + ds
               IF ( s.GE.s1 ) THEN
                  ds = ds - (s-s1)
                  s = s1
                  nstop = nstop + 1
               ENDIF
               ds_fac = HALF*(dsm+ds)
!
! ****** Compute pb and/or b contributions.
!
               IF ( COMpute_pb ) PB(i,j) = PB(i,j) + ds_fac*pb_local
               IF ( COMpute_b ) B(i,j) = B(i,j) + ds_fac*b_local
!
! ****** Compute emissivity weighted integrals if needed.
!
!
               IF ( DO_integral_weighting )                             &
                  & CALL ADD_TO_INTEGRALS(i,j,pb_local,b_local,ds_fac,  &
                  & ex,ey,ez,rv,tv,pv,v_los(1))
!
               dsm = ds
               IF ( nstop.GE.2 ) GOTO 50
            ENDDO
!
 50      ENDDO
      ENDDO
!$omp end parallel do
!
      END
!*==READ_FIELD.spg  processed by SPAG 6.72Dc at 22:10 on 27 May 2023
!#######################################################################
      SUBROUTINE READ_FIELD(Filename,Fieldname,Field)
!
!-----------------------------------------------------------------------
!
! ****** Read a MAS style field.
!
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!
      USE IDENT
      USE PARAMS
      USE GEOMETRY
      USE TYPES
      IMPLICIT NONE
!*--READ_FIELD928
!
!-----------------------------------------------------------------------
!
      CHARACTER(512) :: Filename
      CHARACTER(*) :: Fieldname
      RECORD /MAS_FIELD/ ::Field
!
!-----------------------------------------------------------------------
!
      INTEGER :: ierr
!
!-----------------------------------------------------------------------
!
      IF ( VERbose ) THEN
         WRITE (*,*)
         WRITE (*,*) '### Reading ' , TRIM(Fieldname) , ' from file: ' ,&
                   & TRIM(Filename)
      ENDIF
!
! ****** Read the field from the specified file.
!
      CALL RDHDF(Filename,Field%SDS,ierr)
!
      IF ( ierr.NE.0 ) THEN
         WRITE (*,*)
         WRITE (*,*) '### ERROR in ' , CNAME , ':'
         WRITE (*,*) '### Could not read ' , TRIM(Fieldname) , '.'
         WRITE (*,*) 'IERR (from RDHDF) = ' , ierr
         WRITE (*,*) 'File name: ' , TRIM(Filename)
         CALL EXIT(1)
      ENDIF
!
! ****** Check that it contains a 2D or 3D field.
!
      IF ( Field%SDS%NDIM.NE.2 .AND. Field%SDS%NDIM.NE.3 ) THEN
         WRITE (*,*)
         WRITE (*,*) '### ERROR in ' , CNAME , ':'
         WRITE (*,*) '### ' , TRIM(Fieldname) ,                         &
                    &' file does'//' not contain a 2D or 3D field.'
         WRITE (*,*) 'Number of dimensions = ' , Field%SDS%NDIM
         WRITE (*,*) 'File name: ' , TRIM(Filename)
         CALL EXIT(1)
      ENDIF
!
! ****** Check that it contains scales.
!
      IF ( .NOT.Field%SDS%SCALE ) THEN
         WRITE (*,*)
         WRITE (*,*) '### ERROR in ' , CNAME , ':'
         WRITE (*,*) '### The density file does'//' not contain scales.'
         WRITE (*,*) 'File name: ' , TRIM(Filename)
         CALL EXIT(1)
      ENDIF
!
! ****** Set the flag for an axisymmetric case.
!
      IF ( Field%SDS%NDIM.EQ.2 ) THEN
         TWOdee = .TRUE.
      ELSE
         TWOdee = .FALSE.
      ENDIF
!
! ****** Check that the field array has at least two points in
! ****** each (non-negligible) dimension.
!
      IF ( Field%SDS%DIMS(1).LT.2 ) THEN
         WRITE (*,*)
         WRITE (*,*) '### ERROR in ' , CNAME , ':'
         WRITE (*,*) '### ' , TRIM(Fieldname) ,                         &
                    &' array must have at least'//                      &
                    &' 2 points in the r dimension.'
         WRITE (*,*) 'Number of points in r = ' , Field%SDS%DIMS(1)
         WRITE (*,*) 'File name: ' , TRIM(Filename)
         CALL EXIT(1)
      ENDIF
!
      IF ( Field%SDS%DIMS(2).LT.2 ) THEN
         WRITE (*,*)
         WRITE (*,*) '### ERROR in ' , CNAME , ':'
         WRITE (*,*) '### ' , TRIM(Fieldname) ,                         &
                    &' array must have at least'//                      &
                    &' 2 points in the theta dimension.'
         WRITE (*,*) 'Number of points in theta = ' , Field%SDS%DIMS(2)
         WRITE (*,*) 'File name: ' , TRIM(Filename)
         CALL EXIT(1)
      ENDIF
!
      IF ( .NOT.TWOdee ) THEN
         IF ( Field%SDS%DIMS(3).LT.2 ) THEN
            WRITE (*,*)
            WRITE (*,*) '### ERROR in ' , CNAME , ':'
            WRITE (*,*) '### ' , TRIM(Fieldname) ,                      &
                      & ' array must have at'//                         &
                       &' least 2 points in the phi dimension.'
            WRITE (*,*) 'Number of points in phi = ' , Field%SDS%DIMS(3)
            WRITE (*,*) 'File name: ' , TRIM(Filename)
            CALL EXIT(1)
         ENDIF
      ENDIF
!
! ****** If we are using old MAS files, add a phi point.
!
      IF ( OLDmas ) CALL ADD_PHI_POINT(Field%SDS)
!
! ****** Construct the inverse interpolation table for
! ****** faster interpolation.
!
      CALL BUILD_INVERSE_TABLES(Field%SDS,Field%INVTAB)
!
! ****** Set the maximum radius available in the file but check
! ****** that a smaller value hasn't been set yet (rmax is initialized
! ****** as a negative number).
!
      IF ( RMAx.LT.0. ) THEN
         RMAx = Field%SDS%SCALES(1)%F(Field%SDS%DIMS(1))
      ELSE
         RMAx = MIN(RMAx,Field%SDS%SCALES(1)%F(Field%SDS%DIMS(1)))
      ENDIF
!
! ****** Compute the spline coefficients if needed.
!
      IF ( CUBic ) THEN
         IF ( VERbose ) THEN
            WRITE (*,*)
            WRITE (*,*) '### Computing cubic spline coefficients for ' ,&
                      & TRIM(Fieldname)
         ENDIF
         IF ( Field%SDS%NDIM.EQ.2 ) THEN
            CALL COMPUTE_SPLINE_2D(Field%SDS%DIMS(1),Field%SDS%DIMS(2), &
                                 & Field%SDS%SCALES(1)%F,               &
                                 & Field%SDS%SCALES(2)%F,Field%SDS%F,   &
                                 & Field%SPL2)
         ELSEIF ( Field%SDS%NDIM.EQ.3 ) THEN
            CALL COMPUTE_SPLINE_3D(Field%SDS%DIMS(1),Field%SDS%DIMS(2), &
                                 & Field%SDS%DIMS(3),Field%SDS%SCALES(1)&
                                 & %F,Field%SDS%SCALES(2)%F,            &
                                 & Field%SDS%SCALES(3)%F,Field%SDS%F,   &
                                 & Field%SPL3)
         ENDIF
      ENDIF
!
      END
!*==ADD_PHI_POINT.spg  processed by SPAG 6.72Dc at 22:10 on 27 May 2023
!#######################################################################
      SUBROUTINE ADD_PHI_POINT(S)
!
!-----------------------------------------------------------------------
!
! ****** Add a point in the phi dimension in a periodic manner
! ****** to the SDS in structure S.
!
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!
      USE NUMBER_TYPES
      USE SDS_DEF
      USE CONSTANTS
      IMPLICIT NONE
!*--ADD_PHI_POINT1089
!
!-----------------------------------------------------------------------
!
      RECORD /SDS   / ::S
!
!-----------------------------------------------------------------------
!
      REAL(R_Typ) , DIMENSION(:,:,:) , POINTER :: f
      REAL(R_Typ) , DIMENSION(:) , POINTER :: p
!
!-----------------------------------------------------------------------
!
! ****** For a 2D (axisymmetric) case, this is not needed.
!
      IF ( S%NDIM.NE.3 ) RETURN
!
! ****** Add a point to the phi dimension.
!
      ALLOCATE (f(S%DIMS(1),S%DIMS(2),S%DIMS(3)+1))
      ALLOCATE (p(S%DIMS(3)+1))
!
      f(:,:,1:S%DIMS(3)) = S%f(:,:,:)
      f(:,:,S%DIMS(3)+1) = S%f(:,:,1)
!
      p(1:S%DIMS(3)) = S%SCALES(3)%f(:)
      p(S%DIMS(3)+1) = S%SCALES(3)%f(1) + TWOPI
      S%DIMS(3) = S%DIMS(3) + 1
!
      DEALLOCATE (S%f)
      DEALLOCATE (S%SCALES(3)%f)
!
      S%f => f
      S%SCALES(3)%f => p
!
      END
!*==READ_VIGNETTING_FUNCTION.spg  processed by SPAG 6.72Dc at 22:10 on 27 May 2023
!#######################################################################
      SUBROUTINE READ_VIGNETTING_FUNCTION(Fname)
!
!-----------------------------------------------------------------------
!
! ****** Read the vignetting function from file FNAME.
!
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!
      USE IDENT
      USE SDS_DEF
      USE PARAMS
      USE VIGNETTING_FUNCTION
      IMPLICIT NONE
!*--READ_VIGNETTING_FUNCTION1143
!
!-----------------------------------------------------------------------
!
      CHARACTER(*) :: Fname
!
!-----------------------------------------------------------------------
!
      RECORD /SDS   / ::s
      INTEGER :: ierr , i , j
!
!-----------------------------------------------------------------------
!
      IF ( VERbose ) THEN
         WRITE (*,*)
         WRITE (*,*) '### Reading the vignetting function'//            &
                    &' from file: ' , TRIM(Fname)
      ENDIF
!
! ****** Read the vignetting function from the specified file.
!
      CALL RDHDF(Fname,s,ierr)
!
      IF ( ierr.NE.0 ) THEN
         WRITE (*,*)
         WRITE (*,*) '### ERROR in ' , CNAME , ':'
         WRITE (*,*) '### Could not read the vignetting function.'
         WRITE (*,*) 'IERR (from RDHDF) = ' , ierr
         WRITE (*,*) 'File name: ' , TRIM(Fname)
         CALL EXIT(1)
      ENDIF
!
! ****** Check that it contains a 1D field.
!
      IF ( s%NDIM.NE.1 ) THEN
         WRITE (*,*)
         WRITE (*,*) '### ERROR in ' , CNAME , ':'
         WRITE (*,*) '### The vignetting function file does'//          &
                    &' not contain a 1D field.'
         WRITE (*,*) 'Number of dimensions = ' , s%NDIM
         WRITE (*,*) 'File name: ' , TRIM(Fname)
         CALL EXIT(1)
      ENDIF
!
! ****** Check that it contains scales.
!
      IF ( .NOT.s%SCALE ) THEN
         WRITE (*,*)
         WRITE (*,*) '### ERROR in ' , CNAME , ':'
         WRITE (*,*) '### The vignetting function file does'//          &
                    &' not contain scales.'
         WRITE (*,*) 'File name: ' , TRIM(Fname)
         CALL EXIT(1)
      ENDIF
!
      NVF = s%DIMS(1)
!
! ****** Check that it contains two or more points.
!
      IF ( NVF.LT.2 ) THEN
         WRITE (*,*)
         WRITE (*,*) '### ERROR in ' , CNAME , ':'
         WRITE (*,*) '### The vignetting function file has'//           &
                    &' less than two points.'
         WRITE (*,*) 'Number of points = ' , NVF
         WRITE (*,*) 'File name: ' , TRIM(Fname)
         CALL EXIT(1)
      ENDIF
!
! ****** Transfer the vignetting function to 1D arrays for
! ****** convenience.
!
      ALLOCATE (R_Vf(NVF))
      ALLOCATE (VF(NVF))
!
      R_Vf = s%SCALES(1)%F
      VF = s%F(:,1,1)
!
      CALL DEALLOCATE_SDS(s)
!
! ****** Check that the vignetting function scale is monotonically
! ****** increasing.
!
      DO i = 1 , NVF - 1
         IF ( R_Vf(i+1).LE.R_Vf(i) ) THEN
            WRITE (*,*)
            WRITE (*,*) '### ERROR in ' , CNAME , ':'
            WRITE (*,*) '### The vignetting function scale is not'//    &
                       &' monotonically increasing.'
            WRITE (*,*) 'File name: ' , TRIM(Fname)
            WRITE (*,*) 'Number of points = ' , NVF
            WRITE (*,*) 'Radial locations of points:'
            WRITE (*,*) 'Index' , CHAR(9) , 'Radius'
            DO j = 1 , NVF
               WRITE (*,*) j , CHAR(9) , R_Vf(j)
            ENDDO
            WRITE (*,*) 'Not monotonic at index = ' , i + 1
            CALL EXIT(1)
         ENDIF
      ENDDO
!
! ****** Compute the spline coefficients if needed.
!
      IF ( CUBic ) THEN
         IF ( VERbose ) THEN
            WRITE (*,*)
            WRITE (*,*) '### Computing cubic spline coefficients for vf'
         ENDIF
         CALL COMPUTE_SPLINE_1D(NVF,R_Vf,VF,VF_spl1)
      ENDIF
!
      END
!*==WRITE_IMAGE.spg  processed by SPAG 6.72Dc at 22:10 on 27 May 2023
!#######################################################################
      SUBROUTINE WRITE_IMAGE(Filename,Fieldname,X,Y,F)
!
!-----------------------------------------------------------------------
!
! ****** Write a 2D image (wrapper for wrdhdf_2d and error checking).
!
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!
      USE IDENT
      USE NUMBER_TYPES
      USE PARAMS
      USE IMAGE_REGION
      USE IMAGE_FIELDS
      IMPLICIT NONE
!*--WRITE_IMAGE1274
!
!-----------------------------------------------------------------------
!
      CHARACTER(512) :: Filename
      CHARACTER(*) :: Fieldname
      REAL(R_Typ) , DIMENSION(NX) :: X
      REAL(R_Typ) , DIMENSION(NY) :: Y
      REAL(R_Typ) , DIMENSION(NX,NY) :: F
!
!-----------------------------------------------------------------------
!
      INTEGER :: ierr
!
!-----------------------------------------------------------------------
!
      IF ( VERbose ) THEN
         WRITE (*,*)
         WRITE (*,*) '### Writing ' , TRIM(Fieldname) , ' to file: ' ,  &
                   & TRIM(Filename)
      ENDIF
!
      CALL WRHDF_2D(Filename,.TRUE.,NX,NY,F,X,Y,HDF32,ierr)
!
      IF ( ierr.NE.0 ) THEN
         WRITE (*,*)
         WRITE (*,*) '### ERROR in ' , CNAME , ':'
         WRITE (*,*) '### Could not write the ' , TRIM(Fieldname) ,     &
                    &' output file:'
         WRITE (*,*) 'IERR (from WRHDF_2D) = ' , ierr
         WRITE (*,*) 'File name: ' , TRIM(Fieldname)
         CALL EXIT(1)
      ENDIF
!
      END
!*==VIGNETTE.spg  processed by SPAG 6.72Dc at 22:10 on 27 May 2023
!#######################################################################
      SUBROUTINE VIGNETTE(F)
!
!-----------------------------------------------------------------------
!
! ****** Multiply the field F by the vignetting function, using
! ****** the radius R_MIN.
!
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!
      USE IDENT
      USE NUMBER_TYPES
      USE PARAMS
      USE IMAGE_REGION
      USE IMAGE_FIELDS
      USE VIGNETTING_FUNCTION
      USE EVALUATE_SPLINE_1D_INTERFACE
      IMPLICIT NONE
!*--VIGNETTE1331
!
!-----------------------------------------------------------------------
!
      REAL(R_Typ) , DIMENSION(NX,NY) :: F
!
!-----------------------------------------------------------------------
!
      REAL(R_Typ) , PARAMETER :: ONE = 1._R_TYP
!
!-----------------------------------------------------------------------
!
      INTEGER :: ierr , i , j , ivf , ivfp1
      REAL(R_Typ) :: rv , alpha , vfv
!
!-----------------------------------------------------------------------
!
      DO j = 1 , NY
         DO i = 1 , NX
            rv = R_Min(i,j)
            IF ( rv.LT.ROCc ) GOTO 50
            IF ( rv.LE.R_Vf(1) ) THEN
               vfv = VF(1)
            ELSEIF ( rv.GE.R_Vf(NVF) ) THEN
               vfv = VF(NVF)
            ELSEIF ( CUBic ) THEN
               vfv = EVALUATE_SPLINE_1D(VF_spl1,rv)
            ELSE
               CALL INTERP(NVF,R_Vf,rv,ivf,ivfp1,alpha,ierr)
               IF ( ierr.NE.0 ) THEN
                  WRITE (*,*)
                  WRITE (*,*) '### ERROR in ' , CNAME , ':'
                  WRITE (*,*) '### An error occurred while applying'//  &
                             &' the vignetting function.'
                  WRITE (*,*) 'Please check the vignetting'//           &
                             &' function definition.'
                  CALL EXIT(1)
               ENDIF
               vfv = (ONE-alpha)*VF(ivf) + alpha*VF(ivfp1)
            ENDIF
            F(i,j) = vfv*F(i,j)
 50      ENDDO
      ENDDO
!
      END
!*==MULTIPLY_BY_RADIAL_POWER.spg  processed by SPAG 6.72Dc at 22:10 on 27 May 2023
!#######################################################################
      SUBROUTINE MULTIPLY_BY_RADIAL_POWER(F)
!
!-----------------------------------------------------------------------
!
! ****** Multiply the field F by R_MIN**POWER.
!
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!
      USE IDENT
      USE NUMBER_TYPES
      USE PARAMS
      USE IMAGE_REGION
      USE IMAGE_FIELDS
      IMPLICIT NONE
!*--MULTIPLY_BY_RADIAL_POWER1395
!
!-----------------------------------------------------------------------
!
      REAL(R_Typ) , DIMENSION(NX,NY) :: F
!
!-----------------------------------------------------------------------
!
      INTEGER :: i , j
      REAL(R_Typ) :: rv
!
!-----------------------------------------------------------------------
!
      DO j = 1 , NY
         DO i = 1 , NX
            rv = R_Min(i,j)
            IF ( rv.LT.ROCc ) GOTO 50
            IF ( POWer.LT.0. ) THEN
               F(i,j) = F(i,j)/rv**ABS(POWer)
            ELSE
               F(i,j) = F(i,j)*rv**POWer
            ENDIF
 50      ENDDO
      ENDDO
!
      END
!*==TRANSFORM_POSITION.spg  processed by SPAG 6.72Dc at 22:10 on 27 May 2023
!#######################################################################
      SUBROUTINE TRANSFORM_POSITION(X,Y,Z,R,T,P)
!
!-----------------------------------------------------------------------
!
! ****** Apply the transformations due to the solar P angle and the
! ****** solar B0 angle, and rotate to the specified central meridian.
! ****** Return the transformed position in spherical coordinates
! ****** (R,T,P).
!
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!
      USE NUMBER_TYPES
      USE CONSTANTS
      USE GEOMETRY
      IMPLICIT NONE
!*--TRANSFORM_POSITION1441
!
!-----------------------------------------------------------------------
!
      REAL(R_Typ) :: X , Y , Z
      REAL(R_Typ) :: R , T , P
!
!-----------------------------------------------------------------------
!
      REAL(R_Typ) , PARAMETER :: ZERO = 0.
!
!-----------------------------------------------------------------------
!
      REAL(R_Typ) :: x1 , y1 , z1 , x2 , y2 , z2
!
!-----------------------------------------------------------------------
!
      REAL(R_Typ) , EXTERNAL :: FOLD
!
!-----------------------------------------------------------------------
!
! ****** Apply the solar P angle transformation.
!
      x1 = X
      y1 = COS_p*Y + SIN_p*Z
      z1 = -SIN_p*Y + COS_p*Z
!
! ****** Apply the solar B0 angle transformation.
!
      x2 = COS_b0*x1 - SIN_b0*z1
      y2 = y1
      z2 = SIN_b0*x1 + COS_b0*z1
!
! ****** Convert to spherical coordinates.
!
      CALL C2S(x2,y2,z2,R,T,P)
!
! ****** Transform to the specified central meridian longitude.
!
      P = FOLD(ZERO,TWOPI,P+CML)

      END
!*==FOLD.spg  processed by SPAG 6.72Dc at 22:10 on 27 May 2023
!#######################################################################
      FUNCTION FOLD(X0,X1,X)
!
!-----------------------------------------------------------------------
!
! ****** "Fold" X into the periodic interval [X0,X1].
!
! ****** On return, X is such that X0.le.X.lt.X1.
!
!-----------------------------------------------------------------------
!
! ****** It is assumed that X0 does not equal X1, as is physically
! ****** necessary.  If X0 and X1 are equal, the routine just
! ****** returns with FOLD=X.
!
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!
      USE NUMBER_TYPES
      IMPLICIT NONE
!*--FOLD1506
!
!-----------------------------------------------------------------------
!
      REAL(r_typ) :: FOLD
      REAL(r_typ) :: X0 , X1 , X
!
!-----------------------------------------------------------------------
!
      REAL(r_typ) :: xl
!
!-----------------------------------------------------------------------
!
      FOLD = X
!
      IF ( X0.EQ.X1 ) RETURN
!
      xl = X1 - X0
!
      FOLD = MOD(X-X0,xl) + X0
!
      IF ( FOLD.LT.X0 ) FOLD = FOLD + xl
      IF ( FOLD.GE.X1 ) FOLD = FOLD - xl
!
      END
!*==GET_KERNELS.spg  processed by SPAG 6.72Dc at 22:10 on 27 May 2023
!#######################################################################
      SUBROUTINE GET_KERNELS(R_min,Rv,Tv,Pv,Compute_pb,Compute_b,       &
                           & Pb_local,B_local,Ds)
!
!-----------------------------------------------------------------------
!
! ****** Calculate the contribution to the polarized brightness pB
! ****** and the brightness B at a point (RV,TV,PV) in spherical
! ****** coordinates.
!
!-----------------------------------------------------------------------
!
! ****** This formulation is adapted from Fran Bagenal, after the
! ****** treatment in Billings (1966) and Altschuler and Perry (1972).
!
!-----------------------------------------------------------------------
!
! ****** When COMPUTE_PB=.true., the polarization brightness pB is
! ****** returned in PB_LOCAL.
!
! ****** When COMPUTE_B=.true., the brightness B is returned
! ****** in B_LOCAL.
!
! ****** The smallest mesh spacing at the interpolation point,
! ****** calculated in Cartesian coordinates, is returned in DS.
!
! ****** R_MIN is the closest distance to the Sun along the line
! ****** of sight.  This point defines the Thomson sphere.
! ****** In the case of an infinitely-distant observer, R_MIN
! ****** becomes the plane-of-sky-radius (called rho by Altschuler
! ****** and Perry).
!
! ****** RV and R_MIN are in units of solar radii.
!
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!
      USE NUMBER_TYPES
      USE GEOMETRY
      USE MAS_FIELDS
      USE PARAMS , ONLY:HE_rho
      IMPLICIT NONE
!*--GET_KERNELS1576
!
!-----------------------------------------------------------------------
!
      REAL(R_Typ) :: R_min
      REAL(R_Typ) :: Rv , Tv , Pv
      LOGICAL :: Compute_pb , Compute_b
      REAL(R_Typ) :: Pb_local
      REAL(R_Typ) :: B_local
      REAL(R_Typ) :: Ds
!
!-----------------------------------------------------------------------
!
      REAL(R_Typ) , PARAMETER :: TWO = 2._R_TYP
!
!-----------------------------------------------------------------------
!
      REAL(R_Typ) :: n_e , cf , ef
      REAL(R_Typ) , DIMENSION(3) :: drtp
!
!-----------------------------------------------------------------------
!
      REAL(R_Typ) , EXTERNAL :: DS_S2C
      REAL(R_Typ) , EXTERNAL :: CFUNC_BILLINGS
      REAL(R_Typ) , EXTERNAL :: EFUNC_BILLINGS
!
!-----------------------------------------------------------------------
!
! ****** Interpolate the density.
!
      CALL INTERP_FIELD_DS(RHO%SDS,RHO%INVTAB,Rv,Tv,Pv,n_e,drtp,        &
                         & RHO%SPL2,RHO%SPL3)
!
! ****** Convert N_E to physical units.
!
      n_e = n_e*FN_N/HE_rho
!
! ****** Convert the mesh spacing from spherical coordinates
! ****** to a Cartesian spacing.
!
      Ds = DS_S2C(Rv,Tv,drtp,TWOdee)
!
! ****** Evaluate the "Billings" radial functions.
!
      IF ( Compute_pb .OR. Compute_b ) cf = CFUNC_BILLINGS(Rv)
!
      IF ( Compute_b ) ef = EFUNC_BILLINGS(Rv)
!
! ****** Obtain the local contribution to pB.
!
      IF ( Compute_pb ) THEN
         Pb_local = cf*(R_min/Rv)**2*n_e
      ELSE
         Pb_local = 0.
      ENDIF
!
! ****** Obtain the local contribution to B.
!
      IF ( Compute_b ) THEN
         B_local = (TWO*ef-cf*(R_min/Rv)**2)*n_e
      ELSE
         B_local = 0.
      ENDIF
!
      END
!*==CFUNC_BILLINGS.spg  processed by SPAG 6.72Dc at 22:10 on 27 May 2023
!#######################################################################
      FUNCTION CFUNC_BILLINGS(R)
!
!-----------------------------------------------------------------------
!
! ****** Coronal light scattering function for polarized brightness
! ****** [Billings 1966].
!
!-----------------------------------------------------------------------
!
! ****** R is the radius in units of R_sun.
!
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!
      USE NUMBER_TYPES
      USE PARAMS , ONLY:MU_ld
      IMPLICIT NONE
!*--CFUNC_BILLINGS1662
!
!-----------------------------------------------------------------------
!
      REAL(r_typ) :: R
      REAL(r_typ) :: CFUNC_BILLINGS
!
!-----------------------------------------------------------------------
!
      REAL(r_typ) , PARAMETER :: ONE = 1._R_TYP
      REAL(r_typ) , PARAMETER :: THREE = 3._R_TYP
      REAL(r_typ) , PARAMETER :: EIGHTH = .125_R_TYP
!
!-----------------------------------------------------------------------
!
! ****** Intensity coefficient [cm^3], normalized to I0,
! ****** the central disk brightness, 2.49e10 [erg/cm^2/s/sr].
!
      REAL(r_typ) , PARAMETER :: A = 8.69E-15_R_TYP
!
!-----------------------------------------------------------------------
!
      REAL(r_typ) :: so , co , so_sq , co_sq
      REAL(r_typ) :: log_term
      REAL(r_typ) :: afunc , bfunc
!
!-----------------------------------------------------------------------
!
! ****** Protect against invalid R values.  The physical requirement
! ****** is that R.ge.1.
!
      IF ( R.LE.ONE ) THEN
         so = ONE
         so_sq = ONE
         co = 0.
         co_sq = 0.
         log_term = 0.
      ELSE
         so = ONE/R
         so_sq = so**2
         co_sq = ONE - so_sq
         co = SQRT(ABS(co_sq))
         log_term = LOG((ONE+so)/co)
      ENDIF
!
      afunc = co*so_sq
!
      bfunc = -EIGHTH*(ONE-THREE*so_sq-(co_sq/so)*(ONE+THREE*so_sq)     &
            & *log_term)
!
      CFUNC_BILLINGS = A*((ONE-MU_ld)*afunc+MU_ld*bfunc)
!
      END
!*==EFUNC_BILLINGS.spg  processed by SPAG 6.72Dc at 22:10 on 27 May 2023
!#######################################################################
      FUNCTION EFUNC_BILLINGS(R)
!
!-----------------------------------------------------------------------
!
! ****** Coronal light scattering function for tangential
! ****** polarized brightness [Billings 1966].
!
!-----------------------------------------------------------------------
!
! ****** R is the radius in units of R_sun.
!
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!
      USE NUMBER_TYPES
      USE PARAMS , ONLY:MU_ld
      IMPLICIT NONE
!*--EFUNC_BILLINGS1736
!
!-----------------------------------------------------------------------
!
      REAL(r_typ) :: R
      REAL(r_typ) :: EFUNC_BILLINGS
!
!-----------------------------------------------------------------------
!
      REAL(r_typ) , PARAMETER :: ONE = 1._R_TYP
      REAL(r_typ) , PARAMETER :: THREE = 3._R_TYP
      REAL(r_typ) , PARAMETER :: FOUR = 4._R_TYP
      REAL(r_typ) , PARAMETER :: FIVE = 5._R_TYP
      REAL(r_typ) , PARAMETER :: EIGHTH = .125_R_TYP
      REAL(r_typ) , PARAMETER :: ONE_THIRD = ONE/THREE
      REAL(r_typ) , PARAMETER :: FOUR_THIRDS = FOUR/THREE
!
!-----------------------------------------------------------------------
!
! ****** Intensity coefficient [cm^3], normalized to I0,
! ****** the central disk brightness, 2.49e10 [erg/cm^2/s/sr].
!
      REAL(r_typ) , PARAMETER :: A = 8.69E-15_R_TYP
!
!-----------------------------------------------------------------------
!
      REAL(r_typ) :: so , co , so_sq , co_sq
      REAL(r_typ) :: log_term
      REAL(r_typ) :: cfunc , dfunc
!
!-----------------------------------------------------------------------
!
! ****** Protect against invalid R values.  The physical requirement
! ****** is that R.ge.1.
!
      IF ( R.LE.ONE ) THEN
         so = ONE
         so_sq = ONE
         co = 0.
         co_sq = 0.
         log_term = 0.
      ELSE
         so = ONE/R
         so_sq = so**2
         co_sq = ONE - so_sq
         co = SQRT(ABS(co_sq))
         log_term = LOG((ONE+so)/co)
      ENDIF
!
      cfunc = FOUR_THIRDS - co - ONE_THIRD*co*co_sq
!
      dfunc = EIGHTH*(FIVE+so_sq-(co_sq/so)*(FIVE-so_sq)*log_term)
!
      EFUNC_BILLINGS = A*((ONE-MU_ld)*cfunc+MU_ld*dfunc)
!
      END
!*==C2S.spg  processed by SPAG 6.72Dc at 22:10 on 27 May 2023
!#######################################################################
      SUBROUTINE C2S(X,Y,Z,R,T,P)
!
!-----------------------------------------------------------------------
!
! ****** Convert from Cartesian coordinates (X,Y,Z)
! ****** to spherical coordinates (R,T,P).
!
! ****** This routine returns T and P in radians, in the
! ****** following range:
!
!          0. .le. t .le. pi
!          0. .le. p .lt. 2.*pi
!
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!
      USE NUMBER_TYPES
      USE CONSTANTS
      IMPLICIT NONE
!*--C2S1815
!
!-----------------------------------------------------------------------
!
      REAL(R_Typ) :: X , Y , Z
      REAL(R_Typ) :: R , T , P
!
!-----------------------------------------------------------------------
!
      R = SQRT(X**2+Y**2+Z**2)
!
      IF ( R.EQ.0. ) THEN
         T = 0.
      ELSE
         T = ACOS(Z/R)
      ENDIF
!
      IF ( X.NE.0. ) THEN
         P = ATAN2(Y,X)
      ELSEIF ( Y.GE.0. ) THEN
         P = HALFPI
      ELSE
         P = -HALFPI
      ENDIF
      IF ( P.LT.0. ) P = P + TWOPI
!
      END
!*==S2C.spg  processed by SPAG 6.72Dc at 22:10 on 27 May 2023
!#######################################################################
      SUBROUTINE S2C(R,T,P,X,Y,Z)
!
!-----------------------------------------------------------------------
!
! ****** Convert from spherical coordinates (R,T,P)
! ****** to Cartesian coordinates (X,Y,Z).
!
! ****** This routine assumes that T and P are in radians.
!
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!
      USE NUMBER_TYPES
      IMPLICIT NONE
!*--S2C1860
!
!-----------------------------------------------------------------------
!
      REAL(r_typ) :: R , T , P
      REAL(r_typ) :: X , Y , Z
!
!-----------------------------------------------------------------------
!
      REAL(r_typ) :: st
!
!-----------------------------------------------------------------------
!
      st = SIN(T)
      X = R*st*COS(P)
      Y = R*st*SIN(P)
      Z = R*COS(T)
!
      END
!*==DS_S2C.spg  processed by SPAG 6.72Dc at 22:10 on 27 May 2023
!#######################################################################
      FUNCTION DS_S2C(R,T,Drtp,Axisymmetric)
!
!-----------------------------------------------------------------------
!
! ****** Return the cell size in Cartesian coordinates at the
! ****** spherical position (R,T), calculated as the smallest
! ****** cell size from the mesh spacings in spherical coordinates
! ****** in DRTP.
!
! ****** In the axisymmetric case, AXISYMMETRIC=.true., the phi
! ****** mesh spacing is not present in DRTP, and is not considered.
!
! ****** This routine assumes that all angles are in radians.
!
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!
      USE NUMBER_TYPES
      IMPLICIT NONE
!*--DS_S2C1902
!
!-----------------------------------------------------------------------
!
      REAL(r_typ) :: R , T
      REAL(r_typ) , DIMENSION(3) :: Drtp
      LOGICAL :: Axisymmetric
      REAL(r_typ) :: DS_S2C
!
!-----------------------------------------------------------------------
!
      REAL(r_typ) :: dr , dt , dp , st
!
!-----------------------------------------------------------------------
!
      dr = Drtp(1)
      dt = Drtp(2)
!
      IF ( Axisymmetric ) THEN
         DS_S2C = MIN(dr,R*dt)
      ELSE
         dp = Drtp(3)
         st = MAX(SIN(T),SIN(dt))
         DS_S2C = MIN(dr,R*dt,R*st*dp)
      ENDIF
!
      END
!*==INTERP_FIELD_DS.spg  processed by SPAG 6.72Dc at 22:10 on 27 May 2023
!#######################################################################
      SUBROUTINE INTERP_FIELD_DS(Fld,Invtab,R,T,P,Fv,Ds,Sp2,Sp3)
!
!-----------------------------------------------------------------------
!
! ****** Interpolate the value of the field in SDS FLD at the point
! ****** (R,T,P) in spherical coordinates.
!
! ****** The value of the interpolated field is returned in FV.
!
! ****** The cell spacing at the interpolation point along each
! ****** dimension is returned in array DS.
!
! ****** If cubic spline interpolation is specified, it calls both
! ****** types (linear and cubic) because the spline library does
! ****** not return the info needed for mesh spacing, and I want to keep
! ****** the "outside" logic as is. There is probably a small speed
! ****** penalty, but it is not worth the complexity to speed it up.
! ****** Also cubic can return negative numbers --> floor fv to zero.
!
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!
      USE NUMBER_TYPES
      USE TYPES
      USE PARAMS , ONLY:CUBic
      USE SPLINE_DEF
      USE EVALUATE_SPLINE_2D_INTERFACE
      USE EVALUATE_SPLINE_3D_INTERFACE
      IMPLICIT NONE
!*--INTERP_FIELD_DS1962
!
!-----------------------------------------------------------------------
!
      RECORD /SDS   / ::Fld
      RECORD /VTAB  / ::Invtab
      REAL(r_typ) :: R , T , P
      REAL(r_typ) :: Fv
      REAL(r_typ) , DIMENSION(3) :: Ds
      RECORD /SPL2D / ::Sp2
      RECORD /SPL3D / ::Sp3
!
!-----------------------------------------------------------------------
!
! ****** Note that the interpolated value is set to zero if it is
! ****** outside the bounds of the mesh.
!
      IF ( Fld%NDIM.EQ.2 ) THEN
         CALL INTERP_2D(Fld%DIMS(1),Fld%DIMS(2),Fld%SCALES(1)%F,        &
                      & Fld%SCALES(2)%F,Invtab,Fld%F,R,T,Fv,Ds)
!
         IF ( CUBic .AND. (Fv.NE.0.) ) THEN
            Fv = EVALUATE_SPLINE_2D(Sp2,R,T,Invtab%C(1),Invtab%C(2))
            Fv = MAX(Fv,0.)
         ENDIF
!
      ELSEIF ( Fld%NDIM.EQ.3 ) THEN
         CALL INTERP_3D(Fld%DIMS(1),Fld%DIMS(2),Fld%DIMS(3),            &
                      & Fld%SCALES(1)%F,Fld%SCALES(2)%F,Fld%SCALES(3)%F,&
                      & Invtab,Fld%F,R,T,P,Fv,Ds)
         IF ( CUBic .AND. (Fv.NE.0.) ) THEN
            Fv = EVALUATE_SPLINE_3D(Sp3,R,T,P,Invtab%C(1),Invtab%C(2),  &
               & Invtab%C(3))
            Fv = MAX(Fv,0.)
         ENDIF
      ELSE
         Fv = 0.
         Ds = HUGE(Fv)
      ENDIF
!
      END
!*==INTERP_FIELD.spg  processed by SPAG 6.72Dc at 22:10 on 27 May 2023
!#######################################################################
      SUBROUTINE INTERP_FIELD(Fld,Invtab,R,T,P,Fv,Sp2,Sp3)
!
!-----------------------------------------------------------------------
!
! ****** Interpolate the value of the field in SDS FLD at the point
! ****** (R,T,P) in spherical coordinates.
!
! ****** The value of the interpolated field is returned in FV.
!
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!
      USE NUMBER_TYPES
      USE TYPES
      USE PARAMS , ONLY:CUBic
      USE SPLINE_DEF
      USE EVALUATE_SPLINE_2D_INTERFACE
      USE EVALUATE_SPLINE_3D_INTERFACE
      IMPLICIT NONE
!*--INTERP_FIELD2026
!
!-----------------------------------------------------------------------
!
      RECORD /SDS   / ::Fld
      RECORD /VTAB  / ::Invtab
      REAL(r_typ) :: R , T , P
      REAL(r_typ) :: Fv
      RECORD /SPL2D / ::Sp2
      RECORD /SPL3D / ::Sp3
!
      REAL(r_typ) , DIMENSION(3) :: ds
!
!-----------------------------------------------------------------------
!
      IF ( Fld%NDIM.EQ.2 ) THEN
         IF ( CUBic ) THEN
            Fv = EVALUATE_SPLINE_2D(Sp2,R,T,Invtab%C(1),Invtab%C(2))
         ELSE
            CALL INTERP_2D(Fld%DIMS(1),Fld%DIMS(2),Fld%SCALES(1)%F,     &
                         & Fld%SCALES(2)%F,Invtab,Fld%F,R,T,Fv,ds)
         ENDIF
!
      ELSEIF ( Fld%NDIM.EQ.3 ) THEN
         IF ( CUBic ) THEN
            Fv = EVALUATE_SPLINE_3D(Sp3,R,T,P,Invtab%C(1),Invtab%C(2),  &
               & Invtab%C(3))
         ELSE
            CALL INTERP_3D(Fld%DIMS(1),Fld%DIMS(2),Fld%DIMS(3),         &
                         & Fld%SCALES(1)%F,Fld%SCALES(2)%F,Fld%SCALES(3)&
                         & %F,Invtab,Fld%F,R,T,P,Fv,ds)
         ENDIF
      ENDIF
!
      END
!*==INTERP_2D.spg  processed by SPAG 6.72Dc at 22:10 on 27 May 2023
!
!#######################################################################
      SUBROUTINE INTERP_2D(Nx,Ny,X,Y,Inv,F,Xv,Yv,Fv,Ds)
!
!-----------------------------------------------------------------------
!
! ****** Interpolate the value of the 2D field FV at (XV,YV) from
! ****** array F(NX,NY), defined on the mesh X(NX) x Y(NY).
!
! ****** The structure INV holds the inverse interpolation tables.
!
! ****** If the point (XV,YV) is outside the bounds of the
! ****** X x Y mesh, FV=0. is returned.
!
! ****** The mesh spacing at the interpolation point in each
! ****** dimension is returned in array DS.
!
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!
      USE IDENT
      USE NUMBER_TYPES
      USE TYPES
      IMPLICIT NONE
!*--INTERP_2D2088
!
!-----------------------------------------------------------------------
!
      INTEGER :: Nx , Ny
      REAL(r_typ) , DIMENSION(Nx) :: X
      REAL(r_typ) , DIMENSION(Ny) :: Y
      RECORD /VTAB  / ::Inv
      REAL(r_typ) , DIMENSION(Nx,Ny) :: F
      REAL(r_typ) :: Xv , Yv , Fv
      REAL(r_typ) , DIMENSION(2) :: Ds
!
!-----------------------------------------------------------------------
!
      REAL(r_typ) , PARAMETER :: ONE = 1._R_TYP
!
!-----------------------------------------------------------------------
!
      LOGICAL :: outside
      INTEGER :: i , j , ip1 , jp1
      REAL(r_typ) :: ax , ay
!
!-----------------------------------------------------------------------
!
      REAL(r_typ) , EXTERNAL :: MESH_SPACING
!
!-----------------------------------------------------------------------
!
! ****** Find the cell that contains the interpolation point.
!
      outside = .FALSE.
!
      IF ( Xv.LT.X(1) ) THEN
         outside = .TRUE.
         i = 1
         ip1 = 1
         ax = 0.
      ELSEIF ( Xv.GT.X(Nx) ) THEN
         outside = .TRUE.
         i = Nx
         ip1 = Nx
         ax = 0.
      ELSE
         CALL INTERPI(X,Nx,Inv%C(1),Xv,i,ip1,ax)
      ENDIF
!
      IF ( Yv.LT.Y(1) ) THEN
         outside = .TRUE.
         j = 1
         jp1 = 1
         ay = 0.
      ELSEIF ( Yv.GT.Y(Ny) ) THEN
         outside = .TRUE.
         j = Ny
         jp1 = Ny
         ay = 0.
      ELSE
         CALL INTERPI(Y,Ny,Inv%C(2),Yv,j,jp1,ay)
      ENDIF
!
! ****** Get the mesh spacing at the interpolation point.
!
      Ds(1) = MESH_SPACING(Nx,X,i,ip1,ax)
      Ds(2) = MESH_SPACING(Ny,Y,j,jp1,ay)
!
! ****** If the point is outside the mesh limits, set the
! ****** interpolated value to 0.
!
      IF ( outside ) THEN
         Fv = 0.
         WRITE (*,*)
         WRITE (*,*) '### ERROR in ' , CNAME , ':'
         WRITE (*,*) '### Error in INTERP_2D:'
         WRITE (*,*) '### A requested location is outside the mesh.'
         WRITE (*,*) '### This should not happen anymore!'
         WRITE (*,*)
         WRITE (*,*) 'The requested location was:'
         WRITE (*,*) 'r: ' , Xv
         WRITE (*,*) 't: ' , Yv
         CALL EXIT(1)
      ELSE
         Fv = (ONE-ax)*((ONE-ay)*F(i,j)+ay*F(i,jp1))                    &
            & + ax*((ONE-ay)*F(ip1,j)+ay*F(ip1,jp1))
      ENDIF
!
      END
!*==INTERP_3D.spg  processed by SPAG 6.72Dc at 22:10 on 27 May 2023
!#######################################################################
      SUBROUTINE INTERP_3D(Nx,Ny,Nz,X,Y,Z,Inv,F,Xv,Yv,Zv,Fv,Ds)
!
!-----------------------------------------------------------------------
!
! ****** Interpolate the value of the 3D field FV at (XV,YV,ZV) from
! ****** array F(NX,NY,NZ), defined on the mesh X(NX) x Y(NY) x Z(NZ).
!
! ****** The structure INV holds the inverse interpolation tables.
!
! ****** If the point (XV,YV,ZV) is outside the bounds of the
! ****** X x Y x Z mesh, FV=0. is returned.
!
! ****** The mesh spacing at the interpolation point in each
! ****** dimension is returned in array DS.
!
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!
      USE IDENT
      USE NUMBER_TYPES
      USE TYPES
      IMPLICIT NONE
!*--INTERP_3D2200
!
!-----------------------------------------------------------------------
!
      INTEGER :: Nx , Ny , Nz
      REAL(r_typ) , DIMENSION(Nx) :: X
      REAL(r_typ) , DIMENSION(Ny) :: Y
      REAL(r_typ) , DIMENSION(Nz) :: Z
      RECORD /VTAB  / ::Inv
      REAL(r_typ) , DIMENSION(Nx,Ny,Nz) :: F
      REAL(r_typ) :: Xv , Yv , Zv , Fv
      REAL(r_typ) , DIMENSION(3) :: Ds
!
!-----------------------------------------------------------------------
!
      REAL(r_typ) , PARAMETER :: ONE = 1._R_TYP
!
!-----------------------------------------------------------------------
!
      LOGICAL :: outside
      INTEGER :: i , j , k , ip1 , jp1 , kp1
      REAL(r_typ) :: ax , ay , az
!
!-----------------------------------------------------------------------
!
      REAL(r_typ) , EXTERNAL :: MESH_SPACING
!
!-----------------------------------------------------------------------
!
! ****** Find the cell that contains the interpolation point.
!
      outside = .FALSE.
!
      IF ( Xv.LT.X(1) ) THEN
         outside = .TRUE.
         i = 1
         ip1 = 1
         ax = 0.
      ELSEIF ( Xv.GT.X(Nx) ) THEN
         outside = .TRUE.
         i = Nx
         ip1 = Nx
         ax = 0.
      ELSE
         CALL INTERPI(X,Nx,Inv%C(1),Xv,i,ip1,ax)
      ENDIF
!
      IF ( Yv.LT.Y(1) ) THEN
         outside = .TRUE.
         j = 1
         jp1 = 1
         ay = 0.
      ELSEIF ( Yv.GT.Y(Ny) ) THEN
         outside = .TRUE.
         j = Ny
         jp1 = Ny
         ay = 0.
      ELSE
         CALL INTERPI(Y,Ny,Inv%C(2),Yv,j,jp1,ay)
      ENDIF
!
      IF ( Zv.LT.Z(1) ) THEN
         outside = .TRUE.
         k = 1
         kp1 = 1
         az = 0.
      ELSEIF ( Zv.GT.Z(Nz) ) THEN
         outside = .TRUE.
         k = Nz
         kp1 = Nz
         az = 0.
      ELSE
         CALL INTERPI(Z,Nz,Inv%C(3),Zv,k,kp1,az)
      ENDIF
!
! ****** Get the mesh spacing at the interpolation point.
!
      Ds(1) = MESH_SPACING(Nx,X,i,ip1,ax)
      Ds(2) = MESH_SPACING(Ny,Y,j,jp1,ay)
      Ds(3) = MESH_SPACING(Nz,Z,k,kp1,az)
!
! ****** If the point is outside the mesh limits, set the
! ****** interpolated value to 0.
!
      IF ( outside ) THEN
         Fv = 0.
         WRITE (*,*)
         WRITE (*,*) '### ERROR in ' , CNAME , ':'
         WRITE (*,*) '### Error in INTERP_3D:'
         WRITE (*,*) '### A requested location is outside the mesh.'
         WRITE (*,*) '### This should not happen anymore!'
         WRITE (*,*)
         WRITE (*,*) 'The requested location was:'
         WRITE (*,*) 'r: ' , Xv
         WRITE (*,*) 't: ' , Yv
         WRITE (*,*) 'p: ' , Zv
         CALL EXIT(1)
      ELSE
         Fv = (ONE-ax)*((ONE-ay)*((ONE-az)*F(i,j,k)+az*F(i,j,kp1))      &
            & +ay*((ONE-az)*F(i,jp1,k)+az*F(i,jp1,kp1)))                &
            & + ax*((ONE-ay)*((ONE-az)*F(ip1,j,k)+az*F(ip1,j,kp1))      &
            & +ay*((ONE-az)*F(ip1,jp1,k)+az*F(ip1,jp1,kp1)))
      ENDIF
!
      END
!*==MESH_SPACING.spg  processed by SPAG 6.72Dc at 22:10 on 27 May 2023
!#######################################################################
      FUNCTION MESH_SPACING(Nx,X,I,Ip1,Alpha)
!
!-----------------------------------------------------------------------
!
! ****** Return the mesh spacing in coordinate X at the interpolation
! ****** point specified by I, IP1, and ALPHA.
!
! ****** This is an internal utility routine that does not do
! ****** any error checking.
!
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!
      USE NUMBER_TYPES
      IMPLICIT NONE
!*--MESH_SPACING2324
!
!-----------------------------------------------------------------------
!
      INTEGER :: Nx
      REAL(r_typ) , DIMENSION(Nx) :: X
      INTEGER :: I , Ip1
      REAL(r_typ) :: Alpha
      REAL(r_typ) :: MESH_SPACING
!
!-----------------------------------------------------------------------
!
      REAL(r_typ) , PARAMETER :: ONE = 1._R_TYP
!
!-----------------------------------------------------------------------
!
      INTEGER :: im , ip
      REAL(r_typ) :: d , dxm , dxp
!
!-----------------------------------------------------------------------
!
      im = MAX(1,I-1)
      ip = MIN(Nx,I+1)
      d = REAL(ip-im)
      IF ( d.NE.0. ) THEN
         dxm = (X(ip)-X(im))/d
      ELSE
         dxm = HUGE(ONE)
      ENDIF
!
      im = MAX(1,Ip1-1)
      ip = MIN(Nx,Ip1+1)
      d = REAL(ip-im)
      IF ( d.NE.0. ) THEN
         dxp = (X(ip)-X(im))/d
      ELSE
         dxp = HUGE(ONE)
      ENDIF
!
      MESH_SPACING = (ONE-Alpha)*dxm + Alpha*dxp
!
      END
!*==BUILD_INVERSE_TABLES.spg  processed by SPAG 6.72Dc at 22:10 on 27 May 2023
!#######################################################################
      SUBROUTINE BUILD_INVERSE_TABLES(S,Inv)
!
!-----------------------------------------------------------------------
!
! ****** Build the inverse interpolation tables INV for the SDS
! ****** in structure S.
!
! ****** These arrays are used to to increase the efficiency
! ****** of interpolation lookups.
!
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!
      USE NUMBER_TYPES
      USE TYPES
      IMPLICIT NONE
!*--BUILD_INVERSE_TABLES2386
!
!-----------------------------------------------------------------------
!
      RECORD /SDS   / ::S
      RECORD /VTAB  / ::Inv
!
!-----------------------------------------------------------------------
!
      INTEGER :: i
!
!-----------------------------------------------------------------------
!
! ****** Use a number of points for the inverse interpolation table
! ****** equal to the number in the original scale.
!
      DO i = 1 , S%NDIM
         Inv%C(i)%N = S%DIMS(i)
         ALLOCATE (Inv%C(i)%F(Inv%C(i)%N))
         CALL GETINV(S%SCALES(i)%F,S%DIMS(i),Inv%C(i))
      ENDDO
!
      END
!*==GETINV.spg  processed by SPAG 6.72Dc at 22:10 on 27 May 2023
!#######################################################################
      SUBROUTINE GETINV(X,N,Tab)
!
!-----------------------------------------------------------------------
!
! ****** Build an inverse interpolation table to increase the
! ****** efficiency of table look-up in a nonuniform mesh.
!
! ****** On input, the table X(N) is specified, together with the
! ****** number of points to use in the inverse interpolation
! ****** table, TAB%N.
!
! ****** The output is a structure TAB with the inverse interpolation
! ****** table.  This structure has the following components:
!
! ******    N:  the number of points in the table (input);
! ******    D:  the inverse of the uniform table spacing;
! ******    F:  the inverse interpolation table.
!
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!
      USE IDENT
      USE NUMBER_TYPES
      USE TYPES
      IMPLICIT NONE
!*--GETINV2438
!
!-----------------------------------------------------------------------
!
      INTEGER :: N
      REAL(r_typ) , DIMENSION(N) :: X
      RECORD /ITAB  / ::Tab
!
!-----------------------------------------------------------------------
!
      REAL(r_typ) , PARAMETER :: ONE = 1._R_TYP
!
!-----------------------------------------------------------------------
!
      INTEGER :: i , k , ierr , ip1
      REAL(r_typ) :: dx , xv , alpha , en
!
!-----------------------------------------------------------------------
!
! ****** For the special case when the table X has only one point,
! ****** the inverse table is not used.  It is thus loaded with
! ****** dummy values.
!
      IF ( N.EQ.1 ) THEN
         Tab%D = 0.
         Tab%F = 1.
         RETURN
      ENDIF
!
! ****** Check that the number of points is valid.
!
      IF ( Tab%N.LE.1 ) THEN
         WRITE (*,*)
         WRITE (*,*) '### ERROR in ' , CNAME , ':'
         WRITE (*,*) '### Error in GETINV:'
         WRITE (*,*) '### Invalid number of points specified'//         &
                    &' for the inverse interpolation table.'
         WRITE (*,*)
         WRITE (*,*) 'Number of points = ' , Tab%N
         CALL EXIT(1)
      ENDIF
!
! ****** Set the uniform interval to be used in the inverse
! ****** interpolation.
!
      dx = (X(N)-X(1))/(Tab%N-ONE)
!
      IF ( dx.LE.0. ) THEN
         WRITE (*,*)
         WRITE (*,*) '### ERROR in ' , CNAME , ':'
         WRITE (*,*) '### Error in GETINV:'
         WRITE (*,*) '### Invalid interval for the inverse'//           &
                    &' interpolation table.'
         WRITE (*,*)
         WRITE (*,*) 'Interval = ' , dx
         CALL EXIT(1)
      ENDIF
!
      Tab%D = ONE/dx
!
! ****** Build the inverse interpolation table.
!
      en = N
!
      DO k = 1 , Tab%N
         xv = X(1) + (k-ONE)*dx
         xv = MAX(xv,X(1))
         xv = MIN(xv,X(N))
         CALL INTERP(N,X,xv,i,ip1,alpha,ierr)
         IF ( ierr.NE.0 ) THEN
            WRITE (*,*)
            WRITE (*,*) '### ERROR in ' , CNAME , ':'
            WRITE (*,*) '### Error in GETINV:'
            WRITE (*,*) '### Error in building the inverse'//           &
                       &' interpolation table.'
            CALL EXIT(1)
         ENDIF
         Tab%F(k) = i + alpha
         Tab%F(k) = MAX(Tab%F(k),ONE)
         Tab%F(k) = MIN(Tab%F(k),en)
      ENDDO
!
      END
!*==INTERPI.spg  processed by SPAG 6.72Dc at 22:10 on 27 May 2023
!#######################################################################
      SUBROUTINE INTERPI(X,N,Tab,Xv,I,Ip1,Alpha)
!
!-----------------------------------------------------------------------
!
! ****** Get interpolation factor ALPHA and table indices I and IP1.
!
! ****** This routine does not do the actual interpolation.  Use the
! ****** returned values of I, IP1, and ALPHA to get the
! ****** interpolated value.
!
! ****** This version uses an inverse interpolation table, TAB,
! ****** to improve the efficiency of the search.
!
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!
      USE IDENT
      USE NUMBER_TYPES
      USE TYPES
      IMPLICIT NONE
!*--INTERPI2545
!
!-----------------------------------------------------------------------
!
      INTEGER :: N
      REAL(r_typ) , DIMENSION(N) :: X
      RECORD /ITAB  / ::Tab
      REAL(r_typ) :: Xv
      INTEGER :: I
      INTEGER :: Ip1
      REAL(r_typ) :: Alpha
      INTENT (IN)::X , N , Tab , Xv
      INTENT (OUT)::I , Ip1 , Alpha
!
!-----------------------------------------------------------------------
!
      REAL(r_typ) , PARAMETER :: ONE = 1._R_TYP
!
!-----------------------------------------------------------------------
!
      INTEGER :: ig
      REAL(r_typ) :: xi , fiv
!
!-----------------------------------------------------------------------
!
! ****** For the special case when the table has only one point,
! ****** the inverse table is not used.  In this case it is
! ****** necessary for XV to equal X(1) exactly, otherwise this
! ****** routine exits with an error.
!
      IF ( N.EQ.1 ) THEN
         IF ( Xv.NE.X(1) ) GOTO 100
         I = 1
         Ip1 = 1
         Alpha = 0.
      ENDIF
!
! ****** Get an estimate of the nearest grid point location in
! ****** the (uniform) inverse interpolation table.
!
      xi = ONE + (Xv-X(1))*Tab%D
      I = xi
      I = MAX(I,1)
      I = MIN(I,Tab%N-1)
      Alpha = xi - I
      fiv = (ONE-Alpha)*Tab%F(I) + Alpha*Tab%F(I+1)
!
! ****** Set IG to be the guess for the nearest grid point.
!
      ig = fiv
      ig = MAX(ig,1)
      ig = MIN(ig,N-1)
!
      IF ( Xv.GE.X(ig) ) THEN
!
! ****** Search forwards.
!
         DO I = ig , N - 1
            IF ( Xv.GE.X(I) .AND. Xv.LE.X(I+1) ) THEN
               Alpha = (Xv-X(I))/(X(I+1)-X(I))
               Ip1 = I + 1
               RETURN
            ENDIF
         ENDDO
!
      ELSE
!
! ****** Search backwards.
!
         DO I = ig - 1 , 1 , -1
            IF ( Xv.GE.X(I) .AND. Xv.LE.X(I+1) ) THEN
               Alpha = (Xv-X(I))/(X(I+1)-X(I))
               Ip1 = I + 1
               RETURN
            ENDIF
         ENDDO
!
      ENDIF
!
! ****** ERROR: value not found in table.
!
!
 100  WRITE (*,*)
      WRITE (*,*) '### ERROR in ' , CNAME , ':'
      WRITE (*,*) '### Error in INTERPI:'
      WRITE (*,*) 'Abscissa not found in table.'
      WRITE (*,*)
      WRITE (*,*) 'N = ' , N
      WRITE (*,*) 'X = ' , (X(I),I=1,N)
      WRITE (*,*)
      WRITE (*,*) 'Abscissa requested = ' , Xv
      CALL EXIT(1)
!
      END
!*==INTERP.spg  processed by SPAG 6.72Dc at 22:10 on 27 May 2023
!#######################################################################
      SUBROUTINE INTERP(N,X,Xv,I,Ip1,A,Ierr)
!
!-----------------------------------------------------------------------
!
! ****** Get the interpolation factor at XV from the table X(N).
!
!-----------------------------------------------------------------------
!
! ****** This routine does not do the actual interpolation.  Use the
! ****** returned values of I, IP1, and A to get the interpolated
! ****** value.
!
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!
      USE NUMBER_TYPES
      IMPLICIT NONE
!*--INTERP2660
!
!-----------------------------------------------------------------------
!
      INTEGER :: N
      REAL(r_typ) , DIMENSION(N) :: X
      REAL(r_typ) :: Xv
      INTEGER :: I , Ip1
      REAL(r_typ) :: A
      INTEGER :: Ierr
!
!-----------------------------------------------------------------------
!
      Ierr = 0
!
! ****** For the special case when the table has only one point,
! ****** it is necessary for XV to equal X(1) exactly, otherwise
! ****** this routine exits with an error.
!
      IF ( N.EQ.1 ) THEN
         IF ( Xv.NE.X(1) ) GOTO 100
         I = 1
         Ip1 = 1
         A = 0.
         RETURN
      ENDIF
!
! ****** Find the interval and compute the interpolation factor.
!
      DO I = 1 , N - 1
         IF ( Xv.GE.X(I) .AND. Xv.LE.X(I+1) ) THEN
            Ip1 = I + 1
            IF ( X(I).EQ.X(I+1) ) THEN
               A = 0.
            ELSE
               A = (Xv-X(I))/(X(I+1)-X(I))
            ENDIF
            RETURN
         ENDIF
      ENDDO
!
! ****** ERROR: the value was not found.
!
!
 100  WRITE (*,*)
      WRITE (*,*) '### ERROR in INTERP:'
      WRITE (*,*) '### Value not found in table.'
      WRITE (*,*) 'Value requested = ' , Xv
      WRITE (*,*) 'Min table value = ' , X(1)
      WRITE (*,*) 'Max table value = ' , X(N)
      WRITE (*,*) 'Number of values in table = ' , N
      Ierr = 1
!
      END
!*==RTP2XYZ.spg  processed by SPAG 6.72Dc at 22:10 on 27 May 2023
!#######################################################################
      SUBROUTINE RTP2XYZ(R,T,P,X,Y,Z)
!
!-----------------------------------------------------------------------
!
! ****** Convert spherical position to cartesian position.
!
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!
      USE NUMBER_TYPES
      IMPLICIT NONE
!*--RTP2XYZ2731
!
!-----------------------------------------------------------------------
!
      REAL(r_typ) :: R , T , P
      REAL(r_typ) :: X , Y , Z
!
!-----------------------------------------------------------------------
!
      REAL(r_typ) :: cph , sph , cth , sth
!
!-----------------------------------------------------------------------
!
      cph = COS(P)
      sph = SIN(P)
      cth = COS(T)
      sth = SIN(T)
      X = R*cph*sth
      Y = R*sph*sth
      Z = R*cth
!
      END
!*==SVTOCV.spg  processed by SPAG 6.72Dc at 22:10 on 27 May 2023
!#######################################################################
      SUBROUTINE SVTOCV(T,P,Vr,Vt,Vp,Vx,Vy,Vz)
!
!-----------------------------------------------------------------------
!
! ****** Rotate Spherical Vector Components to Cartesian Components
!
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!
      USE NUMBER_TYPES
      IMPLICIT NONE
!*--SVTOCV2770
!
!-----------------------------------------------------------------------
!
      REAL(r_typ) :: T , P
      REAL(r_typ) :: Vr , Vt , Vp
      REAL(r_typ) :: Vx , Vy , Vz
!
!-----------------------------------------------------------------------
!
      REAL(r_typ) :: cph , sph , cth , sth
!
!-----------------------------------------------------------------------
!
      cph = COS(P)
      sph = SIN(P)
      cth = COS(T)
      sth = SIN(T)
      Vx = Vr*sth*cph + Vt*cth*cph - Vp*sph
      Vy = Vr*sth*sph + Vt*cth*sph + Vp*cph
      Vz = Vr*cth - Vt*sth
!
      END
!*==NORMALIZE.spg  processed by SPAG 6.72Dc at 22:10 on 27 May 2023
!#######################################################################
      SUBROUTINE NORMALIZE(Vec_in,Vec_out)
!
!-----------------------------------------------------------------------
!
! ****** Normalize a vector to unit length.
!
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!
      USE NUMBER_TYPES
      IMPLICIT NONE
!*--NORMALIZE2808
!
!-----------------------------------------------------------------------
!
      REAL(r_typ) , DIMENSION(3) :: Vec_out
      REAL(r_typ) , DIMENSION(3) :: Vec_in
      REAL(r_typ) :: mag
!
!-----------------------------------------------------------------------
!
      mag = SQRT(Vec_in(1)**2+Vec_in(2)**2+Vec_in(3)**2)
      Vec_out(1) = Vec_in(1)/mag
      Vec_out(2) = Vec_in(2)/mag
      Vec_out(3) = Vec_in(3)/mag
!
      END
!*==CROSS.spg  processed by SPAG 6.72Dc at 22:10 on 27 May 2023
!#######################################################################
      SUBROUTINE CROSS(Vec1,Vec2,Vec_out)
!
!-----------------------------------------------------------------------
!
! ****** Compute the cross product of two vectors.
!
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!
      USE NUMBER_TYPES
      IMPLICIT NONE
!*--CROSS2839
!
!-----------------------------------------------------------------------
!
      REAL(r_typ) , DIMENSION(3) :: Vec_out
      REAL(r_typ) , DIMENSION(3) :: Vec1 , Vec2
!
!-----------------------------------------------------------------------
!
      Vec_out(1) = Vec1(2)*Vec2(3) - Vec1(3)*Vec2(2)
      Vec_out(2) = Vec1(3)*Vec2(1) - Vec1(1)*Vec2(3)
      Vec_out(3) = Vec1(1)*Vec2(2) - Vec1(2)*Vec2(1)
!
      END
!*==GET_VECTOR_PROJECTION.spg  processed by SPAG 6.72Dc at 22:10 on 27 May 2023
!#######################################################################
      SUBROUTINE GET_VECTOR_PROJECTION(Ex,Ey,Ez,Rv,Tv,Pv,Vec)
!
!-----------------------------------------------------------------------
!
! ****** Get the vector field at this location and project it into
! ****** the observing frame defined with unit vectors ex, ey, ez.
!
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!
      USE NUMBER_TYPES
      USE MAS_FIELDS
      IMPLICIT NONE
!*--GET_VECTOR_PROJECTION2870
!
!-----------------------------------------------------------------------
!
      REAL(R_Typ) , DIMENSION(3) :: Ex , Ey , Ez
      REAL(R_Typ) :: Rv , Tv , Pv
      REAL(R_Typ) , DIMENSION(3) :: Vec
!
!-----------------------------------------------------------------------
!
      REAL(R_Typ) :: vrv , vtv , vpv
      REAL(R_Typ) :: vx , vy , vz
!
!-----------------------------------------------------------------------
!
! ****** Interpolate the velocity components.
!
      CALL INTERP_FIELD(VR%SDS,VR%INVTAB,Rv,Tv,Pv,vrv,VR%SPL2,VR%SPL3)
      CALL INTERP_FIELD(VT%SDS,VT%INVTAB,Rv,Tv,Pv,vtv,VT%SPL2,VT%SPL3)
      CALL INTERP_FIELD(VP%SDS,VP%INVTAB,Rv,Tv,Pv,vpv,VP%SPL2,VP%SPL3)
!
! ****** Get the cartesian vector components.
!
      CALL SVTOCV(Tv,Pv,vrv,vtv,vpv,vx,vy,vz)
!
! ****** Project the vector into the observer frame.
!
      Vec(1) = Ex(1)*vx + Ex(2)*vy + Ex(3)*vz
      Vec(2) = Ey(1)*vx + Ey(2)*vy + Ey(3)*vz
      Vec(3) = Ez(1)*vx + Ez(2)*vy + Ez(3)*vz
!
      END
!*==GET_UNIT_VECTORS.spg  processed by SPAG 6.72Dc at 22:10 on 27 May 2023
!#######################################################################
      SUBROUTINE GET_UNIT_VECTORS(Losx,Losy,Losz,Ex,Ey,Ez)
!
!-----------------------------------------------------------------------
!
! ****** Compute the unit vectors of the observer frame.
!
! ****** los_vec is a vector pointing to the observer from a point
! ****** along the LOS (doesn't need to be normalized).
!
! ****** ex, ey, ez are the unit vectors of the observer frame.
!
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!
      USE NUMBER_TYPES
      IMPLICIT NONE
!*--GET_UNIT_VECTORS2922
!
!-----------------------------------------------------------------------
!
      REAL(r_typ) , PARAMETER :: ZERO = 0._R_TYP
      REAL(r_typ) , PARAMETER :: ONE = 1._R_TYP
!
!-----------------------------------------------------------------------
!
      REAL(r_typ) , INTENT(IN) :: Losx , Losy , Losz
      REAL(r_typ) , DIMENSION(3) , INTENT(OUT) :: Ex , Ey , Ez
!
!-----------------------------------------------------------------------
!
      REAL(r_typ) , DIMENSION(3) :: vec
!
!-----------------------------------------------------------------------
!
! ****** x-hat is sun-to-observer unit vector --> normalize los_vec.
!
      vec(1) = Losx
      vec(2) = Losy
      vec(3) = Losz
      CALL NORMALIZE(vec,Ex)
!
! ****** y-hat is perp to heliographic z axis and observer x-hat.
!
      vec(1) = ZERO
      vec(2) = ZERO
      vec(3) = ONE
      CALL CROSS(vec,Ex,Ey)
      CALL NORMALIZE(Ey,Ey)
!
! ****** z-hat is perp to heliographic x-hat and y-hat.
!
      CALL CROSS(Ex,Ey,Ez)
      CALL NORMALIZE(Ez,Ez)
!
      END
!*==ADD_TO_INTEGRALS.spg  processed by SPAG 6.72Dc at 22:10 on 27 May 2023
!#######################################################################
      SUBROUTINE ADD_TO_INTEGRALS(I,J,E_pb,E_b,Ds,Ex,Ey,Ez,Rv,Tv,Pv,X)
!
!-----------------------------------------------------------------------
!
! ****** Add the current contribution to the weighted integrals.
!
! ****** Everything is here so it can be used in both integration
! ****** algoriths and skipped if not requested.
!
! ****** i, j are the current pixel indexes.
!
! ****** e_pb and e_b are the local pB and B emissivities.
!
! ****** ds is the path length contribution for the integral.
!
! ****** ex, ey, ez are the unit vectors of the observer frame.
!
! ****** rv, tv, pv are the current coordinates
!
! ****** x is the signed distance from r_min along the LOS (positive is
! ****** towards the observer, negative is away.
!
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!
      USE NUMBER_TYPES
      USE PARAMS
      USE CONSTANTS
      USE IMAGE_FIELDS
      USE MAS_FIELDS
      IMPLICIT NONE
!*--ADD_TO_INTEGRALS2996
!
!-----------------------------------------------------------------------
!
      INTEGER , INTENT(IN) :: I , J
      REAL(R_Typ) , INTENT(IN) :: E_pb , E_b , Ds
      REAL(R_Typ) , DIMENSION(3) , INTENT(IN) :: Ex , Ey , Ez
      REAL(R_Typ) , INTENT(IN) :: Rv , Tv , Pv , X
!
!-----------------------------------------------------------------------
!
      REAL(R_Typ) :: var_now
      REAL(R_Typ) :: e_now
      REAL(R_Typ) :: angle_now
      REAL(R_Typ) , DIMENSION(3) :: vec
!
!-----------------------------------------------------------------------
!
! ****** Select the weighting.
!
      IF ( WEIght_integral_by_b ) THEN
         e_now = E_b
      ELSE
         e_now = E_pb
      ENDIF
!
! ****** Scalar field.
!
      IF ( COMpute_scalar_integration ) THEN
         CALL INTERP_FIELD(SFIeld%SDS,SFIeld%INVTAB,Rv,Tv,Pv,var_now,   &
                         & SFIeld%SPL2,SFIeld%SPL3)
         SF_avg(I,J) = SF_avg(I,J) + Ds*e_now*var_now
      ENDIF
!
! ****** Average LOS angle.
!
      IF ( COMpute_angle_integration ) THEN
         angle_now = ATAN(X/R_Min(I,J))/DEGTORAD
         ANGle_avg(I,J) = ANGle_avg(I,J) + Ds*e_now*angle_now
      ENDIF
!
! ****** Vector field. Mapping from observer frame to image frame is:
!        Observer_x -> Image_LOS
!        Observer_y -> Image_x
!        Observer_z -> Image_y
!
      IF ( COMpute_vector_integration ) THEN
         CALL GET_VECTOR_PROJECTION(Ex,Ey,Ez,Rv,Tv,Pv,vec)
         VLOs_avg(I,J) = VLOs_avg(I,J) + Ds*e_now*vec(1)
         VX_avg(I,J) = VX_avg(I,J) + Ds*e_now*vec(2)
         VY_avg(I,J) = VY_avg(I,J) + Ds*e_now*vec(3)
      ENDIF
!
      END
!*==NORMALIZE_INTEGRALS.spg  processed by SPAG 6.72Dc at 22:10 on 27 May 2023
!#######################################################################
      SUBROUTINE NORMALIZE_INTEGRALS
!
!-----------------------------------------------------------------------
!
! ****** Normalize the weighted integrals to get averages.
!
! ****** Pixels that are occulted or outside the grid get the disk val.
!
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!
      USE NUMBER_TYPES
      USE PARAMS
      USE GEOMETRY
      USE IMAGE_REGION
      USE IMAGE_FIELDS
      IMPLICIT NONE
!*--NORMALIZE_INTEGRALS3071
!
!-----------------------------------------------------------------------
!
      REAL(R_Typ) , PARAMETER :: ONE = 1._R_TYP
!
!-----------------------------------------------------------------------
!
      INTEGER :: i , j
      REAL(R_Typ) :: e_now
!
!-----------------------------------------------------------------------
!
! ****** Loop over pixels, normalize the integrals.
!
      DO j = 1 , NY
         DO i = 1 , NX

            IF ( R_Min(i,j).LE.MAX(ROCc,ONE) .OR. R_Min(i,j).GE.RMAx )  &
               & THEN

               IF ( COMpute_scalar_integration ) SF_avg(i,j)            &
                  & = DISk_value
               IF ( COMpute_angle_integration ) ANGle_avg(i,j)          &
                  & = DISk_value
               IF ( COMpute_vector_integration ) THEN
                  VLOs_avg(i,j) = DISk_value
                  VX_avg(i,j) = DISk_value
                  VY_avg(i,j) = DISk_value
               ENDIF
               GOTO 50
            ENDIF
!
! ****** Divide out emissivity integral to get averages.
!
            IF ( DO_integral_weighting ) THEN
               IF ( WEIght_integral_by_b ) THEN
                  e_now = B(i,j)
               ELSE
                  e_now = PB(i,j)
               ENDIF
!
               IF ( COMpute_scalar_integration ) SF_avg(i,j)            &
                  & = SF_avg(i,j)/e_now
!
               IF ( COMpute_angle_integration ) ANGle_avg(i,j)          &
                  & = ANGle_avg(i,j)/e_now
!
               IF ( COMpute_vector_integration ) THEN
                  VLOs_avg(i,j) = VLOs_avg(i,j)/e_now
                  VX_avg(i,j) = VX_avg(i,j)/e_now
                  VY_avg(i,j) = VY_avg(i,j)/e_now
               ENDIF
            ENDIF
 50      ENDDO
      ENDDO
!
      END
!*==SET_PARAMETERS.spg  processed by SPAG 6.72Dc at 22:10 on 27 May 2023
!#######################################################################
      SUBROUTINE SET_PARAMETERS
!
!-----------------------------------------------------------------------
!
! ****** Set parameters from the command line arguments.
!
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!
      USE IDENT
      USE NUMBER_TYPES
      USE SYNTAX
      USE PARAGRAPH_DEF
      USE GET_USAGE_LINE_INTERFACE
      USE PRINT_PAR_INTERFACE
      USE DELETE_PAR_INTERFACE
      USE PARAMS
      USE IMAGE_REGION
      IMPLICIT NONE
!*--SET_PARAMETERS3152
!
!-----------------------------------------------------------------------
!
! ****** Storage the for usage line.
!
      RECORD /PARAGRAPH/  , POINTER::usage
!
! ****** Storage for the error message.
!
      CHARACTER(72) :: errmsg
!
!-----------------------------------------------------------------------
!
      INTEGER :: ierr , nargs
      CHARACTER(512) :: arg
      LOGICAL :: set
      LOGICAL , SAVE :: print_help = .FALSE.

!
!-----------------------------------------------------------------------
!
      INTEGER , EXTERNAL :: INTVAL
      REAL(R_Typ) , EXTERNAL :: FPVAL
!
!-----------------------------------------------------------------------
!
! ****** Define the syntax.
!
      CALL DEFARG(group_k,'-v',' ',' ')
      CALL DEFARG(group_k,'-help',' ',' ')
      CALL DEFARG(group_k,'-cubic',' ',' ')
      CALL DEFARG(group_ka,'-long','0','<deg>')
      CALL DEFARG(group_ka,'-p','0','<deg>')
      CALL DEFARG(group_ka,'-b0','0','<deg>')
      CALL DEFARG(group_ka,'-r','<none>','<AU>')
      CALL DEFARG(group_ka,'-nx','201','<int>')
      CALL DEFARG(group_ka,'-ny','201','<int>')
      CALL DEFARG(group_ka,'-x0','-3','<val>')
      CALL DEFARG(group_ka,'-x1','3','<val>')
      CALL DEFARG(group_ka,'-y0','-3','<val>')
      CALL DEFARG(group_ka,'-y1','3','<val>')
      CALL DEFARG(group_k,'-wispr1',' ',' ')
      CALL DEFARG(group_k,'-wispr2',' ',' ')
      CALL DEFARG(group_ka,'-dsmult','1','<factor>')
      CALL DEFARG(group_ka,'-rocc','1','<radius>')
      CALL DEFARG(group_ka,'-power','0','<r_power>')
      CALL DEFARG(group_ka,'-disk','-1','<val>')
      CALL DEFARG(group_k,'-oldmas',' ',' ')
      CALL DEFARG(group_ka,'-vf','<none>','<file>')
      CALL DEFARG(group_ka,'-he_frac','0','<val>')
      CALL DEFARG(group_ka,'-mu','0.63','<val>')
      CALL DEFARG(group_ka,'-rho',' ','<file>')
      CALL DEFARG(group_ka,'-scalar','<none>','<file>')
      CALL DEFARG(group_ka,'-vr','<none>','<file>')
      CALL DEFARG(group_ka,'-vt','<none>','<file>')
      CALL DEFARG(group_ka,'-vp','<none>','<file>')
      CALL DEFARG(group_ka,'-avg_scalar','<none>','<file>')
      CALL DEFARG(group_ka,'-avg_vlos','<none>','<file>')
      CALL DEFARG(group_ka,'-avg_vx','<none>','<file>')
      CALL DEFARG(group_ka,'-avg_vy','<none>','<file>')
      CALL DEFARG(group_ka,'-avg_los_angle','<none>','<file>')
      CALL DEFARG(group_k,'-avg_using_b',' ',' ')
      CALL DEFARG(group_ka,'-pb','<none>','<file>')
      CALL DEFARG(group_ka,'-b','<none>','<file>')
!
! ****** Parse the command line.
!
      CALL PARSE(errmsg,ierr)
!
      IF ( ierr.NE.0 ) THEN
!
         WRITE (*,*)
         WRITE (*,*) '### ' , CNAME , ' Version ' , CVERS , ' of ' ,    &
                   & CDATE , '.'
         WRITE (*,*) '### Calculate polarization brightness and'//      &
                    &' brightness.'
!
         IF ( ierr.GT.1 ) THEN
            WRITE (*,*)
            WRITE (*,*) errmsg
         ENDIF
!
! ****** Check if the "-help" option was specified.
!
         CALL NARGS_SPECIFIED(nargs)
         CALL FETCHARG('-help',set,arg)
!
         IF ( nargs.EQ.1 .AND. set ) THEN
            print_help = .TRUE.
         ELSEIF ( ierr.GT.1 ) THEN
            print_help = .FALSE.
            WRITE (*,*)
            WRITE (*,*) TRIM(errmsg)
         ENDIF
!
! ****** Print the usage line.
!
         CALL GET_USAGE_LINE(usage)
!
         WRITE (*,*)
         WRITE (*,*) 'Usage:'
         WRITE (*,*)
!
         CALL PRINT_PAR(usage)
!
         IF ( print_help ) CALL PRINT_HELP_TEXT
!
         WRITE (*,*)
         WRITE (*,*) 'The density is read in from the HDF file'//       &
                    &' specified by -rho <file>.'
         WRITE (*,*) 'To compute pB, specify the output HDF'//          &
                    &' file by setting -pb <file>.'
         WRITE (*,*) 'To compute B, specify the output HDF'//           &
                    &' file by setting -b <file>.'
         WRITE (*,*)
         WRITE (*,*) '-long sets the longitude of the central'//        &
                    &' meridian in degrees [default=0].'
         WRITE (*,*) '-p sets the solar P angle in degrees'//           &
                    &' [default=0].'
         WRITE (*,*) '-b0 sets the solar B0 angle in degrees'//         &
                    &' [default=0].'
         WRITE (*,*)
         WRITE (*,*) '-r sets the distance of the observer from'//      &
                    &' the Sun [AU]. By default, the'
         WRITE (*,*) 'observer is assumed to be infinitely far'//       &
                    &' away from the Sun. Setting this'
         WRITE (*,*) 'distance switches the computation from'//         &
                    &' parallel lines of sight to non-parallel'
         WRITE (*,*) 'lines of sight.'
         WRITE (*,*)
         WRITE (*,*) 'Use -x0, -x1, -y0, and -y1 to set the image'//    &
                    &' region. When the observer distance'
         WRITE (*,*) 'has been specified (non-parallel LOSs),'//        &
                    &' -x0, -x1 set the elongation [deg], and'
         WRITE (*,*) '-y0, -y1 set the altitude from the'//             &
                    &' ecliptic plane [deg]. Otherwise (parallel'
         WRITE (*,*) 'LOSs), -x0, -x1 and -y0, -y1 set the'//           &
                    &' horizontal and vertical plane-of-sky'
         WRITE (*,*) 'coordinates [solar radii]. The default is'//      &
                    &' the region [-3,3] x [-3,3] in both'
         WRITE (*,*) 'cases. The number of mesh points in the'//        &
                    &' image can be set using -nx and -ny'
         WRITE (*,*) '[default=201].'
         WRITE (*,*)
         WRITE (*,*) 'The flag -wispr1 is an abbreviation for'//        &
                    &' the region [13.5,51.4] x [-20,20]'
         WRITE (*,*) 'degrees, appropriate to camera 1 of'//            &
                    &' WISPR on Parker Solar Probe.'
         WRITE (*,*)
         WRITE (*,*) 'The flag -wispr2 is an abbreviation for'//        &
                    &' the region [49.7,104.7] x [-29,29]'
         WRITE (*,*) 'degrees, appropriate to camera 2 of WISPR.'
         WRITE (*,*)
         WRITE (*,*) '-vf specifies a 1D HDF file that defines'//       &
                    &' a radial vignetting function by which'
         WRITE (*,*) 'to multiply pB and B [default=none]. Use'//       &
                    &' -power to set an additional radial'
         WRITE (*,*) 'power by which to multiply pB and B'//            &
                    &' [default=0]. These functions use the'
         WRITE (*,*) 'radius R_MIN, the closest distance to the'//      &
                    &' solar surface along the LOS. (In the'
         WRITE (*,*) 'case of parallel LOSs, this is identical'//       &
                    &' to the plane-of-sky radius.)'
         WRITE (*,*)
         WRITE (*,*) '-rocc sets the occulting disk radius in'//        &
                    &' units of the solar radius [default=1].'
         WRITE (*,*) '-disk sets the value of pB and B for'//           &
                    &' LOSs that intersect the solar disk or'
         WRITE (*,*) 'occulting disk [default=-1]. -dsmult sets'//      &
                    &' the factor by which to multiply the'
         WRITE (*,*) 'LOS integration step size [default=1].'//         &
                    &' The step size is varied to match the'
         WRITE (*,*) 'local mesh in the density file.'
         WRITE (*,*)
         WRITE (*,*) 'Specify -oldmas when using data files from'//     &
                    &' the old MAS code.'
         WRITE (*,*)
         WRITE (*,*) 'Use -he_frac to specify a non-zero helium'//      &
                    &' fraction when getting ne from rho.'
         WRITE (*,*)
         WRITE (*,*) 'Use -mu to specify a custom limb-darkening'//     &
                    &' coefficient (default 0.63).'
         WRITE (*,*)
         WRITE (*,*) 'Specify -cubic to use cubic spline interpolation' &
                   & //' for rho and vf.'
         WRITE (*,*)
         WRITE (*,*) 'Specify -help for information on the special'//   &
                    &' integral averaging options:'
         WRITE (*,*) ' -scalar, -vr, -vt, -vp, -avg_scalar,'//          &
                    &' -avg_los_angle, etc.'
         WRITE (*,*)
!
         CALL DELETE_PAR(usage)
!
         CALL EXIT(1)
!
      ENDIF
!
! ****** Set the parameters.
!
! ****** Verbose flag.
!
      CALL FETCHARG('-v',set,arg)
      VERbose = set
!
! ****** Cubic interpolation flag.
!
      CALL FETCHARG('-cubic',set,arg)
      CUBic = set
!
! ****** Switch to indicate old MAS code data files.
!
      CALL FETCHARG('-oldmas',set,arg)
      OLDmas = set
!
! ****** Carrington rotation longitude.
!
      CALL FETCHARG('-long',set,arg)
      CRLong = FPVAL(arg,'-long')
!
! ****** Solar P angle.
!
      CALL FETCHARG('-p',set,arg)
      PANgle = FPVAL(arg,'-p')
!
! ****** Solar B0 angle.
!
      CALL FETCHARG('-b0',set,arg)
      B0Angle = FPVAL(arg,'-b0')
!
! ****** Distance of the observer from the Sun [AU].
!
      CALL FETCHARG('-r',set,arg)
      IF ( set ) THEN
         R_Obs_set = .TRUE.
         R_Obs = FPVAL(arg,'-r')
      ELSE
         R_Obs_set = .FALSE.
         R_Obs = 1._R_TYP
      ENDIF
!
! ****** Number of points to use for the image.
!
      CALL FETCHARG('-nx',set,arg)
      NX = INTVAL(arg,'-nx')
!
      CALL FETCHARG('-ny',set,arg)
      NY = INTVAL(arg,'-ny')
!
! ****** Image dimensions.
!
      IMAge_limits_set = .FALSE.
!
      CALL FETCHARG('-x0',set,arg)
      IF ( set ) IMAge_limits_set = .TRUE.
      X0 = FPVAL(arg,'-x0')
!
      CALL FETCHARG('-x1',set,arg)
      IF ( set ) IMAge_limits_set = .TRUE.
      X1 = FPVAL(arg,'-x1')
!
      CALL FETCHARG('-y0',set,arg)
      IF ( set ) IMAge_limits_set = .TRUE.
      Y0 = FPVAL(arg,'-y0')
!
      CALL FETCHARG('-y1',set,arg)
      IF ( set ) IMAge_limits_set = .TRUE.
      Y1 = FPVAL(arg,'-y1')
!
! ****** WISPR camera 1 flag.
!
      CALL FETCHARG('-wispr1',set,arg)
      WISpr1 = set
!
! ****** WISPR camera 2 flag.
!
      CALL FETCHARG('-wispr2',set,arg)
      WISpr2 = set
!
! ****** Occulting disk radius.
!
      CALL FETCHARG('-rocc',set,arg)
      ROCc = FPVAL(arg,'-rocc')
!
! ****** Factor by which to multiply the LOS integration step size.
!
      CALL FETCHARG('-dsmult',set,arg)
      DSMult = FPVAL(arg,'-dsmult')
!
! ****** Radial power exponent.
!
      CALL FETCHARG('-power',set,arg)
      POWer_set = set
      POWer = FPVAL(arg,'-power')
!
! ****** Value for pB on the solar disk.
!
      CALL FETCHARG('-disk',set,arg)
      DISk_value = FPVAL(arg,'-disk')
!
! ****** Vignetting function file name.
!
      VF_file = TRIM(py_vf)
      IF ( VF_file.EQ.'<none>' ) VF_file = ' '
!
! ****** Density file name.
!
      RHO_file = TRIM(py_rho)
!
! ****** Vector field r component file name (input).
!
      VR_file = TRIM(py_vr)
      IF ( VF_file.EQ.'<none>' ) THEN
         VR_file = ' '
      ELSE
         COMpute_vector_integration = .TRUE.
         VR_file = TRIM(py_vr)
      ENDIF
!
! ****** Vector field t component file name (input).
!
      VT_file = TRIM(py_vt)
      IF ( VT_file.EQ.'<none>' ) THEN
         VT_file = ' '
      ELSE
         COMpute_vector_integration = .TRUE.
         VT_file = TRIM(py_vt)
      ENDIF

!
! ****** Vector field p component file name (input).
!
      VP_file = TRIM(py_vp)
      IF ( VP_file.EQ.'<none>' ) THEN
         VP_file = ' '
      ELSE
         COMpute_vector_integration = .TRUE.
         VP_file = TRIM(py_vp)
      ENDIF
!
! ****** Scalar field file name (input).
!
      SCAlar_file = TRIM(py_scalar)
      IF ( SCAlar_file.EQ.'<none>' ) THEN
         COMpute_scalar_integration = .FALSE.
         SCAlar_file = ' '
      ELSE
         COMpute_vector_integration = .TRUE.
         SCAlar_file = TRIM(py_scalar)
      ENDIF

!
! ****** Emissivity weighted LOS integral of scalar field.
!
      SCAlar_file_out = TRIM(py_avg_scalar)
      IF ( SCAlar_file_out.EQ.'<none>' ) THEN
         SCAlar_file_out = ' '
      ELSE
         SCAlar_file_out = TRIM(py_avg_scalar)
      ENDIF
!
! ****** Emissivity weighted LOS contribution angle.
!
      ANGle_file_out = TRIM(py_avg_los_angle)
      IF ( ANGle_file_out.EQ.'<none>' ) THEN
         ANGle_file_out = ' '
      ELSE
         COMpute_vector_integration = .TRUE.
         ANGle_file_out = TRIM(py_avg_los_angle)
      ENDIF

!
! ****** Emissivity weighted LOS integral of vector field (vlos).
!
      VLOs_file_out = TRIM(py_avg_vlos)
      IF ( VLOs_file_out.EQ.'<none>' ) THEN
         VLOs_file_out = ' '
      ELSE
         VLOs_file_out = TRIM(py_avg_vlos)
      ENDIF
!
! ****** Emissivity weighted LOS integral of vector field (vx).
!
      VX_file_out = TRIM(py_avg_vx)
      IF ( VX_file_out.EQ.'<none>' ) THEN
         VX_file_out = ' '
      ELSE
         VX_file_out = TRIM(py_avg_vx)
      ENDIF

!
! ****** Emissivity weighted LOS integral of vector field (vy).
!
      VY_file_out = TRIM(py_avg_vy)
      IF ( VY_file_out.EQ.'<none>' ) THEN
         VY_file_out = ' '
      ELSE
         VY_file_out = TRIM(py_avg_vy)
      ENDIF

!
! ****** Flag to weight integral averages with B and not pB.
!
      IF ( py_avg_using_b ) THEN
         WEIght_integral_by_b = .TRUE.
      ELSE
         WEIght_integral_by_b = .FALSE.
      ENDIF
!
! ****** Helium Fraction
!
      HE_frac = py_he_frac
      HE_rho = (1._R_TYP+4._R_TYP*HE_frac)/(1._R_TYP+2._R_TYP*HE_frac)
!
! ****** Value for limb-darkening coefficient.
!
      MU_ld = py_mu
!
! ****** Output pB file name.
!
      IF ( py_pb.EQ.'<none>' ) THEN
         COMpute_pb = .FALSE.
         PB_file = ' '
      ELSE
         COMpute_pb = .TRUE.
         PB_file = TRIM(py_pb)
      ENDIF
!
! ****** Output B file name.
!
      IF ( py_b.EQ.'<none>' ) THEN
         COMpute_b = .FALSE.
         B_File = ' '
      ELSE
         COMpute_b = .TRUE.
         B_File = TRIM(py_p)
      ENDIF

!
      END
!*==PRINT_HELP_TEXT.spg  processed by SPAG 6.72Dc at 22:10 on 27 May 2023
!#######################################################################
      SUBROUTINE PRINT_HELP_TEXT
!
!-----------------------------------------------------------------------
!
! ****** Print the help / syntax information.
!
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
!*--PRINT_HELP_TEXT3604
!
!-----------------------------------------------------------------------
!
      WRITE (*,*)
      WRITE (*,*)
      WRITE (*,*) '### Help for GETPB Special Integration Methods:'
      WRITE (*,*)
      WRITE (*,*) 'GETPB can now compute emissivity weighted'//         &
                 &' averages of other quantities.'
      WRITE (*,*) 'Here are the additions flags/options: '
      WRITE (*,*)
      WRITE (*,*)
      WRITE (*,*) 'SCALAR FIELD AVERAGING:'
      WRITE (*,*)
      WRITE (*,*) ' Use -scalar <file> to specify a 2D or 3D file'//    &
                 &' with a scalar field.'
      WRITE (*,*)
      WRITE (*,*) ' GETPB will integrate the scalar along the LOS'//    &
                 &' to get the weighted average.'
      WRITE (*,*)
      WRITE (*,*) ' Specify the output file with -avg_scalar <file>'
      WRITE (*,*)
      WRITE (*,*)
      WRITE (*,*) 'VECTOR FIELD COMPONENT AVERAGING:'
      WRITE (*,*)
      WRITE (*,*) ' Use -vr <file>, -vt <file>, -vp <file> to'//        &
                 &' specify a 2D or 3D vector field.'
      WRITE (*,*)
      WRITE (*,*) ' GETPB will project the rtp vector field into'//     &
                 &' LOS image coordinates and'
      WRITE (*,*) ' compute their weighted averages.'
      WRITE (*,*) ' '
      WRITE (*,*) ' Specify the output files for each components as:'
      WRITE (*,*) '   X   (+ is right):        -avg_vx <file>'
      WRITE (*,*) '   Y   (+ is up):           -avg_vy <file>'
      WRITE (*,*) '   LOS (+ is towards obs):  -avg_vlos <file>'
      WRITE (*,*)
      WRITE (*,*)
      WRITE (*,*) 'LOS ANGLE:'
      WRITE (*,*)
      WRITE (*,*) ' Use -avg_los_angle <file> to output the'//          &
                 &' weighted average LOS angle.'
      WRITE (*,*)
      WRITE (*,*) ' Here the "LOS angle" is the angle between a'//      &
                 &' position on the LOS and'
      WRITE (*,*) ' the closest approach of the LOS (rmin).'
      WRITE (*,*)
      WRITE (*,*) ' A positive LOS angle is in front of rmin,'//        &
                 &' negative is behind.'
      WRITE (*,*)
      WRITE (*,*) ' For plane-parallel integration, rmin lies'//        &
                 &' in the plane of sky.'
      WRITE (*,*)
      WRITE (*,*)
      WRITE (*,*) 'EMISSIVITY WEIGHTING:'
      WRITE (*,*)
      WRITE (*,*) ' GETPB can weight the LOS averages'//                &
                 &' by b OR pb (not both at once).'
      WRITE (*,*)
      WRITE (*,*) ' The default is to weight by the pb emissivity.'
      WRITE (*,*)
      WRITE (*,*) ' Set the flag -avg_using_b to weight by b'
      WRITE (*,*)
!
      CALL EXIT(1)
!
      END



