c
c-----------------------------------------------------------------------
c
c ****** Calculate the polarization brightness and brightness
c ****** from the electron density (produced by the MAS code).
c
c-----------------------------------------------------------------------
c
c ****** Updates and bug fixes:
c
c        05/19/94, ZM, Version 1.00:
c
c         - Original version of program.
c
c        01/24/95, ZM, Version 1.01:
c
c         - Generalized to handle 3D fields.
c         - Added the capability to enhance the density perturbation.
c         - Added the capability to apply a "vignetting function"
c           to pB.
c
c        02/06/95, ZM, Version 1.02:
c
c         - Added the ability to rotate the output pB image
c           in the theta direction.
c
c        04/28/95, ZM, Version 1.03:
c
c         - Added the ability to multiply pB by a radial vignetting
c           function.
c
c        03/03/97, ZM, Version 1.04:
c
c         - Added the ability to specify the solar B0 and P angles.
c
c        06/15/2004, ZM, Version 1.05:
c
c         - Updated the program to be able to deal with
c           nonuniform meshes in phi.
c
c        05/20/2009, ZM, Version 1.06:
c
c         - Allowed file names to be longer.
c
c        07/06/2010, ZM, Version 1.07:
c
c         - Corrected a bug in routine RDHDF that caused the program
c           to crash if the file being read from is not present.
c           This is a well-know bug in the HDF library that has been
c           corrected in our modern tools libraries.
c
c        08/30/2010, ZM, Version 1.08:
c
c         - Updated the program to use ZM's tools libraries,
c           and to use FORTRAN90.  This was a much-awaited
c           upgrade.  The new program uses a simpler line-of-sight
c           integration scheme, modeled after the one used in
c           GETEIT, with an adaptive integration step size
c           based on the local density mesh.
c         - The computation of the brightness B was added as
c           an option.
c         - Note that the usage line for this version is not
c           compatible with previous versions.
c         - Clarified the normalization of the output pB and B.
c           The output pB and B are in units of I0, the central
c           disk brightness, 2.49e10 [erg/cm^2/s/sr].
c
c        10/25/2018, ZM, Version 1.09:
c
c         - Generalized the program to allow the use of non-parallel
c           lines of sight to compute pB and B.  This allows
c           computations in regions that are large compared to the
c           distance from the Sun to the observer (such as for the
c           STEREO HI imagers and for WISPR on Parker Solar Probe).
c           The previous version assumed that the lines of sight
c           were parallel (i.e., corresponding to an observer at
c           infinity).  This new version has the capability to work
c           as before (i.e., an infinitely distant observer), which
c           is the default behavior.  To use the new feature (i.e.,
c           to calculate along non-parallel lines of sight), specify
c           the distance of the observer from the Sun in AU
c           [astronomical units] using the option -r <dist>.  The
c           previous results can be obtained by omitting this option.
c           When using this new option, for domains that are small
c           compared to the Sun-observer distance, the results
c           ought to be close to those calculated previously.  When
c           using the new option, the image domain [X0,X1] x [Y0,Y1]
c           is specified using elongation and altitude angles
c           [in degrees].
c         - Added shortcuts for the image view for the WISPR cameras
c           on PSP.  Use the flags -wispr1 and -wispr2 to specify
c           the image coordinates for these cameras.
c         - Removed the ability to specify the limits of the
c           integration along the line of sight (-z0 and -z1).
c           These had dubious merit.  The integration is now
c           always performed throughout the whole simulation
c           domain.
c         - This version has OpenMP commands to parallelize the
c           main loop of the computation.  These were adapted from
c           RL's version 1.08b from 04/12/2015.
c
c        08/20/2019, CD, Version 1.10:
c
c         - Updated the integrators to start and end only within the
c           bounds of the density mesh. This improves accuracy.
c         - Added the option for cubic spline interpolation, activated
c           with the -cubic flag.
c         - Cubic interpolation is used for rho and vignette function.
c         - These updates help if one applies special processing
c           to the images, which might otherwise bring out grid and/or
c           interpolation artifacts buried in the image.
c         - getpb now requires the zm_spline library as a dependency.
c
c        02/27/2020, CD, Version 1.11:
c
c         - Bugfix: Integrator update in v1.10 didn't account for LOSs
c           that pass outside the density mesh. Should be fine now.
c
c        11/16/2020, CD, Version 1.12:
c
c         - Added -he_frac flag to get ne from MAS rho if he_frac != 0.
c
c        04/01/2021, CD, Version 1.13:
c
c         - MAJOR UPDATE: Add option to compute emissivity weighted
c           integrals of scalar and vector fields.
c           Set the -help flag for a description of the new options.
c         - Add -mu flag to specify a custom limb-darkening constant.
c         - Rework how files are read and written to be more generic.
c         - SIDE EFFECT: If -scalar or -vr, vt, vp files are read in,
c           then the integration limits are adjusted to be the minimum
c           radius of every file read in (including rho). This can
c           change the integration limits from the half to main mesh,
c           leading to small floating point differences in pb & b
c           compared to previous versions (but that is OK!).
c
c-----------------------------------------------------------------------
c
c#######################################################################
c
c ****** Note that this tool uses modules from ZM's tools library.
c
c#######################################################################
      module ident
c
      character(*), parameter :: cname='GETPB'
      character(*), parameter :: cvers='1.13'
      character(*), parameter :: cdate='04/22/2021'
c
      end module
c#######################################################################
      module types
c
c-----------------------------------------------------------------------
c ****** Definition of data structures.
c-----------------------------------------------------------------------
c
      use number_types
      use sds_def
      use invint_def
      use spline_def
c
      implicit none
c
c ****** Inverse interpolation table structure definitions.
c
      type :: vtab
        type(itab), dimension(3) :: c
      end type
c
c ****** Container for 2D/3D MAS style fields and interpolation tables.
c
      type :: mas_field
        type(sds) :: sds
        type(vtab) :: invtab
        type(spl2d) :: spl2
        type(spl3d) :: spl3
      end type
c
      end module
c#######################################################################
      module constants
c
      use number_types
c
      implicit none
c
      real(r_typ), parameter :: pi=3.14159265358979323846_r_typ
      real(r_typ), parameter :: twopi=pi*2._r_typ
      real(r_typ), parameter :: halfpi=pi*.5_r_typ
      real(r_typ), parameter :: degtorad=pi/180._r_typ
c
c ****** One AU in solar radii.
c
      real(r_typ), parameter :: au_in_rs=214.939469_r_typ
c
c ****** Offset to ensure integration is inside domain.
c
      real(r_typ), parameter :: integration_offset=1.e-10_r_typ
c
      end module
c#######################################################################
      module params
c
c-----------------------------------------------------------------------
c ****** Parameters.
c-----------------------------------------------------------------------
c
      use number_types
c
      implicit none
c
      character(512) :: rho_file
c
      character(512) :: vr_file, vt_file, vp_file
      character(512) :: scalar_file
c
      character(512) :: vlos_file_out, vx_file_out, vy_file_out
      character(512) :: scalar_file_out
      character(512) :: angle_file_out
c
      character(512) :: pb_file
      character(512) :: b_file
c
      character(512) :: vf_file
c
      logical :: verbose
c
      logical :: compute_pb
      logical :: compute_b
c
c ****** Options for emissivity weighted integrals.
c
      logical :: do_integral_weighting=.false.
      logical :: compute_vector_integration=.false.
      logical :: compute_scalar_integration=.false.
      logical :: compute_angle_integration=.false.
c
c ****** Weight the scalar/vector integrals by B instead (default pB).
c
      logical :: weight_integral_by_b=.false.
c
c ****** Flag indicating that the density file is from the
c ****** old MAS code.
c
      logical :: oldmas
c
      real(r_typ) :: crlong
      real(r_typ) :: b0angle
      real(r_typ) :: pangle
      real(r_typ) :: rocc
      real(r_typ) :: power
      real(r_typ) :: disk_value
c
      logical :: power_set
c
c ****** Distance of the observer from the Sun.
c
      logical :: r_obs_set
      real(r_typ) :: r_obs
      real(r_typ) :: r_obs_rs
c
c ****** Step size multiplier.
c
      real(r_typ) :: dsmult
c
c ****** Shorthand flags to set the WISPR camera domains.
c
      logical :: wispr1
      logical :: wispr2
c
c ****** Helium Fraction Parameters (like in MAS).
c ****** HE_RHO is defined by rho(norm)=HE_RHO*n_e(norm),
c
      real(r_typ) :: he_frac = 0.0
      real(r_typ) :: he_rho
c
c ****** Option to use cubic interpolation
c
      logical :: cubic
c
c ****** Limb darkening coefficient at 520 nm. Default is from
c ****** Altschuler & Perry, Solar Physics, 23, 410 (1972).
c
      real(r_typ) :: mu_ld = 0.63_r_typ
c
c ****** Flag for the precision of the output files.
c
      logical :: hdf32
c
      end module
c#######################################################################
      module image_region
c
      use number_types
c
      implicit none
c
c ****** Image dimensions.
c
      integer :: nx
      integer :: ny
c
c ****** Image region limits.
c
      real(r_typ) :: x0
      real(r_typ) :: x1
      real(r_typ) :: y0
      real(r_typ) :: y1
c
      logical :: image_limits_set
c
c ****** Number of OpenMP iterations to do in each thread.
c
      integer, parameter :: iterations_per_thread=500
c
      end module
c#######################################################################
      module geometry
c
      use number_types
c
      implicit none
c
c ****** Sin and cos of the P and B0 angles.
c
      real(r_typ) :: sin_p,cos_p
      real(r_typ) :: sin_b0,cos_b0
c
c ****** Central meridian longitude.
c
      real(r_typ) :: cml
c
c ****** Flag to indicate an axisymmetric density (i.e., 2D).
c
      logical :: twodee
c
c ****** Maximum radius for all of the input meshes.
c
      real(r_typ) :: rmax=-1.0
c
      end module
c#######################################################################
      module mas_fields
c
c-----------------------------------------------------------------------
c ****** Storage for the density.
c-----------------------------------------------------------------------
c
      use number_types
      use types
      use spline_def
c
      implicit none
c
c ****** Density field container.
c
      type(mas_field) :: rho
c
c ****** Containers for arbitrary vector field, v, in spherical coords.
c
      type(mas_field) :: vr
      type(mas_field) :: vt
      type(mas_field) :: vp
c
c ****** Containers for arbitrary scalar field.
c
      type(mas_field) :: sfield
c
c ****** Normalization factor for electron density.
c
c ****** Multiply the normalized density read in [in MAS code units]
c ****** by FN_N to get electron density in [/cm^3].
c
      real(r_typ), parameter :: fn_n=1.e8_r_typ
c
      end module
c#######################################################################
      module image_fields
c
c-----------------------------------------------------------------------
c ****** Storage for the computed pB and B, and their associated
c ****** arrays.
c-----------------------------------------------------------------------
c
      use number_types
      use types
c
      implicit none
c
c ****** Image axes for the case when the observer distance
c ****** has been specified (elongation, altitude).
c
      real(r_typ), dimension(:), pointer :: elong,alt
c
c ****** Image axes for the case when the observer distance
c ****** is infinite (POS x and y).
c
      real(r_typ), dimension(:), pointer :: x_pos,y_pos
c
c ****** Polarization brightness and brightness.
c
      real(r_typ), dimension(:,:), allocatable :: pb
      real(r_typ), dimension(:,:), allocatable :: b
c
c ****** Emissivity weighted averages of the vector field components
c ****** in the frame of the image plane (los is TOWARDS observer).
c
      real(r_typ), dimension(:,:), allocatable :: vlos_avg
      real(r_typ), dimension(:,:), allocatable :: vx_avg
      real(r_typ), dimension(:,:), allocatable :: vy_avg
c
c ****** Emissivity weighted averages of the scalar field.
c
      real(r_typ), dimension(:,:), allocatable :: sf_avg
c
c ****** Emissivity weighted average elongation angle (0 when r=rmin).
c
      real(r_typ), dimension(:,:), allocatable :: angle_avg
c
c ****** Closest distance to the Sun.  This is used for vignetting
c ****** and multiplication by a radial power.
c
      real(r_typ), dimension(:,:), allocatable :: r_min
c
      end module
c#######################################################################
      module vignetting_function
c
c-----------------------------------------------------------------------
c ****** Vignetting function definition.
c-----------------------------------------------------------------------
c
      use number_types
      use spline_def
c
      implicit none
c
c ****** Vignetting function.
c
      integer :: nvf
      real(r_typ), dimension(:), allocatable :: r_vf
      real(r_typ), dimension(:), allocatable :: vf
c
c ****** Spline structure for cubic interpolation.
c
      type(spl1d) :: vf_spl1
c
      end module

c#######################################################################


c#######################################################################
      program GETPB
c
c-----------------------------------------------------------------------
c
      use ident
      use number_types
      use constants
      use params
      use image_region
      use geometry
      use mas_fields
      use image_fields
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      real(r_typ), dimension(:), pointer :: x,y
      REAL, intent(out), dimension(:, :) :: ans
c
c-----------------------------------------------------------------------
c
c ****** Set the parameters.
c
      call python_getpb(py_rho, py_b, py_pb, py_help, py_verbose, py_cubic, py_oldmas, py_long, py_p, py_b0, py_r, py_nx, py_ny, py_x0, py_x1, py_y0, py_y1, py_wispr1, py_wispr2, py_rocc, py_dsmult, py_power, py_disk, py_vf, py_vr, py_vt, py_vp, py_scalar,
     &   py_avg_scalar, py_avg_los_angle, py_avg_vlos, py_avg_vx, py_avg_vy, py_avg_using_b, py_he_frac, py_mu)
c
c ****** Check that at least one of B or pB was requested.
c
      if (.not.(rho_file.eq.'<none>')) then
        write (*,*)
        write (*,*) '### ERROR in ',cname,':'
        write (*,*) '### Rho is required.'
        call exit (1)
      end if

      if (.not.(compute_b.or.compute_pb)) then
        write (*,*)
        write (*,*) '### ERROR in ',cname,':'
        write (*,*) '### At least one of pB or B must be requested.'
        call exit (1)
      end if
c
c ****** Check that the image dimensions are valid.
c
      if (nx.lt.1.or.ny.lt.1) then
        write (*,*)
        write (*,*) '### ERROR in ',cname,':'
        write (*,*) '### The specified number of points for the'//
     &              ' image are invalid:'
        write (*,*) 'NX = ',nx
        write (*,*) 'NY = ',ny
        call exit (1)
      end if
c
c ****** Check that the WISPR camera shorthand flags are
c ****** consistent with the usage of other options.
c
      if (wispr1.and.wispr2) then
        write (*,*)
        write (*,*) '### ERROR in ',cname,':'
        write (*,*) '### You cannot specify both the'//
     &              ' -wispr1 and -wispr2 flags.'
        call exit (1)
      end if
c
      if ((wispr1.or.wispr2).and..not.r_obs_set) then
        write (*,*)
        write (*,*) '### ERROR in ',cname,':'
        write (*,*) '### You must specify the Sun-observer'//
     &              ' distance when using the -wispr1'
        write (*,*) '### or -wispr2 flags.'
        call exit (1)
      end if
c
      if ((wispr1.or.wispr2).and.image_limits_set) then
        write (*,*)
        write (*,*) '### ERROR in ',cname,':'
        write (*,*) '### You cannot set the image coordinate'//
     &              ' limits when using the -wispr1 or'
        write (*,*) '### -wispr2 flags.'
        call exit (1)
      end if
c
c ****** Check that the observer distance is positive (if
c ****** it was specified).
c
      if (r_obs_set) then
        if (r_obs.le.0.) then
          write (*,*)
          write (*,*) '### ERROR in ',cname,':'
          write (*,*) '### The distance of the observer'//
     &                ' from the Sun must be positive:'
          write (*,*) 'Specified distance [AU] = ',r_obs
          call exit (1)
        end if
      end if
c
c ****** Set the image region for WISPR, if requested.
c
      if (wispr1) then
        x0=13.5_r_typ
        x1=51.4_r_typ
        y0=-20._r_typ
        y1=20._r_typ
      else if (wispr2) then
        x0=49.7_r_typ
        x1=104.7_r_typ
        y0=-29._r_typ
        y1=29._r_typ
      end if
c
c ****** Check the specified image plane coordinates.
c
      if (x0.gt.x1.or.y0.gt.y1) then
        write (*,*)
        write (*,*) '### ERROR in ',cname,':'
        write (*,*) '### The specified image plane coordinates'//
     &              ' are invalid:'
        write (*,*) 'X0 = ',x0
        write (*,*) 'X1 = ',x1
        write (*,*) 'Y0 = ',y0
        write (*,*) 'Y1 = ',y1
        call exit (1)
      end if
c
c ****** Read the density.
c
      call read_field( rho_file, 'density', rho)
c
c ****** Read the scalar field (if necessary).
c
      if (compute_scalar_integration) then
        call read_field( scalar_file, 'scalar field', sfield)
      endif
c
c ****** Read the vector fields (if necessary).
c
      if (compute_vector_integration) then
        if (vr_file.eq.' '.or.vt_file.eq.' '.or.vp_file.eq.' ') then
          write (*,*)
          write (*,*) '### ERROR in ',cname,':'
          write (*,*) '### Not all vector files were specified!'
          write (*,*) '  -vr: ', trim(vr_file)
          write (*,*) '  -vt: ', trim(vt_file)
          write (*,*) '  -vp: ', trim(vp_file)
          call exit (1)
        endif
        call read_field( vr_file, 'vector r component', vr)
        call read_field( vt_file, 'vector t component', vt)
        call read_field( vp_file, 'vector p component', vp)
      endif
c
c ****** Read the vignetting function if it was requested.
c
      if (vf_file.ne.' ') then
        call read_vignetting_function (vf_file)
      end if
c
c ****** Convert the observer distance to solar radii.
c
      r_obs_rs=r_obs*au_in_rs
c
      if (verbose) then
        write (*,*)
        write (*,*) 'Longitude [deg] = ',crlong
        write (*,*) 'B0 angle [deg]  = ',b0angle
        write (*,*) 'P angle [deg]   = ',pangle
        write (*,*) 'Occulting disk radius [R_sun] = ',rocc
        if (r_obs_set) then
          write (*,*) 'Observer distance [AU]    = ',r_obs
          write (*,*) 'Observer distance [R_sun] = ',r_obs_rs
        else
          write (*,*) 'Observer distance [AU] = infinite'
        end if
        if (cubic) then
          write (*,*) 'Using cubic spline interpolation.'
        end if
      end if
c
c ****** Convert the Carrington longitude to radians.
c
      cml=crlong*degtorad
c
c ****** Convert the solar P and B0 angles to radians.
c
      pangle=pangle*degtorad
      b0angle=b0angle*degtorad
c
      sin_p=sin(pangle)
      cos_p=cos(pangle)
      sin_b0=sin(b0angle)
      cos_b0=cos(b0angle)
c
      if (verbose) then
        write (*,*)
        write (*,*) '### Image region:'
        if (r_obs_set) then
          write (*,901) 'Elongation [deg]:',
     &                  ' X0 = ',x0,
     &                  ' X1 = ',x1
          write (*,901) 'Altitude   [deg]:',
     &                  ' Y0 = ',y0,
     &                  ' Y1 = ',y1
        else
          write (*,901) 'POS x limits [R_sun]:',
     &                  ' X0 = ',x0,
     &                  ' X1 = ',x1
          write (*,901) 'POS y limits [R_sun]:',
     &                  ' Y0 = ',y0,
     &                  ' Y1 = ',y1
        end if
      end if
  901 format (1x,a,a,f15.6,',',a,f15.6)
c
c ****** Allocate the arrays for the results, initialized with zeros.
c
      allocate (r_min(nx,ny))
      if (compute_pb) then
        allocate (pb(nx,ny))
        pb(:,:)=0.
      endif
      if (compute_b) then
        allocate (b(nx,ny))
        b(:,:)=0.
      endif
      if (compute_scalar_integration) then
        allocate(sf_avg(nx,ny))
        sf_avg(:,:)=0.
      endif
      if (compute_angle_integration) then
        allocate(angle_avg(nx,ny))
        angle_avg(:,:)=0.
      endif
      if (compute_vector_integration) then
        allocate(vlos_avg(nx,ny))
        allocate(vx_avg(nx,ny))
        allocate(vy_avg(nx,ny))
        vlos_avg(:,:)=0.
        vx_avg(:,:)=0.
        vy_avg(:,:)=0.
      endif
c
c ****** Print the helium fraction if it was specified.
c
      if (verbose.and.he_frac.ne.0.) then
        write (*,*)
        write (*,*) 'Using a non-zero helium fraction!'
        write (*,'(A,F10.7)') '  he_frac: ', he_frac
        write (*,'(A,F10.7)') '  he_rho:  ', he_rho
      end if
c
c ****** Print the limb-darkening coefficient.
c
      if (verbose.and.mu_ld.ne.0.63_r_typ) then
        write (*,*)
        write (*,*) 'Using a non-standard limb-darkening coefficient!'
        write (*,'(A,F10.7)') ' mu: ', mu_ld
      endif
c
c ****** Error check the integral weighting option.
c
      if (compute_vector_integration.or.compute_scalar_integration.or.
     &    compute_angle_integration) then
c
        do_integral_weighting=.true.
c
        if (weight_integral_by_b) then
          if (.not.compute_b) then
            write (*,*)
            write (*,*) '### ERROR in ',cname,':'
            write (*,*) '### Integral weighting specified with B but'//
     &                  ' -b was not specified!'
            call exit (1)
          endif
          if (verbose) then
            write (*,*)
            write (*,*) 'Integral averages requested! Weighting: B.'
          endif
        endif
        if (.not.weight_integral_by_b) then
          if(.not.compute_pb) then
            write (*,*)
            write (*,*) '### ERROR in ',cname,':'
            write (*,*) '### Integral weighting specified with pB but'//
     &                  ' -pb was not specified!'
            call exit (1)
          endif
          if (verbose) then
            write (*,*)
            write (*,*) 'Integral averages requested! Weighting: pB.'
          endif
        endif
      endif
c
c ****** Compute pB and/or B.
c
      if (verbose) then
        write (*,*)
        write (*,*) 'Calculating pB and/or B ...'
      end if
c
      if (r_obs_set) then
        call get_image_obs_distance_specified
        x=>elong
        y=>alt
      else
        call get_image_obs_at_infinity
        x=>x_pos
        y=>y_pos
      end if
c
c ****** Normalize wieghted integrals if necessary.
c
      if (do_integral_weighting) then
        call normalize_integrals
      endif
c
c ****** Multiply pB and B by the requested radial power.
c
      if (power_set) then
        if (compute_pb) then
          if (verbose) then
            write (*,*)
            write (*,*) '### Multiplying pB by r_min^p with p = ',
     &                  power
          end if
          call multiply_by_radial_power (pb)
        end if
        if (compute_b) then
          if (verbose) then
            write (*,*)
            write (*,*) '### Multiplying B by r_min^p with p = ',
     &                  power
          end if
          call multiply_by_radial_power (b)
        end if
      end if
c
c ****** Apply the vignetting function (if it was specified).
c
      if (vf_file.ne.' ') then
        if (compute_pb) then
          if (verbose) then
            write (*,*)
            write (*,*) '### Multiplying pB by the vignetting'//
     &                  ' function ... '
          end if
          call vignette (pb)
        end if
        if (compute_b) then
          if (verbose) then
            write (*,*)
            write (*,*) '### Multiplying B by the vignetting'//
     &                  ' function ... '
          end if
          call vignette (b)
        end if
      end if
c
c ****** Set the precision of the output files.
c
      hdf32 = rho%sds%hdf32
c
c ****** Write the pB image.
c
      if (compute_pb) ans = pb
!call write_image(pb_file, 'pB', x, y, pb)
c ****** We want it to return pb/b for the python version

c
c ****** Write the B image.
c
c      if (compute_b) call write_image(b_file, 'B', x, y, b)
       if (compute_b) ans = b
c
c ****** Write the scalar field image.
c
      if (compute_scalar_integration)
     &    call write_image(scalar_file_out, 'LOS averaged scalar field',
     &                     x, y, sf_avg)
c
c ****** Write the LOS contribution angle image.
c
      if (compute_angle_integration)
     &    call write_image(angle_file_out,
     &           'LOS averaged contribution angle', x, y, angle_avg)
c
c ****** Write the vector field images.
c
      if (compute_vector_integration) then
        call write_image(vlos_file_out, 'LOS averaged vlos',
     &                   x,y,vlos_avg)
        call write_image(vx_file_out, 'LOS averaged vx',x,y,vx_avg)
        call write_image(vy_file_out, 'LOS averaged vy',x,y,vy_avg)
      endif
c
      call exit (0)

      end
c#######################################################################
      subroutine get_image_obs_distance_specified
c
c-----------------------------------------------------------------------
c
c ****** Compute pB and/or B for the case when the observer distance
c ****** has been specified.
c
c ****** This uses non-parallel lines of sight for the computation,
c ****** and should be used when the image domain is of significant
c ****** size compared to the Sun-observer distance.
c
c ****** In this case, the image coordinates [X0,X1] x [Y0,Y1] are
c ****** specified in terms of elongation and altitude angles with
c ****** respect to the ecliptic plane.
c
c-----------------------------------------------------------------------
c
      use ident
      use number_types
      use constants
      use params
      use image_region
      use geometry
      use image_fields
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      real(r_typ), parameter :: half=.5_r_typ
      real(r_typ), parameter :: one=1._r_typ
c
c-----------------------------------------------------------------------
c
      integer :: i,j,nstop
      real(r_typ) :: dx,dy
      real(r_typ) :: s0,s1
      real(r_typ) :: dmax
      real(r_typ) :: e_rad,a_rad,cos_a
      real(r_typ) :: s,t,ds,dsm
      real(r_typ) :: t_r_min
      real(r_typ) :: ds_local
      real(r_typ) :: rv,tv,pv
      real(r_typ) :: pb_local,b_local
      real(r_typ), dimension(3) :: v_obs,v_ref,v_los,v_r_min
c
c-----------------------------------------------------------------------
c
      real(r_typ) :: xo,yo,zo,xr,yr,zr,xv,yv,zv
      real(r_typ) :: ds_fac, x_los
c
c ****** Unit vectors of observer frame in heliographic coordinates.
c
      real(r_typ), dimension(3) :: ex, ey, ez
c
c-----------------------------------------------------------------------
c
c ****** Note that the (x,y,z) coordinates are defined such
c ****** that the x-axis points from the Sun to the observer,
c ****** and the z-axis points vertically upward.
c
c ****** Construct the coordinates for the axes of the image.
c ****** These are specified in terms of [elongation,altitude]
c ****** angles, in degrees.
c
c ****** The elongation is the angle between the Sun-observer line
c ****** and the LOS, measured positive towards the West limb
c ****** (with the Sun being at elongation=0).
c
c ****** The altitude is the angle expressing the vertical
c ****** inclination of the LOS with respect to the ecliptic
c ****** plane (with positive angles specifying a LOS above the
c ****** ecliptic plane, and negative angles specifying a LOS
c ****** below the ecliptic plane.
c
      allocate (elong(nx))
      allocate (alt(ny))
c
      if (nx.eq.1) then
        elong(1)=x0
      else
        dx=(x1-x0)/(nx-1)
        do i=1,nx
          elong(i)=x0+(i-1)*dx
        enddo
      end if
c
      if (ny.eq.1) then
        alt(1)=y0
      else
        dy=(y1-y0)/(ny-1)
        do j=1,ny
          alt(j)=y0+(j-1)*dy
        enddo
      end if
c
c ****** Cartesian vector from the Sun to the observer.
c
      v_obs(1)=r_obs_rs
      v_obs(2)=0.
      v_obs(3)=0.
c
c ****** Calculate pB and/or B.
c
c$omp parallel do shared(pb,b,r_min)
c$omp& shared(sf_avg,angle_avg,vlos_avg,vx_avg,vy_avg)
c$omp& shared(v_obs)
c$omp& private(ex,ey,ez)
c$omp& private(i,j,ds_fac,x_los,xo,yo,zo,xr,yr,zr,xv,yv,zv)
c$omp& private(e_rad,a_rad,cos_a,v_ref,t_r_min,v_r_min)
c$omp& private(t,s,dsm,ds,ds_local,nstop,v_los,rv,tv,pv,s0,s1,dmax)
c$omp& private(pb_local,b_local)
c$omp& collapse(2)
c$omp& schedule(dynamic,iterations_per_thread)
      do j=1,ny
        do i=1,nx
c
c ****** Cartesian vector from the Sun to the reference point.
c
c ****** This is a point on a sphere of radius R_OBS centered
c ****** at the observer, where R_OBS is the distance of the
c ****** observer from the Sun.
c
          e_rad=degtorad*elong(i)
          a_rad=degtorad*alt(j)
          cos_a=cos(a_rad)
          v_ref(1)=r_obs_rs*(one+cos(pi-e_rad)*cos_a)
          v_ref(2)=r_obs_rs*sin(pi-e_rad)*cos_a
          v_ref(3)=r_obs_rs*sin(a_rad)
c
c ****** Find the point along this LOS that is closest to the
c ****** Sun.  This is the point on the LOS at which it
c ****** intersects the Thomson sphere.  The distance from this
c ****** point to the Sun is needed in the calculation of pB
c ****** and B.
c
          t_r_min=-( v_obs(1)*(v_ref(1)-v_obs(1))
     &              +v_obs(2)*(v_ref(2)-v_obs(2))
     &              +v_obs(3)*(v_ref(3)-v_obs(3))
     &             )/r_obs_rs**2
c
          v_r_min(1)=t_r_min*v_ref(1)+(one-t_r_min)*v_obs(1)
          v_r_min(2)=t_r_min*v_ref(2)+(one-t_r_min)*v_obs(2)
          v_r_min(3)=t_r_min*v_ref(3)+(one-t_r_min)*v_obs(3)
c
          r_min(i,j)=sqrt( v_r_min(1)**2
     &                    +v_r_min(2)**2
     &                    +v_r_min(3)**2)
c
c ****** Check if this LOS intersects the occulting disk.
c ****** If so, skip the calculation along this LOS.
c
          if (r_min(i,j).le.max(rocc,one)) then
            if (compute_pb) pb(i,j)=disk_value
            if (compute_b) b(i,j)=disk_value
            cycle
          end if
c
c ****** The distance s is along the LOS, starting at the observer
c ****** (s = 0), and increasing toward the reference point
c ****** (s = R_OBS).  The parameter t is proportional to s,
c ****** with t = 0 specifying the observer, and t = 1 specifying
c ****** the reference point.
c
          if (compute_pb) pb(i,j)=0.
          if (compute_b) b(i,j)=0.
c
c ****** Check if the LOS doesn't intersect the grid, leave it zeroed.
c
          if (r_min(i,j).ge.rmax) cycle
c
c ****** Compute unit vectors for LOS if doing vector/angle projection.
c ****** Here the projection changes for each pixel (non-parallel).
c
         if (compute_vector_integration.or.
     &       compute_angle_integration) then
           call transform_position (v_obs(1),v_obs(2),v_obs(3),rv,tv,pv)
           call rtp2xyz(rv,tv,pv,xo,yo,zo)
           call transform_position (v_r_min(1),v_r_min(2),v_r_min(3),
     &                              rv,tv,pv)
           call rtp2xyz(rv,tv,pv,xr,yr,zr)
           xv=xo-xr
           yv=yo-yr
           zv=zo-zr
           call get_unit_vectors(xv,yv,zv,ex,ey,ez)
         endif
c
c ****** Integrate along the LOS, starting at first intersection
c ****** of the LOS with the density mesh.
c
          dmax=sqrt( rmax**2 - r_min(i,j)**2)
          s0 = t_r_min*r_obs_rs - dmax + integration_offset
          s1 = t_r_min*r_obs_rs + dmax - integration_offset
c
          s=s0
          dsm=0.
          nstop=0
          do
c
c ****** Set the position along the LOS.
c
            t=s/r_obs_rs
c
            v_los(1)=t*v_ref(1)+(one-t)*v_obs(1)
            v_los(2)=t*v_ref(2)+(one-t)*v_obs(2)
            v_los(3)=t*v_ref(3)+(one-t)*v_obs(3)
c
            call transform_position (v_los(1),v_los(2),v_los(3),
     &                               rv,tv,pv)
c
            call get_kernels (r_min(i,j),rv,tv,pv,
     &                        compute_pb,compute_b,
     &                        pb_local,b_local,ds_local)
c
            ds=ds_local*dsmult
            s=s+ds
            if (s.ge.s1) then
              ds=ds-(s-s1)
              s=s1
              nstop=nstop+1
            end if
            ds_fac = half*(dsm+ds)
c
            if (compute_pb) then
              pb(i,j)=pb(i,j) + ds_fac*pb_local
            end if
c
            if (compute_b) then
              b(i,j)=b(i,j) + ds_fac*b_local
            end if
c
c ****** Compute emissivity weighted integrals if needed.
c ****** For non-plane-parallel angle need to dot position with x-hat.
c
            if (do_integral_weighting) then
              if (compute_angle_integration) then
                call rtp2xyz(rv,tv,pv,xv,yv,zv)
                x_los = ex(1)*xv + ex(2)*yv + ex(3)*zv
              else
                x_los = 0.
              endif
c
              call add_to_integrals(i,j,pb_local,b_local,ds_fac,
     &               ex,ey,ez,rv,tv,pv,x_los)
            endif
c
            dsm=ds
            if (nstop.ge.2) exit
          enddo
        enddo
      enddo
c$omp end parallel do
c
      end
c#######################################################################
      subroutine get_image_obs_at_infinity
c
c-----------------------------------------------------------------------
c
c ****** Compute pB and/or B for the case when the observer is at
c ****** infinity.
c
c ****** This uses parallel lines of sight for the computation.
c
c ****** In this case, the image coordinates [X0,X1] x [Y0,Y1] are
c ****** specified in terms of horizontal and vertical coordinates
c ****** in the plane of the sky.
c
c-----------------------------------------------------------------------
c
      use ident
      use number_types
      use constants
      use params
      use image_region
      use geometry
      use image_fields
      use mas_fields
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      real(r_typ), parameter :: zero=0._r_typ
      real(r_typ), parameter :: half=.5_r_typ
      real(r_typ), parameter :: one=1._r_typ
c
c-----------------------------------------------------------------------
c
      integer :: i,j,nstop
      real(r_typ) :: dx,dy
      real(r_typ) :: s0,s1
      real(r_typ) :: dmax
      real(r_typ) :: s,ds,dsm
      real(r_typ) :: ds_local
      real(r_typ) :: rv,tv,pv
      real(r_typ) :: pb_local,b_local
      real(r_typ), dimension(3) :: v_los
c
c-----------------------------------------------------------------------
c
      real(r_typ) :: xv,yv,zv
      real(r_typ) :: ds_fac
c
c ****** Unit vectors of observer frame in heliographic coordinates.
c
      real(r_typ), dimension(3) :: ex, ey, ez
c
c-----------------------------------------------------------------------
c
c ****** Note that the (x,y,z) coordinates are defined such
c ****** that the x-axis points from the Sun to the observer,
c ****** and the z-axis points vertically upward.
c
c ****** Construct the coordinates for the axes of the image.
c ****** These are specified in terms of horizontal and vertical
c ****** plane-of-sky positions, in solar radii.
c
      allocate (x_pos(nx))
      allocate (y_pos(ny))
c
      if (nx.eq.1) then
        x_pos(1)=x0
      else
        dx=(x1-x0)/(nx-1)
        do i=1,nx
          x_pos(i)=x0+(i-1)*dx
        enddo
      end if
c
      if (ny.eq.1) then
        y_pos(1)=y0
      else
        dy=(y1-y0)/(ny-1)
        do j=1,ny
          y_pos(j)=y0+(j-1)*dy
        enddo
      end if
c
c ****** Compute unit vectors for LOS if doing vector projection.
c ****** For the plane-parallel assumption, they never change.
c
      if (compute_vector_integration) then
        call transform_position (one,zero,zero,rv,tv,pv)
        call rtp2xyz(rv,tv,pv,xv,yv,zv)
        call get_unit_vectors(xv,yv,zv,ex,ey,ez)
      endif
c
c ****** Calculate pB and/or B.
c
c$omp parallel do shared(pb,b,r_min)
c$omp& shared(sf_avg,angle_avg,vlos_avg,vx_avg,vy_avg)
c$omp& shared(ex,ey,ez)
c$omp& private(i,j,ds_fac)
c$omp& private(s,dsm,ds,ds_local,nstop,v_los,rv,tv,pv,s0,s1,dmax)
c$omp& private(pb_local,b_local)
c$omp& collapse(2)
c$omp& schedule(dynamic,iterations_per_thread)
      do j=1,ny
        do i=1,nx
c
c ****** Calculate the plane-of-sky radius.
c
          r_min(i,j)=sqrt(x_pos(i)**2+y_pos(j)**2)
c
c ****** Check if this LOS intersects the occulting disk.
c ****** If so, skip the calculation along this LOS.
c
          if (r_min(i,j).le.max(rocc,one)) then
            if (compute_pb) pb(i,j)=disk_value
            if (compute_b) b(i,j)=disk_value
            cycle
          end if
c
c ****** The distance s is along the LOS, referenced with respect
c ****** to the solar limb (s = 0), and increasing away from the
c ****** observer (so s = -x == v_los(1) in the observer frame).
c
          if (compute_pb) pb(i,j)=0.
          if (compute_b) b(i,j)=0.
c
c ****** Check if the LOS doesn't intersect the grid, leave it zeroed.
c
          if (r_min(i,j).ge.rmax) cycle
c
c ****** Integrate along the LOS, starting at first intersection
c ****** of the LOS with the density mesh.
c
          dmax=sqrt( rmax**2 - r_min(i,j)**2)
          s0 = -dmax + integration_offset
          s1 =  dmax - integration_offset
c
          s=s0
          dsm=0.
          nstop=0
          do
c
c ****** Set the position along the LOS.
c
            v_los(1)=-s
            v_los(2)=x_pos(i)
            v_los(3)=y_pos(j)
c
            call transform_position (v_los(1),v_los(2),v_los(3),
     &                               rv,tv,pv)
c
            call get_kernels (r_min(i,j),rv,tv,pv,
     &                        compute_pb,compute_b,
     &                        pb_local,b_local,ds_local)
c
            ds=ds_local*dsmult
            s=s+ds
            if (s.ge.s1) then
              ds=ds-(s-s1)
              s=s1
              nstop=nstop+1
            end if
            ds_fac = half*(dsm+ds)
c
c ****** Compute pb and/or b contributions.
c
            if (compute_pb) then
              pb(i,j)=pb(i,j) + ds_fac*pb_local
            end if
            if (compute_b) then
              b(i,j)=b(i,j) + ds_fac*b_local
            end if
c
c ****** Compute emissivity weighted integrals if needed.
c
            if (do_integral_weighting) then
              call add_to_integrals(i,j,pb_local,b_local,ds_fac,
     &               ex,ey,ez,rv,tv,pv,v_los(1))
c
            endif
c
            dsm=ds
            if (nstop.ge.2) exit
          enddo
c
        enddo
      enddo
c$omp end parallel do
c
      end
c#######################################################################
      subroutine read_field( filename, fieldname, field)
c
c-----------------------------------------------------------------------
c
c ****** Read a MAS style field.
c
c-----------------------------------------------------------------------
c
      use ident
      use params
      use geometry
      use types
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      character(512) :: filename
      character(*) :: fieldname
      type(mas_field) :: field
c
c-----------------------------------------------------------------------
c
      integer :: ierr
c
c-----------------------------------------------------------------------
c
      if (verbose) then
        write (*,*)
        write (*,*) '### Reading ',trim(fieldname),
     &              ' from file: ',trim(filename)
      end if
c
c ****** Read the field from the specified file.
c
      call rdhdf (filename,field%sds,ierr)
c
      if (ierr.ne.0) then
        write (*,*)
        write (*,*) '### ERROR in ',cname,':'
        write (*,*) '### Could not read ',trim(fieldname),'.'
        write (*,*) 'IERR (from RDHDF) = ',ierr
        write (*,*) 'File name: ',trim(filename)
        call exit (1)
      end if
c
c ****** Check that it contains a 2D or 3D field.
c
      if (.not.(field%sds%ndim.eq.2.or.field%sds%ndim.eq.3)) then
        write (*,*)
        write (*,*) '### ERROR in ',cname,':'
        write (*,*) '### ',trim(fieldname),' file does'//
     &              ' not contain a 2D or 3D field.'
        write (*,*) 'Number of dimensions = ',field%sds%ndim
        write (*,*) 'File name: ',trim(filename)
        call exit (1)
      end if
c
c ****** Check that it contains scales.
c
      if (.not.field%sds%scale) then
        write (*,*)
        write (*,*) '### ERROR in ',cname,':'
        write (*,*) '### The density file does'//
     &              ' not contain scales.'
        write (*,*) 'File name: ',trim(filename)
        call exit (1)
      end if
c
c ****** Set the flag for an axisymmetric case.
c
      if (field%sds%ndim.eq.2) then
        twodee=.true.
      else
        twodee=.false.
      end if
c
c ****** Check that the field array has at least two points in
c ****** each (non-negligible) dimension.
c
      if (field%sds%dims(1).lt.2) then
        write (*,*)
        write (*,*) '### ERROR in ',cname,':'
        write (*,*) '### ',trim(fieldname),' array must have at least'//
     &              ' 2 points in the r dimension.'
        write (*,*) 'Number of points in r = ',field%sds%dims(1)
        write (*,*) 'File name: ',trim(filename)
        call exit (1)
      end if
c
      if (field%sds%dims(2).lt.2) then
        write (*,*)
        write (*,*) '### ERROR in ',cname,':'
        write (*,*) '### ',trim(fieldname),' array must have at least'//
     &              ' 2 points in the theta dimension.'
        write (*,*) 'Number of points in theta = ',field%sds%dims(2)
        write (*,*) 'File name: ',trim(filename)
        call exit (1)
      end if
c
      if (.not.twodee) then
        if (field%sds%dims(3).lt.2) then
          write (*,*)
          write (*,*) '### ERROR in ',cname,':'
          write (*,*) '### ',trim(fieldname),' array must have at'//
     &                ' least 2 points in the phi dimension.'
          write (*,*) 'Number of points in phi = ',field%sds%dims(3)
          write (*,*) 'File name: ',trim(filename)
          call exit (1)
        end if
      end if
c
c ****** If we are using old MAS files, add a phi point.
c
      if (oldmas) then
        call add_phi_point (field%sds)
      end if
c
c ****** Construct the inverse interpolation table for
c ****** faster interpolation.
c
      call build_inverse_tables (field%sds,field%invtab)
c
c ****** Set the maximum radius available in the file but check
c ****** that a smaller value hasn't been set yet (rmax is initialized
c ****** as a negative number).
c
      if (rmax.lt.0.) then
        rmax=field%sds%scales(1)%f(field%sds%dims(1))
      else
        rmax=min(rmax,field%sds%scales(1)%f(field%sds%dims(1)))
      endif
c
c ****** Compute the spline coefficients if needed.
c
      if (cubic) then
        if (verbose) then
          write (*,*)
          write (*,*) '### Computing cubic spline coefficients for ',
     &                 trim(fieldname)
        end if
        if (field%sds%ndim.eq.2) then
          call compute_spline_2d (field%sds%dims(1),
     &                            field%sds%dims(2),
     &                            field%sds%scales(1)%f,
     &                            field%sds%scales(2)%f,
     &                            field%sds%f,field%spl2)
        else if (field%sds%ndim.eq.3) then
          call compute_spline_3d (field%sds%dims(1),
     &                            field%sds%dims(2),
     &                            field%sds%dims(3),
     &                            field%sds%scales(1)%f,
     &                            field%sds%scales(2)%f,
     &                            field%sds%scales(3)%f,
     &                            field%sds%f,field%spl3)
        endif
      endif
c
      return
      end
c#######################################################################
      subroutine add_phi_point (s)
c
c-----------------------------------------------------------------------
c
c ****** Add a point in the phi dimension in a periodic manner
c ****** to the SDS in structure S.
c
c-----------------------------------------------------------------------
c
      use number_types
      use sds_def
      use constants
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      type(sds) :: s
c
c-----------------------------------------------------------------------
c
      real(r_typ), dimension(:,:,:), pointer :: f
      real(r_typ), dimension(:), pointer :: p
c
c-----------------------------------------------------------------------
c
c ****** For a 2D (axisymmetric) case, this is not needed.
c
      if (s%ndim.ne.3) return
c
c ****** Add a point to the phi dimension.
c
      allocate (f(s%dims(1),s%dims(2),s%dims(3)+1))
      allocate (p(s%dims(3)+1))
c
      f(:,:,1:s%dims(3))=s%f(:,:,:)
      f(:,:,s%dims(3)+1)=s%f(:,:,1)
c
      p(1:s%dims(3))=s%scales(3)%f(:)
      p(s%dims(3)+1)=s%scales(3)%f(1)+twopi
      s%dims(3)=s%dims(3)+1
c
      deallocate (s%f)
      deallocate (s%scales(3)%f)
c
      s%f=>f
      s%scales(3)%f=>p
c
      return
      end
c#######################################################################
      subroutine read_vignetting_function (fname)
c
c-----------------------------------------------------------------------
c
c ****** Read the vignetting function from file FNAME.
c
c-----------------------------------------------------------------------
c
      use ident
      use sds_def
      use params
      use vignetting_function
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      character(*) :: fname
c
c-----------------------------------------------------------------------
c
      type(sds) :: s
      integer :: ierr,i,j
c
c-----------------------------------------------------------------------
c
      if (verbose) then
        write (*,*)
        write (*,*) '### Reading the vignetting function'//
     &              ' from file: ',trim(fname)
      end if
c
c ****** Read the vignetting function from the specified file.
c
      call rdhdf (fname,s,ierr)
c
      if (ierr.ne.0) then
        write (*,*)
        write (*,*) '### ERROR in ',cname,':'
        write (*,*) '### Could not read the vignetting function.'
        write (*,*) 'IERR (from RDHDF) = ',ierr
        write (*,*) 'File name: ',trim(fname)
        call exit (1)
      end if
c
c ****** Check that it contains a 1D field.
c
      if (s%ndim.ne.1) then
        write (*,*)
        write (*,*) '### ERROR in ',cname,':'
        write (*,*) '### The vignetting function file does'//
     &              ' not contain a 1D field.'
        write (*,*) 'Number of dimensions = ',s%ndim
        write (*,*) 'File name: ',trim(fname)
        call exit (1)
      end if
c
c ****** Check that it contains scales.
c
      if (.not.s%scale) then
        write (*,*)
        write (*,*) '### ERROR in ',cname,':'
        write (*,*) '### The vignetting function file does'//
     &              ' not contain scales.'
        write (*,*) 'File name: ',trim(fname)
        call exit (1)
      end if
c
      nvf=s%dims(1)
c
c ****** Check that it contains two or more points.
c
      if (nvf.lt.2) then
        write (*,*)
        write (*,*) '### ERROR in ',cname,':'
        write (*,*) '### The vignetting function file has'//
     &              ' less than two points.'
        write (*,*) 'Number of points = ',nvf
        write (*,*) 'File name: ',trim(fname)
        call exit (1)
      end if
c
c ****** Transfer the vignetting function to 1D arrays for
c ****** convenience.
c
      allocate (r_vf(nvf))
      allocate (vf(nvf))
c
      r_vf=s%scales(1)%f
      vf=s%f(:,1,1)
c
      call deallocate_sds (s)
c
c ****** Check that the vignetting function scale is monotonically
c ****** increasing.
c
      do i=1,nvf-1
        if (r_vf(i+1).le.r_vf(i)) then
          write (*,*)
          write (*,*) '### ERROR in ',cname,':'
          write (*,*) '### The vignetting function scale is not'//
     &                ' monotonically increasing.'
          write (*,*) 'File name: ',trim(fname)
          write (*,*) 'Number of points = ',nvf
          write (*,*) 'Radial locations of points:'
          write (*,*) 'Index',char(9),'Radius'
          do j=1,nvf
            write (*,*) j,char(9),r_vf(j)
          enddo
          write (*,*) 'Not monotonic at index = ',i+1
          call exit (1)
        end if
      enddo
c
c ****** Compute the spline coefficients if needed.
c
      if (cubic) then
        if (verbose) then
          write (*,*)
          write (*,*) '### Computing cubic spline coefficients for vf'
        end if
        call compute_spline_1d (nvf,
     &                          r_vf,
     &                          vf,vf_spl1)
      endif
c
      return
      end
c#######################################################################
      subroutine write_image (filename, fieldname, x, y, f)
c
c-----------------------------------------------------------------------
c
c ****** Write a 2D image (wrapper for wrdhdf_2d and error checking).
c
c-----------------------------------------------------------------------
c
      use ident
      use number_types
      use params
      use image_region
      use image_fields
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      character(512) :: filename
      character(*) :: fieldname
      real(r_typ), dimension(nx) :: x
      real(r_typ), dimension(ny) :: y
      real(r_typ), dimension(nx,ny) :: f
c
c-----------------------------------------------------------------------
c
      integer :: ierr
c
c-----------------------------------------------------------------------
c
      if (verbose) then
        write (*,*)
        write (*,*) '### Writing ',trim(fieldname),' to file: ',
     &               trim(filename)
      end if
c
        call wrhdf_2d (filename,.true.,nx,ny,f,x,y,hdf32,ierr)
c
      if (ierr.ne.0) then
        write (*,*)
        write (*,*) '### ERROR in ',cname,':'
        write (*,*) '### Could not write the ',trim(fieldname),
     &              ' output file:'
        write (*,*) 'IERR (from WRHDF_2D) = ',ierr
        write (*,*) 'File name: ',trim(fieldname)
        call exit (1)
      end if
c
      end
c#######################################################################
      subroutine vignette (f)
c
c-----------------------------------------------------------------------
c
c ****** Multiply the field F by the vignetting function, using
c ****** the radius R_MIN.
c
c-----------------------------------------------------------------------
c
      use ident
      use number_types
      use params
      use image_region
      use image_fields
      use vignetting_function
      use evaluate_spline_1d_interface
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      real(r_typ), dimension(nx,ny) :: f
c
c-----------------------------------------------------------------------
c
      real(r_typ), parameter :: one=1._r_typ
c
c-----------------------------------------------------------------------
c
      integer :: ierr,i,j,ivf,ivfp1
      real(r_typ) :: rv,alpha,vfv
c
c-----------------------------------------------------------------------
c
      do j=1,ny
        do i=1,nx
          rv=r_min(i,j)
          if (rv.lt.rocc) cycle
          if (rv.le.r_vf(1)) then
            vfv=vf(1)
          else if (rv.ge.r_vf(nvf)) then
            vfv=vf(nvf)
          else
            if (cubic) then
              vfv=evaluate_spline_1d(vf_spl1,rv)
            else
              call interp (nvf,r_vf,rv,ivf,ivfp1,alpha,ierr)
              if (ierr.ne.0) then
                write (*,*)
                write (*,*) '### ERROR in ',cname,':'
                write (*,*) '### An error occurred while applying'//
     &                      ' the vignetting function.'
                write (*,*) 'Please check the vignetting'//
     &                      ' function definition.'
                call exit (1)
              end if
              vfv=(one-alpha)*vf(ivf)+alpha*vf(ivfp1)
            endif
          end if
          f(i,j)=vfv*f(i,j)
        enddo
      enddo
c
      return
      end
c#######################################################################
      subroutine multiply_by_radial_power (f)
c
c-----------------------------------------------------------------------
c
c ****** Multiply the field F by R_MIN**POWER.
c
c-----------------------------------------------------------------------
c
      use ident
      use number_types
      use params
      use image_region
      use image_fields
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      real(r_typ), dimension(nx,ny) :: f
c
c-----------------------------------------------------------------------
c
      integer :: i,j
      real(r_typ) :: rv
c
c-----------------------------------------------------------------------
c
      do j=1,ny
        do i=1,nx
          rv=r_min(i,j)
          if (rv.lt.rocc) cycle
          if (power.lt.0.) then
            f(i,j)=f(i,j)/rv**abs(power)
          else
            f(i,j)=f(i,j)*rv**power
          end if
        enddo
      enddo
c
      return
      end
c#######################################################################
      subroutine transform_position (x,y,z,r,t,p)
c
c-----------------------------------------------------------------------
c
c ****** Apply the transformations due to the solar P angle and the
c ****** solar B0 angle, and rotate to the specified central meridian.
c ****** Return the transformed position in spherical coordinates
c ****** (R,T,P).
c
c-----------------------------------------------------------------------
c
      use number_types
      use constants
      use geometry
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      real(r_typ) :: x,y,z
      real(r_typ) :: r,t,p
c
c-----------------------------------------------------------------------
c
      real(r_typ), parameter :: zero=0.
c
c-----------------------------------------------------------------------
c
      real(r_typ) :: x1,y1,z1,x2,y2,z2
c
c-----------------------------------------------------------------------
c
      real(r_typ), external :: fold
c
c-----------------------------------------------------------------------
c
c ****** Apply the solar P angle transformation.
c
      x1=x
      y1= cos_p*y+sin_p*z
      z1=-sin_p*y+cos_p*z
c
c ****** Apply the solar B0 angle transformation.
c
      x2= cos_b0*x1-sin_b0*z1
      y2=y1
      z2= sin_b0*x1+cos_b0*z1
c
c ****** Convert to spherical coordinates.
c
      call c2s (x2,y2,z2,r,t,p)
c
c ****** Transform to the specified central meridian longitude.
c
      p=fold(zero,twopi,p+cml)

      return
      end
c#######################################################################
      function fold (x0,x1,x)
c
c-----------------------------------------------------------------------
c
c ****** "Fold" X into the periodic interval [X0,X1].
c
c ****** On return, X is such that X0.le.X.lt.X1.
c
c-----------------------------------------------------------------------
c
c ****** It is assumed that X0 does not equal X1, as is physically
c ****** necessary.  If X0 and X1 are equal, the routine just
c ****** returns with FOLD=X.
c
c-----------------------------------------------------------------------
c
      use number_types
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      real(r_typ) :: fold
      real(r_typ) :: x0,x1,x
c
c-----------------------------------------------------------------------
c
      real(r_typ) :: xl
c
c-----------------------------------------------------------------------
c
      fold=x
c
      if (x0.eq.x1) return
c
      xl=x1-x0
c
      fold=mod(x-x0,xl)+x0
c
      if (fold.lt.x0) fold=fold+xl
      if (fold.ge.x1) fold=fold-xl
c
      return
      end
c#######################################################################
      subroutine get_kernels (r_min,rv,tv,pv,
     &                        compute_pb,compute_b,
     &                        pb_local,b_local,ds)
c
c-----------------------------------------------------------------------
c
c ****** Calculate the contribution to the polarized brightness pB
c ****** and the brightness B at a point (RV,TV,PV) in spherical
c ****** coordinates.
c
c-----------------------------------------------------------------------
c
c ****** This formulation is adapted from Fran Bagenal, after the
c ****** treatment in Billings (1966) and Altschuler and Perry (1972).
c
c-----------------------------------------------------------------------
c
c ****** When COMPUTE_PB=.true., the polarization brightness pB is
c ****** returned in PB_LOCAL.
c
c ****** When COMPUTE_B=.true., the brightness B is returned
c ****** in B_LOCAL.
c
c ****** The smallest mesh spacing at the interpolation point,
c ****** calculated in Cartesian coordinates, is returned in DS.
c
c ****** R_MIN is the closest distance to the Sun along the line
c ****** of sight.  This point defines the Thomson sphere.
c ****** In the case of an infinitely-distant observer, R_MIN
c ****** becomes the plane-of-sky-radius (called rho by Altschuler
c ****** and Perry).
c
c ****** RV and R_MIN are in units of solar radii.
c
c-----------------------------------------------------------------------
c
      use number_types
      use geometry
      use mas_fields
      use params, only: he_rho
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      real(r_typ) :: r_min
      real(r_typ) :: rv,tv,pv
      logical :: compute_pb,compute_b
      real(r_typ) :: pb_local
      real(r_typ) :: b_local
      real(r_typ) :: ds
c
c-----------------------------------------------------------------------
c
      real(r_typ), parameter :: two=2._r_typ
c
c-----------------------------------------------------------------------
c
      real(r_typ) :: n_e,cf,ef
      real(r_typ), dimension(3) :: drtp
c
c-----------------------------------------------------------------------
c
      real(r_typ), external :: ds_s2c
      real(r_typ), external :: cfunc_billings
      real(r_typ), external :: efunc_billings
c
c-----------------------------------------------------------------------
c
c ****** Interpolate the density.
c
      call interp_field_ds (rho%sds,rho%invtab,rv,tv,pv,n_e,drtp,
     &                   rho%spl2,rho%spl3)
c
c ****** Convert N_E to physical units.
c
      n_e=n_e*fn_n/he_rho
c
c ****** Convert the mesh spacing from spherical coordinates
c ****** to a Cartesian spacing.
c
      ds=ds_s2c(rv,tv,drtp,twodee)
c
c ****** Evaluate the "Billings" radial functions.
c
      if (compute_pb.or.compute_b) then
        cf=cfunc_billings(rv)
      end if
c
      if (compute_b) then
        ef=efunc_billings(rv)
      end if
c
c ****** Obtain the local contribution to pB.
c
      if (compute_pb) then
        pb_local=cf*(r_min/rv)**2*n_e
      else
        pb_local=0.
      end if
c
c ****** Obtain the local contribution to B.
c
      if (compute_b) then
        b_local=(two*ef-cf*(r_min/rv)**2)*n_e
      else
        b_local=0.
      end if
c
      return
      end
c#######################################################################
      function cfunc_billings (r)
c
c-----------------------------------------------------------------------
c
c ****** Coronal light scattering function for polarized brightness
c ****** [Billings 1966].
c
c-----------------------------------------------------------------------
c
c ****** R is the radius in units of R_sun.
c
c-----------------------------------------------------------------------
c
      use number_types
      use params, only: mu_ld
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      real(r_typ) :: r
      real(r_typ) :: cfunc_billings
c
c-----------------------------------------------------------------------
c
      real(r_typ), parameter :: one=1._r_typ
      real(r_typ), parameter :: three=3._r_typ
      real(r_typ), parameter :: eighth=.125_r_typ
c
c-----------------------------------------------------------------------
c
c ****** Intensity coefficient [cm^3], normalized to I0,
c ****** the central disk brightness, 2.49e10 [erg/cm^2/s/sr].
c
      real(r_typ), parameter :: a=8.69e-15_r_typ
c
c-----------------------------------------------------------------------
c
      real(r_typ) :: so,co,so_sq,co_sq
      real(r_typ) :: log_term
      real(r_typ) :: afunc,bfunc
c
c-----------------------------------------------------------------------
c
c ****** Protect against invalid R values.  The physical requirement
c ****** is that R.ge.1.
c
      if (r.le.one) then
        so=one
        so_sq=one
        co=0.
        co_sq=0.
        log_term=0.
      else
        so=one/r
        so_sq=so**2
        co_sq=one-so_sq
        co=sqrt(abs(co_sq))
        log_term=log((one+so)/co)
      end if
c
      afunc=co*so_sq
c
      bfunc=-eighth*( one-three*so_sq
     &               -(co_sq/so)*(one+three*so_sq)*log_term)
c
      cfunc_billings=a*((one-mu_ld)*afunc+mu_ld*bfunc)
c
      return
      end
c#######################################################################
      function efunc_billings (r)
c
c-----------------------------------------------------------------------
c
c ****** Coronal light scattering function for tangential
c ****** polarized brightness [Billings 1966].
c
c-----------------------------------------------------------------------
c
c ****** R is the radius in units of R_sun.
c
c-----------------------------------------------------------------------
c
      use number_types
      use params, only: mu_ld
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      real(r_typ) :: r
      real(r_typ) :: efunc_billings
c
c-----------------------------------------------------------------------
c
      real(r_typ), parameter :: one=1._r_typ
      real(r_typ), parameter :: three=3._r_typ
      real(r_typ), parameter :: four=4._r_typ
      real(r_typ), parameter :: five=5._r_typ
      real(r_typ), parameter :: eighth=.125_r_typ
      real(r_typ), parameter :: one_third=one/three
      real(r_typ), parameter :: four_thirds=four/three
c
c-----------------------------------------------------------------------
c
c ****** Intensity coefficient [cm^3], normalized to I0,
c ****** the central disk brightness, 2.49e10 [erg/cm^2/s/sr].
c
      real(r_typ), parameter :: a=8.69e-15_r_typ
c
c-----------------------------------------------------------------------
c
      real(r_typ) :: so,co,so_sq,co_sq
      real(r_typ) :: log_term
      real(r_typ) :: cfunc,dfunc
c
c-----------------------------------------------------------------------
c
c ****** Protect against invalid R values.  The physical requirement
c ****** is that R.ge.1.
c
      if (r.le.one) then
        so=one
        so_sq=one
        co=0.
        co_sq=0.
        log_term=0.
      else
        so=one/r
        so_sq=so**2
        co_sq=one-so_sq
        co=sqrt(abs(co_sq))
        log_term=log((one+so)/co)
      end if
c
      cfunc=four_thirds-co-one_third*co*co_sq
c
      dfunc=eighth*( five+so_sq
     &              -(co_sq/so)*(five-so_sq)*log_term)
c
      efunc_billings=a*((one-mu_ld)*cfunc+mu_ld*dfunc)
c
      return
      end
c#######################################################################
      subroutine c2s (x,y,z,r,t,p)
c
c-----------------------------------------------------------------------
c
c ****** Convert from Cartesian coordinates (X,Y,Z)
c ****** to spherical coordinates (R,T,P).
c
c ****** This routine returns T and P in radians, in the
c ****** following range:
c
c          0. .le. t .le. pi
c          0. .le. p .lt. 2.*pi
c
c-----------------------------------------------------------------------
c
      use number_types
      use constants
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      real(r_typ) :: x,y,z
      real(r_typ) :: r,t,p
c
c-----------------------------------------------------------------------
c
      r=sqrt(x**2+y**2+z**2)
c
      if (r.eq.0.) then
        t=0.
      else
        t=acos(z/r)
      end if
c
      if (x.eq.0.) then
        if (y.ge.0.) then
          p= halfpi
        else
          p=-halfpi
        end if
      else
        p=atan2(y,x)
      end if
      if (p.lt.0.) p=p+twopi
c
      return
      end
c#######################################################################
      subroutine s2c (r,t,p,x,y,z)
c
c-----------------------------------------------------------------------
c
c ****** Convert from spherical coordinates (R,T,P)
c ****** to Cartesian coordinates (X,Y,Z).
c
c ****** This routine assumes that T and P are in radians.
c
c-----------------------------------------------------------------------
c
      use number_types
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      real(r_typ) :: r,t,p
      real(r_typ) :: x,y,z
c
c-----------------------------------------------------------------------
c
      real(r_typ) :: st
c
c-----------------------------------------------------------------------
c
      st=sin(t)
      x=r*st*cos(p)
      y=r*st*sin(p)
      z=r*cos(t)
c
      return
      end
c#######################################################################
      function ds_s2c (r,t,drtp,axisymmetric)
c
c-----------------------------------------------------------------------
c
c ****** Return the cell size in Cartesian coordinates at the
c ****** spherical position (R,T), calculated as the smallest
c ****** cell size from the mesh spacings in spherical coordinates
c ****** in DRTP.
c
c ****** In the axisymmetric case, AXISYMMETRIC=.true., the phi
c ****** mesh spacing is not present in DRTP, and is not considered.
c
c ****** This routine assumes that all angles are in radians.
c
c-----------------------------------------------------------------------
c
      use number_types
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      real(r_typ) :: r,t
      real(r_typ), dimension(3) :: drtp
      logical :: axisymmetric
      real(r_typ) :: ds_s2c
c
c-----------------------------------------------------------------------
c
      real(r_typ) :: dr,dt,dp,st
c
c-----------------------------------------------------------------------
c
      dr=drtp(1)
      dt=drtp(2)
c
      if (axisymmetric) then
        ds_s2c=min(dr,r*dt)
      else
        dp=drtp(3)
        st=max(sin(t),sin(dt))
        ds_s2c=min(dr,r*dt,r*st*dp)
      end if
c
      return
      end
c#######################################################################
      subroutine interp_field_ds (fld,invtab,r,t,p,fv,ds,sp2,sp3)
c
c-----------------------------------------------------------------------
c
c ****** Interpolate the value of the field in SDS FLD at the point
c ****** (R,T,P) in spherical coordinates.
c
c ****** The value of the interpolated field is returned in FV.
c
c ****** The cell spacing at the interpolation point along each
c ****** dimension is returned in array DS.
c
c ****** If cubic spline interpolation is specified, it calls both
c ****** types (linear and cubic) because the spline library does
c ****** not return the info needed for mesh spacing, and I want to keep
c ****** the "outside" logic as is. There is probably a small speed
c ****** penalty, but it is not worth the complexity to speed it up.
c ****** Also cubic can return negative numbers --> floor fv to zero.
c
c-----------------------------------------------------------------------
c
      use number_types
      use types
      use params, only: cubic
      use spline_def
      use evaluate_spline_2d_interface
      use evaluate_spline_3d_interface
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      type(sds) :: fld
      type(vtab) :: invtab
      real(r_typ) :: r,t,p
      real(r_typ) :: fv
      real(r_typ), dimension(3) :: ds
      type(spl2d) :: sp2
      type(spl3d) :: sp3
c
c-----------------------------------------------------------------------
c
c ****** Note that the interpolated value is set to zero if it is
c ****** outside the bounds of the mesh.
c
      if (fld%ndim.eq.2) then
        call interp_2d (fld%dims(1),
     &                  fld%dims(2),
     &                  fld%scales(1)%f,
     &                  fld%scales(2)%f,
     &                  invtab,
     &                  fld%f,r,t,fv,ds)
c
        if (cubic .and. (fv.ne.0.) ) then
          fv=evaluate_spline_2d(sp2,r,t,invtab%c(1),invtab%c(2))
          fv=max(fv,0.)
        endif
c
      else if (fld%ndim.eq.3) then
        call interp_3d (fld%dims(1),
     &                  fld%dims(2),
     &                  fld%dims(3),
     &                  fld%scales(1)%f,
     &                  fld%scales(2)%f,
     &                  fld%scales(3)%f,
     &                  invtab,
     &                  fld%f,r,t,p,fv,ds)
        if (cubic .and. (fv.ne.0.) ) then
          fv=evaluate_spline_3d(sp3,r,t,p,
     &                          invtab%c(1),invtab%c(2),invtab%c(3))
          fv=max(fv,0.)
        endif
      else
        fv=0.
        ds=huge(fv)
      end if
c
      return
      end
c#######################################################################
      subroutine interp_field (fld,invtab,r,t,p,fv,sp2,sp3)
c
c-----------------------------------------------------------------------
c
c ****** Interpolate the value of the field in SDS FLD at the point
c ****** (R,T,P) in spherical coordinates.
c
c ****** The value of the interpolated field is returned in FV.
c
c-----------------------------------------------------------------------
c
      use number_types
      use types
      use params, only: cubic
      use spline_def
      use evaluate_spline_2d_interface
      use evaluate_spline_3d_interface
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      type(sds) :: fld
      type(vtab) :: invtab
      real(r_typ) :: r,t,p
      real(r_typ) :: fv
      type(spl2d) :: sp2
      type(spl3d) :: sp3
c
      real(r_typ), dimension(3) :: ds
c
c-----------------------------------------------------------------------
c
      if (fld%ndim.eq.2) then
        if (cubic) then
          fv=evaluate_spline_2d(sp2,r,t,invtab%c(1),invtab%c(2))
        else
          call interp_2d (fld%dims(1),
     &                  fld%dims(2),
     &                  fld%scales(1)%f,
     &                  fld%scales(2)%f,
     &                  invtab,
     &                  fld%f,r,t,fv,ds)
        endif
c
      else if (fld%ndim.eq.3) then
        if (cubic) then
          fv=evaluate_spline_3d(sp3,r,t,p,
     &                          invtab%c(1),invtab%c(2),invtab%c(3))
        else
          call interp_3d (fld%dims(1),
     &                  fld%dims(2),
     &                  fld%dims(3),
     &                  fld%scales(1)%f,
     &                  fld%scales(2)%f,
     &                  fld%scales(3)%f,
     &                  invtab,
     &                  fld%f,r,t,p,fv,ds)
        endif
      endif
c
      return
      end
c
c#######################################################################
      subroutine interp_2d (nx,ny,x,y,inv,f,xv,yv,fv,ds)
c
c-----------------------------------------------------------------------
c
c ****** Interpolate the value of the 2D field FV at (XV,YV) from
c ****** array F(NX,NY), defined on the mesh X(NX) x Y(NY).
c
c ****** The structure INV holds the inverse interpolation tables.
c
c ****** If the point (XV,YV) is outside the bounds of the
c ****** X x Y mesh, FV=0. is returned.
c
c ****** The mesh spacing at the interpolation point in each
c ****** dimension is returned in array DS.
c
c-----------------------------------------------------------------------
c
      use ident
      use number_types
      use types
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      integer :: nx,ny
      real(r_typ), dimension(nx) :: x
      real(r_typ), dimension(ny) :: y
      type(vtab) :: inv
      real(r_typ), dimension(nx,ny) :: f
      real(r_typ) :: xv,yv,fv
      real(r_typ), dimension(2) :: ds
c
c-----------------------------------------------------------------------
c
      real(r_typ), parameter :: one=1._r_typ
c
c-----------------------------------------------------------------------
c
      logical :: outside
      integer :: i,j,ip1,jp1
      real(r_typ) :: ax,ay
c
c-----------------------------------------------------------------------
c
      real(r_typ), external :: mesh_spacing
c
c-----------------------------------------------------------------------
c
c ****** Find the cell that contains the interpolation point.
c
      outside=.false.
c
      if (xv.lt.x(1)) then
        outside=.true.
        i=1
        ip1=1
        ax=0.
      else if (xv.gt.x(nx)) then
        outside=.true.
        i=nx
        ip1=nx
        ax=0.
      else
        call interpi (x,nx,inv%c(1),xv,i,ip1,ax)
      end if
c
      if (yv.lt.y(1)) then
        outside=.true.
        j=1
        jp1=1
        ay=0.
      else if (yv.gt.y(ny)) then
        outside=.true.
        j=ny
        jp1=ny
        ay=0.
      else
        call interpi (y,ny,inv%c(2),yv,j,jp1,ay)
      end if
c
c ****** Get the mesh spacing at the interpolation point.
c
      ds(1)=mesh_spacing(nx,x,i,ip1,ax)
      ds(2)=mesh_spacing(ny,y,j,jp1,ay)
c
c ****** If the point is outside the mesh limits, set the
c ****** interpolated value to 0.
c
      if (outside) then
        fv=0.
        write (*,*)
        write (*,*) '### ERROR in ',cname,':'
        write (*,*) '### Error in INTERP_2D:'
        write (*,*) '### A requested location is outside the mesh.'
        write (*,*) '### This should not happen anymore!'
        write (*,*)
        write (*,*) 'The requested location was:'
        write (*,*) 'r: ', xv
        write (*,*) 't: ', yv
        call exit (1)
      else
        fv= (one-ax)*((one-ay)*f(i  ,j  )+ay*f(i  ,jp1))
     &     +     ax *((one-ay)*f(ip1,j  )+ay*f(ip1,jp1))
      end if
c
      return
      end
c#######################################################################
      subroutine interp_3d (nx,ny,nz,x,y,z,inv,f,xv,yv,zv,fv,ds)
c
c-----------------------------------------------------------------------
c
c ****** Interpolate the value of the 3D field FV at (XV,YV,ZV) from
c ****** array F(NX,NY,NZ), defined on the mesh X(NX) x Y(NY) x Z(NZ).
c
c ****** The structure INV holds the inverse interpolation tables.
c
c ****** If the point (XV,YV,ZV) is outside the bounds of the
c ****** X x Y x Z mesh, FV=0. is returned.
c
c ****** The mesh spacing at the interpolation point in each
c ****** dimension is returned in array DS.
c
c-----------------------------------------------------------------------
c
      use ident
      use number_types
      use types
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      integer :: nx,ny,nz
      real(r_typ), dimension(nx) :: x
      real(r_typ), dimension(ny) :: y
      real(r_typ), dimension(nz) :: z
      type(vtab) :: inv
      real(r_typ), dimension(nx,ny,nz) :: f
      real(r_typ) :: xv,yv,zv,fv
      real(r_typ), dimension(3) :: ds
c
c-----------------------------------------------------------------------
c
      real(r_typ), parameter :: one=1._r_typ
c
c-----------------------------------------------------------------------
c
      logical :: outside
      integer :: i,j,k,ip1,jp1,kp1
      real(r_typ) :: ax,ay,az
c
c-----------------------------------------------------------------------
c
      real(r_typ), external :: mesh_spacing
c
c-----------------------------------------------------------------------
c
c ****** Find the cell that contains the interpolation point.
c
      outside=.false.
c
      if (xv.lt.x(1)) then
        outside=.true.
        i=1
        ip1=1
        ax=0.
      else if (xv.gt.x(nx)) then
        outside=.true.
        i=nx
        ip1=nx
        ax=0.
      else
        call interpi (x,nx,inv%c(1),xv,i,ip1,ax)
      end if
c
      if (yv.lt.y(1)) then
        outside=.true.
        j=1
        jp1=1
        ay=0.
      else if (yv.gt.y(ny)) then
        outside=.true.
        j=ny
        jp1=ny
        ay=0.
      else
        call interpi (y,ny,inv%c(2),yv,j,jp1,ay)
      end if
c
      if (zv.lt.z(1)) then
        outside=.true.
        k=1
        kp1=1
        az=0.
      else if (zv.gt.z(nz)) then
        outside=.true.
        k=nz
        kp1=nz
        az=0.
      else
        call interpi (z,nz,inv%c(3),zv,k,kp1,az)
      end if
c
c ****** Get the mesh spacing at the interpolation point.
c
      ds(1)=mesh_spacing(nx,x,i,ip1,ax)
      ds(2)=mesh_spacing(ny,y,j,jp1,ay)
      ds(3)=mesh_spacing(nz,z,k,kp1,az)
c
c ****** If the point is outside the mesh limits, set the
c ****** interpolated value to 0.
c
      if (outside) then
        fv=0.
        write (*,*)
        write (*,*) '### ERROR in ',cname,':'
        write (*,*) '### Error in INTERP_3D:'
        write (*,*) '### A requested location is outside the mesh.'
        write (*,*) '### This should not happen anymore!'
        write (*,*)
        write (*,*) 'The requested location was:'
        write (*,*) 'r: ', xv
        write (*,*) 't: ', yv
        write (*,*) 'p: ', zv
        call exit (1)
      else
        fv= (one-ax)*( (one-ay)*( (one-az)*f(i  ,j  ,k  )
     &                           +     az *f(i  ,j  ,kp1))
     &                +     ay *( (one-az)*f(i  ,jp1,k  )
     &                           +     az *f(i  ,jp1,kp1)))
     &     +     ax *( (one-ay)*( (one-az)*f(ip1,j  ,k  )
     &                           +     az *f(ip1,j  ,kp1))
     &                +     ay *( (one-az)*f(ip1,jp1,k  )
     &                           +     az *f(ip1,jp1,kp1)))
      end if
c
      return
      end
c#######################################################################
      function mesh_spacing (nx,x,i,ip1,alpha)
c
c-----------------------------------------------------------------------
c
c ****** Return the mesh spacing in coordinate X at the interpolation
c ****** point specified by I, IP1, and ALPHA.
c
c ****** This is an internal utility routine that does not do
c ****** any error checking.
c
c-----------------------------------------------------------------------
c
      use number_types
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      integer :: nx
      real(r_typ), dimension(nx) :: x
      integer :: i,ip1
      real(r_typ) :: alpha
      real(r_typ) :: mesh_spacing
c
c-----------------------------------------------------------------------
c
      real(r_typ), parameter :: one=1._r_typ
c
c-----------------------------------------------------------------------
c
      integer :: im,ip
      real(r_typ) :: d,dxm,dxp
c
c-----------------------------------------------------------------------
c
      im=max( 1,i-1)
      ip=min(nx,i+1)
      d=real(ip-im)
      if (d.ne.0.) then
        dxm=(x(ip)-x(im))/d
      else
        dxm=huge(one)
      end if
c
      im=max( 1,ip1-1)
      ip=min(nx,ip1+1)
      d=real(ip-im)
      if (d.ne.0.) then
        dxp=(x(ip)-x(im))/d
      else
        dxp=huge(one)
      end if
c
      mesh_spacing=(one-alpha)*dxm+alpha*dxp
c
      return
      end
c#######################################################################
      subroutine build_inverse_tables (s,inv)
c
c-----------------------------------------------------------------------
c
c ****** Build the inverse interpolation tables INV for the SDS
c ****** in structure S.
c
c ****** These arrays are used to to increase the efficiency
c ****** of interpolation lookups.
c
c-----------------------------------------------------------------------
c
      use number_types
      use types
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      type(sds) :: s
      type(vtab) :: inv
c
c-----------------------------------------------------------------------
c
      integer :: i
c
c-----------------------------------------------------------------------
c
c ****** Use a number of points for the inverse interpolation table
c ****** equal to the number in the original scale.
c
      do i=1,s%ndim
        inv%c(i)%n=s%dims(i)
        allocate (inv%c(i)%f(inv%c(i)%n))
        call getinv (s%scales(i)%f,s%dims(i),inv%c(i))
      enddo
c
      return
      end
c#######################################################################
      subroutine getinv (x,n,tab)
c
c-----------------------------------------------------------------------
c
c ****** Build an inverse interpolation table to increase the
c ****** efficiency of table look-up in a nonuniform mesh.
c
c ****** On input, the table X(N) is specified, together with the
c ****** number of points to use in the inverse interpolation
c ****** table, TAB%N.
c
c ****** The output is a structure TAB with the inverse interpolation
c ****** table.  This structure has the following components:
c
c ******    N:  the number of points in the table (input);
c ******    D:  the inverse of the uniform table spacing;
c ******    F:  the inverse interpolation table.
c
c-----------------------------------------------------------------------
c
      use ident
      use number_types
      use types
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      integer :: n
      real(r_typ), dimension(n) :: x
      type(itab) :: tab
c
c-----------------------------------------------------------------------
c
      real(r_typ), parameter :: one=1._r_typ
c
c-----------------------------------------------------------------------
c
      integer :: i,k,ierr,ip1
      real(r_typ) :: dx,xv,alpha,en
c
c-----------------------------------------------------------------------
c
c ****** For the special case when the table X has only one point,
c ****** the inverse table is not used.  It is thus loaded with
c ****** dummy values.
c
      if (n.eq.1) then
        tab%d=0.
        tab%f=1.
        return
      end if
c
c ****** Check that the number of points is valid.
c
      if (tab%n.le.1) then
        write (*,*)
        write (*,*) '### ERROR in ',cname,':'
        write (*,*) '### Error in GETINV:'
        write (*,*) '### Invalid number of points specified'//
     &              ' for the inverse interpolation table.'
        write (*,*)
        write (*,*) 'Number of points = ',tab%n
        call exit (1)
      end if
c
c ****** Set the uniform interval to be used in the inverse
c ****** interpolation.
c
      dx=(x(n)-x(1))/(tab%n-one)
c
      if (dx.le.0.) then
        write (*,*)
        write (*,*) '### ERROR in ',cname,':'
        write (*,*) '### Error in GETINV:'
        write (*,*) '### Invalid interval for the inverse'//
     &              ' interpolation table.'
        write (*,*)
        write (*,*) 'Interval = ',dx
        call exit (1)
      end if
c
      tab%d=one/dx
c
c ****** Build the inverse interpolation table.
c
      en=n
c
      do k=1,tab%n
        xv=x(1)+(k-one)*dx
        xv=max(xv,x(1))
        xv=min(xv,x(n))
        call interp (n,x,xv,i,ip1,alpha,ierr)
        if (ierr.ne.0) then
          write (*,*)
          write (*,*) '### ERROR in ',cname,':'
          write (*,*) '### Error in GETINV:'
          write (*,*) '### Error in building the inverse'//
     &                ' interpolation table.'
          call exit (1)
        end if
        tab%f(k)=i+alpha
        tab%f(k)=max(tab%f(k),one)
        tab%f(k)=min(tab%f(k),en)
      enddo
c
      return
      end
c#######################################################################
      subroutine interpi (x,n,tab,xv,i,ip1,alpha)
c
c-----------------------------------------------------------------------
c
c ****** Get interpolation factor ALPHA and table indices I and IP1.
c
c ****** This routine does not do the actual interpolation.  Use the
c ****** returned values of I, IP1, and ALPHA to get the
c ****** interpolated value.
c
c ****** This version uses an inverse interpolation table, TAB,
c ****** to improve the efficiency of the search.
c
c-----------------------------------------------------------------------
c
      use ident
      use number_types
      use types
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      integer :: n
      real(r_typ), dimension(n) :: x
      type(itab) :: tab
      real(r_typ) :: xv
      integer :: i
      integer :: ip1
      real(r_typ) :: alpha
      intent(in) :: x,n,tab,xv
      intent(out) :: i,ip1,alpha
c
c-----------------------------------------------------------------------
c
      real(r_typ), parameter :: one=1._r_typ
c
c-----------------------------------------------------------------------
c
      integer :: ig
      real(r_typ) :: xi,fiv
c
c-----------------------------------------------------------------------
c
c ****** For the special case when the table has only one point,
c ****** the inverse table is not used.  In this case it is
c ****** necessary for XV to equal X(1) exactly, otherwise this
c ****** routine exits with an error.
c
      if (n.eq.1) then
        if (xv.eq.x(1)) then
          i=1
          ip1=1
          alpha=0.
        else
          go to 900
        end if
      end if
c
c ****** Get an estimate of the nearest grid point location in
c ****** the (uniform) inverse interpolation table.
c
      xi=one+(xv-x(1))*tab%d
      i=xi
      i=max(i,1)
      i=min(i,tab%n-1)
      alpha=xi-i
      fiv=(one-alpha)*tab%f(i)+alpha*tab%f(i+1)
c
c ****** Set IG to be the guess for the nearest grid point.
c
      ig=fiv
      ig=max(ig,1)
      ig=min(ig,n-1)
c
      if (xv.ge.x(ig)) then
c
c ****** Search forwards.
c
        do i=ig,n-1
          if (xv.ge.x(i).and.xv.le.x(i+1)) then
            alpha=(xv-x(i))/(x(i+1)-x(i))
            ip1=i+1
            return
          end if
        enddo
c
      else
c
c ****** Search backwards.
c
        do i=ig-1,1,-1
          if (xv.ge.x(i).and.xv.le.x(i+1)) then
            alpha=(xv-x(i))/(x(i+1)-x(i))
            ip1=i+1
            return
          end if
        enddo
c
      end if
c
c ****** ERROR: value not found in table.
c
  900 continue
c
      write (*,*)
      write (*,*) '### ERROR in ',cname,':'
      write (*,*) '### Error in INTERPI:'
      write (*,*) 'Abscissa not found in table.'
      write (*,*)
      write (*,*) 'N = ',n
      write (*,*) 'X = ',(x(i),i=1,n)
      write (*,*)
      write (*,*) 'Abscissa requested = ',xv
      call exit (1)
c
      return
      end
c#######################################################################
      subroutine interp (n,x,xv,i,ip1,a,ierr)
c
c-----------------------------------------------------------------------
c
c ****** Get the interpolation factor at XV from the table X(N).
c
c-----------------------------------------------------------------------
c
c ****** This routine does not do the actual interpolation.  Use the
c ****** returned values of I, IP1, and A to get the interpolated
c ****** value.
c
c-----------------------------------------------------------------------
c
      use number_types
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      integer :: n
      real(r_typ), dimension(n) :: x
      real(r_typ) :: xv
      integer :: i,ip1
      real(r_typ) :: a
      integer :: ierr
c
c-----------------------------------------------------------------------
c
      ierr=0
c
c ****** For the special case when the table has only one point,
c ****** it is necessary for XV to equal X(1) exactly, otherwise
c ****** this routine exits with an error.
c
      if (n.eq.1) then
        if (xv.eq.x(1)) then
          i=1
          ip1=1
          a=0.
          return
        else
          go to 900
        end if
      end if
c
c ****** Find the interval and compute the interpolation factor.
c
      do i=1,n-1
        if (xv.ge.x(i).and.xv.le.x(i+1)) then
          ip1=i+1
          if (x(i).eq.x(i+1)) then
            a=0.
          else
            a=(xv-x(i))/(x(i+1)-x(i))
          end if
          return
        end if
      enddo
c
c ****** ERROR: the value was not found.
c
  900 continue
c
      write (*,*)
      write (*,*) '### ERROR in INTERP:'
      write (*,*) '### Value not found in table.'
      write (*,*) 'Value requested = ',xv
      write (*,*) 'Min table value = ',x(1)
      write (*,*) 'Max table value = ',x(n)
      write (*,*) 'Number of values in table = ',n
      ierr=1
c
      return
      end
c#######################################################################
      subroutine rtp2xyz (r,t,p,x,y,z)
c
c-----------------------------------------------------------------------
c
c ****** Convert spherical position to cartesian position.
c
c-----------------------------------------------------------------------
c
c-----------------------------------------------------------------------
c
      use number_types
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      real(r_typ) :: r, t, p
      real(r_typ) :: x, y, z
c
c-----------------------------------------------------------------------
c
      real(r_typ) :: cph, sph, cth, sth
c
c-----------------------------------------------------------------------
c
      cph=cos(p)
      sph=sin(p)
      cth=cos(t)
      sth=sin(t)
      x = r*cph*sth
      y = r*sph*sth
      z = r*cth
c
      return
      end
c#######################################################################
      subroutine svtocv (t,p,vr,vt,vp,vx,vy,vz)
c
c-----------------------------------------------------------------------
c
c ****** Rotate Spherical Vector Components to Cartesian Components
c
c-----------------------------------------------------------------------
c
c-----------------------------------------------------------------------
c
      use number_types
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      real(r_typ) :: t, p
      real(r_typ) :: vr, vt, vp
      real(r_typ) :: vx, vy, vz
c
c-----------------------------------------------------------------------
c
      real(r_typ) :: cph, sph, cth, sth
c
c-----------------------------------------------------------------------
c
      cph=cos(p)
      sph=sin(p)
      cth=cos(t)
      sth=sin(t)
      vx = vr*sth*cph + vt*cth*cph - vp*sph
      vy = vr*sth*sph + vt*cth*sph + vp*cph
      vz = vr*cth - vt*sth
c
      return
      end
c#######################################################################
      subroutine normalize(vec_in, vec_out)
c
c-----------------------------------------------------------------------
c
c ****** Normalize a vector to unit length.
c
c-----------------------------------------------------------------------
c
      use number_types
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      real(r_typ), dimension(3) :: vec_out
      real(r_typ), dimension(3) :: vec_in
      real(r_typ) :: mag
c
c-----------------------------------------------------------------------
c
      mag = sqrt(vec_in(1)**2 + vec_in(2)**2 + vec_in(3)**2)
      vec_out(1) = vec_in(1)/mag
      vec_out(2) = vec_in(2)/mag
      vec_out(3) = vec_in(3)/mag
c
      return
      end
c#######################################################################
      subroutine cross(vec1,vec2,vec_out)
c
c-----------------------------------------------------------------------
c
c ****** Compute the cross product of two vectors.
c
c-----------------------------------------------------------------------
c
      use number_types
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      real(r_typ), dimension(3) :: vec_out
      real(r_typ), dimension(3) :: vec1,vec2
c
c-----------------------------------------------------------------------
c
      vec_out(1) = vec1(2)*vec2(3) - vec1(3)*vec2(2)
      vec_out(2) = vec1(3)*vec2(1) - vec1(1)*vec2(3)
      vec_out(3) = vec1(1)*vec2(2) - vec1(2)*vec2(1)
c
      return
      end
c#######################################################################
      subroutine get_vector_projection(ex,ey,ez,rv,tv,pv,vec)
c
c-----------------------------------------------------------------------
c
c ****** Get the vector field at this location and project it into
c ****** the observing frame defined with unit vectors ex, ey, ez.
c
c-----------------------------------------------------------------------
c
      use number_types
      use mas_fields
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      real(r_typ), dimension(3) :: ex,ey,ez
      real(r_typ) :: rv, tv, pv
      real(r_typ), dimension(3) :: vec
c
c-----------------------------------------------------------------------
c
      real(r_typ) :: vrv, vtv, vpv
      real(r_typ) :: vx, vy, vz
c
c-----------------------------------------------------------------------
c
c ****** Interpolate the velocity components.
c
      call interp_field (vr%sds,vr%invtab,rv,tv,pv,vrv,vr%spl2,vr%spl3)
      call interp_field (vt%sds,vt%invtab,rv,tv,pv,vtv,vt%spl2,vt%spl3)
      call interp_field (vp%sds,vp%invtab,rv,tv,pv,vpv,vp%spl2,vp%spl3)
c
c ****** Get the cartesian vector components.
c
      call svtocv( tv, pv, vrv, vtv, vpv, vx, vy, vz)
c
c ****** Project the vector into the observer frame.
c
      vec(1) = ex(1)*vx + ex(2)*vy + ex(3)*vz
      vec(2) = ey(1)*vx + ey(2)*vy + ey(3)*vz
      vec(3) = ez(1)*vx + ez(2)*vy + ez(3)*vz
c
      return
      end
c#######################################################################
      subroutine get_unit_vectors(losx,losy,losz,ex,ey,ez)
c
c-----------------------------------------------------------------------
c
c ****** Compute the unit vectors of the observer frame.
c
c ****** los_vec is a vector pointing to the observer from a point
c ****** along the LOS (doesn't need to be normalized).
c
c ****** ex, ey, ez are the unit vectors of the observer frame.
c
c-----------------------------------------------------------------------
c
      use number_types
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      real(r_typ), parameter :: zero=0._r_typ
      real(r_typ), parameter :: one=1._r_typ
c
c-----------------------------------------------------------------------
c
      real(r_typ), intent(in) :: losx, losy, losz
      real(r_typ), dimension(3), intent(out) :: ex, ey, ez
c
c-----------------------------------------------------------------------
c
      real(r_typ), dimension(3) :: vec
c
c-----------------------------------------------------------------------
c
c ****** x-hat is sun-to-observer unit vector --> normalize los_vec.
c
      vec(1)=losx
      vec(2)=losy
      vec(3)=losz
      call normalize(vec,ex)
c
c ****** y-hat is perp to heliographic z axis and observer x-hat.
c
      vec(1)=zero
      vec(2)=zero
      vec(3)=one
      call cross(vec,ex,ey)
      call normalize(ey,ey)
c
c ****** z-hat is perp to heliographic x-hat and y-hat.
c
      call cross(ex,ey,ez)
      call normalize(ez,ez)
c
      return
      end
c#######################################################################
      subroutine add_to_integrals(i,j,e_pb,e_b,ds,ex,ey,ez,rv,tv,pv,x)
c
c-----------------------------------------------------------------------
c
c ****** Add the current contribution to the weighted integrals.
c
c ****** Everything is here so it can be used in both integration
c ****** algoriths and skipped if not requested.
c
c ****** i, j are the current pixel indexes.
c
c ****** e_pb and e_b are the local pB and B emissivities.
c
c ****** ds is the path length contribution for the integral.
c
c ****** ex, ey, ez are the unit vectors of the observer frame.
c
c ****** rv, tv, pv are the current coordinates
c
c ****** x is the signed distance from r_min along the LOS (positive is
c ****** towards the observer, negative is away.
c
c-----------------------------------------------------------------------
c
      use number_types
      use params
      use constants
      use image_fields
      use mas_fields
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      integer, intent(in) :: i, j
      real(r_typ), intent(in) :: e_pb, e_b, ds
      real(r_typ), dimension(3), intent(in) :: ex, ey, ez
      real(r_typ), intent(in) :: rv, tv, pv, x
c
c-----------------------------------------------------------------------
c
      real(r_typ) :: var_now
      real(r_typ) :: e_now
      real(r_typ) :: angle_now
      real(r_typ), dimension(3) :: vec
c
c-----------------------------------------------------------------------
c
c ****** Select the weighting.
c
      if (weight_integral_by_b) then
        e_now=e_b
      else
        e_now=e_pb
      endif
c
c ****** Scalar field.
c
      if (compute_scalar_integration) then
        call interp_field (sfield%sds,sfield%invtab,rv,tv,pv,
     &            var_now,sfield%spl2,sfield%spl3)
        sf_avg(i,j)=sf_avg(i,j) + ds*e_now*var_now
      endif
c
c ****** Average LOS angle.
c
      if (compute_angle_integration) then
        angle_now = atan(x/r_min(i,j))/degtorad
        angle_avg(i,j)=angle_avg(i,j) + ds*e_now*angle_now
      endif
c
c ****** Vector field. Mapping from observer frame to image frame is:
c        Observer_x -> Image_LOS
c        Observer_y -> Image_x
c        Observer_z -> Image_y
c
      if (compute_vector_integration) then
         call get_vector_projection(ex,ey,ez,rv,tv,pv,vec)
         vlos_avg(i,j) = vlos_avg(i,j) + ds*e_now*vec(1)
         vx_avg(i,j) = vx_avg(i,j) + ds*e_now*vec(2)
         vy_avg(i,j) = vy_avg(i,j) + ds*e_now*vec(3)
      endif
c
      return
      end
c#######################################################################
      subroutine normalize_integrals
c
c-----------------------------------------------------------------------
c
c ****** Normalize the weighted integrals to get averages.
c
c ****** Pixels that are occulted or outside the grid get the disk val.
c
c-----------------------------------------------------------------------
c
      use number_types
      use params
      use geometry
      use image_region
      use image_fields
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      real(r_typ), parameter :: one=1._r_typ
c
c-----------------------------------------------------------------------
c
      integer :: i, j
      real(r_typ) :: e_now
c
c-----------------------------------------------------------------------
c
c ****** Loop over pixels, normalize the integrals.
c
      do j=1,ny
        do i=1,nx

          if (r_min(i,j).le.max(rocc,one).or.r_min(i,j).ge.rmax) then

            if (compute_scalar_integration) sf_avg(i,j)=disk_value
            if (compute_angle_integration) angle_avg(i,j)=disk_value
            if (compute_vector_integration) then
              vlos_avg(i,j)=disk_value
              vx_avg(i,j)=disk_value
              vy_avg(i,j)=disk_value
            endif
            cycle
          endif
c
c ****** Divide out emissivity integral to get averages.
c
          if (do_integral_weighting) then
            if (weight_integral_by_b) then
              e_now=b(i,j)
            else
              e_now=pb(i,j)
            endif
c
            if (compute_scalar_integration) then
              sf_avg(i,j)=sf_avg(i,j)/e_now
            endif
c
            if (compute_angle_integration) then
              angle_avg(i,j)=angle_avg(i,j)/e_now
            endif
c
            if (compute_vector_integration) then
              vlos_avg(i,j)=vlos_avg(i,j)/e_now
              vx_avg(i,j)=vx_avg(i,j)/e_now
              vy_avg(i,j)=vy_avg(i,j)/e_now
            endif
          endif
        enddo
      enddo
c
      return
      end
c#######################################################################
      subroutine python_getpb(py_rho, py_b, py_pb, py_help, py_verbose, py_cubic, py_oldmas, py_long, py_p, py_b0, py_r, py_nx, py_ny, py_x0, py_x1, py_y0, py_y1, py_wispr1, py_wispr2, py_rocc, py_dsmult, py_power, py_disk, py_vf, py_vr, py_vt, py_vp, py_scalar, py_avg_scalar, py_avg_los_angle, py_avg_vlos, py_avg_vx, py_avg_vy, py_avg_using_b, py_he_frac, py_mu)
c
c-----------------------------------------------------------------------
c
c ****** Set parameters from the command line arguments.
c
c-----------------------------------------------------------------------
c
      use ident
      use number_types
      use syntax
      use paragraph_def
      use get_usage_line_interface
      use print_par_interface
      use delete_par_interface
      use params
      use image_region
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c ****** Storage the for usage line.
c
      type(paragraph), pointer :: usage
c
c ****** Storage for the error message.
c
      character(72) :: errmsg
c
c-----------------------------------------------------------------------
c
      integer :: ierr, nargs
      character(512) :: arg
      logical :: set
      logical, save :: print_help=.false.

c
c-----------------------------------------------------------------------
c
      integer, external :: intval
      real(r_typ), external :: fpval
c
c-----------------------------------------------------------------------
c
c
c
c ****** Parse the command line.
c
c
c ****** Check if the "-help" option was specified.
c
!       call nargs_specified (nargs)
!       call fetcharg ('-help',set,arg)
!
        if (py_help) then
          print_help=.true.
        else
          print_help=.false.
        end if

c ****** Set the parameters.
c
c ****** Verbose flag.
c
      if (py_verbose) then
        verbose=.true.
      else
        verbose=.false.
      end if
c
c ****** Cubic interpolation flag.
c
      if (py_cubic) then
        cubic=.true.
      else
        cubic=.false.
      end if
c
c ****** Switch to indicate old MAS code data files.
c
      if (py_oldmas) then
        oldmas=.true.
      else
        oldmas=.false.
      end if
c
c ****** Carrington rotation longitude.
c

      crlong=py_long
c
c ****** Solar P angle.
c
      pangle=py_p
c
c ****** Solar B0 angle.
c
      b0angle=py_b0
c
c ****** Distance of the observer from the Sun [AU].
c
      if (py_r == 0) then
        r_obs_set=.false.
        r_obs=1._r_typ
      else
         r_obs_set=.true.
        r_obs=py_r
      end if
c
c ****** Number of points to use for the image.
c
      nx=py_nx
c
      ny=py_ny
c
c ****** Image dimensions.
c
      image_limits_set=.false.
c
      if (py_x0 /= -3) image_limits_set=.true.
      x0=py_x0
c
      if (py_x1 /= 3) image_limits_set=.true.
      x1=py_x1
c
      if (py_y0 /= -3) image_limits_set=.true.
      y0=py_y0
c
      if (py_y1 /= 3) image_limits_set=.true.
      y1=py_y1
c
c ****** WISPR camera 1 flag.
c
      if (py_wispr1) then
        wispr1=.true.
      else
         wispr1=.false.
      end if

c
c ****** WISPR camera 2 flag.
c
      if (py_wispr2) then
        wispr2=.true.
      else
         wispr2=.false.
      end if
c
c ****** Occulting disk radius.
c
      rocc=py_rocc
c
c ****** Factor by which to multiply the LOS integration step size.
c
      dsmult=py_dsmult
c
c ****** Radial power exponent.
c
      if (py_power /= 0) then
        power_set=.true.
        power=py_power
      else
        power_set=.false.
        power=0
      endif
c
c ****** Value for pB on the solar disk.
c
      disk_value=py_disk
c
c ****** Vignetting function file name.
c
      vf_file=trim(py_vf)
      if (vf_file.eq.'<none>') then
        vf_file=' '
      end if
c
c ****** Density file name.
c
      rho_file=trim(py_rho)
c
c ****** Vector field r component file name (input).
c
      vr_file=trim(py_vr)
      if (vr_file.eq.'<none>') then
        vr_file=' '
      else
        compute_vector_integration=.true.
        vr_file=trim(py_vr)

      end if
c
c ****** Vector field t component file name (input).
c
      vt_file=trim(py_vt)
      if (vt_file.eq.'<none>') then
        vt_file=' '
      else
        compute_vector_integration=.true.
        vt_file=trim(py_vt)

      end if
c
c ****** Vector field p component file name (input).
c
      vp_file=trim(py_vp)
      if (vp_file.eq.'<none>') then
        vp_file=' '
      else
        compute_vector_integration=.true.
        vp_file=trim(py_vp)

      end if
c
c ****** Scalar field file name (input).
c
      scalar_file=trim(py_scalar)
      if (scalar_file.eq.'<none>') then
        compute_scalar_integration=.false.
        scalar_file=' '
      else
        compute_scalar_integration=.true.
        scalar_file=trim(py_scalar)
      end if
c
c ****** Emissivity weighted LOS integral of scalar field.
c
      scalar_file_out=trim(py_avg_scalar)
      if (scalar_file_out.eq.'<none>') then
        scalar_file_out=' '
      else
        scalar_file_out=trim(py_avg_scalar)
      end if
c
c ****** Emissivity weighted LOS contribution angle.
c
      angle_file_out=trim(py_avg_los_angle)
      if (angle_file_out.eq.'<none>') then
        angle_file_out=' '
      else
        compute_angle_integration=.true.
        angle_file_out=trim(py_avg_los_angle)
      end if
c
c ****** Emissivity weighted LOS integral of vector field (vlos).
c
      vlos_file_out=trim(py_avg_vlos)
      if (vlos_file_out.eq.'<none>') then
        vlos_file_out=' '
      else
        vlos_file_out=trim(py_avg_vlos)
      end if
c
c ****** Emissivity weighted LOS integral of vector field (vx).
c
      vx_file_out=trim(py_avg_vx)
      if (vx_file_out.eq.'<none>') then
        vx_file_out=' '
      else
        vx_file_out=trim(py_avg_vx)
      end if

c
c ****** Emissivity weighted LOS integral of vector field (vy).
c
      vy_file_out=trim(py_avg_vy)
      if (vy_file_out.eq.'<none>') then
        vy_file_out=' '
      else
        vy_file_out=trim(py_avg_vy)
      end if

c
c ****** Flag to weight integral averages with B and not pB.
c
      if (py_avg_using_b) then
        weight_integral_by_b=.true.
      else
        avg_using_b=.false.
      end if
c
c ****** Helium Fraction
c
      he_frac=py_he_frac
      he_rho=(1._r_typ+4._r_typ*he_frac)/(1._r_typ+2._r_typ*he_frac)
c
c ****** Value for limb-darkening coefficient.
c
      mu_ld=py_mu
c
c ****** Output pB file name.
c
      if (py_pb.eq.'<none>') then
        compute_pb=.false.
        pb_file=' '
      else
        compute_pb=.true.
        pb_file=trim(arg)
      end if

c
c ****** Output B file name.
c
      if (py_pb.eq.'<none>') then
        compute_b=.false.
        b_file=' '
      else
        compute_b=.true.
        b_file=trim(arg)
      end if
c
      return
      end
c#######################################################################
      subroutine print_help_text
c
c-----------------------------------------------------------------------
c
c ****** Print the help / syntax information.
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
        write (*,*)
        write (*,*)
        write (*,*) '### Help for GETPB Special Integration Methods:'
        write (*,*)
        write (*,*) 'GETPB can now compute emissivity weighted'//
     &              ' averages of other quantities.'
        write (*,*) 'Here are the additions flags/options: '
        write (*,*)
        write (*,*)
        write (*,*) 'SCALAR FIELD AVERAGING:'
        write (*,*)
        write (*,*) ' Use -scalar <file> to specify a 2D or 3D file'//
     &              ' with a scalar field.'
        write (*,*)
        write (*,*) ' GETPB will integrate the scalar along the LOS'//
     &              ' to get the weighted average.'
        write (*,*)
        write (*,*) ' Specify the output file with -avg_scalar <file>'
        write (*,*)
        write (*,*)
        write (*,*) 'VECTOR FIELD COMPONENT AVERAGING:'
        write (*,*)
        write (*,*) ' Use -vr <file>, -vt <file>, -vp <file> to'//
     &              ' specify a 2D or 3D vector field.'
        write (*,*)
        write (*,*) ' GETPB will project the rtp vector field into'//
     &              ' LOS image coordinates and'
        write (*,*) ' compute their weighted averages.'
        write (*,*) ' '
        write (*,*) ' Specify the output files for each components as:'
        write (*,*) '   X   (+ is right):        -avg_vx <file>'
        write (*,*) '   Y   (+ is up):           -avg_vy <file>'
        write (*,*) '   LOS (+ is towards obs):  -avg_vlos <file>'
        write (*,*)
        write (*,*)
        write (*,*) 'LOS ANGLE:'
        write (*,*)
        write (*,*) ' Use -avg_los_angle <file> to output the'//
     &              ' weighted average LOS angle.'
        write (*,*)
        write (*,*) ' Here the "LOS angle" is the angle between a'//
     &              ' position on the LOS and'
        write (*,*) ' the closest approach of the LOS (rmin).'
        write (*,*)
        write (*,*) ' A positive LOS angle is in front of rmin,'//
     &              ' negative is behind.'
        write (*,*)
        write (*,*) ' For plane-parallel integration, rmin lies'//
     &              ' in the plane of sky.'
        write (*,*)
        write (*,*)
        write (*,*) 'EMISSIVITY WEIGHTING:'
        write (*,*)
        write (*,*) ' GETPB can weight the LOS averages'//
     &              ' by b OR pb (not both at once).'
        write (*,*)
        write (*,*) ' The default is to weight by the pb emissivity.'
        write (*,*)
        write (*,*) ' Set the flag -avg_using_b to weight by b'
        write (*,*)
c
      call exit (1)
c
      end



