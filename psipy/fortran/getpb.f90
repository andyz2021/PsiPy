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

!#######################################################################

      subroutine GETPB(py_rho, py_b, py_pb, py_help, py_verbose, py_cubic, py_oldmas, py_long, &
       py_p, py_b0, py_r, py_nx, py_ny, py_x0, py_x1, py_y0, py_y1, py_wispr1, py_wispr2,&
       py_rocc, py_dsmult, py_power, py_disk, py_vf, py_vr, py_vt, py_vp, py_scalar,&
       py_avg_scalar, py_avg_los_angle, py_avg_vlos, py_avg_vx, py_avg_vy, py_avg_using_b, py_he_frac, py_mu, py_write_to_file, ans)
!
!-----------------------------------------------------------------------
!
      use ident
      use number_types
      use constants
      use params
      use image_region
      use geometry
      use mas_fields
      use image_fields
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      real(r_typ), dimension(:), pointer :: x,y
      REAL, intent(out), dimension(py_nx, py_ny) :: ans

      character(512) :: py_rho
!
      character(512) :: py_vr, py_vt, py_vp
      character(512) :: py_scalar
!
      character(512) :: py_avg_vlos, py_avg_vx, py_avg_vy
      character(512) :: py_avg_scalar
      character(512) :: py_avg_los_angle
!
      character(512) :: py_pb
      character(512) :: py_b
!
      character(512) :: py_vf
!
!
! Options for emissivity weighted integrals.
!
!      logical :: do_integral_weighting=.false.
!      logical :: compute_vector_integration=.false.
!      logical :: compute_scalar_integration=.false.
!      logical :: compute_angle_integration=.false.
!!
!!Weight the scalar/vector integrals by B instead (default pB).
!!
!      logical :: weight_integral_by_b=.false.
!
! ****** Flag indicating that the density file is from the
! ****** old MAS code.
!
      logical :: py_oldmas

      logical :: py_avg_using_b

      logical :: py_verbose

      logical :: py_help

      logical :: py_write_to_file
!
      real :: py_long
      real :: py_b0
      real :: py_p
      real :: py_rocc
      real :: py_power
      real :: py_disk

      integer :: py_y0
      integer :: py_y1
      integer :: py_x0
      integer :: py_x1
      integer :: py_nx
      integer :: py_ny
!
!
! ****** Distance of the observer from the Sun.
!
      real :: py_r
!
! ****** Step size multiplier.
!
      real :: py_dsmult
!
! ****** Shorthand flags to set the WISPR camera domains.
!
      logical :: py_wispr1
      logical :: py_wispr2
!
! ****** Helium Fraction Parameters (like in MAS).
! ****** HE_RHO is defined by rho(norm)=HE_RHO*n_e(norm),
!
      real :: py_he_frac
!
! ****** Option to use cubic interpolation
!
      logical :: py_cubic
!
! ****** Limb darkening coefficient at 520 nm. Default is from
! ****** Altschuler & Perry, Solar Physics, 23, 410 (1972).
!
      real :: py_mu
!
!-----------------------------------------------------------------------
!
! ****** Set the parameters.
!

      call python_getpb(py_rho, py_b, py_pb, py_help, py_verbose, py_cubic, py_oldmas, py_long, &
       py_p, py_b0, py_r, py_nx, py_ny, py_x0, py_x1, py_y0, py_y1, py_wispr1, py_wispr2,&
       py_rocc, py_dsmult, py_power, py_disk, py_vf, py_vr, py_vt, py_vp, py_scalar,&
       py_avg_scalar, py_avg_los_angle, py_avg_vlos, py_avg_vx, py_avg_vy, py_avg_using_b, py_he_frac, py_mu, py_write_to_file)
!
! ****** Check that at least one of B or pB was requested.
!
!      print *, 'rho=', rho_file
!      print *, 'len=', len_trim(rho_file)
      if ((len_trim(rho_file) == 0)) then
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
!
! ****** Check that the image dimensions are valid.
!
      if (nx.lt.1.or.ny.lt.1) then
        write (*,*)
        write (*,*) '### ERROR in ',cname,':'
        write (*,*) '### The specified number of points for the'//&
     &              ' image are invalid:'
        write (*,*) 'NX = ',nx
        write (*,*) 'NY = ',ny
        call exit (1)
      end if
!
! ****** Check that the WISPR camera shorthand flags are
! ****** consistent with the usage of other options.
!
      if (wispr1.and.wispr2) then
        write (*,*)
        write (*,*) '### ERROR in ',cname,':'
        write (*,*) '### You cannot specify both the'//&
     &              ' -wispr1 and -wispr2 flags.'
        call exit (1)
      end if
!
      if ((wispr1.or.wispr2).and..not.r_obs_set) then
        write (*,*)
        write (*,*) '### ERROR in ',cname,':'
        write (*,*) '### You must specify the Sun-observer'//&
     &              ' distance when using the -wispr1'
        write (*,*) '### or -wispr2 flags.'
        call exit (1)
      end if
!
      if ((wispr1.or.wispr2).and.image_limits_set) then
        write (*,*)
        write (*,*) '### ERROR in ',cname,':'
        write (*,*) '### You cannot set the image coordinate'//&
     &              ' limits when using the -wispr1 or'
        write (*,*) '### -wispr2 flags.'
        call exit (1)
      end if
!
! ****** Check that the observer distance is positive (if
! ****** it was specified).
!
      if (r_obs_set) then
        if (r_obs.le.0.) then
          write (*,*)
          write (*,*) '### ERROR in ',cname,':'
          write (*,*) '### The distance of the observer'//&
     &                ' from the Sun must be positive:'
          write (*,*) 'Specified distance [AU] = ',r_obs
          call exit (1)
        end if
      end if
!
! ****** Set the image region for WISPR, if requested.
!
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
!
! ****** Check the specified image plane coordinates.
!
      if (x0.gt.x1.or.y0.gt.y1) then
        write (*,*)
        write (*,*) '### ERROR in ',cname,':'
        write (*,*) '### The specified image plane coordinates'//&
     &              ' are invalid:'
        write (*,*) 'X0 = ',x0
        write (*,*) 'X1 = ',x1
        write (*,*) 'Y0 = ',y0
        write (*,*) 'Y1 = ',y1
        call exit (1)
      end if
!
! ****** Read the density.
!
      call read_field( rho_file, 'density', rho)
!
! ****** Read the scalar field (if necessary).
!
      if (compute_scalar_integration) then
        call read_field( scalar_file, 'scalar field', sfield)
      endif
!
! ****** Read the vector fields (if necessary).
!
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
!
! ****** Read the vignetting function if it was requested.
!
      if (vf_file.ne.' ') then
        call read_vignetting_function (vf_file)
      end if
!
! ****** Convert the observer distance to solar radii.
!
      r_obs_rs=r_obs*au_in_rs
!
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
!
! ****** Convert the Carrington longitude to radians.
!
      cml=crlong*degtorad
!
! ****** Convert the solar P and B0 angles to radians.
!
      pangle=pangle*degtorad
      b0angle=b0angle*degtorad
!
      sin_p=sin(pangle)
      cos_p=cos(pangle)
      sin_b0=sin(b0angle)
      cos_b0=cos(b0angle)
!
      if (verbose) then
        write (*,*)
        write (*,*) '### Image region:'
        if (r_obs_set) then
          write (*,901) 'Elongation [deg]:',&
     &                  ' X0 = ',x0,&
     &                  ' X1 = ',x1
          write (*,901) 'Altitude   [deg]:',&
     &                  ' Y0 = ',y0,&
     &                  ' Y1 = ',y1
        else
          write (*,901) 'POS x limits [R_sun]:',&
     &                  ' X0 = ',x0,&
     &                  ' X1 = ',x1
          write (*,901) 'POS y limits [R_sun]:',&
     &                  ' Y0 = ',y0,&
     &                  ' Y1 = ',y1
        end if
      end if
  901 format (1x,a,a,f15.6,',',a,f15.6)
!
! ****** Allocate the arrays for the results, initialized with zeros.
!
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
!
! ****** Print the helium fraction if it was specified.
!
      if (verbose.and.he_frac.ne.0.) then
        write (*,*)
        write (*,*) 'Using a non-zero helium fraction!'
        write (*,'(A,F10.7)') '  he_frac: ', he_frac
        write (*,'(A,F10.7)') '  he_rho:  ', he_rho
      end if
!
! ****** Print the limb-darkening coefficient.
!
      if (verbose.and.mu_ld.ne.0.63) then
        write (*,*)
        write (*,*) 'Using a non-standard limb-darkening coefficient!'
        write (*,'(A,F10.7)') ' mu: ', mu_ld
      endif
!
! ****** Error check the integral weighting option.
!
      if (compute_vector_integration.or.compute_scalar_integration.or.&
     &    compute_angle_integration) then
!
        do_integral_weighting=.true.
!
        if (weight_integral_by_b) then
          if (.not.compute_b) then
            write (*,*)
            write (*,*) '### ERROR in ',cname,':'
            write (*,*) '### Integral weighting specified with B but'//&
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
            write (*,*) '### Integral weighting specified with pB but'//&
     &                  ' -pb was not specified!'
            call exit (1)
          endif
          if (verbose) then
            write (*,*)
            write (*,*) 'Integral averages requested! Weighting: pB.'
          endif
        endif
      endif
!
! ****** Compute pB and/or B.
!
      if (verbose) then
        write (*,*)
        write (*,*) 'Calculating pB and/or B ...'
      end if
!
      if (r_obs_set) then
        call get_image_obs_distance_specified
        x=>elong
        y=>alt
      else
        call get_image_obs_at_infinity
        x=>x_pos
        y=>y_pos
      end if
!
! ****** Normalize wieghted integrals if necessary.
!
      if (do_integral_weighting) then
        call normalize_integrals
      endif
!
! ****** Multiply pB and B by the requested radial power.
!
      if (power_set) then
        if (compute_pb) then
          if (verbose) then
            write (*,*)
            write (*,*) '### Multiplying pB by r_min^p with p = ',&
     &                  power
          end if
          call multiply_by_radial_power (pb)
        end if
        if (compute_b) then
          if (verbose) then
            write (*,*)
            write (*,*) '### Multiplying B by r_min^p with p = ',&
     &                  power
          end if
          call multiply_by_radial_power (b)
        end if
      end if
!
! ****** Apply the vignetting function (if it was specified).
!
      if (vf_file.ne.' ') then
        if (compute_pb) then
          if (verbose) then
            write (*,*)
            write (*,*) '### Multiplying pB by the vignetting'//&
     &                  ' function ... '
          end if
          call vignette (pb)
        end if
        if (compute_b) then
          if (verbose) then
            write (*,*)
            write (*,*) '### Multiplying B by the vignetting'//&
     &                  ' function ... '
          end if
          call vignette (b)
        end if
      end if
!
! ****** Set the precision of the output files.
!
      hdf32 = rho%sds%hdf32
!
      if (verbose .AND. write_to_file) then
            write (*,*)
            write (*,*) 'Writing to pb/b file.'
          endif
! ****** Write the pB image.
!
      if (compute_pb .AND. write_to_file) call write_image(pb_file, 'pB', x, y, pb)
      if (compute_pb) ans = pb
!call write_image(pb_file, 'pB', x, y, pb)
! ****** We want it to return pb/b for the python version

!
! ****** Write the B image.
!
       if (compute_b .AND. write_to_file) call write_image(b_file, 'B', x, y, b)
       if (compute_b) ans = b
!
! ****** Write the scalar field image.
!

      if (compute_scalar_integration)&
     &    call write_image(scalar_file_out, 'LOS averaged scalar field',&
     &                     x, y, sf_avg)
!
! ****** Write the LOS contribution angle image.
!
      if (compute_angle_integration)&
     &    call write_image(angle_file_out,&
     &           'LOS averaged contribution angle', x, y, angle_avg)
!
! ****** Write the vector field images.
!
      if (compute_vector_integration) then
        call write_image(vlos_file_out, 'LOS averaged vlos',&
     &                   x,y,vlos_avg)
        call write_image(vx_file_out, 'LOS averaged vx',x,y,vx_avg)
        call write_image(vy_file_out, 'LOS averaged vy',x,y,vy_avg)
      endif
!
!      print *, 'hi = ', ans
      end

!#######################################################################
      subroutine python_getpb(py_rho, py_b, py_pb, py_help, py_verbose, py_cubic, py_oldmas, py_long, &
       py_p, py_b0, py_r, py_nx, py_ny, py_x0, py_x1, py_y0, py_y1, py_wispr1, py_wispr2,&
       py_rocc, py_dsmult, py_power, py_disk, py_vf, py_vr, py_vt, py_vp, py_scalar,&
       py_avg_scalar, py_avg_los_angle, py_avg_vlos, py_avg_vx, py_avg_vy, py_avg_using_b, py_he_frac, py_mu, py_write_to_file)

!
!-----------------------------------------------------------------------
!
! ****** Set parameters from the command line arguments.
!
!-----------------------------------------------------------------------
!
      use ident
      use number_types
      use syntax
      use paragraph_def
      use get_usage_line_interface
      use print_par_interface
      use delete_par_interface
      use params
      use image_region
!
!-----------------------------------------------------------------------
!
      implicit none

      character(512) :: py_rho
!
      character(512) :: py_vr, py_vt, py_vp
      character(512) :: py_scalar
!
      character(512) :: py_avg_vlos, py_avg_vx, py_avg_vy
      character(512) :: py_avg_scalar
      character(512) :: py_avg_los_angle
!
      character(512) :: py_pb
      character(512) :: py_b
!
      character(512) :: py_vf
!
!
! Options for emissivity weighted integrals.
!
!      logical :: do_integral_weighting=.false.
!      logical :: compute_vector_integration=.false.
!      logical :: compute_scalar_integration=.false.
!      logical :: compute_angle_integration=.false.
!!
!!Weight the scalar/vector integrals by B instead (default pB).
!!
!      logical :: weight_integral_by_b=.false.
!
! ****** Flag indicating that the density file is from the
! ****** old MAS code.
!
      logical :: py_oldmas

      logical :: py_avg_using_b

      logical :: py_verbose

      logical :: py_help

      logical :: py_write_to_file
!
      real :: py_long
      real :: py_b0
      real :: py_p
      real :: py_rocc
      real :: py_power
      real :: py_disk

      integer :: py_y0
      integer :: py_y1
      integer :: py_x0
      integer :: py_x1
      integer :: py_nx
      integer :: py_ny
!
!
! ****** Distance of the observer from the Sun.
!
      real :: py_r
!
! ****** Step size multiplier.
!
      real :: py_dsmult
!
! ****** Shorthand flags to set the WISPR camera domains.
!
      logical :: py_wispr1
      logical :: py_wispr2
!
! ****** Helium Fraction Parameters (like in MAS).
! ****** HE_RHO is defined by rho(norm)=HE_RHO*n_e(norm),
!
      real :: py_he_frac
!
! ****** Option to use cubic interpolation
!
      logical :: py_cubic
!
! ****** Limb darkening coefficient at 520 nm. Default is from
! ****** Altschuler & Perry, Solar Physics, 23, 410 (1972).
!
      real :: py_mu
!
!-----------------------------------------------------------------------
!
! ****** Storage the for usage line.
!
!
! ****** Storage for the error message.
!
!
!-----------------------------------------------------------------------
!
      logical, save :: print_help=.false.

!
!-----------------------------------------------------------------------
!
      integer, external :: intval
      real, external :: fpval
!
!-----------------------------------------------------------------------
!
!
!
! ****** Parse the command line.
!
!
! ****** Check if the "-help" option was specified.
!
!       call nargs_specified (nargs)
!       call fetcharg ('-help',set,arg)
!
        if (py_help) then
          print_help=.true.
        else
          print_help=.false.
        end if

! ****** Set the parameters.
!
! ****** Verbose flag.
!
      if (py_verbose) then
        verbose=.true.
      else
        verbose=.false.
      end if
!
! ****** Cubic interpolation flag.
!
      if (py_cubic) then
        cubic=.true.
      else
        cubic=.false.
      end if
!
! ****** Switch to indicate old MAS code data files.
!
      if (py_oldmas) then
        oldmas=.true.
      else
        oldmas=.false.
      end if

! ****** Whether to write to file or not

      if (py_write_to_file) then
          write_to_file=.true.
      else
          write_to_file=.false.
      end if


!
! ****** Carrington rotation longitude.
!

      crlong=py_long
!
! ****** Solar P angle.
!
      pangle=py_p
!
! ****** Solar B0 angle.
!
      b0angle=py_b0
!
! ****** Distance of the observer from the Sun [AU].
!
!      print *, 'py_r = ', py_r
      if (py_r == 0) then
        r_obs_set=.false.
        r_obs=1._r_typ
      else
         r_obs_set=.true.
        r_obs=py_r
      end if
!
! ****** Number of points to use for the image.
!
      nx=py_nx
!
      ny=py_ny
!
! ****** Image dimensions.
!
      image_limits_set=.false.
!
      if (py_x0 /= -3) image_limits_set=.true.
      x0=py_x0
!
      if (py_x1 /= 3) image_limits_set=.true.
      x1=py_x1
!
      if (py_y0 /= -3) image_limits_set=.true.
      y0=py_y0
!
      if (py_y1 /= 3) image_limits_set=.true.
      y1=py_y1
!
! ****** WISPR camera 1 flag.
!
      if (py_wispr1) then
        wispr1=.true.
      else
         wispr1=.false.
      end if

!
! ****** WISPR camera 2 flag.
!
      if (py_wispr2) then
        wispr2=.true.
      else
         wispr2=.false.
      end if
!
! ****** Occulting disk radius.
!
      rocc=py_rocc
!
! ****** Factor by which to multiply the LOS integration step size.
!
      dsmult=py_dsmult
!
! ****** Radial power exponent.
!
      if (py_power /= 0) then
        power_set=.true.
        power=py_power
      else
        power_set=.false.
        power=0
      endif
!
! ****** Value for pB on the solar disk.
!
      disk_value=py_disk
!
! ****** Vignetting function file name.
!
      vf_file=trim(py_vf)
      if ((len_trim(vf_file) == 0)) then
        vf_file=' '
      end if
!
! ****** Density file name.
!
      rho_file=trim(py_rho)
!
! ****** Vector field r component file name (input).
!
      vr_file=trim(py_vr)
      if (len_trim(vr_file) == 0) then
        vr_file=' '
      else
        compute_vector_integration=.true.
        vr_file=trim(py_vr)

      end if
!
! ****** Vector field t component file name (input).
!
      vt_file=trim(py_vt)
      if (len_trim(vt_file) == 0) then
        vt_file=' '
      else
        compute_vector_integration=.true.
        vt_file=trim(py_vt)

      end if
!
! ****** Vector field p component file name (input).
!
      vp_file=trim(py_vp)
      if (len_trim(vp_file) == 0) then
        vp_file=' '
      else
        compute_vector_integration=.true.
        vp_file=trim(py_vp)

      end if
!
! ****** Scalar field file name (input).
!
      scalar_file=trim(py_scalar)
      if (len_trim(scalar_file) == 0) then
        compute_scalar_integration=.false.
        scalar_file=' '
      else
        compute_scalar_integration=.true.
        scalar_file=trim(py_scalar)
      end if
!
! ****** Emissivity weighted LOS integral of scalar field.
!
      scalar_file_out=trim(py_avg_scalar)
      if (len_trim(scalar_file_out) == 0) then
        scalar_file_out=' '
      else
        scalar_file_out=trim(py_avg_scalar)
      end if
!
! ****** Emissivity weighted LOS contribution angle.
!
      angle_file_out=trim(py_avg_los_angle)
      if (len_trim(angle_file_out) == 0) then
        angle_file_out=' '
      else
        compute_angle_integration=.true.
        angle_file_out=trim(py_avg_los_angle)
      end if
!
! ****** Emissivity weighted LOS integral of vector field (vlos).
!
      vlos_file_out=trim(py_avg_vlos)
      if (len_trim(vlos_file_out) == 0) then
        vlos_file_out=' '
      else
        vlos_file_out=trim(py_avg_vlos)
      end if
!
! ****** Emissivity weighted LOS integral of vector field (vx).
!
      vx_file_out=trim(py_avg_vx)
      if (len_trim(vx_file_out) == 0) then
        vx_file_out=' '
      else
        vx_file_out=trim(py_avg_vx)
      end if

!
! ****** Emissivity weighted LOS integral of vector field (vy).
!
      vy_file_out=trim(py_avg_vy)
      if (len_trim(vy_file_out) == 0) then
        vy_file_out=' '
      else
        vy_file_out=trim(py_avg_vy)
      end if

!
! ****** Flag to weight integral averages with B and not pB.
!
      if (py_avg_using_b) then
        weight_integral_by_b=.true.
      else
        weight_integral_by_b=.false.
      end if
!
! ****** Helium Fraction
!
      he_frac=py_he_frac
      he_rho=(1._r_typ+4._r_typ*he_frac)/(1._r_typ+2._r_typ*he_frac)
!
! ****** Value for limb-darkening coefficient.
!
      mu_ld=py_mu
!
! ****** Output pB file name.
!
      if (len_trim(py_pb) == 0) then
        compute_pb=.false.
        pb_file=' '
      else
        compute_pb=.true.
        pb_file=trim(py_pb)
      end if

!
! ****** Output B file name.
!
      if (len_trim(py_b) == 0) then
        compute_b=.false.
        b_file=' '
      else
        compute_b=.true.
        b_file=trim(py_b)
      end if
!
      end
!#######################################################################
!      subroutine print_help_text
!!
!!-----------------------------------------------------------------------
!!
!! ****** Print the help / syntax information.
!!
!!-----------------------------------------------------------------------
!!
!      implicit none
!!
!!-----------------------------------------------------------------------
!!
!        write (*,*)
!        write (*,*)
!        write (*,*) '### Help for GETPB Special Integration Methods:'
!        write (*,*)
!        write (*,*) 'GETPB can now compute emissivity weighted'//&
!     &              ' averages of other quantities.'
!        write (*,*) 'Here are the additions flags/options: '
!        write (*,*)
!        write (*,*)
!        write (*,*) 'SCALAR FIELD AVERAGING:'
!        write (*,*)
!        write (*,*) ' Use -scalar <file> to specify a 2D or 3D file'//&
!     &              ' with a scalar field.'
!        write (*,*)
!        write (*,*) ' GETPB will integrate the scalar along the LOS'//&
!     &              ' to get the weighted average.'
!        write (*,*)
!        write (*,*) ' Specify the output file with -avg_scalar <file>'
!        write (*,*)
!        write (*,*)
!        write (*,*) 'VECTOR FIELD COMPONENT AVERAGING:'
!        write (*,*)
!        write (*,*) ' Use -vr <file>, -vt <file>, -vp <file> to'//&
!     &              ' specify a 2D or 3D vector field.'
!        write (*,*)
!        write (*,*) ' GETPB will project the rtp vector field into'//&
!     &              ' LOS image coordinates and'
!        write (*,*) ' compute their weighted averages.'
!        write (*,*) ' '
!        write (*,*) ' Specify the output files for each components as:'
!        write (*,*) '   X   (+ is right):        -avg_vx <file>'
!        write (*,*) '   Y   (+ is up):           -avg_vy <file>'
!        write (*,*) '   LOS (+ is towards obs):  -avg_vlos <file>'
!        write (*,*)
!        write (*,*)
!        write (*,*) 'LOS ANGLE:'
!        write (*,*)
!        write (*,*) ' Use -avg_los_angle <file> to output the'//&
!     &              ' weighted average LOS angle.'
!        write (*,*)
!        write (*,*) ' Here the "LOS angle" is the angle between a'//&
!     &              ' position on the LOS and'
!        write (*,*) ' the closest approach of the LOS (rmin).'
!        write (*,*)
!        write (*,*) ' A positive LOS angle is in front of rmin,'//&
!     &              ' negative is behind.'
!        write (*,*)
!        write (*,*) ' For plane-parallel integration, rmin lies'//&
!     &              ' in the plane of sky.'
!        write (*,*)
!        write (*,*)
!        write (*,*) 'EMISSIVITY WEIGHTING:'
!        write (*,*)
!        write (*,*) ' GETPB can weight the LOS averages'//&
!     &              ' by b OR pb (not both at once).'
!        write (*,*)
!        write (*,*) ' The default is to weight by the pb emissivity.'
!        write (*,*)
!        write (*,*) ' Set the flag -avg_using_b to weight by b'
!        write (*,*)
!!
!      call exit (1)
!!
!      end



