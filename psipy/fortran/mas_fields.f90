      module ident
!
      character(*), parameter :: cname='GETPB'
      character(*), parameter :: cvers='1.13'
      character(*), parameter :: cdate='04/22/2021'
!
      end module
!#######################################################################
      module types
!
!-----------------------------------------------------------------------
! ****** Definition of data structures.
!-----------------------------------------------------------------------
!
      use number_types
      use sds_def
      use invint_def
      use spline_def
!
      implicit none
!
! ****** Inverse interpolation table structure definitions.
!
      type :: vtab
        type(itab), dimension(3) :: c
      end type
!
! ****** Container for 2D/3D MAS style fields and interpolation tables.
!
      !
      type :: mas_field
        type(sds) :: sds
        type(vtab) :: invtab
        type(spl2d) :: spl2
        type(spl3d) :: spl3
      end type
!
      end module
!#######################################################################
      module constants
!
      use number_types
!
      implicit none
!
      real(r_typ), parameter :: pi=3.14159265358979323846_r_typ
      real(r_typ), parameter :: twopi=pi*2._r_typ
      real(r_typ), parameter :: halfpi=pi*.5_r_typ
      real(r_typ), parameter :: degtorad=pi/180._r_typ
!
! ****** One AU in solar radii.
!
      real(r_typ), parameter :: au_in_rs=214.939469_r_typ
!
! ****** Offset to ensure integration is inside domain.
!
      real(r_typ), parameter :: integration_offset=1.e-10_r_typ
!
      end module
!#######################################################################
      module params
!
!-----------------------------------------------------------------------
! ****** Parameters.
!-----------------------------------------------------------------------
!
      use number_types
!
      implicit none
!
      character(512) :: rho_file
!
      character(512) :: vr_file, vt_file, vp_file
      character(512) :: scalar_file
!
      character(512) :: vlos_file_out, vx_file_out, vy_file_out
      character(512) :: scalar_file_out
      character(512) :: angle_file_out
!
      character(512) :: pb_file
      character(512) :: b_file
!
      character(512) :: vf_file
!
      logical :: verbose
!
      logical :: compute_pb
      logical :: compute_b

      logical :: write_to_file
!
! ****** Options for emissivity weighted integrals.
!
      logical :: do_integral_weighting=.false.
      logical :: compute_vector_integration=.false.
      logical :: compute_scalar_integration=.false.
      logical :: compute_angle_integration=.false.
!
! ****** Weight the scalar/vector integrals by B instead (default pB).
!
      logical :: weight_integral_by_b=.false.
!
! ****** Flag indicating that the density file is from the
! ****** old MAS code.
!
      logical :: oldmas
!
      real(r_typ) :: crlong
      real(r_typ) :: b0angle
      real(r_typ) :: pangle
      real(r_typ) :: rocc
      real(r_typ) :: power
      real(r_typ) :: disk_value
!
      logical :: power_set
!
! ****** Distance of the observer from the Sun.
!
      logical :: r_obs_set
      real(r_typ) :: r_obs
      real(r_typ) :: r_obs_rs
!
! ****** Step size multiplier.
!
      real(r_typ) :: dsmult
!
! ****** Shorthand flags to set the WISPR camera domains.
!
      logical :: wispr1
      logical :: wispr2
!
! ****** Helium Fraction Parameters (like in MAS).
! ****** HE_RHO is defined by rho(norm)=HE_RHO*n_e(norm),
!
      real(r_typ) :: he_frac = 0.0
      real(r_typ) :: he_rho
!
! ****** Option to use cubic interpolation
!
      logical :: cubic
!
! ****** Limb darkening coefficient at 520 nm. Default is from
! ****** Altschuler & Perry, Solar Physics, 23, 410 (1972).
!
      real(r_typ) :: mu_ld = 0.63_r_typ
!
! ****** Flag for the precision of the output files.
!
      logical :: hdf32
!
      end module
!#######################################################################
      module image_region
!
      use number_types
!
      implicit none
!
! ****** Image dimensions.
!
      integer :: nx
      integer :: ny
!
! ****** Image region limits.
!
      real(r_typ) :: x0
      real(r_typ) :: x1
      real(r_typ) :: y0
      real(r_typ) :: y1
!
      logical :: image_limits_set
!
! ****** Number of OpenMP iterations to do in each thread.
!
      integer, parameter :: iterations_per_thread=500
!
      end module
!#######################################################################
      module geometry
!
      use number_types
!
      implicit none
!
! ****** Sin and cos of the P and B0 angles.
!
      real(r_typ) :: sin_p,cos_p
      real(r_typ) :: sin_b0,cos_b0
!
! ****** Central meridian longitude.
!
      real(r_typ) :: cml
!
! ****** Flag to indicate an axisymmetric density (i.e., 2D).
!
      logical :: twodee
!
! ****** Maximum radius for all of the input meshes.
!
      real(r_typ) :: rmax=-1.0
!
      end module
!#######################################################################
      module mas_fields
!
!-----------------------------------------------------------------------
! ****** Storage for the density.
!-----------------------------------------------------------------------
!
      use number_types
      use types
      use spline_def
!
      implicit none
!
! ****** Density field container.
!
      type(mas_field) :: rho
!
! ****** Containers for arbitrary vector field, v, in spherical coords.
!
      type(mas_field) :: vr
      type(mas_field) :: vt
      type(mas_field) :: vp
!
! ****** Containers for arbitrary scalar field.
!
      type(mas_field) :: sfield
!
! ****** Normalization factor for electron density.
!
! ****** Multiply the normalized density read in [in MAS code units]
! ****** by FN_N to get electron density in [/cm^3].
!
      real(r_typ), parameter :: fn_n=1.e8_r_typ
!
      end module
!#######################################################################
      module image_fields
!
!-----------------------------------------------------------------------
! ****** Storage for the computed pB and B, and their associated
! ****** arrays.
!-----------------------------------------------------------------------
!
      use number_types
      use types
!
      implicit none
!
! ****** Image axes for the case when the observer distance
! ****** has been specified (elongation, altitude).
!
      real(r_typ), dimension(:), pointer :: elong,alt
!
! ****** Image axes for the case when the observer distance
! ****** is infinite (POS x and y).
!
      real(r_typ), dimension(:), pointer :: x_pos,y_pos
!
! ****** Polarization brightness and brightness.
!
      real(r_typ), dimension(:,:), allocatable :: pb
      real(r_typ), dimension(:,:), allocatable :: b
!
! ****** Emissivity weighted averages of the vector field components
! ****** in the frame of the image plane (los is TOWARDS observer).
!
      real(r_typ), dimension(:,:), allocatable :: vlos_avg
      real(r_typ), dimension(:,:), allocatable :: vx_avg
      real(r_typ), dimension(:,:), allocatable :: vy_avg
!
! ****** Emissivity weighted averages of the scalar field.
!
      real(r_typ), dimension(:,:), allocatable :: sf_avg
!
! ****** Emissivity weighted average elongation angle (0 when r=rmin).
!
      real(r_typ), dimension(:,:), allocatable :: angle_avg
!
! ****** Closest distance to the Sun.  This is used for vignetting
! ****** and multiplication by a radial power.
!
      real(r_typ), dimension(:,:), allocatable :: r_min
!
      end module
!#######################################################################
      module vignetting_function
!
!-----------------------------------------------------------------------
! ****** Vignetting function definition.
!-----------------------------------------------------------------------
!
      use number_types
      use spline_def
!
      implicit none
!
! ****** Vignetting function.
!
      integer :: nvf
      real(r_typ), dimension(:), allocatable :: r_vf
      real(r_typ), dimension(:), allocatable :: vf
!
! ****** Spline structure for cubic interpolation.
!
      type(spl1d) :: vf_spl1
!
      end module

           subroutine get_image_obs_distance_specified
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
      use ident
      use number_types
      use constants
      use params
      use image_region
      use geometry
      use image_fields
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      real(r_typ), parameter :: half=.5_r_typ
      real(r_typ), parameter :: one=1._r_typ
!
!-----------------------------------------------------------------------
!
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
!
!-----------------------------------------------------------------------
!
      real(r_typ) :: xo,yo,zo,xr,yr,zr,xv,yv,zv
      real(r_typ) :: ds_fac, x_los
!
! ****** Unit vectors of observer frame in heliographic coordinates.
!
      real(r_typ), dimension(3) :: ex, ey, ez
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
      allocate (elong(nx))
      allocate (alt(ny))
!
      if (nx.eq.1) then
        elong(1)=x0
      else
        dx=(x1-x0)/(nx-1)
        do i=1,nx
          elong(i)=x0+(i-1)*dx
        enddo
      end if
!
      if (ny.eq.1) then
        alt(1)=y0
      else
        dy=(y1-y0)/(ny-1)
        do j=1,ny
          alt(j)=y0+(j-1)*dy
        enddo
      end if
!
! ****** Cartesian vector from the Sun to the observer.
!
      v_obs(1)=r_obs_rs
      v_obs(2)=0.
      v_obs(3)=0.
!
! ****** Calculate pB and/or B.
!
      do j=1,ny
        do i=1,nx
!
! ****** Cartesian vector from the Sun to the reference point.
!
! ****** This is a point on a sphere of radius R_OBS centered
! ****** at the observer, where R_OBS is the distance of the
! ****** observer from the Sun.
!
          e_rad=degtorad*elong(i)
          a_rad=degtorad*alt(j)
          cos_a=cos(a_rad)
          v_ref(1)=r_obs_rs*(one+cos(pi-e_rad)*cos_a)
          v_ref(2)=r_obs_rs*sin(pi-e_rad)*cos_a
          v_ref(3)=r_obs_rs*sin(a_rad)
!
! ****** Find the point along this LOS that is closest to the
! ****** Sun.  This is the point on the LOS at which it
! ****** intersects the Thomson sphere.  The distance from this
! ****** point to the Sun is needed in the calculation of pB
! ****** and B.
!
          t_r_min=-( v_obs(1)*(v_ref(1)-v_obs(1))&
     &              +v_obs(2)*(v_ref(2)-v_obs(2))&
     &              +v_obs(3)*(v_ref(3)-v_obs(3))&
     &             )/r_obs_rs**2
!
          v_r_min(1)=t_r_min*v_ref(1)+(one-t_r_min)*v_obs(1)
          v_r_min(2)=t_r_min*v_ref(2)+(one-t_r_min)*v_obs(2)
          v_r_min(3)=t_r_min*v_ref(3)+(one-t_r_min)*v_obs(3)
!
          r_min(i,j)=sqrt( v_r_min(1)**2&
     &                    +v_r_min(2)**2&
     &                    +v_r_min(3)**2)
!
! ****** Check if this LOS intersects the occulting disk.
! ****** If so, skip the calculation along this LOS.
!
          if (r_min(i,j).le.max(rocc,one)) then
            if (compute_pb) pb(i,j)=disk_value
            if (compute_b) b(i,j)=disk_value
            cycle
          end if
!
! ****** The distance s is along the LOS, starting at the observer
! ****** (s = 0), and increasing toward the reference point
! ****** (s = R_OBS).  The parameter t is proportional to s,
! ****** with t = 0 specifying the observer, and t = 1 specifying
! ****** the reference point.
!
          if (compute_pb) pb(i,j)=0.
          if (compute_b) b(i,j)=0.
!
! ****** Check if the LOS doesn't intersect the grid, leave it zeroed.
!
          if (r_min(i,j).ge.rmax) cycle
!
! ****** Compute unit vectors for LOS if doing vector/angle projection.
! ****** Here the projection changes for each pixel (non-parallel).
!
         if (compute_vector_integration.or.&
     &       compute_angle_integration) then
           call transform_position (v_obs(1),v_obs(2),v_obs(3),rv,tv,pv)
           call rtp2xyz(rv,tv,pv,xo,yo,zo)
           call transform_position (v_r_min(1),v_r_min(2),v_r_min(3),&
     &                              rv,tv,pv)
           call rtp2xyz(rv,tv,pv,xr,yr,zr)
           xv=xo-xr
           yv=yo-yr
           zv=zo-zr
           call get_unit_vectors(xv,yv,zv,ex,ey,ez)
         endif
!
! ****** Integrate along the LOS, starting at first intersection
! ****** of the LOS with the density mesh.
!
          dmax=sqrt( rmax**2 - r_min(i,j)**2)
          s0 = t_r_min*r_obs_rs - dmax + integration_offset
          s1 = t_r_min*r_obs_rs + dmax - integration_offset
!
          s=s0
          dsm=0.
          nstop=0
          do
!
! ****** Set the position along the LOS.
!
            t=s/r_obs_rs
!
            v_los(1)=t*v_ref(1)+(one-t)*v_obs(1)
            v_los(2)=t*v_ref(2)+(one-t)*v_obs(2)
            v_los(3)=t*v_ref(3)+(one-t)*v_obs(3)
!
            call transform_position (v_los(1),v_los(2),v_los(3),&
     &                               rv,tv,pv)
!
            call get_kernels (r_min(i,j),rv,tv,pv,&
     &                        compute_pb,compute_b,&
     &                        pb_local,b_local,ds_local)
!
            ds=ds_local*dsmult
            s=s+ds
            if (s.ge.s1) then
              ds=ds-(s-s1)
              s=s1
              nstop=nstop+1
            end if
            ds_fac = half*(dsm+ds)
!
            if (compute_pb) then
              pb(i,j)=pb(i,j) + ds_fac*pb_local
            end if
!
            if (compute_b) then
              b(i,j)=b(i,j) + ds_fac*b_local
            end if
!
! ****** Compute emissivity weighted integrals if needed.
! ****** For non-plane-parallel angle need to dot position with x-hat.
!
            if (do_integral_weighting) then
              if (compute_angle_integration) then
                call rtp2xyz(rv,tv,pv,xv,yv,zv)
                x_los = ex(1)*xv + ex(2)*yv + ex(3)*zv
              else
                x_los = 0.
              endif
!
              call add_to_integrals(i,j,pb_local,b_local,ds_fac,&
     &               ex,ey,ez,rv,tv,pv,x_los)
            endif
!
            dsm=ds
            if (nstop.ge.2) exit
          enddo
        enddo
      enddo
!
      end
!#######################################################################
      subroutine get_image_obs_at_infinity
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
      use ident
      use number_types
      use constants
      use params
      use image_region
      use geometry
      use image_fields
      use mas_fields
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      real(r_typ), parameter :: zero=0._r_typ
      real(r_typ), parameter :: half=.5_r_typ
      real(r_typ), parameter :: one=1._r_typ
!
!-----------------------------------------------------------------------
!
      integer :: i,j,nstop
      real(r_typ) :: dx,dy
      real(r_typ) :: s0,s1
      real(r_typ) :: dmax
      real(r_typ) :: s,ds,dsm
      real(r_typ) :: ds_local
      real(r_typ) :: rv,tv,pv
      real(r_typ) :: pb_local,b_local
      real(r_typ), dimension(3) :: v_los
!
!-----------------------------------------------------------------------
!
      real(r_typ) :: xv,yv,zv
      real(r_typ) :: ds_fac
!
! ****** Unit vectors of observer frame in heliographic coordinates.
!
      real(r_typ), dimension(3) :: ex, ey, ez
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
      allocate (x_pos(nx))
      allocate (y_pos(ny))
!
      if (nx.eq.1) then
        x_pos(1)=x0
      else
        dx=(x1-x0)/(nx-1)
        do i=1,nx
          x_pos(i)=x0+(i-1)*dx
        enddo
      end if
!
      if (ny.eq.1) then
        y_pos(1)=y0
      else
        dy=(y1-y0)/(ny-1)
        do j=1,ny
          y_pos(j)=y0+(j-1)*dy
        enddo
      end if
!
! ****** Compute unit vectors for LOS if doing vector projection.
! ****** For the plane-parallel assumption, they never change.
!
      if (compute_vector_integration) then
        call transform_position (one,zero,zero,rv,tv,pv)
        call rtp2xyz(rv,tv,pv,xv,yv,zv)
        call get_unit_vectors(xv,yv,zv,ex,ey,ez)
      endif
!
! ****** Calculate pB and/or B.
!
      do j=1,ny
        do i=1,nx
!
! ****** Calculate the plane-of-sky radius.
!
          r_min(i,j)=sqrt(x_pos(i)**2+y_pos(j)**2)
!
! ****** Check if this LOS intersects the occulting disk.
! ****** If so, skip the calculation along this LOS.
!
          if (r_min(i,j).le.max(rocc,one)) then
            if (compute_pb) pb(i,j)=disk_value
            if (compute_b) b(i,j)=disk_value
            cycle
          end if
!
! ****** The distance s is along the LOS, referenced with respect
! ****** to the solar limb (s = 0), and increasing away from the
! ****** observer (so s = -x == v_los(1) in the observer frame).
!
          if (compute_pb) pb(i,j)=0.
          if (compute_b) b(i,j)=0.
!
! ****** Check if the LOS doesn't intersect the grid, leave it zeroed.
!
          if (r_min(i,j).ge.rmax) cycle
!
! ****** Integrate along the LOS, starting at first intersection
! ****** of the LOS with the density mesh.
!
          dmax=sqrt( rmax**2 - r_min(i,j)**2)
          s0 = -dmax + integration_offset
          s1 =  dmax - integration_offset
!
          s=s0
          dsm=0.
          nstop=0
          do
!
! ****** Set the position along the LOS.
!
            v_los(1)=-s
            v_los(2)=x_pos(i)
            v_los(3)=y_pos(j)
!
            call transform_position (v_los(1),v_los(2),v_los(3),&
     &                               rv,tv,pv)
!
            call get_kernels (r_min(i,j),rv,tv,pv,&
     &                        compute_pb,compute_b,&
     &                        pb_local,b_local,ds_local)
!
            ds=ds_local*dsmult
            s=s+ds
            if (s.ge.s1) then
              ds=ds-(s-s1)
              s=s1
              nstop=nstop+1
            end if
            ds_fac = half*(dsm+ds)
!
! ****** Compute pb and/or b contributions.
!
            if (compute_pb) then
              pb(i,j)=pb(i,j) + ds_fac*pb_local
            end if
            if (compute_b) then
              b(i,j)=b(i,j) + ds_fac*b_local
            end if
!
! ****** Compute emissivity weighted integrals if needed.
!
            if (do_integral_weighting) then
              call add_to_integrals(i,j,pb_local,b_local,ds_fac,&
     &               ex,ey,ez,rv,tv,pv,v_los(1))
!
            endif
!
            dsm=ds
            if (nstop.ge.2) exit
          enddo
!
        enddo
      enddo
!
      end
!#######################################################################
      subroutine read_field( filename, fieldname, field)
!
!-----------------------------------------------------------------------
!
! ****** Read a MAS style field.
!
!-----------------------------------------------------------------------
!
      use ident
      use params
      use geometry
      use types
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      character(512) :: filename
      character(*) :: fieldname
      type(mas_field) :: field
!
!-----------------------------------------------------------------------
!
      integer :: ierr
!
!-----------------------------------------------------------------------
!
      if (verbose) then
        write (*,*)
        write (*,*) '### Reading ',trim(fieldname),&
     &              ' from file: ',trim(filename)
      end if
!
! ****** Read the field from the specified file.
!
      call rdhdf (filename,field%sds,ierr)
!
      if (ierr.ne.0) then
        write (*,*)
        write (*,*) '### ERROR in ',cname,':'
        write (*,*) '### Could not read ',trim(fieldname),'.'
        write (*,*) 'IERR (from RDHDF) = ',ierr
        write (*,*) 'File name: ',trim(filename)
        call exit (1)
      end if
!
! ****** Check that it contains a 2D or 3D field.
!
      if (.not.(field%sds%ndim.eq.2.or.field%sds%ndim.eq.3)) then
        write (*,*)
        write (*,*) '### ERROR in ',cname,':'
        write (*,*) '### ',trim(fieldname),' file does'//&
     &              ' not contain a 2D or 3D field.'
        write (*,*) 'Number of dimensions = ',field%sds%ndim
        write (*,*) 'File name: ',trim(filename)
        call exit (1)
      end if
!
! ****** Check that it contains scales.
!
      if (.not.field%sds%scale) then
        write (*,*)
        write (*,*) '### ERROR in ',cname,':'
        write (*,*) '### The density file does'//&
     &              ' not contain scales.'
        write (*,*) 'File name: ',trim(filename)
        call exit (1)
      end if
!
! ****** Set the flag for an axisymmetric case.
!
      if (field%sds%ndim.eq.2) then
        twodee=.true.
      else
        twodee=.false.
      end if
!
! ****** Check that the field array has at least two points in
! ****** each (non-negligible) dimension.
!
      if (field%sds%dims(1).lt.2) then
        write (*,*)
        write (*,*) '### ERROR in ',cname,':'
        write (*,*) '### ',trim(fieldname),' array must have at least'//&
     &              ' 2 points in the r dimension.'
        write (*,*) 'Number of points in r = ',field%sds%dims(1)
        write (*,*) 'File name: ',trim(filename)
        call exit (1)
      end if
!
      if (field%sds%dims(2).lt.2) then
        write (*,*)
        write (*,*) '### ERROR in ',cname,':'
        write (*,*) '### ',trim(fieldname),' array must have at least'//&
     &              ' 2 points in the theta dimension.'
        write (*,*) 'Number of points in theta = ',field%sds%dims(2)
        write (*,*) 'File name: ',trim(filename)
        call exit (1)
      end if
!
      if (.not.twodee) then
        if (field%sds%dims(3).lt.2) then
          write (*,*)
          write (*,*) '### ERROR in ',cname,':'
          write (*,*) '### ',trim(fieldname),' array must have at'//&
     &                ' least 2 points in the phi dimension.'
          write (*,*) 'Number of points in phi = ',field%sds%dims(3)
          write (*,*) 'File name: ',trim(filename)
          call exit (1)
        end if
      end if
!
! ****** If we are using old MAS files, add a phi point.
!
      if (oldmas) then
        call add_phi_point (field%sds)
      end if
!
! ****** Construct the inverse interpolation table for
! ****** faster interpolation.
!
      call build_inverse_tables (field%sds,field%invtab)
!
! ****** Set the maximum radius available in the file but check
! ****** that a smaller value hasn't been set yet (rmax is initialized
! ****** as a negative number).
!
      if (rmax.lt.0.) then
        rmax=field%sds%scales(1)%f(field%sds%dims(1))
      else
        rmax=min(rmax,field%sds%scales(1)%f(field%sds%dims(1)))
      endif
!
! ****** Compute the spline coefficients if needed.
!
      if (cubic) then
        if (verbose) then
          write (*,*)
          write (*,*) '### Computing cubic spline coefficients for ',&
     &                 trim(fieldname)
        end if
        if (field%sds%ndim.eq.2) then
          call compute_spline_2d (field%sds%dims(1),&
     &                            field%sds%dims(2),&
     &                            field%sds%scales(1)%f,&
     &                            field%sds%scales(2)%f,&
     &                            field%sds%f,field%spl2)
        else if (field%sds%ndim.eq.3) then
          call compute_spline_3d (field%sds%dims(1),&
     &                            field%sds%dims(2),&
     &                            field%sds%dims(3),&
     &                            field%sds%scales(1)%f,&
     &                            field%sds%scales(2)%f,&
     &                            field%sds%scales(3)%f,&
     &                            field%sds%f,field%spl3)
        endif
      endif
!
      return
      end
!#######################################################################
      subroutine add_phi_point (s)
!
!-----------------------------------------------------------------------
!
! ****** Add a point in the phi dimension in a periodic manner
! ****** to the SDS in structure S.
!
!-----------------------------------------------------------------------
!
      use number_types
      use sds_def
      use constants
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      type(sds) :: s
!
!-----------------------------------------------------------------------
!
      real(r_typ), dimension(:,:,:), pointer :: f
      real(r_typ), dimension(:), pointer :: p
!
!-----------------------------------------------------------------------
!
! ****** For a 2D (axisymmetric) case, this is not needed.
!
      if (s%ndim.ne.3) return
!
! ****** Add a point to the phi dimension.
!
      allocate (f(s%dims(1),s%dims(2),s%dims(3)+1))
      allocate (p(s%dims(3)+1))
!
      f(:,:,1:s%dims(3))=s%f(:,:,:)
      f(:,:,s%dims(3)+1)=s%f(:,:,1)
!
      p(1:s%dims(3))=s%scales(3)%f(:)
      p(s%dims(3)+1)=s%scales(3)%f(1)+twopi
      s%dims(3)=s%dims(3)+1
!
      deallocate (s%f)
      deallocate (s%scales(3)%f)
!
      s%f=>f
      s%scales(3)%f=>p
!
      return
      end
!#######################################################################
      subroutine read_vignetting_function (fname)
!
!-----------------------------------------------------------------------
!
! ****** Read the vignetting function from file FNAME.
!
!-----------------------------------------------------------------------
!
      use ident
      use sds_def
      use params
      use vignetting_function
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      character(*) :: fname
!
!-----------------------------------------------------------------------
!
      type(sds) :: s
      integer :: ierr,i,j
!
!-----------------------------------------------------------------------
!
      if (verbose) then
        write (*,*)
        write (*,*) '### Reading the vignetting function'//&
     &              ' from file: ',trim(fname)
      end if
!
! ****** Read the vignetting function from the specified file.
!
      call rdhdf (fname,s,ierr)
!
      if (ierr.ne.0) then
        write (*,*)
        write (*,*) '### ERROR in ',cname,':'
        write (*,*) '### Could not read the vignetting function.'
        write (*,*) 'IERR (from RDHDF) = ',ierr
        write (*,*) 'File name: ',trim(fname)
        call exit (1)
      end if
!
! ****** Check that it contains a 1D field.
!
      if (s%ndim.ne.1) then
        write (*,*)
        write (*,*) '### ERROR in ',cname,':'
        write (*,*) '### The vignetting function file does'//&
     &              ' not contain a 1D field.'
        write (*,*) 'Number of dimensions = ',s%ndim
        write (*,*) 'File name: ',trim(fname)
        call exit (1)
      end if
!
! ****** Check that it contains scales.
!
      if (.not.s%scale) then
        write (*,*)
        write (*,*) '### ERROR in ',cname,':'
        write (*,*) '### The vignetting function file does'//&
     &              ' not contain scales.'
        write (*,*) 'File name: ',trim(fname)
        call exit (1)
      end if
!
      nvf=s%dims(1)
!
! ****** Check that it contains two or more points.
!
      if (nvf.lt.2) then
        write (*,*)
        write (*,*) '### ERROR in ',cname,':'
        write (*,*) '### The vignetting function file has'//&
     &              ' less than two points.'
        write (*,*) 'Number of points = ',nvf
        write (*,*) 'File name: ',trim(fname)
        call exit (1)
      end if
!
! ****** Transfer the vignetting function to 1D arrays for
! ****** convenience.
!
      allocate (r_vf(nvf))
      allocate (vf(nvf))
!
      r_vf=s%scales(1)%f
      vf=s%f(:,1,1)
!
      call deallocate_sds (s)
!
! ****** Check that the vignetting function scale is monotonically
! ****** increasing.
!
      do i=1,nvf-1
        if (r_vf(i+1).le.r_vf(i)) then
          write (*,*)
          write (*,*) '### ERROR in ',cname,':'
          write (*,*) '### The vignetting function scale is not'//&
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
!
! ****** Compute the spline coefficients if needed.
!
      if (cubic) then
        if (verbose) then
          write (*,*)
          write (*,*) '### Computing cubic spline coefficients for vf'
        end if
        call compute_spline_1d (nvf,&
     &                          r_vf,&
     &                          vf,vf_spl1)
      endif
!
      return
      end
!#######################################################################
      subroutine write_image (filename, fieldname, x, y, f)
!
!-----------------------------------------------------------------------
!
! ****** Write a 2D image (wrapper for wrdhdf_2d and error checking).
!
!-----------------------------------------------------------------------
!
      use ident
      use number_types
      use params
      use image_region
      use image_fields
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      character(512) :: filename
      character(*) :: fieldname
      real(r_typ), dimension(nx) :: x
      real(r_typ), dimension(ny) :: y
      real(r_typ), dimension(nx,ny) :: f
!
!-----------------------------------------------------------------------
!
      integer :: ierr
!
!-----------------------------------------------------------------------
!
      if (verbose) then
        write (*,*)
        write (*,*) '### Writing ',trim(fieldname),' to file: ',&
     &               trim(filename)
      end if
!
        call wrhdf_2d (filename,.true.,nx,ny,f,x,y,hdf32,ierr)
!
      if (ierr.ne.0) then
        write (*,*)
        write (*,*) '### ERROR in ',cname,':'
        write (*,*) '### Could not write the ',trim(fieldname),&
     &              ' output file:'
        write (*,*) 'IERR (from WRHDF_2D) = ',ierr
        write (*,*) 'File name: ',trim(fieldname)
        call exit (1)
      end if
!
      end
!#######################################################################
      subroutine vignette (f)
!
!-----------------------------------------------------------------------
!
! ****** Multiply the field F by the vignetting function, using
! ****** the radius R_MIN.
!
!-----------------------------------------------------------------------
!
      use ident
      use number_types
      use params
      use image_region
      use image_fields
      use vignetting_function
      use evaluate_spline_1d_interface
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      real(r_typ), dimension(nx,ny) :: f
!
!-----------------------------------------------------------------------
!
      real(r_typ), parameter :: one=1._r_typ
!
!-----------------------------------------------------------------------
!
      integer :: ierr,i,j,ivf,ivfp1
      real(r_typ) :: rv,alpha,vfv
!
!-----------------------------------------------------------------------
!
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
                write (*,*) '### An error occurred while applying'//&
     &                      ' the vignetting function.'
                write (*,*) 'Please check the vignetting'//&
     &                      ' function definition.'
                call exit (1)
              end if
              vfv=(one-alpha)*vf(ivf)+alpha*vf(ivfp1)
            endif
          end if
          f(i,j)=vfv*f(i,j)
        enddo
      enddo
!
      return
      end
!#######################################################################
      subroutine multiply_by_radial_power (f)
!
!-----------------------------------------------------------------------
!
! ****** Multiply the field F by R_MIN**POWER.
!
!-----------------------------------------------------------------------
!
      use ident
      use number_types
      use params
      use image_region
      use image_fields
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      real(r_typ), dimension(nx,ny) :: f
!
!-----------------------------------------------------------------------
!
      integer :: i,j
      real(r_typ) :: rv
!
!-----------------------------------------------------------------------
!
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
!
      return
      end
!#######################################################################
      subroutine transform_position (x,y,z,r,t,p)
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
      use number_types
      use constants
      use geometry
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      real(r_typ) :: x,y,z
      real(r_typ) :: r,t,p
!
!-----------------------------------------------------------------------
!
      real(r_typ), parameter :: zero=0.
!
!-----------------------------------------------------------------------
!
      real(r_typ) :: x1,y1,z1,x2,y2,z2
!
!-----------------------------------------------------------------------
!
      real(r_typ), external :: fold
!
!-----------------------------------------------------------------------
!
! ****** Apply the solar P angle transformation.
!
      x1=x
      y1= cos_p*y+sin_p*z
      z1=-sin_p*y+cos_p*z
!
! ****** Apply the solar B0 angle transformation.
!
      x2= cos_b0*x1-sin_b0*z1
      y2=y1
      z2= sin_b0*x1+cos_b0*z1
!
! ****** Convert to spherical coordinates.
!
      call c2s (x2,y2,z2,r,t,p)
!
! ****** Transform to the specified central meridian longitude.
!
      p=fold(zero,twopi,p+cml)

      return
      end
!#######################################################################
      function fold (x0,x1,x)
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
      use number_types
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      real(r_typ) :: fold
      real(r_typ) :: x0,x1,x
!
!-----------------------------------------------------------------------
!
      real(r_typ) :: xl
!
!-----------------------------------------------------------------------
!
      fold=x
!
      if (x0.eq.x1) return
!
      xl=x1-x0
!
      fold=mod(x-x0,xl)+x0
!
      if (fold.lt.x0) fold=fold+xl
      if (fold.ge.x1) fold=fold-xl
!
      return
      end
!#######################################################################
      subroutine get_kernels (r_min,rv,tv,pv,&
     &                        compute_pb,compute_b,&
     &                        pb_local,b_local,ds)
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
      use number_types
      use geometry
      use mas_fields
      use params, only: he_rho
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      real(r_typ) :: r_min
      real(r_typ) :: rv,tv,pv
      logical :: compute_pb,compute_b
      real(r_typ) :: pb_local
      real(r_typ) :: b_local
      real(r_typ) :: ds
!
!-----------------------------------------------------------------------
!
      real(r_typ), parameter :: two=2._r_typ
!
!-----------------------------------------------------------------------
!
      real(r_typ) :: n_e,cf,ef
      real(r_typ), dimension(3) :: drtp
!
!-----------------------------------------------------------------------
!
      real(r_typ), external :: ds_s2c
      real(r_typ), external :: cfunc_billings
      real(r_typ), external :: efunc_billings
!
!-----------------------------------------------------------------------
!
! ****** Interpolate the density.
!
      call interp_field_ds (rho%sds,rho%invtab,rv,tv,pv,n_e,drtp,&
     &                   rho%spl2,rho%spl3)
!
! ****** Convert N_E to physical units.
!
      n_e=n_e*fn_n/he_rho
!
! ****** Convert the mesh spacing from spherical coordinates
! ****** to a Cartesian spacing.
!
      ds=ds_s2c(rv,tv,drtp,twodee)
!
! ****** Evaluate the "Billings" radial functions.
!
      if (compute_pb.or.compute_b) then
        cf=cfunc_billings(rv)
      end if
!
      if (compute_b) then
        ef=efunc_billings(rv)
      end if
!
! ****** Obtain the local contribution to pB.
!
      if (compute_pb) then
        pb_local=cf*(r_min/rv)**2*n_e
      else
        pb_local=0.
      end if
!
! ****** Obtain the local contribution to B.
!
      if (compute_b) then
        b_local=(two*ef-cf*(r_min/rv)**2)*n_e
      else
        b_local=0.
      end if
!
      return
      end
!#######################################################################
      function cfunc_billings (r)
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
      use number_types
      use params, only: mu_ld
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      real(r_typ) :: r
      real(r_typ) :: cfunc_billings
!
!-----------------------------------------------------------------------
!
      real(r_typ), parameter :: one=1._r_typ
      real(r_typ), parameter :: three=3._r_typ
      real(r_typ), parameter :: eighth=.125_r_typ
!
!-----------------------------------------------------------------------
!
! ****** Intensity coefficient [cm^3], normalized to I0,
! ****** the central disk brightness, 2.49e10 [erg/cm^2/s/sr].
!
      real(r_typ), parameter :: a=8.69e-15_r_typ
!
!-----------------------------------------------------------------------
!
      real(r_typ) :: so,co,so_sq,co_sq
      real(r_typ) :: log_term
      real(r_typ) :: afunc,bfunc
!
!-----------------------------------------------------------------------
!
! ****** Protect against invalid R values.  The physical requirement
! ****** is that R.ge.1.
!
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
!
      afunc=co*so_sq
!
      bfunc=-eighth*( one-three*so_sq&
     &               -(co_sq/so)*(one+three*so_sq)*log_term)
!
      cfunc_billings=a*((one-mu_ld)*afunc+mu_ld*bfunc)
!
      return
      end
!#######################################################################
      function efunc_billings (r)
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
      use number_types
      use params, only: mu_ld
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      real(r_typ) :: r
      real(r_typ) :: efunc_billings
!
!-----------------------------------------------------------------------
!
      real(r_typ), parameter :: one=1._r_typ
      real(r_typ), parameter :: three=3._r_typ
      real(r_typ), parameter :: four=4._r_typ
      real(r_typ), parameter :: five=5._r_typ
      real(r_typ), parameter :: eighth=.125_r_typ
      real(r_typ), parameter :: one_third=one/three
      real(r_typ), parameter :: four_thirds=four/three
!
!-----------------------------------------------------------------------
!
! ****** Intensity coefficient [cm^3], normalized to I0,
! ****** the central disk brightness, 2.49e10 [erg/cm^2/s/sr].
!
      real(r_typ), parameter :: a=8.69e-15_r_typ
!
!-----------------------------------------------------------------------
!
      real(r_typ) :: so,co,so_sq,co_sq
      real(r_typ) :: log_term
      real(r_typ) :: cfunc,dfunc
!
!-----------------------------------------------------------------------
!
! ****** Protect against invalid R values.  The physical requirement
! ****** is that R.ge.1.
!
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
!
      cfunc=four_thirds-co-one_third*co*co_sq
!
      dfunc=eighth*( five+so_sq&
     &              -(co_sq/so)*(five-so_sq)*log_term)
!
      efunc_billings=a*((one-mu_ld)*cfunc+mu_ld*dfunc)
!
      return
      end
!#######################################################################
      subroutine c2s (x,y,z,r,t,p)
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
      use number_types
      use constants
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      real(r_typ) :: x,y,z
      real(r_typ) :: r,t,p
!
!-----------------------------------------------------------------------
!
      r=sqrt(x**2+y**2+z**2)
!
      if (r.eq.0.) then
        t=0.
      else
        t=acos(z/r)
      end if
!
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
!
      return
      end
!#######################################################################
      subroutine s2c (r,t,p,x,y,z)
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
      use number_types
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      real(r_typ) :: r,t,p
      real(r_typ) :: x,y,z
!
!-----------------------------------------------------------------------
!
      real(r_typ) :: st
!
!-----------------------------------------------------------------------
!
      st=sin(t)
      x=r*st*cos(p)
      y=r*st*sin(p)
      z=r*cos(t)
!
      return
      end
!#######################################################################
      function ds_s2c (r,t,drtp,axisymmetric)
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
      use number_types
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      real(r_typ) :: r,t
      real(r_typ), dimension(3) :: drtp
      logical :: axisymmetric
      real(r_typ) :: ds_s2c
!
!-----------------------------------------------------------------------
!
      real(r_typ) :: dr,dt,dp,st
!
!-----------------------------------------------------------------------
!
      dr=drtp(1)
      dt=drtp(2)
!
      if (axisymmetric) then
        ds_s2c=min(dr,r*dt)
      else
        dp=drtp(3)
        st=max(sin(t),sin(dt))
        ds_s2c=min(dr,r*dt,r*st*dp)
      end if
!
      return
      end
!#######################################################################
      subroutine interp_field_ds (fld,invtab,r,t,p,fv,ds,sp2,sp3)
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
      use number_types
      use types
      use params, only: cubic
      use spline_def
      use evaluate_spline_2d_interface
      use evaluate_spline_3d_interface
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      type(sds) :: fld
      type(vtab) :: invtab
      real(r_typ) :: r,t,p
      real(r_typ) :: fv
      real(r_typ), dimension(3) :: ds
      type(spl2d) :: sp2
      type(spl3d) :: sp3
!
!-----------------------------------------------------------------------
!
! ****** Note that the interpolated value is set to zero if it is
! ****** outside the bounds of the mesh.
!
      if (fld%ndim.eq.2) then
        call interp_2d (fld%dims(1),&
     &                  fld%dims(2),&
     &                  fld%scales(1)%f,&
     &                  fld%scales(2)%f,&
     &                  invtab,&
     &                  fld%f,r,t,fv,ds)
!
        if (cubic .and. (fv.ne.0.) ) then
          fv=evaluate_spline_2d(sp2,r,t,invtab%c(1),invtab%c(2))
          fv=max(fv,0.)
        endif
!
      else if (fld%ndim.eq.3) then
        call interp_3d (fld%dims(1),&
     &                  fld%dims(2),&
     &                  fld%dims(3),&
     &                  fld%scales(1)%f,&
     &                  fld%scales(2)%f,&
     &                  fld%scales(3)%f,&
     &                  invtab,&
     &                  fld%f,r,t,p,fv,ds)
        if (cubic .and. (fv.ne.0.) ) then
          fv=evaluate_spline_3d(sp3,r,t,p,&
     &                          invtab%c(1),invtab%c(2),invtab%c(3))
          fv=max(fv,0.)
        endif
      else
        fv=0.
        ds=huge(fv)
      end if
!
      return
      end
!#######################################################################
      subroutine interp_field (fld,invtab,r,t,p,fv,sp2,sp3)
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
      use number_types
      use types
      use params, only: cubic
      use spline_def
      use evaluate_spline_2d_interface
      use evaluate_spline_3d_interface
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      type(sds) :: fld
      type(vtab) :: invtab
      real(r_typ) :: r,t,p
      real(r_typ) :: fv
      type(spl2d) :: sp2
      type(spl3d) :: sp3
!
      real(r_typ), dimension(3) :: ds
!
!-----------------------------------------------------------------------
!
      if (fld%ndim.eq.2) then
        if (cubic) then
          fv=evaluate_spline_2d(sp2,r,t,invtab%c(1),invtab%c(2))
        else
          call interp_2d (fld%dims(1),&
     &                  fld%dims(2),&
     &                  fld%scales(1)%f,&
     &                  fld%scales(2)%f,&
     &                  invtab,&
     &                  fld%f,r,t,fv,ds)
        endif
!
      else if (fld%ndim.eq.3) then
        if (cubic) then
          fv=evaluate_spline_3d(sp3,r,t,p,&
     &                          invtab%c(1),invtab%c(2),invtab%c(3))
        else
          call interp_3d (fld%dims(1),&
     &                  fld%dims(2),&
     &                  fld%dims(3),&
     &                  fld%scales(1)%f,&
     &                  fld%scales(2)%f,&
     &                  fld%scales(3)%f,&
     &                  invtab,&
     &                  fld%f,r,t,p,fv,ds)
        endif
      endif
!
      return
      end
!
!#######################################################################
      subroutine interp_2d (nx,ny,x,y,inv,f,xv,yv,fv,ds)
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
      use ident
      use number_types
      use types
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      integer :: nx,ny
      real(r_typ), dimension(nx) :: x
      real(r_typ), dimension(ny) :: y
      type(vtab) :: inv
      real(r_typ), dimension(nx,ny) :: f
      real(r_typ) :: xv,yv,fv
      real(r_typ), dimension(2) :: ds
!
!-----------------------------------------------------------------------
!
      real(r_typ), parameter :: one=1._r_typ
!
!-----------------------------------------------------------------------
!
      logical :: outside
      integer :: i,j,ip1,jp1
      real(r_typ) :: ax,ay
!
!-----------------------------------------------------------------------
!
      real(r_typ), external :: mesh_spacing
!
!-----------------------------------------------------------------------
!
! ****** Find the cell that contains the interpolation point.
!
      outside=.false.
!
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
!
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
!
! ****** Get the mesh spacing at the interpolation point.
!
      ds(1)=mesh_spacing(nx,x,i,ip1,ax)
      ds(2)=mesh_spacing(ny,y,j,jp1,ay)
!
! ****** If the point is outside the mesh limits, set the
! ****** interpolated value to 0.
!
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
        fv= (one-ax)*((one-ay)*f(i  ,j  )+ay*f(i  ,jp1))&
     &     +     ax *((one-ay)*f(ip1,j  )+ay*f(ip1,jp1))
      end if
!
      return
      end
!#######################################################################
      subroutine interp_3d (nx,ny,nz,x,y,z,inv,f,xv,yv,zv,fv,ds)
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
      use ident
      use number_types
      use types
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      integer :: nx,ny,nz
      real(r_typ), dimension(nx) :: x
      real(r_typ), dimension(ny) :: y
      real(r_typ), dimension(nz) :: z
      type(vtab) :: inv
      real(r_typ), dimension(nx,ny,nz) :: f
      real(r_typ) :: xv,yv,zv,fv
      real(r_typ), dimension(3) :: ds
!
!-----------------------------------------------------------------------
!
      real(r_typ), parameter :: one=1._r_typ
!
!-----------------------------------------------------------------------
!
      logical :: outside
      integer :: i,j,k,ip1,jp1,kp1
      real(r_typ) :: ax,ay,az
!
!-----------------------------------------------------------------------
!
      real(r_typ), external :: mesh_spacing
!
!-----------------------------------------------------------------------
!
! ****** Find the cell that contains the interpolation point.
!
      outside=.false.
!
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
!
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
!
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
!
! ****** Get the mesh spacing at the interpolation point.
!
      ds(1)=mesh_spacing(nx,x,i,ip1,ax)
      ds(2)=mesh_spacing(ny,y,j,jp1,ay)
      ds(3)=mesh_spacing(nz,z,k,kp1,az)
!
! ****** If the point is outside the mesh limits, set the
! ****** interpolated value to 0.
!
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
        fv= (one-ax)*( (one-ay)*( (one-az)*f(i  ,j  ,k  )&
     &                           +     az *f(i  ,j  ,kp1))&
     &                +     ay *( (one-az)*f(i  ,jp1,k  )&
     &                           +     az *f(i  ,jp1,kp1)))&
     &     +     ax *( (one-ay)*( (one-az)*f(ip1,j  ,k  )&
     &                           +     az *f(ip1,j  ,kp1))&
     &                +     ay *( (one-az)*f(ip1,jp1,k  )&
     &                           +     az *f(ip1,jp1,kp1)))
      end if
!
      return
      end
!#######################################################################
      function mesh_spacing (nx,x,i,ip1,alpha)
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
      use number_types
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      integer :: nx
      real(r_typ), dimension(nx) :: x
      integer :: i,ip1
      real(r_typ) :: alpha
      real(r_typ) :: mesh_spacing
!
!-----------------------------------------------------------------------
!
      real(r_typ), parameter :: one=1._r_typ
!
!-----------------------------------------------------------------------
!
      integer :: im,ip
      real(r_typ) :: d,dxm,dxp
!
!-----------------------------------------------------------------------
!
      im=max( 1,i-1)
      ip=min(nx,i+1)
      d=real(ip-im)
      if (d.ne.0.) then
        dxm=(x(ip)-x(im))/d
      else
        dxm=huge(one)
      end if
!
      im=max( 1,ip1-1)
      ip=min(nx,ip1+1)
      d=real(ip-im)
      if (d.ne.0.) then
        dxp=(x(ip)-x(im))/d
      else
        dxp=huge(one)
      end if
!
      mesh_spacing=(one-alpha)*dxm+alpha*dxp
!
      return
      end
!#######################################################################
      subroutine build_inverse_tables (s,inv)
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
      use number_types
      use types
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      type(sds) :: s
      type(vtab) :: inv
!
!-----------------------------------------------------------------------
!
      integer :: i
!
!-----------------------------------------------------------------------
!
! ****** Use a number of points for the inverse interpolation table
! ****** equal to the number in the original scale.
!
      do i=1,s%ndim
        inv%c(i)%n=s%dims(i)
        allocate (inv%c(i)%f(inv%c(i)%n))
        call getinv (s%scales(i)%f,s%dims(i),inv%c(i))
      enddo
!
      return
      end
!#######################################################################
      subroutine getinv (x,n,tab)
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
      use ident
      use number_types
      use types
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      integer :: n
      real(r_typ), dimension(n) :: x
      type(itab) :: tab
!
!-----------------------------------------------------------------------
!
      real(r_typ), parameter :: one=1._r_typ
!
!-----------------------------------------------------------------------
!
      integer :: i,k,ierr,ip1
      real(r_typ) :: dx,xv,alpha,en
!
!-----------------------------------------------------------------------
!
! ****** For the special case when the table X has only one point,
! ****** the inverse table is not used.  It is thus loaded with
! ****** dummy values.
!
      if (n.eq.1) then
        tab%d=0.
        tab%f=1.
        return
      end if
!
! ****** Check that the number of points is valid.
!
      if (tab%n.le.1) then
        write (*,*)
        write (*,*) '### ERROR in ',cname,':'
        write (*,*) '### Error in GETINV:'
        write (*,*) '### Invalid number of points specified'//&
     &              ' for the inverse interpolation table.'
        write (*,*)
        write (*,*) 'Number of points = ',tab%n
        call exit (1)
      end if
!
! ****** Set the uniform interval to be used in the inverse
! ****** interpolation.
!
      dx=(x(n)-x(1))/(tab%n-one)
!
      if (dx.le.0.) then
        write (*,*)
        write (*,*) '### ERROR in ',cname,':'
        write (*,*) '### Error in GETINV:'
        write (*,*) '### Invalid interval for the inverse'//&
     &              ' interpolation table.'
        write (*,*)
        write (*,*) 'Interval = ',dx
        call exit (1)
      end if
!
      tab%d=one/dx
!
! ****** Build the inverse interpolation table.
!
      en=n
!
      do k=1,tab%n
        xv=x(1)+(k-one)*dx
        xv=max(xv,x(1))
        xv=min(xv,x(n))
        call interp (n,x,xv,i,ip1,alpha,ierr)
        if (ierr.ne.0) then
          write (*,*)
          write (*,*) '### ERROR in ',cname,':'
          write (*,*) '### Error in GETINV:'
          write (*,*) '### Error in building the inverse'//&
     &                ' interpolation table.'
          call exit (1)
        end if
        tab%f(k)=i+alpha
        tab%f(k)=max(tab%f(k),one)
        tab%f(k)=min(tab%f(k),en)
      enddo
!
      return
      end
!#######################################################################
      subroutine interpi (x,n,tab,xv,i,ip1,alpha)
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
      use ident
      use number_types
      use types
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      integer :: n
      real(r_typ), dimension(n) :: x
      type(itab) :: tab
      real(r_typ) :: xv
      integer :: i
      integer :: ip1
      real(r_typ) :: alpha
      intent(in) :: x,n,tab,xv
      intent(out) :: i,ip1,alpha
!
!-----------------------------------------------------------------------
!
      real(r_typ), parameter :: one=1._r_typ
!
!-----------------------------------------------------------------------
!
      integer :: ig
      real(r_typ) :: xi,fiv
!
!-----------------------------------------------------------------------
!
! ****** For the special case when the table has only one point,
! ****** the inverse table is not used.  In this case it is
! ****** necessary for XV to equal X(1) exactly, otherwise this
! ****** routine exits with an error.
!
      if (n.eq.1) then
        if (xv.eq.x(1)) then
          i=1
          ip1=1
          alpha=0.
        else
          go to 900
        end if
      end if
!
! ****** Get an estimate of the nearest grid point location in
! ****** the (uniform) inverse interpolation table.
!
      xi=one+(xv-x(1))*tab%d
      i=xi
      i=max(i,1)
      i=min(i,tab%n-1)
      alpha=xi-i
      fiv=(one-alpha)*tab%f(i)+alpha*tab%f(i+1)
!
! ****** Set IG to be the guess for the nearest grid point.
!
      ig=fiv
      ig=max(ig,1)
      ig=min(ig,n-1)
!
      if (xv.ge.x(ig)) then
!
! ****** Search forwards.
!
        do i=ig,n-1
          if (xv.ge.x(i).and.xv.le.x(i+1)) then
            alpha=(xv-x(i))/(x(i+1)-x(i))
            ip1=i+1
            return
          end if
        enddo
!
      else
!
! ****** Search backwards.
!
        do i=ig-1,1,-1
          if (xv.ge.x(i).and.xv.le.x(i+1)) then
            alpha=(xv-x(i))/(x(i+1)-x(i))
            ip1=i+1
            return
          end if
        enddo
!
      end if
!
! ****** ERROR: value not found in table.
!
  900 continue
!
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
!
      return
      end
!#######################################################################
      subroutine interp (n,x,xv,i,ip1,a,ierr)
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
      use number_types
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      integer :: n
      real(r_typ), dimension(n) :: x
      real(r_typ) :: xv
      integer :: i,ip1
      real(r_typ) :: a
      integer :: ierr
!
!-----------------------------------------------------------------------
!
      ierr=0
!
! ****** For the special case when the table has only one point,
! ****** it is necessary for XV to equal X(1) exactly, otherwise
! ****** this routine exits with an error.
!
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
!
! ****** Find the interval and compute the interpolation factor.
!
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
!
! ****** ERROR: the value was not found.
!
  900 continue
!
      write (*,*)
      write (*,*) '### ERROR in INTERP:'
      write (*,*) '### Value not found in table.'
      write (*,*) 'Value requested = ',xv
      write (*,*) 'Min table value = ',x(1)
      write (*,*) 'Max table value = ',x(n)
      write (*,*) 'Number of values in table = ',n
      ierr=1
!
      return
      end
!#######################################################################
      subroutine rtp2xyz (r,t,p,x,y,z)
!
!-----------------------------------------------------------------------
!
! ****** Convert spherical position to cartesian position.
!
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!
      use number_types
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      real(r_typ) :: r, t, p
      real(r_typ) :: x, y, z
!
!-----------------------------------------------------------------------
!
      real(r_typ) :: cph, sph, cth, sth
!
!-----------------------------------------------------------------------
!
      cph=cos(p)
      sph=sin(p)
      cth=cos(t)
      sth=sin(t)
      x = r*cph*sth
      y = r*sph*sth
      z = r*cth
!
      return
      end
!#######################################################################
      subroutine svtocv (t,p,vr,vt,vp,vx,vy,vz)
!
!-----------------------------------------------------------------------
!
! ****** Rotate Spherical Vector Components to Cartesian Components
!
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!
      use number_types
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      real(r_typ) :: t, p
      real(r_typ) :: vr, vt, vp
      real(r_typ) :: vx, vy, vz
!
!-----------------------------------------------------------------------
!
      real(r_typ) :: cph, sph, cth, sth
!
!-----------------------------------------------------------------------
!
      cph=cos(p)
      sph=sin(p)
      cth=cos(t)
      sth=sin(t)
      vx = vr*sth*cph + vt*cth*cph - vp*sph
      vy = vr*sth*sph + vt*cth*sph + vp*cph
      vz = vr*cth - vt*sth
!
      return
      end
!#######################################################################
      subroutine normalize(vec_in, vec_out)
!
!-----------------------------------------------------------------------
!
! ****** Normalize a vector to unit length.
!
!-----------------------------------------------------------------------
!
      use number_types
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      real(r_typ), dimension(3) :: vec_out
      real(r_typ), dimension(3) :: vec_in
      real(r_typ) :: mag
!
!-----------------------------------------------------------------------
!
      mag = sqrt(vec_in(1)**2 + vec_in(2)**2 + vec_in(3)**2)
      vec_out(1) = vec_in(1)/mag
      vec_out(2) = vec_in(2)/mag
      vec_out(3) = vec_in(3)/mag
!
      return
      end
!#######################################################################
      subroutine cross(vec1,vec2,vec_out)
!
!-----------------------------------------------------------------------
!
! ****** Compute the cross product of two vectors.
!
!-----------------------------------------------------------------------
!
      use number_types
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      real(r_typ), dimension(3) :: vec_out
      real(r_typ), dimension(3) :: vec1,vec2
!
!-----------------------------------------------------------------------
!
      vec_out(1) = vec1(2)*vec2(3) - vec1(3)*vec2(2)
      vec_out(2) = vec1(3)*vec2(1) - vec1(1)*vec2(3)
      vec_out(3) = vec1(1)*vec2(2) - vec1(2)*vec2(1)
!
      return
      end
!#######################################################################
      subroutine get_vector_projection(ex,ey,ez,rv,tv,pv,vec)
!
!-----------------------------------------------------------------------
!
! ****** Get the vector field at this location and project it into
! ****** the observing frame defined with unit vectors ex, ey, ez.
!
!-----------------------------------------------------------------------
!
      use number_types
      use mas_fields
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      real(r_typ), dimension(3) :: ex,ey,ez
      real(r_typ) :: rv, tv, pv
      real(r_typ), dimension(3) :: vec
!
!-----------------------------------------------------------------------
!
      real(r_typ) :: vrv, vtv, vpv
      real(r_typ) :: vx, vy, vz
!
!-----------------------------------------------------------------------
!
! ****** Interpolate the velocity components.
!
      call interp_field (vr%sds,vr%invtab,rv,tv,pv,vrv,vr%spl2,vr%spl3)
      call interp_field (vt%sds,vt%invtab,rv,tv,pv,vtv,vt%spl2,vt%spl3)
      call interp_field (vp%sds,vp%invtab,rv,tv,pv,vpv,vp%spl2,vp%spl3)
!
! ****** Get the cartesian vector components.
!
      call svtocv( tv, pv, vrv, vtv, vpv, vx, vy, vz)
!
! ****** Project the vector into the observer frame.
!
      vec(1) = ex(1)*vx + ex(2)*vy + ex(3)*vz
      vec(2) = ey(1)*vx + ey(2)*vy + ey(3)*vz
      vec(3) = ez(1)*vx + ez(2)*vy + ez(3)*vz
!
      return
      end
!#######################################################################
      subroutine get_unit_vectors(losx,losy,losz,ex,ey,ez)
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
      use number_types
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      real(r_typ), parameter :: zero=0._r_typ
      real(r_typ), parameter :: one=1._r_typ
!
!-----------------------------------------------------------------------
!
      real(r_typ), intent(in) :: losx, losy, losz
      real(r_typ), dimension(3), intent(out) :: ex, ey, ez
!
!-----------------------------------------------------------------------
!
      real(r_typ), dimension(3) :: vec
!
!-----------------------------------------------------------------------
!
! ****** x-hat is sun-to-observer unit vector --> normalize los_vec.
!
      vec(1)=losx
      vec(2)=losy
      vec(3)=losz
      call normalize(vec,ex)
!
! ****** y-hat is perp to heliographic z axis and observer x-hat.
!
      vec(1)=zero
      vec(2)=zero
      vec(3)=one
      call cross(vec,ex,ey)
      call normalize(ey,ey)
!
! ****** z-hat is perp to heliographic x-hat and y-hat.
!
      call cross(ex,ey,ez)
      call normalize(ez,ez)
!
      return
      end
!#######################################################################
      subroutine add_to_integrals(i,j,e_pb,e_b,ds,ex,ey,ez,rv,tv,pv,x)
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
      use number_types
      use params
      use constants
      use image_fields
      use mas_fields
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      integer, intent(in) :: i, j
      real(r_typ), intent(in) :: e_pb, e_b, ds
      real(r_typ), dimension(3), intent(in) :: ex, ey, ez
      real(r_typ), intent(in) :: rv, tv, pv, x
!
!-----------------------------------------------------------------------
!
      real(r_typ) :: var_now
      real(r_typ) :: e_now
      real(r_typ) :: angle_now
      real(r_typ), dimension(3) :: vec
!
!-----------------------------------------------------------------------
!
! ****** Select the weighting.
!
      if (weight_integral_by_b) then
        e_now=e_b
      else
        e_now=e_pb
      endif
!
! ****** Scalar field.
!
      if (compute_scalar_integration) then
        call interp_field (sfield%sds,sfield%invtab,rv,tv,pv,&
     &            var_now,sfield%spl2,sfield%spl3)
        sf_avg(i,j)=sf_avg(i,j) + ds*e_now*var_now
      endif
!
! ****** Average LOS angle.
!
      if (compute_angle_integration) then
        angle_now = atan(x/r_min(i,j))/degtorad
        angle_avg(i,j)=angle_avg(i,j) + ds*e_now*angle_now
      endif
!
! ****** Vector field. Mapping from observer frame to image frame is:
!        Observer_x -> Image_LOS
!        Observer_y -> Image_x
!        Observer_z -> Image_y
!
      if (compute_vector_integration) then
         call get_vector_projection(ex,ey,ez,rv,tv,pv,vec)
         vlos_avg(i,j) = vlos_avg(i,j) + ds*e_now*vec(1)
         vx_avg(i,j) = vx_avg(i,j) + ds*e_now*vec(2)
         vy_avg(i,j) = vy_avg(i,j) + ds*e_now*vec(3)
      endif
!
      return
      end
!#######################################################################
      subroutine normalize_integrals
!
!-----------------------------------------------------------------------
!
! ****** Normalize the weighted integrals to get averages.
!
! ****** Pixels that are occulted or outside the grid get the disk val.
!
!-----------------------------------------------------------------------
!
      use number_types
      use params
      use geometry
      use image_region
      use image_fields
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      real(r_typ), parameter :: one=1._r_typ
!
!-----------------------------------------------------------------------
!
      integer :: i, j
      real(r_typ) :: e_now
!
!-----------------------------------------------------------------------
!
! ****** Loop over pixels, normalize the integrals.
!
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
!
! ****** Divide out emissivity integral to get averages.
!
          if (do_integral_weighting) then
            if (weight_integral_by_b) then
              e_now=b(i,j)
            else
              e_now=pb(i,j)
            endif
!
            if (compute_scalar_integration) then
              sf_avg(i,j)=sf_avg(i,j)/e_now
            endif
!
            if (compute_angle_integration) then
              angle_avg(i,j)=angle_avg(i,j)/e_now
            endif
!
            if (compute_vector_integration) then
              vlos_avg(i,j)=vlos_avg(i,j)/e_now
              vx_avg(i,j)=vx_avg(i,j)/e_now
              vy_avg(i,j)=vy_avg(i,j)/e_now
            endif
          endif
        enddo
      enddo
!
      return
      end

           subroutine print_help_text
!
!-----------------------------------------------------------------------
!
! ****** Print the help / syntax information.
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
        write (*,*)
        write (*,*)
        write (*,*) '### Help for GETPB Special Integration Methods:'
        write (*,*)
        write (*,*) 'GETPB can now compute emissivity weighted'//&
     &              ' averages of other quantities.'
        write (*,*) 'Here are the additions flags/options: '
        write (*,*)
        write (*,*)
        write (*,*) 'SCALAR FIELD AVERAGING:'
        write (*,*)
        write (*,*) ' Use -scalar <file> to specify a 2D or 3D file'//&
     &              ' with a scalar field.'
        write (*,*)
        write (*,*) ' GETPB will integrate the scalar along the LOS'//&
     &              ' to get the weighted average.'
        write (*,*)
        write (*,*) ' Specify the output file with -avg_scalar <file>'
        write (*,*)
        write (*,*)
        write (*,*) 'VECTOR FIELD COMPONENT AVERAGING:'
        write (*,*)
        write (*,*) ' Use -vr <file>, -vt <file>, -vp <file> to'//&
     &              ' specify a 2D or 3D vector field.'
        write (*,*)
        write (*,*) ' GETPB will project the rtp vector field into'//&
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
        write (*,*) ' Use -avg_los_angle <file> to output the'//&
     &              ' weighted average LOS angle.'
        write (*,*)
        write (*,*) ' Here the "LOS angle" is the angle between a'//&
     &              ' position on the LOS and'
        write (*,*) ' the closest approach of the LOS (rmin).'
        write (*,*)
        write (*,*) ' A positive LOS angle is in front of rmin,'//&
     &              ' negative is behind.'
        write (*,*)
        write (*,*) ' For plane-parallel integration, rmin lies'//&
     &              ' in the plane of sky.'
        write (*,*)
        write (*,*)
        write (*,*) 'EMISSIVITY WEIGHTING:'
        write (*,*)
        write (*,*) ' GETPB can weight the LOS averages'//&
     &              ' by b OR pb (not both at once).'
        write (*,*)
        write (*,*) ' The default is to weight by the pb emissivity.'
        write (*,*)
        write (*,*) ' Set the flag -avg_using_b to weight by b'
        write (*,*)
!
      call exit (1)
!
      end