!
!-----------------------------------------------------------------------
!
! ****** Get the field line mapping for MAS code runs, computing
! ****** the structural quantities K and Q.
!
!-----------------------------------------------------------------------
!
! ****** Updates and bug fixes:
!
!        08/30/2005, ZM, Version 1.06:
!
!         - Converted MAPFL v1.05 into a tool.
!         - Improved the tracing accuracy of short field lines
!           (e.g., those near the neutral line).
!         - Improved the field line integrator.
!         - Added the ability to perform a 3D mapping of the
!           field lines.
!
!        06/16/2006, ZM, Version 1.07:
!
!         - Fixed the effect of rondoff errors in generating the
!           meshes for the calculation which caused initial points
!           to lie outside the computation domain.
!
!        06/21/2006, ZM, Version 1.08:
!
!         - Allowed the ability of specifying planes for the
!           field line mapping region and the 3D mapping.
!           This extends the capability of the tool to get
!           Q on 2D slices.
!
!        07/07/2006, ZM, Version 1.09:
!
!         - Added the ability to do new POT3D files,
!           new-old POT3D files, and old POT3D files.
!
!        08/22/2006, ZM, Version 1.10:
!
!         - Corrected a bug in the computation of Q for field
!           lines that wrap periodically at the join between
!           phi=0 and phi=2*pi.
!         - Corrected the computation of Q on open field lines
!           to take into account the difference in the radii
!           of the initial and final field line footpoints.
!         - Added the ability to use cubic spline interpolation.
!           This increases the storage requirements substantially,
!           and makes the computation significantly slower.
!
!        09/17/2006, ZM, Version 1.11:
!
!         - Reverted to the Cartesian field line integrator.
!           The accuracy of the spherical integrator was being
!           cast into doubt, though it was not proven to be
!           bad.
!         - Improved the accuracy with which the final point
!           (i.e., the end point at the r boundaries) was
!           being clipped to the boundary.
!         - Changed the computation of Q to be staggered with
!           respect to the mapping quantities.
!         - Fixed the backward mapping routine to compute Q.
!
!        09/29/2006, ZM, Version 1.12:
!
!         - Improved the field line integrator to use a variable
!           step size.  It is now possible to select either a uniform
!           step size for the field line integration, as before,
!           or to use a variable step size based on the local radius
!           of curvature of the field line, and the local mesh size
!           of the B files along the field line trace.
!         - This has changed the format of the input file.
!
!        01/28/2007, ZM, Version 1.13:
!
!         - Added the ability to compute the mapping on 2D slices
!           through the 3D volume.
!         - This has changed the format of the input file.
!         - Cleaned up some formatting of the output.
!
!        02/19/2008, ZM, Version 1.14:
!
!         - Changed the default behavior that terminated the code
!           when more than 100 bad field line traces were found.
!           By default this is now disabled, but can be put back
!           by setting the variable MAX_BAD_FIELDLINES.
!
!        04/06/2009, ZM, Version 1.15:
!
!         - Added the ability to use an analytic function to define
!           the magnetic field.
!         - Corrected bugs involving how the expansion factor and Q
!           were computed for the backward trace.
!         - Performed a cosmetic clean-up of the code.
!
!        04/21/2009, ZM, Version 1.16:
!
!         - Added the ability to use Slava's new formulation
!           to compute Q directly on slices within the domain
!           by tracing a bundle of field lines from points in
!           the domain forward and backward to the boundaries.
!         - Cleaned up the way mapping along "slices" in the 3D
!           volume is implemented.  These "slices" can now be
!           lines (i.e., 1D files), 2D slices, or 3D volumes.
!           These are all defined by reading in rectilinear
!           files (1D, 2D, or 3D) that contain the (r,t,p) or
!           (x,y,z) starting coordinates for the mapping.
!         - Added a check to make sure that the input file
!           has the correct number of lines.
!
!        02/18/2010, ZM, Version 1.17:
!
!         - Made the ability to trace from a slice more flexible.
!           It is now possible to map from a slice along the
!           forward and backward directions separately.  It is
!           also possible to specify this direction to be either
!           along B or along the direction of increasing radius.
!         - Added the ability to directly compute coronal hole
!           maps at a given radius.  These are coded by magnetic
!           field polarity.
!
!        04/12/2010, ZM, Version 1.18:
!
!         - Added the ability of specifying the r, t, and p meshes
!           to be used for the calculation.  These can be specified
!           as 1D HDF files or as uniform meshes.
!         - Allowed the phi coordinate to be outside the [0,2*pi]
!           range in the input file.  It is now properly wrapped
!           into the main interval during the calculation.
!
!        04/26/2010, ZM, Version 1.19:
!
!         - Added the ability to compute 3D coronal hole maps.
!           These are useful to compute 3D open/closed field regions
!           for use, perhaps, in developing heating masks.
!         - Removed the writing of warning messages about field
!           lines that did not reach the boundaries in the coronal
!           hole map traces.  This was not really necessary since
!           such field lines are already flagged (with values of
!           "-2") in the coronal hole maps.
!
!        02/07/2011, ZM, Version 1.20:
!
!         - Added a multi-threading capability using OpenMP.
!           This version runs in parallel.  This required a
!           restructuring of the code to improve the modularity,
!           to make the code thread-safe, and to improve the
!           parallel performance.
!         - It is necessary to use "thread-safe" programming
!           in the parallel sections from now on.
!
!        05/18/2011, ZM, Version 1.21:
!
!         - Corrected a bug in the periodic wrapping of the phi
!           coordinate when it was outside the range [0,2*pi].
!           It was not being wrapped to the main periodic interval
!           in the GETB routine.
!
!        08/17/2012, ZM, Version 1.22:
!
!         - Fixed a minor bug that was discovered by compiling
!           with GFORTRAN.
!
!        04/29/2013, ZM, Version 1.23:
!
!         - Changed the interpretation of the value of DEBUG_LEVEL
!           for debugging output.  This was done to optimize
!           Slava's output of field line traces.  The functionality
!           of the program was not otherwise modified.
!
!        05/13/2013, ZM, Version 1.24:
!
!         - Added the ability to write field line traces to HDF
!           files.  This option can be selected when tracing from
!           a "slice".  It is not intended to be used when doing
!           very high resolution mapping computations, since it
!           can produce a large amount of data.  To get the HDF
!           files, set DEBUG_LEVEL.ge.2.  The field lines traces
!           for the forward and/or backward trace will be written
!           to HDF files named:
!
!             field_line_traces_forward_rtp.hdf
!             field_line_traces_forward_xyz.hdf
!             field_line_traces_backward_rtp.hdf
!             field_line_traces_backward_xyz.hdf
!
!        09/25/2015, ZM, Version 1.25:
!
!         - Changed the way short field lines are treated.
!           Because of strange behavior noticed by Slava in a
!           certain case, I changed the way the step size for
!           "short" field lines was controlled.  It turned out
!           that a "short" field line could become a very long
!           field line when retraced!  Previously, these field
!           lines were retraced with the miniumum step size, which
!           led to a very long field line with lots of points.
!           I relaxed this requirement on the step size once the
!           number of points in the retraced field line exceeded
!           the minimum number of points.
!
!        01/20/2016, RC, Version 1.26:
!
!         - Modified tracefl() routine to eliminate the use of "goto"s.
!         - Minor speed improvements (~3%).
!         - Manually added spline library into code to allow
!           for the development of future optimizations.
!
!        02/09/2016, ZM, Version 1.26ZM:
!
!         - Fixed the check that adds a phi point to the B arrays
!           for cases when all the B components are on the same
!           mesh.  The code now checks if the point at phi=2*pi is
!           already present, and only adds a point if it is not.
!         - Fixed the storage of field line traces in the field line
!           buffers.  Previously, when "short" field lines were
!           detected, the buffers were not reinitialized when these
!           short field lines were retraced with a smaller step size,
!           leading to incorrect saved field line traces.
!         - Fixed the computation of Q on slices for field line
!           footpoints that approach the poles in routine GETQ.
!           The previous differencing was not accurate for field
!           line footpoints near the poles.  The new scheme switches
!           to Cartesian basis vectors in a small neighborhood of
!           the poles.
!
!        03/25/2016, ZM, Version 1.27ZM:
!
!         - Added code to snap the theta and phi limits of the B
!           meshes to the correct values.  Previously, these could
!           be slightly off due to the precision limitations of
!           32-bit HDF files, leading to possible errors in the
!           tracing of a (very small) subset of field lines.
!
!        04/19/2016, ZM, Version 1.28ZM:
!
!         - Added the ability to write the field line length when
!           mapping forward/backward, and also when computing
!           Q on a slice.
!         - Added the ability to write traces of field lines traced
!           from a slice to individual HDF output files.  These are
!           named:
!
!              <root>_f_####.hdf
!              <root>_b_####.hdf
!
!           where <root> is a string specified in the input file,
!           and #### is a four-digit sequence number (0001, 0002,
!           etc.).  The "f" and "b" labels designate forward/backward
!           traces, as requested in the input file.  There is a flag
!           in the input file to select whether Cartesian (x,y,z)
!           coordinates or spherical (r,t,p) coordinates are
!           written.
!         - Removed the ability to write field lines to a single HDF
!           file when the debug level was set to 2.  This is
!           superseded by the new capability to write the traces to
!           individual HDF files.
!         - COMPTIBILITY CHANGE:  The changes in this version have
!           modified the format
!           of the input file, so unfortunately this version is
!           not backward compatible.
!
!        01/17/2018, ZM, Version 1.29ZM:
!
!         - Changed the way the domain limits are set from the
!           scales in the magnetic field files.  This is an attempt
!           to make sure that roundoff does not cause points that
!           are near theta=pi and phi=2*pi to end up outside the
!           domain.
!
!        01/18/2018, RC, Version 1.30:
!
!         - Merged in 1.26ZM, 1.27ZM, 1.28ZM, and 1.29ZM changes.
!         - Removed the small speed improvement from 1.26 for
!           consistency reasons.
!
!        03/19/2018, RC, Version 1.31:
!
!         - Small change to allow HDF5 trace output files if
!           input b files are HDF5.
!
!        05/10/2018, RC, Version 1.32:
!
!         - Added OpenMP to 3D cubic spline coef calculation.
!
!        05/20/2019, RL, Version 1.33:
!
!         - Added 3D dips map.
!
!        06/28/2019, RL, Version 2.00:
!
!         - COMPATIBILITY CHANGE!!!!
!           Replaced formatted input file with namelist!
!
!        07/01/2019, RC, Version 2.01:
!
!         - Added option to output Slava's signed-log10
!           of Q for forward and backward traces. To use, set
!           SLOGQFFILE and/or SLOGQBFILE to desired output filenames.
!         - Set default values for all namelist parameters and updated
!           mapfl.dat to reflect the default values.
!         - Added error-checking for namelist.
!
!        07/11/2019, RC, Version 2.02:
!
!         - Fixed bug in writing out backward mapping slogq.
!
!        07/12/2019, RC, Version 2.03:
!
!         - Fixed bug in writing out backward mapping slogq (again).
!
!        08/19/2019, RL, Version 2.04:
!
!         - Introduced capability to integrate scalar field along
!           field lines.
!           Speficy INTEGRATE_ALONG_FL=.true. and the name of the
!           file with the field in SCALAR_INPUT_FILE
!           Either TRACE_FWD, or TRACE_BWD, or TRACE_SLICE and
!           COMPUTE_Q_ON_SLICE, must be set true.
!           Specify the file with integrated quantity in LFFILE
!           (for TRACE_FWD or TRACE_BWD true) or in
!           SLICE_LENGTH_OUTPUT_FILE (for TRACE_SLICE and
!           COMPUTE_Q_ON_SLICE true).
!
!        03/02/2021, AP,CD Version 2.0.5:
!
!         - Small modifications for Python/f2py cross compilation.
!         - Debug statements in mesh detection, tweak to OLD MAS check.
!         - Changed version numbering to standard style.
!
!-----------------------------------------------------------------------
!
!#######################################################################
!
! ****** This tool uses modules from ZM's tools library.
!
!#######################################################################
      module ident
!
      character(*), parameter :: cname='MAPFL'
      character(*), parameter :: cvers='2.0.5'
      character(*), parameter :: cdate='03/02/2021'
!
      end module
!#######################################################################
      module spline_def
!
!-----------------------------------------------------------------------
! ****** Definition of cubic spline data structures.
!-----------------------------------------------------------------------
!
      use number_types
!
      implicit none
!
! ***** 1D spline structure.
!
      type :: spl1d
        integer :: nx
        real(r_typ), dimension(:), pointer :: x
        real(r_typ), dimension(:), pointer :: f
        real(r_typ), dimension(:), pointer :: fxx
      end type
!
! ***** 2D spline structure.
!
      type :: spl2d
        integer :: nx
        integer :: ny
        real(r_typ), dimension(:), pointer :: x
        real(r_typ), dimension(:), pointer :: y
        real(r_typ), dimension(:,:), pointer :: f
        real(r_typ), dimension(:,:), pointer :: fxx
        real(r_typ), dimension(:,:), pointer :: fyy
        real(r_typ), dimension(:,:), pointer :: fxxyy
      end type
!
! ***** 3D spline structure.
!
      type :: spl3d
        integer :: nx
        integer :: ny
        integer :: nz
        real(r_typ), dimension(:), pointer :: x
        real(r_typ), dimension(:), pointer :: y
        real(r_typ), dimension(:), pointer :: z
        real(r_typ), dimension(:,:,:), pointer :: f
        real(r_typ), dimension(:,:,:), pointer :: fxx
        real(r_typ), dimension(:,:,:), pointer :: fyy
        real(r_typ), dimension(:,:,:), pointer :: fzz
        real(r_typ), dimension(:,:,:), pointer :: fxxyy
        real(r_typ), dimension(:,:,:), pointer :: fxxzz
        real(r_typ), dimension(:,:,:), pointer :: fyyzz
        real(r_typ), dimension(:,:,:), pointer :: fxxyyzz
      end type
!
      end module
!#######################################################################
      module invint_def
!
!-----------------------------------------------------------------------
! ****** Definition of an inverse interpolation table data structure.
!-----------------------------------------------------------------------
!
      use number_types
!
      implicit none
!
      type :: itab
        integer :: n
        real(r_typ), dimension(:), pointer :: f
        real(r_typ) :: d
      end type
!
      end module
!#######################################################################
      module locate_interval_interface
      interface
        function locate_interval (n,x,xv,tab,ierr)
        use number_types
        use invint_def
        implicit none
        integer :: n
        real(r_typ), dimension(n) :: x
        real(r_typ) :: xv
        type(itab), optional :: tab
        integer, optional :: ierr
        integer :: locate_interval
        end
      end interface
      end module
!#######################################################################
      module evaluate_spline_3d_interface
      interface
        function evaluate_spline_3d (s,x,y,z,tabx,taby,tabz)
        use number_types
        use spline_def
        use invint_def
        use locate_interval_interface
        type(spl3d) :: s
        real(r_typ) :: x,y,z
        type(itab), optional :: tabx,taby,tabz
        real(r_typ) :: evaluate_spline_3d
        end
      end interface
      end module
!#######################################################################
      module debug
!
      implicit none
!
! ****** Debugging level.
!
      integer :: debug_level=0
!
      end module
!#######################################################################
      module constants
!
!-----------------------------------------------------------------------
! ****** Constants.
!-----------------------------------------------------------------------
!
      use number_types
!
      implicit none
!
      real(r_typ), parameter :: pi=3.1415926535897932_r_typ
      real(r_typ), parameter :: halfpi=.5_r_typ*pi
      real(r_typ), parameter :: twopi=2._r_typ*pi
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
! ****** Maximum number of dimensions.
!
      integer, parameter, private :: ndim_max=3
!
! ****** Inverse interpolation table structure for a vector field.
!
      type :: vtab
        type(itab), dimension(ndim_max) :: c
      end type
!
! ****** Vector spline structure.
!
      type :: vspl3d
        type(spl3d) :: r
        type(spl3d) :: t
        type(spl3d) :: p
      end type
!
! ****** Magnetic field vector structure.
!
      type :: vec
        type(sds) :: r
        type(sds) :: t
        type(sds) :: p
        integer :: nrs
        integer :: nts
        integer :: nps
        real(r_typ), dimension(:), pointer :: rs
        real(r_typ), dimension(:), pointer :: ts
        real(r_typ), dimension(:), pointer :: ps
        real(r_typ) :: lim0(ndim_max)
        real(r_typ) :: lim1(ndim_max)
        type(vtab), dimension(ndim_max) :: inv
        real(r_typ), dimension(:), pointer :: drs
        real(r_typ), dimension(:), pointer :: dts
        real(r_typ), dimension(:), pointer :: dps
        real(r_typ), dimension(:), pointer :: sts
        type(itab) :: rs_invtab
        type(itab) :: ts_invtab
        type(itab) :: ps_invtab
        logical :: cubic=.false.
        type(vspl3d) :: spl
        logical :: b_is_32bit
      end type
!
! ****** Data structure to hold file names of a vector component.
!
      type :: vfile
        character(512) :: r=' '
        character(512) :: t=' '
        character(512) :: p=' '
      end type
!
! ****** Initial size for the field line trace buffer.
!
      integer, parameter, private :: fl_buffer_size=1000
!
! ****** Trajectory structure definition.
!
      type :: traj
        integer :: ndim
        integer :: initial_size=fl_buffer_size
        integer :: size
        integer :: npts
        type(rp1d), dimension(:), pointer :: x
      end type
!
! ****** Dual representation Cartesian and spherical position vector.
!
      type :: csvec
        real(r_typ), dimension(3) :: c
        real(r_typ), dimension(3) :: s
      end type
!
! ****** "Inside domain" structure.
!
      type :: inout
        logical :: domain
        logical :: r0
        logical :: r1
        logical :: t0
        logical :: t1
        logical :: p0
        logical :: p1
        logical :: r
        logical :: t
        logical :: p
      end type
!
! ****** Field line integration parameters.
!
      type :: flparam
        logical :: variable=.true.
        integer :: direction
        logical :: direction_is_along_b
        real(r_typ) :: min=0.001_r_typ
        real(r_typ) :: max=0.1_r_typ
        real(r_typ) :: over_rc=0.0025_r_typ
        real(r_typ) :: lmax=100.0_r_typ
        logical :: limit_by_local_mesh=.true.
        real(r_typ) :: local_mesh_factor=1.0_r_typ
        real(r_typ) :: max_increase_factor
        real(r_typ) :: max_decrease_factor
        integer :: short_fl_min_points
        integer :: short_fl_max_tries
        real(r_typ) :: short_fl_shrink_factor
        real(r_typ) :: predictor_min_clip_fraction
      end type
!
! ****** Preferences for Python/f2py (used in other source files).
!
      type :: preferences
        logical :: cubic_vec_field
        logical :: var_dstep
        real(r_typ) :: dstep
        logical :: auto_minmax_dstep
        real(r_typ) :: min_dstep
        real(r_typ) :: max_dstep
        real(r_typ) :: dstep_mult
        logical :: limit_by_local_mesh
        real(r_typ) :: local_mesh_factor
        real(r_typ) :: max_length
        logical :: direction_along_vec_field
        logical :: trace_from_slice_forward
        logical :: trace_from_slice_backward
      end type
!
      end module
!#######################################################################
      module mesh
!
!-----------------------------------------------------------------------
! ****** Meshes.
!-----------------------------------------------------------------------
!
      use number_types
!
      implicit none
!
! ****** Secondary (r,t,p) meshes.
!
      integer :: nrss=54
      integer :: ntss=181
      integer :: npss=361
!
      real(r_typ), dimension(:), pointer :: rss
      real(r_typ), dimension(:), pointer :: tss
      real(r_typ), dimension(:), pointer :: pss
!
      end module
!#######################################################################
      module field
!
!-----------------------------------------------------------------------
! ****** Magnetic field storage.
!-----------------------------------------------------------------------
!
      use number_types
      use types
!
      implicit none
!
! ****** Structure that holds the magnetic field.
!
      type(vec) :: b
!
      end module
!#######################################################################
      module vars
!
!-----------------------------------------------------------------------
! ****** Input variables, switches, etc.
!-----------------------------------------------------------------------
!
      use number_types
      use sds_def
      use types
!
      implicit none
!
! ****** Spherical geometry domain limits.
!
      real(r_typ) :: domain_r_min=1._r_typ
      real(r_typ) :: domain_r_max=300._r_typ
!
! ****** New output mesh flags.
!
      logical :: new_r_mesh=.false.
      logical :: new_t_mesh=.false.
      logical :: new_p_mesh=.false.
!
! ****** New output mesh limits.
!
      real(r_typ) :: r0=0.
      real(r_typ) :: r1=0.
      real(r_typ) :: t0=0.
      real(r_typ) :: t1=0.
      real(r_typ) :: p0=0.
      real(r_typ) :: p1=0.
!
! ****** Flag to use tri-cubic interpolation (when .TRUE.) or
! ****** simple linear interpolation (when .FALSE.) to
! ****** interpolate B between mesh points.
!
      logical :: cubic=.false.
!
! ****** Field line integration.
!
      type(flparam) :: ds
      logical :: set_ds_automatically=.true.
      real(r_typ) :: dsmult=1.0_r_typ
!
! ****** Flag to request a mapping on a slice.
!
      logical :: trace_slice=.false.
!
! ****** Parameters for the slice mapping.
!
      logical :: slice_coords_are_xyz=.false.
!
! ****** Flag to specify whether the tracing direction is along B
! ****** or along the direction of increasing radius.
!
      logical :: trace_slice_direction_is_along_b=.true.
!
! ****** Flags to request tracing in the forward and backward
! ****** directions.
!
      logical :: trace_from_slice_forward=.false.
      logical :: trace_from_slice_backward=.false.
!
! ****** Names of the slice coodrinates.
!
      character, dimension(3) :: slice_coord_name
!
! ****** Structures that hold the slice coordinates.
!
      type (sds) :: slice_c1,slice_c2,slice_c3
!
! ****** Increment to compute Q directly on a slice.
!
      real(r_typ) :: q_increment_h=0.0001_r_typ
!
! ****** Flag to request a 3D mapping.
!
      logical :: trace_3d=.false.
!
! ****** Switch to use an analytic function to define the magnetic
! ****** field.
!
      logical :: use_analytic_function=.false.
!
! ****** Flag to use 32-bit HDF output files.
!
      logical :: hdf32=.true.
!
! ****** Flag to write field line traces originating from a slice
! ****** to individual HDF output files.
!
      logical :: write_traces_to_hdf=.false.
!
! ****** String used for the root file name for the field line
! ****** traces.
!
      character(64) :: write_traces_root='fl'
!
! ****** Flag to write the Cartesian (x,y,z) coordinates for
! ****** field line traces.
!
      logical :: write_traces_as_xyz=.false.
!
! ***** File type for field line tracing output.
!
      character(3) :: fmt='hdf'
!
! ***** Number of segnments along which to search for dips.
!
      integer :: ns_dips=10
!
! ****** Integrate scalar field along field line.
!
      logical :: integrate_along_fl=.false.
!
      end module
!#######################################################################
      module diags
!
!-----------------------------------------------------------------------
! ****** Variables that control diagnostic output.
!-----------------------------------------------------------------------
!
      use number_types
!
      implicit none
!
! ****** Number of iterations between prints of diagnostics
! ****** during execution.
!
      integer :: diagnostic_interval=1000
!
      end module
!#######################################################################
      module files
!
!-----------------------------------------------------------------------
! ****** File names.
!-----------------------------------------------------------------------
!
      use number_types
      use types
!
      implicit none
!
      type(vfile) :: bfile
!
      character(512) :: rffile=' ',tffile=' ',pffile=' ',effile=' '
      character(512) :: kffile=' ',qffile=' ',slogqffile=' '
      character(512) :: rbfile=' ',tbfile=' ',pbfile=' ',ebfile=' '
      character(512) :: kbfile=' ',qbfile=' ',slogqbfile=' '
      character(512) :: lffile=' ',lbfile=' '
!
! ****** File names for the r, t, and p meshes.
!
      character(512) :: mesh_file_r=' '
      character(512) :: mesh_file_t=' '
      character(512) :: mesh_file_p=' '
!
      type(vfile) :: volume3d_output_file
!
      type(vfile) :: slice_input_file
      type(vfile) :: slice_output_file_forward
      type(vfile) :: slice_output_file_backward
!
! ****** File name for the analytic function parameters.
!
      character(512) :: function_params_file&
     &                                  ='magfield_function_params.dat'
!
! ****** File names for slice output quantities.
!
      character(512) :: slice_q_output_file=' '
      character(512) :: slice_length_output_file=' '
!
! ****** File name for the output coronal hole map.
!
      character(512) :: ch_map_output_file='ch.h5'
!
! ****** File name for the output 3D coronal hole map.
!
      character(512) :: ch_map_3d_output_file='ch3d.h5'
!
! ****** File name for the output 3D dips map.
!
      character(512) :: dips_map_3d_output_file='dips3d.h5'
!
! ****** File name for the scalar field to be integrated along fl.
!
      character(512) :: scalar_input_file=' '
!
      end module
!#######################################################################
      module field_line_params
!
!-----------------------------------------------------------------------
! ****** Parameters that control the field line integration.
!-----------------------------------------------------------------------
!
      use number_types
!
      implicit none
!
!-----------------------------------------------------------------------
! ****** Parameters that control variable step-size tracing.
!-----------------------------------------------------------------------
!
! ****** These factors control how much the field line integration
! ****** step size can change from one step to another for the
! ****** case when a variable step size is being used.
!
! ****** MAX_INCREASE_FACTOR should be greater than 1, and
! ****** MAX_DECREASE_FACTOR should be less than 1.

      real(r_typ), parameter :: max_increase_factor=1.5_r_typ
      real(r_typ), parameter :: max_decrease_factor=.1_r_typ
!
!-----------------------------------------------------------------------
! ****** Parameters that control tracing of short field lines.
!-----------------------------------------------------------------------
!
! ****** If a field line trace has a smaller number of points
! ****** than SHORT_FL_MIN_POINTS, the integration step size
! ****** is decreased by the factor SHORT_FL_SHRINK_FACTOR,
! ****** and it is retraced, up to a maximum number of tries
! ****** equal to SHORT_FL_MAX_TRIES.
!
! ****** SHORT_FL_SHRINK_FACTOR should be less than 1.
!
      integer, parameter :: short_fl_min_points=10
      integer, parameter :: short_fl_max_tries=5
      real(r_typ), parameter :: short_fl_shrink_factor=.1_r_typ
!
!-----------------------------------------------------------------------
! ****** Parameters that control clipping to boundaries.
!-----------------------------------------------------------------------
!
! ****** The factor PREDICTOR_MIN_CLIP_FRACTION determines when to
! ****** clip a trace to the radial boundary in the predictor.
! ****** When the normalized distance to the r boundary (as a
! ****** fraction of the current step size) is less than
! ****** PREDICTOR_MIN_CLIP_FRACTION, the field line is clipped
! ****** to the boundary in the predictor without doing a
! ****** corrector step.  This number should be between 0 and 1.
!
      real(r_typ), parameter :: predictor_min_clip_fraction=.1_r_typ
!
!-----------------------------------------------------------------------
! ****** Maximum number of "bad" field line traces after
! ****** which to terminate.  Set to -1 to disable the termination.
!-----------------------------------------------------------------------
!
      integer, parameter :: max_bad_fieldlines=-1
!
      end module
!#######################################################################
      module step_size_stats
!
!-----------------------------------------------------------------------
! ****** Variable step size statistics.
!-----------------------------------------------------------------------
!
      use number_types
!
      implicit none
!
      logical :: gather_stats
!
      integer(8) :: stat_n=0
      real(r_typ) :: stat_ds_sum=0._r_typ
      real(r_typ) :: stat_ds_avg=0._r_typ
      real(r_typ) :: stat_ds_min=huge(0._r_typ)
      real(r_typ) :: stat_ds_max=0._r_typ
!
      end module
!#######################################################################
      module openmp_vars
!
!-----------------------------------------------------------------------
! ****** Variables to control OpenMP parallelization.
!-----------------------------------------------------------------------
!
      implicit none
!
! ****** Number of iterations to do in each thread.
!
      integer :: iterations_per_thread=500
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
      character(512) :: infile
      logical :: verbose
!
      end module
!#######################################################################
      module interp_interface
      interface
        subroutine interp (n,x,xv,i,ip1,alpha,tab)
        use number_types
        use invint_def
        use locate_interval_interface
        integer :: n
        real(r_typ), dimension(n) :: x
        real(r_typ) :: xv
        integer :: i
        integer :: ip1
        real(r_typ) :: alpha
        type(itab), optional :: tab
        end
      end interface
      end module
!#######################################################################
      module tracefl_interface
      interface
        subroutine tracefl (b,ds,s0,s1,bs0,bs1,s,&
     &                      traced_to_r_boundary,xt)
        use number_types
        use types
        type(vec) :: b
        type(flparam) :: ds
        real(r_typ), dimension(3) :: s0,s1
        real(r_typ), dimension(3) :: bs0,bs1
        real(r_typ) :: s
        logical :: traced_to_r_boundary
        type(traj), optional :: xt
        end
      end interface
      end module
!#######################################################################
      module integrate_fl
!
!-----------------------------------------------------------------------
! ****** Internal variables for integration along field line
!-----------------------------------------------------------------------
!
      use number_types
      use types
!
      implicit none
!
      logical :: do_integral_along_fl=.false.
      type(sds) :: scalar_field
      type(vtab) :: inv_sf

!
      end module
!#######################################################################
      program MAPFL
!
!-----------------------------------------------------------------------
!
      use ident
      use params
      use types
      use files
      use mesh
      use field
      use vars
      use field_line_params
      use step_size_stats
      use debug
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      integer :: ierr,i
      real(r_typ) :: ch_map_r=1.0_r_typ
      logical :: trace_fwd=.false.,trace_bwd=.false.
      logical :: compute_q_on_slice=.false.
      logical :: compute_ch_map=.false.
      logical :: compute_ch_map_3d=.false.
      logical :: compute_dips_map_3d=.false.
      character(256) :: errline=' '
!
!-----------------------------------------------------------------------
!
      namelist /datum/&
     &  debug_level, use_analytic_function, function_params_file,&
     &  domain_r_min, domain_r_max, bfile,&
     &  cubic, ds, set_ds_automatically,&
     &  dsmult, trace_fwd, trace_bwd,&
     &  rffile, tffile, pffile, effile, kffile, qffile, lffile,&
     &  rbfile, tbfile, pbfile, ebfile, kbfile, qbfile, lbfile,&
     &  new_r_mesh, mesh_file_r, nrss, r0,r1, new_t_mesh, mesh_file_t,&
     &  ntss, t0,t1, new_p_mesh, mesh_file_p, npss, p0,p1, trace_3d,&
     &  volume3d_output_file, trace_slice, slice_coords_are_xyz,&
     &  trace_slice_direction_is_along_b, compute_q_on_slice,&
     &  q_increment_h, slice_input_file, trace_from_slice_forward,&
     &  slice_output_file_forward, trace_from_slice_backward,&
     &  slice_output_file_backward, slice_q_output_file,&
     &  slice_length_output_file, compute_ch_map, ch_map_r,&
     &  ch_map_output_file, compute_ch_map_3d, ch_map_3d_output_file,&
     &  write_traces_to_hdf, write_traces_root, write_traces_as_xyz,&
     &  compute_dips_map_3d, dips_map_3d_output_file, ns_dips,&
     &  slogqffile,slogqbfile,integrate_along_fl,scalar_input_file
!
!-----------------------------------------------------------------------
!
! ****** Set the parameters.
!
      call set_parameters
!
      call ffopen (1,infile,'r',ierr)
!
      if (ierr.ne.0) then
        write (*,*)
        write (*,*) '### ERROR in MAPFL:'
        write (*,*) '### The input file does not exist'//&
     &              ' or cannot be read.'
        write (*,*) 'File name: ',trim(infile)
        call exit (1)
      end if
!
! ****** Read the input file.
!
      call ffopen (1,trim(infile),'r',ierr)
      read(1,datum,iostat=ierr)
      if (ierr.ne.0) then
        backspace (1)
        read (1,fmt='(A)') errline
        write (*,*)
        write (*,*) '### ERROR reading input file:'
        write (*,*) '### The following line has a problem:'
        write (*,*)
        write (*,*) trim(errline)
        write (*,*)
        write (*,*) '###'
        call exit (1)
      endif
      write (*,*)
      write (*,*) '### Input file contents:'
      write (*,*)
      write(*,datum)
      close (1)
!
      if (verbose) then
        write (*,*)
        write (*,*) '### ',cname,' Version ',cvers,' of ',cdate,'.'
      end if
!
! ****** Read the parameters that define the analytic magnetic
! ****** field function, if requested.
!
      if (use_analytic_function) then
        call read_function_params
      end if
!
! ****** Set the field line integration parameters.
!
      ds%max_increase_factor=max_increase_factor
      ds%max_decrease_factor=max_decrease_factor
      ds%predictor_min_clip_fraction=predictor_min_clip_fraction
      ds%short_fl_min_points=short_fl_min_points
      ds%short_fl_max_tries=short_fl_max_tries
      ds%short_fl_shrink_factor=short_fl_shrink_factor
!
! ****** Read the magnetic field.
!
      if (.not.use_analytic_function) then
        call readb (bfile,b)
      end if
!
! ****** Set the trace output format based on input br
! ****** (for analytic function, sets to hdf)
!
      i=index(bfile%r,'.h');
      if (bfile%r(i+1:i+2).eq.'h5') then
        fmt='h5'
      endif
!
! ****** Set the radial domain limits to those specified.
!
      b%lim0(1)=max(b%lim0(1),domain_r_min)
      b%lim1(1)=min(b%lim1(1),domain_r_max)
!
      if (verbose) then
        write (*,*)
        write (*,*) '### Domain limits:'
        write (*,*) 'Lower boundary value: ',b%lim0(1)
        write (*,*) 'Upper boundary value: ',b%lim1(1)
      end if
!
! ****** Make the new r, t, and p meshes.
!
      call make_new_meshes (b)
!
! ****** Set the default step size.
!
      call set_ds (b)
!
! ****** Set the flag to gather step size statistics.
!
      gather_stats=verbose
!
! ****** Setup the field to integrate along if requested.
!
      if (integrate_along_fl) call set_up_integration
!
! ****** Trace the field lines forward, if requested.
!

      if (trace_fwd) call map_forward
!
! ****** Trace the field lines backward, if requested.
!
      if (trace_bwd) call map_backward
!
! ****** Map the field lines from a 3D rectilinear volume,
! ****** if requested.
!
      if (trace_3d) call map_3d
!
! ****** Map the field lines from a slice, if requested,
! ****** or determine Q on the slice, if requested.
!
      if (trace_slice) then
        call read_slice_coordinates
        if (compute_q_on_slice) then
          call get_q_on_slice
        else
          call map_slice
        end if
        call deallocate_slice_coordinates
      end if
!
! ****** Compute a coronal hole map, if requested.
!
      if (compute_ch_map) then
        call get_ch_map (ch_map_r)
      end if
!
! ****** Compute a 3D coronal hole map, if requested.
!
      if (compute_ch_map_3d) then
        call get_ch_map_3d
      end if
!
! ****** Compute a 3D dips map, if requested.
!
      if (compute_dips_map_3d) then
        call get_dips_map_3d
      end if
!
      if (verbose) then
        stat_ds_avg=0.
        if (stat_n.ne.0) stat_ds_avg=stat_ds_sum/stat_n
        write (*,*)
        write (*,*) '### Field line integration step size statistics:'
        write (*,*) 'Number of field line segments = ',stat_n
        write (*,*) 'Minimum step size used = ',stat_ds_min
        write (*,*) 'Maximum step size used = ',stat_ds_max
        write (*,*) 'Average step size used = ',stat_ds_avg
      end if
!
      call exit (0)
!
      end
!#######################################################################
      subroutine read_function_params
!
!-----------------------------------------------------------------------
!
! ****** Read in the parameters for the analytic magnetic field
! ****** function.
!
!-----------------------------------------------------------------------
!
      use number_types
      use ident
      use params
      use vars
      use files
      use field
      use magfld_func_def
      use magfld_func_index
      use lcase_interface
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      real(r_typ), parameter :: pi=3.1415926535897932_r_typ
      real(r_typ), parameter :: twopi=6.2831853071795864_r_typ
!
!-----------------------------------------------------------------------
!
      if (verbose) then
        write (*,*)
        write (*,*) '### Using an analytic magnetic field'//&
     &              ' function.'
      end if
!
! ****** Read in the parameters (in NAMELIST format) for the
! ****** analytic function.
!
      if (verbose) then
        write (*,*) '### Reading the analytic field'//&
     &              ' function parameters ...'
        write (*,*) '### Parameters file name: ',&
     &              trim(function_params_file)
      end if
!
      call read_magfld_function_params (function_params_file)
!
! ****** Check that the requested function index is valid.
!
      if (function_index.lt.1.or.&
     &    function_index.gt.number_of_functions) then
        write (*,*)
        write (*,*) '### ERROR in ',cname,':'
        write (*,*) '### An invalid function index was requested:'
        write (*,*) 'FUNCTION_INDEX = ',function_index
        call exit (1)
      end if
!
      if (verbose) then
        write (*,*) '### Using the analytic magnetic field'//&
     &              ' function with index = ',function_index
      end if
!
! ****** Check that the selected geometry is consistent with the
! ****** use of an analytic field.
!
      if (.not.(new_r_mesh.and.new_t_mesh.and.new_p_mesh)) then
        write (*,*)
        write (*,*) '### ERROR in ',cname,':'
        write (*,*) '### Inconsistent parameters were specified when'//&
     &              ' requesting an analytic'
        write (*,*) '### magnetic field function:'
        write (*,*)
        write (*,*) 'You must specify the r, t, and p meshes to use.'
        call exit (1)
      end if
!
      if (ds%limit_by_local_mesh) then
        write (*,*)
        write (*,*) '### ERROR in ',cname,':'
        write (*,*) '### Inconsistent parameters were specified when'//&
     &              ' requesting an analytic'
        write (*,*) '### magnetic field function:'
        write (*,*)
        write (*,*) 'You must not attempt to limit the integration'//&
     &              ' step size by the local'
        write (*,*) 'magnetic field mesh.'
        call exit (1)
      end if
!
      if (set_ds_automatically) then
        write (*,*)
        write (*,*) '### ERROR in ',cname,':'
        write (*,*) '### Inconsistent parameters were specified when'//&
     &              ' requesting an analytic'
        write (*,*) '### magnetic field function:'
        write (*,*)
        write (*,*) 'You must not attempt to set the integration'//&
     &              ' step size automatically.'
        call exit (1)
      end if
!
! ****** Load the domain limits into the B structure.
!
      b%lim0(1)=domain_r_min
      b%lim1(1)=domain_r_max
      b%lim0(2)=0.
      b%lim1(2)=pi
      b%lim0(3)=0.
      b%lim1(3)=twopi
!
      return
      end
!#######################################################################
      subroutine read_magfld_function_params (fname)
!
!-----------------------------------------------------------------------
!
! ****** Read the NAMELIST parameters that define the magnetic
! ****** field analytic function.
!
!-----------------------------------------------------------------------
!
      use magfld_func_namelist
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
      integer :: ierr
!
!-----------------------------------------------------------------------
!
      call ffopen (1,fname,'r',ierr)
!
      if (ierr.ne.0) then
        write (*,*)
        write (*,*) '### ERROR in READ_MAGFLD_FUNCTION_PARAMS:'
        write (*,*) '### Could not open the analytic field'//&
     &              ' function parameters file.'
        write (*,*) 'File name: ',trim(fname)
        call exit (1)
      end if
!
! ****** Read the NAMELIST parameters.
!
      read (1,function_parameters)
      close (1)
!
      return
      end
!#######################################################################
      subroutine readb (bfile,b)
!
!-----------------------------------------------------------------------
!
! ****** Read the magnetic field from the files specified by
! ****** BFILE into the magnetic field vector structure B.
!
!-----------------------------------------------------------------------
!
      use number_types
      use types
      use vars
      use params
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      type(vfile) :: bfile
      type(vec) :: b
!
!-----------------------------------------------------------------------
!
      integer :: ierr
!
!-----------------------------------------------------------------------
!
! ****** Read the magnetic field components.
!
! ****** Br.
!
      if (verbose) then
        write (*,*)
        write (*,*) 'Reading data file: ',trim(bfile%r)
      end if
!
      call rdhdf (bfile%r,b%r,ierr)
!
      if (ierr.ne.0) then
        write (*,*)
        write (*,*) '### ERROR in READB:'
        write (*,*) '### Could not read Br.'
        write (*,*) 'IERR (from RDHDF) = ',ierr
        write (*,*) 'File name: ',trim(bfile%r)
        call exit (1)
      end if
!
      if (b%r%ndim.ne.3.or..not.b%r%scale) then
        write (*,*)
        write (*,*) '### ERROR in READB:'
        write (*,*) '### Invalid or missing scales in Br file.'
        write (*,*) 'File name: ',trim(bfile%r)
        call exit (1)
      end if
!
! ****** Bt.
!
      if (verbose) then
        write (*,*) 'Reading data file: ',trim(bfile%t)
      end if
!
      call rdhdf (bfile%t,b%t,ierr)
!
      if (ierr.ne.0) then
        write (*,*)
        write (*,*) '### ERROR in READB:'
        write (*,*) '### Could not read Bt.'
        write (*,*) 'IERR (from RDHDF) = ',ierr
        write (*,*) 'File name: ',trim(bfile%t)
        call exit (1)
      end if
!
      if (b%t%ndim.ne.3.or..not.b%t%scale) then
        write (*,*)
        write (*,*) '### ERROR in READB:'
        write (*,*) '### Invalid or missing scales in Bt file.'
        write (*,*) 'File name: ',trim(bfile%t)
        call exit (1)
      end if
!
! ****** Bp.
!
      if (verbose) then
        write (*,*) 'Reading data file: ',trim(bfile%p)
      end if
!
      call rdhdf (bfile%p,b%p,ierr)
!
      if (ierr.ne.0) then
        write (*,*)
        write (*,*) '### ERROR in READB:'
        write (*,*) '### Could not read Bp.'
        write (*,*) 'IERR (from RDHDF) = ',ierr
        write (*,*) 'File name: ',trim(bfile%p)
        call exit (1)
      end if
!
      if (b%p%ndim.ne.3.or..not.b%p%scale) then
        write (*,*)
        write (*,*) '### ERROR in READB:'
        write (*,*) '### Invalid or missing scales in Bp file.'
        write (*,*) 'File name: ',trim(bfile%p)
        call exit (1)
      end if
!
! ****** Set the type of magnetic field files.
      if (verbose) then
        write (*,*) 'Setting btype.'
      end if
!
      call set_btype (b)
!
! ****** Build the inverse interpolation tables.
!

      if (verbose) then
        write (*,*) 'Building inverse tables.'
      end if
!
      call build_inverse_tables (b%r,b%inv(1))
      call build_inverse_tables (b%t,b%inv(2))
      call build_inverse_tables (b%p,b%inv(3))
!
! ****** If cubic spline interpolation was requested, get the
! ****** spline coefficients.
!
      if (cubic) then
        b%cubic=.true.
        if (verbose) then
          write (*,*)
          write (*,*) 'Computing cubic spline coefficients'//&
     &                ' for Br ...'
        end if
        call compute_spline_3d (b%r%dims(1),b%r%dims(2),b%r%dims(3),&
     &                          b%r%scales(1)%f,&
     &                          b%r%scales(2)%f,&
     &                          b%r%scales(3)%f,&
     &                          b%r%f,b%spl%r)
        if (verbose) then
          write (*,*)
          write (*,*) 'Computing cubic spline coefficients'//&
     &                ' for Bt ...'
        end if
        call compute_spline_3d (b%t%dims(1),b%t%dims(2),b%t%dims(3),&
     &                          b%t%scales(1)%f,&
     &                          b%t%scales(2)%f,&
     &                          b%t%scales(3)%f,&
     &                          b%t%f,b%spl%t)
        if (verbose) then
          write (*,*)
          write (*,*) 'Computing cubic spline coefficients'//&
     &                ' for Bp ...'
        end if
        call compute_spline_3d (b%p%dims(1),b%p%dims(2),b%p%dims(3),&
     &                          b%p%scales(1)%f,&
     &                          b%p%scales(2)%f,&
     &                          b%p%scales(3)%f,&
     &                          b%p%f,b%spl%p)
      else
        b%cubic=.false.
      end if
!
      return
      end
!#######################################################################
      subroutine set_btype (b)
!
!-----------------------------------------------------------------------
!
! ****** Determine the primary (r,t,p) scales and the mesh limits
! ****** from the type of magnetic field in structure B, and
! ****** store them in structure B.
!
!-----------------------------------------------------------------------
!
      use number_types
      use types
      use constants
      use params
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      type(vec) :: b
!
!-----------------------------------------------------------------------
!
! ****** Tolerance for checking the bounds of t and p scales.
! ****** This should be set to several times the roundoff in
! ****** 32-bit representations of pi.
!
      real(r_typ), parameter :: eps=2.e-6_r_typ
!
! ****** The values of pi and 2*pi using 32-bit precision.
!
      real(KIND_REAL_4), parameter :: pi_r4=pi
      real(KIND_REAL_4), parameter :: twopi_r4=twopi
!
!-----------------------------------------------------------------------
!
      real(r_typ), dimension(:,:,:), pointer :: f
      real(r_typ), dimension(:), pointer :: z
      integer :: n1,n2,n3
      logical :: add_phi_point
!
!-----------------------------------------------------------------------
!
! ****** Check the type of magnetic field files read in.
!
      if (verbose) then
        write (*,*) 'SET_BTYPE: Checking type of magnetic field files.'
      end if
!
! ****** Check for new MAS code files.
!
      if (b%r%dims(1).eq.b%t%dims(1)+1.and.&
     &    b%r%dims(1).eq.b%p%dims(1)+1.and.&
     &    b%r%dims(2).eq.b%t%dims(2)-1.and.&
     &    b%r%dims(2).eq.b%p%dims(2)  .and.&
     &    b%r%dims(3).eq.b%t%dims(3)  .and.&
     &    b%r%dims(3).eq.b%p%dims(3)-1) then
!
        if (verbose) then
          write(*,*) 'SET_BTYPE: Detected new MAS code type.'
        end if
!
        b%nrs=b%r%dims(1)-1
        b%nts=b%r%dims(2)
        b%nps=b%r%dims(3)
!
        b%rs=>b%t%scales(1)%f
        b%ts=>b%r%scales(2)%f
        b%ps=>b%r%scales(3)%f
!
        add_phi_point=.false.
!
! ****** Check for old MAS code files.
!
      else if (b%r%dims(1).eq.b%t%dims(1)+1.and.&
     &         b%r%dims(1).eq.b%p%dims(1)+1.and.&
     &         b%r%dims(2).eq.b%t%dims(2)-1.and.&
     &         b%r%dims(2).eq.b%p%dims(2)  .and.&
     &         b%r%dims(3).eq.b%t%dims(3)  .and.&
     &         b%r%dims(3).eq.b%p%dims(3)  ) then
!
        if (verbose) then
          write(*,*) 'SET_BTYPE: Detected old MAS code type.'
        end if
!
        b%nrs=b%r%dims(1)-1
        b%nts=b%r%dims(2)
        b%nps=b%r%dims(3)
!
        b%rs=>b%t%scales(1)%f
        b%ts=>b%r%scales(2)%f
        b%ps=>b%r%scales(3)%f
!
! ****** Do not add a phi point if the phi interval already includes
! ****** the whole interval. Can occur if grid modified outside mapfl.
!
        if (abs((b%ps(b%nps)-b%ps(1))-twopi).lt.eps) then
          if (verbose) then
            write (*,*) 'SET_BTYPE: Phi already wrapped!'
          end if
          add_phi_point=.false.
        else
          add_phi_point=.true.
        end if
!
! ****** Check for new POT3D code files.
!
      else if (b%r%dims(1).eq.b%t%dims(1)-1.and.&
     &         b%r%dims(1).eq.b%p%dims(1)-1.and.&
     &         b%r%dims(2).eq.b%t%dims(2)+1.and.&
     &         b%r%dims(2).eq.b%p%dims(2)  .and.&
     &         b%r%dims(3).eq.b%t%dims(3)  .and.&
     &         b%r%dims(3).eq.b%p%dims(3)+1) then
!
        if (verbose) then
          write(*,*) 'SET_BTYPE: Detected POT3D code type.'
        end if
!
        b%nrs=b%r%dims(1)
        b%nts=b%r%dims(2)-1
        b%nps=b%r%dims(3)-1
!
        b%rs=>b%r%scales(1)%f
        b%ts=>b%t%scales(2)%f
        b%ps=>b%p%scales(3)%f
!
        add_phi_point=.false.
!
! ****** Check for new-old POT3D code files.
!
      else if (b%r%dims(1).eq.b%t%dims(1)-1.and.&
     &         b%r%dims(1).eq.b%p%dims(1)-1.and.&
     &         b%r%dims(2).eq.b%t%dims(2)+1.and.&
     &         b%r%dims(2).eq.b%p%dims(2)  .and.&
     &         b%r%dims(3).eq.b%t%dims(3)  .and.&
     &         b%r%dims(3).eq.b%p%dims(3)-1) then
!
        if (verbose) then
          write(*,*) 'SET_BTYPE: Detected new-old POT3D type.'
        end if
!
        b%nrs=b%r%dims(1)
        b%nts=b%r%dims(2)-1
        b%nps=b%r%dims(3)
!
        b%rs=>b%r%scales(1)%f
        b%ts=>b%t%scales(2)%f
        b%ps=>b%r%scales(3)%f
!
        add_phi_point=.false.
!
! ****** Check for old POT3D code files.
!
      else if (b%r%dims(1).eq.b%t%dims(1)-1.and.&
     &         b%r%dims(1).eq.b%p%dims(1)-1.and.&
     &         b%r%dims(2).eq.b%t%dims(2)+1.and.&
     &         b%r%dims(2).eq.b%p%dims(2)  .and.&
     &         b%r%dims(3).eq.b%t%dims(3)  .and.&
     &         b%r%dims(3).eq.b%p%dims(3)  ) then
!
        if (verbose) then
          write(*,*) 'SET_BTYPE: Detected old POT3D code type.'
        end if
!
        b%nrs=b%r%dims(1)
        b%nts=b%r%dims(2)-1
        b%nps=b%r%dims(3)
!
        b%rs=>b%r%scales(1)%f
        b%ts=>b%t%scales(2)%f
        b%ps=>b%r%scales(3)%f
!
! ****** Do not add a phi point if the phi interval already includes
! ****** the whole interval. Can occur if grid modified outside mapfl.
!
        if (abs((b%ps(b%nps)-b%ps(1))-twopi).lt.eps) then
          if (verbose) then
            write (*,*) 'SET_BTYPE: Phi already wrapped!'
          end if
          add_phi_point=.false.
        else
          add_phi_point=.true.
        end if
!
! ****** Check for non-staggered files.
!
      else if (b%r%dims(1).eq.b%t%dims(1).and.&
     &         b%r%dims(1).eq.b%p%dims(1).and.&
     &         b%r%dims(2).eq.b%t%dims(2).and.&
     &         b%r%dims(2).eq.b%p%dims(2).and.&
     &         b%r%dims(3).eq.b%t%dims(3).and.&
     &         b%r%dims(3).eq.b%p%dims(3)) then
!
        if (verbose) then
          write(*,*) 'SET_BTYPE: Detected non-staggered type.'
        end if
!
        b%nrs=b%r%dims(1)
        b%nts=b%r%dims(2)
        b%nps=b%r%dims(3)
!
        b%rs=>b%r%scales(1)%f
        b%ts=>b%r%scales(2)%f
        b%ps=>b%r%scales(3)%f
!
! ****** Do not add a phi point if the phi interval already includes
! ****** the whole interval. Can occur if grid modified outside mapfl.
!
        if (abs((b%ps(b%nps)-b%ps(1))-twopi).lt.eps) then
          if (verbose) then
            write (*,*) 'SET_BTYPE: Phi already wrapped!'
          end if
          add_phi_point=.false.
        else
          add_phi_point=.true.
        end if
!
      else
!
! ****** Invalid file type.
!
        write (*,*)
        write (*,*) '### ERROR in SET_BTYPE:'
        write (*,*) '### Unrecognized magnetic field file staggering.'
        write (*,*) '  br resolution:', b%r%dims
        write (*,*) '  bt resolution:', b%t%dims
        write (*,*) '  bp resolution:', b%p%dims
        call exit (1)
!
      end if
!
! ****** If appropriate, add a point in the phi dimension to
! ****** take care of periodic wrap-around.
!
      if (add_phi_point) then
!
        if (verbose) then
          write (*,*) 'SET_BTYPE: Adding phi point.'
        end if
!
        n1=b%r%dims(1)
        n2=b%r%dims(2)
        n3=b%r%dims(3)
        allocate (f(n1,n2,n3+1))
        allocate (z(n3+1))
        f(:,:,1:n3)=b%r%f(:,:,:)
        f(:,:,n3+1)=b%r%f(:,:,1)
        z(1:n3)=b%r%scales(3)%f(:)
        z(n3+1)=b%r%scales(3)%f(1)+twopi
        deallocate (b%r%f)
        deallocate (b%r%scales(3)%f)
        b%r%dims(3)=n3+1
        b%r%f=>f
        b%r%scales(3)%f=>z
!
        n1=b%t%dims(1)
        n2=b%t%dims(2)
        n3=b%t%dims(3)
        allocate (f(n1,n2,n3+1))
        allocate (z(n3+1))
        f(:,:,1:n3)=b%t%f(:,:,:)
        f(:,:,n3+1)=b%t%f(:,:,1)
        z(1:n3)=b%t%scales(3)%f(:)
        z(n3+1)=b%t%scales(3)%f(1)+twopi
        deallocate (b%t%f)
        deallocate (b%t%scales(3)%f)
        b%t%dims(3)=n3+1
        b%t%f=>f
        b%t%scales(3)%f=>z
!
        n1=b%p%dims(1)
        n2=b%p%dims(2)
        n3=b%p%dims(3)
        allocate (f(n1,n2,n3+1))
        allocate (z(n3+1))
        f(:,:,1:n3)=b%p%f(:,:,:)
        f(:,:,n3+1)=b%p%f(:,:,1)
        z(1:n3)=b%p%scales(3)%f(:)
        z(n3+1)=b%p%scales(3)%f(1)+twopi
        deallocate (b%p%f)
        deallocate (b%p%scales(3)%f)
        b%p%dims(3)=n3+1
        b%p%f=>f
        b%p%scales(3)%f=>z
!
        b%nps=b%r%dims(3)
        b%ps=>b%r%scales(3)%f
!
      end if
!
! ****** Set the precision of the B that was read in, based
! ****** on the type of the individual HDF files of the components.
!
      if (verbose) then
        write (*,*) 'SET_BTYPE: Setting B precision.'
      end if
!
      if (b%r%hdf32.or.b%t%hdf32.or.b%p%hdf32) then
        b%b_is_32bit=.true.
      else
        b%b_is_32bit=.false.
      end if
!
! ****** Snap the outer (t,p) limits to values that are slightly
! ****** larger than the exact values, if they are close enough.
! ****** This will compensate for the reduced precision inherent in
! ****** magnetic fields read in from 32-bit HDF files, which have
! ****** only ~ 7 digits of accuracy.  In this way, positions
! ****** read from 32-bit HDF files that are near the theta=pi
! ****** and phi=2*pi boundaries will be more likely to end up
! ****** inside the domain.
!
      if (abs(b%ts(b%nts)-pi).lt.eps) then
        if (verbose) then
          write (*,*) 'SET_BTYPE: Snapping t.'
        end if
        b%ts(b%nts)=pi+3._r_typ*spacing(pi_r4)
      end if
!
      if (abs(b%ps(b%nps)-twopi).lt.eps) then
        if (verbose) then
          write (*,*) 'SET_BTYPE: Snapping p.'
        end if
        b%ps(b%nps)=twopi+3._r_typ*spacing(twopi_r4)
      end if
!
! ****** Set the domain limits.
!
      if (verbose) then
        write (*,*) 'SET_BTYPE: Set domain limits.'
      end if
!
      b%lim0(1)=b%rs(1)
      b%lim1(1)=b%rs(b%nrs)
      b%lim0(2)=b%ts(1)
      b%lim1(2)=b%ts(b%nts)
      b%lim0(3)=b%ps(1)
      b%lim1(3)=b%ps(b%nps)
!
! ****** Build the inverse interpolation tables for the
! ****** main mesh.
!
      if (verbose) then
        write (*,*) 'SET_BTYPE: Building inverse interpolation tables.'
      end if
!
      b%rs_invtab%n=b%nrs
      allocate (b%rs_invtab%f(b%rs_invtab%n))
      call getinv (b%rs,b%nrs,b%rs_invtab)
!
      b%ts_invtab%n=b%nts
      allocate (b%ts_invtab%f(b%ts_invtab%n))
      call getinv (b%ts,b%nts,b%ts_invtab)
!
      b%ps_invtab%n=b%nps
      allocate (b%ps_invtab%f(b%ps_invtab%n))
      call getinv (b%ps,b%nps,b%ps_invtab)
!
! ****** Compute the mesh cell dimensions on the main mesh.
!
      if (verbose) then
        write (*,*) 'SET_BTYPE: Computing mesh cell dims.'
      end if
!
! ****** These are used in setting the field line integration
! ****** step size.
!
      allocate (b%drs(b%nrs))
      allocate (b%dts(b%nts))
      allocate (b%dps(b%nps))
      allocate (b%sts(b%nts))
!
      call get_dx (b%nrs,b%rs,b%drs)
      call get_dx (b%nts,b%ts,b%dts)
      call get_dx (b%nps,b%ps,b%dps)
      b%sts=sin(b%ts)
      b%sts(    1)=max(b%sts(    1),sin(b%dts(    1)))
      b%sts(b%nts)=max(b%sts(b%nts),sin(b%dts(b%nts)))
!
      return
      end
!#######################################################################
      subroutine get_dx (n,x,dx)
!
!-----------------------------------------------------------------------
!
! ****** Get the cell size DX(N) from the 1D mesh in X(N).
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
      real(r_typ), dimension(n) :: x,dx
!
!-----------------------------------------------------------------------
!
      real(r_typ), parameter :: half=.5_r_typ
!
!-----------------------------------------------------------------------
!
      integer :: i
!
!-----------------------------------------------------------------------
!
      if (n.le.1) then
        dx(1)=0.
      else if (n.eq.2) then
        dx(1)=x(2)-x(1)
        dx(2)=dx(1)
      else
        do i=2,n-1
          dx(i)=half*(x(i+1)-x(i-1))
        enddo
        dx(1)=dx(2)
        dx(n)=dx(n-1)
      end if
!
      return
      end
!#######################################################################
      subroutine make_new_meshes (b)
!
!-----------------------------------------------------------------------
!
! ****** Make new r, t, and p meshes, if requested, or link them
! ****** to the meshes in the B files.
!
!-----------------------------------------------------------------------
!
      use number_types
      use types
      use mesh
      use params
      use files
      use vars
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      type(vec) :: b
!
!-----------------------------------------------------------------------
!
      type(sds) :: s
      integer :: i,ierr
      real(r_typ) :: d
!
!-----------------------------------------------------------------------
!
! ****** Make the r mesh.
!
      if (new_r_mesh) then
!
! ****** Check if the mesh is to be read from a 1D HDF file.
!
        if (mesh_file_r.ne.' ') then
!
          if (verbose) then
            write (*,*)
            write (*,*) '### Reading the r mesh from file: ',&
     &                  trim(mesh_file_r)
          end if
!
          call rdhdf (mesh_file_r,s,ierr)
!
          if (ierr.ne.0) then
            write (*,*)
            write (*,*) '### ERROR in MAKE_NEW_MESHES:'
            write (*,*) '### Error while reading the r mesh'//&
     &                  ' from a file.'
            write (*,*) '### Could not read the data set.'
            write (*,*) 'File name: ',trim(mesh_file_r)
            call exit (1)
          end if
!
          if (s%ndim.ne.1) then
            write (*,*)
            write (*,*) '### ERROR in MAKE_NEW_MESHES:'
            write (*,*) '### Error while reading the r mesh'//&
     &                  ' from a file.'
            write (*,*) '### The HDF file does not contain a 1D'//&
     &                  ' data set.'
            write (*,*) 'File name: ',trim(mesh_file_r)
            call exit (1)
          end if
!
          nrss=s%dims(1)
          allocate (rss(nrss))
          rss=s%f(:,1,1)
!
          call deallocate_sds (s)
!
          if (verbose) then
            write (*,*)
            write (*,*) '### Mesh read in for the r mesh:'
            write (*,*) 'Number of points = ',nrss
            write (*,*) 'Lower limit = ',rss(1)
            write (*,*) 'Upper limit = ',rss(nrss)
          end if
!
        else
!
! ****** Generate a uniform mesh.
!
          if (nrss.lt.1) then
            write (*,*)
            write (*,*) '### ERROR in MAKE_NEW_MESHES:'
            write (*,*) '### Invalid number of points specified'//&
     &                  ' for the uniform r mesh.'
            write (*,*) 'Number of points specified = ',nrss
            call exit (1)
          end if
!
          allocate (rss(nrss))
!
          if (r0.eq.0..and.r1.eq.0.) then
            r0=b%lim0(1)
            r1=b%lim1(1)
          end if
!
          if (verbose) then
            write (*,*)
            write (*,*) '### Generating a uniform r mesh:'
            write (*,*) 'Number of points = ',nrss
            write (*,*) 'Lower limit = ',r0
            write (*,*) 'Upper limit = ',r1
          end if
!
          if (nrss.ne.1) then
            d=(r1-r0)/(nrss-1)
          else
            d=0.
          end if
!
          do i=1,nrss
            rss(i)=r0+(i-1)*d
          enddo
          rss(1)=r0
          if (nrss.gt.1) then
            rss(nrss)=r1
          end if
!
        end if
!
      else
!
! ****** Use the same mesh as the primary B field mesh.
!
        nrss=b%nrs
        rss=>b%rs
!
      end if
!
! ****** Make the t mesh.
!
      if (new_t_mesh) then
!
! ****** Check if the mesh is to be read from a 1D HDF file.
!
        if (mesh_file_t.ne.' ') then
!
          if (verbose) then
            write (*,*)
            write (*,*) '### Reading the t mesh from file: ',&
     &                  trim(mesh_file_t)
          end if
!
          call rdhdf (mesh_file_t,s,ierr)
!
          if (ierr.ne.0) then
            write (*,*)
            write (*,*) '### ERROR in MAKE_NEW_MESHES:'
            write (*,*) '### Error while reading the t mesh'//&
     &                  ' from a file.'
            write (*,*) '### Could not read the data set.'
            write (*,*) 'File name: ',trim(mesh_file_t)
            call exit (1)
          end if
!
          if (s%ndim.ne.1) then
            write (*,*)
            write (*,*) '### ERROR in MAKE_NEW_MESHES:'
            write (*,*) '### Error while reading the t mesh'//&
     &                  ' from a file.'
            write (*,*) '### The HDF file does not contain a 1D'//&
     &                  ' data set.'
            write (*,*) 'File name: ',trim(mesh_file_t)
            call exit (1)
          end if
!
          ntss=s%dims(1)
          allocate (tss(ntss))
          tss=s%f(:,1,1)
!
          call deallocate_sds (s)
!
          if (verbose) then
            write (*,*)
            write (*,*) '### Mesh read in for the t mesh:'
            write (*,*) 'Number of points = ',ntss
            write (*,*) 'Lower limit = ',tss(1)
            write (*,*) 'Upper limit = ',tss(ntss)
          end if
!
        else
!
! ****** Generate a uniform mesh.
!
          if (ntss.lt.1) then
            write (*,*)
            write (*,*) '### ERROR in MAKE_NEW_MESHES:'
            write (*,*) '### Invalid number of points specified'//&
     &                  ' for the uniform t mesh.'
            write (*,*) 'Number of points specified = ',ntss
            call exit (1)
          end if
!
          allocate (tss(ntss))
!
          if (t0.eq.0..and.t1.eq.0.) then
            t0=b%lim0(2)
            t1=b%lim1(2)
          end if
!
          if (verbose) then
            write (*,*)
            write (*,*) '### Generating a uniform t mesh:'
            write (*,*) 'Number of points = ',ntss
            write (*,*) 'Lower limit = ',t0
            write (*,*) 'Upper limit = ',t1
          end if
!
          if (ntss.ne.1) then
            d=(t1-t0)/(ntss-1)
          else
            d=0.
          end if
!
          do i=1,ntss
            tss(i)=t0+(i-1)*d
          enddo
          tss(1)=t0
          if (ntss.gt.1) then
            tss(ntss)=t1
          end if
!
        end if
!
      else
!
! ****** Use the same mesh as the primary B field mesh.
!
        ntss=b%nts
        tss=>b%ts
!
      end if
!
! ****** Make the p mesh.
!
      if (new_p_mesh) then
!
! ****** Check if the mesh is to be read from a 1D HDF file.
!
        if (mesh_file_p.ne.' ') then
!
          if (verbose) then
            write (*,*)
            write (*,*) '### Reading the p mesh from file: ',&
     &                  trim(mesh_file_p)
          end if
!
          call rdhdf (mesh_file_p,s,ierr)
!
          if (ierr.ne.0) then
            write (*,*)
            write (*,*) '### ERROR in MAKE_NEW_MESHES:'
            write (*,*) '### Error while reading the p mesh'//&
     &                  ' from a file.'
            write (*,*) '### Could not read the data set.'
            write (*,*) 'File name: ',trim(mesh_file_p)
            call exit (1)
          end if
!
          if (s%ndim.ne.1) then
            write (*,*)
            write (*,*) '### ERROR in MAKE_NEW_MESHES:'
            write (*,*) '### Error while reading the p mesh'//&
     &                  ' from a file.'
            write (*,*) '### The HDF file does not contain a 1D'//&
     &                  ' data set.'
            write (*,*) 'File name: ',trim(mesh_file_p)
            call exit (1)
          end if
!
          npss=s%dims(1)
          allocate (pss(npss))
          pss=s%f(:,1,1)
!
          call deallocate_sds (s)
!
          if (verbose) then
            write (*,*)
            write (*,*) '### Mesh read in for the p mesh:'
            write (*,*) 'Number of points = ',npss
            write (*,*) 'Lower limit = ',pss(1)
            write (*,*) 'Upper limit = ',pss(npss)
          end if
!
        else
!
! ****** Generate a uniform mesh.
!
          if (npss.lt.1) then
            write (*,*)
            write (*,*) '### ERROR in MAKE_NEW_MESHES:'
            write (*,*) '### Invalid number of points specified'//&
     &                  ' for the uniform p mesh.'
            write (*,*) 'Number of points specified = ',npss
            call exit (1)
          end if
!
          allocate (pss(npss))
!
          if (p0.eq.0..and.p1.eq.0.) then
            p0=b%lim0(3)
            p1=b%lim1(3)
          end if
!
          if (verbose) then
            write (*,*)
            write (*,*) '### Generating a uniform p mesh:'
            write (*,*) 'Number of points = ',npss
            write (*,*) 'Lower limit = ',p0
            write (*,*) 'Upper limit = ',p1
          end if
!
          if (npss.ne.1) then
            d=(p1-p0)/(npss-1)
          else
            d=0.
          end if
!
          do i=1,npss
            pss(i)=p0+(i-1)*d
          enddo
          pss(1)=p0
          if (npss.gt.1) then
            pss(npss)=p1
          end if
!
        end if
!
      else
!
! ****** Use the same mesh as the primary B field mesh.
!
        npss=b%nps
        pss=>b%ps
!
      end if
!
      return
      end
!#######################################################################
      subroutine set_ds (b)
!
!-----------------------------------------------------------------------
!
! ****** Set the field line integration step size.
!
!-----------------------------------------------------------------------
!
! ****** If SET_DS_AUTOMATICALLY=.T., the miniumum step size is set
! ****** to the minimum of the cell dimensions from the magnetic
! ****** field files, and the maximum step size is set to the
! ****** maximum of the cell dimensions.  Otherwise, the values read
! ****** in for DS%MIN and DS%MAX are used.
!
! ****** After being set in the above way, DS%MIN and DS%MAX are
! ****** multiplied by the factor DSMULT.  Thus, DSMULT provides a
! ****** quick way to change the integration step size.
!
!-----------------------------------------------------------------------
!
      use number_types
      use types
      use constants
      use vars
      use params
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      type(vec) :: b
!
!-----------------------------------------------------------------------
!
      integer :: i,j,k
      real(r_typ) :: dr,dt,dp
      real(r_typ) :: drmin,dtmin,dpmin
      real(r_typ) :: drmax,dtmax,dpmax
!
!-----------------------------------------------------------------------
!
      if (set_ds_automatically) then
!
        drmin=abs(b%lim1(1)-b%lim0(1))
        drmax=0.
        do i=1,b%nrs-1
          dr=abs(b%rs(i+1)-b%rs(i))
          drmin=min(drmin,dr)
          drmax=max(drmax,dr)
        enddo
!
        dtmin=pi
        dtmax=0.
        do j=1,b%nts-1
          dt=abs(b%ts(j+1)-b%ts(j))
          dtmin=min(dtmin,dt)
          dtmax=max(dtmax,dt)
        enddo
!
        dpmin=twopi
        dpmax=0.
        do k=1,b%nps-1
          dp=abs(b%ps(k+1)-b%ps(k))
          dpmin=min(dpmin,dp)
          dpmax=max(dpmax,dp)
        enddo
!
        ds%min=min(drmin,b%lim0(1)*dtmin,b%lim0(1)*dtmin*dpmin)
        ds%max=max(drmax,b%lim1(1)*dtmax,b%lim1(1)*dpmax)
!
      end if
!
      if (dsmult.le.0.) then
        write (*,*)
        write (*,*) '### ERROR in SET_DS:'
        write (*,*) '### DSMULT must be positive.'
        write (*,*) 'DSMULT= ',dsmult
        call exit (1)
      end if
!
      ds%over_rc=ds%over_rc*dsmult
      ds%min=ds%min*dsmult
      ds%max=ds%max*dsmult
!
      if (verbose) then
        write (*,*)
        write (*,*) '### Field line integration parameters:'
        if (ds%variable) then
          write (*,*)
          write (*,*) '### Integration step size control:'//&
     &                ' variable step size'
          write (*,*)
          write (*,*) '### Step size parameters:'
          write (*,*) 'DS%OVER_RC = ',ds%over_rc
          write (*,*) 'DS%MIN = ',ds%min
          write (*,*) 'DS%MAX = ',ds%max
          write (*,*)
          if (ds%limit_by_local_mesh) then
            write (*,*) '### Step size limited by the local'//&
     &                  ' B mesh: yes'
            write (*,*) 'DS%LOCAL_MESH_FACTOR = ',&
     &                  ds%local_mesh_factor
          else
            write (*,*) '### Step size limited by the local'//&
     &                  ' B mesh: no'
          end if
        else
          write (*,*)
          write (*,*) '### Integration step size control:'//&
     &                ' uniform step size'
          write (*,*)
          write (*,*) '### Step size parameters:'
          write (*,*) 'DS = ',ds%min
        end if
      end if
!
      return
      end
!#######################################################################
      subroutine map_forward
!
!-----------------------------------------------------------------------
!
! ****** Trace field lines outward from r=R0.
!
!-----------------------------------------------------------------------
!
      use number_types
      use types
      use mesh
      use field
      use vars
      use files
      use params
      use field_line_params
      use diags
      use openmp_vars
      use tracefl_interface
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      real(r_typ), parameter :: half=.5_r_typ
      real(r_typ), parameter :: quarter=.25_r_typ
!
!-----------------------------------------------------------------------
!
! ****** Storage for the mapping.
!
      real(r_typ), dimension(ntss,npss) :: rfl,tfl,pfl,efl,kfl,length
!
      real(r_typ), dimension(:), allocatable :: tssh,pssh
      real(r_typ), dimension(:,:), allocatable :: qfl,slogqfl
!
!-----------------------------------------------------------------------
!
      integer :: ierr,j,k,nbad
      real(r_typ), dimension(3) :: xfl0,xfl1
      real(r_typ), dimension(3) :: bs0,bs1
      logical :: ttb
      real(r_typ) :: s
      real(r_typ) :: dtdt_m,dtdt_p
      real(r_typ) :: dtdp_m,dtdp_p
      real(r_typ) :: dpdt_m,dpdt_p
      real(r_typ) :: dpdp_m,dpdp_p
      real(r_typ) :: dtdt,dtdp,dpdt,dpdp
      real(r_typ) :: dt,dp,aa,bb,cc,dd,stm,stp,tmav,efav
      logical :: wrote_cr
      integer :: n_completed,n_total,nc,diag_step
      real(r_typ) :: pct_done
!
!-----------------------------------------------------------------------
!
      real(r_typ), external :: modulo_twopi
!
!-----------------------------------------------------------------------
!
! ****** Trace field lines, starting from each (T,P) cell at r=R0,
! ****** until the field line hits r=R1, or goes back to r=R0,
! ****** or exhausts the field line length allowed.
!
      if (verbose) then
        write (*,*)
        write (*,*) '### Computing a forward mapping from R0:'
      end if
!
      if (verbose) then
        write (*,*)
        write (*,*) 'Mapping field lines ...'
        write (*,*)
      end if
!
      nbad=0
!
      ds%direction_is_along_b=.false.
      ds%direction=1
!
      n_total=ntss*npss
      n_completed=0
!
!$omp parallel do
!$omp& private(j,k,xfl0,xfl1,bs0,bs1,s,ttb)
!$omp& private(nc,diag_step,pct_done)
!$omp& collapse(2)
!$omp& schedule(dynamic,iterations_per_thread)
      do k=1,npss
        do j=1,ntss
!
! ****** Update the iteration counter for diagnostic
! ****** purposes.
!
          if (verbose) then
!$omp critical
            n_completed=n_completed+1
            nc=n_completed
!$omp end critical
          end if
!
          xfl0(1)=b%lim0(1)
          xfl0(2)=tss(j)
          xfl0(3)=pss(k)
!
          call tracefl (b,ds,xfl0,xfl1,bs0,bs1,s,ttb)
!
! ****** Check that the field line reached R0 or R1, and set
! ****** the expansion factor.
!
          if (ttb) then
            rfl(j,k)=xfl1(1)
            tfl(j,k)=xfl1(2)
            pfl(j,k)=xfl1(3)
            length(j,k)=s
            if (bs1(1).ne.0.) then
              efl(j,k)=abs((bs0(1)*xfl0(1)**2)/(bs1(1)*xfl1(1)**2))
            else
              efl(j,k)=0.
            end if
            if (bs1(1).ne.0.) then
              kfl(j,k)=log10(max(abs(bs0(1)/bs1(1)),tiny(bs0(1))))
            else
              kfl(j,k)=-50._r_typ
            end if
          else
!$omp critical
            nbad=nbad+1
            write (*,*)
            write (*,*) '### WARNING from MAP_FORWARD:'
            write (*,*) '### A field line did not reach R0 or R1.'
            write (*,*) 'Initial theta = ',xfl0(2)
            write (*,*) 'Initial phi   = ',xfl0(3)
            write (*,*) 'Final field line radius = ',xfl1(1)
!$omp end critical
            rfl(j,k)=-1._r_typ
            tfl(j,k)=-1._r_typ
            pfl(j,k)=-1._r_typ
            length(j,k)=0.
            efl(j,k)=0.
            kfl(j,k)=-50._r_typ
          end if
!
          if (max_bad_fieldlines.gt.0) then
            if (nbad.gt.max_bad_fieldlines) then
!$omp critical
              write (*,*)
              write (*,*) '### ERROR in MAP_FORWARD:'
              write (*,*) '### Too many field lines did not reach'//&
     &                    ' R0 or R1.'
              write (*,*) 'Number of bad traces = ',max_bad_fieldlines
              call exit (1)
!$omp end critical
            end if
          end if
!
! ****** Write progress diagnostics if requested.
!
          if (verbose) then
            diag_step=mod(nc,diagnostic_interval)
            if (diag_step.eq.0) then
              pct_done=100.*nc/n_total
              write (*,910) 'Fraction completed: ',pct_done
  910         format (1x,a,f7.3,'%')
            end if
          end if
!
        enddo
      enddo
!$omp end parallel do
!
! ****** Write the mapping.
!
      wrote_cr=.false.
!
      if (rffile.ne.' ') then
        if (verbose) then
          if (.not.wrote_cr) then
            write (*,*)
            wrote_cr=.true.
          end if
          write (*,*) 'Writing the forward mapping for coordinate '//&
     &                'r to file: ',&
     &                trim(rffile)
        end if
        call wrhdf_2d (rffile,.true.,ntss,npss,rfl,tss,pss,&
     &                 hdf32,ierr)
        if (ierr.ne.0) then
          write (*,*)
          write (*,*) '### ERROR in MAP_FORWARD:'
          write (*,*) '### Could not write the forward mapping'//&
     &                ' file for coordinate r.'
          call exit (1)
        end if
      end if
!
      if (tffile.ne.' ') then
        if (verbose) then
          if (.not.wrote_cr) then
            write (*,*)
            wrote_cr=.true.
          end if
          write (*,*) 'Writing the forward mapping for coordinate '//&
     &                't to file: ',&
     &                trim(tffile)
        end if
        call wrhdf_2d (tffile,.true.,ntss,npss,tfl,tss,pss,&
     &                 hdf32,ierr)
        if (ierr.ne.0) then
          write (*,*)
          write (*,*) '### ERROR in MAP_FORWARD:'
          write (*,*) '### Could not write the forward mapping'//&
     &                ' file for coordinate t.'
          call exit (1)
        end if
      end if
!
      if (pffile.ne.' ') then
        if (verbose) then
          if (.not.wrote_cr) then
            write (*,*)
            wrote_cr=.true.
          end if
          write (*,*) 'Writing the forward mapping for coordinate '//&
     &                'p to file: ',&
     &                trim(pffile)
        end if
        call wrhdf_2d (pffile,.true.,ntss,npss,pfl,tss,pss,&
     &                 hdf32,ierr)
        if (ierr.ne.0) then
          write (*,*)
          write (*,*) '### ERROR in MAP_FORWARD:'
          write (*,*) '### Could not write the forward mapping'//&
     &                ' file for coordinate p.'
          call exit (1)
        end if
      end if
!
      if (effile.ne.' ') then
        if (verbose) then
          if (.not.wrote_cr) then
            write (*,*)
            wrote_cr=.true.
          end if
          write (*,*) 'Writing the forward mapping '//&
     &                'expansion factor to file: ',&
     &                trim(effile)
        end if
        call wrhdf_2d (effile,.true.,ntss,npss,efl,tss,pss,&
     &                 hdf32,ierr)
        if (ierr.ne.0) then
          write (*,*)
          write (*,*) '### ERROR in MAP_FORWARD:'
          write (*,*) '### Could not write the forward mapping'//&
     &                ' file for the expansion factor.'
          call exit (1)
        end if
      end if
!
      if (kffile.ne.' ') then
        if (verbose) then
          if (.not.wrote_cr) then
            write (*,*)
            wrote_cr=.true.
          end if
          write (*,*) 'Writing the forward mapping '//&
     &                'K factor to file: ',&
     &                trim(kffile)
        end if
        call wrhdf_2d (kffile,.true.,ntss,npss,kfl,tss,pss,&
     &                 hdf32,ierr)
        if (ierr.ne.0) then
          write (*,*)
          write (*,*) '### ERROR in MAP_FORWARD:'
          write (*,*) '### Could not write the forward mapping'//&
     &                ' file for the K factor.'
          call exit (1)
        end if
      end if
!
      if (lffile.ne.' ') then
        if (verbose) then
          if (.not.wrote_cr) then
            write (*,*)
            wrote_cr=.true.
          end if
          write (*,*) 'Writing the forward mapping '//&
     &                'length to file: ',&
     &                trim(lffile)
        end if
        call wrhdf_2d (lffile,.true.,ntss,npss,length,tss,pss,&
     &                 hdf32,ierr)
        if (ierr.ne.0) then
          write (*,*)
          write (*,*) '### ERROR in MAP_FORWARD:'
          write (*,*) '### Could not write the forward mapping'//&
     &                'length file.'
          call exit (1)
        end if
      end if
!
! ****** Compute Q (if requested).
!
      if (qffile.eq.' '.and.slogqffile.eq.' ') return
!
! ****** This can only be done if NTSS and NPSS exceed 1.
!
      if (ntss.le.1.or.npss.le.1) then
        write (*,*)
        write (*,*) '### WARNING from MAP_FORWARD:'
        write (*,*) '### Could not compute the Q factor.'
        write (*,*) '### To compute Q, NTSS and NPSS'//&
     &              ' must be greater than 1.'
        return
      end if
!
      allocate (tssh(ntss-1))
      allocate (pssh(npss-1))
      allocate (qfl(ntss-1,npss-1))
      if (slogqffile.ne.' ') then
        allocate (slogqfl(ntss-1,npss-1))
      end if
!
! ****** Define the half-meshes (on which Q is computed).
!
      do j=1,ntss-1
        tssh(j)=half*(tss(j)+tss(j+1))
      enddo
!
      do k=1,npss-1
        pssh(k)=half*(pss(k)+pss(k+1))
      enddo
!
!cc$omp parallel do
!cc$omp& default(private)
!cc$omp& shared(efl,pfl,tfl,qfl,pss,tss,tssh)
!cc$omp& schedule(dynamic,iterations_per_thread)
      do k=1,npss-1
        dp=pss(k+1)-pss(k)
        do j=1,ntss-1
          dt=tss(j+1)-tss(j)
          if (efl(j  ,k  ).eq.0..or.efl(j+1,k  ).eq.0..or.&
             &efl(j  ,k+1).eq.0..or.efl(j+1,k+1).eq.0.) then
            qfl(j,k)=0.
          else
            efav=quarter*(efl(j,k  )+efl(j+1,k  )+&
     &                    efl(j,k+1)+efl(j+1,k+1))
            if (efav.ne.0.) then
              tmav=quarter*(tfl(j,k  )+tfl(j+1,k  )+&
     &                      tfl(j,k+1)+tfl(j+1,k+1))
              stm=sin(tmav)
              stp=sin(tssh(j))
              dtdt_m=(tfl(j+1,k  )-tfl(j  ,k  ))/dt
              dtdt_p=(tfl(j+1,k+1)-tfl(j  ,k+1))/dt
              dtdp_m=(tfl(j  ,k+1)-tfl(j  ,k  ))/dp
              dtdp_p=(tfl(j+1,k+1)-tfl(j+1,k  ))/dp
              dpdt_m=modulo_twopi(pfl(j+1,k  )-pfl(j  ,k  ))/dt
              dpdt_p=modulo_twopi(pfl(j+1,k+1)-pfl(j  ,k+1))/dt
              dpdp_m=modulo_twopi(pfl(j  ,k+1)-pfl(j  ,k  ))/dp
              dpdp_p=modulo_twopi(pfl(j+1,k+1)-pfl(j+1,k  ))/dp
              dtdt=half*(dtdt_m+dtdt_p)
              dtdp=half*(dtdp_m+dtdp_p)
              dpdt=half*(dpdt_m+dpdt_p)
              dpdp=half*(dpdp_m+dpdp_p)
              aa=stm*dpdp/stp
              bb=stm*dpdt
              cc=dtdp/stp
              dd=dtdt
              qfl(j,k)=(aa**2+bb**2+cc**2+dd**2)/efav
            else
              qfl(j,k)=0.
            end if
          end if
        enddo
      enddo
!
      if (qffile.ne.' ') then
        if (verbose) then
          if (.not.wrote_cr) then
            write (*,*)
            wrote_cr=.true.
          end if
          write (*,*) 'Writing the forward mapping '//&
     &                'Q to file: ',&
     &                trim(qffile)
        end if
        call wrhdf_2d (qffile,.true.,ntss-1,npss-1,qfl,tssh,pssh,&
     &                 hdf32,ierr)
        if (ierr.ne.0) then
          write (*,*)
          write (*,*) '### ERROR in MAP_FORWARD:'
          write (*,*) '### Could not write the Q factor file.'
          call exit (1)
        end if
      end if
!
      if (slogqffile.ne.' ') then
!
        call slogq (qfl,slogqfl,b%lim0(1))
!
        if (verbose) then
          write (*,*) 'Writing the forward mapping '//&
     &                'SLOG(Q) to file: ',&
     &                trim(slogqffile)
        end if
        call wrhdf_2d(slogqffile,.true.,ntss-1,npss-1,slogqfl,tssh,pssh,&
     &                hdf32,ierr)
        if (ierr.ne.0) then
          write (*,*)
          write (*,*) '### ERROR in MAP_FORWARD:'
          write (*,*) '### Could not write the SLOG(Q) factor file.'
          call exit (1)
        end if
        deallocate (slogqfl)
      end if
!
      deallocate (tssh)
      deallocate (pssh)
      deallocate (qfl)
!
      return
      end
!#######################################################################
      subroutine slogq (qfl,slogqfl,rlevel)
!
!-----------------------------------------------------------------------
!
! ****** Compute Slava's "signed log" of Q.
!
!-----------------------------------------------------------------------
!
      use number_types
      use field
      use mesh
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      real(r_typ), parameter :: half=.5_r_typ
      real(r_typ), parameter :: one=1.0_r_typ
!
!-----------------------------------------------------------------------
!
      real(r_typ), dimension(ntss-1,npss-1) :: qfl,slogqfl
      real(r_typ) :: rlevel
!
!-----------------------------------------------------------------------
!
      real(r_typ), dimension(ntss-1) :: tssh
      real(r_typ), dimension(npss-1) :: pssh
!
!-----------------------------------------------------------------------
!
      integer :: i,j,k
      real(r_typ) :: q,lq,br
      type(csvec) :: x,bv
!
!-----------------------------------------------------------------------
!
! ****** Define the half-meshes (on which Q was computed).
!
      do j=1,ntss-1
        tssh(j)=half*(tss(j)+tss(j+1))
      enddo
!
      do k=1,npss-1
        pssh(k)=half*(pss(k)+pss(k+1))
      enddo
!
! ***** Calculate slogq
!
!$omp parallel do collapse(2)
!$omp& default(shared) private(i,j,q,lq,x,br,bv)
      do i=1,ntss-1
        do j=1,npss-1
          q=half*qfl(i,j)
          if (q.lt.one) then
            q=one
          else
            q=q+sqrt(q**2-one)
          end if
!
          lq=log10(q)
!
! ****** Get Br at the current location.
!
          x%s(1)=rlevel
          x%s(2)=tssh(i)
          x%s(3)=pssh(j)
          call getb (b,x,bv)
          br=bv%s(1)
          if (br.lt.0) then
            lq=-lq
          end if
          slogqfl(i,j)=lq
        enddo
      enddo
!$omp end parallel do
!
      end subroutine
!#######################################################################
      subroutine map_backward
!
!-----------------------------------------------------------------------
!
! ****** Trace field lines inward from r=R1.
!
!-----------------------------------------------------------------------
!
      use number_types
      use types
      use mesh
      use field
      use vars
      use files
      use params
      use field_line_params
      use diags
      use openmp_vars
      use tracefl_interface
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      real(r_typ), parameter :: half=.5_r_typ
      real(r_typ), parameter :: quarter=.25_r_typ
!
!-----------------------------------------------------------------------
!
! ****** Storage for the mapping.
!
      real(r_typ), dimension(ntss,npss) :: rfl,tfl,pfl,efl,kfl,length
!
      real(r_typ), dimension(:), allocatable :: tssh,pssh
      real(r_typ), dimension(:,:), allocatable :: qfl,slogqfl
!
!-----------------------------------------------------------------------
!
      integer :: ierr,j,k,nbad
      real(r_typ), dimension(3) :: xfl0,xfl1
      real(r_typ), dimension(3) :: bs0,bs1
      logical :: ttb
      real(r_typ) :: s
      real(r_typ) :: dtdt_m,dtdt_p
      real(r_typ) :: dtdp_m,dtdp_p
      real(r_typ) :: dpdt_m,dpdt_p
      real(r_typ) :: dpdp_m,dpdp_p
      real(r_typ) :: dtdt,dtdp,dpdt,dpdp
      real(r_typ) :: dt,dp,aa,bb,cc,dd,stm,stp,tmav,efav
      logical :: wrote_cr
      integer :: n_completed,n_total,nc,diag_step
      real(r_typ) :: pct_done
!
!-----------------------------------------------------------------------
!
      real(r_typ), external :: modulo_twopi
!
!-----------------------------------------------------------------------
!
! ****** Trace field lines, starting from each (T,P) cell at r=R1,
! ****** until the field line hits r=R0, or goes back to r=R1,
! ****** or exhausts the number of segments allowed.
!
      if (verbose) then
        write (*,*)
        write (*,*) '### Computing a backward mapping from R1:'
      end if
!
      if (verbose) then
        write (*,*)
        write (*,*) 'Mapping field lines ...'
        write (*,*)
      end if
!
      nbad=0
!
      ds%direction_is_along_b=.false.
      ds%direction=-1
!
      n_total=ntss*npss
      n_completed=0
!
!$omp parallel do
!$omp& default(shared)
!$omp& private(j,k,xfl0,xfl1,bs0,bs1,s,ttb)
!$omp& private(nc,diag_step,pct_done)
!$omp& collapse(2)
!$omp& schedule(dynamic,iterations_per_thread)
      do k=1,npss
        do j=1,ntss
!
! ****** Update the iteration counter for diagnostic
! ****** purposes.
!
          if (verbose) then
!$omp critical (omp_nc)
            n_completed=n_completed+1
            nc=n_completed
!$omp end critical (omp_nc)
          end if
!
          xfl0(1)=b%lim1(1)
          xfl0(2)=tss(j)
          xfl0(3)=pss(k)
!
          call tracefl (b,ds,xfl0,xfl1,bs0,bs1,s,ttb)
!
          if (ttb) then
            rfl(j,k)=xfl1(1)
            tfl(j,k)=xfl1(2)
            pfl(j,k)=xfl1(3)
            length(j,k)=s
            if (bs0(1).ne.0.) then
              efl(j,k)=abs((bs1(1)*xfl1(1)**2)/(bs0(1)*xfl0(1)**2))
            else
              efl(j,k)=0.
            end if
            if (bs1(1).ne.0.) then
              kfl(j,k)=log10(max(abs(bs0(1)/bs1(1)),tiny(bs0(1))))
            else
              kfl(j,k)=-50._r_typ
            end if
          else
!$omp critical (nbad_count)
            nbad=nbad+1
            write (*,*)
            write (*,*) '### WARNING from MAP_BACKWARD:'
            write (*,*) '### A field line did not reach R0 or R1.'
            write (*,*) 'Initial theta = ',xfl0(2)
            write (*,*) 'Initial phi   = ',xfl0(3)
            write (*,*) 'Final field line radius = ',xfl1(1)
!$omp end critical (nbad_count)
            rfl(j,k)=-1._r_typ
            tfl(j,k)=-1._r_typ
            pfl(j,k)=-1._r_typ
            length(j,k)=0.
            efl(j,k)=0.
            kfl(j,k)=-50._r_typ
          end if
!
          if (max_bad_fieldlines.gt.0) then
            if (nbad.gt.max_bad_fieldlines) then
!$omp critical (nbad2)
              write (*,*)
              write (*,*) '### ERROR in MAP_BACKWARD:'
              write (*,*) '### Too many field lines did not reach'//&
     &                    ' R0 or R1.'
              write (*,*) 'Number of bad traces = ',max_bad_fieldlines
              call exit (1)
!$omp end critical (nbad2)
            end if
          end if
!
! ****** Write progress diagnostics if requested.
!
          if (verbose) then
            diag_step=mod(nc,diagnostic_interval)
            if (diag_step.eq.0) then
              pct_done=100.*nc/n_total
              write (*,910) 'Fraction completed: ',pct_done
  910         format (1x,a,f7.3,'%')
            end if
          end if
!
        enddo
      enddo
!$omp end parallel do
!
! ****** Write the mapping.
!
      wrote_cr=.false.
!
      if (rbfile.ne.' ') then
        if (verbose) then
          if (.not.wrote_cr) then
            write (*,*)
            wrote_cr=.true.
          end if
          write (*,*) 'Writing the backward mapping for coordinate '//&
     &                'r to file: ',&
     &                trim(rbfile)
        end if
        call wrhdf_2d (rbfile,.true.,ntss,npss,rfl,tss,pss,&
     &                 hdf32,ierr)
        if (ierr.ne.0) then
          write (*,*)
          write (*,*) '### ERROR in MAP_BACKWARD:'
          write (*,*) '### Could not write the backward mapping'//&
     &                ' file for coordinate r.'
          call exit (1)
        end if
      end if
!
      if (tbfile.ne.' ') then
        if (verbose) then
          if (.not.wrote_cr) then
            write (*,*)
            wrote_cr=.true.
          end if
          write (*,*) 'Writing the backward mapping for coordinate '//&
     &                't to file: ',&
     &                trim(tbfile)
        end if
        call wrhdf_2d (tbfile,.true.,ntss,npss,tfl,tss,pss,&
     &                 hdf32,ierr)
        if (ierr.ne.0) then
          write (*,*)
          write (*,*) '### ERROR in MAP_BACKWARD:'
          write (*,*) '### Could not write the backward mapping'//&
     &                ' file for coordinate t.'
          call exit (1)
        end if
      end if
!
      if (pbfile.ne.' ') then
        if (verbose) then
          if (.not.wrote_cr) then
            write (*,*)
            wrote_cr=.true.
          end if
          write (*,*) 'Writing the backward mapping for coordinate '//&
     &                'p to file: ',&
     &                trim(pbfile)
        end if
        call wrhdf_2d (pbfile,.true.,ntss,npss,pfl,tss,pss,&
     &                 hdf32,ierr)
        if (ierr.ne.0) then
          write (*,*)
          write (*,*) '### ERROR in MAP_BACKWARD:'
          write (*,*) '### Could not write the backward mapping'//&
     &                ' file for coordinate p.'
          call exit (1)
        end if
      end if
!
      if (ebfile.ne.' ') then
        if (verbose) then
          if (.not.wrote_cr) then
            write (*,*)
            wrote_cr=.true.
          end if
          write (*,*) 'Writing the backward mapping '//&
     &                'expansion factor to file: ',&
     &                trim(ebfile)
        end if
        call wrhdf_2d (ebfile,.true.,ntss,npss,efl,tss,pss,&
     &                 hdf32,ierr)
        if (ierr.ne.0) then
          write (*,*)
          write (*,*) '### ERROR in MAP_BACKWARD:'
          write (*,*) '### Could not write the backward mapping'//&
     &                ' file for the expansion factor.'
          call exit (1)
        end if
      end if
!
      if (kbfile.ne.' ') then
        if (verbose) then
          if (.not.wrote_cr) then
            write (*,*)
            wrote_cr=.true.
          end if
          write (*,*) 'Writing the backward mapping '//&
     &                'K factor to file: ',&
     &                trim(kbfile)
        end if
        call wrhdf_2d (kbfile,.true.,ntss,npss,kfl,tss,pss,&
     &                 hdf32,ierr)
        if (ierr.ne.0) then
          write (*,*)
          write (*,*) '### ERROR in MAP_BACKWARD:'
          write (*,*) '### Could not write the backward mapping'//&
     &                ' file for the K factor.'
          call exit (1)
        end if
      end if
!
      if (lbfile.ne.' ') then
        if (verbose) then
          if (.not.wrote_cr) then
            write (*,*)
            wrote_cr=.true.
          end if
          write (*,*) 'Writing the backward mapping '//&
     &                'length to file: ',&
     &                trim(lbfile)
        end if
        call wrhdf_2d (lbfile,.true.,ntss,npss,length,tss,pss,&
     &                 hdf32,ierr)
        if (ierr.ne.0) then
          write (*,*)
          write (*,*) '### ERROR in MAP_BACKWARD:'
          write (*,*) '### Could not write the backward mapping'//&
     &                'length file.'
          call exit (1)
        end if
      end if
!
! ****** Compute Q (if requested).
!
      if (qbfile.eq.' '.and.slogqbfile.eq.' ') return
!
! ****** This can only be done if NTSS and NPSS exceed 1.
!
      if (ntss.le.1.or.npss.le.1) then
        write (*,*)
        write (*,*) '### WARNING from MAP_BACKWARD:'
        write (*,*) '### Could not compute the Q factor.'
        write (*,*) '### To compute Q, NTSS and NPSS'//&
     &              ' must be greater than 1.'
        return
      end if
!
      allocate (tssh(ntss-1))
      allocate (pssh(npss-1))
      allocate (qfl(ntss-1,npss-1))
      if (slogqbfile.ne.' ') then
        allocate (slogqfl(ntss-1,npss-1))
      end if
!
! ****** Define the half-meshes (on which Q is computed).
!
      do j=1,ntss-1
        tssh(j)=half*(tss(j)+tss(j+1))
      enddo
!
      do k=1,npss-1
        pssh(k)=half*(pss(k)+pss(k+1))
      enddo
!
      do k=1,npss-1
        dp=pss(k+1)-pss(k)
        do j=1,ntss-1
          dt=tss(j+1)-tss(j)
          if (efl(j  ,k  ).eq.0..or.efl(j+1,k  ).eq.0..or.&
             &efl(j  ,k+1).eq.0..or.efl(j+1,k+1).eq.0.) then
            qfl(j,k)=0.
          else
            efav=quarter*(efl(j,k  )+efl(j+1,k  )+&
     &                    efl(j,k+1)+efl(j+1,k+1))
            if (efav.ne.0.) then
              tmav=quarter*(tfl(j,k  )+tfl(j+1,k  )+&
     &                      tfl(j,k+1)+tfl(j+1,k+1))
              stm=sin(tmav)
              stp=sin(tssh(j))
              dtdt_m=(tfl(j+1,k  )-tfl(j  ,k  ))/dt
              dtdt_p=(tfl(j+1,k+1)-tfl(j  ,k+1))/dt
              dtdp_m=(tfl(j  ,k+1)-tfl(j  ,k  ))/dp
              dtdp_p=(tfl(j+1,k+1)-tfl(j+1,k  ))/dp
              dpdt_m=modulo_twopi(pfl(j+1,k  )-pfl(j  ,k  ))/dt
              dpdt_p=modulo_twopi(pfl(j+1,k+1)-pfl(j  ,k+1))/dt
              dpdp_m=modulo_twopi(pfl(j  ,k+1)-pfl(j  ,k  ))/dp
              dpdp_p=modulo_twopi(pfl(j+1,k+1)-pfl(j+1,k  ))/dp
              dtdt=half*(dtdt_m+dtdt_p)
              dtdp=half*(dtdp_m+dtdp_p)
              dpdt=half*(dpdt_m+dpdt_p)
              dpdp=half*(dpdp_m+dpdp_p)
              aa=stm*dpdp/stp
              bb=stm*dpdt
              cc=dtdp/stp
              dd=dtdt
              qfl(j,k)=(aa**2+bb**2+cc**2+dd**2)*efav
            else
              qfl(j,k)=0.
            end if
          end if
        enddo
      enddo
!
      if (qbfile.ne.' ') then
        if (verbose) then
          if (.not.wrote_cr) then
            write (*,*)
            wrote_cr=.true.
          end if
          write (*,*) 'Writing the backward mapping '//&
     &                'Q to file: ',&
     &                trim(qbfile)
        end if
        call wrhdf_2d (qbfile,.true.,ntss-1,npss-1,qfl,tssh,pssh,&
     &                 hdf32,ierr)
        if (ierr.ne.0) then
          write (*,*)
          write (*,*) '### ERROR in MAP_BACKWARD:'
          write (*,*) '### Could not write the Q factor file.'
          call exit (1)
        end if
      end if
!
      if (slogqbfile.ne.' ') then
!
        call slogq (qfl,slogqfl,b%lim1(1))
!
        if (verbose) then
          write (*,*) 'Writing the backward mapping '//&
     &                'SLOG(Q) to file: ',&
     &                trim(slogqbfile)
        end if
        call wrhdf_2d(slogqbfile,.true.,ntss-1,npss-1,slogqfl,tssh,pssh,&
     &                hdf32,ierr)
        if (ierr.ne.0) then
          write (*,*)
          write (*,*) '### ERROR in MAP_BACKWARD:'
          write (*,*) '### Could not write the SLOG(Q) factor file.'
          call exit (1)
        end if
        deallocate (slogqfl)
      end if
!
      deallocate (tssh)
      deallocate (pssh)
      deallocate (qfl)
!
      return
      end
!#######################################################################
      function modulo_twopi (x)
!
!-----------------------------------------------------------------------
!
! ****** Return the smallest value of X, modulo 2*pi.
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
      real(r_typ) :: x
      real(r_typ) :: modulo_twopi
!
!-----------------------------------------------------------------------
!
      real(r_typ) :: xm,xp,x_min
!
!-----------------------------------------------------------------------
!
      xm=abs(x-twopi)
      xp=abs(x+twopi)
      x_min=min(xm,abs(x),xp)
      if (xm.eq.x_min) then
        modulo_twopi=x-twopi
      else if (xp.eq.x_min) then
        modulo_twopi=x+twopi
      else
        modulo_twopi=x
      end if
!
      return
      end
!#######################################################################
      subroutine map_3d
!
!-----------------------------------------------------------------------
!
! ****** Trace field lines from every point on a 3D (r,t,p) mesh.
!
!-----------------------------------------------------------------------
!
      use number_types
      use types
      use mesh
      use field
      use vars
      use files
      use params
      use diags
      use openmp_vars
      use tracefl_interface
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
! ****** Storage for the mapping.
!
      real(r_typ), dimension(nrss,ntss,npss) :: rfl,tfl,pfl
!
!-----------------------------------------------------------------------
!
      integer :: ierr,i,j,k
      real(r_typ), dimension(3) :: xfl0,xfl1
      real(r_typ), dimension(3) :: bs0,bs1
      logical :: ttb
      real(r_typ) :: s
      logical :: wrote_cr
      integer :: n_completed,n_total,nc,diag_step
      real(r_typ) :: pct_done
!
!-----------------------------------------------------------------------
!
! ****** Trace field lines, starting from each (R,T,P) cell,
! ****** until the field line hits the boundaries or exhausts
! ****** the field line length allowed.
!
      if (verbose) then
        write (*,*)
        write (*,*) '### Computing a mapping in 3D:'
      end if
!
      if (verbose) then
        write (*,*)
        write (*,*) 'Mapping field lines ...'
        write (*,*)
      end if
!
      ds%direction_is_along_b=.false.
      ds%direction=-1
!
      n_total=nrss*ntss*npss
      n_completed=0
!
!$omp parallel do
!$omp& private(i,j,k,xfl0,xfl1,bs0,bs1,s,ttb)
!$omp& private(nc,diag_step,pct_done)
!$omp& collapse(3)
!$omp& schedule(dynamic,iterations_per_thread)
      do k=1,npss
        do j=1,ntss
          do i=1,nrss
!
! ****** Update the iteration counter for diagnostic
! ****** purposes.
!
            if (verbose) then
!$omp critical
              n_completed=n_completed+1
              nc=n_completed
!$omp end critical
            end if
!
            xfl0(1)=rss(i)
            xfl0(2)=tss(j)
            xfl0(3)=pss(k)
!
            call tracefl (b,ds,xfl0,xfl1,bs0,bs1,s,ttb)
!
            if (ttb) then
              rfl(i,j,k)=xfl1(1)
            else
              rfl(i,j,k)=-xfl1(1)
            end if
            tfl(i,j,k)=xfl1(2)
            pfl(i,j,k)=xfl1(3)
!
! ****** Write progress diagnostics if requested.
!
            if (verbose) then
              diag_step=mod(nc,diagnostic_interval)
              if (diag_step.eq.0) then
                pct_done=100.*nc/n_total
                write (*,910) 'Fraction completed: ',pct_done
  910           format (1x,a,f7.3,'%')
              end if
            end if
!
          enddo
        enddo
      enddo
!$omp end parallel do
!
! ****** Write the mapping.
!
      wrote_cr=.false.
!
      if (volume3d_output_file%r.ne.' ') then
        if (verbose) then
          if (.not.wrote_cr) then
            write (*,*)
            wrote_cr=.true.
          end if
          write (*,*) 'Writing the 3D mapping for coordinate '//&
     &                'r to file: ',&
     &                trim(volume3d_output_file%r)
        end if
        call wrhdf_3d (volume3d_output_file%r,.true.,&
     &                 nrss,ntss,npss,rfl,rss,tss,pss,&
     &                 hdf32,ierr)
        if (ierr.ne.0) then
          write (*,*)
          write (*,*) '### ERROR in MAP_3D:'
          write (*,*) '### Could not write the 3D mapping'//&
     &                ' file for coordinate r.'
          call exit (1)
        end if
      end if
!
      if (volume3d_output_file%t.ne.' ') then
        if (verbose) then
          if (.not.wrote_cr) then
            write (*,*)
            wrote_cr=.true.
          end if
          write (*,*) 'Writing the 3D mapping for coordinate '//&
     &                't to file: ',&
     &                trim(volume3d_output_file%t)
        end if
        call wrhdf_3d (volume3d_output_file%t,.true.,&
     &                 nrss,ntss,npss,tfl,rss,tss,pss,&
     &                 hdf32,ierr)
        if (ierr.ne.0) then
          write (*,*)
          write (*,*) '### ERROR in MAP_3D:'
          write (*,*) '### Could not write the 3D mapping'//&
     &                ' file for coordinate t.'
          call exit (1)
        end if
      end if
!
      if (volume3d_output_file%p.ne.' ') then
        if (verbose) then
          if (.not.wrote_cr) then
            write (*,*)
            wrote_cr=.true.
          end if
          write (*,*) 'Writing the 3D mapping for coordinate '//&
     &                'p to file: ',&
     &                trim(volume3d_output_file%p)
        end if
        call wrhdf_3d (volume3d_output_file%p,.true.,&
     &                 nrss,ntss,npss,pfl,rss,tss,pss,&
     &                 hdf32,ierr)
        if (ierr.ne.0) then
          write (*,*)
          write (*,*) '### ERROR in MAP_3D:'
          write (*,*) '### Could not write the 3D mapping'//&
     &                ' file for coordinate p.'
          call exit (1)
        end if
      end if
!
      return
      end
!#######################################################################
      subroutine read_slice_coordinates
!
!-----------------------------------------------------------------------
!
! ****** Read the coordinates that define the slice in the 3D volume.
!
!-----------------------------------------------------------------------
!
      use number_types
      use types
      use mesh
      use field
      use vars
      use files
      use params
      use sds_def
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      integer :: ierr
!
!-----------------------------------------------------------------------
!
      logical, external :: same_structure_sds
!
!-----------------------------------------------------------------------
!
! ****** Set the coordinate names.
!
      if (slice_coords_are_xyz) then
        slice_coord_name=(/'x','y','z'/)
      else
        slice_coord_name=(/'r','t','p'/)
      end if
!
! ****** Read the coordinates of the slice.
!
      if (verbose) then
        write (*,*)
        write (*,*) '### Reading the coordinates of the slice ...'
      end if
!
! ****** Read the x/r coordinate file.
!
      if (verbose) then
        write (*,*)
        write (*,*) 'Reading the '//slice_coord_name(1)//&
     &              ' coordinate from file: ',&
     &              trim(slice_input_file%r)
      end if
!
      call rdhdf (slice_input_file%r,slice_c1,ierr)
!
      if (ierr.ne.0) then
        write (*,*)
        write (*,*) '### ERROR in READ_SLICE_COORDINATES:'
        write (*,*) '### Could not read the '//slice_coord_name(1)//&
     &              ' coordinate.'
        write (*,*) 'IERR (from RDHDF) = ',ierr
        write (*,*) 'File name: ',trim(slice_input_file%r)
        call exit (1)
      end if
!
! ****** Read the y/t coordinate file.
!
      if (verbose) then
        write (*,*) 'Reading the '//slice_coord_name(2)//&
     &              ' coordinate from file: ',&
     &              trim(slice_input_file%t)
      end if
!
      call rdhdf (slice_input_file%t,slice_c2,ierr)
!
      if (ierr.ne.0) then
        write (*,*)
        write (*,*) '### ERROR in READ_SLICE_COORDINATES:'
        write (*,*) '### Could not read the '//slice_coord_name(2)//&
     &              ' coordinate.'
        write (*,*) 'IERR (from RDHDF) = ',ierr
        write (*,*) 'File name: ',trim(slice_input_file%t)
        call exit (1)
      end if
!
! ****** Check that the y/t coordinate has the same structure as the
! ****** x/r coordinate.
!
      if (.not.same_structure_sds(slice_c1,slice_c2)) then
        write (*,*)
        write (*,*) '### ERROR in READ_SLICE_COORDINATES:'
        write (*,*) '### The data sets for coordinates '//&
     &              slice_coord_name(1)//' and '//&
     &              slice_coord_name(2)//' do not have'//&
     &              ' the same structure.'
        write (*,*) 'Coordinate '//slice_coord_name(1)//&
     &              ' file name: ',trim(slice_input_file%r)
        write (*,*) 'Coordinate '//slice_coord_name(2)//&
     &              ' file name: ',trim(slice_input_file%t)
        call exit (1)
      end if
!
! ****** Read the z/p coordinate file.
!
      if (verbose) then
        write (*,*) 'Reading the '//slice_coord_name(3)//&
     &              ' coordinate from file: ',&
     &              trim(slice_input_file%p)
      end if
!
      call rdhdf (slice_input_file%p,slice_c3,ierr)
!
      if (ierr.ne.0) then
        write (*,*)
        write (*,*) '### ERROR in READ_SLICE_COORDINATES:'
        write (*,*) '### Could not read the '//slice_coord_name(3)//&
     &              ' coordinate.'
        write (*,*) 'IERR (from RDHDF) = ',ierr
        write (*,*) 'File name: ',trim(slice_input_file%p)
        call exit (1)
      end if
!
! ****** Check that the z/p coordinate has the same structure as the
! ****** x/r coordinate.
!
      if (.not.same_structure_sds(slice_c1,slice_c3)) then
        write (*,*)
        write (*,*) '### ERROR in READ_SLICE_COORDINATES:'
        write (*,*) '### The data sets for coordinates '//&
     &              slice_coord_name(1)//' and '//&
     &              slice_coord_name(3)//' do not have'//&
     &              ' the same structure.'
        write (*,*) 'Coordinate '//slice_coord_name(1)//&
     &              ' file name: ',trim(slice_input_file%r)
        write (*,*) 'Coordinate '//slice_coord_name(3)//&
     &              ' file name: ',trim(slice_input_file%p)
        call exit (1)
      end if
!
      return
      end
!#######################################################################
      subroutine deallocate_slice_coordinates
!
!-----------------------------------------------------------------------
!
      use sds_def
      use vars
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
! ****** Read the coordinates of the slice.
!
      call deallocate_sds (slice_c1)
      call deallocate_sds (slice_c2)
      call deallocate_sds (slice_c3)
!
      return
      end
!#######################################################################
      function same_structure_sds (s1,s2)
!
!-----------------------------------------------------------------------
!
! ****** Check if the two data sets S1 and S2 have the same
! ****** structure.  If they do, return .TRUE; otherwise, return
! ****** .FALSE.
!
!-----------------------------------------------------------------------
!
      use sds_def
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      type(sds) :: s1,s2
      logical :: same_structure_sds
!
!-----------------------------------------------------------------------
!
      integer :: i
!
!-----------------------------------------------------------------------
!
      same_structure_sds=.false.
!
      if (s1%ndim.ne.s2%ndim) return
!
      if (s1%scale.neqv.s2%scale) return
!
      do i=1,s1%ndim
        if (s1%dims(i).ne.s2%dims(i)) return
      enddo
!
      same_structure_sds=.true.
!
      return
      end
!#######################################################################
      subroutine map_slice
!
!-----------------------------------------------------------------------
!
! ****** Trace field lines from points on a slice in the 3D volume.
!
!-----------------------------------------------------------------------
!
      use number_types
      use types
      use mesh
      use field
      use vars
      use files
      use params
      use sds_def
      use diags
      use openmp_vars
      use tracefl_interface
      use debug
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
! ****** Storage for the forward mapping.
!
      real(r_typ), dimension(:,:,:), allocatable, target :: rfl_f
      real(r_typ), dimension(:,:,:), allocatable, target :: tfl_f
      real(r_typ), dimension(:,:,:), allocatable, target :: pfl_f
!
! ****** Storage for the backward mapping.
!
      real(r_typ), dimension(:,:,:), allocatable, target :: rfl_b
      real(r_typ), dimension(:,:,:), allocatable, target :: tfl_b
      real(r_typ), dimension(:,:,:), allocatable, target :: pfl_b
!
!-----------------------------------------------------------------------
!
      type (sds) :: out
      integer :: ierr,i,j,k,n1,n2,n3
      real(r_typ), dimension(3) :: xfl0,xfl1
      real(r_typ), dimension(3) :: bs0,bs1
      real(r_typ), dimension(3) :: c
      logical :: ttb
      real(r_typ) :: s
      type(flparam) :: ds_f,ds_b
      logical :: wrote_cr
      integer :: n_completed,n_total,nc,diag_step
      real(r_typ) :: pct_done
!
!-----------------------------------------------------------------------
!
! ****** Field line trace storage buffers.
!
      type(traj), dimension(:,:,:), allocatable :: xtf
      type(traj), dimension(:,:,:), allocatable :: xtb
      real(r_typ), dimension(:,:), allocatable :: fl
      real(r_typ), dimension(3) :: xyz
      integer :: n,l
      real(r_typ) :: dummy
      character(4) :: ch4
      character(256) :: fname
!
!-----------------------------------------------------------------------
!
! ****** Map the field lines for all points on a slice in 3D.
!
      if (verbose) then
        write (*,*)
        write (*,*) '### Computing a mapping starting on a slice:'
      end if
!
! ****** Check that at least one tracing direction was requested.
!
      if (.not.(trace_from_slice_forward.or.&
     &          trace_from_slice_backward)) then
        write (*,*)
        write (*,*) '### ERROR in MAP_SLICE:'
        write (*,*) '### At least one tracing direction (forward'//&
     &              ' and/or backward) must be requested.'
        call exit (1)
      end if
!
! ****** Set the tracing direction to be either along the direction
! ****** of the magnetic field or along the directon of increasing
! ****** radius.
!
      ds%direction_is_along_b=trace_slice_direction_is_along_b
!
      ds_f=ds
      ds_f%direction=1
!
      ds_b=ds
      ds_b%direction=-1
!
      if (verbose) then
        if (trace_slice_direction_is_along_b) then
          write (*,*)
          write (*,*) '### The forward tracing direction is'//&
     &                ' along the direction of B.'
          write (*,*) '### The backward tracing direction is'//&
     &                ' opposite to the direction of B.'
        else
          write (*,*)
          write (*,*) '### The forward tracing direction is'//&
     &                ' along the direction of increasing radius.'
          write (*,*) '### The backward tracing direction is'//&
     &                ' along the direction of decreasing radius.'
        end if
      end if
!
      n1=slice_c1%dims(1)
      n2=slice_c1%dims(2)
      n3=slice_c1%dims(3)
!
! ****** Allocate storage for the mapping.
!
      if (trace_from_slice_forward) then
        allocate (rfl_f(n1,n2,n3))
        allocate (tfl_f(n1,n2,n3))
        allocate (pfl_f(n1,n2,n3))
      end if
!
      if (trace_from_slice_backward) then
        allocate (rfl_b(n1,n2,n3))
        allocate (tfl_b(n1,n2,n3))
        allocate (pfl_b(n1,n2,n3))
      end if
!
! ****** Allocate the field line storage buffers if the field
! ****** line traces are being written to HDF files.
!
      if (write_traces_to_hdf) then
        if (trace_from_slice_forward) then
          allocate (xtf(n1,n2,n3))
          do k=1,n3
            do j=1,n2
              do i=1,n1
                call allocate_trajectory_buffer (xtf(i,j,k))
              enddo
            enddo
          enddo
        end if
        if (trace_from_slice_backward) then
          allocate (xtb(n1,n2,n3))
          do k=1,n3
            do j=1,n2
              do i=1,n1
                call allocate_trajectory_buffer (xtb(i,j,k))
              enddo
            enddo
          enddo
        end if
      end if
!
      if (verbose) then
        write (*,*)
        write (*,*) 'Mapping field lines ...'
        write (*,*)
      end if
!
! ****** Trace a field line from each point on the slice until the
! ****** field line hits the boundaries or exhausts the length
! ****** allowed.
!
      n_total=n1*n2*n3
      n_completed=0
!
!$omp parallel do
!$omp& private(i,j,k,c,xfl0,xfl1,bs0,bs1,s,ttb)
!$omp& private(nc,diag_step,pct_done)
!$omp& collapse(3)
!$omp& schedule(dynamic,iterations_per_thread)
      do k=1,n3
        do j=1,n2
          do i=1,n1
!
! ****** Update the iteration counter for diagnostic
! ****** purposes.
!
            if (verbose) then
!$omp critical
              n_completed=n_completed+1
              nc=n_completed
!$omp end critical
            end if
!
            if (slice_coords_are_xyz) then
              c=(/slice_c1%f(i,j,k),&
     &            slice_c2%f(i,j,k),&
     &            slice_c3%f(i,j,k)/)
              call c2s (c,xfl0)
            else
              xfl0=(/slice_c1%f(i,j,k),&
     &               slice_c2%f(i,j,k),&
     &               slice_c3%f(i,j,k)/)
            end if
!
! ****** Launch a field line in the positive direction, if requested.
!
            if (trace_from_slice_forward) then
!
              if (write_traces_to_hdf) then
                call tracefl (b,ds_f,xfl0,xfl1,bs0,bs1,s,ttb,&
     &                        xtf(i,j,k))
              else
                call tracefl (b,ds_f,xfl0,xfl1,bs0,bs1,s,ttb)
              end if
!
              if (slice_coords_are_xyz) then
                if (ttb) then
                  call s2c (xfl1,c)
                  rfl_f(i,j,k)=c(1)
                  tfl_f(i,j,k)=c(2)
                  pfl_f(i,j,k)=c(3)
                else
                  rfl_f(i,j,k)=0.
                  tfl_f(i,j,k)=0.
                  pfl_f(i,j,k)=0.
                end if
              else
                if (ttb) then
                  rfl_f(i,j,k)=xfl1(1)
                else
                  rfl_f(i,j,k)=-xfl1(1)
                end if
                tfl_f(i,j,k)=xfl1(2)
                pfl_f(i,j,k)=xfl1(3)
              end if
!
            end if
!
! ****** Launch a field line in the negative direction, if requested.
!
            if (trace_from_slice_backward) then
!
              if (write_traces_to_hdf) then
                call tracefl (b,ds_b,xfl0,xfl1,bs0,bs1,s,ttb,&
     &                        xtb(i,j,k))
              else
                call tracefl (b,ds_b,xfl0,xfl1,bs0,bs1,s,ttb)
              end if
!
              if (slice_coords_are_xyz) then
                if (ttb) then
                  call s2c (xfl1,c)
                  rfl_b(i,j,k)=c(1)
                  tfl_b(i,j,k)=c(2)
                  pfl_b(i,j,k)=c(3)
                else
                  rfl_b(i,j,k)=0.
                  tfl_b(i,j,k)=0.
                  pfl_b(i,j,k)=0.
                end if
              else
                if (ttb) then
                  rfl_b(i,j,k)=xfl1(1)
                else
                  rfl_b(i,j,k)=-xfl1(1)
                end if
                tfl_b(i,j,k)=xfl1(2)
                pfl_b(i,j,k)=xfl1(3)
              end if
!
            end if
!
! ****** Write progress diagnostics if requested.
!
            if (verbose) then
              diag_step=mod(nc,diagnostic_interval)
              if (diag_step.eq.0) then
                pct_done=100.*nc/n_total
                write (*,910) 'Fraction completed: ',pct_done
  910           format (1x,a,f7.3,'%')
              end if
            end if
!
          enddo
        enddo
      enddo
!$omp end parallel do
!
! ****** Write the forward mapping.
!
      wrote_cr=.false.
!
      if (trace_from_slice_forward.and.&
     &    slice_output_file_forward%r.ne.' ') then
        out%ndim=slice_c1%ndim
        out%dims=slice_c1%dims
        out%scale=slice_c1%scale
        out%hdf32=slice_c1%hdf32
        out%scales(1)%f=>slice_c1%scales(1)%f
        out%scales(2)%f=>slice_c1%scales(2)%f
        out%scales(3)%f=>slice_c1%scales(3)%f
        out%f=>rfl_f
        if (verbose) then
          if (.not.wrote_cr) then
            write (*,*)
            wrote_cr=.true.
          end if
          write (*,*) 'Writing the forward mapping for'//&
     &                ' coordinate '//slice_coord_name(1)//&
     &                ' to file: ',&
     &                trim(slice_output_file_forward%r)
        end if
        call wrhdf (slice_output_file_forward%r,out,ierr)
        if (ierr.ne.0) then
          write (*,*)
          write (*,*) '### ERROR in MAP_SLICE:'
          write (*,*) '### Could not write the forward'//&
     &                ' mapping file for coordinate '//&
     &                slice_coord_name(1)//'.'
          call exit (1)
        end if
      end if
!
      if (trace_from_slice_forward.and.&
     &    slice_output_file_forward%t.ne.' ') then
        out%ndim=slice_c2%ndim
        out%dims=slice_c2%dims
        out%scale=slice_c2%scale
        out%hdf32=slice_c2%hdf32
        out%scales(1)%f=>slice_c2%scales(1)%f
        out%scales(2)%f=>slice_c2%scales(2)%f
        out%scales(3)%f=>slice_c2%scales(3)%f
        out%f=>tfl_f
        if (verbose) then
          if (.not.wrote_cr) then
            write (*,*)
            wrote_cr=.true.
          end if
          write (*,*) 'Writing the forward mapping for'//&
     &                ' coordinate '//slice_coord_name(2)//&
     &                ' to file: ',&
     &                trim(slice_output_file_forward%t)
        end if
        call wrhdf (slice_output_file_forward%t,out,ierr)
        if (ierr.ne.0) then
          write (*,*)
          write (*,*) '### ERROR in MAP_SLICE:'
          write (*,*) '### Could not write the forward'//&
     &                ' mapping file for coordinate '//&
     &                slice_coord_name(2)//'.'
          call exit (1)
        end if
      end if
!
      if (trace_from_slice_forward.and.&
     &    slice_output_file_forward%p.ne.' ') then
        out%ndim=slice_c3%ndim
        out%dims=slice_c3%dims
        out%scale=slice_c3%scale
        out%hdf32=slice_c3%hdf32
        out%scales(1)%f=>slice_c3%scales(1)%f
        out%scales(2)%f=>slice_c3%scales(2)%f
        out%scales(3)%f=>slice_c3%scales(3)%f
        out%f=>pfl_f
        if (verbose) then
          if (.not.wrote_cr) then
            write (*,*)
            wrote_cr=.true.
          end if
          write (*,*) 'Writing the forward mapping for'//&
     &                ' coordinate '//slice_coord_name(3)//&
     &                ' to file: ',&
     &                trim(slice_output_file_forward%p)
        end if
        call wrhdf (slice_output_file_forward%p,out,ierr)
        if (ierr.ne.0) then
          write (*,*)
          write (*,*) '### ERROR in MAP_SLICE:'
          write (*,*) '### Could not write the forward'//&
     &                ' mapping file for coordinate '//&
     &                slice_coord_name(3)//'.'
          call exit (1)
        end if
      end if
!
! ****** Write the backward mapping.
!
      wrote_cr=.false.
!
      if (trace_from_slice_backward.and.&
     &    slice_output_file_backward%r.ne.' ') then
        out%ndim=slice_c1%ndim
        out%dims=slice_c1%dims
        out%scale=slice_c1%scale
        out%hdf32=slice_c1%hdf32
        out%scales(1)%f=>slice_c1%scales(1)%f
        out%scales(2)%f=>slice_c1%scales(2)%f
        out%scales(3)%f=>slice_c1%scales(3)%f
        out%f=>rfl_b
        if (verbose) then
          if (.not.wrote_cr) then
            write (*,*)
            wrote_cr=.true.
          end if
          write (*,*) 'Writing the backward mapping for'//&
     &                ' coordinate '//slice_coord_name(1)//&
     &                ' to file: ',&
     &                trim(slice_output_file_backward%r)
        end if
        call wrhdf (slice_output_file_backward%r,out,ierr)
        if (ierr.ne.0) then
          write (*,*)
          write (*,*) '### ERROR in MAP_SLICE:'
          write (*,*) '### Could not write the backward'//&
     &                ' mapping file for coordinate '//&
     &                slice_coord_name(1)//'.'
          call exit (1)
        end if
      end if
!
      if (trace_from_slice_backward.and.&
     &    slice_output_file_backward%t.ne.' ') then
        out%ndim=slice_c2%ndim
        out%dims=slice_c2%dims
        out%scale=slice_c2%scale
        out%hdf32=slice_c2%hdf32
        out%scales(1)%f=>slice_c2%scales(1)%f
        out%scales(2)%f=>slice_c2%scales(2)%f
        out%scales(3)%f=>slice_c2%scales(3)%f
        out%f=>tfl_b
        if (verbose) then
          if (.not.wrote_cr) then
            write (*,*)
            wrote_cr=.true.
          end if
          write (*,*) 'Writing the backward mapping for'//&
     &                ' coordinate '//slice_coord_name(2)//&
     &                ' to file: ',&
     &                trim(slice_output_file_backward%t)
        end if
        call wrhdf (slice_output_file_backward%t,out,ierr)
        if (ierr.ne.0) then
          write (*,*)
          write (*,*) '### ERROR in MAP_SLICE:'
          write (*,*) '### Could not write the backward'//&
     &                ' mapping file for coordinate '//&
     &                slice_coord_name(2)//'.'
          call exit (1)
        end if
      end if
!
      if (trace_from_slice_backward.and.&
     &    slice_output_file_backward%p.ne.' ') then
        out%ndim=slice_c3%ndim
        out%dims=slice_c3%dims
        out%scale=slice_c3%scale
        out%hdf32=slice_c3%hdf32
        out%scales(1)%f=>slice_c3%scales(1)%f
        out%scales(2)%f=>slice_c3%scales(2)%f
        out%scales(3)%f=>slice_c3%scales(3)%f
        out%f=>pfl_b
        if (verbose) then
          if (.not.wrote_cr) then
            write (*,*)
            wrote_cr=.true.
          end if
          write (*,*) 'Writing the backward mapping for'//&
     &                ' coordinate '//slice_coord_name(3)//&
     &                ' to file: ',&
     &                trim(slice_output_file_backward%p)
        end if
        call wrhdf (slice_output_file_backward%p,out,ierr)
        if (ierr.ne.0) then
          write (*,*)
          write (*,*) '### ERROR in MAP_SLICE:'
          write (*,*) '### Could not write the backward'//&
     &                ' mapping file for coordinate '//&
     &                slice_coord_name(3)//'.'
          call exit (1)
        end if
      end if
!
! ****** Write the forward field line traces to individual HDF
! ****** files if requested.
!
      if (write_traces_to_hdf.and.trace_from_slice_forward) then
!
        if (verbose) then
          write (*,*)
        end if
!
! ****** Loop over all points.
!
        n=0
        do k=1,n3
          do j=1,n2
            do i=1,n1
!
! ****** Allocate a temporary array to store the field line trace.
!
              allocate (fl(3,xtf(i,j,k)%npts))
!
! ****** Load the array with the field line coordinates.
!
              do l=1,xtf(i,j,k)%npts
                fl(1,l)=xtf(i,j,k)%x(1)%f(l)
                fl(2,l)=xtf(i,j,k)%x(2)%f(l)
                fl(3,l)=xtf(i,j,k)%x(3)%f(l)
                if (write_traces_as_xyz) then
                  call s2c (fl(1,l),xyz)
                  fl(1:3,l)=xyz
                end if
              enddo
!
! ****** Write the HDF files.
!
              n=n+1
              write (ch4,'(i4.4)') n
              fname=trim(write_traces_root)//'_f_'//ch4//'.'//trim(fmt)
!
              if (verbose) then
                write (*,*) 'Writing a forward field line trace'//&
     &                      ' to file: ',trim(fname)
              end if
!
              call wrhdf_2d (fname,.false.,&
     &                       3,xtf(i,j,k)%npts,fl,&
     &                       dummy,dummy,&
     &                       slice_c1%hdf32,ierr)
!
              if (ierr.ne.0) then
                write (*,*)
                write (*,*) '### ERROR in MAP_SLICE:'
                write (*,*) '### Could not write a forward'//&
     &                      ' field line trace to HDF file:'
                write (*,*) trim(fname)
                call exit (1)
              end if
!
              deallocate (fl)
!
              call deallocate_trajectory_buffer (xtf(i,j,k))
!
            enddo
          enddo
        enddo
!
        deallocate (xtf)
!
      end if
!
! ****** Write the backward field line traces to individual HDF
! ****** files if requested.
!
      if (write_traces_to_hdf.and.trace_from_slice_backward) then
!
        if (verbose) then
          write (*,*)
        end if
!
! ****** Loop over all points.
!
        n=0
        do k=1,n3
          do j=1,n2
            do i=1,n1
!
! ****** Allocate a temporary array to store the field line trace.
!
              allocate (fl(3,xtb(i,j,k)%npts))
!
! ****** Load the array with the field line coordinates.
!
              do l=1,xtb(i,j,k)%npts
                fl(1,l)=xtb(i,j,k)%x(1)%f(l)
                fl(2,l)=xtb(i,j,k)%x(2)%f(l)
                fl(3,l)=xtb(i,j,k)%x(3)%f(l)
                if (write_traces_as_xyz) then
                  call s2c (fl(1,l),xyz)
                  fl(1:3,l)=xyz
                end if
              enddo
!
! ****** Write the HDF files.
!
              n=n+1
              write (ch4,'(i4.4)') n
              fname=trim(write_traces_root)//'_b_'//ch4//'.'//trim(fmt)
!
              if (verbose) then
                write (*,*) 'Writing a backward field line trace'//&
     &                      ' to file: ',trim(fname)
              end if
!
              call wrhdf_2d (fname,.false.,&
     &                       3,xtb(i,j,k)%npts,fl,&
     &                       dummy,dummy,&
     &                       slice_c1%hdf32,ierr)
!
              if (ierr.ne.0) then
                write (*,*)
                write (*,*) '### ERROR in MAP_SLICE:'
                write (*,*) '### Could not write a backward'//&
     &                      ' field line trace to HDF file:'
                write (*,*) trim(fname)
                call exit (1)
              end if
!
              deallocate (fl)
!
              call deallocate_trajectory_buffer (xtb(i,j,k))
!
            enddo
          enddo
        enddo
!
        deallocate (xtb)
!
      end if
!
! ****** Deallocate memory.
!
      if (trace_from_slice_forward) then
        deallocate (rfl_f)
        deallocate (tfl_f)
        deallocate (pfl_f)
      end if
!
      if (trace_from_slice_backward) then
        deallocate (rfl_b)
        deallocate (tfl_b)
        deallocate (pfl_b)
      end if
!
      return
      end
!#######################################################################
      subroutine get_ch_map (rv)
!
!-----------------------------------------------------------------------
!
! ****** Compute a coronal hole map at radius r=RV.
!
!-----------------------------------------------------------------------
!
      use number_types
      use types
      use mesh
      use field
      use vars
      use files
      use params
      use diags
      use openmp_vars
      use tracefl_interface
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      real(r_typ) :: rv
!
!-----------------------------------------------------------------------
!
      real(r_typ), parameter :: one=1._r_typ
      real(r_typ), parameter :: two=2._r_typ
!
!-----------------------------------------------------------------------
!
! ****** Storage for the coronal hole map.
!
      real(r_typ), dimension(npss,ntss) :: ch
!
!-----------------------------------------------------------------------
!
      integer :: ierr,j,k
      real(r_typ), dimension(3) :: xfl0,xfl1
      real(r_typ), dimension(3) :: bs0,bs1
      logical :: ttb
      real(r_typ) :: s
      type(flparam) :: ds_f,ds_b
      logical :: f_trace_reached_boundary
      logical :: b_trace_reached_boundary
      logical :: f_trace_on_r0,f_trace_on_r1
      logical :: b_trace_on_r0,b_trace_on_r1
      logical :: f_br_positive
      logical :: b_br_positive
      integer :: n_completed,n_total,nc,diag_step
      real(r_typ) :: pct_done
!
!-----------------------------------------------------------------------
!
      if (verbose) then
        write (*,*)
        write (*,*) '### Computing a coronal hole map at r = ',rv
      end if
!
! ****** Check that the radius specified is valid.
!
      if (rv.lt.b%lim0(1).or.rv.gt.b%lim1(1)) then
        write (*,*)
        write (*,*) '### ERROR in GET_CH_MAP:'
        write (*,*) '### Invalid radius specified.'
        write (*,*) '### The radius is outside the domain limits:'
        write (*,*) 'Lower radial domain limit = ',b%lim0(1)
        write (*,*) 'Upper radial domain limit = ',b%lim1(1)
        write (*,*) 'Specified radius          = ',rv
        call exit (1)
      end if
!
! ****** Check that the coronal hole map output file name is not
! ****** blank, since this does not make sense.
!
      if (ch_map_output_file.eq.' ') then
        write (*,*)
        write (*,*) '### ERROR in GET_CH_MAP:'
        write (*,*) '### A coronal hole map was requested, yet'//&
     &              ' the output file name is blank.'
        call exit (1)
      end if
!
! ****** Set the tracing direction to be along the direction
! ****** of the magnetic field.
!
      ds%direction_is_along_b=.true.
!
      ds_f=ds
      ds_f%direction=1
!
      ds_b=ds
      ds_b%direction=-1
!
      if (verbose) then
        write (*,*)
        write (*,*) 'Mapping field lines ...'
        write (*,*)
      end if
!
      n_total=ntss*npss
      n_completed=0
!
!$omp parallel do
!$omp& private(j,k,xfl0,xfl1,bs0,bs1,s,ttb)
!$omp& private(f_trace_reached_boundary,f_br_positive)
!$omp& private(f_trace_on_r0,f_trace_on_r1)
!$omp& private(b_trace_reached_boundary,b_br_positive)
!$omp& private(b_trace_on_r0,b_trace_on_r1)
!$omp& private(nc,diag_step,pct_done)
!$omp& collapse(2)
!$omp& schedule(dynamic,iterations_per_thread)
      do k=1,npss
        do j=1,ntss
!
! ****** Update the iteration counter for diagnostic
! ****** purposes.
!
          if (verbose) then
!$omp critical
            n_completed=n_completed+1
            nc=n_completed
!$omp end critical
          end if
!
          xfl0(1)=rv
          xfl0(2)=tss(j)
          xfl0(3)=pss(k)
!
! ****** Trace a field line in the forward direction along B.
!
          call tracefl (b,ds_f,xfl0,xfl1,bs0,bs1,s,ttb)
!
! ****** Check that the field line reached R0 or R1.
!
          if (ttb) then
            f_trace_reached_boundary=.true.
            f_trace_on_r0=xfl1(1).eq.b%lim0(1)
            f_trace_on_r1=xfl1(1).eq.b%lim1(1)
            f_br_positive=bs1(1).ge.0.
          else
            f_trace_reached_boundary=.false.
          end if
!
! ****** Trace a field line in the backward direction along B.
!
          call tracefl (b,ds_b,xfl0,xfl1,bs0,bs1,s,ttb)
!
! ****** Check that the field line reached R0 or R1.
!
          if (ttb) then
            b_trace_reached_boundary=.true.
            b_trace_on_r0=xfl1(1).eq.b%lim0(1)
            b_trace_on_r1=xfl1(1).eq.b%lim1(1)
            b_br_positive=bs1(1).ge.0.
          else
            b_trace_reached_boundary=.false.
          end if
!
! ****** Set the coronal hole map value.
!
! ****** Note that the following values are set in the output
! ****** coronal hole map:
! ******
! ******    -1: open field line with negative polarity
! ******     1: open field line with positive polarity
! ******     0: closed field line with both footpoints
! ******        on r=R0
! ******     2: closed field line with both footpoints
! ******        on r=R1
! ******    -2: field line that does not reach either
! ******        the r=R0 boundary or the r=R1 boundary
!
          if (f_trace_reached_boundary.and.&
     &        b_trace_reached_boundary) then
            if (f_trace_on_r0.and.b_trace_on_r1) then
              if (f_br_positive) then
                ch(k,j)=one
              else
                ch(k,j)=-one
              end if
            else if (f_trace_on_r1.and.b_trace_on_r0) then
              if (b_br_positive) then
                ch(k,j)=one
              else
                ch(k,j)=-one
              end if
            else if (f_trace_on_r0.and.b_trace_on_r0) then
              ch(k,j)=0.
            else if (f_trace_on_r1.and.b_trace_on_r1) then
              ch(k,j)=two
            else
              ch(k,j)=-two
            end if
          else
            ch(k,j)=-two
          end if
!
! ****** Write progress diagnostics if requested.
!
          if (verbose) then
            diag_step=mod(nc,diagnostic_interval)
            if (diag_step.eq.0) then
              pct_done=100.*nc/n_total
              write (*,910) 'Fraction completed: ',pct_done
  910         format (1x,a,f7.3,'%')
            end if
          end if
!
        enddo
      enddo
!$omp end parallel do
!
! ****** Write the coronal hole map.
!
      if (verbose) then
        write (*,*)
        write (*,*) 'Writing the coronal hole map to file: ',&
     &              trim(ch_map_output_file)
      end if
!
      call wrhdf_2d (ch_map_output_file,.true.,&
     &               npss,ntss,ch,pss,tss,&
     &               hdf32,ierr)
!
      if (ierr.ne.0) then
        write (*,*)
        write (*,*) '### ERROR in GET_CH_MAP:'
        write (*,*) '### Could not write the coronal hole map.'
        call exit (1)
      end if
!
      return
      end
!#######################################################################
      subroutine get_ch_map_3d
!
!-----------------------------------------------------------------------
!
! ****** Compute a 3D coronal hole map.
!
!-----------------------------------------------------------------------
!
      use number_types
      use types
      use mesh
      use field
      use vars
      use files
      use params
      use diags
      use openmp_vars
      use tracefl_interface
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      real(r_typ), parameter :: one=1._r_typ
      real(r_typ), parameter :: two=2._r_typ
!
!-----------------------------------------------------------------------
!
! ****** Storage for the coronal hole map.
!
      real(r_typ), dimension(nrss,ntss,npss) :: ch
!
!-----------------------------------------------------------------------
!
      integer :: ierr,i,j,k
      real(r_typ), dimension(3) :: xfl0,xfl1
      real(r_typ), dimension(3) :: bs0,bs1
      logical :: ttb
      real(r_typ) :: s
      type(flparam) :: ds_f,ds_b
      logical :: f_trace_reached_boundary
      logical :: b_trace_reached_boundary
      logical :: f_trace_on_r0,f_trace_on_r1
      logical :: b_trace_on_r0,b_trace_on_r1
      logical :: f_br_positive
      logical :: b_br_positive
      integer :: n_completed,n_total,nc,diag_step
      real(r_typ) :: pct_done
!
!-----------------------------------------------------------------------
!
      if (verbose) then
        write (*,*)
        write (*,*) '### Computing a 3D coronal hole map:'
      end if
!
! ****** Check that the coronal hole map output file name is not
! ****** blank, since this does not make sense.
!
      if (ch_map_3d_output_file.eq.' ') then
        write (*,*)
        write (*,*) '### ERROR in GET_CH_MAP_3D:'
        write (*,*) '### A coronal hole map was requested, yet'//&
     &              ' the output file name is blank.'
        call exit (1)
      end if
!
! ****** Set the tracing direction to be along the direction
! ****** of the magnetic field.
!
      ds%direction_is_along_b=.true.
!
      ds_f=ds
      ds_f%direction=1
!
      ds_b=ds
      ds_b%direction=-1
!
      if (verbose) then
        write (*,*)
        write (*,*) 'Mapping field lines ...'
        write (*,*)
      end if
!
      n_total=nrss*ntss*npss
      n_completed=0
!
!$omp parallel do
!$omp& private(i,j,k,xfl0,xfl1,bs0,bs1,s,ttb)
!$omp& private(f_trace_reached_boundary,f_br_positive)
!$omp& private(f_trace_on_r0,f_trace_on_r1)
!$omp& private(b_trace_reached_boundary,b_br_positive)
!$omp& private(b_trace_on_r0,b_trace_on_r1)
!$omp& private(nc,diag_step,pct_done)
!$omp& collapse(3)
!$omp& schedule(dynamic,iterations_per_thread)
      do i=1,nrss
        do k=1,npss
          do j=1,ntss
!
! ****** Update the iteration counter for diagnostic
! ****** purposes.
!
            if (verbose) then
!$omp critical
              n_completed=n_completed+1
              nc=n_completed
!$omp end critical
            end if
!
            xfl0(1)=rss(i)
            xfl0(2)=tss(j)
            xfl0(3)=pss(k)
!
! ****** Trace a field line in the forward direction along B.
!
            call tracefl (b,ds_f,xfl0,xfl1,bs0,bs1,s,ttb)
!
! ****** Check that the field line reached R0 or R1.
!
            if (ttb) then
              f_trace_reached_boundary=.true.
              f_trace_on_r0=xfl1(1).eq.b%lim0(1)
              f_trace_on_r1=xfl1(1).eq.b%lim1(1)
              f_br_positive=bs1(1).ge.0.
            else
              f_trace_reached_boundary=.false.
            end if
!
! ****** Trace a field line in the backward direction along B.
!
            call tracefl (b,ds_b,xfl0,xfl1,bs0,bs1,s,ttb)
!
! ****** Check that the field line reached R0 or R1.
!
            if (ttb) then
              b_trace_reached_boundary=.true.
              b_trace_on_r0=xfl1(1).eq.b%lim0(1)
              b_trace_on_r1=xfl1(1).eq.b%lim1(1)
              b_br_positive=bs1(1).ge.0.
            else
              b_trace_reached_boundary=.false.
            end if
!
! ****** Set the coronal hole map value.
!
! ****** Note that the following values are set in the output
! ****** coronal hole map:
! ******
! ******    -1: open field line with negative polarity
! ******     1: open field line with positive polarity
! ******     0: closed field line with both footpoints
! ******        on r=R0
! ******     2: closed field line with both footpoints
! ******        on r=R1
! ******    -2: field line that does not reach either
! ******        the r=R0 boundary or the r=R1 boundary
!
            if (f_trace_reached_boundary.and.&
     &        b_trace_reached_boundary) then
              if (f_trace_on_r0.and.b_trace_on_r1) then
                if (f_br_positive) then
                  ch(i,j,k)=one
                else
                  ch(i,j,k)=-one
                end if
              else if (f_trace_on_r1.and.b_trace_on_r0) then
                if (b_br_positive) then
                  ch(i,j,k)=one
                else
                  ch(i,j,k)=-one
                end if
              else if (f_trace_on_r0.and.b_trace_on_r0) then
                ch(i,j,k)=0.
              else if (f_trace_on_r1.and.b_trace_on_r1) then
                ch(i,j,k)=two
              else
                ch(i,j,k)=-two
              end if
            else
              ch(i,j,k)=-two
            end if
!
! ****** Write progress diagnostics if requested.
!
            if (verbose) then
              diag_step=mod(nc,diagnostic_interval)
              if (diag_step.eq.0) then
                pct_done=100.*nc/n_total
                write (*,910) 'Fraction completed: ',pct_done
  910           format (1x,a,f7.3,'%')
              end if
            end if
!
          enddo
        enddo
      enddo
!$omp end parallel do
!
! ****** Write the coronal hole map.
!
      if (verbose) then
        write (*,*)
        write (*,*) 'Writing the 3D coronal hole map to file: ',&
     &              trim(ch_map_3d_output_file)
      end if
!
      call wrhdf_3d (ch_map_3d_output_file,.true.,&
     &               nrss,ntss,npss,ch,rss,tss,pss,&
     &               hdf32,ierr)
!
      if (ierr.ne.0) then
        write (*,*)
        write (*,*) '### ERROR in GET_CH_MAP_3D:'
        write (*,*) '### Could not write the coronal hole map.'
        call exit (1)
      end if
!
      return
      end
!#######################################################################
      subroutine get_q_on_slice
!
!-----------------------------------------------------------------------
!
! ****** Trace field lines from points on a slice in the 3D volume,
! ****** getting Q at each point.
!
!-----------------------------------------------------------------------
!
      use number_types
      use types
      use mesh
      use field
      use vars
      use files
      use params
      use sds_def
      use diags
      use openmp_vars
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
! ****** Storage for the mapping.
!
      real(r_typ), dimension(:,:,:), allocatable, target :: qfl
      real(r_typ), dimension(:,:,:), allocatable, target :: length
!
!-----------------------------------------------------------------------
!
      type (sds) :: out
      integer :: ierr,i,j,k,n1,n2,n3
      real(r_typ), dimension(3) :: s0
      real(r_typ), dimension(3) :: c
      real(r_typ) :: q
      logical :: gotq
      integer :: n_completed,n_total,nc,diag_step
      real(r_typ) :: pct_done
      logical :: save_field_line_length
      real(r_typ) :: lfl
      integer :: iseq
!
!-----------------------------------------------------------------------
!
      if (verbose) then
        write (*,*)
        write (*,*) '### Computing Q on a slice:'
      end if
!
      n1=slice_c1%dims(1)
      n2=slice_c1%dims(2)
      n3=slice_c1%dims(3)
!
! ****** Set a flag to indicate if the field line length
! ****** was requested.
!
      if (slice_length_output_file.ne.' ') then
        save_field_line_length=.true.
      else
        save_field_line_length=.false.
      end if
!
! ****** Allocate the storage for the output Q array.
!
      allocate (qfl(n1,n2,n3))
!
! ****** Allocate the storage for the output field line
! ****** length, if requested.
!
      if (save_field_line_length) then
        allocate (length(n1,n2,n3))
      end if
!
      if (verbose) then
        write (*,*)
        write (*,*) 'Getting Q ...'
        write (*,*)
      end if
!
! ****** Calculate Q at each point on the slice.
!
      ds%direction_is_along_b=.true.
      ds%direction=1
!
      n_total=n1*n2*n3
      n_completed=0
!
!$omp parallel do
!$omp& private(i,j,k,iseq,c,s0,q,gotq,lfl)
!$omp& private(nc,diag_step,pct_done)
!$omp& collapse(3)
!$omp& schedule(dynamic,iterations_per_thread)
      do k=1,n3
        do j=1,n2
          do i=1,n1
!
! ****** Compute a sequence number that is used to construct
! ****** file names when diagnostic field line traces are
! ****** being output to files.
!
            iseq=i+(j-1)*n1+(k-1)*n1*n2
!
! ****** Update the iteration counter for diagnostic
! ****** purposes.
!
            if (verbose) then
!$omp critical
              n_completed=n_completed+1
              nc=n_completed
!$omp end critical
            end if
!
            if (slice_coords_are_xyz) then
              c=(/slice_c1%f(i,j,k),&
     &            slice_c2%f(i,j,k),&
     &            slice_c3%f(i,j,k)/)
              call c2s (c,s0)
            else
              s0=(/slice_c1%f(i,j,k),&
     &             slice_c2%f(i,j,k),&
     &             slice_c3%f(i,j,k)/)
            end if
!
            call getq (iseq,ds,s0,q_increment_h,q,gotq,lfl)
!
            if (gotq) then
              qfl(i,j,k)=q
            else
              qfl(i,j,k)=0.
            end if
!
            if (save_field_line_length) then
              length(i,j,k)=lfl
            end if
!
! ****** Write progress diagnostics if requested.
!
            if (verbose) then
              diag_step=mod(nc,diagnostic_interval)
              if (diag_step.eq.0) then
                pct_done=100.*nc/n_total
                write (*,910) 'Fraction completed: ',pct_done
  910           format (1x,a,f7.3,'%')
              end if
            end if
!
          enddo
        enddo
      enddo
!$omp end parallel do
!
      if (verbose) then
        write (*,*)
      end if
!
! ****** Write the Q slice.
!
      if (slice_q_output_file.ne.' ') then
        out%ndim=slice_c1%ndim
        out%dims=slice_c1%dims
        out%scale=slice_c1%scale
        out%hdf32=slice_c1%hdf32
        out%scales(1)%f=>slice_c1%scales(1)%f
        out%scales(2)%f=>slice_c1%scales(2)%f
        out%scales(3)%f=>slice_c1%scales(3)%f
        out%f=>qfl
        if (verbose) then
          write (*,*) 'Writing Q in the slice to file: ',&
     &                trim(slice_q_output_file)
        end if
        call wrhdf (slice_q_output_file,out,ierr)
        if (ierr.ne.0) then
          write (*,*)
          write (*,*) '### ERROR in GET_Q_ON_SLICE:'
          write (*,*) '### Could not write Q in the slice'//&
     &                ' to file: ',trim(slice_q_output_file)
          call exit (1)
        end if
      end if
!
      deallocate (qfl)
!
! ****** Write the field line length, if requested.
!
      if (save_field_line_length) then
        out%ndim=slice_c1%ndim
        out%dims=slice_c1%dims
        out%scale=slice_c1%scale
        out%hdf32=slice_c1%hdf32
        out%scales(1)%f=>slice_c1%scales(1)%f
        out%scales(2)%f=>slice_c1%scales(2)%f
        out%scales(3)%f=>slice_c1%scales(3)%f
        out%f=>length
        if (verbose) then
          write (*,*) 'Writing field line length in'//&
     &                ' the slice to file: ',&
     &                trim(slice_length_output_file)
        end if
        call wrhdf (slice_length_output_file,out,ierr)
        if (ierr.ne.0) then
          write (*,*)
          write (*,*) '### ERROR in GET_Q_ON_SLICE:'
          write (*,*) '### Could not write the field line'//&
     &                ' length in the slice to file: ',&
     &                trim(slice_length_output_file)
          call exit (1)
        end if
!
        deallocate (length)
!
      end if
!
      return
      end
!#######################################################################
      subroutine getq (iseq,ds,s0,h,q,valid,lfl)
!
!-----------------------------------------------------------------------
!
! ****** Obtain Q at the point given by the spherical coordinates in
! ****** vector S0 by tracing 5 field lines forwards and backwards
! ****** from S0 to the boundaries.
!
! ****** If Q was successfully obtained, return the Q value in
! ****** variable Q, and VALID=.T.; otherwise, return VALID=.F..
!
! ****** The length of the central field line is returned in LFL.
!
! ****** ISEQ is a sequence number that is used to construct the
! ****** file names when field line traces are being written to
! ****** output files.
!
!-----------------------------------------------------------------------
!
      use number_types
      use types
      use field
      use debug
      use diags
      use tracefl_interface
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      integer :: iseq
      type(flparam) :: ds
      real(r_typ), dimension(3) :: s0
      real(r_typ) :: h
      real(r_typ) :: q
      logical :: valid
      real(r_typ) :: lfl
!
!-----------------------------------------------------------------------
!
      real(r_typ), parameter :: one=1._r_typ
!
!-----------------------------------------------------------------------
!
! ****** Set the limit for how to evaluate Q derivatives near the
! ****** poles.
!
      real(r_typ), parameter :: st_pole_limit_max=5.e-3_r_typ
!
! ****** When sin(theta) of the central field line is smaller than
! ****** ST_POLE_LIMIT_MAX, Cartesian basis vectors are used to
! ****** compute derivatives of the field line mapping; otherwise,
! ****** spherical basis vectors are used.
!
!-----------------------------------------------------------------------
!
! ****** Flag to write warning messages.
!
      logical, parameter :: write_warning_messages=.true.
!
!-----------------------------------------------------------------------
!
      logical :: ttb
      type(flparam) :: ds_local
      real(r_typ), dimension(3) :: x0,x00
      real(r_typ), dimension(3) :: s00
      real(r_typ), dimension(3) :: bs0,bc0
      real(r_typ), dimension(3) :: bs1_f,bs1_b
      real(r_typ) :: s
      real(r_typ) :: bmag,br0_f,br0_b
      real(r_typ), dimension(3) :: bhat,bhat_abs
      real(r_typ) :: r_f,r_b,t_f,t_b
      real(r_typ) :: st_f,st_b
      integer, dimension(1) :: index_min
      real(r_typ), dimension(3) :: e,e1,e2
      real(r_typ), dimension(3) :: s1_f,s1_b
      real(r_typ), dimension(3) :: s1_1p_f,s1_1m_f
      real(r_typ), dimension(3) :: s1_1p_b,s1_1m_b
      real(r_typ), dimension(3) :: s1_2p_f,s1_2m_f
      real(r_typ), dimension(3) :: s1_2p_b,s1_2m_b
      real(r_typ) :: a_f,b_f,c_f,d_f
      real(r_typ) :: a_b,b_b,c_b,d_b
      real(r_typ) :: nsq
      type(csvec) :: xcs
      type(inout) :: outside
      logical :: outside_1p,outside_1m
      logical :: outside_2p,outside_2m
      real(r_typ) :: dx_1p,dx_1m,dx_2p,dx_2m
      real(r_typ) :: x_1p_f,x_1m_f,y_1p_f,y_1m_f
      real(r_typ) :: x_2p_f,x_2m_f,y_2p_f,y_2m_f
      real(r_typ) :: x_1p_b,x_1m_b,y_1p_b,y_1m_b
      real(r_typ) :: x_2p_b,x_2m_b,y_2p_b,y_2m_b
      character(6) :: ch_seq
      logical :: save_trace_points
      logical :: forward_ttb,backward_ttb
!
! ****** Field line trace storage buffer.
!
      type(traj) :: xt
!
!-----------------------------------------------------------------------
!
      real(r_typ), external :: modulo_twopi
      logical, external :: outside_domain
!
!-----------------------------------------------------------------------
!
! ****** Set the flag to save field line traces.
!
      if (debug_level.ge.2) then
        save_trace_points=.true.
        call allocate_trajectory_buffer (xt)
        write (ch_seq,'(i6.6)') iseq
      else
        save_trace_points=.false.
      end if
!
! ****** Initialize Q and the validity flag.

      q=0.
      valid=.false.
      lfl=0.
!
! ****** Initilaize the local DS from that supplied in the
! ****** argument list.
!
      ds_local=ds
!
! ****** Set the flag to interpret the tracing direction as the
! ****** direction along the magnetic field line.
!
      ds_local%direction_is_along_b=.true.
!
      if (debug_level.ge.2) then
        write (*,*)
        write (*,*) '### COMMENT from GETQ:'
        write (*,*) '### Diagnostics for computation of Q:'
        write (*,*) 's0=',s0
      end if
!
! ****** Check that the requested position is not outside the domain.
! ****** If it is, return without calculating a valid Q.
!
      xcs%s=s0
      call sph_to_cart (xcs)
!
      if (outside_domain(b,xcs,outside)) return
!
!-----------------------------------------------------------------------
! ****** Trace the central field line.
!-----------------------------------------------------------------------
!
      s00=s0
!
      ds_local%direction=1
      if (save_trace_points) then
        call tracefl (b,ds_local,s00,s1_f,bs0,bs1_f,s,ttb,xt)
      else
        call tracefl (b,ds_local,s00,s1_f,bs0,bs1_f,s,ttb)
      end if
!
      if (debug_level.ge.2) then
!$omp critical
        call write_trace ('fl_00_f_'//ch_seq//'.dat',xt)
!$omp end critical
      end if
!
      if (debug_level.ge.2) then
        write (*,*)
        write (*,*) '### Central field line, forward trace:'
        write (*,*) 's0=',s00
        write (*,*) 's1=',s1_f
        write (*,*) 'b0=',bs0
        write (*,*) 'b1=',bs1_f
        write (*,*) 'ttb=',ttb
        write (*,*) 's=',s
      end if
!
      forward_ttb=ttb
!
      if (.not.ttb) then
        if (write_warning_messages) then
          write (*,*)
          write (*,*) '### WARNING from GETQ:'
          write (*,*) '### The central field line did not reach'//&
     &                ' the domain boundaries'
          write (*,*) '### during the forward trace.'
          write (*,*) 'Central launch point: ',s0
        end if
        return
      else
        lfl=lfl+s
      end if
!
      ds_local%direction=-1
      if (save_trace_points) then
        call tracefl (b,ds_local,s00,s1_b,bs0,bs1_b,s,ttb,xt)
      else
        call tracefl (b,ds_local,s00,s1_b,bs0,bs1_b,s,ttb)
      end if
!
      if (debug_level.ge.2) then
!$omp critical
        call write_trace ('fl_00_b_'//ch_seq//'.dat',xt)
!$omp end critical
      end if
!
      if (debug_level.ge.2) then
        write (*,*)
        write (*,*) '### Central field line, backward trace:'
        write (*,*) 's0=',s00
        write (*,*) 's1=',s1_b
        write (*,*) 'b0=',bs0
        write (*,*) 'b1=',bs1_b
        write (*,*) 'ttb=',ttb
        write (*,*) 's=',s
      end if
!
      backward_ttb=ttb
!
      if (.not.ttb) then
        if (write_warning_messages) then
          write (*,*)
          write (*,*) '### WARNING from GETQ:'
          write (*,*) '### The central field line did not reach'//&
     &                ' the domain boundaries'
          write (*,*) '### during the backward trace.'
          write (*,*) 'Central launch point: ',s0
        end if
        return
      else
        lfl=lfl+s
      end if
!
! ****** Save the central field line endpoints.
!
      r_f=s1_f(1)
      r_b=s1_b(1)
      t_f=s1_f(2)
      t_b=s1_b(2)
      st_f=sin(t_f)
      st_b=sin(t_b)
!
! ****** Get the radial component of the magnetic field at the
! ****** central field line endpoints.
!
      br0_f=bs1_f(1)
      br0_b=bs1_b(1)
!
! ****** If the field line did not reach both boundaries,
! ****** set the field line length to 0.
!
      if (.not.(forward_ttb.and.backward_ttb)) then
        lfl=0.
      end if
!
      if (debug_level.ge.2) then
        write (*,*)
        write (*,*) '### Central field line:'
        write (*,*) 'r_f=',r_f
        write (*,*) 'r_b=',r_b
        write (*,*) 't_f=',t_f
        write (*,*) 't_b=',t_b
        write (*,*) 'br0_f=',br0_f
        write (*,*) 'br0_b=',br0_b
        write (*,*) 'length=',lfl
      end if
!
! ****** Generate the basis vectors that define the plane of
! ****** the Q computation at the starting location
! ****** (in Cartesian coordinates).
!
! ****** Initialize the position vector at the starting point
! ****** in Cartesian coordinates in X0.
!
      call s2c (s0,x0)
!
! ****** Get the Cartesian magnetic field vector at S0.
!
      call sv_to_cv (s0,bs0,bc0)
!
      bmag=sqrt(bc0(1)**2+bc0(2)**2+bc0(3)**2)
!
      if (debug_level.ge.2) then
        write (*,*) 'bs0=',bs0
        write (*,*) 'bc0=',bc0
        write (*,*) 'bmag=',bmag
      end if
!
! ****** If we hit a null point (B=0), exit with an error.
!
      if (bmag.eq.0.) then
        if (debug_level.ge.2) then
          write (*,*)
          write (*,*) '### WARNING in GETQ:'
          write (*,*) 'B = 0 at the launch point.'
          write (*,*) 'Exiting with an error ...'
        end if
        return
      end if
!
      bhat=bc0/bmag
!
      if (debug_level.ge.2) then
        write (*,*) 'bhat=',bhat
      end if
!
! ****** Select the unit vector that is most perpendicular to B.
!
      bhat_abs=abs(bhat)
      index_min=minloc(bhat_abs)
!
      if (debug_level.ge.2) then
        write (*,*) 'index_min=',index_min
      end if
!
      e=0.
      e(index_min(1))=one
!
      if (debug_level.ge.2) then
        write (*,*) 'e=',e
      end if
!
! ****** The triplet E1, E2, and BHAT form an orthogonal basis.
!
      call normalized_cross_product (e,bhat,e1)
      call normalized_cross_product (e1,bhat,e2)
!
      if (debug_level.ge.2) then
        write (*,*) 'e1=',e1
        write (*,*) 'e2=',e2
      end if
!
!-----------------------------------------------------------------------
! ****** Trace the field line at X0+H*E1.
!-----------------------------------------------------------------------
!
      x00=x0+h*e1
!
      call c2s (x00,s00)
!
! ****** If the launch point is outside the domain, use the central
! ****** field line to take a one-sided derivative.
!
      xcs%s=s00
      xcs%c=x00
!
      if (outside_domain(b,xcs,outside)) then
        outside_1p=.true.
        dx_1p=0.
        s1_1p_f=s1_f
        s1_1p_b=s1_b
      else
        outside_1p=.false.
        dx_1p=h
      end if
!
      if (.not.outside_1p) then
!
        ds_local%direction=1
        if (save_trace_points) then
          call tracefl (b,ds_local,s00,s1_1p_f,bs0,bs1_f,s,ttb,xt)
        else
          call tracefl (b,ds_local,s00,s1_1p_f,bs0,bs1_f,s,ttb)
        end if
!
        if (debug_level.ge.3) then
!$omp critical
          call write_trace ('fl_p0_f_'//ch_seq//'.dat',xt)
!$omp end critical
        end if
!
        if (debug_level.ge.3) then
          write (*,*)
          write (*,*) '### x0+h*e1 field line, forward trace:'
          write (*,*) 's0=',s00
          write (*,*) 's1=',s1_1p_f
          write (*,*) 'ttb=',ttb
          write (*,*) 's=',s
        end if
!
        if (.not.ttb) then
          if (write_warning_messages) then
            write (*,*)
            write (*,*) '### WARNING from GETQ:'
            write (*,*) '### The x0+h*e1 field line did not reach'//&
     &                  ' the domain boundaries'
            write (*,*) '### during the forward trace.'
            write (*,*) 'Central launch point: ',s0
            write (*,*) 'Launch point: ',s00
          end if
          return
        end if
!
        ds_local%direction=-1
        if (save_trace_points) then
          call tracefl (b,ds_local,s00,s1_1p_b,bs0,bs1_b,s,ttb,xt)
        else
          call tracefl (b,ds_local,s00,s1_1p_b,bs0,bs1_b,s,ttb)
        end if
!
        if (debug_level.ge.3) then
!$omp critical
          call write_trace ('fl_p0_b_'//ch_seq//'.dat',xt)
!$omp end critical
        end if
!
        if (debug_level.ge.3) then
          write (*,*)
          write (*,*) '### x0+h*e1 field line, backward trace:'
          write (*,*) 's0=',s00
          write (*,*) 's1=',s1_1p_b
          write (*,*) 'ttb=',ttb
          write (*,*) 's=',s
        end if
!
        if (.not.ttb) then
          if (write_warning_messages) then
            write (*,*)
            write (*,*) '### WARNING from GETQ:'
            write (*,*) '### The x0+h*e1 field line did not reach'//&
     &                  ' the domain boundaries'
            write (*,*) '### during the backward trace.'
            write (*,*) 'Central launch point: ',s0
            write (*,*) 'Launch point: ',s00
          end if
          return
        end if
!
      end if
!
!-----------------------------------------------------------------------
! ****** Trace the field line at X0-H*E1.
!-----------------------------------------------------------------------
!
      x00=x0-h*e1
!
      call c2s (x00,s00)
!
! ****** If the launch point is outside the domain, use the central
! ****** field line to take a one-sided derivative.
!
      xcs%s=s00
      xcs%c=x00
!
      if (outside_domain(b,xcs,outside)) then
        outside_1m=.true.
        dx_1m=0.
        s1_1m_f=s1_f
        s1_1m_b=s1_b
      else
        outside_1m=.false.
        dx_1m=h
      end if
!
! ****** If both the plus and minus perturbed launch points
! ****** are outside the domain, return without calculating a
! ****** valid Q.  Write a warning if this happens.
!
      if (outside_1p.and.outside_1m) then
        if (write_warning_messages) then
          write (*,*)
          write (*,*) '### WARNING from GETQ:'
          write (*,*) '### The x0+h*e1 and x0-h*e1 launch points'//&
     &                ' are both outside the domain.'
          write (*,*) 'Central launch point: ',s0
          x00=x0+h*e1
          call c2s (x00,s00)
          write (*,*) 'Launch point x0+h*e1: ',s00
          x00=x0-h*e1
          call c2s (x00,s00)
          write (*,*) 'Launch point x0-h*e1: ',s00
        end if
        return
      end if
!
      if (.not.outside_1m) then
!
        ds_local%direction=1
        if (save_trace_points) then
          call tracefl (b,ds_local,s00,s1_1m_f,bs0,bs1_f,s,ttb,xt)
        else
          call tracefl (b,ds_local,s00,s1_1m_f,bs0,bs1_f,s,ttb)
        end if
!
        if (debug_level.ge.3) then
!$omp critical
          call write_trace ('fl_m0_f_'//ch_seq//'.dat',xt)
!$omp end critical
        end if
!
        if (debug_level.ge.3) then
          write (*,*)
          write (*,*) '### x0-h*e1 field line, forward trace:'
          write (*,*) 's0=',s00
          write (*,*) 's1=',s1_1m_f
          write (*,*) 'ttb=',ttb
          write (*,*) 's=',s
        end if
!
        if (.not.ttb) then
          if (write_warning_messages) then
            write (*,*)
            write (*,*) '### WARNING from GETQ:'
            write (*,*) '### The x0-h*e1 field line did not reach'//&
     &                  ' the domain boundaries'
            write (*,*) '### during the forward trace.'
            write (*,*) 'Central launch point: ',s0
            write (*,*) 'Launch point: ',s00
          end if
          return
        end if
!
        ds_local%direction=-1
        if (save_trace_points) then
          call tracefl (b,ds_local,s00,s1_1m_b,bs0,bs1_b,s,ttb,xt)
        else
          call tracefl (b,ds_local,s00,s1_1m_b,bs0,bs1_b,s,ttb)
        end if
!
        if (debug_level.ge.3) then
!$omp critical
          call write_trace ('fl_m0_b_'//ch_seq//'.dat',xt)
!$omp end critical
        end if
!
        if (debug_level.ge.3) then
          write (*,*)
          write (*,*) '### x0-h*e1 field line, backward trace:'
          write (*,*) 's0=',s00
          write (*,*) 's1=',s1_1m_b
          write (*,*) 'ttb=',ttb
          write (*,*) 's=',s
        end if
!
        if (.not.ttb) then
          if (write_warning_messages) then
            write (*,*)
            write (*,*) '### WARNING from GETQ:'
            write (*,*) '### The x0-h*e1 field line did not reach'//&
     &                  ' the domain boundaries'
            write (*,*) '### during the backward trace.'
            write (*,*) 'Central launch point: ',s0
            write (*,*) 'Launch point: ',s00
          end if
          return
        end if
!
      end if
!
!-----------------------------------------------------------------------
! ****** Trace the field line at X0+H*E2.
!-----------------------------------------------------------------------
!
      x00=x0+h*e2
!
      call c2s (x00,s00)
!
! ****** If the launch point is outside the domain, use the central
! ****** field line to take a one-sided derivative.
!
      xcs%s=s00
      xcs%c=x00
!
      if (outside_domain(b,xcs,outside)) then
        outside_2p=.true.
        dx_2p=0.
        s1_2p_f=s1_f
        s1_2p_b=s1_b
      else
        outside_2p=.false.
        dx_2p=h
      end if
!
      if (.not.outside_2p) then
!
        ds_local%direction=1
        if (save_trace_points) then
          call tracefl (b,ds_local,s00,s1_2p_f,bs0,bs1_f,s,ttb,xt)
        else
          call tracefl (b,ds_local,s00,s1_2p_f,bs0,bs1_f,s,ttb)
        end if
!
        if (debug_level.ge.3) then
!$omp critical
          call write_trace ('fl_0p_f_'//ch_seq//'.dat',xt)
!$omp end critical
        end if
!
        if (debug_level.ge.3) then
          write (*,*)
          write (*,*) '### x0+h*e2 field line, forward trace:'
          write (*,*) 's0=',s00
          write (*,*) 's1=',s1_2p_f
          write (*,*) 'ttb=',ttb
          write (*,*) 's=',s
        end if
!
        if (.not.ttb) then
          if (write_warning_messages) then
            write (*,*)
            write (*,*) '### WARNING from GETQ:'
            write (*,*) '### The x0+h*e2 field line did not reach'//&
     &                  ' the domain boundaries'
            write (*,*) '### during the forward trace.'
            write (*,*) 'Central launch point: ',s0
            write (*,*) 'Launch point: ',s00
          end if
          return
        end if
!
        ds_local%direction=-1
        if (save_trace_points) then
          call tracefl (b,ds_local,s00,s1_2p_b,bs0,bs1_b,s,ttb,xt)
        else
          call tracefl (b,ds_local,s00,s1_2p_b,bs0,bs1_b,s,ttb)
        end if
!
        if (debug_level.ge.3) then
!$omp critical
          call write_trace ('fl_0p_b_'//ch_seq//'.dat',xt)
!$omp end critical
        end if
!
        if (debug_level.ge.3) then
          write (*,*)
          write (*,*) '### x0+h*e2 field line, backward trace:'
          write (*,*) 's0=',s00
          write (*,*) 's1=',s1_2p_b
          write (*,*) 'ttb=',ttb
          write (*,*) 's=',s
        end if
!
        if (.not.ttb) then
          if (write_warning_messages) then
            write (*,*)
            write (*,*) '### WARNING from GETQ:'
            write (*,*) '### The x0+h*e2 field line did not reach'//&
     &                  ' the domain boundaries'
            write (*,*) '### during the backward trace.'
            write (*,*) 'Central launch point: ',s0
            write (*,*) 'Launch point: ',s00
          end if
          return
        end if
!
      end if
!
!-----------------------------------------------------------------------
! ****** Trace the field line at X0-H*E2.
!-----------------------------------------------------------------------
!
      x00=x0-h*e2
!
      call c2s (x00,s00)
!
! ****** If the launch point is outside the domain, use the central
! ****** field line to take a one-sided derivative.
!
      xcs%s=s00
      xcs%c=x00
!
      if (outside_domain(b,xcs,outside)) then
        outside_2m=.true.
        dx_2m=0.
        s1_2m_f=s1_f
        s1_2m_b=s1_b
      else
        outside_2m=.false.
        dx_2m=h
      end if
!
! ****** If both the plus and minus perturbed launch points
! ****** are outside the domain, return without calculating a
! ****** valid Q.  This should never happen.
!
      if (outside_2p.and.outside_2m) then
        if (write_warning_messages) then
          write (*,*)
          write (*,*) '### WARNING from GETQ:'
          write (*,*) '### The x0+h*e2 and x0-h*e2 launch points'//&
     &                ' are both outside the domain.'
          write (*,*) 'Central launch point: ',s0
          x00=x0+h*e2
          call c2s (x00,s00)
          write (*,*) 'Launch point x0+h*e2: ',s00
          x00=x0-h*e2
          call c2s (x00,s00)
          write (*,*) 'Launch point x0-h*e2: ',s00
        end if
        return
      end if
!
      if (.not.outside_2m) then
!
        ds_local%direction=1
        if (save_trace_points) then
          call tracefl (b,ds_local,s00,s1_2m_f,bs0,bs1_f,s,ttb,xt)
        else
          call tracefl (b,ds_local,s00,s1_2m_f,bs0,bs1_f,s,ttb)
        end if
!
        if (debug_level.ge.3) then
!$omp critical
          call write_trace ('fl_0m_f_'//ch_seq//'.dat',xt)
!$omp end critical
        end if
!
        if (debug_level.ge.3) then
          write (*,*)
          write (*,*) '### x0-h*e2 field line, forward trace:'
          write (*,*) 's0=',s00
          write (*,*) 's1=',s1_2m_f
          write (*,*) 'ttb=',ttb
          write (*,*) 's=',s
        end if
!
        if (.not.ttb) then
          if (write_warning_messages) then
            write (*,*)
            write (*,*) '### WARNING from GETQ:'
            write (*,*) '### The x0-h*e2 field line did not reach'//&
     &                  ' the domain boundaries'
            write (*,*) '### during the forward trace.'
            write (*,*) 'Central launch point: Spherical s0 = ',s0
            write (*,*) 'Launch point: Spherical x0-h*e2 = ',s00
          end if
          return
        end if
!
        ds_local%direction=-1
        if (save_trace_points) then
          call tracefl (b,ds_local,s00,s1_2m_b,bs0,bs1_b,s,ttb,xt)
        else
          call tracefl (b,ds_local,s00,s1_2m_b,bs0,bs1_b,s,ttb)
        end if
!
        if (debug_level.ge.3) then
!$omp critical
          call write_trace ('fl_0m_b_'//ch_seq//'.dat',xt)
!$omp end critical
        end if
!
        if (debug_level.ge.3) then
          write (*,*)
          write (*,*) '### x0-h*e2 field line, backward trace:'
          write (*,*) 's0=',s00
          write (*,*) 's1=',s1_2m_b
          write (*,*) 'ttb=',ttb
          write (*,*) 's=',s
        end if
!
        if (.not.ttb) then
          if (write_warning_messages) then
            write (*,*)
            write (*,*) '### WARNING from GETQ:'
            write (*,*) '### The x0-h*e2 field line did not reach'//&
     &                  ' the domain boundaries'
            write (*,*) '### during the backward trace.'
            write (*,*) 'Central launch point: Spherical s0 = ',s0
            write (*,*) 'Launch point: Spherical x0-h*e2 = ',s00
          end if
          return
        end if
!
      end if
!
! ****** Get the value of the Jacobian matrix coefficients.
!
      if (st_f.lt.st_pole_limit_max) then
        x_1p_f=r_f*sin(s1_1p_f(2))*cos(s1_1p_f(3))
        x_1m_f=r_f*sin(s1_1m_f(2))*cos(s1_1m_f(3))
        y_1p_f=r_f*sin(s1_1p_f(2))*sin(s1_1p_f(3))
        y_1m_f=r_f*sin(s1_1m_f(2))*sin(s1_1m_f(3))
        x_2p_f=r_f*sin(s1_2p_f(2))*cos(s1_2p_f(3))
        x_2m_f=r_f*sin(s1_2m_f(2))*cos(s1_2m_f(3))
        y_2p_f=r_f*sin(s1_2p_f(2))*sin(s1_2p_f(3))
        y_2m_f=r_f*sin(s1_2m_f(2))*sin(s1_2m_f(3))
        a_f=(x_1p_f-x_1m_f)/(dx_1p+dx_1m)
        b_f=(x_2p_f-x_2m_f)/(dx_2p+dx_2m)
        c_f=(y_1p_f-y_1m_f)/(dx_1p+dx_1m)
        d_f=(y_2p_f-y_2m_f)/(dx_2p+dx_2m)
      else
        a_f=modulo_twopi(s1_1p_f(3)-s1_1m_f(3))*r_f*st_f/(dx_1p+dx_1m)
        b_f=modulo_twopi(s1_2p_f(3)-s1_2m_f(3))*r_f*st_f/(dx_2p+dx_2m)
        c_f=(s1_1p_f(2)-s1_1m_f(2))*r_f/(dx_1p+dx_1m)
        d_f=(s1_2p_f(2)-s1_2m_f(2))*r_f/(dx_2p+dx_2m)
      end if
!
      if (st_b.lt.st_pole_limit_max) then
        x_1p_b=r_b*sin(s1_1p_b(2))*cos(s1_1p_b(3))
        x_1m_b=r_b*sin(s1_1m_b(2))*cos(s1_1m_b(3))
        y_1p_b=r_b*sin(s1_1p_b(2))*sin(s1_1p_b(3))
        y_1m_b=r_b*sin(s1_1m_b(2))*sin(s1_1m_b(3))
        x_2p_b=r_b*sin(s1_2p_b(2))*cos(s1_2p_b(3))
        x_2m_b=r_b*sin(s1_2m_b(2))*cos(s1_2m_b(3))
        y_2p_b=r_b*sin(s1_2p_b(2))*sin(s1_2p_b(3))
        y_2m_b=r_b*sin(s1_2m_b(2))*sin(s1_2m_b(3))
        a_b=(x_1p_b-x_1m_b)/(dx_1p+dx_1m)
        b_b=(x_2p_b-x_2m_b)/(dx_2p+dx_2m)
        c_b=(y_1p_b-y_1m_b)/(dx_1p+dx_1m)
        d_b=(y_2p_b-y_2m_b)/(dx_2p+dx_2m)
      else
        a_b=modulo_twopi(s1_1p_b(3)-s1_1m_b(3))*r_b*st_b/(dx_1p+dx_1m)
        b_b=modulo_twopi(s1_2p_b(3)-s1_2m_b(3))*r_b*st_b/(dx_2p+dx_2m)
        c_b=(s1_1p_b(2)-s1_1m_b(2))*r_b/(dx_1p+dx_1m)
        d_b=(s1_2p_b(2)-s1_2m_b(2))*r_b/(dx_2p+dx_2m)
      end if
!
      if (debug_level.ge.2) then
        write (*,*)
        write (*,*) '### Jacobian matrix coefficients:'
        write (*,*) 'a_f=',a_f
        write (*,*) 'b_f=',b_f
        write (*,*) 'c_f=',c_f
        write (*,*) 'd_f=',d_f
        write (*,*) 'a_b=',a_b
        write (*,*) 'b_b=',b_b
        write (*,*) 'c_b=',c_b
        write (*,*) 'd_b=',d_b
      end if
!
      if (debug_level.ge.2) then
        write (*,*)
        write (*,*) '### Div B conservation check:'
        write (*,*) 'abs(Det(D_f))*abs(br0_f)/bmag = ',&
     &              abs(a_f*d_f-b_f*c_f)*abs(br0_f)/bmag
        write (*,*) 'abs(Det(D_b))*abs(br0_b)/bmag = ',&
     &              abs(a_b*d_b-b_b*c_b)*abs(br0_b)/bmag
      end if
!
! ****** Get the value of Q.
!
      nsq= (a_f*d_b-b_f*c_b)**2&
     &    +(a_b*b_f-a_f*b_b)**2&
     &    +(d_b*c_f-d_f*c_b)**2&
     &    +(a_b*d_f-b_b*c_f)**2
!
      q=nsq*abs(br0_f*br0_b)/bmag**2
!
      if (debug_level.ge.2) then
        write (*,*)
        write (*,*) '### Final Q value:'
        write (*,*) 'q=',q
      end if
!
      valid=.true.
!
! ****** Deallocate the storage for the field line trace buffer
! ****** if it was used.
!
      if (save_trace_points) then
        call deallocate_trajectory_buffer (xt)
      end if
!
      return
      end
!#######################################################################
      subroutine normalized_cross_product (a,b,c)
!
!-----------------------------------------------------------------------
!
! ****** Return the unit vector C = (A x B)/|A x B|.
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
      real(r_typ), dimension(3) :: a,b,c
!
!-----------------------------------------------------------------------
!
      real(r_typ) :: cnorm
!
!-----------------------------------------------------------------------
!
! ****** Set C to the cross-product of A and B.
!
      c(1)=a(2)*b(3)-a(3)*b(2)
      c(2)=a(3)*b(1)-a(1)*b(3)
      c(3)=a(1)*b(2)-a(2)*b(1)
!
! ****** Normalize C to unit length.
!
      cnorm=sqrt(c(1)**2+c(2)**2+c(3)**2)
!
      if (cnorm.ne.0.) then
        c=c/cnorm
      end if
!
      return
      end
!#######################################################################
      subroutine write_trace (fname,xt)
!
!-----------------------------------------------------------------------
!
! ****** Write the field line trace in structure XT to the
! ****** text file named FNAME.
!
!-----------------------------------------------------------------------
!
      use types
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      character(*) :: fname
      type(traj) :: xt
!
!-----------------------------------------------------------------------
!
      character, parameter :: TAB=achar(9)
!
!-----------------------------------------------------------------------
!
      integer :: i,ierr
!
!-----------------------------------------------------------------------
!
      call ffopen (1,fname,'rw',ierr)
!
      if (ierr.ne.0) then
        write (*,*)
        write (*,*) '### ERROR in WRITE_TRACE:'
        write (*,*) '### Could not create a field line output'//&
     &              ' file.'
        write (*,*) 'File name: ',trim(fname)
        call exit (1)
      end if
!
! ****** Write the output coordinates.
!
      write (1,'(5a)') 'r',TAB,'t',TAB,'p'
      do i=1,xt%npts
        write (1,'(3(1pe23.16,a))') xt%x(1)%f(i),TAB,&
     &                              xt%x(2)%f(i),TAB,&
     &                              xt%x(3)%f(i)
      enddo
!
      close (1)
!
      return
      end
!#######################################################################
      subroutine tracefl (b,ds,s0,s1,bs0,bs1,s,&
     &                    traced_to_r_boundary,xt)
!
!-----------------------------------------------------------------------
!
! ****** Trace a magnetic field line.
!
!-----------------------------------------------------------------------
!
! ****** The 3D magnetic field is specified by structure B.
! ****** The structure DS has the field line integration parameters,
! ****** and S0 contains the spherical coordinates of the launch
! ****** point.
!
! ****** TRACED_TO_R_BOUNDARY=.T. is returned if the field line
! ****** was traced all the way to a radial domain boundary, in
! ****** which case S1 contains the spherical coordinates of the
! ****** final location, and S has the traced field line length.
!
! ****** Otherwise, TRACED_TO_R_BOUNDARY=.F. is returned, and
! ****** S and S1 do not necessarily have valid values.
! ****** This also occurs if |B|=0 is encountered during the trace.
!
! ****** The magnetic field vectors in spherical coordinates
! ****** at the starting and ending footpoints, respectively,
! ****** are returned in BS0 and BS1.
!
! ****** If the field line trace is needed, pass in the optional
! ****** field line trace buffer XT.  The buffer XT needs to be
! ****** allocated prior to being passed in to this routine.
!
!-----------------------------------------------------------------------
!
      use number_types
      use types
      use step_size_stats
      use debug
      use integrate_fl
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      type(vec) :: b
      type(flparam) :: ds
      real(r_typ), dimension(3) :: s0,s1
      real(r_typ), dimension(3) :: bs0,bs1
      real(r_typ) :: s
      logical :: traced_to_r_boundary
      type(traj), optional :: xt
!
      intent(in) :: b,ds,s0
      intent(out) :: s1,bs0,bs1,s,traced_to_r_boundary
!
!-----------------------------------------------------------------------
!
      real(r_typ), parameter :: one=1._r_typ
      real(r_typ), parameter :: half=.5_r_typ
!
!-----------------------------------------------------------------------
!
      logical :: store_trace
      real(r_typ) :: ds0,dss,dsss,frac,dsmult_corrector
      real(r_typ) :: sf=1._r_typ,sf1,sf2
      logical :: done_tracing,first,nullb
      integer :: idir0,n,ntry,max_n,max_ntry
      type (csvec) :: x,xp,xo,bv,bhat1,bhat2
      type (inout) :: outside
      integer :: ierr
      real(r_typ) :: arc_length_remaining,ds_b
      type(flparam) :: current_ds
      logical :: tfc
!
      integer :: local_stat_n
      real(r_typ) :: local_stat_ds_sum,&
     &               local_stat_ds_min,&
     &               local_stat_ds_max
!
!-----------------------------------------------------------------------
!
      logical, external :: outside_domain
!
!-----------------------------------------------------------------------
!
      if (debug_level.ge.4) then
        write (*,*)
        write (*,*) '### COMMENT from TRACEFL:'
        write (*,*) '### Starting a new trace:'
        write (*,*) 'S0 = ',S0
      end if
!
! ****** Initializations.
!
      store_trace=present(xt)
!
      current_ds=ds
      first=.true.
      tfc=.false.
      s1=s0
      max_n=floor(ds%lmax/ds%min)
      max_ntry=ds%short_fl_max_tries
!
      if (gather_stats) then
        local_stat_n=0
        local_stat_ds_sum=0.
        local_stat_ds_min=huge(one)
        local_stat_ds_max=0.
      end if
!
! ****** If storage for the field line trace was requested,
! ****** initialize the field line trace buffer.
!
      if (store_trace) xt%npts=0
!
      do ntry=1,max_ntry
!
!-----------------------------------------------------------------------
! ****** Trace a field line starting at launch point S0.
!-----------------------------------------------------------------------
!
        done_tracing=.false.
!
! ****** If things go wrong, return TRACED_TO_R_BOUNDARY=.F..
!
        traced_to_r_boundary=.false.
!
        x%s=s0
        call sph_to_cart (x)
        s=0.
        if (store_trace) call add_trajectory_point (xt,x%s)
!
! ****** Set the starting step size.
!
        ds0=current_ds%min
!
        if (debug_level.ge.4) then
          write (*,*)
          write (*,*) 'Start of trace:'
          write (*,*) 'X = ',x%s
        end if
!
! ****** If the initial point is outside the domain, stop tracing.
!
        if (outside_domain(b,x,outside)) tfc=.true.

        do n=1,max_n
!
! ****** If file line is done tracing (including restarts), exit loop.
!
          if(tfc) exit
!
!-----------------------------------------------------------------------
! ****** Trace to the next step.
!-----------------------------------------------------------------------
!
! ****** Set the tracing direction.
!
          call getb (b,x,bv)
!
! ****** On the first time in, set the tracing direction in
! ****** variable IDIR0 based on that specified by DS%DIRECTION
! ****** and DS%DIRECTION_IS_ALONG_B.  Store the magnetic field at
! ****** the initial location (in spherical coordinates) in BS0.
!
! ****** When DS%DIRECTION_IS_ALONG_B=.F., DS%DIRECTION=1 traces
! ****** in the direction of increasing r, and DS%DIRECTION=-1
! ****** traces in the opposite direction.
!
! ****** When DS%DIRECTION_IS_ALONG_B=.T., DS%DIRECTION=1 traces
! ****** in the direction of the magnetic field vector, and
! ****** DS%DIRECTION=-1 traces in the opposite direction.
!
          if (first) then
            first=.false.
            bs0=bv%s
            if (ds%direction.gt.0) then
              idir0=1
            else
              idir0=-1
            end if
            if (.not.ds%direction_is_along_b) then
              if (bv%s(1).lt.0.) idir0=-idir0
            end if
          end if
!
          call normalize_v (bv,nullb)
          if (nullb) then
            write (*,*)
            write (*,*) '### WARNING from TRACEFL:'
            write (*,*) '### The trace encountered a null pnt (B = 0).'
            write (*,*) '### This occurred at the start of the trace.'
            write (*,*) '### Abandoning the trace ...'
            write (*,*) 'Location (r,t,p) = ',x%s
            tfc=.true.
            exit
          end if
!
          if (n.gt.1) then
            ds_b=abs(dsss)
            bhat1=bhat2
          end if
          bhat2=bv
!
! ****** If a variable step size is being used, set the step
! ****** size for the next step.
!
! ****** If this is a repeat trace (of a short field line,
! ****** NTRY.gt.1), then leave the step size at the minimum
! ****** value until the number of points exceeds the minimum
! ****** allowed number.
!
          if (ds%variable.and.n.gt.1) then
            if (.not.(ntry.gt.1.and.n.le.ds%short_fl_min_points)) then
              call get_ds (b,x,bhat1,bhat2,ds_b,current_ds,ds0)
            end if
          end if
!
! ****** Check to see if this trace segment ends the trace
! ****** (i.e., exceeds the arc length specified, DS%LMAX).
!
          arc_length_remaining=ds%lmax-s
!
          if (ds0.ge.arc_length_remaining) then
            dss=arc_length_remaining
            done_tracing=.true.
          else
            dss=ds0
          end if
!
          if (dss.le.0.) then
            tfc=.true.
            exit
          end if
!
! ****** Gather step size statistics.
!
          if (gather_stats) then
            local_stat_n=local_stat_n+1
            local_stat_ds_sum=local_stat_ds_sum+abs(ds0)
            local_stat_ds_min=min(local_stat_ds_min,abs(ds0))
            local_stat_ds_max=max(local_stat_ds_max,abs(ds0))
          end if
!
!-----------------------------------------------------------------------
! ****** Predictor.
!-----------------------------------------------------------------------
!
! ****** Advance for a half-step to achieve second-order accuracy.
!
          dsmult_corrector=one
          xo=x
          xp=x
!
          if (debug_level.ge.4) then
            write (*,*)
            write (*,*) 'Predictor:'
            write (*,*) 'N = ',n
            write (*,*) 'BV = ',bv%s
          end if
!
          dsss=half*idir0*dss
          call advance (xp,bv,dsss)
!
          if (debug_level.ge.4) then
            write (*,*) 'DSSS = ',dsss
            write (*,*) 'XP = ',xp%s
          end if
!
! ****** Check if the field line has exited the domain.
!
          if (outside_domain(b,xp,outside)) then
!
            if (outside%t.or.outside%p) then
!
! ****** The field line has crossed the theta or phi boundaries;
! ****** stop tracing.
!
              tfc=.true.
              exit
!
            else if (outside%r) then
!
! ****** The field line has crossed the r boundary.
!
              if (debug_level.ge.4) then
                write (*,*) '### Predictor: Outside r domain:'
              end if
!
! ****** Get the clip fraction, FRAC, that clips the segment
! ****** to the r boundary.
!
              call get_r_clip_fraction (b,xo,xp,outside%r0,frac,ierr)
!
              if (debug_level.ge.4) then
                write (*,*) 'After GET_R_CLIP_FRACTION (predictor):'
                write (*,*) 'IERR = ',ierr
                write (*,*) 'FRAC = ',frac
              end if
!
! ****** If there is an error in getting the clip fraction,
! ****** stop tracing.  This should never happen: write
! ****** detailed debugging information.
!
              if (ierr.ne.0) then
                write (*,*)
                write (*,*) '### ANOMALY in TRACEFL:'
                write (*,*) '### Predictor, 1st step:'
                write (*,*) '### Could not get the clip fraction.'
                write (*,*)
                write (*,*) '### Debugging info:'
                write (*,*) 'S0 = ',s0
                write (*,*) 'DS%DIRECTION_IS_ALONG_B = ',&
     &                      ds%direction_is_along_b
                write (*,*) 'DS%DIRECTION = ',ds%direction
                write (*,*) 'DS%MIN = ',current_ds%min
                write (*,*) 'DS%MAX = ',current_ds%max
                write (*,*) 'DSSS = ',dsss
                write (*,*) 'BV = ',bv
                write (*,*) 'XO = ',xo
                write (*,*) 'XP = ',xp
                tfc=.true.
                exit
              end if
!
              if (frac.eq.0..and.n.eq.1) then
!
! ****** The initial point is exactly on the boundary, and is
! ****** being traced out of the boundary.  We are done.
!
                traced_to_r_boundary=.true.
                done_tracing=.true.
                tfc=.true.
                exit
              end if
!
              if (frac.le.ds%predictor_min_clip_fraction) then
!
! ****** The starting point is close to the boundary.  Clip the
! ****** final point to the radial boundary and stop tracing.
!
                call clip_to_r_boundary (b,xo,xp,outside%r0,frac)
!
                if (debug_level.ge.4) then
                  write (*,*) 'After CLIP_TO_R_BOUNDARY (predictor):'
                  write (*,*) 'XP = ',xp%s
                end if
!
                x=xp
                dsss=frac*dsss
                traced_to_r_boundary=.true.
                done_tracing=.true.
                if (do_integral_along_fl) then
                   call getsf (xo,sf1)
                   call getsf (xp,sf2)
                   sf=half*(sf1+sf2)
                endif
                s=s+abs(dsss)*sf
                exit
!
              end if
!
! ****** If the starting point is not close enough to the radial
! ****** boundary, predict again half way to the radial boundary
! ****** to achieve second-order accuracy in the corrector.
!
              dsss=half*frac*dsss
              xp=x
              call advance (xp,bv,dsss)
!
              if (debug_level.ge.4) then
                write (*,*) 'After 2nd predictor:'
                write (*,*) 'XP = ',xp%s
              end if
!
! ****** Check if the predicted point has exited the domain.
! ****** This should never happen: write detailed debugging
! ****** information.
!
              if (outside_domain(b,xp,outside)) then
                write (*,*)
                write (*,*) '### ANOMALY in TRACEFL:'
                write (*,*) '### Predictor, 2nd step:'
                write (*,*) '### Point is outside the boundary.'
                write (*,*)
                write (*,*) '### Debugging info:'
                write (*,*) 'S0 = ',s0
                write (*,*) 'DS%DIRECTION_IS_ALONG_B = ',&
     &                  ds%direction_is_along_b
                write (*,*) 'DS%DIRECTION = ',ds%direction
                write (*,*) 'DS%MIN = ',current_ds%min
                write (*,*) 'DS%MAX = ',current_ds%max
                write (*,*) 'DSSS = ',dsss
                write (*,*) 'BV = ',bv
                write (*,*) 'XO = ',xo
                write (*,*) 'XP = ',xp
                tfc=.true.
                exit
              end if
!
! ****** Reduce the step size for the corrector.  This minimizes the
! ****** chance of "going through the radial boundary and back into
! ****** the domain again", which can happen when the step size
! ****** is too big.
!
              dsmult_corrector=half*frac
!
            end if
!
          end if
!
!-----------------------------------------------------------------------
! ****** Corrector.
!-----------------------------------------------------------------------
!
          call getb (b,xp,bv)
          call normalize_v (bv,nullb)
          if (nullb) then
            write (*,*)
            write (*,*) '### WARNING from TRACEFL:'
            write (*,*) '### The trace encountered a null pnt (B = 0).'
            write (*,*) '### This occurred during the corrector.'
            write (*,*) '### Abandoning the trace ...'
            write (*,*) 'Location (r,t,p) = ',xp%s
            exit
          end if
!
          dsss=idir0*dsmult_corrector*dss
          call advance (x,bv,dsss)
!
          if (debug_level.ge.4) then
            write (*,*) 'After corrector advance:'
            write (*,*) 'BV = ',bv%s
            write (*,*) 'DSSS = ',dsss
            write (*,*) 'X = ',x%s
          end if
!
          if (outside_domain(b,x,outside)) then
            if (outside%t.or.outside%p) then
              tfc=.true.
              exit
            else if (outside%r) then
!
              if (debug_level.ge.4) then
                write (*,*) '### Corrector: Outside r domain:'
              end if
!
              call get_r_clip_fraction (b,xo,x,outside%r0,frac,ierr)
!
! ****** If there is an error in getting the clip fraction,
! ****** stop tracing.  This should never happen: write
! ****** detailed debugging information.
!
              if (ierr.ne.0) then
                write (*,*)
                write (*,*) '### ANOMALY in TRACEFL:'
                write (*,*) '### Corrector step:'
                write (*,*) '### Could not get the clip fraction.'
                write (*,*)
                write (*,*) '### Debugging info:'
                write (*,*) 'S0 = ',s0
                write (*,*) 'DS%DIRECTION_IS_ALONG_B = ',&
     &                      ds%direction_is_along_b
                write (*,*) 'DS%DIRECTION = ',ds%direction
                write (*,*) 'DS%MIN = ',current_ds%min
                write (*,*) 'DS%MAX = ',current_ds%max
                write (*,*) 'DSSS = ',dsss
                write (*,*) 'BV = ',bv
                write (*,*) 'XO = ',xo
                write (*,*) 'X  = ',x
                tfc=.true.
                exit
              end if
!
              if (debug_level.ge.4) then
                write (*,*) 'After GET_R_CLIP_FRACTION (corrector):'
                write (*,*) 'IERR = ',ierr
                write (*,*) 'FRAC = ',frac
              end if
!
              done_tracing=.true.
              traced_to_r_boundary=.true.
!
              if (frac.eq.0.) then
!
! ****** This should only happen when the previous point was exactly
! ****** on the boundary, and is being traced out of the domain.
! ****** In this case, do not store this point: we are done.
!
                x=xo
!
                if (debug_level.ge.4) then
                  write (*,*) 'The trace went from the boundary to'//&
     &                        ' the outside (corrector):'
                  write (*,*) 'X = ',x%s
                end if
!
                exit
!
              else
!
                call clip_to_r_boundary (b,xo,x,outside%r0,frac)
                dsss=frac*dsss
!
                if (debug_level.ge.4) then
                  write (*,*) 'After CLIP_TO_R_BOUNDARY (corrector):'
                  write (*,*) 'X = ',x%s
                end if
!
              end if
            end if
          end if
!
! ****** Increment the number of points and the arc length.
!
          if (do_integral_along_fl) then
            call getsf (xo,sf1)
            call getsf (x,sf2)
            sf=half*(sf1+sf2)
          endif
          s=s+abs(dsss)*sf
!
! ****** Add the current position to the field line buffer
! ****** if requested.
!
          if (store_trace) call add_trajectory_point (xt,x%s)
!
! ****** Break out of max_n loop if done.
!
          if (done_tracing) exit
!
        enddo !n->max_n
!
! ****** Finished tracing the field line.
!
! ****** Break out of outer loop if field line complete.
!
        if (tfc) exit
!
! ****** Check that the number of points in the field line trace
! ****** exceeds the minimum allowed, DS%SHORT_FL_MIN_POINTS.
! ****** If not, reduce the step size and retrace the field line,
! ****** up to a maximum of DS%SHORT_FL_MAX_TRIES times.
!
        if (n.lt.ds%short_fl_min_points) then
!
          if (debug_level.ge.4) then
            write (*,*) 'Short field line:'
            write (*,*) 'N = ',n
            write (*,*) 'CURRENT_DS%OVER_RC = ',current_ds%over_rc
            write (*,*) 'CURRENT_DS%MIN = ',current_ds%min
            write (*,*) 'CURRENT_DS%MAX = ',current_ds%max
          end if
!
          current_ds%over_rc=current_ds%over_rc*&
     &                       ds%short_fl_shrink_factor
          current_ds%min=current_ds%min*&
     &                   ds%short_fl_shrink_factor
          current_ds%max=current_ds%max*&
     &                   ds%short_fl_shrink_factor
!
        else
          exit !Break out of max_ntry loop if line had enough points.
        end if
!
      end do !ntry->max_ntry
!
      if((ntry.ge.max_ntry) .and. (n.le.ds%short_fl_min_points)) then
        write (*,*)
        write (*,*) '### WARNING from TRACEFL:'
        write (*,*) 'Short field line after max tries:'
        write (*,*) 'Number of points in field line = ',n
        write (*,*) 'Number of tries = ',ntry
        write (*,*) 'CURRENT_DS%OVER_RC = ',current_ds%over_rc
        write (*,*) 'CURRENT_DS%MIN = ',current_ds%min
        write (*,*) 'CURRENT_DS%MAX = ',current_ds%max
        write (*,*) 'S0 = ',s0
        write (*,*) 'S1 = ',x%s
      end if
!
! ****** Update the step size statistics.
!
      if (gather_stats) then
!$omp critical (omp_stat)
        stat_n=stat_n+local_stat_n
        stat_ds_sum=stat_ds_sum+local_stat_ds_sum
        stat_ds_min=min(stat_ds_min,local_stat_ds_min)
        stat_ds_max=max(stat_ds_max,local_stat_ds_max)
!$omp end critical (omp_stat)
      end if
!
! ****** Store the final location in S1.
!
      s1=x%s
!
! ****** If the field line was not traced to the radial boundary,
! ****** do not attempt to return the magnetic field vector
! ****** at the field line endpoint.
!
      if (.not.traced_to_r_boundary) then
        bs1=0.
        return
      end if
!
! ****** Store the magnetic field vector (in spherical coordinates)
! ****** at the field line endpoint in BS1.
!
      call getb (b,x,bv)
      bs1=bv%s
!
      if (debug_level.ge.4) then
        write (*,*)
        write (*,*) '### COMMENT from TRACEFL:'
        write (*,*) '### About to exit:'
        write (*,*) 'N = ',n
        write (*,*) 'S = ',s
        write (*,*) 'S1 = ',s1
        write (*,*) 'BS0 = ',bs0
        write (*,*) 'BS1 = ',bs1
      end if
!
      return
      end
!#######################################################################
      subroutine get_ds (b,x,v0,v1,ds_v,ds,deltas)
!
!-----------------------------------------------------------------------
!
! ****** Set the integration step size based on the radius of
! ****** curvature of the field line, as estimated from the unit
! ****** vectors V0 and V1, and the local mesh cell size of the
! ****** magnetic field B at X (if requested).
!
! ****** The vectors V0 and V1 are assumed to have been evaluated
! ****** along the field line a distance DS_V apart.
!
! ****** On input, DELTAS should have the present value of the
! ****** step size.  On return, DELTAS is overwritten by the new
! ****** step size.
!
! ****** It is assumed that DS_V and DELTAS are positive.
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
      type(vec) :: b
      type(csvec) :: x
      type(csvec) :: v0,v1
      real(r_typ) :: ds_v
      type(flparam) :: ds
      real(r_typ) :: deltas
!
!-----------------------------------------------------------------------
!
      real(r_typ) :: dv,factor,ds_mesh
!
!-----------------------------------------------------------------------
!
! ****** If DELTAS is zero, set it to the minimum value.
!
      if (deltas.eq.0.) then
        deltas=ds%min
        return
      end if
!
! ****** First, set the step size based on the local radius
! ****** of curvature of the field line.
!
! ****** Compute the factor by which DELTAS must be multiplied to
! ****** achieve the specified ratio of step size to radius
! ****** of curvature, DS%OVER_RC.
!
      dv=sqrt( (v0%c(1)-v1%c(1))**2&
     &        +(v0%c(2)-v1%c(2))**2&
     &        +(v0%c(3)-v1%c(3))**2)
!
      if (dv.ne.0.) then
        factor=ds_v*ds%over_rc/(dv*deltas)
      else
        factor=huge(dv)
      end if
!
! ****** Limit the change in DELTAS to the maximum permitted
! ****** change per step.
!
      factor=min(factor,ds%max_increase_factor)
      factor=max(factor,ds%max_decrease_factor)
!
! ****** Set the new step size.
!
      deltas=factor*deltas
!
! ****** Next, if requested, limit DELTAS by the local mesh
! ****** cell size of B at X.
!
! ****** Note that this only needs to be done if DELTAS is bigger
! ****** than DS%MIN, since only then is there a possibility of
! ****** reducing DELTAS further.
!
      if (deltas.gt.ds%min.and.ds%limit_by_local_mesh) then
!
! ****** Get the local mesh cell size at X.
!
        call get_local_mesh_size (b,x%s,ds_mesh)
!
! ****** Set DELTAS so that it does not exceed DS%LOCAL_MESH_FACTOR
! ****** times the local mesh size.
!
        deltas=min(deltas,ds%local_mesh_factor*ds_mesh)
!
      end if
!
! ****** Limit DELTAS by the maximum and mimimum allowed values.
!
      deltas=max(deltas,ds%min)
      deltas=min(deltas,ds%max)
!
      return
      end
!#######################################################################
      subroutine get_local_mesh_size (b,s,ds)
!
!-----------------------------------------------------------------------
!
! ****** Get the local mesh size DS from the magnetic field in
! ****** structure B at the spherical position S.
!
!-----------------------------------------------------------------------
!
      use number_types
      use types
      use interp_interface
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      type(vec) :: b
      real(r_typ), dimension(3) :: s
      real(r_typ) :: ds
!
!-----------------------------------------------------------------------
!
      real(r_typ), parameter :: one=1._r_typ
!
!-----------------------------------------------------------------------
!
      integer :: i,j,k,ip1,jp1,kp1
      real(r_typ) :: ar,at,ap,drv,dtv,dpv,stv
!
!-----------------------------------------------------------------------
!
! ****** Get the local size of the main mesh of B at the
! ****** specified point.
!
      call interp (b%nrs,b%rs,s(1),i,ip1,ar,b%rs_invtab)
      drv=(one-ar)*b%drs(i)+ar*b%drs(ip1)
!
      call interp (b%nts,b%ts,s(2),j,jp1,at,b%ts_invtab)
      dtv=(one-at)*b%dts(j)+at*b%dts(jp1)
      stv=(one-at)*b%sts(j)+at*b%sts(jp1)
!
      call interp (b%nps,b%ps,s(3),k,kp1,ap,b%ps_invtab)
      dpv=(one-ap)*b%dps(k)+ap*b%dps(kp1)
!
      ds=min(drv,s(1)*dtv,s(1)*stv*dpv)
!
      return
      end
!#######################################################################
      subroutine normalize_v (v,null)
!
!-----------------------------------------------------------------------
!
! ****** Normalize the vector V to return a unit vector along V.
! ****** If V has zero length, then set NULL=.T. and leave V
! ****** unchanged; otherwise, set NULL=.F..
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
      type(csvec) :: v
      logical :: null
!
!-----------------------------------------------------------------------
!
      real(r_typ) :: vmag
!
!-----------------------------------------------------------------------
!
! ****** Use the spherical representation to compute the norm.
!
      vmag=sqrt(v%s(1)**2+v%s(2)**2+v%s(3)**2)
!
      if (vmag.eq.0.) then
        null=.true.
      else
        null=.false.
        v%c=v%c/vmag
        v%s=v%s/vmag
      end if
!
      return
      end
!#######################################################################
      subroutine advance (x,v,ds)
!
!-----------------------------------------------------------------------
!
! ****** Advance the position vector X by the step DS using the
! ****** velocity vector V.
!
! ****** This routine updates both Cartesian and spherical
! ****** represenations of X in the dual represenations position
! ****** vector X.
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
      type(csvec) :: x,v
      real(r_typ) :: ds
!
!-----------------------------------------------------------------------
!
! ****** Advance the Cartesian position.
!
      x%c=x%c+ds*v%c
!
! ****** Transform the new Cartesian position to spherical
! ****** coordinates.
!
      call cart_to_sph (x)
!
      return
      end
!#######################################################################
      subroutine get_r_clip_fraction (b,x0,x1,outside_r0,frac,ierr)
!
!-----------------------------------------------------------------------
!
! ****** Get the fraction FRAC (between 0 and 1) that expresses
! ****** the normalized distance between X0 and X1
! ****** corresponding to the location of the radial boundary.
!
! ****** The boundary location is obtained from the 3D magnetic
! ****** field in structure B.
!
! ****** It is assumed that X0 and X1 lie on different sides
! ****** of the r boundary specified by flag OUTSIDE_R0.
!
! ****** FRAC can be used to clip the position to the
! ****** radial boundary.
!
! ****** For a normal return, IERR=0 is returned when it was possible
! ****** to estimate FRAC; in the case of an inconsistency,
! ****** IERR=1 is returned and FRAC is invalid.
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
      type(vec) :: b
      type(csvec) :: x0,x1
      logical :: outside_r0
      real(r_typ) :: frac
      integer :: ierr
!
!-----------------------------------------------------------------------
!
      real(r_typ), parameter :: one=1._r_typ
      real(r_typ), parameter :: two=2._r_typ
      real(r_typ), parameter :: half=.5_r_typ
!
!-----------------------------------------------------------------------
!
      real(r_typ) :: rb,dr0,dr1,dssq,x0dotb
      real(r_typ) :: eps,term1,term2,ratio,disc
      integer :: ipm
!
!-----------------------------------------------------------------------
!
      ierr=0
!
      if (outside_r0) then
!
! ****** The segment crossed the inner radial boundary.
!
        rb=b%lim0(1)
        ipm=-1
!
      else
!
! ****** The segment crossed the outer radial boundary.
!
        rb=b%lim1(1)
        ipm=1
!
      end if
!
! ****** Check that X0 and X1 are consistent with a radial
! ****** boundary crossing.
!
      dr0=rb-x0%s(1)
      dr1=rb-x1%s(1)
!
      if (dr0.eq.0..and.dr1.eq.0.) then
        ierr=1
        return
      end if
!
      if (dr0*dr1.gt.0.) then
        ierr=1
        return
      end if
!
! ****** Treat the special case when X0 is exactly on the boundary.
!
      if (dr0.eq.0.) then
        frac=0.
        return
      end if
!
      dssq= (x1%c(1)-x0%c(1))**2&
     &     +(x1%c(2)-x0%c(2))**2&
     &     +(x1%c(3)-x0%c(3))**2
!
      if (dssq.le.0.) then
        ierr=1
        return
      end if
!
      x0dotb= x0%c(1)*(x1%c(1)-x0%c(1))&
     &       +x0%c(2)*(x1%c(2)-x0%c(2))&
     &       +x0%c(3)*(x1%c(3)-x0%c(3))
!
! ****** Get the square root of the discriminant, expanding small
! ****** arguments for accuracy.
!
      eps=10*sqrt(spacing(one))
!
      term1=x0dotb**2
      term2=dssq*dr0*(two*x0%s(1)+dr0)
      if (term1+term2.lt.0.) then
        ierr=1
        return
      end if
      if (term1.ne.0.) then
        ratio=term2/term1
        if (abs(ratio).lt.eps) then
          disc=sqrt(term1)*(one+half*ratio)
        else
          disc=sqrt(term1+term2)
        end if
      else
        disc=sqrt(term1+term2)
      end if
!
      frac=(-x0dotb+ipm*disc)/dssq
!
      return
      end
!#######################################################################
      subroutine clip_to_r_boundary (b,x0,x1,outside_r0,frac)
!
!-----------------------------------------------------------------------
!
! ****** Reset the position X1 to the normalized distance FRAC
! ****** between X0 and X1.
!
! ****** The boundary location is obtained from the 3D magnetic
! ****** field in structure B.
!
! ****** The flag OUTSIDE_R0=.T. indicates that the lower
! ****** radial boundary was crossed in going from X0 to X1;
! ****** otherwise, the upper radial boundary was crossed.
!
! ****** In conjunction with GET_R_CLIP_FRACTION this routine
! ****** can be used to clip points that cross the radial boundaries
! ****** to the boundary.
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
      type(vec) :: b
      type(csvec) :: x0,x1
      logical :: outside_r0
      real(r_typ) :: frac
!
!-----------------------------------------------------------------------
!
      real(r_typ), parameter :: one=1._r_typ
!
!-----------------------------------------------------------------------
!
      real(r_typ) :: rval
!
!-----------------------------------------------------------------------
!
! ****** Reset the Cartesian position in X1.
!
      x1%c=(one-frac)*x0%c+frac*x1%c
!
! ****** Transform the new Cartesian position to spherical
! ****** coordinates.
!
      call cart_to_sph (x1)
!
! ****** Set the radius to the appropriate radial boundary
! ****** value (exactly).  This is done to take care of roundoff.
!
      if (outside_r0) then
        rval=b%lim0(1)
      else
        rval=b%lim1(1)
      end if
!
      x1%s(1)=rval
!
! ****** Transform the new spherical position to Cartesian
! ****** coordinates.
!
      call sph_to_cart (x1)
!
      return
      end
!#######################################################################
      function outside_domain (b,x,outside)
!
!-----------------------------------------------------------------------
!
! ****** If the spherical position in X lies outside the limits
! ****** of the domain, return a function result of .T.;
! ****** otherwise, return .F..
!
! ****** The detailed in/out location for each coordinate
! ****** is set in structure OUTSIDE.
!
!-----------------------------------------------------------------------
!
      use types
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      type(vec) :: b
      type(csvec) :: x
      type(inout) :: outside
      logical :: outside_domain
!
!-----------------------------------------------------------------------
!
      real(r_typ), parameter :: zero=0.
      real(r_typ), parameter :: twopi=6.2831853071795864_r_typ
!
!-----------------------------------------------------------------------
!
      real(r_typ) :: pv
!
!-----------------------------------------------------------------------
!
      real(r_typ), external :: fold
!
!-----------------------------------------------------------------------
!
! ****** Get the value of phi in the main interval [0,2*pi].
!
      pv=fold(zero,twopi,x%s(3))
!
      outside%r0=x%s(1).lt.b%lim0(1)
      outside%t0=x%s(2).lt.b%lim0(2)
      outside%p0=pv.lt.b%lim0(3)
      outside%r1=x%s(1).gt.b%lim1(1)
      outside%t1=x%s(2).gt.b%lim1(2)
      outside%p1=pv.gt.b%lim1(3)
!
      outside%r=outside%r0.or.outside%r1
      outside%t=outside%t0.or.outside%t1
      outside%p=outside%p0.or.outside%p1
!
      outside%domain=outside%r.or.outside%t.or.outside%p
!
      outside_domain=outside%domain
!
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
      subroutine cart_to_sph (x)
!
!-----------------------------------------------------------------------
!
! ****** Update the dual representation position vector X so that
! ****** the Cartesian position corresponds to the spherical position.
!
!-----------------------------------------------------------------------
!
      use types
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      type(csvec) :: x
!
!-----------------------------------------------------------------------
!
      call c2s (x%c,x%s)
!
      return
      end
!#######################################################################
      subroutine sph_to_cart (x)
!
!-----------------------------------------------------------------------
!
! ****** Update the dual representation position vector X so that
! ****** the spherical position corresponds to the Cartesian position.
!
!-----------------------------------------------------------------------
!
      use types
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      type(csvec) :: x
!
!-----------------------------------------------------------------------
!
      call s2c (x%s,x%c)
!
      return
      end
!#######################################################################
      subroutine c2s (x,s)
!
!-----------------------------------------------------------------------
!
! ****** Convert the vector X = (x,y,z) from Cartesian coordinates
! ****** to spherical coordinates S = (r,t,p).
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
      real(r_typ), dimension(3) :: x
      real(r_typ), dimension(3) :: s
!
!-----------------------------------------------------------------------
!
      real(r_typ) :: r,t,p
!
!-----------------------------------------------------------------------
!
      r=sqrt(x(1)**2+x(2)**2+x(3)**2)
!
      if (r.eq.0.) then
        t=0.
      else
        t=acos(x(3)/r)
      end if
!
      if (x(1).eq.0.) then
        if (x(2).ge.0.) then
          p= halfpi
        else
          p=-halfpi
        end if
      else
        p=atan2(x(2),x(1))
      end if
      if (p.lt.0.) p=p+twopi
!
      s(1)=r
      s(2)=t
      s(3)=p
!
      return
      end
!#######################################################################
      subroutine s2c (s,x)
!
!-----------------------------------------------------------------------
!
! ****** Convert the vector S = (r,t,p) from spherical coordinates
! ****** to Cartesian coordinates X = (x,y,z).
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
      real(r_typ), dimension(3) :: s
      real(r_typ), dimension(3) :: x
!
!-----------------------------------------------------------------------
!
      real(r_typ) :: rst
!
!-----------------------------------------------------------------------
!
      rst=s(1)*sin(s(2))
      x(1)=rst*cos(s(3))
      x(2)=rst*sin(s(3))
      x(3)=s(1)*cos(s(2))
!
      return
      end
!#######################################################################
      subroutine cv_to_sv (s,cv,sv)
!
!-----------------------------------------------------------------------
!
! ****** Convert the Cartesian vector CV to spherical vector SV
! ****** at spherical position S.
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
      real(r_typ), dimension(3) :: s,cv,sv
      intent(in) :: s,cv
      intent(out) :: sv
!
!-----------------------------------------------------------------------
!
      real(r_typ) :: st,ct,sp,cp
!
!-----------------------------------------------------------------------
!
      st=sin(s(2))
      ct=cos(s(2))
      sp=sin(s(3))
      cp=cos(s(3))
!
      sv(1)= cv(1)*st*cp+cv(2)*st*sp+cv(3)*ct
      sv(2)= cv(1)*ct*cp+cv(2)*ct*sp-cv(3)*st
      sv(3)=-cv(1)*sp   +cv(2)*cp
!
      return
      end
!#######################################################################
      subroutine sv_to_cv (s,sv,cv)
!
!-----------------------------------------------------------------------
!
! ****** Convert the spherical vector SV to Cartesian vector CV
! ****** at spherical position S.
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
      real(r_typ), dimension(3) :: s,sv,cv
      intent(in) :: s,sv
      intent(out) :: cv
!
!-----------------------------------------------------------------------
!
      real(r_typ) :: st,ct,sp,cp
!
!-----------------------------------------------------------------------
!
      st=sin(s(2))
      ct=cos(s(2))
      sp=sin(s(3))
      cp=cos(s(3))
!
      cv(1)= sv(1)*st*cp+sv(2)*ct*cp-sv(3)*sp
      cv(2)= sv(1)*st*sp+sv(2)*ct*sp+sv(3)*cp
      cv(3)= sv(1)*ct   -sv(2)*st
!
      return
      end
!#######################################################################
      subroutine getb (b,x,bv)
!
!-----------------------------------------------------------------------
!
! ****** Get the interpolated magnetic field vector BV at the
! ****** dual representation position vector X from the
! ****** magnetic field vector in structure B.
!
! ****** When USE_ANALYTIC_FUNCTION=.TRUE., get the magnetic
! ****** field components by calling a function rather
! ****** than using B.
!
!-----------------------------------------------------------------------
!
      use number_types
      use types
      use evaluate_spline_3d_interface
      use vars
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      type(vec) :: b
      type(csvec) :: x,bv
!
!-----------------------------------------------------------------------
!
      real(r_typ), parameter :: zero=0.
      real(r_typ), parameter :: twopi=6.2831853071795864_r_typ
!
!-----------------------------------------------------------------------
!
      real(r_typ) :: pv
      real(r_typ) :: br,bt,bp,st,ct,sp,cp
      real(r_typ), dimension(3) :: xs
!
!-----------------------------------------------------------------------
!
      real(r_typ), external :: fold
!
!-----------------------------------------------------------------------
!
! ****** Get the value of phi in the main interval [0,2*pi].
!
      pv=fold(zero,twopi,x%s(3))
!
! ****** If we are using an analytic function to define the
! ****** magnetic field, call it.
!
      if (use_analytic_function) then
        xs(1)=x%s(1)
        xs(2)=x%s(2)
        xs(3)=pv
        call magnetic_field_function (0.,.true.,xs,bv%s)
        go to 100
      end if
!
! ****** Get Br, Bt, and Bp.
!
      if (b%cubic) then
!
! ****** Cubic spline interpolation.
!
        br=evaluate_spline_3d(b%spl%r,x%s(1),x%s(2),pv,&
     &                        b%inv(1)%c(1),&
     &                        b%inv(1)%c(2),&
     &                        b%inv(1)%c(3))
!
        bt=evaluate_spline_3d(b%spl%t,x%s(1),x%s(2),pv,&
     &                        b%inv(2)%c(1),&
     &                        b%inv(2)%c(2),&
     &                        b%inv(2)%c(3))
!
        bp=evaluate_spline_3d(b%spl%p,x%s(1),x%s(2),pv,&
     &                        b%inv(3)%c(1),&
     &                        b%inv(3)%c(2),&
     &                        b%inv(3)%c(3))
!
      else
!
! ****** Linear interpolation.
!
        call interp_3d (b%r%dims(1),b%r%dims(2),b%r%dims(3),&
     &                  b%r%scales(1)%f,&
     &                  b%r%scales(2)%f,&
     &                  b%r%scales(3)%f,&
     &                  b%inv(1),&
     &                  b%r%f,x%s(1),x%s(2),pv,br)
!
        call interp_3d (b%t%dims(1),b%t%dims(2),b%t%dims(3),&
     &                  b%t%scales(1)%f,&
     &                  b%t%scales(2)%f,&
     &                  b%t%scales(3)%f,&
     &                  b%inv(2),&
     &                  b%t%f,x%s(1),x%s(2),pv,bt)
!
        call interp_3d (b%p%dims(1),b%p%dims(2),b%p%dims(3),&
     &                  b%p%scales(1)%f,&
     &                  b%p%scales(2)%f,&
     &                  b%p%scales(3)%f,&
     &                  b%inv(3),&
     &                  b%p%f,x%s(1),x%s(2),pv,bp)
!
      end if
!
      bv%s(1)=br
      bv%s(2)=bt
      bv%s(3)=bp
!
  100 continue
!
! ****** Transform the spherical components of the magnetic field
! ****** to the Cartesian components.
!
      st=sin(x%s(2))
      ct=cos(x%s(2))
      sp=sin(pv)
      cp=cos(pv)
      bv%c(1)=bv%s(1)*st*cp+bv%s(2)*ct*cp-bv%s(3)*sp
      bv%c(2)=bv%s(1)*st*sp+bv%s(2)*ct*sp+bv%s(3)*cp
      bv%c(3)=bv%s(1)*ct   -bv%s(2)*st
!
      return
      end
!#######################################################################
      subroutine getsf (x,sf)
!
!-----------------------------------------------------------------------
!
! ****** Get the interpolated scalar field SF at the
! ****** dual representation position vector X from the
! ****** field SCALAR_FIELD.
!
!-----------------------------------------------------------------------
!
      use number_types
      use types
      use vars
      use integrate_fl
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      type(csvec) :: x
      real(r_typ) :: sf
!
!-----------------------------------------------------------------------
!
      real(r_typ), parameter :: zero=0.
      real(r_typ), parameter :: twopi=6.2831853071795864_r_typ
!
!-----------------------------------------------------------------------
!
      real(r_typ) :: pv
      real(r_typ), dimension(3) :: xs
!
!-----------------------------------------------------------------------
!
      real(r_typ), external :: fold
!
!-----------------------------------------------------------------------
!
! ****** Get the value of phi in the main interval [0,2*pi].
!
      pv=fold(zero,twopi,x%s(3))
!
! ****** Linear interpolation.
!
        call interp_3d (scalar_field%dims(1),&
     &                  scalar_field%dims(2),&
     &                  scalar_field%dims(3),&
     &                  scalar_field%scales(1)%f,&
     &                  scalar_field%scales(2)%f,&
     &                  scalar_field%scales(3)%f,&
     &                  inv_sf,&
     &                  scalar_field%f,x%s(1),x%s(2),pv,sf)

!
      return
      end
!#######################################################################
      subroutine interp_3d (nx,ny,nz,x,y,z,inv,f,xv,yv,zv,fv)
!
!-----------------------------------------------------------------------
!
! ****** Interpolate the value of the 3D field FV at (XV,YV,ZV) from
! ****** array F(NX,NY,NZ), defined on the mesh X(NX) x Y(NY) x Z(NZ).
! ****** The structure INV holds the inverse interpolation tables.
!
! ****** Note that if the point (XV,YV,ZV) is outside the bounds of
! ****** the X x Y x Z mesh, FV=0. is returned.
!
!-----------------------------------------------------------------------
!
      use number_types
      use types
      use interp_interface
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
!
!-----------------------------------------------------------------------
!
      real(r_typ), parameter :: one=1._r_typ
!
!-----------------------------------------------------------------------
!
      integer :: i,j,k,ip1,jp1,kp1
      real(r_typ) :: ax,ay,az
!
!-----------------------------------------------------------------------
!
! ****** If the point is outside the data limits, return a
! ****** zero value.
!
      if (xv.lt.x(1).or.xv.gt.x(nx).or.&
     &    yv.lt.y(1).or.yv.gt.y(ny).or.&
     &    zv.lt.z(1).or.zv.gt.z(nz)) then
        fv=0.
        return
      end if
!
      call interp (nx,x,xv,i,ip1,ax,inv%c(1))
      call interp (ny,y,yv,j,jp1,ay,inv%c(2))
      call interp (nz,z,zv,k,kp1,az,inv%c(3))
!
      fv= (one-ax)*( (one-ay)*( (one-az)*f(i  ,j  ,k  )&
     &                         +     az *f(i  ,j  ,kp1))&
     &              +     ay *( (one-az)*f(i  ,jp1,k  )&
     &                         +     az *f(i  ,jp1,kp1)))&
     &   +     ax *( (one-ay)*( (one-az)*f(ip1,j  ,k  )&
     &                         +     az *f(ip1,j  ,kp1))&
     &              +     ay *( (one-az)*f(ip1,jp1,k  )&
     &                         +     az *f(ip1,jp1,kp1)))
!
      return
      end
!#######################################################################
      subroutine allocate_trajectory_buffer (xt)
!
!-----------------------------------------------------------------------
!
! ****** Allocate the trajectory buffer XT.
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
      type(traj) :: xt
!
!-----------------------------------------------------------------------
!
      integer :: i
!
!-----------------------------------------------------------------------
!
      xt%ndim=3
      xt%size=xt%initial_size
!
! ****** Allocate storage for the trajectory buffer.
!
      allocate (xt%x(xt%ndim))
!
      do i=1,xt%ndim
        allocate (xt%x(i)%f(xt%size))
      enddo
!
! ****** Initialize the current number of points in the
! ****** trajectory buffer.
!
      xt%npts=0
!
      return
      end
!#######################################################################
      subroutine deallocate_trajectory_buffer (xt)
!
!-----------------------------------------------------------------------
!
! ****** Deallocate the trajectory buffer XT.
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
      type(traj) :: xt
!
!-----------------------------------------------------------------------
!
      integer :: i
!
!-----------------------------------------------------------------------
!
      if (.not.associated(xt%x)) return
!
      do i=1,xt%ndim
        if (associated(xt%x(i)%f)) then
          deallocate (xt%x(i)%f)
        end if
      enddo
!
      deallocate (xt%x)
!
      xt%npts=0
!
      return
      end
!#######################################################################
      subroutine add_trajectory_point (xt,x)
!
!-----------------------------------------------------------------------
!
! ****** Add the position vector X to the trajectory buffer XT.
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
      type(traj) :: xt
      real(r_typ), dimension(xt%ndim) :: x
!
!-----------------------------------------------------------------------
!
      integer :: i,n
!
!-----------------------------------------------------------------------
!
! ****** Increment the number of points in the trajectory.
! ****** If the buffer is full, expand it.
!
      n=xt%npts
      n=n+1
      if (n.gt.xt%size) call expand_trajectory_buffer (xt)
!
! ****** Add the point to the buffer.
!
      do i=1,xt%ndim
        xt%x(i)%f(n)=x(i)
      enddo
!
      xt%npts=n
!
      return
      end
!#######################################################################
      subroutine expand_trajectory_buffer (xt)
!
!-----------------------------------------------------------------------
!
! ****** Expand the trajectory buffer XT by doubling the number
! ****** of points in the buffer.
!
!-----------------------------------------------------------------------
!
      use number_types
      use types
      use debug
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      type(traj) :: xt
!
!-----------------------------------------------------------------------
!
      real(r_typ), dimension(:), pointer :: f
      integer :: i,n
!
!-----------------------------------------------------------------------
!
! ****** Double the current buffer size, and copy the current
! ****** contents to the expanded buffer.
!
      n=xt%size
      do i=1,xt%ndim
        allocate (f(2*n))
        f(1:n)=xt%x(i)%f(1:n)
        deallocate (xt%x(i)%f)
        xt%x(i)%f=>f
      enddo
      xt%size=2*n
      if (debug_level.ge.5) then
        write (*,*) 'Expanded a trajectory buffer to ',xt%size
      end if
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
! ****** table, NU.
!
! ****** The output is a structure TAB with the inverse interpolation
! ****** table.  This structure has the following components:
!
! ******    N:  the number of points in the table;
! ******    D:  the inverse of the uniform table spacing;
! ******    F:  the inverse interpolation table.
!
!-----------------------------------------------------------------------
!
      use number_types
      use invint_def
      use interp_interface
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
      integer :: i,k,ip1
      real(r_typ) :: dx,xv,alpha,en
!
!-----------------------------------------------------------------------
!
! ****** Check that the number of points is valid.
!
      if (tab%n.le.1) then
        write (*,*)
        write (*,*) '### ERROR in GETINV:'
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
        write (*,*) '### ERROR in GETINV:'
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
        call interp (n,x,xv,i,ip1,alpha)
        tab%f(k)=i+alpha
        tab%f(k)=max(tab%f(k),one)
        tab%f(k)=min(tab%f(k),en)
      enddo
!
      return
      end
!#######################################################################
      subroutine interp (n,x,xv,i,ip1,alpha,tab)
!
!-----------------------------------------------------------------------
!
! ****** Find the interval I in table X(i), i=1,2,...,N, that encloses
! ****** the value XV, such that X(I).le.XV.le.X(I+1).
! ****** For the special case when N=1, XV must equal X(1) exactly.
!
! ****** This routine uses LOCATE_INTERVAL (from the SPLINE library)
! ****** to do the actual work.  If the interval is not found
! ****** LOCATE_INTERVAL terminates with an error.
!
! ****** This routine does not do the actual interpolation.  However,
! ****** the returned values of I, IP1 (which generally equals I+1),
! ****** and ALPHA can be used to get the interpolant.
!
! ****** The optional inverse interpolation table, TAB, can be
! ****** supplied to improve the efficiency of the search.
!
!-----------------------------------------------------------------------
!
      use number_types
      use invint_def
      use locate_interval_interface
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
      integer :: i
      integer :: ip1
      real(r_typ) :: alpha
      type(itab), optional :: tab
      intent(in) :: n,x,xv,tab
      intent(out) :: i,ip1,alpha
!
!-----------------------------------------------------------------------
!
      if (present(tab)) then
        i=locate_interval(n,x,xv,tab)
      else
        i=locate_interval(n,x,xv)
      end if
!
      if (n.eq.1) then
        ip1=1
        alpha=0.
      else
        ip1=i+1
        if (x(i).eq.x(i+1)) then
          alpha=0.
        else
          alpha=(xv-x(i))/(x(i+1)-x(i))
        end if
      end if
!
      return
      end
!#######################################################################
      subroutine set_parameters
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
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
! ****** Storage the for usage line.
!
      type(paragraph), pointer :: usage
!
! ****** Storage for the error message.
!
      character(72) :: errmsg
!
!-----------------------------------------------------------------------
!
      integer :: ierr
      character(256) :: arg
      logical :: set
!
!-----------------------------------------------------------------------
!
! ****** Define the syntax.
!
      call defarg (GROUP_K ,'-v',' ',' ')
      call defarg (GROUP_A ,'infile',' ',' ')
!
! ****** Parse the command line.
!
      call parse (errmsg,ierr)
!
      if (ierr.ne.0) then
!
        write (*,*)
        write (*,*) '### ',cname,' Version ',cvers,' of ',cdate,'.'
        write (*,*) '### Calculate the field line mapping'//&
     &              ' for MAS code runs.'
!
        if (ierr.gt.1) then
          write (*,*)
          write (*,*) errmsg
        end if
!
! ****** Print the usage line.
!
        call get_usage_line (usage)
!
        write (*,*)
        write (*,*) 'Usage:'
        write (*,*)
!
        call print_par (usage)
!
        write (*,*)
        write (*,*) 'Read the parameters from input file <infile>.'
        call delete_par (usage)
!
        call exit (1)
!
      end if
!
! ****** Set the parameters.
!
! ****** Verbose flag.
!
      call fetcharg ('-v',set,arg)
      verbose=set
!
! ****** Input file name.
!
      call fetcharg ('infile',set,arg)
      infile=trim(arg)
!
      return
      end
!#######################################################################
!#######################################################################
      subroutine compute_spline_1d (nx,x,f,s)
!
!-----------------------------------------------------------------------
!
! ****** Find the cubic spline coefficients for the 1D function
! ****** defined by array F(NX) with scale X(NX).
!
! ****** The spline coefficients are returned in structure S.
!
!-----------------------------------------------------------------------
!
      use number_types
      use spline_def
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      integer :: nx
      real(r_typ), dimension(nx) :: x
      real(r_typ), dimension(nx) :: f
      type(spl1d) :: s
!
!-----------------------------------------------------------------------
!
! ****** Allocate storage for the spline coefficients.
!
      s%nx=nx
!
      allocate (s%x(nx))
      allocate (s%f(nx))
      allocate (s%fxx(nx))
!
! ****** Evaluate the spline coefficients.
!
      s%x=x
      s%f=f
!
      call ezspline (nx,x,s%f,s%fxx)
!
      return
      end
!#######################################################################
      subroutine compute_spline_2d (nx,ny,x,y,f,s)
!
!-----------------------------------------------------------------------
!
! ****** Find the cubic spline coefficients for the 2D function
! ****** defined by array F(NX,NY), with scales X(NX) and Y(NY).
!
! ****** The spline coefficients are returned in structure S.
!
!-----------------------------------------------------------------------
!
      use number_types
      use spline_def
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
      real(r_typ), dimension(nx,ny) :: f
      type(spl2d) :: s
!
!-----------------------------------------------------------------------
!
      integer :: i,j
      real(r_typ), dimension(nx) :: gx,gppx
      real(r_typ), dimension(ny) :: gy,gppy
!
!-----------------------------------------------------------------------
!
! ****** Allocate storage for the spline coefficients.
!
      s%nx=nx
      s%ny=ny
!
      allocate (s%x(nx))
      allocate (s%y(ny))
      allocate (s%f(nx,ny))
      allocate (s%fxx(nx,ny))
      allocate (s%fyy(nx,ny))
      allocate (s%fxxyy(nx,ny))
!
! ****** Evaluate the spline coefficients.
!
      s%x=x
      s%y=y
      s%f=f
!
      do j=1,ny
        gx(:)=s%f(:,j)
        call ezspline (nx,x,gx,gppx)
        s%fxx(:,j)=gppx(:)
      enddo
!
      do i=1,nx
        gy(:)=s%f(i,:)
        call ezspline (ny,y,gy,gppy)
        s%fyy(i,:)=gppy(:)
      enddo
!
      do i=1,nx
        gy(:)=s%fxx(i,:)
        call ezspline (ny,y,gy,gppy)
        s%fxxyy(i,:)=gppy(:)
      enddo
!
      return
      end
!#######################################################################
      subroutine compute_spline_3d (nx,ny,nz,x,y,z,f,s)
!
!-----------------------------------------------------------------------
!
! ****** Find the cubic spline coefficients for the 3D function
! ****** defined by array F(NX,NY,NZ), with scales X(NX), Y(NY),
! ****** and Z(NZ).
!
! ****** The spline coefficients are returned in structure S.
!
!-----------------------------------------------------------------------
!
      use number_types
      use spline_def
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
      real(r_typ), dimension(nx,ny,nz) :: f
      type(spl3d) :: s
!
!-----------------------------------------------------------------------
!
      integer :: i,j,k
      real(r_typ), dimension(nx) :: gx,gppx
      real(r_typ), dimension(ny) :: gy,gppy
      real(r_typ), dimension(nz) :: gz,gppz
!
!-----------------------------------------------------------------------
!
! ****** Allocate storage for the spline coefficients.
!
      s%nx=nx
      s%ny=ny
      s%nz=nz
!
      allocate (s%x(nx))
      allocate (s%y(ny))
      allocate (s%z(nz))
      allocate (s%f(nx,ny,nz))
      allocate (s%fxx(nx,ny,nz))
      allocate (s%fyy(nx,ny,nz))
      allocate (s%fzz(nx,ny,nz))
      allocate (s%fxxyy(nx,ny,nz))
      allocate (s%fxxzz(nx,ny,nz))
      allocate (s%fyyzz(nx,ny,nz))
      allocate (s%fxxyyzz(nx,ny,nz))
!
! ****** Evaluate the spline coefficients.
!
      s%x=x
      s%y=y
      s%z=z
      s%f=f
!
!$omp parallel default(shared)
!$omp& private(i,j,k,gx,gppx,gy,gppy,gz,gppz)
!
!$omp do collapse(2) schedule(dynamic)
      do k=1,nz
        do j=1,ny
          gx(:)=s%f(:,j,k)
          call ezspline (nx,x,gx,gppx)
          s%fxx(:,j,k)=gppx(:)
        enddo
      enddo
!$omp end do
!
!$omp do collapse(2) schedule(dynamic)
      do k=1,nz
        do i=1,nx
          gy(:)=s%f(i,:,k)
          call ezspline (ny,y,gy,gppy)
          s%fyy(i,:,k)=gppy(:)
        enddo
      enddo
!$omp end do
!
!$omp do collapse(2) schedule(dynamic)
      do j=1,ny
        do i=1,nx
          gz(:)=s%f(i,j,:)
          call ezspline (nz,z,gz,gppz)
          s%fzz(i,j,:)=gppz(:)
        enddo
      enddo
!$omp end do
!$omp end parallel
!
!$omp parallel default(shared)
!$omp& private(i,j,k,gx,gppx,gy,gppy,gz,gppz)
!$omp do collapse(2) schedule(dynamic)
      do k=1,nz
        do i=1,nx
          gy(:)=s%fxx(i,:,k)
          call ezspline (ny,y,gy,gppy)
          s%fxxyy(i,:,k)=gppy(:)
        enddo
      enddo
!$omp end do
!$omp end parallel
!
!$omp parallel default(shared)
!$omp& private(i,j,k,gx,gppx,gy,gppy,gz,gppz)
!$omp do collapse(2) schedule(dynamic)
      do j=1,ny
        do i=1,nx
          gz(:)=s%fxx(i,j,:)
          call ezspline (nz,z,gz,gppz)
          s%fxxzz(i,j,:)=gppz(:)
          gz(:)=s%fyy(i,j,:)
          call ezspline (nz,z,gz,gppz)
          s%fyyzz(i,j,:)=gppz(:)
          gz(:)=s%fxxyy(i,j,:)
          call ezspline (nz,z,gz,gppz)
          s%fxxyyzz(i,j,:)=gppz(:)
        enddo
      enddo
!$omp end do
!$omp end parallel
!
      return
      end
!#######################################################################
      function evaluate_spline_1d (s,x,tab)
!
!-----------------------------------------------------------------------
!
! ****** Get the value of the 1D spline in structure S at the
! ****** point X.
!
! ****** The optional argument TAB is an inverse interpolation
! ****** table that can be used to speed up the search for the
! ****** interval that contains X.
!
! ****** The cubic spline coefficients for the spline S can
! ****** be obtained using routine COMPUTE_SPLINE_1D.
!
!-----------------------------------------------------------------------
!
      use number_types
      use spline_def
      use invint_def
      use locate_interval_interface
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      type(spl1d) :: s
      real(r_typ) :: x
      type(itab), optional :: tab
      real(r_typ) :: evaluate_spline_1d
!
!-----------------------------------------------------------------------
!
      integer :: i
!
!-----------------------------------------------------------------------
!
      real(r_typ), external :: speval
!
!-----------------------------------------------------------------------
!
! ****** Find the index of the cell enclosing point X.
!
      if (present(tab)) then
        i=locate_interval(s%nx,s%x,x,tab)
      else
        i=locate_interval(s%nx,s%x,x)
      end if
!
! ****** Interpolate in x.
!
      evaluate_spline_1d=speval(s%x(i  ),&
     &                          s%x(i+1),&
     &                          s%f(i  ),&
     &                          s%f(i+1),&
     &                          s%fxx(i  ),&
     &                          s%fxx(i+1),&
     &                          x)
!
      return
      end
!#######################################################################
      function evaluate_spline_2d (s,x,y,tabx,taby)
!
!-----------------------------------------------------------------------
!
! ****** Get the value of the 2D spline in structure S at the
! ****** point (X,Y).
!
! ****** The optional arguments TABX and TABY are inverse
! ****** interpolation tables that can be used to speed up the
! ****** search for the cell that contains (X,Y).
!
! ****** The cubic spline coefficients for the spline S can
! ****** be obtained using routine COMPUTE_SPLINE_2D.
!
!-----------------------------------------------------------------------
!
      use number_types
      use spline_def
      use invint_def
      use locate_interval_interface
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      type(spl2d) :: s
      real(r_typ) :: x,y
      type(itab), optional :: tabx,taby
      real(r_typ) :: evaluate_spline_2d
!
!-----------------------------------------------------------------------
!
      integer :: i,j,ii
      real(r_typ), dimension(0:1) :: f1,fxx1
      real(r_typ) :: a,b,c,d
!
!-----------------------------------------------------------------------
!
      real(r_typ), external :: speval
      real(r_typ), external :: speval_abcd
!
!-----------------------------------------------------------------------
!
! ****** Flag to use the more efficient version of this routine.
!
      logical, parameter :: use_faster=.true.
!
!-----------------------------------------------------------------------
!
! ****** Find the indices of the cell enclosing (X,Y).
!
      if (present(tabx)) then
        i=locate_interval(s%nx,s%x,x,tabx)
      else
        i=locate_interval(s%nx,s%x,x)
      end if
!
      if (present(taby)) then
        j=locate_interval(s%ny,s%y,y,taby)
      else
        j=locate_interval(s%ny,s%y,y)
      end if
!
      if (use_faster) go to 100
!
! ****** Leff efficient version.  This version
! ****** is slower but slightly easier to understand.
!
! ****** Interpolate in y.
!
      do ii=0,1
        f1(ii)=speval(s%y(j  ),&
     &                s%y(j+1),&
     &                s%f(i+ii,j  ),&
     &                s%f(i+ii,j+1),&
     &                s%fyy(i+ii,j  ),&
     &                s%fyy(i+ii,j+1),&
     &                y)
        fxx1(ii)=speval(s%y(j  ),&
     &                  s%y(j+1),&
     &                  s%fxx(i+ii,j  ),&
     &                  s%fxx(i+ii,j+1),&
     &                  s%fxxyy(i+ii,j  ),&
     &                  s%fxxyy(i+ii,j+1),&
     &                  y)
      enddo
!
! ****** Interpolate in x.
!
      evaluate_spline_2d=speval(s%x(i  ),&
     &                          s%x(i+1),&
     &                          f1(0),&
     &                          f1(1),&
     &                          fxx1(0),&
     &                          fxx1(1),&
     &                          x)
!
      return
!
  100 continue
!
! ****** More efficient version.  This version
! ****** is slighty faster than the one above.
!
! ****** Interpolate in y.
!
      call speval_get_abcd (s%y(j),s%y(j+1),y,a,b,c,d)
!
      do ii=0,1
        f1(ii)=speval_abcd(a,b,c,d,&
     &                     s%f(i+ii,j  ),&
     &                     s%f(i+ii,j+1),&
     &                     s%fyy(i+ii,j  ),&
     &                     s%fyy(i+ii,j+1))
        fxx1(ii)=speval_abcd(a,b,c,d,&
     &                       s%fxx(i+ii,j  ),&
     &                       s%fxx(i+ii,j+1),&
     &                       s%fxxyy(i+ii,j  ),&
     &                       s%fxxyy(i+ii,j+1))
      enddo
!
! ****** Interpolate in x.
!
      evaluate_spline_2d=speval(s%x(i  ),&
     &                          s%x(i+1),&
     &                          f1(0),&
     &                          f1(1),&
     &                          fxx1(0),&
     &                          fxx1(1),&
     &                          x)
!
      return
      end
!#######################################################################
      function evaluate_spline_3d (s,x,y,z,tabx,taby,tabz)
!
!-----------------------------------------------------------------------
!
! ****** Get the value of the 3D spline in structure S at the
! ****** point (X,Y,Z).
!
! ****** The optional arguments TABX, TABY, and TABZ are inverse
! ****** interpolation tables that can be used to speed up the
! ****** search for the cell that contains (X,Y,Z).
!
! ****** The cubic spline coefficients for the spline S can
! ****** be obtained using routine COMPUTE_SPLINE_3D.
!
!-----------------------------------------------------------------------
!
      use number_types
      use spline_def
      use invint_def
      use locate_interval_interface
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      type(spl3d) :: s
      real(r_typ) :: x,y,z
      type(itab), optional :: tabx,taby,tabz
      real(r_typ) :: evaluate_spline_3d
!
!-----------------------------------------------------------------------
!
      integer :: i,j,k,ii,jj
      real(r_typ), dimension(0:1) :: f1,fxx1
      real(r_typ), dimension(0:1,0:1) :: f2,fxx2,fyy2,fxxyy2
      real(r_typ) :: a,b,c,d
!
!-----------------------------------------------------------------------
!
      real(r_typ), external :: speval
      real(r_typ), external :: speval_abcd
!
!-----------------------------------------------------------------------
!
! ****** Flag to use the more efficient version of this routine.
!
      logical, parameter :: use_faster=.true.
!
!-----------------------------------------------------------------------
!
! ****** Find the indices of the cell enclosing (X,Y,Z).
!
      if (present(tabx)) then
        i=locate_interval(s%nx,s%x,x,tabx)
      else
        i=locate_interval(s%nx,s%x,x)
      end if
!
      if (present(taby)) then
        j=locate_interval(s%ny,s%y,y,taby)
      else
        j=locate_interval(s%ny,s%y,y)
      end if
!
      if (present(tabz)) then
        k=locate_interval(s%nz,s%z,z,tabz)
      else
        k=locate_interval(s%nz,s%z,z)
      end if
!
      if (use_faster) go to 100
!
! ****** Leff efficient version.  This version
! ****** is slower but slightly easier to understand.
!
! ****** Interpolate in z.
!
      do jj=0,1
        do ii=0,1
          f2(ii,jj)=speval(s%z(k  ),&
     &                     s%z(k+1),&
     &                     s%f(i+ii,j+jj,k  ),&
     &                     s%f(i+ii,j+jj,k+1),&
     &                     s%fzz(i+ii,j+jj,k  ),&
     &                     s%fzz(i+ii,j+jj,k+1),&
     &                     z)
          fxx2(ii,jj)=speval(s%z(k  ),&
     &                       s%z(k+1),&
     &                       s%fxx(i+ii,j+jj,k  ),&
     &                       s%fxx(i+ii,j+jj,k+1),&
     &                       s%fxxzz(i+ii,j+jj,k  ),&
     &                       s%fxxzz(i+ii,j+jj,k+1),&
     &                       z)
          fyy2(ii,jj)=speval(s%z(k  ),&
     &                       s%z(k+1),&
     &                       s%fyy(i+ii,j+jj,k  ),&
     &                       s%fyy(i+ii,j+jj,k+1),&
     &                       s%fyyzz(i+ii,j+jj,k  ),&
     &                       s%fyyzz(i+ii,j+jj,k+1),&
     &                       z)
          fxxyy2(ii,jj)=speval(s%z(k  ),&
     &                         s%z(k+1),&
     &                         s%fxxyy(i+ii,j+jj,k  ),&
     &                         s%fxxyy(i+ii,j+jj,k+1),&
     &                         s%fxxyyzz(i+ii,j+jj,k  ),&
     &                         s%fxxyyzz(i+ii,j+jj,k+1),&
     &                         z)
        enddo
      enddo
!
! ****** Interpolate in y.
!
      do ii=0,1
        f1(ii)=speval(s%y(j  ),&
     &                s%y(j+1),&
     &                f2(ii,0),&
     &                f2(ii,1),&
     &                fyy2(ii,0),&
     &                fyy2(ii,1),&
     &                y)
        fxx1(ii)=speval(s%y(j  ),&
     &                  s%y(j+1),&
     &                  fxx2(ii,0),&
     &                  fxx2(ii,1),&
     &                  fxxyy2(ii,0),&
     &                  fxxyy2(ii,1),&
     &                  y)
      enddo
!
! ****** Interpolate in x.
!
      evaluate_spline_3d=speval(s%x(i  ),&
     &                          s%x(i+1),&
     &                          f1(0),&
     &                          f1(1),&
     &                          fxx1(0),&
     &                          fxx1(1),&
     &                          x)
!
      return
!
  100 continue
!
! ****** More efficient version.  This version
! ****** is about 25% faster than the one above.
!
! ****** Interpolate in z.
!
      call speval_get_abcd (s%z(k),s%z(k+1),z,a,b,c,d)
!
      do jj=0,1
        do ii=0,1
          f2(ii,jj)=speval_abcd(a,b,c,d,&
     &                          s%f(i+ii,j+jj,k  ),&
     &                          s%f(i+ii,j+jj,k+1),&
     &                          s%fzz(i+ii,j+jj,k  ),&
     &                          s%fzz(i+ii,j+jj,k+1))
          fxx2(ii,jj)=speval_abcd(a,b,c,d,&
     &                            s%fxx(i+ii,j+jj,k  ),&
     &                            s%fxx(i+ii,j+jj,k+1),&
     &                            s%fxxzz(i+ii,j+jj,k  ),&
     &                            s%fxxzz(i+ii,j+jj,k+1))
          fyy2(ii,jj)=speval_abcd(a,b,c,d,&
     &                            s%fyy(i+ii,j+jj,k  ),&
     &                            s%fyy(i+ii,j+jj,k+1),&
     &                            s%fyyzz(i+ii,j+jj,k  ),&
     &                            s%fyyzz(i+ii,j+jj,k+1))
          fxxyy2(ii,jj)=speval_abcd(a,b,c,d,&
     &                              s%fxxyy(i+ii,j+jj,k  ),&
     &                              s%fxxyy(i+ii,j+jj,k+1),&
     &                              s%fxxyyzz(i+ii,j+jj,k  ),&
     &                              s%fxxyyzz(i+ii,j+jj,k+1))
        enddo
      enddo
!
! ****** Interpolate in y.
!
      call speval_get_abcd (s%y(j),s%y(j+1),y,a,b,c,d)
!
      do ii=0,1
        f1(ii)=speval_abcd(a,b,c,d,&
     &                     f2(ii,0),&
     &                     f2(ii,1),&
     &                     fyy2(ii,0),&
     &                     fyy2(ii,1))
        fxx1(ii)=speval_abcd(a,b,c,d,&
     &                       fxx2(ii,0),&
     &                       fxx2(ii,1),&
     &                       fxxyy2(ii,0),&
     &                       fxxyy2(ii,1))
      enddo
!
! ****** Interpolate in x.
!
      evaluate_spline_3d=speval(s%x(i  ),&
     &                          s%x(i+1),&
     &                          f1(0),&
     &                          f1(1),&
     &                          fxx1(0),&
     &                          fxx1(1),&
     &                          x)
!
      return
      end
!#######################################################################
      subroutine spline (n,x,f,ibc0,c0,ibc1,c1,fpp)
!
!-----------------------------------------------------------------------
!
! ****** Calculate cubic spline coefficients.
!
!-----------------------------------------------------------------------
!
! ****** Get the coefficients of a cubic spline interpolant to
! ****** the function values F(i) defined at the points X(i),
! ****** i=1,...,N.
!
! ****** The computed coefficients (which are actually the second
! ****** derivatives of F at the mesh points) are returned in the
! ****** array FPP.
!
! ****** Use routine SPLINT to evaluate the spline at a
! ****** particular position.
!
!-----------------------------------------------------------------------
!
! ****** The boundary conditions at the two ends are specified
! ****** by IBC0 and IBC1, and the coefficient arrays C0 and C1.
!
! ****** These are defined at x=X(1) according to the value
! ****** of IBC0 and C0.  (The conditions at x=X(N) are
! ****** specified similarly by IBC1 and C1.)
!
!        IBC0 = 1:  Set the second derivative at one end of the
!                   cell to equal the second derivative at the
!                   other end of the cell. [f''(1)=f''(2)]
!
!        IBC0 = 2:  Set the second derivative to zero.
!                   This corresponds to a "natural spline".
!                   [f''(1)=0.]
!
!        IBC0 = 3:  Set the first derivative to zero.
!                   [f'(1)=0.]
!
!        IBC0 = 4:  Set the second derivative to C0(1).
!                   [f''(1)=C0(1)]
!
!        IBC0 = 5:  Set the first derivative to C0(1).
!                   [f'(1)=C0(1)]
!
!        IBC0 = 6:  Set a linear combination of the first and
!                   second derivatives, according to C0(1),
!                   C0(2), and C0(3).
!                   [C0(2)*f'(1)+C0(3)*f''(1)=C0(1)]
!
!        IBC0 = 7:  Set a linear combination of the second
!                   derivatives at the left and right ends
!                   of the first cell.
!                   [C0(1)*f''(1)+C0(2)*f''(2)=C0(3)]
!
!        IBC0 = 8:  Set a linear combination of the first
!                   derivatives at the left and right ends
!                   of the first cell.
!                   [C0(1)*f'(1)+C0(2)*f'(2)=C0(3)]
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
      real(r_typ), dimension(n) :: x,f,fpp
      integer :: ibc0,ibc1
      real(r_typ), dimension(3) :: c0,c1
!
      intent(in) :: n,x,f,ibc0,c0,ibc1,c1
      intent(out) :: fpp
!
!-----------------------------------------------------------------------
!
      real(r_typ), parameter :: one=1._r_typ
      real(r_typ), parameter :: third=one/3._r_typ
      real(r_typ), parameter :: sixth=one/6._r_typ
!
!-----------------------------------------------------------------------
!
      integer :: i,j,ierr
      real(r_typ) :: dx,dxm,dxp,dxh
      real(r_typ), dimension(n) :: a,b,c
!
!-----------------------------------------------------------------------
!
! ****** Check that there are at least 3 points.
!
      if (n.lt.3) then
        write (*,*)
        write (*,*) '### ERROR in SPLINE:'
        write (*,*) '### Invalid number of points specified.'
        write (*,*) '### At least 3 points must be used.'
        write (*,*) 'Number of points specified = ',n
        call exit (1)
      end if
!
! ****** Check that the mesh is monotonic.
!
      dxm=x(2)-x(1)
      do i=2,n-1
        dx=x(i+1)-x(i)
        if (dx*dxm.le.0.) then
          write (*,*)
          write (*,*) '### ERROR in SPLINE:'
          write (*,*) '### The mesh is not monotonic.'
          write (*,*) 'At mesh point index = ',i
          write (*,*) 'Mesh-point values:'
          do j=1,n
            write (*,*) j,x(j)
          enddo
          call exit (1)
        end if
        dxm=dx
      enddo
!
!-----------------------------------------------------------------------
! ****** Set the coefficients for the tridiagonal solve at the
! ****** internal points.
!-----------------------------------------------------------------------
!
      do i=2,n-1
        dxp=x(i+1)-x(i)
        dxm=x(i)-x(i-1)
        dxh=dxp+dxm
        a(i)=dxh*third
        c(i)=dxm*sixth
        b(i)=dxp*sixth
        fpp(i)=(f(i+1)-f(i))/dxp-(f(i)-f(i-1))/dxm
      enddo
!
!-----------------------------------------------------------------------
! ****** Set the boundary condition at X(1).
!-----------------------------------------------------------------------
!
! ****** Set the unused value to zero.
!
      c(1)=0.
!
      select case (ibc0)
      case (1)
!
! ****** Second derivatives equal at the left and right ends
! ****** of the first cell.
!
        a(1)=one
        b(1)=-one
        fpp(1)=0.
!
      case (2)
!
! ****** Second derivative is zero ("natural splines").
!
        a(1)=one
        b(1)=0.
        fpp(1)=0.
!
      case (3)
!
! ****** First derivative is zero.
!
        dx=x(2)-x(1)
        a(1)=dx*third
        b(1)=dx*sixth
        fpp(1)=(f(2)-f(1))/dx
!
      case (4)
!
! ****** Second derivative is specified.
!
        a(1)=one
        b(1)=0.
        fpp(1)=c0(1)
!
      case (5)
!
! ****** First derivative is specified.
!
        dx=x(2)-x(1)
        a(1)=dx*third
        b(1)=dx*sixth
        fpp(1)=(f(2)-f(1))/dx-c0(1)
!
      case (6)
!
! ****** A combination of the first and second derivatives
! ****** is specified.
!
        if (c0(2).eq.0..and.c0(3).eq.0.) then
          write (*,*)
          write (*,*) '### ERROR in SPLINE:'
          write (*,*) '### Invalid boundary condition specified'//&
     &                ' at X(1).'
          write (*,*) '### Boundary condition type IBC0 = 6.'
          write (*,*) '### It is illegal for both C0(2) and C0(3)'//&
     &                ' to be zero.'
          call exit (1)
        end if
!
        if (c0(2).ne.0.) then
          dx=x(2)-x(1)
          a(1)=dx*third-c0(3)/c0(2)
          b(1)=dx*sixth
          fpp(1)=(f(2)-f(1))/dx-c0(1)/c0(2)
        else
          a(1)=one
          b(1)=0.
          fpp(1)=c0(1)/c0(3)
        end if
!
      case (7)
!
! ****** A linear combination of the second derivatives at the left
! ****** and right ends of the first cell is specified.
!
        a(1)=c0(1)
        b(1)=c0(2)
        fpp(1)=c0(3)
!
      case (8)
!
! ****** A linear combination of the first derivatives at the left
! ****** and right ends of the first cell is specified.
!
        dx=x(2)-x(1)
        a(1)=dx*(c0(1)/3-c0(2)/6)
        b(1)=dx*(c0(1)/6-c0(2)/3)
        fpp(1)=(c0(1)+c0(2))*(f(2)-f(1))/dx-c0(3)
!
      case default
!
        write (*,*)
        write (*,*) '### ERROR in SPLINE:'
        write (*,*) '### Invalid boundary condition specified at X(1).'
        write (*,*) '### IBC0 is invalid.'
        write (*,*) 'IBC0 = ',ibc0
        call exit (1)
!
      end select
!
!-----------------------------------------------------------------------
! ****** Set the boundary condition at X(N).
!-----------------------------------------------------------------------
!
! ****** Set the unused value to zero.
!
      b(n)=0.
!
      select case (ibc1)
      case (1)
!
! ****** Second derivatives equal at the left and right ends
! ****** of the last cell.
!
        a(n)=one
        c(n)=-one
        fpp(n)=0.
!
      case (2)
!
! ****** Second derivative is zero ("natural splines").
!
        a(n)=one
        c(n)=0.
        fpp(n)=0.
!
      case (3)
!
! ****** First derivative is zero.
!
        dx=x(n)-x(n-1)
        a(n)=dx*third
        c(n)=dx*sixth
        fpp(n)=-(f(n)-f(n-1))/dx
!
      case (4)
!
! ****** Second derivative is specified.
!
        a(n)=one
        c(n)=0.
        fpp(n)=c1(1)
!
      case (5)
!
! ****** First derivative is specified.
!
        dx=x(n)-x(n-1)
        a(n)=dx*third
        c(n)=dx*sixth
        fpp(n)=c1(1)-(f(n)-f(n-1))/dx
!
      case (6)
!
! ****** A combination of the first and second derivatives
! ****** is specified.
!
        if (c1(2).eq.0..and.c1(3).eq.0.) then
          write (*,*)
          write (*,*) '### ERROR in SPLINE:'
          write (*,*) '### Invalid boundary condition specified'//&
     &                ' at X(N).'
          write (*,*) '### Boundary condition type IBC1 = 6.'
          write (*,*) '### It is illegal for both C1(2) and C1(3)'//&
     &                ' to be zero.'
          call exit (1)
        end if
!
        if (c1(2).ne.0.) then
          dx=x(n)-x(n-1)
          a(n)=dx*third+c1(3)/c1(2)
          c(n)=dx*sixth
          fpp(n)=c1(1)/c1(2)-(f(n)-f(n-1))/dx
        else
          a(n)=one
          c(n)=0.
          fpp(n)=c1(1)/c1(3)
        end if
!
      case (7)
!
! ****** A linear combination of the second derivatives at the left
! ****** and right ends of the last cell is specified.
!
        a(n)=c1(1)
        c(n)=c1(2)
        fpp(n)=c1(3)
!
      case (8)
!
! ****** A linear combination of the first derivatives at the left
! ****** and right ends of the last cell is specified.
!
        dx=x(n)-x(n-1)
        a(n)=dx*(c1(1)/3-c1(2)/6)
        c(n)=dx*(c1(1)/6-c1(2)/3)
        fpp(n)=c1(3)-(c1(1)+c1(2))*(f(n)-f(n-1))/dx
!
      case default
!
        write (*,*)
        write (*,*) '### ERROR in SPLINE:'
        write (*,*) '### Invalid boundary condition specified at X(N).'
        write (*,*) '### IBC1 is invalid.'
        write (*,*) 'IBC1 = ',ibc1
        call exit (1)
!
      end select
!
!-----------------------------------------------------------------------
! ****** Solve the tridiagonal system for the second derivative.
!-----------------------------------------------------------------------
!
! ****** The second derivative is returned in FPP.
!
      call trid (n,c,a,b,fpp,ierr)
!
      if (ierr.ne.0) then
        write (*,*)
        write (*,*) '### ERROR in SPLINE:'
        write (*,*) '### The tridiagonal matrix relating the'//&
     &              ' spline coefficients was singular.'
        write (*,*) '### The spline could not be computed.'
        write (*,*) 'Number of mesh points = ',n
        write (*,*) 'Mesh-point values:'
        do j=1,n
          write (*,*) j,x(j)
        enddo
        call exit (1)
      end if
!
      return
      end
!#######################################################################
      subroutine spline_periodic_type1 (n,x,f,fpp)
!
!-----------------------------------------------------------------------
!
! ****** Calculate cubic spline coefficients for a periodic function.
!
!-----------------------------------------------------------------------
!
! ****** Get the coefficients of a cubic spline interpolant to
! ****** the function values F(i) defined at the points X(i),
! ****** i=1,...,N.
!
! ****** The computed coefficients (which are actually the second
! ****** derivatives of F at the mesh points) are returned in the
! ****** array FPP.
!
! ****** Use routine SPLINT to evaluate the spline at a
! ****** particular position.
!
!-----------------------------------------------------------------------
!
! ****** This routine assumes that the data in F is periodic, such
! ****** that F(N) = F(1), so that the first and last point are
! ****** repeated and represent the same location.  Thus the
! ****** periodicity length is X(N)-X(1).
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
      real(r_typ), dimension(n) :: x,f,fpp
!
      intent(in) :: n,x,f
      intent(out) :: fpp
!
!-----------------------------------------------------------------------
!
      real(r_typ), parameter :: one=1._r_typ
      real(r_typ), parameter :: third=one/3._r_typ
      real(r_typ), parameter :: sixth=one/6._r_typ
!
!-----------------------------------------------------------------------
!
      integer :: i,j,ierr,m
      real(r_typ) :: dx,dxm,dxp,dxh
      real(r_typ), dimension(n-1) :: a,b,c
!
!-----------------------------------------------------------------------
!
! ****** Check that there are at least 3 points.
!
      if (n.lt.3) then
        write (*,*)
        write (*,*) '### ERROR in SPLINE_PERIODIC_TYPE1:'
        write (*,*) '### Invalid number of points specified.'
        write (*,*) '### At least 3 points must be used.'
        write (*,*) 'Number of points specified = ',n
        call exit (1)
      end if
!
! ****** Check that the mesh is monotonic.
!
      dxm=x(2)-x(1)
      do i=2,n-1
        dx=x(i+1)-x(i)
        if (dx*dxm.le.0.) then
          write (*,*)
          write (*,*) '### ERROR in SPLINE_PERIODIC_TYPE1:'
          write (*,*) '### The mesh is not monotonic.'
          write (*,*) 'At mesh point index = ',i
          write (*,*) 'Mesh-point values:'
          do j=1,n
            write (*,*) j,x(j)
          enddo
          call exit (1)
        end if
        dxm=dx
      enddo
!
!-----------------------------------------------------------------------
! ****** Set the coefficients for the tridiagonal solve at the
! ****** internal points.
!-----------------------------------------------------------------------
!
      do i=2,n-1
        dxp=x(i+1)-x(i)
        dxm=x(i)-x(i-1)
        dxh=dxp+dxm
        a(i)=dxh*third
        c(i)=dxm*sixth
        b(i)=dxp*sixth
        fpp(i)=(f(i+1)-f(i))/dxp-(f(i)-f(i-1))/dxm
      enddo
!
!-----------------------------------------------------------------------
! ****** Set the periodic boundary condition at X(1).
!-----------------------------------------------------------------------
!
      dxp=x(2)-x(1)
      dxm=x(n)-x(n-1)
      dxh=dxp+dxm
      a(1)=dxh*third
      c(1)=dxm*sixth
      b(1)=dxp*sixth
      fpp(1)=(f(2)-f(1))/dxp-(f(n)-f(n-1))/dxm
!
!-----------------------------------------------------------------------
! ****** Solve the (periodic) tridiagonal system for the second
! ****** derivative.
!-----------------------------------------------------------------------
!
! ****** The second derivative is returned in FPP.
!
      m=n-1

      call trid_periodic (m,c,a,b,fpp,ierr)
!
      if (ierr.ne.0) then
        write (*,*)
        write (*,*) '### ERROR in SPLINE_PERIODIC_TYPE1:'
        write (*,*) '### The tridiagonal matrix relating the'//&
     &              ' spline coefficients was singular.'
        write (*,*) '### The spline could not be computed.'
        write (*,*) 'Number of mesh points = ',n
        write (*,*) 'Mesh-point values:'
        do j=1,n
          write (*,*) j,x(j)
        enddo
        call exit (1)
      end if
!
      fpp(n)=fpp(1)
!
      return
      end
!#######################################################################
      subroutine spline_periodic_type2 (n,x,f,fpp)
!
!-----------------------------------------------------------------------
!
! ****** Calculate cubic spline coefficients for a periodic function.
!
!-----------------------------------------------------------------------
!
! ****** Get the coefficients of a cubic spline interpolant to
! ****** the function values F(i) defined at the points X(i),
! ****** i=1,...,N.
!
! ****** The computed coefficients (which are actually the second
! ****** derivatives of F at the mesh points) are returned in the
! ****** array FPP.
!
! ****** Use routine SPLINT to evaluate the spline at a
! ****** particular position.
!
!-----------------------------------------------------------------------
!
! ****** This routine assumes that the data in F is periodic, such
! ****** that F(N-1) = F(1) and F(N)  = F(2), so that two sets of
! ****** points are repeated and represent the same locations.
! ****** Thus the periodicity length is X(N-1)-X(1).
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
      real(r_typ), dimension(n) :: x,f,fpp
!
      intent(in) :: n,x,f
      intent(out) :: fpp
!
!-----------------------------------------------------------------------
!
      real(r_typ), parameter :: one=1._r_typ
      real(r_typ), parameter :: third=one/3._r_typ
      real(r_typ), parameter :: sixth=one/6._r_typ
!
!-----------------------------------------------------------------------
!
      integer :: i,j,ierr,m
      real(r_typ) :: dx,dxm,dxp,dxh
      real(r_typ), dimension(2:n-1) :: a,b,c
!
!-----------------------------------------------------------------------
!
! ****** Check that there are at least 4 points.
!
      if (n.lt.4) then
        write (*,*)
        write (*,*) '### ERROR in SPLINE_PERIODIC_TYPE2:'
        write (*,*) '### Invalid number of points specified.'
        write (*,*) '### At least 4 points must be used.'
        write (*,*) 'Number of points specified = ',n
        call exit (1)
      end if
!
! ****** Check that the mesh is monotonic.
!
      dxm=x(2)-x(1)
      do i=2,n-1
        dx=x(i+1)-x(i)
        if (dx*dxm.le.0.) then
          write (*,*)
          write (*,*) '### ERROR in SPLINE_PERIODIC_TYPE2:'
          write (*,*) '### The mesh is not monotonic.'
          write (*,*) 'At mesh point index = ',i
          write (*,*) 'Mesh-point values:'
          do j=1,n
            write (*,*) j,x(j)
          enddo
          call exit (1)
        end if
        dxm=dx
      enddo
!
!-----------------------------------------------------------------------
! ****** Set the coefficients for the tridiagonal solve at the
! ****** internal points.
!-----------------------------------------------------------------------
!
      do i=2,n-1
        dxp=x(i+1)-x(i)
        dxm=x(i)-x(i-1)
        dxh=dxp+dxm
        a(i)=dxh*third
        c(i)=dxm*sixth
        b(i)=dxp*sixth
        fpp(i)=(f(i+1)-f(i))/dxp-(f(i)-f(i-1))/dxm
      enddo
!
!-----------------------------------------------------------------------
! ****** Solve the (periodic) tridiagonal system for the second
! ****** derivative.
!-----------------------------------------------------------------------
!
! ****** The second derivative is returned in FPP.
!
      m=n-2

      call trid_periodic (m,c,a,b,fpp(2),ierr)
!
      if (ierr.ne.0) then
        write (*,*)
        write (*,*) '### ERROR in SPLINE_PERIODIC_TYPE2:'
        write (*,*) '### The tridiagonal matrix relating the'//&
     &              ' spline coefficients was singular.'
        write (*,*) '### The spline could not be computed.'
        write (*,*) 'Number of mesh points = ',n
        write (*,*) 'Mesh-point values:'
        do j=1,n
          write (*,*) j,x(j)
        enddo
        call exit (1)
      end if
!
      fpp(n)=fpp(2)
      fpp(1)=fpp(n-1)
!
      return
      end
!#######################################################################
      subroutine ezspline (n,x,f,fpp)
!
!-----------------------------------------------------------------------
!
! ****** Calculate cubic spline coefficients.
!
! ****** Easy-to-use version of SPLINE.
!
!-----------------------------------------------------------------------
!
! ****** This routine calls SPLINE with boundary conditions
! ****** IBC0=1 and IBC1=1.  See the comments in SPLINE to
! ****** see what this means.
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
      real(r_typ), dimension(n) :: x,f,fpp
!
      intent(in) :: n,x,f
      intent(out) :: fpp
!
!-----------------------------------------------------------------------
!
      integer :: ibc0,ibc1
      real(r_typ), dimension(3) :: c0,c1
!
!-----------------------------------------------------------------------
!
      ibc0=1
      ibc1=1
!
      c0(:)=0.
      c1(:)=0.
!
      call spline (n,x,f,ibc0,c0,ibc1,c1,fpp)
!
      return
      end
!#######################################################################
      function splint (n,x,f,fpp,xv,tab)
!
!-----------------------------------------------------------------------
!
! ****** Evaluate a 1D cubic spline.
!
!-----------------------------------------------------------------------
!
! ****** Evaluate a cubic spline interpolant at X=XV.
!
! ****** On input, the function values F(i) and the second
! ****** derivatives FPP(i), defined at the points X(i),
! ****** i=1,...,N, need to be specified.
!
! ****** The optional argument TAB is an inverse interpolation
! ****** table that can be used to speed up the search for the
! ****** interval that contains XV.
!
! ****** The routine SPLINE can be used to compute FPP.
!
! ****** This routine uses routine SPEVAL to evaluate the spline.
!
! ****** The value of the spline at XV is returned as the
! ****** function value.
!
!-----------------------------------------------------------------------
!
      use number_types
      use invint_def
      use locate_interval_interface
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      integer :: n
      real(r_typ), dimension(n) :: x,f,fpp
      real(r_typ) :: xv
      type(itab), optional :: tab
      real(r_typ) :: splint
!
      intent(in) :: n,x,f,fpp,xv,tab
!
!-----------------------------------------------------------------------
!
      integer :: i
!
!-----------------------------------------------------------------------
!
      real(r_typ), external :: speval
!
!-----------------------------------------------------------------------
!
! ****** Find the mesh interval that encloses XV.
!
      if (present(tab)) then
        i=locate_interval(n,x,xv,tab)
      else
        i=locate_interval(n,x,xv)
      end if
!
! ****** Evaluate the cubic spline.
!
      splint=speval(x(i),x(i+1),f(i),f(i+1),fpp(i),fpp(i+1),xv)
!
      return
      end
!#######################################################################
      function speval (x1,x2,f1,f2,fpp1,fpp2,xv)
!
!-----------------------------------------------------------------------
!
! ****** Evaluate a cubic spline.
!
!-----------------------------------------------------------------------
!
! ****** Evaluate a cubic spline interpolant at X=XV.
!
! ****** On input, the function values F1 and F2 and the second
! ****** derivatives FPP1 and FPP2, defined at the left and right
! ****** ends of the interval X1 and X2 need to be specified.
!
! ****** The value of the spline at XV is returned as the
! ****** function value.
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
      real(r_typ) :: x1,x2,f1,f2,fpp1,fpp2
      real(r_typ) :: xv
      real(r_typ) :: speval
!
      intent(in) :: x1,x2,f1,f2,fpp1,fpp2,xv
!
!-----------------------------------------------------------------------
!
      real(r_typ), parameter :: one=1._r_typ
      real(r_typ), parameter :: sixth=one/6._r_typ
!
!-----------------------------------------------------------------------
!
      real(r_typ) :: dx,a,b,c,d
!
!-----------------------------------------------------------------------
!
! ****** Evaluate the cubic spline.
!
      dx=x2-x1
!
      b=(xv-x1)/dx
      a=one-b
!
      c=a*(a**2-one)*dx**2*sixth
      d=b*(b**2-one)*dx**2*sixth
!
      speval=a*f1+b*f2+c*fpp1+d*fpp2
!
      return
      end
!#######################################################################
      subroutine speval_get_abcd (x1,x2,xv,a,b,c,d)
!
!-----------------------------------------------------------------------
!
! ****** Evaluate a cubic spline.
!
! ****** This version splits the calculation into two parts for
! ****** efficiency.
!
! ****** First call SPEVAL_GET_ABCD, and then call SPEVAL_ABCD.
!
! ****** Typically SPEVAL_GET_ABCD and SPEVAL_GET_ABCD are used
! ****** when multiple spline evaluations are being performed for
! ****** the same position.
!
!-----------------------------------------------------------------------
!
! ****** Evaluate a cubic spline interpolant at X=XV.
!
! ****** On input, the left and right limits of the spline interval
! ****** X1 and X2, as well as the position at which the spline is
! ****** being evaluated, XV, need to be specified.
!
! ****** The coefficients of the spline A, B, C, and D are
! ****** returned on output.
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
      real(r_typ) :: x1,x2,xv,a,b,c,d
!
      intent(in) :: x1,x2,xv
      intent(out) :: a,b,c,d
!
!-----------------------------------------------------------------------
!
      real(r_typ), parameter :: one=1._r_typ
      real(r_typ), parameter :: sixth=one/6._r_typ
!
!-----------------------------------------------------------------------
!
      real(r_typ) :: dx
!
!-----------------------------------------------------------------------
!
! ****** Evaluate the coefficients of the cubic spline.
!
      dx=x2-x1
!
      b=(xv-x1)/dx
      a=one-b
!
      c=a*(a**2-one)*dx**2*sixth
      d=b*(b**2-one)*dx**2*sixth
!
      return
      end
!#######################################################################
      function speval_abcd (a,b,c,d,f1,f2,fpp1,fpp2)
!
!-----------------------------------------------------------------------
!
! ****** Evaluate a cubic spline.
!
! ****** This version splits the calculation into two parts for
! ****** efficiency.
!
! ****** First call SPEVAL_GET_ABCD, and then call SPEVAL_ABCD.
!
! ****** Typically SPEVAL_GET_ABCD and SPEVAL_GET_ABCD are used
! ****** when multiple spline evaluations are being performed for
! ****** the same position.
!
!-----------------------------------------------------------------------
!
! ****** On input, the coefficients A, B, C, and D, and the
! ****** function values F1 and F2 and the second
! ****** derivatives FPP1 and FPP2, need to be specified.
!
! ****** The value of the spline is returned as the
! ****** function value.
!
! ****** The coefficients A, B, C, and D can be obtained using
! ****** routine SPEVAL_GET_ABCD.
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
      real(r_typ) :: a,b,c,d,f1,f2,fpp1,fpp2
      real(r_typ) :: speval_abcd
!
      intent(in) :: a,b,c,d,f1,f2,fpp1,fpp2
!
!-----------------------------------------------------------------------
!
! ****** Evaluate the cubic spline.
!
      speval_abcd=a*f1+b*f2+c*fpp1+d*fpp2
!
      return
      end
!#######################################################################
      function locate_interval (n,x,xv,tab,ierr)
!
!-----------------------------------------------------------------------
!
! ****** Locate a mesh interval.
!
!-----------------------------------------------------------------------
!
! ****** Find the interval I in table X(i), i=1,2,...,N,
! ****** that encloses the value XV, i.e., such that
! ****** X(I).le.XV.le.X(I+1).
!
! ****** For the special case when N=1, XV must equal X(1)
! ****** exactly, otherwise an error occurs.
!
! ****** The optional argument TAB is an inverse interpolation
! ****** table that can be used to speed up the search for the
! ****** interval.
!
! ****** If the optional argument IERR is specified, then this
! ****** routine will return when an error occurs with IERR=1.
! ****** If no error occurs, IERR=0 is returned.  When IERR is not
! ****** specified, this routine will terminate the program
! ****** with a printed error message.
!
! ****** The mesh interval I is returned as the function value.
!
!-----------------------------------------------------------------------
!
      use number_types
      use invint_def
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
      type(itab), optional :: tab
      integer, optional :: ierr
      integer :: locate_interval
!
      intent(in) :: n,x,xv,tab
!
!-----------------------------------------------------------------------
!
      real(r_typ), parameter :: one=1._r_typ
!
!-----------------------------------------------------------------------
!
      integer :: i,ig
      real(r_typ) :: xi,fiv,alpha
!
!-----------------------------------------------------------------------
!
      if (present(ierr)) then
        ierr=0
      end if
!
! ****** For the special case when the table has only one
! ****** point (N=1), the inverse table is not used.  In this
! ****** case it is necessary for XV to equal X(1) exactly,
! ****** otherwise this routine exits with an error.
!
      if (n.eq.1) then
        if (xv.eq.x(1)) then
          locate_interval=i
          return
        else
          go to 900
        end if
      end if
!
! ****** Search for the interval depending on whether the optional
! ****** inverse interpolation table TAB was specified.
!
      if (.not.present(tab)) then
!
! ****** Search without an inverse interpolation table.
!
        do i=1,n-1
          if (xv.ge.x(i).and.xv.le.x(i+1)) then
            locate_interval=i
            return
          end if
        enddo
!
      else
!
! ****** Search with an inverse interpolation table.
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
              locate_interval=i
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
              locate_interval=i
              return
            end if
          enddo
!
        end if
!
      end if
!
  900 continue
!
! ****** Value not found.
!
! ****** If IERR was passed, set IERR=1 and return; otherwise,
! ****** write an error message and terminate the program.
!
      if (present(ierr)) then
        ierr=1
        return
      else
        write (*,*)
        write (*,*) '### ERROR in LOCATE_INTERVAL:'
        write (*,*) '### The value requested was not found in'//&
     &              ' the table:'
        write (*,*) 'Value requested = ',xv
        write (*,*) 'Minimum table value = ',x(1)
        write (*,*) 'Maximum table value = ',x(n)
        call exit (1)
      end if
!
      end
!#######################################################################
      subroutine trid (n,c,a,b,d,ierr)
!
!-----------------------------------------------------------------------
!
! ****** Solve the tridiagonal system of equations:
!
!         C(i)*X(i-1) + A(i)*X(i) + B(i)*X(i+1) = D(i)
!
!        for i=2,...,N-1, with
!
!           A(1)*X(1) + B(1)*X(2) = D(1)
!
!        and
!
!           C(N)*X(N-1) + A(N)*X(N) = D(N)
!
! ****** Note that C(1) and B(N) are not referenced.
!
! ****** D is overwritten with the solution.
!
! ****** Return IERR=0 for a successful completion.  If the
! ****** matrix is singular, this routine returns IERR=1 and
! ****** D is invalid.
!
! ****** This routine does not do any pivoting, so the solution
! ****** is not guaranteed to be accurate unless the matrix is
! ****** well-conditioned (e.g., diagonally dominant).
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
      real(r_typ), dimension(n) :: c,a,b,d
      integer :: ierr
!
      intent(in) :: n,c,a,b
      intent(inout) :: d
      intent(out) :: ierr
!
!-----------------------------------------------------------------------
!
      real(r_typ), parameter :: one=1._r_typ
!
!-----------------------------------------------------------------------
!
      integer :: i
      real(r_typ) :: denom,ace
      real(r_typ), dimension(n) :: aa
!
!-----------------------------------------------------------------------
!
      ierr=1
!
! ****** Copy A to AA, since it will be overwritten during the
! ****** elimination.  This prevents A from being overwritten.
!
      aa=a
!
! ****** Forward elimination.
!
      if (aa(1).eq.0.) return
      d(1)=d(1)/aa(1)
      aa(1)=b(1)/aa(1)
      do i=2,n
        denom=aa(i)-c(i)*aa(i-1)
        if (denom.eq.0.) return
        ace=one/denom
        if (i.ne.n) aa(i)=ace*b(i)
        d(i)=ace*(d(i)-c(i)*d(i-1))
      enddo
!
! ****** Backward substitution.
!
      do i=n-1,1,-1
        d(i)=d(i)-aa(i)*d(i+1)
      enddo
!
! ****** Set the error return flag to indicate successful completion.
!
      ierr=0
!
      return
      end
!#######################################################################
      subroutine trid_periodic (n,c,a,b,d,ierr)
!
!-----------------------------------------------------------------------
!
! ****** Solve the tridiagonal system of equations:
!
!         C(i)*X(i-1) + A(i)*X(i) + B(i)*X(i+1) = D(i)
!
!        for i=2,...,N-1, with
!
!           C(1)*X(N) + A(1)*X(1) + B(1)*X(2) = D(1)
!
!        and
!
!           C(N)*X(N-1) + A(N)*X(N) + B(N)*X(1) = D(N)
!
! ****** D is overwritten with the solution.
!
! ****** Return IERR=0 for a successful completion.  If the
! ****** matrix is singular, this routine returns IERR=1 and
! ****** D is invalid.
!
! ****** This routine does not do any pivoting, so the solution
! ****** is not guaranteed to be accurate unless the matrix is
! ****** well-conditioned (e.g., diagonally dominant).
!
! ****** This system arises for periodic solutions.
!
! ****** This routine uses the Sherman-Morrison formula for
! ****** updating a matrix inverse with a low-rank modification.
! ****** The modification arises from the changes introduced by the
! ****** periodic boundary conditions.
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
      real(r_typ), dimension(n) :: c,a,b,d
      integer :: ierr
!
      intent(in) :: n,c,a,b
      intent(inout) :: d
      intent(out) :: ierr
!
!-----------------------------------------------------------------------
!
      real(r_typ), parameter :: one=1._r_typ
!
!-----------------------------------------------------------------------
!
      integer :: i
      real(r_typ) :: denom,ace
      real(r_typ), dimension(n) :: aa
      real(r_typ), dimension(n) :: y
      real(r_typ), dimension(n,2) :: z2
      real(r_typ), dimension(2,2) :: t,tinv
      real(r_typ) :: detinv
      real(r_typ), dimension(2) :: tvty
!
!-----------------------------------------------------------------------
!
      ierr=1
!
! ****** First, solve the (non-periodic) system with an RHS
! ****** equal to D.
!
      y=d
!
      call trid (n,c,a,b,y,ierr)
!
      if (ierr.ne.0) return
!
! ****** Next, solve the two systems for the "inhomogenous part".
!
      z2=0.
      z2(1,1)=one
      z2(n,2)=one
!
! ****** Copy A to AA, since it will be overwritten during the
! ****** elimination.  This prevents A from being overwritten.
!
      aa=a
!
! ****** Forward elimination.
!
      if (aa(1).eq.0.) return
      z2(1,:)=z2(1,:)/aa(1)
      aa(1)=b(1)/aa(1)
      do i=2,n
        denom=aa(i)-c(i)*aa(i-1)
        if (denom.eq.0.) return
        ace=one/denom
        if (i.ne.n) aa(i)=ace*b(i)
        z2(i,:)=ace*(z2(i,:)-c(i)*z2(i-1,:))
      enddo
!
! ****** Backward substitution.
!
      do i=n-1,1,-1
        z2(i,:)=z2(i,:)-aa(i)*z2(i+1,:)
      enddo
!
! ****** Invert the 2 x 2 system.
!
      t(1,1)=one+z2(n,1)*b(n)
      t(1,2)=z2(n,2)*b(n)
      t(2,1)=z2(1,1)*c(1)
      t(2,2)=one+z2(1,2)*c(1)
!
      denom=t(1,1)*t(2,2)-t(1,2)*t(2,1)
      if (denom.eq.0.) return
      detinv=one/denom
!
      tinv(1,1)= detinv*t(2,2)
      tinv(2,2)= detinv*t(1,1)
      tinv(1,2)=-detinv*t(1,2)
      tinv(2,1)=-detinv*t(2,1)
!
! ****** Construct the final periodic solution.
!
      tvty(1)=tinv(1,1)*b(n)*y(n)+tinv(1,2)*c(1)*y(1)
      tvty(2)=tinv(2,1)*b(n)*y(n)+tinv(2,2)*c(1)*y(1)
!
      d(:)=y(:)-z2(:,1)*tvty(1)-z2(:,2)*tvty(2)
!
! ****** Set the error return flag to indicate successful completion.
!
      ierr=0
!
      return
      end
!#######################################################################
      subroutine get_dips_map_3d
!
!-----------------------------------------------------------------------
!
! ****** Compute a 3D dips map.
!
!-----------------------------------------------------------------------
!
      use number_types
      use types
      use mesh
      use field
      use vars
      use files
      use params
      use diags
      use openmp_vars
      use tracefl_interface
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
! ****** Storage for the coronal hole map.
!
      real(r_typ), dimension(nrss,ntss,npss) :: dips
!
!-----------------------------------------------------------------------
!
      integer :: ierr,i,j,k
      real(r_typ), dimension(3) :: xfl0
      integer :: n_completed,n_total,nc,diag_step
      real(r_typ) :: pct_done, private_dips
!
!-----------------------------------------------------------------------
!
      if (verbose) then
        write (*,*)
        write (*,*) '### Computing a 3D dips map:'
      end if
!
! ****** Check that the coronal hole map output file name is not
! ****** blank, since this does not make sense.
!
      if (dips_map_3d_output_file.eq.' ') then
        write (*,*)
        write (*,*) '### ERROR in GET_DIPS_MAP_3D:'
        write (*,*) '### A coronal hole map was requested, yet'//&
     &              ' the output file name is blank.'
        call exit (1)
      end if


      if (verbose) then
        write (*,*)
        write (*,*) 'Mapping field lines ...'
        write (*,*)
      end if
!
      n_total=nrss*ntss*npss
      n_completed=0

!
!$omp parallel do
!$omp& private(i,j,k,xfl0,private_dips)
!$omp& private(nc,diag_step,pct_done)
!$omp& collapse(3)
!$omp& schedule(dynamic,iterations_per_thread)
      do i=1,nrss
        do k=1,npss
          do j=1,ntss
!
! ****** Update the iteration counter for diagnostic
! ****** purposes.
!
            if (verbose) then
!$omp critical
              n_completed=n_completed+1
              nc=n_completed
!$omp end critical
            end if
!
            xfl0(1)=rss(i)
            xfl0(2)=tss(j)
            xfl0(3)=pss(k)
            private_dips=0.
            call get_dip(xfl0,private_dips)
            dips(i,j,k)=private_dips
!
! ****** Write progress diagnostics if requested.
!
            if (verbose) then
              diag_step=mod(nc,diagnostic_interval)
              if (diag_step.eq.0) then
                pct_done=100.*nc/n_total
                write (*,910) 'Fraction completed: ',pct_done
  910           format (1x,a,f7.3,'%')
              end if
            end if
!
          enddo
        enddo
      enddo
!$omp end parallel do

!
! ****** Write the coronal hole map.
!
      if (verbose) then
        write (*,*)
        write (*,*) 'Writing the 3D coronal hole map to file: ',&
     &              trim(dips_map_3d_output_file)
      end if
!
      call wrhdf_3d (dips_map_3d_output_file,.true.,&
     &               nrss,ntss,npss,dips,rss,tss,pss,&
     &               hdf32,ierr)
!
      if (ierr.ne.0) then
        write (*,*)
        write (*,*) '### ERROR in GET_DIPS_MAP_3D:'
        write (*,*) '### Could not write the coronal hole map.'
        call exit (1)
      end if
!
      return
      end
!#######################################################################
      subroutine get_dip(xfl0,dips)
!
!-----------------------------------------------------------------------
!
! ****** Compute a 3D dips
!
!-----------------------------------------------------------------------
!
      use number_types
      use types
      use mesh
      use field
      use vars
      use files
      use params
      use diags
      use openmp_vars
      use tracefl_interface
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
!
      real(r_typ) :: dips
!
!-----------------------------------------------------------------------
!
      integer :: is
      real(r_typ), dimension(3) :: xfl0,xfl1
      real(r_typ), dimension(3) :: bs0,bs1
      logical :: ttb
      real(r_typ) :: s
      type(flparam) :: ds_f,ds_b
      logical :: f_trace_reached_boundary
      logical :: b_trace_reached_boundary
      logical :: f_trace_on_r0,f_trace_on_r1
      logical :: b_trace_on_r0,b_trace_on_r1
      logical :: f_br_positive
      logical :: b_br_positive
      logical :: dipb,dipf
!
! ****** Field line trace storage buffers.
!
      type(traj) :: xtf, xtb
!
!-----------------------------------------------------------------------
!
! ****** Set the tracing direction to be either along the direction
! ****** of the magnetic field or along the directon of increasing
! ****** radius.
!
      ds%direction_is_along_b=trace_slice_direction_is_along_b
!
      ds_f=ds
      ds_f%direction=1
!
      ds_b=ds
      ds_b%direction=-1
!
! ****** Trace a field line in the forward direction along B.
!
      call allocate_trajectory_buffer (xtf)
      call allocate_trajectory_buffer (xtb)
      call tracefl (b,ds_f,xfl0,xfl1,bs0,bs1,s,ttb,xtf)
!
! ****** Check that the field line reached R0 or R1.
!
      if (ttb) then
        f_trace_reached_boundary=.true.
        f_trace_on_r0=xfl1(1).eq.b%lim0(1)
        f_trace_on_r1=xfl1(1).eq.b%lim1(1)
        f_br_positive=bs1(1).ge.0.
      else
        f_trace_reached_boundary=.false.
      end if
!
! ****** Trace a field line in the backward direction along B.
!
      call tracefl (b,ds_b,xfl0,xfl1,bs0,bs1,s,ttb,xtb)
!
! ****** Check that the field line reached R0 or R1.
!
      if (ttb) then
        b_trace_reached_boundary=.true.
        b_trace_on_r0=xfl1(1).eq.b%lim0(1)
        b_trace_on_r1=xfl1(1).eq.b%lim1(1)
        b_br_positive=bs1(1).ge.0.
      else
        b_trace_reached_boundary=.false.
      end if
!
      if (f_trace_reached_boundary.and.&
     &      b_trace_reached_boundary) then
        if (f_trace_on_r0.and.b_trace_on_r0) then
          dipf=.false.
          dipb=.false.
          do is=1,min(xtf%npts,ns_dips)
            if(xtf%x(1)%f(is).gt.xfl0(1)) then
              dipf=.true.
              exit
            endif
          enddo
          do is=1,min(xtb%npts,ns_dips)
            if(xtb%x(1)%f(is).gt.xfl0(1)) then
              dipb=.true.
              exit
            endif
          enddo
          if (dipb.and.dipf) then
            dips=one
          endif
        end if
      end if
      call deallocate_trajectory_buffer (xtf)
      call deallocate_trajectory_buffer (xtb)
!
      return
      end
!#######################################################################
      subroutine set_up_integration
!
!-----------------------------------------------------------------------
!
! ****** Set up integration along field lines.
! ****** Read scalar field.
!
!-----------------------------------------------------------------------
!
      use number_types
      use vars
      use params
      use integrate_fl
      use files
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      integer :: ierr
!
!-----------------------------------------------------------------------
!
      do_integral_along_fl=.true.
!
! ****** Read the scalar field
!
      if (verbose) then
        write (*,*)
        write (*,*) 'Reading data file: ',trim(scalar_input_file)
      end if
!
      call rdhdf (scalar_input_file,scalar_field,ierr)
!
      if (ierr.ne.0) then
        write (*,*)
        write (*,*) '### ERROR in SET_UP_INTEGRATION:'
        write (*,*) '### Could not read scalar field.'
        write (*,*) 'IERR (from RDHDF_3D) = ',ierr
        write (*,*) 'File name: ',trim(scalar_input_file)
        call exit (1)
      end if
!
      if (scalar_field%ndim.ne.3.or..not.scalar_field%scale)  then
        write (*,*)
        write (*,*) '### ERROR in SET_UP_INTEGRATION:'
        write (*,*) '### Invalid or missing scales in scalar file.'
        write (*,*) 'File name: ',trim(scalar_input_file)
        call exit (1)
      end if
!
      call build_inverse_tables (scalar_field,inv_sf)
!
      return
      end
