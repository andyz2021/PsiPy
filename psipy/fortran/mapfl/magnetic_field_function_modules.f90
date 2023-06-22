!#######################################################################
      module magfld_func_def
!
!-----------------------------------------------------------------------
! ****** Definition of the analytic function that defines B.
!-----------------------------------------------------------------------
!
      implicit none
!
! ****** Number of defined functions.
!
      integer, parameter :: number_of_functions=2
!
! ****** Mnemonics for defined functions.
! ****** There should be NUMBER_OF_FUNCTIONS of these.
!
      integer, parameter :: FUNC_TYPE_DIPOLE           =1
      integer, parameter :: FUNC_TYPE_PFSS_BKG         =2
!
      end module
!#######################################################################
      module magfld_func_index
!
!-----------------------------------------------------------------------
! ****** Index for the analytic function that defines B.
!-----------------------------------------------------------------------
!
      implicit none
!
! ****** Selected function index.
!
      integer :: function_index=0
!
      end module
!#######################################################################
      module magfld_func_params
!
!-----------------------------------------------------------------------
! ****** Parameters for the analytic function that defines B.
!-----------------------------------------------------------------------
!
      use number_types
!
      implicit none
!
! ****** Parameters for the DIPOLE function.
!
      real(r_typ) :: b0=1._r_typ
!
! ****** Parameters for the PFSS_BKG function.
!
      real(r_typ) :: mu
      real(r_typ) :: rss
!
      end module
!#######################################################################
      module magfld_func_namelist
!
!-----------------------------------------------------------------------
! ****** Namelist parameters for the analytic function that defines B.
!-----------------------------------------------------------------------
!
      use magfld_func_index
      use magfld_func_params
!
      implicit none
!
! ****** NAMELIST to read in the parameters.
!
      namelist /function_parameters/ function_index,&
     &                               b0,&
     &                               mu,rss
!
      end module
