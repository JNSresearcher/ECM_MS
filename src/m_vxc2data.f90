! module with declarations of global variables, arrays and structures
! Fortran code created by J.Sochor   ( https://github.com/JNSresearcher )
! e-mail: JNSresearcher@gmail.com

MODULE m_vxc2data
implicit none

! structures for calculating functions
TYPE tFun
    CHARACTER (LEN=1)               :: ex               ! source direction identifier: X, Y, Z
    REAL(8)                         :: vely             ! function value
    CHARACTER (LEN=50)              :: eqn              ! function expression
    INTEGER                         :: args,nomsch      ! number of arguments and domain number
    CHARACTER (LEN=8),ALLOCATABLE   :: namex(:)         ! argument names
    REAL(8),ALLOCATABLE             :: velx(:)          ! argument values
END TYPE tFun
    TYPE (tFun), ALLOCATABLE :: Fun(:), Vmech(:), Venv(:)    
! structures tFun  for: 
! Fun   - size=numfun,  functions sources of coils; 
! Vmech - size=numMech, functions for coils movement  
! Venv  - size=numVenv, functions for ferromagnetics movement

! structure for cells in which functions are located
TYPE tfun_nod
!  number of cells for vector sources of coils
    INTEGER             :: numnod_Fx,  numnod_Fy,  numnod_Fz 
! arrays with cells numbers for vector  functions
    INTEGER,ALLOCATABLE :: nods_Fx(:), nods_Fy(:), nods_Fz(:)
    REAL(8),ALLOCATABLE :: nods_Fx_cos(:), nods_Fy_cos(:), nods_Fz_cos(:)
    
! speed, distance and relative shift of the source movement along the X Y Z axes
    REAL(8)             :: val_Vmech(3), val_Vmech_old(3), Distance(3), shift(3)  
    INTEGER             :: num_Vmech(3),  &    ! X Y Z axis motion function number
                           move(3),       &    ! sign of motion in input data
                           ism(3),        &    ! sign of change motion 
                           length(3)           ! distance in integers: length = nint(Distance)
END TYPE tfun_nod
TYPE (tfun_nod), ALLOCATABLE :: fun_nod(:)  
! structures tfun_nod - for moving coils
! fun_nod - size=numfun

TYPE tfun_Venv
! speed, distance and relative shift of the source movement along the X Y Z axes
    REAL(8)             :: val_Venv(3), val_Venv_old(3), Distance(3),  shift(3)  
    INTEGER             :: num_Venv(3),   &    ! X Y Z axis motion function number
                           move(3),       &    ! sign of motion in input data
                           ism(3),        &    ! sign of change motion 
                           length(3)           ! distance in integers: length = nint(Distance)
END TYPE tfun_Venv
TYPE (tfun_Venv), ALLOCATABLE :: fun_Venv(:) 
! structures tfun_Venv - for driving ferromagnets
! fun_Venv - size=numVenv

TYPE t_non_lin_MU
    CHARACTER(len=:), ALLOCATABLE :: file_MU
    integer num_domain
    character (LEN=2):: tip_non_lin
END TYPE t_non_lin_MU
TYPE (t_non_lin_MU), ALLOCATABLE :: non_lin_MU(:)
INTEGER siz_non_lin_MU

!  for domain Fe - conduction region
INTEGER ::    sizdom_Fe

! geoPHYS - array with physical domain numbers, numdom_Fe - with ferrum domains numbers
INTEGER(1),ALLOCATABLE:: geoPHYS(:,:,:), numdom_Fe(:)     

! valPHYS - array of domain parameters
! valPHYS(:,1) -  parameters environs
! valPHYS(:,2) - not use
! valPHYS(:,3) - domain speed along the X axis
! valPHYS(:,4) - domain speed along the Y axis
! valPHYS(:,5) - domain speed along the Z axis
REAL(8), ALLOCATABLE  :: valPHYS(:,:) 
CHARACTER (LEN=6), ALLOCATABLE  :: typPHYS(:), &  ! array of symbolic domain types
                                   namePHYS(:)    ! array of domain names

END MODULE m_vxc2data
