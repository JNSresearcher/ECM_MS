program ECM_MS
! Fortran codes created by J.Sochor   ( https://github.com/JNSresearcher )
! e-mail: JNSresearcher@gmail.com
! June 2026

USE m_fparser
USE m_vxc2data
IMPLICIT NONE

type espm
    integer im
    REAL(8) em
    type (espm),pointer ::prec
end type espm
type(espm),pointer :: sp_col, next_col 

!---------------------------------------!
! input data from SUBROUTINE vxc2data   !
!---------------------------------------!
REAL(8) :: delta(3)         ! array grid spacing along X, Y and Z   
REAL(8) :: dt,   &          ! time step
           Time, &          ! stop time
           dtt              ! jump duration

INTEGER :: sdx,sdy,sdz      ! number of cells along X,  Y and  Z   
INTEGER :: nsub,     &      ! number of physical domains
           nsub_air, &      ! number of environment domains
           nsub_glob        ! total number of domains
INTEGER :: numfun,   &      ! number of functions for calculating external sources
           numMech,  &      ! number of functions for calculating movements of external sources
           numVenv          ! number of functions for calculating movements of ferromagnets

!//=======================declare SOR and nonlin-mag
REAL(8) ::  tol, omega, tolmag, omegamag
INTEGER ::  iter, itmax, itmag

! solv - character string corresponding to the name of methods: 'BCG' or 'PRD'
CHARACTER(LEN=3) :: solv 
CHARACTER(16):: files      ! name for output files 

! Ek - vector of sources of contour emf, Ik is the vector of contour currents, Iv is the vector of currents in branches
real(8), allocatable ::  Ek(:), Ik(:),  Iv(:)

integer::  num_k,   &    ! total number of contours
           num_kX,  &    ! number of contours along the axis X 
           num_kY,  &    ! number of contours along the axis Y 
           num_kZ,  &    ! number of contours along the axis Z 
           num_brX, &    ! number of branches along the axis X
           num_brY, &    ! number of branches along the axis Y
           num_brZ, &    ! number of branches along the axis Z
           num_br        ! total number of branches
integer :: sdx2,    &    ! number of nodes along the axis X
           sdy2,    &    ! number of nodes along the axis Y
           sdz2,    &    ! number of nodes along the axis Z
           Nodes,   &    ! total number of nodes
           Cells         ! number of cells

INTEGER :: i,j,k, L, m, n, nl, nn,ios, mm 
integer :: nod(1), kdz2,ip, jp, kp, nc, i_m,  i_p, j_p, k_p
real(8) :: a, Sx, Sy, Sz, MU_new,MU_sum, Z0,Zm,Zp, SmuX, SmuY, SmuZ

! indicators of movement of ferromagnets and windings
integer movestop_Fe(3),flag_move_Fe
integer d_i,d_j,d_k,flag_move, movestop(3), m_new, i_new, j_new, k_new, NcellsX, NcellsY, NcellsZ
! arrays for new cells of field sources as they move 
INTEGER, ALLOCATABLE ::  new_nodesX(:), new_nodesY(:), new_nodesZ(:) 

! output control
integer :: Nout_display,Nprint_display, Nout, Npoint, Nprint, Ntime
integer :: k0,k1,k2
REAL(8) :: T, Tcalc, Tsavedata

! recalculation of cell coordinates into contour coordinates
INTEGER, allocatable ::  cells2cont(:)

integer col_kv(4)
real(8) val_kv(4)
real(8), allocatable :: val_Ckv(:),  val_Cvk(:),  val_CZkv(:),   val_Zkk(:), Zvv(:)
integer, allocatable :: irow_Ckv(:), irow_Cvk(:), irow_CZkv(:) , irow_Zkk(:) 
integer, allocatable :: jcol_Ckv(:), jcol_Cvk(:), jcol_CZkv(:),  jcol_Zkk(:)
integer nz_Ckv, nz_Cvk,  nz_CZkv, nz_Zkk

! domain map with ferromagnets
INTEGER(1) n_Fe
INTEGER(1),ALLOCATABLE:: geoPHYS_Fe(:,:,:)

!  
type t_dom_Fe
    INTEGER, ALLOCATABLE  :: num_br(:,:), &  ! branch numbers in a cell, 1st index cell number, 2st index 1 for X, 2 for Y
                             nodes_cell(:)   ! cell numbers in the domain
    REAL(8) mu_middle, mu_middle_old         ! middle magnetic permeabilities for iterations
    INTEGER n_Fe, num_cells
    REAL(8), ALLOCATABLE :: B_dat(:), H_dat(:), MU_dat(:) ! tabular data B-H B-mu
    integer n_dat, nonlin, tip_nonlin
end type t_dom_Fe
TYPE (t_dom_Fe), ALLOCATABLE :: dom_Fe(:) ! sizdom_Fe transmitted through the module

!========for iterac 
real(8) ::  summ_mu_middle_old, summ_mu_middle

!================================END DECLARE====================


!----- START TIMER
call system_clock(k1,k0);

!------------------- input of initial data
CALL vxc2data ( delta,    dt,        Time,     dtt,                   &
                sdx,    sdy,      sdz,       nsub,     nsub_air, nsub_glob, &
                numfun, numMech,  numVenv,     &
                solv,   files,    tol, omega, itmax, tolmag, itmag, omegamag )  
                

sdx2 = sdx+1                         ! number of nodes along the axis X  
sdy2 = sdy+1                         ! number of nodes along the axis Y  
sdz2 = sdz+1                         ! number of nodes along the axis Z 
Nodes = sdx2*sdy2*sdz2               ! total number of nodes

Cells = sdx*sdy*sdz                  ! number of cells

kdz2 = sdx2 * sdy2 ! auxiliary variable

num_kX = sdx2*sdy*sdz               ! number of contours along the axis X  
num_kY = sdx*sdy2*sdz               ! number of contours along the axis Y  
num_kZ = sdx*sdy*sdz2               ! number of contours along the axis Z  
num_k = num_kx + num_ky + num_kz    ! total number of contours

print '( *(a,i7))', 'num_kx=', num_kx,' num_ky=',num_ky,' num_kz=', num_kz, ' num_k=',num_k

num_brX = sdx*sdy2*sdz2              ! number of branches along the axis X
num_brY = sdx2*sdy*sdz2              ! number of branches along the axis Y
num_brZ = sdx2*sdy2*sdz              ! number of branches along the axis Z
num_br = num_brX + num_brY +num_brZ  ! total number of branches

allocate (Zvv(num_br), source=1.d0) 
allocate (Ek(num_k), Ik(num_k), Iv(num_br), source=0.d0)

print '( *(a,i7))', 'num_brX=',num_brX, ' num_brY=',num_brY,' num_brZ=',num_brZ,' num_br=',num_br, ' num_nodes=',Nodes

if (sizdom_Fe /=0)  then
    ALLOCATE (geoPHYS_Fe(sdx,sdy,sdz),  source=0_1  ) 
    ALLOCATE ( dom_Fe(sizdom_Fe)  )  
    dom_Fe%nonlin = 0
    
    do n = 1, sizdom_Fe 
        n_Fe = numdom_Fe(n)
        dom_Fe(n)%n_Fe = int (n_Fe, 4)
        
        m=0
        do k=1, sdz; do j=1, sdy; do i=1, sdx
            if (geoPHYS(i,j,k) == n_Fe) then
                geoPHYS_Fe(i,j,k) = n_Fe
                ! counting the number of cells in each domain
                m = m+1
            endif
        enddo; enddo; enddo
        ! save 
        dom_Fe(n)%num_cells = m
        ! for each cell Fe the numbers of XYZ branches will be known
        ALLOCATE ( dom_Fe(n)%num_br(m,3), dom_Fe(n)%nodes_cell(m),  source = 0) 
        L = int(n_Fe, 4)
        dom_Fe(n)%mu_middle = valPHYS(L,1)
    enddo

    if (siz_non_lin_MU /= 0 ) then ! if nonlinearities are detected
        do n = 1, siz_non_lin_MU
            do m=1, sizdom_Fe
                if (dom_Fe(m)%n_Fe == non_lin_MU(n)%num_domain) then
                    dom_Fe(m)%nonlin = 1
                    if (non_lin_MU(n)%tip_non_lin == 'BH') then
                        dom_Fe(m)%tip_nonlin = 2
                    ELSEIF (non_lin_MU(n)%tip_non_lin == 'BM') then
                        dom_Fe(m)%tip_nonlin = 1
                    else
                        print*, 'error tip nonlin'
                        stop
                    endif
                    
                    open (1, file = non_lin_MU(n)%file_MU )
                    dom_Fe(m)%B_dat = [real(8):: ]
    
                    if (dom_Fe(m)%tip_nonlin == 2) then
                        dom_Fe(m)%H_dat = [real(8):: ]
                    else 
                        dom_Fe(m)%MU_dat = [real(8):: ]
                    endif
                    L=0
                    do
                        READ(1, *,  iostat=L ) Sx, Sy 
                        IF ( L /= 0 )  then
                            EXIT
                        else
                            dom_Fe(m)%B_dat = [dom_Fe(m)%B_dat, Sx]
                            if (dom_Fe(m)%tip_nonlin == 2) dom_Fe(m)%H_dat  = [dom_Fe(m)%H_dat, Sy]
                            if (dom_Fe(m)%tip_nonlin == 1) dom_Fe(m)%MU_dat = [dom_Fe(m)%MU_dat, Sy]
                        endif
                    END DO
                    dom_Fe(m)%n_dat = size(dom_Fe(m)%B_dat)
                    
                    close(1)
                    print '(a, a, a, i3, a,i3)',  'create tabl nonlin mu from file ',&
                    non_lin_MU(n)%file_MU, ' dom=', m, ' siz=', dom_Fe(m)%n_dat
                    
                    if (dom_Fe(m)%n_dat == 0) then
                        print*, 'ERROR read file ', non_lin_MU(n)%file_MU
                        stop
                    endif
                endif
            enddo
        enddo    
    endif  ! siz_non_lin_MU /= 0
    
!----------------------------------------------------------------------
! generation of information about the numbers of ferromagnetic branches
!----------------------------------------------------------------------
    do nn = 1, sizdom_Fe 
        ip=0; jp=0; kp=0; 
        L=0 ! Fe domain number counter
        do k = 1,sdz2;  do j = 1,sdy2;  do i = 1,sdx2
            m=0
            if ( k<sdz2 .and. j<sdy2 .and. i<sdx2) then
                n = i + sdx * ((j -1) + sdy*(k -1))
                i_new = i; j_new = j; k_new = k   ! coord glob nod
                m=1                               ! signal grid nods
            endif
            
            if ( i < sdx2  ) then
                ip = ip + 1
                i_p = ip                      ! num br X
            else
                i_p = 0
            endif 
            if ( j < sdy2) then
                jp = jp + 1
                j_p = jp + num_brX            ! num br Y
            else
                j_p = 0
            endif 
            if ( k < sdz2) then
                kp = kp + 1
                k_p = kp + num_brX + num_brY  ! num br Z
            else
                k_p = 0
            endif 
            
            if (m == 1) then                  ! we are in a grid of nodes
                if (geoPHYS_Fe (i_new,j_new,k_new) == dom_Fe(nn)%n_Fe ) then
                    if (i_p /=0 .and. j_p /=0 .and. k_p /=0) then
                        L = L+1                                    ! counter cells local
                        dom_Fe(nn)%num_br(L,1) = i_p               ! num br X
                        dom_Fe(nn)%num_br(L,2) = j_p               ! num br Y
                        dom_Fe(nn)%num_br(L,3) = k_p               ! num br Z
                        dom_Fe(nn)%nodes_cell(L) = n
                    endif
                endif
            endif
        end do;end do;end do
    enddo
endif  ! sizdom_Fe /=0

!-----------------------------------------------------------
!                FORMATION OF A CONTOUR MATRIX
!                Contour-branch matrix
!-----------------------------------------------------------
allocate (irow_Ckv(num_k+1), cells2cont(3*Cells),  source=0)
nullify(sp_col, next_col )

irow_Ckv(1) = 1
nz_Ckv = 0 ! non-zero element counter

! X contours in the YZ plane
nn=0
L=0
do i=1,sdx2; do k=1,sdz2; do j=1,sdy2
    if (j==sdy2 .or. k==sdz2) cycle  ! num_kx = sdx2*sdy*sdz   number of contours X 
    nn =  nn + 1                     ! current contour X

    if (i ==sdx2 ) then
    else
       n = i + sdx * ((j -1) + sdy*(k -1)) ! cell number
       ! array for converting cell coordinates to contour coordinates
       cells2cont(n) = nn   !  n - cell number, nn - contour number X
    endif

    col_kv = 0; val_kv = 0.d0
    jp  = i     + sdx2*(j+1-1) + sdx2*sdy2*(k-1) + num_brX + num_brY  !z
    nc  = i     + sdx2*(j-1)   + sdx2*sdy2*(k-1) + num_brX + num_brY  !-z
    kp  = i     + sdx2*(j-1)   + sdx2*sdy*(k+1-1) + num_brX           !-y
    ip  = i     + sdx2*(j-1)   + sdx2*sdy*(k-1)   + num_brX           !y
                !     Iz    -Iz      -Iy       Iy                        !  branches
    col_kv(1:4) =  [ jp,   nc,     kp,     ip]
    val_kv(1:4) =  [ 1.d0, -1.d0,  -1.d0,  1.d0   ]
    !--------------------------X------------------------------------
    ! sort branch numbers in ascending order
    call full_sort(col_kv, val_kv, 4, 1,1) 
    do m = 1, 4
        if (col_kv(m) == 0) then
            print*, 'Error: contour number=',nn, 'm',m,'col_kv(m)',col_kv(m)
            stop
        endif
        allocate(sp_col)
        sp_col=espm(col_kv(m),val_kv(m),next_col) ! COLUMN PACKAGING
        nz_Ckv = nz_Ckv + 1  ! non-zero element counter
        next_col => sp_col
    enddo
    ! pointer to the column index for the next row
    irow_Ckv(nn+1) = irow_Ckv(nn) + 4    
enddo; enddo; enddo

! Y contours in the XZ plane
do j=1,sdy2; do k=1,sdz2; do i=1, sdx2
    if (i==sdx2 .or. k==sdz2) cycle ! num_ky = sdx*sdy2*sdz  number of contours Y 
        nn =  nn + 1                ! current contour  Y
    if (j ==sdy2 ) then
    else
        n = i + sdx * ((j -1) + sdy*(k -1)) + Cells
        cells2cont(n) = nn   ! n - cell number, nn - contour number Y
    endif
    col_kv = 0; val_kv = 0.d0
    jp  = i+1   + sdx2*(j-1) + sdx2*sdy2*(k-1) + num_brX + num_brY !z
    nc  = i     + sdx2*(j-1) + sdx2*sdy2*(k-1) + num_brX + num_brY !z
    kp  = i     + sdx*(j-1)   + sdx*sdy2*(k-1)                     !x
    ip  = i     + sdx*(j-1)   + sdx*sdy2*(k+1-1)                   !x
                !   Iz       -Iz       Ix     -Ix                       !  branches
    col_kv(1:4) = [ nc,     jp ,    ip,    kp]
    val_kv(1:4) = [ 1.d0,  -1.d0,  1.d0,  -1.d0 ]
    !--------------------------Y------------------------------------
    ! sort branch numbers in ascending order
    call full_sort(col_kv, val_kv, 4, 1,1)
    do m = 1, 4
        if (col_kv(m) == 0) then
            print*, 'Error: contour number=',nn, 'm',m,'col_kv(m)',col_kv(m)
            stop
        endif
        allocate(sp_col)
        sp_col=espm(col_kv(m),val_kv(m),next_col)         ! COLUMN PACKAGING
        nz_Ckv = nz_Ckv + 1 ! non-zero element counter
        next_col => sp_col
    enddo
    ! pointer to the column index for the next row
    irow_Ckv(nn+1) = irow_Ckv(nn) + 4                     
enddo; enddo; enddo

!  Z contours in the XY plane
do k=1,sdz2; do j=1,sdy2; do i=1, sdx2
   if (i==sdx2 .or. j==sdy2) cycle ! num_kz = sdx*sdy*sdz2 number of contours Z 
    nn =  nn + 1                   ! current contour Z
    if (k ==sdz2 ) then
    else
       n = i + sdx * ((j -1) + sdy*(k -1)) + 2*Cells
       cells2cont(n) = nn   ! n - cell number, nn - contour number Z
    endif
    col_kv = 0; val_kv = 0.d0
    jp  = i    + sdx*(j+1-1)   + sdx*sdy2*(k-1)           !-x
    nc  = i    + sdx*(j-1)     + sdx*sdy2*(k-1)           !+x
    kp =  i    + sdx2*(j-1)   + sdx2*sdy*(k-1) + num_brX  !-y 
    ip  = i+1  + sdx2*(j-1)   + sdx2*sdy*(k-1) + num_brX  !+y 
                !    -Ix         Ix      -Iy    Iy             ! branches
    col_kv(1:4) = [  jp,       nc,     kp,    ip]
    val_kv(1:4) =    [ -1.d0,    1.d0,   -1.d0,  1.d0   ]
    !-------------------------Z------------------------------------
    ! sort branch numbers in ascending order
    call full_sort(col_kv, val_kv, 4, 1,1) 
    do m = 1, 4
        if (col_kv(m) == 0) then
            print*, 'Error: contour number=',nn, 'm',m,'col_kv(m)',col_kv(m)
            stop
        endif
        allocate(sp_col)
        sp_col=espm(col_kv(m),val_kv(m),next_col) ! COLUMN PACKAGING
        nz_Ckv = nz_Ckv + 1
        next_col => sp_col
    enddo
    ! pointer to the column index for the next row
    irow_Ckv(nn+1) = irow_Ckv(nn) + 4 
enddo; enddo; enddo

! rewriting from stack to matrix
allocate(jcol_Ckv(nz_Ckv), source=0)
allocate(val_Ckv(nz_Ckv),source=0.0d0)
i=nz_Ckv  ! number of nonzero elements
DO WHILE (associated(sp_col))
    jcol_Ckv(i) = sp_col%im  ! column number
    val_Ckv(i) = sp_col%em   ! value
    i=i-1
    sp_col => sp_col%prec
END DO
nullify(sp_col, next_col)

! END FORMATION OF A CONTOUR MATRIX
!-------------------------------------


Z0 = valPHYS(nsub_glob,1)         ! environmental resistance
if ( delta(1) == delta(2)  .and.   delta(1) == delta(3) ) then
    SmuX=delta(1)/(0.12566370964050292d-5*delta(1)*delta(1))
    SmuY=SmuX
    SmuZ=SmuX
    if (sizdom_Fe == 0) Zvv = Z0*SmuX ! if there are no ferromagnets
else
    SmuX=delta(1)/(0.12566370964050292d-5*delta(2)*delta(3))
    SmuY=delta(2)/(0.12566370964050292d-5*delta(1)*delta(3))
    SmuZ=delta(3)/(0.12566370964050292d-5*delta(2)*delta(1)) ! dL/(mu0 dS)
endif

if (sizdom_Fe /=0) then 
    if (siz_non_lin_MU == 0 ) then
        call calc_Zvv
    else
        Zvv=SmuX   ! only preliminary value, will be recalculated 
    endif
endif

! formation of the product matrix Ckv*Zvv = CZkv 
nz_CZkv  = nz_Ckv
nz_Cvk   = nz_Ckv

allocate( val_Cvk(nz_Cvk),                      val_CZkv(nz_CZkv),   source=0.d0)
allocate(jcol_Cvk(nz_Cvk), irow_Cvk(num_br+1), jcol_CZkv(nz_CZkv), irow_CZkv(num_k+1) , source=0)

irow_CZkv = irow_Ckv
jcol_CZkv = jcol_Ckv

! CZkv = Ckv * Zvv 
do j = 1, num_k
    do k = irow_Ckv(j), irow_Ckv(j+1) - 1
        val_CZkv(k) = val_Ckv(k) * Zvv(jcol_Ckv(k))
    end do
end do
            !-------------------------------
            ! transpose of the Ckv matrix
            ! Cvk = transpose (Ckv)
            !--------------------------------
call csr_transpose(num_k, num_br, val_Ckv, jcol_Ckv, irow_Ckv, val_Cvk, jcol_Cvk, irow_Cvk) 
            !----------------------------------------------------
            ! Estimating non-zero values ??for a product CZkv*Cvk
            !----------------------------------------------------
! (Ckv*Z)*Cvk  nrow=num_k   ncol=num_br   ncolb=num_k  ja,ia = jcol_CZkv,irow_CZkv;  jb ib =jcol_Cvk,irow_Cvk iw(num_k) 
allocate( irow_Zkk(num_k+1) , source=0)
call spgemm_symbolic(num_k, num_br,num_k,irow_CZkv,jcol_CZkv,irow_Cvk, jcol_Cvk,irow_Zkk, nz_Zkk)

! memory allocation for Zkk
allocate(val_Zkk(nz_Zkk), source=0.d0)
allocate(jcol_Zkk(nz_Zkk) , source=0)
            !--------------------------------------------
            ! Calculation of the contour matrix Zkk using 
            ! the product of matrices:Zkk = CZkv*Cvk
            !-------------------------------------------
call sp_mm (num_k, num_br, num_k, irow_CZkv,jcol_CZkv,val_CZkv, irow_Cvk,jcol_Cvk,val_Cvk, irow_Zkk,jcol_Zkk,val_Zkk, nz_Zkk )

do  j = 1, num_k
    ip = irow_Zkk(j); jp = irow_Zkk(j+1) - 1
    m = jp - ip + 1
    call full_sort( jcol_Zkk(ip:jp), val_Zkk(ip:jp), m, 1,1) 
enddo

PRINT '(a, i7, a,g10.3,a)', 'sparse matrix Zkk complet! non zero=',  &
                         nz_Zkk, ' density=', 100.0* REAL( nz_Zkk,4)/ (real( num_k,4)*real( num_k,4)),'%'

            !-----------------------------
            ! Formation of contour emf's
            !----------------------------
CALL initf (numfun) 
do i=1,numfun
    do m=1,Fun(i)%args
        if ( trim(Fun(i)%namex(m)) == 'T') then
            Fun(i)%velx(m)=T
        endif
    enddo
    CALL parsef (i, trim(fun(i)%eqn), fun(i)%namex(1:fun(i)%args) )
    fun(i)%vely = evalf (i, fun(i)%velx(1:fun(i)%args)) 
enddo

DO n=1,numfun
    a = Fun(n)%vely
    DO k = 1,fun_nod(n)%numnod_Fx          ! k - local number of cell
        m = fun_nod(n)%nods_Fx(k)          ! m - global number of cell
        i = cells2cont(m)                  ! i - contour number
        Ek(i) = a * fun_nod(n)%nods_Fx_cos(k) 
        
        m = fun_nod(n)%nods_Fy(k)
        i = cells2cont(m) 
        Ek(i) = a * fun_nod(n)%nods_Fy_cos(k) 
        
        m = fun_nod(n)%nods_Fz(k)
        i = cells2cont(m)  
        Ek(i) = a * fun_nod(n)%nods_Fz_cos(k) 
    ENDDO
ENDDO

!  creating a folder for results
call execute_command_line ('mkdir '//trim(files), exitstat=i)
PRINT '( a, i3 )', 'created dir '//trim(files)//': status=', i
    
call system_clock(k2);
Tsavedata =  real(k2-k1)/real(k0)
print '( "Data preparation time= ", g10.3)', Tsavedata

call system_clock(k1,k0);
movestop = 1    ! flag of local start-stop movement in all coordinates coils
movestop_Fe = 1 ! flag of local start-stop movement in all coordinates ferromag

!-------------------------
!  initial iterations SOR
!-------------------------
if (solv == 'SOR') then
    print*, 'start solver SOR'
    call sor_sparse (num_k, val_Zkk, jcol_Zkk, irow_Zkk, Ek, Ik, omega, itmax, tol, iter)
    print*, 'init iter SOR=', iter
    summ_mu_middle = 1.d7 ! init summ_mu_middle
else
    print*,"error: only 'sor'"
    stop
endif

!******** CYCLE PREPARATION by time
T=0.
Ntime=0                             ! counter of each point with step dt
Nout_display=nint(Time/(100.0*DT))  ! jump size for displaying a character '>' on the screen 
Nprint_display=Ntime + Nout_display ! point counter for displaying the '>' symbol with a jump Nout_display

Npoint=0                           ! calculation point counter with step dtt
Nout=nint(DTT/DT)                  ! jump size for output to files
Nprint = Ntime + Nout              ! point counter with step Nout for output to files with jump Nout

!----------------------------------------------------------------!
! preparation of arrays for calculating the movement of ferromag !
!----------------------------------------------------------------!
flag_move_Fe = 0
IF (numVenv /=0 ) THEN
    DO i=1,numVenv  
    ! motion search processing
        fun_Venv(i)%Distance = 0.d0 
        IF (fun_Venv(i)% num_Venv(1) == 0  .and. fun_Venv(i)%move(1) /=0) THEN ! it is not function
        !  distance traveled in time dt with speed Vx in fractions of cell size along dX
            fun_Venv(i)%shift(1) = fun_Venv(i)%val_Venv(1)*dt/delta(1)   !Vx*dt/dX
            flag_move_Fe = 1
        ELSEIF (  fun_Venv(i)%move(1) /=0) THEN 
            flag_move_Fe = 1

        ENDIF
        IF (fun_Venv(i)% num_Venv(2) == 0  .and. fun_Venv(i)%move(2) /=0) THEN ! it is not function
            fun_Venv(i)%shift(2) = fun_Venv(i)%val_Venv(2)*dt/delta(2)
            flag_move_Fe = 1; 
        ELSEIF (  fun_Venv(i)%move(2) /=0) THEN 
            flag_move_Fe = 1
        ENDIF
        IF (fun_Venv(i)% num_Venv(3) == 0 .and. fun_Venv(i)%move(3) /=0 ) THEN ! it is not function, but not a shift either
            fun_Venv(i)%shift(3) = fun_Venv(i)%val_Venv(3)*dt/delta(3)
            flag_move_Fe = 1; 
        ELSEIF (  fun_Venv(i)%move(3) /=0) THEN 
            flag_move_Fe = 1
        ENDIF
    ENDDO
ENDIF !numVenv

!---------------------------------------------------------------!
! preparation of arrays for calculating the movement of sources !
!---------------------------------------------------------------!
flag_move = 0 ! global flag of motion in input data

IF (numfun /=0 ) THEN
    IF (numMech /= 0) THEN  ! for init dispos coil
        CALL initf (numMech) 
        DO i=1,numMech
            DO m=1,Vmech(i)%args
                IF ( trim(Vmech(i)%namex(m)) == 'T') THEN
                    Vmech(i)%velx(m)=T
                ENDIF
            ENDDO
            CALL parsef (i, trim(Vmech(i)%eqn), Vmech(i)%namex(1:Vmech(i)%args) )
            Vmech(i)%vely = evalf (i, Vmech(i)%velx(1:Vmech(i)%args))
        ENDDO
    ENDIF

    DO n=1,numfun
    ! motion search processing
        fun_nod(n)%Distance = 0.d0 
        ! unlock
        do i = 1,3
            IF (fun_nod(n)% num_Vmech(i) == 0  .and. fun_nod(n)%move(i) /=0) THEN ! it is not function
                !  distance traveled in time dt with speed Vx in fractions of cell size along dX
                fun_nod(n)%shift(i) = fun_nod(n)%val_Vmech(i)*dt/delta(i)   !Vx*dt/dX
                flag_move = 1
            ELSEIF (  fun_nod(n)%move(i) /=0) THEN 
                flag_move = 1
                fun_nod(n)%Distance(i) = - (Vmech(fun_nod(n)%num_Vmech(i))%vely * 1.5d0*dt/ delta(i) )
            ENDIF
        enddo
    ENDDO
! initial allocation of array memory for source nodes, regardless of whether there is movement or not
    NcellsX = 0; NcellsY = 0; NcellsZ = 0; 
    DO n=1,numfun
        NcellsX = NcellsX + fun_nod(n)%numnod_Fx
        NcellsY = NcellsY + fun_nod(n)%numnod_Fy
        NcellsZ = NcellsZ + fun_nod(n)%numnod_Fz
    ENDDO
    IF (NcellsX /= 0) ALLOCATE (new_nodesX(NcellsX), source=0)
    IF (NcellsY /= 0) ALLOCATE (new_nodesY(NcellsY), source=0)
    IF (NcellsZ /= 0) ALLOCATE (new_nodesZ(NcellsZ), source=0)
ENDIF !numfun

! filling source node arrays if there is no movement
IF (flag_move == 0 .and. numfun /=0 ) THEN
    NcellsX = 0; NcellsY = 0; NcellsZ = 0; ! counters
    DO n=1,numfun
        DO k = 1,fun_nod(n)%numnod_Fx
            m = fun_nod(n)%nods_Fx(k)    ! k - local number m - global number 
            NcellsX = NcellsX + 1
            new_nodesX(NcellsX) = m
            
            m = fun_nod(n)%nods_Fy(k)
            NcellsY = NcellsY + 1
            new_nodesY(NcellsY) = m 
            
            m = fun_nod(n)%nods_Fz(k)
            NcellsZ = NcellsZ + 1
            new_nodesZ(NcellsZ) = m
        ENDDO
    ENDDO

ENDIF

            !----------------------------------
            !    START OF CALCULATION
            !----------------------------------
2000 CONTINUE   

! calculation of the source function
CALL initf (numfun) 
do i=1,numfun
    do m=1,Fun(i)%args
        if ( trim(Fun(i)%namex(m)) == 'T') then
            Fun(i)%velx(m)=T
        endif
    enddo
    CALL parsef (i, trim(fun(i)%eqn), fun(i)%namex(1:fun(i)%args) )
    fun(i)%vely = evalf (i, fun(i)%velx(1:fun(i)%args)) 
enddo

!-------------------------------------------------!
! Calculation of the speed of movement of sources ! 
!-------------------------------------------------!
! Remember old states
IF (numMech /= 0) THEN  ! Calculating new states
   ! fix fun_nod(n)%ism(i) based on the results of a comparison of new and old speeds
    DO n=1,numfun
        DO i=1,3 
            if (fun_nod(n)%num_Vmech(i) /=0) fun_nod(n)%val_Vmech_old(i) = Vmech(fun_nod(n)%num_Vmech(i))%vely
        enddo
    enddo
    ! calculating the new speed
    CALL initf (numMech) 
    DO n=1,numMech
        DO m=1,Vmech(n)%args
            IF ( trim(Vmech(n)%namex(m)) == 'T') THEN
                Vmech(n)%velx(m)=T
            ENDIF
        ENDDO
        CALL parsef (n, trim(Vmech(n)%eqn), Vmech(n)%namex(1:Vmech(n)%args) )
        Vmech(n)%vely = evalf (n, Vmech(n)%velx(1:Vmech(n)%args)) ! new 
    enddo
    ! remember the new speed
    DO n=1,numfun
        DO i=1,3 
            if (fun_nod(n)%num_Vmech(i) /=0) fun_nod(n)%val_Vmech(i) = Vmech(fun_nod(n)%num_Vmech(i))%vely
        enddo
        ! comparison:
        do i = 1,3
            if ( ( fun_nod(n)%val_Vmech(i) <= 0 .and. fun_nod(n)%val_Vmech_old(i) <= 0) .or. &
                 ( fun_nod(n)%val_Vmech(i) >= 0 .and. fun_nod(n)%val_Vmech_old(i) >= 0) ) then
                fun_nod(n)%ism(i) = 0
            elseif ( fun_nod(n)%val_Vmech(i) < 0 .and. fun_nod(n)%val_Vmech_old(i) > 0) then
                fun_nod(n)%ism(i) = 1
            elseif ( fun_nod(n)%val_Vmech(i) > 0 .and. fun_nod(n)%val_Vmech_old(i) < 0) then
                fun_nod(n)%ism(i) = -1
            endif
        enddo
    enddo
endif

!--------------------------------------------------!
! Distribution of independent sources into domains !
!--------------------------------------------------!
IF (flag_move == 1) THEN 
     Ek=0.d0
    ! calculation of new coordinates of nodes of independent sources 
    IF (numfun /= 0) THEN
        NcellsX=0; NcellsY=0; NcellsZ=0
        DO n=1,numfun
            DO i=1,3
                IF (fun_nod(n)% num_Vmech(i) == 0 ) THEN     !
                    fun_nod(n)%Distance(i) = fun_nod(n)%Distance(i) + movestop(i)*fun_nod(n)%shift(i)
                    fun_nod(n)%length(i) = nint(fun_nod(n)%Distance(i))
                ELSE                                         ! 
                    fun_nod(n)%Distance(i) = fun_nod(n)%Distance(i) +  Vmech(fun_nod(n)%num_Vmech(i))%vely*dt/delta(i)
                    if (fun_nod(n)%ism(i) == 0) then
                        fun_nod(n)%length(i) = nint(fun_nod(n)%Distance(i)  ) 
                    endif
                ENDIF
            ENDDO

            a = Fun(n)%vely
            DO k = 1,fun_nod(n)%numnod_Fx
                ! X
                m = fun_nod(n)%nods_Fx(k)
                CALL new_m (m)
                i = cells2cont(m)   ! m - glob number i - contour number
                Ek(i) = a * fun_nod(n)%nods_Fx_cos(k)
                NcellsX = NcellsX + 1
                new_nodesX(NcellsX) = m
                !  Y
                m = fun_nod(n)%nods_Fy(k)-Cells
                CALL new_m( m ) 
                m = m + Cells
                i = cells2cont(m)
                Ek(i) = a * fun_nod(n)%nods_Fy_cos(k) 
                NcellsY = NcellsY + 1
                new_nodesY(NcellsY) = m 
                ! Z
                m = fun_nod(n)%nods_Fz(k)- 2*Cells
                CALL new_m( m )
                m = m + 2*Cells
                i = cells2cont(m)  
                Ek(i) = a * fun_nod(n)%nods_Fz_cos(k) 
                NcellsZ = NcellsZ + 1
                new_nodesZ(NcellsZ) = m 
            ENDDO                  
        ENDDO
    ENDIF
ELSE
    !--------------------!
    ! no moving sources  !
    !--------------------!
    DO n=1,numfun
        a = Fun(n)%vely
        DO k = 1,fun_nod(n)%numnod_Fx
            m = fun_nod(n)%nods_Fx(k)                ! k - local number m - global number
            i = cells2cont(m)                        ! m - glob number i - contour number
            Ek(i) = a * fun_nod(n)%nods_Fx_cos(k) 
            m = fun_nod(n)%nods_Fy(k)
            i = cells2cont(m) 
            Ek(i) = a * fun_nod(n)%nods_Fy_cos(k) 
            m = fun_nod(n)%nods_Fz(k) 
            i = cells2cont(m)  
            Ek(i) = a * fun_nod(n)%nods_Fz_cos(k) 
        ENDDO
    ENDDO
ENDIF

if (siz_non_lin_MU /= 0 ) then ! only for those domains where there are nonlinear mu
            !------------------------------------------
            ! Start iterations for nonlinear regions
            !------------------------------------------
    do mm = 1, itmag
        summ_mu_middle_old = summ_mu_middle  ! sum of middle mu for each domain
        
        !==================================SOLVE======================
        call sor_sparse (num_k, val_Zkk, jcol_Zkk, irow_Zkk, Ek, Ik, omega, itmax, tol, iter)
        !==================================SOLVE======================
    
        ! Iv = Cvk * Ik; 
        do  i=1, num_br
            ip=irow_Cvk(i); jp=irow_Cvk(i+1)-1
            Iv(i) =  dot_product( val_Cvk(ip:jp), Ik(jcol_Cvk(ip:jp) ) )
        enddo 
    
        summ_mu_middle =0.d0
    !   recount
        do n = 1, sizdom_Fe
            if (dom_Fe(n)%nonlin == 0 ) exit   ! for linear domains skip
            MU_sum = 0.d0
            do m = 1, dom_Fe(n)%num_cells   ! number of cells in a domain
                ! for each m-th cell, branch numbers
                i = dom_Fe(n)%num_br(m,1) ! num br X
                j = dom_Fe(n)%num_br(m,2) ! num br Y
                k = dom_Fe(n)%num_br(m,3) ! num br Z
            
                Sx = Iv(i)
                Sy = Iv(j)
                Sz = Iv(k)
                a = sqrt(Sx*Sx + Sy*Sy + Sz*Sz)  ! amplitude of the induction flux density
                !--------------------------------------------------------
                ! Calculation of magnetic permeability from tabular data 
                if (dom_Fe(n)%tip_nonlin == 1) &
                MU_new = get_mu_from_b_tab(a, dom_Fe(n)%b_dat, dom_Fe(n)%MU_dat, dom_Fe(n)%n_dat) 

                if (dom_Fe(n)%tip_nonlin == 2) &
                MU_new = get_mu_from_bh_tab(a, dom_Fe(n)%b_dat, dom_Fe(n)%h_dat, dom_Fe(n)%n_dat) 
                !--------------------------------------------------------
                MU_sum = MU_sum + MU_new
            enddo
            dom_Fe(n)%mu_middle = real(dom_Fe(n)%num_cells ,8)/ MU_sum    ! save new 1/mu_new
            
            if ( omegamag /= 0.d0) then              ! if the under-relaxation is given
                dom_Fe(n)%mu_middle = (1.d0 - omegamag)*dom_Fe(n)%mu_middle + omegamag*dom_Fe(n)%mu_middle_old
                dom_Fe(n)%mu_middle_old = dom_Fe(n)%mu_middle
            endif
            summ_mu_middle = summ_mu_middle + dom_Fe(n)%mu_middle
        enddo
    
        call calc_Zvv2  ! recalculation of branch resistances 

        ! CZkv = Ckv * Zvv
        val_CZkv = 0.d0
        do j = 1, num_k
            do k = irow_Ckv(j), irow_Ckv(j+1) - 1
                val_CZkv(k) = val_Ckv(k) * Zvv(jcol_Ckv(k))
            end do
        end do

        ! calculation of the matrix of contour resistances
        ! Zkk = CZkv*Cvk
        val_Zkk = 0.d0
        call sp_mm (num_k, num_br, num_k, irow_CZkv,jcol_CZkv,val_CZkv, &
               irow_Cvk,jcol_Cvk,val_Cvk, irow_Zkk,jcol_Zkk,val_Zkk, nz_Zkk )
        do  j = 1, num_k
            ip = irow_Zkk(j); jp = irow_Zkk(j+1) - 1
            m = jp - ip + 1
            call full_sort( jcol_Zkk(ip:jp), val_Zkk(ip:jp), m, 1,1) 
        enddo

        summ_mu_middle = summ_mu_middle/real(sizdom_Fe,8)

        a = abs(summ_mu_middle - summ_mu_middle_old)/ summ_mu_middle_old 

        print  '( a, i4, *(a, e10.3 ))', 'num nonlin iter=', mm,  ' err=', a, ' mu_middle=',1.d0/summ_mu_middle
        if ( a < tolmag ) then  !.and. mm>3
            ! print*,'total iter=', mm
            exit
        else
            if ( mm==itmag ) print  '( a, i4, 2(a, e10.3) )','num nonlin iter=', mm, &
                         'summ_mu_middle', summ_mu_middle , 'summ_mu_middle_old', summ_mu_middle_old 
        endif
    enddo ! end iterax
    
else  ! there are no nonlinear regions
    
    call sor_sparse (num_k, val_Zkk, jcol_Zkk, irow_Zkk, Ek, Ik, omega, itmax, tol, iter)
    ! calculation of currents in branches
    do  i=1, num_br
        ip=irow_Cvk(i); jp=irow_Cvk(i+1)-1
        Iv(i) =  dot_product( val_Cvk(ip:jp), Ik(jcol_Cvk(ip:jp) ) )
    enddo 
endif

!------------------------------------------------------!
! Calculation of the speed of movement of ferromagnets !
!------------------------------------------------------!    
IF (numVenv /= 0) THEN
   ! fix fun_Venv((n)%ism(i) based on the results of a comparison of new and old speeds
    CALL initf (numVenv)
    DO n=1,numVenv
        DO i=1,3 
            if (fun_Venv(n)%num_Venv(i) /=0 ) fun_Venv(n)%val_Venv_old(i) = Venv(fun_Venv(n)%num_Venv(i))%vely
        enddo
        DO m=1,Venv(n)%args
            IF ( trim(Venv(n)%namex(m)) == 'T') THEN
                Venv(n)%velx(m)=T
            ENDIF
        ENDDO
        CALL parsef (n, trim(Venv(n)%eqn), Venv(n)%namex(1:Venv(n)%args) )
        Venv(n)%vely = evalf (n, Venv(n)%velx(1:Venv(n)%args)) ! new
       ! remember the new speed
        DO i=1,3 
            if (fun_Venv(n)%num_Venv(i) /=0 ) fun_Venv(n)%val_Venv(i) = Venv(fun_Venv(n)%num_Venv(i))%vely
        enddo
        !  comparison::
        do i = 1,3
            if ( ( fun_Venv(n)%val_Venv(i) <= 0 .and. fun_Venv(n)%val_Venv_old(i) <= 0) .or. &
                 ( fun_Venv(n)%val_Venv(i) >= 0 .and. fun_Venv(n)%val_Venv_old(i) >= 0) ) then
                fun_Venv(n)%ism(i) = 0
            elseif ( fun_Venv(n)%val_Venv(i) < 0 .and. fun_Venv(n)%val_Venv_old(i) > 0) then
                fun_Venv(n)%ism(i) = 1
            elseif ( fun_Venv(n)%val_Venv(i) > 0 .and. fun_Venv(n)%val_Venv_old(i) < 0) then
                fun_Venv(n)%ism(i) = -1
            endif
        enddo
    enddo  
ENDIF

IF (flag_move_Fe == 1) then   ! The motion of a ferromagnet has been established.

    do n = 1, numVenv  ! in all areas of ferromagnets that are moving
        !-------------
        n_Fe = Venv(n)%nomsch ! moving domain number
        
        ! Before moving, store the numbers of the moving Fe in geoPHYS(i,j,k)
        do k=1,sdz; do j=1,sdy; do i=1,sdx;
            if (geoPHYS_Fe(i,j,k) == n_Fe ) then
                geoPHYS(i,j,k) = n_Fe   !  In geoPHYS (ijk), save the Fe domain number
            ! In geoPHYS_Fe(ijk_n) write air instead of ferromagnet.
                geoPHYS_Fe (i,j,k) = valPHYS(nsub_glob,1) ! 0_1 
            endif
        enddo; enddo; enddo  
        
        !----------calculation of new coordinates
        DO i=1,3
            IF (fun_Venv(n)% num_Venv(i) == 0 ) THEN     !
                fun_Venv(n)%Distance(i) =  movestop_Fe(i)* fun_Venv(n)%shift(i)
                fun_Venv(n)%length(i) = nint(fun_Venv(n)%Distance(i))
            ELSE! 
                fun_Venv(n)%Distance(i) = Venv(fun_Venv(n)%num_Venv(i))%vely * dt/delta(i) ! v *dt/dx
                fun_Venv(n)%length(i) = nint(  fun_Venv(n)%Distance(i)  ) 
            ENDIF
        ENDDO

        if (fun_Venv(n)%ism(1) == 0 ) then
            d_i = fun_Venv(n)%length(1)
        elseif (fun_Venv(n)%ism(1) == 1 ) then  ! X turned back
            d_i = fun_Venv(n)%length(1) - 2
        elseif (fun_Venv(n)%ism(1) == -1 ) then
            d_i = fun_Venv(n)%length(1) + 2
        endif
        
        if (fun_Venv(n)%ism(2) == 0 ) then
            d_j = fun_Venv(n)%length(2) 
        elseif (fun_Venv(n)%ism(2) == 1 ) then  ! Y turned back
            d_j = fun_Venv(n)%length(2) - 2
        elseif (fun_Venv(n)%ism(2) == -1 ) then
            d_j = fun_Venv(n)%length(2) + 2
        endif
        
        if (fun_Venv(n)%ism(3) == 0 ) then 
            d_k = fun_Venv(n)%length(3)
        elseif (fun_Venv(n)%ism(3) == 1 ) then  ! Z turned back
            d_k = fun_Venv(n)%length(3) - 2
        elseif (fun_Venv(n)%ism(3) == -1 ) then
            d_k = fun_Venv(n)%length(3) + 2
        endif
        
        ! The first step is to calculate new coordinates for the area, 
        ! save them to the buffer, and erase the old geoPHYS area.
        do k=1,sdz; do j=1,sdy; do i=1,sdx; 
            ! if geoPHYS (ijk) == Fe region number
            if (geoPHYS (i,j,k) == n_Fe ) then
                ! then find new ijk
                call new_ijk( i,j,k, d_i, d_j, d_k, movestop_Fe,  m_new,i_new,j_new,k_new )
                ! in geoPHYS_Fe(ijk_n) write the Fe region number
                 geoPHYS_Fe(i_new,j_new,k_new) = n_Fe
                ! In geoPHYS (ijk) enter the air number
                geoPHYS (i,j,k) = valPHYS(nsub_glob,1)
            endif
        enddo; enddo; enddo
        
        ! overwrite from buffer in geoPHYS
        do k=1,sdz; do j=1,sdy; do i=1,sdx;
            ! if geoPHYS (ijk) == Fe region number
            if (geoPHYS_Fe (i,j,k) == n_Fe ) then
                geoPHYS (i,j,k) = geoPHYS_Fe (i,j,k)
            !  In geoPHYS_Fe (ijk_new) save the air number
                geoPHYS_Fe (i,j,k) = valPHYS(nsub_glob,1) 
            endif
        enddo; enddo; enddo  
    enddo

    !--------------------------
    ! Recalculation of branch resistances
    !--------------------------
    call calc_Zvv

    !Recalculation of the contours resistance matrix
    ! CZkv = Ckv * Zvv
    val_CZkv(k) = 0.d0
    do j = 1, num_k
        do k = irow_Ckv(j), irow_Ckv(j+1) - 1
            val_CZkv(k) = val_Ckv(k) * Zvv(jcol_Ckv(k))
        end do
    end do
    ! Zkk = CZkv*Cvk
    val_Zkk = 0.d0
    call sp_mm (num_k, num_br, num_k, irow_CZkv,jcol_CZkv,val_CZkv, &
               irow_Cvk,jcol_Cvk,val_Cvk, irow_Zkk,jcol_Zkk,val_Zkk, nz_Zkk )
    do  j = 1, num_k
        ip = irow_Zkk(j); jp = irow_Zkk(j+1) - 1
        m = jp - ip + 1
        call full_sort( jcol_Zkk(ip:jp), val_Zkk(ip:jp), m, 1,1) 
    enddo
    
endif  !flag_move_Fe
! ----------------END MOVING-------------------

            !--------------------------------------------
            !  END OF CALCULATION AND OUTPUT OF RESULTS
            !---------------------------------------------
if ( Ntime >= Nprint .or. Ntime == 0 ) then  
    Nprint = Ntime  + Nout
    Npoint= Npoint + 1                   ! point counter with dtt step
    ios=0
    
    call writeVtk (Npoint, sdx, sdy, sdz, num_brX, num_brY, num_kx, num_ky, delta, Iv,   files) !Ev, 

    IF (flag_move == 1 .or. Npoint == 1) &
        CALL writeVtk_src ( Npoint, numfun,  NcellsX, NcellsY, NcellsZ,  new_nodesX, new_nodesY,  new_nodesZ,  &
                            sdx, sdy, delta, files)  
    if (mod(Npoint,5) == 0) then
        print '( i0,$ )',  Npoint
    else 
        print '( a,$ )', '>'
    endif
endif

If (Ntime >= Nprint_display) then
    Nprint_display = Ntime + Nout_display
endif
Ntime = Ntime + 1   ! this is every point
T = T + DT

! If nonlinearities are present, then the vector 
! of the contour current is zeroed for new iterations
if (siz_non_lin_MU /= 0 ) Ik=0.0

if (T < Time) goto 2000
print*,'|'

call system_clock(k2);
Tcalc =  real(k2-k1)/real(k0)
print '( a, g10.3)', 'solve complet. Tcalc=',  Tcalc 
!=========================================

contains


SUBROUTINE calc_Zvv
! Calculation of branch resistances

nn = 0
! X contours in the YZ plane
do i=1, sdx2; do k=1, sdz2; do j=1, sdy2
    if (j==sdy2 .or. k==sdz2) cycle  ! num_kx = sdx2*sdy*sdz   number of contours X 
    nn =  nn + 1                     ! current contour X
    if (i /= sdx2 ) then
        NL = i + sdx * ( (j -1) + sdy*(k -1) )
        n = int(geoPHYS(i,j,k),4)
        m = int(geoPHYS_Fe (i,j,k),4)
        if (m /= 0 ) then            
            nod = FINDLOC(dom_Fe%n_Fe, m)
            L = nod(1)
            if (L /= 0) then
                if (dom_Fe(L)%nonlin == 1 ) cycle
            endif
        endif

        Z0 = valPHYS(n,1)   ! 1/mu_r

        if (i /= sdx) then
            i_p = int(geoPHYS(i+1,j,k),4)  ! domain number
            Zp = valPHYS(i_p,1)            ! for comparison only
        endif
        if (i == 1  ) then
            Zm = 1.d10
        else
            i_m = int(geoPHYS(i-1,j,k),4)
            Zm = valPHYS(i_m,1) 
        endif
    endif

    if ( Z0 <= Zm ) then 
    ! branches of the current face
        jp  = i + sdx2*(j+1-1) + sdx2*sdy2*(k-1) + num_brX + num_brY  !z
        nc  = i + sdx2*(j-1)   + sdx2*sdy2*(k-1) + num_brX + num_brY  !-z
    
        kp  = i + sdx2*(j-1)   + sdx2*sdy*(k+1-1) + num_brX !-y
        ip  = i + sdx2*(j-1)   + sdx2*sdy*(k-1)   + num_brX !y
               ! Iz    -Iz      -Iy     Iy         !  branches
               ! jp,   nc,     kp,     ip 
        Zvv(jp) = Z0 * SmuZ ! z  1/mu * dL/(mu0 dS)
        Zvv(nc) = Z0 * SmuZ !-z
        Zvv(kp) = Z0 * Smuy !-y
        Zvv(ip) = Z0 * Smuy ! y
    endif
    
    if ( Z0 < Zp .and. i /= sdx) then
        ! on the right there is a different material, which means the same resistances 
        ! for the right edge should be assigned to the pairs of branches
        jp  = i+1 + sdx2*(j+1-1) + sdx2*sdy2*(k-1) + num_brX + num_brY  !z
        nc  = i+1 + sdx2*(j-1)   + sdx2*sdy2*(k-1) + num_brX + num_brY  !-z
    
        kp  = i+1 + sdx2*(j-1)   + sdx2*sdy*(k+1-1) + num_brX           !-y
        ip  = i+1 + sdx2*(j-1)   + sdx2*sdy*(k-1)   + num_brX           !y
            ! Iz    -Iz      -Iy  Iy         !  branches
            ! jp,   nc,     kp,   ip 
        Zvv(jp) = Z0 * SmuZ ! z
        Zvv(nc) = Z0 * SmuZ !-z
        Zvv(kp) = Z0 * Smuy !-y
        Zvv(ip) = Z0 * Smuy ! y
    endif
enddo; enddo; enddo

! Y contours in the XZ plane
do j=1, sdy2; do k=1, sdz2; do i=1, sdx2
    if (i==sdx2 .or. k==sdz2) cycle !  num_ky = sdx*sdy2*sdz ! number of contours Y 
    nn =  nn + 1     ! current contour  Y
    if (j /= sdy2 ) then
        NL = i + sdx * ( (j -1) + sdy*(k -1) )
        n = int(geoPHYS(i,j,k),4)
        
        m = int(geoPHYS_Fe (i,j,k),4)
        if (m /= 0 ) then            
            nod = FINDLOC(dom_Fe%n_Fe, m)
            L = nod(1)
            if (L /= 0) then
                if (dom_Fe(L)%nonlin == 1 ) cycle
            endif
        endif
        
        Z0 = valPHYS(n,1)
        
        if (j /= sdy) then
            i_p = int(geoPHYS(i,j+1,k),4)
            Zp = valPHYS(i_p,1)
        endif
        if (j == 1  ) then
            Zm = 1.d20
        else
            i_m = int(geoPHYS(i,j-1,k),4)
            Zm = valPHYS(i_m,1)
        endif
    endif

    if ( Z0 <= Zm ) then 
        ! branches of the current face
        jp  = i+1   + sdx2*(j-1) + sdx2*sdy2*(k-1) + num_brX + num_brY !z
        nc  = i     + sdx2*(j-1) + sdx2*sdy2*(k-1) + num_brX + num_brY !z
    
        kp  = i     + sdx*(j-1)   + sdx*sdy2*(k-1)   !x
        ip  = i     + sdx*(j-1)   + sdx*sdy2*(k+1-1) !x
        !   Iz       -Iz    Ix    -Ix         !  branches
        !   nc,     jp ,    ip,    kp 
        Zvv(jp) = Z0 * SmuZ! z    z- branches have already been calculated, can  skip them?
        Zvv(nc) = Z0 * SmuZ!-z
        Zvv(kp) = Z0 * SmuX!-x
        Zvv(ip) = Z0 * SmuX! x 
    endif
    
    if ( Z0 < Zp .and. j /= sdy) then
    
        jp  = i+1   + sdx2*(j+1-1) + sdx2*sdy2*(k-1) + num_brX + num_brY !z
        nc  = i     + sdx2*(j+1-1) + sdx2*sdy2*(k-1) + num_brX + num_brY !-z
    
        kp  = i     + sdx*(j+1-1)   + sdx*sdy2*(k-1)   !-x
        ip  = i     + sdx*(j+1-1)   + sdx*sdy2*(k+1-1) !x
        !   Iz       -Iz       Ix  -Ix         !  branches
        !   nc,     jp ,    ip,    kp 
        Zvv(jp) = Z0 * SmuZ ! z
        Zvv(nc) = Z0 * SmuZ !-z
        Zvv(kp) = Z0 * SmuX !-x
        Zvv(ip) = Z0 * SmuX ! x
    endif
enddo;enddo;enddo

end SUBROUTINE calc_Zvv


!----------------------
SUBROUTINE calc_Zvv2
!Calculation of the resistance of ferromagnet branches during iterations
nn = 0
! X contours in the YZ plane
do i=1, sdx2; do k=1, sdz2; do j=1, sdy2
    Z0 = valPHYS(nsub_glob,1)  

    if (j==sdy2 .or. k==sdz2) cycle  !num_kx = sdx2*sdy*sdz   number of contours X 
    nn =  nn + 1                     ! current contour X

    if (i /= sdx2 ) then
        NL = i + sdx * ( (j -1) + sdy*(k -1) ) ! 
        n = int(geoPHYS(i,j,k), 4)
        Z0 = valPHYS(n,1)   ! 1/mu_r
        
        do m = 1, sizdom_Fe 
            if (dom_Fe(m)%nonlin == 0 ) exit  
            if (dom_Fe(m)%n_Fe == n) then
                Z0 = dom_Fe(m)%mu_middle ! domain found
            endif
        enddo

        if (i /= sdx) then
            i_p = int(geoPHYS(i+1,j,k),4)  ! domain number
            Zp = valPHYS(i_p,1)            ! for comparison only
        endif
        if (i == 1  ) then
            Zm = 1.d10
        else
            i_m = int(geoPHYS(i-1,j,k),4)
            Zm = valPHYS(i_m,1) 
        endif
    endif

    if ( Z0 <= Zm ) then 
    ! branches of the current face
        jp  = i + sdx2*(j+1-1) + sdx2*sdy2*(k-1) + num_brX + num_brY  !z
        nc  = i + sdx2*(j-1)   + sdx2*sdy2*(k-1) + num_brX + num_brY  !-z
    
        kp  = i + sdx2*(j-1)   + sdx2*sdy*(k+1-1) + num_brX !-y
        ip  = i + sdx2*(j-1)   + sdx2*sdy*(k-1)   + num_brX !y
               !  Iz    -Iz      -Iy      Iy         !  branches
               !   jp,   nc,     kp,     ip 
        Zvv(jp) = Z0 * SmuZ ! z  1/mu * dL/(mu0 dS) 
        Zvv(nc) = Z0 * SmuZ !-z 
        Zvv(kp) = Z0 * Smuy !-y 
        Zvv(ip) = Z0 * Smuy ! y 
    endif
    
    if ( Z0 < Zp .and. i /= sdx) then
        ! on the right there is a different material, which means the same resistances 
        ! for the right edge should be assigned to the pairs of branches
        jp  = i+1 + sdx2*(j+1-1) + sdx2*sdy2*(k-1) + num_brX + num_brY  !z
        nc  = i+1 + sdx2*(j-1)   + sdx2*sdy2*(k-1) + num_brX + num_brY  !-z
    
        kp  = i+1 + sdx2*(j-1)   + sdx2*sdy*(k+1-1) + num_brX !-y
        ip  = i+1 + sdx2*(j-1)   + sdx2*sdy*(k-1)   + num_brX !y
               !  Iz    -Iz     -Iy    Iy         !  branches
               !  jp,   nc,     kp,     ip 
        Zvv(jp) = Z0 * SmuZ    ! z
        Zvv(nc) = Z0 * SmuZ    !-z
        Zvv(kp) = Z0 * Smuy    !-y
        Zvv(ip) = Z0 * Smuy    ! y
    endif
enddo; enddo; enddo

! Y contours in the XZ plane
do j=1, sdy2; do k=1, sdz2; do i=1, sdx2
    if (i==sdx2 .or. k==sdz2) cycle !  num_ky = sdx*sdy2*sdz ! number of contours Y 
    nn =  nn + 1     ! current contour  Y
    
    if (j /= sdy2 ) then
        NL = i + sdx * ( (j -1) + sdy*(k -1) )
        n = int(geoPHYS(i,j,k),4)
        Z0 = valPHYS(n,1)
        
        do m = 1, sizdom_Fe 
            if (dom_Fe(m)%nonlin == 0 ) exit 
            if (dom_Fe(m)%n_Fe == n) Z0 = dom_Fe(m)%mu_middle !domain found
        enddo
        
        if (j /= sdy) then
            i_p = int(geoPHYS(i,j+1,k),4)
            Zp = valPHYS(i_p,1)
        endif
        if (j == 1  ) then
            Zm = 1.d20
        else
            i_m = int(geoPHYS(i,j-1,k),4)
            Zm = valPHYS(i_m,1)
        endif
    endif

    if ( Z0 <= Zm ) then 
        ! branches of the current face
        jp  = i+1   + sdx2*(j-1) + sdx2*sdy2*(k-1) + num_brX + num_brY !z
        nc  = i     + sdx2*(j-1) + sdx2*sdy2*(k-1) + num_brX + num_brY !z
    
        kp  = i     + sdx*(j-1)   + sdx*sdy2*(k-1)   !-x
        ip  = i     + sdx*(j-1)   + sdx*sdy2*(k+1-1) !x
        !   Iz      -Iz     Ix   -Ix         !  branches
        !   nc,     jp ,    ip,    kp 
        Zvv(jp) = Z0 * SmuZ     ! z   ! already calculated,
        Zvv(nc) = Z0 * SmuZ     !-z   ! If uncomment, the current values will increase
        Zvv(kp) = Z0 * SmuX     !-x
        Zvv(ip) = Z0 * SmuX     ! x
    endif
    
    if ( Z0 < Zp .and. j /= sdy) then
        jp  = i+1   + sdx2*(j+1-1) + sdx2*sdy2*(k-1) + num_brX + num_brY !z
        nc  = i     + sdx2*(j+1-1) + sdx2*sdy2*(k-1) + num_brX + num_brY !z
    
        kp  = i     + sdx*(j+1-1)   + sdx*sdy2*(k-1)   !x
        ip  = i     + sdx*(j+1-1)   + sdx*sdy2*(k+1-1) !x
          ! Iz       -Iz    Ix    -Ix         !  branches
          ! nc,     jp ,    ip,    kp 
        ! Zvv(jp) = Z0 * SmuZ     ! z
        ! Zvv(nc) = Z0 * SmuZ     !-z
        Zvv(kp) = Z0 * SmuX     !-x
        Zvv(ip) = Z0 * SmuX     ! x
    endif
enddo; enddo; enddo

end SUBROUTINE calc_Zvv2


pure function get_mu_from_b_tab(b_val, b_table, mu_table, n)
    real(8), intent(in) :: b_val, b_table(n), mu_table(n)
    integer, intent(in) :: n
    real(8) :: get_mu_from_b_tab
    integer :: i

    ! Simple linear interpolation using the B-mu table
    if (b_val <= b_table(1)) then
        get_mu_from_b_tab =  mu_table(1)
    else if (b_val >= b_table(n)) then
        get_mu_from_b_tab =  mu_table(n)
    else
        do i = 1, n-1
            if (b_val >= b_table(i) .and. b_val < b_table(i+1)) then
                get_mu_from_b_tab = mu_table(i) + (mu_table(i+1)-mu_table(i)) * &
                                   (b_val-b_table(i)) / (b_table(i+1)-b_table(i))
                exit
            end if
        end do
    end if
end function get_mu_from_b_tab


pure function get_mu_from_bh_tab(b_val, b_table, h_table, n)
    real(8), intent(in) :: b_val, b_table(n), h_table(n)
    integer, intent(in) :: n
    real(8) :: get_mu_from_bh_tab,  h_val
    integer :: i

    ! Simple linear interpolation using the B-H table
    if (b_val <= b_table(1)) then
        get_mu_from_bh_tab = b_table(1) / (h_table(1) * 0.12566370964050292d-5 )
    else if (b_val >= b_table(n)) then
        get_mu_from_bh_tab = b_table(n) / (h_table(n) * 0.12566370964050292d-5 )
    else
        do i = 1, n-1
            if (b_val >= b_table(i) .and. b_val < b_table(i+1)) then
                h_val = h_table(i) + (h_table(i+1)-h_table(i)) * &
                        (b_val-b_table(i)) / (b_table(i+1)-b_table(i))
                get_mu_from_bh_tab = b_val / (h_val * 0.12566370964050292d-5 )
                exit
            end if
        end do
    end if
end function get_mu_from_bh_tab
!--------------------------------------------

SUBROUTINE new_ijk  (i_old,j_old,k_old, d_i, d_j, d_k,  movestop,  m_new, i_new,j_new,k_new )
! for the movement of a ferromagnet

    implicit none
    integer, intent(in)  ::  i_old,j_old,k_old, d_i, d_j, d_k
    integer, intent(out) :: movestop(3), m_new, i_new, j_new, k_new
    
    k_new = k_old + d_k  !fun_nod(n)%length(3)
    IF  ( k_new > sdz-2  ) THEN 
        movestop(3) =0; k_new = sdz-2
    ELSEIF(k_new < 2    ) THEN
        movestop(3) =0;  k_new = 2
    ELSEIF ( movestop(3) == 0 .and. (k_new < sdz-2 .or. k_new > 2)  ) THEN
             movestop(3) = 1 
    ENDIF
  
    j_new = j_old + d_j !fun_nod(n)%length(2)
!-------------------------------------------------------------------------------------------
! comment/uncomment to check for out of bounds along the y-axis
!  1 variant
    IF  ( j_new > sdy-2  ) THEN 
        movestop(2) =0; j_new = sdy-2
    ELSEIF(j_new < 2    ) THEN
        movestop(2) =0;  j_new = 2
    ELSEIF ( movestop(2) == 0 .and. (j_new < sdy-2 .or. j_new > 2)  ) THEN
             movestop(2) = 1 
    ENDIF
!======================================================
!  2 variant
    ! IF  ( j_new > sdy  ) THEN 
        ! movestop(2) =0; jnew = sdy
    ! ELSEIF(j_new < 0    ) THEN
        ! movestop(2) =0;  jnew = 1
    ! ELSEIF ( movestop(2) == 0 .and. (j_new < sdy .or. j_new > 0 )  ) THEN
             ! movestop(2) = 1 
    ! ENDIF
!-----------------------------------------------------------------------------------------
    ! i_new = i_old + d_i  !fun_nod(n)%length(1)
    ! ! i_new = d_i

    ! IF  ( i_new > sdx-2  ) THEN  !
        ! movestop(1) =0;  i_new = sdx-2
    ! ELSEIF(i_new < 2    ) THEN !
        ! movestop(1) =0; i_new = 2
    ! ELSEIF ( movestop(1) == 0 .and. (i_new < sdx-2 .or. i_new > 2)  ) THEN
             ! movestop(1) = 1 
    ! ENDIF

    i_new = i_old + d_i  !fun_nod(n)%length(1)
    IF  ( i_new > sdx  ) THEN  !
        movestop(1) =0;  i_new = sdx
    ELSEIF(i_new < 0    ) THEN !
        movestop(1) =0; i_new = 1
    ELSEIF ( movestop(1) == 0 .and. (i_new < sdx .or. i_new > 0)  ) THEN
             movestop(1) = 1 
    ENDIF
    
    m_new  = i_new + sdx*( (j_new-1) + sdy*(k_new-1) )
    
END SUBROUTINE new_ijk

SUBROUTINE new_m (m)
! for the movement of a coils
IMPLICIT none
integer, intent(inout) :: m
integer k, knew, j, jnew, i, inew, nij

    k = ceiling(REAL(m)/( REAL(sdx*sdy) ) ) 
    knew = k + fun_nod(n)%length(3)
    
!    1 var
    ! IF  ( knew > sdz-2  ) THEN 
        ! movestop(3) =0; knew = sdz-2
    ! ELSEIF(knew < 2    ) THEN
        ! movestop(3) =0;  knew = 2
    ! ELSEIF ( movestop(3) == 0 .and. (knew < sdz-2 .or. knew > 2)  ) THEN
             ! movestop(3) = 1 
    ! ENDIF
!-----------------
!    2 var
    IF  ( knew > sdz  ) THEN 
        movestop(3) =0; 
        knew = sdz
    ELSEIF(knew < 0    ) THEN
        movestop(3) =0;  knew = 1
    ELSEIF ( movestop(3) == 0 .and. (knew < sdz .or. knew > 0)  ) THEN
        movestop(3) = 1 
    ENDIF
!--------------

    IF (k == 1) THEN
        nij = m
    ELSE
        nij = m - (k-1)*sdx*sdy
    ENDIF
    j = ceiling( REAL(nij) / REAL(sdx) )  
    jnew = j + fun_nod(n)%length(2)
!-------------------------------------------------------------------------------------------
! comment/uncomment to check for out of bounds along the y-axis
!  1 variant
    ! IF  ( jnew > sdy-2  ) THEN 
        ! movestop(2) =0; jnew = sdy-2
    ! ELSEIF(jnew < 2    ) THEN
        ! movestop(2) =0;  jnew = 2
    ! ELSEIF ( movestop(2) == 0 .and. (jnew < sdy-2 .or. jnew > 2)  ) THEN
             ! movestop(2) = 1 
    ! ENDIF
!======================================================
!  2 variant
    IF  ( jnew > sdy  ) THEN 
        movestop(2) =0; jnew = sdy
    ELSEIF(jnew < 0    ) THEN
        movestop(2) =0;  jnew = 1
    ELSEIF ( movestop(2) == 0 .and. (jnew < sdy .or. jnew > 0 )  ) THEN
             movestop(2) = 1 
    ENDIF
!-----------------------------------------------------------------------------------------

    i = nij - (j - 1) * sdx  
    inew = i + fun_nod(n)%length(1)
!=====================================
    ! IF  ( inew > sdx-2  ) THEN  !
        ! movestop(1) =0;  inew = sdx-2
    ! ELSEIF(inew < 2    ) THEN !
        ! movestop(1) =0; inew = 2
    ! ELSEIF ( movestop(1) == 0 .and. (inew < sdx-2 .or. inew > 2)  ) THEN
             ! movestop(1) = 1 
    ! ENDIF
!=======================
!  2 variant
    IF  ( inew > sdx  ) THEN  !
        movestop(1) =0;  inew = sdx
    ELSEIF(inew < 0    ) THEN !
        movestop(1) =0; inew = 1
    ELSEIF ( movestop(1) == 0 .and. (inew < sdx .or. inew > 0)  ) THEN
        movestop(1) = 1 
    ENDIF
!====================
    m  = inew + sdx*( (jnew-1) + sdy*(knew-1) )
    
END SUBROUTINE new_m

end program ECM_MS


subroutine sor_sparse(n, values, col_ind, row_ptr, b, x, omega, max_iter, tol, iter)
! Fortran code created by J.Sochor   ( https://github.com/JNSresearcher )
! e-mail: JNSresearcher@gmail.com
! Acknowledgment: Developed with architectural assistance from Google AI (Gemini)
! June 2026

implicit none
integer, intent(in) :: n, max_iter
real(8), intent(in) :: values(*), b(n), omega, tol
integer, intent(in) :: col_ind(*), row_ptr(n+1)
integer, intent(out) :: iter
real(8), intent(inout) :: x(n)

integer :: i, j, k
real(8) :: old_x(n), diag, old_xi, diff, max_diff, summ  

do iter = 1, max_iter
    max_diff = 0.0d0
    old_x = x
    do i = 1, n
        summ = 0.d0
        diag = 0.0d0
        old_xi = x(i)
        ! Loop through non-zero elements of row i
        do k = row_ptr(i), row_ptr(i+1)-1
            j = col_ind(k)
            if (j /= i) then
                summ = summ + values(k) * x(j)
            else
                diag = values(k)
            end if
        end do
        x(i) = (1.0d0 - omega) * x(i) + (omega / diag) * (b(i) - summ)
        diff = abs(x(i) - old_xi)
        if (diff > max_diff) max_diff = diff
    end do

    if (mod(iter,10) == 0) then
        diff = norm2(x-old_x)/norm2(x)
        if (diff < tol  ) then
            exit
        endif
    endif
end do

end subroutine sor_sparse

subroutine full_sort(a,b,n,k,m)
implicit none
integer :: n,k,m,a(n,k),i,l
REAL(8)::b(n)
logical IsSorted
do
    IsSorted=.true.
    do i=1,n-1
        do l=1,m!k
            if (a(i,l)/=a(i+1,l)) exit
        end do
        if (l>m) cycle
        if (a(i,l)>a(i+1,l)) then
            call swap_lines(a(i,:),a(i+1,:),b(i),b(i+1),k)
            IsSorted=.false.
        end if
    end do
    if (IsSorted) exit
end do

end subroutine full_sort

subroutine swap_lines(a1,a2,b1,b2,k)
! swap two arrays
integer a1(k),a2(k),k,i,z
REAL(8) b1,b2,x
do i=1,k
z=a1(i); a1(i)=a2(i); a2(i)=z
end do
x=b1; b1=b2; b2=x
end subroutine swap_lines


subroutine spgemm_symbolic(n_rows_a, n_cols_a,   n_cols_b, row_ptr_a, col_ind_a, &
                           row_ptr_b, col_ind_b, row_ptr_c, total_nnz_c)
! Fortran code created by J.Sochor   ( https://github.com/JNSresearcher )
! e-mail: JNSresearcher@gmail.com
! Acknowledgment: Developed with architectural assistance from Google AI (Gemini)
! June 2026

implicit none
! Input parameters (sizes and structures of A and B)
integer, intent(in) :: n_rows_a, n_cols_a, n_cols_b
integer, intent(in) :: row_ptr_a(*), col_ind_a(*)
integer, intent(in) :: row_ptr_b(*), col_ind_b(*)

! Output parameters
integer, intent(out) :: row_ptr_c(n_rows_a + 1)
integer, intent(out) :: total_nnz_c

! Auxiliary variables
integer :: i, j, k, col_a, col_b
integer :: row_start_b, row_end_b
integer :: marker(n_cols_b) ! Stores the current line number i
integer :: row_nnz
    
    marker = 0
    total_nnz_c = 0
    row_ptr_c(1) = 1
    
    ! Iterate through each row of matrix A
    do i = 1, n_rows_a
        row_nnz = 0
        ! For each nonzero element in row A
        do k = row_ptr_a(i), row_ptr_a(i+1) - 1
            col_a = col_ind_a(k)
            ! We take the corresponding row in matrix B
            row_start_b = row_ptr_b(col_a)
            row_end_b   = row_ptr_b(col_a + 1) - 1
            do j = row_start_b, row_end_b
                col_b = col_ind_b(j)
                ! If we haven't encountered this column in the current row C yet
                if (marker(col_b) /= i) then
                    marker(col_b) = i   ! We mark that the column is occupied
                    row_nnz = row_nnz + 1
                end if
            end do
        end do
        ! Filling the string pointer for C
        row_ptr_c(i+1) = row_ptr_c(i) + row_nnz
    end do
    
    ! The total number of elements in the final matrix
    total_nnz_c = row_ptr_c(n_rows_a + 1) - 1
    
end subroutine spgemm_symbolic


subroutine csr_transpose(m, n, a, ja, ia, at, jat, iat)
! Fortran code created by J.Sochor   ( https://github.com/JNSresearcher )
! e-mail: JNSresearcher@gmail.com
! Acknowledgment: Developed with architectural assistance from Google AI (Gemini)
! June 2026

! Input parameters:
! m - number of rows of the original matrix
! n - number of columns of the original matrix
! a - array of nonzero elements of the original matrix (CSR)
! ja - array of column indices for elements of a (CSR)
! ia - array of pointers to the beginning of rows (CSR)
!
! Output parameters:
! at - array of nonzero elements of the transposed matrix (CSR)
! jat - array of row indices for elements of at (CSR)
! iat - array of pointers to the start of columns (CSR)

implicit none
integer, intent(in) :: m, n
real(8), intent(in) :: a(*)
integer, intent(in) :: ja(*), ia(m+1)
real(8), intent(out) :: at(*)
integer, intent(out) :: jat(*), iat(n+1)

integer :: i, j, k, col, p, nz
integer, allocatable :: work(:)
    ! 1. Determine the number of non-zero elements in each column
    allocate(work(n))
    work = 0
    iat=0
    nz = ia(m+1) - 1 ! Total number of non-zero elements
    do k = 1, nz
        col = ja(k)
        work(col) = work(col) + 1
    end do
    ! 2. Let's construct the array iat (analog ia for the transposed matrix)
    iat(1) = 1
    do j = 1, n
        iat(j+1) = iat(j) + work(j)
    end do
    ! 3. Redistribute the elements!
    ! First, copy work to iat to temporarily store the offsets.
    work = iat(1:n)
    do i = 1, m
        do p = ia(i), ia(i+1)-1
            j = ja(p)
            k = work(j)
            at(k) = a(p)
            jat(k) = i  ! Row indices become column indices
            work(j) = work(j) + 1
        end do
    end do
    
    deallocate(work)
end subroutine csr_transpose


subroutine  sp_mm (row_a, col_a, col_b, irow_a,jcol_a,val_a, irow_b,jcol_b,val_b,  irow_c,jcol_c,val_c, curr_c_ptr )
! Matrix multiplication in CSR format

! Fortran code created by J.Sochor   ( https://github.com/JNSresearcher )
! e-mail: JNSresearcher@gmail.com
! Acknowledgment: Developed with architectural assistance from Google AI (Gemini)
! June 2026

implicit none
integer :: row_a
integer :: irow_a(row_a + 1)
integer :: jcol_a(*)
real(8) :: val_a(*)

integer :: col_a
integer :: irow_b(col_a + 1)
integer :: jcol_b(*)
real(8) :: val_b(*)

integer :: irow_c(row_a + 1)
integer :: jcol_c(*) 
real(8) :: val_c(*)

! Auxiliary structures for storing packed arrays (SPA)
integer :: col_b
real(8) :: tmp_val(col_b)
integer :: tmp_occ(col_b)            ! Employment mask
integer :: tmp_idx(col_b), n_idx     ! List of indexes involved

integer :: i, j, k, row_start, row_end, colon_b, curr_c_ptr
real(8) :: a_val
    
    tmp_val   = 0.0
    tmp_occ   = 0  
    curr_c_ptr = 1
    irow_c(1) = 1
    ! The main loop over the rows of matrix A
    do i = 1, row_a
        n_idx = 0
        ! Iterating over non-zero elements of string A
        do k = irow_a(i), irow_a(i+1) - 1
            a_val = val_a(k)
            row_start = irow_b(jcol_a(k))
            row_end   = irow_b(jcol_a(k) + 1) - 1
            ! Traverse the corresponding row of B (linear combination of rows)
            do j = row_start, row_end
                colon_b = jcol_b(j)
                if (tmp_occ(colon_b) /= i) then
                    tmp_occ(colon_b) = i
                    n_idx = n_idx + 1
                    tmp_idx(n_idx) = colon_b
                    tmp_val(colon_b) = a_val * val_b(j)
                else
                    tmp_val(colon_b) = tmp_val(colon_b) + a_val * val_b(j)
                end if
            end do
        end do
        
        ! Transferring accumulated data from SPA to the C structure
        ! For the CSR to be correct, it is advisable to sort the indexes in the row.
        do k = 1, n_idx
            colon_b = tmp_idx(k)
            val_c(curr_c_ptr) = tmp_val(colon_b)
            jcol_c(curr_c_ptr) = colon_b
            curr_c_ptr = curr_c_ptr + 1
            tmp_val(colon_b) = 0.0 ! Cleaning up for the next iteration
        end do
        irow_c(i+1) = curr_c_ptr
    end do

endsubroutine  sp_mm 


SUBROUTINE writeVtk_src ( Npoint, numfun,  NcellsX, NcellsY, NcellsZ, new_nodesX, new_nodesY,  new_nodesZ,  &
                            sdx, sdy,  delta, files)
! Fortran code created by J.Sochor   ( https://github.com/JNSresearcher )
! e-mail: JNSresearcher@gmail.com

! This program generates a .vtk file for sources (vector, moving or stationary)

USE m_vxc2data
IMPLICIT NONE

INTEGER,      INTENT(IN):: Npoint, numfun,  sdx,sdy
INTEGER,      INTENT(IN):: NcellsX, NcellsY, NcellsZ, &
                            new_nodesX(NcellsX),new_nodesY(NcellsY),new_nodesZ(NcellsZ) !(nsub)
REAL(8),      INTENT(IN):: delta(3)
CHARACTER(*), INTENT(IN):: files

CHARACTER (LEN = 1)  ci
CHARACTER (LEN = 4)  ch_sd
CHARACTER (LEN = 30) buf1, buf2
INTEGER:: i,j,k, l, m, n,numcells, nnX, nnY, nnZ, nij,  ios, n_new !nim,njm,nkm,nip,njp,nkp,
REAL(8):: s,sxm,sym,szm

numcells =0
DO n=1,numfun
    numcells = numcells + fun_nod(n)%numnod_Fx 
enddo
WRITE(ch_sd,'(i4)')Npoint
ios=0

OPEN(newunit=n_new,file=trim(files)//'/src_'//trim(adjustl(ch_sd))//'.vtk', form="UNFORMATTED",&
access="STREAM",  convert="big_endian", iostat=ios)  
! BUFFERED='YES', BLOCKSIZE=3270000, 
IF (ios/=0) THEN
    PRINT*,'error! Could not open the outfile result.';   STOP
ENDIF

! header output
WRITE (n_new   ) "# vtk DataFile Version 3.0"//new_line(ci)//"out data result"//new_line(ci)//"BINARY" //new_line(ci)
WRITE(n_new  )"DATASET UNSTRUCTURED_GRID"//new_line(ci)
WRITE(buf1,'(i8)')   numcells*8  !
WRITE(n_new  ) "POINTS "//trim(adjustl(buf1))//" double" //new_line(ci)

! output of coordinates
nnX = 0; nnY = 0;nnZ = 0;   ! vector source node counter
DO n=1,numfun
    DO l = 1, fun_nod(n)%numnod_Fx
        nnX = nnX + 1
        m = new_nodesX(nnX)
        CALL find_coord
        CALL write_coord
    ENDDO
ENDDO
WRITE(n_new)  new_line(ci)

WRITE(buf1,'(i8)')   numcells 
WRITE(buf2,'(i8)') 9*numcells  
WRITE(n_new) "CELLS "//trim(adjustl(buf1))//" "//trim(adjustl(buf2))//new_line(ci)
DO i=0, numcells - 1   !( NcellsX + NcellsY + NcellsZ)-1
    WRITE(n_new) 8, 8*i + [0, 1, 2, 3, 4, 5, 6, 7]
ENDDO
WRITE(n_new)   new_line(ci)

WRITE(n_new ) "CELL_TYPES "//trim(adjustl(buf1))//new_line(ci)
DO i=1,numcells  !( NcellsX + NcellsY + NcellsZ)
    WRITE(n_new) 11
ENDDO
WRITE(n_new)  new_line(ci)

WRITE(n_new) "CELL_DATA "//trim(adjustl(buf1))//new_line(ci)
WRITE(n_new) "VECTORS "//"Vector_field_SRC"//" double"//new_line(ci)
DO n=1,numfun
    s = Fun(n)%vely
    do k=1,fun_nod(n)%numnod_Fx
        WRITE(n_new) s*fun_nod(n)%nods_Fx_cos(k), s*fun_nod(n)%nods_Fy_cos(k), s*fun_nod(n)%nods_Fz_cos(k)
    enddo
ENDDO

CLOSE(n_new)

CONTAINS
    SUBROUTINE find_coord
        k = ceiling(REAL(m)/( REAL(sdx*sdy) ) ) 
        IF (k == 1) THEN
            nij = m
        ELSE
            nij = m - (k-1)*sdx*sdy
        ENDIF
        j = ceiling( REAL(nij) / REAL(sdx) ) 
        i = nij - (j - 1) * sdx  
    END SUBROUTINE find_coord

    SUBROUTINE write_coord
        ! point 0 x = i; y=j ; z=k
        sxm = REAL(i,8)*delta(1) + delta(1)/2.d0
        sym = REAL(j,8)*delta(2) + delta(2)/2.d0 
        szm = REAL(k,8)*delta(3) + delta(3)/2.d0 
        WRITE(n_new ) sxm, sym, szm
        ! point 1 x = i+1; y=j ; z=k
        sxm = REAL(i+1,8)*delta(1) + delta(1)/2.d0
        sym = REAL(j,8)  *delta(2) + delta(2)/2.d0
        szm = REAL(k,8)  *delta(3) + delta(3)/2.d0
        WRITE(n_new) sxm, sym, szm
        ! point 2 x = i; y=j+1 ; z=k
        sxm = REAL(i,8)  *delta(1) + delta(1)/2.d0
        sym = REAL(j+1,8)*delta(2) + delta(2)/2.d0 
        szm = REAL(k,8)  *delta(3) + delta(3)/2.d0
        WRITE(n_new) sxm, sym, szm
        ! point 3 x = i+1; y=j+1 ; z=k
        sxm = REAL(i+1,8)*delta(1) + delta(1)/2.d0
        sym = REAL(j+1,8)*delta(2) + delta(2)/2.d0
        szm = REAL(k,8)  *delta(3) + delta(3)/2.d0
        WRITE(n_new) sxm, sym, szm
        ! point 4 x = i; y=j ; z=k+1
        sxm = REAL(i,8)*delta(1)   + delta(1)/2.d0
        sym = REAL(j,8)*delta(2)   + delta(2)/2.d0
        szm = REAL(k+1,8)*delta(3) + delta(3)/2.d0
        WRITE(n_new ) sxm, sym, szm
        ! point 5 x = i+1; y=j ; z=k+1
        sxm = REAL(i+1,8)*delta(1) + delta(1)/2.d0
        sym = REAL(j,8)  *delta(2) + delta(2)/2.d0
        szm = REAL(k+1,8)*delta(3) + delta(3)/2.d0
        WRITE(n_new ) sxm, sym, szm
        ! point 6 x = i; y=j+1 ; z=k+1
        sxm = REAL(i,8)  *delta(1) + delta(1)/2.d0
        sym = REAL(j+1,8)*delta(2) + delta(2)/2.d0
        szm = REAL(k+1,8)*delta(3) + delta(3)/2.d0
        WRITE(n_new) sxm, sym, szm
        ! point 7 x = i+1; y=j+1 ; z=k+1
        sxm = REAL(i+1,8)*delta(1) + delta(1)/2.d0
        sym = REAL(j+1,8)*delta(2) + delta(2)/2.d0
        szm = REAL(k+1,8)*delta(3) + delta(3)/2.d0
        WRITE(n_new) sxm, sym, szm
    END SUBROUTINE write_coord
END SUBROUTINE writeVtk_src

subroutine writeVtk (Npoint, sdx, sdy, sdz, num_brX, num_brY, num_kx, num_ky, delta, iv, files) !ev, 
! Fortran code created by J.Sochor   ( https://github.com/JNSresearcher )
! e-mail: JNSresearcher@gmail.com

implicit none

integer,      intent(IN):: Npoint, sdx,sdy,sdz, num_brX, num_brY, num_kx, num_ky
real(8),      intent(IN):: iv(*), delta(3)
character(*), intent(IN):: files

character (len = 1)  ci
character (len = 4)  ch_sd
character (len = 30) buf1
integer:: i,j,k, nfx, sdx2,sdy2,sdz2,Nodes, ios, n_new
real(8):: sxm,sym,szm
integer:: ip, jp, kp, Cells
    sdx2=sdx+1; sdy2=sdy+1; sdz2=sdz+1;
    Nodes = sdx2*sdy2*sdz2
    Cells = 3*sdx*sdy*sdz
    write(ch_sd,'(i4)')Npoint
    ios=0
    
    open(newunit=n_new,file=trim(files)//'/field_'//trim(adjustl(ch_sd))//'.vtk',form="UNFORMATTED",  &
       access="STREAM",  convert="big_endian", iostat=ios)  

    if (ios/=0) then
        print *,'fatal error! Could not open the outfile result.';   stop
    end if
    
    write (n_new) "# vtk DataFile Version 3.0"//new_line(ci)//"out data result"//new_line(ci)//"BINARY"//new_line(ci)
    write(buf1,'(i8," ",i8," ",i8)') sdx2, sdy2, sdz2  
    write(n_new)"DATASET STRUCTURED_GRID"//new_line(ci)//"DIMENSIONS "//trim(adjustl(buf1))//new_line(ci)
    write(buf1,'(i8)')   (sdx2) * (sdy2) * (sdz2)
    write(n_new) "POINTS "//trim(adjustl(buf1))//" float"//new_line(ci)

    do k=1,sdz2
        szm = real(k,8)*delta(3)  + delta(3) /2.d0
        do j=1,sdy2
            sym = real(j,8)*delta(2) + delta(2) /2.d0
            do i=1,sdx2
                sxm = real(i,8)*delta(1) + delta(1) /2.d0
                write(n_new) real(sxm,4), real(sym,4), real(szm,4)
            enddo
        enddo
    enddo
    write(n_new)  new_line(ci)

    write(n_new) "POINT_DATA "//trim(adjustl(buf1))//new_line(ci)
    
    write(n_new) "VECTORS "//"Vector_field_B"//" float"//new_line(ci)
    ip=0; jp=0; kp=0
    nfx=num_brX + num_brY
    do k = 1,sdz2;  do j = 1,sdy2;  do i = 1,sdx2
        sxm=0.d0;sym=0.d0;szm=0.d0;
        if ( i < sdx2  ) then
            ip = ip + 1
            sxm =   iv(ip)
        else
            sxm = 0.d0
        endif 
        if ( j < sdy2) then
            jp = jp + 1
            sym =   iv(jp + num_brX )
        else
            sym = 0.d0
        endif 
        if ( k < sdz2) then
            kp = kp + 1
            szm =   iv(kp + num_brX + num_brY )
        else
            szm = 0.d0
        endif 
        write(n_new) real( sxm,4 ), real( sym,4  ), real(szm,4)
    end do;  end do;  end do
    write(n_new)  new_line(ci)

    close(n_new)

end subroutine writeVtk

