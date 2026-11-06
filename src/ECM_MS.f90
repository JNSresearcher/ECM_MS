program ECM_MS
! Fortran codes created by J.Sochor   ( https://github.com/JNSresearcher )
! e-mail: JNSresearcher@gmail.com
! November, 2026

USE m_fparser
USE m_vxc2data

use omp_lib
IMPLICIT NONE

interface
    subroutine permute_csr_split_diag(n, ia, ja, a, colors, num_colors, &
                                        ia_new, ja_new, a_new, diag, old_to_new, color_ptr, ja_diag,ja_nodiag ) 
        integer, intent(in) :: n
        integer, intent(in) :: ia(n+1), ja(*)
        real(8), intent(in) :: a(*)
        integer, intent(in) :: colors(n)
        integer, intent(in) :: num_colors
        
        integer, intent(out) :: ia_new(n+1)              ,ja_diag(n)
        integer, intent(out), allocatable :: ja_new(:)   ,ja_nodiag(:)
        real(8), intent(out), allocatable :: a_new(:)
        real(8), intent(out) :: diag(n)
        integer, intent(out) :: old_to_new(n)
        integer, intent(out) :: color_ptr(num_colors + 1)
    end subroutine permute_csr_split_diag
endinterface

type espm
    integer im
    REAL(8) em
    type (espm),pointer ::prec
end type espm
type(espm),pointer :: sp_col, next_col 

!---------------------------------!
! data from SUBROUTINE vxc2data   !
!---------------------------------!
REAL(8)::   delta(3), &      ! array grid spacing along X, Y and Z   
            dt,       &      ! time step
            Time,     &      ! stop time
            dtt              ! jump duration

INTEGER     sdx,sdy,sdz, &   ! number of cells along X,  Y and  Z   
            nsub,        &   ! number of physical domains
            nsub_air,    &   ! number of environment domains
            nsub_glob,   &   ! total number of domains
            numfun,      &   ! number of functions for calculating external sources
            numMech,     &   ! number of functions for calculating movements of external sources
            numVenv,     &   ! number of functions for calculating movements of ferromagnets
            mode            ! 
!//=======================declare SOR and nonlin-mag
REAL(8)    tol, omega, tolmag, omegamag
INTEGER    iter, itmax, itmag
            
! solv - character string corresponding to the name of methods: only 'SOR'
CHARACTER(LEN=4) :: solv 
CHARACTER(32)    :: files      ! name for output files 

! Ek - vector of sources of contour emf, Ik is vector of contour currents, Iv is the vector of branch currents.
real(8), allocatable::  Ek(:), Ik(:),  Iv(:)

integer, allocatable :: colors(:), old_to_new(:), color_ptr(:)
integer :: num_colors
integer, allocatable :: ia_new(:), ja_new(:), ja_diag(:),ja_nodiag(:)
real(8), allocatable :: a_new(:), diag(:)
real(8), allocatable :: b_perm(:), x_perm(:)
real(8) :: omega_opt, res_norm 

integer num_k,   &    ! total number of contours
        num_kX,  &    ! number of contours along the axis X 
        num_kY,  &    ! number of contours along the axis Y 
        num_kZ,  &    ! number of contours along the axis Z 
        num_brX, &    ! number of branches along the axis X
        num_brY, &    ! number of branches along the axis Y
        num_brZ, &    ! number of branches along the axis Z
        num_br,  &    ! total number of branches
        Cells         !  number of cells
! arrays for converting cell numbers into contour numbers and branch numbers
INTEGER, allocatable::  cells2cont(:), cells2br(:,:) 
! auxiliary variables
real(8) a,s, Sx, Sy, Sz, MU_new,MU_sum, Z0, SmuX, SmuY, SmuZ
INTEGER  i,j,k, L, m, n, nn,ios, mm, iterac 
integer:: nod(1), kdz, ip,jp,kp, nc, i_p,j_p,k_p
! indicators of movement of ferromagnets and windings
integer movestop_Fe(3), flag_move_Fe
integer d_i,d_j,d_k, flag_move, movestop(3), m_new, i_new, j_new, k_new, NcellsX, NcellsY, NcellsZ
! arrays for new cells of field sources as they move
INTEGER, ALLOCATABLE ::  new_nodesX(:), new_nodesY(:), new_nodesZ(:) 

! output control
integer Nout_display,Nprint_display, Nout, Npoint, Nprint, Ntime
integer k0,k1,k2
REAL(8)   T, Tcalc, Tsavedata

! auxiliary arrays
integer col_kv(4)
real(8) val_kv(4)

! topological matrices Ckv, Cvk, CZkv, and contour matrix Zkk in CSR format
real(8), allocatable :: val_Ckv(:),  val_Cvk(:),  val_CZkv(:),   val_Zkk(:),  Zvv(:) 
integer, allocatable :: irow_Ckv(:), irow_Cvk(:), irow_CZkv(:) , irow_Zkk(:)
integer, allocatable :: jcol_Ckv(:), jcol_Cvk(:), jcol_CZkv(:),  jcol_Zkk(:), ipvt(:)
integer nz_Ckv, nz_Cvk,  nz_CZkv, nz_Zkk  ! number of nonzero elements

INTEGER(1) n_Fe
INTEGER(1),ALLOCATABLE:: geoPHYS_Fe(:,:,:) ! map of the location of ferromagnetic domains
integer(1),pointer:: temp_geoPHIS(:)       ! for replacing winding regions with air
INTEGER,ALLOCATABLE:: indexes_geo(:)

type t_dom_Fe  
    INTEGER, ALLOCATABLE:: num_brX(:), num_brY(:),num_brZ(:),  & ! branch numbers
                            cells(:),  &                         ! cell numbers in the center and on the six sides
                            cellsXp(:,:), cellsYp(:,:),cellsZp(:,:), &
                            cellsXm(:,:),   cellsYm(:,:),  cellsZm(:,:)
    REAL(8) mu_mean, mu_mean_old                                ! mean magnetic permeabilities for iterations
    INTEGER n_Fe, num_cells                                    ! domain number and number of cells in a domain
    REAL(8), ALLOCATABLE:: mu_cells(:), B_dat(:), H_dat(:), MU_dat(:)  ! tabular data B-H B-mu
    integer n_dat, nonlin, typ_nonlin
end type t_dom_Fe
TYPE (t_dom_Fe), ALLOCATABLE :: dom_Fe(:)                   ! sizdom_Fe transmitted through the module
!========for iterac 
real(8)  summ_mu_mean_old, summ_mu_mean, gauge
integer  local_iter  
!================================END DECLARE====================

!----- START TIMER
call system_clock(k1,k0);
            !-----------------------
            ! input of initial data
            ! ----------------------
CALL vxc2data ( delta,    dt,        Time,     dtt,                   &
                sdx,    sdy,      sdz,       nsub,     nsub_air, nsub_glob, &
                numfun, numMech,  numVenv,     &
                solv,   files,    tol, omega, itmax, tolmag, itmag, omegamag, mode, gauge )  

num_KX=sdx*(sdy-1)*(sdz-1);      ! number of contours along the axis X 
num_kY=(sdx-1)*sdy*(sdz-1);      ! number of contours along the axis Y 
num_kZ=(sdx-1)*(sdy-1)*(sdz);    ! number of contours along the axis Z 

num_brX = (sdx-1)*sdy*sdz;      ! number of branches along the axis X
num_brY = sdx*(sdy-1)*sdz;      ! number of branches along the axis Y
num_brZ = sdx*sdy*(sdz-1);      ! number of branches along the axis Z

Cells = sdx*sdy*sdz                     ! number of cells
kdz = sdx*sdy                           ! auxiliary variable

num_k = num_kX + num_kY + num_kZ        ! total number of contours
print '( a, *(a,i7))','number of contours along:', ' X-axis=', num_kx,' Y-axis=',num_ky,' Z-axis=', num_kz, ' total number=',num_k

num_br = num_brX + num_brY +num_brZ     ! total number of branches
allocate (Zvv(num_br), source=1.d0) ! 
allocate (Ek(num_k), Ik(num_k), Iv(num_br), source=0.d0)
print '( a, *(a,i7))','number of branches along:', ' X-axis=', &
          num_brX,' Y-axis=',num_brY,' Z-axis=', num_brZ, ' total number=',num_br, ' number of cells', Cells

allocate (cells2br(Cells,3) , source=0)
!---------------------------------------------------------
! generation of information about the numbers of branches
!---------------------------------------------------------
! converting a cell number to a branch number
! OUT: m=cells2br(n,i)  m: branch number  n: cell number; i=1(for X) i=2(for Y) ..
ip=0; jp=0; kp=0;  ! counters
do k = 1,sdz;  do j = 1,sdy;  do i = 1,sdx
    n = i + sdx * ((j -1) + sdy*(k -1))
    if ( i < sdx) then
        ip = ip + 1
        i_p = ip                      ! num br X
    else
        i_p = 0
    endif 
    if ( j < sdy) then
        jp = jp + 1
        j_p = jp + num_brX            ! num br Y
    else
        j_p = 0
    endif 
    if ( k < sdz) then
        kp = kp + 1
        k_p = kp + num_brX + num_brY  ! num br Z
    else
        k_p = 0
    endif 
    cells2br(n,1) = i_p  ! on the border sdx i_p=0
    cells2br(n,2) = j_p  ! on the border sdy j_p=0
    cells2br(n,3) = k_p  ! on the border sdz k_p=0
end do;  end do;  end do

!---------------------------------------------------------
! generation of information about ferromagnetic domains
!---------------------------------------------------------
if (sizdom_Fe /=0)  then
    ALLOCATE (geoPHYS_Fe(sdx,sdy,sdz),  source=int(nsub_glob,1)) 
    ALLOCATE ( dom_Fe(sizdom_Fe) )  
    dom_Fe%nonlin = 0                  ! nonlinear ferromagnet signal
    ! Domain numbers not belonging to the ferromagnetic material (coils) are replaced with air.
    temp_geoPHIS(1:Cells) => geoPHYS 
    indexes_geo =   pack( [(i, i=1, Cells)], [(all(temp_geoPHIS(i) /= numdom_Fe), i=1, Cells)] )
    temp_geoPHIS(indexes_geo) = int(nsub_glob,1)
    nullify(temp_geoPHIS);	deallocate(indexes_geo)
    do n = 1, sizdom_Fe 
        n_Fe = numdom_Fe(n)  ! 
        dom_Fe(n)%n_Fe = int (n_Fe, 4)
        m=0
        do k=1, sdz; do j=1, sdy; do i=1, sdx
            if (geoPHYS(i,j,k) == n_Fe) then
                geoPHYS_Fe(i,j,k) = n_Fe
                ! counting the number of cells in each domain
                m=m+1
            endif
        enddo; enddo; enddo
        ! save 
        dom_Fe(n)%num_cells = m  
        ! for each cell Fe the numbers of XYZ branches will be known
        ALLOCATE(dom_Fe(n)%num_brX(m),dom_Fe(n)%num_brY(m), dom_Fe(n)%num_brZ(m), dom_Fe(n)%cells(m), &
                 dom_Fe(n)%cellsXp(m,2),dom_Fe(n)%cellsYp(m,2),dom_Fe(n)%cellsZp(m,2),  &
                 dom_Fe(n)%cellsXm(m,2),  dom_Fe(n)%cellsYm(m,2),  dom_Fe(n)%cellsZm(m,2), source = 0)
        ALLOCATE(dom_Fe(n)%mu_cells(m), source = 0.d0)
        L = int(n_Fe, 4)
        dom_Fe(n)%mu_mean = valPHYS(L,1)
    enddo
    if (siz_non_lin_MU /= 0 ) then  
    ! siz_non_lin_MU, etc.  - the initial information is generated in the vxc2data subroutine and
    ! passed through the m_vxc2data module.
        do n = 1, siz_non_lin_MU
            do m=1, sizdom_Fe
                if (dom_Fe(m)%n_Fe == non_lin_MU(n)%num_domain) then
                    dom_Fe(m)%nonlin = 1
                    if (non_lin_MU(n)%typ_non_lin == 'BH') then
                        dom_Fe(m)%typ_nonlin = 2
                    ELSEIF (non_lin_MU(n)%typ_non_lin == 'BM') then
                        dom_Fe(m)%typ_nonlin = 1
                    else
                        print*, 'error type nonlin'
                        stop
                    endif
                    open (1, file = non_lin_MU(n)%file_MU )
                    dom_Fe(m)%B_dat = [real(8):: ]
                    if (dom_Fe(m)%typ_nonlin == 2) then
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
                            if (dom_Fe(m)%typ_nonlin == 2) dom_Fe(m)%H_dat = [dom_Fe(m)%H_dat, Sy]
                            if (dom_Fe(m)%typ_nonlin == 1) dom_Fe(m)%MU_dat = [dom_Fe(m)%MU_dat, Sy]
                        endif
                    END DO
                    dom_Fe(m)%n_dat = size(dom_Fe(m)%B_dat)
                    close(1)
                    print '(a, a, a, i3, a,i3)',  'create tabl nonlin mu from file ', non_lin_MU(n)%file_MU, &
                                                  ' dom=', m, ' siz=', dom_Fe(m)%n_dat
                    if (dom_Fe(m)%n_dat == 0) then
                        print*, 'ERROR read file ', non_lin_MU(n)%file_MU
                        stop
                    endif
                    dom_Fe(m)%mu_mean_old = valPHYS(dom_Fe(m)%n_Fe,1)
                endif
            enddo
        enddo    
    endif  ! siz_non_lin_MU /= 0
    ! converting a local cell numbers to a global cell numbers and 
    ! branch numbers for ferromagnetic domains
    do nn = 1, sizdom_Fe
        dom_Fe(nn)%cells = 0
        L=0 ! domain number counter Fe
        do k = 1,sdz;  do j = 1,sdy;  do i = 1,sdx
            m = i + sdx * ((j -1) + sdy*(k -1))    ! global cell number
            if (geoPHYS (i,j,k) == int (dom_Fe(nn)%n_Fe,1)  ) then 
                L = L+1                           ! counter cells local
                dom_Fe(nn)%cells(L) = m
            endif
        end do;  end do;  end do
    enddo
    do nn = 1, sizdom_Fe
        do m = 1, dom_Fe(nn)%num_cells   ! number of cells in a domain
            ! for each m-th cell, branch numbers
            mm = dom_Fe(nn)%cells(m)             ! global cell number Fe 
            dom_Fe(nn)%num_brX(m) = cells2br(mm,1)  
            dom_Fe(nn)%num_brY(m) = cells2br(mm,2)  
            dom_Fe(nn)%num_brZ(m) = cells2br(mm,3)  
        enddo 
    enddo
endif  ! sizdom_Fe /=0
!-----------------------------------------------------------
!                FORMATION OF A CONTOUR MATRIX
!                Contour-branch matrix
!-----------------------------------------------------------
allocate (irow_Ckv(num_k+1), cells2cont(3*Cells),  source=0)
nullify(sp_col, next_col )
irow_Ckv(1)=1
nz_Ckv =0 ! non-zero element counter

! X contours in the YZ plane
nn=0
L=0
do i=1,sdx; do k=1,sdz-1; do j=1,sdy-1
    nn =  nn + 1                     ! current contour X
    n = i + sdx * ((j -1) + sdy*(k -1)) ! cell number
    cells2cont(n) = nn   !  n - cell number, nn - contour number X
    col_kv = 0; val_kv = 0.d0
    jp  = i     + sdx*(j+1-1) + sdx*sdy*(k-1) + num_brX + num_brY    ! z
    nc  = i     + sdx*(j-1)   + sdx*sdy*(k-1) + num_brX + num_brY    !-z
    kp  = i     + sdx*(j-1)   + sdx*(sdy-1)*(k+1-1) + num_brX             !-y
    ip  = i     + sdx*(j-1)   + sdx*(sdy-1)*(k-1)   + num_brX             ! y
                !     Iz    -Iz      -Iy       Iy                       !  branches
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
        nz_Ckv = nz_Ckv + 1                       ! non-zero element counter
        next_col => sp_col
    enddo
    irow_Ckv(nn+1) = irow_Ckv(nn) + 4  ! pointer to the column index for the next row
enddo; enddo; enddo

! Y contours in the XZ plane
do j=1,sdy; do k=1,sdz-1; do i=1, sdx-1
    nn =  nn + 1                ! current contour  Y
    n = i + sdx * ((j -1) + sdy*(k -1)) + Cells
    cells2cont(n) = nn           ! n - cell number, nn - contour number Y
    col_kv = 0; val_kv = 0.d0
    jp  = i+1   + sdx*(j-1) + sdx*sdy*(k-1) + num_brX + num_brY  !-z
    nc  = i     + sdx*(j-1) + sdx*sdy*(k-1) + num_brX + num_brY  ! z
    kp  = i     + (sdx-1)*(j-1)   + (sdx-1)*sdy*(k-1)                      !-x
    ip  = i     + (sdx-1)*(j-1)   + (sdx-1)*sdy*(k+1-1)                    ! x
                !   Iz       -Iz       Ix     -Ix                   !  branches
    col_kv(1:4) = [ nc,     jp ,    ip,    kp]
    val_kv(1:4) = [ 1.d0,  -1.d0,  1.d0,  -1.d0   ]
    !--------------------------Y------------------------------------
    ! sort branch numbers in ascending order
    call full_sort(col_kv, val_kv, 4, 1,1)
    do m = 1, 4
        if (col_kv(m) == 0) then
            print*, 'Error: contour number=',nn, 'm',m,'col_kv(m)',col_kv(m)
            stop
        endif
        allocate(sp_col)
        sp_col=espm(col_kv(m),val_kv(m),next_col)   ! COLUMN PACKAGING
        nz_Ckv = nz_Ckv + 1                         ! non-zero element counter
        next_col => sp_col
    enddo
    ! pointer to the column index for the next row
    irow_Ckv(nn+1) = irow_Ckv(nn) + 4 
enddo; enddo; enddo

!  Z contours in the XY plane
do k=1,sdz; do j=1,sdy-1; do i=1, sdx-1
    nn =  nn + 1                        ! current contour Z
    n = i + sdx * ((j -1) + sdy*(k -1)) + 2*Cells
    cells2cont(n) = nn   ! n - cell number, nn - contour number Z
    col_kv = 0; val_kv = 0.d0
    jp  = i    + (sdx-1)*(j+1-1)   + (sdx-1)*sdy*(k-1)             !-x
    nc  = i    + (sdx-1)*(j-1)     + (sdx-1)*sdy*(k-1)             !+x
    kp =  i    + sdx*(j-1)   + sdx*(sdy-1)*(k-1) + num_brX    !-y 
    ip  = i+1  + sdx*(j-1)   + sdx*(sdy-1)*(k-1) + num_brX    !+y 
                !    -Ix         Ix      -Iy    Iy          !  branches
    col_kv(1:4) = [  jp,       nc,     kp,    ip]
    val_kv(1:4) = [ -1.d0,    1.d0,   -1.d0,  1.d0   ]
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
i=nz_Ckv 
DO WHILE (associated(sp_col))
    jcol_Ckv(i) = sp_col%im  ! column number
    val_Ckv(i) = sp_col%em   ! value
    i=i-1
    sp_col => sp_col%prec
END DO
nullify(sp_col, next_col)
! END FORMATION OF A CONTOUR MATRIX
            !----------------------------
            ! resistance calculation
            !----------------------------
Z0= valPHYS(nsub_glob,1)    ! environmental resistance
if ( delta(1) == delta(2)  .and.   delta(1) == delta(3) ) then
    SmuX=  gauge * delta(1)/(0.12566370964050292d-5*delta(1)*delta(1))
    SmuY=SmuX
    SmuZ=SmuX
    if (sizdom_Fe == 0)  Zvv = Z0*SmuX
else 
    SmuX= gauge * delta(1)/(0.12566370964050292d-5*delta(2)*delta(3))   
    SmuY= gauge * delta(2)/(0.12566370964050292d-5*delta(1)*delta(3))
    SmuZ= gauge * delta(3)/(0.12566370964050292d-5*delta(2)*delta(1)) ! dL/(mu0 dS)
    if (sizdom_Fe == 0) then
        Zvv(1:num_brX) = Z0*SmuX
        Zvv(num_brX+1: num_brX+num_brY) = Z0*SmuY
        Zvv(num_brX+num_brY+1: num_br)  = Z0*SmuZ
    endif
endif

!  mu = (mu1 + mu2)/2  1/mu = 2/(mu1 + mu2)   2/(1/a + 1/s) 
if (sizdom_Fe /= 0) then
    do n = 1, Cells    ! In all cells
        i=mod(n-1,sdx)+1; j=mod((n-1)/sdx,sdy)+1; k=(n-1)/kdz+1
        ip=cells2br(n,1)
        jp=cells2br(n,2)
        kp=cells2br(n,3)
        m = int(geoPHYS(i,j,k),4);
        a = valPHYS(m,1)
        if (ip /=0) then
            mm = int(geoPHYS(i+1,j,k),4); 
            s = valPHYS(mm,1)
            Zvv(ip)  = 2.d0*SmuX/( 1.d0/a + 1.d0/s)   
        endif
        if (jp /=0) then
            mm = int(geoPHYS(i,j+1,k),4); 
            s = valPHYS(mm,1)
            Zvv(jp)  = 2.d0*SmuY/( 1.d0/a + 1.d0/s)
        endif
        if (kp /=0) then
            mm = int(geoPHYS(i,j,k+1),4); 
            s = valPHYS(mm,1)
            Zvv(kp)  = 2.d0*SmuZ/( 1.d0/a + 1.d0/s)
        endif
    enddo
    do n = 1, sizdom_Fe        ! Initialization of mu for all domain cells
        mm =  dom_Fe(n)%n_Fe   ! domain number
        dom_Fe(n)%mu_cells = valPHYS(mm,1)
    enddo
    call bnd_cells  ! Defining the boundary
endif
! formation of the product matrix Ckv*Zvv = CZkv
nz_CZkv  = nz_Ckv
nz_Cvk   = nz_Ckv

allocate( val_Cvk(nz_Cvk),                       val_CZkv(nz_CZkv),   source=0.d0)
allocate(jcol_Cvk(nz_Cvk), irow_Cvk(num_br+1), jcol_CZkv(nz_CZkv), irow_CZkv(num_k+1) , source=0)
irow_CZkv = irow_Ckv
jcol_CZkv = jcol_Ckv
        !--------------------
        ! CZkv = Ckv * Zvv 
        !--------------------
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
            ! Estimating non-zero values for a product CZkv*Cvk
            !---------------------------------------------------- 
allocate( irow_Zkk(num_k+1) , source=0)
call spgemm_symbolic(num_k, num_br,num_k,irow_CZkv,jcol_CZkv,irow_Cvk, jcol_Cvk,irow_Zkk, nz_Zkk)
! memory allocation for Zkk
allocate(val_Zkk(nz_Zkk), source=0.d0)
allocate(jcol_Zkk(nz_Zkk) , source=0)
            !--------------------------------------------
            ! Calculation of the contour matrix Zkk using 
            ! the product of matrices: Zkk = CZkv*Cvk
            !-------------------------------------------
call sp_mm (num_k, num_br, num_k, irow_CZkv,jcol_CZkv,val_CZkv, irow_Cvk,jcol_Cvk,val_Cvk, irow_Zkk,jcol_Zkk,val_Zkk, nz_Zkk )
nz_Zkk = nz_Zkk - 1
do  j = 1, num_k
    ip = irow_Zkk(j); jp = irow_Zkk(j+1) - 1
    m = jp - ip + 1
    n= m/2
    allocate(ipvt(1:n), source=0)
    call FLASH2 (jcol_Zkk(ip:jp), val_Zkk(ip:jp), m, ipvt, n)
    deallocate(ipvt)
enddo
PRINT '(a, i9, a,g12.5,a)', 'sparse contour resistance matrix is ready! non zero=',  &
                        nz_Zkk, ' density=', 100.d0* REAL(nz_Zkk)/ (REAL(num_k)* REAL(num_k)),'%'
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
    DO k = 1,fun_nod(n)%numnod_Fx                 ! k - local number of cell
        m = fun_nod(n)%nods_Fx(k)                 ! m - global number of cell
        i = cells2cont(m)                         ! i - contour number
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
if ( solv == 'SOR2' .or. solv == 'SOR3') then
    print '(a,a,a, $)', 'init the solver multicolor ',solv,' ...'
    !   Launching a multi-color graph coloring process
    allocate(colors(num_k), source=0)
    call greedy_color_csr(num_k, irow_Zkk, jcol_Zkk, colors, num_colors)
    !  Block reordering and diagonal extraction
    allocate(diag(num_k), source=0.d0)
    allocate(ia_new(num_k+1), old_to_new(num_k), color_ptr(num_colors+1), source=0 )
    allocate (ja_diag(num_k), source=0)
    call permute_csr_split_diag(num_k, irow_Zkk, jcol_Zkk, val_Zkk, colors, num_colors, &
                                ia_new, ja_new, a_new, diag, old_to_new, color_ptr, ja_diag,ja_nodiag )

    allocate(b_perm(num_k),x_perm(num_k), source=0.d0 )
    b_perm(old_to_new) = Ek
    ! Matrix analysis for Omega search
    call estimate_optimal_omega(num_k, ia_new, ja_new, a_new, diag, b_perm, &
                                num_colors, color_ptr, omega_opt)
    ! Initialize x_perm for the main calculation.
    x_perm = 0.d0
    ! Main calculation 
    if (solv == 'SOR2') then
        call sor_multicolor(num_k, ia_new, ja_new, a_new, diag, b_perm, x_perm, &
                                omega_opt, itmax, tol, num_colors, color_ptr, res_norm, iter)
    else 
        call sor_multicolor_omp(num_k, ia_new, ja_new, a_new, diag, b_perm, x_perm, &
                                omega_opt, itmax, tol, num_colors, color_ptr, res_norm, iter)
    endif
    Ik = x_perm(old_to_new)
    call system_clock(k2);
    Tsavedata =  real(k2-k1)/real(k0)
ELSEIF (solv == 'SOR1') then
    print '(a, $)', 'init the solver SOR simpl...'
    call sor_sparse (num_k, val_Zkk, jcol_Zkk, irow_Zkk, Ek, Ik, omega, itmax, tol, iter)
    call system_clock(k2);
    Tsavedata =  real(k2-k1)/real(k0)
else
    print*,"error: solver only 'SOR1' or 'SOR2' or 'SOR3' "
    stop
endif

if (iter == itmax  ) then
    print'(a, g10.3)', ' Warning! iterations  reached limit, possible non-convergence. Time=', Tsavedata
    print'(a)', ' The results are written to a file field_0.vtk'
    ! calculation of currents in branches
    do  i=1,num_br
        ip=irow_Cvk(i); jp=irow_Cvk(i+1)-1
        Iv(i) =  dot_product( val_Cvk(ip:jp), Ik(jcol_Cvk(ip:jp) ) )
    enddo 
    call writeVtk (Npoint, sdx, sdy, sdz, Cells, cells2br, delta, iv,  files)
    stop
else
    if (solv == 'SOR2' .or. solv == 'SOR3') then
        print '(a, i5,a, i5, *(a, g10.3))', &
        ' Done, number of init iter for linear equations =', iter, ' num_colors= ', num_colors,' omega=',omega_opt,&
        ' err=',res_norm, ' time=',Tsavedata
    ELSEIF (solv == 'SOR1') then
        print '(a, i5, *(a, g10.3))', &
        ' Done, number of init iter for linear equations =', iter, ' err=',tol,' omega=',omega, ' time=',Tsavedata
    endif
endif
if (siz_non_lin_MU /= 0 ) then
    print '(a, i2,a)', ' start nonlin iterations, mode=', mode,' :'
    summ_mu_mean = 1.d7 ! init summ_mu_mean
endif

!******** CYCLE PREPARATION by time
call system_clock(k1,k0);
T=0.
Ntime=0                              ! counter of each point with step dt
Nout_display=nint(Time/(100.0*DT))   ! jump size for displaying a character '>' on the screen 
Nprint_display=Ntime + Nout_display  ! point counter for displaying the '>' symbol with a jump Nout_display
Npoint=0                ! calculation point counter with step dtt
Nout=nint(DTT/DT)       ! jump size for output to files
Nprint = Ntime + Nout   ! point counter with step Nout for output to files with jump Nout
!----------------------------------------------------------------!
! preparation of arrays for calculating the movement of ferromag !
!----------------------------------------------------------------!
flag_move_Fe = 0
IF (numVenv /=0 ) THEN
    DO i=1,numVenv  
    ! motion search processing
        fun_Venv(i)%Distance = 0.d0 
        IF (fun_Venv(i)% num_Venv(1) == 0  .and. fun_Venv(i)%move(1) /=0) THEN
        !  distance traveled in time dt with speed Vx in fractions of cell size along dX
            fun_Venv(i)%shift(1) = fun_Venv(i)%val_Venv(1)*dt/delta(1)   !Vx*dt/dX
            flag_move_Fe = 1
        ELSEIF (  fun_Venv(i)%move(1) /=0) THEN 
            flag_move_Fe = 1
        ENDIF
        IF (fun_Venv(i)% num_Venv(2) == 0  .and. fun_Venv(i)%move(2) /=0) THEN 
            fun_Venv(i)%shift(2) = fun_Venv(i)%val_Venv(2)*dt/delta(2)
            flag_move_Fe = 1; 
        ELSEIF (  fun_Venv(i)%move(2) /=0) THEN 
            flag_move_Fe = 1
        ENDIF
        IF (fun_Venv(i)% num_Venv(3) == 0 .and. fun_Venv(i)%move(3) /=0 ) THEN 
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
            IF (fun_nod(n)% num_Vmech(i) == 0  .and. fun_nod(n)%move(i) /=0) THEN 
                !  distance traveled in time dt with speed Vx in fractions of cell size along dX, dY, dZ
                fun_nod(n)%shift(i) = fun_nod(n)%val_Vmech(i)*dt/delta(i)   !Vx*dt/dX
                flag_move = 1
            ELSEIF (  fun_nod(n)%move(i) /=0) THEN 
                flag_move = 1
                ! fun_nod(n)%Distance(i) = - (Vmech(fun_nod(n)%num_Vmech(i))%vely * 1.5d0*dt/ delta(i) )
                fun_nod(n)%Distance(i) =  - (Vmech(fun_nod(n)%num_Vmech(i))%vely * dt/ delta(i) )
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
            m = fun_nod(n)%nods_Fx(k)      ! k - local number m - global number
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
                    fun_nod(n)%Distance(i) = fun_nod(n)%Distance(i) + Vmech(fun_nod(n)%num_Vmech(i))%vely*dt/delta(i)
                    fun_nod(n)%length(i) = int(fun_nod(n)%Distance(i)  )
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
                m = fun_nod(n)%nods_Fz(k)-2*Cells
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
        DO k = 1,fun_nod(n)%numnod_Fx                 ! k - local number 
            m = fun_nod(n)%nods_Fx(k)
            i = cells2cont(m)                         ! i - contour number  m - global number
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
!====================== for sor_multicolor
if (solv == 'SOR2' .or. solv == 'SOR3') b_perm(old_to_new) = Ek
! 
if (siz_non_lin_MU /= 0 ) then ! only for those domains where there are nonlinear mu
            !------------------------------------------
            ! Start iterations for nonlinear regions
            !------------------------------------------
    local_iter = 0
    do iterac = 1, itmag
        summ_mu_mean_old = summ_mu_mean  ! sum of mean mu for each domain
        !================================== simple SOR method  ==================
        if (solv == 'SOR1') then
            call sor_sparse (num_k, val_Zkk, jcol_Zkk, irow_Zkk, Ek, Ik, omega, itmax, tol, iter)
        else
        !=======================================colorized SOR method=============
            diag  = val_Zkk(ja_diag)         ! ja_diag(num_k)
            a_new = val_Zkk(ja_nodiag)       ! ja_diag(:)
            if (solv == 'SOR2') then
                call sor_multicolor(num_k, ia_new, ja_new, a_new, diag, b_perm, x_perm, &
                                omega_opt, itmax, tol, num_colors, color_ptr, res_norm, iter)
            else
                call sor_multicolor_omp(num_k, ia_new, ja_new, a_new, diag, b_perm, x_perm, &
                                omega_opt, itmax, tol, num_colors, color_ptr, res_norm, iter)
            endif
            Ik = x_perm(old_to_new)
        endif
        ! Iv = Cvk * Ik
        do  i=1,num_br
            ip=irow_Cvk(i); jp=irow_Cvk(i+1)-1
            Iv(i) =  dot_product( val_Cvk(ip:jp), Ik(jcol_Cvk(ip:jp) ) )
        enddo 
        summ_mu_mean=0.d0   
        do n = 1, sizdom_Fe                                 ! Calculation of mu for all ferromagnetic domains
            if (dom_Fe(n)%nonlin == 0 ) CYCLE               ! for linear domains skip
            MU_sum = 0.d0
            do m = 1, dom_Fe(n)%num_cells                  ! number of cells in a domain
                ! for each m-th cell, branch numbers
                Sx = Iv( dom_Fe(n)%num_brX(m))  !current in X-branches
                Sy = Iv( dom_Fe(n)%num_brY(m))  !current in Y-branches
                Sz = Iv( dom_Fe(n)%num_brZ(m))  !current in Z-branches
                a = sqrt(Sx*Sx + Sy*Sy + Sz*Sz) !amplitude of the induction flux density  
                !--------------------------------------------------------
                ! Calculation of magnetic permeability from tabular data 
                if (dom_Fe(n)%typ_nonlin == 1) &
                MU_new = get_mu_from_b_tab(a, dom_Fe(n)%b_dat, dom_Fe(n)%MU_dat, dom_Fe(n)%n_dat) 
                
                if (dom_Fe(n)%typ_nonlin == 2) &
                MU_new = get_mu_from_bh_tab(a, dom_Fe(n)%b_dat, dom_Fe(n)%h_dat, dom_Fe(n)%n_dat) 
                !--------------------------------------------------------
                if (mode == 1) dom_Fe(n)%mu_cells(m) = MU_new 
                MU_sum = MU_sum + MU_new                            ! sum of the permeabilities of all cells of the domain
            enddo
            dom_Fe(n)%mu_mean = real(dom_Fe(n)%num_cells,8)/ MU_sum     ! mean magnetic reluctivity of the domain
            if (mode == 2 .and. omegamag /= 0.d0) then ! recount dom_Fe(n)%mu_mean
                ! under-relaxation for mean magnetic reluctivity of the domain
                dom_Fe(n)%mu_mean = (1.d0 - omegamag)*dom_Fe(n)%mu_mean + omegamag*dom_Fe(n)%mu_mean_old
                dom_Fe(n)%mu_mean_old = dom_Fe(n)%mu_mean
            endif
            nn = dom_Fe(n)%n_Fe 
            valPHYS(nn,1) = dom_Fe(n)%mu_mean
            summ_mu_mean = summ_mu_mean + dom_Fe(n)%mu_mean  ! accumulation summ_mu_mean
        enddo
        if (mode /= 1) then
            call calc_Zvv23
        else
            call calc_Zvv1
        endif
        ! CZkv = Ckv * Zvv
        val_CZkv = 0.d0
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
            n=m/2
            allocate(ipvt(1:n), source=0)
            call FLASH2 (jcol_Zkk(ip:jp), val_Zkk(ip:jp), m, ipvt, n)
            deallocate(ipvt)
        enddo
        ! recalculation of summ_mu_mean taking into account the number of domains
        summ_mu_mean = summ_mu_mean/ real(siz_non_lin_MU)
        a = abs( (summ_mu_mean - summ_mu_mean_old) / summ_mu_mean_old )
        if ( a /= 0.d0) then
            if (mod(iterac,5) == 0) then
                print '( i0,$ )',  iterac
            else 
                print '( a,$ )', '*'
            endif
        endif
        if (mode == 1) then
            if (  a < tolmag  ) then  
                if (local_iter <= 6 ) then
                    local_iter = local_iter + 1
                    if (time == dt) print  '( a, i4, a, i3, a, e10.3, a )',&  
                    'add iter=', iterac, ' num local_iter=',  local_iter, ' err=', a
                ELSE 
                    if (time == dt) then
                        print  '( a, i4, 3(a, e10.3) )','num nonlin iter=', iterac, &
                    ' Mu_mean=', 1.d0/summ_mu_mean , ' Mu_mean_old=', 1.d0/summ_mu_mean_old , ' err=', a
                    else
                        print '( a )', '|'
                    endif
                    local_iter = 0
                    exit
                endif
            else
                if (  local_iter /= 0 ) local_iter = 0
            endif
            if ( iterac==itmag ) print  '( a, i4, 3(a, e10.3) )','num nonlin iter=', iterac, &
                    ' Mu_mean=', 1.d0/summ_mu_mean , ' Mu_mean_old=', 1.d0/summ_mu_mean_old , ' err=', a
        else
            if (  a < tolmag  ) then 
                if (time == dt) then
                    print  '( a, i4, 3(a, e10.3) )','num nonlin iter=', iterac, &
                    ' Mu_mean=', summ_mu_mean , ' Mean_old=', summ_mu_mean_old , ' err=', a
                else
                    print '( a )', '|'
                endif
                exit
            endif
            if ( iterac==itmag ) print  '( a, i4, 3(a, e10.3) )','num nonlin iter=', iterac, &
                    ' Mu_mean=', summ_mu_mean , ' Mu_mean_old=', summ_mu_mean_old , ' err=', a
        endif
    enddo ! end iterac
else  ! there are no nonlinear regions
    IF (flag_move_Fe == 0 ) then   ! Fe stat and coils stat .and. flag_move == 0
        if (solv == 'SOR1') then
            call sor_sparse (num_k, val_Zkk, jcol_Zkk, irow_Zkk, Ek, Ik, omega, itmax, tol, iter)
        else
            if (solv == 'SOR2') then
                call sor_multicolor(num_k, ia_new, ja_new, a_new, diag, b_perm, x_perm, &
                                omega_opt, itmax, tol, num_colors, color_ptr, res_norm, iter)
            else
                call sor_multicolor_omp(num_k, ia_new, ja_new, a_new, diag, b_perm, x_perm, &
                                omega_opt, itmax, tol, num_colors, color_ptr, res_norm, iter)
            endif
            Ik = x_perm(old_to_new)
        endif
        do  i=1,num_br
            ip=irow_Cvk(i); jp=irow_Cvk(i+1)-1
            Iv(i) =  dot_product( val_Cvk(ip:jp), Ik(jcol_Cvk(ip:jp) ) )
        enddo 
    endif
endif
!------------------------------------------------------!
! Calculation of the speed of movement of ferromagnets !
!------------------------------------------------------!    
IF (numVenv /= 0) THEN
    ! fix fun_Venv((n)%ism(i) based on the results of a comparison of new and old speeds
    CALL initf (numVenv)
    DO n=1,numVenv
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
    enddo  
ENDIF

IF (flag_move_Fe == 1) then  ! The motion of a ferromagnet has been established.
    do n = 1, numVenv        ! in all areas of ferromagnets that are moving
        !-------------
        n_Fe = Venv(n)%nomsch  ! moving domain number
        ! Before moving, store the numbers of the moving Fe in geoPHYS(i,j,k)
        do k=1,sdz; do j=1,sdy; do i=1,sdx
            if (geoPHYS_Fe(i,j,k) == n_Fe ) then
                ! In geoPHYS (ijk) write the number of the Fe region
                geoPHYS(i,j,k) = n_Fe   
                ! In geoPHYS_Fe(ijk_n) write air instead of ferromagnet.
                geoPHYS_Fe (i,j,k) = int(nsub_glob,1) 
            endif
        enddo; enddo; enddo  
        !----------calculation of new coordinates
        DO i=1,3
            IF (fun_Venv(n)% num_Venv(i) == 0 ) THEN     !
                fun_Venv(n)%Distance(i) =  movestop_Fe(i)* fun_Venv(n)%shift(i)
                fun_Venv(n)%length(i) = nint(fun_Venv(n)%Distance(i))
            ELSE! 
                fun_Venv(n)%Distance(i) = Venv(fun_Venv(n)%num_Venv(i))%vely * dt/delta(i) ! v *dt/dx
                ! fun_Venv(n)%length(i) = nint(  fun_Venv(n)%Distance(i)  ) 
                fun_Venv(n)%length(i) = int(  fun_Venv(n)%Distance(i)  )
            ENDIF
        ENDDO
        d_i = fun_Venv(n)%length(1)
        d_j = fun_Venv(n)%length(2) 
        d_k = fun_Venv(n)%length(3)
        ! The first step is to calculate new coordinates for the area, 
        ! save them to the buffer, and erase the old geoPHYS area.
        do k=1,sdz; do j=1,sdy; do i=1,sdx; 
            if (geoPHYS (i,j,k) == n_Fe ) then ! for Fe region number
                ! then find new ijk
                call new_ijk( i,j,k, d_i, d_j, d_k, movestop_Fe,  m_new,i_new,j_new,k_new )
                ! in geoPHYS_Fe(ijk_n) write the Fe region number
                geoPHYS_Fe(i_new,j_new,k_new) = n_Fe
                ! In geoPHYS (ijk) enter the air number
                geoPHYS (i,j,k) = int(nsub_glob,1) ! valPHYS(nsub_glob,1)
            endif
        enddo; enddo; enddo
        ! overwrite from buffer in geoPHYS
        do k=1,sdz; do j=1,sdy; do i=1,sdx;
            if (geoPHYS_Fe (i,j,k) == n_Fe ) then
                geoPHYS (i,j,k) = geoPHYS_Fe (i,j,k)  ! 
            endif
        enddo; enddo; enddo  
    enddo
    ! resistance calculation
    do n = 1, Cells  ! In all cells
        i=mod(n-1,sdx)+1; j=mod((n-1)/sdx,sdy)+1; k=(n-1)/kdz+1
        ip=cells2br(n,1)
        jp=cells2br(n,2)
        kp=cells2br(n,3)
        m = int(geoPHYS(i,j,k),4);
        a = valPHYS(m,1)
        if (ip /=0) then
            mm = int(geoPHYS(i+1,j,k),4); 
            s = valPHYS(mm,1)
            Zvv(ip)  = 2.d0*SmuX/( 1.d0/a + 1.d0/s)
        endif
        if (jp /=0) then
            mm = int(geoPHYS(i,j+1,k),4); 
            s = valPHYS(mm,1)
            Zvv(jp)  = 2.d0*SmuY/( 1.d0/a + 1.d0/s)
        endif
        if (kp /=0) then
            mm = int(geoPHYS(i,j,k+1),4); 
            s = valPHYS(mm,1)
            Zvv(kp)  = 2.d0*SmuZ/( 1.d0/a + 1.d0/s)
        endif
    enddo
    do n = 1, sizdom_Fe        ! Initialization of mu for all domain cells   
        mm =  dom_Fe(n)%n_Fe   ! domain number
        dom_Fe(n)%mu_cells = valPHYS(mm,1)
    enddo
    ! Recalculation of the contours resistance matrix
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
        n=m/2
        allocate(ipvt(1:n), source=0)
        call FLASH2 (jcol_Zkk(ip:jp), val_Zkk(ip:jp), m, ipvt, n)
        deallocate(ipvt)
    enddo
    
    if (solv == 'SOR1') then
        call sor_sparse (num_k, val_Zkk, jcol_Zkk, irow_Zkk, Ek, Ik, omega, itmax, tol, iter)   
    else
        diag  = val_Zkk(ja_diag)
        a_new = val_Zkk(ja_nodiag)  
        if (solv == 'SOR2') then
            call sor_multicolor(num_k, ia_new, ja_new, a_new, diag, b_perm, x_perm, &
                                omega_opt, itmax, tol, num_colors, color_ptr, res_norm, iter)
        else 
            call sor_multicolor_omp(num_k, ia_new, ja_new, a_new, diag, b_perm, x_perm, &
                                omega_opt, itmax, tol, num_colors, color_ptr, res_norm, iter)
        endif
        Ik = x_perm(old_to_new)
    endif
    do  i=1,num_br
        ip=irow_Cvk(i); jp=irow_Cvk(i+1)-1
        Iv(i) =  dot_product( val_Cvk(ip:jp), Ik(jcol_Cvk(ip:jp) ) )
    enddo 

endif  !flag_move_Fe
! ================END MOVING=================
            !--------------------------------------------
            !  END OF CALCULATION AND OUTPUT OF RESULTS
            !---------------------------------------------
if ( Ntime >= Nprint .or. Ntime == 0 ) then  ! output with step dtt skipping 1st point
    Nprint = Ntime  + Nout
    Npoint= Npoint + 1  ! point counter with dtt step
    ios=0
    
    !$OMP barrier
    call writeVtk (Npoint, sdx, sdy, sdz, Cells, cells2br, delta, iv, files)
        
    IF (flag_move == 1 .or. Npoint == 1) &
        CALL writeVtk_src ( Npoint, numfun,  NcellsX, NcellsY, NcellsZ,  new_nodesX, new_nodesY,  new_nodesZ,  &
                            sdx, sdy, delta, files) 
        if (sizdom_Fe /=0) call writeVtk_Fe (Npoint, sdx, sdy, sdz, Cells,  delta, nsub_glob,  geoPHYS_Fe, files)
    if (siz_non_lin_MU /= 0 ) then
        print '( a, i0, a, $ )', 'Point>', Npoint, ' iter:'
    else
        if (mod(Npoint,5) == 0) then
            print '( i0,$ )',  Npoint
        else 
            print '( a,$ )', '>'
        endif
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
Tcalc =  Tsavedata + real(k2-k1)/real(k0)
print '( a, a, a, g10.3)', 'The results in VTK format have been saved in the ', '"'//trim(files)//'"', ' folder. Total time=',  Tcalc 

!=========================================
contains
subroutine bnd_cells  
        do nn = 1, sizdom_Fe 
            dom_Fe(nn)%cells = 0
            L=0                                       ! domain number counter Fe
            do k = 1,sdz;  do j = 1,sdy;  do i = 1,sdx
                m = i + sdx * ((j -1) + sdy*(k -1))    ! global cell number
                if ( geoPHYS (i,j,k) == int (dom_Fe(nn)%n_Fe, 1)  ) then 
                    L = L+1                           ! counter cells local
                    dom_Fe(nn)%cells(L) = m
                endif
            end do;  end do;  end do
        enddo 	
        ! calculation of boundary cells if the main ones are defined
        do n = 1, sizdom_Fe
            ! if (dom_Fe(n)%nonlin == 0 ) cycle               ! for linear domains skip
            do m = 1, dom_Fe(n)%num_cells                     ! number of cells in a domain
                nn = dom_Fe(n)%cells(m)                       ! cell number ferr ijk
                i=mod(nn-1,sdx)+1; j=mod((nn-1)/sdx,sdy)+1; k=(nn-1)/kdz+1
                !============================= BND X
                i_p = i+1 + sdx * ((j -1) + sdy*(k -1))       ! cell number i+1
                nod = FINDLOC (dom_Fe(n)%cells, i_p)          ! local index in a list dom_Fe(n)%cells
                nn = int(geoPHYS(i+1,j,k),4 )         ! global domain number on the right
                if ( nn == dom_Fe(n)%n_Fe ) then      ! current domain? =numdom_Fe(n),1
                    dom_Fe(n)%cellsXp(m,1) = nod(1)   !  if yes, remember the  local cell index.
                    dom_Fe(n)%cellsXp(m,2) = n        !  remember the local  number of domain
                else                                  ! nn -  it is neighboring domain.
                    if (nn == nsub_glob) then              !  then it is air
                        dom_Fe(n)%cellsXp(m,1) = i_p       ! remember the global cell number of air on the right.
                        dom_Fe(n)%cellsXp(m,2) = nsub_glob !   remember domain number on the right
                    else                                   ! right domain is not air. 
                        nod = FINDLOC (dom_Fe%n_Fe, nn)    ! local number of Fe in the list dom_Fe()%n_Fe 
                        mm = nod(1)                        
                        if (mm == 0 ) stop 'not in list Fex'
                        nod = FINDLOC (dom_Fe(mm)%cells, i_p)   ! is this cell i_p in this domain?
                        if (nod(1) == 0 ) stop 'no Fex'  
                        dom_Fe(n)%cellsXp(m,1) = nod(1)         ! remember local number of the neighboring cell .
                        dom_Fe(n)%cellsXp(m,2) = mm             ! local number of the neighboring domain mm
                    endif
                endif
                !============================= BND Xm
                i_p = i-1 + sdx * ((j -1) + sdy*(k -1))
                nod = FINDLOC (dom_Fe(n)%cells, i_p)
                nn = int(geoPHYS(i-1,j,k),4 )
                if ( nn == dom_Fe(n)%n_Fe ) then
                    dom_Fe(n)%cellsXm(m,1) = nod(1)
                    dom_Fe(n)%cellsXm(m,2) = n
                else
                    if (nn == nsub_glob) then
                        dom_Fe(n)%cellsXm(m,1) = i_p
                        dom_Fe(n)%cellsXm(m,2) = nsub_glob
                    else
                        nod = FINDLOC (dom_Fe%n_Fe, nn)
                        mm = nod(1)
                        if (mm == 0 ) stop 'no Fex'
                        nod = FINDLOC (dom_Fe(mm)%cells, i_p)
                        if (nod(1) == 0 ) stop 'no Fex'
                        dom_Fe(n)%cellsXm(m,1) = nod(1)
                        dom_Fe(n)%cellsXm(m,2) = mm
                    endif
                endif
                !============================= BND Yp
                j_p = i + sdx * ((j+1 -1) + sdy*(k -1))
                nod = FINDLOC (dom_Fe(n)%cells, j_p)
                nn = int(geoPHYS(i,j+1,k),4 )

                if ( nn == dom_Fe(n)%n_Fe ) then
                    dom_Fe(n)%cellsYp(m,1) = nod(1)
                    dom_Fe(n)%cellsYp(m,2) = n 
                else
                    if (nn == nsub_glob) then
                        dom_Fe(n)%cellsYp(m,1) = j_p
                        dom_Fe(n)%cellsYp(m,2) = nsub_glob
                    else
                        nod = FINDLOC (dom_Fe%n_Fe, nn)
                        mm = nod(1)
                        if (mm == 0 ) stop 'no Feyp'
                        nod = FINDLOC (dom_Fe(mm)%cells, j_p)
                        if (nod(1) ==0 ) stop 'no Feyp'     
                        dom_Fe(n)%cellsYp(m,1) = nod(1)
                        dom_Fe(n)%cellsYp(m,2) = mm
                    endif
                endif
                !============================= BND Ym
                j_p = i + sdx * ((j-1 -1) + sdy*(k -1))
                nod = FINDLOC (dom_Fe(n)%cells, j_p)
                
                nn = int(geoPHYS(i,j-1,k),4 )
                if ( nn == dom_Fe(n)%n_Fe ) then 
                    dom_Fe(n)%cellsYm(m,1) = nod(1)
                    dom_Fe(n)%cellsYm(m,2) = n
                else    
                    if (nn == nsub_glob) then 
                        dom_Fe(n)%cellsYm(m,1) = j_p 
                        dom_Fe(n)%cellsYm(m,2) = nsub_glob 
                    else    
                        nod = FINDLOC (dom_Fe%n_Fe, nn) 
                        mm = nod(1) 
                        if (mm == 0 ) stop 'no Feym' 
                        nod = FINDLOC (dom_Fe(mm)%cells, j_p) 
                        if (nod(1) ==0) stop 'no Feym'
                        dom_Fe(n)%cellsYm(m,1) = nod(1)
                        dom_Fe(n)%cellsYm(m,2) = mm
                    endif
                endif
                !============================= BND Zp
                k_p = i + sdx * ((j -1) + sdy*(k+1 -1))
                nod = FINDLOC (dom_Fe(n)%cells, k_p)
                nn = int(geoPHYS(i,j,k+1),4 )
                if ( nn == dom_Fe(n)%n_Fe ) then  
                    if (nod(1) == 0) stop ' not find cell in domain Fe Z'
                    dom_Fe(n)%cellsZp(m,1) = nod(1)
                    dom_Fe(n)%cellsZp(m,2) = n  
                else     
                    if (nn == nsub_glob) then   
                        dom_Fe(n)%cellsZp(m,1) = k_p 
                        dom_Fe(n)%cellsZp(m,2) = nsub_glob 
                    else               
                        nod = FINDLOC (dom_Fe%n_Fe, nn)  
                        mm = nod(1)
                        if (mm == 0 ) stop 'no Fez' 
                        nod = FINDLOC (dom_Fe(mm)%cells, k_p)  
                        if (nod(1) ==0 ) stop 'no Fez' 
                        if (nod(1) > dom_Fe(mm)%num_cells ) STOP 'nod(1) > dom_Fe(mm)%num_cells for Z'
                        dom_Fe(n)%cellsZp(m,1) = nod(1) 
                        dom_Fe(n)%cellsZp(m,2) = mm 
                    endif
                endif
!============================= BND Zm
                k_p = i + sdx * ((j -1) + sdy*(k-1 -1)) 
                nod = FINDLOC (dom_Fe(n)%cells, k_p) 
                nn = int(geoPHYS(i,j,k-1),4 ) 
                if ( nn == dom_Fe(n)%n_Fe ) then  
                    if (nod(1) == 0) stop ' not find cell in domain Fe Z'
                    dom_Fe(n)%cellsZm(m,1) = nod(1) 
                    dom_Fe(n)%cellsZm(m,2) = n  
                else      
                    if (nn == nsub_glob) then  
                        dom_Fe(n)%cellsZm(m,1) = k_p 
                        dom_Fe(n)%cellsZm(m,2) = nsub_glob 
                    else      
                        nod = FINDLOC (dom_Fe%n_Fe, nn) 
                        mm = nod(1)
                        if (mm == 0 ) stop 'no Fez'
                        nod = FINDLOC (dom_Fe(mm)%cells, k_p)  
                        if (nod(1) ==0 ) stop 'no Fez' 
                        if (nod(1) > dom_Fe(mm)%num_cells ) STOP 'nod(1) > dom_Fe(mm)%num_cells for Z'
                        dom_Fe(n)%cellsZm(m,1) = nod(1) 
                        dom_Fe(n)%cellsZm(m,2) = mm 
                    endif
                endif
            enddo
        enddo
end subroutine bnd_cells

subroutine calc_Zvv1
    do n = 1, sizdom_Fe
        if (dom_Fe(n)%nonlin == 0 ) CYCLE               ! for linear domains skip
        do m = 1, dom_Fe(n)%num_cells                   ! number of cells in a domain
            a = 1.d0/dom_Fe(n)%mu_cells(m)              ! for each m-th cell 1/mu_fe
            nn = dom_Fe(n)%cells(m)                     ! global cell number
            ! for cells i+1
            if (dom_Fe(n)%cellsXp(m,2) == nsub_glob) then ! it is air
                s = 1.d0
            else                                          ! on the right is either a neighboring domain or own
                if (dom_Fe(n)%cellsXp(m,2) /= n ) then    ! not one's own
                    i_p = dom_Fe(n)%cellsXp(m,2)          ! local number of the neighboring domain Fe
                    mm = dom_Fe(n)%cellsXp(m,1)           ! mm - local cell index at the boundary 
                    s = 1.d0/dom_Fe(i_p)%mu_cells(mm)
                else                                      ! on the right a own domain
                    mm  = dom_Fe(n)%cellsXp(m,1)          ! mm - local cell number to the right of its domain
                    s = 1.d0/dom_Fe(n)%mu_cells(mm)
                endif
            endif
            L = cells2br(nn,1)                            ! branch number X
            Zvv(L) =  2.d0*(1.d0 - omegamag)* SmuX/( 1.d0/a + 1.d0/s) + omegamag* Zvv(L)
            ! for cells j+1 similar to the previous 
            if (dom_Fe(n)%cellsYp(m,2) == nsub_glob) then 
                s = 1.d0
            else   
                if (dom_Fe(n)%cellsYp(m,2) /= n ) then    
                    j_p = dom_Fe(n)%cellsYp(m,2)          
                    mm = dom_Fe(n)%cellsYp(m,1)           
                    s = 1.d0/dom_Fe(j_p)%mu_cells(mm)
                else                                      
                    mm  = dom_Fe(n)%cellsYp(m,1)          
                    s = 1.d0/dom_Fe(n)%mu_cells(mm)
                endif
            endif
            L = cells2br(nn,2)
            Zvv(L) =  2.d0*(1.d0 - omegamag)* SmuY/( 1.d0/a + 1.d0/s) + omegamag* Zvv(L)
            ! for cells k+1 similar to the previous
            if (dom_Fe(n)%cellsZp(m,2) == nsub_glob) then
                s = 1.d0
            else   
                if (dom_Fe(n)%cellsZp(m,2) /= n ) then   
                    k_p = dom_Fe(n)%cellsZp(m,2)         
                    mm  = dom_Fe(n)%cellsZp(m,1)         
                    s = 1.d0/dom_Fe(k_p)%mu_cells(mm)
                else                                     
                    mm  = dom_Fe(n)%cellsZp(m,1)         
                    s = 1.d0/dom_Fe(n)%mu_cells(mm)
                endif
            endif
            L = cells2br(nn,3) 
            Zvv(L) =  2.d0*(1.d0 - omegamag)* SmuZ/( 1.d0/a + 1.d0/s) + omegamag* Zvv(L)
            ! for the air neighboring
            ! for cells i-1  
            if (dom_Fe(n)%cellsXm(m,2) == nsub_glob ) then
                mm = dom_Fe(n)%cellsXm(m,1)      ! mm - cell number in the air i-1 
                if (mm /=0 ) then 
                    L = cells2br(mm,1)           ! number of branch X in the air
                    Zvv(L) =  2.d0*(1.d0 - omegamag)* SmuX/( 1.d0/a + 1.d0) + omegamag* Zvv(L)
                endif
            endif
            ! for cells j-1  similar to the previous
            if (dom_Fe(n)%cellsYm(m,2) == nsub_glob ) then
                mm = dom_Fe(n)%cellsYm(m,1)
                if (mm /=0 ) then   
                    L = cells2br(mm,2) 
                    Zvv(L) =  2.d0*(1.d0 - omegamag)* SmuY/( 1.d0/a + 1.d0) + omegamag* Zvv(L)
                endif
            endif
            ! for cells k-1  similar to the previous
            if (dom_Fe(n)%cellsZm(m,2) == nsub_glob ) then
                mm = dom_Fe(n)%cellsZm(m,1)
                if (mm /=0 ) then  
                    L = cells2br(mm,3) 
                    Zvv(L) =  2.d0*(1.d0 - omegamag)* SmuZ/( 1.d0/a + 1.d0) + omegamag* Zvv(L)
                endif
            endif
        enddo
    enddo
end subroutine calc_Zvv1

subroutine calc_Zvv23
! Similar to Zvv1, but the average mu (dom_Fe(n)%mu_mean) of the domain is used.
    do n = 1, sizdom_Fe
        if (dom_Fe(n)%nonlin == 0 ) CYCLE               ! for linear domains skip
        a = dom_Fe(n)%mu_mean
        do m = 1, dom_Fe(n)%num_cells                   ! number of cells in a domain
            nn = dom_Fe(n)%cells(m)                     ! global cell number
            !for cells i+1
            if (dom_Fe(n)%cellsXp(m,2) == nsub_glob) then ! it is air
                s = 1.d0
            else                   !  on the right is either a neighboring domain or own
                if (dom_Fe(n)%cellsXp(m,2) /= n ) then   ! not one's own
                    i_p = dom_Fe(n)%cellsXp(m,2)    ! local number of the neighboring domain Fe
                    s = dom_Fe(i_p)%mu_mean
                else                              ! on the right a own domain
                    s = dom_Fe(n)%mu_mean
                endif
            endif
            L = cells2br(nn,1)                            ! branch number X
            Zvv(L) =  0.5d0*( a + s) * SmuX
            ! for cells j+1  similar to the previous
            if (dom_Fe(n)%cellsYp(m,2) == nsub_glob) then 
                s = 1.d0
            else
                if (dom_Fe(n)%cellsYp(m,2) /= n ) then 
                    j_p = dom_Fe(n)%cellsYp(m,2)  
                    s = dom_Fe(j_p)%mu_mean
                else
                    s = dom_Fe(n)%mu_mean
                endif
            endif
            L = cells2br(nn,2)
            Zvv(L) =  0.5d0*( a + s) * SmuY
            ! for cells k+1  similar to the previous
            if (dom_Fe(n)%cellsZp(m,2) == nsub_glob) then
                s = 1.d0
            else   
                if (dom_Fe(n)%cellsZp(m,2) /= n ) then 
                    k_p = dom_Fe(n)%cellsZp(m,2)
                    s = dom_Fe(k_p)%mu_mean
                else 
                    s = dom_Fe(n)%mu_mean
                endif
            endif
            L = cells2br(nn,3)   
            Zvv(L) =  0.5d0*( a + s) * SmuZ
                ! for air boundaries
                ! for cells i-1 
            if (dom_Fe(n)%cellsXm(m,2) == nsub_glob) then   ! Xm(m,2) - global domain number
                mm = dom_Fe(n)%cellsXm(m,1)                 ! cell number i-1 
                L = cells2br(mm,1)                          ! number of branch X in the air
                Zvv(L) =  0.5d0*( a + 1.d0) * SmuX
            endif
                            ! for cells  j-1 similar to the previous
            if (dom_Fe(n)%cellsYm(m,2) == nsub_glob) then
                mm = dom_Fe(n)%cellsYm(m,1)
                L = cells2br(mm,2)
                Zvv(L) =  0.5d0*( a + 1.d0) * SmuY
            endif
                            ! for cells k-1  similar to the previous
            if (dom_Fe(n)%cellsZm(m,2) == nsub_glob) then
                mm = dom_Fe(n)%cellsZm(m,1)  
                L = cells2br(mm,3) 
                Zvv(L) =  0.5d0*( a + 1.d0) * SmuZ
            endif
        enddo
    enddo
        
end subroutine calc_Zvv23

!-------------------
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

subroutine greedy_color_csr(n, ia, ja, colors, num_colors)
! Fortran code created by J.Sochor   ( https://github.com/JNSresearcher )
! e-mail: JNSresearcher@gmail.com
! August 2026
    integer, intent(in) :: n
    integer, intent(in) :: ia(n+1), ja(*)
    integer, intent(out) :: colors(n)
    integer, intent(out) :: num_colors
    integer :: i, j, neighbor, c
    logical, allocatable :: available_colors(:)

    allocate(available_colors(n))
    colors = 0
    num_colors = 0
!$omp parallel  private( i,j,c, neighbor)
!$omp do   
    do i = 1, n
        available_colors = .true.
        ! Outdoing the neighbors
        do j = ia(i), ia(i+1) - 1
            neighbor = ja(j)
            ! If the neighbor's index is valid and it has already been colored
            if (neighbor > 0 .and. neighbor <= n) then
                if (neighbor /= i .and. colors(neighbor) > 0) then
                    available_colors(colors(neighbor)) = .false.
                end if
            end if
        end do
        ! Assign the first available color (starting from 1).
        do c = 1, n
            if (available_colors(c)) then
                colors(i) = c
                if (c > num_colors) num_colors = c
                exit
            end if
        end do
    end do
!$omp end do
!$omp end parallel
    deallocate(available_colors)
end subroutine greedy_color_csr

    !------------------------------------------------------------------
    !  Reordering with diagonal extraction
    !------------------------------------------------------------------
subroutine permute_csr_split_diag(n, ia, ja, a, colors, num_colors, &
                                        ia_new, ja_new, a_new, diag, old_to_new, color_ptr, ja_diag,ja_nodiag ) 
! Fortran code created by J.Sochor   ( https://github.com/JNSresearcher )
! e-mail: JNSresearcher@gmail.com
! August 2026
    integer, intent(in) :: n
    integer, intent(in) :: ia(n+1), ja(*)
    real(8), intent(in) :: a(*)
    integer, intent(in) :: colors(n)
    integer, intent(in) :: num_colors
    
    integer, intent(out) :: ia_new(n+1)              ,ja_diag(n)
    integer, intent(out), allocatable :: ja_new(:)   ,ja_nodiag(:)
    real(8), intent(out), allocatable :: a_new(:)
    real(8), intent(out) :: diag(n)
    integer, intent(out) :: old_to_new(n)
    integer, intent(out) :: color_ptr(num_colors + 1)
    ! Auxiliary reverse mapping array: new index -> old index
    integer :: new_to_old(n)
    integer :: i, j, c, curr_row, old_row, nnz_count, count_nodes
    integer :: new_col
    
    old_to_new = 0
    new_to_old = 0
    ! construct a permutation map in both directions (without gaps).
    count_nodes = 0
    do c = 1, num_colors
        color_ptr(c) = count_nodes + 1
        do i = 1, n
            if (colors(i) == c) then
                count_nodes = count_nodes + 1
                old_to_new(i) = count_nodes
                new_to_old(count_nodes) = i
            end if
        end do
    end do
    color_ptr(num_colors + 1) = n + 1
    ! A safeguard in case some nodes did not receive a color in the range from 1 to num_colors.
    if (count_nodes /= n) then
        print *, "CRITICAL ERROR: Nodes colored: ", count_nodes, "from ", n
        print *, "Check the correctness of the color array!"
        stop
    end if
    ! Calculating the size of the new structure excluding the diagonal
    ia_new(1) = 1
    nnz_count = 0
    do curr_row = 1, n
        old_row = new_to_old(curr_row) ! Directly retrieving the old string without a loop-based search
        do j = ia(old_row), ia(old_row+1) - 1
            if (ja(j) /= old_row) then
                nnz_count = nnz_count + 1
            end if
        end do
        ia_new(curr_row+1) = nnz_count + 1
    end do
    allocate(ja_new(nnz_count))
    allocate(a_new(nnz_count))
    allocate(ja_nodiag(nnz_count))
    ! filling the new CSR structure and the diag array.
    nnz_count = 1
    diag = 0.d0
    do curr_row = 1, n
        old_row = new_to_old(curr_row)
        do j = ia(old_row), ia(old_row+1) - 1
            new_col = old_to_new(ja(j))
            if (ja(j) == old_row) then
                ja_diag(curr_row) = j
                diag(curr_row) = a(j)          !  diag = a(ja_diag)
            else
                ja_nodiag(nnz_count) = j
                ja_new(nnz_count) = new_col
                a_new(nnz_count) = a(j)        ! a_new = a(ja_nodiag)
                nnz_count = nnz_count + 1
            end if
        end do
    end do
end subroutine permute_csr_split_diag

subroutine sor_multicolor_omp(n, ia, ja, a, diag, b, x, omega, max_iter, tol, &
                                    num_colors, color_ptr, res_norm, iter)
! Fortran code created by J.Sochor   ( https://github.com/JNSresearcher )
! e-mail: JNSresearcher@gmail.com
! August 2026
    integer, intent(in) :: n
    integer, intent(in) :: ia(n+1), ja(*)
    real(8), intent(in) :: a(*)
    real(8), intent(in) :: diag(n)       ! The diagonal is transmitted separately.
    real(8), intent(in) :: b(n)
    real(8), intent(inout) :: x(n)
    real(8), intent(in) :: omega
    integer, intent(in) :: max_iter
    real(8), intent(in) :: tol
    integer, intent(in) :: num_colors
    integer, intent(in) :: color_ptr(num_colors + 1)
    real(8), intent(out) ::  res_norm
    integer , intent(out) :: iter
    integer :: c, i, j,k 
    real(8) :: sigma
    real(8) :: old_x(n), one_minus_omega, omega_over_diag

    one_minus_omega = 1.d0 - omega
    !$omp parallel private( k,c, i,j, omega_over_diag,sigma) shared(x, old_x, a, ia, ja, b, diag, one_minus_omega)
    do k=1,max_iter
        old_x = x
        do c = 1, num_colors
            !$omp do  schedule(guided) ! reduction(+:sigma) guided static
            do i = color_ptr(c), color_ptr(c+1) - 1
                sigma = 0.d0
                do j = ia(i), ia(i+1) - 1
                    sigma = sigma + a(j) * x(ja(j))
                end do
                ! Calculate the new value using the diag vector.
                omega_over_diag = omega / diag(i)
                x(i) = one_minus_omega * x(i) + omega_over_diag * (b(i) - sigma)
            end do
            !$omp end do nowait
        end do
        !$omp barrier

        !$omp single
        res_norm = norm2(x-old_x)/norm2(x)
        !$omp end single
        if (res_norm < tol  ) exit 
    end do
    iter=k
    !$omp end parallel
    
end subroutine sor_multicolor_omp

subroutine sor_multicolor(n, ia, ja, a, diag, b, x, omega, max_iter, tol, &
                                    num_colors, color_ptr, res_norm, iter)
! Fortran code created by J.Sochor   ( https://github.com/JNSresearcher )
! e-mail: JNSresearcher@gmail.com
! August 2026
    integer, intent(in) :: n
    integer, intent(in) :: ia(n+1), ja(*)
    real(8), intent(in) :: a(*)
    real(8), intent(in) :: diag(n)       ! The diagonal is transmitted separately.
    real(8), intent(in) :: b(n)
    real(8), intent(inout) :: x(n)
    real(8), intent(in) :: omega
    integer, intent(in) :: max_iter
    real(8), intent(in) :: tol
    integer, intent(in) :: num_colors
    integer, intent(in) :: color_ptr(num_colors + 1)
    real(8), intent(out) ::  res_norm
    integer , intent(out) :: iter
    integer :: c, i, j
    real(8) :: sigma
    real(8) :: old_x(n), one_minus_omega, omega_over_diag

    one_minus_omega = 1.d0 - omega
    do iter = 1, max_iter
        old_x = x
        do c = 1, num_colors
            do i = color_ptr(c), color_ptr(c+1) - 1
                sigma = 0.d0
                do j = ia(i), ia(i+1) - 1
                    sigma = sigma + a(j) * x(ja(j))
                end do
                ! Calculate the new value using the diag vector.
                omega_over_diag = omega / diag(i)
                x(i) = one_minus_omega * x(i) + omega_over_diag * (b(i) - sigma)
            end do
        end do
        if (mod(iter,10) == 0) then
            res_norm = norm2(x-old_x)/norm2(x)
            if (res_norm < tol  ) then
                exit
            endif
        endif
    end do
end subroutine sor_multicolor

subroutine sor_sparse(n, values, col_ind, row_ptr, b, x, omega, max_iter, tol, iter)
! Fortran code created by J.Sochor   ( https://github.com/JNSresearcher )
! e-mail: JNSresearcher@gmail.com
! June 2026

implicit none
integer, intent(in) :: n, max_iter
real(8), intent(in) :: values(*), b(n), omega, tol
integer, intent(in) :: col_ind(*), row_ptr(n+1)
integer, intent(out) :: iter
real(8), intent(inout) :: x(n)
integer :: i, j, k
real(8) :: old_x(n), diag, diff,  summ 

    do iter = 1, max_iter
        ! max_diff = 0.0d0
        old_x = x
        do i = 1, n
            summ = 0.d0
            diag = 0.0d0
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
        end do
        if (mod(iter,10) == 0) then
            diff = norm2(x-old_x)/norm2(x)
            if (diff < tol  ) then
                exit
            endif
        endif
    end do
end subroutine sor_sparse


!------------------------------------------------------------------
! PRELIMINARY OMEGA ASSESSMENT
!------------------------------------------------------------------
subroutine estimate_optimal_omega(n, ia, ja, a, diag, b, num_colors, color_ptr, omega_opt)
! Fortran code created by J.Sochor   ( https://github.com/JNSresearcher )
! e-mail: JNSresearcher@gmail.com
! August 2026
integer, intent(in) :: n
integer, intent(in) :: ia(n+1), ja(*)
real(8), intent(in) :: a(*)
real(8), intent(in) :: diag(n)
real(8), intent(in) :: b(n)
integer, intent(in) :: num_colors
integer, intent(in) :: color_ptr(num_colors + 1)
real(8), intent(out) :: omega_opt       ! Found parameter
integer :: iter, c, i, j, c_start, c_end
real(8) :: sigma, res_norm, row_res, rho
real(8) :: res_9, res_14  ! Residuals at the 9th and 14th iterations
real(8), allocatable :: x_test(:)

    allocate(x_test(n))
    x_test = 0.d0
    res_9 = 1.d0
    res_14 = 1.d0
    ! We perform exactly 15 Gauss-Seidel test iterations (omega = 1.0).
    do iter = 1, 15
        do c = 1, num_colors
            c_start = color_ptr(c)
            c_end = color_ptr(c+1) - 1
            !$omp parallel do private(i, j, sigma) shared(x_test, a, ia, ja, b, diag)
            do i = c_start, c_end
                sigma = 0.d0
                do j = ia(i), ia(i+1) - 1
                    sigma = sigma + a(j) * x_test(ja(j))
                end do
                ! Gauss-Seidel formula (omega = 1.0)
                x_test(i) = (b(i) - sigma) / diag(i)
            end do
            !$omp end parallel do
        end do
        !  calculate the residual norm at the control steps.
        if (iter == 9 .or. iter == 14) then
            res_norm = 0.d0
            !$omp parallel do private(i, j, row_res) reduction(+:res_norm)
            do i = 1, n
                row_res = b(i) - diag(i) * x_test(i)
                do j = ia(i), ia(i+1) - 1
                    row_res = row_res - a(j) * x_test(ja(j))
                end do
                res_norm = res_norm + row_res**2
            end do
            !$omp end parallel do
            if (iter == 9)  res_9 = sqrt(res_norm)
            if (iter == 14) res_14 = sqrt(res_norm)
        end if
    end do
    ! ANALYSIS OF CONVERGENCE RATE:
    ! The ratio of residuals over 5 steps (from the 9th to the 14th iteration) yields an approximation of rho^5.
    if (res_9 > 0.d0 .and. res_14 < res_9) then
        rho = (res_14 / res_9)**(1.d0 / 5.d0)
        ! Protection against exceeding mathematical limits (spectral radius < 1)
        if (rho >= 1.d0) then
            rho = 0.99d0
            print*,'Rho protection >= 1.0', rho
        endif
        if (rho <= 0.d0) then
            rho = 0.1d0
            print*,'Rho protection <= 0.0', rho
        endif
        ! Classic formula for converting the spectral radius into the optimal SOR omega.
        omega_opt = 2.d0 / (1.d0 + sqrt(1.d0 - rho**2))
    else
        ! The default safe option if something goes wrong.
        omega_opt = 1.d0 
    end if
    deallocate(x_test)
end subroutine estimate_optimal_omega

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

subroutine spgemm_symbolic(n_rows_a, n_cols_a, n_cols_b, &
                            row_ptr_a, col_ind_a, &
                            row_ptr_b, col_ind_b, &
                            row_ptr_c, total_nnz_c)
! Fortran code created by J.Sochor   ( https://github.com/JNSresearcher )
! e-mail: JNSresearcher@gmail.com
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
integer :: marker(n_cols_b) ! Stores the current line numberi
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
    ! Determine the number of non-zero elements in each column
    allocate(work(n))
    work = 0
    iat=0
    nz = ia(m+1) - 1  ! Total number of non-zero elements
    do k = 1, nz
        col = ja(k)
        work(col) = work(col) + 1
    end do
    ! Let's construct the array iat (analog ia for the transposed matrix)
    iat(1) = 1
    do j = 1, n
        iat(j+1) = iat(j) + work(j)
    end do
    ! Redistribute the elements!
    ! First, copy work to iat to temporarily store the offsets.
    work = iat(1:n)
    do i = 1, m
        do p = ia(i), ia(i+1)-1
            j = ja(p)
            k = work(j)
            at(k) = a(p)
            jat(k) = i   ! Row indices become column indices
            work(j) = work(j) + 1
        end do
    end do
    deallocate(work)
end subroutine csr_transpose


subroutine  sp_mm (row_a, col_a, col_b, irow_a,jcol_a,val_a, irow_b,jcol_b,val_b,  irow_c,jcol_c,val_c, curr_c_ptr )
! Matrix multiplication in CSR format
! Fortran code created by J.Sochor   ( https://github.com/JNSresearcher )
! e-mail: JNSresearcher@gmail.com
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

SUBROUTINE FLASH2 (A, S, N, L, M)
!     SORTS ARRY A WITH N ELEMENTS BY USE OF INDEX VECTOR L 
!     OF DIMENSION M WITH M ABOUT 0.1 N. 
!     Karl-Dietrich Neubert, FlashSort1 Algorithm 
!     in  Dr. Dobbs Journal Feb.1998,p.123

!     Adapted slightly; taken from the website 
!     https://www.neubert.net/FSOIntro.html

integer A(N),  HOLD, ANMIN, C1, FLASH
integer L(M)
real(8) s(N), tmp,  tmp2
!============================ CLASS FORMATION ===== 
ANMIN=A(1)
NMAX=1 
DO I=1,N
    IF( A(I).LT.ANMIN) ANMIN = A(I)
    IF( A(I).GT.A(NMAX)) NMAX = I
END DO
IF (ANMIN.EQ.A(NMAX)) RETURN
C1 = (M - 1) / (A(NMAX) - ANMIN)
DO I=1,N
    K = 1 + INT(C1 * (A(I) - ANMIN))
    L(K)=L(K) + 1
END DO
DO K=2,M
    L(K)=L(K) + L(K - 1)
END DO
HOLD=A(NMAX); A(NMAX) = A(1); A(1)=HOLD
tmp =S(NMAX);       S(NMAX) = S(1);   S(1)=tmp
! =============================== PERMUTATION ===== 
NMOVE=0 
J=1
K=M
DO WHILE (NMOVE .LT. N - 1)
    DO WHILE (J.GT.L(K)) 
        J=J + 1 
        K=1 + INT(C1 * (A(J) - ANMIN)) 
    END DO  
    FLASH=A(J)
    tmp2 =S(J)
    DO WHILE (.NOT.(J.EQ.L(K) + 1)) 
        K = 1 + INT(C1 * (FLASH - ANMIN))
        HOLD=A(L(K)) ; tmp=S(L(K))
        A(L(K))=FLASH; S(L(K)) = tmp2
        FLASH=HOLD;    tmp2 = tmp
        L(K)=L(K) - 1
        NMOVE=NMOVE + 1 
    END DO
END DO
!========================= STRAIGHT INSERTION =====
DO I=N-2,1,-1
    IF  (A(I + 1) .LT. A(I)) THEN
        HOLD = A(I); tmp=S(I); !temp=R(I)
        J=I
        DO WHILE  (A(J + 1) .LT. HOLD)
            A(J) = A(J + 1); S(J)=S(J+1)
            J = J + 1 
        END DO
        A(J) = HOLD; S(J) = tmp
    ENDIF
END DO

END SUBROUTINE FLASH2 


SUBROUTINE writeVtk_src ( Npoint, numfun,  NcellsX, NcellsY, NcellsZ, new_nodesX, new_nodesY,  new_nodesZ,  &
                            sdx, sdy,  delta, files)
! Fortran code created by J.Sochor   ( https://github.com/JNSresearcher )
! e-mail: JNSresearcher@gmail.com
! This program generates a .vtk file for sources (vector, moving or stationary)

USE m_vxc2data
IMPLICIT NONE

INTEGER,      INTENT(IN):: Npoint, numfun,  sdx,sdy
INTEGER,      INTENT(IN):: NcellsX, NcellsY, NcellsZ, &
                            new_nodesX(NcellsX),new_nodesY(NcellsY),new_nodesZ(NcellsZ) 
REAL(8),      INTENT(IN):: delta(3)
CHARACTER(*), INTENT(IN):: files
CHARACTER (LEN = 1)  ci
CHARACTER (LEN = 4)  ch_sd
CHARACTER (LEN = 30) buf1, buf2
INTEGER:: i,j,k, l, m, n,numcells, nnX, nnY, nnZ, nij,  ios, n_new 
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
        sxm = REAL(i,8)*delta(1) - 1.5d0*delta(1)!/2.d0*
        sym = REAL(j,8)*delta(2) - 1.5d0*delta(2)!/2.d0*
        szm = REAL(k,8)*delta(3) - 1.5d0*delta(3)!/2.d0*
        WRITE(n_new ) sxm, sym, szm
        ! point 1 x = i+1; y=j ; z=k
        sxm = REAL(i+1,8)*delta(1) - 1.5d0*delta(1)!/2.d0
        sym = REAL(j,8)  *delta(2) - 1.5d0*delta(2)!/2.d0
        szm = REAL(k,8)  *delta(3) - 1.5d0*delta(3)!/2.d0
        WRITE(n_new) sxm, sym, szm
        ! point 2 x = i; y=j+1 ; z=k
        sxm = REAL(i,8)  *delta(1) - 1.5d0*delta(1)!/2.d0
        sym = REAL(j+1,8)*delta(2) - 1.5d0*delta(2)!/2.d0 
        szm = REAL(k,8)  *delta(3) - 1.5d0*delta(3)!/2.d0
        WRITE(n_new) sxm, sym, szm
        ! point 3 x = i+1; y=j+1 ; z=k
        sxm = REAL(i+1,8)*delta(1) - 1.5d0*delta(1)!/2.d0
        sym = REAL(j+1,8)*delta(2) - 1.5d0*delta(2)!/2.d0
        szm = REAL(k,8)  *delta(3) - 1.5d0*delta(3)!/2.d0
        WRITE(n_new) sxm, sym, szm
        ! point 4 x = i; y=j ; z=k+1
        sxm = REAL(i,8)*delta(1)   - 1.5d0*delta(1)!/2.d0
        sym = REAL(j,8)*delta(2)   - 1.5d0*delta(2)!/2.d0
        szm = REAL(k+1,8)*delta(3) - 1.5d0*delta(3)!/2.d0
        WRITE(n_new ) sxm, sym, szm
        ! point 5 x = i+1; y=j ; z=k+1
        sxm = REAL(i+1,8)*delta(1) - 1.5d0*delta(1)!/2.d0
        sym = REAL(j,8)  *delta(2) - 1.5d0*delta(2)!/2.d0
        szm = REAL(k+1,8)*delta(3) - 1.5d0*delta(3)!/2.d0
        WRITE(n_new ) sxm, sym, szm
        ! point 6 x = i; y=j+1 ; z=k+1
        sxm = REAL(i,8)  *delta(1) - 1.5d0*delta(1)!/2.d0
        sym = REAL(j+1,8)*delta(2) - 1.5d0*delta(2)!/2.d0
        szm = REAL(k+1,8)*delta(3) - 1.5d0*delta(3)!/2.d0
        WRITE(n_new) sxm, sym, szm
        ! point 7 x = i+1; y=j+1 ; z=k+1
        sxm = REAL(i+1,8)*delta(1) - 1.5d0*delta(1)!/2.d0
        sym = REAL(j+1,8)*delta(2) - 1.5d0*delta(2)!/2.d0
        szm = REAL(k+1,8)*delta(3) - 1.5d0*delta(3)!/2.d0
        WRITE(n_new) sxm, sym, szm
    END SUBROUTINE write_coord
END SUBROUTINE writeVtk_src

subroutine writeVtk_Fe (Npoint, sdx, sdy, sdz, Cells,  delta, nsub_glob, geoPHYS_Fe, files)
implicit none
integer,      intent(IN):: Npoint, sdx,sdy,sdz, Cells, nsub_glob
real(8),      intent(IN):: delta(3)
integer(1),   intent(IN):: geoPHYS_Fe(sdx,sdy,sdz)
character(*), intent(IN):: files
character (len = 1)  ci
character (len = 4)  ch_sd
character (len = 50) buf1
integer:: i,j,k,  ios, n_new
integer(1):: temp_geoPHIS_Fe(sdx,sdy,sdz)

    temp_geoPHIS_Fe = 0_1
    do k = 1,sdz;  do j = 1,sdy;  do i = 1,sdx
        if (geoPHYS_Fe(i,j,k) /=  int(nsub_glob, 1) ) temp_geoPHIS_Fe(i,j,k) = 1_1
    enddo; enddo; enddo
    write(ch_sd,'(i4)')Npoint
    ios=0
    open(newunit=n_new,file=trim(files)//'/fe_'//trim(adjustl(ch_sd))//'.vtk',form="UNFORMATTED",  &
        access="STREAM",  convert="big_endian", iostat=ios)  
    if (ios/=0) then
        print *,'fatal error! Could not open the outfile result.';   stop
    end if
    write (n_new) "# vtk DataFile Version 3.0"//new_line(ci)//"Voxel Fe Binary Data"//new_line(ci)//"BINARY"//new_line(ci) 
    write(n_new)"DATASET STRUCTURED_POINTS"//new_line(ci)
	write(buf1,'(i8," ",i8," ",i8)') sdx, sdy, sdz 
	write(n_new) "DIMENSIONS "//trim(adjustl(buf1))
    write(n_new) new_line(ci)
	write(buf1,'(e10.3," ",e10.3," ",e10.3)') 0.0, 0.0, 0.0
	write(n_new) "ORIGIN "//trim(adjustl(buf1))
    write(n_new) new_line(ci)
	write(buf1,'(e10.3," ",e10.3," ",e10.3)') delta(1), delta(2), delta(3)
    write(n_new) "SPACING "//trim(adjustl(buf1))
    write(n_new) new_line(ci)
    write(buf1,'(i8)')   (sdx) * (sdy) * (sdz)
    write(n_new) "POINT_DATA "//trim(adjustl(buf1)) 
    write(n_new) new_line(ci)
    write(n_new) "SCALARS Fe_volume unsigned_char 1" // new_line(ci)
    write(n_new) "LOOKUP_TABLE default" // new_line(ci)
    write(n_new) temp_geoPHIS_Fe
    close(n_new)
end subroutine writeVtk_Fe

subroutine writeVtk (Npoint, sdx, sdy, sdz, Cells, cells2br, delta, iv, files)
implicit none
integer,      intent(IN):: Npoint, sdx,sdy,sdz, Cells, cells2br(Cells,3)
real(8),      intent(IN):: iv(*), delta(3)
character(*), intent(IN):: files
character (len = 1)  ci
character (len = 4)  ch_sd
character (len = 50) buf1
integer:: i,j,k, nfx,  ios, n_new
real(8):: sxm,sym,szm
integer:: ip, jp, kp

    write(ch_sd,'(i4)')Npoint
    ios=0
    open(newunit=n_new,file=trim(files)//'/field_'//trim(adjustl(ch_sd))//'.vtk',form="UNFORMATTED",  &
        access="STREAM",  convert="big_endian", iostat=ios)  
    if (ios/=0) then
        print *,'fatal error! Could not open the outfile result.';   stop
    end if
    write (n_new) "# vtk DataFile Version 3.0"//new_line(ci)//"out data result"//new_line(ci)//"BINARY"//new_line(ci)
    write(buf1,'(i8," ",i8," ",i8)') sdx, sdy, sdz  
    write(n_new)"DATASET STRUCTURED_GRID"//new_line(ci)//"DIMENSIONS "//trim(adjustl(buf1))//new_line(ci)
    write(buf1,'(i8)')   (sdx) * (sdy) * (sdz)
    write(n_new) "POINTS "//trim(adjustl(buf1))//" float"//new_line(ci)
    do k=1,sdz
        szm = real(k,8)*delta(3)  - delta(3) !/2.d0
        do j=1,sdy
            sym = real(j,8)*delta(2) - delta(2) !/2.d0
            do i=1,sdx
                sxm = real(i,8)*delta(1) - delta(1) !/2.d0
                write(n_new) real(sxm,4), real(sym,4), real(szm,4)
            enddo
        enddo
    enddo
    write(n_new)  new_line(ci)
    write(n_new) "POINT_DATA "//trim(adjustl(buf1))//new_line(ci)
    write(n_new) "VECTORS "//"Vector_field_B"//" float"//new_line(ci)
    nfx=0
    do k = 1,sdz;  do j = 1,sdy;  do i = 1,sdx
        nfx=nfx+1
        ip=cells2br(nfx,1)
        jp=cells2br(nfx,2)
        kp=cells2br(nfx,3)
        sxm=0.d0;sym=0.d0;szm=0.d0;
        if ( ip /= 0 ) then
            sxm =   iv(ip)
        else
            sxm = 0.d0
        endif 
        if ( jp /= 0 ) then
            sym =   iv(jp )
        else
            sym = 0.d0
        endif 
        if ( kp /= 0) then
            szm =   iv(kp)
        else
            szm = 0.d0
        endif 
        write(n_new) real( sxm,4 ), real( sym,4  ), real(szm,4)
    end do;  end do;  end do
    write(n_new)  new_line(ci)
    close(n_new)
end subroutine writeVtk
