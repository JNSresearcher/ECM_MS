! This program converts "VoxCad" format for calculating 3D fields.
! Fortran codes created by J.Sochor   ( https://github.com/JNSresearcher )
! e-mail: JNSresearcher@gmail.com

SUBROUTINE vxc2data ( delta,    dt,        Time,     dtt,                   &
                      sdx,    sdy,      sdz,       nsub,     nsub_air, nsub_glob,  & 
                      numfun, numMech,  numVenv,   &
                      solv,   files,    tol, omega, itmax, tolmag, itmag, omegamag  ) !

USE, INTRINSIC :: ISO_C_BINDING

USE m_vxc2data
USE m_fparser

IMPLICIT NONE

character(len=:), allocatable :: str_in

INTERFACE
    SUBROUTINE c_routine (  str_in  ) BIND(C, name="c_routine")
        IMPORT 
        CHARACTER(KIND=C_CHAR), DIMENSION(*), INTENT(IN) :: str_in
        END SUBROUTINE c_routine
END INTERFACE

REAL(8), INTENT(OUT) :: delta(3),  &         ! Grid spacing along X, Y and Z
                        dt,        &         ! time step
                        Time,      &         ! stop time
                        dtt                  ! jump duration
INTEGER, INTENT(OUT) :: sdx,  sdy,  sdz,  &  ! Number of cells along X,  Y and  Z
                        nsub,             &  ! number of physical domains
                        nsub_air,         &  ! number of environment domains
                        nsub_glob            ! total number of domains
INTEGER, INTENT(OUT) :: numfun,       &      ! number of functions for calculating external sources
                        numMech,      &      ! number of functions for calculating movements of external sources
                        numVenv

CHARACTER(LEN=3), INTENT(OUT) :: solv        ! character string corresponding to the name of methods: 'BCG' or 'SOR'
CHARACTER(16),    INTENT(OUT) :: files       ! name for output files
REAL(8),          INTENT(OUT) :: tol, tolmag       ! convergence criterion, 
REAL(8),          INTENT(OUT) :: omega, omegamag   ! coeff relax SOR and for nonlin iterations
INTEGER,          INTENT(OUT) :: itmax, itmag      ! maximum number of iterations, 

!-------------------!
! working variables !
!-------------------!
REAL(8) delta0
CHARACTER (LEN=6), ALLOCATABLE ::  nameFun(:), nameschFun(:), nameVmech(:), nameVenv(:)
INTEGER(1),        ALLOCATABLE::   v(:), v_e(:)
INTEGER,           ALLOCATABLE::   temp(:)

TYPE tmp
    INTEGER jm                   ! function number
    INTEGER im                   ! cell number
    CHARACTER (len = 1)  ch      ! function type: X Y Z or 0
    TYPE (tmp),pointer ::prec
END TYPE tmp
TYPE(tmp),pointer :: sp_F,next_F 

TYPE t_coil
    INTEGER, ALLOCATABLE  :: nods(:), irow(:),jcol(:)
    INTEGER num_coil, num_fun, siz_nods, nnz 
    real(8), ALLOCATABLE  :: valp(:), Uc(:),Sc(:)
    REAL(8) S, parall
    character(len=:), allocatable :: name_coil, name_fun
END TYPE t_coil
TYPE (t_coil), ALLOCATABLE :: dom_coil(:) 
INTEGER sizdom_coil

TYPE t_coil_sec
    INTEGER, ALLOCATABLE  :: nods(:)
    INTEGER num_coil,  num_coil_sec,  siz_nods, orient  
    character(len=:), allocatable :: name_coil, name_sec, connect   ! name_sec - omly for inform
END TYPE t_coil_sec
TYPE (t_coil_sec), ALLOCATABLE :: dom_coil_sec(:)  
INTEGER sizdom_coil_sec

CHARACTER (len = 10) ch_e
CHARACTER (len = 60) words(20) !
CHARACTER (len = 74) letter
CHARACTER (len=:),ALLOCATABLE :: st, ch

INTEGER m1, m2, kdz, Cells, nn, i,j,k,  L, m, n, k0, k1, k2, ios 
INTEGER num_nodx, num_nody, num_nodz, num_nodV                     ! sign of the presence of a vector field
INTEGER lst, neww, np, kp, ier, idev, iter, itmax_bcg, ver_python 
LOGICAL uncompress
CHARACTER (LEN=4)::  varconst(12)
REAL(8)          ::  valconst(12), tol_bcg

REAL(8)  numeric
EXTERNAL numeric

letter='123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz'

! init data DEFAULT for solver  omegamag=0.7
solv='SOR'; tol=0.5d-3; omega=1.92d0;  itmax=10000; itmag=200; tolmag=1.0d-3; omegamag=0.d0
files='out';  ver_python = 0
siz_non_lin_MU = 0

OPEN(1,file='in.vxc')
! input structure
k0=0;k1=0;k2=0; kp=0;  

numfun = 0;  numMech=0; numVenv=0; sizdom_coil = 0; sizdom_coil_sec = 0;

m=0 !  m - it is counter CDATA - not USE in DO WHILE
convert_strings0: &
DO WHILE ( 1==1)
    idev = 1
    CALL readline(st,ch,ier,idev)             !============================= 1 idev= 1

    IF (ier == -1) EXIT convert_strings0
    
    IF (st(1:6)=='</VXC>' ) EXIT
    lst = len_trim(st)
    
    m1=index(st,'<Lattice_Dim>')
    IF (m1>=1) THEN
        m2=index(st,'</Lattice_Dim>')
        ch_e=trim(st(m1+13:m2-1))
        delta0=numeric(ch_e)
        PRINT '( a,g10.3,$ )', 'delta0=',delta0
    ENDIF
    m1=index(st,'<X_Dim_Adj>')
    IF (m1>=1) THEN
        m2=index(st,'</X_Dim_Adj>')
        ch_e=trim(st(m1+11:m2-1))
        delta(1)=numeric(ch_e) * delta0
        WRITE(*,'( a,g10.3,$ )') ' deltaX=',REAL(delta(1))
    ENDIF
    m1=index(st,'<Y_Dim_Adj>')
    IF (m1>=1) THEN
        m2=index(st,'</Y_Dim_Adj>')
        ch_e=trim(st(m1+11:m2-1))
        delta(2)=numeric(ch_e) * delta0
        WRITE(*,'( a,g10.3,$ )') ' deltaY=',REAL(delta(2))
    ENDIF
    m1=index(st,'<Z_Dim_Adj>')
    IF (m1>=1) THEN
        m2=index(st,'</Z_Dim_Adj>')
        ch_e=trim(st(m1+11:m2-1))
        delta(3)=numeric(ch_e) * delta0
        WRITE(*,'( a,g10.3 )') ' deltaZ=',REAL(delta(3))
    ENDIF
    m1=index(st,'Material ID')
    IF (m1>=1) THEN
        kp=kp+1
    ENDIF

    k0=index(st,'<Name>')
    IF (k0>0 ) THEN ! 
        n=index(st,'</')
        st=st(k0+6:n-1)
        i=1
        DO
            k = index( st(i:),'=')
            IF (k/=0) THEN
                i=i+k-1; st(i:i) = ' '
            ELSE
                EXIT
            ENDIF
        ENDDO
        
        CALL Upp (st, st)  
        CALL string2words(trim(st), words, neww)
        
        IF (neww > 20 ) THEN
            PRINT*,'error: numbers words >20  newws=: ', neww,' correct dim words(20) string 60 in vxc2data.f90 '
            STOP
        ENDIF
        !  calc numbers L-domain
        DO i=2,neww
            IF (index( words(i),'MU') /=0 ) THEN

                IF (i+2 > neww) CYCLE ! it is last word
                
                DO j=i+2, neww
                    if ( index( words(j),'VE') /=0) THEN
                        IF ( index( words(j),'VEX') /=0 ) THEN 
                            numVenv = numVenv + 1
                        ELSEIF ( index( words(j),'VEY') /=0 ) THEN 
                            numVenv = numVenv + 1
                        ELSEIF ( index( words(j),'VEZ') /=0 ) THEN 
                            numVenv = numVenv + 1
                        ENDIF
                    endif
                
                    IF ( index( words(j),'FILE') /=0) THEN
                        siz_non_lin_MU = siz_non_lin_MU + 1
                    endif
                
                    IF (  index( words(j),'COIL') /=0  ) THEN
                        sizdom_coil_sec = sizdom_coil_sec + 1
                        if ( index( words(j+3),'XP') /=0  .or. index( words(j+3),'XM') /=0 .or. &
                             index( words(j+3),'YP') /=0  .or. index( words(j+3),'YM') /=0 .or. &
                             index( words(j+3),'ZP') /=0  .or. index( words(j+3),'ZM') /=0 ) then
                             sizdom_coil_sec = sizdom_coil_sec + 1
                        endif
                    endif
                
                    IF ( index( words(j),'SRC') /=0) THEN
                        sizdom_coil = sizdom_coil + 1
                        numfun = numfun + 1
                        DO k=j+2, neww
                            IF     ( index( words(k),'VSY') /=0 ) THEN 
                                numMech = numMech + 1
                            ELSEIF ( index( words(k),'VSX') /=0 ) THEN 
                                numMech = numMech + 1
                            ELSEIF ( index( words(k),'VSZ') /=0) THEN
                                numMech = numMech + 1
                            ELSE
                            ENDIF
                        ENDDO
                    ENDIF  
                enddo
                
            ENDIF
            
            IF (trim(words(i)) == 'TRAN') THEN
                DO j=i+1, neww-1,2
                    IF    ( index( words(j),'STOP') /=0) THEN
                        ch_e=trim(words(j+1) ) !(k+1:))
                        Time=numeric(ch_e)
                        WRITE(*,'( a,g10.3,$ )') ' Time=',REAL(Time)
                    ELSEIF ( index(words(j),'STEP') /=0) THEN
                        ch_e=trim(words(j+1))
                        DT=numeric(ch_e)
                        WRITE(*,'( a,g10.3,$ )') ' DT=',REAL(DT)
                    ELSEIF ( index(words(j),'JUMP') /=0) THEN
                        ch_e=trim(words(j+1))
                        DTT=numeric(ch_e)
                        WRITE(*,'( a,g10.3,$ )') ' DTT=',REAL(DTT)
                    ENDIF
                ENDDO
            ENDIF
            
            IF (trim(words(i)) == 'SOLVER') THEN 
                DO j=i+1,neww-1
                    IF ( index(words(j),'TOLSOR') /=0) THEN
                        ch_e=trim(words(j+1))
                        tol=numeric(ch_e)
                        WRITE(*,'( a,g10.3,$ )') ' tolerance=',REAL(tol)
                    ELSEIF ( index(words(j),'ITSOR') /=0) THEN
                        ch_e=trim(words(j+1))
                        itmax= nint(numeric(ch_e))
                        WRITE(*,'( a,i5,$ )') ' itsor=',itmax
                    ELSEIF ( index(words(j),'ITMAG') /=0) THEN
                        ch_e=trim(words(j+1))
                        itmag= nint(numeric(ch_e))
                        WRITE(*,'( a,i5,$ )') ' itmag=',itmag
                    ELSEIF ( index(words(j),'TOLMAG') /=0) THEN
                        ch_e=trim(words(j+1))
                        tolmag= numeric(ch_e)
                        WRITE(*,'( a,g10.3,$ )') ' tolmag=',tolmag
                    ELSEIF ( index(words(j),'OMESOR') /=0) THEN
                        ch_e=trim(words(j+1))
                        omega= numeric(ch_e)
                        WRITE(*,'( a,g10.3,$ )') ' omega=',omega
                    ELSEIF ( index(words(j),'OMEMAG') /=0) THEN
                        ch_e=trim(words(j+1))
                        omegamag= numeric(ch_e)
                        WRITE(*,'( a,g10.3,$ )') ' omegamag=',omegamag
                    ELSEIF ( index(words(j),'SOLV') /=0) THEN
                        solv = trim(words(j+1))
                        WRITE(*,'( a,a,$ )') ' solver=',solv
                    ELSEIF ( index(words(j),'DIR') /=0) THEN
                        files = trim(words(j+1))

                    ELSEIF ( index(words(j),'PYTHON') /=0) THEN
                        ch_e=trim(words(j+1))
                        ver_python = int(numeric(ch_e))
                        WRITE(*,'( a,i2,$ )') ' version python= ',ver_python
                    ENDIF
                ENDDO
            ENDIF
        ENDDO
    ENDIF !Name
    
    IF (index(st,'Structure Compression="ZLIB"') >= 1) THEN
        uncompress = .true.
    ELSEIF (index(st,'Structure Compression="ASCII_READABLE"') >= 1) THEN
        uncompress = .false.
    ENDIF
    m1=index(st,'<X_Voxels>')
    IF (m1>=1) THEN
        m2=index(st,'</X_Voxels>')
        ch_e=st(m1+10:m2-1)
        READ(ch_e,'(i3)')sdx
        WRITE(*,'( a,i3,$ )') ' sdx=',sdx
    ENDIF
    m1=index(st,'<Y_Voxels>')
    IF (m1>=1) THEN
        m2=index(st,'</Y_Voxels>')
        ch_e=st(m1+10:m2-1)
        READ(ch_e,'(i3)')sdy
        WRITE(*,'( a,i3,$ )') ' sdy=',sdy
    ENDIF
    m1=index(st,'<Z_Voxels>')
    IF (m1>=1 .and. m==0) THEN
        m2=index(st,'</Z_Voxels>')
        ch_e=st(m1+10:m2-1)
        READ(ch_e,'(i3)')sdz
        WRITE(*,'( a,i3 )') ' sdz=',sdz
          
        kdz = sdx*sdy
        ALLOCATE(v(sdx*sdy*sdz), v_e(sdx*sdy*sdz), source=0_1)

        READ(1,'(a)',end=91,iostat=ios) ch_e ! skip line  <Data>
        IF ( uncompress ) THEN
            ! START UNCOMPRESS
            PRINT*,'start uncompress'
            if (ver_python == 0) then
                idev=1
                m=0
                DO j=1,sdz
                    CALL readline(st,ch, ier, idev)             !================== 2 dev=1
                    k=index(st,'CDATA')
                    k2=index(st,']]></')
                    IF (k2 == 0) THEN
                        STOP ' Error: long string CDATA'
                    ENDIF
                    k1=k+6
                    ch=st(k1:k2-1) 
                    str_in = ch//C_NULL_CHAR 
                    CALL c_routine (str_in)
                    open(2, file='tmp.txt')
                    read (2, * ) v(m+1 : m+kdz)
                    ! print '( *(i4) )',  v(m+1 : m+kdz)
                    m = m + kdz
                    close (2)
                enddo
                call execute_command_line ('del tmp.txt', exitstat=i)
                PRINT *, "delete files tmp.txt, status: ", i 
            else
                ios=0
                open(2,file='compress.txt', iostat=ios) 
                idev=1
                DO j=1,sdz
                    CALL readline(st,ch, ier, idev)             !================== 2 dev=1 
                    k=index(st,'CDATA')
                    k2=index(st,']]></')
                    IF (k2 == 0) THEN
                        STOP ' Error: long string CDATA'
                    ENDIF
                    k1=k+6
                    ch=st(k1:k2-1)                        
                    lst = len_trim(ch)
                    WRITE(2, '( a )') ch(1:lst)
                    m = m + lst
                ENDDO
                CLOSE(2)
            
                if (ver_python == 2) then
                    call execute_command_line ('python uncompress_zlib2.py', exitstat=i)
                elseif (ver_python == 3) then
                    call execute_command_line ('python uncompress_zlib3.py', exitstat=i)
                else
                    PRINT *, "not ver_python  required 2 or 3 :", ver_python
                    stop
                endif
                PRINT *, "uncompress.txt created, status: ", i 
                open(2,file='uncompress.txt', iostat=ios)
                st=''; m=0
                idev = 2
                DO j=1,sdz
                    CALL readline(st,ch, ier, idev)               !=======     3 idev= 2
                    lst = len_trim(st)
                    lst=lst-1
                    DO i=1,lst
                        v(m+i) = int( index(letter, st(i:i)), 1  )
                    ENDDO
                    m = m + lst
                ENDDO
                CLOSE(2)                
                call execute_command_line ('del compress.txt; del uncompress.txt ', exitstat=i)
                
                PRINT *, "delete files uncompress.txt and compress.txt, status: ", i 
            endif
        ELSE
            idev=1
            DO j=1,sdz
                CALL readline(st,ch,ier, idev)              !--------------4  idev=1
                k=index(st,'CDATA')
                k2=index(st,']]></')
                IF (k2 == 0)   STOP ' Error: long string CDATA'
                k1=k+6
                ch=st(k1:k2-1)
                lst = len_trim(ch)
                DO i=1,lst
                    v(m+i) = int( index(letter, ch(i:i)), 1  )
                ENDDO
                m = m + lst
            ENDDO
        ENDIF
    ENDIF

91 IF (ios/=0) EXIT

END DO convert_strings0


nsub=maxval(v)        ! only the number of domains
Cells = sdx*sdy*sdz

j=0; k=1
DO i=1,Cells
    IF (v(i) == 0) THEN
        j = j + 1                ! new cell
        IF (j == 500000) THEN    ! 3500 2200  300000  
            j=0                  ! reset cells counter
            k = k + 1            ! new domain
        ENDIF
        v(i) = int( (nsub + k), 1 )         ! another one added for environment
    ENDIF
ENDDO
IF (j==0) k = k - 1              ! there was no enlargement of the cells
nsub_air = k

nsub_glob = nsub + nsub_air

!=======================
IF (sizdom_coil_sec /=0 ) THEN
    ALLOCATE ( dom_coil_sec(sizdom_coil_sec) )
else
    stop 'error: not section!'
ENDIF

IF (sizdom_coil /=0 ) THEN
    ALLOCATE ( dom_coil(sizdom_coil) )
    do i=1,sizdom_coil 
        dom_coil(i)%parall = 1.d0
    enddo
else
    stop 'error: not section!'
ENDIF

IF (numfun /=0 ) THEN
    ALLOCATE ( nameFun(numfun), nameschFun(numfun) ) 
    nameFun=''; nameschFun=''
ENDIF

IF ( numMech /=0 ) THEN
    ALLOCATE (  nameVmech(numMech) )
    nameVmech='';
    ALLOCATE (Vmech(numMech) ) 
    DO i=1,numMech
        Vmech(i)%vely = 0.d0
        Vmech(i)%args = 0
        Vmech(i)%nomsch = 0
        Vmech(i)%ex = ' '
        Vmech(i)%eqn = ' '
    ENDDO
ELSE
ENDIF

IF ( numVenv /=0 ) THEN
    ALLOCATE (  nameVenv(numVenv) )
    nameVenv='';
    ALLOCATE (Venv(numVenv), fun_Venv(numVenv )  ) 
    DO i=1,numVenv
        Venv(i)%vely = 0.d0
        Venv(i)%args = 0
        Venv(i)%nomsch = 0
        Venv(i)%ex = ' '
        Venv(i)%eqn = ' '
        fun_Venv(i)%val_Venv(1:3) = 0.d0; 
        fun_Venv(i)%num_Venv(1:3) = 0;  
        fun_Venv(i)%move(1:3) = 0
    ENDDO
    fun_Venv(:)%length(1)   = 0;    fun_Venv(:)%length(2)   = 0;    fun_Venv(:)%length(3)   = 0;   
    fun_Venv(:)%shift(1)    = 0.d0; fun_Venv(:)%shift(2)    = 0.d0; fun_Venv(:)%shift(3)    = 0.d0;
    fun_Venv(:)%Distance(1) = 0.d0; fun_Venv(:)%Distance(2) = 0.d0; fun_Venv(:)%Distance(3) = 0.d0; 
ELSE
    ! PRINT*,  'var Venv not use'
ENDIF

ALLOCATE ( typPHYS(nsub_glob), namePHYS(nsub_glob ) )
 typPHYS=''; namePHYS='';
ALLOCATE ( valPHYS(nsub_glob,5), source=0.d0 )


IF (siz_non_lin_MU /= 0 ) then
    ALLOCATE (non_lin_MU( siz_non_lin_MU) )
endif


! DEFAULT param for ENVIRONMENT
IF (nsub_air /=0) THEN
    DO j=nsub+1, nsub_glob
        typPHYS(j)   = 'R     '
        namePHYS(j)  = 'AIR   '
        valPHYS(j,1) = 1.d0 ! D
    ENDDO
ENDIF

if(numfun /=0) THEN
    ALLOCATE (Fun(numfun), fun_nod(numfun) )
    DO i=1,numfun
        Fun(i)%vely = 0.d0
        Fun(i)%args = 0
        Fun(i)%nomsch = 0
        Fun(i)%ex = ' '
        Fun(i)%eqn = ' '
        fun_nod(i)%val_Vmech(1:3) = 0.d0; 
        fun_nod(i)%num_Vmech(1:3) = 0;  
        fun_nod(i)%move(1:3) = 0
    ENDDO
    fun_nod(:)%length(1)   = 0;    fun_nod(:)%length(2)   = 0;    fun_nod(:)%length(3)   = 0;   
    fun_nod(:)%shift(1)    = 0.d0; fun_nod(:)%shift(2)    = 0.d0; fun_nod(:)%shift(3)    = 0.d0;
    fun_nod(:)%Distance(1) = 0.d0; fun_nod(:)%Distance(2) = 0.d0; fun_nod(:)%Distance(3) = 0.d0; 
ENDIF
!           1       2        3      4       5         6        7       8          9       10     11       12
varconst=['PI  ','E   ','MU0 ','E0  ','DT  ', 'DX  ', 'DY  ', 'DZ  ', 'TIME', 'NX  ','NY  ','NZ  ']

valconst(1) = 3.1415926535897932384626433832795_8 ! 'PI  '
valconst(2) = 0.27182818284590451e+001_8          ! 'E   '
valconst(3) = 0.12566370964050292e-005_8          ! 'MU0 ' 
valconst(4) = 0.88541878176203908e-011_8          ! 'E0  '
valconst(5) = dt        ! 'DT  '
valconst(6) = delta(1)  ! 'DX  '
valconst(7) = delta(2)  ! 'DY  '
valconst(8) = delta(3)  ! 'DZ  '
valconst(9) = Time      ! 'TIME'
valconst(10)= sdx       ! 'NX  '
valconst(11)= sdy       ! 'NY  '
valconst(12)= sdz       ! 'NZ  '

REWIND(1)

!filling arrays
numdom_Fe = [integer(1):: ]
kp=0; ios=0; numfun=0; numMech=0; numVenv=0;  sizdom_coil = 0; sizdom_coil_sec = 0; 
siz_non_lin_MU = 0

idev = 1
convert_strings1: &
DO WHILE ( 1==1)
    CALL readline(st,ch, ier, idev)              !-------------- 5
    IF (st(1:6)=='</VXC>' ) EXIT
    k0=index(st,'<Name>')
    IF ( k0 > 0 ) THEN 
        kp = kp + 1                            ! count all, including lines with tran param etc
        m2=index(st,'</')
        st=st(k0+6:m2-1)
        i=1
        DO
            k = index( st(i:),'=')
            IF (k/=0) THEN
                i=i+k-1; st(i:i) = ' '
            ELSE
                EXIT
            ENDIF
        ENDDO
    
        CALL Upp (st, st)
        CALL string2words(trim(st),words,neww)
        
        DO i=2,neww
            IF ( trim (words(i)) == 'MU' .and. kp <= nsub ) THEN
                valPHYS(kp,1) = 1.d0/evaluate( words(i+1) ) 
                if ( abs(valPHYS(kp,1) - 1.d0) > 1.d-6  )  then
                    numdom_Fe = [numdom_Fe, int(kp, 1) ]
                    DO j=i+2, neww-1
                        IF ( index( words(j),'VE')  /=0) THEN 
                            numVenv = numVenv + 1
                            Venv(numVenv)% nomsch = kp
                            Venv(numVenv)%ex = '-'!ch 
                            nameVenv(numVenv) = trim(words(j+1)) 
                            CALL calcVenv 
                        elseIF ( index( words(j),'FILE_') /=0 ) THEN 
                            siz_non_lin_MU = siz_non_lin_MU + 1
                            non_lin_MU(siz_non_lin_MU)%num_domain = kp
                            non_lin_MU(siz_non_lin_MU)%file_MU = trim(words(j+1))
                            if ( index( words(j),'_BH') /=0 ) THEN 
                                non_lin_MU(siz_non_lin_MU)%tip_non_lin = 'BH'
                            elseif ( index( words(j),'_BM') /=0 ) THEN 
                                non_lin_MU(siz_non_lin_MU)%tip_non_lin = 'BM'
                            else
                                print*, 'error in tip nonlin FE', non_lin_MU(siz_non_lin_MU)%tip_non_lin
                                stop
                            endif
                        endif 
                    enddo
                endif
                typPHYS(kp) = 'R'
                namePHYS(kp) = trim(words(1))
                                  
                IF (i+1 == neww) THEN
                    CYCLE ! it is last word
                ENDIF
                IF ( index( words(i+2),'COIL') /=0 ) THEN 
                    sizdom_coil_sec = sizdom_coil_sec + 1
                    dom_coil_sec( sizdom_coil_sec)%name_coil = trim(words(i+3))
                    dom_coil_sec( sizdom_coil_sec)%num_coil_sec = kp
                    dom_coil_sec( sizdom_coil_sec)%name_sec = trim(words(1))
                    IF     ( index( words(i+4),'JOIN') /=0 ) THEN 
                        dom_coil_sec( sizdom_coil_sec)%connect = trim(words(i+5))
                        if ( index( words(i+5),'XP') /=0 .or. index( words(i+5),'XM') /=0  .or. &
                             index( words(i+5),'YP') /=0 .or. index( words(i+5),'YM') /=0  .or. &
                             index( words(i+5),'ZP') /=0 .or. index( words(i+5),'ZM') /=0 ) then
                            sizdom_coil_sec = sizdom_coil_sec + 1
                            dom_coil_sec( sizdom_coil_sec)%num_coil_sec = kp  
                            dom_coil_sec( sizdom_coil_sec)%name_coil = trim(words(i+3))
                            dom_coil_sec( sizdom_coil_sec)%name_sec = trim(words(1))//'_D'
                            if ( index( words(i+5),'XP') /=0 )dom_coil_sec( sizdom_coil_sec)%connect = 'XM'
                            if ( index( words(i+5),'XM') /=0 )dom_coil_sec( sizdom_coil_sec)%connect = 'XP'
                            if ( index( words(i+5),'YP') /=0 )dom_coil_sec( sizdom_coil_sec)%connect = 'YM'
                            if ( index( words(i+5),'YM') /=0 )dom_coil_sec( sizdom_coil_sec)%connect = 'YP'
                            if ( index( words(i+5),'ZP') /=0 )dom_coil_sec( sizdom_coil_sec)%connect = 'ZM'
                            if ( index( words(i+5),'ZM') /=0 )dom_coil_sec( sizdom_coil_sec)%connect = 'ZP'
                        endif
                    else
                        stop 'not find key JOIN'
                    ENDIF
                endif
                
                IF ( index( words(i+2),'SRC') /=0 ) THEN 
                    sizdom_coil = sizdom_coil + 1
                    dom_coil(sizdom_coil)%name_coil = trim(words(1))
                    dom_coil(sizdom_coil)%name_fun = trim(words(i+3))
                    dom_coil(sizdom_coil)%num_coil = kp
                    DO j=i+2, neww-1
                        IF     ( index( words(j),'PAR') /=0 ) THEN 
                            dom_coil(sizdom_coil)%parall = evaluate( words(j+1) )   !trim(words(j+1)) 
                            exit
                        ENDIF
                    ENDDO
                    DO j=i+2, neww-1
                        IF     ( index( words(j),'SRC') /=0 ) THEN 
                            numfun = numfun + 1
                            Fun(numfun)% nomsch = kp
                            Fun(numfun)%ex = '-'!ch 
                            nameFun(numfun) = trim(words(j+1)) 
                            nameschFun(numfun) = trim(words(1)) 
                            CALL calcVmech 
                        ENDIF
                    ENDDO
                ENDIF 
            
            ELSEIF ( index( words(i),'FUNC') /=0  )  THEN  !   Key word FUNC
                DO L=1,numfun
                    If (  trim(words(i+1)) == trim(nameFun(L)) ) THEN ! 
                        j=i+1  
                        Fun(L)%eqn = trim(words(j+1))
        
                        DO WHILE (0 ==0  )  ! j - eto namex
                            Fun(L)%args = Fun(L)%args + 1
                            j=j+2 ! next x
                            IF (j+1 > neww) EXIT
                        ENDDO
                        Fun(L)%args = Fun(L)%args - 1
                        
                        ALLOCATE (Fun(L)%velx(1:Fun(L)%args),  Fun(L)%namex(1:Fun(L)%args) )
    
                        k=0; j=i+1 
                        DO WHILE (0 ==0  )  ! j - it is namex
                            k=k+1
                            IF ( k > Fun(L)%args ) EXIT
                            Fun(L)%namex(k) = words(j+2)(1:8)
                            Fun(L)%velx(k) = evaluate(words(j + 3))
                            j=j+2             ! next x
                            IF (j+1 > neww) EXIT
                        ENDDO
                    ENDIF 
                ENDDO ! numfun
            
                if (numMech /=0) then
                    DO L=1,numMech
                        If (  trim(words(i+1)) == trim(nameVmech(L)) ) THEN ! 
                            j=i+1  
                            Vmech(L)%eqn = trim(words(j+1))
            
                            DO WHILE ( 0 ==0  )  ! j - it is  namex
                                Vmech(L)%args = Vmech(L)%args + 1
                                j=j+2            ! next x
                                IF (j+1 > neww) EXIT
                            ENDDO
                            Vmech(L)%args = Vmech(L)%args - 1
                            
                            ALLOCATE (Vmech(L)%velx(1:Vmech(L)%args),  Vmech(L)%namex(1:Vmech(L)%args) )
                            
                            k=0; j=i+1 
                            DO WHILE (0==0)         ! j - it is  namex
                                k=k+1
                                IF ( k > Vmech(L)%args ) EXIT
                                Vmech(L)%namex(k) = words(j+2)(1:8)
                                Vmech(L)%velx(k) = evaluate(words(j + 3))
                                j=j+2                ! next x
                                IF (j+1 > neww) EXIT
                            ENDDO
                        ENDIF 
                    ENDDO ! numMech
                endif
                
                if (numVenv /=0) then
                    DO L=1,numVenv
                        If (  trim(words(i+1)) == trim(nameVenv(L)) ) THEN ! 
                            j=i+1  
                            Venv(L)%eqn = trim(words(j+1))
                            DO WHILE ( 0 ==0  )  ! j - it is  namex
                                Venv(L)%args = Venv(L)%args + 1
                                j=j+2            ! next x
                                IF (j+1 > neww) EXIT
                            ENDDO
                            Venv(L)%args = Venv(L)%args - 1
                            
                            ALLOCATE (Venv(L)%velx(1:Venv(L)%args),  Venv(L)%namex(1:Venv(L)%args) )
                            
                            k=0; j=i+1 
                            DO WHILE (0==0)         ! j - it is  namex
                                k=k+1
                                IF ( k > Venv(L)%args ) EXIT
                                Venv(L)%namex(k) = words(j+2)(1:8)
                                Venv(L)%velx(k) = evaluate(words(j + 3))
                                j=j+2                ! next x
                                IF (j+1 > neww) EXIT
                            ENDDO
                        ENDIF 
                    ENDDO ! numVenv
                endif
            ELSE

            ENDIF  ! key word FUNC
        ENDDO  ! words
    ENDIF ! <Name>

92  IF (ios/=0) EXIT

END DO convert_strings1



if (sizdom_coil == 0 ) stop 'error: not coil'
PRINT '( "Cells= ",i9,  " numfun= ",i3, " numMech= ",i3, " numVenv= ",i3, " nsub= ", i3, " nsub_air= ", i3 ," nsub+nsub_air=", i3 )', &
            Cells ,       numfun,        numMech,       numVenv,             nsub,          nsub_air,          nsub_glob
CLOSE(1)

ALLOCATE (geoPHYS(sdx,sdy,sdz), source=0_1  )
geoPHYS = reshape(v, shape =(/sdx,sdy,sdz/) ,  order = (/ 1, 2, 3/) )
sizdom_Fe = size(numdom_Fe)    ! ordinal numbers of cells in the ferrum region

do L = 1, sizdom_coil
    dom_coil(L)%S=0.d0
    n = count(v == dom_coil(L)%num_coil)
    dom_coil(L)%siz_nods = n 
    if (n == 0) then
        print*,'not nides for coil ', dom_coil(L)%name_coil
        stop
    else
        allocate( dom_coil(L)%nods(n), source = 0 )
        allocate( dom_coil(L)%Uc(n), dom_coil(L)%Sc(n), source = 0.0d0 )
        allocate (dom_coil(L)%irow(n+1), source = 1 )
    endif
enddo

! sections
do L = 1, sizdom_coil_sec
    m = count(v == dom_coil_sec(L)%num_coil_sec)
    dom_coil_sec(L)%siz_nods = m 
    if (m == 0) then
        print*,'not nodes for coil_sec ', dom_coil_sec(L)%name_sec
        stop
    else
        allocate( dom_coil_sec(L)%nods(m), source = 0 )
    endif
enddo

! find  nodes of coils
do L = 1, sizdom_coil
    k=0;
    DO j=1,Cells
        IF (v(j) == dom_coil(L)%num_coil )  THEN
            k=k+1
            dom_coil(L)%nods(k) = j 
        ENDIF
    ENDDO
enddo

! find  nodes of sections
do L = 1, sizdom_coil_sec
    n=0; 
    DO j=1,Cells
        IF (v(j) == dom_coil_sec(L)%num_coil_sec )  THEN
            n=n+1
            dom_coil_sec(L)%nods(n) = j 
        ENDIF
    ENDDO
enddo

! save num coils in sections 
do L = 1, sizdom_coil_sec
    dom_coil_sec(L)%num_coil = 0
    do m = 1,sizdom_coil
        if (dom_coil_sec(L)%name_coil == dom_coil(m)%name_coil ) then
            dom_coil_sec(L)%num_coil = m
            exit
        endif
    enddo
    if (dom_coil_sec(L)%num_coil == 0) then
        print*,'error not find coil',dom_coil_sec(L)%name_coil ,' for section', dom_coil_sec(L)%name_sec
        stop
    endif
enddo

CONNECT: BLOCK
integer sides(6),ind(1),nim, nip, njm, njp, nkm, nkp

kdz=sdx*sdy
do L = 1, sizdom_coil_sec
  m1:do m = 1, dom_coil_sec(L)%siz_nods
        n = dom_coil_sec(L)%nods(m)
        k = dom_coil_sec(L)%num_coil 
        nim = n-1; njm = n-sdx;  nkm = n-kdz; 
        nip = n+1; njp = n+sdx;  nkp = n+kdz;
        sides = [nim, nip, njm, njp, nkm, nkp]
        j=0
        do i = 1,6 
            ind = FINDLOC(dom_coil(k)%nods, sides(i))
            If ( ind(1) /= 0) then 
                if (dom_coil_sec(L)%connect == 'XM' .and. i==1) then
                    j = ind(1)
                    exit 
                elseif (dom_coil_sec(L)%connect == 'XP' .and. i==1  .or. dom_coil_sec(L)%connect == 'XM' .and. i==2 ) then
                    cycle m1
                elseif (dom_coil_sec(L)%connect == 'XP' .and. i==2) then
                    j = ind(1) 
                    exit 
                elseif (dom_coil_sec(L)%connect == 'YM' .and. i==3) then
                    j = ind(1) 
                    exit 
                elseif (dom_coil_sec(L)%connect == 'YP' .and. i==3  .or. dom_coil_sec(L)%connect == 'YM' .and. i==4 ) then
                    cycle m1
                elseif (dom_coil_sec(L)%connect == 'YP' .and. i==4) then
                    j = ind(1)
                    exit
                elseif (dom_coil_sec(L)%connect == 'ZM' .and. i==5) then
                    j = ind(1) 
                    exit 
                elseif (dom_coil_sec(L)%connect == 'ZP' .and. i==5  .or. dom_coil_sec(L)%connect == 'ZM' .and. i==6 ) then
                    cycle m1
                elseif (dom_coil_sec(L)%connect == 'ZP' .and. i==6) then
                    j = ind(1) ! local node
                    exit
                else
                    j = ind(1)
                    exit
                endif
            endif
        enddo
        if (j == 0) then
            print*, 'error not find connect'
            stop 
        endif
        dom_coil_sec(L)%orient = i
        if     (dom_coil_sec(L)%connect == 'N1') then
            dom_coil(k)%Sc(j) = 1.d3
        elseif (dom_coil_sec(L)%connect == 'N2') then
            dom_coil(k)%Sc(j) = -1.d3
        elseif (dom_coil_sec(L)%connect == 'XP' .or. dom_coil_sec(L)%connect == 'XM' .or. &
                dom_coil_sec(L)%connect == 'YP' .or. dom_coil_sec(L)%connect == 'YM' .or. &
                dom_coil_sec(L)%connect == 'ZP' .or. dom_coil_sec(L)%connect == 'ZM' ) then
            if ( index(dom_coil_sec(L)%name_sec, '_D') /=0) then
                dom_coil(k)%Sc(j) = -1.d3
            else
                dom_coil(k)%Sc(j) =  1.d3
            endif 
        else
            print*, 'for section ', dom_coil_sec(L)%name_sec, ' not find N1 or N2'
        endif
    enddo m1
enddo

END BLOCK CONNECT

DEALLOCATE (v)

call gen_sparse_coil_matrix

tol_bcg =5.0d-3 ; itmax_bcg = 10000
do L = 1, sizdom_coil
    CALL sprsBCGstabwr (dom_coil(L)%valp, dom_coil(L)%irow, dom_coil(L)%jcol, dom_coil(L)%siz_nods, &
                                dom_coil(L)%Sc, dom_coil(L)%Uc,  tol_bcg, itmax_bcg, iter)
    IF (iter > itmax_bcg) &
    PRINT*,  'name coil=', dom_coil(L)%name_coil, "Possibly the coil cross-section at 'the bends' is too small"
enddo

IF (numfun /=0 ) THEN
    do i = 1, numfun
         do j = 1, sizdom_coil
            if (dom_coil(j)%name_fun == nameFun(i))  then
                dom_coil(j)%num_fun = i
            endif
        enddo
    enddo

    n=0; m1=0;
    DO k=1,sdz
        DO j=1,sdy
            DO i=1,sdx
                np=geoPHYS(i,j,k); 
                n=n+1
                DO m=1,numfun
                    IF (Fun(m)% nomsch == np) THEN
                        v_e(n) = int( np, 1)
                        m1 = m1+1
                    ENDIF
                ENDDO
            ENDDO
        ENDDO
    ENDDO

    n=0; m2=0
    DO k=1,sdz
        DO j=1,sdy
            DO i=1,sdx
                np=geoPHYS(i,j,k); 
                n=n+1
                DO m=1,sizdom_coil_sec
                    if (dom_coil_sec(m)%num_coil_sec == np) then
                        v_e(n) = int( np, 1)
                        m2 = m2 + 1
                    ENDIF
                ENDDO
            ENDDO
        ENDDO
    ENDDO

    NULLIFY(  sp_F, next_F   )
    fun_nod(:)%numnod_Fx=0;  fun_nod(:)%numnod_Fy=0; fun_nod(:)%numnod_Fz=0; 
    
    nn = 0;
    DO k=1,sdz;DO j=1,sdy;DO i=1,sdx
        nn =  nn + 1
        if(v_e(nn)/=0) THEN 
            DO m=1,numfun
                IF (Fun(m)%nomsch == v_e(nn)) THEN   
                    fun_nod(m)%numnod_Fx = fun_nod(m)%numnod_Fx + 1
                    ALLOCATE(sp_F)
                    sp_F=tmp(m,nn,'X',next_F) 
                    next_F => sp_F
                    
                    fun_nod(m)%numnod_Fy = fun_nod(m)%numnod_Fy + 1
                    ALLOCATE(sp_F)
                    sp_F=tmp(m,nn,'Y',next_F) 
                    next_F => sp_F
                    
                    fun_nod(m)%numnod_Fz = fun_nod(m)%numnod_Fz + 1
                    ALLOCATE(sp_F)
                    sp_F=tmp(m,nn,'Z',next_F) 
                    next_F => sp_F
                ENDIF
            ENDDO
        ENDIF
    ENDDO; ENDDO;ENDDO !i j k

    nn = 0; m1  = 0; st='0'
    DO k=1,sdz;DO j=1,sdy;DO i=1,sdx
        nn =  nn + 1 
        if(v_e(nn)/=0) THEN 
            do m = 1, sizdom_coil_sec
                if (dom_coil_sec(m)%num_coil_sec == v_e(nn) ) THEN  !dom_coil_sec(m)%num_coil_sec
                    L = dom_coil_sec(m)%num_coil 
                    if    (dom_coil_sec(m)%connect == 'XP' .or.  dom_coil_sec(m)%connect == 'XM' .or. &
                           dom_coil_sec(m)%connect == 'YP' .or.  dom_coil_sec(m)%connect == 'YM' .or. & 
                           dom_coil_sec(m)%connect == 'ZP' .or.  dom_coil_sec(m)%connect == 'ZM' ) then
                        if (nn /= m1  .and. trim(st) /= dom_coil_sec(m)%name_sec) then
                            m1 = nn; st = dom_coil_sec(m)%name_sec
                            
                            fun_nod(L)%numnod_Fx = fun_nod(L)%numnod_Fx + 1
                            ALLOCATE(sp_F)
                            sp_F=tmp(L,nn,'X',next_F) 
                            next_F => sp_F
                                
                            fun_nod(L)%numnod_Fy = fun_nod(L)%numnod_Fy + 1
                            ALLOCATE(sp_F)
                            sp_F=tmp(L,nn,'Y',next_F) 
                            next_F => sp_F
                    
                            fun_nod(L)%numnod_Fz = fun_nod(L)%numnod_Fz + 1
                            ALLOCATE(sp_F)
                            sp_F=tmp(L,nn,'Z',next_F) 
                            next_F => sp_F
                        endif
                    else
                        fun_nod(L)%numnod_Fx = fun_nod(L)%numnod_Fx + 1
                        ALLOCATE(sp_F)
                        sp_F=tmp(L,nn,'X',next_F) 
                        next_F => sp_F
                    
                        fun_nod(L)%numnod_Fy = fun_nod(L)%numnod_Fy + 1
                        ALLOCATE(sp_F)
                        sp_F=tmp(L,nn,'Y',next_F) 
                        next_F => sp_F
                    
                        fun_nod(L)%numnod_Fz = fun_nod(L)%numnod_Fz + 1
                        ALLOCATE(sp_F)
                        sp_F=tmp(L,nn,'Z',next_F) 
                        next_F => sp_F
                    endif
                endif
            enddo
        ENDIF
    ENDDO; ENDDO;ENDDO !i j k
    
    num_nodx = 0; num_nody = 0; num_nodz =0 ; 
    DO m=1,numfun 
        IF (fun_nod(m)%numnod_Fx /= 0) THEN
            n = fun_nod(m)%numnod_Fx 
            ALLOCATE(fun_nod(m)%nods_Fx(n), source=0)  
            ALLOCATE(fun_nod(m)%nods_Fx_cos(n),source=0.d0)
            num_nodx = num_nodx + n
        ENDIF
        IF (fun_nod(m)%numnod_Fy /= 0) THEN  
            n = fun_nod(m)%numnod_Fy 
            ALLOCATE(fun_nod(m)%nods_Fy(n), source=0)  
            ALLOCATE(fun_nod(m)%nods_Fy_cos(n),source=0.d0) 
            num_nody = num_nody + n
        ENDIF 
        IF ( fun_nod(m)%numnod_Fz /= 0) THEN 
            n = fun_nod(m)%numnod_Fz 
            ALLOCATE(fun_nod(m)%nods_Fz(n), source=0) 
            ALLOCATE(fun_nod(m)%nods_Fz_cos(n),source=0.d0)
            num_nodz = num_nodz + n
        ENDIF
    ENDDO

    num_nodV = num_nodx + num_nody + num_nodz
    IF (num_nodV ==0 ) print*, 'NOT EXIST SOURCE of COILS!!!'
    fun_nod(:)%numnod_Fx=0; fun_nod(:)%numnod_Fy=0; fun_nod(:)%numnod_Fz=0; 

    DO m=1,num_nodV 
        k  = sp_F%jm 
        n  = sp_F%im 
        ch = sp_F%ch   
        IF       (trim(ch) == 'X') THEN 
            fun_nod(k)%numnod_Fx = fun_nod(k)%numnod_Fx + 1
            k0 = fun_nod(k)%numnod_Fx
            fun_nod(k)%nods_Fx(k0) = n  ! v loc yzel
        ELSEIF (trim(ch) == 'Y') THEN
            fun_nod(k)%numnod_Fy = fun_nod(k)%numnod_Fy + 1
            k0 = fun_nod(k)%numnod_Fy
            fun_nod(k)%nods_Fy(k0) = n + Cells
        ELSEIF (trim(ch) == 'Z') THEN
            fun_nod(k)%numnod_Fz = fun_nod(k)%numnod_Fz + 1
            k0 = fun_nod(k)%numnod_Fz
            fun_nod(k)%nods_Fz(k0) = n + 2*Cells
        ELSE
        ENDIF
        
        sp_F => sp_F%prec
    ENDDO
ENDIF ! numfun

allocate(temp(num_nodV), source = 0)
i = 0 
DO n=1,numfun
    DO k = 1,fun_nod(n)%numnod_Fx
        m = fun_nod(n)%nods_Fx(k) 
        i = i + 1
        temp(i) = m
    ENDDO
    DO k = 1,fun_nod(n)%numnod_Fy
        m = fun_nod(n)%nods_Fy(k)
        i = i + 1
        temp(i) = m
    ENDDO
    DO k = 1,fun_nod(n)%numnod_Fz
        m = fun_nod(n)%nods_Fz(k)
        i = i + 1
        temp(i) = m
    ENDDO
ENDDO

i = num_nodV + 1
DO n=numfun,1,-1
    DO k = 1,fun_nod(n)%numnod_Fz
        i = i - 1
        fun_nod(n)%nods_Fz(k) =  temp(i)
    ENDDO
    DO k = 1,fun_nod(n)%numnod_Fy                
        i = i - 1
        fun_nod(n)%nods_Fy(k) =  temp(i)
    ENDDO
    DO k = 1,fun_nod(n)%numnod_Fx
        i = i - 1
        fun_nod(n)%nods_Fx(k) =  temp(i)
    ENDDO
ENDDO

DEALLOCATE(temp)

COSINs: BLOCK
    INTEGER:: i_row(6), ind(1), i,j, k, L, m, n, kdz, nim, nip, njm, njp, nkm, nkp
    REAL(8):: sxm,sym,szm,sm,  sx,sy,sz
    LOGICAL mask(6)

    kdz=sdx*sdy
do L = 1, sizdom_coil
    dom_coil(L)%num_fun = L
    i = dom_coil(L)%num_fun 
    do n = 1, dom_coil(L)%siz_nods
        m = dom_coil(L)%nods(n)
        nim = m-1;  njm = m-sdx;  nkm = m-kdz;
        nip = m+1;  njp = m+sdx;  nkp = m+kdz;
        ind = FINDLOC(dom_coil(L)%nods, nim) 
        i_row(1) = ind(1); nim = ind(1)
        ind = FINDLOC(dom_coil(L)%nods, nip)
        i_row(2) = ind(1); nip = ind(1)
        ind = FINDLOC(dom_coil(L)%nods, njm)
        i_row(3) = ind(1); njm = ind(1)
        ind = FINDLOC(dom_coil(L)%nods, njp)
        i_row(4) = ind(1); njp = ind(1)
        ind = FINDLOC(dom_coil(L)%nods, nkm)
        i_row(5) = ind(1); nkm = ind(1)
        ind = FINDLOC(dom_coil(L)%nods, nkp)
        i_row(6) = ind(1); nkp = ind(1)
        sx=0.5d0;  sy=0.5d0;  sz=0.5d0;
        mask = .false.
        mask = i_row <= 0
        if  (mask(1))  then; sx=1.d0; nim=n; endif
        if  (mask(2))  then; sx=1.d0; nip=n; endif
        if  (mask(3))  then; sy=1.d0; njm=n; endif
        if  (mask(4))  then; sy=1.d0; njp=n; endif
        if  (mask(5))  then; sz=1.d0; nkm=n; endif
        if  (mask(6))  then; sz=1.d0; nkp=n; endif
        sxm = -sy*(dom_coil(L)%Uc(nip) - dom_coil(L)%Uc(nim))/delta(1) 
        sym = -sz*(dom_coil(L)%Uc(njp) - dom_coil(L)%Uc(njm))/delta(2)
        szm = -sx*(dom_coil(L)%Uc(nkp) - dom_coil(L)%Uc(nkm))/delta(3) 
        sm = sqrt(sxm*sxm + sym*sym + szm*szm)
        sxm = sxm/sm; sym = sym/sm; szm = szm/sm;
        fun_nod(i)%nods_Fx(n) = m
        fun_nod(i)%nods_Fy(n) = m + Cells
        fun_nod(i)%nods_Fz(n) = m + 2*Cells
        fun_nod(i)%nods_Fx_cos(n) = sxm
        fun_nod(i)%nods_Fy_cos(n) = sym
        fun_nod(i)%nods_Fz_cos(n) = szm  
    enddo
enddo

sxm = delta(2)*delta(3)
sym = delta(1)*delta(3)
szm = delta(1)*delta(2)

do L = 1, sizdom_coil_sec
    k = dom_coil_sec(L)%orient
    do m = 1, dom_coil_sec(L)%siz_nods
        n = dom_coil_sec(L)%nods(m)
        DO i=1,numfun
            ind = FINDLOC(fun_nod(i)%nods_Fx, n)
            if (ind(1) /= 0) then
                j=ind(1) 
                fun_nod(i)%nods_Fx(j) = n
                fun_nod(i)%nods_Fy(j) = n + Cells
                fun_nod(i)%nods_Fz(j) = n + 2*Cells
                fun_nod(i)%nods_Fx_cos(j)= 0.d0
                fun_nod(i)%nods_Fy_cos(j)= 0.d0        
                fun_nod(i)%nods_Fz_cos(j)= 0.d0  
                if ( dom_coil_sec(L)%connect  == 'N1' ) then
                    if     (k==1) then; fun_nod(i)%nods_Fx_cos(j)= -1.d0; 
                    elseif (k==2) then; fun_nod(i)%nods_Fx_cos(j)=  1.d0;
                    elseif (k==3) then; fun_nod(i)%nods_Fy_cos(j)= -1.d0;
                    elseif (k==4) then; fun_nod(i)%nods_Fy_cos(j)=  1.d0;
                    elseif (k==5) then; fun_nod(i)%nods_Fz_cos(j)= -1.d0;
                    elseif (k==6) then; fun_nod(i)%nods_Fz_cos(j)=  1.d0;
                    else
                        print*,'error, not fund orient for sec',dom_coil_sec(L)%name_sec
                        stop
                    endif
                elseif ( dom_coil_sec(L)%connect == 'N2') then
                    if     (k==1) then; fun_nod(i)%nods_Fx_cos(j)=  1.d0; dom_coil(i)%S = dom_coil(i)%S + sxm;
                    elseif (k==2) then; fun_nod(i)%nods_Fx_cos(j)= -1.d0; dom_coil(i)%S = dom_coil(i)%S + sxm;
                    elseif (k==3) then; fun_nod(i)%nods_Fy_cos(j)=  1.d0; dom_coil(i)%S = dom_coil(i)%S + sym 
                    elseif (k==4) then; fun_nod(i)%nods_Fy_cos(j)= -1.d0; dom_coil(i)%S = dom_coil(i)%S + sym
                    elseif (k==5) then; fun_nod(i)%nods_Fz_cos(j)=  1.d0; dom_coil(i)%S = dom_coil(i)%S + szm
                    elseif (k==6) then; fun_nod(i)%nods_Fz_cos(j)= -1.d0; dom_coil(i)%S = dom_coil(i)%S + szm
                    else
                        print*,'error, not fund orient for sec',dom_coil_sec(L)%name_sec
                        stop
                    endif
                elseif (dom_coil_sec(L)%connect == 'XP' .or. dom_coil_sec(L)%connect == 'XM' .or. &
                        dom_coil_sec(L)%connect == 'YP' .or. dom_coil_sec(L)%connect == 'YM' .or. &
                        dom_coil_sec(L)%connect == 'ZP' .or. dom_coil_sec(L)%connect == 'ZM' ) then
                    if ( index(dom_coil_sec(L)%name_sec, '_D') /=0) then
                        if     (k==1) then; fun_nod(i)%nods_Fx_cos(j)=  1.d0; 
                        elseif (k==2) then; fun_nod(i)%nods_Fx_cos(j)= -1.d0;
                        elseif (k==3) then; fun_nod(i)%nods_Fy_cos(j)=  1.d0;
                        elseif (k==4) then; fun_nod(i)%nods_Fy_cos(j)= -1.d0;
                        elseif (k==5) then; fun_nod(i)%nods_Fz_cos(j)=  1.d0;
                        elseif (k==6) then; fun_nod(i)%nods_Fz_cos(j)= -1.d0;
                        endif
                    else
                        if     (k==1) then; fun_nod(i)%nods_Fx_cos(j)= -1.d0; dom_coil(i)%S = dom_coil(i)%S + 0.5d0*sxm
                        elseif (k==2) then; fun_nod(i)%nods_Fx_cos(j)=  1.d0; dom_coil(i)%S = dom_coil(i)%S + 0.5d0*sxm
                        elseif (k==3) then; fun_nod(i)%nods_Fy_cos(j)= -1.d0; dom_coil(i)%S = dom_coil(i)%S + 0.5d0*sym
                        elseif (k==4) then; fun_nod(i)%nods_Fy_cos(j)=  1.d0; dom_coil(i)%S = dom_coil(i)%S + 0.5d0*sym
                        elseif (k==5) then; fun_nod(i)%nods_Fz_cos(j)= -1.d0; dom_coil(i)%S = dom_coil(i)%S + 0.5d0*szm
                        elseif (k==6) then; fun_nod(i)%nods_Fz_cos(j)=  1.d0; dom_coil(i)%S = dom_coil(i)%S + 0.5d0*szm
                        endif
                    endif
                endif
            endif
        enddo
    enddo

enddo

print '(a,$ )',  'cross-sectional ' 
DO n=1,numfun
    dom_coil(n)%S = dom_coil(n)%S/dom_coil(n)%parall
    print '(a,a,a, e10.3,a,$)',  ' coil:', dom_coil(n)%name_coil, ' S=', dom_coil(n)%S ,';'
    DO k = 1,fun_nod(n)%numnod_Fx
        fun_nod(n)%nods_Fx_cos(k)= fun_nod(n)%nods_Fx_cos(k) /dom_coil(n)%S
        fun_nod(n)%nods_Fy_cos(k)= fun_nod(n)%nods_Fy_cos(k) /dom_coil(n)%S
        fun_nod(n)%nods_Fz_cos(k)= fun_nod(n)%nods_Fz_cos(k) /dom_coil(n)%S
    ENDDO
ENDDO
print*
end BLOCK COSINs

IF (numfun /=0) THEN
    NULLIFY(sp_F, next_F )
    DEALLOCATE(nameFun)
ENDIF
IF (numMech /=0 ) DEALLOCATE( nameVmech)
IF (numVenv /=0 ) DEALLOCATE( nameVenv)

DEALLOCATE(v_e, st, ch)
DEALLOCATE(dom_coil, dom_coil_sec)

WRITE(*,'( a)') 'input data complet!'

CONTAINS

SUBROUTINE gen_sparse_coil_matrix

logical mask(7) !nij,
integer i_row(7),coltmp(7), ind(1), countC, Ltmp, i,j, m,n,L, LL,  nim, nip, njm, njp, nkm, nkp,k0, k1,k2,nFix, nFiy, nFiz,  num_nztmp
integer, allocatable :: temp(:) 
real(8) v_row(7), valtmp(7), sx, sy, sz, s
real(8), allocatable :: v_row_add(:)

type espm
    integer im 
    DOUBLE PRECISION em
    type (espm),pointer ::prec
end type espm
type(espm),pointer :: sp_coltmp,next_tmp

sz = 1.d0/(delta(3)*delta(3))
sy = 1.d0/(delta(2)*delta(2))
sx = 1.d0/(delta(1)*delta(1)) 
s = 2.d0*(sx + sy + sz);
nullify(sp_coltmp,next_tmp)
do LL = 1, sizdom_coil
    dom_coil(LL)%nnz = 0;
    countC = 0
    num_nztmp=0
    do n = 1, dom_coil(LL)%siz_nods   
        L = dom_coil(LL)%nods(n)
        valtmp = 0.0d0; coltmp = 0; Ltmp=0; 
        nFix=0; nFiy=0; nFiz=0  
        countC = countC + 1 
        nim = L-1;  njm = L-sdx;  nkm = L-kdz;
        nip = L+1;  njp = L+sdx;  nkp = L+kdz;
        i_row(1:7) =   [nim, nip, njm, njp, nkm, nkp, L ] 
        v_row(1:7) =  [-sx, -sx, -sy, -sy, -sz, -sz,  s ]
        
        ind = FINDLOC(dom_coil(LL)%nods, nim) 
        i_row(1) = ind(1)
        ind = FINDLOC(dom_coil(LL)%nods, nip)
        i_row(2) = ind(1)
        ind = FINDLOC(dom_coil(LL)%nods, njm)
        i_row(3) = ind(1)
        ind = FINDLOC(dom_coil(LL)%nods, njp)
        i_row(4) = ind(1)
        ind = FINDLOC(dom_coil(LL)%nods, nkm)
        i_row(5) = ind(1)
        ind = FINDLOC(dom_coil(LL)%nods, nkp)
        i_row(6) = ind(1)
        ind = FINDLOC(dom_coil(LL)%nods, L)
        i_row(7) = ind(1)
        mask = .false.
        mask = i_row > 0
        temp = pack([(i, i=1, 7)], mask)
        Ltmp = size(temp)
        coltmp(1:Ltmp) = i_row(temp)
        valtmp(1:Ltmp) = v_row(temp)
        ! 8 corners, 12 edges, 6 faces   
        if ( Ltmp /= 7) then 
            allocate( v_row_add(Ltmp) ,source=0.d0)
            j=0
            do i = 1,6,2
                if    (mask(i) .and. mask(i+1)) then
                    j=j+1; v_row_add(j) = 1.d0
                    j=j+1; v_row_add(j) = 1.d0
                ELSEIF (mask(i) ) then
                    j=j+1; v_row_add(j) = 2.d0
                ELSEIF (mask(i+1)) then
                    j=j+1; v_row_add(j) = 2.d0
                endif
            enddo
            valtmp(1:Ltmp-1) = v_row(temp(1:Ltmp-1))*v_row_add(1:Ltmp-1);
            valtmp(Ltmp) = -sum (valtmp(1:Ltmp-1))
        endif
        k0=0 
        mm2: do k1=1,Ltmp-1  
            do k2=k1+1, Ltmp
                if  (coltmp(k1) == coltmp(k2) ) then 
                    k0 = coltmp(k2)
                    exit mm2
                endif
            enddo
        enddo mm2
        if (k0 /=0 ) then
            print*,'node dom_coil double', k0 !, 'i=',i, 'j=',j, 'k=',k
            stop
        endif
        call full_sort(coltmp, valtmp, Ltmp, 1,1)
        do m = 1, Ltmp 
            if (coltmp(m) <= 0) then
                print*, 'cell=',n, 'm',m,'coltmp(m)',coltmp(m)
                stop
            endif
            allocate(sp_coltmp)
            sp_coltmp=espm(coltmp(m), valtmp(m), next_tmp) 
            num_nztmp = num_nztmp + 1 
            next_tmp => sp_coltmp
        enddo
        
        dom_coil(LL)%irow(countC + 1) = dom_coil(LL)%irow( countC) + Ltmp
        deallocate(temp)
        if (allocated (v_row_add)) deallocate ( v_row_add)
    enddo

    if (num_nztmp /= 0) then
        dom_coil(LL)%nnz = num_nztmp
        if (countC /= dom_coil(LL)%siz_nods) stop 'countC /= div_C(LL)%siz_nods'
        allocate(dom_coil(LL)%jcol(num_nztmp), source=0)
        allocate(dom_coil(LL)%valp(num_nztmp), source=0.d0)
    endif
    i= num_nztmp
    DO WHILE (associated(sp_coltmp))
        dom_coil(LL)%jcol(i) = sp_coltmp%im  
        dom_coil(LL)%valp(i) = sp_coltmp%em
        i=i-1
        sp_coltmp => sp_coltmp%prec  
    END DO
    
    nullify(sp_coltmp, next_tmp)
    PRINT '(a,a,  a,  i7, a,g12.5,a)', 'sparse matrix for coil:', dom_coil(LL)%name_coil,  ' non zero=',  &
                        dom_coil(LL)%nnz, ' density=', 100.0* REAL(num_nztmp)/REAL(countC)/REAL(countC),'%'
enddo

end SUBROUTINE gen_sparse_coil_matrix
!====================

SUBROUTINE readline(line,ch,ier, idev)
    CHARACTER(len=:),ALLOCATABLE,INTENT(out) :: line ,ch
    INTEGER,INTENT(out)                      :: ier
    INTEGER                                  :: idev
    CHARACTER :: buffer
    line=''
    ier=0
    do
        READ( idev, '(a)', advance='no', iostat=ier ) buffer
        IF ( ier /= 0 ) THEN
            EXIT
        ELSE 
            line = line // buffer
        ENDIF
    END DO
    ch = line
END SUBROUTINE readline

FUNCTION evaluate (string) 
    REAL(8)  :: evaluate
    CHARACTER(*) string

    SELECT CASE (string(1:1))
        CASE ('"',"'",'`')          ! function specified in quotes
            m=LEN_TRIM(string)
            CALL initf(1); CALL parsef(1,trim(string(2:m-1)),varconst(1:12));
            evaluate = evalf(1,valconst(1:12))
        CASE DEFAULT
            evaluate = numeric( string )
    END SELECT
END FUNCTION evaluate 

SUBROUTINE calcVmech 
    DO n=1,6  
        IF ( j+1 + n+1 <= neww) THEN                        !   number for Vsx or Vsy Vsz
            IF    ( index( words(j+1 + n ),'VSX') /=0 ) THEN 
                fun_nod(numfun)%move(1) = 1                 ! it is move X
                SELECT CASE (words(j+1 + n+1)(1:1) )        ! analysis of the character after the = sign
                    CASE('A':'Z')                           ! function name specified
                        numMech = numMech + 1               ! this is for the calculation of the mechanical function itself
                        Vmech(numMech)% nomsch = kp
                        Vmech(numMech)%ex = 'X'
                        nameVmech(numMech) = trim(words(j+1+ n+1) )
                        fun_nod(numfun)%num_Vmech(1) = numMech  ! the number of the mechanical func is recorded
                        fun_nod(numfun)%val_Vmech(1) = 0.d0
                    CASE DEFAULT                                ! the rest are either a number or quotes
                        fun_nod(numfun)%num_Vmech(1) = 0        ! writing the number of the mechanical function = 0
                        fun_nod(numfun)%val_Vmech(1) = evaluate( words(j+1 + n+1) )
                END SELECT
            ELSEIF ( index( words(j+1 + n ),'VSY') /=0 ) THEN
                fun_nod(numfun)%move(2) = 1                 ! it is move Y
                SELECT CASE (words(j+1 + n+1)(1:1) )        ! analysis of the character after the = sign
                    CASE('A':'Z') ! function name specified
                        numMech = numMech + 1               ! this is for the calculation of the mechanical function itself
                        Vmech(numMech)% nomsch = kp
                        Vmech(numMech)%ex = 'Y'
                        nameVmech(numMech) = trim(words(j+1+ n+1) )
                        fun_nod(numfun)%num_Vmech(2) = numMech   ! the number of the mechanical func is recorded
                        fun_nod(numfun)%val_Vmech(2) = 0.d0 
                    CASE DEFAULT                                 ! the rest are either a number or quotes
                        fun_nod(numfun)%num_Vmech(2) = 0         ! writing the number of the mechanical function = 0
                        fun_nod(numfun)%val_Vmech(2) = evaluate( words(j+1 + n+1) )
                END SELECT
            ELSEIF ( index( words(j+1 + n ),'VSZ') /=0 ) THEN
                fun_nod(numfun)%move(3) = 1 ! it is move Z
                SELECT CASE (words(j+1 + n+1)(1:1) ) 
                    CASE('A':'Z')
                        numMech = numMech + 1 
                        Vmech(numMech)% nomsch = kp
                        Vmech(numMech)%ex = 'Z'
                        nameVmech(numMech) = trim(words(j+1+ n+1) )
                        fun_nod(numfun)%num_Vmech(3) = numMech 
                        fun_nod(numfun)%val_Vmech(3) = 0.d0 
                    CASE DEFAULT 
                        fun_nod(numfun)%num_Vmech(3) = 0 
                        fun_nod(numfun)%val_Vmech(3) = evaluate( words(j+1 + n+1) )
                END SELECT
            ENDIF
        ENDIF
    ENDDO
END SUBROUTINE calcVmech

SUBROUTINE calcVenv
    DO n=0,5 
        IF ( j + n+1 <= neww) THEN                          !   number for Vsx or Vsy Vsz
            IF    ( index( words(j + n ),'VEX') /=0 ) THEN 
                fun_Venv(numVenv)%move(1) = 1
                SELECT CASE (words(j + n+1)(1:1) )          ! analysis of the character after the = sign
                    CASE('A':'Z')                           ! function name specified
                        ! numVenv = numVenv + 1               ! this is for the calculation of the mechanical function itself
                        Venv(numVenv)% nomsch = kp
                        Venv(numVenv)%ex = 'X'
                        nameVenv(numVenv) = trim(words(j+ n+1) )
                        fun_Venv(numVenv)%num_Venv(1) = numVenv 
                        fun_Venv(numVenv)%val_Venv(1) = 0.d0 
                    CASE DEFAULT                                ! the rest are either a number or quotes
                        valPHYS(kp,3) = evaluate( words(j + n+1) ) * valPHYS(kp,2)  !Vx * mu0*sigma
                        fun_Venv(numVenv)%num_Venv(1) = 0 
                        fun_Venv(numVenv)%val_Venv(1) = evaluate( words(j+1 + n+1) )
                END SELECT
                exit
            ELSEIF ( index( words(j+ n ),'VEY') /=0 ) THEN
                fun_Venv(numVenv)%move(2) = 1
                SELECT CASE (words(j + n+1)(1:1) )           ! analysis of the character after the = sign
                    CASE('A':'Z')                            ! function name specified
                        ! numVenv = numVenv + 1                ! this is for the calculation of the mechanical function itself
                        Venv(numVenv)% nomsch = kp
                        Venv(numVenv)%ex = 'Y'
                        nameVenv(numVenv) = trim(words(j+ n+1) )
                        fun_Venv(numVenv)%num_Venv(2) = numVenv 
                        fun_Venv(numVenv)%val_Venv(2) = 0.d0 
                    CASE DEFAULT                                 ! the rest are either a number or quotes
                        valPHYS(kp,4) = evaluate( words(j + n+1) ) * valPHYS(kp,2)
                        fun_Venv(numVenv)%num_Venv(2) = 0 
                        fun_Venv(numVenv)%val_Venv(2) = evaluate( words(j+1 + n+1) )
                END SELECT
                exit
            ELSEIF ( index( words(j + n ),'VEZ') /=0 ) THEN
                fun_Venv(numVenv)%move(3) = 1
                SELECT CASE (words(j + n+1)(1:1) ) 
                    CASE('A':'Z')
                        ! numVenv = numVenv + 1                ! this is for the calculation of the mechanical function itself
                        Venv(numVenv)% nomsch = kp
                        Venv(numVenv)%ex = 'Z'
                        nameVenv(numVenv) = trim(words(j+ n+1) )
                        fun_Venv(numVenv)%num_Venv(3) = numVenv 
                        fun_Venv(numVenv)%val_Venv(3) = 0.d0 
                    CASE DEFAULT 
                        valPHYS(kp,5) = evaluate( words(j + n+1) ) * valPHYS(kp,2)
                        fun_Venv(numVenv)%num_Venv(3) = 0 
                        fun_Venv(numVenv)%val_Venv(3) = evaluate( words(j+1 + n+1) )
                END SELECT
            ENDIF
        ENDIF
    ENDDO
END SUBROUTINE calcVenv

END SUBROUTINE vxc2data

SUBROUTINE sprsBCGstabWR (valA, irow, jcol, n, b, x, tolerance, itmax, iter)
! Fortran code created by J.Sochor   ( https://github.com/JNSresearcher )
! e-mail: JNSresearcher@gmail.com

INTEGER jcol(*), irow(n+1)
REAL(8) valA(*)
INTEGER iter, itmax,  n 
REAL(8):: alpha, beta, omega, tolerance, rr0, rr0_new, Bnorm, &
                    R(n), R0(n), P(n), AP(n), S(n), AS(n), B(n), X(n)
    iter=0 
    R = sprsAx(X)
    do j=1,n
        R(j) = B(j) - R(j)
    enddo 
    R0 = R
    P = R

    Bnorm=norm2(b)
    if (Bnorm == 0.d0) return

    do
        if (iter > itmax) then
            print*, norm2(R)
            exit
        endif
        iter = iter + 1
        AP = sprsAx(P)
        rr0 = DOT_PRODUCT (R,R0)
        alpha = rr0/dot_product(AP,R0) 
        S = R - alpha*AP
        if (norm2(S) < tolerance) then
            X = X + alpha*P;
            exit
        endif
        AS = sprsAx(S)
        omega = dot_product(AS,S)/dot_product(AS,AS)
        X = X + alpha*P + omega*S
        R = S - omega*AS
        if ( (norm2(R)/Bnorm) < tolerance) exit
        rr0_new = dot_product(R,R0)
        beta = (alpha/omega)*rr0_new/rr0
        P = R + beta*(P - omega*AP)
        if ( (abs(rr0_new)/Bnorm) < tolerance) then
            R0 = R; P = R
        endif
    enddo

contains

function sprsAx  (  V) 
    REAL(8)  :: sprsAx(n), V(n)
    INTEGER i1,i2
    do  i=1,n
        i1=irow(i); i2=irow(i+1)-1
        sprsAx(i) =  dot_product(valA(i1:i2), V(jcol(i1:i2) ) )
    enddo 
END function sprsAx 

END SUBROUTINE sprsBCGstabWR


SUBROUTINE string2words(Ligne,words,num_words)
! converting string to array of words
CHARACTER (LEN = *), INTENT(IN) :: Ligne
CHARACTER (LEN = *), INTENT(OUT):: words(*) 
INTEGER Start
INTEGER Finish
INTEGER num_words
INTEGER length 
    length= len_trim(Ligne)
    Start = 0;Finish=0;num_words=0
    DO j = 1, length
        IF (Ligne(j:j) == ' ' .or. Ligne(j:j)== '	') THEN
            IF (Start > 0) THEN
                num_words = num_words + 1
                words(num_words) = Ligne(Start : Finish)
                Start = 0
            ENDIF
        ELSE IF (Start == 0) THEN
            Start = j;  Finish = j
        ELSE
            Finish = Finish + 1
        ENDIF
    ENDDO
    IF (Start > 0) THEN
        num_words = num_words + 1
        words(num_words) = Ligne(Start : Finish)
        Start = 0
    ENDIF

END SUBROUTINE string2words

SUBROUTINE UPP (st1, st2)
    CHARACTER (LEN=*),  INTENT(in) :: st1
    CHARACTER (LEN=*), INTENT(out) :: st2
    CHARACTER (LEN=*),   PARAMETER :: lsym = 'abcdefghijklmnopqrstuvwxyz'
    CHARACTER (LEN=*),   PARAMETER :: usym = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    st2 = st1
    DO j=1,LEN_TRIM(st1)
        k = INDEX(lsym,st1(j:j))
        IF (k > 0) st2(j:j) = usym(k:k)
    END DO
END SUBROUTINE UPP

REAL(8) function numeric(SA)  !numerus
! Fortran code created by J.Sochor   ( https://github.com/JNSresearcher )
! e-mail: JNSresearcher@gmail.com
!                                       1   2   3    4    5   6   7    8     9
    character*1 :: ch,prefixes(9)=(/'M','K','U','N','P','G','T','F','S'/) 
    character (len = *) sa
    character*20 inter
    REAL(8) multiplier   !factor
    multiplier = 1.0_8;lp=0;lm=0
            m=index(SA,',')
            if (m.NE.0) SA(m:m)='.'
    do i=1,9
        l=index(SA,prefixes(i))
        if (l.NE.0) then
            m=1
            if ( i.EQ.1 ) then
                lm=index(SA, 'MEG')
                if (lm /= 0) then
                    multiplier=1.E+06_8
                else
                    multiplier=1.E-3_8
                endif
                exit
            elseif ( i.EQ.2 ) then
                multiplier=1.E+03_8
                exit
            elseif ( i.EQ.3 ) then
                multiplier=1.E-06_8
                exit
            elseif ( i.EQ.4 ) then
                multiplier=1.E-09_8
                exit
            elseif ( i.EQ.5 ) then
                lm=index(SA, 'PET')
                if (lm /= 0) then
                    multiplier=1.E+15_8
                else
                    multiplier=1.E-12_8
                endif
                exit
            elseif ( i.EQ.6 ) then
                multiplier=1.E+09_8
                exit
            elseif ( i.EQ.7 ) then
                multiplier=1.E+12_8
                exit
            elseif ( i.EQ.8 ) then
                multiplier=1.E-15_8
                exit
            elseif ( i.EQ.9 ) then
                multiplier=1.E-2_8
                exit
            endif
        endif
    enddo

    k=index(SA,'.');j=len_trim(SA)
    if ( (m.NE.0).AND.(k.EQ.0) )  then
        SA(l:l)='.'
        if (lm /= 0) then
            SA(lm+1:lm+2)='  '
            sa=sa(1:lm)//sa(lm+3:)
        endif
    else
        if (lm /= 0) then
            SA(lm:lm+2)='   '
            sa=sa(1:lm)//sa(lm+3:)
        endif
    endif
    j=len_trim(SA)
    i=index(SA,'E')
    if (i.eq.0) then
        do i=1,j
            ch=SA(i:i)
            n=ichar(ch)
            if ((n.GT.57).or.(n.LT.48)) then
                if (ch.EQ.'.') then
                    cycle
                elseif (ch == '-') then
                    cycle
                else
                    SA(i:i)=' '
                endif
            else
            endif
        enddo
    endif

    write (inter,'(A20)') SA(1:j)
    read (inter,'(BN,G20.0)') numeric
    numeric = multiplier * numeric

end function numeric