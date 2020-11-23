























!
!     C H I R A L    T W O B O D Y    P O T E N T I A L S
!
!**********************************************************************
!
! PROGRAM PACKAGE TO COMPUTE ALL POSSIBLE CHIRAL TWOBODY POTENTIALS
! AT LEADING-ORDER (LO), NEXT-TO LEADING-ORDER (NLO)
! AND NEXT-TO-NEXT-TO LEADING-ORDER (NNLO)
!
! THIS CODE IS BASED ON THE WORK AND HELP OF RUPRECHT MACHLEIDT et al.
! AT MOSCOW UNIVERSITY, IDAHO.
!
! AUTHOR: Andreas Ekström
! ADDRESS:  Dept. Physics, University Oslo, N-0316 OSLO
! EMAIL: jaeks@fys.uio.no
! LAST UPDATE: 2015-05-07
!
! N3LO-addition:
!   AUTHOR: Boris D. Carlsson
!   EMAIL: borisc@chalmers.se
!
! IMPORTANT: COMPILE WITH INTEL FORTRAN
! OTHER COMPILERS, SUCH AS PORTLAND GIVE
! WRONG RESULTS. WILL CORRECT THIS IN THE
! FUTURE.
!
!**********************************************************************
MODULE chp_aux


  IMPLICIT NONE
  
! basic mesh info
  TYPE, PUBLIC :: chp_mesh_info
     INTEGER  :: amount
     REAL(8) :: xmin, xmax
  END TYPE chp_mesh_info
  
! GAUSS-LEGENDRE MESH POINTS
! Gauss-Legendre mesh point x, corresponding integration weight w and corresponding x*x*w-value
  TYPE, PUBLIC :: chp_gauleg_mesh_point
     SEQUENCE
     REAL(8) :: x, w, xxw
  END TYPE chp_gauleg_mesh_point
  
! mesh points and weights in momentum space
  TYPE, PUBLIC :: chp_gauleg_mesh
     TYPE(chp_mesh_info)                                    :: info
     TYPE(chp_gauleg_mesh_point), DIMENSION(:), ALLOCATABLE :: pnt
  END TYPE chp_gauleg_mesh

  TYPE, PUBLIC :: chp_real_type_RA
     REAL(8)       :: val
     CHARACTER(LEN=12)   :: name
     LOGICAL             :: set
  END TYPE chp_real_type_RA

  TYPE, PUBLIC :: chp_real_type
     REAL(8)           :: val
     CHARACTER(LEN=12):: name
     LOGICAL          :: set
  END TYPE chp_real_type
  
  TYPE, PUBLIC :: chp_char2_type
     CHARACTER(LEN=2) :: val
     CHARACTER(LEN=12):: name
     LOGICAL          :: set
  END TYPE chp_char2_type
  
  TYPE, PUBLIC :: chp_int_type
     INTEGER          :: val
     CHARACTER(LEN=12):: name
     LOGICAL          :: set
  END TYPE chp_int_type

  TYPE, PUBLIC :: chp_chn_type
     INTEGER          :: S,j,T,tz
     LOGICAL          :: coup
  END TYPE chp_chn_type
  
! mathematical constants
  TYPE(chp_real_type), PARAMETER :: chp_pi=chp_real_type(ACOS(-1.0D0),'pi',.TRUE.)
  REAL(8)             , PARAMETER :: fourpi=4.0D0*chp_pi%val
  REAL(8)             , PARAMETER :: twopi =2.0D0*chp_pi%val
  REAL(8)             , PARAMETER :: pi2   = chp_pi%val*chp_pi%val
  REAL(8)             , PARAMETER :: pi_inv= 1D0 / chp_pi%val
  
!--- standard formatted contacts to pw contacts conversion matrix
  REAL(8), PUBLIC :: LOST2PW(2,2)
  REAL(8), PUBLIC :: LOPW2ST(2,2)
  REAL(8), PUBLIC :: NLOST2PW(7,7)
  REAL(8), PUBLIC :: NLOPW2ST(7,7)
  REAL(8), PUBLIC :: N3LOST2PW(15,15)
  REAL(8), PUBLIC :: N3LOPW2ST(15,15)

CONTAINS
!
!  This function sets up the recursive relation
!  for the associated Legendre polynomials
!
  REAL(8) FUNCTION chp_legendre_polynomials(l, m, x)
    IMPLICIT NONE
    INTEGER,  INTENT(IN) :: l, m
    REAL(8), INTENT(IN) :: x
    REAL(8)             :: fact,pll,pmm,pmmp1,somx2
    INTEGER              :: i,ll
    

!  check whether m, l and x are ok
    IF((M < 0).OR.(M > L).OR.(ABS(X) > 1.0D0)) THEN
       WRITE(*,*) 'legendre_polynomials: bad arguments', m, l, x
       chp_legendre_polynomials = 0.0D0
       RETURN
    ENDIF

!  calculate now pmm as starting point for iterations
    pmm=1.0
    IF (m > 0) THEN
       somx2=SQRT((1.0-x)*(1.0+x))
       fact=1.0D0;
       DO i=1, m
          pmm = -fact*somx2*pmm
          fact = fact+2.0D0
       ENDDO
    ENDIF

!  if l == m we do not need to use recursion relation
    IF (l == m) THEN
       chp_legendre_polynomials=pmm

!  recursive relation for associated Legendre polynomials
    ELSE
       pmmp1=x*(2*m+1)*pmm

!  analytical formula for the case l == m+1
       IF (l == (m+1)) THEN
          chp_legendre_polynomials=pmmp1
       ELSE
          DO ll=m+2, l
             pll=(x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m)
             pmm=pmmp1
             pmmp1=pll
          ENDDO
          chp_legendre_polynomials= pll
       ENDIF
    ENDIF
    
  END FUNCTION chp_legendre_polynomials
  
!
!      This routine calculates gauss-legendre mesh points and weights
!      INPUT:
!      mesh%info%xmin     : lower limit of the integration interval
!      mesh%info%xmax     : upper limit ---------- "" -------------
!      mesh%info%amount   : the desired number of mesh points
!      OUTPUT:
!      mesh%pnt(:)%x      : gauss-legendre mesh points on the interval (x1,x2)
!      mesh%pnt(:)%w      : the corresponding weights
!      FROM               : Numerical recipes
!      F90 version        : M. Hjorth-Jensen
!      Object interface   : M. Kartamyshev
!
  SUBROUTINE chp_setup_gauleg_mesh (mesh)


    TYPE(chp_gauleg_mesh), INTENT(INOUT)        :: mesh
    INTEGER                                 :: i, j, m, n
    REAL(8)                                :: x1, x2
    REAL(8), DIMENSION(:), ALLOCATABLE     :: x, w
    REAL(8)                                :: p1,p2,p3,pp,xl,xm,z,z1
    REAL(8), PARAMETER                     :: EPS = 3.D-14
    
    ALLOCATE(x(mesh%info%amount))
    ALLOCATE(w(mesh%info%amount))

! allocate points and weights storages
    CALL chp_destroy_gauleg_mesh (mesh)
    
    IF (mesh%info%amount <= 0 .OR. mesh%info%xmin > mesh%info%xmax) THEN
       WRITE(*,*) ': incorrect mesh info', mesh%info ; STOP
    ENDIF
    
    ALLOCATE( mesh%pnt( 1:mesh%info%amount ) )
    mesh%pnt(:) = chp_gauleg_mesh_point(0.0D0, 0.0D0, 0.0D0)
    
! set values of local variables
    x1 = mesh%info%xmin ; x2 = mesh%info%xmax; n = mesh%info%amount
    
    m=(n+1)/2
    xm=0.5D0*(x2+x1)
    xl=0.5D0*(x2-x1)
    DO i=1,m
       z1=0.0D0
       z=COS(chp_pi%val*(i - 0.25D0)/(n + 0.5D0))
       DO WHILE ( ABS(z-z1) > EPS)
          p1=1.0D0
          p2=0.0D0
          DO j=1,n
             p3=p2
             p2=p1
             p1=((2.0D0*j-1.0D0)*z*p2-(j-1.0D0)*p3)/j
          END DO
          pp=n*(z*p1-p2)/(z*z-1.0D0)
          z1=z
          z=z-p1/pp
       END DO
       x(i)=xm-xl*z
       x(n+1-i)=xm+xl*z
       w(i)=2.0D0*xl/((1.0D0-z*z)*pp*pp)
       w(n+1-i)=w(i)
    ENDDO
    
! set return values
    mesh%pnt(:)%x = x(:) ; mesh%pnt(:)%w = w(:) ; mesh%pnt(:)%xxw = x(:) * x(:) * w(:)

    DEALLOCATE(w)
    DEALLOCATE(x)

  END SUBROUTINE chp_setup_gauleg_mesh
  
  SUBROUTINE chp_destroy_gauleg_mesh (mesh)
    TYPE(chp_gauleg_mesh), INTENT(INOUT) :: mesh
    IF (ALLOCATED(mesh%pnt) ) DEALLOCATE (mesh%pnt)
  END SUBROUTINE chp_destroy_gauleg_mesh

! phase factor of type (-1)**arg
  FUNCTION chp_minus_power (arg) RESULT (res)
    INTEGER              :: res
    INTEGER, INTENT (IN) :: arg
    INTEGER :: exponent
    
    exponent = ABS(arg) ! (-1)**N = (-1)**(-N)
    
    SELECT CASE (MOD(exponent,2))
    CASE (0)
       res =  1
    CASE (1)
       res = -1
    END SELECT

  END FUNCTION chp_minus_power

!
!     Given an NxN matrix A(N,N), this routine replaces it by the LU
!     decomposed one, where the matrix elements are stored in the same
!     matrix A. The array indx is  an output vector which records the row
!     permutation effected by the partial pivoting. d is the determinant
!
  SUBROUTINE chp_lu_decompose(a,n,indx,d)
    
    IMPLICIT NONE
    INTEGER :: n, i, j, k, imax
    REAL(8) :: sum , tiny, aamax, dum, d
    REAL(8), DIMENSION(n,n) :: a
    INTEGER, DIMENSION(n) :: indx
    REAL(8), ALLOCATABLE :: vv(:)
    
    tiny=1.0e-20
    ALLOCATE ( vv(n) )
    D=1.
    DO i=1,n
       aamax=0.
       DO j=1,n
          IF (ABS(a(i,j)) > aamax) aamax=ABS(a(i,j))
       ENDDO
!     Zero is the largest element
       IF (aamax == 0.) STOP 'Singular matrix.'
!     No nonzero largest element
       vv(i)=1./aamax
    ENDDO
!     loop over columns
    DO j=1,n
!     solves equation 2.3.12 except for i=j of Numerical Recipes
       IF (j > 1) THEN
          DO i=1,j-1
             sum=a(i,j)
             IF (i > 1)THEN
                DO k=1,i-1
                   sum=sum-a(i,k)*a(k,j)
                ENDDO
                a(i,j)=sum
             ENDIF
          ENDDO
       ENDIF
!    start searching for largest pivot element
       aamax=0.
       DO i=j,n
          sum=a(i,j)
          IF (j > 1)THEN
             DO k=1,j-1
                sum=sum-a(i,k)*a(k,j)
             ENDDO
             a(i,j)=sum
          ENDIF
          dum=vv(i)*ABS(sum)
          IF (dum >= aamax) THEN
             imax=i
             aamax=dum
          ENDIF
       ENDDO
!    interchange of rows
       IF (j /= imax)THEN
          DO k=1,n
             dum=a(imax,k)
             a(imax,k)=a(j,k)
             a(j,k)=dum
          ENDDO
!    change of parity for determinant
          d=-d
          vv(imax)=vv(j)
       ENDIF
       indx(j)=imax
       IF(j /= n) THEN
          IF(a(j,j) == 0.) a(j,j)=tiny
          dum=1./a(j,j)
          DO i=j+1,n
             a(i,j)=a(i,j)*dum
          ENDDO
       ENDIF
!    set up determinant
       d=d*a(j,j)
    ENDDO
    IF(a(n,n) == 0.)  a(n,n)=tiny
    DEALLOCATE ( vv)
    
  END SUBROUTINE chp_lu_decompose
  
!     Solves set of linear equations Ax=b, A is input as an LU decompomsed
!     matrix and indx keeps track of the permutations of the rows. b is input
!     as the right-hand side vector b and returns the solution x. A, n and indx
!     are not modified by this routine. This function takes into that b can contain
!     many zeros and is therefore suitable for matrix inversion
  
  
  SUBROUTINE chp_lu_linear_equation(a,n,indx,b)
    
    IMPLICIT NONE
    INTEGER :: n, ii, ll, i, j
    REAL(8) :: sum
    REAL(8), DIMENSION(n,n) :: a
    REAL(8), DIMENSION(n) :: b
    INTEGER, DIMENSION(n) :: indx
    
    ii=0
!     First we solve equation 2.3.6 of numerical recipes
    DO i=1,n
       ll=indx(i)
       sum=b(ll)
       b(ll)=b(i)
       IF (ii /= 0)THEN
          DO j=ii,i-1
             sum=sum-a(i,j)*b(j)
          ENDDO
       ELSEIF (sum /= 0.) THEN
          ii=i
       ENDIF
       b(i)=sum
    ENDDO
!     then we solve equation 2.3.7
    DO i=n,1,-1
       sum=b(i)
       IF (i < n) THEN
          DO j=i+1,n
             sum=sum-a(i,j)*b(j)
          ENDDO
       ENDIF
!     store a component of the solution x in the same place as b
       b(i)=sum/a(i,i)
    ENDDO
    
  END SUBROUTINE chp_lu_linear_equation
  
!            Routines to do mtx inversion, from Numerical
!            Recepies, Teukolsky et al. Routines included
!            below are MATINV, LUDCMP and LUBKSB. See chap 2
!            of Numerical Recipes for further details
!            Recoded in FORTRAN 90 by M. Hjorth-Jensen
!
  SUBROUTINE chp_matinv(a,n)
    
    IMPLICIT NONE
    
    INTEGER, INTENT(IN) :: n
    INTEGER :: i, j
    REAL(8), DIMENSION(n,n), INTENT(INOUT)  :: a
    REAL(8), ALLOCATABLE :: y(:,:)
    REAL(8) :: d
    INTEGER, ALLOCATABLE :: indx(:)
    
    ALLOCATE (y( n, n))  ; ALLOCATE ( indx (n))
    y=0.
!     setup identity matrix
    DO i=1,n
       y(i,i)=1.
    ENDDO
!     LU decompose the matrix just once
    CALL  chp_lu_decompose(a,n,indx,d)
    
!     Find inverse by columns
    DO j=1,n
       CALL chp_lu_linear_equation(a,n,indx,y(:,j))
    ENDDO
!     The original matrix a was destroyed, now we equate it with the inverse y
    a=y
    
    DEALLOCATE ( y ); DEALLOCATE ( indx )
    
  END SUBROUTINE chp_matinv
  
  SUBROUTINE set_contact_conversion_matrices
    
    LOST2PW = transpose(reshape((/ 1.0D0, -3.0D0,  &
         1.0D0,  1.0D0 /), shape(LOST2PW)))
    
    LOST2PW = fourpi * LOST2PW

    LOPW2ST = LOST2PW
    
    CALL chp_matinv(LOPW2ST,2)

    NLOST2PW = transpose(reshape((/ 1.0D0, 0.25D0, -3.0D0, -0.75D0, 0.0D0, -1.0D0, -0.25D0, &
                                    -2.0D0/3.0D0, 1.0D0/6.0D0, -2.0D0/3.0D0, 1.0D0/6.0D0, -2.0D0/3.0D0, 2.0D0, - 0.5D0, &
                                    -2.0D0/3.0D0, 1.0D0/6.0D0, 2.0D0, -0.5D0, 0.0D0, 2.0D0/3.0D0, -1.0D0/6.0D0, &
                                    -2.0D0/3.0D0, 1.0D0/6.0D0, -2.0D0/3.0D0, 1.0D0/6.0D0, -1.0D0/3.0D0, -4.0D0/3.0D0, 1.0D0/3.0D0, &
                                    1.0D0, 0.25D0, 1.0D0, 0.25D0, 0.0D0, 1.0D0/3.0D0, 1.0D0/12.0D0, &
                                    0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, -2.0D0*DSQRT(2.0D0)/3.0D0, -DSQRT(2.0D0)/6.0D0, &
                                    -2.0D0/3.0D0, 1.0D0/6.0D0, -2.0D0/3.0D0, 1.0D0/6.0D0, +1.0D0/3.0D0, 0.0D0, 0.0D0 /), shape(NLOST2PW)))
    NLOST2PW = fourpi * NLOST2PW

    NLOPW2ST = NLOST2PW
    
    CALL chp_matinv(NLOPW2ST,7)
    
    N3LOST2PW = transpose(reshape((/ &
                    1D0, 1D0/16D0, 0.25D0, 0D0, -3D0, -3D0/16D0, -0.75D0, 0D0, 0D0, 0D0, -1D0, -0.25D0, -0.25D0, -1D0/16D0, 0D0, &
                    10D0/3D0, 5D0/24D0, 1D0/6D0, 2D0/3D0, -10D0, -0.625D0, -0.5D0, -2D0, 0D0, 0D0, -10D0/3D0, -1D0/6D0, -1D0/6D0, -5D0/24D0, -2D0/3D0, &
                    -4D0/3D0, 1D0/12D0, 0D0, 0D0, -4D0/3D0, 1D0/12D0, 0D0, 0D0, -2D0/3D0, -1D0/6D0, 8D0/3D0, 1D0/3D0, -1D0/3D0, -1D0/6D0, 0D0, &
                    -4D0/3D0, 1D0/12D0, 0D0, 0D0, 4D0, -0.25D0, 0D0, 0D0, 0D0, 0D0, 4D0/3D0, 0D0, 0D0, -1D0/12D0, 0D0, &
                    -4D0/3D0, 1D0/12D0, 0D0, 0D0, -4D0/3D0, 1D0/12D0, 0D0, 0D0, -1D0/3D0, -1D0/12D0, -2D0, -1D0/6D0, 1D0/6D0, 1D0/8D0, 0D0, &
                    1D0, 1D0/16D0, 0.25D0, 0D0, 1D0, 1D0/16D0, 0.25D0, 0D0, 0D0, 0D0, 1D0/3D0, 1D0/12D0, 1D0/12D0, 1D0/48D0, 0D0, &
                    10D0/3D0, 5D0/24D0, 1D0/6D0, 2D0/3D0, 10D0/3D0, 5D0/24D0, 1D0/6D0, 2D0/3D0, 0D0, 0D0, 10D0/9D0, 1D0/18D0, 1D0/18D0, 5D0/72D0, 2D0/9D0, &
                    8D0/15D0, 1D0/30D0, -2D0/15D0, -2D0/15D0, 8D0/15D0, 1D0/30D0, -2D0/15D0, -2D0/15D0, 2D0/5D0, -0.1D0, -4D0/9D0, 1D0/9D0, 1D0/9D0, -1D0/36D0, -16D0/45D0, &
                    0D0, 0D0, 0D0, 0D0, 0D0, 0D0, 0D0, 0D0, 0D0, 0D0, -DSQRT(2D0)*2D0/3D0, -DSQRT(2D0)/6D0, -DSQRT(2D0)/6D0, -DSQRT(2D0)/24D0, 0D0, &
                    0D0, 0D0, 0D0, 0D0, 0D0, 0D0, 0D0, 0D0, 0D0, 0D0, -DSQRT(2D0)*14D0/9D0, DSQRT(2D0)/18D0, DSQRT(2D0)/18D0, -DSQRT(2D0)*7D0/72D0, DSQRT(2D0)*2D0/9D0, &
                    8D0/15D0, 1D0/30D0, -2D0/15D0, -2D0/15D0, -8D0/5D0, -0.1D0, 0.4D0, 0.4D0, 0D0, 0D0, -8D0/15D0, 2D0/15D0, 2D0/15D0, -1D0/30D0, 2D0/15D0, &
                    8D0/15D0, 1D0/30D0, -2D0/15D0, -2D0/15D0, 8D0/15D0, 1D0/30D0, -2D0/15D0, -2D0/15D0, 2D0/15D0, -1D0/30D0, 0.8D0, -0.2D0, -0.2D0, 0.05D0, 4D0/15D0, &
                    -4D0/3D0, 1D0/12D0, 0D0, 0D0, -4D0/3D0, 1D0/12D0, 0D0, 0D0, 1D0/3D0, 1D0/12D0, -2D0/15D0, 1D0/30D0, -1D0/30D0, 1D0/120D0, 0D0, &
                    0D0, 0D0, 0D0, 0D0, 0D0, 0D0, 0D0, 0D0, 0D0, 0D0, DSQRT(6D0)*4D0/15D0, -DSQRT(6D0)/15D0, DSQRT(6D0)/15D0, -DSQRT(6D0)/60D0, 0D0, &
                    8D0/15D0, 1D0/30D0, -2D0/15D0, -2D0/15D0, 8D0/15D0, 1D0/30D0, -2D0/15D0, -2D0/15D0, -4D0/15D0, 1D0/15D0, 0D0, 0D0, 0D0, 0D0, -2D0/15D0 &
                /), shape(N3LOST2PW)))
    N3LOST2PW = fourpi * N3LOST2PW
    N3LOPW2ST = N3LOST2PW
    CALL chp_matinv(N3LOPW2ST, 15)

  END SUBROUTINE set_contact_conversion_matrices
     
END MODULE chp_aux

! Ref: K. Erkelenz et al., NPA176, 413-432 (1971)
MODULE twobody_pwd
  
  USE chp_aux

  IMPLICIT NONE
  
  INTEGER, PRIVATE :: NZ, JMAX

  TYPE(chp_gauleg_mesh) :: zmesh
!                 z,J
  REAL(8), DIMENSION(:,:), ALLOCATABLE, PRIVATE :: legP
!            z**l z,l
  REAL(8), DIMENSION(:,:), ALLOCATABLE, PRIVATE :: zl
CONTAINS
  
  SUBROUTINE setup_twobody_pwd(set_NZ, set_JMAX)
    
    INTEGER :: set_NZ, set_JMAX
        
    NZ   = set_NZ
    JMAX = set_JMAX

    CALL setup_zmesh
    CALL setup_legendre_polynomials
    
  END SUBROUTINE setup_twobody_pwd
  
  SUBROUTINE print_pwd_numerics(unit)
    
    INTEGER, INTENT(IN) :: unit

    WRITE(unit,"(A8,I6)") 'JMAX:', JMAX
    WRITE(unit,"(A8,I6)") 'NZ  :', NZ

  END SUBROUTINE print_pwd_numerics

  SUBROUTINE setup_zmesh
    
    zmesh%info = chp_mesh_info(NZ, -1.0D0, +1.0D0)
    CALL chp_setup_gauleg_mesh(zmesh)
    
  END SUBROUTINE setup_zmesh

  SUBROUTINE setup_legendre_polynomials
    
    INTEGER :: J, iz
    REAL(8)  :: z
    REAL(8)  :: val, err

    IF (ALLOCATED(legP)) DEALLOCATE(legP)
    IF (ALLOCATED(zl)) DEALLOCATE(zl)
    ALLOCATE(legP(1:zmesh%info%amount,0:JMAX))
    ALLOCATE(zl(1:zmesh%info%amount,0:JMAX))
    DO J=0, JMAX
       DO iz=1, zmesh%info%amount
          
          z = zmesh%pnt(iz)%x
          legP(iz,J) = chp_legendre_polynomials(J,0,z)
          zl(iz,J) = z**J
       END DO

       val = sum(legP(:,J) * zmesh%pnt(:)%w)
       if(J == 0) then
           err = abs(val - 2)
       else
           err = abs(val)
       end if
       if(err > 1.e-13) then
           write(*,*) 'Legendre polynomial ', J, ' bad, error is ', err
           stop
       end if
    END DO
    
  END SUBROUTINE setup_legendre_polynomials
  
!   c     pwd of central force
! ref: K. Erkelenz et al. , NPA 176, (1971), 413
  SUBROUTINE pwd_c(W, CHN, pfinal, pinit, POT)
    
    REAL(8)     ,INTENT(INOUT) :: W(1:NZ)
    TYPE(chp_chn_type),INTENT(IN)    :: CHN
    REAL(8)            ,INTENT(IN)    :: pfinal, pinit
    REAL(8)     ,INTENT(INOUT) :: POT(1:6)
    INTEGER :: j
    
    j = CHN%j
! uncoupled singlet
    IF (.NOT. CHN%coup .AND. CHN%S==0) THEN
       
       POT(1) = POT(1) + 2.0D0*pwd_integral(W,0,j)
       
    END IF
    
    IF (.NOT. CHN%coup .AND. CHN%S==1) THEN
       
       POT(2) = POT(2) + 2.0D0*pwd_integral(W,0,j)
       
! once here, return
       RETURN
    END IF
    
! coupled channels
    IF (CHN%coup) THEN
       
!++Vj
       POT(3) = POT(3) + 2.0D0*pwd_integral(W,0,j+1)

!--Vj
       IF (j==0) RETURN
       POT(4) = POT(4) + 2.0D0*pwd_integral(W,0,j-1)
       
       POT(5) = POT(5) + 0.0D0
       
       POT(6) = POT(6) + 0.0D0

    END IF
    
  END SUBROUTINE pwd_c

!   c     pwd of central force
! ref: E. Epelbaum et al., NPA 747 (2005), 362
!SUBROUTINE pwd_c_E(W, CHN, pfinal, pinit, POT)
!
!  REAL(8)            ,INTENT(INOUT) :: W(1:NZ)
!  TYPE(chp_chn_type),INTENT(IN)    :: CHN
!  REAL(8)            ,INTENT(IN)    :: pfinal, pinit
!  REAL(8)            ,INTENT(INOUT) :: POT(1:6)
!  INTEGER :: j
!
!  j = CHN%j
!  ! uncoupled singlet
!  IF (.NOT. CHN%coup .AND. CHN%S==0) THEN
!
!     POT(1) = POT(1) + 2.0D0*pwd_integral(W,0,j)
!
!  END IF
!
!  IF (.NOT. CHN%coup .AND. CHN%S==1) THEN
!
!     POT(2) = POT(2) + 2.0D0*pwd_integral(W,0,j)
!
!     ! once here, return
!     RETURN
!  END IF
!
!  ! coupled channels
!  IF (CHN%coup) THEN
!
!     POT(3) = POT(3) + 2.0D0*pwd_integral(W,0,j+1)
!     IF (j==0) RETURN
!     POT(4) = POT(4) + 2.0D0*pwd_integral(W,0,j-1)
!
!     POT(5) = POT(5) + 0.0D0
!
!     POT(6) = POT(6) + 0.0D0
!
!  END IF
!
!END SUBROUTINE pwd_c_E

!   s     spin-spin
! ref: K. Erkelenz et al. , NPA 176, (1971), 413
  SUBROUTINE pwd_s(W, CHN, pfinal, pinit, POT)
    
    REAL(8)     ,INTENT(INOUT) :: W(1:NZ)
    TYPE(chp_chn_type),INTENT(IN)    :: CHN
    REAL(8)            ,INTENT(IN)    :: pfinal, pinit
    REAL(8)     ,INTENT(INOUT) :: POT(1:6)
    INTEGER :: j
    
    j = CHN%j
! uncoupled singlet
    IF (.NOT. CHN%coup .AND. CHN%S==0) THEN
       
       POT(1) = POT(1) - 6.0D0*pwd_integral(W,0,j)
       
    ELSE IF (.NOT. CHN%coup .AND. CHN%S==1) THEN
       
       POT(2) = POT(2) + 2.0D0*pwd_integral(W,0,j)
       
! coupled channels
    ELSE IF (CHN%coup) THEN
       
       POT(3) = POT(3) + 2.0D0*pwd_integral(W,0,j+1)
       IF (j==0) RETURN
       POT(4) = POT(4) + 2.0D0*pwd_integral(W,0,j-1)
       
       POT(5) = POT(5) + 0.0D0
       
       POT(6) = POT(6) + 0.0D0

    END IF
    
  END SUBROUTINE pwd_s

!   s     spin-spin
! ref: E. Epelbaum et al., NPA 747 (2005), 362
!SUBROUTINE pwd_s_E(W, CHN, pfinal, pinit, POT)
!
!  REAL(8)            ,INTENT(INOUT) :: W(1:NZ)
!  TYPE(chp_chn_type),INTENT(IN)    :: CHN
!  REAL(8)            ,INTENT(IN)    :: pfinal, pinit
!  REAL(8)            ,INTENT(INOUT) :: POT(1:6)
!  INTEGER :: j
!
!  j = CHN%j
!  ! uncoupled singlet
!  IF (.NOT. CHN%coup .AND. CHN%S==0) THEN
!
!     POT(1) = POT(1) - 6.0D0*pwd_integral(W,0,j)
!
!  ELSE IF (.NOT. CHN%coup .AND. CHN%S==1) THEN
!
!     POT(2) = POT(2) + 2.0D0*pwd_integral(W,0,j)
!
!  ! coupled channels
!  ELSE IF (CHN%coup) THEN
!
!     POT(3) = POT(3) + 2.0D0*pwd_integral(W,0,j+1)
!     IF (j==0) RETURN
!     POT(4) = POT(4) + 2.0D0*pwd_integral(W,0,j-1)
!
!     POT(5) = POT(5) + 0.0D0
!
!     POT(6) = POT(5) + 0.0D0
!
!  END IF
!
!END SUBROUTINE pwd_s_E

!   LS    spin-orbit
! ref: K. Erkelenz et al. , NPA 176, (1971), 413
  SUBROUTINE pwd_LS(W, CHN, pfinal, pinit, POT)
    
    REAL(8)     ,INTENT(IN)    :: W(1:NZ)
    TYPE(chp_chn_type),INTENT(IN)    :: CHN
    REAL(8)            ,INTENT(IN)    :: pfinal, pinit
    REAL(8)     ,INTENT(INOUT) :: POT(1:6)
    REAL(8)  :: jj, p, pp, p2, pp2, rj
    INTEGER :: j
    
    rj = REAL(CHN%j,kind=8)
    jj = 2.0D0*rj+1.0D0
    j = CHN%j
    
    pp  = pfinal; p  = pinit
    pp2 = pp*pp ; p2 = p*p
! uncoupled singlet
    IF (.NOT. CHN%coup .AND. CHN%S==0) THEN
       
       POT(1) = POT(1) + 0.0D0
       
    END IF
    
    IF (.NOT. CHN%coup .AND. CHN%S==1) THEN
       
       POT(2) = POT(2) + 2.0D0*pp*p*( pwd_integral(W,0,j+1) - pwd_integral(W,0,j-1))/jj
       
! once here, return
       RETURN
    END IF
    
! coupled channels
    IF (CHN%coup) THEN
       
       POT(3) = POT(3) + 2.0D0*pp*p*(rj+2.0D0)*( pwd_integral(W,0,j+2) - pwd_integral(W,0,j))/(2.0D0*rj+3.0D0)
       IF (j==0) RETURN
       POT(4) = POT(4) + 2.0D0*pp*p*(rj-1.0D0)*( pwd_integral(W,0,j-2) - pwd_integral(W,0,j))/(2.0D0*rj-1.0D0)
       
       POT(5) = POT(5) + 0.0D0
       
       POT(6) = POT(6) + 0.0D0

    END IF
    
  END SUBROUTINE pwd_LS

!  LS     spin-orbit
! ref: E. Epelbaum et al., NPA 747 (2005), 362
!SUBROUTINE pwd_LS_E(W, CHN, pfinal, pinit, POT)
!
!  REAL(8)            ,INTENT(IN)    :: W(1:NZ)
!  TYPE(chp_chn_type),INTENT(IN)    :: CHN
!  REAL(8)            ,INTENT(IN)    :: pfinal, pinit
!  REAL(8)            ,INTENT(INOUT) :: POT(1:6)
!  REAL(8)  :: jj, p, pp, p2, pp2, rj
!  INTEGER :: j
!
!  rj = REAL(CHN%j,kind=8)
!  jj = 2.0D0*rj+1.0D0
!  j = CHN%j
!
!  pp  = pfinal; p  = pinit
!  pp2 = pp*pp ; p2 = p*p
!
!  ! uncoupled singlet
!  IF (.NOT. CHN%coup .AND. CHN%S==0) THEN
!
!     POT(1) = POT(1) + 0.0D0
!
!  END IF
!
!  IF (.NOT. CHN%coup .AND. CHN%S==1) THEN
!
!     POT(2) = POT(2) + 4.0D0*pp*p*pwd_integral(W,1,j) - 2.0D0*pp*p*(pwd_integral(W,0,j-1) + &
!          pwd_integral(W,0,j+1))
!
!     ! once here, return
!     RETURN
!  END IF
!
!  ! coupled channels
!  IF (CHN%coup) THEN
!
!     POT(3) = POT(3) - 2.0D0*pp*p*pwd_integral(W,0,j)+2.0D0*pp*p*pwd_integral(W,1,j+1)
!     IF (j==0) RETURN
!     POT(4) = POT(4) - 2.0D0*pp*p*pwd_integral(W,0,j)+2.0D0*pp*p*pwd_integral(W,1,j-1)
!
!     POT(5) = POT(5) + 0.0D0
!
!     POT(6) = POT(6) + 0.0D0
!
!  END IF
!
!END SUBROUTINE pwd_LS_E
  
!   sigL  sigma-L
! ref: K. Erkelenz et al. , NPA 176, (1971), 413
  SUBROUTINE pwd_sigL(W, CHN, pfinal, pinit, POT)
    
    REAL(8)     ,INTENT(IN)    :: W(1:NZ)
    TYPE(chp_chn_type),INTENT(IN)    :: CHN
    REAL(8)            ,INTENT(IN)    :: pfinal, pinit
    REAL(8)     ,INTENT(INOUT) :: POT(1:6)
    REAL(8)  :: jj, jj1, p, pp, p2, pp2, rj
    INTEGER :: j
    
    rj = REAL(CHN%j,kind=8)
    jj = 2.0D0*rj+1.0D0
    jj1 = rj*(rj+1.0D0)
    j = CHN%j
    
    pp  = pfinal; p  = pinit
    pp2 = pp*pp ; p2 = p*p
! uncoupled singlet
    IF (.NOT. CHN%coup .AND. CHN%S==0) THEN
       
       POT(1) = POT(1) + 2.0D0*pp2*p2*(pwd_integral(W,2,j) - pwd_integral(W,0,j))
       
    END IF
    
    IF (.NOT. CHN%coup .AND. CHN%S==1) THEN
       
       POT(2) = POT(2) + 2.0D0*pp2*p2*( -1.0D0*pwd_integral(W,0,j)  + &
            ((rj-1.0D0)*pwd_integral(W,1,j+1) + (rj+2.0D0)*pwd_integral(W,1,j-1))/jj)
       
! once here, return
       RETURN
    END IF
    
! coupled channels
    IF (CHN%coup) THEN
       
       POT(3) = POT(3) + 2.0D0*pp2*p2*( -1.0D0*pwd_integral(W,2,j+1)  + &
            ((2.0D0*rj+3.0D0)*pwd_integral(W,0,j+1) - (2.0D0)*pwd_integral(W,1,j))/jj)
       IF (j==0) RETURN
       POT(4) = POT(4) + 2.0D0*pp2*p2*( -1.0D0*pwd_integral(W,2,j-1)  + &
            ((2.0D0*rj-1.0D0)*pwd_integral(W,0,j-1) + (2.0D0)*pwd_integral(W,1,j))/jj)

       POT(5) = POT(5) - DSQRT(jj1)*4.0D0*pp2*p2*(pwd_integral(W,0,j+1) - pwd_integral(W,0,j-1))/jj**2

       POT(6) = POT(6) - DSQRT(jj1)*4.0D0*pp2*p2*(pwd_integral(W,0,j+1) - pwd_integral(W,0,j-1))/jj**2

    END IF
    
  END SUBROUTINE pwd_sigL

!   sigL  sigma-L
! ref: E. Epelbaum et al., NPA 747 (2005), 362
!SUBROUTINE pwd_sigL_E(W, CHN, pfinal, pinit, POT)
!
!  REAL(8)            ,INTENT(IN)    :: W(1:NZ)
!  TYPE(chp_chn_type),INTENT(IN)    :: CHN
!  REAL(8)            ,INTENT(IN)    :: pfinal, pinit
!  REAL(8)            ,INTENT(INOUT) :: POT(1:6)
!  REAL(8)  :: jj, jj1, p, pp, p2, pp2, rj
!  INTEGER :: j
!
!  rj = REAL(CHN%j,kind=8)
!  jj = 2.0D0*rj+1.0D0
!  jj1 = rj*(rj+1.0D0)
!  j = CHN%j
!
!  pp  = pfinal; p  = pinit
!  pp2 = pp*pp ; p2 = p*p
!  ! uncoupled singlet
!  IF (.NOT. CHN%coup .AND. CHN%S==0) THEN
!
!     POT(1) = POT(1) + 2.0D0*pp2*p2*(pwd_integral(W,2,j) -pwd_integral(W,0,j))
!
!  END IF
!
!  IF (.NOT. CHN%coup .AND. CHN%S==1) THEN
!
!     POT(2) = POT(2) - 2.0D0*pp2*p2*(pwd_integral(W,0,j) +3.0D0*pwd_integral(W,2,j) &
!          -2.0D0*(pwd_integral(W,1,j-1) + pwd_integral(W,1,j+1)))
!
!     ! once here, return
!     RETURN
!  END IF
!
!  ! coupled channels
!  IF (CHN%coup) THEN
!
!     POT(3) = POT(3) + 4.0D0*pp2*p2*(pwd_integral(W,0,j+1) - pwd_integral(W,1,j))/jj + &
!          2.0D0*pp2*p2*(pwd_integral(W,0,j+1) - pwd_integral(W,2,j+1))
!     IF (j==0) RETURN
!     POT(4) = POT(4) + 4.0D0*pp2*p2*(pwd_integral(W,1,j) - pwd_integral(W,0,j-1))/jj + &
!          2.0D0*pp2*p2*(pwd_integral(W,0,j-1) - pwd_integral(W,2,j-1))
!
!     POT(5) = POT(5) - 4.0D0*DSQRT(jj1)*pp2*p2*(pwd_integral(W,0,j+1) - pwd_integral(W,0,j-1))/jj**2
!
!     POT(6) = POT(6) - 4.0D0*DSQRT(jj1)*pp2*p2*(pwd_integral(W,0,j+1) - pwd_integral(W,0,j-1))/jj**2
!
!  END IF
!
!END SUBROUTINE pwd_sigL_E

!   T     tensor
! ref: K. Erkelenz et al. , NPA 176, (1971), 413
  SUBROUTINE pwd_T(W, CHN, pfinal, pinit, POT)
    
    REAL(8)     ,INTENT(IN)    :: W(1:NZ)
    TYPE(chp_chn_type),INTENT(IN)    :: CHN
    REAL(8)            ,INTENT(IN)    :: pfinal, pinit
    REAL(8)     ,INTENT(INOUT) :: POT(1:6)
    REAL(8)  :: jj, jj1, p, pp, p2, pp2, rj
    INTEGER :: j
  
    rj  = REAL(CHN%j,kind=8)
    jj  = 2.0D0*rj+1.0D0
    jj1 = rj*(rj+1.0D0)
    j   = CHN%j
    
    pp  = pfinal; p  = pinit
    pp2 = pp*pp ; p2 = p*p
! uncoupled singlet
    IF (.NOT. CHN%coup .AND. CHN%S==0) THEN
       
       POT(1) = POT(1) + 2.0D0*(-(pp2 + p2) * pwd_integral(W,0,j) + &
            2.0D0*pp*p*pwd_integral(W,1,j))
       
    ELSE IF (.NOT. CHN%coup .AND. CHN%S==1) THEN
       
       POT(2) = POT(2) + 2.0D0*( (pp2 + p2)*pwd_integral(W,0,j) - &
            2.0D0*pp*p*(rj*pwd_integral(W,0,j+1) + &
            (rj+1.0D0)*pwd_integral(W,0,j-1))/jj)
       
! coupled channels
    ELSE IF(CHN%coup) THEN
       
       POT(3) = POT(3) + 2.0D0*(-(pp2 + p2) * pwd_integral(W,0,j+1) + &
            2.0D0*pp*p*pwd_integral(W,0,j))/jj
       IF (j==0) RETURN
       POT(4) = POT(4) + 2.0D0*((pp2 + p2) * pwd_integral(W,0,j-1) - &
            2.0D0*pp*p*pwd_integral(W,0,j))/jj
       
       POT(5) = POT(5) - 4.0D0*DSQRT(jj1)*(p2*pwd_integral(W,0,j+1) + &
            pp2*pwd_integral(W,0,j-1) - 2.0D0*pp*p* &
            pwd_integral(W,0,j) )/jj
       
       POT(6) = POT(6) - 4.0D0*DSQRT(jj1)*(p2*pwd_integral(W,0,j-1) + &
            pp2*pwd_integral(W,0,j+1) - 2.0D0*pp*p* &
            pwd_integral(W,0,j) )/jj
  
    END IF
    
  END SUBROUTINE pwd_T

!   T     tensor
! ref: E. Epelbaum et al., NPA 747 (2005), 362
!SUBROUTINE pwd_T_E(W, CHN, pfinal, pinit, POT)
!
!  REAL(8)            ,INTENT(IN)    :: W(1:NZ)
!  TYPE(chp_chn_type),INTENT(IN)    :: CHN
!  REAL(8)            ,INTENT(IN)    :: pfinal, pinit
!  REAL(8)            ,INTENT(INOUT) :: POT(1:6)
!  REAL(8)  :: jj, jj1, p, pp, p2, pp2, rj
!  INTEGER :: j
!
!  rj  = REAL(CHN%j,kind=8)
!  jj  = 2.0D0*rj+1.0D0
!  jj1 = rj*(rj+1.0D0)
!  j   = CHN%j
!
!  pp  = pfinal; p  = pinit
!  pp2 = pp*pp ; p2 = p*p
!
!  ! uncoupled singlet
!  IF (.NOT. CHN%coup .AND. CHN%S==0) THEN
!
!     POT(1) = POT(1)-2.0D0*(pp2+p2)*pwd_integral(W,0,j) + 4.0D0*pp*p*pwd_integral(W,1,j)
!
!  ELSE IF (.NOT. CHN%coup .AND. CHN%S==1) THEN
!
!     POT(2) = POT(2) + 2.0D0*((pp2+p2)*pwd_integral(W,0,j) + 2.0D0*pp*p*pwd_integral(W,1,j) - &
!          2.0D0*pp*p*(pwd_integral(W,0,j-1)+pwd_integral(W,0,j+1)))
!
!     ! coupled channels
!  ELSE IF (CHN%coup) THEN
!
!     POT(3) = POT(3) + (4.0D0*pp*p*pwd_integral(W,0,j)-2.0D0*(pp2+p2)*pwd_integral(W,0,j+1))/jj
!     IF (j==0) RETURN
!     POT(4) = POT(4) + (-4.0D0*pp*p*pwd_integral(W,0,j)-2.0D0*(pp2+p2)*pwd_integral(W,0,j-1))/jj
!
!     POT(5) = POT(5) - DSQRT(jj1)*(-8.0D0*pp*p*pwd_integral(W,0,j) + 4.0D0*pp2*pwd_integral(W,0,j-1) + &
!          4.0D0*p2*pwd_integral(W,0,j+1) )/jj
!
!     POT(6) = POT(6) - DSQRT(jj1)*(-8.0D0*pp*p*pwd_integral(W,0,j) + 4.0D0*pp2*pwd_integral(W,0,j+1) + &
!          4.0D0*p2*pwd_integral(W,0,j-1) )/jj
!
!  END IF
!
!END SUBROUTINE pwd_T_E
  
!   sigk momentum-tensor
! ref: K. Erkelenz et al. , NPA 176, (1971), 413
  SUBROUTINE pwd_sigk(W, CHN, pfinal, pinit, POT)
    
    REAL(8)     ,INTENT(IN)    :: W(1:NZ)
    TYPE(chp_chn_type),INTENT(IN)    :: CHN
    REAL(8)            ,INTENT(IN)    :: pfinal, pinit
    REAL(8)     ,INTENT(INOUT) :: POT(1:6)
    REAL(8)  :: jj, jj1, p, pp, p2, pp2, rj
    INTEGER :: j
    
    rj = REAL(CHN%j,kind=8)
    jj = 2.0D0*rj+1.0D0
    jj1 = rj*(rj+1.0D0)
    j = CHN%j
    
    pp  = pfinal; p  = pinit
    pp2 = pp*pp ; p2 = p*p
! uncoupled singlet
    IF (.NOT. CHN%coup .AND. CHN%S==0) THEN
       
       POT(1) = POT(1) + 0.5D0*(-(pp2 + p2) * pwd_integral(W,0,j) - &
            2.0D0*pp*p*pwd_integral(W,1,j))
       
    END IF
    
    IF (.NOT. CHN%coup .AND. CHN%S==1) THEN
       
       POT(2) = POT(2) + 0.5D0*( (pp2 + p2)*pwd_integral(W,0,j) + &
            2.0D0*pp*p*(rj*pwd_integral(W,0,j+1) + &
            (rj+1.0D0)*pwd_integral(W,0,j-1))/jj)
       
! once here, return
       RETURN
    END IF
    
! coupled channels
    IF (CHN%coup) THEN
       
       POT(3) = POT(3) + 0.5D0*(-(pp2 + p2) * pwd_integral(W,0,j+1) - &
            2.0D0*pp*p*pwd_integral(W,0,j))/jj
       IF (j==0) RETURN
       POT(4) = POT(4) + 0.5D0*((pp2 + p2) * pwd_integral(W,0,j-1) + &
            2.0D0*pp*p*pwd_integral(W,0,j))/jj
       
       POT(5) = POT(5) - DSQRT(jj1)*(p2*pwd_integral(W,0,j+1) + &
            pp2*pwd_integral(W,0,j-1) + 2.0D0*pp*p* &
            pwd_integral(W,0,j))/jj
       
       POT(6) = POT(6) - DSQRT(jj1)*(p2*pwd_integral(W,0,j-1) + &
            pp2*pwd_integral(W,0,j+1) + 2.0D0*pp*p* &
            pwd_integral(W,0,j))/jj

    END IF
    
  END SUBROUTINE pwd_sigk

!   sigk     momentum-tensor
! ref: E. Epelbaum et al., NPA 747 (2005), 362
!SUBROUTINE pwd_sigk_E(W, CHN, pfinal, pinit, POT)
!
!  REAL(8)            ,INTENT(IN)    :: W(1:NZ)
!  TYPE(chp_chn_type),INTENT(IN)    :: CHN
!  REAL(8)            ,INTENT(IN)    :: pfinal, pinit
!  REAL(8)            ,INTENT(INOUT) :: POT(1:6)
!  REAL(8)  :: jj, jj1, p, pp, p2, pp2, rj
!  INTEGER :: j
!
!  rj = REAL(CHN%j,kind=8)
!  jj = 2.0D0*rj+1.0D0
!  jj1 = rj*(rj+1.0D0)
!  j = CHN%j
!
!  pp  = pfinal; p  = pinit
!  pp2 = pp*pp ; p2 = p*p
!
!  ! uncoupled singlet
!  IF (.NOT. CHN%coup .AND. CHN%S==0) THEN
!
!     POT(1) = POT(1)-0.5D0*(pp2+p2)*pwd_integral(W,0,j) - pp*p*pwd_integral(W,1,j)
!
!  END IF
!
!  IF (.NOT. CHN%coup .AND. CHN%S==1) THEN
!
!     POT(2) = POT(2) + 0.5D0*(pp2+p2)*pwd_integral(W,0,j) - pp*p*pwd_integral(W,1,j) + &
!          pp*p*(pwd_integral(W,0,j-1)+pwd_integral(W,0,j+1))
!
!     ! once here, return
!     RETURN
!  END IF
!
!  ! coupled channels
!  IF (CHN%coup) THEN
!
!     POT(3) = POT(3) + (-1.0D0*pp*p*pwd_integral(W,0,j)-0.5D0*(pp2+p2)*pwd_integral(W,0,j+1))/jj
!     IF (j==0) RETURN
!     POT(4) = POT(4) + (+1.0D0*pp*p*pwd_integral(W,0,j)+0.5D0*(pp2+p2)*pwd_integral(W,0,j-1))/jj
!
!     POT(5) = POT(5) - (DSQRT(jj1)*(2.0D0*pp*p*pwd_integral(W,0,j) + pp2*pwd_integral(W,0,j-1) + &
!          p2*pwd_integral(W,0,j+1)))/jj
!
!     POT(6) = POT(6) - (DSQRT(jj1)*(2.0D0*pp*p*pwd_integral(W,0,j) + pp2*pwd_integral(W,0,j+1) + &
!          p2*pwd_integral(W,0,j-1)))/jj
!
!  END IF
!
!END SUBROUTINE pwd_sigk_E
  
  FUNCTION pwd_integral(W,l,j) RESULT(res)
    
    REAL(8), INTENT(IN)  :: W(1:NZ)
    INTEGER, INTENT(IN) :: l, j
    REAL(8) :: res
    
    res = 0.0D0

    IF (j<0) RETURN
    
! not using precalculated z**l
!res = chp_pi%val*SUM(W(:)*zmesh%pnt(:)%x**l*legP(:,j)*zmesh%pnt(:)%w)
! using precalculated z**l (atleast a factor of 2 in speedup)
    res = chp_pi%val*SUM(W(:)*zl(:,l)*legP(:,j)*zmesh%pnt(:)%w)
    
  END FUNCTION pwd_integral
  
END MODULE twobody_pwd

![1] = Machleidts Phys Rep 503 (2011) 1-75 AKA EM2011
![2] = Epelbaum NPA 747 (2005) 362-424     AKA Ep2005
![3] = Machleidts PRC 91 014002 (2015)     AKA EM2015
![4] = Epelbaum EPJA 51, 53 (2015)         AKA Ep2015
MODULE idaho_chiral_potential
  
  USE chp_aux
  USE twobody_pwd

  IMPLICIT NONE
  
!MAXIMUM J
  INTEGER, PARAMETER, PRIVATE :: maximum_angular_momentum = 40 ! maximum angular momentum.
! Z=cos(theta), theta angle between p' and p, i.e., final and initial momenta
! q = (p'-p)
! q^2 = p'^2 + p^2 - 2p'pZ
  INTEGER, PARAMETER, PRIVATE :: nof_theta_int_points   = 96
! x and z are integration variables. Some of the expressions for the
! 2PE 2-loop diagrams contains integrals that must be solved numerically,
! The integration variables are x and z = 2*m_\pi / \mu using the notation
! from [1] and [2]. The number of integration points used here are
! enough to get highly accurate results.
! The integral over x is independent of q, and thus need only be done once
! for each system (np, pp, nn).
  INTEGER, PARAMETER, PRIVATE :: nof_2PE_2loop_int_x_points = 50
  INTEGER, PARAMETER, PRIVATE :: nof_2PE_2loop_int_z_points = 30
  
! FORTRAN UNITS
  INTEGER, PARAMETER, PRIVATE :: CHP_ERR = 6
  INTEGER, PARAMETER, PRIVATE :: CHP_SCR = 6
  
! chiral order
  TYPE(chp_int_type), PRIVATE :: chp_chiral_order ! LO, NLO, NNLO, N3LO
  TYPE(chp_int_type), PRIVATE :: chp_chiral_mode  ! chiral_mode_{EM2011,EM2015,Ep2005,Ep2015}
! chiral order definition
  INTEGER, PARAMETER :: NO_ORDER = -1
  INTEGER, PARAMETER :: LO   = 0
! all contributions vanish at order 1
! due to parity and time-reversal  invariance
  INTEGER, PARAMETER :: NLO  = 2
  INTEGER, PARAMETER :: NNLO = 3
  INTEGER, PARAMETER :: N3LO = 4
  INTEGER, PARAMETER :: chiral_mode_EM2011 = 100 ! Mimic [1]
  INTEGER, PARAMETER :: chiral_mode_EM2015 = 101 ! Mimic [3]
  INTEGER, PARAMETER :: chiral_mode_Ep2005 = 102 ! Mimic [2]
  INTEGER, PARAMETER :: chiral_mode_Ep2015 = 103 ! Mimic [4]
! Summary:
!  base LO  : +1PE +NN.Ct*
!             +minimal_relativity
!  base NLO : +2PE.1-loop_0 +NN.C*
!             +1PE.CIB +NN.Ct.CIB
!  base N2LO: +2PE.1-loop_d
!  base N3LO: +2PE.1loop_r +2PE.1-loop_dd +2PE.2-loop +2PE.2-loop.int +NN.D*
!             +2PE.1-loop_r.mode = "EM2015"
!
!  EM2011 N2LO: +2PE.1-loop_r
!               +2PE.1-loop_r.mode = "EM2011"
!  EM2011 N3LO: +2PE.1-loop_dr +2PE.1-loop_rr -2PE.2-loop.int +1PE.gamma
!               +2PE.CSB.correct_mass
!
!  Ep2005 LO  : -minimal_relativity +kamada_glockle_transform
!  Ep2005 N3LO: +1PE.relcorr +2PE.1-loop_r.mode = "Ep2005"
!
!  Ep2005 LO  : -minimal_relativity +kamada_glockle_transform
!  Ep2015 N3LO: +1PE.relcorr +2PE.1-loop_r.mode = "Ep2015"
  
  REAL(8), PRIVATE  :: rirrel = -99.99D0
  REAL(8),PRIVATE    :: RAirrel
  INTEGER, PRIVATE :: iirrel = -99
!CHARACTER(LEN=2), PRIVATE :: cirrel = 'XX'
! masses, MeV
  TYPE(chp_real_type), PRIVATE :: chp_mnuc(-1:1) ! -1: pp , 0: pn , +1: nn
  TYPE(chp_real_type), PRIVATE :: chp_mpi(-1:2)  ! -1: pi-, 0: pi0, +1: pi+, 2: average
  
! renormalization stuff
  TYPE(chp_real_type), PRIVATE :: chp_ren ! SF or DR, with SF-cut in MeV
! regularization stuff
  TYPE(chp_real_type), PRIVATE :: chp_lambda ! regularization cutoff in MeV, typically ~ 500
  TYPE(chp_real_type), PRIVATE :: chp_regcut_1PE, chp_regcut_2PE
! potential coupling constants
  TYPE(chp_real_type_RA), PRIVATE :: chp_gA  ! dimensionless
  TYPE(chp_real_type_RA), PRIVATE :: chp_fpi ! MeV
  TYPE(chp_real_type), PRIVATE :: chp_fine_structure ! dimensionless

! optional features. All of these have default values depending on chiral_order and chiral_mode
! Thus, these need not be explicitly set if the defaults are to be used.
  TYPE(chp_int_type) , PRIVATE :: chiral_Ct ! 1 or 0 -- NN contact LECs with 0 derivatives
  TYPE(chp_int_type) , PRIVATE :: chiral_Ct_CIB ! 1 or 0 -- Splitting of Ct_1S0 into Ct_1S0{nn,np,pp}
  TYPE(chp_int_type) , PRIVATE :: chiral_C ! 1 or 0 -- NN contact LECs with 1 derivative
  TYPE(chp_int_type) , PRIVATE :: chiral_D ! 1 or 0 -- NN contact LECs with 2 derivatives
  TYPE(chp_int_type) , PRIVATE :: chiral_1PE ! 1 or 0 -- one pion exchange
  TYPE(chp_int_type) , PRIVATE :: chiral_1PE_CIB ! 1 or 0 -- 1PE correct for pion mass differences
  TYPE(chp_int_type) , PRIVATE :: chiral_1PE_gamma ! 1 or 0 -- pion-gamma exchange
  TYPE(chp_int_type) , PRIVATE :: chiral_1PE_relcorr ! 1 or 0 -- relativistic correction to 1PE
  TYPE(chp_int_type) , PRIVATE :: chiral_2PE_1loop_0 ! 1 or 0 -- leading two pion exchange 1-loop diagrams
  TYPE(chp_int_type) , PRIVATE :: chiral_2PE_1loop_d ! 1 or 0 -- 2PE 1-loop prop. to ci
  TYPE(chp_int_type) , PRIVATE :: chiral_2PE_1loop_r ! 1 or 0 -- 2PE 1-loop prop. to 1 / M_N
! chiral_2PE_1loop_r_mode determines what expressions to use for the 2PE 1-loop relativistic diagrams.
! EM2011, EM2015, Ep2005, Ep1015 will use the expressions from [1], [3], [2] and [4] respectively
  TYPE(chp_int_type) , PRIVATE :: chiral_2PE_1loop_r_mode ! chiral_2PE_1loop_r_mode_{EM2011,EM2015,Ep2005,Ep2015}
  INTEGER, PARAMETER :: chiral_2PE_1loop_r_mode_EM2011 = 201
  INTEGER, PARAMETER :: chiral_2PE_1loop_r_mode_EM2015 = 202
  INTEGER, PARAMETER :: chiral_2PE_1loop_r_mode_Ep2005 = 203
  INTEGER, PARAMETER :: chiral_2PE_1loop_r_mode_Ep2015 = 204
  TYPE(chp_int_type) , PRIVATE :: chiral_2PE_1loop_dd ! 1 or 0 -- 2PE 1-loop prop. to ci*cj
  TYPE(chp_int_type) , PRIVATE :: chiral_2PE_1loop_dr ! 1 or 0 -- 2PE 1-loop prop. to ci / M_N
  TYPE(chp_int_type) , PRIVATE :: chiral_2PE_1loop_rr ! 1 or 0 -- 2PE 1-loop prop. to 1 / M_N^2
  TYPE(chp_int_type) , PRIVATE :: chiral_2PE_2loop ! 1 or 0 -- 2PE 2-loop diagrams
  TYPE(chp_int_type) , PRIVATE :: chiral_2PE_2loop_int ! 1 or 0 -- 2PE 2-loop non-analytical terms
  TYPE(chp_int_type) , PRIVATE :: chiral_2PE_CSB_correct_mass ! 1 or 0 -- use correct nucleon mass in 2PE terms
  TYPE(chp_int_type) , PRIVATE :: chiral_minimal_relativity ! 1 or 0 -- use the minimal relativity prescription
  TYPE(chp_int_type) , PRIVATE :: chiral_kamada_glockle_transform ! 1 or 0 -- use Kamada-Glöckle transformation
  
  REAL(8), PRIVATE :: MR_factor ! minimal relativity / kamada-glockle factor that the potential is multiplied with

! variables for calculating optional 2PE 2loop diagrams
  LOGICAL, PRIVATE :: chp_2PE_2loop_int_VTS_data_set
  LOGICAL, PRIVATE :: chp_2PE_2loop_int_WC_DR_data_set
  LOGICAL, PRIVATE :: chp_2PE_2loop_int_WC_SFR_data_set
! Integration meshs for non-analytical integrals in some of the 2PE 2-loop diagrams
  TYPE(chp_gauleg_mesh), PRIVATE :: chp_2PE_2loop_int_mesh_x
  TYPE(chp_gauleg_mesh), PRIVATE :: chp_2PE_2loop_int_mesh_z
  TYPE(chp_gauleg_mesh), PRIVATE :: chp_2PE_2loop_int_mesh_z_VC_SFR
! Factors in the x-integrations that need only be calculated once
  REAL(8) , PRIVATE :: chp_2PE_2loop_int_VTS_x_fact(nof_2PE_2loop_int_z_points)
  REAL(8) , PRIVATE :: chp_2PE_2loop_int_WC_DR_x_fact(nof_2PE_2loop_int_z_points)
  REAL(8) , PRIVATE :: chp_2PE_2loop_int_WC_SFR_x_fact(nof_2PE_2loop_int_z_points)

! LECs ci
  TYPE(chp_real_type_RA), PRIVATE :: chp_c1
  TYPE(chp_real_type_RA), PRIVATE :: chp_c3
  TYPE(chp_real_type_RA), PRIVATE :: chp_c4

! New at N3LO
  TYPE(chp_real_type_RA), PRIVATE :: chp_c2
  TYPE(chp_real_type_RA), PRIVATE :: chp_d1_plus_d2
  TYPE(chp_real_type_RA), PRIVATE :: chp_d3
  TYPE(chp_real_type_RA), PRIVATE :: chp_d5
  TYPE(chp_real_type_RA), PRIVATE :: chp_d14_minus_d15

  TYPE(chp_real_type_RA), PRIVATE :: chp_e14
  TYPE(chp_real_type_RA), PRIVATE :: chp_e15
  TYPE(chp_real_type_RA), PRIVATE :: chp_e16
  TYPE(chp_real_type_RA), PRIVATE :: chp_e17
  TYPE(chp_real_type_RA), PRIVATE :: chp_e18

! CONTACT TERMS
! THEY CAN BE INPUT IN PWD FORMAT
! OR STANDARD FORMAT. THE CODE USES
! THE PWD FORMAT. WHEN PRINTED, THEY
! ARE CONVERTED AND PRINTED IN BOTH
  TYPE(chp_int_type), PRIVATE :: chp_contact_format! PW (partial-wave) or ST (standard)
!
! @LO      PW
! CLO(1) : Ctilde 1S0  CS
! CLO(2) : Ctilde 3S1  CT
  TYPE(chp_real_type_RA), PRIVATE :: chp_CLO(1:2)
  TYPE(chp_real_type), PRIVATE :: chp_CLO_n(1:2)
! @NLO     PW
! CLO(1) : C1S0     C1
! CLO(2) : C3P0     C2
! CLO(3) : C1P1     C3
! CLO(4) : C3P1     C4
! CLO(5) : C3S1     C5
! CLO(6) : C3S1-3D1 C6
! CLO(7) : C3P1     C7
  TYPE(chp_real_type_RA), PRIVATE :: chp_CNLO(1:7)
  TYPE(chp_real_type), PRIVATE :: chp_CNLO_n(1:7)
!  at the NLO order we also have CIB LO contacts
  TYPE(chp_real_type_RA), PRIVATE :: chp_CIB_CLO(-1:1,1:2)
! NOTE: Ctilde pp 3S1 = Ctilde pn 3S1 = Ctilde nn 3S1
! i.e., all CIB in 1S0 contact
! CIB_CLO(-1,1) : Ctilde pp 1S0  CS pp
! CIB_CLO(-1,2) : Ctilde pp 3S1  CT pp
! CIB_CLO( 0,1) : Ctilde pn 1S0  CS pn
! CIB_CLO( 0,2) : Ctilde pn 3S1  CT pn
! CIB_CLO(+1,1) : Ctilde nn 1S0  CS nn
! CIB_CLO(+1,2) : Ctilde nn 3S1  CT nn
  
! @N3LO       PW
! DN3LO( 1) : D1S0t     D1
! DN3LO( 2) : D1S0      D2
! DN3LO( 3) : D3P0      D3
! DN3LO( 4) : D1P1      D4
! DN3LO( 5) : D3P1      D5
! DN3LO( 6) : D3S1t     D6
! DN3LO( 7) : D3S1      D7
! DN3LO( 8) : D3D1      D8
! DN3LO( 9) : D3S1-3D1t D9
! DN3LO(10) : D3S1-3D1  D10
! DN3LO(11) : D1D2      D11
! DN3LO(12) : D3D2      D12
! DN3LO(13) : D3P2      D13
! DN3LO(14) : D3P2-3F2  D14
! DN3LO(15) : D3D3      D15
  TYPE(chp_real_type_RA), PRIVATE :: chp_DN3LO(1:15)
  TYPE(chp_real_type), PRIVATE :: chp_DN3LO_n(1:15)


! LIST OF OPERATOR STRUCTURES
  
  REAL(8), PRIVATE :: Vc(nof_theta_int_points)
  REAL(8), PRIVATE :: Wc(nof_theta_int_points)
  REAL(8), PRIVATE :: Vs(nof_theta_int_points)
  REAL(8), PRIVATE :: Ws(nof_theta_int_points)
  REAL(8), PRIVATE :: VLS(nof_theta_int_points)
  REAL(8), PRIVATE :: WLS(nof_theta_int_points)
  REAL(8), PRIVATE :: VsigL(nof_theta_int_points)
  REAL(8), PRIVATE :: WsigL(nof_theta_int_points)
  REAL(8), PRIVATE :: VT(nof_theta_int_points)
  REAL(8), PRIVATE :: WT(nof_theta_int_points)
  REAL(8), PRIVATE :: Vsigk(nof_theta_int_points)
  REAL(8), PRIVATE :: Wsigk(nof_theta_int_points)
  
! DERIVED CONSTANTS THAT ARE REPEATEDLY USED, ONCE THE MASSES AND COUPLINGS ARE SET
! THESE ARE COMPUTED IN THE ROUTINE set_units_and_derive_constants.
  
  REAL(8), PRIVATE :: c1,c3,c4
  REAL(8), PRIVATE :: c2, d1_plus_d2, d3, d5, d14_minus_d15
  REAL(8), PRIVATE :: mnuc2(-1:1)   ! nucleon mass squared
  REAL(8), PRIVATE :: mnuc_inv(-1:1)! nucleon mass inversed
  REAL(8), PRIVATE :: mpi2(-1:2)    ! pion mass squared
  REAL(8), PRIVATE :: mpi3(-1:2)    ! pion mass cubed
  REAL(8), PRIVATE :: mpi4(-1:2)    ! mpi2 squared
  REAL(8), PRIVATE :: mpi5(-1:2)    ! mpi^5
  REAL(8), PRIVATE :: twompi(-1:2)  ! two times pion mass
  REAL(8), PRIVATE :: fourmpi2(-1:2)! four times pion mass squared
  
  REAL(8), PRIVATE :: sfr                 ! sfr cutoff
  REAL(8), PRIVATE :: sfr2                ! sfr cutoff squared
  REAL(8), PRIVATE :: sfr_heavyside(-1:2) ! THETA(sfr-twompi)

  REAL(8), PRIVATE :: iso(0:1)      ! expectation value of tau_1 \cdot tau_2
  REAL(8), PRIVATE :: fpi2          ! fpi squared
  REAL(8), PRIVATE :: fpi4          ! fpi2 squared
  REAL(8), PRIVATE :: fpi_inv       ! fpi inverse
  REAL(8), PRIVATE :: gA2           ! gA squared
  REAL(8), PRIVATE :: gA4           ! gA2 squared
  
  REAL(8), PRIVATE :: const(50)      

 INTERFACE chp_print_rconst_if
         MODULE PROCEDURE chp_print_rconst, chp_print_rconst_RA
 END INTERFACE

 CONTAINS
! uncoupled singlet
! POT(1) = <L=J  ,S=0,J|V|L=J  ,S=0,J>
! uncoupled triplet
! POT(2) = <L=J  ,S=1,J|V|L=J  ,S=1,J>
! coupled triplets
! POT(3) = <L=J+1,S=1,J|V|L=J+1,S=0,J>
! POT(4) = <L=J-1,S=1,J|V|L=J-1,S=1,J>
! POT(5) = <L=J+1,S=1,J|V|L=J-1,S=1,J>
! POT(6) = <L=J-1,S=1,J|V|L=J+1,S=1,J>
!
  SUBROUTINE chp(pout, pin, coup, S, j, T, tz, POT)
    
    REAL(8) ,INTENT(IN)    :: pout   ! initial relative momentum
    REAL(8) ,INTENT(IN)    :: pin    ! final relative momentum
    LOGICAL,INTENT(IN)    :: coup   ! IF true, coupled channel
    INTEGER,INTENT(IN)    :: S      ! total spin: 0 or 1
    INTEGER,INTENT(IN)    :: j      ! relative angular momentum
    INTEGER,INTENT(IN)    :: T      ! total isospin: 0 or 1
    INTEGER,INTENT(IN)    :: tz     ! isospin projection: -1, 0, 1 ! pp pn nn
    REAL(8),INTENT(INOUT) :: POT(6) ! potential in LSJ formalism ! MeV
    TYPE(chp_chn_type)    :: CHN    ! coup S, j, T, tz
!REAL(8)            :: vlo(6), vnlo(6), vnnlo(6), vn3lo(6), v_contact(6)
    REAL(8)            :: v_contact(6), v_1PE(6), v_2PE(6)

    
!
! Z=cos(theta), theta angle between p' and p, i.e., final and initial momenta
! q = (p'-p)
! q^2 = p'^2 + p^2 - 2p'pZ
    REAL(8) :: pfinal, pinit
    REAL(8) :: pinit2 ! pinit^2
    REAL(8) :: pfinal2! pfinal^2
    REAL(8) :: q2(1:nof_theta_int_points) ! q^2
    REAL(8) :: q(1:nof_theta_int_points)  ! q
    REAL(8) :: k2(1:nof_theta_int_points) ! k^2
! loop functions
    REAL(8) :: loop_w(nof_theta_int_points)
    REAL(8) :: loop_s
    REAL(8) :: loop_L(nof_theta_int_points)
    REAL(8) :: loop_wtilde2(nof_theta_int_points)
    REAL(8) :: loop_A(nof_theta_int_points)
    REAL(8) :: relcorr_1PE
    INTEGER :: imnuc_2PE
    LOGICAL :: pwd_Vc_2PE_needed, pwd_Wc_2PE_needed, pwd_VLS_2PE_needed, pwd_WLS_2PE_needed
    LOGICAL :: pwd_VsT_2PE_needed, pwd_WsT_2PE_needed, pwd_VsigL_2PE_needed
    LOGICAL :: loop_w_and_L_needed, loop_wtilde_and_A_needed
    LOGICAL :: k2_needed
!

    POT     = 0.0D0
    Vc      = 0.0D0
    Wc      = 0.0D0
    Vs      = 0.0D0
    Ws      = 0.0D0
    VLS     = 0.0D0
    WLS     = 0.0D0
    VsigL   = 0.0D0
    WsigL   = 0.0D0
    VT      = 0.0D0
    WT      = 0.0D0
    Vsigk   = 0.0D0
    Wsigk   = 0.0D0
    loop_w_and_L_needed = .FALSE.
    loop_wtilde_and_A_needed = .FALSE.
    pwd_Vc_2PE_needed = .FALSE.
    pwd_Wc_2PE_needed = .FALSE.
    pwd_VLS_2PE_needed = .FALSE.
    pwd_WLS_2PE_needed = .FALSE.
    pwd_VsT_2PE_needed = .FALSE.
    pwd_WsT_2PE_needed = .FALSE.
    pwd_VsigL_2PE_needed = .FALSE.
    k2_needed = .FALSE.
    
    v_contact = 0.0D0
    v_1PE = 0.0D0
    v_2PE = 0.0D0

    
    loop_w       = 0.0D0
    loop_s       = 0.0D0
    loop_L       = 0.0D0
    loop_wtilde2 = 0.0D0
    loop_A       = 0.0D0

    IF (S==0 .AND. coup) THEN
       WRITE(CHP_ERR,"(A)")
       WRITE(CHP_ERR,"(A)") '   ****************************'
       WRITE(CHP_ERR,"(A)") '   error(chp): illegal channel'
       WRITE(CHP_ERR,"(A)") '   ****************************'
       WRITE(CHP_ERR,"(A)")
       RETURN
    END IF
!
! set the channel quantum numbers
    CHN%j    = j
    CHN%S    = S
    CHN%T    = T
    CHN%tz   = tz
    CHN%coup = coup
!
!DERIVED QUANTITIES
! initial and final momenta squared
    pinit = pin
    pfinal = pout
    pinit2  = pinit*pinit
    pfinal2 = pfinal*pfinal
    
    !IF(pin > 2.0D0*chp_lambda%val .OR. pout > 2.0D0*chp_lambda%val) THEN
    !   POT = 0.0D0
    !   RETURN
    !END IF

! squared momentum transfer
    q2 = pfinal2 + pinit2 - 2.0D0*pinit*pfinal*zmesh%pnt(:)%x
! momentum transfer
    q  = DSQRT(q2)

! Nucleon mass to use for 2PE diagrams
    if(chiral_2PE_CSB_correct_mass%val == 1) then
       imnuc_2PE = tz
    else
       imnuc_2PE = 0
    end if

! LO contact terms
    IF (chiral_Ct_CIB%val == 1) THEN
       CALL chp_CIB_LO_contact_terms(CHN, pinit, pfinal, v_contact)
    ELSE IF (chiral_Ct%val == 1) THEN
       CALL chp_LO_contact_terms(CHN, pinit, pfinal, v_contact)
    END IF

! NLO contact terms
    IF (chiral_C%val == 1) THEN
       CALL chp_NLO_contact_terms(CHN, pinit, pfinal, v_contact)
    END IF

! N3LO contact terms
    IF (chiral_D%val == 1) THEN
       CALL chp_N3LO_contact_terms(CHN, pinit, pfinal, v_contact)
    END IF

! one-pion exchange
    IF (chiral_1PE%val == 1) THEN
       
! relativistic 1/m^2 correction to static 1PE
! (1-(p^2+pp^2)/2mN^2 + ...)
! multiplied inside chp_one_pion_exchange
! Formally, this correction should only enter
! the static 1PE at N3LO.
       IF(chiral_1PE_relcorr%val == 1) THEN
! FIXME: if we use the 2PE isospin breaking mass,
!        then we should use the same for 1PE also.
          relcorr_1PE = (1.0D0-((pfinal2 + pinit2)/(2.0D0*chp_mnuc(0)%val**2)))
       ELSE
          relcorr_1PE=1.0D0
       END IF
    
       IF (chiral_1PE_CIB%val == 1) THEN
! LO CIB one-pion exchanges
! [1] Eq. 4.77-4.79
          IF (tz == 0) THEN
              WT = -1.0D0*chp_one_pion_exchange(q2, relcorr_1PE, 0) + &
                   chp_minus_power(CHN%T+1)*2.0D0*chp_one_pion_exchange(q2, relcorr_1PE, 1)
             if(chiral_1PE_gamma%val == 1) then
                WT = WT + chp_minus_power(CHN%T+1)*2.0D0*chp_one_pion_exchange_gamma(q2, 1)
             end if
          ELSE
             WT = chp_one_pion_exchange(q2, relcorr_1PE, 0)
          END IF
       ELSE
! static OPE contribution
! using average pion mass only
!
          WT = iso(T)*chp_one_pion_exchange(q2, relcorr_1PE, 2)
       END IF
       
       CALL pwd_T(WT, CHN, pfinal, pinit, v_1PE)
       WT = 0.0D0
    END IF

    IF (chiral_2PE_1loop_0%val == 1) THEN
       loop_w_and_L_needed = .TRUE.
    END IF
    IF (chiral_2PE_1loop_d%val == 1  .or. chiral_2PE_1loop_r%val == 1 &
            .or. chiral_2PE_1loop_dd%val == 1 .or. chiral_2PE_1loop_dr%val == 1 &
            .or. chiral_2PE_1loop_rr%val == 1 .or. chiral_2PE_2loop%val == 1) THEN
       loop_w_and_L_needed = .TRUE.
       loop_wtilde_and_A_needed = .TRUE.
    END IF

    IF (chiral_2PE_1loop_rr%val == 1) THEN
       k2_needed = .TRUE.
    END IF

    IF (loop_w_and_L_needed) THEN
       loop_w = chp_NLO_two_pion_exchange_loop_w(q2,2)
       IF (chp_ren%name == 'SF') THEN
          loop_s = chp_NLO_sfr_two_pion_exchange_loop_s(2)
          loop_L = chp_NLO_sfr_two_pion_exchange_loop_L(q, q2, loop_w, loop_s, 2)
       ELSE IF(chp_ren%name == 'DR') THEN
          loop_L = chp_NLO_dr_two_pion_exchange_loop_L(q,loop_w,2)
       END IF
    END IF

    IF (loop_wtilde_and_A_needed) THEN
! NNLO loop functions
       IF (chp_ren%name == 'SF') THEN
          loop_wtilde2 = chp_NNLO_two_pion_exchange_loop_wtilde2(q2,2)
          loop_A       = chp_NNLO_sfr_two_pion_exchange_loop_A(q, q2, 2)
       ELSE IF (chp_ren%name == 'DR') THEN
          loop_wtilde2 = chp_NNLO_two_pion_exchange_loop_wtilde2(q2,2)
          loop_A       = chp_NNLO_dr_two_pion_exchange_loop_A(q,2)
       END IF
    END IF

    IF (k2_needed) THEN
       k2 = 0.25D0 * (pfinal2 + pinit2) + 0.5D0*pinit*pfinal*zmesh%pnt(:)%x
    END IF
          
! All the different 2PE contributions
    IF (chiral_2PE_1loop_0%val == 1) THEN
! irreducible, or non-polynomial, NLO two-pion exchanges
       Wc = Wc + iso(T) * chp_two_pion_exchange_1loop_0_Wc(q2, loop_L, loop_w, 2)
       Vs = Vs +          chp_two_pion_exchange_1loop_0_Vs(q2, loop_L, 2)
       VT = VT +          chp_two_pion_exchange_1loop_0_VT(loop_L, 2)
       pwd_Wc_2PE_needed = .TRUE.
       pwd_VsT_2PE_needed = .TRUE.
    END IF

    IF (chiral_2PE_1loop_d%val == 1) THEN
       Vc  = Vc +          chp_two_pion_exchange_1loop_d_Vc (q2,         loop_A, loop_wtilde2, 2)
       WT  = WT + iso(T) * chp_two_pion_exchange_1loop_d_WT (    loop_w, loop_A                 )
       Ws  = Ws + iso(T) * chp_two_pion_exchange_1loop_d_Ws (q2, loop_w, loop_A                 )
       pwd_Vc_2PE_needed = .TRUE.
       pwd_WsT_2PE_needed = .TRUE.
    END IF
       
    IF (chiral_2PE_1loop_r%val == 1) THEN
       Vc  = Vc  +          chp_two_pion_exchange_1loop_r_Vc (q2, loop_w, loop_A, loop_wtilde2, 2, imnuc_2PE)
       Wc  = Wc  + iso(T) * chp_two_pion_exchange_1loop_r_Wc (q2, loop_w, loop_A, loop_wtilde2, 2, imnuc_2PE)
       VLS = VLS +          chp_two_pion_exchange_1loop_r_VLS(            loop_A, loop_wtilde2,    imnuc_2PE)
       WLS = WLS + iso(T) * chp_two_pion_exchange_1loop_r_WLS(    loop_w, loop_A              ,    imnuc_2PE)
       VT  = VT  +          chp_two_pion_exchange_1loop_r_VT (q2, loop_w, loop_A, loop_wtilde2, 2, imnuc_2PE)
       Vs  = Vs  +          chp_two_pion_exchange_1loop_r_Vs (q2, loop_w, loop_A, loop_wtilde2, 2, imnuc_2PE)
       WT  = WT  + iso(T) * chp_two_pion_exchange_1loop_r_WT (q2, loop_w, loop_A              , 2, imnuc_2PE)
       Ws  = Ws  + iso(T) * chp_two_pion_exchange_1loop_r_Ws (q2, loop_w, loop_A              , 2, imnuc_2PE)
       pwd_Vc_2PE_needed = .TRUE.
       pwd_Wc_2PE_needed = .TRUE.
       pwd_VLS_2PE_needed = .TRUE.
       pwd_WLS_2PE_needed = .TRUE.
       pwd_VsT_2PE_needed = .TRUE.
       pwd_WsT_2PE_needed = .TRUE.
    END IF
       
    IF (chiral_2PE_1loop_dd%val == 1) THEN
       Vc = Vc +          chp_N3LO_2PE_Vc_1loop_ci2(        loop_w, loop_L, loop_wtilde2,         2)
       WT = WT + iso(T) * chp_N3LO_2PE_WT_1loop_ci2(        loop_w, loop_L,                       2)
       Ws = Ws + iso(T) * chp_N3LO_2PE_Ws_1loop_ci2(    q2, loop_w, loop_L,                       2)
       pwd_Vc_2PE_needed = .TRUE.
       pwd_WsT_2PE_needed = .TRUE.
    END IF
       
    IF (chiral_2PE_1loop_dr%val == 1) THEN
       Vc  = Vc  +          chp_N3LO_2PE_Vc_1loop_ci1 (    q2, loop_w, loop_L,                       2, imnuc_2PE)
       Wc  = Wc  + iso(T) * chp_N3LO_2PE_Wc_1loop_ci1 (    q2, loop_w, loop_L,                       2, imnuc_2PE)
       VLS = VLS +          chp_N3LO_2PE_VLS_1loop_ci1(        loop_w, loop_L,                       2, imnuc_2PE)
       WLS = WLS + iso(T) * chp_N3LO_2PE_WLS_1loop_ci1(    q2, loop_w, loop_L,                       2, imnuc_2PE)
       WT  = WT  + iso(T) * chp_N3LO_2PE_WT_1loop_ci1 (    q2, loop_w, loop_L,                       2, imnuc_2PE)
       Ws  = Ws  + iso(T) * chp_N3LO_2PE_Ws_1loop_ci1 (    q2, loop_w, loop_L,                       2, imnuc_2PE)
       pwd_Vc_2PE_needed = .TRUE.
       pwd_Wc_2PE_needed = .TRUE.
       pwd_VLS_2PE_needed = .TRUE.
       pwd_WLS_2PE_needed = .TRUE.
       pwd_WsT_2PE_needed = .TRUE.
    END IF
       
    IF (chiral_2PE_1loop_rr%val == 1) THEN
       Vc  = Vc  +          chp_N3LO_2PE_Vc_1loop_ci0 (    q2, loop_w, loop_L,                       2, imnuc_2PE)
       Wc  = Wc  + iso(T) * chp_N3LO_2PE_Wc_1loop_ci0 (k2, q2, loop_w, loop_L,                       2, imnuc_2PE)
       VLS = VLS +          chp_N3LO_2PE_VLS_1loop_ci0(    q2, loop_w, loop_L,                       2, imnuc_2PE)
       WLS = WLS + iso(T) * chp_N3LO_2PE_WLS_1loop_ci0(    q2, loop_w, loop_L,                       2, imnuc_2PE)
       VT  = VT  +          chp_N3LO_2PE_VT_1loop_ci0 (k2, q2, loop_w, loop_L,                       2, imnuc_2PE)
       Vs  = Vs  +          chp_N3LO_2PE_Vs_1loop_ci0 (k2, q2, loop_w, loop_L,                       2, imnuc_2PE)
       WT  = WT  + iso(T) * chp_N3LO_2PE_WT_1loop_ci0 (    q2, loop_w, loop_L,                       2, imnuc_2PE)
       Ws  = Ws  + iso(T) * chp_N3LO_2PE_Ws_1loop_ci0 (    q2, loop_w, loop_L,                       2, imnuc_2PE)
       VsigL = VsigL +      chp_N3LO_2PE_VsigL_1loop_ci0(              loop_L,                       2, imnuc_2PE)
       pwd_Vc_2PE_needed = .TRUE.
       pwd_Wc_2PE_needed = .TRUE.
       pwd_VLS_2PE_needed = .TRUE.
       pwd_WLS_2PE_needed = .TRUE.
       pwd_VsT_2PE_needed = .TRUE.
       pwd_WsT_2PE_needed = .TRUE.
       pwd_VsigL_2PE_needed = .TRUE.
    END IF
       
    IF (chiral_2PE_2loop%val == 1) THEN
       IF (chp_ren%name == 'DR') THEN
          Vc = Vc +          chp_N3LO_2PE_Vc_2loop_DR  (    q2,                 loop_wtilde2, loop_A, 2)
          Wc = Wc + iso(T) * chp_N3LO_2PE_Wc_2loop_DR_a(    q2, loop_w, loop_L, loop_wtilde2,         2)
          Wc = Wc + iso(T) * chp_N3LO_2PE_Wc_2loop_DR_b(    q2,                                       2)
          VT = VT +          chp_N3LO_2PE_VT_2loop_DR_a(        loop_w, loop_L,                       2)
          VT = VT +          chp_N3LO_2PE_VT_2loop_b   (    q2,                                       2)
          Vs = Vs +          chp_N3LO_2PE_Vs_2loop_DR_a(    q2, loop_w, loop_L,                       2)
          Vs = Vs +          chp_N3LO_2PE_Vs_2loop_b   (    q2,                                       2)
          WT = WT + iso(T) * chp_N3LO_2PE_WT_2loop_DR  (        loop_w,                       loop_A, 2)
          Ws = Ws + iso(T) * chp_N3LO_2PE_Ws_2loop_DR  (    q2, loop_w,                       loop_A, 2)
       ELSE IF (chp_ren%name == 'SF') THEN
          Vc = Vc +          chp_N3LO_2PE_Vc_2loop_SFR  (q2, 2)
          Wc = Wc + iso(T) * chp_N3LO_2PE_Wc_2loop_SFR_a(q2, 2)
          Wc = Wc + iso(T) * chp_N3LO_2PE_Wc_2loop_SFR_b(q2, 2)
          VT = VT +          chp_N3LO_2PE_VT_2loop_SFR_a(q2, 2)
          VT = VT +          chp_N3LO_2PE_VT_2loop_b    (q2, 2)
          Vs = Vs +          chp_N3LO_2PE_Vs_2loop_SFR_a(q2, 2)
          Vs = Vs +          chp_N3LO_2PE_Vs_2loop_b    (q2, 2)
          WT = WT + iso(T) * chp_N3LO_2PE_WT_2loop_SFR  (q2, 2)
          Ws = Ws + iso(T) * chp_N3LO_2PE_Ws_2loop_SFR  (q2, 2)
       ELSE
           STOP 'Unknown regulator form, need DR or SF'
       END IF
       pwd_Vc_2PE_needed = .TRUE.
       pwd_Wc_2PE_needed = .TRUE.
       pwd_VsT_2PE_needed = .TRUE.
       pwd_WsT_2PE_needed = .TRUE.
    END IF
       
! project to LSJ basis
    IF (pwd_Vc_2PE_needed) THEN
       CALL pwd_c(Vc, CHN, pfinal, pinit, v_2PE)
    END IF
    IF (pwd_Wc_2PE_needed) THEN
       CALL pwd_c(Wc, CHN, pfinal, pinit, v_2PE)
    END IF
    IF (pwd_VLS_2PE_needed) THEN
       CALL pwd_LS(VLS, CHN, pfinal, pinit, v_2PE)
    END IF
    IF (pwd_WLS_2PE_needed) THEN
       CALL pwd_LS(WLS, CHN, pfinal, pinit, v_2PE)
    END IF
    IF (pwd_VsT_2PE_needed) THEN
       CALL pwd_T(VT, CHN, pfinal, pinit, v_2PE)
       CALL pwd_s(Vs, CHN, pfinal, pinit, v_2PE)
    END IF
    IF (pwd_WsT_2PE_needed) THEN
       CALL pwd_T(WT, CHN, pfinal, pinit, v_2PE)
       CALL pwd_s(Ws, CHN, pfinal, pinit, v_2PE)
    END IF
    IF (pwd_VsigL_2PE_needed) THEN
       CALL pwd_sigL(VsigL, CHN, pfinal, pinit, v_2PE)
    END IF


    
    IF (chiral_minimal_relativity%val == 1) then
! kinematical considerations
! - adopt minimal relativity by default
       MR_factor = chp_mnuc(CHN%tz)%val/DSQRT(cm_energy(pinit,CHN%tz)*cm_energy(pfinal,CHN%tz))
! - unless it is explicitly switched off
    ELSE
       MR_factor = 1.0D0
    END IF
    
    if(chiral_kamada_glockle_transform%val == 1) then
       if (chiral_minimal_relativity%val == 1) stop 'cannot use minimal relativity and kamada-glockle transform simultaneously'
       MR_factor = SQRT((1.0D0+pinit2/(2.0D0*chp_mnuc(CHN%tz)%val**2))*SQRT(1.0D0+pinit2**2/(4.0D0*chp_mnuc(CHN%tz)%val**2))) * &
                   SQRT((1.0D0+pfinal2/(2.0D0*chp_mnuc(CHN%tz)%val**2))*SQRT(1.0D0+pfinal2/(4.0D0*chp_mnuc(CHN%tz)%val**2)))
    end if
    
    v_1PE = v_1PE * const(1) * MR_factor
    v_2PE = v_2PE * const(1) * MR_factor
    v_contact = v_contact * const(1) * MR_factor
    
! add the contributions computed above
! and regulate them according to their chiral order
    POT = v_1PE*freg(pfinal,pinit,chp_regcut_1PE%val) + &
          v_2PE*freg(pfinal,pinit,chp_regcut_2PE%val) + &
          v_contact

    
  END SUBROUTINE chp

  SUBROUTINE chp_orig(pout, pin, coup, S, j, T, tz, POT)
    
    REAL(8),INTENT(IN)    :: pout   ! initial relative momentum
    REAL(8),INTENT(IN)    :: pin    ! final relative momentum
    LOGICAL,INTENT(IN)    :: coup   ! IF true, coupled channel
    INTEGER,INTENT(IN)    :: S      ! total spin: 0 or 1
    INTEGER,INTENT(IN)    :: j      ! relative angular momentum
    INTEGER,INTENT(IN)    :: T      ! total isospin: 0 or 1
    INTEGER,INTENT(IN)    :: tz     ! isospin projection: -1, 0, 1 ! pp pn nn
    REAL(8),INTENT(INOUT) :: POT(6) ! potential in LSJ formalism ! MeV

    CALL chp(pout, pin, coup, S, j, T, tz, POT)

  END SUBROUTINE chp_orig
   
! z = 2*m_\pi / \mu, using the notation of [1] and [2] for the 2PE 2-loop diagrams
! In DR, the integration is done from z_low = 0 and in SFR from z_low = 2 * m_\pi / \Lambda_{SFR}
! The expressions for some of the 2PE 2-loop diagrams contains integrals that must be solved numerically
! This function sets up the integration meshs used in the numerical integrations
  SUBROUTINE chp_setup_2PE_2loop_int_data(nof_x_points, nof_z_points, z_low)
    INTEGER, INTENT(IN) :: nof_x_points
    INTEGER, INTENT(IN) :: nof_z_points
    REAL(8), INTENT(IN) :: z_low

    chp_2PE_2loop_int_mesh_x%info = chp_mesh_info(nof_x_points, 0.0D0, 1.0D0)
    CALL chp_setup_gauleg_mesh(chp_2PE_2loop_int_mesh_x)

    chp_2PE_2loop_int_mesh_z%info = chp_mesh_info(nof_z_points, z_low, 1.0D0)
    CALL chp_setup_gauleg_mesh(chp_2PE_2loop_int_mesh_z)

! TODO: Magic number: 0.95. The reason for using 0.95 is that the z integral is a bit
! tricky near z = 1, therefore I do numerical integration up to 0.95, then I use the
! integral of the taylor expansion around z=1 from 0.95 to 1.0 to get the last part.
! This is needed to get accurate results.
    chp_2PE_2loop_int_mesh_z_VC_SFR%info = chp_mesh_info(nof_z_points, z_low, 0.95D0)
    CALL chp_setup_gauleg_mesh(chp_2PE_2loop_int_mesh_z_VC_SFR)

    chp_2PE_2loop_int_VTS_data_set     = .false.
    chp_2PE_2loop_int_WC_DR_data_set   = .false.
    chp_2PE_2loop_int_WC_SFR_data_set  = .false.
  END SUBROUTINE

! static one pion exchange, [1] 4.5
! without isospin structure
! q2  : momentum transfer squared
! impi: determines which mpi2 to use,
  FUNCTION chp_one_pion_exchange(q2, relcorr_1PE, impi) RESULT(res)
    
    REAL(8) , INTENT(IN) :: q2(1:nof_theta_int_points)
    REAL(8) , INTENT(IN) :: relcorr_1PE
    INTEGER, INTENT(IN) :: impi
    REAL(8) :: res(1:nof_theta_int_points)
    
    res = -1.0D0 * const(2)/ (q2 + mpi2(impi))
    
! This adds the relativistic corrections if they are included, otherwise relcorr_1PE = 1
    res = res*relcorr_1PE
    
  END FUNCTION chp_one_pion_exchange
  
! FIXME: Add reference for this expression
  FUNCTION chp_one_pion_exchange_gamma(q2, impi) RESULT(res)
    
    REAL(8) , INTENT(IN) :: q2(1:nof_theta_int_points)
    REAL(8) :: res(1:nof_theta_int_points)
    REAL(8) :: b2(1:nof_theta_int_points)
    INTEGER, INTENT(IN) :: impi
    
    b2 = q2 / mpi2(impi)
    res = -(0.25d0 * chp_fine_structure%val * gA2 / (chp_pi%val * fpi2 * mpi2(impi))) * (-((1-b2)**2 / (2.0d0 * b2**2 * (1 + b2))) * dlog(b2+1) + 1.0d0 / (2*b2))

  END FUNCTION chp_one_pion_exchange_gamma
  
! NLO loop function w [1] Eq 4.12 (DR and SFR)
! q2  : momentum transfer squared
! impi: determines which mpi2 to use,
  FUNCTION chp_NLO_two_pion_exchange_loop_w(q2, impi) RESULT(res)
    
    REAL(8) , INTENT(IN) :: q2(1:nof_theta_int_points)
    INTEGER, INTENT(IN) :: impi
    REAL(8) :: res(1:nof_theta_int_points)
    
    res = DSQRT(fourmpi2(impi) + q2)
    
  END FUNCTION chp_NLO_two_pion_exchange_loop_w

! NLO SFR loop function s [2] Eq 2.16
! impi: determines which mpi2 to use,
  FUNCTION chp_NLO_sfr_two_pion_exchange_loop_s(impi) RESULT(res)
    
    INTEGER, INTENT(IN) :: impi
    REAL(8) :: res
    
    res = DSQRT(sfr2 - fourmpi2(impi))
    
  END FUNCTION chp_NLO_sfr_two_pion_exchange_loop_s

! NLO dr loop function L [1] Eq 4.11
! with isospin structure
! q2  : momentum transfer squared
! impi: determines which mpi2 to use,
  FUNCTION chp_NLO_dr_two_pion_exchange_loop_L(q, w, impi) RESULT(res)
    
    REAL(8) , INTENT(IN) :: q(1:nof_theta_int_points)
    REAL(8) , INTENT(IN) :: w(1:nof_theta_int_points)
    INTEGER, INTENT(IN) :: impi
    REAL(8) :: res(1:nof_theta_int_points)
    
    res = w * dlog((w+q)/(twompi(impi) ) )/q
    
  END FUNCTION chp_NLO_dr_two_pion_exchange_loop_L

! NLO SFR loop function L [2] Eq 2.16
! q   : momentum transfer
! q2  : momentum transfer squared
! w,s : SFR loop functions
! impi: determines which mpi2 to use,
  FUNCTION chp_NLO_sfr_two_pion_exchange_loop_L(q, q2, w, s, impi) RESULT(res)
    
    REAL(8) , INTENT(IN) :: q(1:nof_theta_int_points)
    REAL(8) , INTENT(IN) :: q2(1:nof_theta_int_points)
    REAL(8) , INTENT(IN) :: w(1:nof_theta_int_points)
    REAL(8) , INTENT(IN) :: s
    INTEGER, INTENT(IN) :: impi
    REAL(8) :: res(1:nof_theta_int_points)
    
    res = 0.0D0
    IF (sfr_heavyside(impi) == 0.0D0) return
    
    res = w * dlog( (sfr2*w*w+q2*s*s+2.0D0*sfr*q*w*s)/( fourmpi2(impi)*(sfr2+q2) ) )/(2.0D0*q)
    
  END FUNCTION chp_NLO_sfr_two_pion_exchange_loop_L
  
! NLO function [1] 4.9 OR [2] 2.14 (W_C part) OR [3] B1 OR [4] 6
! q2  : momentum transfer squared
! impi: determines which mpi2 to use,
  FUNCTION chp_two_pion_exchange_1loop_0_Wc(q2, L, w, impi) RESULT(res)
    
    REAL(8) , INTENT(IN) :: q2(1:nof_theta_int_points)
    REAL(8) , INTENT(IN) :: L(1:nof_theta_int_points)
    REAL(8) , INTENT(IN) :: w(1:nof_theta_int_points)
    INTEGER, INTENT(IN) :: impi
    REAL(8) :: res(1:nof_theta_int_points)
    
    res = -L *(fourmpi2(impi)*const(4) + q2*const(5) + const(6)*mpi4(impi)/(w*w))/const(3)
    
  END FUNCTION chp_two_pion_exchange_1loop_0_Wc
  
! NLO function [1] 4.10 OR [2] 2.14 (V_S part) OR [3] B2 OR [4] 6
! q2  : momentum transfer squared
! impi: determines which mpi2 to use,
  FUNCTION chp_two_pion_exchange_1loop_0_Vs(q2,L,impi) RESULT(res)
    
    REAL(8) , INTENT(IN) :: q2(1:nof_theta_int_points)
    REAL(8) , INTENT(IN) :: L(1:nof_theta_int_points)
    INTEGER, INTENT(IN) :: impi
    REAL(8) :: res(1:nof_theta_int_points)
    
    res = const(7)*L*q2/const(8)
    
  END FUNCTION chp_two_pion_exchange_1loop_0_Vs

! NLO function [1] 4.10 OR [2] 2.14 (V_T part) OR [3] B2 OR [4] 6
! q2  : momentum transfer squared
! impi: determines which mpi2 to use,
  FUNCTION chp_two_pion_exchange_1loop_0_VT(L,impi) RESULT(res)
    
    REAL(8) , INTENT(IN) :: L(1:nof_theta_int_points)
    INTEGER, INTENT(IN) :: impi
    REAL(8) :: res(1:nof_theta_int_points)
    
    res = -const(7)*L/const(8)
    
  END FUNCTION chp_two_pion_exchange_1loop_0_VT

! NNLO loop function wtilde SQUARED [1] Eq 4.20 (DR)
! q2  : momentum transfer squared
! impi: determines which mpi2 to use,
  FUNCTION chp_NNLO_two_pion_exchange_loop_wtilde2(q2, impi) RESULT(res)
    
    REAL(8) , INTENT(IN) :: q2(1:nof_theta_int_points)
    INTEGER, INTENT(IN) :: impi
    REAL(8) :: res(1:nof_theta_int_points)
    
    res = 0.0D0
    res = 2.0D0*mpi2(impi) + q2
    
  END FUNCTION chp_NNLO_two_pion_exchange_loop_wtilde2

! NNLO loop function wtilde [1] Eq 4.19 (DR)
! q   : momentum transfer
! impi: determines which mpi to use,
  FUNCTION chp_NNLO_dr_two_pion_exchange_loop_A(q, impi) RESULT(res)
    
    REAL(8) , INTENT(IN) :: q(1:nof_theta_int_points)
    INTEGER, INTENT(IN) :: impi
    REAL(8) :: res(1:nof_theta_int_points)
    
    res = 0.0D0
    res = datan(q/(twompi(impi)))/(2.0D0*q)
    
  END FUNCTION chp_NNLO_dr_two_pion_exchange_loop_A

! NNLO loop function wtilde [2] Eq 2.17 (SFR)
! q   : momentum transfer
! q2  : momentum transfer squared
! impi: determines which mpi to use,
  FUNCTION chp_NNLO_sfr_two_pion_exchange_loop_A(q, q2, impi) RESULT(res)
    
    REAL(8) , INTENT(IN) :: q(1:nof_theta_int_points)
    REAL(8) , INTENT(IN) :: q2(1:nof_theta_int_points)
    INTEGER, INTENT(IN) :: impi
    REAL(8) :: res(1:nof_theta_int_points)
    
    res = 0.0D0
    IF (sfr_heavyside(impi) == 0.0D0) return

    res = datan( q*(sfr-twompi(impi) )/(q2 + sfr*twompi(impi) ) )/(2.0D0*q)
    
  END FUNCTION chp_NNLO_sfr_two_pion_exchange_loop_A

! NNLO function [1] Eq 4.13 (ci part) OR [2] Eq 2.15 (V_C part) OR [3] Eq C1
! q2    : momentum transfer squared
! A     : NNLO loop function A [1]Eq. 4.19
! wt2   : NNLO loop function wtilde Squared(!) ([1] Eq. 4.20)
! impi  : determines which mpi2 to use
  FUNCTION chp_two_pion_exchange_1loop_d_Vc(q2, A, wt2, impi) RESULT(res)
    
    REAL(8) , INTENT(IN) :: q2(1:nof_theta_int_points)
    REAL(8) , INTENT(IN) :: A(1:nof_theta_int_points)
    REAL(8) , INTENT(IN) :: wt2(1:nof_theta_int_points)
    INTEGER, INTENT(IN) :: impi
    REAL(8) :: res(1:nof_theta_int_points)
    
    res = const(9)*(- &
            (mpi2(impi)*const(11) - q2*c3) *wt2*A)

  END FUNCTION chp_two_pion_exchange_1loop_d_Vc

! NNLO/N3LO function [1] Eq 4.13 (M_N part) and Eq 4.21 OR [2] Eq 2.23 (V_C part) OR [3] Eq D7
! q2    : momentum transfer squared
! w     : Eq 4.12
! A     : NNLO loop function A Eq. 4.19
! wt2   : NNLO loop function wtilde Squared(!) (Eq. 4.20)
! impi  : determines which mpi2 to use
! chp_itope;  EM or KW, choose which iterated 2pe 'model' to use
! this is set in the module header
  FUNCTION chp_two_pion_exchange_1loop_r_Vc(q2, w, A, wt2, impi, imnuc) RESULT(res)
    
    REAL(8) , INTENT(IN) :: q2(1:nof_theta_int_points)
    REAL(8) , INTENT(IN) :: w(1:nof_theta_int_points)
    REAL(8) , INTENT(IN) :: A(1:nof_theta_int_points)
    REAL(8) , INTENT(IN) :: wt2(1:nof_theta_int_points)
    INTEGER, INTENT(IN) :: impi
    INTEGER, INTENT(IN) :: imnuc
    REAL(8) :: res(1:nof_theta_int_points)
    
    if(chiral_2PE_1loop_r_mode%val == chiral_2PE_1loop_r_mode_EM2011) then
! [1] 4.13
       res = const(9)*const(10)*(mpi5(impi)/(w*w) + q2*3.0D0*wt2*A) / chp_mnuc(imnuc)%val
! [1] 4.21
       res = res - const(12)*(chp_mpi(impi)%val*w*w+wt2*wt2*A)/chp_mnuc(imnuc)%val
    elseif(chiral_2PE_1loop_r_mode%val == chiral_2PE_1loop_r_mode_Ep2005) then
! [2] 2.23 OR [4] 19
       res = const(15)*(2*mpi5(impi)/(w*w) - 3*(4*mpi4(impi) - q2*q2)*A) / chp_mnuc(imnuc)%val
    elseif(chiral_2PE_1loop_r_mode%val == chiral_2PE_1loop_r_mode_EM2015) then
! [3] D7 OR [4] 19+22
       res = 3*const(14)*(mpi5(impi)/(2*w*w) + (2*mpi2(impi) + q2) * (q2 - mpi2(impi))*A) / chp_mnuc(imnuc)%val
    elseif(chiral_2PE_1loop_r_mode%val == chiral_2PE_1loop_r_mode_Ep2015) then
! [4] 19
       res = const(15)*((2.0D0*mpi5(impi)/(w**2)) -3.0D0*(4.0D0*mpi4(impi) - q2**2)*A)/chp_mnuc(imnuc)%val
! [4] 22
       res = res + const(15)*(2.0D0*mpi2(impi) + q2)**2*A/chp_mnuc(imnuc)%val
    else
       error stop 'Unknown chiral_2PE_1loop_r_mode'
    end if

  END FUNCTION chp_two_pion_exchange_1loop_r_Vc

! NNLO function Eq [1] 4.14 and 4.22
! q2    : momentum transfer squared
! w     : Eq 4.12
! A     : NNLO loop function A Eq. 4.19
! wt2   : NNLO loop function wtilde Squared(!) (Eq. 4.20)
! impi  : determines which mpi2 to use
! chp_itope;  EM or KW, choose which iterated 2pe 'model' to use
! this is set in the module header
  FUNCTION chp_two_pion_exchange_1loop_r_Wc(q2, w, A, wt2, impi, imnuc) RESULT(res)
    
    REAL(8) , INTENT(IN) :: q2(1:nof_theta_int_points)
    REAL(8) , INTENT(IN) :: w(1:nof_theta_int_points)
    REAL(8) , INTENT(IN) :: A(1:nof_theta_int_points)
    REAL(8) , INTENT(IN) :: wt2(1:nof_theta_int_points)
    INTEGER, INTENT(IN) :: impi
    INTEGER, INTENT(IN) :: imnuc
    REAL(8) :: res(1:nof_theta_int_points)
    
    if(chiral_2PE_1loop_r_mode%val == chiral_2PE_1loop_r_mode_EM2011) then
! [1] 4.14
       res = const(13) * (3.0D0*gA2*mpi5(impi)/(w*w) - &
             (fourmpi2(impi) + 2.0D0*q2 - gA2*(fourmpi2(impi)+3.0D0*q2))*wt2*A)/chp_mnuc(imnuc)%val
! [1] 4.22
       res = res + const(14)*(chp_mpi(impi)%val*w*w + wt2*wt2*A)/chp_mnuc(imnuc)%val
    elseif(chiral_2PE_1loop_r_mode%val == chiral_2PE_1loop_r_mode_Ep2005) then
        error stop 'Not implemented'
    elseif(chiral_2PE_1loop_r_mode%val == chiral_2PE_1loop_r_mode_EM2015) then
![3] D8 = [4] 19 + 22
       res = 2.0D0*const(13) * (3.0D0*gA2*mpi5(impi)/(2.0D0*w*w) + &
             (gA2*(3.0D0*mpi2(impi) + 2.0D0*q2) - 2.0D0*mpi2(impi) - q2)*&
             (2.0D0*mpi2(impi) + q2)*A)/chp_mnuc(imnuc)%val
    elseif(chiral_2PE_1loop_r_mode%val == chiral_2PE_1loop_r_mode_Ep2015) then
![4] 19
       res =  const(13) * ((3.0D0*gA2*mpi5(impi)/(w*w)) - &
            (4.0D0*mpi2(impi) + 2.0D0*q2 - gA2*(7.0D0*mpi2(impi) + 4.5D0*q2))*&
            (2.0D0*mpi2(impi) + q2)*A)/chp_mnuc(imnuc)%val
       
![4] 22
       res = res - const(18)*(2.0D0*mpi2(impi) + q2)**2*A/chp_mnuc(imnuc)%val
    else
       error stop 'Unknown chiral_2PE_1loop_r_mode'
    end if
    
  END FUNCTION chp_two_pion_exchange_1loop_r_Wc
  
! NNLO function Eq [1] 4.15 and 4.23
! w     : Eq 4.12
! A     : NNLO loop function A Eq. 4.19
! wt2   : NNLO loop function wtilde Squared(!) (Eq. 4.20)
! impi  : determines which mpi2 to use
! chp_itope;  EM or KW, choose which iterated 2pe 'model' to use
! this is set in the module header
  FUNCTION chp_two_pion_exchange_1loop_r_VT(q2, w, A, wt2, impi, imnuc) RESULT(res)

    REAL(8) , INTENT(IN) :: q2(1:nof_theta_int_points)
    REAL(8) , INTENT(IN) :: w(1:nof_theta_int_points)
    REAL(8) , INTENT(IN) :: A(1:nof_theta_int_points)
    REAL(8) , INTENT(IN) :: wt2(1:nof_theta_int_points)
    INTEGER, INTENT(IN) :: impi
    INTEGER, INTENT(IN) :: imnuc
    REAL(8) :: res(1:nof_theta_int_points)
    
    if(chiral_2PE_1loop_r_mode%val == chiral_2PE_1loop_r_mode_EM2011) then
! [1] 4.15
       res = 3.0D0*const(15)*wt2*A/chp_mnuc(imnuc)%val
! [1] 4.23
       res = res + const(15)*(chp_mpi(impi)%val + w*w*A )/chp_mnuc(imnuc)%val
    elseif(chiral_2PE_1loop_r_mode%val == chiral_2PE_1loop_r_mode_Ep2005) then
        error stop 'Not implemented'
    elseif(chiral_2PE_1loop_r_mode%val == chiral_2PE_1loop_r_mode_EM2015) then
![3] D9 = [4] 19 + 22
       res = 3.0D0*const(18)*(5.0D0*mpi2(impi) + 2.0D0*q2)*A/chp_mnuc(imnuc)%val
    elseif(chiral_2PE_1loop_r_mode%val == chiral_2PE_1loop_r_mode_Ep2015) then
![4] 19
       res = 3.0D0*const(15)*(4.0D0*mpi2(impi) + 1.5D0*q2)*A/chp_mnuc(imnuc)%val
![4] 22
       res = res -0.5D0*const(15)*(4.0D0*mpi2(impi) + q2)*A/chp_mnuc(imnuc)%val
    else
       error stop 'Unknown chiral_2PE_1loop_r_mode'
    end if
    
  END FUNCTION chp_two_pion_exchange_1loop_r_VT

! NNLO function Eq [1] 4.15 and 4.23
! q2    : momentum transfer squared
! w     : Eq 4.12
! A     : NNLO loop function A Eq. 4.19
! wt2   : NNLO loop function wtilde Squared(!) (Eq. 4.20)
! impi  : determines which mpi2 to use
! chp_itope;  EM or KW, choose which iterated 2pe 'model' to use
! this is set in the module header
  FUNCTION chp_two_pion_exchange_1loop_r_Vs(q2, w, A, wt2, impi, imnuc) RESULT(res)
    
    REAL(8) , INTENT(IN) :: q2(1:nof_theta_int_points)
    REAL(8) , INTENT(IN) :: w(1:nof_theta_int_points)
    REAL(8) , INTENT(IN) :: A(1:nof_theta_int_points)
    REAL(8) , INTENT(IN) :: wt2(1:nof_theta_int_points)
    INTEGER, INTENT(IN) :: impi
    INTEGER, INTENT(IN) :: imnuc
    REAL(8) :: res(1:nof_theta_int_points)
    
    res = -q2 * chp_two_pion_exchange_1loop_r_VT(q2, w, A, wt2, impi, imnuc)

  END FUNCTION chp_two_pion_exchange_1loop_r_Vs
  
! NNLO function Eq [1] 4.16
! q2    : momentum transfer squared
! w     : Eq 4.12
! A     : NNLO loop function A Eq. 4.19
! impi  : determines which mpi2 to use
! chp_itope;  EM or KW, choose which iterated 2pe 'model' to use
! this is set in the module header
  FUNCTION chp_two_pion_exchange_1loop_d_WT(w, A) RESULT(res)
    
    REAL(8) , INTENT(IN) :: w(1:nof_theta_int_points)
    REAL(8) , INTENT(IN) :: A(1:nof_theta_int_points)
    REAL(8) :: res(1:nof_theta_int_points)
    
! [1] 4.16
    res = -1.0D0*const(16)*A*( c4*w*w )
    
  END FUNCTION chp_two_pion_exchange_1loop_d_WT

! NNLO function Eq [1] 4.16 and 4.24
! q2    : momentum transfer squared
! w     : Eq 4.12
! A     : NNLO loop function A Eq. 4.19
! impi  : determines which mpi2 to use
! chp_itope;  EM or KW, choose which iterated 2pe 'model' to use
! this is set in the module header
  FUNCTION chp_two_pion_exchange_1loop_r_WT(q2, w, A, impi, imnuc) RESULT(res)
    
    REAL(8) , INTENT(IN) :: q2(1:nof_theta_int_points)
    REAL(8) , INTENT(IN) :: w(1:nof_theta_int_points)
    REAL(8) , INTENT(IN) :: A(1:nof_theta_int_points)
    INTEGER, INTENT(IN) :: impi
    INTEGER, INTENT(IN) :: imnuc
    REAL(8) :: res(1:nof_theta_int_points)
    
    if(chiral_2PE_1loop_r_mode%val == chiral_2PE_1loop_r_mode_EM2011) then
! [1] 4.16
       res = -1.0D0*const(16)*A*( (1.0D0/(4.0D0*chp_mnuc(imnuc)%val))*w*w - &
             const(17)*(10.0D0*mpi2(impi) + 3.0D0*q2)/chp_mnuc(imnuc)%val)
! [1] 4.24
       res = res - const(18)*(chp_mpi(impi)%val + w*w*A)/chp_mnuc(imnuc)%val
    elseif(chiral_2PE_1loop_r_mode%val == chiral_2PE_1loop_r_mode_Ep2005) then
        error stop 'Not implemented'
    elseif(chiral_2PE_1loop_r_mode%val == chiral_2PE_1loop_r_mode_EM2015) then
! [3] D10 = [4] 19 + 22
       res = const(13)*(gA2*(3.0D0*mpi2(impi) + q2) - w**2)*A/chp_mnuc(imnuc)%val
    elseif(chiral_2PE_1loop_r_mode%val == chiral_2PE_1loop_r_mode_Ep2015) then
![4] 19
       res = -0.5D0*const(13)*(8.0D0*mpi2(impi) + 2.0D0*q2 - &
            gA2*(4.0D0*mpi2(impi) + 1.5D0*q2))*A/chp_mnuc(imnuc)%val
![4] 22
       res = res + 0.5D0*const(18)*(4.0D0*mpi2(impi) + q2)*A/chp_mnuc(imnuc)%val
    else
       error stop 'Unknown chiral_2PE_1loop_r_mode'
    end if
    
  END FUNCTION chp_two_pion_exchange_1loop_r_WT

! NNLO function Eq [1] 4.16
! q2    : momentum transfer squared
! w     : Eq 4.12
! A     : NNLO loop function A Eq. 4.19
! impi  : determines which mpi2 to use
! chp_itope;  EM or KW, choose which iterated 2pe 'model' to use
! this is set in the module header
  FUNCTION chp_two_pion_exchange_1loop_d_Ws(q2, w, A) RESULT(res)
    
    REAL(8) , INTENT(IN) :: q2(1:nof_theta_int_points)
    REAL(8) , INTENT(IN) :: w(1:nof_theta_int_points)
    REAL(8) , INTENT(IN) :: A(1:nof_theta_int_points)
    REAL(8) :: res(1:nof_theta_int_points)
    
    res = -q2 * chp_two_pion_exchange_1loop_d_WT(w, A)
    
  END FUNCTION chp_two_pion_exchange_1loop_d_Ws

! NNLO function Eq [1] 4.16 and 4.24
! q2    : momentum transfer squared
! w     : Eq 4.12
! A     : NNLO loop function A Eq. 4.19
! impi  : determines which mpi2 to use
! chp_itope;  EM or KW, choose which iterated 2pe 'model' to use
! this is set in the module header
  FUNCTION chp_two_pion_exchange_1loop_r_Ws(q2, w, A, impi, imnuc) RESULT(res)
    
    REAL(8) , INTENT(IN) :: q2(1:nof_theta_int_points)
    REAL(8) , INTENT(IN) :: w(1:nof_theta_int_points)
    REAL(8) , INTENT(IN) :: A(1:nof_theta_int_points)
    INTEGER, INTENT(IN) :: impi
    INTEGER, INTENT(IN) :: imnuc
    REAL(8) :: res(1:nof_theta_int_points)
    
    res = -q2 * chp_two_pion_exchange_1loop_r_WT(q2, w, A, impi, imnuc)
    
  END FUNCTION chp_two_pion_exchange_1loop_r_Ws

! NNLO function Eq [1] 4.17
! A     : NNLO loop function A Eq. 4.19
! wt2   : NNLO loop function wtilde Squared(!) (Eq. 4.20)
! impi  : determines which mpi2 to use
  FUNCTION chp_two_pion_exchange_1loop_r_VLS(A, wt2, imnuc) RESULT(res)
    
    REAL(8) , INTENT(IN) :: A(1:nof_theta_int_points)
    REAL(8) , INTENT(IN) :: wt2(1:nof_theta_int_points)
    INTEGER, INTENT(IN) :: imnuc
    REAL(8) :: res(1:nof_theta_int_points)
    
    res = const(19) * wt2*A/chp_mnuc(imnuc)%val
             
  END FUNCTION chp_two_pion_exchange_1loop_r_VLS

! NNLO function Eq [1] 4.18
! w   : NNLO loop function wtilde Squared(!) (Eq. 4.20)
! A     : NNLO loop function A Eq. 4.19
! impi  : determines which mpi2 to use
  FUNCTION chp_two_pion_exchange_1loop_r_WLS(w, A, imnuc) RESULT(res)
    
    REAL(8) , INTENT(IN) :: w(1:nof_theta_int_points)
    REAL(8) , INTENT(IN) :: A(1:nof_theta_int_points)
    INTEGER, INTENT(IN) :: imnuc
    REAL(8) :: res(1:nof_theta_int_points)
    
    res = const(20)*w*w*A/chp_mnuc(imnuc)%val
             
  END FUNCTION chp_two_pion_exchange_1loop_r_WLS

! N3LO two pion exchange functions
! k2    : mean momentum squared
! q2    : momentum transfer squared
! w     : Eq 4.12
! L     : L loop function (depends on what regularization is used)
! wt2   : NNLO loop function wtilde Squared(!) (Eq. 4.20)
! A     : NNLO loop function A Eq. 4.19
! impi  : determines which mpi2 to use

! [1]Eq. D.1
  FUNCTION chp_N3LO_2PE_Vc_1loop_ci2 (        w, L, wt2,    impi) RESULT(res)
    REAL(8) , INTENT(IN) :: w(1:nof_theta_int_points)
    REAL(8) , INTENT(IN) :: L(1:nof_theta_int_points)
    REAL(8) , INTENT(IN) :: wt2(1:nof_theta_int_points)
    INTEGER, INTENT(IN) :: impi
    REAL(8) :: res(1:nof_theta_int_points)
    
    res = c2 * (1D0/6D0) * w * w + c3 * wt2 - c1 * fourmpi2(impi)
    res = ((3D0 / 16D0) * const(21)) * L * (res*res + c2*c2*w*w*w*w*(1D0/45D0))
  END FUNCTION chp_N3LO_2PE_Vc_1loop_ci2

! [1]Eq. D.4
  FUNCTION chp_N3LO_2PE_Vc_1loop_ci1 (    q2, w, L,             impi, imnuc) RESULT(res)
    REAL(8) , INTENT(IN) :: q2(1:nof_theta_int_points)
    REAL(8) , INTENT(IN) :: w(1:nof_theta_int_points)
    REAL(8) , INTENT(IN) :: L(1:nof_theta_int_points)
    INTEGER, INTENT(IN) :: impi
    INTEGER, INTENT(IN) :: imnuc
    REAL(8) :: res(1:nof_theta_int_points)
    
    res = (-gA2 * const(21) * mnuc_inv(imnuc) * (1D0 / 32D0)) * L * ( &
              (c2 - 6D0*c3) * q2 * q2 + 4D0 * (6D0*c1 + c2 - 3D0*c3)*q2*mpi2(impi) + &
              6D0 * (c2 - 2D0*c3)*mpi4(impi) + &
              24D0*(2D0*c1 + c3)*mpi4(impi)*mpi2(impi)/(w*w) )

  END FUNCTION chp_N3LO_2PE_Vc_1loop_ci1

! [1]Eq. D.9
  FUNCTION chp_N3LO_2PE_Vc_1loop_ci0 (    q2, w, L,             impi, imnuc) RESULT(res)
    REAL(8) , INTENT(IN) :: q2(1:nof_theta_int_points)
    REAL(8) , INTENT(IN) :: w(1:nof_theta_int_points)
    REAL(8) , INTENT(IN) :: L(1:nof_theta_int_points)
    INTEGER, INTENT(IN) :: impi
    INTEGER, INTENT(IN) :: imnuc
    REAL(8) :: res(1:nof_theta_int_points)
    
    res = -gA2*gA2*const(21)*mnuc_inv(imnuc)*mnuc_inv(imnuc)*(1D0/32D0) * ( &
              L*(2D0*mpi4(impi)*mpi4(impi)/(w*w*w*w) + 8D0*mpi3(impi)*mpi3(impi)/(w*w) - &
                  q2*q2 - 2D0*mpi4(impi) ) + mpi3(impi)*mpi3(impi)*0.5D0/(w*w) )
  END FUNCTION chp_N3LO_2PE_Vc_1loop_ci0

! [1]Eq. D.18
  FUNCTION chp_N3LO_2PE_Vc_2loop_DR(    q2,       wt2, A, impi) RESULT(res)
    REAL(8) , INTENT(IN) :: q2(1:nof_theta_int_points)
    REAL(8) , INTENT(IN) :: wt2(1:nof_theta_int_points)
    REAL(8) , INTENT(IN) :: A(1:nof_theta_int_points)
    INTEGER, INTENT(IN) :: impi
    REAL(8) :: res(1:nof_theta_int_points)
    
    res = 3D0*gA2*wt2*A*const(21)*const(2)*(1D0/256D0)*( &
              (mpi2(impi) + 2D0*q2)*(2D0*chp_mpi(impi)%val + wt2*A) + &
              4D0*gA2*chp_mpi(impi)%val*wt2)
  END FUNCTION chp_N3LO_2PE_Vc_2loop_DR

! [2]Eq. 2.19 and 2.20, the part for V_c(q)
! Due to the large derivative of the z-integrand w.r.t z close to z=1, the integral
! is done numerically up to 0.95, then the contribution from 0.95 to 1.0 as a function
! of q^2 / m_\pi^2 is calculated from a taylor expansion around z = 1.
! The relative error of the approximate expression for the integral from 0.95 to 1.0
! is smaller than 2.6e-10 for all values of q^2 / m_\pi^2
  FUNCTION chp_N3LO_2PE_Vc_2loop_SFR(    q2,                   impi) RESULT(res)
    REAL(8) , INTENT(IN) :: q2(1:nof_theta_int_points)
    INTEGER, INTENT(IN) :: impi
    REAL(8) :: res(1:nof_theta_int_points)
    INTEGER :: i, j
    INTEGER :: z_nr
    REAL(8) :: fact
    REAL(8) :: a1, s
    REAL(8) :: s2
    REAL(8) :: z, z2, zt
! Magic numbers, I got these from a Taylor expansion calculated in Mathematica
    INTEGER, PARAMETER :: appr_order = 6
    REAL(8), PARAMETER, DIMENSION(0:appr_order) :: c = (/ 3593.92733718761d0, 5431.46429291511d0, 3420.6675650193d0, 1149.1163208075d0, 217.17174525573d0, 21.89317093303d0, 0.9197592690185d0 /)
    REAL(8) :: L
    REAL(8) :: x, m_pi_inv
    REAL(8) :: d1, d2, d3, d4, d5
    
    if(chiral_2PE_2loop_int%val == 0) then
      res = 0
      return
    end if

    z_nr = chp_2PE_2loop_int_mesh_z_VC_SFR%info%amount
    L = chp_2PE_2loop_int_mesh_z_VC_SFR%info%xmin

    fact = fpi_inv**6 * gA2*gA2 * 3.0d0 / (pi2*2.0d0**15)
    a1 = 4.0D0 * mpi2(impi)
    m_pi_inv = 1 / chp_mpi(impi)%val
    d1 = (L-1)*32*(1+4*gA2)*mpi2(impi)*mpi2(impi);
    d2 = -8*(L-1)*(-29 + L + L*L + 4*(-11 + L + L*L)*gA2)*mpi2(impi)/3;
    d3 = 2*(L-1)*(193 + 172*gA2 + L*(L+1)*(-47 - 68*gA2 + 3*L*L*(1 + 4*gA2)))/15;
    d4 = 32*chp_mpi(impi)%val*(1 + 4*gA2)*mpi2(impi);
    d5 = 64*chp_mpi(impi)%val * (1 + gA2);

!$OMP parallel do schedule(static) default(shared) shared(q2, impi, res, fact, a1, chp_2PE_2loop_int_mesh_z_VC_SFR, z_nr, gA2, m_pi_inv, d1, d2, d3, d4, d5, L) private(i, j, s, s2, z, z2, zt, x)
    DO i = 1, nof_theta_int_points
! Calculate the part from 0.95 to 1.0
      x = q2(i) * m_pi_inv**2
      s = c(appr_order)
      DO j = appr_order - 1, 0, -1
        s = x * s + c(j)
      END DO
      s = s * m_pi_inv**2 * (1/(4+x))**(appr_order+1)
! Calculate the rest as a numerical integration from z_low to 0.95
      DO j = 1, z_nr
        z = chp_2PE_2loop_int_mesh_z_VC_SFR%pnt(j)%x
        z2 = z*z
        zt = sqrt(1 - z2)
        s = s + chp_2PE_2loop_int_mesh_z_VC_SFR%pnt(j)%w * (2-z2)*(z2-8)*(z2-2)*z*0.5d0*log((1+z)/(1-z)) / (a1 + z2 * q2(i))
      END DO
      s = s * q2(i)**3
      x = 0.5*sqrt(x)
      s2 = s + d1 + q2(i)*(d2 + q2(i)*d3) + (0.5*a1 + q2(i)) * (d4 + d5*q2(i)) * (atan(x) - atan(L*x)) / sqrt(q2(i))
      res(i) = fact * s2
    END DO
!$OMP end parallel do
  END FUNCTION chp_N3LO_2PE_Vc_2loop_SFR

! [1]Eq. D.5
  FUNCTION chp_N3LO_2PE_Wc_1loop_ci1 (    q2, w, L,             impi, imnuc) RESULT(res)
    REAL(8) , INTENT(IN) :: q2(1:nof_theta_int_points)
    REAL(8) , INTENT(IN) :: w(1:nof_theta_int_points)
    REAL(8) , INTENT(IN) :: L(1:nof_theta_int_points)
    INTEGER, INTENT(IN) :: impi
    INTEGER, INTENT(IN) :: imnuc
    REAL(8) :: res(1:nof_theta_int_points)

    res = -c4*q2*L*const(21)*mnuc_inv(imnuc)*(1D0/192D0) * ( &
              gA2*(8D0*mpi2(impi) + 5D0*q2) + w*w)
  END FUNCTION chp_N3LO_2PE_Wc_1loop_ci1

! [1]Eq. D.10
  FUNCTION chp_N3LO_2PE_Wc_1loop_ci0 (k2, q2, w, L,             impi, imnuc) RESULT(res)
    REAL(8) , INTENT(IN) :: k2(1:nof_theta_int_points)
    REAL(8) , INTENT(IN) :: q2(1:nof_theta_int_points)
    REAL(8) , INTENT(IN) :: w(1:nof_theta_int_points)
    REAL(8) , INTENT(IN) :: L(1:nof_theta_int_points)
    INTEGER, INTENT(IN) :: impi
    INTEGER, INTENT(IN) :: imnuc
    REAL(8) :: res(1:nof_theta_int_points)
    
    res = -const(21)*mnuc_inv(imnuc)*mnuc_inv(imnuc)*(1D0/768D0)*( &
              L*(8D0*gA2*(1.5D0*q2*q2 + 3D0*mpi2(impi)*q2 + 3D0*mpi4(impi) - &
                      6D0*mpi3(impi)*mpi3(impi)/(w*w) - k2*(8D0*mpi2(impi) + 5D0*q2)) + &
                  4D0*gA2*gA2*(k2*(20D0*mpi2(impi) + 7D0*q2 - 16D0*mpi4(impi)/(w*w)) + &
                      16D0*mpi4(impi)*mpi4(impi)/(w*w*w*w) + &
                      12D0*mpi3(impi)*mpi3(impi)/(w*w) - 4D0*mpi4(impi)*q2/(w*w) - &
                      5D0*q2*q2 - 6D0*mpi2(impi)*q2 - 6D0*mpi4(impi)) - 4D0*k2*w*w) + &
              16D0*gA2*gA2*mpi3(impi)*mpi3(impi)/(w*w))
  END FUNCTION chp_N3LO_2PE_Wc_1loop_ci0

! [1]Eq. D.20
  FUNCTION chp_N3LO_2PE_Wc_2loop_DR_a(    q2, w, L, wt2,    impi) RESULT(res)
    REAL(8) , INTENT(IN) :: q2(1:nof_theta_int_points)
    REAL(8) , INTENT(IN) :: w(1:nof_theta_int_points)
    REAL(8) , INTENT(IN) :: L(1:nof_theta_int_points)
    REAL(8) , INTENT(IN) :: wt2(1:nof_theta_int_points)
    INTEGER, INTENT(IN) :: impi
    REAL(8) :: res(1:nof_theta_int_points)
    
! const(21) = 1 / (pi^2 * f_pi^4)
    res = L*const(21)*fpi_inv*fpi_inv*pi_inv*pi_inv*(1D0/18432D0)*( &
              192D0*pi2*fpi2*w*w*d3*(2D0*gA2*wt2 - 0.6D0*(gA2 - 1D0)*w*w) + &
              (6D0*gA2*wt2 - (gA2 - 1D0)*w*w)*(384D0*pi2*fpi2*(wt2*d1_plus_d2 + &
                      4D0*mpi2(impi)*d5) + L*(4D0*mpi2(impi)*(1D0+2D0*gA2) + &
                          q2*(1D0 + 5D0*gA2)) - (q2*(1D0/3D0)*(5D0 + 13D0*gA2) + &
                              8D0*mpi2(impi)*(1D0 + 2D0*gA2))))
  END FUNCTION chp_N3LO_2PE_Wc_2loop_DR_a

! Helper function needed by chp_calculate_2PE_2loop_int_*_data functions
  FUNCTION chp_N3LO_2PE_2loop_int_helper(x) RESULT(res)
    REAL(8), INTENT(IN) :: x
    REAL(8) :: x2, x4, x2_inv, A
    REAL(8) :: res

    IF(x < 0.05) THEN
! Use Taylor expansion around 0 with terms up to x^6 for small x
      x2 = x*x
      x4 = x2*x2

      res = (((-8.0D0/315.0D0)*x4*x2 + (2.0D0/35.0D0)*x4) - 0.2*x2) - 4.0D0/3.0D0
    ELSE
      x2_inv = 1.0D0 / (x*x)
      A = sqrt(1 + x2_inv)

      res = x2_inv - (1+x2_inv)*A * log(x * (1 + A))
    END IF
  END FUNCTION chp_N3LO_2PE_2loop_int_helper

! [1]Eq. D.21 and D.22
! Calculates all q-independent stuff in the z-integral (z = 2*m_\pi / \mu), i.e. it need only be done once.
  SUBROUTINE chp_calculate_2PE_2loop_int_WC_DR_data
    INTEGER :: z_nr, x_nr, i, j
    REAL(8) :: fact
    REAL(8) :: C1
    REAL(8) :: z, z2, zt
    REAL(8) :: res
    REAL(8) :: D1
    REAL(8) :: w, x, y, y2, y_z

    IF(chp_2PE_2loop_int_WC_DR_data_set) return

    z_nr = chp_2PE_2loop_int_mesh_z%info%amount
    x_nr = chp_2PE_2loop_int_mesh_x%info%amount
    fact = 1.0D0 / (2048.0D0 * pi2 * pi2)
    C1 = 2.0D0 * gA2 * gA2 / 3.0D0

!$OMP parallel do schedule(static) default(none) shared(z_nr, x_nr, fact, C1, chp_2PE_2loop_int_mesh_z, gA2, chp_2PE_2loop_int_mesh_x, chp_2PE_2loop_int_WC_DR_x_fact) private(i, j, res, z, z2, zt, D1, w, x, y, y2, y_z)
    DO i = 1, z_nr
      res = 0
      z = chp_2PE_2loop_int_mesh_z%pnt(i)%x
      z2 = z * z
      zt = sqrt(1 - z2)
      D1 = gA2 * (2 - z2)

      DO j = 1, x_nr
        w = chp_2PE_2loop_int_mesh_x%pnt(j)%w
        x = chp_2PE_2loop_int_mesh_x%pnt(j)%x
        y = x * zt
        y2 = y * y
        y_z = y / z

        res = res + fact * w * zt * z * ((gA2 - 1) * y2 - D1) * (-y2 + 2*y * sqrt(z2 + y2) * log(y_z + sqrt(1 + y_z*y_z)) + C1 * (2 - y2 - z2) * (chp_N3LO_2PE_2loop_int_helper(y_z) + 5.0D0/6.0D0))
      END DO

      chp_2PE_2loop_int_WC_DR_x_fact(i) = res
    END DO
!$OMP end parallel do

    chp_2PE_2loop_int_WC_DR_data_set = .TRUE.
  END SUBROUTINE chp_calculate_2PE_2loop_int_WC_DR_data

! [2]Eq. 2.19 and 2.20, the part for W_C(q)
! Calculates all q-independent stuff in the z-integral (z = 2*m_\pi / \mu), i.e. it need only be done once.
  SUBROUTINE chp_calculate_2PE_2loop_int_WC_SFR_data
    INTEGER :: z_nr, x_nr, i, j
    REAL(8) :: c1, c2, cz1, cz2
    REAL(8) :: z, z2, zt
    REAL(8) :: res
    REAL(8) :: w, x

    IF(chp_2PE_2loop_int_WC_SFR_data_set) return

    z_nr = chp_2PE_2loop_int_mesh_z%info%amount
    x_nr = chp_2PE_2loop_int_mesh_x%info%amount
    c1 = -gA2**3
    c2 = gA2*gA2*(2*gA2-1)

!$OMP parallel do schedule(static) default(none) shared(z_nr, x_nr, c1, c2, chp_2PE_2loop_int_mesh_z, gA2, chp_2PE_2loop_int_mesh_x, chp_2PE_2loop_int_WC_SFR_x_fact) private(i, j, res, z, z2, zt, cz1, cz2, w, x)
    DO i = 1, z_nr
      res = 0
      z = chp_2PE_2loop_int_mesh_z%pnt(i)%x
      z2 = z * z
      zt = sqrt(1 - z2)
      cz1 = c1 * (2-z2)**2
      cz2 = c2 * (2-z2)*(1-z2)

      DO j = 1, x_nr
        w = chp_2PE_2loop_int_mesh_x%pnt(j)%w
        x = chp_2PE_2loop_int_mesh_x%pnt(j)%x

        res = res + w * chp_N3LO_2PE_2loop_int_helper(x*zt/z) * (cz1 + cz2*x*x)
      END DO

      chp_2PE_2loop_int_WC_SFR_x_fact(i) = res
    END DO
!$OMP end parallel do

    chp_2PE_2loop_int_WC_SFR_data_set = .TRUE.
  END SUBROUTINE chp_calculate_2PE_2loop_int_WC_SFR_data

! [1]Eq. D.21 and D.22
  FUNCTION chp_N3LO_2PE_Wc_2loop_DR_b (    q2,                   impi) RESULT(res)
    REAL(8) , INTENT(IN) :: q2(1:nof_theta_int_points)
    INTEGER, INTENT(IN) :: impi
    REAL(8) :: res(1:nof_theta_int_points)
    INTEGER :: i, j
    REAL(8) :: fact
    REAL(8) :: a1
    REAL(8) :: s
    REAL(8) :: tmp
    INTEGER :: z_nr
    
    if(chiral_2PE_2loop_int%val == 0) then
      res = 0
      return
    end if

    call chp_calculate_2PE_2loop_int_WC_DR_data

    z_nr = chp_2PE_2loop_int_mesh_z%info%amount

    fact = fpi_inv**6
    a1 = 4.0D0 * mpi2(impi)

!$OMP parallel do schedule(static) default(none) shared(q2, impi, res, fact, a1, z_nr, chp_2PE_2loop_int_mesh_z, chp_2PE_2loop_int_WC_DR_x_fact) private(i, j, s, tmp)
    DO i = 1, nof_theta_int_points
      s = 0
      tmp = q2(i)**3
      DO j = 1, z_nr
        s = s + tmp * (chp_2PE_2loop_int_mesh_z%pnt(j)%w * chp_2PE_2loop_int_WC_DR_x_fact(j) / (a1 + chp_2PE_2loop_int_mesh_z%pnt(j)%x**2 * q2(i)))
      END DO
      res(i) = fact * s
    END DO
!$OMP end parallel do
  END FUNCTION chp_N3LO_2PE_Wc_2loop_DR_b

! [2]Eq. 2.19 and 2.20, the part for W_c(q), LEC independent part
  FUNCTION chp_N3LO_2PE_Wc_2loop_SFR_a(    q2,                   impi) RESULT(res)
    REAL(8) , INTENT(IN) :: q2(1:nof_theta_int_points)
    INTEGER, INTENT(IN) :: impi
    REAL(8) :: res(1:nof_theta_int_points)
    INTEGER :: i, j
    INTEGER :: z_nr
    REAL(8) :: fact
    REAL(8) :: a1, a2, a3, a4, b1, b2, b3
    REAL(8) :: d1
    REAL(8) :: z, z2, z4, zt
    REAL(8) :: s
    
    if(chiral_2PE_2loop_int%val == 0) then
      res = 0
      return
    end if

    call chp_calculate_2PE_2loop_int_WC_SFR_data

    z_nr = chp_2PE_2loop_int_mesh_z%info%amount

    fact = fpi_inv**6 / (pi2*pi2*3.0d0*2.0d0**10)
    d1 = 4.0D0 * mpi2(impi)
    a1 = -(1+2*gA2)**2/3.0d0
    a2 = 1 + 6*gA2 + 8*gA2*gA2
    a3 = -gA2*(8+15*gA2)
    a4 = (-4 + 29*gA2 + 122*gA2*gA2 + 3*gA2**3) / 15.0d0
    b1 = -(171 + 2*gA2*(1+gA2)*(327+49*gA2)) / 450.0d0
    b2 = (-73 + 1748*gA2 + 2549*gA2*gA2 + 726*gA2**3) / 450.0d0
    b3 = -(-64 + 389*gA2 + 1782*gA2*gA2 + 1093*gA2**3) / 450.0d0

!$OMP parallel do schedule(static) default(none) shared(q2, impi, res, fact, a1, a2, a3, a4, b1, b2, b3, d1, chp_2PE_2loop_int_mesh_z, chp_2PE_2loop_int_WC_SFR_x_fact, z_nr) private(i, j, s, z, z2, z4, zt)
    DO i = 1, nof_theta_int_points
      s = 0
      DO j = 1, z_nr
        z = chp_2PE_2loop_int_mesh_z%pnt(j)%x
        z2 = z**2
        z4 = z2**2
        zt = sqrt(1-z2)
        s = s + chp_2PE_2loop_int_mesh_z%pnt(j)%w * z / (d1 + z2 * q2(i)) * &
                ((a1*z4*z2 + a2*z4 + a3*z2 + a4)*log((1+zt)/z) + &
                (chp_2PE_2loop_int_WC_SFR_x_fact(j) + b1*z4 + b2*z2 + b3)*zt)
      END DO
      res(i) = fact * q2(i)**3 * s
    END DO
!$OMP end parallel do
  END FUNCTION chp_N3LO_2PE_Wc_2loop_SFR_a

! [2]Eq. 2.19 and 2.20, the part for W_c(q), LEC dependent part
  FUNCTION chp_N3LO_2PE_Wc_2loop_SFR_b(    q2,                   impi) RESULT(res)
    REAL(8) , INTENT(IN) :: q2(1:nof_theta_int_points)
    INTEGER, INTENT(IN) :: impi
    REAL(8) :: res(1:nof_theta_int_points)
    INTEGER :: i, j
    INTEGER :: z_nr
    REAL(8) :: fact1, fact2, fact3
    REAL(8) :: z, z2, zt
    REAL(8) :: s1, s2, s3
    REAL(8) :: d1, tmp
    
    if(chiral_2PE_2loop_int%val == 0) then
      res = 0
      return
    end if

    call chp_calculate_2PE_2loop_int_WC_SFR_data

    z_nr = chp_2PE_2loop_int_mesh_z%info%amount

    d1 = 4.0D0 * mpi2(impi)
    fact1 = -d1_plus_d2 * fpi_inv**4 / (pi2 * 96 )
    fact2 = -d3         * fpi_inv**4 / (pi2 * 480)
    fact3 =  d5         * fpi_inv**4 / (pi2 * 48 )

!$OMP parallel do schedule(static) default(none) shared(q2, impi, res, fact1, fact2, fact3, d1, chp_2PE_2loop_int_mesh_z, z_nr, gA2) private(i, j, s1, s2, s3, z, z2, zt, tmp)
    DO i = 1, nof_theta_int_points
      s1 = 0
      s2 = 0
      s3 = 0
      DO j = 1, z_nr
        z = chp_2PE_2loop_int_mesh_z%pnt(j)%x
        z2 = z**2
        zt = sqrt(1-z2)
        tmp = chp_2PE_2loop_int_mesh_z%pnt(j)%w * z * zt / (d1 + z2 * q2(i))
        s1 = s1 + tmp * (2-z2)*(  (gA2-1)*(1-z2) - 3*gA2*(2-z2))
        s2 = s2 + tmp * (1-z2)*(3*(gA2-1)*(1-z2) - 5*gA2*(2-z2))
        s3 = s3 + tmp * z2    *(  (gA2-1)*(1-z2) - 3*gA2*(2-z2))
      END DO
      res(i) = q2(i)**3 * (fact1 * s1 + fact2 * s2 + fact3 * s3)
    END DO
!$OMP end parallel do
  END FUNCTION chp_N3LO_2PE_Wc_2loop_SFR_b

! [1]Eq. D.7
  FUNCTION chp_N3LO_2PE_VLS_1loop_ci1(        w, L,             impi, imnuc) RESULT(res)
    REAL(8) , INTENT(IN) :: w(1:nof_theta_int_points)
    REAL(8) , INTENT(IN) :: L(1:nof_theta_int_points)
    INTEGER, INTENT(IN) :: impi
    INTEGER, INTENT(IN) :: imnuc
    REAL(8) :: res(1:nof_theta_int_points)
    
    res = (c2*gA2*const(21)*mnuc_inv(imnuc)*(1D0/8D0))*w*w*L
  END FUNCTION chp_N3LO_2PE_VLS_1loop_ci1

! [1]Eq. D.13
  FUNCTION chp_N3LO_2PE_VLS_1loop_ci0(    q2, w, L,             impi, imnuc) RESULT(res)
    REAL(8) , INTENT(IN) :: q2(1:nof_theta_int_points)
    REAL(8) , INTENT(IN) :: w(1:nof_theta_int_points)
    REAL(8) , INTENT(IN) :: L(1:nof_theta_int_points)
    INTEGER, INTENT(IN) :: impi
    INTEGER, INTENT(IN) :: imnuc
    REAL(8) :: res(1:nof_theta_int_points)
    
    res = (gA2*gA2*const(21)*0.25D0*mnuc_inv(imnuc)*mnuc_inv(imnuc))*L*( &
              (11D0/32D0)*q2 + mpi4(impi)/(w*w))
  END FUNCTION chp_N3LO_2PE_VLS_1loop_ci0

! [1]Eq. D.8
  FUNCTION chp_N3LO_2PE_WLS_1loop_ci1(    q2, w, L,             impi, imnuc) RESULT(res)
    REAL(8) , INTENT(IN) :: q2(1:nof_theta_int_points)
    REAL(8) , INTENT(IN) :: w(1:nof_theta_int_points)
    REAL(8) , INTENT(IN) :: L(1:nof_theta_int_points)
    INTEGER, INTENT(IN) :: impi
    INTEGER, INTENT(IN) :: imnuc
    REAL(8) :: res(1:nof_theta_int_points)
    
    res = (-c4*const(21)*mnuc_inv(imnuc)*(1D0/48D0))*L*( &
              gA2*(8D0*mpi2(impi) + 5D0*q2) + w*w)
  END FUNCTION chp_N3LO_2PE_WLS_1loop_ci1

! [1]Eq. D.14
  FUNCTION chp_N3LO_2PE_WLS_1loop_ci0(    q2, w, L,             impi, imnuc) RESULT(res)
    REAL(8) , INTENT(IN) :: q2(1:nof_theta_int_points)
    REAL(8) , INTENT(IN) :: w(1:nof_theta_int_points)
    REAL(8) , INTENT(IN) :: L(1:nof_theta_int_points)
    INTEGER, INTENT(IN) :: impi
    INTEGER, INTENT(IN) :: imnuc
    REAL(8) :: res(1:nof_theta_int_points)
    
    res = (const(21)*mnuc_inv(imnuc)*mnuc_inv(imnuc)*(1D0/256D0))*L*( &
              16D0*gA2*(mpi2(impi) + 0.375D0*q2) + (4D0/3D0)*gA2*gA2*( &
                  4D0*mpi4(impi)/(w*w) - (11D0/4D0)*q2 - 9D0*mpi2(impi)) - w*w)
  END FUNCTION chp_N3LO_2PE_WLS_1loop_ci0

! [1]Eq. D.11
  FUNCTION chp_N3LO_2PE_VT_1loop_ci0 (k2, q2, w, L,             impi, imnuc) RESULT(res)
    REAL(8) , INTENT(IN) :: k2(1:nof_theta_int_points)
    REAL(8) , INTENT(IN) :: q2(1:nof_theta_int_points)
    REAL(8) , INTENT(IN) :: w(1:nof_theta_int_points)
    REAL(8) , INTENT(IN) :: L(1:nof_theta_int_points)
    INTEGER, INTENT(IN) :: impi
    INTEGER, INTENT(IN) :: imnuc
    REAL(8) :: res(1:nof_theta_int_points)
    
    res = (gA2*gA2*const(21)*mnuc_inv(imnuc)*mnuc_inv(imnuc)*(1D0/32D0))*L*( &
              k2 + 0.625D0*q2 + mpi4(impi)/(w*w))
  END FUNCTION chp_N3LO_2PE_VT_1loop_ci0

! [1]Eq. D.24
  FUNCTION chp_N3LO_2PE_VT_2loop_DR_a   (        w, L,             impi) RESULT(res)
    REAL(8) , INTENT(IN) :: w(1:nof_theta_int_points)
    REAL(8) , INTENT(IN) :: L(1:nof_theta_int_points)
    INTEGER, INTENT(IN) :: impi
    REAL(8) :: res(1:nof_theta_int_points)
    
    res = (-gA2*const(21)*(1D0/32D0)*d14_minus_d15)*w*w*L
  END FUNCTION chp_N3LO_2PE_VT_2loop_DR_a   

! [2]Eq. 2.19 and 2.20, the part for V_T(q), LEC dependent part
  FUNCTION chp_N3LO_2PE_VT_2loop_SFR_a (    q2,                   impi) RESULT(res)
    REAL(8) , INTENT(IN) :: q2(1:nof_theta_int_points)
    INTEGER, INTENT(IN) :: impi
    REAL(8) :: res(1:nof_theta_int_points)
    INTEGER :: i, j
    INTEGER :: z_nr
    REAL(8) :: fact, s
    REAL(8) :: a1
    REAL(8) :: z, z2, zt
    
    if(chiral_2PE_2loop_int%val == 0) then
      res = 0
      return
    end if

    z_nr = chp_2PE_2loop_int_mesh_z%info%amount

    fact = -d14_minus_d15 * fpi_inv**4 * gA2 / (pi2*32.0d0)
    a1 = 4.0D0 * mpi2(impi)

!$OMP parallel do schedule(static) default(none) shared(q2, impi, res, fact, a1, z_nr, chp_2PE_2loop_int_mesh_z) private(i, j, s, z, z2, zt)
    DO i = 1, nof_theta_int_points
      s = 0
      DO j = 1, z_nr
        z = chp_2PE_2loop_int_mesh_z%pnt(j)%x
        z2 = z*z
        zt = sqrt(1 - z2)
        s = s + chp_2PE_2loop_int_mesh_z%pnt(j)%w * z*zt*(1-z2) / (a1 + z2 * q2(i))
      END DO
      res(i) = fact * q2(i)**2 * s
    END DO
!$OMP end parallel do
  END FUNCTION chp_N3LO_2PE_VT_2loop_SFR_a

! [1]Eq. D.25 (LEC independent part), identical with [2]Eq. 2.19 and 2.20, part for V_TS (LEC independent part).
! Only the integration limits differ in [1] (DR) and [2] (SFR)
! Calculates all q-independent stuff in the z-integral (z = 2*m_\pi / \mu) for the LEC independent part, i.e. it need only be done once.
  SUBROUTINE chp_calculate_2PE_2loop_int_VTS_data
    INTEGER :: z_nr, x_nr, i, j
    REAL(8) :: fact, z_fact
    REAL(8) :: res, z, zt
    REAL(8) :: w, x

    IF(chp_2PE_2loop_int_VTS_data_set) return

    z_nr = chp_2PE_2loop_int_mesh_z%info%amount
    x_nr = chp_2PE_2loop_int_mesh_x%info%amount
    fact = -1.0D0 / (1024.0D0 * pi2 * pi2)

!$OMP parallel do schedule(static) default(none) shared(z_nr, x_nr, fact, chp_2PE_2loop_int_mesh_z, chp_2PE_2loop_int_mesh_x, chp_2PE_2loop_int_VTS_x_fact) private(i, j, res, z, zt, w, x, z_fact)
    DO i = 1, z_nr
      res = 0
      z = chp_2PE_2loop_int_mesh_z%pnt(i)%x
      zt = sqrt(1 - z*z)
      z_fact = fact * chp_2PE_2loop_int_mesh_z%pnt(i)%w * z * zt * (1-z*z)

      DO j = 1, x_nr
        w = chp_2PE_2loop_int_mesh_x%pnt(j)%w
        x = chp_2PE_2loop_int_mesh_x%pnt(j)%x

        res = res + z_fact * w * (1 - x*x) * (chp_N3LO_2PE_2loop_int_helper(x*zt/z) - 1.0D0/6.0D0)
      END DO

      chp_2PE_2loop_int_VTS_x_fact(i) = res
    END DO
!$OMP end parallel do

    chp_2PE_2loop_int_VTS_data_set = .TRUE.
  END SUBROUTINE chp_calculate_2PE_2loop_int_VTS_data

! [1]Eq. D.25 (LEC independent part), identical with [2]Eq. 2.19 and 2.20, part for V_T(q) (LEC independent part).
  FUNCTION chp_N3LO_2PE_VT_2loop_b   (    q2,                   impi) RESULT(res)
    REAL(8) , INTENT(IN) :: q2(1:nof_theta_int_points)
    INTEGER, INTENT(IN) :: impi
    REAL(8) :: res(1:nof_theta_int_points)
    INTEGER :: i, j
    REAL(8) :: fact
    REAL(8) :: a1
    INTEGER :: z_nr
    
    if(chiral_2PE_2loop_int%val == 0) then
      res = 0
      return
    end if

    call chp_calculate_2PE_2loop_int_VTS_data

    z_nr = chp_2PE_2loop_int_mesh_z%info%amount
    fact = (chp_gA%val * fpi_inv)**6
    a1 = 4.0D0 * mpi2(impi)

!$OMP parallel do schedule(static) default(none) shared(z_nr, q2, impi, res, fact, a1, chp_2PE_2loop_int_mesh_z, chp_2PE_2loop_int_VTS_x_fact) private(i, j)
    DO i = 1, nof_theta_int_points
      res(i) = 0
      DO j = 1, z_nr
        res(i) = res(i) + chp_2PE_2loop_int_VTS_x_fact(j) / (a1 + chp_2PE_2loop_int_mesh_z%pnt(j)%x**2 * q2(i))
      END DO
      res(i) = res(i) * fact * q2(i)**2
    END DO
!$OMP end parallel do
  END FUNCTION chp_N3LO_2PE_VT_2loop_b   

! [1]Eq. D.2
  FUNCTION chp_N3LO_2PE_WT_1loop_ci2 (        w, L,             impi) RESULT(res)
    REAL(8) , INTENT(IN) :: w(1:nof_theta_int_points)
    REAL(8) , INTENT(IN) :: L(1:nof_theta_int_points)
    INTEGER, INTENT(IN) :: impi
    REAL(8) :: res(1:nof_theta_int_points)
    
    res = (c4*c4*const(21)*(1D0/96D0))*w*w*L
  END FUNCTION chp_N3LO_2PE_WT_1loop_ci2

! [1]Eq. D.6
  FUNCTION chp_N3LO_2PE_WT_1loop_ci1 (q2,     w, L,             impi, imnuc) RESULT(res)
    REAL(8) , INTENT(IN) :: q2(1:nof_theta_int_points)
    REAL(8) , INTENT(IN) :: w(1:nof_theta_int_points)
    REAL(8) , INTENT(IN) :: L(1:nof_theta_int_points)
    INTEGER, INTENT(IN) :: impi
    INTEGER, INTENT(IN) :: imnuc
    REAL(8) :: res(1:nof_theta_int_points)
    
    res = (-c4*const(21)*mnuc_inv(imnuc)*(1D0/192D0))*L*( &
              gA2*(16D0*mpi2(impi) + 7D0*q2) - w*w)
  END FUNCTION chp_N3LO_2PE_WT_1loop_ci1

! [1]Eq. D.12
  FUNCTION chp_N3LO_2PE_WT_1loop_ci0 (    q2, w, L,             impi, imnuc) RESULT(res)
    REAL(8) , INTENT(IN) :: q2(1:nof_theta_int_points)
    REAL(8) , INTENT(IN) :: w(1:nof_theta_int_points)
    REAL(8) , INTENT(IN) :: L(1:nof_theta_int_points)
    INTEGER, INTENT(IN) :: impi
    INTEGER, INTENT(IN) :: imnuc
    REAL(8) :: res(1:nof_theta_int_points)
    
    res = (mnuc_inv(imnuc)*mnuc_inv(imnuc)*const(21)*(1D0/1536D0))*L*( &
              4D0*gA2*gA2*(7D0*mpi2(impi) + (17D0/4D0)*q2 + 4D0*mpi4(impi)/(w*w)) - &
              32D0*gA2*(mpi2(impi) + (7D0/16D0)*q2) + w*w)
  END FUNCTION chp_N3LO_2PE_WT_1loop_ci0

! [1]Eq. D.27
  FUNCTION chp_N3LO_2PE_WT_2loop_DR(        w,             A, impi) RESULT(res)
    REAL(8) , INTENT(IN) :: w(1:nof_theta_int_points)
    REAL(8) , INTENT(IN) :: A(1:nof_theta_int_points)
    INTEGER, INTENT(IN) :: impi
    REAL(8) :: res(1:nof_theta_int_points)
    
    res = (gA2*gA2*const(21)*fpi_inv*fpi_inv*(1D0/2048D0))*w*w*A*( &
              w*w*A + 2D0*chp_mpi(impi)%val*(1D0 + 2D0*gA2))
  END FUNCTION chp_N3LO_2PE_WT_2loop_DR

! [2]Eq. 2.19 and 2.20, par for W_T(q)
  FUNCTION chp_N3LO_2PE_WT_2loop_SFR(    q2,                   impi) RESULT(res)
    REAL(8) , INTENT(IN) :: q2(1:nof_theta_int_points)
    INTEGER, INTENT(IN) :: impi
    REAL(8) :: res(1:nof_theta_int_points)
    INTEGER :: i, j
    INTEGER :: z_nr
    REAL(8) :: fact
    REAL(8) :: a1
    REAL(8) :: s
    REAL(8) :: z, z2, zt
    
    if(chiral_2PE_2loop_int%val == 0) then
      res = 0
      return
    end if

    z_nr = chp_2PE_2loop_int_mesh_z%info%amount

    fact = -fpi_inv**6 * gA2 * gA2 / (pi2*4096.0d0)
    a1 = 4.0D0 * mpi2(impi)

!$OMP parallel do schedule(static) default(none) shared(q2, impi, z_nr, gA2, res, fact, a1, chp_2PE_2loop_int_mesh_z) private(i, j, s, z, z2, zt)
    DO i = 1, nof_theta_int_points
      s = 0
      DO j = 1, z_nr
        z = chp_2PE_2loop_int_mesh_z%pnt(j)%x
        z2 = z*z
        zt = sqrt(1 - z2)
        s = s + chp_2PE_2loop_int_mesh_z%pnt(j)%w * (1-z2) * ((z2-1)*z*0.5d0*log((1+z)/(1-z)) + (1+2*gA2)*z2) / (a1 + z2 * q2(i))
      END DO
      res(i) = fact * q2(i)**2 * s
    END DO
!$OMP end parallel do
  END FUNCTION chp_N3LO_2PE_WT_2loop_SFR

! [1]Eq. D.11
  FUNCTION chp_N3LO_2PE_Vs_1loop_ci0 (k2, q2, w, L,             impi, imnuc) RESULT(res)
    REAL(8) , INTENT(IN) :: k2(1:nof_theta_int_points)
    REAL(8) , INTENT(IN) :: q2(1:nof_theta_int_points)
    REAL(8) , INTENT(IN) :: w(1:nof_theta_int_points)
    REAL(8) , INTENT(IN) :: L(1:nof_theta_int_points)
    INTEGER, INTENT(IN) :: impi
    INTEGER, INTENT(IN) :: imnuc
    REAL(8) :: res(1:nof_theta_int_points)
    
    res = -q2 * chp_N3LO_2PE_VT_1loop_ci0(k2, q2, w, L, impi, imnuc)
  END FUNCTION chp_N3LO_2PE_Vs_1loop_ci0

! [1]Eq. D.24
  FUNCTION chp_N3LO_2PE_Vs_2loop_DR_a   (    q2, w, L,             impi) RESULT(res)
    REAL(8) , INTENT(IN) :: q2(1:nof_theta_int_points)
    REAL(8) , INTENT(IN) :: w(1:nof_theta_int_points)
    REAL(8) , INTENT(IN) :: L(1:nof_theta_int_points)
    INTEGER, INTENT(IN) :: impi
    REAL(8) :: res(1:nof_theta_int_points)
    
    res = -q2 * chp_N3LO_2PE_VT_2loop_DR_a(w, L, impi)
  END FUNCTION chp_N3LO_2PE_Vs_2loop_DR_a

! [2]Eq. 2.19 and 2.20, the part for V_S(q), LEC dependent part
  FUNCTION chp_N3LO_2PE_Vs_2loop_SFR_a   (    q2,            impi) RESULT(res)
    REAL(8) , INTENT(IN) :: q2(1:nof_theta_int_points)
    INTEGER, INTENT(IN) :: impi
    REAL(8) :: res(1:nof_theta_int_points)
    
    res = -q2 * chp_N3LO_2PE_VT_2loop_SFR_a(q2, impi)
  END FUNCTION chp_N3LO_2PE_Vs_2loop_SFR_a

! [1]Eq. D.25 (LEC independent part), identical with [2]Eq. 2.19 and 2.20, part for V_S(q) (LEC independent part).
  FUNCTION chp_N3LO_2PE_Vs_2loop_b   (    q2,                   impi) RESULT(res)
    REAL(8) , INTENT(IN) :: q2(1:nof_theta_int_points)
    INTEGER, INTENT(IN) :: impi
    REAL(8) :: res(1:nof_theta_int_points)
    
    res = -q2 * chp_N3LO_2PE_VT_2loop_b(q2, impi)
  END FUNCTION chp_N3LO_2PE_Vs_2loop_b   

! [1]Eq. D.2
  FUNCTION chp_N3LO_2PE_Ws_1loop_ci2 (    q2, w, L,             impi) RESULT(res)
    REAL(8) , INTENT(IN) :: q2(1:nof_theta_int_points)
    REAL(8) , INTENT(IN) :: w(1:nof_theta_int_points)
    REAL(8) , INTENT(IN) :: L(1:nof_theta_int_points)
    INTEGER, INTENT(IN) :: impi
    REAL(8) :: res(1:nof_theta_int_points)
    
    res = -q2 * chp_N3LO_2PE_WT_1loop_ci2(w, L, impi)
  END FUNCTION chp_N3LO_2PE_Ws_1loop_ci2

! [1]Eq. D.6
  FUNCTION chp_N3LO_2PE_Ws_1loop_ci1 (    q2, w, L,             impi, imnuc) RESULT(res)
    REAL(8) , INTENT(IN) :: q2(1:nof_theta_int_points)
    REAL(8) , INTENT(IN) :: w(1:nof_theta_int_points)
    REAL(8) , INTENT(IN) :: L(1:nof_theta_int_points)
    INTEGER, INTENT(IN) :: impi
    INTEGER, INTENT(IN) :: imnuc
    REAL(8) :: res(1:nof_theta_int_points)
    
    res = -q2 * chp_N3LO_2PE_WT_1loop_ci1(q2, w, L, impi, imnuc)
  END FUNCTION chp_N3LO_2PE_Ws_1loop_ci1

! [1]Eq. D.12
  FUNCTION chp_N3LO_2PE_Ws_1loop_ci0 (    q2, w, L,             impi, imnuc) RESULT(res)
    REAL(8) , INTENT(IN) :: q2(1:nof_theta_int_points)
    REAL(8) , INTENT(IN) :: w(1:nof_theta_int_points)
    REAL(8) , INTENT(IN) :: L(1:nof_theta_int_points)
    INTEGER, INTENT(IN) :: impi
    INTEGER, INTENT(IN) :: imnuc
    REAL(8) :: res(1:nof_theta_int_points)
    
    res = -q2 * chp_N3LO_2PE_WT_1loop_ci0(q2, w, L, impi, imnuc)
  END FUNCTION chp_N3LO_2PE_Ws_1loop_ci0

! [1]Eq. D.27
  FUNCTION chp_N3LO_2PE_Ws_2loop_DR  (    q2, w,             A, impi) RESULT(res)
    REAL(8) , INTENT(IN) :: q2(1:nof_theta_int_points)
    REAL(8) , INTENT(IN) :: w(1:nof_theta_int_points)
    REAL(8) , INTENT(IN) :: A(1:nof_theta_int_points)
    INTEGER, INTENT(IN) :: impi
    REAL(8) :: res(1:nof_theta_int_points)
    
    res = -q2 * chp_N3LO_2PE_WT_2loop_DR(w, A, impi)
  END FUNCTION chp_N3LO_2PE_Ws_2loop_DR

! [2]Eq. 2.19 and 2.20, par for W_S(q)
  FUNCTION chp_N3LO_2PE_Ws_2loop_SFR(    q2, impi) RESULT(res)
    REAL(8) , INTENT(IN) :: q2(1:nof_theta_int_points)
    INTEGER, INTENT(IN) :: impi
    REAL(8) :: res(1:nof_theta_int_points)
    
    res = -q2 * chp_N3LO_2PE_WT_2loop_SFR(q2, impi)
  END FUNCTION chp_N3LO_2PE_Ws_2loop_SFR

! [1]Eq. D.15
  FUNCTION chp_N3LO_2PE_VsigL_1loop_ci0(         L,             impi, imnuc) RESULT(res)
    REAL(8) , INTENT(IN) :: L(1:nof_theta_int_points)
    INTEGER, INTENT(IN) :: impi
    INTEGER, INTENT(IN) :: imnuc
    REAL(8) :: res(1:nof_theta_int_points)

    res = (gA2*gA2*const(21)*mnuc_inv(imnuc)*mnuc_inv(imnuc)*(1D0/32D0))*L
  END FUNCTION chp_N3LO_2PE_VsigL_1loop_ci0



  SUBROUTINE chp_LO_contact_terms(chn, pfinal, pinit, POT)
    
    IMPLICIT NONE

    TYPE(chp_chn_type), INTENT(IN)    :: chn
    REAL(8), INTENT(INOUT) :: POT(6)
    REAL(8), intent(in) :: pfinal, pinit
    REAL(8) :: p, pp

    p = pinit
    pp = pfinal

!1S0
    IF (CHN%J == 0 .AND. CHN%S == 0 .AND. .NOT. CHN%coup) THEN
       POT(1) = POT(1) + chp_CLO(1)%val * freg(pp,p,chp_CLO_n(1)%val)
    ELSE IF (CHN%J == 1 .AND. CHN%S == 1 .AND. CHN%coup) THEN
!3S1
       POT(4)  = POT(4) + chp_CLO(2)%val * freg(pp,p,chp_CLO_n(2)%val)
    END IF

  END SUBROUTINE chp_LO_contact_terms
  
  SUBROUTINE chp_CIB_LO_contact_terms(chn, pfinal, pinit, POT)
    
    IMPLICIT NONE

    TYPE(chp_chn_type), INTENT(IN)    :: chn
    REAL(8), INTENT(INOUT) :: POT(6)
    REAL(8), intent(in) :: pfinal, pinit
    REAL(8) :: p, pp
    
    p = pinit
    pp = pfinal

! CIB contacts
!1S0
    IF (CHN%J == 0 .AND. CHN%S == 0 .AND. .NOT. CHN%coup) THEN
       POT(1) = POT(1) + chp_CIB_CLO(CHN%tz,1)%val * freg(pp,p,chp_CLO_n(1)%val)
!3S1
    ELSE IF (CHN%J == 1 .AND. CHN%S == 1 .AND. CHN%coup) THEN
       POT(4)  = POT(4) + chp_CIB_CLO(CHN%tz,2)%val * freg(pp,p,chp_CLO_n(2)%val)
    END IF
    
  END SUBROUTINE chp_CIB_LO_contact_terms

  SUBROUTINE chp_NLO_contact_terms(chn, pfinal, pinit, POT)
    
    IMPLICIT NONE

    TYPE(chp_chn_type), INTENT(IN)    :: chn
    REAL(8), INTENT(IN)    :: pfinal, pinit
    REAL(8), INTENT(INOUT) :: POT(6)
    REAL(8) :: p, pp, ppp
    
    p  = pinit
    pp = pfinal
    ppp=p*pp

! NLO contacts
!1S0
    IF (CHN%J == 0 .AND. CHN%S == 0 .AND. .NOT. CHN%coup) THEN
       POT(1) = POT(1) + chp_CNLO(1)%val*(p*p+pp*pp) * freg(pp,p,chp_CNLO_n(1)%val)
    ELSE IF (CHN%J == 0 .AND. CHN%S == 1 .AND. CHN%coup) THEN
!3P0
       POT(3)  = POT(3) + chp_CNLO(2)%val*ppp * freg(pp,p,chp_CNLO_n(2)%val)
    ELSE IF (CHN%J == 1 .AND. CHN%S == 0 .AND. .NOT. CHN%coup) THEN
!1P1
       POT(1) = POT(1) + chp_CNLO(3)%val*ppp * freg(pp,p,chp_CNLO_n(3)%val)
    ELSE IF (CHN%J == 1 .AND. CHN%S == 1 .AND. .NOT. CHN%coup) THEN
!3P1
       POT(2) = POT(2) + chp_CNLO(4)%val*ppp * freg(pp,p,chp_CNLO_n(4)%val)
    ELSE IF (CHN%J == 1 .AND. CHN%S == 1 .AND. CHN%coup) THEN
!3S1
       POT(4)  = POT(4) + chp_CNLO(5)%val*(p*p+pp*pp) * freg(pp,p,chp_CNLO_n(5)%val)
!3S1-3D1
       POT(6)  = POT(6) + chp_CNLO(6)%val*pp*pp * freg(pp,p,chp_CNLO_n(6)%val)
!3D1-3S1
       POT(5)  = POT(5) + chp_CNLO(6)%val*p*p * freg(pp,p,chp_CNLO_n(6)%val)
    ELSE IF (CHN%J == 2 .AND. CHN%S == 1 .AND. CHN%coup) THEN
!3P2
       POT(4)  = POT(4) + chp_CNLO(7)%val*p*pp * freg(pp,p,chp_CNLO_n(7)%val)
    END IF

  END SUBROUTINE chp_NLO_contact_terms

  SUBROUTINE chp_N3LO_contact_terms(chn, pfinal, pinit, POT)
    
    IMPLICIT NONE

    TYPE(chp_chn_type), INTENT(IN)    :: chn
    REAL(8)           , INTENT(IN)    :: pfinal, pinit
    REAL(8)        , INTENT(INOUT) :: POT(6)
    REAL(8) :: p, pp, ppp, p2, pp2
    
    p  = pinit
    p2 = p * p
    pp = pfinal
    pp2 = pp * pp
    ppp=p*pp

! N3LO contacts
! Eq [1] E.1, the contacts are numbered in the order they appear
! in the equation [1] E.1.
! FIXME: borisc: p and pp seems to be switched compared with what [1] uses. I do it this
!                way to be consistent with chp_NLO_contact_terms
!                UPDATE: On closer inspection, it seems that the function is called with
!                        reversed arguments, so that pinit <-> pfinal, so everything is OK
    
!1S0
    IF (CHN%J == 0 .AND. CHN%S == 0 .AND. .NOT. CHN%coup) THEN
       POT(1) = POT(1) + chp_DN3LO(1)%val*(p2 * p2 + pp2 * pp2) * freg(pp,p,chp_DN3LO_n(1)%val) + chp_DN3LO(2)%val*(p2 * pp2) * freg(pp,p,chp_DN3LO_n(2)%val)
    END IF
    
!3P0
    IF (CHN%J == 0 .AND. CHN%S == 1 .AND. CHN%coup) THEN
       POT(3) = POT(3) + chp_DN3LO(3)%val*ppp*(p2 + pp2) * freg(pp,p,chp_DN3LO_n(3)%val)
    END IF
    
!1P1
    IF (CHN%J == 1 .AND. CHN%S == 0 .AND. .NOT. CHN%coup) THEN
       POT(1) = POT(1) + chp_DN3LO(4)%val*ppp*(p2 + pp2) * freg(pp,p,chp_DN3LO_n(4)%val)
    END IF
    
!3P1
    IF (CHN%J == 1 .AND. CHN%S == 1 .AND. .NOT. CHN%coup) THEN
       POT(2) = POT(2) + chp_DN3LO(5)%val*ppp*(p2 + pp2) * freg(pp,p,chp_DN3LO_n(5)%val)
    END IF
    
!3S1
    IF (CHN%J == 1 .AND. CHN%S == 1 .AND. CHN%coup) THEN
       POT(4) = POT(4) + chp_DN3LO(6)%val*(p2 * p2 + pp2 * pp2) * freg(pp,p,chp_DN3LO_n(6)%val) + chp_DN3LO(7)%val*(p2 * pp2) * freg(pp,p,chp_DN3LO_n(7)%val)
    END IF
    
!3D1
    IF (CHN%J == 1 .AND. CHN%S == 1 .AND. CHN%coup) THEN
       POT(3) = POT(3) + chp_DN3LO(8)%val*(p2 * pp2) * freg(pp,p,chp_DN3LO_n(8)%val)
    END IF
    
!3S1-3D1
    IF (CHN%J == 1 .AND. CHN%S == 1 .AND. CHN%coup) THEN
       POT(6) = POT(6) + chp_DN3LO(9)%val*(pp2 * pp2) * freg(pp,p,chp_DN3LO_n(9)%val) + chp_DN3LO(10)%val*(p2 * pp2) * freg(pp,p,chp_DN3LO_n(10)%val)
    END IF
    
!3D1-3S1
    IF (CHN%J == 1 .AND. CHN%S == 1 .AND. CHN%coup) THEN
       POT(5) = POT(5) + chp_DN3LO(9)%val*(p2 * p2) * freg(pp,p,chp_DN3LO_n(9)%val) + chp_DN3LO(10)%val*(p2 * pp2) * freg(pp,p,chp_DN3LO_n(10)%val)
    END IF
    
!1D2
    IF (CHN%J == 2 .AND. CHN%S == 0 .AND. .NOT. CHN%coup) THEN
       POT(1) = POT(1) + chp_DN3LO(11)%val*(p2 * pp2) * freg(pp,p,chp_DN3LO_n(11)%val)
    END IF
    
!3D2
    IF (CHN%J == 2 .AND. CHN%S == 1 .AND. .NOT. CHN%coup) THEN
       POT(2) = POT(2) + chp_DN3LO(12)%val*(p2 * pp2) * freg(pp,p,chp_DN3LO_n(12)%val)
    END IF
    
!3P2
    IF (CHN%J == 2 .AND. CHN%S == 1 .AND. CHN%coup) THEN
       POT(4) = POT(4) + chp_DN3LO(13)%val*ppp*(p2 + pp2) * freg(pp,p,chp_DN3LO_n(13)%val)
    END IF
    
!3P2-3F2
    IF (CHN%J == 2 .AND. CHN%S == 1 .AND. CHN%coup) THEN
       POT(6) = POT(6) + chp_DN3LO(14)%val*ppp*(pp2) * freg(pp,p,chp_DN3LO_n(14)%val)
    END IF
    
!3F2-3P2
    IF (CHN%J == 2 .AND. CHN%S == 1 .AND. CHN%coup) THEN
       POT(5) = POT(5) + chp_DN3LO(14)%val*ppp*(p2) * freg(pp,p,chp_DN3LO_n(14)%val)
    END IF
    
!3D3
    IF (CHN%J == 3 .AND. CHN%S == 1 .AND. CHN%coup) THEN
       IF (chp_DN3LO_n(15)%val < 0.0D0) THEN
          POT(4) = POT(4) + chp_DN3LO(15)%val*(p2 * pp2) * 0.5D0 * (freg(pp,p,2.0D0) + freg(pp,p,3.0D0))
       ELSE
          POT(4) = POT(4) + chp_DN3LO(15)%val*(p2 * pp2) * freg(pp,p,chp_DN3LO_n(15)%val)
       END IF
    END IF

  END SUBROUTINE chp_N3LO_contact_terms


! cm energy defined on p.28
  FUNCTION cm_energy(p,tz) RESULT(res)

    REAL(8), INTENT(IN)  :: p
    INTEGER, INTENT(IN) :: tz
    REAL(8) :: res

    res = DSQRT(mnuc2(tz) + p*p)

  END FUNCTION cm_energy

! regulator, eq 4.63
! pfinal : final momentum
! pinit  : initial momentum
! n      : cutoff order
! LAMBDA is accessed from the chp constant chp_lambda
  FUNCTION freg(pfinal, pinit, n) RESULT(res)
    
    REAL(8) , INTENT(IN) :: pfinal, pinit
    REAL(8) , INTENT(IN) :: n
    REAL(8) :: res,exponent,lambda

    lambda = chp_lambda%val
    
    exponent = (pfinal/lambda)**(2.0D0*n) + &
               (pinit/lambda)**(2.0D0*n)
    
    res = dexp(-exponent)
    
  END FUNCTION freg

  SUBROUTINE initialize_chiral_potential
    
    CALL set_contact_conversion_matrices
    CALL setup_twobody_pwd(nof_theta_int_points, maximum_angular_momentum)
    RAirrel          = rirrel
    chp_mnuc(-1)     = chp_real_type(rirrel,'m_prot' ,.FALSE.)
    chp_mnuc(0)      = chp_real_type(rirrel,'m_nucl' ,.FALSE.)
    chp_mnuc(+1)     = chp_real_type(rirrel,'m_neut' ,.FALSE.)
    chp_mpi(-1)      = chp_real_type(rirrel,'m_pi-'  ,.FALSE.)
    chp_mpi(0)       = chp_real_type(rirrel,'m_pi0'  ,.FALSE.)
    chp_mpi(+1)      = chp_real_type(rirrel,'m_pi+'  ,.FALSE.)
    chp_mpi(+2)      = chp_real_type(rirrel,'m_pi '  ,.FALSE.)
    chp_ren          = chp_real_type(rirrel, 'REN'   , .FALSE.)
    chp_lambda       = chp_real_type(rirrel, 'LAMBDA', .FALSE.)
    chp_regcut_1PE   = chp_real_type(rirrel, 'rcut 1', .FALSE.)
    chp_regcut_2PE   = chp_real_type(rirrel, 'rcut 2', .FALSE.)
    chp_gA           = chp_real_type_RA(RAirrel, 'gA'    , .FALSE.)
    chp_fpi          = chp_real_type_RA(RAirrel, 'fpi'   , .FALSE.)
    chp_fine_structure = chp_real_type(rirrel, 'alpha'   , .FALSE.)

    chp_chiral_order = chp_int_type(iirrel, 'nu'     , .FALSE.)
    chp_chiral_mode  = chp_int_type(iirrel, 'nu_mode', .FALSE.)

    chiral_Ct = chp_int_type(iirrel, 'NN LO LEC', .FALSE.)
    chiral_Ct_CIB = chp_int_type(iirrel, 'NN LO CIB', .FALSE.)
    chiral_C = chp_int_type(iirrel, 'NN NLO LEC', .FALSE.)
    chiral_D = chp_int_type(iirrel, 'NN N3LO LEC', .FALSE.)
    chiral_1PE = chp_int_type(iirrel, '1PE', .FALSE.)
    chiral_1PE_CIB = chp_int_type(iirrel, '1PE CIB', .FALSE.)
    chiral_1PE_gamma = chp_int_type(iirrel, '1PE gamma', .FALSE.)
    chiral_1PE_relcorr = chp_int_type(iirrel, '1PE relcorr', .FALSE.)
    chiral_2PE_1loop_0 = chp_int_type(iirrel, '2PE 1lp 0', .FALSE.)
    chiral_2PE_1loop_d = chp_int_type(iirrel, '2PE 1lp d', .FALSE.)
    chiral_2PE_1loop_r = chp_int_type(iirrel, '2PE 1lp r', .FALSE.)
    chiral_2PE_1loop_r_mode = chp_int_type(iirrel, '2PE 1lp r M', .FALSE.)
    chiral_2PE_1loop_dd = chp_int_type(iirrel, '2PE 1lp dd', .FALSE.)
    chiral_2PE_1loop_dr = chp_int_type(iirrel, '2PE 1lp dr', .FALSE.)
    chiral_2PE_1loop_rr = chp_int_type(iirrel, '2PE 1lp rr', .FALSE.)
    chiral_2PE_2loop = chp_int_type(iirrel, '2PE 2lp', .FALSE.)
    chiral_2PE_2loop_int = chp_int_type(iirrel, '2PE 2lp int', .FALSE.)
    chiral_2PE_CSB_correct_mass = chp_int_type(iirrel, '2PE CSB mass', .FALSE.)
    chiral_minimal_relativity = chp_int_type(iirrel, 'Minimal Rel.', .FALSE.)
    chiral_kamada_glockle_transform = chp_int_type(0, 'Kam-Glo', .FALSE.)

    chp_2PE_2loop_int_VTS_data_set = .FALSE.
    chp_2PE_2loop_int_WC_DR_data_set  = .FALSE.
    chp_2PE_2loop_int_WC_SFR_data_set  = .FALSE.
    
! N2LO LECs c1 c3 c4
    chp_c1           = chp_real_type_RA(RAirrel, 'c1'   , .FALSE.)
    chp_c3           = chp_real_type_RA(RAirrel, 'c3'   , .FALSE.)
    chp_c4           = chp_real_type_RA(RAirrel, 'c4'   , .FALSE.)

! N3LO LECs c2, d1+d2, d3, d5, d14-d15
    chp_c2            = chp_real_type_RA(RAirrel, 'c2'     , .FALSE.)
    chp_d1_plus_d2    = chp_real_type_RA(RAirrel, 'd1+d2'  , .FALSE.)
    chp_d3            = chp_real_type_RA(RAirrel, 'd3'     , .FALSE.)
    chp_d5            = chp_real_type_RA(RAirrel, 'd5'     , .FALSE.)
    chp_d14_minus_d15 = chp_real_type_RA(RAirrel, 'd14-d15', .FALSE.)

    chp_e14           = chp_real_type_RA(RAirrel, 'e14'    , .FALSE.)
    chp_e15           = chp_real_type_RA(RAirrel, 'e15'    , .FALSE.)
    chp_e16           = chp_real_type_RA(RAirrel, 'e16'    , .FALSE.)
    chp_e17           = chp_real_type_RA(RAirrel, 'e17'    , .FALSE.)
    chp_e18           = chp_real_type_RA(RAirrel, 'e18'    , .FALSE.)

! LO contacts
    chp_CLO(1)       = chp_real_type_RA(RAirrel, 'x', .FALSE.)
    chp_CLO(2)       = chp_real_type_RA(RAirrel, 'x', .FALSE.)
    chp_CLO_n(1)     = chp_real_type(rirrel, 'x', .FALSE.)
    chp_CLO_n(2)     = chp_real_type(rirrel, 'x', .FALSE.)
    
! NLO contacts
    chp_CIB_CLO(-1,1)       = chp_real_type_RA(RAirrel, 'x', .FALSE.)
    chp_CIB_CLO(-1,2)       = chp_real_type_RA(RAirrel, 'x', .FALSE.)
    chp_CIB_CLO( 0,1)       = chp_real_type_RA(RAirrel, 'x', .FALSE.)
    chp_CIB_CLO( 0,2)       = chp_real_type_RA(RAirrel, 'x', .FALSE.)
    chp_CIB_CLO(+1,1)       = chp_real_type_RA(RAirrel, 'x', .FALSE.)
    chp_CIB_CLO(+1,2)       = chp_real_type_RA(RAirrel, 'x', .FALSE.)
    
    chp_CNLO(1)       = chp_real_type_RA(RAirrel, 'x', .FALSE.)
    chp_CNLO(2)       = chp_real_type_RA(RAirrel, 'x', .FALSE.)
    chp_CNLO(3)       = chp_real_type_RA(RAirrel, 'x', .FALSE.)
    chp_CNLO(4)       = chp_real_type_RA(RAirrel, 'x', .FALSE.)
    chp_CNLO(5)       = chp_real_type_RA(RAirrel, 'x', .FALSE.)
    chp_CNLO(6)       = chp_real_type_RA(RAirrel, 'x', .FALSE.)
    chp_CNLO(7)       = chp_real_type_RA(RAirrel, 'x', .FALSE.)
    chp_CNLO_n(1)     = chp_real_type(rirrel, 'x', .FALSE.)
    chp_CNLO_n(2)     = chp_real_type(rirrel, 'x', .FALSE.)
    chp_CNLO_n(3)     = chp_real_type(rirrel, 'x', .FALSE.)
    chp_CNLO_n(4)     = chp_real_type(rirrel, 'x', .FALSE.)
    chp_CNLO_n(5)     = chp_real_type(rirrel, 'x', .FALSE.)
    chp_CNLO_n(6)     = chp_real_type(rirrel, 'x', .FALSE.)
    chp_CNLO_n(7)     = chp_real_type(rirrel, 'x', .FALSE.)

! N3LO contacts
    chp_DN3LO( 1) = chp_real_type_RA(RAirrel, 'x', .FALSE.)
    chp_DN3LO( 2) = chp_real_type_RA(RAirrel, 'x', .FALSE.)
    chp_DN3LO( 3) = chp_real_type_RA(RAirrel, 'x', .FALSE.)
    chp_DN3LO( 4) = chp_real_type_RA(RAirrel, 'x', .FALSE.)
    chp_DN3LO( 5) = chp_real_type_RA(RAirrel, 'x', .FALSE.)
    chp_DN3LO( 6) = chp_real_type_RA(RAirrel, 'x', .FALSE.)
    chp_DN3LO( 7) = chp_real_type_RA(RAirrel, 'x', .FALSE.)
    chp_DN3LO( 8) = chp_real_type_RA(RAirrel, 'x', .FALSE.)
    chp_DN3LO( 9) = chp_real_type_RA(RAirrel, 'x', .FALSE.)
    chp_DN3LO(10) = chp_real_type_RA(RAirrel, 'x', .FALSE.)
    chp_DN3LO(11) = chp_real_type_RA(RAirrel, 'x', .FALSE.)
    chp_DN3LO(12) = chp_real_type_RA(RAirrel, 'x', .FALSE.)
    chp_DN3LO(13) = chp_real_type_RA(RAirrel, 'x', .FALSE.)
    chp_DN3LO(14) = chp_real_type_RA(RAirrel, 'x', .FALSE.)
    chp_DN3LO(15) = chp_real_type_RA(RAirrel, 'x', .FALSE.)
    chp_DN3LO_n( 1) = chp_real_type(rirrel, 'x', .FALSE.)
    chp_DN3LO_n( 2) = chp_real_type(rirrel, 'x', .FALSE.)
    chp_DN3LO_n( 3) = chp_real_type(rirrel, 'x', .FALSE.)
    chp_DN3LO_n( 4) = chp_real_type(rirrel, 'x', .FALSE.)
    chp_DN3LO_n( 5) = chp_real_type(rirrel, 'x', .FALSE.)
    chp_DN3LO_n( 6) = chp_real_type(rirrel, 'x', .FALSE.)
    chp_DN3LO_n( 7) = chp_real_type(rirrel, 'x', .FALSE.)
    chp_DN3LO_n( 8) = chp_real_type(rirrel, 'x', .FALSE.)
    chp_DN3LO_n( 9) = chp_real_type(rirrel, 'x', .FALSE.)
    chp_DN3LO_n(10) = chp_real_type(rirrel, 'x', .FALSE.)
    chp_DN3LO_n(11) = chp_real_type(rirrel, 'x', .FALSE.)
    chp_DN3LO_n(12) = chp_real_type(rirrel, 'x', .FALSE.)
    chp_DN3LO_n(13) = chp_real_type(rirrel, 'x', .FALSE.)
    chp_DN3LO_n(14) = chp_real_type(rirrel, 'x', .FALSE.)
    chp_DN3LO_n(15) = chp_real_type(rirrel, 'x', .FALSE.)



  END SUBROUTINE initialize_chiral_potential
  
  SUBROUTINE chp_set_units_and_derive_constants
    
    USE chp_aux
    
    implicit none
    INTEGER :: i, row
    REAL(8), ALLOCATABLE, DIMENSION(:) :: tmp_contacts
    REAL(8) :: reg_val
    REAL(8) :: z_low

    IF (.NOT. chp_chiral_order%set) THEN
       ERROR STOP 'Must set chiral_order'
    END IF
    IF (.NOT. chp_chiral_mode%set .AND. chp_chiral_order%val >= LO) THEN
       ERROR STOP 'Must set chiral_mode'
    END IF

    IF (chp_chiral_order%val >= LO) THEN
       IF (.NOT. chiral_1PE%set) call chp_set_chiral_1PE(1)
       IF (.NOT. chiral_Ct%set) call chp_set_chiral_Ct(1)
       IF (chp_chiral_mode%val == chiral_mode_Ep2005 .OR. chp_chiral_mode%val == chiral_mode_Ep2015) THEN
          IF (.NOT. chiral_kamada_glockle_transform%set) call chp_set_chiral_kamada_glockle_transform(1)
       ELSE
          IF (.NOT. chiral_minimal_relativity%set) call chp_set_chiral_minimal_relativity(1)
       END IF
    END IF
    IF (chp_chiral_order%val >= NLO) THEN
       IF (.NOT. chiral_2PE_1loop_0%set) call chp_set_chiral_2PE_1loop_0(1)
       IF (.NOT. chiral_C%set) call chp_set_chiral_C(1)
       IF (.NOT. chiral_1PE_CIB%set) call chp_set_chiral_1PE_CIB(1)
       IF (.NOT. chiral_Ct_CIB%set) call chp_set_chiral_Ct_CIB(1)
    END IF
    IF (chp_chiral_order%val >= NNLO) THEN
       IF (.NOT. chiral_2PE_1loop_d%set) call chp_set_chiral_2PE_1loop_d(1)
       IF (chp_chiral_mode%val == chiral_mode_EM2011) THEN
          IF (.NOT. chiral_2PE_1loop_r%set) call chp_set_chiral_2PE_1loop_r(1)
          IF (.NOT. chiral_2PE_1loop_r_mode%set) call chp_set_chiral_2PE_1loop_r_mode(chiral_2PE_1loop_r_mode_EM2011)
       END IF
    END IF
    IF (chp_chiral_order%val >= N3LO) THEN
       IF (.NOT. chiral_2PE_1loop_r%set) call chp_set_chiral_2PE_1loop_r(1)
       IF (.NOT. chiral_2PE_1loop_dd%set) call chp_set_chiral_2PE_1loop_dd(1)
       IF (.NOT. chiral_2PE_2loop%set) call chp_set_chiral_2PE_2loop(1)
       IF (chp_chiral_mode%val /= chiral_mode_EM2011) THEN
          IF (.NOT. chiral_2PE_2loop_int%set) call chp_set_chiral_2PE_2loop_int(1)
       END IF
       IF (.NOT. chiral_D%set) call chp_set_chiral_D(1)
       IF (chp_chiral_mode%val == chiral_mode_Ep2005) THEN
          IF (.NOT. chiral_2PE_1loop_r_mode%set) call chp_set_chiral_2PE_1loop_r_mode(chiral_2PE_1loop_r_mode_Ep2005)
       ELSE IF (chp_chiral_mode%val == chiral_mode_Ep2015) THEN
          IF (.NOT. chiral_2PE_1loop_r_mode%set) call chp_set_chiral_2PE_1loop_r_mode(chiral_2PE_1loop_r_mode_Ep2015)
       ELSE
          IF (.NOT. chiral_2PE_1loop_r_mode%set) call chp_set_chiral_2PE_1loop_r_mode(chiral_2PE_1loop_r_mode_EM2015)
       END IF
       IF (chp_chiral_mode%val == chiral_mode_EM2011) THEN
          IF (.NOT. chiral_2PE_1loop_dr%set) call chp_set_chiral_2PE_1loop_dr(1)
          IF (.NOT. chiral_2PE_1loop_rr%set) call chp_set_chiral_2PE_1loop_rr(1)
          IF (.NOT. chiral_1PE_gamma%set) call chp_set_chiral_1PE_gamma(1)
          IF (.NOT. chiral_2PE_CSB_correct_mass%set) call chp_set_chiral_2PE_CSB_correct_mass(1)
       END IF
       IF (chp_chiral_mode%val == chiral_mode_Ep2005 .OR. chp_chiral_mode%val == chiral_mode_Ep2015) THEN
          IF (.NOT. chiral_1PE_relcorr%set) call chp_set_chiral_1PE_relcorr(1)
       END IF
    END IF

! set the regulator cutoffs depending on the chiral order
    IF (.NOT. chp_regcut_1PE%set) THEN
       IF (chp_chiral_order%val == N3LO .AND. chp_chiral_mode%val == chiral_mode_EM2011) THEN
          CALL chp_set_1PE_reg_par(4.0D0)
       ELSE
          CALL chp_set_1PE_reg_par(3.0D0)
       END IF
    END IF
    IF (.NOT. chp_regcut_2PE%set) THEN
       IF (chp_chiral_order%val == N3LO .AND. chp_chiral_mode%val == chiral_mode_EM2011) THEN
          CALL chp_set_2PE_reg_par(2.0D0)
       ELSE
          CALL chp_set_2PE_reg_par(3.0D0)
       END IF
    END IF

! Set regulator parameter for LO contacts if not set
    DO i = 1, 2
       IF(chp_CLO_n(i)%set) CYCLE
       CALL chp_set_LO_contact_reg_par(i, 3.0D0)
    END DO

! Set regulator parameter for NLO contacts if not set
    DO i = 1, 7
       IF(chp_CNLO_n(i)%set) CYCLE
       IF(chp_chiral_order%val == N3LO .AND. chp_chiral_mode%val == chiral_mode_EM2011) THEN
          reg_val = 2.0D0
       ELSE
          reg_val = 3.0D0
       END IF
       CALL chp_set_NLO_contact_reg_par(i, reg_val)
    END DO

! Set regulator parameter for N3LO contacts if not set
    DO i = 1, 15
       IF(chp_DN3LO_n(i)%set) CYCLE
       IF (chp_chiral_order%val == N3LO .AND. chp_chiral_mode%val == chiral_mode_EM2011) THEN
          SELECT CASE(i)
          CASE(1)
             reg_val = 2.0D0
          CASE(2)
             reg_val = 2.0D0
          CASE(3)
             reg_val = 3.0D0
          CASE(4)
             reg_val = 2.0D0
          CASE(5)
             reg_val = 4.0D0
          CASE(6)
             reg_val = 2.0D0
          CASE(7)
             reg_val = 2.0D0
          CASE(8)
             reg_val = 2.0D0
          CASE(9)
             reg_val = 2.0D0
          CASE(10)
             reg_val = 2.0D0
          CASE(11)
             reg_val = 4.0D0
          CASE(12)
             reg_val = 2.0D0
          CASE(13)
             reg_val = 2.0D0
          CASE(14)
             reg_val = 4.0D0
          CASE(15)
             reg_val = -1.0D0
          END SELECT
       ELSE
          reg_val = 3.0D0
       END IF
       CALL chp_set_N3LO_contact_reg_par(i, reg_val)
    END DO



    c1 = chp_c1%val
    c3 = chp_c3%val
    c4 = chp_c4%val
    c2 = chp_c2%val
    d1_plus_d2 = chp_d1_plus_d2%val
    d3 = chp_d3%val
    d5 = chp_d5%val
    d14_minus_d15 = chp_d14_minus_d15%val

    mnuc2(:)   = chp_mnuc(:)%val*chp_mnuc(:)%val
    mnuc_inv(:)= 1D0 / chp_mnuc(:)%val
    mpi2(:)    = chp_mpi(:)%val*chp_mpi(:)%val
    mpi3(:)    = chp_mpi(:)%val*mpi2(:)
    mpi4(:)    = mpi2*mpi2
    mpi5(:)    = mpi4*chp_mpi(:)%val
    twompi(:)  = 2.0D0*chp_mpi(:)%val
    fourmpi2(:)= 4.0D0*mpi2
    
    IF (chp_ren%name == 'SF') THEN
       
       sfr  = chp_ren%val
       sfr2 = sfr*sfr

       DO i=-1,2
          IF ( (sfr-twompi(i)) <  0.0D0 ) sfr_heavyside(i) = 0.0D0
          IF ( (sfr-twompi(i)) >= 0.0D0 ) sfr_heavyside(i) = 1.0D0
       END DO

       z_low = 2 * chp_mpi(2)%val / sfr
       IF (z_low > 1) z_low = 1

    ELSE
       z_low = 0
    END IF
! z = 2*m_\pi / \mu, using the notation of [1] and [2] for the 2PE 2-loop diagrams
! In DR, the integration is done from z = 0 and in SFR from z = 2 * m_\pi / \Lambda_{SFR}
    CALL chp_setup_2PE_2loop_int_data(nof_2PE_2loop_int_x_points, nof_2PE_2loop_int_z_points, z_low)
    
    gA2        = chp_gA%val*chp_gA%val
    gA4        = gA2*gA2
    
    fpi2       = chp_fpi%val*chp_fpi%val
    fpi4       = fpi2*fpi2
    fpi_inv    = 1D0 / chp_fpi%val

    iso(0) = -3.0D0
    iso(1) = +1.0D0
    
    const = 0.0D0
!1: 1/(2pi)^3
    const(1) = 1.0D0/(twopi**3)
!2: gA^2/(4*fpi^2)
    const(2) = gA2/(4.0D0*fpi2)
!3: 384pi^2*fpi^4
    const(3) = 384.0D0*pi2*fpi4
!4: 5gA^4-4gA^2-1
    const(4) = 5.0D0*gA4-4.0D0*gA2-1.0D0
!5: 23gA^4-10gA^2-1
    const(5) = 23.0D0*gA4-10.0D0*gA2 - 1.0D0
!6: 48gA^4
    const(6) = 48.0D0*gA4
!7: 3gA^4
    const(7) = 3.0D0*gA4
!8: 64pi^2fpi^4
    const(8) = 64.0D0*pi2*fpi4
!9: 3gA^2/16pifpi^4
    const(9) = 3.0D0*gA2/(16.0D0*chp_pi%val*fpi4)
!10: ga^2/16
    const(10) = ga2/16.0D0
!11: 2.0D0*(2c1-c3)
    const(11) = 2.0D0*(2.0D0*c1-c3)
!12 : const(7)/256pifpi^4
    const(12) = const(7)/(256.0D0*chp_pi%val*fpi4)
!13: gA^2/(128pifpi^4)
    const(13) = gA2/(128.0D0*chp_pi%val*fpi4)
!14: gA^4/(128pifpi^4)
    const(14) = gA4/(128.0D0*chp_pi%val*fpi4)
!15: 3gA^4/(512pifpi^4)
    const(15) = 3.0D0*gA4/(512.0D0*chp_pi%val*fpi4)
!16: gA2/(32pifpi^4)
    const(16) = gA2/(32.0D0*chp_pi%val*fpi4)
!17: gA2/8
    const(17) = gA2/8.0D0
!18: gA4/(256pifpi^4)
    const(18) = gA4/(256.0D0*chp_pi%val*fpi4)
!19: 3gA4/(32pifpi^4)
    const(19) = 3.0D0*gA4/(32.0D0*chp_pi%val*fpi4)
!20: const(16)*(1-gA2)
    const(20) = const(16)*(1.0D0-gA2)
!21: 1 / (pi^2fpi^4)
    const(21) = 1D0 / (pi2 * fpi4)
! -----------------------------------------------------
    
! if the contacts were input in PW, transform to ST and print
    IF (chp_contact_format%name == 'PW') THEN
       
!WRITE(CHP_SCR,*)
!WRITE(CHP_SCR,"(A)") '   *************************'
!WRITE(CHP_SCR,"(A)") '   CONTACTS INPUT FORMAT: PW'
!WRITE(CHP_SCR,"(A)") '   IN THE ST FORMAT THEY ARE:'
!WRITE(CHP_SCR,"(A)") '   *************************'
!WRITE(CHP_SCR,*)
       
       IF (chiral_Ct_CIB%val == 1) THEN
! do the CIB LO contacts
          DO i=-1,1
             
             IF (ALLOCATED(tmp_contacts)) DEALLOCATE(tmp_contacts)
             ALLOCATE(tmp_contacts(1:2))
             tmp_contacts = 0.0D0
             
!simple matrix vector multiplication
             
             DO row=1, 2
                tmp_contacts(row) = SUM(LOPW2ST(row,:)*chp_CIB_CLO(i,:)%val)
             END DO

             SELECT CASE (i)
             CASE(-1)
!WRITE(CHP_SCR,"(A12,F30.16)") 'CSpp', tmp_contacts(1)
!WRITE(CHP_SCR,"(A12,F30.16)") 'CTpp', tmp_contacts(2)
             CASE( 0)
!WRITE(CHP_SCR,"(A12,F30.16)") 'CSpn', tmp_contacts(1)
!WRITE(CHP_SCR,"(A12,F30.16)") 'CTpn', tmp_contacts(2)
             CASE(+1)
!WRITE(CHP_SCR,"(A12,F30.16)") 'CSnn', tmp_contacts(1)
!WRITE(CHP_SCR,"(A12,F30.16)") 'CTnn', tmp_contacts(2)
             END SELECT
          END DO
       ELSE IF (chiral_Ct%val == 1) THEN
! LO CONTACTS
          IF (ALLOCATED(tmp_contacts)) DEALLOCATE(tmp_contacts)
          ALLOCATE(tmp_contacts(1:2))
          tmp_contacts = 0.0D0

!simple matrix vector multiplication
          
          DO row=1, 2
             tmp_contacts(row) = SUM(LOPW2ST(row,:)*chp_CLO(:)%val)
          END DO
          
!WRITE(CHP_SCR,"(A12,F30.16)") 'CS', tmp_contacts(1)
!WRITE(CHP_SCR,"(A12,F30.16)") 'CT', tmp_contacts(2)
       END IF

       IF (chiral_C%val == 1) THEN
! do the NLO contacts
          IF (ALLOCATED(tmp_contacts)) DEALLOCATE(tmp_contacts)
          ALLOCATE(tmp_contacts(1:7))
          tmp_contacts = 0.0D0
          
!simple matrix vector multiplication
          
          DO row=1,7
             tmp_contacts(row) = SUM(NLOPW2ST(row,:)*chp_CNLO(:)%val)
          END DO
          
!WRITE(CHP_SCR,"(A12,F30.16)") 'C1', tmp_contacts(1)
!WRITE(CHP_SCR,"(A12,F30.16)") 'C2', tmp_contacts(2)
!WRITE(CHP_SCR,"(A12,F30.16)") 'C3', tmp_contacts(3)
!WRITE(CHP_SCR,"(A12,F30.16)") 'C4', tmp_contacts(4)
!WRITE(CHP_SCR,"(A12,F30.16)") 'C5', tmp_contacts(5)
!WRITE(CHP_SCR,"(A12,F30.16)") 'C6', tmp_contacts(6)
!WRITE(CHP_SCR,"(A12,F30.16)") 'C7', tmp_contacts(7)
       END IF
       
       IF (chiral_D%val == 1) THEN
! do the N3LO contacts
          
          IF (ALLOCATED(tmp_contacts)) DEALLOCATE(tmp_contacts)
          ALLOCATE(tmp_contacts(1:15))
          tmp_contacts = 0.0D0
          
!simple matrix vector multiplication
          
          DO row=1,15
             tmp_contacts(row) = SUM(N3LOPW2ST(row,:)*chp_DN3LO(:)%val)
          END DO
          
!WRITE(CHP_SCR,"(A12,F30.16)") 'D1' , tmp_contacts( 1)
!WRITE(CHP_SCR,"(A12,F30.16)") 'D2' , tmp_contacts( 2)
!WRITE(CHP_SCR,"(A12,F30.16)") 'D3' , tmp_contacts( 3)
!WRITE(CHP_SCR,"(A12,F30.16)") 'D4' , tmp_contacts( 4)
!WRITE(CHP_SCR,"(A12,F30.16)") 'D5' , tmp_contacts( 5)
!WRITE(CHP_SCR,"(A12,F30.16)") 'D6' , tmp_contacts( 6)
!WRITE(CHP_SCR,"(A12,F30.16)") 'D7' , tmp_contacts( 7)
!WRITE(CHP_SCR,"(A12,F30.16)") 'D8' , tmp_contacts( 8)
!WRITE(CHP_SCR,"(A12,F30.16)") 'D9' , tmp_contacts( 9)
!WRITE(CHP_SCR,"(A12,F30.16)") 'D10', tmp_contacts(10)
!WRITE(CHP_SCR,"(A12,F30.16)") 'D11', tmp_contacts(11)
!WRITE(CHP_SCR,"(A12,F30.16)") 'D12', tmp_contacts(12)
!WRITE(CHP_SCR,"(A12,F30.16)") 'D13', tmp_contacts(13)
!WRITE(CHP_SCR,"(A12,F30.16)") 'D14', tmp_contacts(14)
!WRITE(CHP_SCR,"(A12,F30.16)") 'D15', tmp_contacts(15)
        END IF
       
    END IF
       
! if the contacts were input in ST, transform to PW replace, and print
    IF (chp_contact_format%name == 'ST') THEN
       
       WRITE(CHP_SCR,*)
       WRITE(CHP_SCR,"(A)") '   *************************'
       WRITE(CHP_SCR,"(A)") '   CONTACTS INPUT FORMAT: ST'
       WRITE(CHP_SCR,"(A)") '   IN THE PW FORMAT THEY ARE:'
       WRITE(CHP_SCR,"(A)") '   *************************'
       WRITE(CHP_SCR,*)

! LO
! FIXME: borisc: What about higher orders? They don't seem to
!                be converted to PW format.
       IF (chiral_Ct%val == 1) THEN
          
! matrix vector multiplication
          
          IF (ALLOCATED(tmp_contacts)) DEALLOCATE(tmp_contacts)
          ALLOCATE(tmp_contacts(1:2))
          tmp_contacts = 0.0D0

!simple matrix vector multiplication
          
          DO row=1, 2
             tmp_contacts(row) = SUM(LOST2PW(row,:)*chp_CLO(:)%val)
          END DO
          
          WRITE(CHP_SCR,"(A12,F30.16)") '1S0', tmp_contacts(1)
          WRITE(CHP_SCR,"(A12,F30.16)") '3S1', tmp_contacts(2)
          
! replace the ST matrix elements with st ones
          
          DO i=1,2
             chp_CLO(i)%val = tmp_contacts(i)
          END DO
          
! change the format identifier
          chp_contact_format%name = 'PW'
       END IF

    END IF
    
! convert the contact parameters to units of [MeV]
    
! LO
! the LO contacts are input in units of 10^4/GeV^2
! this = 10^-2/MeV^2
!
!IF (chp_chiral_order%val == LO) THEN
!
!   DO i=1, 2
!      chp_CLO(i)%val = chp_CLO(i)%val * 0.01D0
!   END DO
!
!END IF

! NLO
! the LO CIBcontacts are input in units of 10^4/GeV^2
! this = 10^-2/MeV^2
! the NLO contacts are input in units of 10^4/GeV^4
! this = 10^-8/MeV^4
! the N3LO contacts are input in units of 10^4/GeV^6
! this = 10^(4-3*6)/MeV^6 = 10^-14/MeV^6
    
!IF (chp_chiral_order%val > LO) THEN
!
!   DO i=1, 2
!      DO j=-1,1
!         chp_CIB_CLO(j,i)%val = chp_CIB_CLO(j,i)%val * 0.01D0
!      END DO
!   END DO
!
!   DO i=1,7
!      chp_CNLO(i)%val = chp_CNLO(i)%val * 1.E-08
!   END DO
!
!END IF
    
!IF (chp_chiral_order%val > NNLO) THEN
!   DO i=1,15
!      chp_DN3LO(i)%val = chp_DN3LO(i)%val * 1.D-14
!   END DO
!END IF


! -----------------------------------------------------
    
! DO SOME CONSISTENCY CHECKS

!IF (chp_gA%val == 1.29D0 .AND. chp_chiral_order%val == LO) THEN
!   WRITE(CHP_ERR,"(A)")
!   WRITE(CHP_ERR,"(A)") '   ****************************************'
!   WRITE(CHP_ERR,"(A)") '   AT LO, YOU DONT HAVE TO CORRECT FOR THE '
!   WRITE(CHP_ERR,"(A)") '   GOLDBERGER-TREIMAN DISCREPANCY. CHANGE'
!   WRITE(CHP_ERR,"(A)") '   gA FROM 1.29 TO 1.276'
!   WRITE(CHP_ERR,"(A)") '   ****************************************'
!   WRITE(CHP_ERR,"(A)")
!END IF
    IF (ALLOCATED(tmp_contacts)) DEALLOCATE(tmp_contacts)
    
  END SUBROUTINE chp_set_units_and_derive_constants

  SUBROUTINE chp_print_constants(title,unit)
    
    INTEGER :: unit,i,j
    CHARACTER(LEN=42) :: title

    WRITE(unit,"(A42)") '------------------------------------------'
    WRITE(unit,"(A42)") title
    WRITE(unit,"(A42)") '------------------------------------------'

    WRITE(unit,"(A)") 'CHIRAL POTENTIAL: NUMERICS'
    CALL print_pwd_numerics(unit)
    WRITE(unit,"(A)") 'CHIRAL POTENTIAL: CONSTANTS'
    CALL chp_print_rconst_if(unit, chp_mnuc(-1))
    CALL chp_print_rconst_if(unit, chp_mnuc(0))
    CALL chp_print_rconst_if(unit, chp_mnuc(+1))
    CALL chp_print_rconst_if(unit, chp_mpi(-1))
    CALL chp_print_rconst_if(unit, chp_mpi(0))
    CALL chp_print_rconst_if(unit, chp_mpi(+1))
    CALL chp_print_rconst_if(unit, chp_mpi(+2))
    WRITE(unit,"(A42)") '------------------------------------------'
    CALL chp_print_iconst(unit, chp_chiral_order)
    CALL chp_print_iconst(unit, chp_chiral_mode)
    CALL chp_print_rconst_if(unit, chp_ren)
    CALL chp_print_rconst_if(unit, chp_lambda)
    CALL chp_print_rconst_if(unit, chp_regcut_1PE)
    CALL chp_print_rconst_if(unit, chp_regcut_2PE)
    CALL chp_print_iconst(unit, chiral_Ct)
    CALL chp_print_iconst(unit, chiral_Ct_CIB)
    CALL chp_print_iconst(unit, chiral_C)
    CALL chp_print_iconst(unit, chiral_D)
    CALL chp_print_iconst(unit, chiral_1PE)
    CALL chp_print_iconst(unit, chiral_1PE_CIB)
    CALL chp_print_iconst(unit, chiral_1PE_gamma)
    CALL chp_print_iconst(unit, chiral_1PE_relcorr)
    CALL chp_print_iconst(unit, chiral_2PE_1loop_0)
    CALL chp_print_iconst(unit, chiral_2PE_1loop_d)
    CALL chp_print_iconst(unit, chiral_2PE_1loop_r)
    CALL chp_print_iconst(unit, chiral_2PE_1loop_r_mode)
    CALL chp_print_iconst(unit, chiral_2PE_1loop_dd)
    CALL chp_print_iconst(unit, chiral_2PE_1loop_dr)
    CALL chp_print_iconst(unit, chiral_2PE_1loop_rr)
    CALL chp_print_iconst(unit, chiral_2PE_2loop)
    CALL chp_print_iconst(unit, chiral_2PE_2loop_int)
    CALL chp_print_iconst(unit, chiral_2PE_CSB_correct_mass)
    CALL chp_print_iconst(unit, chiral_minimal_relativity)
    CALL chp_print_iconst(unit, chiral_kamada_glockle_transform)
    CALL chp_print_rconst_if(unit, chp_gA)
    CALL chp_print_rconst_if(unit, chp_fpi)
    CALL chp_print_rconst_if(unit, chp_fine_structure)
    WRITE(unit,"(A42)") '------------------------------------------'
    IF (chiral_Ct_CIB%val == 1) THEN
       DO i=1,2
          DO j=-1,1
             CALL chp_print_rconst_if(unit, chp_CIB_CLO(j,i))
          END DO
       END DO
    ELSEIF (chiral_Ct%val == 1) THEN
       CALL chp_print_rconst_if(unit, chp_CLO(1))
       CALL chp_print_rconst_if(unit, chp_CLO(2))
    END IF
    IF (chiral_C%val == 1) THEN
       DO i=1,7
          CALL chp_print_rconst_if(unit, chp_CNLO(i))
       END DO
    END IF
    IF (chiral_D%val == 1) THEN
       DO i=1,15
          CALL chp_print_rconst_if(unit, chp_DN3LO(i))
       END DO
    END IF
    IF (chp_chiral_order%val >= NNLO) THEN
       CALL chp_print_rconst_if(unit, chp_c1)
       CALL chp_print_rconst_if(unit, chp_c3)
       CALL chp_print_rconst_if(unit, chp_c4)
    END IF
    IF (chp_chiral_order%val >= N3LO) THEN
       CALL chp_print_rconst_if(unit, chp_c2)
       CALL chp_print_rconst_if(unit, chp_d1_plus_d2)
       CALL chp_print_rconst_if(unit, chp_d3)
       CALL chp_print_rconst_if(unit, chp_d5)
       CALL chp_print_rconst_if(unit, chp_d14_minus_d15)
    END IF

    
  END SUBROUTINE chp_print_constants
  
  SUBROUTINE chp_print_rconst(unit, chp_const)
    
    INTEGER, INTENT(IN) :: unit
    TYPE(chp_real_type) :: chp_const
    
    IF (chp_const%set)  THEN
       WRITE(unit,"(A12,F30.16)") chp_const%name, chp_const%val
    ELSE
       WRITE(unit,"(A12,A30)") chp_const%name, '----------------'
    END IF
    
  END SUBROUTINE chp_print_rconst

! MODIFIED, NEEDED TO PRINT RARealD-consts
  SUBROUTINE chp_print_rconst_RA(unit, chp_const)
    
    INTEGER, INTENT(IN) :: unit
    TYPE(chp_real_type_RA) :: chp_const
    
    IF (chp_const%set)  THEN
       WRITE(unit,"(A12,F30.16)") chp_const%name, chp_const%val
    ELSE
       WRITE(unit,"(A12,A30)") chp_const%name, '----------------'
    END IF
  END SUBROUTINE chp_print_rconst_RA
    
  SUBROUTINE chp_print_cconst(unit, chp_const)
    
    INTEGER, INTENT(IN) :: unit
    TYPE(chp_char2_type) :: chp_const
    
    IF (chp_const%set)  THEN
       WRITE(unit,"(A12,A30)") chp_const%name, chp_const%val
    ELSE
       WRITE(unit,"(A12,A30)") chp_const%name, '----------------'
    END IF
    
  
END SUBROUTINE chp_print_cconst

  SUBROUTINE chp_print_iconst(unit, chp_const)
    
    INTEGER, INTENT(IN) :: unit
    TYPE(chp_int_type) :: chp_const
    
    IF (chp_const%set)  THEN
       WRITE(unit,"(A12,I30)") chp_const%name, chp_const%val
    ELSE
       WRITE(unit,"(A12,A30)") chp_const%name, '----------------'
    END IF
    
  END SUBROUTINE chp_print_iconst

  SUBROUTINE chp_set_mass_nucleon(set_mnuc)
    
    REAL(8), INTENT(IN) :: set_mnuc(-1:1)
    
    chp_mnuc(-1)%val  = set_mnuc(-1) ! proton mass
    chp_mnuc(-1)%set  = .TRUE.
    chp_mnuc(0)%val   = set_mnuc(0) ! nucleon mass
    chp_mnuc(0)%set   = .TRUE.
    chp_mnuc(+1)%val  = set_mnuc(+1) ! neutron mass
    chp_mnuc(+1)%set  = .TRUE.
    
  END SUBROUTINE chp_set_mass_nucleon
  
  SUBROUTINE chp_get_mass_nucleon(get_mnuc)
    
    REAL(8), INTENT(INOUT) :: get_mnuc(-1:1)
    
    IF (.NOT. chp_mnuc(-1)%set) THEN
       WRITE(CHP_ERR,"(A)") 'mnuc(-1) NOT SET'
       STOP
    END IF
    
    get_mnuc(-1) = chp_mnuc(-1)%val ! proton mass

    IF (.NOT. chp_mnuc(0)%set) THEN
       WRITE(CHP_ERR,"(A)") 'mnuc(0) NOT SET'
       STOP
    END IF
    
    get_mnuc(0) = chp_mnuc(0)%val ! nucleon mass

    IF (.NOT. chp_mnuc(+1)%set) THEN
       WRITE(CHP_ERR,"(A)") 'mnuc(+1) NOT SET'
       STOP
    END IF
    
    get_mnuc(+1) = chp_mnuc(+1)%val ! neutron mass
        
  END SUBROUTINE chp_get_mass_nucleon

  SUBROUTINE chp_set_mass_pion(set_mpi)
    
    REAL(8), INTENT(IN) :: set_mpi(-1:1)
    
    chp_mpi(-1)%val  = set_mpi(-1) ! pi-
    chp_mpi(-1)%set  = .TRUE.
    chp_mpi(0)%val   = set_mpi(0)  ! pi0
    chp_mpi(0)%set   = .TRUE.
    chp_mpi(+1)%val  = set_mpi(+1) ! pi+
    chp_mpi(+1)%set  = .TRUE.
    
    chp_mpi(+2)%val = (chp_mpi(+1)%val + chp_mpi(-1)%val + chp_mpi(0)%val)/3.0D0
    chp_mpi(+2)%set = .TRUE.
    
  END SUBROUTINE chp_set_mass_pion

  SUBROUTINE chp_get_mass_pion(get_mpi)
    
    REAL(8), INTENT(INOUT) :: get_mpi(-1:2)
    
    IF (.NOT. chp_mpi(-1)%set) THEN
       WRITE(CHP_ERR,"(A)") 'chp_mpi(-1) NOT SET'
       STOP
    END IF
    
    get_mpi(-1) = chp_mpi(-1)%val! pi-

    IF (.NOT. chp_mpi(01)%set) THEN
       WRITE(CHP_ERR,"(A)") 'chp_mpi(01) NOT SET'
       STOP
    END IF
    
    get_mpi(0) = chp_mpi(0)%val! pi0

    IF (.NOT. chp_mpi(+1)%set) THEN
       WRITE(CHP_ERR,"(A)") 'chp_mpi(+1) NOT SET'
       STOP
    END IF
    
    get_mpi(+1) = chp_mpi(+1)%val! pi+

    IF (.NOT. chp_mpi(+2)%set) THEN
       WRITE(CHP_ERR,"(A)") 'chp_mpi(+2) NOT SET'
       STOP
    END IF
    
    get_mpi(+2) = chp_mpi(+2)%val! average pion mass

  END SUBROUTINE chp_get_mass_pion
  
  SUBROUTINE chp_set_chiral_order(set_chiral_order)
    
    INTEGER, INTENT(IN) :: set_chiral_order
    
    IF (set_chiral_order == 1 .OR. set_chiral_order>4 .OR. set_chiral_order <-1) THEN
       WRITE(CHP_ERR,"(A,I5)") 'error(set_chiral_order): illegal chiral order', set_chiral_order
       STOP
    END IF
    chp_chiral_order%val = set_chiral_order
    chp_chiral_order%set = .TRUE.

  END SUBROUTINE chp_set_chiral_order

  SUBROUTINE chp_get_chiral_order(get_chiral_order)
    
    INTEGER, INTENT(INOUT) :: get_chiral_order
    
    IF (.NOT. chp_chiral_order%set) THEN
       WRITE(CHP_ERR,"(A)") 'chp_chiral_order NOT SET'
       STOP
    END IF
    
    get_chiral_order = chp_chiral_order%val
    
  END SUBROUTINE chp_get_chiral_order
  
  SUBROUTINE chp_set_chiral_mode(set_chiral_mode)
    
    INTEGER, INTENT(IN) :: set_chiral_mode
    
    IF (set_chiral_mode < chiral_mode_EM2011 .OR. set_chiral_mode > chiral_mode_Ep2015) THEN
       WRITE(CHP_ERR,"(A,I5)") 'error(set_chiral_mode): illegal chiral mode', set_chiral_mode
       STOP
    END IF
    chp_chiral_mode%val = set_chiral_mode
    chp_chiral_mode%set = .TRUE.

  END SUBROUTINE chp_set_chiral_mode

  SUBROUTINE chp_get_chiral_mode(get_chiral_mode)
    
    INTEGER, INTENT(INOUT) :: get_chiral_mode
    
    IF (.NOT. chp_chiral_mode%set) THEN
       WRITE(CHP_ERR,"(A)") 'chp_chiral_mode NOT SET'
       STOP
    END IF
    
    get_chiral_mode = chp_chiral_mode%val
    
  END SUBROUTINE chp_get_chiral_mode
  
  SUBROUTINE chp_set_1PE_reg_par(val)
    REAL(8), INTENT(IN) :: val

    chp_regcut_1PE%set = .TRUE.
    chp_regcut_1PE%val = val
  END SUBROUTINE chp_set_1PE_reg_par

  SUBROUTINE chp_set_2PE_reg_par(val)
    REAL(8), INTENT(IN) :: val

    chp_regcut_2PE%set = .TRUE.
    chp_regcut_2PE%val = val
  END SUBROUTINE chp_set_2PE_reg_par

  SUBROUTINE chp_set_reg(set_reg, set_cutoff)
    
    CHARACTER(LEN=2), INTENT(IN) :: set_reg
    REAL(8), INTENT(IN) :: set_cutoff
    
    IF (set_reg /= 'SF' .AND. set_reg/= 'DR') THEN
       WRITE(CHP_ERR, "(A,A)") 'error(chp_set_reg): illegal regularization choice', set_reg
       STOP
    END IF
    
    chp_ren%val  = set_cutoff
    chp_ren%name = set_reg
    chp_ren%set  = .TRUE.
    
  END SUBROUTINE chp_set_reg
  
  SUBROUTINE chp_get_reg(get_reg_type,get_cutoff)
    
    CHARACTER(LEN=2), INTENT(INOUT) :: get_reg_type
    REAL(8), INTENT(INOUT) :: get_cutoff
    
    IF (.NOT. chp_ren%set) THEN
       WRITE(CHP_ERR,"(A)") 'chp_ren NOT SET'
       STOP
    END IF
    
    get_cutoff = chp_ren%val
    get_reg_type = chp_ren%name(1:2)
    
  END SUBROUTINE chp_get_reg
  
  SUBROUTINE chp_set_contact_format(set_cnt_fmt)
    
    CHARACTER(LEN=2), INTENT(IN) :: set_cnt_fmt
    
    IF (set_cnt_fmt /= 'PW' .AND. set_cnt_fmt/= 'ST') THEN
       WRITE(CHP_ERR, "(A,A)") 'error(chp_contact_format): illegal contact format', set_cnt_fmt
       STOP
    END IF
    
    chp_contact_format%name = set_cnt_fmt
    chp_contact_format%set  = .TRUE.
    
  END SUBROUTINE chp_set_contact_format

  SUBROUTINE chp_get_contact_format(get_cnt_fmt)
    
    CHARACTER(LEN=2), INTENT(INOUT) :: get_cnt_fmt
    
    IF (.NOT. chp_contact_format%set) THEN
       WRITE(CHP_ERR,"(A)") 'chp_contact_format NOT SET'
       STOP
    END IF
    
    get_cnt_fmt = chp_contact_format%name(1:2)
    
  END SUBROUTINE chp_get_contact_format
  
  SUBROUTINE chp_set_LO_contact_reg_par(contact_no, val)
    INTEGER, INTENT(IN) :: contact_no
    REAL(8), INTENT(IN) :: val

    chp_CLO_n(contact_no)%set = .TRUE.
    chp_CLO_n(contact_no)%val = val
  END SUBROUTINE chp_set_LO_contact_reg_par

  SUBROUTINE chp_set_LO_contact(contact_no, val)    
    INTEGER, INTENT(IN) :: contact_no
    REAL(8), INTENT(IN) :: val

    IF (.NOT. chp_contact_format%set) THEN
       WRITE(CHP_ERR,"(A)") 'error(chp_set_LO_contact): you must set the format before the contact value'
       STOP
    END IF
! LO
! the LO contacts are input in units of 10^4/GeV^2
! this = 10^-2/MeV^2
!

    chp_CLO(contact_no)%val  = val*0.01D0
    chp_CLO(contact_no)%set  = .TRUE.
   
    SELECT CASE(contact_no)
    CASE(1)
       IF (chp_contact_format%name == 'PW') chp_CLO(1)%name = '1S0'
       IF (chp_contact_format%name == 'ST') chp_CLO(1)%name = 'CS->1S0'
    CASE(2)
       IF (chp_contact_format%name == 'PW') chp_CLO(2)%name = '3S1'
       IF (chp_contact_format%name == 'ST') chp_CLO(2)%name = 'CT->3S1'
    CASE DEFAULT
       WRITE(CHP_ERR,"(A,I5)") 'error(chp_set_LO_contact): illegal contact nummber', contact_no
       STOP
    END SELECT

  END SUBROUTINE chp_set_LO_contact



  SUBROUTINE chp_get_LO_contact(contact_no, val, name)    
    INTEGER, INTENT(IN) :: contact_no
    REAL(8), INTENT(INOUT) :: val
    CHARACTER(LEN=12), INTENT(INOUT) :: name
    
! LO
! the LO contacts are input in units of 10^4/GeV^2
! this = 10^-2/MeV^2
!
    
    IF (.NOT. chp_CLO(contact_no)%set) THEN
       WRITE(CHP_ERR,"(A,I5,A)") 'chp_CLO',contact_no, ' NOT SET'
       STOP
    END IF
    
    val = chp_CLO(contact_no)%val * 100.0D0
    name = chp_CLO(contact_no)%name

  END SUBROUTINE chp_get_LO_contact
  
  SUBROUTINE chp_set_CIB_LO_contact(contact_no, tz, val)    
    INTEGER, INTENT(IN) :: contact_no, tz
    REAL(8), INTENT(IN) :: val

   IF (.NOT. chp_contact_format%set) THEN
       WRITE(CHP_ERR,"(A)") 'error(chp_set_LO_contact): you must set the format before the contact value'
       STOP
    END IF
    
    chp_CIB_CLO(tz,contact_no)%val  = val*0.01D0
    chp_CIB_CLO(tz,contact_no)%set  = .TRUE.
    
    IF (tz == -1) THEN
       SELECT CASE(contact_no)
       CASE(1)
          IF (chp_contact_format%name == 'PW') chp_CIB_CLO(tz,1)%name = '1S0pp'
          IF (chp_contact_format%name == 'ST') chp_CIB_CLO(tz,1)%name = 'CS->1S0pp'
       CASE(2)
          IF (chp_contact_format%name == 'PW') chp_CIB_CLO(tz,2)%name = '3S1pp'
          IF (chp_contact_format%name == 'ST') chp_CIB_CLO(tz,2)%name = 'CT->3S1pp'
       CASE DEFAULT
          WRITE(CHP_ERR,"(A,I5)") 'error(chp_set_LO_contact): illegal contact nummber', contact_no
          STOP
       END SELECT
    END IF

    IF (tz == 0) THEN
       SELECT CASE(contact_no)
       CASE(1)
          IF (chp_contact_format%name == 'PW') chp_CIB_CLO(tz,1)%name = '1S0pn'
          IF (chp_contact_format%name == 'ST') chp_CIB_CLO(tz,1)%name = 'CS->1S0pn'
       CASE(2)
          IF (chp_contact_format%name == 'PW') chp_CIB_CLO(tz,2)%name = '3S1pn'
          IF (chp_contact_format%name == 'ST') chp_CIB_CLO(tz,2)%name = 'CT->3S1pn'
       CASE DEFAULT
          WRITE(CHP_ERR,"(A,I5)") 'error(chp_set_LO_contact): illegal contact nummber', contact_no
          STOP
       END SELECT
    END IF

    IF (tz == +1) THEN
       SELECT CASE(contact_no)
       CASE(1)
          IF (chp_contact_format%name == 'PW') chp_CIB_CLO(tz,1)%name = '1S0nn'
          IF (chp_contact_format%name == 'ST') chp_CIB_CLO(tz,1)%name = 'CS->1S0nn'
       CASE(2)
          IF (chp_contact_format%name == 'PW') chp_CIB_CLO(tz,2)%name = '3S1nn'
          IF (chp_contact_format%name == 'ST') chp_CIB_CLO(tz,2)%name = 'CT->3S1nn'
       CASE DEFAULT
          WRITE(CHP_ERR,"(A,I5)") 'error(chp_set_LO_contact): illegal contact nummber', contact_no
          STOP
       END SELECT
    END IF
  END SUBROUTINE chp_set_CIB_LO_contact



  SUBROUTINE chp_get_CIB_LO_contact(contact_no, tz, val, name)    
    INTEGER, INTENT(IN) :: contact_no, tz
    REAL(8), INTENT(INOUT) :: val
    CHARACTER(LEN=12), INTENT(INOUT) :: name
    
    IF (.NOT. chp_CIB_CLO(tz,contact_no)%set) THEN
       WRITE(CHP_ERR,"(A,2I5,A)") 'chp_CIB_CLO ', tz, contact_no, ' NOT SET'
       STOP
    END IF

! the LO CIBcontacts are input in units of 10^4/GeV^2
! this = 10^-2/MeV^2
    
    val = chp_CIB_CLO(tz,contact_no)%val*100.0D0
    name = chp_CIB_CLO(tz,contact_no)%name
    
  END SUBROUTINE chp_get_CIB_LO_contact
  
  SUBROUTINE chp_set_NLO_contact_reg_par(contact_no, val)
    INTEGER, INTENT(IN) :: contact_no
    REAL(8), INTENT(IN) :: val

    chp_CNLO_n(contact_no)%set = .TRUE.
    chp_CNLO_n(contact_no)%val = val
  END SUBROUTINE chp_set_NLO_contact_reg_par

  SUBROUTINE chp_set_NLO_contact(contact_no, val)    
    INTEGER, INTENT(IN) :: contact_no
    REAL(8), INTENT(IN) :: val

   IF (.NOT. chp_contact_format%set) THEN
       WRITE(CHP_ERR,"(A)") 'error(chp_set_NLO_contact): you must set the format before the contact value'
       STOP
    END IF
    
! the NLO contacts are input in units of 10^4/GeV^4
! this = 10^-8/MeV^4
    chp_CNLO(contact_no)%val  = val*1.E-08
    chp_CNLO(contact_no)%set  = .TRUE.
    
    SELECT CASE(contact_no)
    CASE(1)
       IF (chp_contact_format%name == 'PW') chp_CNLO(1)%name = '1S0'
       IF (chp_contact_format%name == 'ST') chp_CNLO(1)%name = 'C1->1S0'
    CASE(2)
       IF (chp_contact_format%name == 'PW') chp_CNLO(2)%name = '3P0'
       IF (chp_contact_format%name == 'ST') chp_CNLO(2)%name = 'C2->3P0'
    CASE(3)
       IF (chp_contact_format%name == 'PW') chp_CNLO(3)%name = '1P1'
       IF (chp_contact_format%name == 'ST') chp_CNLO(3)%name = 'C3->1P1'
    CASE(4)
       IF (chp_contact_format%name == 'PW') chp_CNLO(4)%name = '3P1'
       IF (chp_contact_format%name == 'ST') chp_CNLO(4)%name = 'C4->3P1'
    CASE(5)
       IF (chp_contact_format%name == 'PW') chp_CNLO(5)%name = '3S1'
       IF (chp_contact_format%name == 'ST') chp_CNLO(5)%name = 'C5->3S1'
    CASE(6)
       IF (chp_contact_format%name == 'PW') chp_CNLO(6)%name = '3S1-3D1'
       IF (chp_contact_format%name == 'ST') chp_CNLO(6)%name = 'C6->3S1-3D1'
    CASE(7)
       IF (chp_contact_format%name == 'PW') chp_CNLO(7)%name = '3P2'
       IF (chp_contact_format%name == 'ST') chp_CNLO(7)%name = 'C7->3P2'
    CASE DEFAULT
       WRITE(CHP_ERR,"(A,I5)") 'error(chp_set_LO_contact): illegal contact nummber', contact_no
       STOP
    END SELECT
    
  END SUBROUTINE chp_set_NLO_contact



  SUBROUTINE chp_get_NLO_contact(contact_no, val, name)    
    
    INTEGER, INTENT(IN) :: contact_no
    REAL(8), INTENT(INOUT) :: val
    CHARACTER(LEN=12), INTENT(INOUT) :: name
    
    IF (.NOT. chp_CNLO(contact_no)%set) THEN
       WRITE(CHP_ERR,"(A,I5,A)") 'chp_CNLO ', contact_no, ' NOT SET'
       STOP
    END IF
! the NLO contacts are input in units of 10^4/GeV^4
! this = 10^-8/MeV^4
    val = chp_CNLO(contact_no)%val*1.0E+08
    name = chp_CNLO(contact_no)%name
    
  END SUBROUTINE chp_get_NLO_contact

  SUBROUTINE chp_set_N3LO_contact_reg_par(contact_no, val)
    INTEGER, INTENT(IN) :: contact_no
    REAL(8), INTENT(IN) :: val

    chp_DN3LO_n(contact_no)%set = .TRUE.
    chp_DN3LO_n(contact_no)%val = val
  END SUBROUTINE chp_set_N3LO_contact_reg_par

  SUBROUTINE chp_set_N3LO_contact(contact_no, val)    
    
    INTEGER, INTENT(IN) :: contact_no
    REAL(8), INTENT(IN) :: val
    
    IF (.NOT. chp_contact_format%set) THEN
       WRITE(CHP_ERR,"(A)") 'error(chp_set_N3LO_contact): you must set the format before the contact value'
       STOP
    END IF

! the N3LO contacts are input in units of 10^4/GeV^6
! this = 10^(4-3*6)/MeV^6 = 10^-14/MeV^6
    chp_DN3LO(contact_no)%val  = val * 1.D-14
    chp_DN3LO(contact_no)%set  = .TRUE.
    
    SELECT CASE(contact_no)
    CASE(1)
       IF (chp_contact_format%name == 'PW') chp_DN3LO(1)%name = 't1S0'
       IF (chp_contact_format%name == 'ST') chp_DN3LO(1)%name = 'D1->t1S0'
    CASE(2)
       IF (chp_contact_format%name == 'PW') chp_DN3LO(2)%name = '1S0'
       IF (chp_contact_format%name == 'ST') chp_DN3LO(2)%name = 'D2->1S0'
    CASE(3)
       IF (chp_contact_format%name == 'PW') chp_DN3LO(3)%name = '3P0'
       IF (chp_contact_format%name == 'ST') chp_DN3LO(3)%name = 'D3->3P0'
    CASE(4)
       IF (chp_contact_format%name == 'PW') chp_DN3LO(4)%name = '1P1'
       IF (chp_contact_format%name == 'ST') chp_DN3LO(4)%name = 'D4->1P1'
    CASE(5)
       IF (chp_contact_format%name == 'PW') chp_DN3LO(5)%name = '3P1'
       IF (chp_contact_format%name == 'ST') chp_DN3LO(5)%name = 'D5->3P1'
    CASE(6)
       IF (chp_contact_format%name == 'PW') chp_DN3LO(6)%name = 't3S1'
       IF (chp_contact_format%name == 'ST') chp_DN3LO(6)%name = 'D6->t3S1'
    CASE(7)
       IF (chp_contact_format%name == 'PW') chp_DN3LO(7)%name = '3S1'
       IF (chp_contact_format%name == 'ST') chp_DN3LO(7)%name = 'D7->3S1'
    CASE(8)
       IF (chp_contact_format%name == 'PW') chp_DN3LO(8)%name = '3D1'
       IF (chp_contact_format%name == 'ST') chp_DN3LO(8)%name = 'D8->3D1'
    CASE(9)
       IF (chp_contact_format%name == 'PW') chp_DN3LO(9)%name = 't3S1-3D1'
       IF (chp_contact_format%name == 'ST') chp_DN3LO(9)%name = 'D9->t3S1-3D1'
    CASE(10)
       IF (chp_contact_format%name == 'PW') chp_DN3LO(10)%name = '3S1-3D1'
       IF (chp_contact_format%name == 'ST') chp_DN3LO(10)%name = 'D10->3S1-3D1'
    CASE(11)
       IF (chp_contact_format%name == 'PW') chp_DN3LO(11)%name = '1D2'
       IF (chp_contact_format%name == 'ST') chp_DN3LO(11)%name = 'D11->1D2'
    CASE(12)
       IF (chp_contact_format%name == 'PW') chp_DN3LO(12)%name = '3D2'
       IF (chp_contact_format%name == 'ST') chp_DN3LO(12)%name = 'D12->3D2'
    CASE(13)
       IF (chp_contact_format%name == 'PW') chp_DN3LO(13)%name = '3P2'
       IF (chp_contact_format%name == 'ST') chp_DN3LO(13)%name = 'D13->3P2'
    CASE(14)
       IF (chp_contact_format%name == 'PW') chp_DN3LO(14)%name = '3P2-3F2'
       IF (chp_contact_format%name == 'ST') chp_DN3LO(14)%name = 'D14->3P2-3F2'
    CASE(15)
       IF (chp_contact_format%name == 'PW') chp_DN3LO(15)%name = '3D3'
       IF (chp_contact_format%name == 'ST') chp_DN3LO(15)%name = 'D15->3D3'
    CASE DEFAULT
       WRITE(CHP_ERR,"(A,I5)") 'error(chp_set_N3LO_contact): illegal contact nummber', contact_no
       STOP
    END SELECT
    
  END SUBROUTINE chp_set_N3LO_contact



  SUBROUTINE chp_get_N3LO_contact(contact_no, val, name)    
    
    INTEGER, INTENT(IN) :: contact_no
    REAL(8), INTENT(INOUT) :: val
    CHARACTER(LEN=12), INTENT(INOUT) :: name
    
    
    IF (.NOT. chp_DN3LO(contact_no)%set) THEN
       WRITE(CHP_ERR,"(A,I5,A)") 'chp_DN3LO ', contact_no, ' NOT SET'
       STOP
    END IF
! the N3LO contacts are input in units of 10^4/GeV^6
! this = 10^(4-3*6)/MeV^6 = 10^-14/MeV^6
    val = chp_DN3LO(contact_no)%val*1.0D+14
    name = chp_DN3LO(contact_no)%name
    
  END SUBROUTINE chp_get_N3LO_contact



  SUBROUTINE chp_set_c1(set_c1)
    
    REAL(8), INTENT(IN) :: set_c1
    
    chp_c1%val  = set_c1 * 1.0E-03 ! transform to MeV^-1
    chp_c1%set  = .TRUE.

  END SUBROUTINE chp_set_c1
  


  SUBROUTINE chp_get_c1(get_c1, name)
    
    REAL(8), INTENT(INOUT) :: get_c1
    CHARACTER(LEN=12), INTENT(INOUT) :: name
    
    IF (.NOT. chp_c1%set) THEN
       WRITE(CHP_ERR,"(A)") 'c1 NOT SET'
       STOP
    END IF
    
    get_c1 = chp_c1%val * 1000.0D0 ! transform to  GeV^-1
    name = chp_c1%name

  END SUBROUTINE chp_get_c1

  SUBROUTINE chp_set_c3(set_c3)
    
    REAL(8), INTENT(IN) :: set_c3
    
    chp_c3%val  = set_c3 * 1.0E-03 ! transform to MeV^-1
    chp_c3%set  = .TRUE.

  END SUBROUTINE chp_set_c3
  


  SUBROUTINE chp_get_c3(get_c3,name)
    
    REAL(8), INTENT(INOUT) :: get_c3
    CHARACTER(LEN=12), INTENT(INOUT) :: name
    
    IF (.NOT. chp_c3%set) THEN
       WRITE(CHP_ERR,"(A)") 'c3 NOT SET'
       STOP
    END IF
    get_c3 = chp_c3%val * 1000.0D0 ! transform to GeV^-1
    name = chp_c3%name
    
  END SUBROUTINE chp_get_c3
  
  SUBROUTINE chp_set_c4(set_c4)
    
    REAL(8), INTENT(IN) :: set_c4
    
    chp_c4%val  = set_c4 * 1.0E-03 ! transform to MeV^-1
    chp_c4%set  = .TRUE.

  END SUBROUTINE chp_set_c4
  


  SUBROUTINE chp_get_c4(get_c4, name)
    
    REAL(8), INTENT(INOUT) :: get_c4
    CHARACTER(LEN=12), INTENT(INOUT) :: name
    
    IF (.NOT. chp_c4%set) THEN
       WRITE(CHP_ERR,"(A)") 'c4 NOT SET'
       STOP
    END IF
    
    get_c4 = chp_c4%val * 1000.0D0 ! transform to GeV^-1
    name = chp_c4%name

  END SUBROUTINE chp_get_c4
  
  SUBROUTINE chp_set_c2(set_c2)
    
    REAL(8), INTENT(IN) :: set_c2
    
    chp_c2%val  = set_c2 * 1.0E-03 ! transform to MeV^-1
    chp_c2%set  = .TRUE.

  END SUBROUTINE chp_set_c2
  


  SUBROUTINE chp_get_c2(get_c2,name)
    
    REAL(8), INTENT(INOUT) :: get_c2
    CHARACTER(LEN=12), INTENT(INOUT) :: name
    
    IF (.NOT. chp_c2%set) THEN
       WRITE(CHP_ERR,"(A)") 'c2 NOT SET'
       STOP
    END IF
    
    get_c2 = chp_c2%val * 1000.0D0 ! transform to GeV^-1
    name = chp_c2%name

  END SUBROUTINE chp_get_c2

  SUBROUTINE chp_set_d1_plus_d2(set_d1_plus_d2)
    
    REAL(8), INTENT(IN) :: set_d1_plus_d2
    
    chp_d1_plus_d2%val  = set_d1_plus_d2 * 1.0D-6
    chp_d1_plus_d2%set  = .TRUE.
    
  END SUBROUTINE chp_set_d1_plus_d2
  


  SUBROUTINE chp_get_d1_plus_d2(get_d1_plus_d2,name)
    
    REAL(8), INTENT(INOUT) :: get_d1_plus_d2
    CHARACTER(LEN=12), INTENT(INOUT) :: name
    
    IF (.NOT. chp_d1_plus_d2%set) THEN
       WRITE(CHP_ERR,"(A)") 'd1+d2 NOT SET'
       STOP
    END IF
    
    get_d1_plus_d2 = chp_d1_plus_d2%val * 1.0D6 ! transform to GeV^-2
    name = chp_d1_plus_d2%name

  END SUBROUTINE chp_get_d1_plus_d2

  SUBROUTINE chp_set_d3(set_d3)
    
    REAL(8), INTENT(IN) :: set_d3
    
    chp_d3%val  = set_d3 * 1.0D-6
    chp_d3%set  = .TRUE.
    
  END SUBROUTINE chp_set_d3
  


  SUBROUTINE chp_get_d3(get_d3,name)
    
    REAL(8), INTENT(INOUT) :: get_d3
    CHARACTER(LEN=12), INTENT(INOUT) :: name
    
    IF (.NOT. chp_d3%set) THEN
       WRITE(CHP_ERR,"(A)") 'd3 NOT SET'
       STOP
    END IF
    
    get_d3 = chp_d3%val * 1.0D6 ! transform to GeV^-2
    name = chp_d3%name

  END SUBROUTINE chp_get_d3

  SUBROUTINE chp_set_d5(set_d5)
    
    REAL(8), INTENT(IN) :: set_d5
    
    chp_d5%val  = set_d5 * 1.0D-6
    chp_d5%set  = .TRUE.
    
  END SUBROUTINE chp_set_d5
  


  SUBROUTINE chp_get_d5(get_d5,name)
    
    REAL(8), INTENT(INOUT) :: get_d5
    CHARACTER(LEN=12), INTENT(INOUT) :: name
    
    IF (.NOT. chp_d5%set) THEN
       WRITE(CHP_ERR,"(A)") 'd5 NOT SET'
       STOP
    END IF
    
    get_d5 = chp_d5%val * 1.0D6 ! transform to GeV^-2
    name = chp_d5%name

  END SUBROUTINE chp_get_d5

  SUBROUTINE chp_set_d14_minus_d15(set_d14_minus_d15)
    
    REAL(8), INTENT(IN) :: set_d14_minus_d15
    
    chp_d14_minus_d15%val  = set_d14_minus_d15 * 1.0D-6
    chp_d14_minus_d15%set  = .TRUE.
    
  END SUBROUTINE chp_set_d14_minus_d15



  SUBROUTINE chp_get_d14_minus_d15(get_d14_minus_d15,name)
    
    REAL(8), INTENT(INOUT) :: get_d14_minus_d15
    CHARACTER(LEN=12), INTENT(INOUT) :: name
    
    IF (.NOT. chp_d14_minus_d15%set) THEN
       WRITE(CHP_ERR,"(A)") 'd14_minus_d15 NOT SET'
       STOP
    END IF
    
    get_d14_minus_d15 = chp_d14_minus_d15%val * 1.0D6 ! transform to GeV^-2
    name = chp_d14_minus_d15%name

  END SUBROUTINE chp_get_d14_minus_d15

  SUBROUTINE chp_set_e14(set_e14)
    
    REAL(8), INTENT(IN) :: set_e14
    
    chp_e14%val  = set_e14 * 1.0E-09 ! transform to MeV^-1
    chp_e14%set  = .TRUE.

  END SUBROUTINE chp_set_e14
  


  SUBROUTINE chp_get_e14(get_e14, name)
    
    REAL(8), INTENT(INOUT) :: get_e14
    CHARACTER(LEN=12), INTENT(INOUT) :: name
    
    IF (.NOT. chp_e14%set) THEN
       WRITE(CHP_ERR,"(A)") 'e14 NOT SET'
       STOP
    END IF
    
    get_e14 = chp_e14%val * 1.0D9 ! transform to  GeV^-3
    name = chp_e14%name

  END SUBROUTINE chp_get_e14

  SUBROUTINE chp_set_e15(set_e15)
    
    REAL(8), INTENT(IN) :: set_e15
    
    chp_e15%val  = set_e15 * 1.0E-09 ! transform to MeV^-1
    chp_e15%set  = .TRUE.

  END SUBROUTINE chp_set_e15
  


  SUBROUTINE chp_get_e15(get_e15, name)
    
    REAL(8), INTENT(INOUT) :: get_e15
    CHARACTER(LEN=12), INTENT(INOUT) :: name
    
    IF (.NOT. chp_e15%set) THEN
       WRITE(CHP_ERR,"(A)") 'e15 NOT SET'
       STOP
    END IF
    
    get_e15 = chp_e15%val * 1.0D9 ! transform to  GeV^-3
    name = chp_e15%name

  END SUBROUTINE chp_get_e15

  SUBROUTINE chp_set_e16(set_e16)
    
    REAL(8), INTENT(IN) :: set_e16
    
    chp_e16%val  = set_e16 * 1.0E-09 ! transform to MeV^-1
    chp_e16%set  = .TRUE.

  END SUBROUTINE chp_set_e16
  


  SUBROUTINE chp_get_e16(get_e16, name)
    
    REAL(8), INTENT(INOUT) :: get_e16
    CHARACTER(LEN=12), INTENT(INOUT) :: name
    
    IF (.NOT. chp_e16%set) THEN
       WRITE(CHP_ERR,"(A)") 'e16 NOT SET'
       STOP
    END IF
    
    get_e16 = chp_e16%val * 1.0D9 ! transform to  GeV^-3
    name = chp_e16%name

  END SUBROUTINE chp_get_e16

  SUBROUTINE chp_set_e17(set_e17)
    
    REAL(8), INTENT(IN) :: set_e17
    
    chp_e17%val  = set_e17 * 1.0E-09 ! transform to MeV^-1
    chp_e17%set  = .TRUE.

  END SUBROUTINE chp_set_e17
  


  SUBROUTINE chp_get_e17(get_e17, name)
    
    REAL(8), INTENT(INOUT) :: get_e17
    CHARACTER(LEN=12), INTENT(INOUT) :: name
    
    IF (.NOT. chp_e17%set) THEN
       WRITE(CHP_ERR,"(A)") 'e17 NOT SET'
       STOP
    END IF
    
    get_e17 = chp_e17%val * 1.0D9 ! transform to  GeV^-3
    name = chp_e17%name

  END SUBROUTINE chp_get_e17

  SUBROUTINE chp_set_e18(set_e18)
    
    REAL(8), INTENT(IN) :: set_e18
    
    chp_e18%val  = set_e18 * 1.0E-09 ! transform to MeV^-1
    chp_e18%set  = .TRUE.

  END SUBROUTINE chp_set_e18
  


  SUBROUTINE chp_get_e18(get_e18, name)
    
    REAL(8), INTENT(INOUT) :: get_e18
    CHARACTER(LEN=12), INTENT(INOUT) :: name
    
    IF (.NOT. chp_e18%set) THEN
       WRITE(CHP_ERR,"(A)") 'e18 NOT SET'
       STOP
    END IF
    
    get_e18 = chp_e18%val * 1.0D9 ! transform to  GeV^-3
    name = chp_e18%name

  END SUBROUTINE chp_get_e18

  SUBROUTINE chp_set_lambda(set_lambda)
    
    REAL(8), INTENT(IN) :: set_lambda

    chp_lambda%val  = set_lambda
    chp_lambda%set  = .TRUE.

  END SUBROUTINE chp_set_lambda

  SUBROUTINE chp_get_lambda(get_lambda)
    
    REAL(8), INTENT(INOUT) :: get_lambda

    IF (.NOT. chp_lambda%set) THEN
       WRITE(CHP_ERR,"(A)") 'lambda NOT SET'
       STOP
    END IF
    
    get_lambda = chp_lambda%val
  
  END SUBROUTINE chp_get_lambda

  SUBROUTINE chp_set_gA(set_gA)
    
    REAL(8), INTENT(IN) :: set_gA

! This depends on gA, so make it unset when gA is changed
    IF(.NOT. chp_gA%set .OR. chp_gA%val /= set_gA) THEN
       chp_2PE_2loop_int_WC_DR_data_set  = .FALSE.
       chp_2PE_2loop_int_WC_SFR_data_set  = .FALSE.
    END IF
    chp_gA%val  = set_gA
    chp_gA%set  = .TRUE.

  END SUBROUTINE chp_set_gA



  SUBROUTINE chp_get_gA(get_gA)
    
    REAL(8), INTENT(INOUT) :: get_gA
    
    IF (.NOT. chp_gA%set) THEN
       WRITE(CHP_ERR,"(A)") 'gA NOT SET'
       STOP
    END IF
    
    get_gA = chp_gA%val
  
  END SUBROUTINE chp_get_gA

  SUBROUTINE chp_set_fpi(set_fpi)
    
    REAL(8), INTENT(IN) :: set_fpi

    chp_fpi%val  = set_fpi
    chp_fpi%set  = .TRUE.
    
  END SUBROUTINE chp_set_fpi



  SUBROUTINE chp_get_fpi(get_fpi)
    
    REAL(8), INTENT(INOUT) :: get_fpi

    IF (.NOT. chp_fpi%set) THEN
       WRITE(CHP_ERR,"(A)") 'fpi NOT SET'
       STOP
    END IF
    
    get_fpi = chp_fpi%val
    
  END SUBROUTINE chp_get_fpi

  SUBROUTINE chp_set_fine_structure(set_fine_structure)
    
    REAL(8), INTENT(IN) :: set_fine_structure

    chp_fine_structure%val  = set_fine_structure
    chp_fine_structure%set  = .TRUE.
    
  END SUBROUTINE chp_set_fine_structure

  SUBROUTINE chp_get_fine_structure(get_fine_structure)
    
    REAL(8), INTENT(INOUT) :: get_fine_structure

    IF (.NOT. chp_fine_structure%set) THEN
       WRITE(CHP_ERR,"(A)") 'fine_structure NOT SET'
       STOP
    END IF
    
    get_fine_structure = chp_fine_structure%val
    
  END SUBROUTINE chp_get_fine_structure



  SUBROUTINE chp_set_chiral_Ct(set)
    
    INTEGER, INTENT(IN) :: set

    IF (set < 0 .OR. set > 1) THEN
       ERROR STOP 'Expected 0 or 1 in chp_set_chiral_Ct'
    END IF
    chiral_Ct%val  = set
    chiral_Ct%set  = .TRUE.
    
  END SUBROUTINE chp_set_chiral_Ct

  SUBROUTINE chp_get_chiral_Ct(get)
    
    INTEGER, INTENT(INOUT) :: get

    IF (.NOT. chiral_Ct%set) THEN
       WRITE(CHP_ERR,"(A)") 'chiral_Ct NOT SET'
       STOP
    END IF
    
    get = chiral_Ct%val
    
  END SUBROUTINE chp_get_chiral_Ct

  SUBROUTINE chp_set_chiral_Ct_CIB(set)
    
    INTEGER, INTENT(IN) :: set

    IF (set < 0 .OR. set > 1) THEN
       ERROR STOP 'Expected 0 or 1 in chp_set_chiral_Ct_CIB'
    END IF
    chiral_Ct_CIB%val  = set
    chiral_Ct_CIB%set  = .TRUE.
    
  END SUBROUTINE chp_set_chiral_Ct_CIB

  SUBROUTINE chp_get_chiral_Ct_CIB(get)
    
    INTEGER, INTENT(INOUT) :: get

    IF (.NOT. chiral_Ct_CIB%set) THEN
       WRITE(CHP_ERR,"(A)") 'chiral_Ct_CIB NOT SET'
       STOP
    END IF
    
    get = chiral_Ct_CIB%val
    
  END SUBROUTINE chp_get_chiral_Ct_CIB

  SUBROUTINE chp_set_chiral_C(set)
    
    INTEGER, INTENT(IN) :: set

    IF (set < 0 .OR. set > 1) THEN
       ERROR STOP 'Expected 0 or 1 in chp_set_chiral_C'
    END IF
    chiral_C%val  = set
    chiral_C%set  = .TRUE.
    
  END SUBROUTINE chp_set_chiral_C

  SUBROUTINE chp_get_chiral_C(get)
    
    INTEGER, INTENT(INOUT) :: get

    IF (.NOT. chiral_C%set) THEN
       WRITE(CHP_ERR,"(A)") 'chiral_C NOT SET'
       STOP
    END IF
    
    get = chiral_C%val
    
  END SUBROUTINE chp_get_chiral_C

  SUBROUTINE chp_set_chiral_D(set)
    
    INTEGER, INTENT(IN) :: set

    IF (set < 0 .OR. set > 1) THEN
       ERROR STOP 'Expected 0 or 1 in chp_set_chiral_D'
    END IF
    chiral_D%val  = set
    chiral_D%set  = .TRUE.
    
  END SUBROUTINE chp_set_chiral_D

  SUBROUTINE chp_get_chiral_D(get)
    
    INTEGER, INTENT(INOUT) :: get

    IF (.NOT. chiral_D%set) THEN
       WRITE(CHP_ERR,"(A)") 'chiral_D NOT SET'
       STOP
    END IF
    
    get = chiral_D%val
    
  END SUBROUTINE chp_get_chiral_D

  SUBROUTINE chp_set_chiral_1PE(set)
    
    INTEGER, INTENT(IN) :: set

    IF (set < 0 .OR. set > 1) THEN
       ERROR STOP 'Expected 0 or 1 in chp_set_chiral_1PE'
    END IF
    chiral_1PE%val  = set
    chiral_1PE%set  = .TRUE.
    
  END SUBROUTINE chp_set_chiral_1PE

  SUBROUTINE chp_get_chiral_1PE(get)
    
    INTEGER, INTENT(INOUT) :: get

    IF (.NOT. chiral_1PE%set) THEN
       WRITE(CHP_ERR,"(A)") 'chiral_1PE NOT SET'
       STOP
    END IF
    
    get = chiral_1PE%val
    
  END SUBROUTINE chp_get_chiral_1PE

  SUBROUTINE chp_set_chiral_1PE_CIB(set)
    
    INTEGER, INTENT(IN) :: set

    IF (set < 0 .OR. set > 1) THEN
       ERROR STOP 'Expected 0 or 1 in chp_set_chiral_1PE_CIB'
    END IF
    chiral_1PE_CIB%val  = set
    chiral_1PE_CIB%set  = .TRUE.
    
  END SUBROUTINE chp_set_chiral_1PE_CIB

  SUBROUTINE chp_get_chiral_1PE_CIB(get)
    
    INTEGER, INTENT(INOUT) :: get

    IF (.NOT. chiral_1PE_CIB%set) THEN
       WRITE(CHP_ERR,"(A)") 'chiral_1PE_CIB NOT SET'
       STOP
    END IF
    
    get = chiral_1PE_CIB%val
    
  END SUBROUTINE chp_get_chiral_1PE_CIB

  SUBROUTINE chp_set_chiral_1PE_gamma(set)
    
    INTEGER, INTENT(IN) :: set

    IF (set < 0 .OR. set > 1) THEN
       ERROR STOP 'Expected 0 or 1 in chp_set_chiral_1PE_gamma'
    END IF
    chiral_1PE_gamma%val  = set
    chiral_1PE_gamma%set  = .TRUE.
    
  END SUBROUTINE chp_set_chiral_1PE_gamma

  SUBROUTINE chp_get_chiral_1PE_gamma(get)
    
    INTEGER, INTENT(INOUT) :: get

    IF (.NOT. chiral_1PE_gamma%set) THEN
       WRITE(CHP_ERR,"(A)") 'chiral_1PE_gamma NOT SET'
       STOP
    END IF
    
    get = chiral_1PE_gamma%val
    
  END SUBROUTINE chp_get_chiral_1PE_gamma

  SUBROUTINE chp_set_chiral_1PE_relcorr(set)
    
    INTEGER, INTENT(IN) :: set

    IF (set < 0 .OR. set > 1) THEN
       ERROR STOP 'Expected 0 or 1 in chp_set_chiral_1PE_relcorr'
    END IF
    chiral_1PE_relcorr%val  = set
    chiral_1PE_relcorr%set  = .TRUE.
    
  END SUBROUTINE chp_set_chiral_1PE_relcorr

  SUBROUTINE chp_get_chiral_1PE_relcorr(get)
    
    INTEGER, INTENT(INOUT) :: get

    IF (.NOT. chiral_1PE_relcorr%set) THEN
       WRITE(CHP_ERR,"(A)") 'chiral_1PE_relcorr NOT SET'
       STOP
    END IF
    
    get = chiral_1PE_relcorr%val
    
  END SUBROUTINE chp_get_chiral_1PE_relcorr

  SUBROUTINE chp_set_chiral_2PE_1loop_0(set)
    
    INTEGER, INTENT(IN) :: set

    IF (set < 0 .OR. set > 1) THEN
       ERROR STOP 'Expected 0 or 1 in chp_set_chiral_2PE_1loop_0'
    END IF
    chiral_2PE_1loop_0%val  = set
    chiral_2PE_1loop_0%set  = .TRUE.
    
  END SUBROUTINE chp_set_chiral_2PE_1loop_0

  SUBROUTINE chp_get_chiral_2PE_1loop_0(get)
    
    INTEGER, INTENT(INOUT) :: get

    IF (.NOT. chiral_2PE_1loop_0%set) THEN
       WRITE(CHP_ERR,"(A)") 'chiral_2PE_1loop_0 NOT SET'
       STOP
    END IF
    
    get = chiral_2PE_1loop_0%val
    
  END SUBROUTINE chp_get_chiral_2PE_1loop_0

  SUBROUTINE chp_set_chiral_2PE_1loop_d(set)
    
    INTEGER, INTENT(IN) :: set

    IF (set < 0 .OR. set > 1) THEN
       ERROR STOP 'Expected 0 or 1 in chp_set_chiral_2PE_1loop_d'
    END IF
    chiral_2PE_1loop_d%val  = set
    chiral_2PE_1loop_d%set  = .TRUE.
    
  END SUBROUTINE chp_set_chiral_2PE_1loop_d

  SUBROUTINE chp_get_chiral_2PE_1loop_d(get)
    
    INTEGER, INTENT(INOUT) :: get

    IF (.NOT. chiral_2PE_1loop_d%set) THEN
       WRITE(CHP_ERR,"(A)") 'chiral_2PE_1loop_d NOT SET'
       STOP
    END IF
    
    get = chiral_2PE_1loop_d%val
    
  END SUBROUTINE chp_get_chiral_2PE_1loop_d

  SUBROUTINE chp_set_chiral_2PE_1loop_r(set)
    
    INTEGER, INTENT(IN) :: set

    IF (set < 0 .OR. set > 1) THEN
       ERROR STOP 'Expected 0 or 1 in chp_set_chiral_2PE_1loop_r'
    END IF
    chiral_2PE_1loop_r%val  = set
    chiral_2PE_1loop_r%set  = .TRUE.
    
  END SUBROUTINE chp_set_chiral_2PE_1loop_r

  SUBROUTINE chp_get_chiral_2PE_1loop_r(get)
    
    INTEGER, INTENT(INOUT) :: get

    IF (.NOT. chiral_2PE_1loop_r%set) THEN
       WRITE(CHP_ERR,"(A)") 'chiral_2PE_1loop_r NOT SET'
       STOP
    END IF
    
    get = chiral_2PE_1loop_r%val
    
  END SUBROUTINE chp_get_chiral_2PE_1loop_r

  SUBROUTINE chp_set_chiral_2PE_1loop_r_mode(set)
    
    INTEGER, INTENT(IN) :: set

    IF (set < chiral_2PE_1loop_r_mode_EM2011 .OR. set > chiral_2PE_1loop_r_mode_Ep2015) THEN
       ERROR STOP 'Bad value in chp_set_chiral_2PE_1loop_r_mode'
    END IF
    chiral_2PE_1loop_r_mode%val  = set
    chiral_2PE_1loop_r_mode%set  = .TRUE.
    
  END SUBROUTINE chp_set_chiral_2PE_1loop_r_mode

  SUBROUTINE chp_get_chiral_2PE_1loop_r_mode(get)
    
    INTEGER, INTENT(INOUT) :: get

    IF (.NOT. chiral_2PE_1loop_r_mode%set) THEN
       WRITE(CHP_ERR,"(A)") 'chiral_2PE_1loop_r_mode NOT SET'
       STOP
    END IF
    
    get = chiral_2PE_1loop_r_mode%val
    
  END SUBROUTINE chp_get_chiral_2PE_1loop_r_mode

  SUBROUTINE chp_set_chiral_2PE_1loop_dd(set)
    
    INTEGER, INTENT(IN) :: set

    IF (set < 0 .OR. set > 1) THEN
       ERROR STOP 'Expected 0 or 1 in chp_set_chiral_2PE_1loop_dd'
    END IF
    chiral_2PE_1loop_dd%val  = set
    chiral_2PE_1loop_dd%set  = .TRUE.
    
  END SUBROUTINE chp_set_chiral_2PE_1loop_dd

  SUBROUTINE chp_get_chiral_2PE_1loop_dd(get)
    
    INTEGER, INTENT(INOUT) :: get

    IF (.NOT. chiral_2PE_1loop_dd%set) THEN
       WRITE(CHP_ERR,"(A)") 'chiral_2PE_1loop_dd NOT SET'
       STOP
    END IF
    
    get = chiral_2PE_1loop_dd%val
    
  END SUBROUTINE chp_get_chiral_2PE_1loop_dd

  SUBROUTINE chp_set_chiral_2PE_1loop_dr(set)
    
    INTEGER, INTENT(IN) :: set

    IF (set < 0 .OR. set > 1) THEN
       ERROR STOP 'Expected 0 or 1 in chp_set_chiral_2PE_1loop_dr'
    END IF
    chiral_2PE_1loop_dr%val  = set
    chiral_2PE_1loop_dr%set  = .TRUE.
    
  END SUBROUTINE chp_set_chiral_2PE_1loop_dr

  SUBROUTINE chp_get_chiral_2PE_1loop_dr(get)
    
    INTEGER, INTENT(INOUT) :: get

    IF (.NOT. chiral_2PE_1loop_dr%set) THEN
       WRITE(CHP_ERR,"(A)") 'chiral_2PE_1loop_dr NOT SET'
       STOP
    END IF
    
    get = chiral_2PE_1loop_dr%val
    
  END SUBROUTINE chp_get_chiral_2PE_1loop_dr

  SUBROUTINE chp_set_chiral_2PE_1loop_rr(set)
    
    INTEGER, INTENT(IN) :: set

    IF (set < 0 .OR. set > 1) THEN
       ERROR STOP 'Expected 0 or 1 in chp_set_chiral_2PE_1loop_rr'
    END IF
    chiral_2PE_1loop_rr%val  = set
    chiral_2PE_1loop_rr%set  = .TRUE.
    
  END SUBROUTINE chp_set_chiral_2PE_1loop_rr

  SUBROUTINE chp_get_chiral_2PE_1loop_rr(get)
    
    INTEGER, INTENT(INOUT) :: get

    IF (.NOT. chiral_2PE_1loop_rr%set) THEN
       WRITE(CHP_ERR,"(A)") 'chiral_2PE_1loop_rr NOT SET'
       STOP
    END IF
    
    get = chiral_2PE_1loop_rr%val
    
  END SUBROUTINE chp_get_chiral_2PE_1loop_rr

  SUBROUTINE chp_set_chiral_2PE_2loop(set)
    
    INTEGER, INTENT(IN) :: set

    IF (set < 0 .OR. set > 1) THEN
       ERROR STOP 'Expected 0 or 1 in chp_set_chiral_2PE_2loop'
    END IF
    chiral_2PE_2loop%val  = set
    chiral_2PE_2loop%set  = .TRUE.
    
  END SUBROUTINE chp_set_chiral_2PE_2loop

  SUBROUTINE chp_get_chiral_2PE_2loop(get)
    
    INTEGER, INTENT(INOUT) :: get

    IF (.NOT. chiral_2PE_2loop%set) THEN
       WRITE(CHP_ERR,"(A)") 'chiral_2PE_2loop NOT SET'
       STOP
    END IF
    
    get = chiral_2PE_2loop%val
    
  END SUBROUTINE chp_get_chiral_2PE_2loop

  SUBROUTINE chp_set_chiral_2PE_2loop_int(set)
    
    INTEGER, INTENT(IN) :: set

    IF (set < 0 .OR. set > 1) THEN
       ERROR STOP 'Expected 0 or 1 in chp_set_chiral_2PE_2loop_int'
    END IF
    chiral_2PE_2loop_int%val  = set
    chiral_2PE_2loop_int%set  = .TRUE.
    
  END SUBROUTINE chp_set_chiral_2PE_2loop_int

  SUBROUTINE chp_get_chiral_2PE_2loop_int(get)
    
    INTEGER, INTENT(INOUT) :: get

    IF (.NOT. chiral_2PE_2loop_int%set) THEN
       WRITE(CHP_ERR,"(A)") 'chiral_2PE_2loop_int NOT SET'
       STOP
    END IF
    
    get = chiral_2PE_2loop_int%val
    
  END SUBROUTINE chp_get_chiral_2PE_2loop_int

  SUBROUTINE chp_set_chiral_2PE_CSB_correct_mass(set)
    
    INTEGER, INTENT(IN) :: set

    IF (set < 0 .OR. set > 1) THEN
       ERROR STOP 'Expected 0 or 1 in chp_set_chiral_2PE_CSB_correct_mass'
    END IF
    chiral_2PE_CSB_correct_mass%val  = set
    chiral_2PE_CSB_correct_mass%set  = .TRUE.
    
  END SUBROUTINE chp_set_chiral_2PE_CSB_correct_mass

  SUBROUTINE chp_get_chiral_2PE_CSB_correct_mass(get)
    
    INTEGER, INTENT(INOUT) :: get

    IF (.NOT. chiral_2PE_CSB_correct_mass%set) THEN
       WRITE(CHP_ERR,"(A)") 'chiral_2PE_CSB_correct_mass NOT SET'
       STOP
    END IF
    
    get = chiral_2PE_CSB_correct_mass%val
    
  END SUBROUTINE chp_get_chiral_2PE_CSB_correct_mass

  SUBROUTINE chp_set_chiral_minimal_relativity(set)
    
    INTEGER, INTENT(IN) :: set

    IF (set < 0 .OR. set > 1) THEN
       ERROR STOP 'Expected 0 or 1 in chp_set_chiral_minimal_relativity'
    END IF
    chiral_minimal_relativity%val  = set
    chiral_minimal_relativity%set  = .TRUE.
    
  END SUBROUTINE chp_set_chiral_minimal_relativity

  SUBROUTINE chp_get_chiral_minimal_relativity(get)
    
    INTEGER, INTENT(INOUT) :: get

    IF (.NOT. chiral_minimal_relativity%set) THEN
       WRITE(CHP_ERR,"(A)") 'chiral_minimal_relativity NOT SET'
       STOP
    END IF
    
    get = chiral_minimal_relativity%val
    
  END SUBROUTINE chp_get_chiral_minimal_relativity

  SUBROUTINE chp_set_chiral_kamada_glockle_transform(set)
    
    INTEGER, INTENT(IN) :: set

    IF (set < 0 .OR. set > 1) THEN
       ERROR STOP 'Expected 0 or 1 in chp_set_chiral_kamada_glockle_transform'
    END IF
    chiral_kamada_glockle_transform%val  = set
    chiral_kamada_glockle_transform%set  = .TRUE.
    
  END SUBROUTINE chp_set_chiral_kamada_glockle_transform

  SUBROUTINE chp_get_chiral_kamada_glockle_transform(get)
    
    INTEGER, INTENT(INOUT) :: get

    IF (.NOT. chiral_kamada_glockle_transform%set) THEN
       WRITE(CHP_ERR,"(A)") 'chiral_kamada_glockle_transform NOT SET'
       STOP
    END IF
    
    get = chiral_kamada_glockle_transform%val
    
  END SUBROUTINE chp_get_chiral_kamada_glockle_transform

END MODULE idaho_chiral_potential

