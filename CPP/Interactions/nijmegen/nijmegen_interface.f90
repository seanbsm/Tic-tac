module nijmegen

  implicit none
  
contains

  subroutine vnijm(pout, pin, coup_, S, j_, T, tz, POT) bind(C,name="nijmegen_fort_interface")

    REAL(8) ,INTENT(IN)    :: pout 
    REAL(8) ,INTENT(IN)    :: pin  
    LOGICAL,INTENT(IN)    :: coup_  
    INTEGER,INTENT(IN)    :: S      
    INTEGER,INTENT(IN)    :: j_     
    INTEGER,INTENT(IN)    :: T   ! not used/tested in interface. only in pnijm   
    INTEGER,INTENT(IN)    :: tz     
    REAL(8),INTENT(INOUT) :: POT(6) 

    ! interaction common block 
    INTEGER  :: inn, j, IDPAR
    REAL(8)  :: v,xmev,ymev
    CHARACTER (LEN=4) :: label
    LOGICAL :: coup, sing, trip, heform, endep

    COMMON/CHOICE/IDPAR
    COMMON /cnn/ inn
    COMMON /cpot/ v(6),xmev,ymev
    COMMON /cstate/ j,heform, sing, trip, coup, endep, label

    coup = .FALSE. ; sing = .FALSE. ; trip = .FALSE. ; heform = .FALSE.
    ! ignore label
    !IDPAR = 0, 1, 2 for Nijm 93, I, II(local), respectively.
    IDPAR = 1
    
    SELECT CASE(tz)
    CASE(-1)
       ! pp
       inn = 1
    CASE(0)
       ! pn
       inn = 2
    CASE(+1)
       ! nn
       inn = 3
    END SELECT

    IF (coup_)  coup = .TRUE.
    IF (S == 1) trip = .TRUE.
    IF (S == 0) sing = .TRUE.

    j = j_
    
    xmev = pout
    ymev = pin 
    
    call nijm
    POT = v
    
  end subroutine vnijm
 
end module nijmegen
