        !COMPILER-GENERATED INTERFACE MODULE: Fri Jan 09 15:25:07 2015
        MODULE CALCOLOCALL__genmod
          INTERFACE 
            SUBROUTINE CALCOLOCALL(NVV,NPP,NPPMAT,NPPTIME,IDXX,PAR,     &
     &TMP_C_P,FUT,KKSTRIKE,OPTC,OPTP)
              INTEGER(KIND=4) :: NPPTIME
              INTEGER(KIND=4) :: NPPMAT
              INTEGER(KIND=4) :: NPP
              INTEGER(KIND=4) :: NVV
              INTEGER(KIND=4) :: IDXX(1:NPPTIME)
              REAL(KIND=8) :: PAR(NVV)
              REAL(KIND=8) :: TMP_C_P(1:NPPMAT,1:NPP)
              REAL(KIND=8) :: FUT(1:NPPMAT)
              REAL(KIND=8) :: KKSTRIKE(1:NPP)
              REAL(KIND=8) :: OPTC(1:NPPMAT,1:NPP)
              REAL(KIND=8) :: OPTP(1:NPPMAT,1:NPP)
            END SUBROUTINE CALCOLOCALL
          END INTERFACE 
        END MODULE CALCOLOCALL__genmod
