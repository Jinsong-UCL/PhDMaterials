        !COMPILER-GENERATED INTERFACE MODULE: Fri Jan 09 09:55:38 2015
        MODULE FUNZHESTON__genmod
          INTERFACE 
            SUBROUTINE FUNZHESTON(TAU,K,EPSIL,THETA,RHO12,K1,V0,DELTAA, &
     &NN,A)
              INTEGER(KIND=4) :: NN
              REAL(KIND=8) :: TAU
              REAL(KIND=8) :: K(1:NN+1)
              REAL(KIND=8) :: EPSIL(2)
              REAL(KIND=8) :: THETA(2)
              REAL(KIND=8) :: RHO12(2)
              REAL(KIND=8) :: K1(2)
              REAL(KIND=8) :: V0(2)
              REAL(KIND=8) :: DELTAA
              COMPLEX(KIND=8) :: A(1:NN+1)
            END SUBROUTINE FUNZHESTON
          END INTERFACE 
        END MODULE FUNZHESTON__genmod
