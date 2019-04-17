module constants

    implicit none
    INTEGER(KIND=8),PARAMETER  :: NDIM = 20000 
    COMPLEX(KIND=8),PARAMETER  :: Imu =(0.0D0,1.0D0)
    REAL(KIND=8), PARAMETER    :: hbar=1.0D0, kb=1.0D0
    REAL(KIND=8), PARAMETER    :: A2au=1.889726, au2A=0.529177
    REAL(KIND=8), PARAMETER    :: eV2au=0.036749, au2eV=27.211383,invcm2au = 4.559489400000001e-06
    Character(LEN=40), PARAMETER:: sep='----------------------------------------'

    save  

end module constants
