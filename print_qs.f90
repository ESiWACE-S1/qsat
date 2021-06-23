PROGRAM PRINTQS
IMPLICIT NONE

INTERFACE
  FUNCTION dales_formula(T)
    REAL :: dales_formula
    REAL, INTENT(IN) :: T
  END FUNCTION dales_formula

  FUNCTION dales_formula_ice(T)
    REAL :: dales_formula_ice
    REAL, INTENT(IN) :: T
  END FUNCTION dales_formula_ice

  FUNCTION Huang(T)
    REAL :: Huang
    REAL, INTENT(IN) :: T
  END FUNCTION Huang

  FUNCTION Huang_Ice(T)
    REAL :: Huang_Ice
    REAL, INTENT(IN) :: T
  END FUNCTION Huang_Ice

  FUNCTION WagnerPruss(T)
    REAL ::  WagnerPruss
    REAL, INTENT(IN) :: T
  END FUNCTION WagnerPruss

  FUNCTION WagnerPruss_Ice(T)
    REAL ::  WagnerPruss_Ice
    REAL, INTENT(IN) :: T
  END FUNCTION WagnerPruss_Ice

  FUNCTION Improved_Magnus(T)
    REAL ::  Improved_Magnus
    REAL, INTENT(IN) :: T
  END FUNCTION Improved_Magnus

  FUNCTION Improved_Magnus_Ice(T)
    REAL ::  Improved_Magnus_Ice
    REAL, INTENT(IN) :: T
  END FUNCTION Improved_Magnus_Ice
END INTERFACE

INTEGER :: n
REAL :: T  ! temperature [K]

WRITE (*,*) "Over water"
WRITE (*,*) "n, WagnerPruss(T), Improved_Magnus(T), Huang(T), dales_formula"
DO n=-40,40
  T = 273.15 - n
  WRITE (*,*) n, WagnerPruss(T), Improved_Magnus(T), Huang(T), dales_formula(T) !, WagnerPruss_Ice(T), Improved_Magnus_Ice(T)
END DO

WRITE (*,*) "Over Ice"
WRITE (*,*) "n, WagnerPruss(T), Improved_Magnus(T), Huang(T), dales_formula"
DO n=-40,40
  T = 273.15 - n
  WRITE (*,*) n, WagnerPruss_Ice(T), Improved_Magnus_Ice(T), Huang_Ice(T), dales_formula_ice(T)
END DO

END PROGRAM

! Improved_Magnus [Pa] saturation vapor pressure over water
! T [K] temperature
FUNCTION Improved_Magnus(T)

IMPLICIT NONE
REAL, INTENT(IN) :: T
REAL  :: Improved_Magnus

Improved_Magnus = 610.94 * EXP((17.625 * (T-273.15)) / (T + 243.04 - 273.15))
END FUNCTION


! Improved_Magnus [Pa] saturation vapor pressure over ice
! T [K] temperature
FUNCTION Improved_Magnus_Ice(T)

IMPLICIT NONE
REAL, INTENT(IN) :: T
REAL  :: Improved_Magnus_Ice

Improved_Magnus_Ice = 611.21 * EXP((22.587 * (T-273.15)) / (T + 273.86 - 273.15))
END FUNCTION Improved_Magnus_Ice


! WagnerPruss [Pa] saturation vapor pressure of water
! T [K] temperature
FUNCTION WagnerPruss(T)

IMPLICIT NONE
REAL, INTENT(IN) :: T
REAL :: WagnerPruss

REAL, PARAMETER :: TC = 647.096 ! critical point temperature [K]

REAL :: x

x = (1 - (T/TC))
WagnerPruss = 22064000.0 * EXP( (TC/T) * ( &
  - 7.85951783 * x &
  + 1.84408259 * x**1.5 &
  - 11.7866497 * x**3 &
  + 22.6807411 * x**3.5 &
  - 15.9618719 * x**4 &
  + 1.80122502 * x**7.5 &
  ))

END FUNCTION WagnerPruss


! WagnerPrussIce [Pa] saturation vapor pressure over ice 
! T [K] temperature
! NOTE: those values are all wrong.. check formula!
FUNCTION WagnerPruss_Ice(T)

IMPLICIT NONE
REAL, INTENT(IN) :: T
REAL :: WagnerPruss_Ice

REAL, PARAMETER :: TT = 273.16 ! Triple point temperature of water [K]

WagnerPruss_Ice = 611.657 * EXP( &
  (TT/T) * ( &
    + 21.2144006 * (T/TT)**0.00333333333 &
    + 27.3203819 * (T/TT)**1.20666667 &
    +  6.1059813 * (T/TT)**1.70333333 &
  ))

END FUNCTION WagnerPruss_Ice


! Huang [Pa] saturation vapor pressure over water
! T [K] temperature
FUNCTION Huang(T)

IMPLICIT NONE
REAL, INTENT(IN) :: T
REAL  :: Huang

REAL :: TC

TC = T - 273.15  !  TC [degC]

Huang = EXP(34.494 - 4924.99 / (TC + 237.1)) /  (TC + 105)**1.57

END FUNCTION Huang


! Huang [Pa] saturation vapor pressure over water
! T [K] temperature
FUNCTION Huang_Ice(T)

IMPLICIT NONE
REAL, INTENT(IN) :: T
REAL  :: Huang_Ice

REAL :: TC

TC = T - 273.15  !  TC [degC]

Huang_Ice = EXP(43.494 -  6545.8 / (TC + 278)) / (TC + 868)**2

END FUNCTION Huang_Ice

! Current implementation in DALES
FUNCTION dales_formula(T)
  IMPLICIT NONE
  INTERFACE
    FUNCTION sattab_l(m)
      REAL :: sattab_l
      INTEGER, INTENT(IN) :: m
    END FUNCTION sattab_l

    FUNCTION sattab_i(m)
      REAL :: sattab_i
      INTEGER, INTENT(IN) :: m
    END FUNCTION sattab_i
  END INTERFACE

  real,parameter :: rd      = 287.04           !<    * gas constant for dry air.
  real,parameter :: rv      = 461.5            !<    * gas constant for water vapor.
  real,parameter :: tup     = 268.             !<    * Temperature range over which mixed phase occurs (high)
  real,parameter :: tdn     = 253.             !<    * Temperature range over which mixed phase occurs (low)

  real, intent(in) :: T
  real :: dales_formula
  real :: esl, qvsl, qvsi

  real :: ilratio
  real :: esi, thi, tlo
  integer :: tlonr, thinr

  tlonr=int(T*5.)
  tlo = tlonr - T*5.
  esl=(tlo + 1.0)*sattab_l(tlonr - 750)-tlo*sattab_l(tlonr - 750 + 1)
  esi=(tlo + 1.0)*sattab_i(tlonr - 750)-tlo*sattab_i(tlonr - 750 + 1)
  dales_formula = esl
end function dales_formula

FUNCTION dales_formula_ice(T)
  IMPLICIT NONE
  INTERFACE
    FUNCTION sattab_l(m)
      REAL :: sattab_l
      INTEGER, INTENT(IN) :: m
    END FUNCTION sattab_l

    FUNCTION sattab_i(m)
      REAL :: sattab_i
      INTEGER, INTENT(IN) :: m
    END FUNCTION sattab_i
  END INTERFACE

  real,parameter :: rd      = 287.04           !<    * gas constant for dry air.
  real,parameter :: rv      = 461.5            !<    * gas constant for water vapor.
  real,parameter :: tup     = 268.             !<    * Temperature range over which mixed phase occurs (high)
  real,parameter :: tdn     = 253.             !<    * Temperature range over which mixed phase occurs (low)

  real, intent(in) :: T
  real :: dales_formula_ice
  real :: esl, qvsl, qvsi

  real :: ilratio
  real :: esi, thi, tlo
  integer :: tlonr, thinr

  tlonr=int(T*5.)
  tlo = tlonr - T*5.
  esl=(tlo + 1.0)*sattab_l(tlonr - 750)-tlo*sattab_l(tlonr - 750 + 1)
  esi=(tlo + 1.0)*sattab_i(tlonr - 750)-tlo*sattab_i(tlonr - 750 + 1)
  dales_formula_ice = esl
END FUNCTION dales_formula_ice

FUNCTION sattab_l(m)
  implicit none
  integer, intent(in) :: m
  real :: sattab_l, t
  
  t = 150.+0.2*m
  sattab_l = exp(54.842763-6763.22/t-4.21*log(t)+0.000367*t+&
       tanh(0.0415*(t-218.8))*(53.878-1331.22/t-9.44523*log(t)+ 0.014025*t))
END FUNCTION sattab_l

FUNCTION sattab_i(m)
  implicit none
  integer, intent(in) :: m
  real :: sattab_i, t

  t = 150.+0.2*m
  sattab_i = exp(9.550426-5723.265/t+3.53068*log(t)-0.00728332*t)
END FUNCTION sattab_i

