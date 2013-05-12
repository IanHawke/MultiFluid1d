module variables_mod

  implicit none

  integer, parameter :: &
       NF  = 2, & ! Number of fluids
       NFF = 3, & ! Number of fluid pairs, = NF (NF + 1) / 2
       ND  = 1, & ! Number of spatial dimensions
       NV  = 2, & ! Number of variables, = ND + 1
       NEVOLVE = 4, & ! Number of evolved variables, = NF * NV
       LWA = 33   ! Size of the work array for minpack; depends on NFF as
                  ! LWA >= (NFF (3 NFF + 13)) / 2

  integer, parameter :: &
       N_UP_FILE = 10, &
       N_DOWN_FILE = 11, &
       MU_UP_FILE = 12, &
       MU_DOWN_FILE = 13, &
       N_SQD_FILE = 14, &
       V_UP_FILE = 15, &
       C_CC_FILE = 16, &
       C_ENT_FILE = 17, &
       CS2_FILE = 18, &
       CONSTRAINTS_FILE = 19, &
       FFT_FILE = 20

  integer, parameter :: &
       M_FN_MULTI_POLY = 1, &
       M_FN_TWO_STREAM = 2, &
       M_FN_COSMOLOGICAL_POLY = 3, &
       M_FN_GENERIC_PAPER = 4
  
end module variables_mod
