! --------------------------------------------------------------------------------- !
!                                                                                   !
!     Copyright (C) 2002, 2006 Dr. Jon E. Wallevik (jon.wallevik@vvpf.net),         !
!     Software homepage: http://www.vvpf.net                                        !
!     The Norwegian University of Science and Technology (NTNU), 1998 - 2003.       !
!     Icelandic Building Research Institute, 2003 - 2004.                           !
!     The Norwegian University of Science and Technology (NTNU), 2004 - 2006.       !
!     Reykjavik University, 2006.                                                   !
!                                                                                   !
!     This file is part of Viscometric-ViscoPlastic-Flow, version 2.0.              !
!                                                                                   !
!     Viscometric-ViscoPlastic-Flow, is free software; you can redistribute it      !
!     and/or modify it under the terms of the GNU General Public License as         !
!     published by the Free Software Foundation; either version 2 of the            !
!     License, or (at your option) any later version.                               !
!                                                                                   !
!     Viscometric-ViscoPlastic-Flow, is distributed in the hope that it will be     !
!     useful, but WITHOUT ANY WARRANTY; without even the implied warranty of        !
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General      !
!     Public License for more details.                                              !
!                                                                                   !
!     You should have received a copy of the GNU General Public License             !
!     along with Viscometric-ViscoPlastic-Flow; if not, write to the Free           !
!     Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,                !
!     MA  02111-1307  USA                                                           !
!                                                                                   !
! --------------------------------------------------------------------------------- !
! File name: param.f90 (MODULE)                                                     !
! This code defines and sets all variables of relevance, like R_i, R_o, h, dr, dz,  !
! dt, tol, tol_RMS and so forth.                                                    !
! --------------------------------------------------------------------------------- !
MODULE CONSTANTS_AND_PARAMETERS
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: CONSTANTS
CONTAINS
! ================================================================================= !
SUBROUTINE CONSTANTS(rho,CALCULATE_TIME_DEPENDENT_PROBL,ZERO_TIME,tol_Newton,&
                     tol_Plastic,tol_U_3,tol_Usb_3,tol_Usb_4,tol_RMS,dt_Plastic,&
                     dt_Newton,count_max,dr,R_i,R_o,h,k_OUTPUT_rms,dt_OUTPUT,&
                     N_history,OUTER_CYLINDER_ROTATES)
! --------------------------------------------------------------------------------- !
DOUBLE PRECISION,INTENT(OUT) :: rho,ZERO_TIME,tol_Newton,tol_Plastic,&
                                tol_U_3,tol_Usb_3,tol_Usb_4,tol_RMS,&
                                dt_Plastic,dt_Newton,dr,R_i,R_o,h,dt_OUTPUT

LOGICAL,INTENT(OUT)          :: CALCULATE_TIME_DEPENDENT_PROBL,&
                                OUTER_CYLINDER_ROTATES
INTEGER,INTENT(OUT)          :: count_max,k_OUTPUT_rms,N_history
! --------------------------------------------------------------------------------- !
! Informing the main routine about if it is the parabolic- or the elliptic problem
! that is to be solved:
! If parabolic problem (time dependence), then:
CALCULATE_TIME_DEPENDENT_PROBL = .TRUE.
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
! If elliptic problem (time independence), then:
! CALCULATE_TIME_DEPENDENT_PROBL = .FALSE.
! --------------------------------------------------------------------------------- !
! Informing the main routine if it is the outer- or the inner cylinder that rotates:
! If it is the outer cylinder that rotates, then:
OUTER_CYLINDER_ROTATES = .TRUE.
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
! If it is the inner cylinder that rotates, then:
! OUTER_CYLINDER_ROTATES = .FALSE.
! --------------------------------------------------------------------------------- !
rho          = 2090.0D0 ! Density of the test material in kg/m^3.                   !
! --------------------------------------------------------------------------------- !
tol_Newton   = 0.5D-4   ! Used as condition for time-independence in the Newtonian  !
                        ! case. Also used as tolerance for the successive           !
                        ! substitution in the Newtonian case.                       !
tol_RMS      = 1.0D-20  ! Condition for time-independence in the viscoplastic case  !
                        ! (see also Equation 7.75).                                 !
tol_Plastic  = 1.0D-12  ! The successive substitution tolerance for both the time-  !
                        ! dependent and the time-independent case (Equation 7.73).  !
tol_U_3      = 1.0D-12  ! The successive substitution tolerance for U_3 (only used  !
                        ! for the time-dependent case).                             !
tol_Usb_3    = 1.0D-12  ! The successive substitution tolerance for Usb_3 (only     !
                        ! used for the time-dependent case).                        !
tol_Usb_4    = 1.0D-12  ! The successive substitution tolerance for Usb_4 (only     !
                        ! used for the time-dependent case).                        !
! --------------------------------------------------------------------------------- !
dt_Newton    = 0.1D-2   ! Only used for the Newtonian case.                         !
dt_Plastic   = 0.5D-4   ! Minimum value for dt_Plastic is 0.1D-6. If a lower value  !
                        ! is preferred, then adjustments has to be made to:         !
                        ! NUMBER_OF_TIME_ITERATIONS = IDNINT(REAL_TIME/dt_Plastic)  !
                        ! in "main.f90". Possibly, one can then use instead:        !
                        ! NUMBER_OF_TIME_ITERATIONS = DINT(REAL_TIME/dt_Plastic)    !
! --------------------------------------------------------------------------------- !
count_max    = 500      ! Maximum number of successive substitution iterations, for !
                        ! each time step k.                                         !
! --------------------------------------------------------------------------------- !
! ZERO_TIME (Section 7.11): Maximum pseudotransient iteration time in seconds, when !
! solving the time-independent problem. It will either be this time (i.e. ZERO_TIME)!
! or the condition "tol_RMS" (i.e. "tol_RMS_active"; see main.f90) that will term-  !
! inate the time-independent (i.e. the elliptical) iteration.                       !
ZERO_TIME    = 5.0D0    ! Example: 5.0D0 => 5 seconds of steady state calculation.  !
! --------------------------------------------------------------------------------- !
k_OUTPUT_rms = 10000    ! Information output every k_OUTPUT_rms times, i.e. every   !
                        ! "dt*k_OUTPUT_rms" second to console. During the time-     !
                        ! independent calculation, k_OUTPUT_rms also activates the  !
                        ! the subroutine "WRITE2FILE_rms" to write the RMS data     !
                        ! into the file "log.dat".                                  !
! --------------------------------------------------------------------------------- !
dt_OUTPUT    = 0.1D0    ! vel, eta, shear stress and etc. written to files every    !
                        ! dt_OUTPUT second [sec].                                   !
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
! NOTE: If a strange or inconsistent results is gained in the files dSR_dt.dat and  !
! ddSR_dtdt.dat, put dt_OUTPUT = 0.01D0 or lower and reanalyze.                     !
! --------------------------------------------------------------------------------- !
! NOTE [IMPORTANT]: If the variable "ddSR_dtdt" is used in the material functions   !
! CDF3, SBF3 or SBF4, defined in the file "viscous.f90", then it is recommended     !
! that the following applies: N_history * dt_Plastic = 0.05 sec (or so).            !
! For example:                                                                      !
!   (1) N_history * dt_Plastic = 5000 * 1.0D-5 sec = 0.05 sec.                      !
!   (2) N_history * dt_Plastic =  500 * 1.0D-4 sec = 0.05 sec.                      !
!   (3) N_history * dt_Plastic =  100 * 5.0D-4 sec = 0.05 sec.                      !
! This time of 0.05 sec, represents the time period which regression is made when   !
! calculating the time derivative of the shear rate (dSR_dt; i.e. gamma double dot) !
! which is again used in calculating the double time derivative of the shear rate   !
! (ddSR_dtdt; i.e. gamma triple dot). Longer time than 0.05 will result in more     !
! incorrect result for ddSR_dtdt (and to somewhat lesser extent for dSR_dt), and    !
! shorter time will make numerical oscillations (errors) more pronounced.           !
! NOTE that the above is only a recommendation and can be incorrect in some cases.  !
! USE "ddSR_dtdt" WITH CAUTION!!!                                                   !
! - - - - - - - - - - - ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
N_history = 20          ! Length of the SR_history array, which is the registry of  !
                        ! current (k+1) and previous (..., k-4, k-3, k-2, k-1, k)   !
                        ! shear rate values. The SR_history array is used in the    !
                        ! calculation of dSR_dt (time derivative of the shear rate; !
                        ! i.e. gamma double dot). Note that the minimum value of    !
                        ! N_history is 3; No maximum value.                         !
! --------------------------------------------------------------------------------- !
dr  = 0.25D-3           ! =>  0.5 mm = Spacing between grid points in r-direction.  !
R_i = 0.085D0           ! =>  8.5 cm = Inner radius of viscometer.                  !
R_o = 0.101D0           ! => 10.1 cm = Outer radius of viscometer.                  !
h   = 0.116D0           ! => 11.6 cm = Height where torque is measured.             !
! --------------------------------------------------------------------------------- !
RETURN
END SUBROUTINE CONSTANTS
! ================================================================================= !
END MODULE CONSTANTS_AND_PARAMETERS
! --------------------------------------------------------------------------------- !
