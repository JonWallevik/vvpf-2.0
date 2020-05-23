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
! No official documentation exist for this version. However, algorithms and other   !
! features of this program are well documented for the older version 1.0, which is  !
! in the following Ph.D.-thesis (Chapters 7,8 and 9; Appendix A):                   !
! Wallevik, J. E. (2003); Rheology of Particle Suspensions - Fresh Concrete, Mortar !
! and Cement Paste with Various Types of Lignosulfonates (Ph.D.-thesis); Department !
! of Structural Engineering, The Norwegian University of Science and Technology,    !
! ISBN 82-471-5566-4, ISSN 0809-103X.                                               !
!    ** Download: http://www.diva-portal.org/ntnu/theses/abstract.xsql?dbid=319 **  !
! NOTE: The comments present in all the source codes of this software are also most !
! valuable as documentation.                                                        !
! --------------------------------------------------------------------------------- !
! The files that make up the software VVPF 2.0 are "main.f90", "param.f90",         !
! "motion.f90", "viscous.f90", "write2f.f90", "shear.f90" and "update.f90".         !
! See the top part of each corresponding source code about its purpose.             !
! --------------------------------------------------------------------------------- !
! File name: main.f90 (PROGRAM)                                                     !
! This is the center of the whole software, holding and passing information to and  !
! from the different subroutines. Some subroutines interact directly with each      !
! other without going through the channels defined by main.f90 (this applies mostly !
! for the subroutines in the files update.f90, shear.f90 and viscous.f90).          !
! --------------------------------------------------------------------------------- !
PROGRAM MAIN_ROUTINE

USE CONSTANTS_AND_PARAMETERS
USE SHEAR_VISCOSITY
USE RATE_OF_SHEAR
USE ROTATION
USE MATRIX
USE WRITE_INFORMATION
IMPLICIT NONE
! --------------------------------------------------------------------------------- !
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:) :: ARRAY

DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)   :: VECTOR,v_kp1_new

DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:) :: SR_history,      SR_history_tmp,&
                                               dSR_dt_history,  dSR_dt_history_tmp

DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)   :: time_history,time_history_tmp

DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)   :: &
                                   U3_o,    U_3_k,    U_3_kp1_new,    U_3_kp1,&
                                   Usb3_o,  Usb_3_k,  Usb_3_kp1_new,  Usb_3_kp1,&
                                   Usb4_o,  Usb_4_k,  Usb_4_kp1_new,  Usb_4_kp1

DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)   :: &
                                   SR_kp1,      dSR_dt_kp1,    ddSR_dtdt_kp1,&
                                   VELOCITY_k,  VELOCITY_kp1,  VELOCITY_kp1_new

DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)   :: v_k,v_kp1
! --------------------------------------------------------------------------------- !
DOUBLE PRECISION :: & 
              tol_U_3,    RMS_U_3,    U_3_norm,    U_3,    U3_o_scalar,    CDF3_kp1,&
              tol_Usb_3,  RMS_Usb_3,  Usb_3_norm,  Usb_3,  Usb3_o_scalar,  SBF3_kp1,&
              tol_Usb_4,  RMS_Usb_4,  Usb_4_norm,  Usb_4,  Usb4_o_scalar,  SBF4_kp1

DOUBLE PRECISION :: time,dr,dt,R_i,R_o,h,omega,tol,a,b,SR,dSR_dt,ddSR_dtdt,&
                    dt_Newton,dt_Plastic,tol_Newton,tol_Plastic,tol_RMS_active,&
                    tol_RMS,RMS_vel,vel_norm,small_zero,EPS,nn_min,nn_max,ZERO_TIME,&
                    REAL_TIME,dt_OUTPUT,rho,Lambda

INTEGER          :: i,k,problem,NR,N_Lambda,N_Lambda_MAX,N_min,N_max,&
                    N_dt,N_read,N_history,count_n,MAX_NUMBER_OF_ITERATIONS,&
                    NUMBER_OF_TIME_ITERATIONS,count_rms,count_max,&
                    k_OUTPUT_rms,k_OUTPUT

LOGICAL          :: CONVERGENCE,TIME_INDEPENDENCE,WARNING_SIGN,FALSE_CONVERGENCE,&
                    CALCULATE_TIME_DEPENDENT_PROBL,OUTER_CYLINDER_ROTATES
                               
CHARACTER        :: IGNORED_INPUT
! --------------------------------------------------------------------------------- !
! The value "1501" used below must always be the same as the value of "N_read".     !
! The corresponding value should also be changed in the file "motion.f90".          !
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
DOUBLE PRECISION,DIMENSION(1501) :: raw_time_data,raw_omega_data
COMMON /RAW_DATA/ raw_time_data,raw_omega_data
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
N_read         = 1501    ! => The length of raw data points to be read. If "N_read" !
                         ! is not equal to the length of the input files (which are !
                         ! two, "time_input.dat" and "omega_input.dat"), then it    !
                         ! must be corrected manually in the files "main.f90" and   !
                         ! "motion.f90". Search for the string "N_read" and give it !
                         ! the correct value (for example 2000, if each input file  !
                         ! has 2000 data points). If such correction is not made, a !
                         ! "run time error" might occur during the computation and/ !
                         ! or strange numerical outcome might result from this.     !
! --------------------------------------------------------------------------------- !
raw_time_data  = 0.0D0   ! => Initializing raw_time_data.                           !
raw_omega_data = 0.0D0   ! => Initializing raw_omega_data.                          !
! --------------------------------------------------------------------------------- !
small_zero     = 0.1D-7  ! => Used in relation to screen and file output. This      !
                         !    value does generally not have to be changed.          !
EPS            = 1.0D-15 ! => Used in relation to vel_norm.                         !
! --------------------------------------------------------------------------------- !
PRINT *, " "
PRINT *, " "
PRINT *, " __________________________________________________________ "
PRINT *, "                                                            "
PRINT *, "         Viscometric-ViscoPlastic-Flow, version 2.0         "
PRINT *, "                     http://www.vvpf.net                    "
PRINT *, "                                                            "
PRINT *, "        Copyright (C) 2002, 2006 Dr. Jon E. Wallevik        "
PRINT *, "                    jon.wallevik@vvpf.net                   "
PRINT *, " __________________________________________________________ "
PRINT *, "                                                            "
PRINT *, " This software is free software; you can redistribute it    "
PRINT *, " and/or modify it under the terms of the GNU General Public "
PRINT *, " License as published by the Free Software Foundation;      "
PRINT *, " either version 2 of the License, or (at your option) any   "
PRINT *, " later version. This software is distributed in the hope    "
PRINT *, " that it will be useful, but WITHOUT ANY WARRANTY; without  "
PRINT *, " even the implied warranty of MERCHANTABILITY or FITNESS    "
PRINT *, " FOR A PARTICULAR PURPOSE.                                  "
PRINT *, " See the GNU General Public License for more details.       "
PRINT *, " __________________________________________________________ "
PRINT *, "                                                            "
WRITE(*,"(A)",ADVANCE="NO") "                PRESS 'ENTER' TO CONTINUE"
PRINT *, " "
READ (*,"(A)") IGNORED_INPUT
PRINT *, " "
! --------------------------------------------------------------------------------- !
! Retrieving variables like R_i, R_o, h, dr, dz, dt, tol, tol_RMS and so forth:
CALL CONSTANTS(rho,CALCULATE_TIME_DEPENDENT_PROBL,ZERO_TIME,tol_Newton,tol_Plastic,&
               tol_U_3,tol_Usb_3,tol_Usb_4,tol_RMS,dt_Plastic,dt_Newton,count_max,&
               dr,R_i,R_o,h,k_OUTPUT_rms,dt_OUTPUT,N_history,OUTER_CYLINDER_ROTATES)
! --------------------------------------------------------------------------------- !
! Retrieving the coagulated state at the time t = 0 sec for particle size number 3 
! (i.e. U3_o_scalar):
CALL COAG_DISP_3(0.0D0,U3_o_scalar,1.0D0,1.0D0,1.0D0,1.0D0,1.0D0,CDF3_kp1)

IF ((U3_o_scalar < 0.0D0).OR.(U3_o_scalar > 1.0D0)) THEN
  PRINT *, " ------------------------------------------------------------------- "
  PRINT *, "  ERROR: U3_o must be in the domain [0,1]; U3_o =                    "
  PRINT *,    U3_o_scalar
  PRINT *, "  Make correction in the file viscous.f90 (SUBROUTINE COAG_DISP_3)!  "
  PRINT *, " ------------------------------------------------------------------- "
  STOP
END IF
! --------------------------------------------------------------------------------- !
! For the structural-breakdown (sb = structural-breakdown; particle size number 3):
! Retrieving the linked state at the time t = 0 sec (i.e. Usb3_o_scalar):
CALL STRUCT_BREAK_3(0.0D0,Usb3_o_scalar,1.0D0,1.0D0,1.0D0,1.0D0,SBF3_kp1)

IF ((Usb3_o_scalar < 0.0D0).OR.(Usb3_o_scalar > 1.0D0)) THEN
  PRINT *, " ---------------------------------------------------------------------- "
  PRINT *, "  ERROR: Usb3_o must be in the domain [0,1]; Usb3_o =                   "
  PRINT *,    Usb3_o_scalar
  PRINT *, "  Make correction in the file viscous.f90 (SUBROUTINE STRUCT_BREAK_3)!  "
  PRINT *, " ---------------------------------------------------------------------- "
  STOP
END IF
! --------------------------------------------------------------------------------- !
! For the structural-breakdown (sb = structural-breakdown; particle size number 4):
! Retrieving the linked state at the time t = 0 sec (i.e. Usb4_o_scalar):
CALL STRUCT_BREAK_4(0.0D0,Usb4_o_scalar,1.0D0,1.0D0,1.0D0,1.0D0,SBF4_kp1)

IF ((Usb4_o_scalar < 0.0D0).OR.(Usb4_o_scalar > 1.0D0)) THEN
  PRINT *, " ---------------------------------------------------------------------- "
  PRINT *, "  ERROR: Usb4_o must be in the domain [0,1]; Usb4_o =                   "
  PRINT *,    Usb4_o_scalar
  PRINT *, "  Make correction in the file viscous.f90 (SUBROUTINE STRUCT_BREAK_4)!  "
  PRINT *, " ---------------------------------------------------------------------- "
  STOP
END IF
! --------------------------------------------------------------------------------- !
! Reading the time values (in seconds) from file:
PRINT *, " Reading time values from the file 'time_input.dat'..."
OPEN(UNIT=8,file="time_input.dat",status="old",action="read",&
     position="rewind",iostat=problem)
IF (problem /= 0) THEN
  PRINT *, " ------------------------------------------------- "
  PRINT *, "  ERROR: Could not open the file: time_input.dat!  "
  PRINT *, "  Make sure that this file exists!                 "
  PRINT *, " ------------------------------------------------- "
  STOP
ELSE
  DO i = 1,N_read
    READ (UNIT=8,FMT=*) raw_time_data(i)
  END DO
END IF

CLOSE (UNIT=8)
PRINT *, " ...done! "
! --------------------------------------------------------------------------------- !
! Reading the angular velocity (in rad/s) from file:
PRINT *, " Reading omega values from the file 'omega_input.dat'..."
OPEN(UNIT=8,file="omega_input.dat",status="old",action="read",&
     position="rewind",iostat=problem)
IF (problem /= 0) THEN
  PRINT *, " -------------------------------------------------- "
  PRINT *, "  ERROR: Could not open the file: omega_input.dat!  "
  PRINT *, "  Make sure that this file exists!                  "
  PRINT *, " -------------------------------------------------- "
  STOP
ELSE
  DO i = 1,N_read
    READ (UNIT=8,FMT=*) raw_omega_data(i)
  END DO
END IF

CLOSE (UNIT=8)
PRINT *, " ...done!"
! --------------------------------------------------------------------------------- !
! Calculating maximum and minimum values, and also determine their locations:
nn_min = raw_time_data(1)
nn_max = raw_time_data(1)
N_min  = 1
N_max  = 1
DO i = 2,N_read
  IF (raw_time_data(i) < nn_min) THEN
    nn_min = raw_time_data(i)
    N_min  = i
  ELSE IF (raw_time_data(i) > nn_min) THEN
    nn_max = raw_time_data(i)
    N_max  = i
  END IF
END DO
! ================================================================================= !
! ============================= Begin info to screen ============================== !
! ================================================================================= !
PRINT *, " __________________________________________________________ "
PRINT *, "                                                            "
PRINT *, "  The number of raw time-omega data points that was read is:"
PRINT *, "  N_read = ", N_read
PRINT *, "  Make sure that this value correspond to the length of     "
PRINT *, "  'time_input.dat' and 'omega_input.dat'!                   "
PRINT *, "  To change this value, edit the files 'main.f90' and       "
PRINT *, "  'motion.f90'. Note that the first entry of the file       "
PRINT *, "  'omega_input.dat' will be used as the angular velocity    "
PRINT *, "  for the time-independent calculation.                     "
PRINT *, "                                                            "

IF (OUTER_CYLINDER_ROTATES) THEN
  PRINT *, "  It is the outer cylinder that rotates! (see param.f90)  "
  PRINT *, "  The first entry in the file 'omega_input.dat' is...     "
  PRINT *, "  omega_[t=0](R_o) = ", raw_omega_data(1), " rad/s"
  PRINT *, "  or equally...                                           "
  PRINT *, "      f_[t=0](R_o) = ", raw_omega_data(1)/(2.0D0*ACOS(-1.0D0)), " rps"
ELSE IF (.NOT.OUTER_CYLINDER_ROTATES) THEN
  PRINT *, "  It is the inner cylinder that rotates! (see param.f90)  "
  PRINT *, "  The first entry in the file 'omega_input.dat' is...     "
  PRINT *, "  omega_[t=0](R_i) = ", raw_omega_data(1), " rad/s"
  PRINT *, "  or equally...                                           "
  PRINT *, "      f_[t=0](R_i) = ", raw_omega_data(1)/(2.0D0*ACOS(-1.0D0)), " rps"
ELSE
  PRINT *, "  The variable OUTER_CYLINDER_ROTATES is not set!         "
  PRINT *, "  Execution terminated!                                   "
  STOP
END IF

PRINT *, " "
PRINT *, "  Minimum time is ", nn_min," seconds"
PRINT *, "  ... and is at i = ",N_min
PRINT *, "  Maximum time is ", nn_max," seconds"
PRINT *, "  ... and is at i = ",N_max
PRINT *, " __________________________________________________________ "
PRINT *, "                                                            "
PRINT *, "              ___________________________                   "
WRITE(*,"(A)",ADVANCE="NO") "                PRESS 'ENTER' TO CONTINUE"
PRINT *, " "
READ (*,"(A)") IGNORED_INPUT
PRINT *, " "

IF (.NOT.OUTER_CYLINDER_ROTATES) THEN
  PRINT *, " __________________________________________________________ "
  PRINT *, "                                                            "
  PRINT *, " ATTENTION: Since it is the inner cylinder that is set to   "
  PRINT *, " rotate and the outer cylinder is set stationary, the       "
  PRINT *, " plugged state will now be in the same area as zero         "
  PRINT *, " velocity. This will make the system of algebraic equations "
  PRINT *, " more unstable, which means more successive substitution    "
  PRINT *, " iterations to compensate and hence LONGER COMPUTATION TIME."
  PRINT *, " The subroutines in this software are intelligent in the    "
  PRINT *, " sense that the next time step k + 2 is never calculated    "
  PRINT *, " until the program is satisfied with the accuracy of the    "
  PRINT *, " current time step k + 1. Hence, the longer calculation     "
  PRINT *, " time mentioned above, is just a consequence of the same    "
  PRINT *, " accuracy present in the calculations for either case of    "
  PRINT *, " rotating outer cylinder, or rotating inner cylinder.       "
  PRINT *, " __________________________________________________________ "
  PRINT *, "                                                            "
  WRITE(*,"(A)",ADVANCE="NO") "                PRESS 'ENTER' TO CONTINUE"
  PRINT *, " "
  READ (*,"(A)") IGNORED_INPUT
  PRINT *, " "
END IF
! --------------------------------------------------------------------------------- !
! PRINT *, raw_time_data
! PRINT *, raw_omega_data
! STOP
! --------------------------------------------------------------------------------- !
REAL_TIME = nn_max
! --------------------------------------------------------------------------------- !
IF (DABS(IDNINT(REAL_TIME/dt_Plastic)*dt_Plastic - REAL_TIME) > 1.0D-14) THEN
  PRINT *, "                                                        "
  PRINT *, " ------------------------------------------------------ "
  PRINT *, "  TERMINAL ERROR (main.f90):                            "
  PRINT *, "  The term 'IDNINT(REAL_TIME/dt_Plastic)*dt_Plastic'    "
  PRINT *, "  is not equal to 'REAL_TIME'!                          "
  PRINT *, " ------------------------------------------------------ "
  PRINT *, "  IDNINT(REAL_TIME/dt_Plastic)*dt_Plastic (seconds):    "
  PRINT *,    IDNINT(REAL_TIME/dt_Plastic)*dt_Plastic
  PRINT *, " ------------------------------------------------------ "
  PRINT *, "  REAL_TIME (seconds):                                  "
  PRINT *,    REAL_TIME
  PRINT *, " ------------------------------------------------------ "
  PRINT *, "  Select another time step for 'dt_Plastic'. This is    "
  PRINT *, "  done in the file 'param.f90'!                         "
  PRINT *, " ------------------------------------------------------ "
  PRINT *, "                                                        "
  STOP
END IF
! --------------------------------------------------------------------------------- !
NR  = IDNINT((R_o - R_i)/dr) + 1  ! (0.101 - 0.085)/0.0005 + 1 = 33 grid points     !
R_o = dr*DBLE(NR-1) + R_i         !  0.0005*(33-1) + 0.085 = 0.101                  !
! --------------------------------------------------------------------------------- !
! The variable "omega" is the angular velocity [rad/s] (here, at t = 0.0 sec).
CALL ANGULAR_VELOCITY(0.0D0,dt_Plastic,omega)
! --------------------------------------------------------------------------------- !
14 FORMAT(7X,"NR = ",(I3,2X),"; dr = ",F7.5," m;  rho = ",F6.1," kg/m^3")
15 FORMAT(7X,"R_i = ",F6.4," m  ; R_o = ",F6.4," m  ; h = ",F6.4," m")
16 FORMAT(7X,"dt_Plastic = ",E9.3," s")
17 FORMAT(7X,"U3_o = ",F6.4,"; Usb3_o = ",F6.4,"; Usb4_o = ",F6.4)
18 FORMAT(7X,"Regression time = ",F7.5," s")
PRINT '(5X,A40)',"Geometric values, time values and etc.:"
PRINT 14, NR,dr,rho
PRINT 15, R_i,R_o,h
PRINT 16, dt_Plastic
PRINT 17, U3_o_scalar,Usb3_o_scalar,Usb4_o_scalar
PRINT 18, dt_Plastic*N_history
PRINT *, "    "
! --------------------------------------------------------------------------------- !
CALL WARNING_FOR_WRITING(NR)
! --------------------------------------------------------------------------------- !
MAX_NUMBER_OF_ITERATIONS  = IDNINT(ZERO_TIME/dt_Plastic)   
NUMBER_OF_TIME_ITERATIONS = IDNINT(REAL_TIME/dt_Plastic) 

PRINT *, "      -------------------------------------------  "
PRINT "( '         MAX_NUMBER_OF_ITERATIONS: ', I10 )        ",&
                   MAX_NUMBER_OF_ITERATIONS
PRINT *, "      -------------------------------------------  "
PRINT *, "                                                   "

IF (CALCULATE_TIME_DEPENDENT_PROBL) THEN
  PRINT *, "      -------------------------------------------  "
  PRINT "( '         NUMBER_OF_TIME_ITERATIONS: ', I10 )       ",&
                     NUMBER_OF_TIME_ITERATIONS                
  PRINT *, "      -------------------------------------------  "
  PRINT *, "                                                   "
END IF

PRINT *, "              ___________________________                   "
WRITE(*,"(A)",ADVANCE="NO") "                PRESS 'ENTER' TO CONTINUE"
PRINT *, " "
READ (*,"(A)") IGNORED_INPUT
PRINT *, " "
! ================================================================================= !
! ============================== End info to screen =============================== !
! ================================================================================= !
! Creating log file and making the first RMS entry:
OPEN(UNIT=8,file="log.dat",status="replace",action="write",&
     position="rewind",iostat=problem)
IF (problem /= 0) THEN
  PRINT *, " Could not create the file: log.dat! "
  STOP
ELSE
  WRITE (UNIT=8,FMT=*) 0,0.0D0
END IF

CLOSE (UNIT=8)
! --------------------------------------------------------------------------------- !
ALLOCATE(ARRAY(NR-2,NR-2),VECTOR(NR-2),v_kp1_new(NR-2),stat=problem)
IF (problem /= 0) THEN 
  PRINT *, " MAIN_ROUTINE: The program could not allocate space!  "
  PRINT *, " Error code 1 in 'main.f90' and execution terminated! "
  STOP
END IF

ALLOCATE(SR_history(N_history,NR),     SR_history_tmp(N_history,NR),&
         dSR_dt_history(N_history,NR), dSR_dt_history_tmp(N_history,NR),&
         stat=problem)
IF (problem /= 0) THEN 
  PRINT *, " MAIN_ROUTINE: The program could not allocate space!  "
  PRINT *, " Error code 2 in 'main.f90' and execution terminated! "
  STOP
END IF

ALLOCATE(time_history(N_history),time_history_tmp(N_history),stat=problem)
IF (problem /= 0) THEN 
  PRINT *, " MAIN_ROUTINE: The program could not allocate space!  "
  PRINT *, " Error code 3 in 'main.f90' and execution terminated! "
  STOP
END IF

ALLOCATE(U3_o(NR),    U_3_k(NR),    U_3_kp1(NR),    U_3_kp1_new(NR),&
         Usb3_o(NR),  Usb_3_k(NR),  Usb_3_kp1(NR),  Usb_3_kp1_new(NR),&
         Usb4_o(NR),  Usb_4_k(NR),  Usb_4_kp1(NR),  Usb_4_kp1_new(NR),&
         stat=problem)
IF (problem /= 0) THEN 
  PRINT *, " MAIN_ROUTINE: The program could not allocate space!  "
  PRINT *, " Error code 4 in 'main.f90' and execution terminated! "
  STOP
END IF

ALLOCATE(SR_kp1(NR),      dSR_dt_kp1(NR),    ddSR_dtdt_kp1(NR),&
         VELOCITY_k(NR),  VELOCITY_kp1(NR),  VELOCITY_kp1_new(NR),&
         v_k(NR),         v_kp1(NR),         stat=problem)
IF (problem /= 0) THEN 
  PRINT *, " MAIN_ROUTINE: The program could not allocate space!  "
  PRINT *, " Error code 5 in 'main.f90' and execution terminated! "
  STOP
END IF
! --------------------------------------------------------------------------------- !
! U3_o is the coagulated state at time t = 0 sec, for the size 3 cement particles.  !
! If the coagulated state at time t = 0 sec is known as a function of radius r,     !
! one can assign different values to U3_o at different grid points. This means that !
! the value U3_o_scalar is overwritten.                                             !
! Example: U3_o(i) = DSQRT(DBLE(i)*U3_o_scalar/(DBLE(NR)))                          !
! --------------------------------------------------------------------------------- !
DO i = 1,NR
  U3_o(i) = U3_o_scalar
END DO
! --------------------------------------------------------------------------------- !
! Usb3_o is the linked state at time t = 0 sec (a "structural-breakdown" term), for !
! the size 3 cement particles:                                                      !
! Example: Usb3_o(i) = DSQRT(DBLE(i)*Usb3_o_scalar/(DBLE(NR)))                      !
! --------------------------------------------------------------------------------- !
DO i = 1,NR
  Usb3_o(i) = Usb3_o_scalar
END DO
! --------------------------------------------------------------------------------- !
! Usb4_o is the linked state at time t = 0 sec (a "structural-breakdown" term), for !
! the size 4 cement particles:                                                      !
! Example: Usb4_o(i) = DSQRT(DBLE(i)*Usb4_o_scalar/(DBLE(NR)))                      !
! --------------------------------------------------------------------------------- !
DO i = 1,NR
  Usb4_o(i) = Usb4_o_scalar
END DO
! --------------------------------------------------------------------------------- !
! Initialization: 
U_3_k              = U3_o
U_3_kp1            = U3_o
U_3_kp1_new        = U3_o

Usb_3_k            = Usb3_o
Usb_3_kp1          = Usb3_o
Usb_3_kp1_new      = Usb3_o

Usb_4_k            = Usb4_o
Usb_4_kp1          = Usb4_o
Usb_4_kp1_new      = Usb4_o

v_k                = 0.0D0
v_kp1              = 0.0D0
v_kp1_new          = 0.0D0

VELOCITY_k         = 0.0D0
VELOCITY_kp1       = 0.0D0
VELOCITY_kp1_new   = 0.0D0

SR_kp1             = 0.0D0

SR_history         = 0.0D0
SR_history_tmp     = 0.0D0
dSR_dt_kp1         = 0.0D0

dSR_dt_history     = 0.0D0
dSR_dt_history_tmp = 0.0D0
ddSR_dtdt_kp1      = 0.0D0

time_history       = 0.0D0
time_history_tmp   = 0.0D0

ARRAY              = 0.0D0
VECTOR             = 0.0D0

WARNING_SIGN       = .FALSE. 
FALSE_CONVERGENCE  = .FALSE.
! --------------------------------------------------------------------------------- !
! The variable "omega" is the angular velocity [rad/s] (here, at t = 0.0 sec).      !
CALL ANGULAR_VELOCITY(0.0D0,dt_Plastic,omega)
! --------------------------------------------------------------------------------- !
! Initialization of the Dirichlet boundary condition at t = 0.0 sec.                !
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
IF (OUTER_CYLINDER_ROTATES) THEN
! If it is the outer cylinder (R_o) that rotates:
  VELOCITY_k(1)  = 0.0D0
  VELOCITY_k(NR) = R_o*omega
ELSE
! If it is the inner cylinder (R_i) that rotates:
  VELOCITY_k(1)  = R_i*omega
  VELOCITY_k(NR) = 0.0D0
END IF
! --------------------------------------------------------------------------------- !
! Linear approximation for the velocity, to speed up the calculations:
a = VELOCITY_k(1)
b = VELOCITY_k(NR)
DO i = 2,NR-1
  VELOCITY_k(i) = a - (a - b)*DBLE(i-1)/DBLE(NR-1)
END DO
! CHECK OUT IF VELOCITY_k IS OK:
! CALL WRITE2FILE_VEL_DEBUG(VELOCITY_k)
! STOP
! ================================================================================= !
! ============================ Begin of CONTINUATION ============================== !
! ================================================================================= !
N_Lambda_MAX = 1
CONTINUATION: DO N_Lambda = 0,N_Lambda_MAX
! Lambda => The Continuation Method (see Section 7.8 in the Ph.D.-thesis).
Lambda = DBLE(N_Lambda)/DBLE(N_Lambda_MAX)
PRINT *, "________________________________________________________"
PRINT *, "CONTINUATION:",Lambda

! --------------------------------------------------------------------------------- !
IF (N_Lambda == 0) THEN       ! "N_Lambda = 0" => NEWTONIAN FLUID ::                !
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
  dt  = dt_Newton             ! Time step for the Newtonian case.                   !
  tol = tol_Newton            ! Tolerance for the successive substitution algorithm.!
  tol_RMS_active = tol_Newton ! Tolerance for time independence (see Equation 7.75).!
! --------------------------------------------------------------------------------- !
ELSE                          ! "N_Lambda > 0" => VISCOPLASTIC FLUID ::             !
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
  dt  = dt_Plastic            ! Time step for the viscoplastic case.                !
  tol = tol_Plastic           ! Tolerance for the successive substitution algorithm.!
  tol_RMS_active = tol_RMS    ! Tolerance for time independence (see Equation 7.75).!
! --------------------------------------------------------------------------------- !
END IF                                                                              !
! --------------------------------------------------------------------------------- !

TIME_INDEPENDENCE = .FALSE.
! Initializing time step k, for each CONTINUATION step:
k = 0
! ================================================================================= !
! =========================== Begin of the time loop ============================== !
! ================================================================================= !
! The iteration scheme used, when solving the time independent problem, is called   !
! "the pseudotransient method". It consists of gaining solution for the steady      !
! state problem, by time marching the transient governing equation (then with a     !
! constant boundary condition), until a steady state is reached. As such, the time  !
! (i.e. the time step k) plays the role of iteration parameter. For the full        !
! transient problem (shown further below), the same approach is used, however with  !
! continuous changes in the boundary condition. For further information, see        !
! Section 7.8, Page 167 in the pre-mentioned Ph.D.-thesis.                          !
! ================================================================================= !
! In Section 7.11.1 (Page 174) is a similar algorithm shown, to the one used below: !
! ================================================================================= !
ZERO_TIME_LOOP: DO WHILE (.NOT.TIME_INDEPENDENCE)
CONVERGENCE = .FALSE.
IF (ABS(MOD((k+1),k_OUTPUT_rms)) < small_zero) THEN
  PRINT *, "______________________________________________________"
  PRINT *, "                                                      "
  PRINT *, "  PSEUDO-TRANSIENT time step: k+1 = ",k+1
  PRINT *, "------------------------------------------------------"
END IF
! --------------------------------------------------------------------------------- !
! -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*- !
! The following routines are to update the Dirichlet boundary condition for the     !
! time step k+1. Most of these routines are redundant since the boundary conditions !
! are not changing here. However it can be a good practice to include them, while   !
! debugging the source code, if by some unfortunate accident some of the boundary   !
! values are overwritten.                                                           !
! --------------------------------------------------------------------------------- !
! The variable "omega" is the angular velocity [rad/s] (here, at t = 0.0 sec).      !
! CALL ANGULAR_VELOCITY(0.0D0,dt,omega)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
! IF (OUTER_CYLINDER_ROTATES) THEN
! ! If it is the outer cylinder (R_o) that rotates:
!   VELOCITY_k(1)  = 0.0D0
!   VELOCITY_k(NR) = R_o*omega
! ELSE
! ! If it is the inner cylinder (R_i) that rotates:
!   VELOCITY_k(1)  = R_i*omega
!   VELOCITY_k(NR) = 0.0D0
! END IF
! -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*- !
! --------------------------------------------------------------------------------- !
! A guess for the time step k+1:
VELOCITY_kp1     = VELOCITY_k
! -----------------------------
! Updating boundary condition for "VELOCITY_kp1_new":
VELOCITY_kp1_new = VELOCITY_kp1
! --------------------------------------------------------------------------------- !
count_n = 0
! ================================================================================= !
! ====================== BEGIN OF SUCCESSIVE SUBSTITUTION ========================= !
! ================================================================================= !
! The iteration loop here is because of the non-linearity of the governing equation !
! (i.e. of Equation 2.18, Page 16 in the pre-mentioned Ph.D.-thesis) and to come    !
! around this problem, the successive substitution approach is used (see Section    !
! 7.8). Within this loop, the value of "k" and "k+1" are constants. More precisely, !
! the time = k*dt and time = (k+1)*dt are frozen here:                              !
! --------------------------------------------------------------------------------- !
! (a) VELOCITY_k (i.e. "v_k") is a constant =>                                      !
! --------------------------------------------------------------------------------- !
v_k = VELOCITY_k(:)
! --------------------------------------------------------------------------------- !
! (b) VELOCITY_kp1_new is updated through "MATRIX_UPDATE" and "MATRIX_SOLVER";      !
! (c) VELOCITY_kp1 (i.e. "v_kp1") <- VELOCITY_kp1_new;                              !
! ================================================================================= !
CONVERGE: DO WHILE (.NOT.CONVERGENCE)

count_n = count_n + 1
10 FORMAT(4X,"Successive substitution number = ",1(I4,1X))

IF (ABS(MOD((k+1),k_OUTPUT_rms)) < small_zero) THEN
  PRINT 10, count_n
END IF

!=======================================================================
!========================= BEGIN OF ITERATION ==========================
v_kp1 = VELOCITY_kp1(:)

CALL MATRIX_UPDATE(U_3_kp1,Usb_3_kp1,Usb_4_kp1,rho,k,dt,Lambda,&
                   dr,R_i,NR,v_k,v_kp1,ARRAY,VECTOR)
CALL MATRIX_SOLVER(ARRAY,VECTOR,v_kp1_new,NR-2)

VELOCITY_kp1_new(2:NR-1) = v_kp1_new
! ----------------------------------------------------------------------
! It should be noted that no PDE's for the coagulated- and linked states
! are needed here in the time independent part, since only their 
! (constant) initial values are required, c.f. ...
! U_3_k       = U3_o;   Usb_3_k       = Usb3_o;   Usb_4_k       = Usb4_o
! U_3_kp1     = U3_o;   Usb_3_kp1     = Usb3_o;   Usb_4_kp1     = Usb4_o
! U_3_kp1_new = U3_o;   Usb_3_kp1_new = Usb3_o;   Usb_4_kp1_new = Usb4_o
!========================== END OF ITERATION ===========================
!=======================================================================
! --------- PAUSE FOR DEBUGGING --------- 
! CALL WRITE2FILE_VEL_DEBUG(VELOCITY_kp1_new)
! WRITE(*,"(A)",ADVANCE="NO") " PRESS 'ENTER' TO CONTINUE "
! PRINT *, " "
! READ (*,"(A)") IGNORED_INPUT
! PRINT *, " "
!=======================================================================

! Settings for testing of convergence (or rather stability):
CONVERGENCE = .TRUE.
RMS_vel     = 0.0D0
vel_norm    = 1.0D0

! ...checking the difference at each grid point, one by one...
DO i = 2,NR-1
  vel_norm =  (VELOCITY_kp1_new(i) + VELOCITY_kp1(i))/2.0D0 + EPS
  RMS_vel  = ((VELOCITY_kp1_new(i) - VELOCITY_kp1(i))/vel_norm)**2.0D0
  IF (RMS_vel > tol) THEN
    CONVERGENCE = .FALSE.
  END IF
END DO

! Updating ..._kp1
VELOCITY_kp1 = VELOCITY_kp1_new

IF (count_n == count_max) THEN
  CONVERGENCE       = .TRUE.
  FALSE_CONVERGENCE = .TRUE.
  PRINT *, " WARNING: FALSE CONVERGENCE! TIME STEP k = ", k
  PRINT *, " Maximum amount of successive substitutions is = ", count_max
  PRINT *, " ----------------------------------------------------------- "
  PRINT *, " RECOMMENDATION: Kill this application and reduce the time   "
  PRINT *, " step by an order of magnitude: dt -> dt/10                  "
  PRINT *, " ----------------------------------------------------------- "
END IF

END DO CONVERGE 
! ================================================================================= !
! ======================= END OF SUCCESSIVE SUBSTITUTION ========================== !
! ================================================================================= !
! Checking for time independence:
count_rms = 0
RMS_vel   = 0.0D0
vel_norm  = 1.0D0

! ...checking the overall difference between kp1 and k...
DO i = 2,NR-1
  count_rms = count_rms + 1
  vel_norm  =  (VELOCITY_kp1(i) + VELOCITY_k(i))/2.0D0 + EPS
  RMS_vel   = ((VELOCITY_kp1(i) - VELOCITY_k(i))/vel_norm)**2.0D0 + RMS_vel
END DO

RMS_vel = DSQRT(RMS_vel/count_rms)

IF (ABS(MOD((k+1),k_OUTPUT_rms)) < small_zero) THEN
  PRINT *, "   RMS_vel   =", RMS_vel
  CALL WRITE2FILE_rms(k+1,RMS_vel)
END IF

IF (RMS_vel < tol_RMS_active) THEN
  TIME_INDEPENDENCE = .TRUE.
ELSE 
  TIME_INDEPENDENCE = .FALSE.
END IF

IF (k == MAX_NUMBER_OF_ITERATIONS) THEN
  TIME_INDEPENDENCE = .TRUE.
  WARNING_SIGN      = .TRUE.
END IF

! Updating velocity and k, for the next pseudotransient time step...
VELOCITY_k = VELOCITY_kp1_new
k = k + 1

END DO ZERO_TIME_LOOP
! ================================================================================= !
! ============================ End of the time loop =============================== !
! ================================================================================= !
END DO CONTINUATION
! ================================================================================= !
! ============================= End of CONTINUATION =============================== !
! ================================================================================= !
IF (WARNING_SIGN) THEN 
  PRINT *, " --------------------------------------------------------- "
  PRINT *, " WARNING: k = MAX_NUMBER_OF_ITERATIONS; See log.dat        "
  PRINT *, " RECOMMENDATIONS:                                          "
  PRINT *, " I)  Rerun this application with reduced time step.        "
  PRINT *, "     Try order of magnitude less: dt -> dt/10.             "
  PRINT *, " II) Increase ZERO_TIME in the file param.f90.             "
  PRINT *, " --------------------------------------------------------- "
END IF

IF (FALSE_CONVERGENCE) THEN
  PRINT *, " --------------------------------------------------------- "
  PRINT *, " WARNING: FALSE CONVERGENCE WAS ACHIEVED. Reduce the time  "
  PRINT *, " step by order of magnitude: dt -> dt/10 and then rerun    "
  PRINT *, " the application.                                          "
  PRINT *, " --------------------------------------------------------- "
END IF

PRINT *, " PSEUDO-TRANSIENT CALCULATION FINISHED! "
PRINT *, " -------------------------------------- "
! --------------------------------------------------------------------------------- !
IF (.NOT.CALCULATE_TIME_DEPENDENT_PROBL) THEN
  PRINT *, " Number of grid points (not including Dirichlet " 
  PRINT *, " boundary points) = ", count_rms
  PRINT *, "                                                "
  PRINT *, " NO TIME DEPENDENT CALCULATION IS DONE SINCE    "
  PRINT *, " CALCULATE_TIME_DEPENDENT_PROBL = .FALSE.       "
  PRINT *, "                                                "
  PRINT *, " WRITING INFORMATION TO FILE, STAND BY...       "
  CALL WRITE2FILE_ZERO(U3_o,Usb3_o,Usb4_o,VELOCITY_k,NR,0,dt,Lambda,R_i,dr,h,omega)
  PRINT *, " ...DONE! "
  PRINT *, " EXECUTION FINISHED! "
  STOP
END IF
! ================================================================================= !
! ================================================================================= !
! ================================================================================= !
dt  = dt_Plastic
tol = tol_Plastic
! --------------------------------------------------------------------------------- !
! The flag "Lambda = -1" is used to inform all relevant subroutines that the time   !
! dependent calculation has begun, with time a dependent shear viscosity ETA =      !
! ETA(SR,t,...). More precisely, the shear viscosity is no longer frozen to its     !
! initial value ETA(SR,0,...) (which of course applies only at time t = 0.0 sec).   !
Lambda = - 1.0D0
! --------------------------------------------------------------------------------- !
! Note that at the moment, then VELOCITY_k = VELOCITY_kp1_new. 
PRINT *, "                                                "
PRINT *, " WRITING 't=0' INFORMATION TO FILE, STAND BY... "
CALL WRITE2FILE_ZERO(U3_o,Usb3_o,Usb4_o,VELOCITY_k,NR,0,dt,Lambda,R_i,dr,h,omega)
PRINT *, " ...DONE! "
PRINT *, "          "
PRINT *, " BEGINNING WITH THE TIME DEPENDENT PROBLEM... "
! --------------------------------------------------------------------------------- !
! Initializing the SR_history and time_history:
CALL ROS_PROFILE(VELOCITY_kp1_new,NR,R_i,dr,SR_kp1)
DO k = 1,N_history
  time_history(k)     = 0.0D0
  SR_history(k,:)     = SR_kp1(:)
  dSR_dt_history(k,:) = 0.0D0
END DO
! --------------------------------------------------------------------------------- !
! Initializing the coagulated- and linked state again, just in case it was 
! overwritten in the above. These routines are really redundant for stable versions
! of VVPF and should generally be commented:
! -------------------------
! U_3_k         = U3_o
! U_3_kp1       = U3_o
! U_3_kp1_new   = U3_o
! -------------------------
! Usb_3_k       = Usb3_o
! Usb_3_kp1     = Usb3_o
! Usb_3_kp1_new = Usb3_o
! -------------------------
! Usb_4_k       = Usb4_o
! Usb_4_kp1     = Usb4_o
! Usb_4_kp1_new = Usb4_o
! --------------------------------------------------------------------------------- !
N_dt = (NUMBER_OF_TIME_ITERATIONS - 1)
! ================================================================================= !
! ======================= Begin of the time loop (TIME) =========================== !
! ================================================================================= !
! In Section 7.11.2 (Page 175) is a similar algorithm shown, to the one used below: !
! ================================================================================= !
TIME_LOOP: DO k = 0,N_dt
! --------------------------------------------------------------------------------- !
! The software has already the solution for the time t = k*dt. Now, it is starting  !
! to calculate the solution for the time t=(k+1)*dt, i.e. the new velocity is:      !
! 'CALL ANGULAR_VELOCITY(DBLE(k)+1.0D0,dt,omega)'. For example, "k = 0", means      !
! "t = 0 sec", and that solution is already obtained from the time-independent      !
! calculation done above (that solution is at the moment present in "VELOCITY_k").  !
! Hence, for this particular example, the software is starting to calculate the     !
! solution for the time "t = 1*dt".                                                 !
! --------------------------------------------------------------------------------- !
CONVERGENCE = .FALSE.

IF (ABS(MOD((k+1),k_OUTPUT_rms)) < small_zero) THEN
  PRINT *, "______________________________________________________________________"
  PRINT *, "______________________________________________________________________"
  PRINT *, "                                                                      "
  PRINT *, "  Calculating now for the time step k+1 = ",k+1," -> ->               "
  PRINT *, "  time = (k+1)dt   = ",(DBLE(k)+1.0D0)*dt," SEC"
  PRINT *, "----------------------------------------------------------------------"
END IF
! --------------------------------------------------------------------------------- !
! A guess for the time step k+1:
VELOCITY_kp1 = VELOCITY_k
U_3_kp1      = U_3_k
Usb_3_kp1    = Usb_3_k
Usb_4_kp1    = Usb_4_k
! --------------------------------------------------------------------------------- !
! The variable "omega" is the angular velocity [rad/s] (here, at t = (k+1)*dt sec). !
CALL ANGULAR_VELOCITY(DBLE(k)+1.0D0,dt,omega)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
! Updating Dirichlet boundary condition at k+1:                                     !
IF (OUTER_CYLINDER_ROTATES) THEN
! If it is the outer cylinder (R_o) that rotates:
  VELOCITY_kp1(1)  = 0.0D0
  VELOCITY_kp1(NR) = R_o*omega
ELSE
! If it is the inner cylinder (R_i) that rotates:
  VELOCITY_kp1(1)  = R_i*omega
  VELOCITY_kp1(NR) = 0.0D0
END IF
! --------------------------------------------------------------------------------- !
! Updating boundary condition for "VELOCITY_kp1_new":
VELOCITY_kp1_new = VELOCITY_kp1
! --------------------------------------------------------------------------------- !
count_n = 0
! ================================================================================= !
! =================== BEGIN OF SUCCESSIVE SUBSTITUTION (TIME) ===================== !
! ================================================================================= !
! The iteration loop here is because of the non-linearity of the governing equation !
! (i.e. of Equation 2.18, Page 16 in the pre-mentioned Ph.D.-thesis) and to come    !
! around this problem, the successive substitution approach is used (see Section    !
! 7.8). Within this loop, the value of "k" and "k+1" are constants. More precisely, !
! the time = k*dt and time = (k+1)*dt are frozen here.                              !
! What is occurring here is listed with the points (a) to (l):                       !
! --------------------------------------------------------------------------------- !
! (a) VELOCITY_k (i.e. "v_k") is a constant =>                                      !
! --------------------------------------------------------------------------------- !
v_k = VELOCITY_k(:)
! --------------------------------------------------------------------------------- !
! (b) VELOCITY_kp1_new is updated through "MATRIX_UPDATE" and "MATRIX_SOLVER";      !
! (c) VELOCITY_kp1 (i.e. "v_kp1") <- VELOCITY_kp1_new;                              !
! --------------------------------------------------------------------------------- !
! (d) U_3_k = constant;                                                             !
! (e) U_3_kp1_new is updated through the successive substitution.                   !
! (f) U_3_kp1 = U_3_kp1_new;                                                        !
! --------------------------------------------------------------------------------- !
! (g) Usb_3_k = constant;                                                           !
! (h) Usb_3_kp1_new is updated through the successive substitution.                 !
! (i) Usb_3_kp1 = Usb_3_kp1_new;                                                    !
! --------------------------------------------------------------------------------- !
! (j) Usb_4_k = constant;                                                           !
! (k) Usb_4_kp1_new is updated through the successive substitution.                 !
! (l) Usb_4_kp1 = Usb_4_kp1_new;                                                    !
! ================================================================================= !
TIME_CONVERGE: DO WHILE (.NOT.CONVERGENCE)

count_n = count_n + 1
12 FORMAT(4X,"Successive substitution (time) number = ",1(I4,1X))

IF (ABS(MOD((k+1),k_OUTPUT_rms)) < small_zero) THEN
  PRINT 12, count_n
END IF

!=======================================================================
!=============== (TIME) == BEGIN OF ITERATION == (TIME) ================
v_kp1 = VELOCITY_kp1(:)

CALL MATRIX_UPDATE(U_3_kp1,Usb_3_kp1,Usb_4_kp1,rho,k,dt,Lambda,&
                   dr,R_i,NR,v_k,v_kp1,ARRAY,VECTOR)
CALL MATRIX_SOLVER(ARRAY,VECTOR,v_kp1_new,NR-2)

VELOCITY_kp1_new(2:NR-1) = v_kp1_new
!=======================================================================
! Here is the calculation of the different coagulated- and linked states
! U_3, Usb_3 and Usb_4.
! ----------------------------------------------------------------------
time_history(N_history) = DBLE(k+1)*dt
time                    = DBLE(k+1)*dt
! ----------------------------------------------------------------------
! Calculating the shear rate at k + 1 (i.e. at kp1):
CALL ROS_PROFILE(VELOCITY_kp1_new,NR,R_i,dr,SR_kp1)

! ...and updating the newest SR info (kp1_new) into the SR_history:
SR_history(N_history,:) = SR_kp1
! ----------------------------------------------------------------------
! Calculating the newest time derivative of shear rate (from the "_new"):
CALL dROS_dt(k,NR,N_history,time_history,SR_history,dSR_dt_kp1)

! ...and updating the newest dSR_dt info into the dSR_dt_history:
dSR_dt_history(N_history,:) = dSR_dt_kp1
! ----------------------------------------------------------------------
! Calculating the newest double time derivative of shear rate 
! (from the "_new"):
CALL ddROS_dtdt(k,NR,N_history,time_history,dSR_dt_history,ddSR_dtdt_kp1)
! ----------------------------------------------------------------------
DO i = 1,NR
  SR        = SR_kp1(i)
  dSR_dt    = dSR_dt_kp1(i)
  ddSR_dtdt = ddSR_dtdt_kp1(i)

  U_3       = U_3_kp1(i)
  Usb_3     = Usb_3_kp1(i)
  Usb_4     = Usb_4_kp1(i)

  CALL COAG_DISP_3(time,U3_o_scalar,U_3,Usb_3,SR,dSR_dt,ddSR_dtdt,CDF3_kp1)
  U_3_kp1_new(i) = U_3_k(i) + dt*CDF3_kp1

  CALL STRUCT_BREAK_3(time,Usb3_o_scalar,Usb_3,SR,dSR_dt,ddSR_dtdt,SBF3_kp1)
  Usb_3_kp1_new(i) = Usb_3_k(i) + dt*SBF3_kp1

  CALL STRUCT_BREAK_4(time,Usb4_o_scalar,Usb_4,SR,dSR_dt,ddSR_dtdt,SBF4_kp1)
  Usb_4_kp1_new(i) = Usb_4_k(i) + dt*SBF4_kp1

  ! --------------------------------------------------------------------------------!
  ! Because of the truncation error, the U values can accidentally become less than !
  ! zero. This is usually related to the problem when for example somewhere in the  !
  ! calculations, a very large number is multiplied with a very small number and    !
  ! the result is then subtracted from a medium sized number, resulting in that     !
  ! relevant information is lost on the way. This type of problem arises because a  !
  ! computer can only work with limited number of significant figures (see also the !
  ! footnote on Page 159 in the pre-mentioned Ph.D.-thesis). Since the software     !
  ! works with DOUBLE PRECISION variables, only small negative values initially     !
  ! results for U, when occurring. However, it can (and will) increase (i.e. become !
  ! more negative), causing substantial problem on the way, during the simulation.  !
  ! The following routine is to eliminate any negative values at its birth.         !
  ! --------------------------------------------------------------------------------!
  IF (U_3_kp1_new(i) < 0.0D0) THEN
    U_3_kp1_new(i) = 0.0D0
  END IF

  IF (Usb_3_kp1_new(i) < 0.0D0) THEN
    Usb_3_kp1_new(i) = 0.0D0
  END IF

  IF (Usb_4_kp1_new(i) < 0.0D0) THEN
    Usb_4_kp1_new(i) = 0.0D0
  END IF
END DO
!=============== (TIME) === END OF ITERATION === (TIME) ================
!=======================================================================
! Settings for testing of convergence (or rather stability):
CONVERGENCE = .TRUE.

RMS_vel     = 0.0D0
vel_norm    = 1.0D0

RMS_U_3     = 0.0D0
U_3_norm    = 1.0D0

RMS_Usb_3   = 0.0D0
Usb_3_norm  = 1.0D0

RMS_Usb_4   = 0.0D0
Usb_4_norm  = 1.0D0

! ...checking the difference at each grid point, one by one...
DO i = 2,NR-1
  vel_norm   =  (VELOCITY_kp1_new(i) + VELOCITY_kp1(i))/2.0D0 + EPS
  RMS_vel    = ((VELOCITY_kp1_new(i) - VELOCITY_kp1(i))/vel_norm)**2.0D0

  U_3_norm   =  (U_3_kp1_new(i) + U_3_kp1(i))/2.0D0 + EPS
  RMS_U_3    = ((U_3_kp1_new(i) - U_3_kp1(i))/U_3_norm)**2.0D0

  Usb_3_norm =  (Usb_3_kp1_new(i) + Usb_3_kp1(i))/2.0D0 + EPS
  RMS_Usb_3  = ((Usb_3_kp1_new(i) - Usb_3_kp1(i))/Usb_3_norm)**2.0D0

  Usb_4_norm =  (Usb_4_kp1_new(i) + Usb_4_kp1(i))/2.0D0 + EPS
  RMS_Usb_4  = ((Usb_4_kp1_new(i) - Usb_4_kp1(i))/Usb_4_norm)**2.0D0

  IF (RMS_vel > tol) THEN
    CONVERGENCE = .FALSE.
  ELSE IF ((RMS_U_3 > tol_U_3).OR.(RMS_Usb_3 > tol_Usb_3)) THEN
    CONVERGENCE = .FALSE.
  ELSE IF (RMS_Usb_4 > tol_Usb_4) THEN
    CONVERGENCE = .FALSE.
  END IF
END DO

! Updating ..._kp1
VELOCITY_kp1 = VELOCITY_kp1_new
U_3_kp1      = U_3_kp1_new
Usb_3_kp1    = Usb_3_kp1_new
Usb_4_kp1    = Usb_4_kp1_new

IF (count_n == count_max) THEN
  CONVERGENCE       = .TRUE.
  FALSE_CONVERGENCE = .TRUE.
  PRINT *, " WARNING: FALSE CONVERGENCE! TIME STEP k= ", k
  PRINT *, " Maximum amount of successive substitutions is = ", count_max
  PRINT *, " ----------------------------------------------------------- "
  PRINT *, " RECOMMENDATION: Kill this application and reduce the time   "
  PRINT *, " step by an order of magnitude: dt -> dt/10                  "
  PRINT *, " ----------------------------------------------------------- "
END IF

END DO TIME_CONVERGE
! ================================================================================= !
! ==================== END OF SUCCESSIVE SUBSTITUTION (TIME) ====================== !
! ================================================================================= !
! Updating the velocity VELOCITY_k, U_3, Usb_3, Usb_4 and SR_history: kp1 -> k, and
! the next time step starts...
VELOCITY_k = VELOCITY_kp1_new
U_3_k      = U_3_kp1_new
Usb_3_k    = Usb_3_kp1_new
Usb_4_k    = Usb_4_kp1_new
! --------------------------------------------------------------------------------- !
! Moving the stack and making room for new information from the coming step k+2,    !
! i.e. the freed stack will not be updated until in the next iteration (kp1 -> k):  !
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
time_history_tmp                = time_history
time_history(1:N_history-1)     = time_history_tmp(2:N_history)
SR_history_tmp                  = SR_history
SR_history(1:N_history-1,:)     = SR_history_tmp(2:N_history,:)
dSR_dt_history_tmp              = dSR_dt_history
dSR_dt_history(1:N_history-1,:) = dSR_dt_history_tmp(2:N_history,:)
! --------------------------------------------------------------------------------- !
k_OUTPUT = IDNINT(dt_OUTPUT/dt)
IF (ABS(MOD((k+1),k_OUTPUT)) < small_zero) THEN
  PRINT *, " =================================================================== "
  PRINT *, " WRITING DATA TO FILE AT TIME ", DBLE(k+1)*dt,"SECONDS"
  CALL WRITE2FILE(U_3_kp1,Usb_3_kp1,Usb_4_kp1,VELOCITY_kp1,dSR_dt_kp1,&
                  ddSR_dtdt_kp1,NR,k+1,dt,Lambda,R_i,dr,h,omega)
  PRINT *, " =================================================================== "
END IF
! --------------------------------------------------------------------------------- !
END DO TIME_LOOP
! ================================================================================= !
! ======================== End of the time loop (TIME) ============================ !
! ================================================================================= !
IF (WARNING_SIGN) THEN
  PRINT *, " --------------------------------------------------------- "
  PRINT *, " WARNING: For the time independent case, then:             "
  PRINT *, " k = MAX_NUMBER_OF_ITERATIONS; See log.dat                 "
  PRINT *, " --------------------------------------------------------- "
END IF

IF (FALSE_CONVERGENCE) THEN
  PRINT *, " --------------------------------------------------------- "
  PRINT *, " WARNING: For the time independent case, then:             "
  PRINT *, " FALSE CONVERGENCE WAS ACHIEVED. Reduce the time           "
  PRINT *, " step by order of magnitude: dt -> dt/10 and then rerun    "
  PRINT *, " the application.                                          "
  PRINT *, " --------------------------------------------------------- "
END IF

IF (.NOT.OUTER_CYLINDER_ROTATES) THEN
  PRINT *, " For this calculation, it is the inner cylinder that rotates. "
  PRINT *, " In this case, torque is here calculated as FROM the test     "
  PRINT *, " material ON the inner cylinder. Hence, the torque should be  "
  PRINT *, " negative (meaning that the test material is trying to slow   "
  PRINT *, " down the rotation of the inner cylinder).                    "
  PRINT *, " The torque applied FROM the inner cylinder ON the test       "
  PRINT *, " material, is simply gained by putting a minus sign in front  "
  PRINT *, " of the current calculated values, and would then result in   "
  PRINT *, " positive values. This operation is reserved for the user.    "
END IF

PRINT *, " Number of grid points (not including the Dirichlet "
PRINT *, " boundary points) =   ", count_rms
PRINT *, "                      "
PRINT *, " EXECUTION FINISHED!  "
! --------------------------------------------------------------------------------- !
END PROGRAM MAIN_ROUTINE
! --------------------------------------------------------------------------------- !
