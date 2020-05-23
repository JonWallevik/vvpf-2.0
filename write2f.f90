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
! File name: write2f.f90 (MODULE)                                                   !
! This file takes care of writing all computed data into the different files. It is !
! only the source main.f90 that makes such request.                                 !
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
! USE "ddSR_dtdt" WITH CAUTION!!!                                                   !
! NOTE that the above is only a recommendation and can be incorrect in some cases.  !
! --------------------------------------------------------------------------------- !
! ADDITIONAL NOTE: If a strange or inconsistent results are gained, put dt_OUTPUT = !
! 0.01D0 or lower, and reanalyze.                                                   !
! --------------------------------------------------------------------------------- !
MODULE WRITE_INFORMATION
  USE SHEAR_VISCOSITY  
  USE RATE_OF_SHEAR
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: WARNING_FOR_WRITING,  WRITE2FILE_VEL_DEBUG,  ARRAY_VECTOR_debug,&
            WRITE2FILE_ZERO,      WRITE2FILE,            WRITE2FILE_rms
CONTAINS
! ================================================================================= !
SUBROUTINE WARNING_FOR_WRITING(NR)
! --------------------------------------------------------------------------------- !
INTEGER,INTENT(IN) :: NR
! --------------------------------------------------------------------------------- !
IF (NR > 300) THEN
  PRINT *, " ERROR: NR > 300 ( NR = ",NR,")"
  PRINT *, " FORMAT STATEMENT IN THE FILE 'write2f.f90' IS TO SHORT:     "
  PRINT *, " ERROR -> 10 FORMAT(1X,300(F7.4,1X))                         "
  PRINT *, " PLEASE MAKE THE NECESSARY ADJUSTMENT IN ALL THE SUBROUTINES "
  PRINT *, " OF THIS FILE. TERMINAL ERROR!                               "
  STOP
END IF
! --------------------------------------------------------------------------------- !
RETURN
END SUBROUTINE WARNING_FOR_WRITING
! ================================================================================= !
SUBROUTINE WRITE2FILE_VEL_DEBUG(VELOCITY)
! --------------------------------------------------------------------------------- !
DOUBLE PRECISION,DIMENSION(:),INTENT(IN) :: VELOCITY
INTEGER                                  :: problem
! --------------------------------------------------------------------------------- !
10 FORMAT(1X,300(F10.7,1X))
! --------------------------------------------------------------------------------- !
OPEN(UNIT=8,file="vel_debug.dat",status="replace",action="write",&
     position="rewind",iostat=problem)
IF (problem /= 0) THEN
  PRINT *, " TERMINAL ERROR (write2f.f90):             "
  PRINT *, " Could not create the file: vel_debug.dat! "
  STOP
ELSE
  WRITE (UNIT=8,FMT=10) VELOCITY(:)
END IF

CLOSE (UNIT=8)
! --------------------------------------------------------------------------------- !
RETURN
END SUBROUTINE WRITE2FILE_VEL_DEBUG
! ================================================================================= !
SUBROUTINE WRITE2FILE_rms(k,rms) 
! --------------------------------------------------------------------------------- !
DOUBLE PRECISION,INTENT(IN) :: rms
INTEGER,INTENT(IN)          :: k
INTEGER                     :: problem
! --------------------------------------------------------------------------------- !
OPEN(UNIT=8,file="log.dat",status="old",action="write",&
     position="append",iostat=problem)
IF (problem /= 0) THEN
  PRINT *, " TERMINAL ERROR (write2f.f90):                    "
  PRINT *, " Could not write into the existing file: log.dat! "
  STOP
ELSE
  WRITE (UNIT=8,FMT=*) k,rms
END IF

CLOSE (UNIT=8)
! --------------------------------------------------------------------------------- !
RETURN
END SUBROUTINE WRITE2FILE_rms
! ================================================================================= !
SUBROUTINE WRITE2FILE_ZERO(U_3_kp1,Usb_3_kp1,Usb_4_kp1,V,NR,kp1,&
                           dt,Lambda,R_i,dr,h,omega)
! --------------------------------------------------------------------------------- !
DOUBLE PRECISION,DIMENSION(:),INTENT(IN)  :: U_3_kp1,Usb_3_kp1,Usb_4_kp1,V
DOUBLE PRECISION,INTENT(IN)               :: dt,Lambda,R_i,dr,h,omega
INTEGER,INTENT(IN)                        :: NR,kp1

DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: SR,ETA,TORQUE_r,Shear_Stress
! --------------------------------------------------------------------------------- !
DOUBLE PRECISION :: U_3,Usb_3,Usb_4,TIME,ETA_tmp,SR_tmp
INTEGER          :: problem,i
LOGICAL          :: WRITE2FILE_ERROR
! --------------------------------------------------------------------------------- !
WRITE2FILE_ERROR = .FALSE.

TIME = DBLE(kp1)*dt

ALLOCATE(SR(NR),ETA(NR),TORQUE_r(NR),Shear_Stress(NR),stat=problem)
IF (problem /= 0) THEN 
  PRINT *, " WRITE2FILE_ZERO: The program could not allocate space!  "
  PRINT *, " Error code 1 in 'write2f.f90' and execution terminated! "
  STOP
END IF

SR           = 0.0D0
ETA          = 0.0D0
TORQUE_r     = 0.0D0
Shear_Stress = 0.0D0
! --------------------------------------------------------------------------------- !
! CALCULATING THE SHEAR RATE AT ALL POINTS (shear.f90): SR 
CALL ROS_PROFILE(V,NR,R_i,dr,SR)
! --------------------------------------------------------------------------------- !
! CALCULATING THE SHEAR VISCOSITY (i.e. apparent viscosity) AT ALL POINTS:
ETA_tmp = 0.0D0
DO i = 1,NR
  SR_tmp = SR(i)
  U_3    = U_3_kp1(i)
  Usb_3  = Usb_3_kp1(i)
  Usb_4  = Usb_4_kp1(i)
  CALL VISCOSITY(U_3,Usb_3,Usb_4,TIME,Lambda,SR_tmp,ETA_tmp) 
  ETA(i) = ETA_tmp
END DO
! --------------------------------------------------------------------------------- !
! TORQUE_r(1)  IS THE TORQUE APPLIED ON THE INNER CYLINDER FROM THE TEST MATERIAL.  !
! TORQUE_r(NR) IS THE TORQUE APPLIED FROM THE OUTER CYLINDER ON THE TEST MATERIAL.  !
! TORQUE_r(i)  IS THE TORQUE APPLIED FROM THE OUTER MATERIAL ON THE INNER MATERIAL. !
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
CALL TORQUE(V,ETA,NR,R_i,dr,h,TORQUE_r,Shear_Stress)
! --------------------------------------------------------------------------------- !
! WRITING INFORMATION TO FILES:
! --------------------------------------------------------------------------------- !
10 FORMAT(1X,300(F9.5,1X))
11 FORMAT(1X,300(F13.5,1X))
12 FORMAT(1X,2(F9.5,1X))
! --------------------------------------------------------------------------------- !
! WRITING THE RADIUS OF INNER- AND OUTER CYLINDER TO FILE:
OPEN(UNIT=8,file="RiRo_dr.dat",status="replace",action="write",&
     position="rewind",iostat=problem)
IF (problem /= 0) THEN
  PRINT *, " Could not create the file: RiRo_dr.dat! "
  WRITE2FILE_ERROR = .TRUE.
ELSE
  WRITE (UNIT=8,FMT=12) R_i
  WRITE (UNIT=8,FMT=12) DBLE(NR-1)*dr + R_i
  WRITE (UNIT=8,FMT=12) dr
END IF

CLOSE (UNIT=8)
! --------------------------------------------------------------------------------- !
! WRITING ANGULAR VELOCITY TO FILE:
OPEN(UNIT=8,file="time_omega.dat",status="replace",action="write",&
     position="rewind",iostat=problem)
IF (problem /= 0) THEN
  PRINT *, " Could not create the file: time_omega.dat! "
  WRITE2FILE_ERROR = .TRUE.
ELSE
  WRITE (UNIT=8,FMT=12) DBLE(kp1)*dt,omega
END IF

CLOSE (UNIT=8)
! --------------------------------------------------------------------------------- !
! WRITING COAGULATED STATE OF SIZE 3 CEMENT PARTICLES U_3 TO FILE:
OPEN(UNIT=8,file="U_3.dat",status="replace",action="write",&
     position="rewind",iostat=problem)
IF (problem /= 0) THEN
  PRINT *, " Could not create the file: U_3.dat! "
  WRITE2FILE_ERROR = .TRUE.
ELSE
  WRITE (UNIT=8,FMT=10) U_3_kp1(:)
END IF

CLOSE (UNIT=8)
! --------------------------------------------------------------------------------- !
! WRITING LINKED STATE OF SIZE 3 CEMENT PARTICLES Usb_3 TO FILE:
OPEN(UNIT=8,file="Usb_3.dat",status="replace",action="write",&
     position="rewind",iostat=problem)
IF (problem /= 0) THEN
  PRINT *, " Could not create the file: Usb_3.dat! "
  WRITE2FILE_ERROR = .TRUE.
ELSE
  WRITE (UNIT=8,FMT=10) Usb_3_kp1(:)
END IF

CLOSE (UNIT=8)
! --------------------------------------------------------------------------------- !
! WRITING LINKED STATE OF SIZE 4 CEMENT PARTICLES Usb_4 TO FILE:
OPEN(UNIT=8,file="Usb_4.dat",status="replace",action="write",&
     position="rewind",iostat=problem)
IF (problem /= 0) THEN
  PRINT *, " Could not create the file: Usb_4.dat! "
  WRITE2FILE_ERROR = .TRUE.
ELSE
  WRITE (UNIT=8,FMT=10) Usb_4_kp1(:)
END IF

CLOSE (UNIT=8)
! --------------------------------------------------------------------------------- !
! WRITING VELOCITY TO FILE:
OPEN(UNIT=8,file="velocity.dat",status="replace",action="write",&
     position="rewind",iostat=problem)
IF (problem /= 0) THEN
  PRINT *, " Could not create the file: velocity.dat! "
  WRITE2FILE_ERROR = .TRUE.
ELSE
  WRITE (UNIT=8,FMT=10) V(:)
END IF

CLOSE (UNIT=8)
! --------------------------------------------------------------------------------- !
! WRITING SHEAR RATE (i.e. gamma dot) TO FILE:
OPEN(UNIT=8,file="shear_rate.dat",status="replace",action="write",&
     position="rewind",iostat=problem)
IF (problem /= 0) THEN
  PRINT *, " Could not create the file: shear_rate.dat! "
  WRITE2FILE_ERROR = .TRUE.
ELSE
  WRITE (UNIT=8,FMT=10) SR(:)
END IF

CLOSE (UNIT=8)
! --------------------------------------------------------------------------------- !
! WRITING THE TIME DERIVATIVE OF SHEAR RATE (i.e. gamma double dot) TO FILE:        !
! A fiction value of 0.0D0 is put into the file "dSR_dt.dat" at t = 0 sec. This is  !
! OK, since dSR_dt at t = 0 sec is newer used (because the coagulated- and linked   !
! states U_3, Usb_3 and Usb_4 are always known at t = 0 sec; => i.e. [U_3 = U3_o at !
! t = 0 sec], [Usb_3 = Usb3_o at t = 0 sec], and [Usb_4 = Usb4_o at t = 0 sec]).    !
! The initial values U3_o, Usb3_o and Usb4_o are set in the file "viscous.f90".     !
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
OPEN(UNIT=8,file="dSR_dt.dat",status="replace",action="write",&
     position="rewind",iostat=problem)
IF (problem /= 0) THEN
  PRINT *, " Could not create the file: dSR_dt.dat! "
  WRITE2FILE_ERROR = .TRUE.
ELSE
  WRITE (UNIT=8,FMT=11) 0.0D0*V(:)
END IF

CLOSE (UNIT=8)
! --------------------------------------------------------------------------------- !
! WRITING THE DOUBLE TIME DERIVATIVE OF SHEAR RATE (i.e. gamma triple dot) TO FILE: !
! A fiction value of 0.0D0 is put into the file "ddSR_dtdt.dat" at t = 0 sec. This  !
! is OK, since ddSR_dtdt at t = 0 sec is newer used (because the coagulated- and    !
! linked states U_3, Usb_3 and Usb_4 are always known at t = 0 sec; => i.e. [U_3 =  !
! U3_o at t = 0 sec], [Usb_3 = Usb3_o at t = 0 sec], and [Usb_4 = Usb4_o at t = 0   !
! sec]). The initial values U3_o, Usb3_o and Usb4_o are set in "viscous.f90".       !
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
OPEN(UNIT=8,file="ddSR_dtdt.dat",status="replace",action="write",&
     position="rewind",iostat=problem)
IF (problem /= 0) THEN
  PRINT *, " Could not create the file: ddSR_dtdt.dat! "
  WRITE2FILE_ERROR = .TRUE.
ELSE
  WRITE (UNIT=8,FMT=11) 0.0D0*V(:)
END IF

CLOSE (UNIT=8)
! --------------------------------------------------------------------------------- !
! WRITING SHEAR VISCOSITY (i.e. apparent viscosity) TO FILE:
OPEN(UNIT=8,file="eta.dat",status="replace",action="write",&
     position="rewind",iostat=problem)
IF (problem /= 0) THEN
  PRINT *, " Could not create the file: eta.dat! "
  WRITE2FILE_ERROR = .TRUE.
ELSE
  WRITE (UNIT=8,FMT=11) ETA(:)
END IF

CLOSE (UNIT=8)
! --------------------------------------------------------------------------------- !
! WRITING SHEAR STRESS TO FILE:
OPEN(UNIT=8,file="shear_stress.dat",status="replace",action="write",&
     position="rewind",iostat=problem)
IF (problem /= 0) THEN
  PRINT *, " Could not create the file: shear_stress.dat! "
  WRITE2FILE_ERROR = .TRUE.
ELSE
  WRITE (UNIT=8,FMT=11) Shear_Stress(:)
END IF

CLOSE (UNIT=8)
! --------------------------------------------------------------------------------- !
! WRITING TORQUE TO FILE:
OPEN(UNIT=8,file="torque.dat",status="replace",action="write",&
     position="rewind",iostat=problem)
IF (problem /= 0) THEN
  PRINT *, " Could not create the file: torque.dat! "
  WRITE2FILE_ERROR = .TRUE.
ELSE
  WRITE (UNIT=8,FMT=12) DBLE(kp1)*dt,TORQUE_r(1)
END IF

CLOSE (UNIT=8)
! --------------------------------------------------------------------------------- !
! WRITING TORQUE(r) TO FILE:
OPEN(UNIT=8,file="torque_r.dat",status="replace",action="write",&
     position="rewind",iostat=problem)
IF (problem /= 0) THEN
  PRINT *, " Could not create the file: torque_r.dat! "
  WRITE2FILE_ERROR = .TRUE.
ELSE
  WRITE (UNIT=8,FMT=10) TORQUE_r(:)
END IF

CLOSE (UNIT=8)
! --------------------------------------------------------------------------------- !
IF (WRITE2FILE_ERROR) THEN
  PRINT *, " WRITE2FILE_ZERO: Could not create file(s)!               "
  PRINT *, " Error code 2 in 'write2f.f90' and execution terminated!  "
  STOP
END IF
! --------------------------------------------------------------------------------- !
DEALLOCATE(SR,ETA,TORQUE_r,Shear_Stress,stat=problem)
IF (problem /= 0) THEN 
  PRINT *, " WRITE2FILE_ZERO: The program could not deallocate space! "
  PRINT *, " Error code 3 in 'write2f.f90' and execution terminated!  "
  STOP
END IF
! --------------------------------------------------------------------------------- !
RETURN
END SUBROUTINE WRITE2FILE_ZERO
! ================================================================================= !
SUBROUTINE WRITE2FILE(U_3_kp1,Usb_3_kp1,Usb_4_kp1,V,dSR_dt,ddSR_dtdt,&
                      NR,kp1,dt,Lambda,R_i,dr,h,omega)
! --------------------------------------------------------------------------------- !
DOUBLE PRECISION,DIMENSION(:),INTENT(IN)  :: U_3_kp1,Usb_3_kp1,Usb_4_kp1,V,&
                                             dSR_dt,ddSR_dtdt
DOUBLE PRECISION,INTENT(IN)               :: dt,Lambda,R_i,dr,h,omega
INTEGER,INTENT(IN)                        :: NR,kp1

DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: SR,ETA,TORQUE_r,Shear_Stress
! --------------------------------------------------------------------------------- !
DOUBLE PRECISION :: U_3,Usb_3,Usb_4,TIME,ETA_tmp,SR_tmp
INTEGER          :: problem,i
LOGICAL          :: WRITE2FILE_ERROR
! --------------------------------------------------------------------------------- !
WRITE2FILE_ERROR = .FALSE.

TIME = DBLE(kp1)*dt

ALLOCATE(SR(NR),ETA(NR),TORQUE_r(NR),Shear_Stress(NR),stat=problem)
IF (problem /= 0) THEN 
  PRINT *, " WRITE2FILE: The program could not allocate space!       "
  PRINT *, " Error code 4 in 'write2f.f90' and execution terminated! "
  STOP
END IF

SR           = 0.0D0
ETA          = 0.0D0
TORQUE_r     = 0.0D0
Shear_Stress = 0.0D0
! --------------------------------------------------------------------------------- !
! CALCULATING THE SHEAR RATE AT ALL POINTS (shear.f90): SR 
CALL ROS_PROFILE(V,NR,R_i,dr,SR)
! --------------------------------------------------------------------------------- !
! CALCULATING THE SHEAR VISCOSITY (i.e. apparent viscosity) AT ALL POINTS:
ETA_tmp = 0.0D0
DO i = 1,NR
  SR_tmp = SR(i)
  U_3    = U_3_kp1(i)
  Usb_3  = Usb_3_kp1(i)
  Usb_4  = Usb_4_kp1(i)
  CALL VISCOSITY(U_3,Usb_3,Usb_4,TIME,Lambda,SR_tmp,ETA_tmp) 
  ETA(i) = ETA_tmp
END DO
! --------------------------------------------------------------------------------- !
! TORQUE_r(1)  IS THE TORQUE APPLIED ON THE INNER CYLINDER FROM THE TEST MATERIAL.  !
! TORQUE_r(NR) IS THE TORQUE APPLIED FROM THE OUTER CYLINDER ON THE TEST MATERIAL.  !
! TORQUE_r(i)  IS THE TORQUE APPLIED FROM THE OUTER MATERIAL ON THE INNER MATERIAL. !
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
CALL TORQUE(V,ETA,NR,R_i,dr,h,TORQUE_r,Shear_Stress)
! --------------------------------------------------------------------------------- !
! WRITING INFORMATION TO FILES:
! --------------------------------------------------------------------------------- !
10 FORMAT(1X,300(F9.5,1X))
11 FORMAT(1X,300(F13.5,1X))
12 FORMAT(1X,2(F9.5,1X))
! --------------------------------------------------------------------------------- !
! WRITING ANGULAR VELOCITY TO FILE:
OPEN(UNIT=8,file="time_omega.dat",status="old",action="write",&
     position="append",iostat=problem)
IF (problem /= 0) THEN
  PRINT *, " Could not open the file: time_omega.dat! "
  WRITE2FILE_ERROR = .TRUE.
ELSE
  WRITE (UNIT=8,FMT=12) DBLE(kp1)*dt,omega
END IF

CLOSE (UNIT=8)
! --------------------------------------------------------------------------------- !
! WRITING COAGULATED STATE OF SIZE 3 CEMENT PARTICLES U_3 TO FILE:
OPEN(UNIT=8,file="U_3.dat",status="old",action="write",&
     position="append",iostat=problem)
IF (problem /= 0) THEN
  PRINT *, " Could not open the file: U_3.dat! "
  WRITE2FILE_ERROR = .TRUE.
ELSE
  WRITE (UNIT=8,FMT=10) U_3_kp1(:)
END IF

CLOSE (UNIT=8)
! --------------------------------------------------------------------------------- !
! WRITING LINKED STATE OF SIZE 3 CEMENT PARTICLES Usb_3 TO FILE:
OPEN(UNIT=8,file="Usb_3.dat",status="old",action="write",&
     position="append",iostat=problem)
IF (problem /= 0) THEN
  PRINT *, " Could not open the file: Usb_3.dat! "
  WRITE2FILE_ERROR = .TRUE.
ELSE
  WRITE (UNIT=8,FMT=10) Usb_3_kp1(:)
END IF

CLOSE (UNIT=8)
! --------------------------------------------------------------------------------- !
! WRITING LINKED STATE OF SIZE 4 CEMENT PARTICLES Usb_4 TO FILE:
OPEN(UNIT=8,file="Usb_4.dat",status="old",action="write",&
     position="append",iostat=problem)
IF (problem /= 0) THEN
  PRINT *, " Could not open the file: Usb_4.dat! "
  WRITE2FILE_ERROR = .TRUE.
ELSE
  WRITE (UNIT=8,FMT=10) Usb_4_kp1(:)
END IF

CLOSE (UNIT=8)
! --------------------------------------------------------------------------------- !
! WRITING VELOCITY TO FILE:
OPEN(UNIT=8,file="velocity.dat",status="old",action="write",&
     position="append",iostat=problem)
IF (problem /= 0) THEN
  PRINT *, " Could not open the file: velocity.dat! "
  WRITE2FILE_ERROR = .TRUE.
ELSE
  WRITE (UNIT=8,FMT=10) V(:)
END IF

CLOSE (UNIT=8)
! --------------------------------------------------------------------------------- !
! WRITING SHEAR RATE (i.e. gamma dot) TO FILE:
OPEN(UNIT=8,file="shear_rate.dat",status="old",action="write",&
     position="append",iostat=problem)
IF (problem /= 0) THEN
  PRINT *, " Could not open the file: shear_rate.dat! "
  WRITE2FILE_ERROR = .TRUE.
ELSE
  WRITE (UNIT=8,FMT=10) SR(:)
END IF

CLOSE (UNIT=8)
! --------------------------------------------------------------------------------- !
! WRITING THE TIME DERIVATIVE OF SHEAR RATE (i.e. gamma double dot) TO FILE:
OPEN(UNIT=8,file="dSR_dt.dat",status="old",action="write",&
     position="append",iostat=problem)
IF (problem /= 0) THEN
  PRINT *, " Could not open the file: dSR_dt.dat! "
  WRITE2FILE_ERROR = .TRUE.
ELSE
  WRITE (UNIT=8,FMT=11) dSR_dt(:)
END IF

CLOSE (UNIT=8)
! --------------------------------------------------------------------------------- !
! WRITING THE DOUBLE TIME DERIVATIVE OF SHEAR RATE (i.e. gamma triple dot) TO FILE:
OPEN(UNIT=8,file="ddSR_dtdt.dat",status="old",action="write",&
     position="append",iostat=problem)
IF (problem /= 0) THEN
  PRINT *, " Could not open the file: ddSR_dtdt.dat! "
  WRITE2FILE_ERROR = .TRUE.
ELSE
  WRITE (UNIT=8,FMT=11) ddSR_dtdt(:)
END IF

CLOSE (UNIT=8)
! --------------------------------------------------------------------------------- !
! WRITING SHEAR VISCOSITY (i.e. apparent viscosity) TO FILE:
OPEN(UNIT=8,file="eta.dat",status="old",action="write",&
     position="append",iostat=problem)
IF (problem /= 0) THEN
  PRINT *, " Could not open the file: eta.dat! "
  WRITE2FILE_ERROR = .TRUE.
ELSE
  WRITE (UNIT=8,FMT=11) ETA(:)
END IF

CLOSE (UNIT=8)
! --------------------------------------------------------------------------------- !
! WRITING SHEAR STRESS TO FILE:
OPEN(UNIT=8,file="shear_stress.dat",status="old",action="write",&
     position="append",iostat=problem)
IF (problem /= 0) THEN
  PRINT *, " Could not open the file: shear_stress.dat! "
  WRITE2FILE_ERROR = .TRUE.
ELSE
  WRITE (UNIT=8,FMT=11) Shear_Stress(:)
END IF

CLOSE (UNIT=8)
! --------------------------------------------------------------------------------- !
! WRITING TORQUE TO FILE:
OPEN(UNIT=8,file="torque.dat",status="old",action="write",&
     position="append",iostat=problem)
IF (problem /= 0) THEN
  PRINT *, " Could not open the file: torque.dat! "
  WRITE2FILE_ERROR = .TRUE.
ELSE
  WRITE (UNIT=8,FMT=12) DBLE(kp1)*dt,TORQUE_r(1)
END IF

CLOSE (UNIT=8)
! --------------------------------------------------------------------------------- !
! WRITING TORQUE(r) TO FILE:
OPEN(UNIT=8,file="torque_r.dat",status="old",action="write",&
     position="append",iostat=problem)
IF (problem /= 0) THEN
  PRINT *, " Could not open the file: torque_r.dat! "
  WRITE2FILE_ERROR = .TRUE.
ELSE
  WRITE (UNIT=8,FMT=10) TORQUE_r(:)
END IF

CLOSE (UNIT=8)
! --------------------------------------------------------------------------------- !
IF (WRITE2FILE_ERROR) THEN
  PRINT *, " WRITE2FILE: Could not create file(s)!                   "
  PRINT *, " Error code 5 in 'write2f.f90' and execution terminated! "
  STOP
END IF
! --------------------------------------------------------------------------------- !
DEALLOCATE(SR,ETA,TORQUE_r,Shear_Stress,stat=problem)
IF (problem /= 0) THEN 
  PRINT *, " WRITE2FILE: The program could not deallocate space!     "
  PRINT *, " Error code 6 in 'write2f.f90' and execution terminated! "
  STOP
END IF
! --------------------------------------------------------------------------------- !
RETURN
END SUBROUTINE WRITE2FILE
! ================================================================================= !
SUBROUTINE ARRAY_VECTOR_debug(M,D,v,dim_n)
! --------------------------------------------------------------------------------- !
DOUBLE PRECISION,DIMENSION(:,:),INTENT(IN) :: M
DOUBLE PRECISION,DIMENSION(:),INTENT(IN)   :: D,v
INTEGER,INTENT(IN)                         :: dim_n
! --------------------------------------------------------------------------------- !
INTEGER :: problem,i
LOGICAL :: WRITE2FILE_ERROR
! --------------------------------------------------------------------------------- !
WRITE2FILE_ERROR = .FALSE.

10 FORMAT(1X,300(F10.4,1X))
11 FORMAT(1X,F10.4)
! --------------------------------------------------------------------------------- !
OPEN(UNIT=8,file="ARRAY_debug.dat",status="replace",action="write",&
     position="rewind",iostat=problem)
IF (problem /= 0) THEN
  PRINT *, " Could not create the file: ARRAY_debug.dat! "
  WRITE2FILE_ERROR = .TRUE.
ELSE
  DO i = 1,dim_n
    WRITE (UNIT=8,FMT=10) M(i,:)
  END DO
END IF

CLOSE (UNIT=8)
! --------------------------------------------------------------------------------- !
OPEN(UNIT=8,file="VECTOR_debug.dat",status="replace",action="write",&
     position="rewind",iostat=problem)
IF (problem /= 0) THEN
  PRINT *, " Could not create the file: VECTOR_debug.dat! "
  WRITE2FILE_ERROR = .TRUE.
ELSE
  DO i = 1,dim_n
    WRITE (UNIT=8,FMT=11) D(i)
  END DO
END IF

CLOSE (UNIT=8)
! --------------------------------------------------------------------------------- !
OPEN(UNIT=8,file="vel_debug.dat",status="replace",action="write",&
     position="rewind",iostat=problem)
IF (problem /= 0) THEN
  PRINT *, " Could not create the file: vel_debug.dat! "
  WRITE2FILE_ERROR = .TRUE.
ELSE
  DO i = 1,dim_n
    WRITE (UNIT=8,FMT=11) v(i)
  END DO
END IF

CLOSE (UNIT=8)
! --------------------------------------------------------------------------------- !
IF (WRITE2FILE_ERROR) THEN
  PRINT *, " ARRAY_VECTOR_debug: Could not create file(s)!           "
  PRINT *, " Error code 7 in 'write2f.f90' and execution terminated! "
  STOP
END IF
! --------------------------------------------------------------------------------- !
RETURN
END SUBROUTINE ARRAY_VECTOR_debug
! ================================================================================= !
SUBROUTINE TORQUE(V,ETA,NR,R_i,dr,h,TORQUE_r,Shear_Stress)
! --------------------------------------------------------------------------------- !
DOUBLE PRECISION,DIMENSION(:),INTENT(IN)  :: V,ETA
DOUBLE PRECISION,DIMENSION(:),INTENT(OUT) :: TORQUE_r,Shear_Stress

INTEGER,INTENT(IN)                        :: NR
DOUBLE PRECISION,INTENT(IN)               :: R_i,dr,h
! --------------------------------------------------------------------------------- !
INTEGER          :: i
DOUBLE PRECISION :: r,PI,pmSR
! --------------------------------------------------------------------------------- !
PI = DACOS(-1.0D0)
! --------------------------------------------------------------------------------- !
! In this software (i.e. in VVPF ver. 2.0), it is assumed that the cylinder (i.e.   !
! the outer- or the inner cylinder) rotates counterclockwise. Hence, when the outer !
! cylinder rotates, the torque values (i.e. "TORQUE_r(i)") should be positive,      !
! meaning that the test material is trying to rotate the stationary inner cylinder  !
! counterclockwise. However, if it is the inner cylinder that rotates, the torque   !
! should be negative, meaning that the test material is trying to slow down the     !
! rotation of the inner cylinder (i.e. this is the torque applied from the test     !
! material on the inner cylinder [however, the torque applied FROM the inner        !
! cylinder ON the test material should be positive, which is simply gained by       !
! putting a minus sign in front of the current calculated values, and would then    !
! result in positive values]). See Equations 3.12 and 3.16 (Pages 57 and 58) in     !
! the Ph.D.-thesis (mentioned in the beginning of the file "main.f90") about the    !
! calculation of torque. In either case of rotating outer- or of rotating inner     !
! cylinder, these equations are valid and are used in this software. When plotting  !
! the torque as a function of time for the latter case, it is more convenient to    !
! use its absolute value.                                                           !
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
! Rotation is assumed counterclockwise in this software, but a clockwise rotation   !
! in the viscometer is also valid. That is, this software can still be used to      !
! analyze the outcome when a clockwise rotation is present (the test material does  !
! not care if it is rotated clockwise or counterclockwise). See Page 57 (and 58) in !
! the Ph.D.-thesis mentioned in the beginning of the file 'main.f90', about this    !
! subject. The only change is mathematical, in terms of +/- sign in torque values.  !
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
! TORQUE_r(1)  IS THE TORQUE APPLIED ON THE INNER CYLINDER FROM THE TEST MATERIAL.  !
! TORQUE_r(NR) IS THE TORQUE APPLIED FROM THE OUTER CYLINDER ON THE TEST MATERIAL.  !
! TORQUE_r(i)  IS THE TORQUE APPLIED FROM THE OUTER MATERIAL ON THE INNER MATERIAL. !
! --------------------------------------------------------------------------------- !
! Shear stress (i.e. "Shear_Stress(i)") is calculated by Equation 7.5 (Page 156)    !
! and is always a positive value. See also the last paragraph in Section 3.2.1      !
! (Page 53) about the shear stress and why it must always be a positive value.      !
! --------------------------------------------------------------------------------- !
! CALCULATING THE SHEAR STRESS AND TORQUE AT THE INNER CYLINDER (r = R_i):
i    = 1
r    = DBLE(i-1)*dr + R_i
pmSR = ((4.0D0*V(i+1) - V(i+2) - 3.0D0*V(i))/(2.0D0*dr) - V(i)/r)

Shear_Stress(i) = ETA(i)*DABS(pmSR)
TORQUE_r(i)     = r*((ETA(i)*pmSR)*h*(2*PI*r))
! --------------------------------------------------------------------------------- !
! CALCULATING THE SHEAR STRESS AND TORQUE IN BULK (i.e. within the test material):
DO i = 2,NR-1
  r    = DBLE(i-1)*dr + R_i
  pmSR = ((V(i+1) - V(i-1))/(2.0D0*dr) - V(i)/r)

  Shear_Stress(i) = ETA(i)*DABS(pmSR)
  TORQUE_r(i)     = r*((ETA(i)*pmSR)*h*(2*PI*r))
END DO
! --------------------------------------------------------------------------------- !
! CALCULATING THE SHEAR STRESS AND TORQUE AT THE OUTER CYLINDER (r = R_o):
i    = NR
r    = DBLE(i-1)*dr + R_i
pmSR = ((-4.0D0*V(i-1) + V(i-2) + 3.0D0*V(i))/(2.0D0*dr) - V(i)/r)

Shear_Stress(i) = ETA(i)*DABS(pmSR)
TORQUE_r(i)     = r*((ETA(i)*pmSR)*h*(2*PI*r))
! --------------------------------------------------------------------------------- !
RETURN
END SUBROUTINE TORQUE
! ================================================================================= !
END MODULE WRITE_INFORMATION
! --------------------------------------------------------------------------------- !
