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
! File name: motion.f90 (MODULE)                                                    !
! This file reads the basic information from the "COMMON /RAW_DATA/" to produce the !
! angular velocity "omega". The information about the angular velocity is requested !
! by the routine main.f90.                                                          !
! --------------------------------------------------------------------------------- !
MODULE ROTATION
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: ANGULAR_VELOCITY
CONTAINS
! ================================================================================= !
SUBROUTINE ANGULAR_VELOCITY(dbl_kp1,dt,omega)
! --------------------------------------------------------------------------------- !
DOUBLE PRECISION,INTENT(IN)  :: dbl_kp1,dt
DOUBLE PRECISION,INTENT(OUT) :: omega
! --------------------------------------------------------------------------------- !
DOUBLE PRECISION :: time,omega_a,omega_b,time_a,time_b
INTEGER          :: N_read,NT_LOCK,i
! --------------------------------------------------------------------------------- !
! The value "1501" used below must always be the same as the value of "N_read".     !
! The corresponding value should also be changed in the file "main.f90".            !
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
DOUBLE PRECISION,DIMENSION(1501) :: raw_time_data,raw_omega_data
COMMON /RAW_DATA/ raw_time_data,raw_omega_data
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
N_read  = 1501        ! -> The length of raw data points to be read.                !
! --------------------------------------------------------------------------------- !
time    = dbl_kp1*dt  ! => time = (k+1)*dt, c.f. main.f90                           !
! --------------------------------------------------------------------------------- !
NT_LOCK = N_read - 1  ! => This statement is necessary for the last time step, i.e. !
                      !    at k+1 = N_dt + 1 = (NUMBER_OF_TIME_ITERATIONS - 1) + 1  !
                      !    = NUMBER_OF_TIME_ITERATIONS = REAL_TIME/dt_Plastic       !
                      !    (e.g. at (k+1)*dt = REAL_TIME = 50 s). This is because   !
                      !    the "IF"-sentence below is never executed for the last   !
                      !    time iteration, as its condition never becomes TRUE.     !
! --------------------------------------------------------------------------------- !
! "raw_time_data(i) <= time < raw_time_data(i+1)"
DO i = 1,N_read-1
  IF ((raw_time_data(i) <= time).AND.(time < raw_time_data(i+1))) THEN
    NT_LOCK = i
    EXIT
  END IF
END DO

omega_a = raw_omega_data(NT_LOCK)
omega_b = raw_omega_data(NT_LOCK + 1)
time_a  = raw_time_data(NT_LOCK)
time_b  = raw_time_data(NT_LOCK + 1)

omega   = omega_a + ((omega_b - omega_a)/(time_b - time_a))*(time - time_a)
! --------------------------------------------------------------------------------- !
RETURN
END SUBROUTINE ANGULAR_VELOCITY
! ================================================================================= !
END MODULE ROTATION
! --------------------------------------------------------------------------------- !
