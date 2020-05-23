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
! File name: shear.f90 (MODULE)                                                     !
! This routine calculates the shear rate SR(i), its time derivative dSR_dt(i) and   !
! its double time derivative ddSR_dtdt(i), from the computed velocity profile       !
! VELOCITY_k(i). This information is requested by update.f90 and main.f90. The      !
! shear rate is calculated by Equation 3.10, Page 57. See Section 7.5 (Page 162)    !
! about an example of numerical implementation. Note that "ROS" and "SR" means the  !
! same thing; i.e. ROS = rate of shear = SR = shear rate. General definition of     !
! shear rate (i.e. rate of shear) is given by Equation 2.24, Page 18.               !
! --------------------------------------------------------------------------------- !
MODULE RATE_OF_SHEAR
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: ROS,ROS_PROFILE,dROS_dt,ddROS_dtdt
CONTAINS
! ================================================================================= !
! Calculation of the shear rate at a single grid point (gamma dot):
SUBROUTINE ROS(rp1,r,rm1,dr,V_i,V_ip1,V_im1,SR_i,SR_ip12,SR_im12)
! --------------------------------------------------------------------------------- !
DOUBLE PRECISION,INTENT(IN)  :: rp1,r,rm1,dr,V_i,V_ip1,V_im1
DOUBLE PRECISION,INTENT(OUT) :: SR_i,SR_ip12,SR_im12
! --------------------------------------------------------------------------------- !
SR_i    = DABS((V_ip1 - V_im1)/(2.0D0*dr) - V_i/r)
SR_ip12 = DABS((V_ip1 - V_i)/dr - (V_ip1 + V_i)/(rp1 + r))
SR_im12 = DABS((V_i - V_im1)/dr - (V_i + V_im1)/(r + rm1))
! --------------------------------------------------------------------------------- !
RETURN
END SUBROUTINE ROS
! ================================================================================= !
! Calculation of the shear rate at all grid points (gamma dot):
SUBROUTINE ROS_PROFILE(V,NR,R_i,dr,SR)
! --------------------------------------------------------------------------------- !
DOUBLE PRECISION,DIMENSION(:),INTENT(IN)  :: V
DOUBLE PRECISION,DIMENSION(:),INTENT(OUT) :: SR
INTEGER,INTENT(IN)                        :: NR
DOUBLE PRECISION,INTENT(IN)               :: R_i,dr
! --------------------------------------------------------------------------------- !
INTEGER          :: i
DOUBLE PRECISION :: r,EPS
! --------------------------------------------------------------------------------- !
EPS = 1.0D-15
! --------------------------------------------------------------------------------- !
! CALCULATING THE SHEAR RATE IN THE BULK:                                           !
! See Section 7.5 about the formulas for the shear rate (SR). Note that ROS and SR  !
! means the same thing: ROS = rate of shear = SR = shear rate.                      !
! --------------------------------------------------------------------------------- !
DO i = 2,NR-1
  r     = DBLE(i-1)*dr + R_i
  SR(i) = DABS((V(i+1) - V(i-1))/(2.0D0*dr) - V(i)/r)
END DO
! ---------------------------------------------------------
! CALCULATING THE SHEAR RATE ON THE "LEFT" WALL (r=R_i):
i     = 1
r     = DBLE(i-1)*dr + R_i
SR(i) = DABS((4.0D0*V(i+1) - V(i+2) - 3.0D0*V(i))/(2.0D0*dr) - V(i)/r)
! ---------------------------------------------------------
! CALCULATING THE SHEAR RATE ON THE "RIGHT" WALL (r=R_o):
i     = NR
r     = DBLE(i-1)*dr + R_i
SR(i) = DABS((-4.0D0*V(i-1) + V(i-2) + 3.0D0*V(i))/(2.0D0*dr) - V(i)/r)
! --------------------------------------------------------------------------------- !
RETURN
END SUBROUTINE ROS_PROFILE
! ================================================================================= !
! Calculation of the time derivative of the shear rate at all grid points,
! (gamma double dot):
SUBROUTINE dROS_dt(k,NR,N_history,time_history,SR_history,dSR_dt)
! --------------------------------------------------------------------------------- !
DOUBLE PRECISION,DIMENSION(:,:),INTENT(IN) :: SR_history
DOUBLE PRECISION,DIMENSION(:),INTENT(IN)   :: time_history
DOUBLE PRECISION,DIMENSION(:),INTENT(OUT)  :: dSR_dt
INTEGER,INTENT(IN)                         :: k,NR,N_history
! --------------------------------------------------------------------------------- !
INTEGER          :: i,N,N_start,N_stop
DOUBLE PRECISION :: slope
! --------------------------------------------------------------------------------- !
! If everything else fails, one can use...                                          !
! dSR_dt = (SR_history(N_history,:) - SR_history(N_history-1,:))/dt                 !
! ...however, this approach is more sensitive to truncation and instability errors. !
! --------------------------------------------------------------------------------- !
IF ((k+2) >= N_history) THEN
  N = N_history
  DO i = 1,NR
    CALL SLOPE_CALC(N,time_history,SR_history(:,i),slope)
    dSR_dt(i) = slope
  END DO
ELSE
  N_start = N_history - (k + 1)
  N_stop  = N_history
  N       = N_stop - N_start + 1
  DO i = 1,NR
    CALL SLOPE_CALC(N,time_history(N_start:N_stop),&
                    SR_history(N_start:N_stop,i),slope)
    dSR_dt(i) = slope
  END DO
END IF
! --------------------------------------------------------------------------------- !
RETURN
END SUBROUTINE dROS_dt
! ================================================================================= !
! Calculation of the double time derivative of the shear rate at all grid points,
! (gamma triple dot):
SUBROUTINE ddROS_dtdt(k,NR,N_history,time_history,dSR_dt_history,ddSR_dtdt)
! --------------------------------------------------------------------------------- !
DOUBLE PRECISION,DIMENSION(:,:),INTENT(IN) :: dSR_dt_history
DOUBLE PRECISION,DIMENSION(:),INTENT(IN)   :: time_history
DOUBLE PRECISION,DIMENSION(:),INTENT(OUT)  :: ddSR_dtdt
INTEGER,INTENT(IN)                         :: k,NR,N_history
! --------------------------------------------------------------------------------- !
INTEGER          :: i,N,N_start,N_stop
DOUBLE PRECISION :: slope,dt
! --------------------------------------------------------------------------------- !
dt        = time_history(N_history) - time_history(N_history-1)
ddSR_dtdt = (dSR_dt_history(N_history,:) - dSR_dt_history(N_history-1,:))/dt
! --------------------------------------------------------------------------------- !
! If everything else fails, one can use the routine below. This routine is more     !
! stable (meaning less oscilations in ddSR_dtdt), but it will give more incorrect   !
! values. In addition, this routine gives a much longer computation time.           !
! --------------------------------------------------------------------------------- !
! IF ((k+2) >= N_history) THEN
!   N = N_history
!   DO i = 1,NR
!     CALL SLOPE_CALC(N,time_history,dSR_dt_history(:,i),slope)
!     ddSR_dtdt(i) = slope
!   END DO
! ELSE
!   N_start = N_history - (k + 1)
!   N_stop  = N_history
!   N       = N_stop - N_start + 1
!   DO i = 1,NR
!     CALL SLOPE_CALC(N,time_history(N_start:N_stop),&
!                     dSR_dt_history(N_start:N_stop,i),slope)
!     ddSR_dtdt(i) = slope
!   END DO
! END IF
! --------------------------------------------------------------------------------- !
RETURN
END SUBROUTINE ddROS_dtdt
! ================================================================================= !
! Slope calculation of the data points (x_i,y_i) by the method of least squares.
! The number of data points is N (that is, i = 1,2,3,... N).
SUBROUTINE SLOPE_CALC(N,x,y,slope)
! --------------------------------------------------------------------------------- !
DOUBLE PRECISION,DIMENSION(:),INTENT(IN) :: x,y
DOUBLE PRECISION,INTENT(OUT)             :: slope
INTEGER,INTENT(IN)                       :: N
! --------------------------------------------------------------------------------- !
DOUBLE PRECISION :: SUM_x,SUM_y,SUM_xy,SUM_xx
! --------------------------------------------------------------------------------- !
SUM_x  = SUM(x)
SUM_y  = SUM(y)
SUM_xy = DOT_PRODUCT(x,y)
SUM_xx = DOT_PRODUCT(x,x)

slope = (SUM_x*SUM_y - N*SUM_xy)/(SUM_x*SUM_x - N*SUM_xx)
! --------------------------------------------------------------------------------- !
RETURN
END SUBROUTINE SLOPE_CALC
! ================================================================================= !
END MODULE RATE_OF_SHEAR
! --------------------------------------------------------------------------------- !
