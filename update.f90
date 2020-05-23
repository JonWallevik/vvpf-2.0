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
! File name: update.f90 (MODULE)                                                    !
! This file sets up the system of algebraic equations (see for example Equations    !
! 7.28 to 7.31 in the Ph.D.-thesis mentioned in the file "main.f90"). This file     !
! also contains the Thomas algorithm that is used in solving this system.           !
! --------------------------------------------------------------------------------- !
MODULE MATRIX
  USE RATE_OF_SHEAR
  USE SHEAR_VISCOSITY
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: MATRIX_UPDATE,MATRIX_SOLVER
CONTAINS
! ================================================================================= !
SUBROUTINE MATRIX_UPDATE(U_3_kp1,Usb_3_kp1,Usb_4_kp1,rho,k,dt,Lambda,&
                         dr,R_i,NR,v_k,v_kp1,ARRAY,VECTOR)
! --------------------------------------------------------------------------------- !
DOUBLE PRECISION,DIMENSION(:,:),INTENT(OUT) :: ARRAY
DOUBLE PRECISION,DIMENSION(:),INTENT(OUT)   :: VECTOR

DOUBLE PRECISION,DIMENSION(:),INTENT(IN)    :: U_3_kp1,Usb_3_kp1,Usb_4_kp1,v_k,v_kp1

DOUBLE PRECISION,INTENT(IN)                 :: dt,dr,R_i,Lambda,rho
INTEGER,INTENT(IN)                          :: k,NR
! --------------------------------------------------------------------------------- !
DOUBLE PRECISION :: &
           U_3_im1kp1,    U_3_ikp1,    U_3_ip1kp1,    U_3_im12kp1,    U_3_ip12kp1,&
           Usb_3_im1kp1,  Usb_3_ikp1,  Usb_3_ip1kp1,  Usb_3_im12kp1,  Usb_3_ip12kp1,&
           Usb_4_im1kp1,  Usb_4_ikp1,  Usb_4_ip1kp1,  Usb_4_im12kp1,  Usb_4_ip12kp1

DOUBLE PRECISION :: SR_ikp1,   SR_ip12kp1,   SR_im12kp1,&
                    ETA_ikp1,  ETA_ip12kp1,  ETA_im12kp1

DOUBLE PRECISION :: V_im1k,    V_ik,    V_ip1k,&
                    V_im1kp1,  V_ikp1,  V_ip1kp1,&
                    A_kp1,     B_kp1,   C_kp1

DOUBLE PRECISION :: BETA,THETA,rp1,r,rm1,TIME_k,TIME_kp1

INTEGER          :: I
! --------------------------------------------------------------------------------- !
BETA     =  dt/(rho*dr)          !                                                  !
TIME_k   =  DBLE(k)*dt           ! <= The solution for this time is known.          !
TIME_kp1 = (DBLE(k) + 1.0D0)*dt  ! <= This is the time, which is solved for.        !
! --------------------------------------------------------------------------------- !
ARRAY    =  0.0D0
VECTOR   =  0.0D0
! --------------------------------------------------------------------------------- !
! I = 1    => i = 2 in main.f90, i.e. near (and not at) the inner cylinder R_i.     !
! I = NR-2 => i = NR-1 in main.f90, i.e. near (and not at) the outer cylinder R_o.  !
! An array of size "(NR-2)x(NR-2)" is generated in the following routine.           !
! --------------------------------------------------------------------------------- !
DO I = 1,NR-2

  rp1           = DBLE(I+1)*dr + R_i
  r             = DBLE(I)*dr   + R_i
  rm1           = DBLE(I-1)*dr + R_i
  THETA         = (2.0D0*dt)/(rho*r)

  V_im1k        = v_k(I)
  V_ik          = v_k(I+1)
  V_ip1k        = v_k(I+2)

  V_im1kp1      = v_kp1(I)
  V_ikp1        = v_kp1(I+1)
  V_ip1kp1      = v_kp1(I+2)

  U_3_im1kp1    = U_3_kp1(I)
  U_3_ikp1      = U_3_kp1(I+1)
  U_3_ip1kp1    = U_3_kp1(I+2)

  U_3_im12kp1   = (U_3_ikp1 + U_3_im1kp1)/2.0D0
  U_3_ip12kp1   = (U_3_ip1kp1 + U_3_ikp1)/2.0D0

  Usb_3_im1kp1  = Usb_3_kp1(I)
  Usb_3_ikp1    = Usb_3_kp1(I+1)
  Usb_3_ip1kp1  = Usb_3_kp1(I+2)

  Usb_3_im12kp1 = (Usb_3_ikp1 + Usb_3_im1kp1)/2.0D0
  Usb_3_ip12kp1 = (Usb_3_ip1kp1 + Usb_3_ikp1)/2.0D0

  Usb_4_im1kp1  = Usb_4_kp1(I)
  Usb_4_ikp1    = Usb_4_kp1(I+1)
  Usb_4_ip1kp1  = Usb_4_kp1(I+2)

  Usb_4_im12kp1 = (Usb_4_ikp1 + Usb_4_im1kp1)/2.0D0
  Usb_4_ip12kp1 = (Usb_4_ip1kp1 + Usb_4_ikp1)/2.0D0

  CALL ROS(rp1,r,rm1,dr,V_ikp1,V_ip1kp1,V_im1kp1,SR_ikp1,SR_ip12kp1,SR_im12kp1)

  CALL ETA(TIME_kp1,    Lambda,&
           SR_ikp1,     SR_ip12kp1,     SR_im12kp1,&
           U_3_ikp1,    U_3_ip12kp1,    U_3_im12kp1,&
           Usb_3_ikp1,  Usb_3_ip12kp1,  Usb_3_im12kp1,&
           Usb_4_ikp1,  Usb_4_ip12kp1,  Usb_4_im12kp1,&
           ETA_ikp1,    ETA_ip12kp1,    ETA_im12kp1)

  A_kp1 = (BETA*ETA_ip12kp1)*(1.0D0/dr - 1.0D0/(rp1+r)) + (THETA*ETA_ikp1)/(2.0D0*dr)
  B_kp1 = (BETA*ETA_ip12kp1)*(1.0D0/dr + 1.0D0/(rp1+r)) + &
          (BETA*ETA_im12kp1)*(1.0D0/dr - 1.0D0/(r+rm1)) + (THETA*ETA_ikp1)/r 
  C_kp1 = (BETA*ETA_im12kp1)*(1.0D0/dr + 1.0D0/(r+rm1)) - (THETA*ETA_ikp1)/(2.0D0*dr)

  ! ------------------------------------------------------------------------------- !
  ! Note, since "I = 1" is equivalent to "i = 2" then "im1 = i-1 = 1" (rm1 = R_i).  !
  ! Since "I = NR-2" is equivalent to "i = NR-1" then "ip1 = i+1 = NR" (rp1 = R_o). !
  ! ------------------------------------------------------------------------------- !
  IF (I == 1) THEN
    ARRAY(I,I)   = -(1.0D0 + B_kp1)
    ARRAY(I,I+1) =   A_kp1
    VECTOR(I)    = - V_ik - C_kp1*V_im1kp1
  ELSE IF (I == NR-2) THEN
    ARRAY(I,I-1) =   C_kp1
    ARRAY(I,I)   = -(1.0D0 + B_kp1)
    VECTOR(I)    = - V_ik - A_kp1*V_ip1kp1
  ELSE 
    ARRAY(I,I-1) =   C_kp1
    ARRAY(I,I)   = -(1.0D0 + B_kp1)
    ARRAY(I,I+1) =   A_kp1
    VECTOR(I)    = - V_ik
  END IF

END DO
! --------------------------------------------------------------------------------- !
RETURN
END SUBROUTINE MATRIX_UPDATE
! ================================================================================= !
SUBROUTINE MATRIX_SOLVER(M,D,v,dim_n)
! --------------------------------------------------------------------------------- !
! Subroutine that solves the trigonal system "M * v = D" with the Thomas algorithm. !
! This routine is also known as the: "Crout reduction for tridiagonal linear        !
! systems" - algorithm.                                                             !
! M ->  Left side of the linear system (a tridiagonal array).                       !
! D ->  Right side of the linear system (a vector).                                 !
! v ->  The variable to be solved: v = (M)^(-1) * D.                                !
! --------------------------------------------------------------------------------- !
DOUBLE PRECISION,DIMENSION(:,:),INTENT(INOUT) :: M
DOUBLE PRECISION,DIMENSION(:),INTENT(INOUT)   :: D
DOUBLE PRECISION,DIMENSION(:),INTENT(OUT)     :: v
INTEGER,INTENT(IN)                            :: dim_n
! --------------------------------------------------------------------------------- !
DOUBLE PRECISION :: coef
INTEGER          :: I
! --------------------------------------------------------------------------------- !
v = -1.0D0
! --------------------------------------------------------------------------------- !
! Forward and then back substitution:
DO I = 1,dim_n-1
  coef       = M(I+1,(I-1)+1)/M(I,I)
  M(I+1,I+1) = M(I+1,I+1) - coef*M(I,I+1)
  D(I+1)     = D(I+1) - coef*D(I)
END DO
v(dim_n) = D(dim_n)/M(dim_n,dim_n)
DO I = dim_n-1,1,-1
  v(I) = (D(I) - M(I,I+1)*v(I+1))/M(I,I)
END DO
! --------------------------------------------------------------------------------- !
RETURN
END SUBROUTINE MATRIX_SOLVER
! ================================================================================= !
END MODULE MATRIX
! --------------------------------------------------------------------------------- !
