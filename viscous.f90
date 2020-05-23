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
! File name: viscous.f90 (MODULE)                                                   !
! In this file, the shear viscosity function ETA = ETA(SR,t,...) is defined and     !
! calculated. This information is requested by update.f90 and write2f.f90.          !
! --------------------------------------------------------------------------------- !
! NOTE [IMPORTANT]: If the variable "ddSR_dtdt" is used in the material functions   !
! CDF3, SBF3 or SBF4, defined in this file (viscous.f90), then it is recommended    !
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
MODULE SHEAR_VISCOSITY
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: ETA,VISCOSITY,COAG_DISP_3,STRUCT_BREAK_3,STRUCT_BREAK_4
CONTAINS
! ================================================================================= !
! -*-*-*-*-*-*-*-*- [The ETA subroutine is called by update.f90 ] -*-*-*-*-*-*-*-*- !
SUBROUTINE ETA(time,     Lambda,&
               SR_i,     SR_ip12,     SR_im12,&
               U_3_i,    U_3_ip12,    U_3_im12,&
               Usb_3_i,  Usb_3_ip12,  Usb_3_im12,&
               Usb_4_i,  Usb_4_ip12,  Usb_4_im12,&
               ETA_i,    ETA_ip12,    ETA_im12)
! --------------------------------------------------------------------------------- !
DOUBLE PRECISION,INTENT(IN)  :: time,     Lambda,&
                                SR_i,     SR_ip12,     SR_im12,&
                                U_3_i,    U_3_ip12,    U_3_im12,&
                                Usb_3_i,  Usb_3_ip12,  Usb_3_im12,&
                                Usb_4_i,  Usb_4_ip12,  Usb_4_im12

DOUBLE PRECISION,INTENT(OUT) :: ETA_i,    ETA_ip12,    ETA_im12
! --------------------------------------------------------------------------------- !
CALL VISCOSITY(U_3_i,    Usb_3_i,    Usb_4_i,    time, Lambda, SR_i,    ETA_i)
CALL VISCOSITY(U_3_ip12, Usb_3_ip12, Usb_4_ip12, time, Lambda, SR_ip12, ETA_ip12)
CALL VISCOSITY(U_3_im12, Usb_3_im12, Usb_4_im12, time, Lambda, SR_im12, ETA_im12)
! --------------------------------------------------------------------------------- !
RETURN
END SUBROUTINE ETA
! ================================================================================= !
! [The VISCOSITY subroutine is called by the ETA subroutine above, and write2f.f90] !
SUBROUTINE VISCOSITY(U_3,Usb_3,Usb_4,time_tmp,Lambda_tmp,SR,ETA)
! --------------------------------------------------------------------------------- !
DOUBLE PRECISION,INTENT(IN)  :: U_3,Usb_3,Usb_4,time_tmp,Lambda_tmp,SR
DOUBLE PRECISION,INTENT(OUT) :: ETA
! --------------------------------------------------------------------------------- !
DOUBLE PRECISION :: mu,       mu_tot,   mu_thix,  mu_sb3,   mu_sb4,&
                    tau,      tau_tot,  tau_thix, tau_sb3,  tau_sb4

DOUBLE PRECISION :: xi_1cd3,  xi_2cd3,  xi_1sb3,  xi_2sb3,&
                    xi_1sb4,  xi_2sb4,  delta,    Lambda,   time
! --------------------------------------------------------------------------------- !
! "Lambda >= 0" (Lambda_tmp > -0.5D0) means that the time independent calculation   !
! is active, with a time independent shear viscosity ETA. That is, ETA is frozen    !
! to its initial value that applies at time = t = 0 sec. This means that ETA only   !
! depends on radius r; i.e. ETA = ETA(r,t) = ETA(r,0) [ETA is still a function of r,!
! because of its dependency on the shear rate SR = SR(r,t) = SR(r,0)]. In this case,!
! "Lambda" is related to the Continuation Method (see Section 7.8 the Ph.D.-thesis  !
! mentioned in the file "main.f90"). As explained in "main.f90", then U_3, Usb_3,   !
! and Usb_4 are frozen to their initial values during the time independent          !
! calculations. This is implemented in the file "main.f90".                         !
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
! "Lambda = -1" (Lambda_tmp < -0.5D0) means that the time dependent calculation is  !
! active with a time dependent shear viscosity ETA = ETA(r,t).                      !
! --------------------------------------------------------------------------------- !
IF (Lambda_tmp > -0.5D0) THEN        ! => Time independent calculations.            !
  Lambda = Lambda_tmp                !    Lambda = continuation parameter.          !
  time   = 0.0D0                     ! => Time is frozen to zero.                   !
ELSE IF (Lambda_tmp < -0.5D0) THEN   ! => Time dependent calculations.              !
  Lambda = 1.0D0                     !    Lambda = dummy variable.                  !
  time   = time_tmp                  ! => Time is now an increasing value.          !
END IF                               !    [Time is mostly used for older models].   !
! --------------------------------------------------------------------------------- !
delta    = 0.004D0  ! <- The regularization parameter (see Section 7.9).            !
! --------------------------------------------------------------------------------- !
!   cd = "COAGULATION-DISPERSION";           sb = "STRUCTURAL-BREAKDOWN";           !
!    3 = "particle size number 3 (U_3)";      4 = "particle size number 4 (U_4)";   !
! --------------------------------------------------------------------------------- !
xi_1cd3  =  10.70D0     ! (xi_1cd3 = xi_I) 
xi_2cd3  =  35.90D0     ! (xi_2cd3 = xi_II)
! -------------------------------------------
xi_1sb3  =   3.50D0     ! (xi_1sb3 = xi_III)
xi_2sb3  =  14.10D0     ! (xi_2sb3 = xi_IV)
! -------------------------------------------
xi_1sb4  =   2.50D0     ! (xi_1sb4 = xi_V)
xi_2sb4  =  12.10D0     ! (xi_2sb4 = xi_VI)
! -------------------------------------------
mu       =   0.80D0
tau      =  42.50D0
! ---------------------------------------------------
mu_thix  = xi_1cd3*(U_3**(2.0D0/3.0D0))
tau_thix = xi_2cd3*(U_3**(2.0D0/3.0D0))
! -------------------------------------------
mu_sb3   = xi_1sb3*(Usb_3**(2.0D0/3.0D0))
tau_sb3  = xi_2sb3*(Usb_3**(2.0D0/3.0D0))
! -------------------------------------------
mu_sb4   = xi_1sb4*Usb_4**1.0D0
tau_sb4  = xi_2sb4*Usb_4**1.0D0
! -------------------------------------------
mu_tot   =  mu +  mu_thix +  mu_sb3 +  mu_sb4
tau_tot  = tau + tau_thix + tau_sb3 + tau_sb4
! ---------------------------------------------------
ETA      = mu_tot + (tau_tot*Lambda)/(SR + delta)
! --------------------------------------------------------------------------------- !
RETURN
END SUBROUTINE VISCOSITY
! ================================================================================= !
! -*-*-*-*-*-*- [The COAG_DISP_3 subroutine is called by main.f90 ] -*-*-*-*-*-*-*- ! 
! CDF3 = COAGULATION-DISPERSION-FUNCTION for cement particle size number 3:
SUBROUTINE COAG_DISP_3(time,U3_o,U_3,Usb_3,SR,dSR_dt,ddSR_dtdt,CDF3)
! --------------------------------------------------------------------------------- !
DOUBLE PRECISION,INTENT(IN)  :: time,U_3,Usb_3,SR,dSR_dt,ddSR_dtdt
DOUBLE PRECISION,INTENT(OUT) :: U3_o,CDF3
! --------------------------------------------------------------------------------- !
DOUBLE PRECISION :: H_3,I_3,H_param
! --------------------------------------------------------------------------------- !
U3_o = 0.650D0 ! => This value is assigned to U3_o_scalar in main.f90.              !
               !    Note that the maximum value of U3_o + Usb3_o is one, i.e. 1.0D0 !
! --------------------------------------------------------------------------------- !
! Dispersion-rate function I_3:
I_3 = 0.40D0*SR**2.3D0
! ---------------------------------
! Coagulation-rate function H_3:
H_param = 1.0D0*DEXP(-0.0004D0*(dSR_dt - 0.0D0)**2.0D0) + 0.00D0
! - - - - - - - - - - - - - - - - -
H_3 = ((0.40D0*SR**0.01D0 + 0.15D0*SR**0.1D0)*H_param + 1.0D-4)/(SR**2.0D0 + 10.0D0)
! --------------------------------- - - - - - - - - - - - - - - - - - - - - - - - - !
! If linked cement particles that are broken a part, are allowed to participate in
! the process of coagulation, dispersion and re-coagulation, then use the following:
CDF3 = H_3*((1.0D0 - Usb_3) - U_3)**2.0D0 - I_3*U_3**2.0D0
! --------------------------------------------------------------------------------- !
! If linked cement particles that are broken a part, are NOT allowed to participate
! in the process of coagulation, dispersion and re-coagulation, then use the 
! following (if used, remember to define Usb3_o above => DOUBLE PRECISION :: Usb3_o):
! - - - - - - - - - - - - - - - - -
! Usb3_o = 0.350D0 ! => Remember to update this value in "SUBROUTINE STRUCT_BREAK_3".
! CDF3   = H_3*((1.0D0 - Usb3_o) - U_3)**2.0D0 - I_3*U_3**2.0D0
! --------------------------------------------------------------------------------- !
RETURN
END SUBROUTINE COAG_DISP_3
! ================================================================================= !
! -*-*-*-*-*-*- [The STRUCT_BREAK_3 subroutine is called by main.f90] -*-*-*-*-*-*- !
! SBF3 = STRUCTURAL-BREAKDOWN-FUNCTION for cement particle size number 3 (based on
! the concept of structural-breakdown by Tattersall and Banfill):
SUBROUTINE STRUCT_BREAK_3(time,Usb3_o,Usb_3,SR,dSR_dt,ddSR_dtdt,SBF3)
! --------------------------------------------------------------------------------- !
DOUBLE PRECISION,INTENT(IN)  :: time,Usb_3,SR,dSR_dt,ddSR_dtdt
DOUBLE PRECISION,INTENT(OUT) :: Usb3_o,SBF3
! --------------------------------------------------------------------------------- !
DOUBLE PRECISION :: Hsb_3,Isb_3
! --------------------------------------------------------------------------------- !
Usb3_o = 0.350D0 ! The value shown here (i.e. in this particular subroutine) is as- !
                 ! signed to Usb3_o_scalar in main.f90. Note that the maximum value !
                 ! of U3_o + Usb3_o is one, i.e. 1.0D0.                             !
                 ! ---------------------------------------------------------------- !
                 ! If linked cement particles that are broken a part, are NOT al-   !
                 ! lowed to re-coagulate, then remember to update this value in     !
                 ! SUBROUTINE COAG_DISP_3 (above). I.e. this value must be the same !
                 ! as the Usb3_o in the "SUBROUTINE  COAG_DISP_3".                  !
! --------------------------------------------------------------------------------- !
! The "break-apart" function, which describes the rate which the linkages between 
! the size 3 cement particles are broken apart (also designated as dispersion rate). 
! This is the structural-breakdown part for size 3 cement particles.
Isb_3 = 0.045D0*SR**0.4D0
! ---------------------------------
! Extremely slow (or no) reformation of linkages (i.e. hydrate membrane) should
! exist, relative to the experiment of 50 seconds:
Hsb_3 = 0.0D0
! ---------------------------------
SBF3  = - Isb_3*Usb_3**1.0D0
! ---------------------------------
! SBF3 = Hsb_3*(1.0D0 - Usb_3)**2.0D0 - Isb_3*Usb_3**2.0D0
! --------------------------------------------------------------------------------- !
RETURN
END SUBROUTINE STRUCT_BREAK_3
! ================================================================================= !
! -*-*-*-*-*-*- [The STRUCT_BREAK_4 subroutine is called by main.f90] -*-*-*-*-*-*- !
! SBF4 = STRUCTURAL-BREAKDOWN-FUNCTION for cement particle size number 4 (based on
! the concept of structural-breakdown by Tattersall and Banfill):
SUBROUTINE STRUCT_BREAK_4(time,Usb4_o,Usb_4,SR,dSR_dt,ddSR_dtdt,SBF4)
! --------------------------------------------------------------------------------- !
DOUBLE PRECISION,INTENT(IN)  :: time,Usb_4,SR,dSR_dt,ddSR_dtdt
DOUBLE PRECISION,INTENT(OUT) :: Usb4_o,SBF4
! --------------------------------------------------------------------------------- !
DOUBLE PRECISION :: Hsb_4,Isb_4
! --------------------------------------------------------------------------------- !
Usb4_o = 1.0D0 ! => This value is assigned to Usb4_o_scalar in main.f90.
! ---------------------------------
! The "break-apart" function, which describes the rate which the linkages between 
! the size 4 cement particles are broken apart (also designated as dispersion rate). 
! This is the structural-breakdown part for size 4 cement particles.
Isb_4 = 1.20D0*SR**0.1D0
! ---------------------------------
! Extremely slow (or no) reformation of linkages (i.e. hydrate membrane) should
! exist, relative to the experiment of 50 seconds:
Hsb_4 = 0.0D0
! ---------------------------------
SBF4  = - Isb_4*Usb_4**1.0D0
! ---------------------------------
! SBF4 = Hsb_4*(1.0D0 - Usb_4)**2.0D0 - Isb_4*Usb_4**2.0D0
! --------------------------------------------------------------------------------- !
RETURN
END SUBROUTINE STRUCT_BREAK_4
! ================================================================================= !
END MODULE SHEAR_VISCOSITY
! --------------------------------------------------------------------------------- !
