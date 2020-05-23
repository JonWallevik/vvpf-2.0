# Makefile for vvpf 2.0.
# Written by Dr. Jon E. Wallevik (last updated 17.12.07).
# Example of usage: "make gnu_opt", or just "make".
# --------------------------------------------------------
GNUOPT = -O3 -funroll-loops -D__NO_MATH_INLINES -ffast-math -march=nocona -ftree-vectorize
INTOPT = -O3 -xP -static -align -Zp16 -ipo
SOURCE = param.f90 shear.f90 viscous.f90 update.f90 motion.f90 write2f.f90 main.f90
NAME   = -o vvpf20
# --------------------------------------------------------
gnu_opt:
	gfortran $(GNUOPT) $(NAME) $(SOURCE)
g95_no_opt:
	g95 $(NAME) $(SOURCE)
intel_no_opt:
	ifort $(NAME) $(SOURCE)
intel_opt:
	ifort $(INTOPT) $(NAME) $(SOURCE)
# --------------------------------------------------------
