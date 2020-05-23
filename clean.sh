#!/bin/sh
# -------------------------------------------------
# rm *.mod
# rm *.dat
# rm first
# -------------------------------------------------
echo "---------------------------------------------"
# =================================================
select item in clean_mod clea_dat l_count quit
do
if [ $item = "clean_mod" ]; then
  if [ -s matrix.mod ]; then
    echo "Removing mod files..."
    echo "   "
    rm -v *.mod
    echo "   "
    echo " ...done!"
  else
    echo "The .mod files are not present (nothing to remove!)"
  fi
  break # haettir i thessa do luppu.
elif [ $item = "clea_dat" ]; then
  if [ -s log.dat ]; then
    echo "Removing all data files (except for input files)..."
    echo "   "
    mv omega_input.dat omega_input.tmp
    mv time_input.dat time_input.tmp
    rm -v *.dat
    mv omega_input.tmp omega_input.dat
    mv time_input.tmp time_input.dat
    echo "   "
    echo " ...done!"
  else
    echo "The data files are not present (nothing to remove!)"
  fi
  break # haettir i thessa do luppu.
elif [ $item = "l_count" ]; then
  echo "---------------------------------------------"
  echo "Line count of the sources:"
  echo "   "
  echo `wc -l main.f90`
  echo `wc -l motion.f90`
  echo `wc -l shear.f90`
  echo `wc -l viscous.f90`
  echo `wc -l param.f90`
  echo `wc -l update.f90`
  echo `wc -l write2f.f90` 
  echo "   "
  echo "Quitting (nothing is changed!)"
  echo "---------------------------------------------"
  exit  # haettir alfarid i thessu forriti.
elif [ $item = "quit" ]; then
  echo "---------------------------------------------"
  echo "Quitting (nothing is changed!)"
  echo "---------------------------------------------"
  exit  # haettir alfarid i thessu forriti.
fi
done
echo "---------------------------------------------"
# =================================================
