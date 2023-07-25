#! /bin/bash
#
gfortran -c -Wall sphere_delaunay.f90
if [ $? -ne 0 ]; then
  echo "sphere_delaunay compile error."
  exit
fi
#
gfortran -c -Wall stripack.f90
if [ $? -ne 0 ]; then
  echo "stripack compile error."
  exit
fi

gfortran sphere_delaunay.o stripack.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm sphere_delaunay.o
#
mv a.out sphere_delaunay_bonds
#
echo "Normal end of execution."
