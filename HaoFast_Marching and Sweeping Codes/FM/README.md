Some notes on running fast_marching and fast_sweeping
=====================================================

1. Both fast_marching.f90 and fast_sweeping.f90 can be compiled using gfortran 
   or ifort, for example:

   gfortran -cpp fast_sweeping.f90 -Ofast -o fs

2. input file is very simple, i.e. a bunch of 2D coordinates. First line needs 
   to be: # xxx, where xxx is the number points. See b2.dat for example.

3. Test_cases_pnt.tgz contains other files with points and elements. Note that
   only points are needed for fm or fs.

4. In the FM folder, b3.dat and b4.dat are now added.

5. Output files are: fm.g and fm.q in PLOT3D format and bcurve.dat in Tecplot
   format.

17-Oct-23
Hao
