# Requirements
===============
1. The C++ BOOST library from http://boost.org -- particularly boost/foreach.hpp

If on a mac, you can install these using 'macports' or 'homebrew'.


# Compiling
===============
Run the command 'recompile', which generates an executable 'consciousness'.  I.e.,

$ ./recompile
$ ./consciousness filename.txt

Computes the measures for the system specified in filename.txt





# Some notes
=======================

The directory "e/" contains example systems to compute the phi/psi.  For example, to compute the measures for the system "transitions/3RN.txt", you'd do:

$ ./consciousness e/transitions/3RN.txt

-----

The directory 'balduzzi_python' contains the original python code from David Balduzzi to compute the phi in the 2008 paper, "Integrated Information in Discrete Dynamical Systems"

-----

The directory 'tests' is a series of simple programs that spit out system diagnostic information.  It's unlikely you'll ever need them.  You can safely ignore this directory.

-----

The directory 'pics' contains two pretty pictures for the Google Code site.  You can ignore this directory.

-----

