# A root finder algorithm (companion matrix)

  This algorithm uses a companion matrix to find roots of a polynomial.

  It is necessary to install lapack and blas to run the algorithm. For Ubuntu users, just type
  the command line below in the terminal,

  > #### In the terminal:
  >
  > sudo apt-get install liblapack-dev libblas-dev


  To compile the root_finder_companion_matrix.f90

  > #### Compilation:
  >
  > gfortran -llapack -lblas root_finder_companion_matrix.f90 -o root_finder && ./root_finder
