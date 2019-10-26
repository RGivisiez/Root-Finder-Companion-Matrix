# Root finder algorithm (Companion matrix)

  This algorithm uses a companion matrix to find roots of a polynomial. What the algorithm
does is change the problem of find roots to a problem of diagonalize a matrice. The program
uses subroutines from LAPACK and BLAS to diagonalize the matrice, for Ubuntu users
it is easy to install these packages using the command line below.

  > #### In the terminal:
  >
  > sudo apt-get install liblapack-dev libblas-dev


 To run root_finder_companion_matrix.f90 using gfortran, type the command line below. 

  > #### Compilation:
  >
  > gfortran -llapack -lblas root_finder_companion_matrix.f90 -o root_finder && ./root_finder
