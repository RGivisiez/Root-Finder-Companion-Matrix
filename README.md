# Root finder algorithm (Companion matrix)

 This algorithm uses a companion matrix to find the roots of a polynomial.
 
 The companion matrix technic changes the problem of finding the polynomial roots to an eigenvalue problem where we have to diagonalize a matrix.

 The program uses subroutines from LAPACK and BLAS to diagonalize the matrice. For Ubuntu users, it is easy to install these packages using the command line below.

  > #### In the terminal:
  >
  > sudo apt-get install liblapack-dev libblas-dev


 To run root_finder_companion_matrix.f90 using gfortran, type the command line below. 

  > #### Compilation:
  >
  > gfortran -llapack -lblas root_finder_companion_matrix.f90 -o root_finder && ./root_finder

## License

See the [LICENSE](LICENSE.md) file for license rights and limitations (MIT).
