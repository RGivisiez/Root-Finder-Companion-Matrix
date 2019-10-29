# Root finder algorithm (Companion matrix)

 This algorithm uses a companion matrix to find the roots of a polynomial.
 
 In resume, the companion matrix technic changes the problem of finding polynomial roots to a problem of matrix
diagonalization, where the eigenvalues of this matrix are the polynomial roots. To diagonalize the matrix I used
subroutines from LAPACK and BLAS, so it is necessary to have these installed on your computer. For Ubuntu users,
it is easy to install those packages using the command line below.

  #### Terminal:
 ```bash
 $ sudo apt-get install liblapack-dev libblas-dev
 ```


 To run root_finder_companion_matrix.f90 using gfortran, type the command line below. 

  #### Compilation:
  ```bash 
  $ gfortran -llapack -lblas root_finder_companion_matrix.f90 -o root_finder && ./root_finder
  ```
## License

See the [LICENSE](LICENSE.md) file for license rights and limitations (MIT).
