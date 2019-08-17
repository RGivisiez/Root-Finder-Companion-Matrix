Program Main

  Implicit None

  Real*8, Allocatable, Dimension (:) :: coefficients                    ! Polynomial coefficients.
  Real*8, Allocatable, Dimension (:) :: img_part_root, real_part_root   ! Imaginary and real part of the roots.
  Integer*4 pol_degree                                                  ! Polynomial degree.

  Character*60 arq_name   ! File name where roots values will be saved.

  Allocate(coefficients(pol_degree + 1)) ! The +1 is to count the constant term c. (c + x + x^2)

  Allocate(img_part_root(pol_degree + 1), real_part_root(pol_degree + 1)) ! The +1 is to count the constant term c. (c + x + x^2)

  real_part_root = 0
  img_part_root = 0

  !Example: 4 - x^2 = 0

  pol_degree = 2        ! Polynomial degree

  coefficients(1) = 4    ! Constant (c)
  coefficients(2) = 0    ! First term (x^1)
  coefficients(3) = -1   ! Second term (x^2)

  !Root finder subroutine

  Call lapack_root_finder(coefficients, pol_degree, real_part_root, img_part_root)

  arq_name = 'roots.dat'

  Call save_roots_arq(pol_degree, real_part_root, img_part_root, arq_name)

  print*, 'Real part of the roots:', real_part_root(1:pol_degree)
  print*, 'Imaginary part of the roots:', img_part_root(1:pol_degree)

End Program Main

Subroutine lapack_root_finder(coefficients, pol_degree, real_part_root, img_part_root)

  Implicit None
  
  Integer*4, Parameter :: Nb = 64
  Integer*4 :: LWORK
    
  Integer*4 LDA, Info
  
  Real*8 Dummy(1,1)
  Real*8, Allocatable, Dimension (:,:) :: A
  Real*8, Allocatable, Dimension (:) :: WORK

  Integer*4 j, pol_degree
  Real*8 normalization
  Real*8 real_part_root(pol_degree + 1), img_part_root(pol_degree + 1)
  Real*8 coefficients(pol_degree + 1)
 
  LWORK = (2 + Nb) * 3 * (pol_degree + 1)
  
  LDA = pol_degree

  Allocate(A(LDA, pol_degree), WORK(LWORK))
  
  A = 0.0d0
  
  Do j = 2, pol_degree

     A(j - 1, j) = 1.0d0

  end do  
  
  normalization = coefficients(pol_degree + 1)

  Do j = 1, pol_degree

    A(pol_degree, j) = ( -1.0d0 * coefficients(j) ) / normalization

  end do
  
  Call DGEEV('No left vectors', 'No right vectors', pol_degree, A, LDA, real_part_root, &
              & img_part_root, Dummy, 1, Dummy, 1, WORK, -1, Info)
 
  LWORK = min(int(LWORK),INT(WORK(1)));

  Deallocate(WORK)
  Allocate(WORK(LWORK))

  Call DGEEV('No left vectors', 'No right vectors', pol_degree, A, LDA, real_part_root, &
              & img_part_root,Dummy, 1, Dummy, 1, WORK, LWORK, Info)
  
  if( Info /= 0 )then 
   print*,' '
   write(*,*)'INFO:', Info, LWORK
   write(*,*)'The algorithm had some problems finding the roots.'
   print*,' '
  end if
  
End Subroutine lapack_root_finder

Subroutine save_roots_arq(pol_degree, real_part_root, img_part_root, arq_name)
  
  Integer*4 pol_degree
  Real*8 real_part_root(pol_degree + 1), img_part_root(pol_degree + 1)

  Character*60 arq_name

  Open(Unit = 99, file = arq_name)
  
  write(99,'("#",6X,"Real Part",15X,"Imaginary Part")')  
  
  Do j = 1, pol_degree

    If(img_part_root(j) .eq. 0.0d0)then

      write(99,*) real_part_root(j), 0.0d0            ! Real roots.

    else

      write(99,*) real_part_root(j), img_part_root(j) ! Imaginary roots.

    end if
  end do

  close(99)

End Subroutine save_roots_arq