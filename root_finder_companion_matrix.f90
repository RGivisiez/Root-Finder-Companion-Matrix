Module root_finder

  Implicit None

  Contains

  Subroutine lapack_root_finder(coefficients, pol_degree, real_part_root, img_part_root)

    Implicit None
    
    Real*8, Allocatable, Dimension(:), intent(out) :: real_part_root, img_part_root
    Real*8, intent(in) :: coefficients(pol_degree + 1)
    Integer*4, intent(in) ::  pol_degree
    
    Integer*4, Parameter :: Nb = 64
    Integer*4 :: LWORK
      
    Integer*4 LDA, Info
    
    Real*8 Dummy(1,1)
    Real*8, Allocatable, Dimension (:,:) :: Companion_matrix
    Real*8, Allocatable, Dimension (:) :: WORK

    Integer*4 j
    Real*8 normalization
   
   Allocate(real_part_root(pol_degree), img_part_root(pol_degree))

    LWORK = (2 + Nb) * 3 * (pol_degree + 1)
    
    LDA = pol_degree

    Allocate(Companion_matrix(LDA, pol_degree), WORK(LWORK))
    
    Companion_matrix = 0.0d0
    
    Do j = 2, pol_degree

       Companion_matrix(j - 1, j) = 1.0d0

    end do  
    
    normalization = coefficients(pol_degree + 1)

    Do j = 1, pol_degree

      Companion_matrix(pol_degree, j) = ( -1.0d0 * coefficients(j) ) / normalization

    end do
    
    ! First the program find the best array size to work with.
    Call DGEEV('No left vectors', 'No right vectors', pol_degree, &
               Companion_matrix, LDA, real_part_root, &
               img_part_root, Dummy, 1, Dummy, 1, WORK, -1, Info)
   
    LWORK = min(int(LWORK), int(WORK(1)));

    Deallocate(WORK)

    Allocate(WORK(LWORK))

    ! Now it diagonalize the matrice.
    Call DGEEV('No left vectors', 'No right vectors', pol_degree, &
               Companion_matrix, LDA, real_part_root, &
               img_part_root, Dummy, 1, Dummy, 1, WORK, LWORK, Info)
    
    if( Info /= 0 )then 
     print*,' '
     write(*,*)'INFO:', Info, LWORK
     write(*,*)'The algorithm had some problems finding the roots.'
     print*,' '
    end if
    
  End Subroutine lapack_root_finder

  Subroutine save_roots_arq(pol_degree, real_part_root, img_part_root, arq_name)
    
    Integer*4, intent(in) :: pol_degree
    Real*8, Dimension (:), intent(in) :: real_part_root, img_part_root
    Character(len=*), intent(in) :: arq_name

    Integer*4 j

    Open(Unit = 100, file = arq_name)

    write(100,'("#",6X,"Real Part",15X,"Imaginary Part")')  
    
    Do j = 1, pol_degree

      If(img_part_root(j) .eq. 0.0d0)then

        write(100,*) real_part_root(j), 0.0d0            ! Write real roots.

      else

        write(100,*) real_part_root(j), img_part_root(j) ! Write imaginary roots.

      end if
    end do

    close(100)

  End Subroutine save_roots_arq

  Subroutine number_of_rows_in_file(arq_name, number_of_rows)
    
    Implicit None

    Character(len=*), intent(in) :: arq_name
    Integer*4, intent(out) :: number_of_rows

    Integer*4 end_of_file, err

    Open(100, file=Trim(arq_name), status='old', iostat=err)

    Call open_file_error(arq_name, err)

    number_of_rows = 0

    Do

      Read(100, *, IOSTAT = end_of_file)

      if( IS_IOSTAT_END(end_of_file) ) exit

      number_of_rows = number_of_rows + 1

    end do

    Close(100)

  End Subroutine number_of_rows_in_file

  Subroutine open_file_error(arq_name, err)

    Implicit None

    Integer*4, intent(in) :: err
    Character(len=*), intent(in) :: arq_name

    if (err /= 0) then
      print*,'File (', Trim(arq_name), ') does not exist.'
      stop
    end if

  End Subroutine open_file_error

End Module root_finder

Program Main

  Implicit None

  Call system('mkdir roots')

  Call simple_polynomial

  Call high_degree_polynomial

End Program Main

Subroutine simple_polynomial

  Use root_finder
  Implicit None

  Real*8, Allocatable, Dimension (:) :: coefficients                    ! Polynomial coefficients.
  Real*8, Allocatable, Dimension (:) :: img_part_root, real_part_root   ! Imaginary and real part of the roots.
  Integer*4 pol_degree                                                  ! Polynomial degree.

  Character*60 arq_name   ! File name where roots values will be saved.

  !Example: 4 - x^2 = 0

  pol_degree = 2        ! Polynomial degree

  Allocate(coefficients(pol_degree + 1)) ! The +1 is to count the constant term c. (c + x + x^2)

  coefficients(1) = 4    ! Constant (c)
  coefficients(2) = 0    ! First term (x^1)
  coefficients(3) = -1   ! Second term (x^2)

  !Root finder subroutine

  Call lapack_root_finder(coefficients, pol_degree, real_part_root, img_part_root)

  arq_name = 'roots/simple_roots.dat'

  Call save_roots_arq(pol_degree, real_part_root, img_part_root, arq_name)

  print*, ''
  print*, 'Solving: 4 - x^2 = 0'
  print*, 'Roots:', real_part_root(1), img_part_root(1)
  print*, 'Roots:', real_part_root(2), img_part_root(2)
  print*, ''

End Subroutine simple_polynomial

Subroutine high_degree_polynomial

  Use root_finder
  Implicit None

  Real*8, Allocatable, Dimension (:) :: coefficients                    ! Polynomial coefficients.
  Real*8, Allocatable, Dimension (:) :: img_part_root, real_part_root   ! Imaginary and real part of the roots.
  Integer*4 pol_degree, number_of_rows, i

  Character(len=60) arq_name

  ! File name where the coefficients are stored.
  arq_name = 'HD_coefficients/HD_coefficients.dat'

  Call number_of_rows_in_file(arq_name, number_of_rows)

  pol_degree = number_of_rows - 1

  Allocate(coefficients(pol_degree + 1)) ! The +1 is to count the constant term c. (c + x + x^2)

  Open(100, file=arq_name)

  Do i = 1, number_of_rows
    read(100, *) coefficients(i)
  end do 

  print*, ''
  print*, 'Solving a high degree polynomial.'

  Call lapack_root_finder(coefficients, pol_degree, real_part_root, img_part_root)

  arq_name = 'roots/HD_polynomial_roots.dat'

  Call save_roots_arq(pol_degree, real_part_root, img_part_root, arq_name)

  print*,
  print*,'Done.'
  print*,
  
  print*,'See the roots in file: ', Trim(arq_name)
  print*, ''

End Subroutine high_degree_polynomial
