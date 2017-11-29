module Common_Functions
    ! Author: Callum Schofield - 2017
    ! Collection of common functionality
    implicit none

    public init_random_seed, timestamp, permutation, mat_inv, file_dimensions, identity, Matusita_dist

    contains

    ! Initialises a random seed using the system clock
    subroutine init_random_seed()

      integer :: i, n, clock
      integer, dimension(:), allocatable :: seed

      call random_seed(size = n)
      allocate(seed(n))

      call system_clock(COUNT=clock)

      seed = clock + 37 * (/ (i - 1, i = 1, n) /)
      call random_seed(PUT = seed)

      deallocate(seed)
    end subroutine

    ! Returns a string which contains the current time
    function timestamp() result(str)
      implicit none
      character(len=20) :: str
      integer :: values(8)
      character(len=4) :: year
      character(len=2) :: month
      character(len=2) :: day, hour, minute, second
      character(len=5) :: zone

      ! Get current time
      call date_and_time(VALUES=values, ZONE=zone)  
      write(year,'(i4.4)')    values(1)
      write(month,'(i2.2)')   values(2)
      write(day,'(i2.2)')     values(3)
      write(hour,'(i2.2)')    values(5)
      write(minute,'(i2.2)')  values(6)
      write(second,'(i2.2)')  values(7)

      str = year//'-'//month//'-'//day//'_'//hour//':'//minute//':'//second
    end function

    ! Returns the identity matrix of size n
    subroutine identity(I, n)
        implicit none
        real(8), dimension(:,:), intent(out) :: I
        integer, intent(in) :: n
        integer :: j

        I(:,:) = 0
        forall(j = 1:n) I(j,j) = 1
        
    end subroutine

    ! Returns a vector which contains a random permutation
    subroutine permutation( n, perm)
      integer, intent(in) :: n
      integer, dimension(:), intent(out) :: perm
      integer :: i, j, tmp
      real(8) :: rand

      ! Initialise the permutation vector 
      perm = (/ (i, i=1,n) /)

      ! Choose an index value between 1, i. Replace the value at current index with this value, repeat.
      do i = n, 2, -1
        call random_number(rand)
        j = floor( i*rand)+1
        tmp = perm(j)
        perm(j) = perm(i)
        perm(i) = tmp
      end do

    end

    ! Returns the inverse of the input matrix
    subroutine mat_inv( A, Ainv)
      real(8), intent(in), dimension(:,:) :: A
      real(8), intent(inout), dimension (:,:) :: Ainv
      real(8), allocatable :: WORK(:)
      integer :: n , info, LWORK
      integer, allocatable :: IPIV(:) 

      n = size(A,1)
      Ainv = real( A)
      LWORK = n**2
      allocate(IPIV(n), WORK(LWORK))
      call dgetrf( n, n, Ainv, n, IPIV, info)
      if(info /= 0) print*,'Error factorising matrix for inversion. Info = ', info
      call dgetri( n, Ainv, n, IPIV, WORK, LWORK, info)
      if(info /= 0) print*,'Error inverting matrix. Info = ', info
      deallocate(IPIV, WORK)

    end
    
    ! Returns the number of lines in the input file
    ! id_in - option to include an id to ensure unique file id
    subroutine file_dimensions(file, n, id_in)
        implicit none
        integer, intent(out) :: n 
        character(len=*), intent(in) :: file
        integer :: readState, id
        integer, optional :: id_in

        id = 1000
        if(present(id_in)) id = id_in

        open(id, file=file, action='read')
        ! Note that I have found some architecture requires n to start at 1 to return correct result
        n = 0
        ! n = 1
        ! Loop until end of file
        do 
            read(id,*, iostat=readState)
            ! If iostat is <0 then end of file reached
            if( readState < 0) exit
            n = n + 1
        end do
        close(id)
    end

    ! Returns Matusita distance between two vectors A and B
    function Matusita_dist( A, B) result(d)
      implicit none
      integer :: i, n
      real(8) :: d
      real(8), intent(inout), dimension(:) :: A, B

      d = 0
      ! Ensure real numbers
      A = abs(A)
      B = abs(B)
      n = size(A)
      do i = 1, n
          d = d + (sqrt(A(i)) - sqrt(B(i)) )**2
      end do
      d = sqrt(d)

    end 

end module