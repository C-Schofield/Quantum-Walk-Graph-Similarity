module NodeAffinity_Functions
    ! Author: Callum Schofield - 2017
    ! Functions to determine similarity using the DeltaCon algorithm developed by Koutra, Vogelstein and Faloutsos
    use m_refsor
    use Common_Functions

    implicit none

    public NodeAffinity_Sim, NodeAffinity_Compare, belief

    ! Interface to allow double and integer matrices
    interface NodeAffinity_Sim 
        module procedure D_NodeAffinity_Sim, I_NodeAffinity_Sim 
    end interface

    contains

    ! Takes two double adjacency matrices to determine the similarity index using the node affinity method
    function D_NodeAffinity_Sim( A, B) result(sim)
        implicit none
        real(8) :: d, sim
        integer :: i, j, n
        real(8), dimension(:,:), intent(inout) :: A, B
        real(8), dimension(:,:), allocatable :: S1, S2

        d = 0
        n = size(A, 1)

        allocate(S1(n,n), S2(n,n))

        ! Find the node affinities
        call belief(A, S1)
        call belief(B, S2)

        ! Determine the difference between the node affinities
        do i = 1, n 
            do j = 1, n
                d = d + (sqrt( abs(S1(i,j)) ) - sqrt( abs(S2(i,j)) ) )**2
            end do 
        end do

        sim = 1.d0 / (1.d0 + sqrt(d))

        deallocate(S1, S2)

    end function

    ! Takes two integer adjacency matrices to determine the similarity index using the node affinity method
    function I_NodeAffinity_Sim( A, B) result(sim)
        implicit none
        real(8) :: d, sim
        integer :: i, j, n
        integer, dimension(:,:), intent(in) :: A, B
        real(8), dimension(:,:), allocatable :: S1, S2, A_real, B_real
        d = 0
        n = size(A, 1)

        allocate(S1(n,n), S2(n,n), A_real(n,n), B_real(n,n))
        ! belief takes real matrices for A,B
        A_real = real(A,8)
        B_real = real(B,8)

        ! Find the node affinities
        call belief(A_real, S1)
        call belief(B_real, S2)

        ! Determine the difference between the node affinities
        do i = 1, n 
            do j = 1, n
                d = d + (sqrt( abs(S1(i,j)) ) - sqrt( abs(S2(i,j)) ) )**2
            end do 
        end do

        sim = 1.d0 / (1.d0 + sqrt(d))

        deallocate(S1, S2, A_real, B_real)

    end function

    ! Calculates the similarity score using two already calculated node affinity matrices S1, S2
    function NodeAffinity_Compare( S1, S2) result(sim)
        real(8) :: d, sim
        integer :: i, j, n 
        real(8), dimension(:,:) :: S1, S2

        d = 0
        n = size(S1, 1)

        do i = 1, n 
            do j = 1, n
                d = d + (sqrt( abs(S1(i,j)) ) - sqrt( abs(S2(i,j)) ) )**2
            end do 
        end do

        sim = 1.d0 / (1.d0 + sqrt(d))
    end function

    ! Calculates node affinity matrix S for input adjacency matrix A
    subroutine belief( A, S)
        implicit none
        real(8), dimension(:,:), intent(inout) :: A
        real(8), dimension(:,:), intent(out) :: S
        real(8), dimension(:,:), allocatable :: Id, D, WORK(:)
        real(8) :: eps, total
        integer :: i, j, n, info, LWORK
        integer, allocatable :: IPIV(:)
        
        n = size(A,1)
        do i = 1, n
            total = sum( A(:,i)) 
            !if(total == 0) then
            if(total <= 1.d-10) then
                A(i,i) = 1          ! If no connections, self loop
            end if
        end do
        
        LWORK = n**2
        allocate(Id(n,n), D(n,n), IPIV(n), WORK(LWORK))
        call identity(Id, n)
        D = 0

        do i = 1, n 
            do j = 1, n
                D(i,i) = D(i,i) + A(j, i) ! Sum over rows in column i
            end do 
        end do

        eps = 1 / (1 + maxval(D) )

        S = Id + eps**2*D - eps*A

        ! Invert S
        call dgetrf( n, n, S, n, IPIV, info)
        if(info /= 0) print*,'Error factorising matrix for inversion. Info = ', info
        call dgetri( n, S, n, IPIV, WORK, LWORK, info)
        if(info /= 0) print*,'Error inverting matrix. Info = ', info
        deallocate(Id, D, IPIV, WORK)
    end subroutine

end module