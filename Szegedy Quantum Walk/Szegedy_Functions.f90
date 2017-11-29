module Szegedy_Functions
    ! Author: Callum Schofield - 2017
    ! Subroutines used to model the Szegedy Quantum Walk
    use Common_Functions
    
    implicit none

    private :: Uop, swapOp, colNormalise
    public :: SzegStep, SzegInit, kron_1D, kron_2D, MapProb

    contains

    ! Applies a single step of the Szegedy quantum walk
    ! U - Szegedy step operator 
    ! psi - system state 
    ! n - size of graph
    subroutine SzegStep( U, psi, n)
        implicit none
        real(8), dimension(:,:), intent(in) :: U
        real(8), dimension(:), intent(inout) :: psi
        real(8), allocatable, dimension(:) :: tmp
        integer, intent(in) :: n

        allocate( tmp(n**2))
        call dgemv( 'N', n**2, n**2, 1.d0, U, n**2, psi, 1, 0.d0, tmp, 1)
        psi = tmp
        deallocate(tmp)

    end subroutine

    ! Generate the U operator and psi
    subroutine SzegInit( A, psi, U)
        implicit none
        real(8), dimension(:,:), intent(in) :: A
        real(8), dimension(:,:), intent(out) :: U
        real(8), dimension(:), intent(out) :: psi
        real(8), allocatable, dimension(:,:) :: proj, states, P
        integer :: n, i

        n = size(A,1)

        allocate( proj(n**2, n**2), states(n,n**2), P(n,n))
        P = A

        ! Ensure input matrix is stochastic / column normalised
        ! Will turn an adjacency matrix into a stochastic transition matrix
        call colNormalise( P)
        call genStates( P, states)
        ! Calculate the projection operator
        call projOp( states, proj, n)
        ! Calculate the Szegedy step operator
        call UOp( proj, U, n)
        ! Determine psi0
        do i = 1, n**2
            psi(i) = sum( states(:,i)) / sqrt(real(n))
        end do
        deallocate(proj, states, P)
    end subroutine

    ! Normalise the columns of input  matrix
    subroutine colNormalise( P)
        implicit none
        real(8), dimension(:,:), intent(inout) :: P
        real(8) :: total
        integer :: n, i

        n = size(P, 1)
        ! If no entries in column, set index (i,i) to 1
        do i = 1, n
            total = sum( P(:,i)) 
            if(total == 0) then
                P(i,i) = 1
            else 
                P(:,i) = P(:,i) / total
            end if
            
        end do

    end subroutine

    ! Compute the U operator
    subroutine UOp( proj, U, n)
        implicit none
        real(8), dimension(:,:), intent(in) :: proj
        real(8), dimension(:,:), intent(out) :: U
        real(8), allocatable, dimension(:,:) :: Id, tmp, swap
        integer, intent(in) :: n
        
        allocate(Id(n**2,n**2), swap(n**2,n**2), tmp(n**2,n**2))

        call identity( Id, n**2)
        call swapOp( swap, n)

        tmp = 2*proj - Id

        call dgemm( 'N', 'N', n**2, n**2, n**2, 1.d0, swap, n**2, tmp, n**2, 0.d0, U, n**2)

        deallocate(Id, tmp, swap)

    end subroutine

    subroutine genStates( P, states)
        implicit none
        real(8), allocatable, dimension(:) :: ketJ, psiJ
        real(8), dimension(:,:), intent(in) :: P
        real(8), dimension(:,:), intent(out) :: states
        integer :: n, i

        n = size(P,1)

        ! Psi is set of n x n**2 psiJ states
        allocate( psiJ(n**2), ketJ(n) )

        do i = 1, n
            ketJ = 0
            ketJ(i) = 1
            call kron_1D( ketJ, sqrt(P(:,i)), psiJ)
            states(i,:) = psiJ
        end do

        deallocate( psiJ, ketJ )
    end subroutine

    ! Compute n**2 by n**2 the swap operator
    subroutine swapOp( swap, n)
        implicit none
        integer :: j, k
        integer, intent(in) :: n
        real(8), dimension(:,:), intent(out) :: swap

        swap(:,:) = 0
        do j = 1, n
            do k = 1, n
                swap( (j-1)*n + k, (k-1)*n + j) = 1
            end do                 
        end do
    end subroutine

    ! Compute the projection operator
    subroutine projOp( states, proj, n)
        implicit none
        real(8), dimension(:,:), intent(in) :: states
        real(8), dimension(:,:), intent(out) :: proj
        integer, intent(in) :: n

        ! states is n*n**2
        ! Computes states^T.states
        ! Note states matrix is defined by rows
        call dgemm('T', 'N', n**2, n**2, n, 1.d0, states, n, states, n, 0.d0, proj, n**2)

    end subroutine

    ! Maps the state vector to the probability vector
    subroutine MapProb( psi, prob, n)
        implicit none
        real(8), dimension(:), intent(in) :: psi
        real(8), dimension(:), intent(out) :: prob
        real(8), allocatable, dimension(:) :: tmp
        integer, intent(in) :: n
        integer :: i

        allocate( tmp(n**2))

        tmp = abs(psi)**2
        
        do i = 1, n
            ! Sum n lots of elements
            prob(i) = sum( tmp(n*(i-1)+1:n*i))
        end do

        deallocate(tmp)

    end subroutine MapProb

    ! Determine the kronecker product of two 1D matrices 
    ! Input A, B 
    ! Output C
    subroutine kron_1D( A, B, C)
        implicit none 
        real(8), dimension(:), intent(in) :: A, B
        real(8), allocatable, dimension(:), intent(out) :: C
        integer :: nA, nB, i, j, x
        ! A is 1 x nA matrix
        nA = size(A)
        ! B is 1 x nB matrix
        nB = size(B)

        if(allocated(C)) deallocate(C)
        allocate( C(nA*nB) )

        ! C is a mA*mB x nA*nB matrix
        do i = 1, nA
            do j = 1, nB
                x = nB*(i-1) + j
                C(x) = A(i)*B(j)
            end do    
        end do
    end subroutine

    ! Determine the kronecker product of two 2D matrices 
    ! Input A, B 
    ! Output C
    subroutine kron_2D( A, B, C)
        implicit none  
        real(8), dimension(:,:), intent(in) :: A, B
        real(8), allocatable, dimension(:,:), intent(out) :: C
        integer :: mA, nA, mB, nB, i, j, k, l, x, y
        ! A is mA x nA matrix
        mA = size(A,1)
        nA = size(A,2)
        ! B is mB x nB matrix
        mB = size(B,1)
        nB = size(B,2)

        if(allocated(C)) deallocate(C)
        allocate( C(mA*mB, nA*nB) )

        ! C is a mA*mB x nA*nB matrix
        do i = 1, mA
            do j = 1, nA
                do k = 1, mB
                    do l = 1, nB
                        x = mB*(i-1) + k
                        y = nB*(j-1) + l
                        C(x,y) = A(i,j)*B(k,l)
                    end do    
                end do
            end do     
        end do
    end subroutine

end module Szegedy_Functions