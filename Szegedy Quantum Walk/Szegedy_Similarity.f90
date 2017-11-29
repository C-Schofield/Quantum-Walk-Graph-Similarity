module Szegedy_Similarity
    ! Author: Callum Schofield - 2017
    ! Provides functions which calculate the similarity between two graphs using Szegedy walks
    use Szegedy_Functions
    use Common_Functions
    use m_refsor

    implicit none

    private :: genProb
    public :: Szegedy_Compare, Szegedy_Sim, genProbAll

    contains

    ! Compare a graph with the probability from another graph already determined
    function Szegedy_Compare( A, probOrig, nSteps, errMult) result(sim)
        implicit none
        integer, intent(in)    :: nSteps
        real(8), intent(inout) :: A(:,:), probOrig(:,:)
        real(8), intent(in)    :: errMult
        real(8), allocatable, dimension(:,:) :: U
        real(8), allocatable, dimension(:) :: psi, prob
        real(8) :: err, totalErr(3), sim(3)
        integer :: n, i, j

        n = size(A,1)
        allocate( U(n**2,n**2), psi(n**2), prob(n) )

        call SzegInit( A, psi, U)
        totalErr(1) = 0.d0
        totalErr(2) = 0.d0
        totalErr(3) = 0.d0
        do i = 1, nSteps
            err = 0
            call SzegStep( U, psi, n)
            call MapProb(psi, prob, n)
            call refsor(prob)
            do j = 1, n
                err = err + abs( probOrig(i,j) - prob(j) )
            end do
            ! Absolute difference between probability vectors
            totalErr(1) = totalErr(1) + err
            ! Euclidean distance between probability vectors
            totalErr(2) = totalErr(2) + Norm2(probOrig(i,:) - prob)
            ! Matusita distance between probability vectors
            totalErr(3) = totalErr(3) + Matusita_dist(probOrig(i,:), prob)
        end do 
        do i = 1, 3
            ! Average error over all steps in Szegedy walk
            totalErr(i) = totalErr(i) / nSteps
            ! If errMult is not > 0, return only d. This way various weightings can be applied post completion.
            if( errMult > 0 )then 
                sim(i) = 1.d0 / ( 1.d0 + totalErr(i)*errMult) 
            else 
                sim(i) = totalErr(i)
            end if    
        end do
        deallocate( U, psi, prob )
    end

    ! Find the Similarity between two graphs using Szegedy walk
    function Szegedy_Sim( A, B, nSteps, errMult) result(sim)
        implicit none
        integer, intent(in)    :: nSteps
        real(8), intent(in)    :: errMult
        real(8), intent(inout), dimension(:,:)  :: A, B
        real(8), allocatable, dimension(:,:)    :: U_A, U_B
        real(8), allocatable, dimension(:)      :: psiA, psiB, probA, probB
        real(8) :: err, totalErr(3), sim(3)
        integer :: n, i, j

        n = size(A,1)
        allocate( U_A(n**2,n**2), U_B(n**2,n**2), psiA(n**2), probA(n), psiB(n**2), probB(n) )
        call SzegInit( A, psiA, U_A)
        call SzegInit( B, psiB, U_B)
        totalErr(1) = 0.d0
        totalErr(2) = 0.d0
        totalErr(3) = 0.d0

        do i = 1, nSteps
            err = 0
            call SzegStep( U_A, psiA, n)
            call SzegStep( U_B, psiB, n)
            call MapProb(psiA, probA, n)
            call MapProb(psiB, probB, n)
            call refsor(probA)
            call refsor(probB)

            do j = 1, n
                 err = err + abs( probA(j) - probB(j) )
            end do
            ! Absolute difference between probability vectors
            totalErr(1) = totalErr(1) + err
            ! Euclidean distance between probability vectors
            totalErr(2) = totalErr(2) + Norm2(probA - probB)
            ! Matusita distance between probability vectors
            totalErr(3) = totalErr(3) + Matusita_dist(probA, probB)
        end do
        ! Average error over all steps in Szegedy walk
        do i = 1, 3
            totalErr(i) = totalErr(i)/nSteps
            ! If errMult is not > 0, return only d. This way various weightings can be applied post completion.
            if( errMult > 0 )then 
                sim(i) = 1.d0 / ( 1.d0 + totalErr(i)*errMult) 
            else 
                sim(i) = totalErr(i)
            end if    
        end do
        deallocate( U_A, U_B, psiA, probA, psiB, probB)
    end function

    ! Generates a 1 x n matrix to hold the probabilities of each state for the final step
    subroutine genProb( A, probStates, nSteps)
        implicit none
        real(8), intent(in) :: A(:,:)
        real(8), intent(out) :: probStates(:)
        integer, intent(in) :: nSteps
        real(8), allocatable, dimension(:,:) :: U
        real(8), allocatable, dimension(:) :: psi, prob
        integer :: i ,n

        n = size(A, 1)

        allocate( U(n**2,n**2), psi(n**2), prob(n))
        call SzegInit( A, psi, U)
        
        ! Apply Szegedy step operator
        do i = 1, nSteps
            call SzegStep( U, psi, n)
        end do
        ! Only store the final set of probabilities
        call MapProb(psi, prob, n)
        ! Sort the probability vector
        call refsor(prob)
        probStates(:) = prob
        deallocate( U, psi, prob)
    end subroutine
    
    ! Generates a nSteps x n matrix to hold the probabilities of each state
    subroutine genProbAll( A, probStates, nSteps)
        implicit none
        real(8), intent(in) :: A(:,:)
        real(8), intent(out) :: probStates(:,:)
        integer, intent(in) :: nSteps
        real(8), allocatable, dimension(:,:) :: U
        real(8), allocatable, dimension(:) :: psi, prob
        integer :: i ,n

        n = size(A, 1)

        allocate( U(n**2,n**2), psi(n**2), prob(n))
        call SzegInit( A, psi, U)
        
        do i = 1, nSteps
            ! Apply Szegedy step operator
            call SzegStep( U, psi, n)
            ! Store the probability vector for each time step
            call MapProb(psi, prob, n)
            ! Sort the probability vector
            call refsor(prob)
            probStates(i,:) = prob
        end do
        deallocate( U, psi, prob)
    end subroutine

end module
