module QW_Similarity
    ! Author: Callum Schofield - 2017
    ! Determine the similarity value for a comparison between 2 graphs, A and B 
    ! sim output:
    !   sim(1) - prob1Norm + prob2Norm > eps = 1
    !   sim(2) - prob1Norm/nSteps + prob2Norm/nSteps
    !   sim(3) - (prob1Norm/nSteps)^2 + (prob2Norm/nSteps)^2

    use Walk_Functions
    use data_types
    use Common_Functions
    use omp_lib

    implicit none
    public  QW_Calc_Sim
    private Sum_S_2node, Sum_S_1node, MatrixPermutation, Apply_Walk

    contains

    subroutine QW_Calc_Sim( AdjA, AdjB, sim, nPhase, nSteps, phi, eps)
        integer, intent(inout) :: AdjA(:,:), AdjB(:,:)
        integer, allocatable, dimension(:,:) :: AdjA2, AdjB2, Atmp, Btmp, perm(:)
        real(8), intent(out) :: sim(:)
        real(8), intent(inout) :: phi(:,:), eps
        integer, intent(inout) :: nPhase, nSteps
        real(8) :: d, sumAB(4), sumAA(4), sumBB(4)
        integer :: i, n

        n = size(AdjA,1)

        allocate( AdjA2(n,n), AdjB2(n,n), Atmp(n,n), Btmp(n,n), perm(n))
        ! Permutation vector used to permute input graphs (same vector used for each input)
        call permutation(n, perm)
        call MatrixPermutation(AdjA, AdjA2, perm)
        call MatrixPermutation(AdjB, AdjB2, perm)

        ! Store input graphs, to ensure they remain untampered 
        Atmp = AdjA
        Btmp = AdjB

        deallocate(perm)

        ! Determine the comparison score
        ! Note that the code allows for 0 or 1 reference node phase addition schemes, but these do not successfully determine isomorphisms
        if( nPhase == 2) then
            ! 2 reference nodes
            call Sum_S_2node(AdjA, AdjA2, sumAA, nPhase, nSteps, phi, eps)
            call Sum_S_2node(AdjB, AdjB2, sumBB, nPhase, nSteps, phi, eps)
            AdjA = Atmp
            AdjB = Btmp
            call Sum_S_2node(AdjA, AdjB, sumAB, nPhase, nSteps, phi, eps)
        else
            ! 0 or 1 reference nodes, note phi = 0 for 0 reference nodes
            call Sum_S_1node(AdjA, AdjA2, sumAA, nPhase, nSteps, phi, eps)
            call Sum_S_1node(AdjB, AdjB2, sumBB, nPhase, nSteps, phi, eps)
            AdjA = Atmp
            AdjB = Btmp
            call Sum_S_1node(AdjA, AdjB, sumAB, nPhase, nSteps, phi, eps)
        end if

        deallocate( AdjA2, AdjB2, Atmp, Btmp)

        ! Calculate the similarity score
        do i = 1, size(sumAA)
            if( sumAA(i) > 0) then 
                d = Abs(sumAA(i) - sumAB(i))/sumAA(i)
            else 
                d = sumAB(i)
            end if

            if( sumBB(i) > 0) then 
                d = d + Abs(sumBB(i) - sumAB(i))/sumBB(i)
            else 
                d = d + sumAB(i)
            end if

            sim(i) = 1 / (1 + d)
        end do
    end 

    ! Find the sum of the matrix S using 2 reference nodes
    ! n - number of nodes in graph
    ! AdjA, AdjB - A and B adjacency graph matrices
    ! sum - output sum of matrix SAB
    subroutine Sum_S_2node( AdjA, AdjB, result_sum, nPhase, nSteps, phi, eps)
        type(node), allocatable, dimension(:) :: state_Aorig, state_A, state_Atmp, &
                                               & state_Borig, state_B, state_Btmp
        integer, dimension(:,:), intent(inout) :: AdjA, AdjB
        integer, intent(in) :: nPhase, nSteps
        ! integer, intent(out) :: result_sum
        real(8), intent(inout) :: phi(:,:), eps
        real(8), intent(out) :: result_sum(:)
        integer :: i, j, k, l, m, n
        real(8) :: prob_A(nPhase, nSteps), prob_B(nPhase, nSteps), d(3), tmp

        n = size(AdjA,1)

        allocate( state_Aorig(n), state_Borig(n) )
        call init_state( AdjA, state_Aorig)
        call init_state( AdjB, state_Borig)
        result_sum = 0      
        
        ! -- Shared in loop --
        ! state_Xorig   - original state of graph X (A,B)
        ! -- Private to each thread -- 
        ! state_X       - current state of graph X (A,B)
        ! state_Xtmp    - temporary flip state of current thread iteration for graph X (A,B)
    
        !$omp parallel private(state_A, state_B, state_Atmp, state_Btmp, prob_A, prob_B, d, tmp) &
        !$omp&         shared(AdjA, AdjB, state_Aorig, state_Borig, nPhase, nSteps, phi, n, eps)
            allocate( state_A(n), state_Atmp(n), state_B(n), state_Btmp(n)) 
            ! Use copy state to allocate temporary state variables. Only need to allocate space for flip operation
            call copy_state( state_Aorig, state_Atmp, 1)    
            call copy_state( state_Borig, state_Btmp, 1)
            !$omp do reduction(+: result_sum)
            do i = 1, n - 1
                do j = i + 1, n
                    call copy_state( state_Aorig, state_A)
                    !call copy_state( state_Aorig, state_Atmp)
                    call Apply_Walk( state_A, state_Atmp, AdjA, prob_A, (/ i, j /), nPhase, nSteps, phi)
                    do k = 1, n
                        do l = 1, n
                            call copy_state( state_Borig, state_B)
                            !call copy_state( state_Borig, state_Btmp)
                            call Apply_Walk( state_B, state_Btmp, AdjB, prob_B, (/ k, l /), nPhase, nSteps, phi)
                            d = 0
                            do m = 1, nPhase
                                ! Find the euclidean distance between the probability vectors
                                tmp = Norm2(prob_A(m,:) - prob_B(m,:))/nSteps
                                d(1) = d(1) + tmp
                                d(2) = d(2) + tmp**2
                                d(3) = Matusita_dist(prob_A(m,:), prob_B(m,:)) / nSteps
                            end do
                            ! Add sum by 1 if probability matches, based on threshold
                            if( d(1)/nPhase > eps) then           
                                result_sum(1) = result_sum(1) + 1
                            end if
                            do m = 1, 3
                                result_sum(m+1) = result_sum(m+1) + d(m)
                            end do
                        end do    
                    end do
                end do 
            end do
            !$omp end do
            call clear_state(state_A)
            call clear_state(state_B)
            call clear_state(state_Atmp)
            call clear_state(state_Btmp)
            deallocate( state_A, state_Atmp, state_B, state_Btmp)
        !$omp end parallel

        ! Clear state variables
        call clear_state(state_Aorig)
        call clear_state(state_Borig)
        !deallocate(state_Aorig, state_A, state_Atmp, state_A_Flip, state_Borig, state_B, state_Btmp, state_B_Flip)
        deallocate(state_Aorig, state_Borig)

    end

    ! Find the sum of the matrix S using 1 reference node
    ! n - number of nodes in graph
    ! AdjA, AdjB - A and B adjacency graph matrices
    ! sum - output sum of matrix SAB
    subroutine Sum_S_1node( AdjA, AdjB, result_sum, nPhase, nSteps, phi, eps)
        type(node), allocatable, dimension(:) :: state_Aorig, state_A, state_Atmp, state_Borig, state_B, state_Btmp
        integer, dimension(:,:), intent(inout) :: AdjA, AdjB
        integer, intent(in) :: nPhase, nSteps
        real(8), intent(out) :: result_sum(:)
        real(8), intent(inout) :: phi(:,:), eps
        integer :: i, j, m, n
        real(8) :: prob_A(nPhase, nSteps), prob_B(nPhase, nSteps), d(3), tmp

        n = size(AdjA,1)

        allocate( state_Aorig(n), state_Borig(n) )
        call init_state( AdjA, state_Aorig)
        call init_state( AdjB, state_Borig)
        result_sum = 0      
        
        ! -- Shared in loop --
        ! state_Xorig   - original state of graph X (A,B)
        ! -- Private to each thread -- 
        ! state_X       - current state of graph X (A,B)
        ! state_Xtmp    - temporary flip state of current thread iteration for graph X (A,B)        

        !$omp parallel private(state_A, state_B, state_Atmp, state_Btmp, prob_A, prob_B, d, tmp) &
        !$omp& shared(AdjA, AdjB, state_Aorig, state_Borig, nPhase, nSteps, phi, n, eps)
            allocate( state_A(n), state_Atmp(n), state_B(n), state_Btmp(n)) 
            ! Use copy state to allocate temporary state variables. Only need to allocate space for flip operation
            call copy_state( state_Aorig, state_Atmp, 1)    
            call copy_state( state_Borig, state_Btmp, 1)
            !$omp do reduction(+: result_sum)
            do i = 1, n
                call copy_state( state_Aorig, state_A)
                call Apply_Walk( state_A, state_Atmp, AdjA, prob_A, (/ i /), nPhase, nSteps, phi)
                do j = 1, n
                    call copy_state( state_Borig, state_B)
                    call Apply_Walk( state_B, state_Btmp, AdjB, prob_B, (/ j /), nPhase, nSteps, phi)
                    d = 0
                    do m = 1, nPhase
                        ! Find the euclidean distance between the probability vectors
                        tmp = Norm2(prob_A(m,:) - prob_B(m,:))/nSteps
                        d(1) = d(1) + tmp
                        d(2) = d(2) + tmp**2
                        d(3) = Matusita_dist(prob_A(m,:), prob_B(m,:)) / nSteps
                    end do
                    ! Add sum by 1 if probability matches, based on threshold
                    if( d(1) > eps) then           
                        result_sum(1) = result_sum(1) + 1
                    end if
                    do m = 1, 3
                        result_sum(m+1) = result_sum(m+1) + d(m)
                    end do
                end do    
            end do
            !$omp end do
            call clear_state(state_A)
            call clear_state(state_B)
            call clear_state(state_Atmp)
            call clear_state(state_Btmp)
            deallocate( state_A, state_Atmp, state_B, state_Btmp)
        !$omp end parallel

        ! Clear state variables
        call clear_state(state_Aorig)
        call clear_state(state_Borig)
        deallocate(state_Aorig, state_Borig)

    end

    ! Applies the coined quantum walk to the initial state
    ! state - initial input state
    ! state_tmp - temporary state vector, initialised with correct memory allocation for flip operation
    ! Adj - Graph adjacency matrix
    ! nodes - vector containing reference nodes
    ! nSteps - number of steps the walker must take
    ! prob - output probability vector
    subroutine Apply_Walk( state, state_tmp, Adj, prob, nodes, nPhase, nSteps, phi)
        type(node), dimension(:) :: state, state_tmp
        integer, intent(in) :: nodes(:), nPhase, nSteps
        integer, intent(in), dimension(:,:) :: Adj
        real(8), intent(inout), dimension(:,:) :: prob, phi
        real(8) :: probTmp(size(state))
        integer :: m, n

        ! n = size(prob,1)
        do m = 1, nSteps
            call Walk( state, state_tmp, 1, Adj, nodes, phi(m,:))
            call ProbDist(state,probTmp)
            do n = 1, nPhase
                prob(n,m) = probTmp(nodes(n))
            end do
        end do
    end

    ! Determines a pemutation matrix copy of the input matrix, using a permutation vector 
    ! A - input adjacency matrix 
    ! A2 - output permutated copy of A 
    ! perm - input permutation vector
    subroutine MatrixPermutation(A, A2, perm)
        integer, dimension(:,:), intent(inout) :: A, A2
        integer, dimension(:), intent(in) :: perm
        real(8), allocatable, dimension(:,:) :: Pinv
        integer, allocatable :: P(:,:)
        integer :: n, i

        n = size(A,1)

        allocate( P(n,n), Pinv(n,n))

        P = 0
        do i = 1, n 
            P(i, perm(i)) = 1
        end do

        call mat_inv(real(P,8), Pinv) 
      
        ! P.A.Inv(P)
        A2 = matmul( matmul(P,A), int(Pinv))

        deallocate(P, Pinv)
    end

end module