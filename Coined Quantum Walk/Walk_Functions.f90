module Walk_Functions
    ! Author: Callum Schofield - 2017
    ! As the space state and coin operator dimensions depend on each nodes degree we 
    ! cannot declare the state space / coin operators using a single matrix of regular size.
    ! I have defined a data type which contains the space and coin operator for each
    ! individual node, which then allows the state space information to be stored in an array 
    ! of (type) node of length n. 
    use data_types

    implicit none

    public Walk, ProbDist, init_state, clear_state, copy_state
    private FlipCoin, Step, G, KroneckerDelta
    
    interface copy_state 
        module procedure copy_state_all, copy_state_single
    end interface

    contains

    ! --- PUBLIC FUNCTIONS --- !

    ! Performs a walk on the graph. Usage:
    !   Use init_state to initialise state and state_tmp variables
    !   t - number of steps
    !   Adj - Adjacency matrix of graph to be walked on 
    !   nodes- reference odes to apply phase to
    !   phi - phase to apply to reference nodes
    subroutine Walk( state, state_tmp, t, Adj, nodes, phi)
        type(node), dimension(:), intent(inout) :: state, state_tmp
        real(8), intent(in) :: phi(:)
        integer, intent(in) :: t, nodes(:), Adj(:,:)
        integer :: i
        
        do i = 1, t
            call Step(state, state_tmp, Adj, nodes, phi)
        end do

    end

    ! Gives the probability to be in a certain node for a given state space
    subroutine ProbDist( state, prob)
        type(node), dimension(:), intent(in) :: state
        real(8), dimension(:), intent(out) :: prob
        integer :: i, n

        n = size(state)
        do i = 1, n
            prob(i) = sum(abs(state(i)%space)**2)
        end do
    end

    ! Initialises the space state and coin operators. Usage:
    !   A - Adjacency matrix of graph
    !   state - allocates and initialises the state
    !   state_tmp - allocates the temporary state - for use in FlipCoin 
    !   Note: States MUST be allocated to length n prior to calling. This only allocates the variables contained within the state
    subroutine init_state(A, state)
        type(node), dimension(:), intent(inout) :: state
        integer, dimension(:,:), intent(inout) :: A
        integer, allocatable :: d(:)
        integer :: n, i, j, tmp

        n = size(state)
        ! Get the degree of each node (sum rows)
        allocate(d(n))
        d = sum(A,1)
        ! Create self loops for isolated nodes
        do i = 1, n
            if(d(i) == 0)then
                A(i,i) = 1
                d(i) = 1
            end if
        end do
        ! Relabel edges in adjacency matrix, for mapping to space stats 
        do i = 1, n 
            tmp = 0
            do j = 1,n
                if(A(i,j) /= 0)then
                    tmp = tmp + 1
                    A(i,j) = tmp
                end if
            end do
        end do
        ! Initialise the space state to equal superposition 
        do i = 1, n 
            ! allocate space vector i in state space to length degree(i)
            allocate( state(i)%space(d(i)) )   
            state(i)%space = 1 / sqrt(real(n*d(i)))
        end do
        ! Initialise the coin operator for each node
        do i = 1, n 
            allocate( state(i)%coin(d(i), d(i)))
            call G(d(i), state(i)%coin)
        end do 

        deallocate(d)
        
    end

    ! Allocate a new state and copy data into it from another state 
    subroutine copy_state_all( state_old, state_new)
        implicit none 
        type(node), dimension(:), intent(in)    :: state_old 
        type(node), dimension(:), intent(inout) :: state_new
        integer :: n, i, d

        n = size(state_old) 
        
        do i = 1, n 
            d = size( state_old(i)%space )
            if(.not. allocated(state_new(i)%space) ) allocate( state_new(i)%space(d) )
            if(.not. allocated(state_new(i)%coin) ) allocate( state_new(i)%coin(d, d) )
            state_new(i)%space = state_old(i)%space 
            state_new(i)%coin = state_old(i)%coin
        end do 
    end
    
    ! Allocate a new state and copy data into it from another state 
    ! This form allows to specifiy whether to copy space/coin only
    !   mode 1 - copy space
    !        2 - copy coins
    subroutine copy_state_single( state_old, state_new, mode)
        implicit none 
        type(node), dimension(:), intent(in)    :: state_old 
        type(node), dimension(:), intent(inout) :: state_new
        integer, intent(in) :: mode
        integer :: n, i, d

        n = size(state_old) 
        
        do i = 1, n 
            d = size( state_old(i)%space )
            if( mode == 1) then 
                if(.not. allocated(state_new(i)%space) ) allocate( state_new(i)%space(d) )
                state_new(i)%space = state_old(i)%space 
            else 
                if(.not. allocated(state_new(i)%coin) ) allocate( state_new(i)%coin(d, d) )
                state_new(i)%coin = state_old(i)%coin
            end if 
        end do 
    end

    ! Deallocate state memory
    subroutine clear_state( state)
        type(node), dimension(:), intent(inout) :: state
        integer :: n, i

        n = size(state)
        do i = 1, n 
            if(allocated(state(i)%space)) deallocate(state(i)%space)
            if(allocated(state(i)%coin)) deallocate(state(i)%coin)
        end do

    end

    ! --- PRIVATE FUNCTIONS --- !

    ! Applies the coin flip operation on input state and stores it in the temporary state
    ! input - state_x 
    ! output - state_y (temporary state)
    subroutine FlipCoin( state_x, state_y)
        type(node), dimension(:), intent(in)  :: state_x
        type(node), dimension(:), intent(inout) :: state_y

        integer :: i, j, n, m 

        n = size(state_x)
        do i = 1, n
            m = size(state_x(i)%space)
            do j = 1, m
                state_y(i)%space(j) = dot_product(state_x(i)%coin(j,:),state_x(i)%space)
            end do
        end do

    end

    ! Applies a single step of the coined quantum walk
    ! state_x - initial state 
    ! state_y - allocated temporary state 
    ! Adj - Adjacency matrix 
    ! nodes - reference nodes to apply phase to 
    ! phi - phase to apply to reference nodes
    subroutine Step(state_x, state_y, Adj, nodes, phi)
        type(node), dimension(:), intent(inout) :: state_x, state_y
        integer, intent(in) :: nodes(:), Adj(:,:)
        real(8), intent(in) :: phi(:)
        integer :: n, i, j

        ! Initialise temporary state, y
        n = size(state_x)

        ! Apply phase to each reference node
        do i = 1, size(nodes)
            state_x(nodes(i))%space = state_x(nodes(i))%space * exp(ii*phi(i))
        end do

        ! Apply the coin operator 
        call FlipCoin( state_x, state_y)

        ! Apply translation operator
        do i = 1, n
            do j = 1, n
                if( Adj(j, i) /= 0)then 
                    state_x(i)%space(Adj(i,j)) = state_y(j)%space(Adj(j,i))
                end if
            end do 
        end do
    end

    ! Returns the kronecker delta of inputs i,j
    function KroneckerDelta(i, j) result(output) 
        integer :: i, j, output 
        
        output = 0
        if( i == j)output = 1
    end 

    ! Returns the grover coin of dimension s
    subroutine G(s, result)
        real(8), dimension(:,:), intent(out) :: result
        integer, intent(in) :: s
        integer :: i, j
        
        result = 2.d0 / s
        do i = 1, s 
            do j = 1, s 
                result(i,j) = result(i,j) - KroneckerDelta(i,j)
            end do 
        end do
    end

end module