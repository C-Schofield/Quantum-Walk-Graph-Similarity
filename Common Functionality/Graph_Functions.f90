module Graph_Functions
    ! Author: Callum Schofield - 2017
    ! Collection of functions which are used to gain information or modify graphs
    implicit none

    public Remove_Rand_Edge, Count_Edges, Add_Rand_Edge, Change_Graph, Modify_Weight
    private D_Remove_Rand_Edge, I_Remove_Rand_Edge
    private D_Count_Edges, I_Count_Edges
    private D_Add_Rand_Edge, I_Add_Rand_Edge
    private D_Change_Graph, I_Change_Graph

    ! Interfaces to allow double or integer input graphs
    interface Add_Rand_Edge 
        module procedure D_Add_Rand_Edge, I_Add_Rand_Edge
    end interface

    interface Remove_Rand_Edge
        module procedure D_Remove_Rand_Edge, I_Remove_Rand_Edge
    end interface Remove_Rand_Edge

    interface Count_Edges
        module procedure D_Count_Edges, I_Count_Edges
    end interface Count_Edges

    interface Change_Graph 
        module procedure D_Change_Graph, I_Change_Graph
    end interface

    contains

    ! -- Change_Graph -- !
    ! Adj - input graph 
    ! nChanges - number of modifications to make 
    ! directed - boolean to declare if the input graph is to be treated as directed or not 
    ! mode - type of random modifications to be made 
    !   1 - Add edges 
    !   2 - Remove edges 
    !   3 - Add and/or remove edges 
    !   4 - Change edge weights (without adding or removing edges) 

    ! Make modifications to real input graph
    subroutine D_Change_Graph( Adj, nChanges, directed, mode)
        implicit none 
        real(8), intent(inout), dimension(:,:) :: Adj
        integer, intent(inout) :: nChanges
        integer, intent(in) :: mode
        integer :: i, ierr
        real(8) :: tmp
        logical, intent(in) :: directed
        
        select case(mode)
            case(1)
                ! Add edges
                do i = 1, nChanges
                    call Add_Rand_Edge(Adj, directed, ierr)
                    if(ierr == 1) then 
                        ! Cannot add more edges, stop and return actual number of changes made
                        nChanges = i - 1
                        return
                    end if
                end do  
            case(2)
                ! Remove edges
                do i = 1, nChanges
                    call Remove_Rand_Edge(Adj, directed, ierr)
                    if(ierr == 1) then 
                        ! Cannot remove more edges, stop and return actual number of changes made
                        nChanges = i - 1
                        return
                    end if
                end do 
            case(3)
                ! Add and remove edges (random 50/50)
                do i = 1, nChanges
                    ierr = 1
                    ! If edge cannot be added/removed, try again
                    do while( ierr == 1)
                        call random_number(tmp)
                        ! Round tmp to 0 or 1
                        tmp = floor( tmp + 0.5d0)
                        if( tmp == 0 )then 
                            call Add_Rand_Edge(Adj, directed, ierr)
                        else
                            call Remove_Rand_Edge(Adj, directed, ierr)
                        end if
                    end do
                end do
            case(4)
                ! Change weight of a selected edge 
                do i = 1, nChanges 
                    ierr = 0
                    call Modify_Weight( Adj, directed, ierr)
                    if( ierr == 1)then 
                        write(*,*)'Cannot modify edge weight if graph has no edges.'
                        exit
                    end if
                end do
            case default 
                write(*,*)'Invalid mode for changing graph, valid modes are (1,2,3,4). &
                            & Using mode 2 (removing edges only).'
                ! Remove edges
                do i = 1, nChanges
                    call Remove_Rand_Edge(Adj, directed, ierr)
                    if(ierr == 1) then 
                        ! Cannot remove more edges, stop and return actual number of changes made
                        nChanges = i - 1
                        return
                    end if
                end do
        end select      
    end

    ! Make modifications to integer input graph
    subroutine I_Change_Graph( Adj, nChanges, directed, mode)
        implicit none 
        integer, intent(inout), dimension(:,:) :: Adj
        integer, intent(inout) :: nChanges
        integer, intent(in) :: mode
        integer :: i, ierr
        real(8) :: tmp
        logical :: directed
        
        select case(mode)
            case(1)
                ! Add edges
                do i = 1, nChanges
                    call Add_Rand_Edge(Adj, directed, ierr)
                    if(ierr == 1) then 
                        ! Cannot add more edges, stop and return actual number of changes made
                        nChanges = i - 1
                        return
                    end if
                end do  
            case(2)
                ! Remove edges
                do i = 1, nChanges
                    call Remove_Rand_Edge(Adj, directed, ierr)
                    if(ierr == 1) then 
                        ! Cannot remove more edges, stop and return actual number of changes made
                        nChanges = i - 1
                        return
                    end if
                end do 
            case(3)
                ! Add and remove edges (random 50/50)
                do i = 1, nChanges
                    ierr = 1
                    ! If edge cannot be added/removed, try again
                    do while( ierr == 1)
                        call random_number(tmp)
                        ! Round tmp to 0 or 1
                        tmp = floor( tmp + 0.5d0)
                        if( tmp == 0 )then 
                            call Add_Rand_Edge(Adj, directed, ierr)
                        else
                            call Remove_Rand_Edge(Adj, directed, ierr)
                        end if
                    end do
                end do
            case default
                write(*,*)'Invalid mode for changing graph, valid modes are (1,2,3). &
                            & Using mode 2 (removing edges only).'
                ! Remove edges
                do i = 1, nChanges
                    call Remove_Rand_Edge(Adj, directed, ierr)
                    if(ierr == 1) then 
                        ! Cannot remove more edges, stop and return actual number of changes made
                        nChanges = i - 1
                        return
                    end if
                end do 
        end select        
    end

    ! -- Remove_Rand_Edge -- !
    ! Remove an edge from a random location in the graph A 
    ! A - input adjacency matrix 
    ! directed - boolean to declare if the input graph is to be treated as directed or not 
    ! ierr - error tag

    ! Remove random edge in input double precision matrix
    subroutine D_Remove_Rand_Edge( A, directed, ierr)
        implicit none
        real(8), dimension(:,:), intent(inout) :: A
        logical, intent(in) :: directed
        integer, intent(out) :: ierr
        real(8) :: rand
        integer :: nEdges, edgeToDel, edgeCount, i, j, n
        logical :: isDeleted

        ierr = 0
        n = size(A, 1)
        ! count the number of edges
        nEdges = D_Count_Edges(A, directed)

        ! Return/don't remove an edge if graph has no edges
        if(nEdges == 0) then 
            ierr = 1
            return
        end if 

        edgeCount = 0
        ! Randomly determine the edge to be deleted, between 1 and current number of edges
        call random_number(rand)
        edgeToDel = floor(nEdges*rand) + 1
        isDeleted = .false.
        if( directed)then
            do i = 1, n
                do j = 1, n
                    if( A(i,j) > 0 ) edgeCount = edgeCount + 1
                    if( edgeCount == edgeToDel) then
                        A(i,j) = 0
                        isDeleted = .true.
                    end if
                    if( isDeleted) exit
                end do
                if(  isDeleted) exit
            end do
        else 
            ! Only count upper diagonal, then set lower diagonal accordingly
            do i = 1, n
                do j = i, n
                    if( A(i,j) > 0 ) edgeCount = edgeCount + 1
                    if( edgeCount == edgeToDel) then
                        A(i,j) = 0
                        A(j,i) = 0
                        isDeleted = .true.
                    end if
                    if( isDeleted) exit
                end do
                if(  isDeleted) exit
            end do
        end if
    end subroutine

    ! Remove random edge in input integer precision matrix
    subroutine I_Remove_Rand_Edge( A, directed, ierr)
        implicit none
        integer, dimension(:,:), intent(inout) :: A
        logical, intent(in) :: directed
        integer, intent(out) :: ierr
        real(8) :: rand
        integer :: nEdges, edgeToDel, edgeCount, i, j, n
        logical :: isDeleted

        ierr = 0
        n = size(A, 1)
        ! count the number of edges
        nEdges = I_Count_Edges(A, directed)

        ! Return/don't remove an edge if graph has no edges
        if(nEdges == 0) then 
            ierr = 1
            return
        end if 

        edgeCount = 0
        ! Randomly determine the edge to be deleted, between 1 and current number of edges
        call random_number(rand)
        edgeToDel = floor(nEdges*rand) + 1
        isDeleted = .false.
        if( directed)then
            do i = 1, n
                do j = 1, n
                    if( A(i,j) == 1 ) edgeCount = edgeCount + 1
                    if( edgeCount == edgeToDel) then
                        A(i,j) = 0
                        isDeleted = .true.
                    end if
                    if( isDeleted) exit
                end do
                if(  isDeleted) exit
            end do
        else 
            ! Only count upper diagonal, then set lower diagonal accordingly (include diagonal)
            do i = 1, n
                do j = i, n
                    if( A(i,j) == 1 ) edgeCount = edgeCount + 1
                    if( edgeCount == edgeToDel) then
                        A(i,j) = 0
                        A(j,i) = 0
                        isDeleted = .true.
                    end if
                    if( isDeleted) exit
                end do
                if(  isDeleted) exit
            end do
        end if
    end subroutine

    ! -- Add_Rand_Edge -- !
    ! Add an edge to a random location in the graph A (where there is no edge)
    ! A - input adjacency matrix 
    ! directed - boolean to declare if the input graph is to be treated as directed or not 
    ! ierr - error tag

    ! Add random edge to double precision matrix
    ! weight_tmp - weight to set new edge to (optional)
    subroutine D_Add_Rand_Edge(A, directed, ierr, weight_tmp)
        implicit none   
        real(8), intent(inout) :: A(:,:)
        real(8), intent(in), optional :: weight_tmp
        logical, intent(in) :: directed
        integer, intent(out) :: ierr
        real(8) :: rand, weight
        integer :: maxEdges, nEdges, edgeToAdd, edgeCount, i, j, n
        logical :: isAdded

        ierr = 0
        weight = 1.d0

        ! Update weight if it is included
        if(present(weight_tmp))weight = weight_tmp

        n = size(A, 1)
        ! count the number of edges
        nEdges = D_Count_Edges(A, directed)

        ! Return/don't add an edge if graph is complete
        if( directed)then
            maxEdges = n**2
        else 
            maxEdges = (n**2-n)/2 + n
        end if
        if(nEdges == maxEdges) then 
            ierr = 1
            return
        end if 
        edgeCount = 0
        ! Randomly determine the edge to be added, between 1 and current number of "empty" edges
        call random_number(rand)
        edgeToAdd = floor((maxEdges - nEdges)*rand) + 1
        isAdded = .false.
        if( directed)then
            do i = 1, n
                do j = 1, n
                    ! Cycle through until enough valid "non-egdes" have been found
                    if( A(i,j) < 0.0000001d0 ) edgeCount = edgeCount + 1
                    if( edgeCount == edgeToAdd) then
                        A(i,j) = weight
                        isAdded = .true.
                    end if
                    if( isAdded) exit
                end do
                if(  isAdded) exit
            end do
        else 
            ! Only count upper diagonal, then set lower diagonal accordingly
            do i = 1, n
                do j = i, n
                    if( A(i,j)  < 0.0000001d0 ) edgeCount = edgeCount + 1
                    if( edgeCount == edgeToAdd) then
                        A(i,j) = weight
                        A(j,i) = weight
                        isAdded = .true.
                    end if
                    if( isAdded) exit
                end do
                if(  isAdded) exit
            end do
        end if
    end 

    ! Add random edge to integer precision matrix
    subroutine I_Add_Rand_Edge(A, directed, ierr)
        implicit none   
        Integer, intent(inout) :: A(:,:)
        logical, intent(in) :: directed
        integer, intent(out) :: ierr
        real(8) :: rand
        integer :: maxEdges, nEdges, edgeToAdd, edgeCount, i, j, n
        logical :: isAdded

        ierr = 0
        n = size(A, 1)
        ! count the number of edges
        nEdges = I_Count_Edges(A, directed)

        ! Return/don't add an edge if graph is complete
        if( directed)then
            maxEdges = n**2
        else 
            maxEdges = (n**2-n)/2 + n
        end if
        if(nEdges == maxEdges) then 
            ierr = 1
            return
        end if 
        edgeCount = 0
        ! Randomly determine the edge to be added, between 1 and current number of "empty" edges
        call random_number(rand)
        edgeToAdd = floor((maxEdges - nEdges)*rand) + 1
        isAdded = .false.
        if( directed)then
            do i = 1, n
                do j = 1, n
                    ! Cycle through until enough valid "non-egdes" have been found
                    if( A(i,j) == 0 ) edgeCount = edgeCount + 1
                    if( edgeCount == edgeToAdd) then
                        A(i,j) = 1
                        isAdded = .true.
                    end if
                    if( isAdded) exit
                end do
                if(  isAdded) exit
            end do
        else 
            ! Only count upper diagonal, then set lower diagonal accordingly
            do i = 1, n
                do j = i, n
                    if( A(i,j)  == 0 ) edgeCount = edgeCount + 1
                    if( edgeCount == edgeToAdd) then
                        A(i,j) = 1
                        A(j,i) = 1
                        isAdded = .true.
                    end if
                    if( isAdded) exit
                end do
                if(  isAdded) exit
            end do
        end if

    end 
    
    ! -- Modify_Weight -- !
    ! Changes weight of edges in input double precision graph, note this does not add/remove (eps <= edge weight <= 1)
    ! A - input graph adjacency matrix 
    ! directed - boolean to declare if the input graph is to be treated as directed or not 
    ! ierr - error tag
    subroutine Modify_Weight( A, directed, ierr)
        implicit none
        real(8), dimension(:,:), intent(inout) :: A
        logical, intent(in) :: directed
        integer, intent(out) :: ierr
        real(8) :: rand, eps
        integer :: nEdges, edgeToMod, edgeCount, i, j, n
        logical :: isChanged

        eps = 0.001d0
        ierr = 0
        n = size(A, 1)
        ! count the number of edges
        nEdges = D_Count_Edges(A, directed)

        ! Can't change edge weight if graph has no edges
        if(nEdges == 0) then 
            ierr = 1
            return
        end if 

        edgeCount = 0
        ! Randomly determine the edge to be deleted, between 1 and current number of edges
        call random_number(rand)
        edgeToMod = floor(nEdges*rand) + 1
        isChanged = .false.
        if( directed)then
            do i = 1, n
                do j = 1, n
                    if( A(i,j) > 0 ) edgeCount = edgeCount + 1
                    if( edgeCount == edgeToMod) then
                        ! Ensure eps < rand
                        do while( rand < eps)
                            call random_number(rand)
                        end do 
                        A(i,j) = rand
                        isChanged = .true.
                    end if
                    if( isChanged) exit
                end do
                if(  isChanged) exit
            end do
        else 
            ! Only count upper diagonal, then set lower diagonal accordingly
            do i = 1, n
                do j = i, n
                    if( A(i,j) > 0 ) edgeCount = edgeCount + 1
                    if( edgeCount == edgeToMod) then
                        ! Ensure eps > rand
                        do while( rand < eps)
                            call random_number(rand)
                        end do 
                        A(i,j) = rand
                        A(j,i) = rand
                        isChanged = .true.
                    end if
                    if( isChanged) exit
                end do
                if(  isChanged) exit
            end do
        end if
    end subroutine

    ! -- Count_Edges -- !
    ! Counts the number of edges in the input graph
    ! A - input adjacency matrix
    ! directed - boolean to declare if the input graph is to be treated as directed or not 
    ! nEdges - resulting number of edges

    ! Double matrix
    function D_Count_Edges( A, directed)  result(nEdges)
        implicit none 
        real(8), dimension(:,:), intent(in) :: A
        integer :: nEdges
        logical, intent(in) :: directed
        integer :: i, n, nDiag
        logical, allocatable :: mask(:,:), diag(:)

        n = size(A,1)
        allocate( mask(n,n), diag(n))

        ! Count the number of elements that are greater than 0
        mask = A > 0
        nEdges = count(mask)

        if(.not. directed) then 
            ! If the graph is undirected, then total number of edges is sum of off diagonal + diagonal
            ! i.e. (total sum - diagonal / 2) + diagonal
            do i = 1, n 
               diag(i) = A(i,i) > 0     ! Sum of the diagonal egdes
            end do
            nDiag = count(diag)
            nEdges = nDiag + (nEdges - nDiag)/2 
        end if

        deallocate( mask)
    end 

    ! Integer matrix
    function I_Count_Edges( A, directed)  result(nEdges)
        implicit none 
        integer, dimension(:,:), intent(in) :: A
        integer :: nEdges
        logical, intent(in) :: directed
        integer :: i, n 

        n = size(A,1)

        if(directed) then
            ! Number of edges is given by sum of adjacency matrix  
            nEdges = sum(A) 
        else
            ! If the graph is undirected, then total number of edges is sum of off diagonal + diagonal
            ! i.e. (total sum - diagonal / 2) + diagonal
            nEdges = 0
            ! Sum up the diagonal egdes
            do i = 1, n 
                nEdges = nEdges + A(i,i)
            end do
            nEdges = nEdges + (sum(A) - nEdges)/2 
        end if
    end 

end module