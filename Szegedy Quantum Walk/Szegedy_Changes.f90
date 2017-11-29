program Szegedy_Changes
    ! Author: Callum Schofield - 2017
    ! This code makes a number of changes to graphs to make comparisons
    !   Input 1 file with the filename/path to each graph to be tested 
    !   Each graph will then be tested against itself with a number changes made 
    !   Second input file can be used to specify parameters, else default values file used
    !   Mode specifies how to change the input graph:
    !       1 - Add edges to the graph
    !       2 - Remove edges from the graph 
    !       3 - Add and remove edges to/from the graph, each iteration has a 50% chance of either
    !
    ! Inputs  - nThreads:   number of threads to be used (must be specified in slurm file for Pawsey, not here)          
    !           graphFile:    list containing adjacency matrices to be tested (all matrices must be in the same folder as program)
    !           nSteps:     number of steps to be done
    !           errMult:    error multiplier to be used in Szegedy similarity index classification.
    !           begin:      Number of changes to begin with 
    !           end:        Number of changes to end with. 
    !                        If negative, assume to change as many edges as in the graph and step via percentage given in step
    !           step:       Number of changes to skip by each iteration. Denotes percentage to skip if end < 0
    !                        If negative, set step size to be 1, regardless of end
    !           mode:       Type of changes to be made (1 - Add, 2 - Remove, 3 - Add and Remove)
    !           
    ! Outputs - Output.dat: data file containing similarity index (both classical and quantum) for each graph compared 
    !
    ! Optional command line arguments (note not all arguments contained in parameters file can be entered on command line 
    !   -p, -param:     Parameter file location. NOTE use this first if entering other parameters, as it will overwrite any already loaded parameters
    !   -g, -graphs:    Graph file location
    !   -o, -out:       Output file location 
    !   -t, -threads:   Number of threads to use
    !   -r, -errmult:   Error multiplier to be used in Szegedy similarity index classification. Input a value <= 0 if you want the algorithm to return
    !                    only d, that is 1 / (1+ w*d) can be determined post completion for varying w.
    !   -b, -begin:     Number of changes to start with
    !   -e, -end:       Number of changes to end with
    !   -s, -step:      Number of changes to step with between iteration 
    !   -m, -mode:      Type of changes to be made (1 - Add, 2 - Remove, 3 - Add and Remove)
    !	-n, -trials:    Number of trials to be conducted at each number of changes
    !   -d, -directed:  Determines whether to treat edges as directed or undirected when making changes

    use omp_lib
    use mpi 
    use Szegedy_Similarity
    use Common_Functions
    use Graph_Functions
    use Szegedy_Functions
    use NodeAffinity_Functions

    implicit none
    integer, parameter :: master = 0        ! Master node ID
    integer, parameter :: from_master = 1   ! MPI tag sending to workers
    integer, parameter :: to_master = 2     ! MPI tag receiving from workers

    integer :: nSteps, nThreads, cr, cm, numGraphs, ierr, task_id, num_tasks
    integer :: startChanges, endChanges, stepchanges, nTrials, mode
    real(8) :: errMult, start, finish
    character(len = 100) :: arg, graphFile, OutFile
    logical :: directed

    ! Use MPI to compare multiple graphs across multiple nodes
    call MPI_INIT( ierr)

    call MPI_COMM_RANK (MPI_COMM_WORLD, task_id, ierr)
    call MPI_COMM_SIZE (MPI_COMM_WORLD, num_tasks, ierr) 
    ! Set timer function
    ! call system_clock(start)
    ! call system_clock(count_rate=cr)
    ! call system_clock(count_max=cm)
    ! rate = real(cr)
    start = MPI_Wtime()

    ! Initialise parameters to default
    nSteps = 25                             ! Number of steps of Szegedy to do per iteration
    graphFile = 'Parameters/GraphList.dat'  ! Adjacency matrix input file
    OutFile = 'Output/output.dat'           ! Output file
    nThreads = 12                           ! Number of threads to use
    errMult = 1.d0                          ! Error multiplier 
    startChanges = 0                        ! Number of changes to start with
    endChanges = 200                        ! Number of changes to end with
    stepChanges = 50                        ! Number of changes to step with
    nTrials = 5                             ! Number of trials to do at each change number
    mode = 2                                ! Graph change mode to use 
    directed = .true.                       ! Determines whether to treat edges as directed or undirected when making changes

    call init()

    if( task_id == master )then
        call file_dimensions( graphFile, numGraphs) 
    end if

    call MPI_BCAST( numGraphs, 1, MPI_INT, master, MPI_COMM_WORLD, ierr)

    call run_prog( )

    if( task_id == master )then
        ! call system_clock(finish)
        finish = MPI_Wtime()
    end if
    
    call mpi_finalize( ierr)

contains

    ! Subroutine to initialise program 
    subroutine init() 
        implicit none 
        integer :: i, tmp
        character(len=100) :: arg, file_parameters

       ! Check command line for inputs or parameter file
        do i = 1, iargc(), 2
            call getarg(i,arg)
            select case(arg)
                case ("-p","-param")
                    call getarg(i+1,arg)
                    read(arg,'(A)')file_parameters
                    ! Initialise parameters using input file
                    call parse_input(file_parameters)
                case ("-g", "-graphs")
                    call getarg(i+1,arg)
                    read(arg,'(A)')graphFile
                case ("-o", "-out")
                    call getarg(i+1,arg)
                    read(arg,'(A)')outFile
                case ("-t", "-threads")
                    call getarg(i+1,arg)
                    read(arg,*)nThreads  
                case ("-w", "-weight")
                    call getarg(i+1,arg)
                    read(arg,*)errMult 
                case ("-b", "-begin")
                    call getarg(i+1,arg)
                    read(arg,*)startChanges 
                case ("-e", "-end")
                    call getarg(i+1,arg)
                    read(arg,*)endChanges 
                case ("-s", "-step")
                    call getarg(i+1,arg)
                    read(arg,*)stepChanges 
                case ("-m", "-mode")
                    call getarg(i+1,arg)
                    read(arg,*)mode   
                case ("-n", "-trials")
                    call getarg(i+1,arg)
                    read(arg,*)nTrials 
                case ("-d", "-directed")
                    call getarg(i+1,arg)
                    read(arg,*)tmp
                    if(tmp == 1) then 
                        directed = .true.
                    else if(tmp == 0) then
                        directed = .false.
                    else 
                        directed = .true.
                        write(*,*)'Incorrect input for directed. &
                        & Enter 0 for undirected and 1 for directed. Program will continue using directed.' 
                    end if 
            end select
        end do

        call init_random_seed()
        call omp_set_num_threads(nThreads)

    end 

    subroutine run_prog()
        implicit none
        real(8), allocatable, dimension(:,:) :: A1, A2, probOrig, SOrig, S2
        real(8), allocatable, dimension(:)   :: simIndex, simIndex2, simIndex3, simBelief, simIndex_all, &
                                              & simIndex_all2, simIndex_all3, simBelief_all, tmpGraph
        integer, allocatable, dimension(:)   :: nChanges, displ, graph_split_all, edges_all, edges
        real(8) :: totalTime, time1, time2, simTmp(3)
        integer :: n, i, j, k, nTest, change_err, status(MPI_STATUS_SIZE), nEdges
        integer :: graph_extra, graph_splitTmp, graph_split, tmp, stepChanges_tmp
        logical :: outFile_exists, changePerc
        character(len = 50) :: strTmp
        character(len=100)  :: graph
        character(len=200)  :: line
        
        if( task_id == master )then
            changePerc = .false.
            if(endChanges < 0)changePerc = .true.   ! Make changes based on percentage of edges
            ! Open list of adjacency matrices
            open(unit=10, file=graphFile, action='read')
            ! Open output file
            inquire(file = outFile, exist = outFile_exists)
            if( .not. outFile_exists) then 
                open(20, file=outFile, action='write')
                write(20,*)'The similarity values are output as such:'
                write(20,*)'  Changes made in second graph'
                write(20,*)'  Percentage change in edges between the two graphs'
                write(20,*)'  Similarity index d value using Szegedy: abs(prob1 - prob2), need to  1 / (1 + d*ErrMult)'
                write(20,*)'  Similarity index d value using Szegedy: 2Norm(prob1, prob2), need to  1 / (1 + d*ErrMult)'
                write(20,*)'  Similarity index d value using Szegedy: Matus_dist(prob1, prob2), need to  1 / (1 + d*ErrMult)'
                write(20,*)'  Similarity index using Node Affinity'
            else 
                open(20, file=outFile, position='append', action='write')
            end if 
            
            allocate( graph_split_all(num_tasks), displ(num_tasks))
        end if
        do i = 1, numGraphs
            ! call system_clock(time1)
            if( task_id == master )then
                time1 = MPI_Wtime()
                graph_split_all = 0
                read(10,'(A)')line
                read(line,'(A)')graph          ! input for graph A location 
                call file_dimensions( graph, n)
                print*, 'Beginning graph ', graph
            end if

            call MPI_BCAST( n, 1, MPI_INT, master, MPI_COMM_WORLD, ierr)
            allocate( A1(n,n), A2(n,n), tmpGraph(n), probOrig(nSteps,n), SOrig(n,n), S2(n,n))

            if( task_id == master )then
                open(unit=40, file=graph, action='read')
                read(40, *)A1
                close(40)

                if( size(A1,1) /= size(A1,2) ) then
                    write(*,*)'Adjacency matrix ',graph,' is not n x n'
                    deallocate( A1 , probOrig, SOrig, S2)
                    cycle
                end if

                ! Fortran writes column major - so transpose matrix
                A1 = transpose(A1)
                nEdges = Count_Edges(A1, directed)
                stepChanges_tmp = stepChanges
                ! Change based on percentage of edges
                if(changePerc)then 
                    endChanges = nEdges
                    ! stepChanges is given as a %
                    stepChanges_tmp = endChanges*stepChanges / 100.d0
                    if(stepChanges_tmp == 0)stepChanges_tmp = 1     ! Set step to 1 if too small
                end if
                
                if( stepChanges < 0)stepChanges_tmp = 1 !stepChanges is set to 1
                        
                nTest = (endChanges - startChanges) / stepChanges_tmp + 1
                allocate(nChanges(nTest*nTrials), edges_all(nTest*nTrials), simIndex_all(nTest*nTrials), &
                       & simIndex_all2(nTest*nTrials), simIndex_all3(nTest*nTrials), simBelief_all(nTest*nTrials) )
                ! Create an array containing all the tests to be made (i.e. number of changes) 
                do j = 1, nTest
                    do k = 1, nTrials
                        nChanges(k + (j - 1)*nTrials) = stepChanges_tmp * (j - 1) + startChanges 
                    end do 
                end do
                ! Determine how many tests must be completed by each node, and send this number to them
                graph_splitTmp = nTest*nTrials / num_tasks
                graph_extra = mod(nTest*nTrials, num_tasks)
                tmp = graph_splitTmp + 1
                do j = 1, num_tasks - 1
                    if( j <= graph_extra) then
                        graph_split = graph_splitTmp + 1
                    else
                        graph_split = graph_splitTmp
                    end if
                    graph_split_all(j+1) = graph_split
                    call MPI_Send( graph_split, 1, MPI_INT, j, from_master, MPI_COMM_WORLD, ierr)
                    ! Send the list of changes to be made and compared to each node as well
                    call MPI_Send( nChanges(tmp:tmp+graph_split-1), graph_split, MPI_INT, j, from_master, MPI_COMM_WORLD, ierr)
                    tmp = tmp + graph_split
                end do
                graph_split = graph_splitTmp
                graph_split_all(1) = graph_split
                ! Displ gives the displacement for the gathering. 
                displ(1) = 0! nTest - graph_split_all(1) 
                if( num_tasks > 1)then
                    do j = 2, num_tasks 
                        displ(j) = displ(j-1) + graph_split_all(j-1)        
                    end do
                end if
                
                ! Generate the probability matrix for the original matrix so it is not determined each iteration
                ! Use A2 as a tmp variable, in case algorithm modifies input matrix
                A2 = A1
                call genProbAll( A2, probOrig, nSteps)      
                call belief( A2, SOrig) 
                do j = 1, num_tasks - 1
                    do k = 1, n
                        tmpGraph = A1(:,k)
                        call MPI_Send( tmpGraph, n, MPI_DOUBLE, j, from_master, MPI_COMM_WORLD, ierr)
                    end do
                end do
            end if
            
            if( task_id /= master )then
                call MPI_Recv( graph_split, 1, MPI_INT, master, from_master, MPI_COMM_WORLD, status, ierr)
                allocate( nChanges(graph_split))
                call MPI_Recv( nChanges, graph_split, MPI_INT, master, from_master, MPI_COMM_WORLD, status, ierr)
                do j = 1, n
                    call MPI_Recv( tmpGraph, n, MPI_DOUBLE, master, from_master, MPI_COMM_WORLD, status, ierr)
                    A1(:,j) = tmpGraph
                end do
            end if
            deallocate(tmpGraph)
            allocate( simIndex(graph_split), simIndex2(graph_split), simIndex3(graph_split), &
                    & simBelief(graph_split), edges(graph_split) )
            if( num_tasks > 1)then
                ! Broadcast the original probabilities and classical results for comparison
                call MPI_BCAST( probOrig, nSteps*n, MPI_DOUBLE, master, MPI_COMM_WORLD, ierr)
                call MPI_BCAST( SOrig, n**2, MPI_DOUBLE, master, MPI_COMM_WORLD, ierr)
            end if

            !$omp parallel private( A2, S2, simTmp) &
            !$omp&         shared( simIndex, simIndex2, simIndex3, simBelief, SOrig, probOrig, nChanges, mode, edges, nSteps)
                !$omp do schedule(static)
                    do j = 1, graph_split
                        A2 = A1
                        call change_graph( A2, nChanges(j), directed, mode)
                        edges(j) = Count_Edges(A2, directed)
                        simTmp = Szegedy_Compare( A2, probOrig, nSteps, errMult)
                        simIndex(j) = simTmp(1)
                        simIndex2(j) = simTmp(2)
                        simIndex3(j) = simTmp(3)
                        call belief(A2, S2)
                        simBelief(j) = NodeAffinity_Compare( SOrig, S2)
                    end do
                !$omp end do
            !$omp end parallel

            deallocate( A1, A2)

            ! Collect the results
            if( num_tasks > 0)then
                call MPI_GatherV( simIndex, graph_split, MPI_DOUBLE, simIndex_all, graph_split_all, displ, &
                                    & MPI_DOUBLE, master, MPI_COMM_WORLD, ierr)
                call MPI_GatherV( simIndex2, graph_split, MPI_DOUBLE, simIndex_all2, graph_split_all, displ, &
                                    & MPI_DOUBLE, master, MPI_COMM_WORLD, ierr)
                call MPI_GatherV( simIndex3, graph_split, MPI_DOUBLE, simIndex_all3, graph_split_all, displ, &
                                    & MPI_DOUBLE, master, MPI_COMM_WORLD, ierr)
                call MPI_GatherV( simBelief, graph_split, MPI_DOUBLE, simBelief_all, graph_split_all, displ, &
                                    & MPI_DOUBLE, master, MPI_COMM_WORLD, ierr)
                call MPI_GatherV( edges, graph_split, MPI_INT, edges_all, graph_split_all, displ, &
                                    & MPI_INT, master, MPI_COMM_WORLD, ierr)
            else
                ! No collecting to do if only 1 node
                simIndex_all = simIndex
                simIndex_all2 = simIndex2
                simIndex_all3 = simIndex3
                simBelief_all = simBelief
                edges_all = edges
            end if
            if( task_id == master )then
                
                time2 = MPI_Wtime()
                ! totalTime = (time2 - time1)/rate
                totalTime = time2 - time1
                write(strTmp,'(A,F10.4)')'Execution time: ', totalTime
                write(*,'(A,A,F10.4)')trim(adjustl(graph)), ' completed. Time taken: ', totalTime
                write(20,*)'-----------------------------------'
                write(20,'(A)')trim(adjustl(graph))
                write(20,'(A, I3)')'Number of vertices: ', n
                write(20,'(A, I5)')'Number of edges: ', nEdges
                write(20,'(A, I3)')'Number of steps: ', nSteps
                write(20,'(A, F10.4)')'Error multiplier: ', errMult
                write(20,'(A, I3)')'Threads: ', nThreads
                write(20,'(A, I3)')'Nodes: ', num_tasks
                write(20,'(A)')trim(adjustl(strTmp))
                write(20,*)nChanges
                write(20,*)100.d0*(edges_all - nEdges)/nEdges
                write(20,*)simIndex_all 
                write(20,*)simIndex_all2
                write(20,*)simIndex_all3
                write(20,*)simBelief_all
                flush(20)
                deallocate( simIndex_all, simIndex_all2, simIndex_all3, simBelief_all, edges_all )
            end if
            deallocate(simIndex, simIndex2, simIndex3, simBelief, probOrig, SOrig, edges, nChanges, S2 )
        end do
        if( task_id == master )then
            close(10)
            close(20)
            deallocate(graph_split_all, displ)
        end if
    end subroutine

    subroutine parse_input( inFile)
        implicit none 
        character(len=*), intent(in) :: inFile 
        character(len=100) :: line
        integer :: i, readState, tmp

        open(50, file=inFile)
        ! Index current input
        i = 1

        ! Loop until end of file
        do
            read(50, '(A)', iostat=readState)line
            ! If iostat is <0 then end of file reached
            if( readState < 0) exit
            ! If line starts with comment character, cycle to next line
            if ( index(line, "#") == 1 )  cycle
            select case (i)
                case (1)
                    read(line,'(A)')graphFile 
                case (2)
                    read(line,'(A)')outFile
                case (3)
                    read(line,*)nSteps
                case (4)
                    read(line,*)nThreads
                case (5)
                    read(line,*)startChanges
                case(6)
                    read(line,*)endChanges
                case(7)
                    read(line,*)stepChanges
                case(8)
                    read(line,*)nTrials
                case(9)
                    read(line,*)mode
                case(10)
                    read(line,*)tmp
                    if(tmp == 1) then 
                        directed = .true.
                    else if(tmp == 0) then
                        directed = .false.
                    else 
                        directed = .true.
                        write(*,*)'Incorrect input for directed. &
                        & Enter 0 for undirected and 1 for directed. Program will continue using directed.' 
                    end if
                case(11)
                    read(line,*)errMult
            end select
            i = i + 1
        end do
        close(50)
    end 

end program Szegedy_Changes