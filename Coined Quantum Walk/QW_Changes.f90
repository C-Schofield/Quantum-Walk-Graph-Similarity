program QW_Changes
    ! Author: Callum Schofield - 2017
    ! This code makes a number of changes to graphs to make comparisons
    !   Input 1 file with the filename/path to each graph to be tested 
    !   Each graph will then be tested against itself with a number changes made 
    !   Second input file can be used to specify parameters, else parameters in default parameter file used
    !   Mode specifies how to change the input graph:
    !       1 - Add edges to the graph
    !       2 - Remove edges from the graph 
    !       3 - Add and remove edges to/from the graph, each iteration has a 50% chance of either
    
    ! Optional command line arguments (note not all arguments contained in parameters file can be entered on command line 
    !   -p, -param:     Parameter file location. NOTE use this first if entering other parameters, as it will overwrite any already loaded parameters
    !   -g, -graphs:    Graph file location
    !   -o, -out:       Output file location 
    !   -t, -threads:   Number of threads to use 
    !   -b, -begin:     Number of changes to start with
    !   -e, -end:       Number of changes to end with
    !   -s, -step:      Number of changes to step with between iteration 
    !   -m, -mode:      Type of changes to be made (1 - Add, 2 - Remove, 3 - Add and Remove)
    !   -n, -trials:    Number of trials to be conducted at each number of changes

    use data_types
    use NodeAffinity_Functions
    use Common_Functions
    use QW_Similarity 
    use Graph_Functions
    use omp_lib

    implicit none
    character(len=100) :: graphFile, outFile
    integer :: nSteps, nPhase, nThreads, startChanges, endChanges, stepChanges, nTrials, mode
    real(8) :: eps
    real(8), allocatable :: phi(:,:)
    
    ! Default values (change all via parameter file, or some via command line)       
    outFile = "Output/Out_Change.dat"
    graphFile = "Parameters/GraphList.dat"
    nSteps = 10
    nPhase = 2
    nThreads = 12
    startChanges = 50
    endChanges = 200
    stepChanges = 50
    nTrials = 5
    mode = 2
    eps = 0.001d0
    ! Initialise program before running
    call init()
    call run()

    deallocate(phi)

contains

    ! Subroutine to initialise program to localise variables
    subroutine init() 
        implicit none 
        integer :: i
        real(8) :: offset
        character(len=100) :: arg, file_parameters
        
        ! Check command line for inputs or parameter file
        do i = 1, iargc(), 2
            call getarg(i,arg)
            select case (arg)
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
                case ("-b", "-begin")
                    call getarg(i+1,arg)
                    read(arg,*)startChanges 
                case ("-e", "-end")
                    call getarg(i+1,arg)
                    read(arg,*)endChanges 
                case ("-n", "-trials")
                    call getarg(i+1,arg)
                    read(arg,*)nTrials 
                case ("-m", "-mode")
                    call getarg(i+1,arg)
                    read(arg,*)mode   
                case ("-s", "-step")
                    call getarg(i+1,arg)
                    read(arg,*)stepChanges 
            end select
        end do

        allocate(phi(nSteps,2))
        call init_random_seed()
        call omp_set_num_threads(nThreads)
        
        offset = 0.1d0

        ! Create a list of random phases to use in each step
        do i = 1, nSteps
            call random_number(phi(i,:))
            ! Scale random number between 0 and 2pi, to some offset (want to avoid 0, 2pi, ideally also multiples of pi/2 and pi)
            phi(i,:) = phi(i,:)*((2-offset)*pi - offset) + offset
        end do    

        ! Code is currently only set up to work with 0, 1 or 2 phase additions, only '2' detects isomorphism 
        if( nPhase < 0 .or. nPhase > 2 ) then 
            write(*,*)'Can only add phase to 0, 1 or 2 nodes. Setting number of nodes to 2'
            nPhase = 2
        else if( nPhase == 0) then
            ! Treat 0 phase additions as phi = 0 applied to 1 reference node
            ! This way, all combinations of nodes is tested, but without any phase added
            phi = 0     
            nPhase = 1
        end if

    end 

    subroutine run()
        implicit none
        integer, allocatable, dimension(:,:) :: AdjA, AdjA_Orig, AdjB
        integer :: cr, cm, start, finish, start2, finish2, n, readState, i, j
        integer :: numTest, nChanges, nEdges_Orig, nEdges!, sumAB, sumAA, sumBB
        real(8) :: rate, sim(4), NA_sim, total_time
        logical :: outFile_exists
        character(len=100) :: graph, strTmp
        character(len=200) :: line
    
        ! Number of graphs to be tested in total
        numTest = (endChanges - startChanges) / stepChanges + 1

        ! Initialise timer variables 
        call system_clock(count_rate=cr)
        call system_clock(count_max=cm)
        rate = real(cr)
        ! If the output file doesn't exist, include a header to show output format
        inquire(file = outFile, exist = outFile_exists)
        if( .not. outFile_exists) then 
            open(10, file=outFile, action='write')
            write(10,*)'The similarity values are output as such:'
            write(10,*)'1st column - number of changes'
            write(10,*)'2nd column - % difference'
            write(10,*)'3rd column - NA Similarity'
            write(10,*)'4th-7th column - QW using sim = 1 / (1+dx) for d1-d4 below'
            write(10,*)'  d1 - prob1Norm + prob2Norm > eps = 1'
            write(10,*)'  d2 - (prob1Norm + prob2Norm)/nSteps'
            write(10,*)'  d3 - (prob1Norm/nSteps)^2 + (prob2Norm/nSteps)^2'
            write(10,*)'  d4 - (prob1Mat_d/nSteps) + (prob2Mat_d/nSteps)'
        else 
            open(10, file=outFile, position='append', action='write')
        end if 

        open(20, file=graphFile, action='read')

        do 
            read(20,'(A)',iostat=readState)line
            if(readstate /= 0)exit         ! End of graph input file reached
            read(line,'(A)')graph          ! input for graph A location 
            print*,'Testing ',graph
            ! Find number of nodes by counting rows in file        
            call file_dimensions(graph, n)

            open(30, file=graph, iostat=readState, action='read')

            if( readstate /= 0 ) then ! Error occured when loading in graph
                write(*,*)'Cannot open ', trim(adjustl(graph)), ': ', readstate
                ! Close opened files 
                close(30)
                cycle   
            end if

            allocate( AdjA(n,n), AdjA_Orig(n,n), AdjB(n,n))

            ! Quantum walk code may modify adjacency matrix - thus keep a copy of the original
            read(30, *)AdjA_Orig
            close(30)

            ! sim(1) - prob1Norm + prob2Norm > eps = 1
            ! sim(2) - (prob1Norm + prob2Norm )/nSteps
            ! sim(3) - (prob1Norm/nSteps)^2 + (prob2Norm/nSteps)^2
            ! sim(4) - MatusitaDistance1 + MatusitaDistance2

            ! Start timing
            call system_clock(start)
            
            nEdges_Orig = Count_Edges( AdjA_Orig, .false.)
            
            write(10,*)'--------------------------------------------'
            write(10,*)trim(adjustl(graph)), ' vs changes'
            write(10,'(A, 3I3)')'Number of steps; Number of phase additions; Number of threads: ', nSteps, nPhase, nThreads 
            write(10,'(A,F5.3)')'eps: ', eps 
            
            do i = 1, numTest
                
                ! Loop over the number of trials to do at current edge number
                do j = 1, nTrials
                    call system_clock(start2)
                    ! Set graphs to be tested
                    AdjA = AdjA_Orig
                    AdjB = AdjA_Orig 
                    ! Change matrix B
                    nChanges = stepChanges * (i - 1) + startChanges 
                    call change_graph( AdjB, nChanges, .false., mode)
                    nEdges = Count_Edges( AdjB, .false.)
                    ! Make similarity measure
                    call QW_Calc_Sim( AdjA, AdjB, sim, nPhase, nSteps, phi, eps)
                    NA_sim =  NodeAffinity_Sim( AdjA, AdjB)
                    call system_clock(finish2)
                    total_time = (finish2 - start2)/rate
                    write(strTmp,'(A,I4,A,I3,A,F20.4)')'Time taken for ',nChanges,' changes trial ',j,': ', total_time
                    write(*,*)trim(adjustl(strTmp))
                    write(10,*)nChanges, 100.d0*(nEdges - nEdges_Orig)/nEdges_Orig, NA_sim, sim           
                    flush(10)   ! Flush to output results, in case of algorithm prematurely ending
                end do 
                call system_clock(finish)
            end do 
            ! Finish timing
            call system_clock(finish)
            total_time = (finish - start)/rate
            write(strTmp,'(A,F20.4)')'Time taken: ', total_time
            write(10,*)trim(adjustl(strTmp))
            deallocate(AdjA, AdjB, AdjA_Orig)

        end do

        close(10)
        close(20)

    end
    
    ! Read input from parameter file
    subroutine parse_input( inFile)
        implicit none 
        character(len=*), intent(in) :: inFile 
        character(len=100) :: line
        integer :: i, readState

        open(30, file=inFile)
        ! Index current input
        i = 1

        ! Loop until end of file
        do
            read(30, '(A)', iostat=readState)line
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
                    read(line,*)nPhase
                case (5)
                    read(line,*)nThreads
                case (6)
                    read(line,*)startChanges
                case (7)
                    read(line,*)endChanges
                case (8)
                    read(line,*)stepChanges
                case (9)
                    read(line,*)nTrials
                case (10)
                    read(line,*)mode
                case (11)
                    read(line,*)eps
            end select
            i = i + 1
        end do
        close(30)
    end 
end