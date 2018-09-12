program QW_Multiple
    ! Author: Callum Schofield - 2017
    ! This code makes straight comparisons between the graphs given in 1 input file 
    !   Each line of the graph input file contains the filepath to the graphs to be tested separated by '<>'
    !   Parameters can also be specified via an input file, if a file is not specified then default values are used
    
    ! Optional command line arguments (note not all arguments contained in parameters file can be entered on command line 
    !   -p, -param:   Parameter file location. NOTE use this first if entering other parameters, as it will overwrite any already loaded parameters
    !   -g, -graphs:  Graph file location
    !   -o, -out:     Output file location 
    !   -t, -threads: Number of threads to use 

    ! use Walk_Functions
    use data_types
    use NodeAffinity_Functions
    use Common_Functions
    use QW_Similarity 
    use omp_lib
    
    implicit none
    character(len=100) :: graphFile, outFile
    integer :: nSteps, nPhase, nThreads
    real(8) :: eps
    real(8), allocatable :: phi(:,:)
 
    ! Default values (change all via parameter file, or some via command line)   
    outFile = "Output/Out_Multiple.dat"
    graphFile = "Parameters/GraphListComparison.dat"
    nSteps = 20
    nPhase = 2
    nThreads = 18
    eps = 0.000001d0

    call init()
    call run()

    deallocate(phi)

contains

    ! Subroutine to initialise program to localise variables
    subroutine init() 
        implicit none 
        integer :: i 
        ! real(8) :: offset
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
            end select
        end do

        allocate(phi(nSteps,2))
        call init_random_seed()
        call omp_set_num_threads(nThreads)
        
        !offset = 0.1d0

        ! Create a list of random phases to use in each step
        !do i = 1, nSteps
        !    call random_number(phi(i,:))
        !    ! Scale random number between 0 and 2pi, to some offset (want to avoid 0, 2pi, ideally also multiples of pi/2 and pi)
        !    phi(i,:) = phi(i,:)*((2-offset)*pi - offset) + offset
        !end do 
           
        ! Phi is still set-up to allow varying phases for each step, however for consistency between iterations a constant phase should be applied
        ! For optimisation sake, this should be reduced to a single value if the above functionality is deemed unsuitable
        phi(:,1) = pi/2
        phi(:,2) = 3*pi/2

        ! Code is currently only set up to work with 0, 1 or 2 phase additions, only '2' detects isomorphism 
        if( nPhase < 0 .or. nPhase > 2 ) then 
            ! Code is currently only set up to work with 0, 1 or 2 phase additions
            write(*,*)'Can only add phase to 0, 1 or 2 nodes. Setting number of reference nodes to 2'
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
        integer, allocatable, dimension(:,:) :: AdjA, AdjB
        integer :: cr, cm, start, finish, n1, n2, i, readStateA, readStateB!, sumAB, sumAA, sumBB
        real(8) :: rate, sim(4), NA_sim, total_time
        logical :: outFile_exists
        character(len=50) :: strTmp
        character(len=100) :: graphA, graphB
        character(len=200) :: line

        ! Initialise timer variables 
        call system_clock(count_rate=cr)
        call system_clock(count_max=cm)
        rate = real(cr)
        ! If the output file doesn't exist, include a header to show output format
        inquire(file = outFile, exist = outFile_exists)
        if( .not. outFile_exists) then 
            open(10, file=outFile, action='write')
            write(10,*)'The similarity values are output as such:'
            write(10,*)'First line - node affinity similarity.'
            write(10,*)'Second line - QW using sim = 1 / (1+dx) where each column dx is given by the sum of:'
            write(10,*)'  d1 - prob1Norm + prob2Norm > eps = 1'
            write(10,*)'  d2 - (prob1Norm + prob2Norm)/nSteps'
            write(10,*)'  d3 - (prob1Norm/nSteps)^2 + (prob2Norm/nSteps)^2'
            write(10,*)'  d4 - (prob1Mat_d/nSteps) + (prob2Mat_d/nSteps)'
        else 
            open(10, file=outFile, position='append', action='write')
        end if 

        open(20, file=graphFile, action='read')

        do 
            read(20,'(A)',iostat=readStateA)line
            if(readstateA /= 0)exit         ! End of graph input file reach
            i = index(line, "<>")           ! '<>' acts as a separater between graph pathsed
            if(i == 0 )cycle                ! Incorrect input/no pair on this line  
            read(line(1:i-1),'(A)')graphA   ! input for graph A location 
            read(line(i+2:),'(A)')graphB    ! input for graph B location

            print*,'Comparing ', trim(adjustl(graphA)), ' vs ', trim(adjustl(graphB))

            ! Allow same graph to be entered for self comparison
            if(graphA == graphB)then
                ! Find number of nodes by counting rows in file        
                call file_dimensions(graphA, n1)
                open(30, file=graphA, iostat=readStateA, action='read')

                if( readstateA /= 0 ) then ! Error occured when loading in graph
                    write(*,*)'Cannot open ', trim(adjustl(graphA)), ': ', readstateA
                    ! Close opened files 
                    close(30)
                    cycle   
                end if
                allocate( AdjA(n1,n1), AdjB(n1,n1))
                read(30, *)AdjA
                AdjB = AdjA
                close(30)
            else
                ! Find number of nodes by counting rows in file        
                call file_dimensions(graphA, n1)
                call file_dimensions(graphB, n2)
                
                open(30, file=graphA, iostat=readStateA, action='read')
                open(31, file=graphB, iostat=readStateB, action='read')

                if( readstateA /= 0 .or. readStateB /= 0) then ! Error occured when loading in either graph
                    write(*,*)'Cannot open ', trim(adjustl(graphA)), ': ', readstateA, ' or ', &
                        & trim(adjustl(graphB)), ': ', readStateB
                    ! Close opened files 
                    close(30)
                    close(31)
                    cycle   
                end if

                if(n1 /= n2)then 
                    write(*,*)'Graphs have mismatched number of nodes. &
                    & Code is currently only setup to compare graphs with the same number of nodes.'
                    cycle
                end if

                allocate( AdjA(n1,n1), AdjB(n2,n2))
                read(30, *)AdjA
                read(31, *)AdjB
                close(30)
                close(31)

            end if
     
            ! sim(1) - prob1Norm + prob2Norm > eps = 1
            ! sim(2) - (prob1Norm + prob2Norm)/nSteps
            ! sim(3) - (prob1Norm/nSteps)^2 + (prob2Norm/nSteps)^2
            ! sim(4) - MatusitaDistance1 + MatusitaDistance2

            ! Start timing
            call system_clock(start)

            !call QW_Calc_Sim( AdjA, AdjB, sim, nPhase, nSteps, phi, eps)

            ! Finish timing
            call system_clock(finish)
            total_time = (finish - start)/rate
            write(strTmp,'(A,F10.4)')'Time taken: ', total_time

            NA_sim =  NodeAffinity_Sim( AdjA, AdjB)
            deallocate(AdjA, AdjB)

            write(10,*)'--------------------------------------------'
            write(10,*)trim(adjustl(graphA)), ' vs ', trim(adjustl(graphB))
            write(10,*)trim(adjustl(strTmp))
            write(10,'(A, 3I3)')'Number of steps; Number of phase additions; Number of threads: ', nSteps, nPhase, nThreads 
            write(10,'(A,F10.7,F10.7)')'theta, phi: ', phi(1,1), phi(1,2) 
            write(10,'(A,F15.10)')'eps: ', eps 
            write(10,*)NA_sim 
            ! write(10,*)sim
            flush(10)

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
                    read(line,*)nPhase
                case (5)
                    read(line,*)nThreads
                case (6)
                    read(line,*)eps
            end select
            i = i + 1
        end do
        close(50)
    end 
end