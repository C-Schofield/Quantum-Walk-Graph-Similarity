program Szegedy_Multiple
    ! Author: Callum Schofield - 2017
    ! This code makes comparisons between given input graph pairs
    ! The code allows for MPI - however this would only be beneficial for a large number of graphs
    !   note that the first node only allocates computation to other nodes + receives results
    !
    ! Inputs  - nThreads:   number of threads to be used (must be specified in slurm file for Pawsey, not here)          
    !           graphFile:    list containing adjacency matrices to be tested (all matrices must be in the same folder as program)
    !           nSteps:     number of steps to be done
    !           errMult:    error multiplier to be used in Szegedy similarity index classification.
    !           
    ! Outputs - Output.dat: data file containing similarity index (both classical and quantum) for each graph compared 
    !
    ! Optional command line arguments (note not all arguments contained in parameters file can be entered on command line 
    !   -p, -param:   Parameter file location. NOTE use this first if entering other parameters, as it will overwrite any already loaded parameters
    !   -g, -graphs:  Graph file location
    !   -o, -out:     Output file location 
    !   -t, -threads: Number of threads to use
    !   -r, -errmult: Error multiplier to be used in Szegedy similarity index classification. Input a value <= 0 if you want the algorithm to return
    !                 only d, that is 1 / (1+ w*d) can be determined post completion for varying w.

    use m_refsor
    use omp_lib
    use Szegedy_Similarity
    use Common_Functions
    use Graph_Functions
    use Szegedy_Functions
    use NodeAffinity_Functions

    implicit none

    integer :: nSteps, nThreads, cr, cm, numGraphs, start, finish!, ierr, task_id, num_tasks
    real(8) :: errMult, rate
    character(len = 100) :: arg, graphFile, OutFile

    ! Set timer function
    call system_clock(start)
    call system_clock(count_rate=cr)
    call system_clock(count_max=cm)
    rate = real(cr)

    ! Initialise parameters to default
    nSteps = 25                                         ! Number of steps of Szegedy to do per iteration
    graphFile = 'Parameters/GraphListComparison.dat'    ! Adjacency matrix input file
    OutFile = 'Output/output.dat'                       ! Output file
    nThreads = 12                                       ! Number of threads to use
    errMult = 1.d0                                      ! Error multiplier 

    call init()

    call file_dimensions( graphFile, numGraphs) 

    call run_prog()

    call system_clock(finish)

    write(*,*)'The program took ', (finish - start), 'seconds.'

contains

    ! Subroutine to initialise program 
    subroutine init() 
        implicit none 
        integer :: i 
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
                case ("-w", "-weight")
                    call getarg(i+1,arg)
                    read(arg,*)errMult 
            end select
        end do

        call init_random_seed()
        call omp_set_num_threads(nThreads)

    end 

    ! Run code 
    subroutine run_prog()
        implicit none
        real(8), allocatable, dimension(:,:) :: A1, A2
        real(8) :: totalTime, simIndex(3), simBelief
        integer :: n, n2, i, j, k, nTest, tmp, change_err, readStateA, readStateB, id, time1, time2
        logical :: outFile_exists
        integer(kind=omp_lock_kind) :: lck1, lck2
        character(len = 50) :: strTmp
        character(len=100)  :: graphA, graphB
        character(len=200)  :: line

        call omp_init_lock(lck1)
        call omp_init_lock(lck2)
        ! Open list of adjacency matrices
        open(unit=10, file=graphFile, action='read')
        ! Open output file
        inquire(file = outFile, exist = outFile_exists)
        if( .not. outFile_exists) then 
            open(20, file=outFile, action='write')
            write(20,*)'The similarity values are output as such:'
            write(20,*)'  Similarity index using Szegedy: abs(prob1 - prob2), d; sim = 1 / (1+d)'
            write(20,*)'  Similarity index using Szegedy: 2Norm(prob1, prob2), d; sim = 1 / (1+d)'
            write(20,*)'  Similarity index using Szegedy: Matusita(prob1, prob2), d; sim = 1 / (1+d)'
            write(20,*)'  Similarity index using Node Affinity'
        else 
            open(20, file=outFile, position='append', action='write')
        end if 

        !$omp parallel private(time1, time2, totalTime, line, graphA, graphB, A1, A2, n, n2, & 
        !$omp&                 readStateA, readStateB, simIndex, simBelief) shared(lck1, lck2, numGraphs)
        !$omp do !schedule(static)
            do i = 1, numGraphs
                call system_clock(time1)
                ! time1 = MPI_Wtime()
                call omp_set_lock(lck1)
                read(10,'(A)',iostat=readStateA)line
                ! if(readstateA /= 0)exit          ! End of graph input file reach
                id = index(line, "<>")           ! '<>' acts as a separater between graph paths
                if(id == 0 ) then                ! Incorrect input/no pair on this line  
                    write(*,*)'Incorrect input in: ', line
                    call omp_unset_lock(lck1)
                    cycle                        
                end if
                read(line(1:id-1),'(A)')graphA   ! input for graph A location 
                read(line(id+2:),'(A)')graphB    ! input for graph B location
                !call omp_unset_lock(lck1)
                write(*,*)'Comparing ', trim(adjustl(graphA)), ' vs ', trim(adjustl(graphB))
                ! Lock opening/reading graphs in case same graph file is being read
                ! Allow same graph to be entered for self comparison
                if(graphA == graphB)then
                    ! Find number of nodes by counting rows in file        
                    call file_dimensions(graphA, n, i)
                    open(30+i, file=graphA, iostat=readStateA, action='read')

                    if( readstateA /= 0 ) then ! Error occured when loading in graph
                        write(*,*)'Cannot open ', trim(adjustl(graphA)), ': ', readstateA
                        ! Close opened files 
                        close(30+i)
                        call omp_unset_lock(lck1)
                        cycle   
                    end if
                    allocate( A1(n,n), A2(n,n))
                    read(30+i, *)A1
                    A2 = A1
                    close(30+i)
                else
                    ! Find number of nodes by counting rows in file        
                    call file_dimensions(graphA, n, 60 + i)
                    call file_dimensions(graphB, n2, 70 + i)
                    
                    open(30+i, file=graphA, iostat=readStateA, action='read')
                    open(40+i, file=graphB, iostat=readStateB, action='read')

                    if( readstateA /= 0 .or. readStateB /= 0) then ! Error occured when loading in either graph
                        write(*,*)'Cannot open ', trim(adjustl(graphA)), ': ', readstateA, ' or ', &
                            & trim(adjustl(graphB)), ': ', readStateB
                        ! Close opened files 
                        close(30+i)
                        close(40+i)
                        call omp_unset_lock(lck1)
                        cycle   
                    end if

                    if(n /= n2)then 
                        write(*,*)trim(adjustl(graphA)), ' and ', trim(adjustl(graphB)), &
                        & ': Graphs have mismatched number of nodes. '!, &
                                ! & 'Code is currently only setup to compare graphs with the same number of nodes.'
                        call omp_unset_lock(lck1)
                        cycle
                    end if
                    allocate( A1(n,n), A2(n2,n2))
                    read(30+i, *)A1
                    read(40+i, *)A2
                    close(30+i)
                    close(40+i)
                end if
                call omp_unset_lock(lck1)

                ! Fortran writes column major - so transpose matrix
                A1 = transpose(A1)
                A2 = transpose(A2)
                ! nEdges = Count_Edges(A1, .true.)
                simIndex = Szegedy_Sim( A1, A2, nSteps, errMult)
                simBelief = NodeAffinity_Sim( A1, A2)
                ! time2 = MPI_Wtime()
                ! totalTime = time2 - time1
                call system_clock(time2)
                totalTime = (time2 - time1)/rate
                call omp_set_lock(lck2)
                write(strTmp,'(A,F10.4)')'Execution time: ', totalTime
                write(20,*)'-----------------------------------'
                write(20,*)trim(adjustl(graphA)), ' vs ', trim(adjustl(graphB))
                write(20,'(A, I3)')'Number of vertices: ', n
                write(20,'(A, I3)')'Number of steps: ', nSteps
                write(20,'(A, F10.4)')'Error multiplier: ', errMult
                write(20,'(A, I3)')'Threads: ', nThreads
                write(20,'(A)')trim(adjustl(strTmp))
                write(20,*)simIndex(1)
                write(20,*)simIndex(2)
                write(20,*)simIndex(3)
                write(20,*)simBelief
                flush(20)
                call omp_unset_lock(lck2)
                deallocate( A1, A2)
            end do
        !$omp end do 
        !$omp end parallel 
        close(10)
        close(20)   
        call omp_destroy_lock(lck1)
        call omp_destroy_lock(lck2)
    end subroutine

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
                    read(line,*)nThreads
                case (5)
                    read(line,*)errMult
            end select
            i = i + 1
        end do
        close(50)
    end 

end program Szegedy_Multiple
