program QW_Single

    ! use Walk_Functions
    use data_types
    use NodeAffinity_Functions
    use Common_Functions
    use QW_Similarity 

    implicit none
    character(len=100) :: fileA, fileB, outFile
    integer :: nSteps, nPhase, nThreads
    real(8) :: eps
    real(8), allocatable :: phi(:,:)

    nSteps = -1
    ! Initialise program before running
    call init()
    if( nSteps == -1) then 
        write(*,*)'No parameter file found. Program stopping.'
        stop
    end if
    call run()

    deallocate(phi)

contains

    ! Subroutine to initialise program to localise variables
    subroutine init() 
        implicit none 
        integer :: i 
        ! real(8) :: offset
        character(len=100) :: arg, file_parameters

        ! Default parameters file 
        file_parameters = 'Parameters/default_parameters_single.dat'

        ! Check to see if user inputs a parameters file
        if( iargc() > 0) then
            call getarg(1,arg)
            read(arg,*)file_parameters 
        end if

        ! Initialise parameters using input file
        call parse_input(file_parameters)

        allocate(phi(nSteps,2))
        call init_random_seed()
        call omp_set_num_threads(nThreads)

        offset = 0.1d0

        ! Create a list of random phases to use in each step
        ! do i = 1, nSteps
        !     call random_number(phi(i,:))
        !     ! Scale random number between 0 and 2pi, to some offset (want to avoid 0, 2pi, ideally also multiples of pi/2 and pi)
        !     phi(i,:) = phi(i,:)*((2-offset)*pi - offset) + offset
        ! end do    

        ! Phi is still set-up to allow varying phases for each step, however for consistency between iterations a constant phase should be applied
        ! For optimisation sake, this should be reduced to a single value if the above functionality is deemed unsuitable
        phi(:,1) = pi/2
        phi(:,2) = 3*pi/2

        ! Code is currently only set up to work with 0, 1 or 2 phase additions
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
        integer, allocatable, dimension(:,:) :: AdjA, AdjB
        integer :: cr, cm, start, finish, n1, n2, i!, sumAB, sumAA, sumBB
        real(8) :: rate, sim(5), NA_sim, total_time
        character(len=100) :: strTmp
        logical :: outFile_exists

        ! Find number of nodes by counting rows in file        
        call file_dimensions(fileA, n1)
        call file_dimensions(fileB, n2)

        if(n1 /= n2)then 
            write(*,*)'Graphs have mismatched number of nodes. &
            & Code is currently only setup to compare graphs with the same number of nodes.'
            stop
        end if
        allocate( AdjA(n1,n1), AdjB(n2,n2))
        open(20, file=fileA, action='read')
        open(21, file=fileB, action='read')
        read(20, *)AdjA
        read(21, *)AdjB
        close(20)
        close(21)

        ! sim(1) - prob1Norm + prob2Norm > eps = 1
        ! sim(2) - prob1Norm + prob2Norm
        ! sim(3) - prob1Norm^2 + prob2Norm^2
        ! sim(4) - bin prob1Norm/nSteps, prob2Norm/nSteps to nearest 0.1, then add
        ! sim(5) - bin prob1Norm/nSteps, prob2Norm/nSteps to nearest 0.1, square, then add

        ! Currently sim is returned as d, actual sim is 1 / (1+d) or some variation thereof

        ! Initialise timer variables and start timing
        call system_clock(count_rate=cr)
        call system_clock(count_max=cm)
        rate = real(cr)
        call system_clock(start)

        call QW_Find_Sim( AdjA, AdjB, sim, nPhase, nSteps, phi, eps)

        ! Finish timing
        call system_clock(finish)

        total_time = (finish - start)/rate
        write(strTmp,'(A,F10.4)')'Time taken: ', total_time

        NA_sim =  NodeAffinity_Sim( AdjA, AdjB)

        ! If the output file doesn't exist, include a header to show output format
        inquire(file = outFile, exist = outFile_exists)
        if( .not. outFile_exists) then 
            open(40, file=outFile, action='write')
            write(40,*)'The similarity values are output as such:'
            write(40,*)'First line - node affinity similarity.'
            write(40,*)'Second line - QW using sim = 1 / (1+d) where each column d is given by the sum of:'
            write(40,*)'  1 - prob1Norm + prob2Norm > eps = 1'
            write(40,*)'  2 - prob1Norm + prob2Norm'
            write(40,*)'  3 - prob1Norm^2 + prob2Norm^2'
            write(40,*)'  4 - bin prob1Norm/nSteps, prob2Norm/nSteps to nearest 0.1'
            write(40,*)'  5 - bin prob1Norm/nSteps, prob2Norm/nSteps to nearest 0.1, square'
            write(40,*)'Third line - QW using sim = 1 / (1+d^2) where each column d is as above'
        else 
            open(40, file=outFile, position='append', action='write')
        end if 
        
        write(40,*)'--------------------------------------------'
        write(40,*)trim(adjustl(fileA)), ' vs ', trim(adjustl(fileB))
        write(40,*)trim(adjustl(strTmp))
        write(40,'(A, 3I3)')'Number of steps; Number of phase additions; Number of threads: ', nSteps, nPhase, nThreads 
        write(40,'(A,F5.3)')'eps: ', eps 
        write(40,*)NA_sim 
        write(40,*)1 / (1 + sim)
        write(40,*)1 / (1 + sim**2)

        close(40)
        deallocate( AdjA, AdjB)

    end

    subroutine file_dimensions(file, n)
        implicit none
        integer, intent(out) :: n 
        character(len=*), intent(in) :: file
        integer :: readState

        open(10, file=file, action='read')
        n = 0
        ! Loop until end of file
        do 
            read(10,*, iostat=readState)
            ! If iostat is <0 then end of file reached
            if( readState < 0) exit
            n = n + 1
        end do
        close(10)
    end

    subroutine parse_input( inFile)
        implicit none 
        character(len=*), intent(in) :: inFile 
        character(len=100) :: data
        integer :: i, readState

        open(30, file=inFile)
        ! Index current input
        i = 1

        ! Loop until end of file
        do
            read(30, '(A)', iostat=readState)data
            ! If iostat is <0 then end of file reached
            if( readState < 0) exit
            ! If line starts with comment character, cycle to next line
            if ( index(data, "#") == 1 )  cycle
            select case (i)
                case (1)
                    read(data,'(A)')fileA 
                case (2)
                    read(data,'(A)')fileB  
                case (3)
                    read(data,'(A)')outFile
                case (4)
                    read(data,*)nSteps
                case (5)
                    read(data,*)nPhase
                case (6)
                    read(data,*)nThreads
                case (7)
                    read(data,*)eps
            end select
            i = i + 1
        end do
        close(30)
    end 
end