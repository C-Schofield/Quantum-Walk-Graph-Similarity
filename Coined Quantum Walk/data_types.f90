module data_types
    ! Author: Callum Schofield - 2017
    ! Defines parameters and types for coined quantum walk.
    implicit none

    real(8), parameter      :: pi = 4.d0*ATan(1.d0) ! Approximation for pi
    complex(8), parameter   :: ii = (0.d0,1.d0)     ! Complex i

    private 
    public :: pi, ii, node

    type node
        real(8), allocatable, dimension(:,:)    :: coin     ! Contains coin operator for each node, variable length depends on degree 
        complex(8), allocatable, dimension(:)   :: space    ! Contains state space information for each node, variable length depends on degree
    end type 

end module data_types