! =================================================================================================
!                                         TIME MODULE
! =================================================================================================


module time_module

    implicit none
    private
    
    integer,    parameter   :: dp = kind(1.0d0) !   
    integer                 :: clock0, clock1, clockrate
    
    public  :: start_clock, read_clock
    
    contains
    
    subroutine start_clock()
    
        implicit none
        
        call system_clock(clock0, count_rate = clockrate)        
        
    end subroutine start_clock
    
    subroutine read_clock(time)
    
        implicit none
        real(dp), intent(out) :: time
        
        call system_clock(clock1, count_rate = clockrate)
        time = real(clock1-clock0,dp) / real(clockrate,dp)
        
        write(*,'(A,ES10.2)') 'Seconds elapsed: ', time
    
    end subroutine read_clock

end module time_module