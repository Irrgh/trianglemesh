module surface_functions
    implicit none
    contains
    
    pure elemental function test1 (x,y) result(z)
        real(8), intent(in) :: x,y    
        real(8) :: z
        z = sin(0.5*(x**2+y**2)**0.4) * sin(ATAN2(x,y)*5.0) * 2.0 + x / 2.0
    end function
    
    pure elemental function test2(x,y) result(z)
        real(8), intent(in) :: x,y
        real(8) :: z
        z = sin(x*0.2+y*0.2)
    end function
    
end module