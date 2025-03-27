module convolution
    implicit none
    
    real(8), dimension(2,2),parameter :: roberts_operator_x = TRANSPOSE(RESHAPE((/1,0,0,-1/),(/2,2/)))
    real(8), dimension(2,2),parameter :: roberts_operator_y = TRANSPOSE(RESHAPE((/0,1,-1,0/),(/2,2/)))
    
    real(8), dimension(3,3), parameter :: sobel_operator_x = RESHAPE((/-1,0,1, -2,0,2, -1,0,1/), (/3,3/))
    real(8), dimension(3,3), parameter :: sobel_operator_y = RESHAPE((/-1,-2,-1, 0, 0, 0, 1,2, 1/), (/3,3/))
    
    real(8), dimension(3,3), parameter :: scharr_operator_x = RESHAPE((/3,0,-3, 10,0,-10, 3,0,-3/), (/3,3/))
    real(8), dimension(3,3), parameter :: scharr_operator_y = RESHAPE((/-3,-10, -3, 0, 0, 0, 3,10, 3/), (/3,3/))
    
    public convolute_real
    !public roberts_operator_x robertsoperator_y
    !public sobel_operator_x, sobel_operator_y
    !public scharr_operator_x, scharr_operator_y
    
    contains 
    
    
    
    
    
    
    subroutine convolute_real (input, kernel, output)
        real(8), intent(in) :: input(:,:), kernel(:,:)
        real(8), intent(inout) :: output(:,:)
        integer :: ix, iy, kx, ky, ox, oy, i, j, m, n
        ix = SIZE(input,1)
        iy = SIZE(input,2)
        ox = SIZE(output,1)
        oy = SIZE(output,2)
        kx = SIZE(kernel, 1)
        ky = SIZE(kernel,2)
        
        if (.NOT. (ix == ox .AND. iy == oy)) then
            print *, "Input and Output matrix dont have the same dimension: input(",ix,",",iy,"), output(",ox,",",oy,")"
            error stop
        end if
        
        do i = 2, ix-1
            do j = 2, iy-1
                do m = -1, 1
                    do n = -1, 1
                        output(i, j) = output(i, j) + input(i+m, j+n) * kernel(m+2, n+2)
                    end do
                end do
            end do
        end do
    end subroutine
    
    
    
    
    
    
    
    
    
    
    
    
    
end module
    
    