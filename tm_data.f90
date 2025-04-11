module tm_data
    use vector
    implicit none
    
    
   
    public randvec2, encodePos, decodePos, openBin, readBin, generate_faces
    
    contains
    
    
    function randvec2 () result (res)
        type (vec2_f64) :: res
        call RANDOM_NUMBER(res%arr(1))
        call RANDOM_NUMBER(res%arr(2))
    end function
    
    
    
    ! compresses 3 20 bit reals in a 64 bit integer
    function encodePos (pos) result (enc)
        type (vec3_f64) :: pos
        integer(8) :: enc

        enc = INT(pos%arr(1) * 1000,8)
        enc = ISHFT(enc,20)
        enc = OR(enc, INT(pos%arr(2) * 1000,8))
        enc = ISHFT(enc,20)
        enc = OR(enc, INT(pos%arr(3) * 1000,8))     
    end function
    
    ! decompresses 3 20 bit reals from a 64 bit integer
    function decodePos (enc) result (pos)
        type (vec3_f64) :: pos
        integer(8) :: enc
        integer(4), parameter :: shiftConst = ISHFT(1,20) - 1
        
        pos%arr(3) = real(AND(enc,shiftConst),8) / 1000
        enc = ISHFT(enc,-20)
        pos%arr(2) = real(AND(enc,shiftConst),8) / 1000
        enc = ISHFT(enc,-20)
        pos%arr(1) = real(AND(enc,shiftConst),8) / 1000
    end function
    
    
    subroutine readBin (path,arr,iunit)
        real(8), allocatable, intent(inout) :: arr(:)
        character(len = *) :: path
        integer, intent(in) :: iunit
        integer :: ios, n
    
        call openBin(path,iunit)
        
        inquire(unit=iunit, size=n)
        n = n / 8  ! Since REAL*8 is 8 bytes
        allocate(arr(n))

        rewind(iunit)
        read(iunit) arr
        close(iunit)
    end subroutine
    
    subroutine openBin (path,iunit)
        integer, intent(in) :: iunit
        character(len = *), intent(in) :: path
        integer :: ios
        
        open (unit=iunit, file=path, form="unformatted", access="stream", status="unknown", iostat=ios)
        if (ios /= 0) then
            print *, 'Error opening file: ', path, ios
            stop
        end if
    end subroutine
    
    subroutine generate_faces (arr,nsize)
        integer(4), allocatable, intent(inout) :: arr(:)
        integer, intent(in) :: nsize
        integer :: b, arr_size, i
        
        arr_size = nsize*(nsize-1)-1
        allocate(arr(arr_size*6))
        
        b = 1
        do i = 1, nsize*(nsize-1)-1
            if (mod(i,nsize) /= nsize-1) then
                arr((b-1)*6+1) = i
                arr((b-1)*6+2) = i+1
                arr((b-1)*6+3) = i+nsize+1
                
                arr((b-1)*6+4) = i
                arr((b-1)*6+5) = i+nsize
                arr((b-1)*6+6) = i+nsize + 1
                b = b + 1
            end if
        end do
    end subroutine
    
    
    
end module tm_data