program main
    use tm_data
    use convolution
    use delaunay
    use, intrinsic :: iso_c_binding
    implicit none
    character(len = *), parameter :: in_file = "C:\Temp\Mediafunction_hpt_400_400.bin"
    integer :: i, j
    integer, parameter :: nsize = 401
    real(8), allocatable :: data(:), datad1(:,:), datad1x(:,:), datad1y(:,:)
    integer(4), allocatable :: face_data(:)
    type(mesh) :: del
    
    
    allocate(datad1x(nsize,nsize), datad1y(nsize,nsize), datad1(nsize,nsize))
    
    call readBin(in_file, data,9)
    
    do i = 1, nsize
        do j = 1, nsize
            !datad1(i,j) = data((nsize*(i-1)+j)*3)
        end do
    end do
    
    call convolute_real(datad1,scharr_operator_x,datad1x)
    call convolute_real(datad1,scharr_operator_y,datad1y)
    
   
    do i = 1, nsize
        do j = 1, nsize
            data((nsize*(i-1)+j)*3) = sqrt(datad1y(i,j)**2 + datad1y(i,j)**2)
        end do
    end do
    
    call init(del, vec3f(0,0,0,"BL"), vec3f(0,1,0,"TL"), vec3f(1,1,0,"TR"), vec3f(1,0,0,"BR"))
    
    
    
    
    call insert_site(del,vec3f(0.4,0.6,1,"AH"))
    
    call quad_enum(del%root, v_proc, c_null_ptr)
    
    
    
    
    
    
    call generate_faces(face_data,nsize)
    
    call openBin("faces_400_400.bin",10)
    write(10) face_data
    close(10)
    deallocate(face_data)
    
    call openBin("vertex_400_400.bin",11)
    write(11) data(1:nsize*nsize*3+1)
    close(11)
    deallocate(data)
end program