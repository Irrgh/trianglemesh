program main
    use tm_data
    use convolution
    use delaunay
    use mesh_helper
    use, intrinsic :: iso_c_binding
    implicit none
    character(len = *), parameter :: in_file = "C:\Temp\Mediafunction_hpt_400_400.bin"
    integer :: i, j
    integer, parameter :: nsize = 401
    real(8), allocatable :: data(:), datad1(:,:), datad1x(:,:), datad1y(:,:)
    integer(4), allocatable :: face_data(:)
    type(tm_del) :: del
    type(mesh) :: m
    real(8) :: x,y
    
    !allocate(datad1x(nsize,nsize), datad1y(nsize,nsize), datad1(nsize,nsize))
    !
    !call readBin(in_file, data,9)
    !
    !print *, loc(data), loc(data(1)), loc(data(2))
    !
    !do i = 1, nsize
    !    do j = 1, nsize
    !        !datad1(i,j) = data((nsize*(i-1)+j)*3)
    !    end do
    !end do
    !
    !call convolute_real(datad1,scharr_operator_x,datad1x)
    !call convolute_real(datad1,scharr_operator_y,datad1y)
    !
    !
    !do i = 1, nsize
    !    do j = 1, nsize
    !        data((nsize*(i-1)+j)*3) = sqrt(datad1y(i,j)**2 + datad1y(i,j)**2)
    !    end do
    !end do
    
    call init(del, 100.0,1100)
    
    !call insert_site(del,vec3f(0,-25,10,"01"))
    !call insert_site(del,vec3f(-35,0,-5,"02"))
    !call insert_site(del,vec3f(0,-74.4,16,"03"))
    !call insert_site(del,vec3f(-34.25,54.34,64.5,"04"))
    
    call RANDOM_SEED()
    do i = 1, 1000
        call RANDOM_NUMBER(x)
        call RANDOM_NUMBER(y)
        x = x*120-60
        y = y*120-60
        if (MOD(i,10) == 0) then
            print *, i
        end if
        
        
        call insert_site(del, vec3f(x,y,(SIN(x*0.2) + COS(y*0.2))*2,"__"))
    end do
    
    
    !call quad_enum(del%root,adjacency_list,c_null_ptr)
    
    print *, del%vc, del%ec, euler_faces(del%vc, del%ec)
    
    m = get_mesh(del)
    
    !allocate(data(3*del%vc))
    !
    !!call generate_faces(face_data,nsize)
    !
    !
    call openBin("delaunay_faces.bin",10)
    write(10) m%faces
    close(10)
    !!deallocate(face_data)
    !
    call openBin("delaunay_vertex.bin",11)
    write(11) m%vertices
    close(11)
    !deallocate(data)
end program