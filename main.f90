program main
    use tm_data
    use convolution
    use delaunay
    use vector
    use euler
    use, intrinsic :: iso_c_binding
    implicit none
    character(len = *), parameter :: in_file = "C:\Temp\Mediafunction_hpt_50_50.bin"
    integer :: i, j
    type(tm_del) :: del
    type(mesh) :: m
    real(8) :: x,y
    
    call init(del, 100.0,50000)
    
   
    call RANDOM_SEED()
    do i = 1, 50000
        call RANDOM_NUMBER(x)
        call RANDOM_NUMBER(y)
        x = x*120-60
        y = y*120-60
        if (MOD(i,100) == 0) then
            print *, i
        end if
        
        call insert_site(del, vec3_f64(x,y,(SIN(x*0.2) + COS(y*0.2))*4))
    end do
    
    call finalize(del)
    
    
    print *, del%vc, del%ec, euler_faces(del%vc, del%ec)
    
    
    m = get_mesh(del)
    print *, "Mesh created successfully"
 
    call openBin("delaunay_faces.bin",10)
    write(10) m%faces
    
    close(10)
    
    call openBin("delaunay_vertex.bin",11)
    
    write(11) m%vertices
    close(11)
end program