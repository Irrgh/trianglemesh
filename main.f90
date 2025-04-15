program main
    use tm_data
    use convolution
    use delaunay
    use vector
    use euler
    use quadtree
    use, intrinsic :: iso_c_binding
    implicit none
    character(len = *), parameter :: in_file = "C:\Temp\Mediafunction_hpt_50_50.bin"
    integer :: i, j
    type(tm_del) :: del
    type(mesh) :: m
    type(quad_bvh) :: bvh
    real(8) :: x,y,z
    
    call init(del, 100.0_8,1000)
    
   
    call RANDOM_SEED()
    do i = 1, 1000
        call RANDOM_NUMBER(x)
        call RANDOM_NUMBER(y)
        x = x*120-60
        y = y*120-60
        if (MOD(i,100) == 0) then
            print *, i
        end if
        z = sin(0.5*(x**2+y**2)**0.4) * sin(ATAN2(x,y)*5) * 2 + x / 10
        call insert_site(del, vec3_f64_create(x,y,z))
    end do
    
    call finalize(del)
    
    
    print *, del%vc, del%ec, euler_faces(del%vc, del%ec)
    
    
    m = get_mesh(del)
    
    call quad_bvh_create(bvh,m)
    
    j = 0
    do i=1, bvh%nodes_used
        
        if (bvh%nodes(i)%prim_count > 0) then
            print *, i, bvh%nodes(i)%first_pc, bvh%nodes(i)%prim_count
            j = j + bvh%nodes(i)%prim_count
        end if
            
    end do
    
    
    
    print *, "Mesh created successfully"
 
    call openBin("delaunay_faces.bin",10)
    write(10) m%faces
    
    close(10)
    
    call openBin("delaunay_vertex.bin",11)
    
    write(11) m%vertices
    close(11)
end program