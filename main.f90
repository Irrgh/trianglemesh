program main
    use tm_data
    use convolution
    use delaunay
    use vector
    use euler
    use bvh_acceleration
    use surface_functions
    use, intrinsic :: iso_c_binding
    implicit none
    character(len = *), parameter :: in_file = "C:\Temp\Mediafunction_hpt_50_50.bin"
    integer :: i, j
    type(tm_del) :: del
    type(tm_mesh) :: m
    type(tm_bvh) :: bvh
    real(8) :: x,y,z, sum, res, approx
    
    call init(del, 100.0_8,500000)
    
   
    call RANDOM_SEED()
    do i = 1, 500000
        call RANDOM_NUMBER(x)
        call RANDOM_NUMBER(y)
        x = x*120.0-60.0
        y = y*120.0-60.0
        if (MOD(i,1000) == 0) then
            print *, i
        end if
        z = test1(x,y)
        call insert_site(del, vec3_f64_create(x,y,z))
    end do
    
    call finalize(del)
    print *, del%vc, del%ec, euler_faces(del%vc, del%ec)
    
    m = get_mesh(del)
    call quad_bvh_create(bvh,m)
    
    
    
    
    sum = 0.0
    do i = 1, 1e7
    
        call RANDOM_NUMBER(x)
        call RANDOM_NUMBER(y)
        x = x*100.0-50.0
        y = y*100.0-50.0
        
        approx = intersection_bvh(bvh,ray(vec3_f64((/x,y, 100.0_8/)), vec3_f64((/0.0_8,0.0_8,-1.0_8/))))
        res = 100.0 - test1(x,y)
        
        
        
        if (approx < 1e30) then
            sum = sum + abs(res - approx)
        else
            print *, "No intersection at: ", x, y
        end if
        
        if (mod(i,100000) == 0) print *, i
    end do
    
    
    print *, "Average error for 10m sampled points with a base resolution of 500k points: ", sum / 1e7
    
    
    print *, "Mesh created successfully"
 
    call openBin("delaunay_faces.bin",10)
    write(10) m%faces
    
    close(10)
    
    call openBin("delaunay_vertex.bin",11)
    write(11) m%vertices
    close(11)
    
      
    
end program