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
    type(vec3f), allocatable :: point_list(:)
    
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
    
    call init(del, 100.0,5050)
         
    allocate(point_list(5000))
    
    !call insert_site(del, vec3f(0.398373349505228,-2.61666070782412,1.89149817489868,"__"))
    !call insert_site(del, vec3f(57.2057069477088,-57.9972388175155,-0.669062087691127,"__"))
    !call insert_site(del, vec3f(-16.2936447211355,-12.8600001303014,-1.45050420812468,"__"))
    !call insert_site(del, vec3f(-18.2886812531101,-24.0702092768475,1.19001561342115, "__"))
    !call insert_site(del, vec3f(5.61248064835583,51.2934302817756,0.458397139901375,"__"))  ! 5
    !
    !call insert_site(del, vec3f(-30.7604368378582,-58.6998469054173,1.61649960826247,"__"))
    !call insert_site(del, vec3f(-56.1496004614588,25.3230530640386,2.63530078075907, "__"))
    !call insert_site(del, vec3f(27.4302501937239,-54.7048173425298,-1.53988146175693,"__"))
    !call insert_site(del, vec3f(-24.9515137918661,59.9528869967884,3.60078124630293,"__"))
    !call insert_site(del, vec3f(-53.5794500700446,38.8879494021997,2.07492009407233,"__"))  ! 10
    !
    !call insert_site(del, vec3f(46.7265555224180,-53.0622703071185,-0.588866101697988,"__"))
    ! 
    !call insert_site(del, vec3f(32.1277398759767,16.8599324920652,-1.66338830163992,"__"))
    !
    ! 
    !call insert_site(del, vec3f(50.9507589558205,3.51003633735368,0.141481578104608,"__"))
    !call insert_site(del, vec3f(5.92569695957189,-7.55782528892875,1.97149628022757,"__"))
    !call insert_site(del, vec3f(40.9010723496744,-59.6307644101861,3.49844425644001,"__"))  ! 15
    !
    !call insert_site(del, vec3f(-41.5401400210857,-24.3603232738691,-1.47936312577435,"__"))
    !call insert_site(del, vec3f(-26.3083986268444,45.3135573638845,-0.164543836076551,"__"))
    !call insert_site(del, vec3f(23.4833919331841,-0.427601773453134,-7.062459335288818D-003,"__"))
    !call insert_site(del, vec3f(1.67326780139818,-59.3887544367668,2.20113192244139,"__")) ! 19
    
    call RANDOM_SEED()
    do i = 1, 5000
        call RANDOM_NUMBER(x)
        call RANDOM_NUMBER(y)
        x = x*120-60
        y = y*120-60
        if (MOD(i,1) == 0) then
            print *, i
        end if
        
        point_list(i) = vec3f(x,y,(SIN(x*0.2) + COS(y*0.2))*2,"__")
        
        call insert_site(del, point_list(i))
    end do
    
    call finalize(del)
    
    
    print *, del%vc, del%ec, euler_faces(del%vc, del%ec)
    
    m = get_mesh(del)
 
    call openBin("delaunay_faces.bin",10)
    write(10) m%faces
    close(10)
    
    call openBin("delaunay_vertex.bin",11)
    write(11) m%vertices
    close(11)
end program