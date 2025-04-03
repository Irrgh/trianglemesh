module mesh_helper
    use delaunay
    use quad_edge
    implicit none
        
    type edge_acumulator
        integer(c_intptr_t), allocatable :: edges(:)
        integer :: index
    end type
    
    type vec3f_
        real(8) :: x,y,z
    end type
    
    type vec3i_
        integer(4) :: x,y,z
    end type
    
    type vec2i_
        integer(4) :: x,y
    end type
    
    type mesh 
        type(vec3f_), allocatable :: vertices(:)
        type(vec2i_), allocatable :: edges(:)
        type(vec3i_), allocatable :: faces(:)
    end type
    
    private edge_acumulator, edge_reduce
    
    contains
    
    function euler_vertices (edges, faces) result (vertices)
        integer :: vertices, faces, edges
        ! 1 = v - e + f
        ! 1 - f = v - e
        ! 1 + e - f = v
        vertices = edges - faces + 1
    end function
    
    function euler_edges (vertices, faces) result (edges)
        integer :: vertices, faces, edges
        ! 1 = v - e + f
        ! 1 + e = v + f
        ! e = v + f - 1
        edges = vertices + faces - 1
    end function
    
    function euler_faces (vertices,edges) result (faces)
        integer :: vertices, faces, edges
        ! 1 = v - e + f
        ! 1 - v = -e + f
        ! 1 + e - v = f
        faces = edges - vertices + 1
    end function
    
    ! euler characteristic for planar graphs 
    function euler_char (vertices, edges, faces) result (char)
        integer :: vertices, faces, edges
        logical :: char
        char = 1 == vertices - edges + faces
    end function
    
    !--------------------!
    !   tm_del operators !
    !--------------------!
    
    subroutine edge_reduce(edge,closure) 
        integer(c_intptr_t), intent(in) :: edge
        type(c_ptr), intent(in) :: closure
        type(edge_acumulator), pointer :: acc
        call c_f_pointer(closure,acc)
        
        acc%edges(acc%index) = edge
        acc%index = acc%index + 1
    end subroutine
    
    function list_edges (del) result(edges)
        type(tm_del) :: del
        integer(c_intptr_t), allocatable :: edges(:)
        type(edge_acumulator), target :: acc
        type(c_ptr) :: ptr
        allocate(acc%edges(del%ec))
        acc%edges = 0
        acc%index = 1
        ptr = c_loc(acc)
        
        call quad_enum(del%root, edge_reduce, ptr) 
        edges = acc%edges
    end function
    
    function get_mesh (del) result (m)
        type(tm_del) :: del
        integer(c_intptr_t), allocatable :: list(:)
        type(mesh) :: m
        integer(c_intptr_t) :: e0, e1, e2
        type(edge_struct), pointer :: p0, p1, p2
        type(vec3f) :: tmp
        integer(c_intptr_t) :: t0, t1, t2, h0, h1, h2, i, f_idx, v_loc
        integer, parameter :: stride = 32    ! size of type would be 26 but alignment is 32
        
        v_loc = loc(del%vertices)
        
        
        list = list_edges(del)
        allocate(m%vertices(del%vc), m%edges(del%ec), m%faces(euler_faces(del%vc,del%ec)))
        
        

        do i = 1, del%vc            ! maps the internal delaunay vertices
            tmp = del%vertices(i)
            m%vertices(i) = vec3f_(tmp%x, tmp%y, tmp%z)    
        end do
        
        do i = 1, del%ec
            p0 => deref(list(i))
            p0%mark = ISHFT(i,2)
        end do
        
        f_idx = 1
        do i = 1, del%ec
            e0 = list(i)
            e1 = LNEXT(e0)
            e2 = LNEXT(e1)
            
            if (LNEXT(e2) /= e0) then   ! cw triangle -> ccw triangle
                e0 = SYM(e0)
                e1 = LNEXT(e0)
                e2 = LNEXT(e1)
            end if
            
            p0 => deref(e0)
            p1 => deref(e1)
            p2 => deref(e2)
            
            t0 = (e0 .AND. 2) / 2
            t1 = (e1 .AND. 2) / 2
            t2 = (e2 .AND. 2) / 2
            
            h0 = p0%mark .AND. ISHFT(1,t0)
            h1 = p1%mark .AND. ISHFT(1,t1)
            h2 = p2%mark .AND. ISHFT(1,t2)
            
            
            if (h0 == 0 .AND. h1 == 0 .AND. h2 == 0) then       ! has this orientation of the edges been used before
                p0%mark = p0%mark .OR. ISHFT(1,t0)
                p1%mark = p1%mark .OR. ISHFT(1,t1)
                p2%mark = p2%mark .OR. ISHFT(1,t2)
                
                t0 = TRANSFER(ODATA(e0),t0)
                t1 = TRANSFER(ODATA(e1),t1)
                t2 = TRANSFER(ODATA(e2),t2)
                
                t0 = (t0 - v_loc) / stride
                t1 = (t1 - v_loc) / stride
                t2 = (t2 - v_loc) / stride    
                
                m%faces(f_idx) = vec3i_(t0,t1,t2)
                f_idx = f_idx + 1
            end if
            
            h0 = TRANSFER(DDATA(e0),h0)
            h0 = (h0 - v_loc) / stride
            m%edges(i)   = vec2i_(t0,h0)
        end do
    end function
        
        
        
    
    
    
    
    
    function extract_tm_del (del) result(faces)
        type(tm_del) :: del
        type(vec3i), allocatable :: faces(:)
        integer, allocatable :: visited(:)
        integer :: face_count, i
        integer(c_intptr_t) :: a,b,c
        face_count = euler_faces(del%vc, del%ec)
        allocate(faces(face_count), visited(del%ec)) 
        
    end function
    
    subroutine adjacency_list (e,c)
        integer(c_intptr_t), intent(in) :: e
        type(c_ptr), intent(in) :: c
        integer(c_intptr_t) :: s
        type(edge_struct), pointer :: p
        character(len =2) :: start_org
        print *, org(e)%n," => ", dest(e)%n, ":"
        start_org = org(e)%n
        
        
        s = ONEXT(e)
        do while (e /= s)
            if (org(s)%n /= start_org) then
                print *, "--> ", org(s)%n, " => ", dest(s)%n, "THIS SHOULD NOT BE POSSIBLE" 
            else 
                print *, "--> ", org(s)%n, " => ", dest(s)%n 
            end if
            
            s = ONEXT(s)
        end do
    end subroutine 
    
end module