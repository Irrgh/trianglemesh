module mesh_helper
    use delaunay
    use quad_edge
    implicit none
        
    type edge_acumulator
        integer(c_intptr_t), allocatable :: edges(:)
        integer :: index
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
    !   Mesh operators   !
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
        type(mesh) :: del
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
    
    function adjacency (del) result (list)
        type(mesh) :: del
        integer(c_intptr_t), allocatable :: list(:)
        integer, allocatable :: visit_count
        integer :: i
        integer(c_intptr_t) :: tmp_edge
        type(edge_struct), pointer :: edge_ptr
        integer, pointer :: visit_ptr
        
        
        allocate(visit_count(del%ec))
        list = list_edges(del)
        visit_count = 0
        
        do i = 1, del%ec
            edge_ptr => deref(list(i))
            visit_ptr => visit_count(i)
            edge_ptr%tmp = c_loc(visit_count(i))
        end do        
    end function
        
        
        
    
    
    
    
    
    function extract_mesh (del) result(faces)
        type(mesh) :: del
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