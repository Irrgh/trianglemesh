module mesh_helper
    use delaunay
    use quad_edge
    implicit none
        
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
    
    function extract_mesh (del) result(faces)
        type(mesh) :: del
        type(vec3i), allocatable :: faces(:)
        allocate(faces(100))
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