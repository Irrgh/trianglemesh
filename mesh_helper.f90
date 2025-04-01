module mesh_helper
    use delaunay
    implicit none
        
    contains
    
    function euler_vertices (edges, faces) result (vertices)
        integer :: vertices, faces, edges
        ! 2 = v - e + f
        ! 2 - f = v - e
        ! 2 + e - f = v
        vertices = edges - faces + 2
    end function
    
    function euler_edges (vertices, faces) result (edges)
        integer :: vertices, faces, edges
        ! 2 = v - e + f
        ! 2 + e = v + f
        ! e = v + f - 2
        edges = vertices + faces - 2
    end function
    
    function euler_faces (vertices,edges) result (faces)
        integer :: vertices, faces, edges
        ! 2 = v - e + f
        ! 2 - v = -e +f
        ! 2 + e - v = f
        faces = edges - vertices + 2
    end function
    
    function euler_char (vertices, edges, faces) result (char)
        integer :: vertices, faces, edges
        logical :: char
        char = 2 == vertices - edges + faces
    end function
    
    function extract_mesh (del) result(faces)
        type(mesh) :: del
        type(vec3i), allocatable :: faces(:)
        allocate(faces(100))
    end function
    
end module