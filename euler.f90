module euler
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
end module