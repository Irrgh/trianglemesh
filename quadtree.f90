module quadtree
    use delaunay
    use vector
    
    integer(4), parameter :: root_index = 1
    integer(4), parameter :: max_depth = 64
    
    type triangle
        integer(4) :: p0,p1,p2
        type(vec2_f64) :: c
    end type
    
    
    type quad_bvh
        type(node), allocatable :: nodes(:) 
        type(triangle), allocatable :: tris(:)
        type(vec3_f64), allocatable :: vertices(:)
    end type
    
    type node
        type(vec2_i32) :: min, max  ! 16 bit
        integer(4) :: first_pc       ! first prim or child
        integer(4) :: prim_count
    end type
    
    contains
    
    function is_leaf(n) result (leaf) 
        type(node) :: n
        logical :: leaf
        leaf = n%prim_count > 0
    end function
    
    function centroid(a,b,c) result (d)
        type(vec3_f64) :: a,b,c
        type(vec2_f64) :: d
        a = vec3_f64_add(a,b)
        a = vec3_f64_add(a,c)
        a = vec3_f64_scale(a, 0.3333_8)
        d = vec2_f64(a%x,b%x)
    end function
    
    subroutine quad_bvh_create(bvh, m)
        type(quad_bvh), intent(inout) :: bvh
        type(mesh), intent(in) :: m
        type(node), root
        
        allocate(bvh%vertices(SIZE(m%vertices)))
        allocate(bvh%nodes(SIZE(m%faces)))
        allocate(bvh%tris(SIZE(m%faces)))
        bvh%vertices = m%vertices
        
        call init_tris(bvh,m%faces)
        
        root%first_pc = 1
        root%prim_count = SIZE(m%faces)
        
        bvh%nodes(root_index) = root
        
        call update_bounds(bvh,root_index)
        call subdivide(bvh,root_index)
        
    end subroutine
    
    subroutine init_tris (bvh,faces) 
        type(quad_bvh), intent(inout) :: bvh
        type(vec3_i32), intent(in), allocatable :: faces(:)
        type(triangle) :: t
        type(vec3_i32) :: f
        integer(4) :: i
        
        !OMP$ SIMD                
        do i = i, SIZE(faces)
            f = faces(i)
            t%p0 = f%x-1
            t%p1 = f%y-1
            t%p2 = f%z-1
            t%c = centroid(bvh%vertices(f%x-1),bvh%vertices(f%y-1),bvh%vertices(f%z-1))
            bvh%tris(i) = t
        end do
        
    end subroutine
    
    
    
    
end module