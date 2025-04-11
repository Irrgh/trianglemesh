module quadtree
    use delaunay
    use vector
    implicit none
    
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
        integer(4) :: nodes_used
    end type
    
    type node
        type(vec3_f64) :: min, max  ! 16 bit
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
        d%arr = (/a%arr(1),b%arr(1)/)
    end function
    
    subroutine quad_bvh_create(bvh, m)
        type(quad_bvh), intent(inout) :: bvh
        type(mesh), intent(in) :: m
        type(node) :: root
        
        allocate(bvh%vertices(SIZE(m%vertices)))
        allocate(bvh%nodes(SIZE(m%faces)))
        allocate(bvh%tris(SIZE(m%faces)))
        bvh%vertices = m%vertices
        
        call init_tris(bvh,m%faces)
        
        root%first_pc = 1
        root%prim_count = SIZE(m%faces)
        
        bvh%nodes_used = bvh%nodes_used + 1
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
        
        !!OMP$ SIMD                
        do i = 1, SIZE(faces)
            f = faces(i)
            t%p0 = f%arr(1)+1
            t%p1 = f%arr(2)+1
            t%p2 = f%arr(3)+1
            t%c = centroid(bvh%vertices(f%arr(1)+1),bvh%vertices(f%arr(2)+1),bvh%vertices(f%arr(3)+1))
            bvh%tris(i) = t
        end do
        
    end subroutine
    
    subroutine update_bounds (bvh,index)
        type(quad_bvh), intent(inout) :: bvh
        integer(4), intent(in) :: index
        real(8) :: min, max
        type(node) :: n
        type(triangle) :: tri
        integer(4) :: i
        
        min = TINY(1.0_8)
        max = HUGE(1.0_8)
        n = bvh%nodes(index)
        n%min%arr = (/max,max,max/)
        n%max%arr = (/min,min,min/)
        
        do i=1, n%prim_count
            
            tri = bvh%tris(n%first_pc+i-1)
            
            call vec3_f64_min(n%min,n%min,bvh%vertices(tri%p0))
            call vec3_f64_min(n%min,n%min,bvh%vertices(tri%p1))
            call vec3_f64_min(n%min,n%min,bvh%vertices(tri%p2))
            
            call vec3_f64_max(n%max,n%max,bvh%vertices(tri%p0))
            call vec3_f64_max(n%max,n%max,bvh%vertices(tri%p1))
            call vec3_f64_max(n%max,n%max,bvh%vertices(tri%p2))
            
        end do
    end subroutine
    
    recursive subroutine subdivide(bvh,index)
        type(quad_bvh), intent(inout) :: bvh
        integer(4), intent(in) :: index
        type(node) :: n
        type(vec3_f64) :: dim
        type(triangle) :: tmp
        real(8) :: split_pos
        integer(4) :: axis, i, j, left_count, left_index, right_index
        
        n = bvh%nodes(index)
        if (n%prim_count <= 2) return
        
        
        dim = vec3_f64_sub(n%max,n%min)
        axis = 1
        if (dim%arr(2) > dim%arr(1)) axis = 2
        !if (dim%arr(3) > dim%arr(axis)) axis = 2
        split_pos = n%min%arr(axis) + dim%arr(axis) * 0.5
        
        
        i = n%first_pc
        j = i + n%prim_count - 1
        
        do while (i <= j)
            if (bvh%tris(i)%c%arr(axis) < split_pos) then
                i = i + 1
            else
                tmp = bvh%tris(i)
                bvh%tris(i) =  bvh%tris(j)
                bvh%tris(j) = tmp
                j = j - 1
            end if
        end do
        
        left_count = i - n%first_pc
        if (left_count == 0 .OR. left_count == n%prim_count) return
    
        left_index = bvh%nodes_used + 1
        right_index = bvh%nodes_used + 2
        bvh%nodes_used = right_index
        
        bvh%nodes(left_index)%first_pc = n%first_pc
        bvh%nodes(left_index)%prim_count = left_count
        bvh%nodes(right_index)%first_pc = i
        bvh%nodes(right_index)%prim_count = n%prim_count - left_count
        n%first_pc = left_index
        n%prim_count = 0
        
        call update_bounds(bvh,left_index)
        call update_bounds(bvh,right_index)
        
        call subdivide(bvh,left_index)
        call subdivide(bvh, right_index)
        
    end subroutine
    
    
    
end module