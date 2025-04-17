module bvh_acceleration
    use vector
    use delaunay
    implicit none
    
    integer(4), parameter :: root_index = 1
    integer(4), parameter :: max_depth = 32
    
    type triangle
        type(vec3_i32) :: t     ! triangle indices
        type(vec3_f64) :: c     ! triangle aabb center
    end type
    
    type ray
        type(vec3_f64) :: org, dir
    end type
    
    type tm_bvh
        type(node), allocatable :: nodes(:) 
        type(triangle), allocatable :: tris(:)
        type(vec3_f64), allocatable :: vertices(:)
        integer(4) :: nodes_used
    end type
    
    type node
        type(vec3_f64) :: min, max  ! 16 bit
        integer(4) :: first_pc  ! first prim or child
        integer(4) :: prim_count
    end type
    
    contains
    
    pure function is_leaf(n) result (leaf) 
        type(node), intent(in) :: n
        logical :: leaf
        leaf = n%prim_count > 0
    end function
    
    function centroid(a,b,c) result (d)
        type(vec3_f64) :: a,b,c,d
        d = vec3_f64_add(a,b)
        d = vec3_f64_add(d,c)
        d = vec3_f64_scale(d, 0.3333_8)
    end function
    
    
    function aabb_center(a,b,c) result (d)
        type(vec3_f64) :: a,b,c
        type(vec3_f64) :: d, min, max
        
        max%arr = (/-1e30,-1e30,-1e30/)
        min%arr = (/1e30, 1e30, 1e30/)

        min = vec3_f64_min(min,a)
        min = vec3_f64_min(min,b)
        min = vec3_f64_min(min,c)
        
        max = vec3_f64_max(max,c)
        max = vec3_f64_max(max,a)
        max = vec3_f64_max(max,b)
        
        d%arr = (max%arr - min%arr) / 2 + min%arr
    end function
    
    subroutine quad_bvh_create(bvh, m)
        type(tm_bvh), intent(inout), target :: bvh
        type(mesh), intent(in) :: m
        type(node), pointer :: root
        
        print *, SIZE(m%vertices), SIZE(m%edges), SIZE(m%faces)
        
        allocate(bvh%vertices(SIZE(m%vertices)))
        allocate(bvh%nodes(SIZE(m%faces)))
        allocate(bvh%tris(SIZE(m%faces)))
        bvh%vertices = m%vertices
        
        call init_tris(bvh,m%faces)
        
        root => bvh%nodes(root_index)
        root%first_pc = 1
        root%prim_count = SIZE(m%faces)
        bvh%nodes_used = 1
        
        
        call update_bounds(bvh,root_index)
        
        
        call subdivide(bvh,root_index)
        
    end subroutine
    
    subroutine init_tris (bvh,faces) 
        type(tm_bvh), intent(inout) :: bvh
        type(vec3_i32), intent(in), allocatable :: faces(:)
        type(triangle) :: t
        type(vec3_i32) :: f
        integer(4) :: i
        
        do i = 1, SIZE(faces)
            f = faces(i)
            t%t%arr = f%arr+1
            t%c = aabb_center(bvh%vertices(f%arr(1)+1),bvh%vertices(f%arr(2)+1),bvh%vertices(f%arr(3)+1))
            bvh%tris(i) = t
        end do
        
    end subroutine
    
    subroutine update_bounds (bvh,index)
        type(tm_bvh), intent(inout), target :: bvh
        integer(4), intent(in) :: index
        real(8), parameter :: min = -1e30
        real(8), parameter :: max = 1e30
        type(node), pointer :: n
        type(vec3_f64) :: tri(3)
        integer(4) :: i
        
        n => bvh%nodes(index)
        n%min%arr = (/max,max,max/)
        n%max%arr = (/min,min,min/)
        
        do i=1, n%prim_count
            
            tri = bvh%vertices(bvh%tris(n%first_pc+i-1)%t%arr)
            
            n%min = vec3_f64_min(n%min,tri(1))
            n%min = vec3_f64_min(n%min,tri(2))
            n%min = vec3_f64_min(n%min,tri(3))
            
            n%max = vec3_f64_max(n%max,tri(1))
            n%max = vec3_f64_max(n%max,tri(2))
            n%max = vec3_f64_max(n%max,tri(3))
            
        end do  
        n => NULL()
    end subroutine
    
    recursive subroutine subdivide(bvh,index)
        type(tm_bvh), intent(inout), target :: bvh
        integer(4), intent(in) :: index
        type(node), pointer :: n
        type(vec3_f64) :: dim
        type(triangle) :: tmp
        real(8) :: split_pos
        integer(4) :: axis, i, j, left_count, left_index, right_index
        
        n => bvh%nodes(index)
        if (n%prim_count <= 4) then
            return
        end if
        
        
        dim = vec3_f64_sub(n%max,n%min)
        axis = 1
        if (dim%arr(2) > dim%arr(1)) axis = 2
        if (dim%arr(2) > dim%arr(axis)) axis = 3
        split_pos = n%min%arr(axis) + dim%arr(axis) * 0.5
        
        i = n%first_pc
        j = i + n%prim_count - 1
            
        do while (i <= j)
            if (bvh%tris(i)%c%arr(axis) < split_pos) then
                i = i + 1
            else
                tmp = bvh%tris(i)
                bvh%tris(i) = bvh%tris(j)
                bvh%tris(j) = tmp
                j = j - 1
            end if
        end do
        
        left_count = i - n%first_pc
        if (left_count == 0 .OR. left_count == n%prim_count) then
            return
        end if
    
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
        
        !if (vec3_f64_equals(bvh%nodes(left_index)%min,n%min) .and. vec3_f64_equals(bvh%nodes(left_index)%max,n%max)) then
        !    print *, "Left child has same aabb as parent", left_index, left_count
        !end if
        !
        !if (vec3_f64_equals(bvh%nodes(right_index)%min,n%min) .and. vec3_f64_equals(bvh%nodes(right_index)%max,n%max)) then
        !    print *, "Right child has same aabb as parent", right_index, bvh%nodes(right_index)%prim_count
        !end if
        
        
        call subdivide(bvh,left_index)
        call subdivide(bvh, right_index)
        
    end subroutine
    
    pure elemental function valid_aabb_intersection(tmin) result (l)
        type(vec2_f64), intent(in) :: tmin
        logical :: l
        l = tmin%arr(2) >= 0 .and. tmin%arr(2) >= tmin%arr(1)
    end function
    
    pure elemental function intersection_aabb(aabb_min,aabb_max,r) result (tmin)
        type(vec3_f64), intent(in) :: aabb_min,aabb_max
        type(ray), intent(in) :: r
        type(vec2_f64) :: tmin      ! tmin%arr(1) == entry      tmin%arr(2) == exit
        type(vec3_f64) :: inv
        real(8) :: tx1, tx2, ty1, ty2, tz1, tz2
        inv%arr = 1.0 / r%dir%arr
        
        tx1 = (aabb_min%arr(1) - r%org%arr(1)) * inv%arr(1)
        tx2 = (aabb_max%arr(1) - r%org%arr(1)) * inv%arr(1)
        
        tmin%arr(1) = MIN(tx1,tx2)
        tmin%arr(2) = MAX(tx1,tx2)
        
        ty1 = (aabb_min%arr(2) - r%org%arr(2)) * inv%arr(2)
        ty2 = (aabb_max%arr(2) - r%org%arr(2)) * inv%arr(2)
        
        tmin%arr(1) = MAX(tmin%arr(1),MIN(ty1,ty2))
        tmin%arr(2) = MIN(tmin%arr(2),MAX(ty1,ty2))
        
        tz1 = (aabb_min%arr(3) - r%org%arr(3)) * inv%arr(3)
        tz2 = (aabb_max%arr(3) - r%org%arr(3)) * inv%arr(3)
        
        tmin%arr(1) = MAX(tmin%arr(1),MIN(tz1,tz2))
        tmin%arr(2) = MIN(tmin%arr(2),MAX(tz1,tz2))
    end function
    
    pure elemental function intersection_triangle(bvh,index,r) result(tmin)
        type(tm_bvh), intent(in) :: bvh
        integer(4), intent(in) :: index
        type(ray), intent(in) :: r
        type(vec3_f64) :: v0, v1, v2, ab, ac, p, t, q
        real(8) :: tmin, det, inv_det, u, v, dist
        
        tmin = -1.0
        
        v0 = bvh%vertices(bvh%tris(index)%t%arr(1))
        v1 = bvh%vertices(bvh%tris(index)%t%arr(2))
        v2 = bvh%vertices(bvh%tris(index)%t%arr(3))
        
        ab = vec3_f64_sub(v1,v0)
        ac = vec3_f64_sub(v2,v0)
        p = vec3_f64_cross(r%dir,ac)
        det = vec3_f64_dot(ab,p)
        
        if (abs(det) < EPS) then
            return
        end if
    
        inv_det = 1 / det
        t = vec3_f64_sub(r%org,v0)
        u = vec3_f64_dot(t,p) * inv_det
        
        if (u < 0.0 .OR. u > 1.0) then
            return
        end if
        
        q = vec3_f64_cross(t,ab)
        v = vec3_f64_dot(r%dir,q) * inv_det
        
        if (v < 0.0 .OR. u + v > 1.0) then
            return
        end if
        
        dist = vec3_f64_dot(ac,q) * inv_det
        if (dist < 0.0) then
            return
        end if
        
        tmin = dist        
    end function
    
    
    pure elemental function intersection_bvh(bvh,r) result(tmin)
        type(tm_bvh), intent(in) :: bvh
        type(ray), intent(in) :: r
        integer(4) :: intersection_stack(max_depth), s_idx, i
        real(8) :: tmin, tcurr, distance_stack(max_depth) 
        type(node) :: root, n, child1, child2
        type(vec2_f64) :: root_intersection, t1, t2
        
        tmin = 1e30
        s_idx = 1
        root = bvh%nodes(root_index)
        root_intersection = intersection_aabb(root%min,root%max,r)
        
        if (valid_aabb_intersection(root_intersection)) then
            intersection_stack(s_idx) = root_index
            distance_stack(s_idx) = root_intersection%arr(1)
        end if
        
        do while (s_idx > 0)
            
            if (tmin < distance_stack(s_idx)) then
                return
            end if
            
            
            n = bvh%nodes(intersection_stack(s_idx))
            
            if (is_leaf(n)) then
               
                do i = 1, n%prim_count
                    tcurr = intersection_triangle(bvh, n%first_pc+i-1, r)
                    if (tcurr >= 0 .and. tcurr < tmin) tmin = tcurr
                end do
            else
                
                child1 = bvh%nodes(n%first_pc)
                child2 = bvh%nodes(n%first_pc+1)
 
                t1 = intersection_aabb(child1%min, child1%max, r)
                t2 = intersection_aabb(child2%min, child2%max, r)
                
                if (valid_aabb_intersection(t1) .and. valid_aabb_intersection(t2)) then
                    if (t1%arr(1) >= 0 .and. t2%arr(1) >= 0) then       ! r%org outside both children
                        if (t1%arr(1) > t2%arr(1)) then                 ! t2 is closer -> on top of stack
                            intersection_stack(s_idx) = n%first_pc
                            intersection_stack(s_idx+1) = n%first_pc + 1
                            distance_stack(s_idx) = t1%arr(1)
                            distance_stack(s_idx+1) = t2%arr(1)
                        else                                            ! t1 is closer -> on top of stack
                            intersection_stack(s_idx) = n%first_pc + 1
                            intersection_stack(s_idx+1) = n%first_pc
                            distance_stack(s_idx) = t2%arr(1)
                            distance_stack(s_idx+1) = t1%arr(1)
                        end if
                    else if (t1%arr(1) < 0 .and. t2%arr(1) < 0) then    ! r%org inside both children
                        if (t1%arr(2) > t2%arr(2)) then                 ! t2 is closer -> on top of stack
                            intersection_stack(s_idx) = n%first_pc
                            intersection_stack(s_idx+1) = n%first_pc + 1
                            distance_stack(s_idx) = t1%arr(1)
                            distance_stack(s_idx+1) = t2%arr(1)
                        else                                            ! t1 is closer -> on top of stack
                            intersection_stack(s_idx) = n%first_pc + 1
                            intersection_stack(s_idx+1) = n%first_pc
                            distance_stack(s_idx) = t2%arr(1)
                            distance_stack(s_idx+1) = t1%arr(1)
                        end if
                    else if (t1%arr(1) < 0) then                        ! r%org inside child1
                        intersection_stack(s_idx) = n%first_pc + 1
                        intersection_stack(s_idx+1) = n%first_pc
                        distance_stack(s_idx) = t2%arr(1)
                        distance_stack(s_idx+1) = t1%arr(1)
                    else if (t2%arr(1) < 0) then                        ! r%org inside child2
                        intersection_stack(s_idx) = n%first_pc
                        intersection_stack(s_idx+1) = n%first_pc + 1
                        distance_stack(s_idx) = t1%arr(1)
                        distance_stack(s_idx+1) = t2%arr(1)
                    end if
                    s_idx = s_idx + 2
                    
                else
                    if (valid_aabb_intersection(t1)) then
                        intersection_stack(s_idx) = n%first_pc
                        distance_stack(s_idx) = t1%arr(1)
                        s_idx = s_idx + 1
                    end if
                    if (valid_aabb_intersection(t2)) then
                        intersection_stack(s_idx) = n%first_pc + 1
                        distance_stack(s_idx) = t2%arr(1)
                        s_idx = s_idx + 1
                    end if
                end if
            end if
            s_idx = s_idx - 1
        end do 
    end function
    
    
    
    
    
end module