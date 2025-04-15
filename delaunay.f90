module delaunay
    use quad_edge
    use iso_c_binding
    use vector
    use euler
    
    type edge_acumulator
        integer(c_intptr_t), allocatable :: edges(:)
        integer :: index
    end type
    
    type mesh 
        type(vec3_f64), allocatable :: vertices(:)
        type(vec2_i32), allocatable :: edges(:)
        type(vec3_i32), allocatable :: faces(:)
    end type
    
    type tm_del
        type(vec3_f64), allocatable :: vertices(:)
        integer(c_intptr_t) :: root         ! edge to start navigating at
        integer(4) :: vc, ec, capacity
        logical :: finalized
        real(8) :: max_radius
    end type
    
    private
    public tm_del, mesh, init, insert_site, get_mesh, list_edges, finalize
    
    contains
    
    subroutine add_edge (del)
        type(tm_del), intent(inout) :: del
        del%ec = del%ec + 1
    end subroutine
    
    subroutine remove_edge(del)
        type(tm_del), intent(inout) :: del
        del%ec = del%ec - 1
    end subroutine
    
    function add_vertex (del,vertex) result (addr)
        type(tm_del), target :: del
        type(vec3_f64) :: vertex
        type(c_ptr) :: addr
        del%vc = del%vc + 1
        del%vertices(del%vc) = vertex
        addr = c_loc(del%vertices(del%vc))
    end function
    
    subroutine remove_vertex (del)
        type(tm_del), intent(inout) :: del
        del%vc = del%vc - 1    
    end subroutine
    
    function get_vertices(del) result (vertices)
        type(tm_del) :: del
        type(vec3_f64), allocatable :: vertices(:)
        allocate(vertices(del%vc))
        vertices = del%vertices(1:del%vc)
    end function
        
    !---------------------------------------!
    !    Auxilliary quad edge operations    !
    !---------------------------------------!
    
    
    subroutine end_points(e,org,dest)
        integer(c_intptr_t), intent(in) :: e
        type(c_ptr), intent(in) :: dest, org
        type(edge_struct), pointer :: tmp
        integer :: io,id
        io = (3 .and. (e .and. 2))+1
        id = (3 .and. (e .and. 2) + 2)+1
        
        tmp => deref(e)
        tmp%data(io) = org
        tmp%data(id) = dest
    end subroutine
    
    ! Add a new connecting edge between a and b, so
    ! that all three edges have the same left face.
    ! Data pointers are also set.
    function connect (a,b) result (e)
        integer(c_intptr_t) :: a,b,e
        e = make_edge()
        call splice(e,LNEXT(a))
        call splice(SYM(e), b)
        call end_points(e,DDATA(a), ODATA(b))
    end function
    
    subroutine swap (e)
        integer(c_intptr_t), intent(in) :: e
        integer(c_intptr_t) :: a,b
        a = OPREV(e)
        b = OPREV(SYM(e))                   ! SYM(OPREV(e)) == RNEXT(e)
        call splice (e, a)
        call splice (SYM(e), b)             
        call splice (e, LNEXT(a))           ! LNEXT(e) == DNEXT(e)
        call splice (SYM(e), LNEXT(b))      ! LNEXT(b) == ONEXT(e)
        call end_points(e, DDATA(a), DDATA(b))
    end subroutine
    
    !--------------------------!
    !                          !
    !            tm            ! 
    !           /  \           !
    !       a  /    \ c        !
    !         /      \         !
    !        /        \        !
    !       bl ------ br       !
    !             b            !
    !--------------------------!
    
    ! Initialises a delaunay triangulation with a bounding triangle.
    ! It is garanteed that all points p with length_vec3f(p) <= max_radius
    ! are valid point to add to the triangulation. 
    subroutine init (del,max_radius,capacity)
        type(tm_del), intent(inout) :: del
        real(8), intent(in) :: max_radius
        integer, intent(in) :: capacity
        type(c_ptr) :: t1,t2,t3
        integer(c_intptr_t) :: a, b, c
        type(vec3_f64), target :: bl, br, tm
        real(8) :: length, height
        real(8), parameter :: pi = 3.14159265358979323846
          
        length = 2 * max_radius * SIN((60.0 / 180.0) * pi)  ! half the side length
        height = 2 * max_radius                             ! from center to tm        SQRT(length**2 + max_radius**2)    
        
        bl%arr(1) = -length
        bl%arr(2) = -max_radius
        bl%arr(3) = 0
        
        br%arr(1) = length
        br%arr(2) = -max_radius
        br%arr(3) = 0
        
        tm%arr(1) = 0
        tm%arr(2) = height
        tm%arr(3) = 0 
        
        !bl = vec3_f64_create(-length, -max_radius, 0.0_8)
        !br = vec3_f64_create(length, -max_radius, 0.0_8)
        !tm = vec3_f64_create(0.0_8, height, 0.0_8)
        
        del%max_radius = max_radius
        del%capacity = capacity
        
        a = make_edge()
        b = make_edge()
        c = make_edge()
        
        
        allocate(del%vertices(capacity+3))
        del%vc = 0
        del%ec = 3
        
        t1 = add_vertex(del,tm)
        t2 = add_vertex(del,bl)
        t3 = add_vertex(del,br)
        call end_points(a,t1,t2)    
        call end_points(b,t2,t3)
        call end_points(c,t3,t1)
        
        call splice(SYM(a),b)
        call splice(SYM(b),c)
        call splice(SYM(c),a)
        
        del%root = a
        del%finalized = .FALSE.
    end subroutine
            

    
    function org(e) result (vec)
        type(vec3_f64), pointer :: vec
        integer(c_intptr_t) :: e
       
        call c_f_pointer(ODATA(e),vec)
    end function
    
    function dest(e) result (vec)
        type(vec3_f64), pointer :: vec
        integer(c_intptr_t) :: e
        
        call c_f_pointer(DDATA(e),vec)
    end function
    
    subroutine debug (e)
        integer(c_intptr_t), intent(in) :: e
        print *, "----------------------------"
        print *, "e:"
        print *, org(e), dest(e), e
        print *, "SYM(e):"
        print *, org(SYM(e)), dest(SYM(e)), SYM(e)
        print *, "----------------------------"
        print *, "LPREV(e):"
        print *, org(LPREV(e)), dest(LPREV(e)), LPREV(e)
        print *, "ONEXT(e):"
        print *, org(ONEXT(e)), dest(ONEXT(e)), ONEXT(e)
        print *, "DPREV(e):"
        print *, org(DPREV(e)), dest(DPREV(e)), DPREV(e)
        print *, "LNEXT(e):"
        print *, org(LNEXT(e)), dest(LNEXT(e)), LNEXT(e)
        print *, "-----------------------------"
        print *, "OPREV(e):"
        print *, org(OPREV(e)), dest(OPREV(e)), OPREV(e)
        print *, "RNEXT(e):"
        print *, org(RNEXT(e)), dest(RNEXT(e)), RNEXT(e)
        print *, "RPREV(e):"
        print *, org(RPREV(e)), dest(RPREV(e)), RPREV(e)
        print *, "DNEXT(e):"
        print *, org(DNEXT(e)), dest(DNEXT(e)), DNEXT(e)
    end subroutine
    
    
    
    
    ! Twice the area of the triangle, positive if ccw
    function tri_area(a,b,c) result (area)
        type(vec3_f64) :: a,b,c
        real(8) :: area
        area = (b%arr(1) - a%arr(1)) * (c%arr(2) - a%arr(2)) - (b%arr(2) - a%arr(2)) * (c%arr(1) - a%arr(1))
    end function
    
    function in_circle(a,b,c,d) result (l)
        type(vec3_f64) :: a,b,c,d
        logical :: l
        l = (a%arr(1)**2 + a%arr(2)**2) * tri_area(b,c,d) - &
          & (b%arr(1)**2 + b%arr(2)**2) * tri_area(a,c,d) + &
          & (c%arr(1)**2 + c%arr(2)**2) * tri_area(a,b,d) - &
          & (d%arr(1)**2 + d%arr(2)**2) * tri_area(a,b,c) > 0
    end function
    
    function ccw (a,b,c) result(l)
        type(vec3_f64) :: a,b,c
        logical :: l
        l = tri_area(a,b,c) > 0
    end function
    
    function right_of (p,e) result (l)
        type(vec3_f64) :: p
        integer(c_intptr_t) :: e
        logical :: l
        l = ccw(p,dest(e), org(e))
    end function
    
    function left_of (p,e) result (l)
        type(vec3_f64) :: p
        integer(c_intptr_t) :: e
        logical :: l
        l = ccw(p,org(e), dest(e))
    end function
    
    function on_edge(p,e) result (l)
        type(vec3_f64) :: p, o, d
        integer(c_intptr_t) :: e
        logical :: l
        real(8) :: t1,t2,t3,m
        
        o = org(e)
        d = dest(e)
        t1 = vec3_f64_length(vec3_f64_sub(p,o))
        t2 = vec3_f64_length(vec3_f64_sub(p,d))
        
        if (t1 < EPS .OR. t2 < EPS) then
            l = .TRUE.
            return
        end if
        t3 = vec3_f64_length(vec3_f64_sub(o, d))
        
        l = abs(tri_area(p,o,d)) < EPS * 2 * t3
    end function
    
    ! Finds one edge of the triangle of the delaunay triangulation that contains the point.
    function locate(del, p) result (e)
        type(tm_del) :: del
        type(vec3_f64) :: p
        integer(c_intptr_t) :: e,i
        type(mesh) :: m
        e = del%root
        i = 1
   
        do while (.TRUE.)            
            if (i > del%ec * 2) then 
                m = get_mesh(del)
                print *, m%vertices
                print *, m%faces
                print *, m%edges
            end if
            
            if (vec3_f64_equals(p,org(e)) .OR. vec3_f64_equals(p,dest(e))) then
                return
            else if (right_of(p,e)) then
                e = SYM(e)
            else if (.NOT. right_of(p,ONEXT(e))) then
                e = ONEXT(e)
            else if (.NOT. right_of(p,DPREV(e))) then
                e = DPREV(e)
            else
                return
            end if
            i = i + 1
        end do
        
        
    end function
    
    
    ! Inserts a point p into a existing delaunay triangulation.
    ! Valid if length_vec3f(p) <= del%max_radius
    ! If point is already contained in the triangulation do nothing.
    subroutine insert_site (del,p)
        type(tm_del), intent(inout) :: del
        type(vec3_f64), intent(in), target :: p
        type(c_ptr) :: tmp
        integer(c_intptr_t) :: e,b,s,t
        
        if (del%finalized) then
            print *,"Point can not be inserted: Triangulation is finalized."
            return
        end if
        
        if (vec3_f64_length(p) > del%max_radius) then
            print *,"Point can not be inserted: p(",p%arr(1),",",p%arr(2), ") is outside of specified max_radius of ", del%max_radius
            return
        end if
        
        if (del%vc - 3 >= del%capacity) then
            print *, "Point can not be inserted: Maximum capacity of ", del%capacity, " vertices reached."
            return
        end if
        
        e = locate(del,p)       ! one edge of the containing triangle
        
        
        if (vec3_f64_equals(p,org(e)) .OR. vec3_f64_equals(p,dest(e))) then
            return
        else if (on_edge(p,e)) then
            e = OPREV(e)
            call destroy_edge(ONEXT(e))
            call remove_edge(del)
        end if
        
        tmp = add_vertex(del,p)
        b = make_edge()
        call add_edge(del)
    
        s = b
        call end_points(b, ODATA(e), tmp)
        call splice(b,e)
        
        
        do while (LNEXT(e) /= s)
            b = connect(e,SYM(b))
            e = OPREV(b)
            call add_edge(del)
        end do
        
        do while (.TRUE.)
            t = OPREV(e)
            if (right_of(dest(t),e) .AND. in_circle(org(e),dest(t), dest(e), p)) then
                call swap(e)
                e = t
            else if (ONEXT(e) == s) then
                return
            else
                e = LPREV(ONEXT(e))
            end if
        end do
    end subroutine
    
    subroutine finalize (del)
        type(tm_del), intent(inout) :: del
        integer(c_intptr_t) :: a,b,c,t0,t1
        
        if (del%finalized) then
            print *, "Triangulation cannot be finalized: Triangulation is already finalized"
            return
        end if
        
        a = del%root
        b = RPREV(a)
        c = RPREV(b)
        
        
        if (RPREV(c) == a) then
            
            del%root = RPREV(LNEXT(a))
             
            t0 = ONEXT(b)
     
            do while (t0 /= SYM(a))     ! destroys all ccw edges between b and a
                t1 = ONEXT(t0)
                call destroy_edge(t0)
                call remove_edge(del)
                t0 = t1
            end do
                
            t0 = ONEXT(a)
            
            do while (t0 /= SYM(c))     ! destorys all ccw edges between c and b
                t1 = ONEXT(t0)
                call destroy_edge(t0)
                call remove_edge(del)
                t0 = t1
            end do
            
            t0 = ONEXT(c)               
            
            do while (t0 /= SYM(b))     ! destroys all ccw edges between a and c
                t1 = ONEXT(t0)
                call destroy_edge(t0)
                call remove_edge(del)
                t0 = t1
            end do
            
            call destroy_edge(a)
            call destroy_edge(b)
            call destroy_edge(c)
            call remove_edge(del)
            call remove_edge(del)
            call remove_edge(del)
            
            call remove_vertex(del)
            call remove_vertex(del)
            call remove_vertex(del)
            
            del%finalized = .TRUE.
        else 
            print *, "Triangulation cannot be finalized: RPREV() LOOP of del%root does not exist. Ensure del%root is either tm, bl or br?"
            return
        end if
    end subroutine
    
    
    subroutine edge_reduce(edge,closure) 
        integer(c_intptr_t), intent(in) :: edge
        type(c_ptr), intent(in) :: closure
        type(edge_acumulator), pointer :: acc
        call c_f_pointer(closure,acc)
        
        acc%edges(acc%index) = edge
        acc%index = acc%index + 1
    end subroutine
    
    function list_edges_ (del) result(edges)
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
    
    
    
    function list_edges (del) result(edges)
        type(tm_del) :: del
        integer(c_intptr_t), allocatable :: edges(:), stack(:)
        integer(c_intptr_t) :: e
        type(edge_struct), pointer :: edge
        integer :: i,j, mark
        
        allocate(edges(del%ec), stack(del%ec))
        i = 1
        j = 1
        mark = next_mark
        next_mark = next_mark + 1
        edges = 0
        stack = 0
        stack(1) = del%root
        
        do while (j - i + 1 > 0)
            
            e = stack(j-i+1)
            stack(j-i+1) = 0
            edge => deref(e)
            i = i + 1
            do while (edge%mark /= mark)
                edges(j) = e 
                edge%mark = mark
                j = j + 1
                stack(j-i+1) = ONEXT(e)
                e = ONEXT(SYM(e))
                edge => deref(e)
            end do
        end do    
        deallocate(stack)
    end function
    
    function get_mesh (del) result (m)
        type(tm_del) :: del
        integer(c_intptr_t), allocatable :: list(:)
        type(mesh) :: m
        integer(c_intptr_t) :: e0, e1, e2
        type(edge_struct), pointer :: p0, p1, p2
        integer(c_intptr_t) :: t0, t1, t2, h0, h1, h2, i, f_idx, v_loc, v_offset, inv_c
        integer :: stride
        
        stride = SIZEOF(del%vertices(1))
        v_loc = loc(del%vertices)
        
        list = list_edges(del)
        
        allocate(m%edges(del%ec), m%faces(euler_faces(del%vc,del%ec)), m%vertices(del%vc))
        
        v_offset = 0
        if (del%finalized) v_offset = 3

        m%vertices = del%vertices(1+v_offset:del%vc+v_offset)
        
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
                if (LNEXT(e2) == e0) then
                    inv_c = inv_c + 1
                    print *, "Inverted Triangle:", i
                else
                    print *, "wtf", i, f_idx
                end if
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
                
                t0 = ((t0 - v_loc) / stride) - v_offset
                t1 = ((t1 - v_loc) / stride) - v_offset
                t2 = ((t2 - v_loc) / stride) - v_offset
                
                m%faces(f_idx)%arr = (/t0,t1,t2/)
                f_idx = f_idx + 1
            end if
            
            h0 = TRANSFER(DDATA(e0),h0)
            h0 = ((h0 - v_loc) / stride) - v_offset
            m%edges(i)%arr = (/t0,h0/)
        end do      
    end function
    
    
end module