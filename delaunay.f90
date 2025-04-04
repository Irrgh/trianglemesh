module delaunay
    use quad_edge
    use iso_c_binding
    
    type vec3f
        real(8) :: x,y,z
        character(len = 2) :: n
    end type
    
    
    type vec3i
        integer(4) :: x,y,z
    end type
    
    
    type tm_del
        type(vec3f), allocatable :: vertices(:)
        integer(c_intptr_t) :: root         ! edge to start navigating at
        integer(4) :: vc, ec, capacity
        real(8) :: max_radius
    end type
    
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
        type(vec3f) :: vertex
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
        type(vec3f), allocatable :: vertices(:)
        allocate(vertices(del%vc))
        vertices = del%vertices(1:del%vc)
    end function
        
    
    
    
    !-----------------------------------!
    !          Bad vector math          !
    !-----------------------------------!
    
    function add_vec3f(a,b) result (c)
        type(vec3f) :: a,b,c
        c%x = a%x + b%x
        c%y = a%y + b%y
        c%z = a%z + b%z
    end function
    
    function sub_vec3f(a,b) result (c)
        type(vec3f) :: a,b,c
        c%x = a%x - b%x
        c%y = a%y - b%y
        c%z = a%z - b%z
    end function
    
    function length_vec3f(a) result (l)
        type(vec3f) :: a
        real(8) :: l
        l = sqrt(a%x**2 + a%y**2)
    end function
    
    function equals_vec3f(a,b) result (l)
        type(vec3f) :: a,b
        logical :: l
        real(8), parameter :: EPS = 0.001
        l = abs(a%x - b%x) < EPS .AND. abs(a%y - b%y) < EPS
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
        real, intent(in) :: max_radius
        integer, intent(in) :: capacity
        type(c_ptr) :: t1,t2,t3
        integer(c_intptr_t) :: a, b, c
        type(vec3f), target :: bl, br, tm
        real(8) :: length, height
        real(8), parameter :: pi = 3.14159265358979323846
          
        length = 2 * max_radius * SIN((60.0 / 180.0) * pi)  ! half the side length
        height = 2 * max_radius                             ! from center to tm        SQRT(length**2 + max_radius**2)    
        
        bl = vec3f(-length, -max_radius, 0, "bl")
        br = vec3f(length, -max_radius, 0, "br")
        tm = vec3f(0, height, 0, "tm")
        
        del%max_radius = max_radius
        del%capacity = capacity
        
        a = make_edge()
        b = make_edge()
        c = make_edge()
        
        
        allocate(del%vertices(capacity))
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
    end subroutine
            
    subroutine print_edge (edge,closure)
        integer(c_intptr_t), intent(in) :: edge
        type(c_ptr), intent(in) :: closure
        print *, "Traversed edge at: ", edge, ", org: ", org(edge)%n, ", dest: ", dest(edge)%n
    end subroutine
    
    function org(e) result (vec)
        type(vec3f), pointer :: vec
        integer(c_intptr_t) :: e
       
        call c_f_pointer(ODATA(e),vec)
    end function
    
    function dest(e) result (vec)
        type(vec3f), pointer :: vec
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
        type(vec3f) :: a,b,c
        real(8) :: area
        area = (b%x - a%x) * (c%y - a%y) - (b%y - a%y) * (c%x - a%x)
    end function
    
    function in_circle(a,b,c,d) result (l)
        type(vec3f) :: a,b,c,d
        logical :: l
        
        l = (a%x**2 + a%y**2) * tri_area(b,c,d) - &
          & (b%x**2 + b%y**2) * tri_area(a,c,d) + &
          & (c%x**2 + c%y**2) * tri_area(a,b,d) - &
          & (d%x**2 + d%y**2) * tri_area(a,b,c) > 0
        
    end function
    
    function ccw (a,b,c) result(l)
        type(vec3f) :: a,b,c
        logical :: l
        l = tri_area(a,b,c) > 0
    end function
    
    function right_of (p,e) result (l)
        type(vec3f) :: p
        integer(c_intptr_t) :: e
        logical :: l
        l = ccw(p,dest(e), org(e))
    end function
    
    function left_of (p,e) result (l)
        type(vec3f) :: p
        integer(c_intptr_t) :: e
        logical :: l
        l = ccw(p,org(e), dest(e))
    end function
    
    function on_edge(p,e) result (l)
        type(vec3f) :: p, o, d
        integer(c_intptr_t) :: e
        logical :: l
        real(8) :: t1,t2,t3,m
        real(8), parameter :: EPS = 1e-3
        
        o = org(e)
        d = dest(e)
        t1 = length_vec3f(sub_vec3f(p,o))
        t2 = length_vec3f(sub_vec3f(p,d))
        
        if (t1 < EPS .OR. t2 < EPS) then
            l = .TRUE.
            return
        end if
        t3 = length_vec3f(sub_vec3f(o, d))
        
        if (t1 > t3 .or. t2 > t3) then
            l = .TRUE.
            return
        end if
        
        l = abs(tri_area(p,o,d)) < EPS * 2 * t3
    end function
    
    ! Finds one edge of the triangle of the delaunay triangulation that contains the point.
    function locate(del, p) result (e)
        type(tm_del) :: del
        type(vec3f) :: p
        integer(c_intptr_t) :: e, i
        e = del%root
        
        i = 1
        do while (.TRUE.)
            
            if (del%ec * 4 < i) then
                !call print_mesh_info(del)
            end if
    
            i = i + i
            
            if (equals_vec3f(p,org(e)) .OR. equals_vec3f(p,dest(e))) then
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
        end do
        
        
    end function
    
    
    ! Inserts a point p into a existing delaunay triangulation.
    ! Valid if length_vec3f(p) <= del%max_radius
    ! If point is already contained in the triangulation do nothing.
    subroutine insert_site (del,p)
        type(tm_del), intent(inout) :: del
        type(vec3f), intent(in), target :: p
        type(c_ptr) :: tmp
        integer(c_intptr_t) :: e,b,s,t
        
        if (length_vec3f(p) > del%max_radius) then
            print *,"Cannot add point to delaunay triangulation: p(",p%x,",",p%y, ") is outside of specified max_radius of ", del%max_radius
            stop
        end if
        
        if (del%vc >= del%capacity) then
            print *, "Maximum capacity of ", del%capacity, " vertices reached."
            stop
        end if
        
        
        e = locate(del,p)       ! one edge of the containing triangle
        
        
        if (equals_vec3f(p,org(e)) .OR. equals_vec3f(p,dest(e))) then
            return
        else if (on_edge(p,e)) then
            e = OPREV(e)
            call destroy_edge(ONEXT(e))
            call remove_edge(del)
        end if
        
        tmp = add_vertex(del,p)
        b = make_edge()
    
        call end_points(b, ODATA(e), tmp)
        call splice(b,e)
        
        s = b
        call add_edge(del)
        
        do while (LNEXT(e) /= s)
            b = connect(e,SYM(b))
            e = OPREV(b)
            call add_edge(del)
        end do
        
        do while (.TRUE.)
            t = OPREV(e)
            if (right_of(dest(t),e) .AND. in_circle(org(e),dest(t), dest(e), p)) then
                call swap(e)
                e = OPREV(e)
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
        
        a = del%root
        b = RPREV(a)
        c = RPREV(b)
        
        if (RPREV(c) == a) then
            
            del%root = RPREV(LNEXT(a))
             
            !t0 = ONEXT(b)
            !
            !do while (t0 /= SYM(a))
            !    t1 = ONEXT(t0)
            !    call destroy_edge(t0)
            !    call remove_edge(del)
            !    t0 = t1
            !end do
            !
            !t0 = ONEXT(a)
            
            !do while (t0 /= SYM(c))
            !    t1 = ONEXT(t0)
            !    call destroy_edge(t0)
            !    call remove_edge(del)
            !    t0 = t1
            !end do
            !
            !t0 = ONEXT(c)
            !
            !do while (t0 /= SYM(b))
            !    t1 = ONEXT(t0)
            !    call destroy_edge(t0)
            !    call remove_edge(del)
            !    t0 = t1
            !end do
            
            call destroy_edge(a)
            call destroy_edge(b)
            !call destroy_edge(c)
            call remove_edge(del)
            call remove_edge(del)
            !call remove_edge(del)
            !
            !call remove_vertex(del)
            !call remove_vertex(del)
            !call remove_vertex(del)
            
        end if
    end subroutine
    
    
    
end module