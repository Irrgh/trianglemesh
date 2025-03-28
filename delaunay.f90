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
    
    type mesh
        type(vec3f), allocatable :: vert(:) ! position of vertices in 3d space
        type(vec3i), allocatable :: tris(:) ! indices of the 3 edges constructing a triangle
        integer(c_intptr_t) :: root         ! edge to start navigating at
        integer(4) :: vc, tc
    end type
    
    contains
    
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
        l = (a%x == b%x) .AND. (a%y == b%y)
    end function
    
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
        call end_points(e,DDATA(a), ODATA(b))
        call splice(e,LNEXT(a))
        call splice(SYM(e), b)
    end function
    
    subroutine swap (e)
        integer(c_intptr_t), intent(in) :: e
        integer(c_intptr_t) :: a,b
        a = OPREV(e)
        b = SYM(OPREV(e))                   ! SYM(OPREV(e)) == RNEXT(e)
        call splice (e, a)
        call splice (SYM(e), b)             
        call splice (e, LNEXT(a))           ! LNEXT(e) == DNEXT(e)
        call splice (SYM(e), LNEXT(b))      ! LNEXT(b) == ONEXT(e)
        call end_points(e, DDATA(a), DDATA(b))
    end subroutine
    
    !-------------------------------!
    !               c               !
    !       tl ----------- tr       !
    !       |       b   -   |       !
    !     a |       -       | d     !
    !       |   -           |       !
    !       bl ----------- br       !
    !               e               !
    !-------------------------------!
    
    subroutine init (m,bl,tl,tr,br)
        type(vec3f), target, intent(in) :: bl,tl,tr,br
        type(mesh), intent(inout) :: m
        type(vec3f), pointer :: t1,t2
        type(edge_struct), pointer :: p1, p2, p3
        integer(c_intptr_t) :: a, b, c, d, e
        
        a = make_edge()
        b = make_edge()
        !c = make_edge()
        
        p1 => deref(a)
        p2 => deref(b)
        p3 => deref(c)
        print *, "a,b,c:", a,b,c
        !print *, "correct_init([a,b,c]):", correct_init(a), correct_init(b), correct_init(c)
        !print *, "I wonder if this works aswell:", verify(a), verify(b), verify(c)
        print *, ccw(tl,bl,tr)
        
        
        t1 => tl
        t2 => bl
        call end_points(a,c_loc(t1),c_loc(t2))    
        t1 => bl
        t2 => tr
        call end_points(b,c_loc(t1),c_loc(t2))
        
        
        
        call splice(SYM(a),b)
        c = connect(b,a)
            
        t1 => tr
        t2 => tl
        call end_points(c,c_loc(t1), c_loc(t2))
        
        
        
        !d = make_edge()
        !call splice(c, d)
        !t1 => tr
        !t2 => br
        !call end_points(d, c_loc(t1), c_loc(t2))
        !
        !e = make_edge()
        !call splice(SYM(d),e)
        !t1 => br
        !t2 => bl
        !call end_points(e, c_loc(t1),c_loc(t2))
        
        call debug(a)
        
        
        
        m%root = a
    end subroutine
        
    subroutine v_proc (edge,closure)
        integer(c_intptr_t), intent(in) :: edge
        type(c_ptr), intent(in) :: closure
        print *, "Traversed edge at: ", edge, ", org: ", org(edge), ", dest: ", dest(edge)
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
        call debug(e)
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
    
    function locate(del, p) result (e)
        type(mesh) :: del
        type(vec3f) :: p
        integer(c_intptr_t) :: e
        e = del%root
        
        do while (.TRUE.)
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
    
    
    subroutine insert_site (del,p)
        type(mesh), intent(inout) :: del
        type(vec3f), intent(in), target :: p
        integer(c_intptr_t) :: e,b,s,t
        type(vec3f), pointer :: tmp
        
        e = locate(del,p)
        
        if (equals_vec3f(p,org(e)) .OR. equals_vec3f(p,dest(e))) then
            return
        else if (on_edge(p,e)) then
            e = OPREV(e)
            call destroy_edge(e)
        end if
        
        b = make_edge()
        tmp => p
        
        call end_points(b, ODATA(e), c_loc(tmp))
        call splice(b,e)
        s = b
             
        do while (LNEXT(e) /= s)
            b = connect(e,SYM(b))
            e = OPREV(b)
            
            call debug(e)
        end do
    
        
        print *, "-------------------------------------"
        print *, org(e), dest(e), org(ONEXT(e)), dest(ONEXT(e)), org(ONEXT(e)), dest(ONEXT(ONEXT(e)))
        print *, "- - - - - - - - - - - - - - - - - - -"
        print *,  org(e), dest(e), org(OPREV(e)), dest(OPREV(e)), org(OPREV(e)), dest(OPREV(OPREV(e)))
        e = LNEXT(e)
        print *, "-------------------------------------"
        print *, org(e), dest(e), org(ONEXT(e)), dest(ONEXT(e)), org(ONEXT(e)), dest(ONEXT(ONEXT(e)))
        print *, "- - - - - - - - - - - - - - - - - - -"
        print *,  org(e), dest(e), org(OPREV(e)), dest(OPREV(e)), org(OPREV(e)), dest(OPREV(OPREV(e)))
        e = LNEXT(e)
        print *, "-------------------------------------"
        print *, org(e), dest(e), org(ONEXT(e)), dest(ONEXT(e)), org(ONEXT(e)), dest(ONEXT(ONEXT(e)))
        print *, "- - - - - - - - - - - - - - - - - - -"
        print *,  org(e), dest(e), org(OPREV(e)), dest(OPREV(e)), org(OPREV(e)), dest(OPREV(OPREV(e)))
        e = LNEXT(e)
        
        
        do while (.TRUE.)
            t = OPREV(e)
            if (right_of(dest(t),e) .AND. in_circle(org(e),dest(t), dest(e), p)) then
                call swap(e)
                e = OPREV(e)
            else if (ONEXT(e) == s) then
                return
            else
                print *, "------------------------------"
                print *, "e:", org(e), dest(e)
                print *, "ONEXT(e):", org(ONEXT(e)), dest(ONEXT(e))
                e = SYM(OPREV(ONEXT(e)))
                
                print *, "LPREV(ONEXT(e):", org(e), dest(e)
            end if
        end do
    end subroutine
    
    
end module