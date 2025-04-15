module quad_edge
    use iso_c_binding
    implicit none
    
    !---------------------------------------!
    !   Interface for quad_enum callback    !
    !---------------------------------------!
    abstract interface      
        subroutine visit_proc(edge,closure)
            use iso_c_binding, only: c_intptr_t, c_ptr
            integer(c_intptr_t), intent(in) :: edge
            type(c_ptr), intent(in) :: closure
        end subroutine
    end interface
    
    type edge_struct
        integer(c_intptr_t) :: next(4)
        type(c_ptr) :: data(4)
        integer(c_int64_t) :: mark
    end type
    
    integer(c_int64_t) :: next_mark = 1
    
    private quad_do_enum
    
    contains
    
    function deref (e) result (edge)
        integer(c_intptr_t), intent(in) :: e
        type(c_ptr) :: p
        type(edge_struct), pointer :: edge
        p = TRANSFER(ISHFT(ISHFT(e,-2),2),p)
        call c_f_pointer(p,edge)
    end function
    
    !-----------------------------------------------!
    !           Edge orientation operators          !
    !-----------------------------------------------!
    
    pure function ROT (e) result (rot_e)
        integer(c_intptr_t), intent(in) :: e
        integer(c_intptr_t) :: rot_e
        rot_e = ISHFT(ISHFT(e,-2),2) + ((e+1) .AND. 3)
    end function
        
    pure function SYM (e) result (sym_e)
        integer(c_intptr_t), intent(in) :: e
        integer(c_intptr_t) :: sym_e
        sym_e = ISHFT(ISHFT(e,-2),2) + ((e+2) .AND. 3)
    end function
    
    pure function TOR (e) result (tor_e)
        integer(c_intptr_t), intent(in) :: e
        integer(c_intptr_t) :: tor_e
        tor_e = ISHFT(ISHFT(e,-2),2) + ((e+3) .AND. 3)
    end function
    
    !-------------------------------------------!
    !       Vertex/face walking operators       !
    !-------------------------------------------!
    
    function ONEXT (e) result (onext_e)
        integer(c_intptr_t), intent(in) :: e
        integer(c_intptr_t) :: onext_e
        type(edge_struct), pointer :: edge
        edge => deref(e)
        onext_e = edge%next((e .AND. 3)+1)
    end function
    
    function ROTRNEXT (e) result (rotrnext_e)
        integer(c_intptr_t), intent(in) :: e
        integer(c_intptr_t) :: rotrnext_e
        type(edge_struct), pointer :: edge
        edge => deref(e)
        rotrnext_e = edge%next(((e+1) .AND. 3)+1)
    end function
    
    function SYMDNEXT (e) result (symdnext_e)
        integer(c_intptr_t), intent(in) :: e
        integer(c_intptr_t) :: symdnext_e
        type(edge_struct), pointer :: edge
        edge => deref(e)
        symdnext_e = edge%next(((e+2) .AND. 3)+1)
    end function
    
    function TORLNEXT (e) result (torlnext_e)
        integer(c_intptr_t), intent(in) :: e
        integer(c_intptr_t) :: torlnext_e
        type(edge_struct), pointer :: edge
        integer :: idx
        idx = ((e+3) .AND. 3)
        edge => deref(e)
        torlnext_e = edge%next(((e+3) .AND. 3)+1)
    end function

    
    function RNEXT (e) result (rnext_e)
        integer(c_intptr_t), intent(in) :: e
        integer(c_intptr_t) :: rnext_e
        rnext_e = TOR(ROTRNEXT(e))
    end function
    
    function DNEXT (e) result (dnext_e)
        integer(c_intptr_t), intent(in) :: e
        integer(c_intptr_t) :: dnext_e
        dnext_e = SYM(SYMDNEXT(e))
    end function
    
    function LNEXT (e) result (lnext_e)
        integer(c_intptr_t), intent(in) :: e
        integer(c_intptr_t) :: lnext_e
        lnext_e = ROT(TORLNEXT(e))
    end function
    
   
    
    function OPREV (e) result (oprev_e)
        integer(c_intptr_t), intent(in) :: e
        integer(c_intptr_t) :: oprev_e
        oprev_e = ROT(ROTRNEXT(e))
    end function
    
    function DPREV (e) result (dprev_e)
        integer(c_intptr_t), intent(in) :: e
        integer(c_intptr_t) :: dprev_e
        dprev_e = TOR(TORLNEXT(e))
    end function
    
    function RPREV (e) result (rprev_e)
        integer(c_intptr_t), intent(in) :: e
        integer(c_intptr_t) :: rprev_e
        rprev_e = SYMDNEXT(e)
    end function
    
    function LPREV (e) result (lprev_e)
        integer(c_intptr_t), intent(in) :: e
        integer(c_intptr_t) :: lprev_e
        lprev_e = SYM(ONEXT(e))
    end function
    
    
    !---------------------------------------!
    !             Data Pointers             !
    !---------------------------------------!
    
    
    function ODATA (e) result (odata_v)
        integer(c_intptr_t), intent(in) :: e
        type(c_ptr) :: odata_v
        type(edge_struct), pointer :: edge
        edge => deref(e)
        odata_v = edge%data((e .AND. 3)+1)
    end function
    
    function RDATA (e) result (rdata_v)
        integer(c_intptr_t), intent(in) :: e
        type(c_ptr) :: rdata_v
        type(edge_struct), pointer :: edge
        edge => deref(e)
        rdata_v = edge%data(((e+1) .AND. 3)+1)
    end function
    
    function DDATA (e) result (ddata_v)
        integer(c_intptr_t), intent(in) :: e
        type(c_ptr) :: ddata_v
        type(edge_struct), pointer :: edge
        edge => deref(e)
        ddata_v = edge%data(((e+2) .AND. 3)+1)
    end function
    
    function LDATA (e) result (ldata_v)
        integer(c_intptr_t), intent(in) :: e
        type(c_ptr) :: ldata_v
        type(edge_struct), pointer :: edge
        edge => deref(e)
        ldata_v = edge%data(((e+3) .AND. 3)+1)
    end function
    
    !-------------------------------!
    !       Make a new edge         !
    !-------------------------------!
    
    pure function make_edge () result (e)
        type(edge_struct), pointer :: edge
        type(c_ptr) :: p
        integer(c_intptr_t) :: e
        allocate(edge)
        p = c_loc(edge)
        e = TRANSFER(p,e)

        edge%next(1) = e
        edge%next(3) = SYM(e)   ! SYMDNEXT
        edge%next(2) = TOR(e)   ! ROTRNEXT
        edge%next(4) = ROT(e)   ! TORLNEXT
        edge%data = c_null_ptr
        edge%mark = 0
    end function
    
    function correct_init (e) result (l)
        logical :: l, c1, c2, c3, c4
        integer(c_intptr_t), intent(in) :: e
        integer(c_intptr_t) :: f
        c1 = (LNEXT(e) == RNEXT(e)) .and. (RNEXT(e) == SYM(e))
        c2 = (ONEXT(e) == OPREV(e)) .and. (OPREV(e) == e)
        f = ROT(e)
        c3 = (LNEXT(f) == RNEXT(f)) .and. (RNEXt(f) == f)
        c4 = (ONEXT(f) == OPREV(f)) .and. (OPREV(f) == SYM(f)) 
        l = (c1 .and. c2) .and. (c3 .and. c4)
    end function
    
    function verify (e) result (l)
        logical :: l, c1, c2, c3, c4
        integer(c_intptr_t), intent(in) :: e
        c1 = ONEXT(e) == e
        c2 = ONEXT(SYM(e)) == SYM(e)
        c3 = ONEXT(ROT(e)) == TOR(e)
        c4 = ONEXT(TOR(e)) == ROT(e)
        l = (c1 .and. c2) .and. (c3 .and. c4)
    end function
    
    
    !---------------------------!
    !       Delete an edge      !
    !---------------------------!
    
    subroutine destroy_edge (e)
        integer(c_intptr_t), intent(in) :: e
        integer(c_intptr_t) :: f
        type(edge_struct), pointer :: edge
        
        f = SYM(e)
    
        if (ONEXT(e) /= e) then
            call splice(e,OPREV(e))
        end if
        if (ONEXT(f) /= f) then
            call splice(f,OPREV(f))
        end if
    
        edge => deref(e)
        edge%next = 0
        edge%data = c_null_ptr
        deallocate(edge)
        edge => NULL()
    end subroutine
    
    !-------------------------------!
    !       Splice primitive        !
    !-------------------------------!
    subroutine splice(a,b)
        integer(c_intptr_t), intent(in) :: a, b    
        integer(c_intptr_t), target :: alpha, beta, t1,t2,t3,t4
        type(edge_struct), pointer :: tmp
        
        alpha = ROT(ONEXT(a))
        beta = ROT(ONEXT(b))
        t1 = ONEXT(a)
        t2 = ONEXT(b)
        t3 = ONEXT(ROT(ONEXT(a)))
        t4 = ONEXT(ROT(ONEXT(b)))
        
        tmp => deref(a)
        tmp%next((a .AND. 3)+1) = t2
        
        tmp => deref(alpha)
        tmp%next((alpha .AND. 3)+1) = t4
        
        tmp => deref(b)
        tmp%next((b .AND. 3)+1) = t1
        
        tmp => deref(beta)
        tmp%next((beta .AND. 3)+1) = t3
        
        tmp => NULL()
    end subroutine
    
    subroutine splice_check (a,b) 
        integer(c_intptr_t), intent(in) :: a,b
        integer(c_intptr_t) :: an(4), bn(4)
        type(edge_struct), pointer :: tmp
        
        tmp => deref(a)
        an = tmp%next
        tmp => deref(b)
        bn = tmp%next
        
        call splice(a,b)
        call splice(a,b)
        
        tmp => deref(a)
        
        print *, an(1) == tmp%next(1), an(2) == tmp%next(2), an(3) == tmp%next(3), an(4) == tmp%next(4)
        tmp%next = an
        tmp => deref(b)
        print *, bn(1) == tmp%next(1), bn(2) == tmp%next(2), bn(3) == tmp%next(3), bn(4) == tmp%next(4)
        tmp%next = bn
    end subroutine       
    
    !-------------------------------------------------------------------------------!
    !   Enumerates undirected primal edges reachable from $a$.                      !
    !                                                                               !
    !   Calls visit_proc(e, closure) for every edge $e$ that can be reached from    !
    !   edge $a$ by a chain of SYM and ONEXT calls; except that exactly one         !
    !   of $e$ and SYM(e) is visited.                                               !
    !-------------------------------------------------------------------------------!
    subroutine quad_enum(a,v_proc,closure)
        integer(c_intptr_t), intent(in) :: a
        procedure(visit_proc) :: v_proc
        type(c_ptr),intent(in) :: closure
        integer(c_int64_t) :: mark
        mark = next_mark
        next_mark = next_mark + 1
        
        if (next_mark == 0) then 
            next_mark = 1
        end if
        
        call quad_do_enum(a,v_proc,closure,mark)
    end subroutine
    
    recursive subroutine quad_do_enum(a,v_proc,closure,mark)
        integer(c_intptr_t), intent(in) :: a
        procedure(visit_proc) :: v_proc
        type(c_ptr),intent(in) :: closure
        type(edge_struct), pointer :: edge
        integer(c_int64_t), intent(in) :: mark
        integer(c_intptr_t) :: t
        t = a
        do while (deref(t)%mark /= mark)
            call v_proc(t,closure)
            edge => deref(t)
            edge%mark = mark
            call quad_do_enum(ONEXT(SYM(t)), v_proc, closure, mark)
            t = ONEXT(t)
        end do
        edge => NULL()
    end subroutine
    

    subroutine test_print(e,c)
        integer(c_intptr_t), intent(in) :: e
        type(c_ptr), intent(in) :: c
        type(edge_struct), pointer :: edge
        edge => deref(e)
        print *, "Pass through edge at: " ,e
        edge => NULL()
    end subroutine 
    
    
    
    
    ! Ported from C to Fortran by Erik Schellenberger @ Siemens Energy Global Gmbh, Goerlitz
    
    
    ! Copyright notice:
    !
    ! Copyright 1996 Institute of Computing, Unicamp.
    !
    ! Permission to use this software for any purpose is hereby granted,
    ! provided that any substantial copy or mechanically derived version
    ! of this file that is made available to other parties is accompanied
    ! by this copyright notice in full, and is distributed under these same
    ! terms. 
    !
    ! DISCLAIMER: This software is provided "as is" with no explicit or
    ! implicit warranty of any kind.  Neither the authors nor their
    ! employers can be held responsible for any losses or damages
    ! that might be attributed to its use.
    !
    ! End of copyright notice.
    
end module