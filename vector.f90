module vector
    implicit none
    
    
    type vec2_f64
        real(8) :: arr(2)
    end type
    
    type vec3_f64
        real(8) :: arr(3)
    end type
    
    type vec2_i32
        integer(4) :: arr(2)
    end type
    
    type vec3_i32
        integer(4) :: arr(3)
    end type
        
    real(8), parameter :: EPS = 1e-9
    
    contains 
    
    type(vec2_f64) function vec2_f64_create(x,y) result(vec)
        real(8) :: x,y
        vec%arr(1) = x
        vec%arr(2) = y
    end function
    
    
    type(vec3_f64) function vec3_f64_create(x,y,z) result(vec)
        real(8) :: x,y,z
        vec%arr(1) = x
        vec%arr(2) = y
        vec%arr(3) = z
    end function
    
    type(vec2_i32) function vec2_i32_create(x,y) result(vec)
        integer(4) :: x,y
        vec%arr(1) = x
        vec%arr(2) = y
    end function
    
    type(vec3_i32) function vec3_i32_create(x,y,z) result(vec)
        integer(4) :: x,y,z
        vec%arr(1) = x
        vec%arr(2) = y
        vec%arr(3) = z
    end function
    
    
    
    
    function vec2_f64_add(a,b) result(c)
        type(vec2_f64) :: a,b,c
        c%arr(1) = a%arr(1) + b%arr(1)
        c%arr(2) = a%arr(2) + b%arr(2)
    end function
    
    function vec3_f64_add(a,b) result(c)
        type(vec3_f64) :: a,b,c
        c%arr(1) = a%arr(1) + b%arr(1)
        c%arr(2) = a%arr(2) + b%arr(2)
        c%arr(3) = a%arr(3) + b%arr(3)
    end function
    
    function vec2_i32_add(a,b) result(c)
        type(vec2_i32) :: a,b,c
        c%arr(1) = a%arr(1) + b%arr(1)
        c%arr(2) = a%arr(2) + b%arr(2)
    end function
    
    function vec3_i32_add(a,b) result(c)
        type(vec3_i32) :: a,b,c
        c%arr(1) = a%arr(1) + b%arr(1)
        c%arr(2) = a%arr(2) + b%arr(2)
        c%arr(3) = a%arr(3) + b%arr(3)
    end function
    
    function vec2_f64_sub(a,b) result(c)
        type(vec2_f64) :: a,b,c
        c%arr(1) = a%arr(1) - b%arr(1)
        c%arr(2) = a%arr(2) - b%arr(2)
    end function
    
    function vec3_f64_sub(a,b) result(c)
        type(vec3_f64) :: a,b,c
        c%arr(1) = a%arr(1) - b%arr(1)
        c%arr(2) = a%arr(2) - b%arr(2)
        c%arr(3) = a%arr(3) - b%arr(3)
    end function
    
    function vec2_i32_sub(a,b) result(c)
        type(vec2_i32) :: a,b,c
        c%arr(1) = a%arr(1) - b%arr(1)
        c%arr(2) = a%arr(2) - b%arr(2)
    end function
    
    function vec3_i32_sub(a,b) result(c)
        type(vec3_i32) :: a,b,c
        c%arr(1) = a%arr(1) - b%arr(1)
        c%arr(2) = a%arr(2) - b%arr(2)
        c%arr(3) = a%arr(3) - b%arr(3)
    end function
    
    function vec2_f64_scale(v,s) result (v_)
        type(vec2_f64) :: v, v_
        real(8) :: s
        v_%arr(1) = v%arr(1) * s
        v_%arr(2) = v%arr(2) * s
    end function
    
    function vec3_f64_scale(v,s) result (v_)
        type(vec3_f64) :: v, v_
        real(8) :: s
        v_%arr(1) = v%arr(1) * s
        v_%arr(2) = v%arr(2) * s
        v_%arr(3) = v%arr(3) * s
    end function
    
    
    function vec2_f64_length(v) result (l)
        type(vec2_f64) :: v
        real(8) :: l
        l = SQRT(v%arr(1)**2 + v%arr(2)**2)
    end function
    
    function vec3_f64_length(v) result (l)
        type(vec3_f64) :: v
        real(8) :: l
        l = SQRT(v%arr(1)**2 + v%arr(2)**2 + v%arr(3)**2)
    end function
    
    function vec2_f64_dist(a,b) result (d)
        type(vec2_f64) :: a,b
        real(8) :: d
        d = vec2_f64_length(vec2_f64_sub(a,b))
    end function
    
    function vec3_f64_dist(a,b) result (d)
        type(vec3_f64) :: a,b
        real(8) :: d
        d = vec3_f64_length(vec3_f64_sub(a,b))
    end function
    
    function vec2_f64_norm(v) result (v_)
        type(vec2_f64) :: v, v_
        v_ = vec2_f64_scale(v,1 / vec2_f64_length(v))
    end function
    
    function vec3_f64_norm(v) result (v_)
        type(vec3_f64) :: v, v_
        v_ = vec3_f64_scale(v,1 / vec3_f64_length(v))
    end function
    
    function vec2_f64_dot(a,b) result(d)
        type(vec2_f64) :: a,b
        real(8) :: d
        d = a%arr(1) * b%arr(1) + a%arr(2) * b%arr(2)
    end function
    
    function vec3_f64_dot(a,b) result(d)
        type(vec3_f64) :: a,b
        real(8) :: d
        d = a%arr(1) * b%arr(1) + a%arr(2) * b%arr(2) + a%arr(3) * b%arr(3)
    end function
    
    function vec2_f64_cross(a,b) result (c)
        type(vec2_f64) :: a,b
        type(vec3_f64) :: c
        c%arr(1) = 0
        c%arr(2) = 0
        c%arr(3) = a%arr(1) * b%arr(2) - a%arr(2) * b%arr(1)
    end function
    
    function vec3_f64_cross(a,b) result (c)
        type(vec3_f64) :: a,b,c
        c%arr(1) = a%arr(2) * b%arr(3) - a%arr(3) * b%arr(2)
        c%arr(2) = a%arr(3) * b%arr(1) - a%arr(1) * b%arr(3)
        c%arr(3) = a%arr(1) * b%arr(2) - a%arr(2) * b%arr(1)
    end function
    
    
    
    function vec2_f64_angle(a,b) result (alpha)
        type(vec2_f64) :: a,b
        real(8) :: alpha
        real(8), parameter :: pi = ACOS(REAL(-1.0,8))
        alpha = ACOS(vec2_f64_dot(a,b) / (vec2_f64_length(a) * vec2_f64_length(b)))
    end function
    
    function vec3_f64_angle(a,b) result (alpha)
        type(vec3_f64) :: a,b
        real(8) :: alpha
        real(8), parameter :: pi = ACOS(REAL(-1.0,8))
        alpha = ACOS(vec3_f64_dot(a,b) / (vec3_f64_length(a) * vec3_f64_length(b)))
    end function
    
    function vec2_f64_equals(a,b) result (l)
        type(vec2_f64) :: a,b
        logical :: l
        l = ABS(a%arr(1) - b%arr(1)) < EPS .AND. ABS(a%arr(2) - b%arr(2)) < EPS
    end function
    
    function vec3_f64_equals(a,b) result (l)
        type(vec3_f64) :: a,b
        logical :: l
        l = ABS(a%arr(1) - b%arr(1)) < EPS .AND. ABS(a%arr(2) - b%arr(2)) < EPS .AND. ABS(a%arr(3) - b%arr(3)) < EPS
    end function
    
    function vec3_f64_in_sphere(a,b,c,d) result (l)
        type(vec3_f64) :: a,b,c,d,bc,ca,ab,o,bary
        logical :: l
        real(8) :: l0, l1, l2, r, sum
        
        bc = vec3_f64_sub(b,c)
        ca = vec3_f64_sub(c,a)
        ab = vec3_f64_sub(a,b)
        
        l0 = vec3_f64_length(bc)
        l1 = vec3_f64_length(ca)
        l2 = vec3_f64_length(ab)
        
        bary%arr(1) = (l1**2 + l2**2 - l0**2)
        bary%arr(2) = (l2**2 + l0**2 - l1**2)
        bary%arr(3) = (l0**2 + l1**2 - l2**2)
        
        sum = bary%arr(1) + bary%arr(2) + bary%arr(3)
        
        bary%arr(1) = bary%arr(1) / sum
        bary%arr(2) = bary%arr(2) / sum
        bary%arr(3) = bary%arr(3) / sum
        
        o%arr(1) = bary%arr(1) * a%arr(1) + bary%arr(2) * b%arr(1) + bary%arr(3) * c%arr(1)
        o%arr(1) = bary%arr(1) * a%arr(2) + bary%arr(2) * b%arr(2) + bary%arr(3) * c%arr(2)
        o%arr(1) = bary%arr(1) * a%arr(3) + bary%arr(2) * b%arr(3) + bary%arr(3) * c%arr(3)
        
        r = vec3_f64_dist(o,a)
        
        l = vec3_f64_dist(o,d) < r
    end function
    
    subroutine vec2_f64_min(out,a,b)
        type(vec2_f64), intent(inout) :: out
        type(vec2_f64), intent(in) :: a,b
        out%arr(1) = MIN(a%arr(1), b%arr(1))
        out%arr(2) = MIN(a%arr(2), b%arr(2))
    end subroutine
    
    subroutine vec2_f64_max(out,a,b)
        type(vec2_f64), intent(inout) :: out
        type(vec2_f64), intent(in) :: a,b
        out%arr(1) = MAX(a%arr(1), b%arr(1))
        out%arr(2) = MAX(a%arr(2), b%arr(2))
    end subroutine
    
    
    
    subroutine vec3_f64_min(out,a,b)
        type(vec3_f64), intent(inout) :: out
        type(vec3_f64), intent(in) :: a,b
        out%arr(1) = MIN(a%arr(1), b%arr(1))
        out%arr(2) = MIN(a%arr(2), b%arr(2))
        out%arr(3) = MIN(a%arr(3), b%arr(3))
    end subroutine
    
    subroutine vec3_f64_max(out,a,b)
        type(vec3_f64), intent(inout) :: out
        type(vec3_f64), intent(in) :: a,b
        out%arr(1) = MAX(a%arr(1), b%arr(1))
        out%arr(2) = MAX(a%arr(2), b%arr(2))
        out%arr(3) = MAX(a%arr(3), b%arr(3))
    end subroutine
    
    
end module