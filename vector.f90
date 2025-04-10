module vector
    implicit none
    
    
    type vec2_f64
        real(8) :: x,y
    end type
    
    type vec3_f64
        real(8) :: x,y,z
    end type
    
    type vec2_i32
        integer(4) :: x,y
    end type
    
    type vec3_i32
        integer(4) :: x,y,z
    end type
        
    real(8), parameter :: EPS = 1e-9
    
    contains 
    
    function vec2_f64_add(a,b) result(c)
        type(vec2_f64) :: a,b,c
        c%x = a%x + b%x
        c%y = a%y + b%y
    end function
    
    function vec3_f64_add(a,b) result(c)
        type(vec3_f64) :: a,b,c
        c%x = a%x + b%x
        c%y = a%y + b%y
        c%z = a%z + b%z
    end function
    
    function vec2_i32_add(a,b) result(c)
        type(vec2_i32) :: a,b,c
        c%x = a%x + b%x
        c%y = a%y + b%y
    end function
    
    function vec3_i32_add(a,b) result(c)
        type(vec3_i32) :: a,b,c
        c%x = a%x + b%x
        c%y = a%y + b%y
        c%z = a%z + b%z
    end function
    
    function vec2_f64_sub(a,b) result(c)
        type(vec2_f64) :: a,b,c
        c%x = a%x - b%x
        c%y = a%y - b%y
    end function
    
    function vec3_f64_sub(a,b) result(c)
        type(vec3_f64) :: a,b,c
        c%x = a%x - b%x
        c%y = a%y - b%y
        c%z = a%z - b%z
    end function
    
    function vec2_i32_sub(a,b) result(c)
        type(vec2_i32) :: a,b,c
        c%x = a%x - b%x
        c%y = a%y - b%y
    end function
    
    function vec3_i32_sub(a,b) result(c)
        type(vec3_i32) :: a,b,c
        c%x = a%x - b%x
        c%y = a%y - b%y
        c%z = a%z - b%z
    end function
    
    function vec2_f64_scale(v,s) result (v_)
        type(vec2_f64) :: v, v_
        real(8) :: s
        v_%x = v%x * s
        v_%y = v%y * s
    end function
    
    function vec3_f64_scale(v,s) result (v_)
        type(vec3_f64) :: v, v_
        real(8) :: s
        v_%x = v%x * s
        v_%y = v%y * s
        v_%z = v%z * s
    end function
    
    
    function vec2_f64_length(v) result (l)
        type(vec2_f64) :: v
        real(8) :: l
        l = SQRT(v%x**2 + v%y**2)
    end function
    
    function vec3_f64_length(v) result (l)
        type(vec3_f64) :: v
        real(8) :: l
        l = SQRT(v%x**2 + v%y**2 + v%z**2)
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
        d = a%x * b%x + a%y * b%y
    end function
    
    function vec3_f64_dot(a,b) result(d)
        type(vec3_f64) :: a,b
        real(8) :: d
        d = a%x * b%x + a%y * b%y + a%z * b%z
    end function
    
    function vec2_f64_cross(a,b) result (c)
        type(vec2_f64) :: a,b
        type(vec3_f64) :: c
        c%x = 0
        c%y = 0
        c%z = a%x * b%y - a%y * b%x
    end function
    
    function vec3_f64_cross(a,b) result (c)
        type(vec3_f64) :: a,b,c
        c%x = a%y * b%z - a%z * b%y
        c%y = a%z * b%x - a%x * b%z
        c%z = a%x * b%y - a%y * b%x
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
        l = ABS(a%x - b%x) < EPS .AND. ABS(a%y - b%y) < EPS
    end function
    
    function vec3_f64_equals(a,b) result (l)
        type(vec3_f64) :: a,b
        logical :: l
        l = ABS(a%x - b%x) < EPS .AND. ABS(a%y - b%y) < EPS .AND. ABS(a%z - b%z) < EPS
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
        
        bary%x = (l1**2 + l2**2 - l0**2)
        bary%y = (l2**2 + l0**2 - l1**2)
        bary%z = (l0**2 + l1**2 - l2**2)
        
        sum = bary%x + bary%y + bary%z
        
        bary%x = bary%x / sum
        bary%y = bary%y / sum
        bary%z = bary%z / sum
        
        o%x = bary%x * a%x + bary%y * b%x + bary%z * c%x
        o%x = bary%x * a%y + bary%y * b%y + bary%z * c%y
        o%x = bary%x * a%z + bary%y * b%z + bary%z * c%z
        
        r = vec3_f64_dist(o,a)
        
        l = vec3_f64_dist(o,d) < r
    end function
    
    
end module