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
    
    
    
end module