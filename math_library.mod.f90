module Math_Library

implicit none

contains

subroutine cross_product(v1,v2,v3,v4,v5,v6,rv1,rv2,rv3)
    integer(kind=8) :: v1,v2,v3,v4,v5,v6
    integer(kind=8) :: rv1,rv2,rv3
    rv1 = v2 * v6 - v3 * v5
    rv2 = v3 * v4 - v1 * v6
    rv3 = v1 * v5 - v2 * v4
    
end subroutine cross_product

subroutine find_common_factor(a,b,c,common_factor)

integer(kind=8) :: a,b,c, common_factor,i

common_factor = 1

if ( a .ne. 0 .and. b .ne. 0 .and. c .ne. 0 ) then
    
    common_factor = gcd(gcd(abs(a), abs(b)), abs(c))
    
else if ( a .eq. 0 .and. b .ne. 0 .and. c .ne. 0) then

    do i=1,min(abs(b),abs(c))
        if ( ( mod(abs(b), i) .eq. 0 ) .and. ( mod(abs(c), i) .eq. 0 ) ) then
            common_factor = i
        end if
    end do
    
else if ( a .ne. 0 .and. b .eq. 0 .and. c .ne. 0) then

    do i=1,min(abs(a),abs(c))
        if ( ( mod(abs(a), i) .eq. 0 ) .and. ( mod(abs(c), i) .eq. 0 ) ) then
            common_factor = i
        end if
    end do

else if ( a .ne. 0 .and. b .ne. 0 .and. c .eq. 0) then

    do i=1,min(abs(a),abs(b))
        if ( ( mod(abs(a), i) .eq. 0 ) .and. ( mod(abs(b), i) .eq. 0 ) ) then
            common_factor = i
        end if
    end do
    
else if ( a .ne. 0 .and. b .eq. 0 .and. c .eq. 0) then
            common_factor = abs(a)
else if ( a .eq. 0 .and. b .ne. 0 .and. c .eq. 0) then
            common_factor = abs(b)
else if ( a .eq. 0 .and. b .eq. 0 .and. c .ne. 0) then
            common_factor = abs(c)
else
    !print *, "unexpected vector input for common factor search, will use conventional algorithm!"
    do i=1,min(abs(a),abs(b),abs(c))
        if ( (mod(abs(a), i) .eq. 0 ) .and. ( mod(abs(b), i) .eq. 0 ) .and. ( mod(abs(c), i) .eq. 0 ) ) then
            common_factor = i
        end if
    end do
end if

contains
    function gcd(x, y) result(g)
        integer(kind=8) :: x, y, g
        integer(kind=8) :: r, tx, ty
        tx = x ; ty = y
        do while (ty /= 0)
            r = mod(tx, ty)
            tx = ty
            ty = r
        end do
        g = tx
    end function gcd
end subroutine find_common_factor

subroutine find_least_common_multiple(a,b,least_common_multiple)
    ! input: integers a and b
    implicit none
    integer(kind=8), intent(in)  :: a, b
    integer(kind=8), intent(out) :: least_common_multiple
    integer(kind=8) :: gcd, x, y, r
    
    ! input = 0
    if (a == 0 .or. b == 0) then
        least_common_multiple = 0
        return
    end if
    
    ! 
    x = abs(a)
    y = abs(b)
    do while (y /= 0)
        r = mod(x, y)
        x = y
        y = r
    end do
    gcd = x
    
    ! gets LCM = |a * b| / GCD
    ! 注意：先除后乘，避免整数溢出
    least_common_multiple = abs(a / gcd) * abs(b)
    
end subroutine find_least_common_multiple


subroutine find_vector_pair_minimum_angle(v1,v2,v3,v4,v5,v6,minimum_angle_deg)
    
integer(kind=8) :: v1,v2,v3,v4,v5,v6
real(kind=8) :: minimum_angle_deg,minimum_angle_rad
real(kind=8) :: cos_value
real(kind=8) :: tol = 1e-7

cos_value = 1.0d0*(v1*v4+v2*v5+v3*v6)/sqrt(1.0d0*v1*v1+v2*v2+v3*v3)/sqrt(1.0d0*v4*v4+v5*v5+v6*v6)
cos_value = max(-1.0d0, min(1.0d0, cos_value))
minimum_angle_rad = acos( cos_value )

if (abs(minimum_angle_rad) .lt. tol) then
    minimum_angle_rad = 0.0d0
end if

minimum_angle_deg = abs(minimum_angle_rad/3.1415926535d0*180.0d0)
minimum_angle_deg = min(minimum_angle_deg,180.0d0-minimum_angle_deg)

end subroutine find_vector_pair_minimum_angle

subroutine find_vector_pair_angle(v1,v2,v3,v4,v5,v6,angle_deg)
    
integer(kind=8) :: v1,v2,v3,v4,v5,v6
real(kind=8) :: angle_deg,angle_rad
real(kind=8) :: cos_value
real(kind=8) :: tol = 1e-7

cos_value = dble(v1*v4+v2*v5+v3*v6)/sqrt(dble(v1*v1+v2*v2+v3*v3))/sqrt(dble(v4*v4+v5*v5+v6*v6))
cos_value = max(-1.0d0, min(1.0d0, cos_value))
angle_rad = acos( cos_value )

if (abs(angle_rad) .lt. tol) then
    angle_rad = 0.0d0
end if

angle_deg = angle_rad/3.1415926535d0*180.0d0

end subroutine find_vector_pair_angle



subroutine reduce_vector(v1,v2,v3,v_common_factor)
integer(kind=8) :: v1,v2,v3,v_common_factor
!get common factor
call find_common_factor(v1,v2,v3,v_common_factor)
!reduce vector
v1 = v1/v_common_factor ; v2 = v2/v_common_factor ; v3 = v3/v_common_factor 
end subroutine reduce_vector

!#available instruction
!"P"/"p" ONLY reduce vector p if p = kq
!"Q"/"q" ONLY reduce vector q if q = kp
!Other instructions, reduce both vectors by their common factors
subroutine reduce_vectors(p1,p2,p3,q1,q2,q3,p_common_factor,q_common_factor,instruction)
integer(kind=8) :: p1,p2,p3,q1,q2,q3,p_common_factor,q_common_factor,pq_common_factor
character(len=1):: instruction
!get common factor
call find_common_factor(p1,p2,p3,p_common_factor)
call find_common_factor(q1,q2,q3,q_common_factor)
call find_common_factor(p_common_factor,q_common_factor,p_common_factor,pq_common_factor)

!reduce vector
if (instruction .eq. "P" .OR. instruction .eq. "p") then !ONLY reduce P
    
    if (p_common_factor .ge. q_common_factor .AND. mod(p_common_factor,q_common_factor) .eq. 0) then
        pq_common_factor = p_common_factor/q_common_factor
        p1 = p1/pq_common_factor ; p2 = p2/pq_common_factor ; p3 = p3/pq_common_factor
        p_common_factor = pq_common_factor ; q_common_factor = 1
    else
        p_common_factor = 1 ; q_common_factor = 1
    end if
 
else if (instruction .eq. "Q" .OR. instruction .eq. "q") then !ONLY reduce Q
    
    if (q_common_factor .ge. p_common_factor .AND. mod(q_common_factor,p_common_factor) .eq. 0) then
        pq_common_factor = q_common_factor/p_common_factor
        q1 = q1/pq_common_factor ; q2 = q2/pq_common_factor ; q3 = q3/pq_common_factor
        p_common_factor = 1 ; q_common_factor = pq_common_factor
    else
        p_common_factor = 1 ; q_common_factor = 1
    end if
 
else
    p1 = p1/pq_common_factor ; p2 = p2/pq_common_factor ; p3 = p3/pq_common_factor ; p_common_factor = pq_common_factor
    q1 = q1/pq_common_factor ; q2 = q2/pq_common_factor ; q3 = q3/pq_common_factor ; q_common_factor = pq_common_factor
end if
end subroutine reduce_vectors


!check if the matrix is orthogonal and right-handed

!#available instruction
!if the matrix is NOT orthogonal and right-handed
!"E"/"S" stop if the matrix 
!"X"/"x" inverse x to -x and recheck
!"Y"/"y" inverse y to -y and recheck
!"Z"/"z" inverse z to -z and recheck
subroutine check_matrix(x1,x2,x3,y1,y2,y3,z1,z2,z3,location,instruction)
  
    integer(kind=8) :: x1,x2,x3,y1,y2,y3,z1,z2,z3
    character(len=6) :: location
    character(len=1) :: instruction
    !local variables
    integer(kind=8) :: x(3),y(3),z(3)
    integer(kind=8) :: dot_xy,dot_xz,dot_yz
    integer(kind=16) :: det
    integer(kind=8) :: cross(3)
  
    x = (/ x1, x2, x3 /) ; y = (/ y1, y2, y3 /) ; z = (/ z1, z2, z3 /)
    dot_xy = x(1)*y(1) + x(2)*y(2) + x(3)*y(3)
    dot_xz = x(1)*z(1) + x(2)*z(2) + x(3)*z(3)
    dot_yz = y(1)*z(1) + y(2)*z(2) + y(3)*z(3)

!check if any two are orthogonal
    if (dot_xy .ne. 0) then
        print '(A,A)','[E] Input xy is not equal to 0 at ',location
        print '(A,3(I8,X),A,3(I8,X),I8,X)',"x",x(:),"y",y(:),dot_xy
        stop
    end if
  
    if (dot_xz .ne. 0) then
        print '(A,A)','[E] Input xz is not equal to 0 at ',location
        print '(A,3(I8,X),A,3(I8,X),I8,X)',"x",x(:),"z",z(:),dot_xz
        stop
    end if
  
    if (dot_yz .ne. 0) then
        print '(A,A)','[E] Input yz is not equal to 0 at ',location
        print '(A,3(I8,X),A,3(I8,X),I8,X)',"y",y(:),"z",z(:),dot_yz
        stop
    end if

!check if the three are right-handed
    det = x(1)* (y(2)*z(3) - y(3)*z(2)) &
        - x(2)* (y(1)*z(3) - y(3)*z(1)) &
        + x(3)* (y(1)*z(2) - y(2)*z(1))
  
    if (det .le. 0) then
        print '(A,A)','[W] input x y z is LEFT-handed or coplanar x × y ≠ z, at ',location
        print '(A,3(I12,X),A,3(I12,X),A,3(I12,X))',"x",x(:),"y",y(:),"z",z(:)
        print '(A,A)','The check subroutine will comply with your instruction ',instruction
        call sleep(1)

        if (instruction .eq. "E" .OR. instruction .eq. "S") then
            print '(A)','[E] Stop to check your input'
            stop
        else if (instruction .eq. "X" .OR. instruction .eq. "x") then
            x1 = -1 * x1 ; x2 = -1 * x2 ; x3 = -1 * x3 
        else if (instruction .eq. "Y" .OR. instruction .eq. "y") then
            y1 = -1 * y1 ; y2 = -1 * y2 ; y3 = -1 * y3 
        else if (instruction .eq. "Z" .OR. instruction .eq. "z") then
            z1 = -1 * z1 ; z2 = -1 * z2 ; z3 = -1 * z3 
        else
            stop
        end if
        x = (/ x1, x2, x3 /) ; y = (/ y1, y2, y3 /) ; z = (/ z1, z2, z3 /)
        cross(1) = x(2)*y(3) - x(3)*y(2)
        cross(2) = x(3)*y(1) - x(1)*y(3)
        cross(3) = x(1)*y(2) - x(2)*y(1)
!check if the three are right-handed
        det = x(1)* (y(2)*z(3) - y(3)*z(2)) &
            - x(2)* (y(1)*z(3) - y(3)*z(1)) &
            + x(3)* (y(1)*z(2) - y(2)*z(1))
  
        if (det .le. 0) then
            print '(A,3(I12,X),A,3(I12,X),A,3(I12,X))',"x",x(:),"y",y(:),"z",z(:)
            print '(A,X,I12)',"determinant",det
            print '(A,A)','[E] input x y z is not right-handed even after your instruction, the program will stop'
            stop
        end if
    end if
    
    !Passed check
    print '(A,A)',"[N] Input x y z is orthogonal and right-handed at ",location
  
end subroutine check_matrix



end module Math_Library