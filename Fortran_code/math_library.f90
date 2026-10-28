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



subroutine reduce_matrix(rx1,rx2,rx3,ry1,ry2,ry3,rz1,rz2,rz3,x_common_factor,y_common_factor,z_common_factor)
                         
integer(kind=8) :: rx1,rx2,rx3,ry1,ry2,ry3,rz1,rz2,rz3
integer(kind=8) :: x_common_factor,y_common_factor,z_common_factor
    !get common factors
    call find_common_factor(rx1,rx2,rx3,x_common_factor)
    call find_common_factor(ry1,ry2,ry3,y_common_factor)
    call find_common_factor(rz1,rz2,rz3,z_common_factor)
    !reduce rotation matrix
    rx1 = rx1/x_common_factor ; rx2 = rx2/x_common_factor ; rx3 = rx3/x_common_factor
    ry1 = ry1/y_common_factor ; ry2 = ry2/y_common_factor ; ry3 = ry3/y_common_factor
    rz1 = rz1/z_common_factor ; rz2 = rz2/z_common_factor ; rz3 = rz3/z_common_factor
end subroutine reduce_matrix

!#available instruction
!"P"/"p" ONLY reduce matrix P if Px = kQx / Py = kQy / Pz = kQz
!"Q"/"q" ONLY reduce matrix Q if Qx = kPx / Qy = kPy / Qz = kPz
!Other instructions, reduce both matrices by their common factors
subroutine reduce_matrices(P_x1,P_x2,P_x3,P_y1,P_y2,P_y3,P_z1,P_z2,P_z3,&
                           Q_x1,Q_x2,Q_x3,Q_y1,Q_y2,Q_y3,Q_z1,Q_z2,Q_z3,&
                           P_x_common_factor,P_y_common_factor,P_z_common_factor,&
                           Q_x_common_factor,Q_y_common_factor,Q_z_common_factor,&
                           instruction)

integer(kind=8) :: P_x1,P_x2,P_x3,P_y1,P_y2,P_y3,P_z1,P_z2,P_z3
integer(kind=8) :: Q_x1,Q_x2,Q_x3,Q_y1,Q_y2,Q_y3,Q_z1,Q_z2,Q_z3
integer(kind=8) :: P_x_common_factor,P_y_common_factor,P_z_common_factor
integer(kind=8) :: Q_x_common_factor,Q_y_common_factor,Q_z_common_factor
character(len=1):: instruction
!temporary
integer(kind=8) :: P_common_factor,Q_common_factor,common_factor

if (instruction .eq. "P" .OR. instruction .eq. "p") then !ONLY reduce P
    
    !get X common factors
    call find_common_factor(P_x1,P_x2,P_x3,P_common_factor)
    call find_common_factor(Q_x1,Q_x2,Q_x3,Q_common_factor)
    if (P_common_factor .ge. Q_common_factor .AND. mod(P_common_factor,Q_common_factor) .eq. 0) then
        P_x_common_factor = P_common_factor/Q_common_factor
        Q_x_common_factor = 1
        P_x1 = P_x1/P_x_common_factor ; P_x2 = P_x2/P_x_common_factor ; P_x3 = P_x3/P_x_common_factor
    else
        P_x_common_factor = 1 ; Q_x_common_factor = 1
    end if
    !get Y common factors
    call find_common_factor(P_y1,P_y2,P_y3,P_common_factor)
    call find_common_factor(Q_y1,Q_y2,Q_y3,Q_common_factor)
    if (P_common_factor .ge. Q_common_factor .AND. mod(P_common_factor,Q_common_factor) .eq. 0) then
        P_y_common_factor = P_common_factor/Q_common_factor
        Q_y_common_factor = 1
        P_y1 = P_y1/P_y_common_factor ; P_y2 = P_y2/P_y_common_factor ; P_y3 = P_y3/P_y_common_factor
    else
        P_y_common_factor = 1 ; Q_y_common_factor = 1
    end if
    !get Z common factors
    call find_common_factor(P_z1,P_z2,P_z3,P_common_factor)
    call find_common_factor(Q_z1,Q_z2,Q_z3,Q_common_factor)
    if (P_common_factor .ge. Q_common_factor .AND. mod(P_common_factor,Q_common_factor) .eq. 0) then
        P_z_common_factor = P_common_factor/Q_common_factor
        Q_z_common_factor = 1
        P_z1 = P_z1/P_z_common_factor ; P_z2 = P_z2/P_z_common_factor ; P_z3 = P_z3/P_z_common_factor
    else
        P_z_common_factor = 1 ; Q_z_common_factor = 1
    end if
    
else if (instruction .eq. "Q" .OR. instruction .eq. "q") then !ONLY reduce Q
    
    !get X common factors
    call find_common_factor(P_x1,P_x2,P_x3,P_common_factor)
    call find_common_factor(Q_x1,Q_x2,Q_x3,Q_common_factor)
    if (Q_common_factor .ge. P_common_factor .AND. mod(Q_common_factor,P_common_factor) .eq. 0) then
        P_x_common_factor = 1
        Q_x_common_factor = Q_common_factor/P_common_factor
        Q_x1 = Q_x1/Q_x_common_factor ; Q_x2 = Q_x2/Q_x_common_factor ; Q_x3 = Q_x3/Q_x_common_factor
    else
        P_x_common_factor = 1 ; Q_x_common_factor = 1
    end if
    !get Y common factors
    call find_common_factor(P_y1,P_y2,P_y3,P_common_factor)
    call find_common_factor(Q_y1,Q_y2,Q_y3,Q_common_factor)
    if (Q_common_factor .ge. P_common_factor .AND. mod(Q_common_factor,P_common_factor) .eq. 0) then
        P_y_common_factor = 1
        Q_y_common_factor = Q_common_factor/P_common_factor
        Q_y1 = Q_y1/Q_y_common_factor ; Q_y2 = Q_y2/Q_y_common_factor ; Q_y3 = Q_y3/Q_y_common_factor
    else
        P_y_common_factor = 1 ; Q_y_common_factor = 1
    end if
    !get Z common factors
    call find_common_factor(P_z1,P_z2,P_z3,P_common_factor)
    call find_common_factor(Q_z1,Q_z2,Q_z3,Q_common_factor)
    if (Q_common_factor .ge. P_common_factor .AND. mod(Q_common_factor,P_common_factor) .eq. 0) then
        Q_z_common_factor = 1
        Q_z_common_factor = Q_common_factor/P_common_factor
        Q_z1 = Q_z1/Q_z_common_factor ; Q_z2 = Q_z2/Q_z_common_factor ; Q_z3 = Q_z3/Q_z_common_factor
    else
        P_z_common_factor = 1 ; Q_z_common_factor = 1
    end if
    
else
    !get X common factors
    call find_common_factor(P_x1,P_x2,P_x3,P_common_factor)
    call find_common_factor(Q_x1,Q_x2,Q_x3,Q_common_factor)
    call find_common_factor(P_common_factor,Q_common_factor,P_common_factor,common_factor)
    !reduce X orientations
    P_x1 = P_x1/common_factor ; P_x2 = P_x2/common_factor ; P_x3 = P_x3/common_factor ; P_x_common_factor = common_factor
    Q_x1 = Q_x1/common_factor ; Q_x2 = Q_x2/common_factor ; Q_x3 = Q_x3/common_factor ; Q_x_common_factor = common_factor
    !get Y common factors
    call find_common_factor(P_y1,P_y2,P_y3,P_common_factor)
    call find_common_factor(Q_y1,Q_y2,Q_y3,Q_common_factor)
    call find_common_factor(P_common_factor,Q_common_factor,P_common_factor,common_factor)
    !reduce Y orientations
    P_y1 = P_y1/common_factor ; P_y2 = P_y2/common_factor ; P_y3 = P_y3/common_factor ; P_y_common_factor = common_factor
    Q_y1 = Q_y1/common_factor ; Q_y2 = Q_y2/common_factor ; Q_y3 = Q_y3/common_factor ; Q_y_common_factor = common_factor
    !get Z common factors
    call find_common_factor(P_z1,P_z2,P_z3,P_common_factor)
    call find_common_factor(Q_z1,Q_z2,Q_z3,Q_common_factor)
    call find_common_factor(P_common_factor,Q_common_factor,P_common_factor,common_factor)
    !reduce Z orientations
    P_z1 = P_z1/common_factor ; P_z2 = P_z2/common_factor ; P_z3 = P_z3/common_factor ; P_z_common_factor = common_factor
    Q_z1 = Q_z1/common_factor ; Q_z2 = Q_z2/common_factor ; Q_z3 = Q_z3/common_factor ; Q_z_common_factor = common_factor
end if
end subroutine reduce_matrices

subroutine rescale_matrices(IG1_x1,IG1_x2,IG1_x3,IG1_y1,IG1_y2,IG1_y3,IG1_z1,IG1_z2,IG1_z3,&
                            IG2_x1,IG2_x2,IG2_x3,IG2_y1,IG2_y2,IG2_y3,IG2_z1,IG2_z2,IG2_z3,&
                            RG1_x1,RG1_x2,RG1_x3,RG1_y1,RG1_y2,RG1_y3,RG1_z1,RG1_z2,RG1_z3,&
                            RG2_x1,RG2_x2,RG2_x3,RG2_y1,RG2_y2,RG2_y3,RG2_z1,RG2_z2,RG2_z3)
                            
integer(kind=8) :: IG1_x1,IG1_x2,IG1_x3,IG1_y1,IG1_y2,IG1_y3,IG1_z1,IG1_z2,IG1_z3
integer(kind=8) :: IG2_x1,IG2_x2,IG2_x3,IG2_y1,IG2_y2,IG2_y3,IG2_z1,IG2_z2,IG2_z3

integer(kind=8) :: RG1_x1,RG1_x2,RG1_x3,RG1_y1,RG1_y2,RG1_y3,RG1_z1,RG1_z2,RG1_z3
integer(kind=8) :: RG2_x1,RG2_x2,RG2_x3,RG2_y1,RG2_y2,RG2_y3,RG2_z1,RG2_z2,RG2_z3

integer(kind=8) :: square_G1,square_G2,least_common_multiple
integer(kind=8) :: common_factor_1,common_factor_2,common_factor

RG1_x1 = IG1_x1 ; RG1_x2 = IG1_x2 ; RG1_x3 = IG1_x3
RG1_y1 = IG1_y1 ; RG1_y2 = IG1_y2 ; RG1_y3 = IG1_y3
RG1_z1 = IG1_z1 ; RG1_z2 = IG1_z2 ; RG1_z3 = IG1_z3

RG2_x1 = IG2_x1 ; RG2_x2 = IG2_x2 ; RG2_x3 = IG2_x3
RG2_y1 = IG2_y1 ; RG2_y2 = IG2_y2 ; RG2_y3 = IG2_y3
RG2_z1 = IG2_z1 ; RG2_z2 = IG2_z2 ; RG2_z3 = IG2_z3

square_G1 = RG1_x1*RG1_x1 + RG1_x2*RG1_x2 + RG1_x3*RG1_x3
square_G2 = RG2_x1*RG2_x1 + RG2_x2*RG2_x2 + RG2_x3*RG2_x3

!1st evaluation X
if (square_G1 .eq. square_G2) then
else 
    call find_least_common_multiple(square_G1,square_G2,least_common_multiple)
    RG1_x1 = RG1_x1*nint( sqrt(1.0*least_common_multiple/square_G1) )
    RG1_x2 = RG1_x2*nint( sqrt(1.0*least_common_multiple/square_G1) )
    RG1_x3 = RG1_x3*nint( sqrt(1.0*least_common_multiple/square_G1) )
    
    RG2_x1 = RG2_x1*nint( sqrt(1.0*least_common_multiple/square_G2) )
    RG2_x2 = RG2_x2*nint( sqrt(1.0*least_common_multiple/square_G2) )
    RG2_x3 = RG2_x3*nint( sqrt(1.0*least_common_multiple/square_G2) )
    
    print '(A,I12,X,A,I12,X,I12)', "rescaling X orientations by the least common multiple: ",least_common_multiple,&
                                   "square ratio: ",square_G1,square_G2
    call find_common_factor(RG1_x1,RG1_x2,RG1_x3,common_factor_1)
    call find_common_factor(RG2_x1,RG2_x2,RG2_x3,common_factor_2)
    call find_common_factor(common_factor_1,common_factor_2,common_factor_1,common_factor)

    RG1_x1 = RG1_x1 / common_factor ; RG1_x2 = RG1_x2 / common_factor ; RG1_x3 = RG1_x3 / common_factor
    RG2_x1 = RG2_x1 / common_factor ; RG2_x2 = RG2_x2 / common_factor ; RG2_x3 = RG2_x3 / common_factor
    
    print '(A,X,3(I8,X),A,X,3(I8,X))', "rescaled X orientations",RG1_x1,RG1_x2,RG1_x3,"&",RG2_x1,RG2_x2,RG2_x3
end if
!2nd evaluation X
square_G1 = RG1_x1*RG1_x1 + RG1_x2*RG1_x2 + RG1_x3*RG1_x3
square_G2 = RG2_x1*RG2_x1 + RG2_x2*RG2_x2 + RG2_x3*RG2_x3
if (square_G1 .ne. square_G2) then
    print *, "G1 x size & G2 x size NOT match",square_G1," != ",square_G2
    print '(A)', "[E] Stop at <rescale_matrices>"
    stop
end if
!1st evaluation Y
square_G1 = RG1_y1*RG1_y1 + RG1_y2*RG1_y2 + RG1_y3*RG1_y3
square_G2 = RG2_y1*RG2_y1 + RG2_y2*RG2_y2 + RG2_y3*RG2_y3
if (square_G1 .eq. square_G2) then
else 
    call find_least_common_multiple(square_G1,square_G2,least_common_multiple)
    RG1_y1 = RG1_y1*nint( sqrt(1.0*least_common_multiple/square_G1) )
    RG1_y2 = RG1_y2*nint( sqrt(1.0*least_common_multiple/square_G1) )
    RG1_y3 = RG1_y3*nint( sqrt(1.0*least_common_multiple/square_G1) )
    
    RG2_y1 = RG2_y1*nint( sqrt(1.0*least_common_multiple/square_G2) )
    RG2_y2 = RG2_y2*nint( sqrt(1.0*least_common_multiple/square_G2) )
    RG2_y3 = RG2_y3*nint( sqrt(1.0*least_common_multiple/square_G2) )
    
    print '(A,I12,X,A,I12,X,I12)', "rescaling Y orientations by the least common multiple: ",least_common_multiple,&
                                   "square ratio: ",square_G1,square_G2
    call find_common_factor(RG1_y1,RG1_y2,RG1_y3,common_factor_1)
    call find_common_factor(RG2_y1,RG2_y2,RG2_y3,common_factor_2)
    call find_common_factor(common_factor_1,common_factor_2,common_factor_1,common_factor)

    RG1_y1 = RG1_y1 / common_factor ; RG1_y2 = RG1_y2 / common_factor ; RG1_y3 = RG1_y3 / common_factor
    RG2_y1 = RG2_y1 / common_factor ; RG2_y2 = RG2_y2 / common_factor ; RG2_y3 = RG2_y3 / common_factor
    
    print '(A,X,3(I8,X),A,X,3(I8,X))', "rescaled Y orientations",RG1_y1,RG1_y2,RG1_y3,"&",RG2_y1,RG2_y2,RG2_y3
end if
!2nd evaluation Y
square_G1 = RG1_y1*RG1_y1 + RG1_y2*RG1_y2 + RG1_y3*RG1_y3
square_G2 = RG2_y1*RG2_y1 + RG2_y2*RG2_y2 + RG2_y3*RG2_y3
if (square_G1 .ne. square_G2) then
    print *, "G1 y size & G2 y size NOT match",square_G1," != ",square_G2
    print '(A)', "[E] Stop at <rescale_matrices>"
    stop
end if
!1st evaluation Z
square_G1 = RG1_z1*RG1_z1 + RG1_z2*RG1_z2 + RG1_z3*RG1_z3
square_G2 = RG2_z1*RG2_z1 + RG2_z2*RG2_z2 + RG2_z3*RG2_z3
if (square_G1 .eq. square_G2) then
else 
    call find_least_common_multiple(square_G1,square_G2,least_common_multiple)
    RG1_z1 = RG1_z1*nint( sqrt(1.0*least_common_multiple/square_G1) )
    RG1_z2 = RG1_z2*nint( sqrt(1.0*least_common_multiple/square_G1) )
    RG1_z3 = RG1_z3*nint( sqrt(1.0*least_common_multiple/square_G1) )
    
    RG2_z1 = RG2_z1*nint( sqrt(1.0*least_common_multiple/square_G2) )
    RG2_z2 = RG2_z2*nint( sqrt(1.0*least_common_multiple/square_G2) )
    RG2_z3 = RG2_z3*nint( sqrt(1.0*least_common_multiple/square_G2) )
    
    print '(A,I12,X,A,I12,X,I12)', "rescaling Z orientations by the least common multiple: ",least_common_multiple,&
                                   "square ratio: ",square_G1,square_G2
    call find_common_factor(RG1_z1,RG1_z2,RG1_z3,common_factor_1)
    call find_common_factor(RG2_z1,RG2_z2,RG2_z3,common_factor_2)
    call find_common_factor(common_factor_1,common_factor_2,common_factor_1,common_factor)

    RG1_z1 = RG1_z1 / common_factor ; RG1_z2 = RG1_z2 / common_factor ; RG1_z3 = RG1_z3 / common_factor
    RG2_z1 = RG2_z1 / common_factor ; RG2_z2 = RG2_z2 / common_factor ; RG2_z3 = RG2_z3 / common_factor
    
    print '(A,X,3(I8,X),A,X,3(I8,X))', "rescaled Z orientations",RG1_z1,RG1_z2,RG1_z3,"&",RG2_z1,RG2_z2,RG2_z3
end if
!2nd evaluation Z
square_G1 = RG1_z1*RG1_z1 + RG1_z2*RG1_z2 + RG1_z3*RG1_z3
square_G2 = RG2_z1*RG2_z1 + RG2_z2*RG2_z2 + RG2_z3*RG2_z3
if (square_G1 .ne. square_G2) then
    print *, "G1 z size & G2 z size NOT match",square_G1," != ",square_G2
    print '(A)', "[E] Stop at <rescale_matrices>"
    stop
end if

end subroutine rescale_matrices


subroutine three_dimensional_vector_projection(x1,x2,x3,y1,y2,y3,z1,z2,z3,v1,v2,v3,v_length,px,py,pz)

    integer(kind=8) :: x1,x2,x3,y1,y2,y3,z1,z2,z3,v1,v2,v3
    
    real(kind=8) :: px,py,pz
    real(kind=8) :: v(3),x(3),y(3),z(3)
    real(kind=8) :: x_length,y_length,z_length,v_length
    real(kind=8) :: dot_vx,dot_vy,dot_vz

    !convert to real
    v = (/ dble(v1), dble(v2), dble(v3) /)
    x = (/ dble(x1), dble(x2), dble(x3) /)
    y = (/ dble(y1), dble(y2), dble(y3) /)
    z = (/ dble(z1), dble(z2), dble(z3) /)
    
    x_length = sqrt(dot_product(x,x))
    y_length = sqrt(dot_product(y,y))
    z_length = sqrt(dot_product(z,z))
    
    !x = x / x_length
    !y = y / y_length
    !z = z / z_length
    
    !dot product
    dot_vx = dot_product(v,x)
    dot_vy = dot_product(v,y)
    dot_vz = dot_product(v,z)
    
    if (sqrt(dot_product(v,v)) <= 0.0d0) then
        px = 0.0d0; py = 0.0d0; pz = 0.0d0
        return
    end if
    
    !vector projection
    if (x_length > 0.0d0) then
        !px = v_length * dot_vx / (sqrt(dot_product(v,v)) * x_length) * x / x_length
        px = v_length * dot_vx / (sqrt(dot_product(v,v)) * x_length)
    else
        px = 0.0d0
    end if

    if (y_length > 0.0d0) then
        !py = v_length * dot_vy / (sqrt(dot_product(v,v)) * y_length) * y / y_length
        py = v_length * dot_vy / (sqrt(dot_product(v,v)) * y_length)
    else
        py = 0.0d0
    end if

    if (z_length > 0.0d0) then
        !pz = v_length * dot_vz / (sqrt(dot_product(v,v)) * z_length) * z / z_length
        pz = v_length * dot_vz / (sqrt(dot_product(v,v)) * z_length)
    else
        pz = 0.0d0
    end if

end subroutine three_dimensional_vector_projection







subroutine inverse_unimodular_matrix(M, M_inv)
    implicit none
    integer(kind=8), intent(in)  :: M(3,3)
    integer(kind=8), intent(out) :: M_inv(3,3)

    integer(kind=8) :: C(3,3)
    integer(kind=8) :: det

    !------------------------------------------------------------------
    ! 计算代数余子式矩阵 C
    !------------------------------------------------------------------
    C(1,1) =   M(2,2)*M(3,3) - M(2,3)*M(3,2)
    C(1,2) = -(M(2,1)*M(3,3) - M(2,3)*M(3,1))
    C(1,3) =   M(2,1)*M(3,2) - M(2,2)*M(3,1)
    
    C(2,1) = -(M(1,2)*M(3,3) - M(1,3)*M(3,2))
    C(2,2) =   M(1,1)*M(3,3) - M(1,3)*M(3,1)
    C(2,3) = -(M(1,1)*M(3,2) - M(1,2)*M(3,1))

    C(3,1) =   M(1,2)*M(2,3) - M(1,3)*M(2,2)
    C(3,2) = -(M(1,1)*M(2,3) - M(1,3)*M(2,1))
    C(3,3) =   M(1,1)*M(2,2) - M(1,2)*M(2,1)

    !------------------------------------------------------------------
    ! 计算 determinant
    !
    ! det(U) = U(1,1)*C(1,1)
    !        + U(1,2)*C(1,2)
    !        + U(1,3)*C(1,3)
    !------------------------------------------------------------------
    det = M(1,1)*C(1,1) + &
          M(1,2)*C(1,2) + &
          M(1,3)*C(1,3)
    !------------------------------------------------------------------
    ! 检查是否为奇异矩阵
    !------------------------------------------------------------------

    if (det == 0) then
        print *, "ERROR: matrix is singular!"
        print *, "det(U) = 0"
        M_inv = 0
        return
    end if

    !------------------------------------------------------------------
    ! 检查是否为幺模矩阵
    !
    ! 幺模矩阵必须满足:
    !
    !       det(U) = +1 或 -1
    !------------------------------------------------------------------
    if (abs(det) /= 1) then
        print *, "ERROR: matrix is not unimodular!"
        print *, "det(U) = ", det
        M_inv = 0
        return
    end if

    !------------------------------------------------------------------
    ! adj(U) = transpose(C)
    !
    ! U^(-1) = adj(U) / det(U)
    !------------------------------------------------------------------
    M_inv = transpose(C) / det
end subroutine inverse_unimodular_matrix

subroutine inverse_diagonal_matrix(M, M_inv)
    implicit none
    real(kind=8), intent(in)  :: M(3,3)
    real(kind=8), intent(out) :: M_inv(3,3)
    integer(kind=8) :: i
    M_inv = 0.0d0
    do i = 1, 3
        M_inv(i,i) = 1.0d0 / M(i,i)
    end do
end subroutine inverse_diagonal_matrix

subroutine adjugate_diagonal_matrix(M, M_adj)
    implicit none
    real(kind=8), intent(in)  :: M(3,3)
    real(kind=8), intent(out) :: M_adj(3,3)
    M_adj = 0.0d0
    M_adj(1,1) = M(2,2) * M(3,3)
    M_adj(2,2) = M(1,1) * M(3,3)
    M_adj(3,3) = M(1,1) * M(2,2)
end subroutine adjugate_diagonal_matrix




subroutine adjugate_matrix(M, M_adj, ierr)
    implicit none
    integer(kind=8), intent(in)  :: M(3,3)
    integer(kind=8), intent(out) :: M_adj(3,3)
    integer(kind=8), intent(out) :: ierr

    integer(kind=8) :: C(3,3)   ! 代数余子式矩阵
    integer(kind=8) :: minor(2,2)
    integer(kind=8) :: i, j, r, d, rr, dd

    ierr = 0

    ! ============================================================
    ! 1. 计算代数余子式矩阵 C
    !    C(i,j) = (-1)^(i+j) * 删掉第 i 行第 j 列后的 2x2 行列式
    ! ============================================================
    do i = 1, 3
        do j = 1, 3

            ! ---- 构造删掉第 i 行第 j 列后的 2x2 子矩阵 ----
            rr = 0
            do r = 1, 3
                if (r == i) cycle
                rr = rr + 1
                dd = 0
                do d = 1, 3
                    if (d == j) cycle
                    dd = dd + 1
                    minor(rr, dd) = M(r, d)
                end do
            end do

            ! ---- 2x2 行列式 ----
            C(i,j) = minor(1,1)*minor(2,2) - minor(1,2)*minor(2,1)

            ! ---- 符号 (-1)^(i+j) ----
            if (mod(i + j, 2) == 1) C(i,j) = -C(i,j)

        end do
    end do

    ! ============================================================
    ! 2. 伴随矩阵 = 代数余子式矩阵的转置
    !    adj(M)(i,j) = C(j,i)
    ! ============================================================
    do i = 1, 3
        do j = 1, 3
            M_adj(i,j) = C(j,i)
        end do
    end do

end subroutine adjugate_matrix



subroutine inverse_matrix(M, M_inv, ierr)
    implicit none
    real(kind=8), intent(in)  :: M(3,3)
    real(kind=8), intent(out) :: M_inv(3,3)
    integer(kind=8), intent(out)      :: ierr
    
    real(kind=8) :: Aug(3,6)   ! 增广矩阵 [M | I]
    real(kind=8) :: pivot, factor
    integer(kind=8) :: i, j, k, pivot_row
    
    ierr = 0
    ! ============================================================
    ! 3×3 矩阵求逆（高斯-约当消元法）
    ! 输入: M(3,3) - 待求逆矩阵
    ! 输出: M_inv(3,3) - 逆矩阵
    !       ierr - 错误码: 0=成功, 1=奇异矩阵
    ! ============================================================

    ! ============================================================
    ! 1. 构造增广矩阵 [M | I]
    ! ============================================================
    Aug(:, 1:3) = M
    Aug(:, 4:6) = 0.0d0
    do i = 1, 3
        Aug(i, 3+i) = 1.0d0
    end do
    
    ! ============================================================
    ! 2. 高斯-约当消元
    ! ============================================================
    do i = 1, 3
        ! 选主元
        pivot_row = i
        do j = i+1, 3
            if (abs(Aug(j, i)) > abs(Aug(pivot_row, i))) then
                pivot_row = j
            end if
        end do
        
        ! 检查主元是否为零
        if (abs(Aug(pivot_row, i)) < 1.0d-12) then
            ierr = 1
            return
        end if
        
        ! 交换行
        if (pivot_row /= i) then
            do j = 1, 6
                pivot = Aug(i, j)
                Aug(i, j) = Aug(pivot_row, j)
                Aug(pivot_row, j) = pivot
            end do
        end if
        
        ! 归一化主元行
        pivot = Aug(i, i)
        do j = 1, 6
            Aug(i, j) = Aug(i, j) / pivot
        end do
        
        ! 消去其他行
        do j = 1, 3
            if (j /= i) then
                factor = Aug(j, i)
                do k = 1, 6
                    Aug(j, k) = Aug(j, k) - factor * Aug(i, k)
                end do
            end if
        end do
    end do
    ! ============================================================
    ! 3. 提取逆矩阵
    ! ============================================================
    M_inv = Aug(:, 4:6)
end subroutine inverse_matrix


subroutine inverse_transformation_matrix(M, M_inv, det)
    implicit none
    integer(kind=8), dimension(3,3), intent(in)  :: M
    integer(kind=8), intent(out)                 :: det
    integer(kind=8), dimension(3,3), intent(out) :: M_inv
    M_inv = 0
    det = M(1,1) * (M(2,2)*M(3,3) - M(2,3)*M(3,2)) &
        - M(1,2) * (M(2,1)*M(3,3) - M(2,3)*M(3,1)) &
        + M(1,3) * (M(2,1)*M(3,2) - M(2,2)*M(3,1))
    
    ! 代数余子式
    ! C11 = +(a22*a33 - a23*a32)
    ! C12 = -(a21*a33 - a23*a31)
    ! C13 = +(a21*a32 - a22*a31)
    ! C21 = -(a12*a33 - a13*a32)
    ! C22 = +(a11*a33 - a13*a31)
    ! C23 = -(a11*a32 - a12*a31)
    ! C31 = +(a12*a23 - a13*a22)
    ! C32 = -(a11*a23 - a13*a21)
    ! C33 = +(a11*a22 - a12*a21)

    M_inv(1,1) =  (M(2,2)*M(3,3) - M(2,3)*M(3,2))
    M_inv(1,2) = -(M(1,2)*M(3,3) - M(1,3)*M(3,2))
    M_inv(1,3) =  (M(1,2)*M(2,3) - M(1,3)*M(2,2))
    
    M_inv(2,1) = -(M(2,1)*M(3,3) - M(2,3)*M(3,1))
    M_inv(2,2) =  (M(1,1)*M(3,3) - M(1,3)*M(3,1))
    M_inv(2,3) = -(M(1,1)*M(2,3) - M(1,3)*M(2,1))
    
    M_inv(3,1) =  (M(2,1)*M(3,2) - M(2,2)*M(3,1))
    M_inv(3,2) = -(M(1,1)*M(3,2) - M(1,2)*M(3,1))
    M_inv(3,3) =  (M(1,1)*M(2,2) - M(1,2)*M(2,1))
end subroutine inverse_transformation_matrix


!#available instruction
!"T"/"t" output the transpose of the input matrix
subroutine calculate_transformation_matrix(x1,x2,x3,y1,y2,y3,z1,z2,z3,TM,location,instruction)

integer(kind=8) :: x1,x2,x3,y1,y2,y3,z1,z2,z3
integer(kind=8),dimension(3,3) :: TM
character(len=6) location
character(len=1) instruction
real(kind=8) :: tol = 1E-8

!examine the input transformation matrix
if (sqrt( dble( x1*x1 + x2*x2 + x3*x3 ) ) < tol .OR. &
    sqrt( dble( y1*y1 + y2*y2 + y3*y3 ) ) < tol .OR. &
    sqrt( dble( z1*z1 + z2*z2 + z3*z3 ) ) < tol) then
    
    print '(A,A)', "[E] One of the basis vectors of the rotation matrix is 0, at ",location
    print '(3(I8,1X))',x1,x2,x3
    print '(3(I8,1X))',y1,y2,y3
    print '(3(I8,1X))',z1,z2,z3
    stop
end if

TM(1,1) = x1
TM(1,2) = x2
TM(1,3) = x3
TM(2,1) = y1
TM(2,2) = y2
TM(2,3) = y3
TM(3,1) = z1
TM(3,2) = z2
TM(3,3) = z3

if (instruction .eq. "T" .OR. instruction .eq. "t") then
    TM = transpose(TM)
else

end if

end subroutine calculate_transformation_matrix
!#available instruction
!"T"/"t" output the transpose of the input matrix
subroutine calculate_rotation_matrix(x1,x2,x3,y1,y2,y3,z1,z2,z3,RM,location,instruction)

integer(kind=8) :: x1,x2,x3,y1,y2,y3,z1,z2,z3
real(kind=8),dimension(3,3) :: RM
character(len=6) location
character(len=1) instruction
real(kind=8) :: tol = 1E-8

!examine the input transformation matrix
if (sqrt( dble( x1*x1 + x2*x2 + x3*x3 ) ) < tol .OR. &
    sqrt( dble( y1*y1 + y2*y2 + y3*y3 ) ) < tol .OR. &
    sqrt( dble( z1*z1 + z2*z2 + z3*z3 ) ) < tol) then
    
    print '(A,A)', "[E] One of the basis vectors of the rotation matrix is 0, at ",location
    print '(3(I8,1X))',x1,x2,x3
    print '(3(I8,1X))',y1,y2,y3
    print '(3(I8,1X))',z1,z2,z3
    stop
end if

RM(1,1) = dble(x1) / sqrt( dble( x1*x1 + x2*x2 + x3*x3 ) )
RM(1,2) = dble(x2) / sqrt( dble( x1*x1 + x2*x2 + x3*x3 ) )
RM(1,3) = dble(x3) / sqrt( dble( x1*x1 + x2*x2 + x3*x3 ) )
RM(2,1) = dble(y1) / sqrt( dble( y1*y1 + y2*y2 + y3*y3 ) )
RM(2,2) = dble(y2) / sqrt( dble( y1*y1 + y2*y2 + y3*y3 ) )
RM(2,3) = dble(y3) / sqrt( dble( y1*y1 + y2*y2 + y3*y3 ) )
RM(3,1) = dble(z1) / sqrt( dble( z1*z1 + z2*z2 + z3*z3 ) )
RM(3,2) = dble(z2) / sqrt( dble( z1*z1 + z2*z2 + z3*z3 ) )
RM(3,3) = dble(z3) / sqrt( dble( z1*z1 + z2*z2 + z3*z3 ) )

if (instruction .eq. "T" .OR. instruction .eq. "t") then
    RM = transpose(RM)
else

end if

end subroutine calculate_rotation_matrix

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