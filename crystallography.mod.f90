module Crystallography

use Math_Library

contains


!find the orientation periodicity of [u v w]
subroutine find_orientation_periodicity(lattice_type,lattice_constant,u,v,w,n,periodicity)
character(len=3) :: lattice_type
real(kind=8) :: lattice_constant
integer(kind=8) :: u,v,w,n
real(kind=8) :: periodicity
!temporary
integer(kind=8) :: u1,v1,w1,common_factor,n_odd_indices

!assign values
u1 = u ; v1 = v ; w1 = w
!reduce orientation to the simplest ratio
call reduce_vector(u1,v1,w1,common_factor)
!number of the odd indices
n_odd_indices = mod(abs(u1),2) + mod(abs(v1),2) + mod(abs(w1),2)

    if ( lattice_type .eq. "FCC" .OR. lattice_type .eq. "fcc" ) then
        if ( n_odd_indices .eq. 2) then !two odd indices
            n = 2
        else
            n = 1
        end if
    else if ( lattice_type .eq. "BCC" .OR. lattice_type .eq. "bcc" ) then
        if ( n_odd_indices .eq. 3) then !three odd indices
            n = 2
        else
            n = 1
        end if
    else if ( lattice_type .eq. "DIA" .OR. lattice_type .eq. "dia" ) then
        if ( n_odd_indices .eq. 3) then !three odd indices
            n = 4
        else if ( n_odd_indices .eq. 2) then !two odd indices  
            n = 2
        else
            n = 1 
        end if
    else
        n = 1
    end if
    
    !do not use reduced orientation here
    periodicity = lattice_constant * sqrt( dble(u*u + v*v + w*w) ) / dble(n)
    
end subroutine find_orientation_periodicity

!find the orientation periodicity of orientation pair [u1 v1 w1] and [u2 v2 w2]
subroutine find_orientation_pair_periodicity(lattice_type,lattice_constant,u1,v1,w1,n1,u2,v2,w2,n2,periodicity)
character(len=3) :: lattice_type
real(kind=8) :: lattice_constant
integer(kind=8) :: u1,v1,w1,n1,u2,v2,w2,n2
real(kind=8) :: periodicity,periodicity1,periodicity2

call find_orientation_periodicity(lattice_type,lattice_constant,u1,v1,w1,n1,periodicity1)
call find_orientation_periodicity(lattice_type,lattice_constant,u2,v2,w2,n2,periodicity2)

    if (abs( mod( max(periodicity1, periodicity2), min(periodicity1, periodicity2)) ) .lt. 1E-6) then    
        periodicity = max(periodicity1, periodicity2)
    else 
        print '(A)', "[E] Fractional coordinates exist at <find_orientation_pair_periodicity>"
        print '(2(A,X,3(I8,X)))', "[E] Input orientation pair is",u1,v1,w1,"&",u2,v2,w2
        stop
    end if
    
end subroutine find_orientation_pair_periodicity


!find the plane distance of (h k l)
subroutine find_plane_distance(lattice_type,lattice_constant,h,k,l,m,distance)
character(len=3) :: lattice_type
real(kind=8) :: lattice_constant
integer(kind=8) :: h,k,l,m
real(kind=8) :: distance
!temporary
integer(kind=8) :: h1,k1,l1,common_factor,n_odd_indices,sum_hkl

!assign values
h1 = h ; k1 = k ; l1 = l
!reduce orientation to the simplest ratio
call reduce_vector(h1,k1,l1,common_factor)

    if ( lattice_type .eq. "FCC" .OR. lattice_type .eq. "fcc" ) then
        !number of the odd indices
        n_odd_indices = mod(abs(h1),2) + mod(abs(k1),2) + mod(abs(l1),2)
        if ( n_odd_indices .eq. 2 .OR. n_odd_indices .eq. 1) then !two odd indices
            m = 2
        else
            m = 1
        end if
    end if
    
    
    if ( lattice_type .eq. "BCC" .OR. lattice_type .eq. "bcc" ) then
        !sum of (h k l)
        sum_hkl = h1 + k1 + l1
        if ( mod(abs(sum_hkl),2) .eq. 1) then !h + k + l = odd
            m = 2
        else
            m = 1
        end if
    end if
    
    if ( lattice_type .eq. "DIA" .OR. lattice_type .eq. "dia" ) then
        !number of the odd indices
        n_odd_indices = mod(abs(h1),2) + mod(abs(k1),2) + mod(abs(l1),2)
        if ( n_odd_indices .eq. 3) then !three odd indices
            m = 1
        else if ( n_odd_indices .eq. 2) then !two odd indices  
            m = 2
        else if ( n_odd_indices .eq. 1) then !one odd indices 
            m = 2
        else
            !sum of (h k l)
            sum_hkl = h1 + k1 + l1
            if ( mod(abs(sum_hkl),4) .eq. 0) then !h + k + l = 4 * n
                m = 1
            else
                m = 2
            end if
        end if
    end if
    
    !do not use reduced plane here
    distance = lattice_constant / sqrt( dble(h1*h1 + k1*k1 + l1*l1) ) / dble(m)
    
end subroutine find_plane_distance


!find the plane distance of plane pair (h1 k1 l1) and (h2 k2 l2)
subroutine find_plane_pair_distance(lattice_type,lattice_constant,h1,k1,l1,m1,h2,k2,l2,m2,distance)
character(len=3) :: lattice_type
real(kind=8) :: lattice_constant
integer(kind=8) :: h1,k1,l1,m1,h2,k2,l2,m2
real(kind=8) :: distance,distance1,distance2

call find_plane_distance(lattice_type,lattice_constant,h1,k1,l1,m1,distance1)
call find_plane_distance(lattice_type,lattice_constant,h2,k2,l2,m2,distance2)

if (abs( mod( max(distance1, distance2), min(distance1, distance2) ) ).lt. 1E-6) then
    distance = min(distance1,distance2)
else if (abs( mod( max(distance1, distance2), min(distance1, distance2) ) - min(distance1,distance2)) .lt. 1E-6) then
    distance = min(distance1,distance2)
else 
    print '(A)', "[E] Fractional coordinates exist at <find_plane_pair_distance>"
    print '(2(A,X,3(I8,X)))', "[E] Input plane pair is",h1,k1,l1,"&",h2,k2,l2
    print '(A,X,F8.4,X,F8.4)', "[E] Plane distances are",distance1,distance2
    stop
end if

end subroutine find_plane_pair_distance



!check if [u1 v1 w1] and [u2 v2 w2] belongs to the same crystallographic orientation family
subroutine check_orientation_family(u1,v1,w1,u2,v2,w2,is_orientation_family)
!input orientation
integer(kind=8) :: u1,v1,w1,u2,v2,w2
real(kind=8) :: deviation_angle,min_deviation_angle

integer(kind=8) :: n_symmetry_operator
integer(kind=8),dimension(24,3,3) :: seTM
character(len=18),dimension(24) :: symmetry_operator
character(len=18) :: symmetrically_optimal_operator

integer(kind=8),dimension(3) :: input_orientation,symmetrically_equivalent_orientation
logical :: is_orientation_family

integer(kind=8) :: i
!initialize
call construct_symmetrically_equivalent_transformation_matrices(n_symmetry_operator,seTM,symmetry_operator)

symmetrically_optimal_operator=symmetry_operator(1)
input_orientation(1) = u1 ; input_orientation(2) = v1 ; input_orientation(3) = w1

call find_vector_pair_minimum_angle(u1,v1,w1,u2,v2,w2,min_deviation_angle)
deviation_angle = min_deviation_angle


do i = 1, n_symmetry_operator
    symmetrically_equivalent_orientation = matmul(input_orientation,seTM(i,:,:))

    call find_vector_pair_minimum_angle(symmetrically_equivalent_orientation(1),&
                                        symmetrically_equivalent_orientation(2),&
                                        symmetrically_equivalent_orientation(3),&
                                        u2,v2,w2,deviation_angle)

    if (deviation_angle .lt. 1E-6 .OR. abs(deviation_angle-180) .lt. 1E-6) then
        is_orientation_family = .TRUE.
    else
        is_orientation_family = .FALSE.
    end if
end do

end subroutine check_orientation_family






subroutine construct_symmetrically_equivalent_transformation_matrices(n_symmetry_operator,seTM,symmetry_operator)
!symmetrically equivalent transformation matrices
integer(kind=8) :: n_symmetry_operator
integer(kind=8),dimension(24,3,3) :: seTM
character(len=18),dimension(24) :: symmetry_operator

!number of symmetry operator is 24 for FCC metals
n_symmetry_operator = 24

!set all elements to 0
seTM(:,:,:) = 0

!list & complie 24 transformation matrices below

!1 identity to itself
symmetry_operator( 1) = "no_transformation_"
seTM( 1,:,:) = reshape([ 1,  0,  0, &
                         0,  1,  0, &
                         0,  0,  1  ],[3,3])

!### rotation by <100> orientation family ###
!2 rotate  90 degree by [100]
symmetry_operator( 2) = "_90deg+[100]axis__"
seTM( 2,:,:) = reshape([ 1,  0,  0, &
                         0,  0, -1, &
                         0,  1,  0  ],[3,3])

!3 rotate 180 degree by [100]
symmetry_operator( 3) = "180deg+[100]axis__"
seTM( 3,:,:) = reshape([ 1,  0,  0, &
                         0, -1,  0, &
                         0,  0, -1  ],[3,3])
!4 rotate 270 degree by [100]
symmetry_operator( 4) = "270deg+[100]axis__"
seTM( 4,:,:) = reshape([ 1,  0,  0, &
                         0,  0,  1, &
                         0, -1,  0  ],[3,3])

!5 rotate  90 degree by [010]
symmetry_operator( 5) = "_90deg+[010]axis__"
seTM( 5,:,:) = reshape([ 0,  0,  1, &
                         0,  1,  0, &
                        -1,  0,  0  ],[3,3])
!6 rotate 180 degree by [010]
symmetry_operator( 6) = "180deg+[010]axis__"
seTM( 6,:,:) = reshape([-1,  0,  0, &
                         0,  1,  0, &
                         0,  0, -1  ],[3,3])
!7 rotate 270 degree by [010]
symmetry_operator( 7) = "270deg+[010]axis__"
seTM( 7,:,:) = reshape([ 0,  0, -1, &
                         0,  1,  0, &
                         1,  0,  0  ],[3,3])

!8  rotate  90 degree by [001]
symmetry_operator( 8) = "_90deg+[001]axis__"
seTM( 8,:,:) = reshape([ 0, -1,  0, &
                         1,  0,  0, &
                         0,  0,  1  ],[3,3])
!9  rotate 180 degree by [001]
symmetry_operator( 9) = "180deg+[001]axis__"
seTM( 9,:,:) = reshape([-1,  0,  0, &
                         0, -1,  0, &
                         0,  0,  1  ],[3,3])
!10 rotate 270 degree by [001]
symmetry_operator(10) = "270deg+[001]axis__"
seTM(10,:,:) = reshape([ 0,  1,  0, &
                        -1,  0,  0, &
                         0,  0,  1  ],[3,3])
                          
!### rotation by <110> orientation family ###
!11 rotate 180 degree by [110]
symmetry_operator(11) = "180deg+[110]axis__"
seTM(11,:,:) = reshape([ 0,  1,  0, &
                         1,  0,  0, &
                         0,  0, -1  ],[3,3])
!                          _
!12 rotate 180 degree by [110]
symmetry_operator(12) = "180deg+[1-10]axis_"
seTM(12,:,:) = reshape([ 0, -1,  0, &
                        -1,  0,  0, &
                         0,  0, -1  ],[3,3])
!13 rotate 180 degree by [101]
symmetry_operator(13) = "180deg+[101]axis__"
seTM(13,:,:) = reshape([ 0,  0,  1, &
                         0, -1,  0, &
                         1,  0,  0  ],[3,3])
!                           _
!14 rotate 180 degree by [101]
symmetry_operator(14) = "180deg+[10-1]axis_"
seTM(14,:,:) = reshape([ 0,  0, -1, &
                         0, -1,  0, &
                        -1,  0,  0  ],[3,3])
!15 rotate 180 degree by [011]
symmetry_operator(15) = "180deg+[011]axis__"
seTM(15,:,:) = reshape([-1,  0,  0, &
                         0,  0,  1, &
                         0,  1,  0  ],[3,3])
!                           _
!16 rotate 180 degree by [011]
symmetry_operator(16) = "180deg+[01-1]axis_"
seTM(16,:,:) = reshape([-1,  0,  0, &
                         0,  0, -1, &
                         0, -1,  0  ],[3,3])

!### rotation by <111> orientation family ###
!17 rotate 120 degree by [111]
symmetry_operator(17) = "120deg+[111]axis__"
seTM(17,:,:) = reshape([ 0,  0,  1, &
                         1,  0,  0, &
                         0,  1,  0  ],[3,3])
!18 rotate 240 degree by [111]
symmetry_operator(18) = "240deg+[111]axis__"
seTM(18,:,:) = reshape([ 0,  1,  0, &
                         0,  0,  1, &
                         1,  0,  0  ],[3,3])
!                          __
!19 rotate 120 degree by [111]
symmetry_operator(19) = "120deg+[1-1-1]axis"
seTM(19,:,:) = reshape([ 0, -1,  0, &
                         0,  0, -1, &
                         1,  0,  0  ],[3,3])
!                          __                      
!20 rotate 240 degree by [111]
symmetry_operator(20) = "240deg+[1-1-1]axis"
seTM(20,:,:) = reshape([ 0,  0,  1, &
                        -1,  0,  0, &
                         0, -1,  0  ],[3,3])        
!                         _ _
!21 rotate 120 degree by [111]
symmetry_operator(21) = "120deg+[-11-1]axis"
seTM(21,:,:) = reshape([ 0,  0, -1, &
                         1,  0,  0, &
                         0, -1,  0  ],[3,3])
!                         _ _                    
!22 rotate 240 degree by [111]
symmetry_operator(22) = "240deg+[-11-1}axis"
seTM(22,:,:) = reshape([ 0,  1,  0, &
                         0,  0,  1, &
                        -1,  0,  0  ],[3,3])
!                         __  
!23 rotate 120 degree by [111]
symmetry_operator(23) = "120deg+[-1-11]axis"
seTM(23,:,:) = reshape([ 0, -1,  0, &
                         0,  0,  1, &
                        -1,  0,  0  ],[3,3])
!                         __  
!24 rotate 240 degree by [111]
symmetry_operator(24) = "240deg+[-1-11]axis"
seTM(24,:,:) = reshape([ 0,  0, -1, &
                        -1,  0,  0, &
                         0,  1,  0  ],[3,3])

end subroutine construct_symmetrically_equivalent_transformation_matrices


!find the symmetrical equivalent orientation [su sv sw] of [uvw] that has the minimum deviation angle with [u0,v0,w0]
subroutine find_symmetrically_optimal_orientation(u0,v0,w0,u,v,w,su,sv,sw,min_deviation_angle,symmetrically_optimal_operator)
!input orientation
integer(kind=8) :: u0,v0,w0,u,v,w,su,sv,sw
real(kind=8) :: deviation_angle,min_deviation_angle

integer(kind=8) :: n_symmetry_operator
integer(kind=8),dimension(24,3,3) :: seTM
character(len=18),dimension(24) :: symmetry_operator
character(len=18) :: symmetrically_optimal_operator

integer(kind=8),dimension(3) :: symmetrically_equivalent_orientation,input_orientation
integer(kind=8) :: i

!initialize
call construct_symmetrically_equivalent_transformation_matrices(n_symmetry_operator,seTM,symmetry_operator)

symmetrically_optimal_operator=symmetry_operator(1)
input_orientation(1) = u ; input_orientation(2) = v ; input_orientation(3) = w
!set input as default
su = u ; sv = v ; sw = w
call find_vector_pair_minimum_angle(u,v,w,u0,v0,w0,min_deviation_angle)
deviation_angle = min_deviation_angle



do i = 1, n_symmetry_operator
    symmetrically_equivalent_orientation = matmul(input_orientation,seTM(i,:,:))

    call find_vector_pair_minimum_angle(symmetrically_equivalent_orientation(1),&
                                        symmetrically_equivalent_orientation(2),&
                                        symmetrically_equivalent_orientation(3),&
                                        u0,v0,w0,deviation_angle)

    if (deviation_angle .lt. min_deviation_angle) then
        min_deviation_angle = deviation_angle
        su = symmetrically_equivalent_orientation(1)
        sv = symmetrically_equivalent_orientation(2)
        sw = symmetrically_equivalent_orientation(3)
        symmetrically_optimal_operator=symmetry_operator(i)
    end if
end do

end subroutine find_symmetrically_optimal_orientation

!#available instruction
!"M"/"m" output the true disorientation_angle & disorientation_axis, but the output grain transformation matrices will be those with minimal X Y Z 
!"C"/"c" check if the input grain transformation matrices are orthogonal and right-handed
subroutine find_symmetrically_optimal_transformation_matrices(P_x1,P_x2,P_x3,P_y1,P_y2,P_y3,P_z1,P_z2,P_z3,&
                                                              Q_x1,Q_x2,Q_x3,Q_y1,Q_y2,Q_y3,Q_z1,Q_z2,Q_z3,&
                                                              SP_x1,SP_x2,SP_x3,SP_y1,SP_y2,SP_y3,SP_z1,SP_z2,SP_z3,&
                                                              SQ_x1,SQ_x2,SQ_x3,SQ_y1,SQ_y2,SQ_y3,SQ_z1,SQ_z2,SQ_z3,&
                                                              disorientation_angle,disorientation_axis,&
                                                              symmetrically_optimal_operator,instruction)
                                               
!input transformation matrices
integer(kind=8) :: P_x1,P_x2,P_x3,P_y1,P_y2,P_y3,P_z1,P_z2,P_z3
integer(kind=8) :: Q_x1,Q_x2,Q_x3,Q_y1,Q_y2,Q_y3,Q_z1,Q_z2,Q_z3
!output symmetrically optimal transformation matrices
integer(kind=8) :: SP_x1,SP_x2,SP_x3,SP_y1,SP_y2,SP_y3,SP_z1,SP_z2,SP_z3
integer(kind=8) :: SQ_x1,SQ_x2,SQ_x3,SQ_y1,SQ_y2,SQ_y3,SQ_z1,SQ_z2,SQ_z3
!output disorientation angle calculated from the optimal matrices
real(kind=8) :: disorientation_angle,local_disorientation_angle,min_disorientation_angle
real(kind=8),dimension(3) :: disorientation_axis,min_disorientation_axis
real(kind=8) :: x_deviation_angle,min_x_deviation_angle
real(kind=8) :: y_deviation_angle,min_y_deviation_angle
real(kind=8) :: z_deviation_angle,min_z_deviation_angle
real(kind=8) :: min_deviation_angle
!
integer(kind=8) :: n_symmetry_operator
!symmetrically optimal transformation matrices
integer(kind=8),dimension(24,3,3) :: seTM
character(len=18),dimension(24) :: symmetry_operator
character(len=18) :: symmetrically_optimal_operator

character(len=1) :: instruction
!independent transformation matrices for symmetrical operations (transpose)
integer(kind=8),dimension(3,3) :: iseTM,iseTMt
!transformation matrices for grains p & q, r = rotated, d = difference
integer(kind=8),dimension(3,3) :: TMp,TMq,TMpr,TMqr,TMd
!rotation matrices for grains p & q, r = rotated, d = difference
real(kind=8),dimension(3,3) :: RMp,RMq,RMpr,RMqr,RMd
real(kind=8) :: delta_trace,min_delta_trace,temp
integer(kind=8) :: i

logical :: is_minimal_deviation_angle
real(kind=8) :: pi=3.1415926535,tol=1E-6


SP_x1 = P_x1 ; SP_x2 = P_x2 ; SP_x3 = P_x3
SP_y1 = P_y1 ; SP_y2 = P_y2 ; SP_y3 = P_y3
SP_z1 = P_z1 ; SP_z2 = P_z2 ; SP_z3 = P_z3
        
SQ_x1 = Q_x1 ; SQ_x2 = Q_x2 ; SQ_x3 = Q_x3
SQ_y1 = Q_y1 ; SQ_y2 = Q_y2 ; SQ_y3 = Q_y3
SQ_z1 = Q_z1 ; SQ_z2 = Q_z2 ; SQ_z3 = Q_z3

!check if the input is a legal transformation/rotation matrix
if (instruction .eq. "C" .OR. instruction .eq. "c") then
    call check_matrix(P_x1,P_x2,P_x3,P_y1,P_y2,P_y3,P_z1,P_z2,P_z3,"SYM_OP","E")
    call check_matrix(Q_x1,Q_x2,Q_x3,Q_y1,Q_y2,Q_y3,Q_z1,Q_z2,Q_z3,"SYM_OP","E")
end if

!complie rotation matrix
call calculate_rotation_matrix(P_x1,P_x2,P_x3,P_y1,P_y2,P_y3,P_z1,P_z2,P_z3,RMp,"SYM_OP","N")
call calculate_rotation_matrix(Q_x1,Q_x2,Q_x3,Q_y1,Q_y2,Q_y3,Q_z1,Q_z2,Q_z3,RMq,"SYM_OP","N")
             
!complie transformation matrix                            
call calculate_transformation_matrix(P_x1,P_x2,P_x3,P_y1,P_y2,P_y3,P_z1,P_z2,P_z3,TMp,"SYM_OP","N")
call calculate_transformation_matrix(Q_x1,Q_x2,Q_x3,Q_y1,Q_y2,Q_y3,Q_z1,Q_z2,Q_z3,TMq,"SYM_OP","N")

call construct_symmetrically_equivalent_transformation_matrices(n_symmetry_operator,seTM,symmetry_operator)

disorientation_angle = 999999999.0
min_disorientation_angle = 999999999.0
disorientation_axis(:) = 999999999.0
min_disorientation_axis(:) = 999999999.0
min_x_deviation_angle = 999999999.0
min_y_deviation_angle = 999999999.0
min_z_deviation_angle = 999999999.0
min_deviation_angle = 999999999.0

min_delta_trace = 999999999.0
symmetrically_optimal_operator = symmetry_operator(1)

do i = 1, n_symmetry_operator
    iseTM (:,:) =           seTM(i,:,:)
    iseTMt(:,:) = transpose(seTM(i,:,:))
    TMpr = matmul(TMp,iseTMt)
    TMqr = matmul(TMq,iseTMt)
    call calculate_rotation_matrix(TMpr(1,1),TMpr(1,2),TMpr(1,3),&
                                   TMpr(2,1),TMpr(2,2),TMpr(2,3),&
                                   TMpr(3,1),TMpr(3,2),TMpr(3,3),RMpr,"SYM_OP","N")
                                         
    call calculate_rotation_matrix(TMqr(1,1),TMqr(1,2),TMqr(1,3),&
                                   TMqr(2,1),TMqr(2,2),TMqr(2,3),&
                                   TMqr(3,1),TMqr(3,2),TMqr(3,3),RMqr,"SYM_OP","N")

    !calculate the matrix from p to qs, i.e., p-1qs
    RMd = matmul(transpose(RMqr),RMp)
    
    delta_trace = RMd(1,1) + RMd(2,2) + RMd(3,3)
    
    temp = (delta_trace - 1.0) / 2.0
    temp = max(-1.0, min(1.0, temp))   !limit to  [-1, 1]
    local_disorientation_angle = acos(temp) * 180.0 / 3.1415926
    
    call find_vector_pair_angle(TMp(1,1),TMp(1,2),TMp(1,3),TMqr(1,1),TMqr(1,2),TMqr(1,3),x_deviation_angle)
    call find_vector_pair_angle(TMp(2,1),TMp(2,2),TMp(2,3),TMqr(2,1),TMqr(2,2),TMqr(2,3),y_deviation_angle)
    call find_vector_pair_angle(TMp(3,1),TMp(3,2),TMp(3,3),TMqr(3,1),TMqr(3,2),TMqr(3,3),z_deviation_angle)
    
    if (local_disorientation_angle .le. min_disorientation_angle) then
        min_delta_trace = delta_trace
        min_disorientation_angle = local_disorientation_angle
        min_disorientation_axis(1) = sqrt( (RMd(1,1)-cos(min_disorientation_angle/180*pi)) / &
                                           (1-cos(min_disorientation_angle/180*pi)))
        min_disorientation_axis(2) = sqrt( (RMd(2,2)-cos(min_disorientation_angle/180*pi)) / &
                                           (1-cos(min_disorientation_angle/180*pi)))
        min_disorientation_axis(3) = sqrt( (RMd(3,3)-cos(min_disorientation_angle/180*pi)) / &
                                           (1-cos(min_disorientation_angle/180*pi)))
        if (instruction .ne. "M" .AND. instruction .ne. "m") then
            symmetrically_optimal_operator = symmetry_operator(i)
            !extract the symmetrically optimal transformation matrices
            SP_x1 = TMp(1,1) ; SP_x2 = TMp(1,2) ; SP_x3 = TMp(1,3)
            SP_y1 = TMp(2,1) ; SP_y2 = TMp(2,2) ; SP_y3 = TMp(2,3)
            SP_z1 = TMp(3,1) ; SP_z2 = TMp(3,2) ; SP_z3 = TMp(3,3)
            SQ_x1 = TMqr(1,1) ; SQ_x2 = TMqr(1,2) ; SQ_x3 = TMqr(1,3)
            SQ_y1 = TMqr(2,1) ; SQ_y2 = TMqr(2,2) ; SQ_y3 = TMqr(2,3)
            SQ_z1 = TMqr(3,1) ; SQ_z2 = TMqr(3,2) ; SQ_z3 = TMqr(3,3)
        end if
    end if
    

    if (instruction .eq. "M" .OR. instruction .eq. "m") then
        is_minimal_deviation_angle = .FALSE.
        if (x_deviation_angle - min_deviation_angle .le. tol) then
            is_minimal_deviation_angle = .TRUE.
            min_deviation_angle = x_deviation_angle
        end if
        if (y_deviation_angle - min_deviation_angle .le. tol) then
            is_minimal_deviation_angle = .TRUE.
            min_deviation_angle = y_deviation_angle
        end if
        if (z_deviation_angle - min_deviation_angle .le. tol) then
            is_minimal_deviation_angle = .TRUE.
            min_deviation_angle = z_deviation_angle
        end if
        
        if (is_minimal_deviation_angle .eqv. .TRUE.) then
            symmetrically_optimal_operator = symmetry_operator(i)
            !extract the symmetrically optimal transformation matrices
            SP_x1 = TMp(1,1) ; SP_x2 = TMp(1,2) ; SP_x3 = TMp(1,3)
            SP_y1 = TMp(2,1) ; SP_y2 = TMp(2,2) ; SP_y3 = TMp(2,3)
            SP_z1 = TMp(3,1) ; SP_z2 = TMp(3,2) ; SP_z3 = TMp(3,3)
            SQ_x1 = TMqr(1,1) ; SQ_x2 = TMqr(1,2) ; SQ_x3 = TMqr(1,3)
            SQ_y1 = TMqr(2,1) ; SQ_y2 = TMqr(2,2) ; SQ_y3 = TMqr(2,3)
            SQ_z1 = TMqr(3,1) ; SQ_z2 = TMqr(3,2) ; SQ_z3 = TMqr(3,3)
        end if
    end if
        
end do

!assign to output
disorientation_angle = min_disorientation_angle
disorientation_axis(:) = min_disorientation_axis(:)

end subroutine find_symmetrically_optimal_transformation_matrices
    

end module Crystallography
