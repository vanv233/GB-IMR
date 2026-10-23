module Supercell
    
    use Math_Library
    use Utility
    use Crystallography
    
contains

!resolve the bicrystal supercells for two standard transformation matrices of a coincident-site-lattice grain boundary character

!input:
!P_x1,P_x2,P_x3,P_y1,P_y2,P_y3,P_z1,P_z2,P_z3 (#grain transformation matrix)
!Q_x1,Q_x2,Q_x3,Q_y1,Q_y2,Q_y3,Q_z1,Q_z2,Q_z3,(#grain transformation matrix)
!vector_search_range,n_vector_ratio,vector_ratios (#integer pairs of unique ratios for search)
!lattice_type,lattice_constant,sigma
!n_supercell_orientation (#number of in-boundary-plane lattice orientations with the minimal orientation lengths)
!n_supercell_solution (#number of output supercells)

!output:
!supercell_solution_list (#two lattice orientation pairs and their common orientation reduction factors)
!supercell_parameter_list (#minimal orientation lengths, deviation_angles to input X lattice orientations, relative angles and supercell size parameters)

!#available instruction
!"Q"/"q" quick computation, skip some unnecessary steps
subroutine resolve_bicrystal_supercells(P_x1,P_x2,P_x3,P_y1,P_y2,P_y3,P_z1,P_z2,P_z3,&
                                        Q_x1,Q_x2,Q_x3,Q_y1,Q_y2,Q_y3,Q_z1,Q_z2,Q_z3,&
                                        vector_search_range,n_vector_ratio,vector_ratios,&
                                        lattice_type,lattice_constant,sigma,n_supercell_orientation,&
                                        n_supercell_solution,supercell_solution_list,supercell_parameter_list,&
                                        instruction)

!input, grain transformation matrices P & Q(standard)
integer(kind=8) :: P_x1,P_x2,P_x3,P_y1,P_y2,P_y3,P_z1,P_z2,P_z3
integer(kind=8) :: Q_x1,Q_x2,Q_x3,Q_y1,Q_y2,Q_y3,Q_z1,Q_z2,Q_z3
!input, for search propose
integer(kind=8) :: vector_search_range
integer(kind=8) :: n_vector_ratio
integer(kind=8),dimension(n_vector_ratio,2)::vector_ratios
!input, material parameters
character(len=3) :: lattice_type
real(kind=8) :: lattice_constant
integer(kind=8) :: sigma

!number of the orientation list, set to 10
integer(kind=8) :: n_supercell_orientation
!number of the solution list, set to 20
integer(kind=8) :: n_supercell_solution

character(len=1) :: instruction !ignore unnecessary output for fast computation

!output, the assembled supercell solutions & their parameters, -3=fixed y orientation supercell -2=fixed x orientation supercell -1=orthogonal supercell 0=recommend irreducible supercell 1-n_supercell_solution=all computed irreducible supercells
integer(kind=8),dimension(-3:n_supercell_solution,2,8) :: supercell_solution_list !1,1-4=Pw_orientation + Pw_factor 1,5-8=Qw_orientation + Qw_factor 2,1-4=Pv_orientation + Pv_factor 2,5-8=Qv_orientation + Qv_factor
real(kind=8),dimension(-3:n_supercell_solution,8) :: supercell_parameter_list !1=w_periodicity 2=v_periodicity 3=xw_angle 4=xv_angle 5=wv_angle 6=GB_size 7=miGB_size 8=siGB_size

!temporary
integer(kind=8) :: P_w1,P_w2,P_w3,P_wn,Q_w1,Q_w2,Q_w3,Q_wn,P_v1,P_v2,P_v3,P_vn,Q_v1,Q_v2,Q_v3,Q_vn
real(kind=8)    :: w_periodicity,v_periodicity,xw_angle,xv_angle,wv_angle,GB_size,miGB_size,siGB_size


!temporary, for orthogonal supercell orientation
integer(kind=8) :: P_o1,P_o2,P_o3,P_on,Q_o1,Q_o2,Q_o3,Q_on
real(kind=8)    :: o_periodicity

!temporary, for irreducible supercell orientation
integer(kind=8) :: P_i1,P_i2,P_i3,P_in,Q_i1,Q_i2,Q_i3,Q_in
real(kind=8)    :: i_periodicity

!store the temporarily obtained supercell edge orientations
integer(kind=8),dimension(8) :: wv_orientation !1,2,3=orientation 4=factor
real(kind=8),dimension(5) :: wv_property !1=periodicity 2=angle to x 3=i 4=j 5=k

!store the minimal supercell edge orientations
integer(kind=8),dimension(n_supercell_orientation,8) :: wv_orientation_list
real(kind=8),dimension(n_supercell_orientation,5) :: wv_property_list

integer(kind=8) :: i,j,k,cross_dot,insert_index,recommend_index
integer(kind=8) :: P_common_factor,Q_common_factor
real(kind=8) :: temp_angle1,temp_angle2
logical :: is_duplicated


if (instruction .eq. "Q" .OR. instruction .eq. "q") then !instruction
!do not check
!do nothing
else
    !check the input
    if ( n_supercell_solution .gt. n_supercell_orientation * (n_supercell_orientation-1) / 2 ) then
        print *, "[E] Input n_supercell_solution is greater than the free combination of the input n_supercell_orientation!"
        print *, "[E] Stop at <resolve_bicrystal_supercells>"
        stop
    end if
    
    if (P_x1*P_x1 + P_x2*P_x2 + P_x3*P_x3 .ne. Q_x1*Q_x1 + Q_x2*Q_x2 + Q_x3*Q_x3) then
        print *, "[E] Input |Px| != |Qx|"
        print *, "Px",P_x1,P_x2,P_x3,"Qx",Q_x1,Q_x2,Q_x3
        print *, "[E] Stop at <resolve_bicrystal_supercells>"
        stop
    end if
    
    if (P_y1*P_y1 + P_y2*P_y2 + P_y3*P_y3 .ne. Q_y1*Q_y1 + Q_y2*Q_y2 + Q_y3*Q_y3) then
        print *, "[E] Input |Py| != |Qy|"
        print *, "Py",P_y1,P_y2,P_y3,"Qy",Q_y1,Q_y2,Q_y3
        print *, "[E] Stop at <resolve_bicrystal_supercells>"
        stop
    end if
end if

!initialize the list
wv_orientation_list = huge(1)
wv_property_list = huge(1.0)
!initialize the list
do k = -3, n_supercell_solution
supercell_solution_list(k,:,:) = 99999999 !huge(1)
supercell_parameter_list(k,:) = 99999999.9 !huge(1.0)
end do

do k = 1,n_vector_ratio
    i = vector_ratios(k,1) ; j = vector_ratios(k,2)
    
    !simultaneously search two quadrants
    !first quadrant
    P_i1 =  i * P_x1 + j * P_y1
    P_i2 =  i * P_x2 + j * P_y2
    P_i3 =  i * P_x3 + j * P_y3
    Q_i1 =  i * Q_x1 + j * Q_y1
    Q_i2 =  i * Q_x2 + j * Q_y2
    Q_i3 =  i * Q_x3 + j * Q_y3

    !reduce the orientation
    call reduce_vectors(P_i1,P_i2,P_i3,Q_i1,Q_i2,Q_i3,P_common_factor,Q_common_factor,"N")
        
    !atomic periodicity along the orientation pair
    call find_orientation_pair_periodicity(lattice_type,lattice_constant,P_i1,P_i2,P_i3,P_in,Q_i1,Q_i2,Q_i3,Q_in,i_periodicity)
    
    !complie data required to construct the list
    wv_orientation(1) = P_i1 ; wv_orientation(2) = P_i2 ; wv_orientation(3) = P_i3 ; wv_orientation(4) = P_in
    wv_orientation(5) = Q_i1 ; wv_orientation(6) = Q_i2 ; wv_orientation(7) = Q_i3 ; wv_orientation(8) = Q_in
    
    wv_property(1) = i_periodicity
    call find_vector_pair_angle(P_i1,P_i2,P_i3,P_x1,P_x2,P_x3,wv_property(2)) !deviation_angle relative to P_x
    wv_property(3) = dble(i) ; wv_property(4) = dble(j) ; wv_property(5) = 0.0

    !input true length
    call construct_orientation_list(n_supercell_orientation,8_8,5_8,wv_orientation_list,wv_property_list,&
                                    wv_orientation,wv_property,1_8,2_8)
                                    

if (instruction .eq. "Q" .OR. instruction .eq. "q") then !instruction
!do not compute orthogonal supercell
!do nothing
else
    !processing the orthogonal supercell now
    call cross_product(P_i1,P_i2,P_i3,P_z1,P_z2,P_z3,P_o1,P_o2,P_o3)
    call cross_product(Q_i1,Q_i2,Q_i3,Q_z1,Q_z2,Q_z3,Q_o1,Q_o2,Q_o3)

    call reduce_vectors(P_o1,P_o2,P_o3,Q_o1,Q_o2,Q_o3,P_common_factor,Q_common_factor,"N")
    
    !atomic periodicity along the orientation pair
    call find_orientation_pair_periodicity(lattice_type,lattice_constant,P_o1,P_o2,P_o3,P_on,Q_o1,Q_o2,Q_o3,Q_on,o_periodicity)


    P_w1 = P_i1 ; P_w2 = P_i2 ; P_w3 = P_i3 ; P_wn = P_in
    Q_w1 = Q_i1 ; Q_w2 = Q_i2 ; Q_w3 = Q_i3 ; Q_wn = Q_in
        
    w_periodicity = i_periodicity
    call find_vector_pair_angle(P_i1,P_i2,P_i3,P_x1,P_x2,P_x3,xw_angle)
        
    P_v1 = P_o1 ; P_v2 = P_o2 ; P_v3 = P_o3 ; P_vn = P_on
    Q_v1 = Q_o1 ; Q_v2 = Q_o2 ; Q_v3 = Q_o3 ; Q_vn = Q_on
        
    v_periodicity = o_periodicity
    xv_angle = xw_angle + 90 ; wv_angle = xv_angle - xw_angle
    GB_size = w_periodicity * v_periodicity
    miGB_size = GB_size / lattice_constant / lattice_constant
    siGB_size = GB_size / lattice_constant / lattice_constant / sigma

    if (miGB_size - supercell_parameter_list(-1,7) .lt. -1E-4) then
        call assign_values(-1_8)
    else if (abs(miGB_size - supercell_parameter_list(-1,7)) .lt. 1E-4) then
        if (xw_angle - supercell_parameter_list(-1,3) .lt. -1E-4) then
            call assign_values(-1_8)
        end if
    else
        !do nothing
    end if
    
end if !instruction
end do







if (instruction .eq. "Q" .OR. instruction .eq. "q") then !instruction
!do not print
!do nothing
else

    print  *, " "
    print '(A,X,I6,X,A)', "[N] Find",n_supercell_orientation,"candidate orientation pairs ..."
    do i = 1,n_supercell_orientation
        print '(I4,A)', i,"-th candidate orientation pair"
        print '(2(A,X,3(I8,X),A,X,I2,X))', &
              "Pw/Pv",wv_orientation_list(i,1:3),"reduction factor",wv_orientation_list(i,4),&
              "Qw/Qv",wv_orientation_list(i,5:7),"reduction factor",wv_orientation_list(i,8)
              
        print '(2(A,X,F10.2,X),A,X,2(F10.2,X))', &
              "periodicity/lattice_constant",wv_property_list(i,1)/lattice_constant,&
              "angle to Px",wv_property_list(i,2),&
              "combination i & j",wv_property_list(i,3),wv_property_list(i,4)
    end do
    print  *, " "
    print '(A)', "[N] Assembling and evaluating irreducible supercells now"


    !assemble fixed x orientation supercell
    do i = 1,n_supercell_orientation
        P_w1 = P_x1 ; P_w2 = P_x2 ; P_w3 = P_x3
        Q_w1 = Q_x1 ; Q_w2 = Q_x2 ; Q_w3 = Q_x3
        call find_orientation_pair_periodicity(lattice_type,lattice_constant,P_w1,P_w2,P_w3,P_wn,Q_w1,Q_w2,Q_w3,Q_wn,w_periodicity)
        
        xw_angle = 0.0d0
        
        P_v1 = wv_orientation_list(i,1) ; P_v2 = wv_orientation_list(i,2) ; P_v3 = wv_orientation_list(i,3) 
        P_vn = wv_orientation_list(i,4)
        Q_v1 = wv_orientation_list(i,5) ; Q_v2 = wv_orientation_list(i,6) ; Q_v3 = wv_orientation_list(i,7)
        Q_vn = wv_orientation_list(i,8)
        
        v_periodicity = wv_property_list(i,1)
        xv_angle = wv_property_list(i,2)
        
        wv_angle = xv_angle - xw_angle
        GB_size = w_periodicity * v_periodicity * sin(wv_angle/180*3.1415926535)
        miGB_size = GB_size / lattice_constant / lattice_constant
        siGB_size = GB_size / lattice_constant / lattice_constant / sigma
        
        !discard the linearly dependent combinations
        if (abs(wv_angle - 0.0) .lt. 1E-4 .OR. abs(wv_angle - 180.0) .lt. 1E-4 .OR. abs(wv_angle + 180.0) .lt. 1E-4) then
            cycle
        end if
        
        if (miGB_size - supercell_parameter_list(-2,7) .lt. -1E-4) then
            call assign_values(-2_8)
        else if (abs(miGB_size - supercell_parameter_list(-2,7)) .lt. 1E-4) then
            if (xv_angle - supercell_parameter_list(-2,4) .lt. -1E-4) then
                call assign_values(-2_8)
            end if
        else
            !do nothing
        end if
        
    end do

    !assemble fixed y orientation supercell
    do i = 1,n_supercell_orientation
        P_w1 = P_y1 ; P_w2 = P_y2 ; P_w3 = P_y3
        Q_w1 = Q_y1 ; Q_w2 = Q_y2 ; Q_w3 = Q_y3
        call find_orientation_pair_periodicity(lattice_type,lattice_constant,P_w1,P_w2,P_w3,P_wn,Q_w1,Q_w2,Q_w3,Q_wn,w_periodicity)

        xw_angle = 90.0d0
        
        P_v1 = wv_orientation_list(i,1) ; P_v2 = wv_orientation_list(i,2) ; P_v3 = wv_orientation_list(i,3) 
        P_vn = wv_orientation_list(i,4)
        Q_v1 = wv_orientation_list(i,5) ; Q_v2 = wv_orientation_list(i,6) ; Q_v3 = wv_orientation_list(i,7)
        Q_vn = wv_orientation_list(i,8)

        v_periodicity = wv_property_list(i,1)
        xv_angle = wv_property_list(i,2)
        
        if (xv_angle .lt. xw_angle) then
            P_v1 = -1*P_v1 ; P_v2 = -1*P_v2 ; P_v3 = -1*P_v3 ; Q_v1 = -1*Q_v1 ; Q_v2 = -1*Q_v2 ; Q_v3 = -1*Q_v3
            xv_angle = xv_angle + 180
        end if
        
        call find_vector_pair_angle(P_w1,P_w2,P_w3,P_v1,P_v2,P_v3,wv_angle)
        GB_size = w_periodicity * v_periodicity * sin(wv_angle/180*3.1415926535)
        miGB_size = GB_size / lattice_constant / lattice_constant
        siGB_size = GB_size / lattice_constant / lattice_constant / sigma
        
        !discard the linearly dependent combinations
        if (abs(wv_angle - 0.0) .lt. 1E-4 .OR. abs(wv_angle - 180.0) .lt. 1E-4 .OR. abs(wv_angle + 180.0) .lt. 1E-4) then
            cycle
        end if
        
        if (miGB_size - supercell_parameter_list(-3,7) .lt. -1E-4) then
            call assign_values(-3_8)
        else if (abs(miGB_size - supercell_parameter_list(-3,7)) .lt. 1E-4) then
            if (xv_angle - supercell_parameter_list(-3,4) .lt. -1E-4) then
                call assign_values(-3_8)
            end if
        else
            !do nothing
        end if
        
    end do

end if !instruction




!assemble the supercell
!get the n_supercell minimal w & v pairs by the free combination of n_orientation * (n_orientation - 1)/2 u_list
do i = 1, n_supercell_orientation - 1
    do j = i + 1, n_supercell_orientation
    
        P_w1 = wv_orientation_list(i,1) ; P_w2 = wv_orientation_list(i,2) ; P_w3 = wv_orientation_list(i,3)
        P_wn = wv_orientation_list(i,4)
        Q_w1 = wv_orientation_list(i,5) ; Q_w2 = wv_orientation_list(i,6) ; Q_w3 = wv_orientation_list(i,7)
        Q_wn = wv_orientation_list(i,8)
        w_periodicity = wv_property_list(i,1)
        xw_angle = wv_property_list(i,2)
        
        P_v1 = wv_orientation_list(j,1) ; P_v2 = wv_orientation_list(j,2) ; P_v3 = wv_orientation_list(j,3)
        P_vn = wv_orientation_list(j,4)
        Q_v1 = wv_orientation_list(j,5) ; Q_v2 = wv_orientation_list(j,6) ; Q_v3 = wv_orientation_list(j,7)
        Q_vn = wv_orientation_list(j,8)
        v_periodicity = wv_property_list(j,1)
        xv_angle = wv_property_list(j,2)
        
        !w must always in the right hand of v. 
        cross_dot = (P_w2*P_v3 - P_w3*P_v2) * P_z1 + &
                    (P_w3*P_v1 - P_w1*P_v3) * P_z2 + &
                    (P_w1*P_v2 - P_w2*P_v1) * P_z3
        !force v to located in the first & second quadrants
        if (cross_dot .lt. 0) then
            P_v1 = -1*P_v1 ; P_v2 = -1*P_v2 ; P_v3 = -1*P_v3 ; Q_v1 = -1*Q_v1 ; Q_v2 = -1*Q_v2 ; Q_v3 = -1*Q_v3
            xv_angle = xv_angle + 180
        end if
        
        !calculate wv_angle
        call find_vector_pair_angle(P_w1,P_w2,P_w3,P_v1,P_v2,P_v3,wv_angle)
        
        is_duplicated = .FALSE.
        !discard the linearly dependent combinations
        if (abs(wv_angle - 0.0) .lt. 1E-4 .OR. abs(wv_angle - 180.0) .lt. 1E-4 .OR. abs(wv_angle + 180.0) .lt. 1E-4) then
            is_duplicated = .TRUE.
        end if
        !discard the duplicated combinations
        do k = 1, n_supercell_solution
            call find_vector_pair_angle(P_w1,P_w2,P_w3,&
                 supercell_solution_list(k,1,1),supercell_solution_list(k,1,2),supercell_solution_list(k,1,3),&
                 temp_angle1)
            call find_vector_pair_angle(P_v1,P_v2,P_v3,&
                 supercell_solution_list(k,2,1),supercell_solution_list(k,2,2),supercell_solution_list(k,2,3),&
                 temp_angle2)
                 
            if ( (abs(temp_angle1 - 0.0  ) .lt. 1E-4  .OR. &
                  abs(temp_angle1 - 180.0) .lt. 1E-4  .OR. &
                  abs(temp_angle1 + 180.0) .lt. 1E-4) .AND. &
                 (abs(temp_angle2 - 0.0  ) .lt. 1E-4  .OR. &
                  abs(temp_angle2 - 180.0) .lt. 1E-4  .OR. &
                  abs(temp_angle2 + 180.0) .lt. 1E-4) ) then
                is_duplicated = .TRUE.
            end if
        end do

        if (is_duplicated .eqv. .TRUE.) then
            cycle
        end if
        
        
        GB_size = w_periodicity * v_periodicity * sin(wv_angle/180*3.1415926535)
        miGB_size = GB_size / lattice_constant / lattice_constant
        siGB_size = GB_size / lattice_constant / lattice_constant / sigma
        
        
        !determine the insert_index if this is acceptable supercell solution in the list
        insert_index = n_supercell_solution + 1
    
        do k = 1, n_supercell_solution
            if (miGB_size - supercell_parameter_list(k,7) .lt. -1E-4) then
                insert_index = k
                exit
            else if (abs(miGB_size - supercell_parameter_list(k,7)) .lt. 1E-4) then
                if (xw_angle - supercell_parameter_list(k,3) .lt. -1E-4) then
                    insert_index = k
                    exit
                else if (abs(xw_angle - supercell_parameter_list(k,3)) .lt. 1E-4) then
                    if (xv_angle .lt. supercell_parameter_list(k,4)) then
                        insert_index = k
                        exit
                    end if
                end if
            end if
        end do
        
        !insert and remove
        if (insert_index .le. n_supercell_solution) then
            do k = n_supercell_solution, insert_index + 1, -1
                supercell_solution_list(k,:,:) = supercell_solution_list(k-1,:,:)
                supercell_parameter_list(k,:) = supercell_parameter_list(k-1,:)
            end do

            call assign_values(insert_index)
            
        end if
        
    end do
end do




if (instruction .eq. "Q" .OR. instruction .eq. "q") then !instruction
    !do not print
    !do not recommend supercell
    supercell_solution_list(0,:,:) = supercell_solution_list(1,:,:)
    supercell_parameter_list(0,:) = supercell_parameter_list(1,:)
    !do not output orthogonal supercell
    supercell_solution_list(-1,:,:) = supercell_solution_list(1,:,:)
    supercell_parameter_list(-1,:) = supercell_parameter_list(1,:)
    !do not output fixed x orientation supercell
    supercell_solution_list(-2,:,:) = supercell_solution_list(1,:,:)
    supercell_parameter_list(-2,:) = supercell_parameter_list(1,:)
    !do not output fixed y orientation supercell
    supercell_solution_list(-3,:,:) = supercell_solution_list(1,:,:)
    supercell_parameter_list(-3,:) = supercell_parameter_list(1,:)
else

    !recommend a supercell for simulation (wv_angle = 90degree, if any)
    recommend_index = 1
    do k = n_supercell_solution,-2,-1 !do not recommend -3
        if ( abs(supercell_parameter_list(k,7)-supercell_parameter_list(1,7)) .lt. 1E-4 .AND. &
             abs(supercell_parameter_list(k,5)-90.0) .lt. 1E-4 ) then
            recommend_index = k
        end if
    end do
    supercell_solution_list(0,:,:) = supercell_solution_list(recommend_index,:,:)
    supercell_parameter_list(0,:) = supercell_parameter_list(recommend_index,:) 
    
    
    print  *, " "
    print '(A,X,I6,X,A,X,I6,A)', "[N] Find",n_supercell_solution,"supercell solutions, and recommend the",recommend_index,"-th"
    print  *, " "
    do k = -3, n_supercell_solution
        if ( k .eq. 0) then
            print '(A)', "[N] Recommended supercell solution"
        else if ( k .eq. -1) then
            print '(I4,A)', k,"-th supercell solution (orthogonal)"
        else if ( k .eq. -2) then
            print '(I4,A)', k,"-th supercell solution (fixed x orientation)"    
        else if ( k .eq. -3) then
            print '(I4,A)', k,"-th supercell solution (fixed y orientation)"    
        else
            print '(I4,A)', k,"-th supercell solution (irreducible)"
        end if
        
        print '(2(A,X,3(I8,X),A,X,I2,X),A,X,F8.2)', &
              "Pw",supercell_solution_list(k,1,1:3),"reduction factor",supercell_solution_list(k,1,4),&
              "Qw",supercell_solution_list(k,1,5:7),"reduction factor",supercell_solution_list(k,1,8),&
              "periodicity/lattice_constant",supercell_parameter_list(k,1)/lattice_constant
        
        print '(2(A,X,3(I8,X),A,X,I2,X),A,X,F8.2)', &
              "Pv",supercell_solution_list(k,2,1:3),"reduction factor",supercell_solution_list(k,2,4),&
              "Qv",supercell_solution_list(k,2,5:7),"reduction factor",supercell_solution_list(k,2,8),&
              "periodicity/lattice_constant",supercell_parameter_list(k,2)/lattice_constant
              
        print '(3(A,X,F6.2,X),A,X,F8.2,X,F8.2,X,F8.4,X)', &
              "xw_angle",supercell_parameter_list(k,3),&
              "xv_angle",supercell_parameter_list(k,4),&
              "wv_angle",supercell_parameter_list(k,5),&
              "GB_size/miGB_size/siGB_size",supercell_parameter_list(k,6:8)
    end do
    print  *, " "
    
end if !instruction




    contains 
    
    subroutine assign_values(insert_index)
        implicit none
        integer(kind=8), intent(in) :: insert_index
        supercell_solution_list(insert_index,1,1) = P_w1
        supercell_solution_list(insert_index,1,2) = P_w2
        supercell_solution_list(insert_index,1,3) = P_w3
        supercell_solution_list(insert_index,1,4) = P_wn
        supercell_solution_list(insert_index,1,5) = Q_w1
        supercell_solution_list(insert_index,1,6) = Q_w2
        supercell_solution_list(insert_index,1,7) = Q_w3 
        supercell_solution_list(insert_index,1,8) = Q_wn
        supercell_solution_list(insert_index,2,1) = P_v1
        supercell_solution_list(insert_index,2,2) = P_v2
        supercell_solution_list(insert_index,2,3) = P_v3
        supercell_solution_list(insert_index,2,4) = P_vn
        supercell_solution_list(insert_index,2,5) = Q_v1
        supercell_solution_list(insert_index,2,6) = Q_v2
        supercell_solution_list(insert_index,2,7) = Q_v3
        supercell_solution_list(insert_index,2,8) = Q_vn
        !1=w_length 2=v_length 3=xw_angle 4=xv_angle 5=wv_angle 6=GB_size 7=miGB_size 8=siGB_size
        supercell_parameter_list(insert_index,1) = w_periodicity
        supercell_parameter_list(insert_index,2) = v_periodicity
        supercell_parameter_list(insert_index,3) = xw_angle
        supercell_parameter_list(insert_index,4) = xv_angle
        supercell_parameter_list(insert_index,5) = wv_angle
        supercell_parameter_list(insert_index,6) = GB_size
        supercell_parameter_list(insert_index,7) = miGB_size
        supercell_parameter_list(insert_index,8) = siGB_size
    end subroutine assign_values
    
end subroutine resolve_bicrystal_supercells

end module Supercell