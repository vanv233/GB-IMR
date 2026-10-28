
program main
    use Math_Library
    use Utility
    use Crystallography
    use Supercell
implicit none
!input, grain transformation matrices P & Q(standard)
integer(kind=8) :: P_x1,P_x2,P_x3,P_y1,P_y2,P_y3,P_z1,P_z2,P_z3
integer(kind=8) :: Q_x1,Q_x2,Q_x3,Q_y1,Q_y2,Q_y3,Q_z1,Q_z2,Q_z3
!input, for search propose
integer(kind=8) :: vector_search_range
integer(kind=8) :: n_vector_ratio
integer(kind=8),allocatable,dimension(:,:) :: vector_ratios
!input, material parameters
character(len=3) :: lattice_type
real(kind=8) :: lattice_constant
integer(kind=8) :: sigma

!number of the orientation list, set to 10
integer(kind=8),parameter :: n_supercell_orientation = 4
!number of the solution list, set to 20
integer(kind=8),parameter :: n_supercell_solution = 6 

character(len=1) :: instruction !ignore unnecessary output for fast computation

!output, the assembled supercell solutions & their parameters, -3=fixed y orientation supercell -2=fixed x orientation supercell -1=orthogonal supercell 0=recommend irreducible supercell 1-n_supercell_solution=all computed irreducible supercells
integer(kind=8),dimension(-3:n_supercell_solution,2,8) :: supercell_solution_list !1,1-4=Pw_orientation + Pw_factor 1,5-8=Qw_orientation + Qw_factor 2,1-4=Pv_orientation + Pv_factor 2,5-8=Qv_orientation + Qv_factor
real(kind=8),dimension(-3:n_supercell_solution,8) :: supercell_parameter_list !1=w_periodicity 2=v_periodicity 3=xw_angle 4=xv_angle 5=wv_angle 6=GB_size 7=miGB_size 8=siGB_size
    
    lattice_type = "FCC"
    lattice_constant = 4.05
    
    !determine it by yourself
    sigma = 5
    
    P_x1 = 3
    P_x2 = 0
    P_x3 = 1
    
    P_y1 = 0
    P_y2 = 1
    P_y3 = 0
    
    P_z1 = -1
    P_z2 = 0
    P_z3 = 3
    
    Q_x1 = 3
    Q_x2 = 0
    Q_x3 = -1
    
    Q_y1 = 0
    Q_y2 = 1
    Q_y3 = 0
    
    Q_z1 = 1
    Q_z2 = 0
    Q_z3 = 3
    
    vector_search_range = 200
    
    call find_unique_ratios(vector_search_range,vector_ratios,n_vector_ratio)
    
    
    call resolve_bicrystal_supercells(P_x1,P_x2,P_x3,P_y1,P_y2,P_y3,P_z1,P_z2,P_z3,&
                                      Q_x1,Q_x2,Q_x3,Q_y1,Q_y2,Q_y3,Q_z1,Q_z2,Q_z3,&
                                      n_vector_ratio,vector_ratios,&
                                      lattice_type,lattice_constant,sigma,n_supercell_orientation,&
                                      n_supercell_solution,supercell_solution_list,supercell_parameter_list,&
                                      "N")
    
end program main

