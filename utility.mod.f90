module Utility

contains

!return greatest common divisor of a & b
function gcd2(a, b)
    integer(kind=8) :: gcd2
    integer(kind=8) :: a, b
    integer(kind=8) :: x, y, r

    x = abs(a)
    y = abs(b)
    do while (y /= 0)
        r = mod(x, y)
        x = y
        y = r
    end do
    gcd2 = x
end function gcd2

!return greatest common divisor of a & b & c
function gcd3(a, b, c)
    integer(kind=8) :: gcd3
    integer(kind=8) :: a, b, c
    gcd3 = gcd2(gcd2(a,b),c)
end function gcd3

!prepare the 2D integer list for linear combination trials 
subroutine find_unique_ratios(search_range,ratios,n_ratio)
    integer(kind=8) :: search_range,n_ratio
    integer(kind=8),allocatable,dimension(:,:) :: ratios
    integer(kind=8) :: i,j,counter,counter_first_quadrant

    counter = 0
    do i = 1, search_range
        do j = 1, search_range
            if (gcd2(i,j) == 1) then
                counter = counter + 1
            end if
        end do
    end do

    print '(A,X,I10,X,A)', "Find",counter,"+2 irreducible ratios"
    counter_first_quadrant = counter + 2

    n_ratio = 2 * counter_first_quadrant
    allocate(ratios(n_ratio,2))

    counter = 0
    ratios = 0

    do i = 1, search_range
        do j = 1, search_range
            if (gcd2(i,j) == 1) then
                counter = counter + 1
                ratios(counter,1) = i
                ratios(counter,2) = j
            end if
        end do
    end do

    !copy irreducible ratios to output
    ratios(counter+1,1) = 0 ; ratios(counter+1,2) = 1
    ratios(counter+2,1) = 1 ; ratios(counter+2,2) = 0

    j = 1
    do i = 1,counter_first_quadrant
        ratios(i + j * counter_first_quadrant, 1) = -1 * ratios(i,1)
        ratios(i + j * counter_first_quadrant, 2) =      ratios(i,2)
    end do
    
    print '(A,X,I10)', "Total size is",n_ratio

end subroutine find_unique_ratios

!prepare the 3D integer list for linear combination trials 
subroutine find_unique_combinations(search_range,combinations,n_combination)
    integer(kind=8) :: search_range,n_combination
    integer(kind=8),allocatable,dimension(:,:) :: combinations!4*search_range*search_range*search_range+4*1
    integer(kind=8) :: i,j,k,counter,counter_first_quadrant

    counter = 0
    do i = 0, search_range
        do j = 0, search_range
            do k = 1, search_range / 10
                if (i .eq. 0 .AND. j .eq. 0) cycle
                if (gcd3(i, j, k) .eq. 1) then
                    counter = counter + 1
                end if
            end do
        end do
    end do

    print '(A,X,I10,X,A)', "Find",counter,"+4 irreducible ratios"
    counter_first_quadrant = counter + 1

    n_combination = 4 * counter_first_quadrant
    allocate(combinations(n_combination,3))

    counter = 0
    combinations = 0

    do i = 0, search_range
        do j = 0, search_range
            do k = 1, search_range / 10
                if (i .eq. 0 .AND. j .eq. 0) cycle
                if (gcd3(i, j, k) == 1) then
                    counter = counter + 1
                    combinations(counter,1) = i
                    combinations(counter,2) = j
                    combinations(counter,3) = k
                end if
            end do
        end do
    end do
    
    combinations(counter+1,1) = 0 ; combinations(counter+1,2) = 0 ; combinations(counter+1,3) = 1
   !combinations(counter+2,1) = 0 ; combinations(counter+2,2) = 1 ; combinations(counter+2,3) = 0
   !combinations(counter+3,1) = 1 ; combinations(counter+3,2) = 0 ; combinations(counter+3,3) = 0
    
    j = 1
    do i = 1,counter_first_quadrant
        combinations(i + j * counter_first_quadrant,1) = -1 * combinations(i,1)
        combinations(i + j * counter_first_quadrant,2) =      combinations(i,2)
        combinations(i + j * counter_first_quadrant,3) =      combinations(i,3)
    end do
    
    j = 2
    do i = 1,counter_first_quadrant
        combinations(i + j * counter_first_quadrant,1) =      combinations(i,1)
        combinations(i + j * counter_first_quadrant,2) = -1 * combinations(i,2)
        combinations(i + j * counter_first_quadrant,3) =      combinations(i,3)
    end do
    
    j = 3
    do i = 1,counter_first_quadrant
        combinations(i + j * counter_first_quadrant,1) = -1 * combinations(i,1)
        combinations(i + j * counter_first_quadrant,2) = -1 * combinations(i,2)
        combinations(i + j * counter_first_quadrant,3) =      combinations(i,3)
    end do
    
    print '(A,X,I10)', "Total size is",n_combination

end subroutine find_unique_combinations



!maintain a list of the unique candidate orientation pairs with the minimal periodicity
subroutine construct_orientation_list(n_list,list_size,property_size,vector_list,property_list,vector_input,property_input,&
                                      sort_property1_index,sort_property2_index)

integer(kind=8) :: n_list,list_size,property_size
integer(kind=8),dimension(n_list,list_size) :: vector_list !update
real(kind=8),dimension   (n_list,property_size) :: property_list !update
integer(kind=8),dimension(list_size) :: vector_input
real(kind=8),dimension   (property_size) :: property_input !1 is always the length
integer(kind=8) :: sort_property1_index,sort_property2_index !index of the property list to sort

!temporary
integer(kind=8) :: i,insert_index
integer(kind=8),dimension(n_list,3) :: vector_list_short !update

do i = 1,n_list
    vector_list_short(i,1) = vector_list(i,1)
    vector_list_short(i,2) = vector_list(i,2)
    vector_list_short(i,3) = vector_list(i,3)
end do

!sort parameter 1 is too large or two small
if (property_input(1) .lt. 1E-6 .OR. property_input(1) .gt. 1E8) return

!unduplicate, or linearly independent with the existing list, 
do i = 1, n_list
    if (abs(property_list(i,1) - property_input(1)) < 1E-6 .AND. &
        abs(dble(vector_list_short(i,2)*vector_input(3) - vector_list_short(i,3)*vector_input(2))) < 1E-6 .AND. &
        abs(dble(vector_list_short(i,3)*vector_input(1) - vector_list_short(i,1)*vector_input(3))) < 1E-6 .AND. &
        abs(dble(vector_list_short(i,1)*vector_input(2) - vector_list_short(i,2)*vector_input(1))) < 1E-6 ) then
        return
    end if
end do

insert_index = n_list + 1

!determine the insert index in the list
do i = 1, n_list
    if (property_input(sort_property1_index) - property_list(i,sort_property1_index) .lt. -1E-6) then !the minimum property 1 in the list
        insert_index = i
        exit
    else if (abs(property_input(sort_property1_index) - property_list(i,sort_property1_index)) .lt. 1E-6) then !property 1 is equal, compare property 2
        if (property_input(sort_property2_index) .lt. property_list(i,sort_property2_index)) then
            insert_index = i
            exit
        end if
    else
    !do nothing
    end if
end do


!do not insert if out of the list
if (insert_index .ge. n_list + 1) return

!move the element after insert_index
do i = n_list, insert_index + 1, -1
    vector_list(i,:) = vector_list(i-1, :)
    property_list(i,:) = property_list(i-1,:)
end do

!insert
vector_list(insert_index,:) = vector_input(:)
property_list(insert_index,:) = property_input(:)

end subroutine construct_orientation_list



end module Utility
