# GB-IMR
A Fortran module code to calculate the irreducible macroscopic representation of coincident-site-lattice grain boundaries and construct simulation supercell.

Please integrate the module to your own codes! Brief introductions are given in the notes of each module file. 

An example of the program entry to use these modules is given here.

!!!!!!!!!!!!
Program main

!define variables

!read integer transformation matrices for a given grain boundary (2 * 9 integers)

!call find_unique_ratios (generate the list of linear combination trials)

!call resolve_bicrystal_supercells (get the irreducible macroscopic representation)

!extract the supercell data from array (supercell_solution_list and supercell_parameter_list)

!pass the data to LAMMPS, typically by generating a LAMMPS script with the variables calculated by the Fortran code.

end program main
!!!!!!!!!!!!

Code structure may be rewritten in the future.
