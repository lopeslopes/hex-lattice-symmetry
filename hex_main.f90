program hex_lattices
use hex_utils
implicit none

real*16, dimension(:,:), allocatable :: latA1, latA2, latB1, latB2
real*16, dimension(:,:), allocatable :: latAA, latAB, lat_AB_aux, lat_AA_aux
real*16, dimension(3)                :: d1, origin1, origin2, cur_point, eqv_point
real*16, dimension(:), allocatable   :: distsA2, distsB2, distsB1
real*16                              :: angle, z, a, d, cur_dist, eqv_angle
integer                              :: i, j, k, ind_angle, cur_loc
integer                              :: n_points, half_n, same_distA, same_distB
logical                              :: condAB, condBA, condAA, condBB
logical                              :: AB_stacking, hex_center_pivot

! INITIAL DEFINITIONS
n_points = 20000
half_n = n_points/2
a = 2.46e0_16
z = 3.35e0_16
hex_center_pivot = .false.
AB_stacking = .true.

! STACKING AND ORIGIN DEFINITION
if (AB_stacking) then
    write(*,*) "Stacking mode: AB (Bernal stacking)"
else
    write(*,*) "Stacking mode: AA (no displacement)"
endif

if (hex_center_pivot) then
    write(*,*) "Pivot point: empty center of hexagonal cell"
    d = sqrt((a**2)/(2.e0_16*(1.e0_16-cos(2.e0_16*pi/3.e0_16))))
    d1 = [d*cos(pi/6.e0_16), d*sin(pi/6.e0_16), z]
    origin1 = d1 - [0.e0_16, 0.e0_16, z]
    origin2 = d1
else
    write(*,*) "Pivot point: node at origin"
    origin1 = [0.e0_16, 0.e0_16, 0.e0_16]
    origin2 = [0.e0_16, 0.e0_16, z]
endif

! ALLOCATION OF LATTICES AND FIRST CREATION
allocate(latA1(half_n,3))
allocate(latB1(half_n,3))
allocate(latA2(half_n,3))
allocate(latB2(half_n,3))

call create_lattice_eh(latA1, latB1, 0.e0_16, a, .false.)
call create_lattice_eh(latA2, latB2, z      , a, AB_stacking)

! OBTAINING POSSIBLE ANGLES
allocate(distsA2(half_n))
allocate(distsB2(half_n))
allocate(distsB1(half_n))

call distance_from_origin(latA2, origin2, distsA2)
call distance_from_origin(latB2, origin2, distsB2)
call distance_from_origin(latB1, origin1, distsB1)

open(30, file="anglesA.dat")
open(31, file="anglesB.dat")

do i=1, half_n
    cur_dist = distsB1(i)
    cur_point = [latB1(i,1), latB1(i,2), z]
    same_distA = count(abs(distsA2-cur_dist) .lt. tol)
    same_distB = count(abs(distsB2-cur_dist) .lt. tol)
    
    do while (same_distA .gt. 0)
    cur_loc = findloc(distsA2, cur_dist, 1)
    if (cur_loc == 0) then
        exit
    endif
    eqv_point = latA2(cur_loc,:)
    eqv_angle = acos((1.e0_16/(cur_dist**2)) * dot_product(eqv_point, cur_point))
    if (isnan(eqv_angle) .eqv. .false.) then
        write(30,*) eqv_angle
    endif
    same_distA = same_distA - 1
    distsA2(cur_loc) = 5000.e0_16
    enddo
    
    do while (same_distB .gt. 0)
    cur_loc = findloc(distsB2, cur_dist, 1)
    if (cur_loc == 0) then
        exit
    endif
    eqv_point = latB2(cur_loc,:)
    eqv_angle = acos((1.e0_16/(cur_dist**2)) * dot_product(eqv_point, cur_point))
    if (isnan(eqv_angle) .eqv. .false.) then
        write(31,*) eqv_angle
    endif
    same_distB = same_distB - 1
    distsB2(cur_loc) = 5000.e0_16
    enddo
enddo

close(30)
close(31)

deallocate(distsA2, distsB2, distsB1)

! TESTING OBTAINED ANGLES
!do ind_angle=1, 30
    ! ANGLE DEFINITION
    angle = magic_angle(1)
    write(*,*) "Angle in degrees: ", (angle*180.e0_16)/pi 

    ! ROTATE LATTICE 2
    call rotate_lattice(latA2, angle, origin2)
    call rotate_lattice(latB2, angle, origin2)

    ! FIND OVERLAPPING POINTS FROM LATTICES 1 AND 2, SEPARATING A AND B
    allocate(lat_AA_aux(n_points,3))
    allocate(lat_AB_aux(n_points,3))
    j = 0
    k = 0
    do i=1, half_n
        condAB = exists_in_lattice([latB2(i,1),latB2(i,2),0.e0_16], latA1)
        condBA = exists_in_lattice([latA2(i,1),latA2(i,2),0.e0_16], latB1)
        condAA = exists_in_lattice([latA2(i,1),latA2(i,2),0.e0_16], latA1)
        condBB = exists_in_lattice([latB2(i,1),latB2(i,2),0.e0_16], latB1)
        if (condAB) then
            j = j + 1
            lat_AB_aux(j,:) = latB2(i,:)
        endif
        if (condBA) then
            j = j + 1
            lat_AB_aux(j,:) = latA2(i,:)
        endif
        if (condAA) then
            k = k + 1
            lat_AA_aux(k,:) = latA2(i,:)
        endif
        if (condBB) then
            k = k + 1
            lat_AA_aux(k,:) = latB2(i,:)
        endif
    enddo

    allocate(latAA(k,3))
    allocate(latAB(j,3))
    latAA = lat_AA_aux(1:k,:)
    latAB = lat_AB_aux(1:j,:)
    deallocate(lat_AB_aux, lat_AA_aux)

    write(*,*) "  Number of AB points: ", j
    write(*,*) "  Number of AA points: ", k
!enddo

! WRITE LATTICES TO FILE
call write_lattice(latA1, "latticeA1.dat")
call write_lattice(latB1, "latticeB1.dat")
call write_lattice(latA2, "latticeA2.dat")
call write_lattice(latB2, "latticeB2.dat")
call write_lattice(latAA, "latticeAA.dat")
call write_lattice(latAB, "latticeAB.dat")

! DEALLOCATE
deallocate(latA1, latB1, latA2, latB2, latAA, latAB)

end program
