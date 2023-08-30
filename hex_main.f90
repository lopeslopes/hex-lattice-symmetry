program hex_lattices
use hex_utils
implicit none

real*16, dimension(:,:), allocatable :: latA1, latA2, latB1, latB2
real*16, dimension(3)                :: d1, origin1, origin2
real*16, dimension(:), allocatable   :: distsA, distsB
real*16                              :: angle, z, a, d, cur_dist
integer                              :: i, j, k, max_ind, min_ind, ind_angle
integer                              :: n_points, half_n, same_dist
logical                              :: condAB, condBA, condAA, condBB
logical                              :: AB_stacking, hex_center_pivot

! INITIAL DEFINITIONS
n_points = 24
half_n = n_points/2
a = 2.46e0_16
z = 3.35e0_16
hex_center_pivot = .false.
AB_stacking = .false.

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

ind_angle = 30
! LOOP OVER ANGLES CAN START HERE
!do ind_angle=1, 30
    ! ANGLE DEFINITION
    angle = magic_angle(ind_angle)
    write(*,*) "Angle in degrees: ", (angle*180.e0_16)/pi 

    ! ALLOCATE AND CREATE LATTICES
    allocate(latA1(half_n,3))
    allocate(latB1(half_n,3))
    allocate(latA2(half_n,3))
    allocate(latB2(half_n,3))
    call create_lattice_eh(latA1, latB1, 0.e0_16, a, .false.)
    call create_lattice_eh(latA2, latB2, z      , a, AB_stacking)


    ! TEST SECTION --------------------
    allocate(distsA(half_n))
    allocate(distsB(half_n))

    call distance_from_origin(latA1, origin1, distsA)
    call distance_from_origin(latA2, origin2, distsB)

    same_dist = 0
    do i=1, half_n
        write(*,"(F10.5, 2X, F10.5)") distsA(i), distsB(i)
    enddo

    deallocate(distsA, distsB)
    ! ---------------------------------

    ! ROTATE LATTICE 2
    !call rotate_lattice(latA2, angle, origin2)
    !call rotate_lattice(latB2, angle, origin2)

    !! FIND OVERLAPPING POINTS FROM LATTICES 1 AND 2, SEPARATING A AND B
    !j = 0
    !k = 0
    !max_ind = 20
    !min_ind = 0
    !do i=1, size(latB2,1)
    !    condAB = exists_in_lattice([latB2(i,1), latB2(i,2), 0.e0_16], latA1(i-min_ind:i+max_ind,:))
    !    condBA = exists_in_lattice([latA2(i,1), latA2(i,2), 0.e0_16], latB1(i-min_ind:i+max_ind,:))
    !    condAA = exists_in_lattice([latA2(i,1), latA2(i,2), 0.e0_16], latA1(i-min_ind:i+max_ind,:))
    !    condBB = exists_in_lattice([latB2(i,1), latB2(i,2), 0.e0_16], latB1(i-min_ind:i+max_ind,:))
    !    if (condAB) then
    !        j = j + 1
    !    endif
    !    if (condBA) then
    !        j = j + 1
    !    endif
    !    if (condAA) then
    !        k = k + 1
    !    endif
    !    if (condBB) then
    !        k = k + 1
    !    endif
    !    if (i < 20) then
    !        min_ind = min_ind + 1
    !    else
    !        if (i > n_points - 20) then
    !            max_ind = max_ind - 1
    !        endif
    !    endif
    !enddo

    !write(*,*) "  Number of AB points: ", j
    !write(*,*) "  Number of AA points: ", k

    !! WRITE LATTICES TO FILE
    call write_lattice(latA1, "latticeA1.dat")
    call write_lattice(latB1, "latticeB1.dat")
    call write_lattice(latA2, "latticeA2.dat")
    call write_lattice(latB2, "latticeB2.dat")

    ! DEALLOCATE
    deallocate(latA1, latB1, latA2, latB2)
!enddo

end program
