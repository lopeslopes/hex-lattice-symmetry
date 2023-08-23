program hex_lattices
use hex_utils
implicit none

real*16, dimension(:,:), allocatable :: lat1, lat2, lat1f, lat2f, lat3f, lat4f
real*16, dimension(:,:), allocatable :: lat1_aux, lat2_aux, lat3_aux, lat4_aux
real*16, dimension(:,:), allocatable :: latA1, latA2, latB1, latB2, latAB
real*16, dimension(3)                :: d1, pivot
real*16                              :: angle, z, a, d
integer                              :: i, j, k, n_points, max_ind, min_ind, ind_angle
logical                              :: cond3, condAB, condBA, condAA, condBB
logical                              :: AB_stacking, hex_center_pivot

! INITIAL DEFINITIONS
n_points = 50000
a = 2.46e0_16
z = 3.35e0_16
hex_center_pivot = .false.
AB_stacking = .true.

if (AB_stacking) then
    write(*,*) "Stacking mode: AB (Bernal stacking)"
else
    write(*,*) "Stacking mode: AA (no displacement)"
endif
if (hex_center_pivot) then
    write(*,*) "Pivot point: empty center of hexagonal cell"
else
    write(*,*) "Pivot point: node at origin"
endif

! LOOP OVER ANGLES CAN START HERE
do ind_angle=25, 35
    ! ANGLE DEFINITION
    angle = magic_angle(ind_angle)
    write(*,*) "Angle in degrees: ", (angle*180.e0_16)/pi 

    ! ALLOCATE AND CREATE LATTICES
    allocate(lat1(n_points,3))
    allocate(lat2(n_points,3))
    call create_lattice_eh(lat1, 0.e0_16, a, .false.)
    call create_lattice_eh(lat2, z      , a, AB_stacking)

    ! ROTATE LATTICE 2
    d = sqrt((a**2)/(2.e0_16*(1.e0_16-cos(2.e0_16*pi/3.e0_16))))
    d1 = [d*cos(pi/6.e0_16), d*sin(pi/6.e0_16), z]
    if (hex_center_pivot) then
        pivot = d1
    else
        pivot = [0.e0_16, 0.e0_16, z]
    endif

    call rotate_lattice(lat2, angle, pivot)

    ! TRIM THE LATTICE 1 ------------------
    !allocate(lat1_aux(n_points, 3))
    !j = 0
    !call trim_lattice(lat1, -10000.e0_16, 10000.e0_16, -10000.e0_16, 10000.e0_16, j, lat1_aux)
    !allocate(lat1f(j,3))
    !lat1f = lat1_aux(1:j,:)
    !deallocate(lat1_aux)
    !deallocate(lat1)

    ! OR DON'T TRIM
    j = size(lat1,1)
    allocate(lat1f(j,3))
    lat1f = lat1
    deallocate(lat1)
    ! -------------------------------------

    ! TRIM THE LATTICE 2 ------------------
    !allocate(lat2_aux(n_points, 3))
    !j = 0
    !call trim_lattice(lat2, -10000.e0_16, 10000.e0_16, -10000.e0_16, 10000.e0_16, j, lat2_aux)
    !allocate(lat2f(j,3))
    !lat2f = lat2_aux(1:j,:)
    !deallocate(lat2_aux)
    !deallocate(lat2)

    ! OR DON'T TRIM
    j = size(lat2,1)
    allocate(lat2f(j,3))
    lat2f = lat2
    deallocate(lat2)
    ! -------------------------------------

    ! SEPARATE LATTICE INTO A AND B
    allocate(latA1(size(lat1f,1)/2, 3))
    allocate(latA2(size(lat1f,1)/2, 3))
    allocate(latB1(size(lat2f,1)/2, 3))
    allocate(latB2(size(lat2f,1)/2, 3))
    call separate_lattice(lat1f, latA1, latB1)
    call separate_lattice(lat2f, latA2, latB2)

    ! FIND OVERLAPPING POINTS FROM LATTICES 1 AND 2, SEPARATING A AND B
    !allocate(lat3_aux(n_points,3))
    !allocate(lat4_aux(n_points,3))
    j = 0
    k = 0
    max_ind = 2000
    min_ind = 0
    do i=1, size(latB2,1)
        condAB = exists_in_lattice([latB2(i,1), latB2(i,2), 0.e0_16], latA1(i-min_ind:i+max_ind,:))
        condBA = exists_in_lattice([latA2(i,1), latA2(i,2), 0.e0_16], latB1(i-min_ind:i+max_ind,:))
        condAA = exists_in_lattice([latA2(i,1), latA2(i,2), 0.e0_16], latA1(i-min_ind:i+max_ind,:))
        condBB = exists_in_lattice([latB2(i,1), latB2(i,2), 0.e0_16], latB1(i-min_ind:i+max_ind,:))
        
        if (condAB) then
            j = j + 1
            !lat3_aux(j,:) = latB2(i,:)
        endif

        if (condBA) then
            j = j + 1
            !lat3_aux(j,:) = latA2(i,:)
        endif
        
        if (condAA) then
            k = k + 1
            !lat4_aux(k,:) = latA2(i,:)
        endif
        
        if (condBB) then
            k = k + 1
            !lat4_aux(k,:) = latB2(i,:)
        endif

        if (i < 2000) then
            min_ind = min_ind + 1
        else
            if (i > n_points - 2000) then
                max_ind = max_ind - 1
            endif
        endif
    enddo
    !allocate(lat3f(j,3))
    !lat3f = lat3_aux(1:j,:)
    !allocate(lat4f(k,3))
    !lat4f = lat4_aux(1:k,:)
    !deallocate(lat3_aux)
    !deallocate(lat4_aux)

    write(*,*) "  Number of AB points: ", j
    write(*,*) "  Number of AA points: ", k

    ! WRITE LATTICES TO FILE
    !call write_lattice(lat1f, "lattice1t.dat")
    !call write_lattice(lat2f, "lattice2t.dat")
    !call write_lattice(lat3f, "lattice3t.dat")
    !call write_lattice(lat4f, "lattice4t.dat")

    ! DEALLOCATE
    deallocate(lat1f)
    deallocate(lat2f)
    !deallocate(lat3f)
    !deallocate(lat4f)
    deallocate(latA1)
    deallocate(latB1)
    deallocate(latA2)
    deallocate(latB2)
enddo

end program
