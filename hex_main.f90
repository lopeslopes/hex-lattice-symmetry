program hex_lattices
use hex_utils
implicit none

real*16, dimension(:,:), allocatable :: lat1, lat2, lat1f, lat2f, lat3f
real*16, dimension(:,:), allocatable :: lat1_aux, lat2_aux, lat3_aux
real*16, dimension(:,:), allocatable :: latAA, latAB
real*16                              :: angle, z, a
integer                              :: i, j, n_points, max_ind, min_ind
logical                              :: cond3

n_points = 100000
angle = magic_angle(30)!(1.08455e0_16*pi)/180.e0_16!magic_angle(29)
write(*,*) "angle in radians: ", angle
write(*,*) "angle in degrees: ", (angle*180.e0_16)/pi
z = 3.35e0_16
a = 2.46e0_16

! ALLOCATE AND CREATE LATTICES
allocate(lat1(n_points,3))
allocate(lat2(n_points,3))
call create_lattice_eh(lat1, 0.e0_16, a, .false.)
call create_lattice_eh(lat2, z      , a, .false.)

! ROTATE LATTICE 2
call rotate_lattice(lat2, angle)

! TRIM THE LATTICE 1
allocate(lat1_aux(n_points, 3))
j = 0
call trim_lattice(lat1, -10000.e0_16, 10000.e0_16, -10000.e0_16, 10000.e0_16, j, lat1_aux)
allocate(lat1f(j,3))
lat1f = lat1_aux(1:j,:)
deallocate(lat1_aux)
deallocate(lat1)

! TRIM THE LATTICE 2
allocate(lat2_aux(n_points, 3))
j = 0
call trim_lattice(lat2, -10000.e0_16, 10000.e0_16, -10000.e0_16, 10000.e0_16, j, lat2_aux)
allocate(lat2f(j,3))
lat2f = lat2_aux(1:j,:)
deallocate(lat2_aux)
deallocate(lat2)

! FIND POINTS FROM LAT2 WHERE THERE'S A LAT1 POINT RIGHT ABOVE/BELOW
! TODO: SEPARATE BETWEEN A AND B LATTICES TO CHECK WHICH POINTS ARE OVERLAPPING
allocate(lat3_aux(n_points,3))
j = 0
max_ind = 10000
min_ind = 0
do i=1, size(lat2f,1)
    cond3 = exists_in_lattice([lat2f(i,1), lat2f(i,2), 0.e0_16], lat1f(i-min_ind:i+max_ind,:))
    if (cond3) then
        j = j + 1
        lat3_aux(j,:) = lat2f(i,:)
    endif
    
    if (i < 10000) then
        min_ind = min_ind + 1
    else
        if (i > n_points - 10000) then
            max_ind = max_ind - 1
        endif
    endif
enddo
allocate(lat3f(j,3))
lat3f = lat3_aux(1:j,:)

! WRITE LATTICES TO FILE
call write_lattice(lat1f, "lattice1t.dat")
call write_lattice(lat2f, "lattice2t.dat")
call write_lattice(lat3f, "lattice3t.dat")

! DEALLOCATE
deallocate(lat1f)
deallocate(lat2f)
deallocate(lat3f)

end program
