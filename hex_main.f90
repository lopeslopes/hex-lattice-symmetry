program hex_lattices
use hex_utils
implicit none

real*16, dimension(:,:), allocatable :: lat1, lat2, lat1f, lat2f
real*16, dimension(:,:), allocatable :: lat1_aux, lat2_aux
real*16                              :: angle, z, a
integer                              :: i, j, n_points

n_points = 100
angle = (10.e0_16*pi)/180.e0_16
z = 2.e0_16
a = 2.42e0_16

! ALLOCATE AND CREATE LATTICES
allocate(lat1(n_points,3))
allocate(lat2(n_points,3))
call create_lattice_eh(lat1, 0.e0_16, a)
call create_lattice_eh(lat2, z      , a)

! ROTATE LATTICE 2
call rotate_lattice(lat2, angle) 

! TRIM THE LATTICE 1
allocate(lat1_aux(n_points, 3))
j = 0
call trim_lattice(lat1, -5.e0_16, 5.e0_16, -2.e0_16, 2.e0_16, j, lat1_aux)
allocate(lat1f(j,3))
lat1f = lat1_aux(1:j,:)
deallocate(lat1_aux)
deallocate(lat1)

! TRIM THE LATTICE 2
allocate(lat2_aux(n_points, 3))
j = 0
call trim_lattice(lat2, -5.e0_16, 5.e0_16, -2.e0_16, 2.e0_16, j, lat2_aux)
allocate(lat2f(j,3))
lat2f = lat2_aux(1:j,:)
deallocate(lat2_aux)
deallocate(lat2)

! FIND POINTS FROM LAT2 WHERE THERE'S A LAT1 POINT RIGHT ABOVE/BELOW
! PART OF OLD VERSION
!allocate(lat3_aux(num_points,3))
!lat3_aux = 0.e0_16
!l = 0
!do i=1, min(size(lat1f,1), size(lat2f,1))
!    cond3 = exists_in_lattice([lat2f(i,1), lat2f(i,2), 0.e0_16], lat1f(i-10000:i+10000,:))
!    if (cond3) then
!        l = l + 1
!        lat3_aux(l,:) = lat2f(i,:)
!    endif
!enddo


! WRITE LATTICES TO FILE
call write_lattice(lat1f, "lattice1t.dat")
call write_lattice(lat2f, "lattice2t.dat")

! DEALLOCATE
deallocate(lat1f)
deallocate(lat2f)

end program
