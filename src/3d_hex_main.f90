program hex_lattices
use hex_utils
implicit none

real(kind=16), dimension(:,:), allocatable :: latA1, latA2, latB1, latB2
real(kind=16), dimension(:,:), allocatable :: latAA, latAB, lat_AA_aux, lat_AB_aux
real(kind=16), dimension(:,:), allocatable :: latBB, latBA, lat_BB_aux, lat_BA_aux
real(kind=16), dimension(3)                :: d1, origin1, origin2
real(kind=16)                              :: z, a, d, angle
integer                                    :: i, j, k, l, m, n
logical, dimension(:), allocatable         :: condAB, condBA, condAA, condBB
logical                                    :: AB_stacking, hex_center_pivot

! INITIAL DEFINITIONS
n = 1000000
a = 2.46e0_16
z = 3.35e0_16
hex_center_pivot = .false.
AB_stacking = .false.

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
    origin2 = [0.e0_16, 0.e0_16,    z]
endif

! ALLOCATION OF LATTICES AND FIRST CREATION
write(*,*) "Creating lattices..."
allocate(latA1(n/2,3))
allocate(latB1(n/2,3))
allocate(latA2(n/2,3))
allocate(latB2(n/2,3))

call create_honeycomb_lattice(latA1, latB1, 0.e0_16, a,     .false.)
call create_honeycomb_lattice(latA2, latB2,    z, a, AB_stacking)

angle = 1.91601993415282663305148436041635378e-2_16!0.019166696333079717e0_16
write(*,*) "Angle in radians: ", angle
write(*,*) "Angle in degrees: ", (angle*180.e0_16)/pi

! ROTATE SECOND LATTICE BY THE ANGLE
call rotate_lattice(latA2, angle, origin2, 3)
call rotate_lattice(latB2, angle, origin2, 3)
call write_lattice(latA1, "latticeA1.dat")
call write_lattice(latB1, "latticeB1.dat")
call write_lattice(latA2, "latticeA2.dat")
call write_lattice(latB2, "latticeB2.dat")

! CHECK IF THERE'S ANY POINTS OVERLAPING (CHECK TOLERANCE)
write(*,*) "Checking overlapping points after rotation..."
allocate(lat_AA_aux(n,3))
allocate(lat_AB_aux(n,3))
allocate(lat_BB_aux(n,3))
allocate(lat_BA_aux(n,3))
j = 0
k = 0
l = 0
m = 0

! CHECKING EXISTENCE IN LATTICE USING OPENMP
allocate(condAA(n/2))
allocate(condAB(n/2))
allocate(condBA(n/2))
allocate(condBB(n/2))
condAA = .false.
condAB = .false.
condBA = .false.
condBB = .false.

!$omp parallel do private(i) shared(latA1, latB1, latA2, latB2) num_threads(6)
do i=1, n/2
    condAB(i) = any(all(abs(latA1 - spread([latB2(i,1),latB2(i,2),0.e0_16], dim=1, ncopies=n/2)) .lt. tol2, dim=2))
    condBA(i) = any(all(abs(latB1 - spread([latA2(i,1),latA2(i,2),0.e0_16], dim=1, ncopies=n/2)) .lt. tol2, dim=2))
    condAA(i) = any(all(abs(latA1 - spread([latA2(i,1),latA2(i,2),0.e0_16], dim=1, ncopies=n/2)) .lt. tol2, dim=2))
    condBB(i) = any(all(abs(latB1 - spread([latB2(i,1),latB2(i,2),0.e0_16], dim=1, ncopies=n/2)) .lt. tol2, dim=2))
enddo
!$omp end parallel do

do i=1, n/2
    if (condAB(i)) then
        j = j + 1
        lat_AB_aux(j,:) = latB2(i,:)
    endif
    if (condBA(i)) then
        k = k + 1
        lat_BA_aux(k,:) = latA2(i,:)
    endif
    if (condAA(i)) then
        l = l + 1
        lat_AA_aux(l,:) = latA2(i,:)
    endif
    if (condBB(i)) then
        m = m + 1
        lat_BB_aux(m,:) = latB2(i,:)
    endif
enddo

allocate(latAB(j,3))
allocate(latBA(k,3))
allocate(latAA(l,3))
allocate(latBB(m,3))
latAA = lat_AA_aux(1:l,:)
latAB = lat_AB_aux(1:j,:)
latBA = lat_BA_aux(1:k,:)
latBB = lat_BB_aux(1:m,:)
deallocate(lat_AB_aux, lat_AA_aux, lat_BA_aux, lat_BB_aux)

! WRITE LATTICES TO FILE
call write_lattice(latAA, "latticeAA.dat")
call write_lattice(latAB, "latticeAB.dat")
call write_lattice(latBA, "latticeBA.dat")
call write_lattice(latBB, "latticeBB.dat")

! DEALLOCATE
deallocate(latA1, latB1, latA2, latB2)
deallocate(latAA, latAB, latBA, latBB)

write(*,*) "  Number of AB points: ", j
write(*,*) "  Number of BA points: ", k
write(*,*) "  Number of AA points: ", l
write(*,*) "  Number of BB points: ", m

end program
