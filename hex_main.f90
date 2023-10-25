program hex_lattices
use hex_utils
implicit none

real*16, dimension(:,:), allocatable :: latA1, latA2, latB1, latB2
!real*16, dimension(:,:), allocatable :: latAA, latAB, lat_AB_aux, lat_AA_aux
real*16, dimension(3)                :: d1, origin1, origin2
real*16, dimension(:), allocatable   :: distsA, distsB
real*16                              :: angle, z, a, d, diff
integer                              :: i, j, k, l, m, n, i1, j1
logical                              :: condAB, condBA, condAA, condBB
logical                              :: AB_stacking, hex_center_pivot

! INITIAL DEFINITIONS
n = 50000
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
    origin2 = [0.e0_16, 0.e0_16, z]
endif

! ALLOCATION OF LATTICES AND FIRST CREATION
allocate(latA1(n/2,3))
allocate(latB1(n/2,3))
allocate(latA2(n/2,3))
allocate(latB2(n/2,3))

call create_honeycomb_lattice(latA1, latB1, 0.e0_16, a, .false.)
call create_honeycomb_lattice(latA2, latB2,       z, a, AB_stacking)

allocate(distsA(n/2))
allocate(distsB(n/2))
call distance_from_origin(latA1, origin1, distsA)
call distance_from_origin(latB1, origin1, distsB)

do i1=1, n/2
    do j1=1, n/2
        diff = abs(distsA(i1)-distsB(j1))
        if (diff .lt. 5.e-3_16) then
            write(*,*) diff, i1, j1
        endif
    enddo
enddo

!call sym_operation(11, latA2, origin2)
!call sym_operation(11, latB2, origin2)

!! TESTING ANGLES
!angle = pi/3.e0_16
!write(*,*) "Angle in radians: ", angle 
!
!! ! ROTATE LATTICE 2
!! call rotate_lattice_2(latA2, angle, origin2, [0.e0_16, 0.e0_16, 1.e0_16])
!! call rotate_lattice_2(latB2, angle, origin2, [0.e0_16, 0.e0_16, 1.e0_16])
!
!call rotate_lattice(latA2, angle, origin2, 3)
!call rotate_lattice(latB2, angle, origin2, 3)

!! FIND OVERLAPPING POINTS FROM LATTICES 1 AND 2, SEPARATING A AND B
!allocate(lat_AA_aux(n,3))
!allocate(lat_AB_aux(n,3))
!j = 0
!k = 0
!l = 0
!m = 0
!do i=1, n/2
!    condAB = exists_in_lattice([latB2(i,1),latB2(i,2),0.e0_16], latA1)
!    condBA = exists_in_lattice([latA2(i,1),latA2(i,2),0.e0_16], latB1)
!    condAA = exists_in_lattice([latA2(i,1),latA2(i,2),0.e0_16], latA1)
!    condBB = exists_in_lattice([latB2(i,1),latB2(i,2),0.e0_16], latB1)
!    if (condAB) then
!        j = j + 1
!        lat_AB_aux(j+k,:) = latB2(i,:)
!    endif
!    if (condBA) then
!        k = k + 1
!        lat_AB_aux(k+j,:) = latA2(i,:)
!    endif
!    if (condAA) then
!        l = l + 1
!        lat_AA_aux(l+m,:) = latA2(i,:)
!    endif
!    if (condBB) then
!        m = m + 1
!        lat_AA_aux(m+l,:) = latB2(i,:)
!    endif
!enddo
!
!allocate(latAA(m+l,3))
!allocate(latAB(j+k,3))
!latAA = lat_AA_aux(1:l+m,:)
!latAB = lat_AB_aux(1:j+k,:)
!deallocate(lat_AB_aux, lat_AA_aux)
!
!write(*,*) "  Number of AB points: ", j
!write(*,*) "  Number of BA points: ", k
!write(*,*) "  Number of AA points: ", l
!write(*,*) "  Number of BB points: ", m

! WRITE LATTICES TO FILE
call write_lattice(latA1, "latticeA1.dat")
call write_lattice(latB1, "latticeB1.dat")
call write_lattice(latA2, "latticeA2.dat")
call write_lattice(latB2, "latticeB2.dat")
!call write_lattice(latAA, "latticeAA.dat")
!call write_lattice(latAB, "latticeAB.dat")

! DEALLOCATE
!deallocate(latA1, latB1, latA2, latB2, latAA, latAB)
deallocate(latA1, latB1, latA2, latB2)

end program
