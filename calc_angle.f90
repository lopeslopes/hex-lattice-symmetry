program hex_lattices
use hex_utils
implicit none

real(kind=8), dimension(:,:), allocatable   :: latA1, latA2, latB1, latB2
real(kind=8), dimension(:,:), allocatable   :: latA_ol, latB_ol, latA_ol_aux, latB_ol_aux
real(kind=8), dimension(:,:), allocatable   :: latAq1, latBq1, latauxq1a, latauxq1b
real(kind=8), dimension(3)                  :: d1, origin1, origin2
real(kind=8), dimension(:), allocatable     :: distsA, distsB
real(kind=8)                                :: angle, z, a, d, diff, min_diff, min_angle
integer                                     :: i, j, k, l, m, n, i1, j1, k1, l1
logical                                     :: AB_stacking, hex_center_pivot, q1a, q1b

! INITIAL DEFINITIONS
n = 100000
a = 2.46d0
z = 3.35d0
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
    d = sqrt((a**2)/(2.d0*(1.d0-cos(2.d0*pi/3.d0))))
    d1 = [d*cos(pi/6.d0), d*sin(pi/6.d0), z]
    origin1 = d1 - [0.d0, 0.d0, z]
    origin2 = d1
else
    write(*,*) "Pivot point: node at origin"
    origin1 = [0.d0, 0.d0, 0.d0]
    origin2 = [0.d0, 0.d0,    z]
endif

! ALLOCATION OF LATTICES AND FIRST CREATION
write(*,*) "Creating lattices..."
allocate(latA1(n/2,3))
allocate(latB1(n/2,3))
allocate(latA2(n/2,3))
allocate(latB2(n/2,3))

call create_honeycomb_lattice(latA1, latB1, 0.d0, a,     .false.)
call create_honeycomb_lattice(latA2, latB2,    z, a, AB_stacking)

! TRIM THE LATTICE TO SEARCH ONLY FIRST SEXTANT
k1 = 0
l1 = 0
allocate(latauxq1a(n/2,3))
allocate(latauxq1b(n/2,3))
do i1=1, n/2
    q1a = ((latA1(i1,1) .gt. 0.d0) .and. (latA1(i1,2) .gt. 0.d0) .and. (latA1(i1,2) .lt. sqrt(3.d0)*latA1(i1,1)))
    q1b = ((latB1(i1,1) .gt. 0.d0) .and. (latB1(i1,2) .gt. 0.d0) .and. (latB1(i1,2) .lt. sqrt(3.d0)*latB1(i1,1)))
    if (q1a) then 
        k1 = k1 + 1
        latauxq1a(k1,:) = latA1(i1,:)
    endif
    if (q1b) then
        l1 = l1 + 1
        latauxq1b(l1,:) = latB1(i1,:)
    endif
enddo
allocate(latAq1(k1,3))
allocate(latBq1(l1,3))
latAq1 = latauxq1a(1:k1,:)
latBq1 = latauxq1b(1:l1,:)
deallocate(latauxq1a, latauxq1b)

! CALCULATE DISTANCES FROM ORIGIN
write(*,*) "Calculating distances..."
allocate(distsA(k1))
allocate(distsB(l1))
call distance_from_origin(latAq1, origin1, distsA)
call distance_from_origin(latBq1, origin1, distsB)

! SEARCH FOR POINTS THAT ARE CLOSE TO EACH OTHER IN DISTANCE
write(*,*) "Searching for similar distances..."
write(*,*) "tolerance: ", tol
open(25, file="angles.dat")
m = 0
min_diff = 5.d0
allocate(latA_ol_aux(5*n,3))
allocate(latB_ol_aux(5*n,3))
do i1=1, k1
    do j1=1, l1
        diff = abs(distsA(i1)-distsB(j1))
        if (diff .lt. tol) then
            angle = find_angle(latAq1(i1,:), latBq1(j1,:))
            write(25,*) angle
            m = m + 1
            if (diff .lt. min_diff) then
                min_diff = diff
                min_angle = angle
            endif
            latA_ol_aux(m,:) = latAq1(i1,:)
            latB_ol_aux(m,:) = latBq1(j1,:)
        endif
    enddo
enddo
close(25)
allocate(latA_ol(m,3))
allocate(latB_ol(m,3))
latA_ol = latA_ol_aux(1:m,:)
latB_ol = latB_ol_aux(1:m,:)

call write_lattice(latA_ol, "latticeOA.dat")
call write_lattice(latB_ol, "latticeOB.dat")

deallocate(latA_ol_aux, latB_ol_aux, latAq1, latBq1)
deallocate(latA_ol, latB_ol)
write(*,*) "Minimum distance: ", min_diff

end program
