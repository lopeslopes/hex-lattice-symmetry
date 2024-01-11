program calc_angle
use hex_utils
implicit none

real(kind=16), dimension(:,:), allocatable   :: latA1, latA2, latB1, latB2
real(kind=16), dimension(:,:), allocatable   :: latA_ol, latB_ol, latA_ol_aux, latB_ol_aux
real(kind=16), dimension(:,:), allocatable   :: latAq1, latBq1, latauxq1a, latauxq1b
real(kind=16), dimension(2)                  :: d1, origin1, origin2
real(kind=16), dimension(:), allocatable     :: distsA, distsB
real(kind=16)                                :: angle, z, a, d, diff, min_diff, min_angle
integer                                      :: i, j, k, l, m, n, i1, j1, k1, l1
logical                                      :: AB_stacking, hex_center_pivot, q1a, q1b

! INITIAL DEFINITIONS
n = 8000000
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
    d1 = [d*cos(pi/6.e0_16), d*sin(pi/6.e0_16)]
    origin1 = d1 - [0.e0_16, 0.e0_16]
    origin2 = d1
else
    write(*,*) "Pivot point: node at origin"
    origin1 = [0.e0_16, 0.e0_16]
    origin2 = [0.e0_16, 0.e0_16]
endif

! ALLOCATION OF LATTICES AND FIRST CREATION
write(*,*) "Creating lattices..."
allocate(latA1(n/2,2))
allocate(latB1(n/2,2))
allocate(latA2(n/2,2))
allocate(latB2(n/2,2))

call create_honeycomb_lattice(latA1, latB1, a,     .false.)
call create_honeycomb_lattice(latA2, latB2, a, AB_stacking)

! TRIM THE LATTICE TO SEARCH ONLY FIRST SEXTANT
k1 = 0
l1 = 0
allocate(latauxq1a(n/2,2))
allocate(latauxq1b(n/2,2))
do i1=1, n/2
    q1a = ((latA1(i1,1) .gt. 0.e0_16) .and. (latA1(i1,2) .gt. 0.e0_16) .and. (latA1(i1,2) .lt. sqrt(3.e0_16)*latA1(i1,1)))
    q1b = ((latB1(i1,1) .gt. 0.e0_16) .and. (latB1(i1,2) .gt. 0.e0_16) .and. (latB1(i1,2) .lt. sqrt(3.e0_16)*latB1(i1,1)))
    if (q1a) then 
        k1 = k1 + 1
        latauxq1a(k1,:) = latA1(i1,:)
    endif
    if (q1b) then
        l1 = l1 + 1
        latauxq1b(l1,:) = latB1(i1,:)
    endif
enddo
allocate(latAq1(k1,2))
allocate(latBq1(l1,2))
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
min_diff = 5.e0_16
allocate(latA_ol_aux(5*n,2))
allocate(latB_ol_aux(5*n,2))
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
allocate(latA_ol(m,2))
allocate(latB_ol(m,2))
latA_ol = latA_ol_aux(1:m,:)
latB_ol = latB_ol_aux(1:m,:)

call write_lattice(latA_ol, "latticeOA.dat")
call write_lattice(latB_ol, "latticeOB.dat")

deallocate(latA_ol_aux, latB_ol_aux, latAq1, latBq1)
deallocate(latA_ol, latB_ol)
write(*,*) "Minimum distance: ", min_diff

end program
