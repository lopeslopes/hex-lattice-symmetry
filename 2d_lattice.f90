program symmetry
implicit none

real*16, dimension(:,:), allocatable :: lat1A, lat1B, lat2A, lat2B, lat3
real*16, dimension(:,:), allocatable :: lat1Af, lat1Bf, lat2Af, lat2Bf
real*16, dimension(:,:), allocatable :: lat1A_aux, lat1B_aux, lat2A_aux, lat2B_aux, lat3_aux
real*16                              :: rot_angle, lat_separation, lat_const
real*16                              :: pi
integer                              :: i, j, k, l, m, num_points
logical                              :: cond1A, cond2A, cond1B, cond2B
logical                              :: cond3AA, cond3AB, cond3BA, cond3BB

! INITIAL DEFINITIONS
pi = 4.e0_16*atan(1.e0_16)

num_points      = 4000000
lat_separation  = 10.e0_16
lat_const       = 2.42d0
rot_angle       = (1.08445d0/180.e0_16) * pi

! ALLOCATE THE TWO LATTICES
allocate(lat1A(num_points/2, 3))
allocate(lat1B(num_points/2, 3))
allocate(lat2A(num_points/2, 3))
allocate(lat2B(num_points/2, 3))
lat1A = 0.e0_16
lat1B = 0.e0_16
lat2A = 0.e0_16
lat2B = 0.e0_16

! GENERATE THE LATTICE POINTS
call create_lattice_eh_AB(lat1A, lat1B, 0.e0_16       , lat_const, 1000)
call create_lattice_eh_AB(lat2A, lat2B, lat_separation, lat_const, 1000)

! ROTATE SECOND LATTICE
do i=1, num_points/2
    lat2A(i,:) = rotate_point(lat2A(i,:), [0.e0_16, 0.e0_16, lat_separation], rot_angle)
    lat2B(i,:) = rotate_point(lat2B(i,:), [0.e0_16, 0.e0_16, lat_separation], rot_angle)
enddo

! TRIM THE LATTICES (ONLY USE WHAT HAS -50 < x < 50 AND y > 0) FOR PERFORMANCE
allocate(lat1A_aux(num_points/2,3))
allocate(lat2A_aux(num_points/2,3))
allocate(lat1B_aux(num_points/2,3))
allocate(lat2B_aux(num_points/2,3))
lat1A_aux = 0.e0_16
lat2A_aux = 0.e0_16
lat1B_aux = 0.e0_16
lat2B_aux = 0.e0_16

j = 0
k = 0
l = 0
m = 0
do i=1, num_points/2
    cond1A = (lat1A(i,1) > -50.e0_16) .and. (lat1A(i,1) < 50.e0_16) .and. (lat1A(i,2) >= 0.e0_16)
    cond1B = (lat1B(i,1) > -50.e0_16) .and. (lat1B(i,1) < 50.e0_16) .and. (lat1B(i,2) >= 0.e0_16)
    cond2A = (lat2A(i,1) > -50.e0_16) .and. (lat2A(i,1) < 50.e0_16) .and. (lat2A(i,2) >= 0.e0_16)
    cond2B = (lat2B(i,1) > -50.e0_16) .and. (lat2B(i,1) < 50.e0_16) .and. (lat2B(i,2) >= 0.e0_16)    
    if (cond1A) then
        j = j + 1
        lat1A_aux(j,:) = lat1A(i,:)
    endif
    if (cond2A) then
        k = k + 1
        lat2A_aux(k,:) = lat2A(i,:)
    endif
    if (cond1B) then
        l = l + 1
        lat1B_aux(j,:) = lat1B(i,:)
    endif
    if (cond2B) then
        m = m + 1
        lat2B_aux(k,:) = lat2B(i,:)
    endif
enddo

! STORE TRIMMED LATTICES IN SMALLER ARRAYS
allocate(lat1Af(j,3))
allocate(lat2Af(k,3))
allocate(lat1Bf(l,3))
allocate(lat2Bf(m,3))
lat1Af = lat1A_aux(1:j,:)
lat2Af = lat2A_aux(1:k,:)
lat1Bf = lat1B_aux(1:l,:)
lat2Bf = lat2B_aux(1:m,:)


! WRITE THEM TO FILES
call write_lattice(lat1Af, "lattice1A.dat")
call write_lattice(lat2Af, "lattice2A.dat")
call write_lattice(lat1Bf, "lattice1B.dat")
call write_lattice(lat2Bf, "lattice2B.dat")

! GENERATE THIRD LATTICE TO STORE WHICH POINTS OF THE FIRST ONE
! ARE EXACTLY BELOW A POINT FROM THE SECOND ONE
allocate(lat3_aux(num_points,3))
lat3_aux = 0.e0_16
l = 0
do i=1, min(j,k)
    cond3AA = exists_in_lattice([lat2Af(i,1), lat2Af(i,2), 0.e0_16], lat1Af(i-2000:i+2000,:))
    cond3AB = exists_in_lattice([lat2Af(i,1), lat2Af(i,2), 0.e0_16], lat1Bf(i-2000:i+2000,:))
    cond3BA = exists_in_lattice([lat2Bf(i,1), lat2Bf(i,2), 0.e0_16], lat1Af(i-2000:i+2000,:))
    cond3BB = exists_in_lattice([lat2Bf(i,1), lat2Bf(i,2), 0.e0_16], lat1Bf(i-2000:i+2000,:))

    if (cond3AA .or. cond3AB) then
        l = l + 1
        lat3_aux(l,:) = lat2Af(i,:)
    endif
    if (cond3BA .or. cond3BB) then
        l = l + 1
        lat3_aux(l,:) = lat2Bf(i,:)
    endif
enddo

! ALLOCATE ARRAY OF THE RIGHT SIZE TO STORE THIRD LATTICE AND WRITE
allocate(lat3(l,3))
lat3 = lat3_aux(1:l,:)
call write_lattice(lat3, "lattice03.dat")

! DEALLOCATE EVERYTHING
deallocate(lat1A)
deallocate(lat2A)
deallocate(lat1B)
deallocate(lat2B)
deallocate(lat3)
deallocate(lat1Af)
deallocate(lat2Af)
deallocate(lat1Bf)
deallocate(lat2Bf)
deallocate(lat1A_aux)
deallocate(lat2A_aux)
!deallocate(lat1B_aux)
!deallocate(lat2B_aux)
!deallocate(lat3_aux)

contains

function exists_in_lattice(point, lattice) result(flag)
    real*16, dimension(3), intent(in)     :: point
    real*16, dimension(:,:), intent(in)   :: lattice
    real*16                               :: c1, c2, c3, tol
    logical, dimension(:), allocatable   :: condition_list
    logical                              :: flag
    integer                              :: i, cont

    allocate(condition_list(size(lattice,1)))

    tol = 1.d-3
    do i=1, size(lattice, 1)
        c1 = abs(lattice(i,1) - point(1))
        c2 = abs(lattice(i,2) - point(2))
        c3 = abs(lattice(i,3) - point(3))
        condition_list(i) = (c1 < tol) .and. (c2 < tol) .and. (c3 < tol)
    enddo

    flag = any(condition_list .eqv. .true.)
end function

! CREATES LATTICE OF BODY-CENTERED HEXAGONS
subroutine create_lattice_bch(lattice, z, a, num_columns)
    real*16, dimension(:,:), intent(inout) :: lattice
    real*16, intent(in)                    :: z, a
    real*16, dimension(3)                  :: v1, v2, origin_gen, lat_origin
    real*16                                :: angle, pi, size_x, size_y
    integer                               :: i, j, num_columns, row, i0

    pi = 4.e0_16 * atan(1.e0_16)
    angle = (120.e0_16/180.e0_16) * pi
    v1 = [a, 0.e0_16, 0.e0_16]
    v2 = [a*cos(angle), a*sin(angle), 0.e0_16]

    origin_gen = [0.e0_16, 0.e0_16, z]

    row = 1
    lattice(1,:) = origin_gen
    i = 2
    do while(i < size(lattice,1))
        do j=1, num_columns
            lattice(i,:) = origin_gen + real(j)*v1
            i = i + 1
        enddo
        row  = row + 1
        if (modulo(row,2) == 1) then
            origin_gen = origin_gen + v1
        endif
        origin_gen = origin_gen + v2
    enddo
    lattice(1,:) = [0.e0_16, 0.e0_16, z]

    i0 = ((row/2) * num_columns) + (num_columns/2)
    lat_origin = lattice(i0,:)
    lattice(:,1) = lattice(:,1) - lat_origin(1)
    lattice(:,2) = lattice(:,2) - lat_origin(2)
    lattice(:,3) = lattice(:,3) - lat_origin(3)
end subroutine

! CREATES LATTICE OF EMPTY HEXAGONS (HONEYCOMB)
subroutine create_lattice_eh(lattice, z, a, num_columns)
    real*16, dimension(:,:), intent(inout) :: lattice
    real*16, intent(in)                    :: z, a
    real*16, dimension(3)                  :: v1, v2, v3, v4, origin_a, origin_b, lat_origin
    real*16                                :: angle, pi, d1
    integer                                :: i, j, num_columns, row, i0

    pi = 4.e0_16 * atan(1.e0_16)
    angle = (30.e0_16/180.e0_16) * pi
    d1 = a/sqrt(3.e0_16)

    v1 = [1.5d0*d1, 0.e0_16, 0.e0_16]
    v2 = [a*cos(angle), a*sin(angle), 0.e0_16]
    v3 = [0.e0_16, a, 0.e0_16]

    origin_a = [0.e0_16, 0.e0_16, z]
    origin_b = origin_a + [d1, 0.e0_16, 0.e0_16]

    row = 1
    lattice(1,:) = origin_a
    lattice(2,:) = origin_b
    i = 3
    do while(i < size(lattice,1))
        do j=1, num_columns, 2
            lattice(i,:) = origin_a + real(j)*v1
            i = i + 1
            lattice(i,:) = origin_b + real(j)*v1
            i = i + 1
        enddo
        row  = row + 1
        if (modulo(row,2) == 1) then
            origin_a = origin_a + v2
            origin_b = origin_b + v2
        else
            origin_a = origin_a - v2 + v3
            origin_b = origin_b - v2 + v3
        endif
    enddo
    lattice(1,:) = [0.e0_16, 0.e0_16, z]

    i0 = ((row/2) * num_columns) + (num_columns/2)
    lat_origin = lattice(i0,:)
    lattice(:,1) = lattice(:,1) - lat_origin(1)
    lattice(:,2) = lattice(:,2) - lat_origin(2)
    lattice(:,3) = lattice(:,3) - lat_origin(3)
end subroutine

subroutine create_lattice_eh_AB(latticeA, latticeB, z, a, num_columns)
    real*16, dimension(:,:), intent(inout) :: latticeA, latticeB
    real*16, intent(in)                    :: z, a
    real*16, dimension(3)                  :: v1, v2, v3, d1, origin_a, origin_b, lat_origin
    real*16                                :: angle, pi, d
    integer                                :: i, j, num_columns, row, i0

    pi = 4.e0_16 * atan(1.e0_16)
    angle = (60.e0_16/180.e0_16) * pi

    v1 = [a, 0.e0_16, 0.e0_16]
    v2 = [a*cos(angle), a*sin(angle), 0.e0_16]
    v3 = [-a*cos(angle), a*sin(angle), 0.e0_16]

    d = sqrt((a**2)/(2.e0_16*(1.e0_16-cos(2.e0_16*angle))))
    d1 = [d*cos(angle/2.e0_16), d*sin(angle/2.e0_16), 0.e0_16]

    origin_a = [0.e0_16, 0.e0_16, z]
    origin_b = origin_a + d1

    row = 1
    latticeA(1,:) = origin_a
    latticeB(1,:) = origin_b
    i = 2
    do while(i < size(latticeA,1))
        do j=1, num_columns
            latticeA(i,:) = origin_a + real(j,16)*v1
            latticeB(i,:) = origin_b + real(j,16)*v1
            i = i + 1
        enddo
        row  = row + 1
        if (modulo(row,2) == 1) then
            origin_a = origin_a + v2
            origin_b = origin_b + v2
        else
            origin_a = origin_a + v3
            origin_b = origin_b + v3
        endif
    enddo
    latticeA(1,:) = [0.e0_16, 0.e0_16, z]
    latticeB(1,:) = [d*cos(angle/2.e0_16), d*sin(angle/2.e0_16), z]

    i0 = ((row/2) * num_columns) + (num_columns/2)
    lat_origin = latticeA(i0,:)
    latticeA(:,1) = latticeA(:,1) - lat_origin(1)
    latticeA(:,2) = latticeA(:,2) - lat_origin(2)
    latticeA(:,3) = latticeA(:,3) - lat_origin(3)
    latticeB(:,1) = latticeB(:,1) - lat_origin(1)
    latticeB(:,2) = latticeB(:,2) - lat_origin(2)
    latticeB(:,3) = latticeB(:,3) - lat_origin(3)
end subroutine

function rotate_point(point, pivot, angle) result(new_point)
    real*16, dimension(3), intent(in) :: point
    real*16, intent(in)               :: angle
    real*16, dimension(3)             :: aux1, aux2, pivot, new_point

    aux1 = point - pivot
    
    aux2(1) = aux1(1)*cos(angle) - aux1(2)*sin(angle)
    aux2(2) = aux1(1)*sin(angle) + aux1(2)*cos(angle)
    aux2(3) = aux1(3)

    new_point = aux2 + pivot
end function

subroutine write_lattice(lattice, filename)
    real*16, dimension(:,:), intent(in) :: lattice
    character(13), intent(in)          :: filename
    integer                            :: i

    open(23, file=filename)
    do i=1, size(lattice, 1)
        write(23,*) lattice(i,1), ";", lattice(i,2), ";", lattice(i,3)
    enddo
    close(23)
end subroutine

end program symmetry
