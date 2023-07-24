program symmetry
implicit none

real*16, dimension(:,:), allocatable :: lat1, lat1f, lat1_aux, lat2, lat2f, lat2_aux, lat3, lat3_aux
real*16                              :: rot_angle, lat_separation, lat_const
real*16                              :: pi
integer                              :: i, j, k, l, num_points
logical                              :: cond1A, cond2A, cond1B, cond2B, cond3

! INITIAL DEFINITIONS
pi = 4.e0_16*atan(1.e0_16)

num_points      = 4000000
lat_separation  = 10.e0_16
lat_const       = 2.42e0_16
rot_angle       = (1.08445e0_16/180.e0_16) * pi

! ALLOCATE THE TWO LATTICES
allocate(lat1(num_points, 3))
allocate(lat2(num_points, 3))
lat1 = 0.e0_16
lat2 = 0.e0_16

! GENERATE THE LATTICE POINTS
call create_lattice_eh(lat1, 0.e0_16       , lat_const, 2000)
call create_lattice_eh(lat2, lat_separation, lat_const, 2000)

! ROTATE SECOND LATTICE
do i=1, num_points
    lat2(i,:) = rotate_point(lat2(i,:), [0.e0_16, 0.e0_16, lat_separation], rot_angle)
enddo

!! TRIM THE LATTICES (ONLY USE WHAT HAS -50 < x < 50 AND y > 0) FOR PERFORMANCE
allocate(lat1_aux(num_points,3))
allocate(lat2_aux(num_points,3))
lat1_aux = 0.e0_16
lat2_aux = 0.e0_16

j = 0
k = 0
do i=1, num_points
    cond1A = (lat1(i,1) > -100.e0_16) .and. (lat1(i,1) < 100.e0_16)
    cond1B = (lat1(i,2) <= 1000.e0_16) .and. (lat1(i,2) >= 0.e0_16)
    cond2A = (lat2(i,1) > -100.e0_16) .and. (lat2(i,1) < 100.e0_16)
    cond2B = (lat2(i,2) <= 1000.e0_16) .and. (lat2(i,2) >= 0.e0_16)

    if (cond1A .and. cond1B) then
        j = j + 1
        lat1_aux(j,:) = lat1(i,:)
    endif
    if (cond2A .and. cond2B) then
        k = k + 1
        lat2_aux(k,:) = lat2(i,:)
    endif
enddo

! STORE TRIMMED LATTICES IN SMALLER ARRAYS
allocate(lat1f(j,3))
allocate(lat2f(k,3))
lat1f = lat1_aux(1:j,:)
lat2f = lat2_aux(1:k,:)

!! SECTION TO USE IF NOT TRIMMING THE LATTICES
!allocate(lat1f(size(lat1,1), size(lat1,2)))
!allocate(lat2f(size(lat2,1), size(lat2,2)))
!lat1f = lat1
!lat2f = lat2

! WRITE THEM TO FILES
call write_lattice(lat1f, "lattice01.dat")
call write_lattice(lat2f, "lattice02.dat")

! GENERATE THIRD LATTICE TO STORE WHICH POINTS OF THE FIRST ONE
! ARE EXACTLY BELOW A POINT FROM THE SECOND ONE
allocate(lat3_aux(num_points,3))
lat3_aux = 0.e0_16
l = 0
do i=1, min(size(lat1f,1), size(lat2f,1))
    cond3 = exists_in_lattice([lat2f(i,1), lat2f(i,2), 0.e0_16], lat1f(i-10000:i+10000,:))
    if (cond3) then
        l = l + 1
        lat3_aux(l,:) = lat2f(i,:)
    endif
enddo

! ALLOCATE ARRAY OF THE RIGHT SIZE TO STORE THIRD LATTICE AND WRITE
allocate(lat3(l,3))
lat3 = lat3_aux(1:l,:)
call write_lattice(lat3, "lattice03.dat")

! DEALLOCATE EVERYTHING
deallocate(lat1)
deallocate(lat2)
deallocate(lat3)
deallocate(lat1f)
deallocate(lat2f)
!deallocate(lat1_aux)
!deallocate(lat2_aux)
!deallocate(lat3_aux)

contains

function exists_in_lattice(point, lattice) result(flag)
    real*16, dimension(3), intent(in)     :: point
    real*16, dimension(:,:), intent(in)   :: lattice
    real*16                               :: tol
    logical                               :: flag
    real*16, dimension(3)                 :: vec1
    logical, dimension(:), allocatable    :: condition_list
    integer                               :: i

    tol = 1.d-3
    !flag = any(all(abs(spread(point, dim=1, ncopies=size(lattice,1)) - lattice) .lt. tol, dim=2))

    allocate(condition_list(size(lattice,1)))
    do i=1, size(lattice, 1)
        vec1 = abs(lattice(i,:) - point)
        condition_list(i) = all(vec1 .lt. tol)
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
    lattice(1,:) = origin_a
    lattice(2,:) = origin_b
    i = 3
    do while(i < size(lattice,1))
        do j=1, num_columns
            lattice(i,:) = origin_a + real(j,16)*v1
            i = i + 1
            lattice(i,:) = origin_b + real(j,16)*v1
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
    lattice(1,:) = [0.e0_16, 0.e0_16, z]
    lattice(2,:) = [d*cos(angle/2.e0_16), d*sin(angle/2.e0_16), z]

    i0 = ((row/2) * num_columns) + (num_columns/2)
    lat_origin = lattice(i0,:)
    lattice(:,1) = lattice(:,1) - lat_origin(1)
    lattice(:,2) = lattice(:,2) - lat_origin(2)
    lattice(:,3) = lattice(:,3) - lat_origin(3)
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
