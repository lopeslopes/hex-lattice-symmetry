program symmetry
implicit none

real*8, dimension(:,:), allocatable :: lat1, lat2, lat3, lat1f, lat2f
real*8, dimension(:,:), allocatable :: lat1_aux, lat2_aux, lat3_aux
real*8                              :: rot_angle, lat_separation, lat_const
real*8                              :: pi
integer                             :: i, j, k, l, num_points
logical                             :: cond1, cond2

pi = 4.d0*atan(1.d0)

num_points      = 4000000
lat_separation  = 10.d0
lat_const       = 2.42d0
rot_angle       = (1.08445d0/180.d0) * pi

allocate(lat1(num_points, 3))
allocate(lat2(num_points, 3))
lat1 = 0.d0
lat2 = 0.d0

call create_lattice_eh(lat1, 0.d0          , lat_const, 2000)
call create_lattice_eh(lat2, lat_separation, lat_const, 2000)

allocate(lat1_aux(num_points,3))
allocate(lat2_aux(num_points,3))
lat1_aux = 0.d0
lat2_aux = 0.d0
j = 0
k = 0
do i=1, num_points
    lat2(i,:) = rotate_point(lat2(i,:), [0.d0, 0.d0, lat_separation], rot_angle)
enddo

do i=1, num_points
    cond1 = (lat1(i,1) > -50.d0) .and. (lat1(i,1) < 50.d0) .and. (lat1(i,2) >= 0.d0)
    cond2 = (lat2(i,1) > -50.d0) .and. (lat2(i,1) < 50.d0) .and. (lat2(i,2) >= 0.d0)
    if (cond1) then
        j = j + 1
        lat1_aux(j,:) = lat1(i,:)
    endif
    if (cond2) then
        k = k + 1
        lat2_aux(k,:) = lat2(i,:)
    endif
enddo

allocate(lat1f(j,3))
allocate(lat2f(k,3))
lat1f = lat1_aux(1:j,:)
lat2f = lat2_aux(1:k,:)

call write_lattice(lat1f, "lattice01.dat")
call write_lattice(lat2f, "lattice02.dat")

allocate(lat3_aux(num_points,3))
lat3_aux = 0.d0
l = 0
do i=1, min(j,k)
    if (exists_in_lattice([lat2f(i,1), lat2f(i,2), 0.d0], lat1f(i-1000:i+1000,:))) then
        l = l + 1
        lat3_aux(l,:) = lat2f(i,:)
    endif
enddo

allocate(lat3(l,3))
lat3 = lat3_aux(1:l,:)
call write_lattice(lat3, "lattice03.dat")

deallocate(lat1)
deallocate(lat2)
deallocate(lat3)
deallocate(lat1f)
deallocate(lat2f)
deallocate(lat1_aux)
deallocate(lat2_aux)
deallocate(lat3_aux)

contains

function exists_in_lattice(point, lattice) result(flag)
    real*8, dimension(3), intent(in)     :: point
    real*8, dimension(:,:), intent(in)   :: lattice
    real*8                               :: c1, c2, c3, tol
    logical, dimension(:), allocatable   :: condition_list
    logical                              :: flag
    integer                              :: i, cont

    allocate(condition_list(size(lattice,1)))

    tol = 1e-8
    do i=1, size(lattice, 1)
        c1 = abs(lattice(i,1) - point(1))
        c2 = abs(lattice(i,2) - point(2))
        c3 = abs(lattice(i,3) - point(3))
        condition_list(i) = (c1 < tol) .and. (c2 < tol) .and. (c3 < tol)
    enddo

    cont = count(condition_list .eqv. .true.)
    flag = any(condition_list .eqv. .true.)
end function

subroutine create_lattice_bch(lattice, z, a, num_columns) ! body-centered hexagon
    real*8, dimension(:,:), intent(inout) :: lattice
    real*8, intent(in)                    :: z, a
    real*8, dimension(3)                  :: v1, v2, origin_gen, lat_origin
    real*8                                :: angle, pi, size_x, size_y
    integer                               :: i, j, num_columns, row, i0

    pi = 4.d0 * atan(1.d0)
    angle = (120.d0/180.d0) * pi
    v1 = [a, 0.d0, 0.d0]
    v2 = [a*cos(angle), a*sin(angle), 0.d0]

    origin_gen = [0.d0, 0.d0, z]

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
    lattice(1,:) = [0.d0, 0.d0, z]

    i0 = ((row/2) * num_columns) + (num_columns/2)
    lat_origin = lattice(i0,:)
    lattice(:,1) = lattice(:,1) - lat_origin(1)
    lattice(:,2) = lattice(:,2) - lat_origin(2)
    lattice(:,3) = lattice(:,3) - lat_origin(3)
end subroutine

subroutine create_lattice_eh(lattice, z, a, num_columns) ! empty hexagon
    real*8, dimension(:,:), intent(inout) :: lattice
    real*8, intent(in)                    :: z, a
    real*8, dimension(3)                  :: v1, v2, v3, v4, origin_a, origin_b, lat_origin
    real*8                                :: angle, pi, d1
    integer                               :: i, j, num_columns, row, i0

    pi = 4.d0 * atan(1.d0)
    angle = (30.d0/180.d0) * pi
    d1 = a/sqrt(3.d0)

    v1 = [1.5d0*d1, 0.d0, 0.d0]
    v2 = [a*cos(angle), a*sin(angle), 0.d0]
    v3 = [0.d0, a, 0.d0]

    origin_a = [0.d0, 0.d0, z]
    origin_b = origin_a + [d1, 0.d0, 0.d0]

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
    lattice(1,:) = [0.d0, 0.d0, z]

    i0 = ((row/2) * num_columns) + (num_columns/2)
    lat_origin = lattice(i0,:)
    lattice(:,1) = lattice(:,1) - lat_origin(1)
    lattice(:,2) = lattice(:,2) - lat_origin(2)
    lattice(:,3) = lattice(:,3) - lat_origin(3)
end subroutine

function rotate_point(point, pivot, angle) result(new_point)
    real*8, dimension(3), intent(in) :: point
    real*8, intent(in)               :: angle
    real*8, dimension(3)             :: aux1, aux2, pivot, new_point

    aux1 = point - pivot
    
    aux2(1) = aux1(1)*cos(angle) - aux1(2)*sin(angle)
    aux2(2) = aux1(1)*sin(angle) + aux1(2)*cos(angle)
    aux2(3) = aux1(3)

    new_point = aux2 + pivot
end function

subroutine write_lattice(lattice, filename)
    real*8, dimension(:,:), intent(in) :: lattice
    character(13), intent(in)          :: filename
    integer                            :: i

    open(23, file=filename)
    do i=1, size(lattice, 1)
        write(23,*) lattice(i,1), ";", lattice(i,2), ";", lattice(i,3)
    enddo
    close(23)
end subroutine

end program symmetry
