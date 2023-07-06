program symmetry
implicit none

real*8, dimension(:,:), allocatable :: lat1, lat2
real*8                              :: rot_angle, lat_separation, lat_const
real*8                              :: pi
integer                             :: i, num_points

pi = 4.d0*atan(1.d0)

num_points      = 10000
lat_separation  = 10.d0
lat_const       = 1.d0
rot_angle       = (10.d0/180.d0) * pi

allocate(lat1(num_points, 3))
allocate(lat2(num_points, 3))
lat1 = 0.d0
lat2 = 0.d0

call create_square_lattice(lat1, 0.d0          , lat_const, 100)
call create_square_lattice(lat2, lat_separation, lat_const, 100)

do i=1, num_points
    lat2(i,:) = rotate_point(lat2(i,:), [0.d0, 0.d0, lat_separation], rot_angle)
enddo

call write_lattice(lat1, "lattice01.dat")
call write_lattice(lat2, "lattice02.dat")

contains

function exists_in_lattice(point, lattice) result(flag)
    real*8, dimension(3), intent(in)     :: point
    real*8, dimension(:,:), intent(in)   :: lattice
    real*8                               :: c1, c2, c3, tol
    logical, dimension(:), allocatable   :: condition_list
    logical                              :: flag
    integer                              :: i, cont

    allocate(condition_list(size(lattice,1)))

    tol = 1e-12
    do i=1, size(lattice, 1)
        c1 = abs(lattice(i,1) - point(1))
        c2 = abs(lattice(i,2) - point(2))
        c3 = abs(lattice(i,3) - point(3))
        condition_list(i) = (c1 < tol) .and. (c2 < tol) .and. (c3 < tol)
    enddo

    cont = count(condition_list .eqv. .true.)
    flag = any(condition_list .eqv. .true.)
end function

subroutine create_lattice(lattice, z, a)
    real*8, dimension(:,:), intent(inout) :: lattice
    real*8, intent(in)                    :: z, a
    real*8, dimension(3)                  :: p_aux
    real*8                                :: angle_rad, pi
    integer                               :: i, ind

    print*, "Creating lattice"

    pi = 4.d0 * atan(1.d0)
    lattice(:,3) = z
    lattice(1,:) = [0.d0, 0.d0, z]
    i = 2
    ind = 2

    do while(ind < size(lattice, 1))
        p_aux = lattice(i-1,:) + [a, 0.d0, 0.d0]
        if (.not. exists_in_lattice(p_aux, lattice)) then
            lattice(ind,:) = p_aux
            ind = ind + 1
        endif
        
        angle_rad = (60.d0/180.d0) * pi
        p_aux = rotate_point(lattice(i,:), lattice(i-1,:), angle_rad)
        if (.not. exists_in_lattice(p_aux, lattice)) then
            lattice(ind,:) = p_aux
            ind = ind + 1
        endif
        
        angle_rad = (120.d0/180.d0) * pi
        p_aux = rotate_point(lattice(i,:), lattice(i-1,:), angle_rad)
        if (.not. exists_in_lattice(p_aux, lattice)) then
            lattice(ind,:) = p_aux
            ind = ind + 1
        endif

        angle_rad = (180.d0/180.d0) * pi
        p_aux = rotate_point(lattice(i,:), lattice(i-1,:), angle_rad)
        if (.not. exists_in_lattice(p_aux, lattice)) then
            lattice(ind,:) = p_aux
            ind = ind + 1
        endif

        angle_rad = (240.d0/180.d0) * pi
        p_aux = rotate_point(lattice(i,:), lattice(i-1,:), angle_rad)
        if (.not. exists_in_lattice(p_aux, lattice)) then
            lattice(ind,:) = p_aux
            ind = ind + 1
        endif

        angle_rad = (300.d0/180.d0) * pi
        p_aux = rotate_point(lattice(i,:), lattice(i-1,:), angle_rad)
        if (.not. exists_in_lattice(p_aux, lattice)) then
            lattice(ind,:) = p_aux
            ind = ind + 1
        endif

        i = i + 1
    enddo
    lattice(1,:) = [0.d0, 0.d0, z]
end subroutine

subroutine create_square_lattice(lattice, z, a, num_columns)
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
