module hex_utils
implicit none

real*16, parameter, public :: pi = 4.e0_16*atan(1.e0_16)

contains

function exists_in_lattice(point, lattice) result(flag)
    real*16, dimension(3), intent(in)     :: point
    real*16, dimension(:,:), intent(in)   :: lattice
    real*16                               :: tol
    logical                               :: flag
    real*16, dimension(3)                 :: vec1
    logical, dimension(:), allocatable    :: condition_list
    integer                               :: i

    tol = 1.e-3_16
    !flag = any(all(abs(spread(point, dim=1, ncopies=size(lattice,1)) - lattice) .lt. tol, dim=2))

    allocate(condition_list(size(lattice,1)))
    do i=1, size(lattice, 1)
        vec1 = abs(lattice(i,:) - point)
        condition_list(i) = all(vec1 .lt. tol)
    enddo
    flag = any(condition_list .eqv. .true.)
end function

! CREATES LATTICE OF BODY-CENTERED HEXAGONS
subroutine create_lattice_bch(lattice, z, a)
    real*16, dimension(:,:), intent(inout) :: lattice
    real*16, intent(in)                    :: z, a
    real*16, dimension(3)                  :: v1, v2, origin_gen, lat_origin
    real*16                                :: angle, size_x, size_y
    integer                                :: i, j, num_columns, row, i0

    num_columns = int(sqrt(real(size(lattice,1),16)))
    angle = (120.e0_16/180.e0_16) * pi
    v1 = [a, 0.e0_16, 0.e0_16]
    v2 = [a*cos(angle), a*sin(angle), 0.e0_16]

    origin_gen = [0.e0_16, 0.e0_16, z]

    row = 1
    lattice(1,:) = origin_gen
    i = 2
    do while(i < size(lattice,1))
        do j=1, num_columns
            if (i >= size(lattice,1)) then
                exit
            endif
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
subroutine create_lattice_eh(lattice, z, a)
    real*16, dimension(:,:), intent(inout) :: lattice
    real*16, intent(in)                    :: z, a
    real*16, dimension(3)                  :: v1, v2, v3, d1, origin_a, origin_b, lat_origin
    real*16                                :: angle, d
    integer                                :: i, j, num_columns, row, i0

    num_columns = int(sqrt(real(size(lattice,1),16)))
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
    do while(i < size(lattice,1)-1)
        do j=1, num_columns
            if (i >= size(lattice,1)-1) then
                exit
            endif
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

subroutine rotate_lattice(lattice, angle)
    real*16, dimension(:,:), intent(inout)  :: lattice
    real*16, intent(in)                     :: angle
    real*16, dimension(3)                   :: pivot
    integer                                 :: i

    pivot = [0.e0_16, 0.e0_16, lattice(1,3)]
    do i=1, size(lattice,1)
        lattice(i,:) = rotate_point(lattice(i,:), pivot, angle)
    enddo
end subroutine

subroutine trim_lattice(lattice, x_min, x_max, y_min, y_max, trimmed_size, aux_lattice)
    real*16, dimension(:,:), intent(in)    :: lattice
    real*16, dimension(:,:), intent(inout) :: aux_lattice
    real*16, dimension(:,:), allocatable   :: trimmed_lattice
    real*16, intent(in)                    :: x_min, x_max, y_min, y_max
    logical                                :: cond_x, cond_y
    integer, intent(inout)                 :: trimmed_size
    integer                                :: i, j

    j = 0
    do i=1, size(lattice,1)
        cond_x = (lattice(i,1) >= x_min) .and. (lattice(i,1) <= x_max)
        cond_y = (lattice(i,2) >= y_min) .and. (lattice(i,2) <= y_max)

        if (cond_x .and. cond_y) then
            j = j + 1
            aux_lattice(j,:) = lattice(i,:)
        endif
    enddo

    trimmed_size = j
end subroutine

subroutine write_lattice(lattice, filename)
    real*16, dimension(:,:), intent(in) :: lattice
    character(len=13), intent(in)       :: filename
    integer                             :: i

    open(23, file=filename)
    do i=1, size(lattice, 1)
        write(23,*) lattice(i,1), ";", lattice(i,2), ";", lattice(i,3)
    enddo
    close(23)
end subroutine

end module
