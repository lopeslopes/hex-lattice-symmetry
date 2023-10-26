module hex_utils
implicit none

real*16, parameter, public :: pi = 4.e0_16*atan(1.e0_16)
real*16, parameter, public :: tol = 2.4e-3_16

contains

function magic_angle(ind) result(angle)
    integer, intent(in) :: ind
    real*16             :: angle

    angle = acos((real(3*ind*ind + 3*ind, 16) + 0.5e0_16)/real(3*ind*ind + 3*ind + 1, 16))
end function

function exists_in_lattice(point, lattice) result(flag)
    real*16, dimension(3), intent(in)     :: point
    real*16, dimension(:,:), intent(in)   :: lattice
    logical                               :: flag
    real*16, dimension(3)                 :: vec1
    logical, dimension(:), allocatable    :: condition_list
    integer                               :: i

    allocate(condition_list(size(lattice,1)))
    do i=1, size(lattice, 1)
        vec1 = abs(lattice(i,:) - point)
        condition_list(i) = all(vec1 .lt. tol)
    enddo
    flag = any(condition_list .eqv. .true.)
end function

subroutine create_honeycomb_lattice(latticeA, latticeB, z, a, ab_stacking)
    real*16, dimension(:,:), intent(inout) :: latticeA, latticeB
    real*16, intent(in)                    :: z, a
    real*16, dimension(3)                  :: v1, v2, v3, d1, origin_a, origin_b, lat_origin
    real*16                                :: angle, d
    integer                                :: i, j, num_columns, row, i0
    logical, intent(in)                    :: ab_stacking

    num_columns = 2*int(sqrt(real(size(latticeA,1)*2,16)))/3
    angle = (60.e0_16/180.e0_16) * pi

    v1 = [a, 0.e0_16, 0.e0_16]
    v2 = [a*cos(angle), a*sin(angle), 0.e0_16]
    v3 = [-a*cos(angle), a*sin(angle), 0.e0_16]

    d = sqrt((a**2)/(2.e0_16*(1.e0_16-cos(2.e0_16*angle))))
    d1 = [d*cos(angle/2.e0_16), d*sin(angle/2.e0_16), 0.e0_16]

    origin_a = [0.e0_16, 0.e0_16, z]
    origin_b = origin_a + d1

    row = 1
    i = 1
    do while(i < size(latticeA,1))
        do j=1, num_columns
            if (i > size(latticeA,1)) then
                exit
            endif
            latticeA(i,:) = origin_a + real(j-1,16)*v1
            latticeB(i,:) = origin_b + real(j-1,16)*v1
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

    i0 = ((row/2) * num_columns) + (num_columns/2)
    lat_origin = [latticeA(i0,1), latticeA(i0,2), 0.e0_16]
    
    latticeA = latticeA - spread(lat_origin, dim=1, ncopies=size(latticeA,1))
    latticeB = latticeB - spread(lat_origin, dim=1, ncopies=size(latticeB,1))

    if (ab_stacking) then
        latticeA = latticeA + spread(d1, dim=1, ncopies=size(latticeA,1))
        latticeB = latticeB + spread(d1, dim=1, ncopies=size(latticeB,1))
    endif
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

subroutine distance_from_origin(lattice, origin, distances)
    real*16, dimension(:,:), intent(in)  :: lattice
    real*16, dimension(3), intent(in)    :: origin
    real*16, dimension(3)                :: aux_vec
    real*16, dimension(:), intent(inout) :: distances
    integer                              :: i

    do i=1, size(lattice,1)
        aux_vec = lattice(i,:) - origin
        if (sum(abs(aux_vec)) .lt. tol) then
            distances(i) = 1000.e0_16
        else
            distances(i) = sqrt(dot_product(aux_vec, aux_vec))
        endif
    enddo
end subroutine

function rotate_point(point, pivot, angle, axis) result(new_point)
    real*16, dimension(3), intent(in) :: point
    real*16, intent(in)               :: angle
    real*16, dimension(3)             :: aux1, aux2, pivot, new_point
    real*16, dimension(3,3)           :: rot_matrix
    integer, intent(in)               :: axis

    ! AXIS: 1=x, 2=y, 3=z

    select case (axis)
        case(1)
            rot_matrix(1,:) = [1.e0_16, 0.e0_16, 0.e0_16]
            rot_matrix(2,:) = [0.e0_16, cos(angle), -sin(angle)]
            rot_matrix(3,:) = [0.e0_16, sin(angle), cos(angle)]
        case(2)
            rot_matrix(1,:) = [cos(angle), 0.e0_16, sin(angle)]
            rot_matrix(2,:) = [0.e0_16, 1.e0_16, 0.e0_16]
            rot_matrix(3,:) = [-sin(angle), 0.e0_16, cos(angle)]
        case(3)
            rot_matrix(1,:) = [cos(angle), -sin(angle), 0.e0_16]
            rot_matrix(2,:) = [sin(angle), cos(angle), 0.e0_16]
            rot_matrix(3,:) = [0.e0_16, 0.e0_16, 1.e0_16]
    end select

    aux1 = point - pivot
    aux2 = matmul(rot_matrix, aux1)
    new_point = aux2 + pivot
end function

subroutine rotate_lattice(lattice, angle, pivot, axis)
    real*16, dimension(:,:), intent(inout)  :: lattice
    real*16, intent(in)                     :: angle
    real*16, dimension(3), intent(in)       :: pivot
    integer                                 :: i
    integer, intent(in)                     :: axis


    do i=1, size(lattice,1)
        lattice(i,:) = rotate_point(lattice(i,:), pivot, angle, axis)
    enddo
end subroutine

function find_angle(point1, point2) result(angle)
    real*16, dimension(3), intent(in)   :: point1, point2
    real*16                             :: angle
    real*16                             :: x1, x2, y1, y2, l1, l2, l3

    x1 = max(point1(1), point2(1))
    x2 = min(point1(1), point2(1))
    y1 = min(point1(2), point2(2))
    y2 = max(point1(2), point2(2))

    l1 = sqrt(dot_product(point1-point2, point1-point2))
    l2 = sqrt(dot_product(point1, point1))
    l3 = sqrt(dot_product(point2, point2))

    angle = acos((l2**2 + l3**2 - l1**2)/(2.e0_16*l2*l3))
end function

end module
