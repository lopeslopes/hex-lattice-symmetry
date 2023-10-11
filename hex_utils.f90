module hex_utils
implicit none

real*16, parameter, public :: pi = 4.e0_16*atan(1.e0_16)
real*16, parameter, public :: tol = 1.e-8_16

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

subroutine make_std_rotation_matrix(matrix, axis_index, angle)
    real*16, dimension(3,3), intent(inout)  :: matrix
    real*16, intent(in)                     :: angle
    integer                                 :: axis_index
    
    select case (axis_index)
        case(1)
            matrix(1,:) = [1.e0_16,    0.e0_16,     0.e0_16]
            matrix(2,:) = [0.e0_16, cos(angle), -sin(angle)]
            matrix(3,:) = [0.e0_16, sin(angle),  cos(angle)]
        case(2)
            matrix(1,:) = [ cos(angle), 0.e0_16, sin(angle)]
            matrix(2,:) = [    0.e0_16, 1.e0_16,    0.e0_16]
            matrix(3,:) = [-sin(angle), 0.e0_16, cos(angle)]
        case(3)
            matrix(1,:) = [cos(angle), -sin(angle), 0.e0_16]
            matrix(2,:) = [sin(angle),  cos(angle), 0.e0_16]
            matrix(3,:) = [   0.e0_16,     0.e0_16, 1.e0_16]
    end select
end subroutine

subroutine rotate_lattice_2(lattice, angle, pivot, axis)
    real*16, dimension(:,:), intent(inout)  :: lattice
    real*16, dimension(3), intent(in)       :: pivot, axis
    real*16, intent(in)                     :: angle
    real*16, dimension(3,3)                 :: uxu, skew_u, rotation, identity
    integer                                 :: i

    identity(1,:) = [1.e0_16, 0.e0_16, 0.e0_16]
    identity(2,:) = [0.e0_16, 1.e0_16, 0.e0_16]
    identity(3,:) = [0.e0_16, 0.e0_16, 1.e0_16]

    uxu(1,:) = [axis(1)*axis(1), axis(1)*axis(2), axis(1)*axis(3)]
    uxu(2,:) = [axis(2)*axis(1), axis(2)*axis(2), axis(2)*axis(3)]
    uxu(3,:) = [axis(3)*axis(1), axis(3)*axis(2), axis(3)*axis(3)]

    skew_u(1,:) = [ 0.e0_16, -axis(3),  axis(2)]
    skew_u(2,:) = [ axis(3),  0.e0_16, -axis(2)]
    skew_u(3,:) = [-axis(2),  axis(1),  0.e0_16]

    rotation = cos(angle/2.e0_16)*identity + &
             & (1.e0_16 - sin(angle/2.e0_16))*uxu + &
             & sin(angle/2.e0_16)*skew_u

    lattice = lattice - spread(pivot, dim=1, ncopies=size(lattice,1))
    do i=1, size(lattice,1)
        lattice(i,:) = matmul(rotation,lattice(i,:))
    enddo
    lattice = lattice + spread(pivot, dim=1, ncopies=size(lattice,1))
end subroutine

subroutine sym_operation(index, lattice, pivot)
    real*16, dimension(:,:), intent(inout)  :: lattice
    real*16, dimension(3), intent(in)       :: pivot
    integer, intent(in)                     :: index
    real*16                                 :: rot_angle
    integer                                 :: axis

    ! TABLE OF CORRESPONDENCE OF INDEXES TO SYM. OP.
    !  1 --> E
    !  2 --> C6_1
    !  3 --> C6_2
    !  4 --> C3_1
    !  5 --> C3_2
    !  6 --> C2
    !  7 --> C'2_1
    !  8 --> C'2_2
    !  9 --> C'2_3
    ! 10 --> C''2_1
    ! 11 --> C''2_2
    ! 12 --> C''2_3

    select case (index)
        case(1) ! E
            rot_angle = 0.e0_16
            axis = 1
        case(2) ! C6_1
            rot_angle = pi/3.e0_16
            axis = 3
        case(3) ! C6_2
            rot_angle = (5.e0_16*pi)/3.e0_16
            axis = 3
        case(4) ! C3_1
            rot_angle = (2.e0_16*pi)/3.e0_16
            axis = 3
        case(5) ! C3_2
            rot_angle = (4.e0_16*pi)/3.e0_16
            axis = 3
        case(6) ! C2
            rot_angle = pi
            axis = 3
        case(7) ! C'2_1
            rot_angle = pi
            axis = 1
        case(8) ! C'2_2 FIX AXIS TO USE THIS
            rot_angle = 0.e0_16
            axis = 1
        case(9) ! C'2_3 FIX AXIS TO USE THIS
            rot_angle = 0.e0_16
            axis = 1
        case(10) ! C''2_1 FIX AXIS TO USE THIS
            rot_angle = 0.e0_16
            axis = 1
        case(11) ! C''2_2
            rot_angle = pi
            axis = 2
        case(12) ! C''2_3 FIX AXIS TO USE THIS
            rot_angle = 0.e0_16
            axis = 1
    end select

    call rotate_lattice(lattice, rot_angle, pivot, axis)
end subroutine

end module
