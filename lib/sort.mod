module sort
    implicit none
    contains

    recursive subroutine sort_int8(a, n)
        integer(8), intent(inout) :: a(:)
        integer, intent(in) :: n
        integer(8) :: pivot, tmp
        integer :: i, j

        if (n < 2) return
        pivot = a(n/2 + 1)
        i = 1
        j = n
        do
            do while (a(i) < pivot); i = i + 1; end do
            do while (a(j) > pivot); j = j - 1; end do
            if (i >= j) exit
            tmp = a(i); a(i) = a(j); a(j) = tmp
            i = i + 1
            j = j - 1
        end do
        call sort_int8(a(1:i-1), i-1)
        call sort_int8(a(j+1:n), n-j)
    end subroutine

    ! quicksort.f -*-f90-*-
    ! Author: t-nissie
    ! License: GPLv3
    ! Gist: https://gist.github.com/t-nissie/479f0f16966925fa29ea
    !

    recursive subroutine quicksort(a, first, last)
    implicit none
    real*8  a(*), x, t
    integer first, last
    integer i, j

    x = a( (first+last) / 2 )
    i = first
    j = last
    do
        do while (a(i) < x)
            i=i+1
        end do
        do while (x < a(j))
            j=j-1
        end do
        if (i >= j) exit
        t = a(i);  a(i) = a(j);  a(j) = t
        i=i+1
        j=j-1
    end do
    if (first < i-1) call quicksort(a, first, i-1)
    if (j+1 < last)  call quicksort(a, j+1, last)
    end subroutine quicksort

    

end module