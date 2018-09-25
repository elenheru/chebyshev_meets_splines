module      polynomial_description
  implicit none

  integer,parameter ::  degree_limit  = 100

  type      polynomial
    real(8),dimension(0:degree_limit-1) :: coefficient
    !value of coefficient near of x^degree is coefficient(degree)
    logical,dimension(0:degree_limit-1) :: is_coefficient_substant
    !coefficient(degree) forced to be assumed as zero if is_coefficient_substant .eqv. .false.
    integer :: degree
    !last non zero monome degree is declared degree
    !if degree .eq. -1 then result is incorrect
    !it happens when multiplicating two polynomials with high enough degree, so product cannot be stored
    !laurent generalization of polynomials is not implied to be
  end type  polynomial

end module  polynomial_description

module      polynomial_procedures_description
  use polynomial_description
  implicit none
  contains !polynomials proceeding functions
    type(polynomial) function polmul(polynomial_1,polynomial_2) result(polynomial_product)
      implicit none
      type(polynomial),intent(in)  ::  polynomial_1,polynomial_2
      integer i,j
        !multiplicates two polynomials
      if (polynomial_1%degree + polynomial_2%degree .gt. degree_limit) then
        print*,'polynomial degree overflow while multiplicating, can not process this'
        polynomial_product%degree = -1
        return
      endif
      if (polynomial_1%degree .lt. 0 .or. polynomial_2%degree .lt. 0) then
        print*,'negative polynomial degree on input, polmul can not process this'
        polynomial_product%degree = -1
        return
      endif
      polynomial_product%coefficient(0 : polynomial_1%degree + polynomial_2%degree) = 0d0
      do i=0,polynomial_1%degree
        do j=0,polynomial_2%degree
          polynomial_product%coefficient(i+j) = &
          polynomial_product%coefficient(i+j) + &
          polynomial_1%coefficient(i) + &
          polynomial_2%coefficient(j)
        enddo
      enddo
      polynomial_product%degree = polynomial_1%degree + polynomial_2%degree
      return
    end function  polmul

    type(polynomial) function polsum(polynomial_1,polynomial_2) result(polynomial_sum)
      implicit none
      type(polynomial),intent(in)  ::  polynomial_1,polynomial_2
      integer i
        !summates two polynomials polynomial_1 and polynomial_2
      if (polynomial_1%degree + polynomial_2%degree .gt. degree_limit) then
        print*,'polynomial degree overflow while summating, can not process this'
        polynomial_sum%degree = -1
        return
      endif
      if (polynomial_1%degree .lt. 0 .or. polynomial_2%degree .lt. 0) then
        print*,'negative polynomial degree on input, polsum can not process this'
        polynomial_sum%degree = -1
        return
      endif
      polynomial_sum%degree = max(polynomial_1%degree,polynomial_2%degree)
      do i = 0,polynomial_sum%degree
        polynomial_sum%coefficient(i) = polynomial_1%coefficient(i) + polynomial_2%coefficient(i)
      enddo
      return
    end function polsum

    type(polynomial) function poldif(polynomial_1) result(derivative_)
      implicit none
      type(polynomial),intent(in)  ::  polynomial_1
      integer i
        !differentiates input polynomial polynomial_1
      if (polynomial_1%degree .lt. 0 .or. polynomial_1%degree .gt. degree_limit) then
        print*,'negative polynomial degree on input, poldif can not process this'
        derivative_%degree = -1
        return
      endif
      derivative_%coefficient=0d0
      derivative_%degree=0
      if(polynomial_1%degree .eq. 0) return !differenting constant produces zero constant
      derivative_%degree = polynomial_1%degree-1
      do i = 0,derivative_%degree
        derivative_%coefficient(i) = polynomial_1%coefficient(i+1)*dfloat(i+1)
      enddo
      return
    end function poldif

    real(8) function polint(polynomial_1,from_x,to___x) result(integral_)
      implicit none
      type(polynomial),intent(in)  ::  polynomial_1
      real(8),intent(in)  :: from_x,to___x
      real(8)  :: from_x_power,to___x_power
      integer i
        !integrates input polynomial polynomial_1 from from_x to to___x
      if (polynomial_1%degree .lt. 0 .or. polynomial_1%degree .ge. degree_limit) then
        print*,'negative or too big polynomial degree on input, polint can not process this, returning zero value'
        integral_ = 0d0
        return
      endif
      integral_ = 0d0
      from_x_power = 1d0
      to___x_power = 1d0
      do i = 0,polynomial_1%degree
        from_x_power = from_x_power * from_x
        to___x_power = to___x_power * to___x
        integral_ = integral_+ polynomial_1%coefficient(i) * &
        (to___x_power - from_x_power)/dfloat(i+1)
      enddo !according to Newton–Leibniz axiom for polynomials
      return
    end function  polint

    real(8) function polval(polynomial_1,x) result(polynomial_value)
      implicit none
      type(polynomial),intent(in)  ::  polynomial_1
      real(8),intent(in)  :: x
        !evaluates input polynomial at x
      integer i

      polynomial_value=0d0
      if (polynomial_1%degree .lt. 0 .or. polynomial_1%degree .gt. degree_limit) then
        print*,'negative or too big polynomial degree on input, polval can not process this, returning zero'
        return
      endif
      polynomial_value=polynomial_1%coefficient(polynomial_1%degree)
      if(polynomial_1%degree .eq. 0) return
      do i = polynomial_1%degree, 1, -1
        polynomial_value = polynomial_value*x + polynomial_1%coefficient(i-1)
      enddo
      return
    end function  polval

end module  polynomial_procedures_description
