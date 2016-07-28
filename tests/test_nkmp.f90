module test_nkmp
  use, intrinsic :: iso_fortran_env, only: wp => real64

  use fruit
  implicit none

contains

  
  subroutine test_nkmp_1
    
    use model_nkmp_t, only: model

    type(model) :: nkmp

    real(wp) :: p0(13)

    nkmp = model()

    p0 = [2.09_wp,0.98_wp,2.25_wp,0.65_wp,0.81_wp,0.98_wp,0.93_wp,0.34_wp,3.16_wp,0.51_wp,0.19_wp,0.65_wp,0.24_wp]

    call assert_equals(-306.20_wp, nkmp%lik(p0), 0.01_wp)

  end subroutine test_nkmp_1

end module test_nkmp
