module random 
  use global_var, only: wp
  implicit none
  private
  public :: RandomInt

contains
  function RandomInt(lowLimit, upLimit) result(res)
    integer, intent(in) :: upLimit, lowLimit  !< 分别对于想要生成的整数的上限和下限
    integer :: up, low
    integer :: res(2)
    real(wp) :: tempI, tempJ

    up = max(upLimit, lowLimit)
    low = min(upLimit, lowLimit)
    do !! 随机生成2个原子ID, 用于交换原子
      call Random_Number(tempI)
      call Random_Number(tempJ)
      res(1) = nint(tempI * up)
      res(2) = nint(tempJ * up)
      if (res(1) /= res(2) .and. All(res > low)) exit
    end do
  end function RandomInt
end module random