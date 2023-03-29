module model
  !! 此模组包含用于建立FCC, BCC单晶模型的子程序, 目前还没有来得及完成
  use global_var, only: wp
  implicit none

  private
  public :: Atom
  type Atom
    real(wp) :: x   !! 原子的x坐标
    real(wp) :: y   !! 原子的y坐标
    real(wp) :: z   !! 原子的z坐标
    integer, allocatable :: slave(:)    !! 原子的所有在截断半径以内的邻居
  end type

contains

end module model
