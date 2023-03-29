module global_var
  use iso_fortran_env, only: real32, real64, input_unit, output_unit, int32, int64
  implicit none
  private
  integer, public, parameter :: wp = real64             !< 工作的精度(working precision), 默认是双精度, 如果把real64改成real32就是单精度
  integer, public, parameter :: stdin = input_unit      !< 标准输入通道的整数id
  integer, public, parameter :: stdout = output_unit    !< 标准输出通道的整数id
  integer, public, parameter :: si = int32              !< 短整型所占用的精度
  integer, public, parameter :: li = int64              !< 长整型所占用的精度
  integer, public, parameter :: NDIM = 3                !< 仿真维数
  integer, public, parameter :: THREADS = 4             !< openMP的并行线程数
  real(wp), public, parameter :: KB = 1.d0              !< 波兹曼常数(Bolzman constant)
  real(wp), public, parameter :: PI = 4.d0*atan(1.d0)   !< 圆周率π
  real(wp), public, parameter :: Hartree = sqrt(27.2d0) !< 单位ev
  real(wp), public, parameter :: Bohr = sqrt(0.529d0)   !< 单位埃米(Angstroms)
  real(wp), public, parameter :: Bolz =  0.0000861727708792736d0       !< 波兹曼常数: 1.380649e-23 J/K = 0.0000861727708792736 ev/k
  character(len=*), public, parameter :: filename = 'parameter.nml'    !< 输入文件的名称, 按照namelist的格式撰写

end module global_var
