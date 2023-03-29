module unittest
  use global_var, only: wp, stdout
  use fileIO, only: eamFile, ReadEamfile, eamAlloyFile, ReadEamAlloyFile
  use computeUE, only: interpolate, Embed, Prof
  implicit none

  private
  public :: TestRho, TestF, TestZ

contains

  subroutine TestRho(filename)
    !! 测试势函数中对数组ρ的插值
    character(len=*), intent(in) :: filename    !< 势函数文件名称
    type(eamFile) :: eampot                     !< 包含了势函数所有信息和处理之后的数据的集合
    integer :: i                                !< 循环变量
    real(wp) :: rho                             !< 插值计算出来的ρ

    call ReadEamfile(filename, eampot)             !! 读取势函数文件

    do i = 0, eampot%Nr - 1
      rho = Prof(eampot, i * eampot%dr)         !! 调用Prof函数, 根据不同的dr进行插值计算并与文件中的同样位置的数值进行对比
      !! 如果插值和数组对应位置的误差过大, 说明插值不正确
      if (abs(rho - eampot%Rho(i+1)) > 1e-10) stop 'Interpolation error when interpolate array of rho.'
    end do

    write(stdout, '(A)') 'Interpolation of Rho successful!'

  end subroutine TestRho

  subroutine TestF(filename)
    !! 测试势函数中对数组F的插值
    character(len=*), intent(in) :: filename  !< 势函数文件名称
    type(eamFile) :: eampot                   !< 包含了势函数所有信息和处理之后的数据的集合
    integer :: i                              !< 循环变量
    real(wp) :: F                             !< 插值计算出来的ρ

    call ReadEamfile(filename, eampot)

    do i = 0, eampot%Nrho - 1
      F = Embed(eampot, i * eampot%drho)      !! 调用Embed函数, 根据不同的dr进行插值计算并与文件中的同样位置的数值进行对比
      !! 如果插值和数组对应位置的误差过大, 说明插值不正确
      if (abs(F - eampot%F(i+1)) > 1e-10) stop 'Interpolation error when interpolate array of F.'
    end do

    write(stdout, '(A)') 'Interpolation of F successful!'
  end subroutine TestF

  subroutine TestZ(filename)
    !! 测试势函数中对数组Z的插值
    character(len=*), intent(in) :: filename    !< 势函数文件名称
    type(eamFile) :: eampot                     !< 包含了势函数所有信息和处理之后的数据的集合
    integer :: i                                !< 循环变量
    real(wp) :: Z                               !< 插值计算出来的ρ

    call ReadEamfile(filename, eampot)

    do i = 0, eampot%Nr - 1
      Z = interpolate(eampot%dri, eampot%rx, eampot%Z, i * eampot%dr)   !! 调用interpolate函数, 根据不同的dr进行插值计算并与文件中的同样位置的数值进行对比
      !! 如果插值和数组对应位置的误差过大, 说明插值不正确
      if (abs(z - eampot%z(i+1)) > 1e-10) stop 'Interpolation error when interpolate array of Z.'
    end do

    write(stdout, '(A)') 'Interpolation of Z successful!'

  end subroutine TestZ

end module unittest
