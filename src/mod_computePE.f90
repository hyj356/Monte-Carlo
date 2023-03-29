module computeUE
  use global_var, only: wp, stdout
  use fileIO, only: eamFile, Atom, region, AllocMem, eamAlloyFile
  implicit none
  private
  public :: Interpolate, Energy, Distance, SwapType
  logical, private, save :: isNeighbored = .false.    !! 用来判断当前模型的近邻列表是否已经建立, 避免反复建立邻居列表增加运行时间

  interface Pair
    module procedure :: PairEam, PairEamAlloy
  end interface Pair

  interface Prof
    module procedure :: ProfEam, ProfEamAlloy
  end interface Prof

  interface Embed
    module procedure :: EmbedEam, EmbedEamAlloy
  end interface Embed

  interface Energy
    module procedure :: EamEnergy, EamAlloyEnergy
  end interface Energy

contains

  pure function Interpolate( dri, x, y, xi ) result(res)
    !! 这里因为是只需要对势函数进行插值, 而势函数的数组是确定的, 所以传入的数据有所不同
    real(wp), intent(in) ::  dri                  !< dri等于1/dr, 为了避免频繁的计算除法, 且dr可以早早从势函数文件中读取, 因此这里提前计算好
    real(wp), dimension(:), intent(in) :: x       !< 用于插值的数组, x取值
    real(wp), dimension(:), intent(in) :: y       !< 用于插值的数组, y取值
    real(wp), intent(in)               :: xi      !< 被插的点
    real(wp) :: res

    integer                            :: prev    !< 数组中位于xi之前的下标
    integer                            :: next    !< 数组中位于xi之后的下标
    integer                            :: xdim    !< x数组的维数
    integer                            :: ydim    !< y数组的维数

    !! 在实际运算的时候, 必须确保原子之间距离小于截断半径才会进行运算, 所以无需担心待插的xi超出范围
    !! 寻找传入距离对应的2个下标
    prev = int(xi * dri) + 1      !! 如果太靠近起点, 这里是0, 但是Fortran数组从1开始
    next = prev + 1               !! 所以这里要注意

    if (xi < maxval(x)) then      !! 如果在范围以内就进行线性内插
      res = y(prev) + ((xi - x(prev)) * (y(next) - y(prev))) * dri        !! 进行内插
    else                          !! 在范围以外就利用最后2个点进行线性外插, lammps中也是线性外插处理, 但不知道算法是否相同
      xdim = size(x); ydim = size(y)
      res = y(ydim) + ((xi - x(xdim)) * (y(ydim) - y(ydim - 1))) * dri    !! 进行线性外推, 在lammps的源代码中也是一样的做法
    end if

  end function Interpolate

  pure subroutine SwapType(Atoms, i, j)
    !! 交换第i个和第j个原子
    type(Atom), intent(inout), dimension(:) :: Atoms  !< 存储原子坐标和邻居列表的ATOM类数组
    integer, intent(in) :: i, j                       !< 需要交换原子类型的第i个和第j个原子
    integer :: tempID                                 !< 用于存储的临时的类型id

    tempID = Atoms(i)%typeId
    Atoms(i)%typeId = Atoms(j)%typeId
    Atoms(j)%typeId = tempID

  end subroutine SwapType

  pure function PairEam(eamPot, r) result(res)
    !! 根据传入的2个原子之间的距离线性内插计算2体势
    type(eamFile), intent(in) :: eamPot  !< 一个包含了势函数所有数据和处理好数据之后的类
    real(wp), intent(in):: r             !< 2个原子之间的距离
    real(wp) :: EffectiveZ               !< 用于计算2体势φ的有效电荷值
    real(wp) :: res                      !< 计算对势值, 单位为ev

    !! 根据从势函数读取到的Z数组, 内插获得在r距离下的有效电荷值, 注意这里应该使用Z2R数组
    EffectiveZ = interpolate(eamPot%dri, eamPot%rx, eamPot%Z2R, r)  !! 这里用一个中间变量是为了避免平方和2次计算interpolate函数的额外运算
    res = EffectiveZ * EffectiveZ / r    !! 计算真正的对势值

  end function PairEam

  pure function PairEamAlloy(eampot, itype, jtype, r) result(res)
    !! 根据eam/alloy势函数文件和2个原子的类型和之间的距离线性内插计算2个原子之间的两体势
    !! 注意: 这里的ϕ(r)不像EAM那样是有效电荷Z(r)，也不是真的能量ϕ(r)，而是r×ϕ(r)。
    type(eamAlloyFile), intent(in) :: eamPot  !< 一个包含了势函数所有数据和处理好数据之后的类
    real(wp), intent(in):: r             !< 2个原子之间的距离
    integer, intent(in) :: itype         !< 主原子的原子类型id
    integer, intent(in) :: jtype         !< 近邻原子的原子类型id                    
    real(wp) :: res                      !< 计算对势值, 单位为ev
    real(wp) :: phiR                     !< 这个是内插结果, 还需要除以原子之间的距离才是对体势
    
    phiR = interpolate(eamPot%dri, eamPot%rx, eamPot%Z(:, itype, jtype), r)   !! 内插计算phiR
    res = phiR / r    !! 计算真正的两体势

  end function PairEamAlloy

  pure function ProfEam(eampot, r) result(res)
    !! 此程序用于计算电子密度, 根据2个原子之间的距离线性内计算电子密度
    type(eamFile), intent(in) :: eamPot  !< 一个包含了势函数所有数据和处理好数据之后的类
    real(wp), intent(in):: r             !< 2个原子之间的距离
    real(wp) :: res

    res = interpolate(eamPot%dri, eamPot%rx, eamPot%Rho, r)

  end function ProfEam

  pure function ProfEamAlloy(eampot, jtype, r) result(res)
    !! 根据eam/alloy势函数文件和近邻原子的类型和之间的距离线性内插计算2个原子之间的两体势
    type(eamAlloyFile), intent(in) :: eamPot  !< 一个包含了势函数所有数据和处理好数据之后的类
    real(wp), intent(in):: r             !< 2个原子之间的距离
    integer, intent(in) :: jtype         !< 近邻原子的原子类型id                    
    real(wp) :: res                      !< 计算对势值, 单位为ev

    res = interpolate(eamPot%dri, eamPot%rx, eamPot%Rho(:, jtype), r)

  end function ProfEamAlloy

  pure function EmbedEam(eampot, rho) result(res)
    !! 根据传入的电子密度ρ(rho), 线性内插计算嵌入能F
    type(eamFile), intent(in) :: eamPot  !< 一个包含了势函数所有数据和处理好数据之后的类
    real(wp), intent(in):: rho           !< 中心原子与近邻原子形成的电子密度
    real(wp) :: res

    res = interpolate(eamPot%drhoi, eamPot%rhox, eamPot%F, rho)

  end function EmbedEam

  pure function EmbedEamAlloy(eampot, rho, itype) result(res)
    !! 根据传入的电子密度ρ(rho)和势函数文件, 线性内插计算嵌入能F
    type(eamAlloyFile), intent(in) :: eamPot  !< 一个包含了势函数所有数据和处理好数据之后的类
    real(wp), intent(in):: rho           !< 中心原子与近邻原子形成的电子密度
    integer, intent(in) :: itype         !< 主原子的原子类型id
    real(wp) :: res

    res = interpolate(eamPot%drhoi, eamPot%rhox, eamPot%F(:, itype), rho)

  end function EmbedEamAlloy

  pure function Distance(Atoms, i, j, box) result(res)
    !! 考虑周期性边界条件计算2个原子之间的距离
    type(Atom), intent(in), dimension(:) :: Atoms     !< 存储原子坐标和邻居列表的ATOM类数组
    integer, intent(in) :: i, j                       !< 2个原子之间的下标
    type(region), intent(in)  :: box                  !< 仿真盒子, 包含了坐标上下限和三个维度长度的信息
    real(wp) :: dx, dy, dz                            !< xyz方向上的差距, 考虑周期性边界条件
    real(wp) :: res                                   !< 考虑周期性边界条件下第i个与第j个原子之间的距离


    dx = Atoms(i)%x - Atoms(j)%x
    dx = dx - nint(dx / box%lx) * box%lx
    dy = Atoms(i)%y - Atoms(j)%y
    dy = dy - nint(dy / box%ly) * box%ly
    dz = Atoms(i)%z - Atoms(j)%z
    dz = dz - nint(dz / box%lz) * box%lz

    res = dx * dx + dy * dy + dz * dz

  end function Distance

  pure subroutine BuildList(Atoms, cutoff2, box)
    !! 构建每个原子的邻居列表, 使用verlet方法
    type(Atom), intent(inout), allocatable :: Atoms(:) !< 存储原子坐标和邻居列表的ATOM类数组
    real(wp), intent(in) :: cutoff2                    !< 建立邻居列表的截断半径的平方
    type(region), intent(in)  :: box                   !< 仿真盒子, 包含了坐标上下限和三个维度长度的信息
    real(wp) :: dR                                     !< 两个原子之间的距离
    integer :: nAtoms                                  !< 原子个数
    integer :: i, j                                    !< 循环变量

    !! 初始化变量
    nAtoms = size(Atoms)

    !! 开始verlet循环构建邻居列表, 时间复杂度为N(N - 1)/2, N为原子数量
    do i = 1, nAtoms - 1
      do j = i + 1, nAtoms
        dR = Distance(Atoms, i, j, box)
        if (dR < cutoff2) then
          Atoms(i)%slave = [Atoms(i)%slave, j]    !! Fortran2003特性, 类似于python的list类型的append, 但因涉及反复申请内存, 效率极低
          Atoms(j)%slave = [Atoms(j)%slave, i]    !! 比较好的折中办法是使用指针和派生数据类型搭建链表,
        end if
      end do
    end do

    !! 邻居列表构建完毕
  end subroutine BuildList

  pure subroutine EamEnergy(eampot, Atoms, power, box)
    !! 根据eam势函数文件与模型文件计算模型的势能
    type(eamFile), intent(in) :: eampot                 !< 一个包含了势函数所有数据和处理好数据之后的类
    type(Atom), intent(inout), allocatable :: Atoms(:)  !< 存储原子坐标和邻居列表的ATOM类数组
    real(wp), intent(out) :: power                      !< 输出的模型的能量
    type(region), intent(in)  :: box                    !< 仿真盒子, 包含了坐标上下限和三个维度长度的信息
    real(wp) :: phi     !< 两体势
    real(wp) :: rho     !< 电子密度
    !real(wp) :: emb    !< 嵌入势, 其实计算用不太上, 但还是放到这里提醒一下读者
    real(wp) :: dR      !< 两个原子之间的距离
    integer :: i, j

    !! 首先构建邻居列表
    call BuildList(Atoms, eampot%cutoff2, box)

    !! 初始化每一个变量
    phi = 0.d0; rho = 0.d0;
    power = 0.d0
    !! 第一层循环遍历每一个原子, 第二层循环遍历每一个原子的近邻
    do i = 1, size(Atoms)
      do j = 1, size(Atoms(i)%slave)
        dR = sqrt(Distance(Atoms, i, Atoms(i)%slave(j), box))
        phi = phi + pair(eampot, dR)
        rho = rho + prof(eampot, dR)
      end do
      power = power + Embed(eampot, rho)
      rho = 0.d0
    end do
    power = power +  0.5d0 * phi

  end subroutine EamEnergy

  subroutine EamAlloyEnergy(eampot, Atoms, power, box)
    !! 根据eam势函数文件与模型文件计算模型的势能
    type(eamAlloyFile), intent(in) :: eampot            !< 一个包含了势函数所有数据和处理好数据之后的类
    type(Atom), intent(inout), allocatable :: Atoms(:)  !< 存储原子坐标和邻居列表的ATOM类数组
    real(wp), intent(out) :: power                      !< 输出的模型的能量
    type(region), intent(in)  :: box                    !< 仿真盒子, 包含了坐标上下限和三个维度长度的信息
    real(wp) :: phi     !< 两体势
    real(wp) :: rho     !< 电子密度
    !real(wp) :: emb    !< 嵌入势, 其实计算用不上, 但还是放到这里提醒一下读者
    real(wp) :: dR      !< 两个原子之间的距离
    integer :: i, j

    !! 首先检查有没有构建邻居列表, 没有的话就构建一次, 一般来说只需要构建一次即可
    if (.not. isNeighbored) then
      call BuildList(Atoms, eampot%cutoff2, box)
      isNeighbored = .true.
    end if 

    !! 初始化每一个变量
    phi = 0.d0; rho = 0.d0;
    power = 0.d0; dR = 0.d0
    !! 第一层循环遍历每一个原子, 第二层循环遍历每一个原子的近邻
    !$OMP parallel do default(firstprivate) shared(Atoms) reduction(+:power, phi)
    do i = 1, size(Atoms)
      do j = 1, size(Atoms(i)%slave)
        dR = sqrt(Distance(Atoms, i, Atoms(i)%slave(j), box)) !! 考虑周期性边界条件计算原子之间的距离
        phi = phi + pair(eampot, Atoms(i)%typeId, Atoms(Atoms(i)%slave(j))%typeId, dR)  !! 计算2个原子之间的两体势
        rho = rho + prof(eampot, Atoms(Atoms(i)%slave(j))%typeId, dR)   !! 计算2个原子之间的电子密度
      end do
      !! 根据计算出来的rho计算原子的嵌入能, 并记得清0第i个原子的电子密度rho
      power = power + Embed(eampot, rho, Atoms(i)%typeId)
      rho = 0.d0
    end do
    !$OMP end parallel do
    power = power +  0.5d0 * phi

  end subroutine EamAlloyEnergy
end module
