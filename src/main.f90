program main
  use random, only: RandomInt
  use global_var, only: wp, stdout, Bolz
  use fileIO, only:  Atom, ReadData, region, eamAlloyFile, &
      ReadEamAlloyFile, WriteData, InitialParameter,       &
      WriteLog  
  use computeUE, only: Energy, SwapType
  implicit none
  type(eamAlloyFile) :: Potential        !< 势函数文件中的内容
  type(Atom), allocatable :: Model(:)    !< 模型文件
  type(region)  :: box                   !< 仿真盒子
  character(len=128) :: potentialName    !< 势函数文件名称
  character(len=128) :: modelName        !< 原子模型文件名称
  character(len=3)   :: cTHREADS         !< 字符串形式的线程数
  real(wp) :: Ep                         !< 迭代之前的能量(previous)
  real(wp) :: En                         !< 迭代之后的能量(next)
  real(wp) :: temprature                 !< MC迭代的特征温度
  real(wp) :: start, end                 !< 用于统计程序的运行时间
  real(wp) :: P                          !< Metropolis判据对应的可能性
  real(wp) :: posibility                 !< 一个在0-1之间的均匀分布的随机数, 与Metropolis判据对比判断是否接受交换
  real(wp) :: kT                         !< 物理意义是Metropolis指数部分的分母, 为了减小计算量值等于1/kT  
  integer :: iterSteps                   !< MC迭代步数
  integer :: i                           !< 循环变量
  integer :: atomID(2)                   !< 需要交换的2个原子的ID
  integer :: nAtoms                      !< 原子数量
  integer :: success = 0                 !< 交换成功的次数
  integer :: failed = 0                  !< 交换没成功的次数
  integer :: logID                       !< log文件的ID
  logical :: accepted = .false.          !< 判断交换原子前后的模型是否接受
  integer :: THREADS                     !< 并行线程数
  
  ! 结果: CoCuFeNiPd-0.txt: -35939.7479428878 ev
  call random_seed    !! 初始化随机种子
  call InitialParameter(potentialName, modelName, temprature, iterSteps)  !! 从命令行读取含有namelist参数的文件
  call ReadEamAlloyFile(potentialName, Potential)     !! 读取势函数文件里面的全部信息
  call ReadData(modelName, box, Model)                !! 读取模型文件里面的全部信息
  call Energy(Potential, Model, Ep, box)              !! 计算模型的初始势能
  write(stdout, '(T8, *(A, 8X))') 'Steps', 'Potential Energy(ev)', 'time(s)', 'Success', 'Failed'
  write(stdout, *)  0, Ep, 0.0000000000000d0, success, failed
  open(newunit=logID, file='log.mc', action='write')
  call WriteLog(logID, 0, Ep, 0.d0)
  
  nAtoms = size(Model)    !! 记录模型中含有几个原子
  kT = 1.d0 / Bolz * temprature
  call getenv('OMP_NUM_THREADS', cTHREADS)    !! 从环境变量中获取有几核并行
  read(cTHREADS, *) THREADS                   !! 将字符串变量转成整型变量
  call cpu_time(start)                        !! 开始计时

  !! 正式开始MC迭代, 最多迭代iterSteps次
  do i = 1, iterSteps
    atomID = RandomInt(1, nAtoms)   !! 随机生成2个原子ID用于交换
    call SwapType(Model, atomID(1), atomID(2))  !! 交换原子ID
    call Energy(Potential, Model, En, box)      !! 计算交换原子之后的势能
    accepted = En < Ep            !! 对比交换前后两次模型的势能, 如果减小了就接受
    if (.not. accepted) then      !! 如果发现模型势能没有变小, 进行Metropolis判别
      P = exp((Ep - En) * kT)     !! 计算Metropolis判据
      call random_number(posibility)  !! 生成0-1之间的均匀分布的随机数
      if (posibility <= P) then   !! 如果发现当前模型符合Metropolis判据, 也接受交换
        Ep = En                   !! 减小了之后, 记录这一次的交换原子之后的能量
        success = success + 1     !! 成功次数加一
      else                        !! 如果势能没有减小且不符合Metropolis判据, 交换回去
        call SwapType(Model, atomID(2), atomID(1))
        failed = failed + 1
      end if
    else
      Ep = En                   !! 减小了之后, 记录这一次的交换原子之后的能量
      success = success + 1     !! 成功次数加一
    end if

    if (mod(i, 1000) == 0) then !! 每1000步输出一次相关的信息
      call cpu_time(end)
      write(stdout, *) i, Ep, (end - start) / THREADS, success, failed
      call WriteLog(logID, i, Ep, (end - start) / THREADS)
      call WriteData('../model/MC_Model.lmp', Potential, box, Model)
    end if
  end do
  close(logID)
  
end program main
 !! 参考教程: 1. https://zhuanlan.zhihu.com/p/512798840?utm_source=wechat_session&utm_medium=social&s_r=0
 !!          2. https://docs.lammps.org/pair_eam.html
