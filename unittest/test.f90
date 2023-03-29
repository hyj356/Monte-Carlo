program test
  use unittest, only: TestRho, TestF, TestZ
  implicit none

  call TestF('Cu_u3.eam')
  CALL TestZ('Cu_u3.eam')
  call TestRho('Cu_u3.eam')

end program test
