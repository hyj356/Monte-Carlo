FC = gfortran
FFLAGS = -O3 -Wall -Wno-unused-function #-g
FILES = global_var.f90 mod_fileIO.f90 mod_computePE.f90 main.f90
TESTS = global_var.f90 mod_fileIO.f90 mod_computePE.f90 mod_unittest.f90 test.f90

.PHONY: target clean

target:
	${FC} ${FILES} -o $@ ${FFLAGS}	# 计算模型的势能

test:
	${FC} ${TESTS} -o $@ ${FFLAGS}	# 用来做单元测试, 看看线性插值是否起效的程序

clean:
	rm -rf *.mod