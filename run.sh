#!/bin/bash

# 此脚本的用法: 首先赋予其可直接执行的权限
# chmod 755 run.sh
# 然后我们输入2个参数, 一个是并行核数, 一个是输入参数的文件名称, 比如
# ./run.sh 40 parameter.nml
# 意思就是使用40核来运行程序, 包含MC迭代信息文件的名字是./src/parameter.nml

export OMP_NUM_THREADS=$1
./src/mc $2