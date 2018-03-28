cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c                   源程序结构说明 
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c	该程序有两个子程序，sub用来得到观测方程，sub1用来解法方程
c	编译时先gfortran -c sub.f;gfortran -c sub1.f;gfortran -c main.f
c	然后将其链接在一起gfortran -o main sub.o sub1.o main.o
c	便编译得到可执行程序main
c	角度直接使用角度制，如直角为90
c	注意，该程序中对于观测方程或法方程的Ax=b中的x即g,h排列是
c	按g(1),g(2),...h(1),h(2)....的顺序排列的，即先排完g,再排h,且对于恒0项的g,h不算在x中
c	这对于确定A的排列很重要
c	由于试算的时候并没有用估值，直接计算就可解出x，故没理解为什么要先设x估值
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c                    文件说明 
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c	文件说明:
c	 1.点的坐标是放在data_position文件中的，
c	 2.正演算得的x分量值是放在data_lb中的
c	反演至今只用到x分量，以后有需要时再扩展
c	声明！！！！！！！！！！！！！！！！！！！！！！！_2014_01_25
c	之前的文件说明失效了，即data_position和data_lb已经不是程序
c	的输入文件了，统一整合到data文件中，data文件的结构如下：
c	lon1 lat1 X_magnitude1
c	lon2 lat2 X_magnitude2
c	...
c	即每行分别是每个测点的经度、纬度、X分量强度。之间用空格格开。
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c                    使用说明   _2014_01_25添加 
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c	首先，由于为方便使用，之前便将源代码进行了稍微修改，所需文件也有
c	所改动，但没有细说，于是重新进行了说明（具体见上文件说明）
c	使用说明:
c	 1.在编译好的执行文件目录下按要求格式放好data文件(规范见文件说明)
c	 2.直接运行执行文件，默认将结果显示到屏幕上
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c                     变量
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c	变量			含义
c	r			地球平均半径
c	pi			圆周率
c	N_max			模型的阶数
c	data_N			观测点个数
c	LA			观测方程系数矩阵
c	Lb			观测方程右边观测值
c	NA			法方程系数阵
c	Nb			法方程右边值
c	colat			余纬(0:180)
c	colon			经度(0:360)
c	lat			纬度(-90:90)
c	lon			经度(-180:180)
c	a			观测点高度(地表为r)
c	p			勒让德函数
c	grad_p			勒让德函数梯度
c	i,j			程序循环变量
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


	program refutations
	use sub1
	use sub
	implicit none
	real*8,parameter::r=6371200d0,pi=atan(1d0)*4d0
	real*8,parameter::zero=0.0d0,e=0.00000000000000000000000001
	integer,parameter::N_max=13,data_N=60000
c	LA表示观测方程中的系数矩阵，LB为观测值,即(LA)x=Lb
	real*8 LA(data_N,N_max*(N_max+2)),Lb(data_N),x(N_max*(N_max+2))
c	NA表示法方程中系数矩阵，Nb为观测值，即(NA)x=Nb
	real*8 NA(N_max*(N_max+2),N_max*(N_max+2)),Nb(N_max*(N_max+2))
	real*8 colat,colon,a,lat,lon
	real*8 p(N_max*(N_max+3)/2+1)
	real*8 grad_p(N_max*(N_max+3)/2+1)
c	real*8 gh_gu(N_max*(N_max+2))
c	real*8 gh_v(N_max*(N_max+2)),gh(N_max*(N_max+2))
c	real*8 w(data_N)
	integer i,j,k,m,n

c	Following is to remind you how to use this program
	print *,"This program is to get coefficient of IGRF"
	print *,"PS: at most 13 orders of coefficient intotal"
	print *,"second release, by don on 2014/01/25"
	print *,"Any problem,contact dpdeng@whu.edu.cn"
	print *,"Be careful there must be a file named data,"
	print *,"in the current directory "
	print *,"file data is the only one we need."
	print *,"The content of data is like following"
	write(*,"(2A6,A26)")"lon1","lat1","magnitude_of_Xcomponent1"
	write(*,"(2A6,A26)")"lon2","lat2","magnitude_of_Xcomponent2"
	print *,"..."
	print *,""
	print *,"result is coming.Please wait"

c	open(14,file="data_lb")
c	open(12,file="data_position")
c	do i=1,data_N
c	read(14,*)lb(i)
c	read(12,*)lon,lat
c	colat=90-lat
c	colon=lon
c	colat=colat*pi/180d0
c	colon=colon*pi/180d0
c	call getLAofx(colat,colon,r,N_max,LA,data_N,i)
c	end do
c	close(12)
c	close(14)

	open(12,file="data")
	do i=1,data_N
	read(12,*)lon,lat,lb(i)
	colat=90-lat
	colon=lon
	colat=colat*pi/180d0
	colon=colon*pi/180d0
	call getLAofx(colat,colon,r,N_max,LA,data_N,i)
	end do
	close(12)

	NA=matmul(transpose(LA),LA)
c	w=uz_dim-matmul(la,gh_gu)
	Nb=matmul(transpose(LA),lb)
	call cg(NA,Nb,x,N_max*(N_max+2),e)

c	Following is to output
	print *,"coefficient file"
	print *,"D means degree,O means order"
	write(*,"(2A3,2A13)")"D","O","g","h"
	i=1
	k=1
	do n=1,N_max
	do m=0,n
	if(m>0) then
	write(*,"(2I3,2F13.1)")n,m,x(i),x(k+N_max*(N_max+3)/2)
	k=k+1
	else
	write(*,"(2I3,2F13.1)")n,m,x(i),zero
	end if
	i=i+1
	end do
	end do
	
	end program refutations
