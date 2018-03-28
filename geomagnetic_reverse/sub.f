	module sub
	contains
	
c	This subroutine is used to calculate p and grad_p which use the data of 
c	latitude and longtitude 
c	x means lon,y means lat,N_max means Nmax here
	subroutine point(N_max,x,y,p,grad_p)
	implicit none
	real*8 x
	real*8 y
	real*8 p(N_max*(N_max+3)/2+1)
	real*8 grad_p(N_max*(N_max+3)/2+1)
	integer N_max,n,m,i
	
	p(1)=1
	p(2)=cos(y)
	p(3)=sin(y)
	grad_p(1)=0
	grad_p(2)=(cos(y)*p(2)-p(1))/sin(y)
	grad_p(3)=cos(y)*p(3)/sin(y)
	if (N_max>=2) then
	i=4
	do n=2,N_max
	do m=0,n
	if(m==n)then
	p(i)=p(i-n-1)*sin(y)*((2.0*m-1)/(2.0*m))**0.5
	grad_p(i)=n*cos(y)*p(i-n-1)*((2.0*m-1)/(2.0*m))**0.5
	else
	p(i)=((2*n-1)*cos(y)*p(i-n)-((n-1)**2-m**2)** 
     &	0.5*p(i-n-n+1))/(n**2-m**2)**0.5
	grad_p(i)=(n*cos(y)*p(i)-(n**2-m**2)** 
     &	0.5*p(i-n))/sin(y)
	end if                            
	i=i+1
	end do
	end do
	end if
	end subroutine point




c	得到观测方程的系数矩阵LA，j表示计算LA第几行，即第j个观测点的系数行。
	subroutine getLAofx(lat,lon,r,N_max,LA,data_N,j)
	implicit none
	integer N_max,j,data_N
	real*8 lat,lon,r,LA(data_N,N_max*(N_max+2))
	real*8 p(N_max*(N_max+3)/2+1),grad_p(N_max*(N_max+3)/2+1)
	integer i,n,m,k
	call point(N_max,lon,lat,p,grad_p)
c	ux=ux+(g(i)*cos(m*x)+h(i)*sin(m*x))*grad_p(i+1)
	i=1
	k=1
	do n=1,N_max
	do m=0,n
	LA(j,i)=cos(m*lon)*grad_p(i+1)
	if(m>0) then
	LA(j,k+N_max*(N_max+3)/2)=sin(m*lon)*grad_p(i+1)
	k=k+1
	end if
	i=i+1
	end do
	end do
	end subroutine
	


	end module
