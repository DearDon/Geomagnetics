!	This is a module file to help to compute 
!	It include 3 subroutines to make the program function well

       module sub
       contains
!       

!	This subroutine is used to calculate p and grad_p which use the data of 
!	latitude and longtitude from user
       subroutine point(N0,x,y,p,grad_p)
       real*8 x
       real*8 y
       real*8 p(N0*(N0+3)/2+1) !to get p,we need p(0),so total number+1
       real*8 sy,cy
       real*8 grad_p(N0*(N0+3)/2+1)!to fit with p,total number also plus one
       integer N0,n,m,i !we use n,m,i in loop
       
	sy=sin(y)
	cy=cos(y)
       p(1)=1
       p(2)=cy
       p(3)=sy
       grad_p(1)=0!we don't use it,give it whatever you want
       grad_p(2)=(cos(y)*p(2)-p(1))/sin(y)
       grad_p(3)=cos(y)*p(3)/sin(y)
       i=4
       do n=2,N0
       	do m=0,n
       		if(m==n	)then
       			p(i)=p(i-n-1)*sy*((2.0d0*m-1)/(2.0d0*m))**0.5
       			grad_p(i)=n*cy*p(i-n-1)*((2.0d0*m-1)/(2.0d0*m))**0.5
       		else
       			p(i)=((2.0d0*n-1)*cy*p(i-n)-((n-1)**2-m**2)** &
       		&	0.5*p(i-n-n+1))/(n**2-m**2)**0.5
       			grad_p(i)=(n*cy*p(i)-(n**2-m**2)**&
       		&	0.5*p(i-n))/sy
       		end if                            
       		i=i+1
       	end do
       end do			
       end subroutine point



!	得到观测方程的系数矩阵NA，j表示计算NA第几行，即第j个观测点的系数行。
!	subroutine getNAofx(lat,lon,r,Nmax,NA,j)
!	real*8 lat,lon,r,NA(Nmax*(Nmax+3),Nmax*(Nmax+3))
!	integer Nmax,j,i
!	real*8 p(Nmax*(Nmax+3)/2+1),grad_p(Nmax*(Nmax+3)/2+1)
!	call subroutine point(Nmax,lon,lat,p,grad_p)
!		!ux=ux+(g(i)*cos(m*x)+h(i)*sin(m*x))*grad_p(i+1)
!	i=1
!	do n=1,N0
!		do m=0,n
!		NA(j,i)=cos(m*x)*grad_p(i+1)
!		NA(j,i+Nmax*(Nmax+3)/2)=sin(m*x)*grad_p(i+1)
!		i=i+1
!		end do
!	end do
!	end subroutine


	subroutine getgh(N0,g,h)!to calculate g,h
	integer N0
	real*8 g(N0*(N0+3)/2)
	real*8 h(N0*(N0+3)/2)
11	format(f12.2,f12.2)
	open(10,file='new.cof')
	do i=1,65
		read(10,11)g(i),h(i)		
	
	end do
	close(10)
	end subroutine getgh

	subroutine Utotal(N0,r,x,y,g,h,p,grad_p,u,ux,uy,uz)!to calculate u
	
	real*8 g(N0*(N0+3)/2)
	real*8 h(N0*(N0+3)/2)
	real*8 u,r,uz,ux,uy
	real*8 x
	real*8 y
	real*8 p(N0*(N0+3)/2+1) !to get p,we need p(0,0),so total number+1
	real*8 grad_p(N0*(N0+3)/2+1)!to fit with p,total number also plus one
	integer N0,n,m,i !we use n,m,i in loop
	real*8 date,alt,colat,elong,f
	isv=0;date=1900.0;itype=2;alt=6371.2;colat=60;elong=30
	u=0d0
	uy=0d0
	ux=0d0
	uz=0d0
	i=1
	do n=1,N0
		do m=0,n
		u=u+(g(i)*cos(m*x)+h(i)*sin(m*x))*p(i+1)	
		uy=uy+(g(i)*sin(m*x)-h(i)*cos(m*x))*m*p(i+1)
		uz=uz-(n+1)*(g(i)*cos(m*x)+h(i)*sin(m*x))*p(i+1)
		ux=ux+(g(i)*cos(m*x)+h(i)*sin(m*x))*grad_p(i+1)
		i=i+1
		
		end do
	end do
	u=u*r
	uy=uy/sin(y)
	end subroutine utotal
end
