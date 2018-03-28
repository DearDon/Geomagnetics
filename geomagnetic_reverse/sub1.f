	module sub1
	contains 
	subroutine cg(A,b,x,N,e)
c	variable	meaning
c	A		a matrix,from main program,is coefficient matrix of Ax=b
c	b		a vector,from main program,is righthand of Ax=b
c	x		a vector,the answer of Ax=b,is what we need,our goal
c	r		a vector,minus grads of 0.5xAx-bx at x point,says b-Ax
c	p		a vector,the direction of iteration better than r
c	w		a vector,value is A*p,is useful to simplify the process
c	q0		a number,value is r0*r0,is standard of loop times
c	q1		a number,value is rk-1*rk-1
c	q2		a number, value is rk*rk
c	ba,ar		a number,named by their pronounciation
c	e		a number,standard of loop times,input by client
c	test		a number,value is matmul(r,w)
c	pw		a number,value is matmul(p,w)
c	i		a number,count variable
c	N		a number,the degree of A
	real*8 A(N,N)
	real*8 b(N),x(N),r(N),p(N),w(N)
c	real*8 A(2,2),b(2),x(2),r(2),p(2),w(2)
	real*8 q0,q1,q2,ba,ar,e,test,pw
	integer i,N
c	write(*,*)"you want the x_error less than?"
c	read(*,*)e
	x=0.0d0
	r=b
	call onedimenmul(r,r,N,q0)
	q2=q0
	p=r

c	w=matmul(A,p)
c	call onedimenmul(r,w,2,test)
c	ar=q2/test
c	x=x+ar*p
c	r=r-ar*w
c	q1=q2;call onedimenmul(r,r,2,q2)
c	r=r-a*w
	i=1
	do while(q2>=q0*e)
	q1=q2
c	ba=q2/q1
c	p=r+ba*p
	
	w=matmul(A,p)
	call onedimenmul(p,w,N,pw)
c	pw is p*w
	ar=q1/pw
	x=x+ar*p
	r=r-ar*w
	call onedimenmul(r,r,N,q2)
c	r=r-a*w
	ba=q2/q1
	i=i+1
	p=r+ba*p
	end do
	end subroutine cg

c	This subroutine is to solve one dimention's multiplication
	subroutine onedimenmul(m1,m2,n,ans)
	integer n
	real*8 m1(n),m2(n),ans
	ans=0
	do i=1,n
	ans=m1(i)*m2(i)+ans
	end do
	
	end subroutine onedimenmul

	end module sub1

