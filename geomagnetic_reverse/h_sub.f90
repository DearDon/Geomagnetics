
module sub
	contains

	subroutine lerangde(Nmax,lat,lon,p,grad_p)!lon is lamata,lat is theta
	real*8 lon
	real*8 lat
	real*8 p(0:Nmax,0:Nmax) 
	real*8 grad_p(1:Nmax,0:Nmax)
	integer Nmax,n,m,i,j

	
	p(0,0)=1
	p(1,0)=cos(lat)
	p(1,1)=sin(lat)
	print *, p(1,1)
	grad_p(1,0)=(cos(lat)*p(1,0)-p(0,0))/sin(lat)
	grad_p(1,1)=cos(lat)*p(1,1)/sin(lat)



	do i=0,Nmax
	  do j=0,Nmax
	    if(j<=i) then
     	  do n=2,Nmax
		   do m=0,n
			 if(m==n	)then
				p(m,m)=p(m-1,m-1)*sin(lat)*((2*m-1)/(2*m))**0.5
				grad_p(m,m)=n*cos(lat)*p(m,m)*((2*m-1)/(2*m))**0.5
			 else
				p(n,m)=((2*n-1)*cos(lat)*p(n-1,m)-((n-1)**2-m**2)** &
			 &	0.5*p(n-2,m))/(n**2-m**2)**0.5
				grad_p(n,m)=(n*cos(lat)*p(m,n)-(n**2-m**2)**&
			 &	0.5*p(n-1,m))/sin(y)
			 end if                            
			end do
		   end do
		 else
		   p(i,j)=0
		 end if
	  end do 
	end do			
	end subroutine lerangde

	
	subroutine cosine(Nmax,lon,lat,a,r,NA,p)!to calculate g,h
	real*8  lon
	real*8  lat
	real*8 a,r
	real*8 p(0:Nmax,0:Nmax) 
	real*8 :: NA(Nmax*(Nmax+2))
	integer i,m,n

	i=1
	 do n=1,Nmax 
	     do m=0,n
		   NA(i)=-cos(m*lon)*p(n,m)*(n+1)*(a/r)**(n+2)
		   i=i+1
		 end do
		 
	     do m=1,n
		   NA(i)=-sin(m*lon)*p(n,m)*(n+1)*(a/r)**(n+2)
		   i=i+1
		 end do
	 end do

	end subroutine cosine

end module sub
