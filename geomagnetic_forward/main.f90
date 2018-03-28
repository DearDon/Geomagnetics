! PROGRAM:
! This is a Geomagnetic Field Modeling software by Don!
! this program calculates the geomagnetic field values from a spherical harmonic! model.The coefficients are all from IGRF model
! this program can get x,y,z,f of magnet and U only if we input latitude and 
! longitude .In this simple program we only compute the value on surface of 
! earth.
! HISTORY:
! 2011/12/5
! First release by Don





!This main program call three subroutines to finish the calculation
!variable claim
!in program   		real meaning
!g(N0*(N0+3)/2)		spherical harmonic coefficients
!h(N0*(N0+3)/2)		spherical harmonic coefficients
!p(N0*(N0+3)/2)		lerande function
!grad_p			the grad of p
!x			"longitude" range from 0 to 360
!y			"latitude"  range from 0 to 180
!u			U --u
!uz			vertically-downward component --z
!ux			northward component --x
!uy			eastward component --y
!T			total intensity --f
!r			earth's radius --a
!N0			the max n,in other word, degree


program main
use sub


integer,parameter::N0=10 !degree
real*8,parameter::r=6371200d0,pi=4d0*atan(1d0)!pi=3.14159265358979... !r is radii
real*8 g(N0*(N0+3)/2)
real*8 h(N0*(N0+3)/2)
real*8 p(N0*(N0+3)/2+1)
real*8 grad_p(N0*(N0+3)/2+1)
real*8 x,y,u,uz,ux,uy,T
!x y used to give the position

!write(*,*)"please write x,y(angel,x from -180 to 180,y from -90 to 90 )"
!read(*,*)x,y
do i=-179,179
do j=-89,89
x=i
y=j
y=90-y
x=x*pi/180d0
y=y*pi/180d0
call point(N0,x,y,p,grad_p)!to use x,y to calculate p and grad_p
call getgh(N0,g,h)!to get g and h from a file
call Utotal(N0,r,x,y,g,h,p,grad_p,u,ux,uy,uz)!to calculate magnet
T=(ux**2+uy**2+uz**2)**0.5 !to get f,total intensity
!write(*,*)"u",u
write(*,*)"x",i,j,ux
!write(*,*)"y",uy
!write(*,*)"z",uz
!write(*,*)"T",T
end do
end do
end program
