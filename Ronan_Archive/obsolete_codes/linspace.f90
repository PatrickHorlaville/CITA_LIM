program linspace 

integer :: n=10000
real, dimension(10000) :: rcens
real, dimension(10001) :: redge

call linspc(rcens, redge) 

print *,rcens

contains
subroutine linspc(rcens, redge)
 real, dimension(10001) :: redge
 real, dimension(10000) :: rcens
 real a
 integer i
 real s
 real b 
 a = 5994.489132365306 ! Initial value
 b = 7089.244513924822 ! end value
 s = (b-a)/n
 redge = (/((a + i*s),i=0,n+1)/) ! assign to y in increments of 1.5 starting at 1.5
 rcens = (/((a + (s/2) + i*s),i=0,n)/) ! assign to y in increments of 1.5 starting at 1.5
end subroutine linspc

end program linspace
