module var
	integer,parameter :: realKind=selected_real_kind(5),n=1000000
	real(kind=realKind), dimension(2) :: x,y,dist,pot,delta
	real,parameter :: max_delta=2,min_delta=0.001

	real, dimension(n) :: random_table,result_table
	real(kind=realKind) :: t1,t2
	integer :: criteria=10000
end module

subroutine run()
	use var
	implicit none
	real :: j
	integer :: i=0,t=0,k
	i=0
	t=0
	do while( i < criteria )

		t=t+1

		call random_number(delta)
		delta = max_delta*(1-2*delta)


		x(2) = x(1) + delta(1)
		y(2) = y(2) + delta(2)

		dist(1) = sqrt(x(1)**2 + y(1)**2)
		dist(2) = sqrt(x(2)**2 + y(2)**2)
		
		pot(1) = 2*((1/dist(1)**12) - (1/dist(1)**6))
		pot(2) = 2*((1/dist(2)**12) - (1/dist(2)**6))
		
		i=i+1

		if( (pot(2)-pot(1))<0 ) then
			x(1) = x(2)
			y(1) = y(2)
			i=0
		end if
		

		if(mod(t,500) == 0) then
			call random_number(j)
			j = j * 400
			k = floor(j)
			call random_seed(k)
		end if


	end do
end subroutine


program part_inter
	use var
	implicit none
	integer :: l,j
	real(kind=realKind) :: s
	open(1,file="results/potXtemp.dat")
	open(2,file="results/time.dat")

	open(3,file="results/nTempts.dat")


	do l=0,n
		j = irand()
		call random_seed(j)

		call random_number(x)
		call random_number(y)

		x = x*10
		y = y*10

		call run()
		
		write(1,*) dist(1)
		result_table(l)= dist(1)

	end do
	do l=0,n
		s = result_table(l) + s
	end do
	s = s/n
	write(*,*) 'Media final:', s
	do l=0,n
		result_table(l) = abs(result_table(l) - s)
	end do
	result_table = result_table * result_table
	do l=0,n
		s = result_table(l) + s
	end do
	write(*,*) "Desvio Pad:", sqrt(s)

end program