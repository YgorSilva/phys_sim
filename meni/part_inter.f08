module var
	integer,parameter :: realKind=selected_real_kind(4),criteria=1000,n=10000
	real(kind=realKind), dimension(2) :: dist,pot
	real :: delta
	real,parameter :: max_delta=1

	real, dimension(n) :: random_table,result_table
end module

subroutine run()
	use var
	implicit none
	integer :: i=0,t=0
	i=0
	do
		t = t+1
		call random_number(delta)
		delta = max_delta*(1-2*delta)

		!write(*,*) dist(1)
		
		dist(2) = dist(1) + delta
		
		pot(1) = 2*((1/dist(1)**12) - (1/dist(1)**6))
		pot(2) = 2*((1/dist(2)**12) - (1/dist(2)**6))
		
		!write(*,*)t, dist(1)
		
		if( (pot(2)-pot(1))<0 ) then
			dist(1) = dist(2)
			i=0
		elseif(i > criteria) then
			exit
		end if
		i=i+1
	end do
end subroutine


program part_inter
	use var
	implicit none
	integer :: l,j
	real :: s
	open(1,file="results/potXtemp.dat")

	call random_number(random_table)

	do l=0,n
		j = irand()
		call random_seed(j)
		call random_number(dist)
		
		dist(1) = dist(1)*10

		call run()
		!write(*,*) dist(1)
		result_table(l)= dist(1)

	end do
	do l=0,n
		s = result_table(l) + s
	end do
	write(*,*) s/n

end program