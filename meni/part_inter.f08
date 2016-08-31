module var
	integer,parameter :: realKind=selected_real_kind(5),n=1000 !Definindo parametros como precisao do real e numero de simulacoes
	real(kind=realKind), dimension(2) :: x,y,dist,pot,delta !Vetores da simulacao, (1) é o valor real e (2) é o suposto valor exceto em delta
	real,parameter :: max_delta=2 !Valor maximo para uma tentativa

	real, dimension(n) :: result_table,err !Tabela que armazena o resultado final de cada simulacao e erro (desvio pad)
	real(kind=realKind) :: t1,t2 !Variaveis para comparar tempo
	integer :: criteria=10000
end module

subroutine run()
	use var
	implicit none
	real :: j
	integer :: i=0,t=0,k,g
	i=0
	t=0

	!Cria uma seed inteira aleatoria
	g = irand()
	call random_seed(g)
	!Chama valores aleatorios para x,y
	call random_number(x)
	call random_number(y)
	!Define os valores maximos para x e y
	x = x*10
	y = y*10

	do while( i < criteria )
		!Contadores
		t=t+1
		i=i+1

		call random_number(delta)
		delta = max_delta*(1-2*delta)

		!Adiciona a variacao aleatoria nas cordenadas futuras
		x(2) = x(1) + delta(1)
		y(2) = y(1) + delta(2)
		!Calcula a distancia
		dist(1) = sqrt(x(1)**2 + y(1)**2)
		dist(2) = sqrt(x(2)**2 + y(2)**2)
		!Calcula o potencial
		pot(1) = 2*((1/dist(1)**12) - (1/dist(1)**6))
		pot(2) = 2*((1/dist(2)**12) - (1/dist(2)**6))
		!Compara se o potencial melhorou
		if( pot(2)<pot(1) ) then
			x(1) = x(2)
			y(1) = y(2)
			i=0
		!Reseta o contador e usa as condicoes futuras como condicoes iniciais
		end if

		!De 500 em 500 gera uma seed nova para maior aleatoriedade
		if(mod(t,500) == 0) then
			call random_number(j)
			j = j * 400
			k = floor(j)
			call random_seed(k)
		end if

	end do
end subroutine

subroutine desvio_pad()
	use var
	implicit none
	integer :: l
	real(kind=realKind) :: s=0
	!Calcula o desvio padrao do array result_table
	do l=0,n
		s = result_table(l) + s !Soma a media
	end do
	s = s/(n+1) !Realiza a media
	write(*,"(A33,1X,F15.4)") 'Media final: ', s
	do l=0,n
		result_table(l) = abs(result_table(l) - s) !Pega a diferenca em relaçãoo a media
	end do
	result_table = result_table * result_table !Eleva ao quadrado
	do l=0,n
		s = result_table(l) + s !Soma denovo
	end do
	err = sqrt(s/(n+1))
	write(*,"(A34,1X,E15.3)") "Desvio Padrão: ", sqrt(s/(n+1)) !Tira a raiz da variância
end subroutine
subroutine time_write()
	use var

	write(*,"(A35,1X,F15.2,A)") "Operação realizada em: ",(t2-t1),"s"
	write(*,"(A35,1X,I15)") "Numero de operações: ", n
	write(*,"(A35,1X,E15.3,A)") "Tempo medio de cada simulação: " ,(t2-t1)/(n+1),"s"

end subroutine
program part_inter
	use var
	implicit none
	integer :: l,j
	open(1,file="results/sim_vals.dat")

	call CPU_TIME(t1) !Diferença de tempo
	do l=0,n
		!Roda a simulacao
		call run()
		!Escreve o resultado em um arquivo
		write(1,*) dist(1)
		!Salva o resultado para futuro calculo de erro
		result_table(l)= dist(1)
	end do
	call CPU_TIME(t2) !Diferença de tempo

	call time_write()
	!Calcula o desvio padrão e a media final
	call desvio_pad()

end program