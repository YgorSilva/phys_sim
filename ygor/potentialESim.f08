module randModule
    implicit none
    contains
        real function rand(min, max)
            implicit none
            real:: min, max, x

            call random_number(x)
            rand = (max-min)*x + min
        end function
        function randArrUD(max, length) !função que retorna um vetor de números aleatórios entre dois números (min e max)
            implicit none
            real:: max
            real, dimension(length):: x, randArrUD
            integer:: length

            call random_number(x)
            randArrUD = max*(1-(2*x))
        end function randArrUD

        real function randUD(max) !função que retorna um número aleatório entre dois números (min e max)
            implicit none
            real:: max, x

            call random_number(x)
            randUD = max*(1-(2*x))
        end function randUD
end module randModule
program monteCarloIntegration
    use randModule
    implicit none
    integer:: s = 1
    real:: x, t = 0.3, e
    e = getEuler(1, 100)

    call regErrorVariation()
    contains
    subroutine getEqSpot(dlt)
            implicit none
            real, dimension(s):: results
            real:: x0, x, dlt, err, dltx, eqSpot
            integer::i, j, seed, nonMoves
            logical:: continuing

            eqSpot = 0.0
            results = 0.0
            seed = irand()
            call random_seed(seed)
            do i = 1, s
                nonMoves = 0
                x0 = rand(0.0, 10.0)
                do while(nonMoves < 10000)
                    dltx = randUD(dlt)
                    x = x0 + dltx
                    if(potentialFunc(x) < potentialFunc(x0))then
                        x0 = x
                        nonMoves = 0
                    else
                        nonMoves = nonMoves + 1
                    end if
                end do
                write(*,*) s, i, x0, potentialFunc(x0)
                results(i) = x0
                eqSpot = eqSpot + x0
            end do
            eqSpot = eqSpot/s

            !   do i = 1, s
            !end do
            err = getError(eqSpot, results)
            write(1,*) s, eqSpot, err
        end subroutine getEqSpot

        real function potentialFunc(x)
            implicit none
            real:: x

            potentialFunc = 2.0*((1.0/(x)**12)- (1.0/(x)**6))
        end function potentialFunc

        real function getError(md, samples)
            implicit none
            real:: md
            real, dimension(s):: samples
            integer:: i
            getError = 0.0
            do i = 1, s
                getError = getError + (abs(samples(i)-md)*abs(samples(i)-md))
            end do
            getError = sqrt((getError/(s-1)))
        end function getError

        subroutine regErrorVariation()
            implicit none
            integer:: k, max = 100, dltS = 1, seed
            real:: delta = 0.01

            open(1, file='errorVariation.dat')
            do k=1, max
                call getEqSpot(delta)
                s = s+dltS
            end do
            close(1)
        end subroutine regErrorVariation

        real function getNSProb(dltE)
            implicit none
            real:: dltE

            getNSProb= e**(-dltE/t)
        end function getNSProb

         recursive function getEuler(indx, mx) result(e) !retorna e
            implicit none
            integer, parameter::  r = selected_real_kind(24)
            real(kind=r):: e
            integer:: indx, x = 0, mx, nxtindx

            nxtindx = indx+1
            if(indx > mx)then
                e = 1
            else if(indx == 1)then
                e = 2.0 + 1.0/getEuler(nxtindx, mx)
            else if(indx == 2)then
                e = 1.0 + 1.0/getEuler(nxtindx, mx)
            elseif((mod(x,3) == 0).and.(indx > 2).and.(indx <= mx))then
                x = x+1
                e = ((real(x-1)+3.0)/3.0)*2.0 + 1.0/getEuler(nxtindx, mx)
            elseif((mod(x,3) /= 0).and.(indx > 2).and.(indx <= mx))then
                x = x+1
                e = 1.0 + 1.0/getEuler(nxtindx, mx)
            end if
        end function getEuler
end program