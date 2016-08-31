module integration
    implicit none
    integer:: n = 3000000, s = 10
    contains
        function randArr(min, max, length) !função que retorna um vetor de números aleatórios entre dois números (min e max)
            implicit none
            real:: min, max
            real, dimension(length):: x, randArr
            integer:: length

            call random_number(x)
            randArr = ((max-min)*x)+min
        end function randArr

        real function rand(min, max) !função que retorna um número aleatório entre dois números (min e max)
            implicit none
            real:: min, max, x

            call random_number(x)
            rand = ((max-min)*x)+min
        end function rand

        real function getExtreme(f, xmin, xmax, ismax)!retorna o máximo ou o mínimo de uma função f(x), x->[xmin,xmax]
            implicit none
            real, external:: f
            real, intent(in):: xmin, xmax
            real::y, x, dlt, offset, xofmax, xmx, xmn
            integer:: d = 10000, i
            logical:: ismax

            xmn = xmin
            xmx = xmax
            x = xmn
            dlt = xmx - xmn
            getExtreme = f(x)
            xofmax = x
            do while(dlt > 0.000009)
                do i = 1, d
                    x = x + dlt/d
                    y = f(x)
                    if(ismax)then
                        if(y > getExtreme)then
                        getExtreme = y
                        xofmax = x
                        end if
                    else
                        if(y < getExtreme)then
                            getExtreme = y
                            xofmax = x
                        end if
                    end if
                end do
                if((xofmax - dlt/3) < xmn)then
                    offset = xmn - (xofmax - dlt/3)
                    xmx = (xofmax + dlt/3) + offset
                else if((xofmax + dlt/3) > xmx)then
                    offset = (xofmax + dlt/3) - xmx
                    xmn = (xofmax - dlt/3) - offset
                else
                    xmn = xofmax - dlt/3
                    xmx = xofmax + dlt/3
                end if
                dlt = xmx - xmn
                x = xmn
            end do
        end function getExtreme

        real function getRoot(f, xmin, xmax)!retorna uma raíz função f(x), x->[xmin,xmax]
            implicit none
            real, external:: f
            real, allocatable, dimension(:):: rts
            real:: xmin, xmax
            real::y, ly, x, dlt, offset
            integer:: d = 1000, i, j = 1

            x = xmin
            dlt = xmax - xmin
            ly = f(x)
            allocate(rts(j))
            do while(dlt > 0.0001)
                do i = 1, d
                    x = x + dlt/d
                    y = f(x)
                    if(abs(y) < abs(ly))then
                        ly = y
                        rts = x
                    end if
                    print *, dlt
                end do
                if((x - dlt/3) < xmin)then
                    offset = xmin - (x - dlt/3)
                    xmax = (x + dlt/3) + offset
                else if((x + dlt/3) > xmax)then
                    offset = (x + dlt/3) - xmax
                    xmin = (x - dlt/3) - offset
                else
                    xmin = x - dlt/3
                    xmax = x + dlt/3
                end if
                    dlt = xmax - xmin
                    x = xmin
            end do
        end function getRoot

        real function integrate(f, xmin, xmax)!integra uma função f(x), x->[xmin,xmax]
            implicit none
            real, intent(in):: xmin, xmax
            real:: ymax, ymin, switch, err
            real, dimension(n):: x, y
            real, external:: f
            real, dimension(s):: results
            integer :: i, j, counter

            ymax = getExtreme(f, xmin, xmax, .true.)
            ymin = getExtreme(f, xmin, xmax, .false.)

            if(ymax > 0 .and. ymin > 0)then
                ymin = 0
            else if(ymin < 0 .and. ymax < 0)then
                ymax = 0
            end if
            !print *, ymax, ymin, xmax, xmin

            x = randArr(xmin, xmax, n)
            y = randArr(ymin, ymax, n)

            do i = 1, s
                counter = 0.0
                do j = 1, n
                    if(f(x(j)) > 0.0 .and. f(x(j)) >= y(j) .and. y(j) >= 0.0)then
                        counter = counter + 1.0
                    else if(f(x(j)) < 0.0 .and. f(x(j)) <= y(j) .and. y(j) <= 0.0)then
                        counter = counter - 1.0
                    end if
                end do
                results(i) = (real(counter)/real(n))*((xmax-xmin)*(ymax-ymin))
                integrate = integrate + results(i)
            end do
            !print *,real(counter), real(n)
            integrate = integrate/real(s)

            !print *, integrate
            err = getError(integrate, results)
        end function

        real function integrateToInf(f, dlt)!integra uma função f(x), x->[xmin,inf)
            implicit none
            real, external:: f
            real, allocatable, dimension(:):: results, sr
            real:: dlt, xmin = 0.0, xmax, dif = 0.0, pinf = (1.0/1e-29), switch, rmd = 0
            integer:: i=1, j
            logical:: convergence = .false., divergence = .false.

            xmax = dlt
            !print *, pinf
            allocate(results(i))
            allocate(sr(i))
            do while(.not.convergence .and. .not.divergence)
                results(i) = integrate(f, xmin, xmax)
                if(xmax .ne. pinf)then
                    !print 1000, xmax, xmin
                    xmin = xmin + dlt
                    xmax = xmax + dlt
                else
                    divergence = .true.
                end if
                if(i >= 5)then
                    do j = 1, i-1
                        rmd = rmd+((results(i+1) - results(i))*(results(i+1) - results(i)))
                    end do
                    rmd =(rmd/(i-1))**(0.5)
                    print *, rmd
                    if(rmd < 1.0e-6)then
                        convergence = .true.
                    else if(integrateToInf == pinf)then
                        divergence = .true.
                    end if
                end if
                integrateToInf = integrateToInf + results(i)
                sr = results
                deallocate(results)
                allocate(results(i+1))
                results(:i) = sr(:i)
                deallocate(sr)
                allocate(sr(i+1))
                i = i+1
            end do
            if(abs(nint(integrateToInf) - integrateToInf) < 0.001)then
                integrateToInf = nint(integrateToInf)
            end if
            1000 format(/f10.6, 1x, f10.6/)
        end function integrateToInf

        real function getPI()!retorna PI
            implicit none
            real:: err = 0.0
            real, dimension(n):: x, y
            real, dimension(s):: results
            integer:: i, j, counter

            do i = 1, s
                x = randArr(-1.0, 1.0, n)
                y = randArr(-1.0, 1.0, n)
                counter = 0
                do j = 1, n
                    if((x(j)*x(j))+(y(j)*y(j)) <= 1)then
                        counter = counter + 1
                    end if
                end do
                results(i) = ((real(counter)/real(n))*4.0)
                getPI = getPI + results(i)
            end do
            getPI = getPI/real(s)

            do i = 1, s
                err = err+((results(i) - getPI)*(results(i) - getPI))
            end do
            err =(err/(s-1))**(0.5)

            print *,err
        end function getPI
        ! 3.1415915 0.00000506769

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

            1000 format(/i6,1x,i1,1x,a3/)
        end function getEuler

        real function square(x)!f(x) = x^2
            implicit none
            real:: x

            square = (x*x*x)-(6*(x*x))+(4*x)+12
        end function square

        real function seno(x)
            implicit none
            real:: x

            seno = sin(x)
        end function seno

        real function etonx(x)!f(x) = e^-x
            implicit none
            real:: x
            etonx = exp(-x)
        end function etonx

        real function getEqSpot(dlt)
            implicit none
            real, dimension(s):: results
            real:: x0, x, dlt, err, moves = 0.0, total = 0.0, dltx
            integer::i, j
            logical:: continuing

            getEqSpot = 0.0
            do i = 1, s
                x0 = rand(0.0, 10.0)
                continuing = .true.
                moves = 0.0
                total = 0.0
                do while(continuing)
                    dltx = rand((-1.0)*dlt, dlt)
                    x = x0 + dltx
                    if(potentialFunc(x) < potentialFunc(x0))then
                        x0 = x
                        moves = moves + 1.0
                    end if
                    total = total + 1.0
                    if(moves/total < 0.05 .and. total >= 100.0)then
                        continuing = .false.
                    end if
                end do
                results(i) = x0
                getEqSpot = getEqSpot + x0
            end do
            getEqSpot = getEqSpot/real(s)
            err = getError(getEqSpot, results)
            write(1,*) s, getEqSpot, err
        end function getEqSpot

        real function simParticleOscilation(x1)
            implicit none
            real:: x1, y, x2, eqx, x, dlt = 0.2, ctrl
            integer:: i
            logical:: continuing = .true.

            y = potentialFunc(x1)
            eqx = getEqSpot(0.5)
            x = x1

            if(x1 < eqx)then
                ctrl = 1.0
            else
                ctrl = 0.0
            end if

            do while(continuing)
                if(potentialFunc(x) <= y)then
                    x = x + rand((ctrl-1.0)*dlt, ctrl*dlt)
                else
                    x = x + rand(-(ctrl)*dlt, (-1.0*(ctrl-1.0))*dlt)
                end if

                if(abs(potentialFunc(x) - y) <= 1e-6)then
                    x2 = x
                    continuing = .false.
                end if
            end do

            dlt = (x2-x1)/20
            x = x1
            do i = 1, 40
                write(*,1000) x, potentialFunc(x)
                x = x+dlt
                if(i == 20)then
                    dlt = -1.0*dlt
                end if
            end do
            1000 format(f10.6, 1x, f10.6)
        end function simParticleOscilation

        real function potentialFunc(x)
            implicit none
            real:: x

            potentialFunc = 2.0*((1.0/(x)**12)- (1.0/(x)**6))
        end function potentialFunc

        real function getError(x, results)
            implicit none
            real:: x
            real, dimension(s):: results
            integer:: i
            do i = 1, s
                getError = getError + ((results(i) - x)**2)
            end do
            getError = (getError/(s))**(0.5)
        end function getError

        real function regErrorVariation()
            implicit none
            real:: x
            integer:: i, max = 990, dlt = 1, seed

            open(1, file='errorVariation.dat')
            do i=1, max
                seed = irand()
                call random_seed()
                x = getEqSpot(0.1)
                s = s+dlt
            end do
            close(1)
        end function regErrorVariation
end module integration
program monteCarloIntegration
    use integration
    implicit none
    integer:: i, max = 1000
    real:: x, eqSpot

    x = regErrorVariation()
end program