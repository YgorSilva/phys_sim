reset
set term gif animate delay 15
set output 'animate.gif'
set xrange[0:10]
set yrange[-2:1]
n = "`awk 'END {print NR}' < moves.dat`"
i=0; while i<n{plot 2*(1/x**12 - 1/x**6), 'moves.dat' every ::i::i using 1:2 with circles; i=i+1}
set output