#set out "RE_cont.eps"
#set term post eps color

set size square

set title "Continuation of Lagrange's equilateral RE \n(color code = curvature K)"
set ylabel "y"
set xlabel "x"

plot \
"Lagr.res" u 2:3:1 not palette z, \
"Lagr.res" u 4:5:1 not palette z, \
"Lagr.res" u 6:7:1 not palette z

set term pop
set out
