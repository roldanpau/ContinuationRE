#set out "RE_1_2_3.eps"
#set term post eps color

set size square

set grid

set title "Lagrange's equilateral RE with masses 1,2,3\n (color code = curvature K)"
set ylabel "y"
set xlabel "x"

plot \
"RE_1_2_3.res" u 2:3:1 not palette z, \
"RE_1_2_3.res" u 4:5:1 not palette z, \
"RE_1_2_3.res" u 6:7:1 not palette z

set term pop
set out
