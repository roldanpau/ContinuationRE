#set out "RE_0.5_1_1_1_pos.eps"
#set term post eps color

set size ratio -1

set grid

#set title "Continuation of Lagrange's equilateral RE \n(color code = curvature K)"
set ylabel "y"
set xlabel "x"

set arrow 1 from -0.00000000000000,0.56884869024927,0 \
to 0.42801097159491,0.08939696006771,0 nohead palette z

set arrow 2 from 0.42801097159491,0.08939696006771,0 \
to 0.00000000000000,-0.46321826526006,0 nohead palette z

set arrow 3 from 0.00000000000000,-0.46321826526006,0 \
to -0.42801097159491,0.08939696006771,0 nohead palette z

set arrow 4 from -0.42801097159491,0.08939696006771,0 \
to -0.00000000000000,0.56884869024927,0 nohead palette z

plot \
"RE_0.5_1_1_1_pos.res" u 2:3:1 not palette z pt 1, \
"RE_0.5_1_1_1_pos.res" u 4:5:1 not palette z pt 1, \
"RE_0.5_1_1_1_pos.res" u 6:7:1 not palette z pt 1, \
"RE_0.5_1_1_1_pos.res" u 8:9:1 not palette z pt 1

set term pop
set out
