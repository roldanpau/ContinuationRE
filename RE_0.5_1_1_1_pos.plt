#set out "RE_0.5_1_1_1_pos.eps"
#set term post eps color

set multiplot

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

set label "  m1" at -0.00000000000000,0.56884869024927 left
set label "m2  " at 0.42801097159491,0.08939696006771 right
set label "  m3" at 0.00000000000000,-0.46321826526006 left
set label "  m4" at -0.42801097159491,0.08939696006771 left

set object 1 ellipse center 0.42801097159491,0.08939696006771 \
size 0.2,0.1 behind

set arrow 5 from 0.52801097159491,0.08939696006771 \
to 0.7,0.08939696006771 front lt 3

set xrange [-0.7:1.8]
#set yrange [-0.6:0.6]

plot \
"RE_0.5_1_1_1_pos.res" u 2:3:1 not palette z pt 1, \
"RE_0.5_1_1_1_pos.res" u 4:5:1 not palette z pt 1, \
"RE_0.5_1_1_1_pos.res" u 6:7:1 not palette z pt 1, \
"RE_0.5_1_1_1_pos.res" u 8:9:1 not palette z pt 1

set origin .53, .3
set size .30, .30 
clear
unset key
set grid
#unset arrow 1
#unset arrow 2
unset arrow 3
unset arrow 4
unset arrow 5
unset object 1

set xtics 0.05
set ytics 0.03

set arrow 1 from 0.375,0.15,0 \
to 0.42801097159491,0.08939696006771,0 nohead palette z

set arrow 2 from 0.42801097159491,0.08939696006771,0 \
to 0.4,0.05 nohead palette z

#set bmargin 1.5
#set tmargin 1.5
#set lmargin 1.5
#set rmargin 1.5

unset xlabel
unset ylabel
unset colorbox

set xrange [0.35:0.5]
set yrange [0.05:0.15]

plot \
"RE_0.5_1_1_1_pos.res" u 8:9:1 not palette z pt 1

unset multiplot

set term pop
set out
