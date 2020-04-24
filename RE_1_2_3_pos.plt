#set out "RE_1_2_3_pos.eps"
#set term post eps color

set multiplot

set size ratio -1

set grid

#set title "Lagrange's equilateral RE with masses 1,2,3\n (color code = curvature K)"
set ylabel "y"
set xlabel "x"

set arrow 1 from 0.65569691463853,0.07571335803467,0 \
to -0.13113938292771,0.52999350624271,0 nohead palette z

set arrow 2 from -0.13113938292771,0.52999350624271,0 \
to -0.13113938292771,-0.37856679017336,0 nohead palette z

set arrow 3 from -0.13113938292771,-0.37856679017336,0 \
to 0.65569691463853,0.07571335803467,0 nohead palette z

set label "m1  " at 0.65569691463853,0.07571335803467 right
set label "  m2" at -0.13113938292771,0.52999350624271 left
set label "  m3" at -0.13113938292771,-0.37856679017336 left

set object 1 ellipse center -0.13113938292771,-0.37856679017336 \
size 0.2,0.1 behind

set arrow 4 from -0.02,-0.37856679017336 \
to 0.3,-0.37856679017336 front lt 3

set xrange [-0.3:1.4]
set yrange [-0.8:]

plot \
"RE_1_2_3_pos.res" u 2:3:1 not palette z pt 1, \
"RE_1_2_3_pos.res" u 4:5:1 not palette z pt 1, \
"RE_1_2_3_pos.res" u 6:7:1 not palette z pt 1

set origin .4, .17
set size .30, .30 
clear
unset key
set grid
unset arrow 1
#unset arrow 2
#unset arrow 3
unset arrow 4
unset object 1

set xtics 0.04
set ytics 0.02

set arrow 2 from -0.13113938292771,-0.36,0 \
to -0.13113938292771,-0.37856679017336,0 nohead palette z

set arrow 3 from -0.13113938292771,-0.37856679017336,0 \
to -0.10,-0.36,0 nohead palette z

#set bmargin 1.5
#set tmargin 1.5
#set lmargin 1.5
#set rmargin 1.5

unset xlabel
unset ylabel
unset colorbox

set xrange [-0.15:-0.1]
set yrange [-0.4:-0.36]

plot \
"RE_1_2_3_pos.res" u 6:7:1 not palette z pt 1

unset multiplot


set term pop
unset out
