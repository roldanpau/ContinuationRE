#set out "RE_1_2_3.eps"
#set term post eps color

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

plot \
"RE_1_2_3.res" u 2:3:1 not palette z pt 1, \
"RE_1_2_3.res" u 4:5:1 not palette z pt 1, \
"RE_1_2_3.res" u 6:7:1 not palette z pt 1

set term pop
set out
