#set out "RE_1_1_1_pos.eps"
#set term post eps color

set size ratio -1

set grid

#set title "Continuation of Lagrange's equilateral RE \n(color code = curvature K)"
set ylabel "y"
set xlabel "x"

set arrow 1 from 0.41634158882780,-0.00000000000000,0 \
to -0.20817079441390,0.36056239257685,0 nohead palette z

set arrow 2 from -0.20817079441390,0.36056239257685,0 \
to -0.20817079441390,-0.36056239257685,0 nohead palette z

set arrow 3 from -0.20817079441390,-0.36056239257685,0 \
to 0.41634158882780,-0.00000000000000,0 nohead palette z

plot \
"RE_1_1_1_pos.res" u 2:3:1 not palette z pt 1, \
"RE_1_1_1_pos.res" u 4:5:1 not palette z pt 1, \
"RE_1_1_1_pos.res" u 6:7:1 not palette z pt 1

set term pop
set out
