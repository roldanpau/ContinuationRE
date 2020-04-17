#set out "L_masses_1.eps"
#set term post eps 

set title "L as a function of radius, equal masses m=sqrt(3)"
set ylabel "L(u)"
set xlabel "radius"
plot "Lagr.res" notitle w linespoints

set term pop
set out
