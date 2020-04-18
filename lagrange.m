clc
clear
format long

l = cbrt(6.0)/2.0;
m1 = 1;
m2 = 2;
m3 = 3;
m = m1+m2+m3;

x1 = (l/sqrt(3))*(1-m1/m) + (l/(2*sqrt(3)))*(m2/m) + (l/(2*sqrt(3)))*(m3/m)
y1 = (l/(2*m))*(m3-m2)

x2 = -(l*m1)/(sqrt(3)*m) - l*(1-m2/m)/(2*sqrt(3)) + (l*m3)/(2*sqrt(3)*m)
y2 = l*(1-m2/m)/2 + (l*m3)/(2*m)

x3 = -(m1*l)/(m*sqrt(3)) + (m2*l)/(2*m*sqrt(3)) - l*(1-m3/m)/(2*sqrt(3))
y3 = -(m2*l)/(2*m) - l*(1-m3/m)/2