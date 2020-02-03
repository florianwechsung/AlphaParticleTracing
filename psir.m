function z=psir(x,y,A,C,c1,c2,c3)
z=C*x.^3/2+A*(x.*log(x)+x/2-x.^3/2)+2*c2*x+4*c3*(x.^3-2*x.*y.^2);