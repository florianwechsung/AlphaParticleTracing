function z=psirr(x,y,A,C,c1,c2,c3)
z=3/2*C*x.^2+A*(log(x)-3/2*x.^2+3/2)+2*c2+4*c3*(3*x.^2-2*y.^2);