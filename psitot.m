function z=psitot(x,y,A,C,c1,c2,c3)
z=C*x.^4/8+A*(x.^2.*log(x)/2-x.^4/8)+c1+c2*x.^2+c3*(x.^4-4*x.^2.*y.^2);