function dy = MyFieldLines(phi,y,C,A,c1,c2,c3,Bt,xb)

dy = [-y(1)*psiz(y(1),y(2),A,C,c1,c2,c3)/sqrt(Bt^2*xb^2-2*A*psitot(y(1),y(2),A,C,c1,c2,c3))
    y(1)*psir(y(1),y(2),A,C,c1,c2,c3)/sqrt(Bt^2*xb^2-2*A*psitot(y(1),y(2),A,C,c1,c2,c3))];
