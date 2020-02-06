function dy = MyParticleTrajectory(tpart,y,C,A,c1,c2,c3,Bt,xb,omega_c)

dy = [y(2)
    omega_c*(y(1)*y(4)*Bz(y(1),y(5),A,C,c1,c2,c3)...
    -y(6)*Bphi(y(1),y(5),A,C,c1,c2,c3,Bt,xb))+y(1)*y(4)^2
    y(4)
    omega_c/y(1)*(y(6)*Br(y(1),y(5),A,C,c1,c2,c3)...
    -y(2)*Bz(y(1),y(5),A,C,c1,c2,c3))-2*y(2)/y(1)*y(4)
    y(6)
    omega_c*(y(2)*Bphi(y(1),y(5),A,C,c1,c2,c3,Bt,xb)...
    -y(1)*y(4)*Br(y(1),y(5),A,C,c1,c2,c3))];
