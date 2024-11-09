function Jp = SCARAjacP(Q,Qp,L)
% time derivative of the Jacobian matrix: SCARA robot
l1=L(1);
l2=L(2);
Jp= zeros(2,2);
Jp(1,1)= -l1*cos(Q(1))*Qp(1)- l2*cos(Q(1)+Q(2))*(Qp(1)+Qp(2));
Jp(1,2)= -l2*cos(Q(1)+Q(2))*(Qp(1)+Qp(2));
Jp(2,1)= -l1*sin(Q(1))*Qp(1)-l2*sin(Q(1)+Q (2))*(Qp(1)+Qp(2));
Jp(2,2)= -l2*sin(Q(1)+Q(2))*(Qp(1)+Qp(2));
end
