clear all;
close all; 
clc;
% robot parameters
% l1=100; %[cm]
% l2=100; %[cm]
l1=1; %[cm]
l2=1; %[cm]
g1=0.28; %[m]
g2=0.283; %[m]
L=[l1;l2;g1;g2]; %[m]
M=diag([0.5,0.5,7,7,0.0565,7.8,7.8,0.0565]);
Fse=zeros(8,1);
Fse(1)=1;
Fse(2)=1;

Q=[-0.244;2.375];
Qp=[-2.66;-2.2];
Qpp=[26.63;-12.36];

%position of the extended 2R(operational space)
Sd=SCARAdir_din(Q,L);
%Jacobians matrix and its time derivative
J=SCARAjacdin(Q,L);
Jp=SCARAjacPdin(Q,Qp,L);
% acceleration of the extended 2R
Spp=Jp*Qp+J*Qpp;
%inertial forces
Fsi=-M*Spp;
% forces in the operational space
Fs=(Fse+Fsi);
% forces in the joint space
Fcq=-J'*Fs;
display(Fcq);
S0=[0;0]; % gripper position for links equal to 1 cm
Sf=[ 2; 0];
% S0=[0;0]; % gripper position for links equal to 1 meter
% Sf=[ 200; 0];
dS=Sf-S0;
dT =1/100;
T=5;
tt =[0:dT:T];
n= length(tt);
t3 = 5;

% Initialize matrices to collect Manipulability Indices
M_I_matrix = [];
M_I_d_matrix = [];
D_I_matrix = [];
G_I_matrix = [];
detJ_matrix = [];
MAP_Axis_matrix=[];
MIP_Axis_matrix=[];
%% Achive the position, velocity and acceleration with respect to the motion curve
for i=1:n
 resx=cycloidal(tt(i),t3,S0(1),dS(1));
 resy=cycloidal(tt(i),t3,S0(2),dS(2));


S=[resx.pos;resy.pos];
Sp=[resx.vel;resy.vel];
Spp=[resx.acc;resy.acc];

Q=SCARAinv(S,L,1);
J=SCARAjac(Q,L);
Qp=inv(J)*Sp;
Jp=SCARAjacP(Q,Qp,L);
Qpp=inv(J)*(Spp-Jp*Qp);

Sd=SCARAdir_din(Q,L);
Je=SCARAjacdin(Q,L);
M = Je * Je';
[ev,ei]=eig(M);
Jep=SCARAjacPdin(Q,Qp,L);

Sepp=Jep*Qp+Je*Qpp;
Fsi=-M*Sepp;
Fs=(Fse+Fsi);
Fcq=-Je'*Fs;
Fq1(i) = Fcq(1);
Fq2(i) = Fcq(2);

if i == 1
    Q1s = Q;%storing the first value for ploting in the plot scara function
    L1s = L;
end

if i == n
    Q1f = Q;%storing the last value for ploting in the plot scara function
    L1f = L;
end
grid on;
PlotScara(Q,L,'k',1);

[MAP_Axis, MIP_Axis ,M_I, M_I_d, D_I,G_I,detJ] = plotEllipsoid(L, S);

 % Collecting the delta values in matrices
 M_I_matrix = [M_I_matrix; M_I];
 M_I_d_matrix = [M_I_d_matrix; M_I_d];
 D_I_matrix = [D_I_matrix ; D_I];
 G_I_matrix = [ G_I_matrix; G_I];
 detJ_matrix = [detJ_matrix ;detJ];
 MAP_Axis_matrix=[MAP_Axis_matrix;MAP_Axis];
 MIP_Axis_matrix=[MIP_Axis_matrix;MIP_Axis];

hold off;
end 
%% Plotting the first and the last position of the Scara and the Ellipsoid
figure();
hold on;
grid on;
Q1=SCARAinv(S,L ,1); % first solution (beta >0)
Q2=SCARAinv(S,L,-1); % second solution( beta <0)i
PlotAreaSCARA(L,1); % draw the work area
PlotScara(Q1s,L1s,'r',1);
PlotScara(Q1f,L1f,'b',1);
plotEllipsoid(L,S);
hold on
plot(Sf,S0);
hold off
figure(2); 
%grid;
plot(tt,Fq1,tt,Fq2);

%% Plotting the results related to the indices of the ellipsoids
figure;
subplot(4,1,1);plot(tt,M_I_matrix);grid on;
xlabel('time') ;
title('Manipulability Index')
subplot(4,1,2);plot(tt,M_I_d_matrix,'color','r');grid on;
xlabel('time') ;
title('Determinant of Manipulability Index')
subplot(4,1,3);plot(tt,detJ_matrix);grid on;
xlabel('time') ;
title('Determinant of J')
subplot(4,1,4);plot(tt,G_I_matrix);grid on;
xlabel('Cycloidal Motion Curve') ;
ylabel(' GII ');
title('Global Isotropic Index');
figure();grid on;
plot(tt,M_I_matrix);grid on;
hold on;
plot(tt,M_I_d_matrix,'color','r');grid on;
xlabel('time') ;
ylabel(' Manipulability Index ');
hold off;
figure();grid on;
plot(tt,M_I_matrix);grid on;
xlabel('time') ;
ylabel(' W  ');
title('Manipulability Index links equal to 1 cm')
hold off;
figure();grid on;
plot(tt,D_I_matrix,'color','g');grid on;
hold on;
xlabel('time') ;
ylabel(' K ');
title('Dextrity Index');
hold off;
figure();
plot(tt,detJ_matrix,'color','g');grid on;
xlabel('time') ;
ylabel(' Determinant of J ');
figure();
plot(tt,G_I_matrix,'color','g');grid on;
xlabel('time') ;
title(' Global Isotropic Index or Condition Number ');
figure();
plot(tt,MAP_Axis_matrix,'color','g');grid on;
hold on
plot(tt,MIP_Axis_matrix,'color','r');grid on;
xlabel('time') ;
ylabel(' Length of principle Axis ');
title('Length of principle Axis')
if D_I == inf
    disp('Dextrity index is equal to INFINITY so the robot is in a singularity Position')
end