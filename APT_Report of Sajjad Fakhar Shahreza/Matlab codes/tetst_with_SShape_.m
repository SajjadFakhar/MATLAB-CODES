clear all; clc;close all;

l1=100;l2=70;
L=[l1;l2];

t11=1;t21=4;t31=5;%t1,2,3(first)motion curve
t12=1;t22=4;t32=5;%t1,2,3(second)motion curve
T=t31;
Fs(1) = 10;
Fs(2)= 10;

S1=[30;0];
S2=[170;0];
dS=S2-S1;
i=1; % conunter in For loop

%[MAP_Axis, MIP_Axis ,M_I, M_I_d, D_I,G_I,detJ]
% Initialize matrices to collect delta values
M_I_matrix = [];
M_I_d_matrix = [];
D_I_matrix = [];
G_I_matrix = [];
detJ_matrix = [];
MAP_Axis_matrix=[];
MIP_Axis_matrix=[];

%% Achive the position, velocity and acceleration with respect to the motion curve
for t=0:0.01:T
    resx=Sshape(t,S1(1),dS(1),t11,t21,t31);
    resy=Sshape(t,S1(2),dS(2),t12,t22,t32);
    %% position
    px(i)=resx.pos; py(i)=resy.pos;
    
    %% velocity
    vx(i)=resx.vel; vy(i)=resy.vel; 
    
    %% acceleration
    ax(i)=resx.acc; ay(i)=resy.acc; 
    
    %% storing the data in matrices 
    SS=[px(i); py(i)]; %end effector positipon
    SSp=[vx(i),vy(i)]; %end effector velocity
    SSa=[ax(i),ay(i)]; %end effector acceleration
    Q1=SCARAinv(SS,L,1);
    q1(i)=Q1(1); %alpha
    q2(i)=Q1(2); %beta
    J = SCARAjac(Q1,L); 
    Q1p = J^(-1).*SSp; %(Q' = j^-1  * S') , S' = [X',Y']TRANSPOSE (VELOCITY OF ENDEFFECTOR)
    q1p(i) =Q1p(1); %alpha'
    q2p(i) =Q1p(2); %beta'
    Fq(i,:)=-(J')*Fs';%this line allows us to solve the kinetik static problem
    grid on;
    PlotScara(Q1,L,'r',1); % here it will show and give the movement of joints to the final point
    hold on;
    
    [MAP_Axis, MIP_Axis ,M_I, M_I_d, D_I,G_I,detJ] = plotEllipsoid(L, SS);
    % Collecting the delta values in matrices
    M_I_matrix = [M_I_matrix; M_I];
    M_I_d_matrix = [M_I_d_matrix; M_I_d];
    D_I_matrix = [D_I_matrix ; D_I];
    G_I_matrix = [ G_I_matrix; G_I];
    detJ_matrix = [detJ_matrix ;detJ];
    MAP_Axis_matrix=[MAP_Axis_matrix;MAP_Axis];
    MIP_Axis_matrix=[MIP_Axis_matrix;MIP_Axis];
    hold off;
    
    time(i)=t;
    i=i+1;
end
%% Plotting the data that we achieved based on the motion curve 
QT = [q1(T*100),q2(T*100)];
PlotAreaSCARA(L,1);
hold on
PlotScara(QT,L,'r',1);grid; 
hold on
plotEllipsoid(L,S2)
hold on;
plot(px,py,'LineWidth',1,'Color','b');grid; % this line will show the linear movement between two joints and it shows on this plot PlotScara(Q1,L,'b',1);grid;
title('work space of this scara robot')
xlabel('x') 
ylabel('y')
hold off
figure;plot(px,py,'LineWidth',2,'Color','b');grid;
xlabel('-2\pi < x < 2\pi') 
ylabel('Sine and Cosine Values')
figure;plot(q1,q2);grid;
figure;
subplot(2,2,1);plot(time,vx);grid;
subplot(2,2,2);plot(time,vy);grid;
subplot(2,2,3);plot(time,q1p,time(1:end-1),diff(q1)/0.01);grid;
subplot(2,2,4);plot(time,q2p);grid;
figure;
plot(time,Fq(:,1),time,Fq(:,2));grid;
%% Plotting what we achieved for the Maipulability Ellipsoids
figure;
subplot(3,2,1);plot(time,M_I_matrix);grid on;
xlabel('time') ;
title('Manipulability Index (Area of the ellipse)');
subplot(3,2,2);plot(time,M_I_d_matrix,'color','r');grid on;
xlabel('time') ;
title('Manipulability Index (Area of the ellipse)');
subplot(3,2,4);plot(time,D_I_matrix);grid on;
xlabel('time') ;
title('Dextrity Index');
subplot(3,2,3);plot(time,detJ_matrix);grid on;
xlabel('time') ;
ylabel(' Determinant of J ');
title('Determinant of J');
subplot(3,2,5);plot(time,G_I_matrix);grid on;
xlabel('time') ;
ylabel(' GII ');
title('Global Isotropic Index or Condition Number');
figure();grid on;
plot(time,M_I_matrix);grid on;
hold on;
plot(time,M_I_d_matrix,'color','r');grid on;
xlabel('time') ;
ylabel('Manipulability Index (Area of the ellipse)');
hold off;
figure();grid on;
plot(time,M_I_matrix);grid on;
xlabel('time') ;
ylabel(' W  ');
title('Manipulability Index (Area of the ellipse)')
hold off;
figure();grid on;
plot(time,D_I_matrix,'color','g');grid on;
hold on;
title('Dextrity Index');
hold off;
figure();
plot(time,detJ_matrix,'color','g');grid on;
xlabel('time') ;
ylabel(' Determinant of J ');
figure();
plot(time,G_I_matrix);grid on;
xlabel('time') ;
ylabel(' GII ');
title('Global Isotropic Index or Condition Number');
figure();
plot(time,MAP_Axis_matrix,'color','g');grid on;
hold on
plot(time,MIP_Axis_matrix,'color','r');grid on;
xlabel('time') ;
title(' Length of principle Axis ');
if D_I == inf 
    disp('Dextrity index is equal to INFINITY so the robot is in a singularity Position')
end