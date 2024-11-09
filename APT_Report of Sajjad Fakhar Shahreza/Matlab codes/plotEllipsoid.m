
function [MAP_Axis, MIP_Axis ,M_I, M_I_d, D_I,G_I,detJ] = plotEllipsoid(L, S)
%[MAP_Axis, MIP_Axis ,M_I, M_I_d, D_I,G_I,detJ] these values are taken to
%show the indices with respect to the time based on the motion curve.

%% Calculation of the Jacobian Matrix
Q = SCARAinv(S, L, 1);
J = SCARAjac(Q, L);
detJ = abs(det(J)); 

%% Compute Singular Values of J that are equal to Dimensions of the Axes
S_V = svd(J);
Maximum_SV = max(S_V);
Minimum_SV= min(S_V);

%% Compute JJT and its Eigenvalues and Eigenvectors
M = (J* J');
[EV, EI] = eig(M);

%% Plot the ellipsoid centered at the TCP position
ponit_angles = [0 : 0.01 : 2*pi] ; 
Ellipsoid_points_Force = EV * [Maximum_SV * cos(ponit_angles); Minimum_SV * sin(ponit_angles)];
Ellipsoid_points_Velocity = EV * [Minimum_SV * cos(ponit_angles);Maximum_SV * sin(ponit_angles)];

%% All the points that we achieved for ellipsoids must be located in the TCP(S = [X Y])
Ellipsoid_Velocity = Ellipsoid_points_Velocity + S; 
Force_Ellipsoid= Ellipsoid_points_Force + S;

hold on;
plot(Ellipsoid_Velocity(1, :), Ellipsoid_Velocity(2, :), 'k', 'LineWidth', 0.1);
plot(Force_Ellipsoid(1, :), Force_Ellipsoid(2, :), 'b', 'LineWidth', 0.1);

%% Plot eigenvectors with annotations
S_f = 50;  % Adjust the scale factor for better visualization
% S_f = 3; % This scale factor used for the time that links are equal to 1 cm
quiver(S(1,1), S(2,1), EV(1, 1) * S_f, EV(2, 1) * S_f, 'r', 'LineWidth', 1);
quiver(S(1,1), S(2,1), EV(1, 2) * S_f, EV(2, 2) * S_f, 'b', 'LineWidth', 1);

%% Calculating the indices of the Manipulability Ellipsoids
M_I = (sqrt((max(eig(M)))*(min(eig(M))))); %Manipulability Index = (sqrt(det(M))) Area of the ellipse
M_I_d = (sqrt(det(J*J')));% W in the time that the robot is in nonredundent situation( area of the ellipse)
D_I = (Maximum_SV/Minimum_SV);%Dextrity Index
G_I = min(S_V) / max(S_V);%Global Isotropy Index or Condition Number
MAP_Axis=Maximum_SV;
MIP_Axis=Minimum_SV;

end
