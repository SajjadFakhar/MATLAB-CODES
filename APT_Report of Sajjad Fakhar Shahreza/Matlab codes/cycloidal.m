function res=cycloidal(t,t3,S0,dS)
% cycloidal motion curve
res.pos=S0+dS*(t/t3-sin(2*pi*t/t3)/(2*pi));
res.vel=dS*(1-cos(2*pi*t/t3))/t3;
res.acc=2*pi*dS/t3^2*sin(2*pi*t/t3);
end