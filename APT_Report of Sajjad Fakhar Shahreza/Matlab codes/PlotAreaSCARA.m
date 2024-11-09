
function PlotAreaSCARA(L,fig)
figure(fig); 
hold off; 
plot ([0],[0]);
hold on; 
grid on;
l1=L(1); l2=L(2);
d1=l1+l2; d2=abs(l1 -l2);
nn =50;
for i=0 : nn
fi=(2*pi/nn)*i;
x1=l1*cos(fi); y1=l1*sin(fi);
xa=d1*cos(fi); ya=d1*sin(fi);
xb=d2*cos(fi); yb=d2*sin(fi);
if i~=0
plot([x1old x1],[y1old y1],'LineWidth',1,'color','green');
plot([xaold xa],[yaold ya],'LineWidth',2,'color','green');
plot([xbold xb],[ybold yb],'LineWidth',2,'color','green');
end
x1old =x1; y1old =y1;
xaold =xa; yaold =ya;
xbold =xb; ybold =yb;
end
hold on
end
