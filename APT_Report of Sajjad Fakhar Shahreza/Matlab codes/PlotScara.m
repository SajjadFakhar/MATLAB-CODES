function PlotScara (Q,L,colore,fig)
figure(fig)
q1=Q(1);
q2=Q(2);
l1=L(1);
l2=L(2);
x1=l1*cos(q1);
y1=l1*sin(q1);
x2=x1+l2*cos(q1+q2);
y2=y1+l2*sin(q1+q2);
plot([0 x1 x2],[0 y1 y2],'LineWidth',1.1,'color',colore);
axis equal;
app=(l1+l2 )*2.2;
xlim([-app app])
ylim([-app app]);
end
