function status = dblPendPlot(t,u,flag)

if strcmp(flag,'done')
    status = 0;
    return;
end

th1 = u(1,end);
th2 = u(2,end);
x1 = 2*sin(th1);
y1 = -2*cos(th1);
x2 = 2*sin(th1) + sin(th2);
y2 = -(2*cos(th1) + cos(th2));
x = [0 x1 x2];
y = [0 y1 y2];
plot(x,y,'k-o');
xlim([-3 3]);
ylim([-3 3]);
daspect([1 1 1]);
pause(0.05)
drawnow;
status = 0;
