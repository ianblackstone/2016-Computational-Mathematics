opt = odeset('OutputFcn',@dblPendPlot);
[t,u] = ode45(@dblPend,[0 100],[0 pi-0.0001 0 0],opt);
