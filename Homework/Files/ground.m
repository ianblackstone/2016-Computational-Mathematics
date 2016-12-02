function [gstop,isstop,direction] = ground(t,u)

gstop = u(2);
isstop = 1; %Flag that tells matlab wiether or not to stop the simulation
direction = -1;%telling matlab wheither or not is can approach to find t 
%values wether going right or left (negative or positive values)
end
