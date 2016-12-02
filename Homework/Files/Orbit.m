function f = Orbit(t,pos)

mu = 0.012277471;
mus = 1 - mu;

D1 = ((pos(1)+mu)^2 + pos(2)^2)^(3/2);
D2 = ((pos(1)-mus)^2 + pos(2)^2)^(3/2);

f = [pos(3);pos(4);pos(1)+2*pos(4)-mus*(pos(1)+mu)/D1-mu*(pos(1)-mus)/D2;pos(2)-2*pos(3)-mus*pos(2)/D1-mu*pos(2)/D2];

end