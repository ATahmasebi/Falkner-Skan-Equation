function y_out = falkner_skan(x,y_in)

global beta

y_out(1,1) = y_in(2);
y_out(2,1) = y_in(3);
y_out(3,1) = -y_in(1)*y_in(3) - beta*(1-y_in(2)^2);

end