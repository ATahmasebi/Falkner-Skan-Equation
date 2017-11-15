function y = coupled_equations(eta,f) 
% f() is a vector containing f and its derivatives. S s=is also in f form:  
% f = f(1)  f' = f(2) f" = f(3) S=f(4)  S'=f(5)
% y() is a vector in which derivatives of f() are assigned 
y = zeros(5,1); 
y(1) = f(2);           
y(2) = f(3); 
y(3) = -f(1)*f(3)-0.1*((f(2))^2-1-f(4)); 
y(4) = f(5);
y(5) = -f(1)*f(5);
end
