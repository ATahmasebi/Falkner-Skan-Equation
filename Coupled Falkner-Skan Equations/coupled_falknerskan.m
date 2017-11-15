%Col1--Col--Col3--Col4--Col5
% f     f'   f''   S     S'
clear;clf;clc
tolerance=1e-6;
SW=1 ;% this is Sw, for boundary condition S(0)=S_W
eta_infinity=8;
eta_span = 0:0.01:eta_infinity;

est_f(1)=0.3; %For falkner equation
est_s(1)=0.5; % For S equation

% 1st iteration
[eta,f]=ode45(@BlasiusFunction,eta_span,[0 0 est_f(1) SW est_s(1)]);
targ_f(1)=f(length(f),2);% last value of column 2 (f')
targ_s(1)=f(length(f),4);% last value of column 4 (S)

% 2nd iteration
est_f(2)=0.7;
est_s(2)=0.6;

[eta,f]=ode45(@BlasiusFunction,eta_span,[0 0 est_f(2) SW est_s(2)]);
targ_f(2)=f(length(f),2);% last value of column 2 f'
targ_s(2)=f(length(f),4);% last value of column 4 s'

% Work parameters for bisection
A=(targ_f(1)-1)/(1-targ_f(2)); % for boundary condition f'(infinity)=1
B=(targ_s(1)-0)/(0-targ_s(2)); % for boundary condition S(infinity)=0

% Estimate new wall shear
estimate1=(est_f(1)+A*est_f(2))/(1+A);
estimate2=(est_s(1)+B*est_s(2))/(1+B);
[eta,f]=ode45(@BlasiusFunction,eta_span,[0 0 estimate1 SW estimate2]);

% Perform iterations until wall shear converges below a prescribed tolerance
it=0;
while abs(targ_f(2)-1)>tolerance && abs(targ_s(2)-0)>tolerance
    
    it=it+1;
    est_f(1)=est_f(2); 
    est_f(2)=estimate1;
    est_s(1)=est_s(2); 
    est_s(2)=estimate2; 
    
    targ_f(1)=targ_f(2);
    targ_f(2)=f(length(f),2);
    targ_s(1)=targ_s(2);
    targ_s(2)=f(length(f),4);
    
    A=(targ_f(1)-1)/(1-targ_f(2));
    B=(targ_s(1)-0)/(0-targ_s(2));
    
    estimate1=(est_f(1)+A*est_f(2))/(1+A); % for falkner function boundary condition
    estimate2=(est_s(1)+B*est_s(2))/(1+B); % estimate S function boundary condition
    [eta,f]=ode45(@BlasiusFunction,eta_span,[0 0 estimate1 SW estimate2]);
end

figure(1)
plot(eta,f(:,1),'k','linewidth',2);grid, ylabel('f'); xlabel('eta');
figure(2)
plot(eta,f(:,2),'k','linewidth',2);grid, ylabel('f^1'); xlabel('eta');
figure(3)
plot(eta,f(:,3),'k','linewidth',2);grid, ylabel('f^2'); xlabel('eta')
figure(4)
plot(eta,f(:,4),'k','linewidth',2);grid, ylabel('S'); xlabel('eta')
figure(5)
plot(eta,f(:,5),'k','linewidth',2);grid, ylabel('S^1'); xlabel('eta')

% Print wall shear
W1 = sprintf('Wall Shear is %d.',f(1,3));
W2 = sprintf('Wall Shear is %d.',f(1,4));
disp(W1)
disp(W2)

X = sprintf('F''(inf) tends to %d.',round(f(length(f),2)));
Y = sprintf('S(inf) tends to %d.',round(f(length(f),4)));
disp(X)
disp(Y)

