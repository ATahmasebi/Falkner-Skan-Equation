clear

global beta
beta = 1;
right_hand_endpoint =5;
stopping_distance = 0.0001;

shoot_parameter(1) = 1;
initial_value = [0,0,shoot_parameter];
[x_out, y_out] = ode45(@falkner_skan,[0, right_hand_endpoint],initial_value);
f_dash_endpoint(1) = y_out(end,2);

shoot_parameter(2) = 1.1;

n = 2;
while abs(shoot_parameter(n) - shoot_parameter(n-1)) > stopping_distance 
    initial_value = [0,0,shoot_parameter(n)];
    [x_out, y_out] = ode45(@falkner_skan,[0, right_hand_endpoint],initial_value);
    f_dash_endpoint(n) = y_out(end,2);
    plot(x_out,y_out','LineWidth',2)
    grid on
    xlabel('Time') % x-axis label
    ylabel('') % y-axis label
    LG=legend('f','f1','f2')
    pause(1) % pausing execution
    shoot_parameter(n+1) = shoot_parameter(n) - (f_dash_endpoint(n)-1)*(shoot_parameter(n)-shoot_parameter(n-1))/(f_dash_endpoint(n)-f_dash_endpoint(n-1));
    n = n+1;
end
fprintf('Done')