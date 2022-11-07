function dydt = funcBlock(t,y)
%% function defines ODE governing block sliding around loop

global R mu

% Redefine input variables
s = y(1);
theta = s/R; % angle [rad] defined positive in clockwise direction from 
     % standard pos-x axis
v = y(2);

% Create output array and fill with zeros
dydt = zeros(size(y)); % don't change

dydt(1) = y(2); % don't change
dydt(2) = 9.81*cos(theta) - mu*9.81*sin(theta)*sign(v) - mu*v^2/R*sign(v); 
     % this last line is the expression for the tangential component of
     % acceleration in your problem