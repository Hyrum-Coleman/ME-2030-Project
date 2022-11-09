%% Initialize workspace
clear all, close all, clc
% This code, along with funcBlock.m, simulates the frictional sliding of a 
% block down a ramp and around a circular (slotted) loop. The loop is said
% to be slotted because normal force between the block and the loop is 
% allowed to take on a negative value, instead of just disappearing (going
% to zero), as would happen in an actual open loop.

%% Define project constants and initial conditions
global R mu %% Define R and mu as global variables
R = 5 * 0.0254; %% semi-circle radius [in], converted to [m]
angRamp = 50; %% [deg], ramp angle relative to horizontal
mu = .1; %% coefficent of friction
mass = 1; %% mass of block [kg]; this is arbitrary in this problem
H = .6061569;

t = 0:0.001:2; %% time paratmeter (start time:time step:end time) for solving
     % ODE [s]; you may want to adjust this to shorten or lengthen the simulation

iters = 0;
flag = true;
normFlag = false;
angFlag = false;

while flag
    % Determine velocity of block when it enters loop
    angInit = 90 - angRamp; %% initial loop angle, rel to pos-x [deg]
    angInitRad = angInit*pi/180; %% convert deg to rad
    hLoop = H - R*(1-sin(angInitRad)); %% height change between drop and loop entry
    vLoop = sqrt(2.*9.81.*hLoop.*(1-mu.*tan(angInitRad))); %% velocity entering loop, 
         % from work-energy
    sLoop = R*angInitRad; %% initial loop position (defined from initial angle) [m]
    
    % Solve ODE
    % organize initial conditions into array
    y0 = [sLoop vLoop]; %% [<initial pos> <initial velocity>]
    
    % solve ODE at each specified time
    [t,y] = ode45(@funcBlock,t,y0); %%solve ODE defined in "funcBlock"
    
    % Evaluate results
    % relabel output
    pos = y(:,1); % pos [m]
    vel = y(:,2); % velocity [m/s]
    ang = pos/R; % loop angle [rad]
    angDeg = ang*180/pi;
    fNorm = mass*(9.81*sin(ang)+vel.^2/R); % normal force
    index = find((angDeg >= 270 - .1) & (angDeg <= 270 + .1));
    fNorm(index)
    avg = mean(fNorm(index))
    if (avg >= 0 - .001) && (avg <= 0 + .001)
        normFlag = true
    end

    if any(angDeg > 270)
        angFlag = true
    end

    if angFlag && normFlag
        flag = false;
    end
    H = H + .00000001;
    iters = iters + 1
end
fprintf("The minimum height is %f inches\n", H * 39.37)
%angAdjust = -(90 - pos/R*180/pi); % redefines angle to be zero at bottom of loop

% prepare figure 1 (angle and velocity as a function of time)
figure(1);
hold on;
grid on; box on;
xlabel('Time (s)')

% plot position results
yyaxis left
posPlot = plot(t,angDeg,'-','LineWidth',2);
ylabel("Ang (deg)");

% plot velocity results
yyaxis right
velPlot = plot(t,vel,'-','LineWidth',2);
ylabel("Vel (m/s)");

legend([posPlot,velPlot],'Ang','Vel');
hold off;

% prepare figure 2 (normal force and velocity as a function of loop angle,
     % where the angle is zero at the bottom of the loop)
figure(2);
hold on;
grid on; box on;
xlabel('Angle (deg)')

% plot velocity results
yyaxis right
velPlot = plot(angDeg, vel,'-','LineWidth',2);
ylabel("Vel (m/s)");

% plot normal force results
yyaxis left
fNormPlot = plot(angDeg,fNorm,'-','LineWidth',2);
ylabel("Force (N)");

legend([fNormPlot,velPlot],'Normal Force','Vel');
hold off;
