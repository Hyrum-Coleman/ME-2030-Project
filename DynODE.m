%% Initialize workspace
clear
clc
close all
% This code, along with funcBlock.m, simulates the frictional sliding of a 
% block down a ramp and around a circular (slotted) loop. The loop is said
% to be slotted because normal force between the block and the loop is 
% allowed to take on a negative value, instead of just disappearing (going
% to zero), as would happen in an actual open loop.

%% Calculation for mu = .1
tic
global R mu %% Define R and mu as global variables
R = 5 * 0.0254; %% semi-circle radius [in], converted to [m]
angRamp = 50; %% [deg], ramp angle relative to horizontal
mu = .1; %% coefficent of friction
mass = 1; %% mass of block [kg]; this is arbitrary in this problem
 H = .60611;

t = 0:0.0001:2; %% time paratmeter (start time:time step:end time) for solving
     % ODE [s]; you may want to adjust this to shorten or lengthen the simulation

iters = 0;
flag = true;
normFlag = false;
angFlag = false;

fprintf("Starting loop for mu = %.2f\n", mu)
tic
while flag
    % Setup
    angInit = 90 - angRamp;
    angInitRad = angInit*pi/180;
    hLoop = H - R*(1-sin(angInitRad));
    vLoop = sqrt(2.*9.81.*hLoop.*(1-mu.*tan(angInitRad)));
    sLoop = R*angInitRad; %% initial loop position (defined from initial angle) [m]
    
    % Solve ODE
    y0 = [sLoop vLoop];
    [t,y] = ode45(@funcBlock,t,y0);
    
    % Evaluate results
    pos = y(:,1); % position
    vel = y(:,2); % velocity
    ang = pos/R;  % angles
    angDeg = ang*180/pi; % degree versions of angles
    fNorm = mass*(9.81*sin(ang)+vel.^2/R); % normal force

    % Check flags
    index = find((angDeg >= 270 - 1) & (angDeg <= 270 + 1)); % where the angle is ~ 270 (at the top of the loop)
    avg = mean(fNorm(index)); % average normal force at the top of the loop

    if any(angDeg > 270)
        angFlag = true; % Set to true if we have enough height to reach the top of the loop
    end
    if (avg >= 0 - .001) && (avg <= 0 + .001)
        normFlag = true; % Set to true if the average NF at the top is ~0
    end

    if angFlag && normFlag
        flag = false; % exit loop if other flags are true
    end

    % Iterate
    H = H + .00000001;
    iters = iters + 1;
end
fprintf("The minimum height for mu = %.2f is %f inches\n", mu, H * 39.37);
fprintf("It took %d iterations\n", iters);
toc
fprintf("\n")

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

% prepare figure 2
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

%% Calculation for mu = .2
mu = .2;
H = 1.16634;
iters = 0;
flag = true;
normFlag = false;
angFlag = false;

fprintf('Starting loop for mu = %.2f\n', mu)
tic
while flag
    % Setup 
    angInit = 90 - angRamp;
    angInitRad = angInit*pi/180;
    hLoop = H - R*(1-sin(angInitRad));
    vLoop = sqrt(2.*9.81.*hLoop.*(1-mu.*tan(angInitRad)));
    sLoop = R*angInitRad; 
    
    % Solve ODE
    y0 = [sLoop vLoop];
    [t,y] = ode45(@funcBlock,t,y0);
    
    % Evaluate results
    pos = y(:,1);
    vel = y(:,2);
    ang = pos/R;
    angDeg = ang*180/pi;
    fNorm = mass*(9.81*sin(ang)+vel.^2/R);

    % Check flags
    index = find((angDeg >= 270 - 1) & (angDeg <= 270 + 1));
    avg = mean(fNorm(index));
    if (avg >= 0 - .001) && (avg <= 0 + .001)
        normFlag = true;
    end

    if any(angDeg > 270)
        angFlag = true;
    end

    if angFlag && normFlag
        flag = false;
    end

    % Increment
    H = H + .00000001;
    iters = iters + 1;
end
fprintf("The minimum height for mu = %.2f is %f inches\n", mu, H * 39.37);
fprintf("It took %d iterations\n", iters);
toc
fprintf("\n");

% prepare figure 3
figure(3);
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
muPlot = plot(0, 0);

legend([posPlot,velPlot, muPlot],'Ang','Vel', 'Mu = 0.2');
hold off;

% prepare figure 4
figure(4);
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
muPlot = plot(0, 0);

legend([fNormPlot,velPlot, muPlot],'Normal Force','Vel', 'Mu = 0.2');
hold off;

%% Calculation for mu = .5
mu = .5;
H = 9.26;
iters = 0;
flag = true;
normFlag = false;
angFlag = false;

fprintf("Starting loop for mu = %.2f\n", mu)
tic
while flag
    % Setup 
    angInit = 90 - angRamp;
    angInitRad = angInit*pi/180;
    hLoop = H - R*(1-sin(angInitRad));
    vLoop = sqrt(2.*9.81.*hLoop.*(1-mu.*tan(angInitRad)));
    sLoop = R*angInitRad; 
    
    % Solve ODE
    y0 = [sLoop vLoop];
    [t,y] = ode45(@funcBlock,t,y0);
    
    % Evaluate results
    pos = y(:,1);
    vel = y(:,2);
    ang = pos/R;
    angDeg = ang*180/pi;
    fNorm = mass*(9.81*sin(ang)+vel.^2/R);

    % Check flags
    index = find((angDeg >= 270 - .1) & (angDeg <= 270 + .1));
    avg = mean(fNorm(index));
    if (avg >= 0 - .001) && (avg <= 0 + .001)
        normFlag = true;
    end

    if any(angDeg > 270)
        angFlag = true;
    end

    if angFlag && normFlag
        flag = false;
    end

    % Increment
    H = H + .00000001;
    iters = iters + 1;
end
fprintf("The minimum height for mu = %.2f is %f inches\n", mu, H * 39.37)
fprintf("It took %d iterations\n", iters);
toc
toc

% prepare figure 5
figure(5);
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
muPlot = plot(0, 0);

legend([posPlot,velPlot, muPlot],'Ang','Vel', 'Mu = 0.5');
hold off;

% prepare figure 6
figure(6);
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
muPlot = plot(0, 0);

legend([fNormPlot,velPlot, muPlot],'Normal Force','Vel', 'Mu = 0.5');
hold off;
