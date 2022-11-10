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
global R mu %% Define R and mu as global variables
R = 5 * 0.0254; %% semi-circle radius [in], converted to [m]
angRamp = 50; %% [deg], ramp angle relative to horizontal
mu = .1; %% coefficent of friction
mass = 1; %% mass of block [kg]; this is arbitrary in this problem
A = .6;
B = .61;
t = 0:0.0001:2; %% time paratmeter (start time:time step:end time) for solving
     % ODE [s]; you may want to adjust this to shorten or lengthen the simulation
H = (A + B) / 2;
iters = 0;
flag = true;

fprintf("Starting loop for mu = %.1f\n", mu)
tic
while flag
    iters = iters + 1;
    hOld = H;
    H = (A + B) / 2;
    % Setup
    angInit = 90 - angRamp;
    angInitRad = angInit*pi/180;
    hLoop = H - R*(1-sin(angInitRad));
    vHLoop = sqrt(2.*9.81.*hLoop.*(1-mu.*tan(angInitRad)));
    sLoop = R*angInitRad; %% initial loop position (defined from initial angle) [m]

    ALoop = A - R*(1-sin(angInitRad));
    vALoop = sqrt(2.*9.81.*ALoop.*(1-mu.*tan(angInitRad)));
    
    % Solve ODE
    yH0 = [sLoop vHLoop];
    yA0 = [sLoop vALoop];
    [t,yH] = ode45(@funcBlock,t,yH0);
    [~, yA] = ode45(@funcBlock, t, yA0);
    
    % Evaluate results
    posH = yH(:,1); % position
    velH = yH(:,2); % velocity
    angH = posH/R;  % angles
    angDegH = angH*180/pi; % degree versions of angles
    fNormH = mass*(9.81*sin(angH)+velH.^2/R); % normal force

    posA = yA(:,1); % position
    velA = yA(:,2); % velocity
    angA = posA/R;  % angles
    angDegA = angA*180/pi; % degree versions of angles
    fNormA = mass*(9.81*sin(angA)+velA.^2/R); % normal force
    
    % Pull out useful values
    indexH = find((angDegH >= 270 - .1) & (angDegH <= 270 + .1)); % where the angle is ~ 270 (at the top of the loop)
    avgH = mean(fNormH(indexH)); % average normal force at the top of the loop

    indexA = find((angDegA >= 270 - .1) & (angDegA <= 270 + .1)); % where the angle is ~ 270 (at the top of the loop)
    avgA = mean(fNormA(indexA)); % average normal force at the top of the loop

    ea = abs((H - hOld) / H);
    if iters > 2
        if ea < 1e-16
            break
        end
    end

    if avgA * avgH < 0
        B = H;
    else
        A = H;
    end

    if mod(iters, 10) == 0
        fprintf("Average: %.16f", avgH)
        fprintf("  Iters: %d\n", iters)
    end

end
pointOneH = H * 39.37;
fprintf("The minimum height for mu = %.1f is %.16f inches\n", mu, pointOneH);
fprintf("It took %d iterations\n", iters);
time1 = toc;
fprintf("%f seconds elapsed\n", time1)
fprintf("\n")

% prepare figure 1 (angle and velocity as a function of time)
figure(1);
hold on;
grid on; box on;
xlabel('Time (s)')

% plot position results
yyaxis left
posPlot = plot(t,angDegH,'-','LineWidth',2);
ylabel("Ang (deg)");

% plot velocity results
yyaxis right
velPlot = plot(t,velH,'-','LineWidth',2);
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
velPlot = plot(angDegH, velH,'-','LineWidth',2);
ylabel("Vel (m/s)");

% plot normal force results
yyaxis left
fNormPlot = plot(angDegH,fNormH,'-','LineWidth',2);
ylabel("Force (N)");

legend([fNormPlot,velPlot],'Normal Force','Vel');
hold off;

%% Calculation for mu = .2
mu = .2;
iters = 0;
flag = true;
A = 1;
B = 2;

fprintf('Starting loop for mu = %.1f\n', mu)
tic
while flag
    iters = iters + 1;
    hOld = H;
    H = (A + B) / 2;
    % Setup
    angInit = 90 - angRamp;
    angInitRad = angInit*pi/180;
    hLoop = H - R*(1-sin(angInitRad));
    vHLoop = sqrt(2.*9.81.*hLoop.*(1-mu.*tan(angInitRad)));
    sLoop = R*angInitRad; %% initial loop position (defined from initial angle) [m]

    ALoop = A - R*(1-sin(angInitRad));
    vALoop = sqrt(2.*9.81.*ALoop.*(1-mu.*tan(angInitRad)));
    
    % Solve ODE
    yH0 = [sLoop vHLoop];
    yA0 = [sLoop vALoop];
    [t,yH] = ode45(@funcBlock,t,yH0);
    [~, yA] = ode45(@funcBlock, t, yA0);
    
    % Evaluate results
    posH = yH(:,1); % position
    velH = yH(:,2); % velocity
    angH = posH/R;  % angles
    angDegH = angH*180/pi; % degree versions of angles
    fNormH = mass*(9.81*sin(angH)+velH.^2/R); % normal force

    posA = yA(:,1); % position
    velA = yA(:,2); % velocity
    angA = posA/R;  % angles
    angDegA = angA*180/pi; % degree versions of angles
    fNormA = mass*(9.81*sin(angA)+velA.^2/R); % normal force
    
    % Pull out useful values
    indexH = find((angDegH >= 270 - .1) & (angDegH <= 270 + .1)); % where the angle is ~ 270 (at the top of the loop)
    avgH = mean(fNormH(indexH)); % average normal force at the top of the loop

    indexA = find((angDegA >= 270 - .1) & (angDegA <= 270 + .1)); % where the angle is ~ 270 (at the top of the loop)
    avgA = mean(fNormA(indexA)); % average normal force at the top of the loop

    ea = abs((H - hOld) / H);
    if iters > 2
        if ea < 1e-16
            break
        end
    end

    if avgA * avgH < 0
        B = H;
    else
        A = H;
    end

    if mod(iters, 10) == 0
        fprintf("Average: %.16f", avgH)
        fprintf("  Iters: %d\n", iters)
    end

end
pointTwoH = H * 39.37;
fprintf("The minimum height for mu = %.1f is %.16f inches\n", mu, pointTwoH);
fprintf("It took %d iterations\n", iters);
time2 = toc;
fprintf("%f minutes elapsed\n", (time2 / 60))
fprintf("\n");

% prepare figure 3
figure(3);
hold on;
grid on; box on;
xlabel('Time (s)')

% plot position results
yyaxis left
posPlot = plot(t,angDegH,'-','LineWidth',2);
ylabel("Ang (deg)");

% plot velocity results
yyaxis right
velPlot = plot(t,velH,'-','LineWidth',2);
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
velPlot = plot(angDegH, velH,'-','LineWidth',2);
ylabel("Vel (m/s)");

% plot normal force results
yyaxis left
fNormPlot = plot(angDegH,fNormH,'-','LineWidth',2);
ylabel("Force (N)");
muPlot = plot(0, 0);

legend([fNormPlot,velPlot, muPlot],'Normal Force','Vel', 'Mu = 0.2');
hold off;

%% Calculation for mu = .5
mu = .5;
A = 9;
B = 9.5;
iters = 0;
flag = true;

fprintf("Starting loop for mu = %.1f\n", mu)
tic
while flag
    iters = iters + 1;
    hOld = H;
    H = (A + B) / 2;
    % Setup
    angInit = 90 - angRamp;
    angInitRad = angInit*pi/180;
    hLoop = H - R*(1-sin(angInitRad));
    vHLoop = sqrt(2.*9.81.*hLoop.*(1-mu.*tan(angInitRad)));
    sLoop = R*angInitRad; %% initial loop position (defined from initial angle) [m]

    ALoop = A - R*(1-sin(angInitRad));
    vALoop = sqrt(2.*9.81.*ALoop.*(1-mu.*tan(angInitRad)));
    
    % Solve ODE
    yH0 = [sLoop vHLoop];
    yA0 = [sLoop vALoop];
    [t,yH] = ode45(@funcBlock,t,yH0);
    [~, yA] = ode45(@funcBlock, t, yA0);
    
    % Evaluate results
    posH = yH(:,1); % position
    velH = yH(:,2); % velocity
    angH = posH/R;  % angles
    angDegH = angH*180/pi; % degree versions of angles
    fNormH = mass*(9.81*sin(angH)+velH.^2/R); % normal force

    posA = yA(:,1); % position
    velA = yA(:,2); % velocity
    angA = posA/R;  % angles
    angDegA = angA*180/pi; % degree versions of angles
    fNormA = mass*(9.81*sin(angA)+velA.^2/R); % normal force
    
    % Pull out useful values
    indexH = find((angDegH >= 270 - .1) & (angDegH <= 270 + .1)); % where the angle is ~ 270 (at the top of the loop)
    avgH = mean(fNormH(indexH)); % average normal force at the top of the loop

    indexA = find((angDegA >= 270 - .1) & (angDegA <= 270 + .1)); % where the angle is ~ 270 (at the top of the loop)
    avgA = mean(fNormA(indexA)); % average normal force at the top of the loop

    ea = abs((H - hOld) / H);
    if iters > 2
        if ea < 1e-16
            break
        end
    end

    if avgA * avgH < 0
        B = H;
    else
        A = H;
    end

    if mod(iters, 10) == 0
        fprintf("Average: %.16f", avgH)
        fprintf("  Iters: %d\n", iters)
    end

end
pointFiveH = H * 39.37;
fprintf("The minimum height for mu = %.1f is %.16f inches\n", mu, pointFiveH)
fprintf("It took %d iterations\n", iters);
time3 = toc;
fprintf("%f minutes elapsed\n", (time3 / 60))
% prepare figure 5
figure(5);
hold on;
grid on; box on;
xlabel('Time (s)')

% plot position results
yyaxis left
posPlot = plot(t,angDegH,'-','LineWidth',2);
ylabel("Ang (deg)");

% plot velocity results
yyaxis right
velPlot = plot(t,velH,'-','LineWidth',2);
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
velPlot = plot(angDegH, velH,'-','LineWidth',2);
ylabel("Vel (m/s)");

% plot normal force results
yyaxis left
fNormPlot = plot(angDegH,fNormH,'-','LineWidth',2);
ylabel("Force (N)");
muPlot = plot(0, 0);

legend([fNormPlot,velPlot, muPlot],'Normal Force','Vel', 'Mu = 0.5');
hold off;

fprintf("\nmu = .1  H = %.10f inches\n", pointOneH)
fprintf("mu = .2  H = %.10f inches\n", pointTwoH)
fprintf("mu = .5  H = %.10f inches\n", pointFiveH)
fprintf("Total time elapsed: %f minutes\n", (time1 + time2 + time3) / 60)
