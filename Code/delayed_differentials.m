%% Code to simulate the effect of delays on neural network activity
% Using parameters from: https://doi.org/10.1098/rsta.2008.0256

function delayed_differentials( ...
    thetaU, ...
    thetaV, ...
    a, ...
    b, ...
    c, ...
    d, ...
    tau1Range, ...
    tau2Range, ...
    %{
    minU, ...
    maxU, ...
    uStep, ...
    minV, ...
    maxV, ...
    vStep, ...
    %}
    alpha, ...
    beta, ...
    time)

arguments
    %% BACKGROUND DRIVE/BIAS
    thetaU (1,1) {mustBeNumeric} = -2;
    thetaV (1,1) {mustBeNumeric} = -4;

    %% WEIGHTS
    a (1,1) {mustBeNumeric} = 10; % u self-feedback
    b (1,1) {mustBeNumeric} = -10; % v on u
    c (1,1) {mustBeNumeric} = 10; % u on v
    d (1,1) {mustBeNumeric} = 2; % v self-feedback

    %% TIME DELAY
    tau1Range (1,:) {mustBeVector} = [0,1]; % delay for self-feedback
    tau2Range (1,:) {mustBeVector} = [0,1]; % delay between populations

    %% INITIAL VALUES (NOT USED)
    %{
    minU (1,1) {mustBeNumeric} = 0;
    maxU (1,1) {mustBeGreaterThanOrEqual(maxU,minU)} = 1;
    uStep (1,1) {mustBePositive} = 0.1;
    
    minV (1,1) {mustBeNumeric} = 0;
    maxV (1,1) {mustBeGreaterThanOrEqual(maxV,minV)} = 1;
    vStep (1,1) {mustBePositive} = 0.1;

    %}

    %% PARAMETERS
    alpha (1,1) {mustBeNumeric} = 1;
    beta (1,1) {mustBeNumeric} = 1;

    %% SIMULATION OPTIONS
    time (2,1) {mustBeNumeric} = [0,600]; % time to simulate for
end
uSpan = linspace(0,1,10);
vSpan = linspace(0,1,10);
close all;
clf;
colours = hsv(length(tau1Range));

for tau2 = tau2Range
    for tau1 = tau1Range
        if(tau1 ~= tau2)
            continue;
        end
        figure(1);
        hold on;
        axis([0 1 0 1]);
        position = find(tau1==tau1Range);
        subplot(1,length(tau1Range),position);
        initialise(thetaU,thetaV,a,b,c,d,alpha,beta,tau1,tau2);

        if(tau1 ~= 0 && tau2 ~=0)
            solutions = dde23(@(t,y,Z) getDifferential(t,y,Z,thetaU,thetaV,a,b,c,d,alpha,beta),[tau1, tau2],@ddex1hist,time);
            plot(solutions.x,solutions.y(1,:),'DisplayName','u population');
            plot(solutions.x,solutions.y(2,:),'DisplayName','v population');

        else
            allDUs = zeros(length(vSpan),length(uSpan));
            allDVs = zeros(length(vSpan),length(uSpan));
            vIndex = 0;

            for initialV = vSpan
                vIndex = vIndex + 1;
                uIndex = 0;
                for initialU = uSpan
                    uIndex = uIndex + 1;
                    [t,y] = ode45(@(t,y) getDifferentialWithoutLags(t,y,thetaU,thetaV,a,b,c,d,alpha,beta),time,[initialU,initialV]);
                    u = y(:,1);
                    v = y(:,2);
                    du = gradient(u);
                    dv = gradient(v);
                    allDUs(uIndex,vIndex) = du(1);
                    allDVs(uIndex,vIndex) = dv(1);
                end
            end

            hold on;

            [t,y] = ode45(@(t,y) getDifferentialWithoutLags(t,y,thetaU,thetaV,a,b,c,d,alpha,beta),time,[0,0]);
            plot(t,y(:,1));
            plot(t,y(:,2));
        end

    end
    for tau1 = tau1Range
        if(tau1 ~= tau2)
            continue;
        end
        position = find(tau1==tau1Range);

        figure(5);

        leg = legend('show');
        title(leg,'Simulations')
        hold on;
        title('(U,V) Phase Plane');
        xlabel('v population');
        ylabel('u population');
        axis([0 0.6 0 0.9]);
        if(tau1 ~= 0 && tau2 ~=0)
            solutions = dde23(@(t,y,Z) getDifferential(t,y,Z,thetaU,thetaV,a,b,c,d,alpha,beta),[tau1, tau2],@ddex1hist,time);
            u = solutions.y(1,:);
            v = solutions.y(2,:);
            label = sprintf('\\tau_1=%g,  \\tau_2=%g',tau1,tau2);
        else
            [t,y] = ode45(@(t,y) getDifferentialWithoutLags(t,y,thetaU,thetaV,a,b,c,d,alpha,beta),time,[0,0]);
            u = y(:,1);
            v = y(:,2);
            label = sprintf('\\tau_1=%g,  \\tau_2=%g',tau1,tau2);
        end
        du = gradient(u);
        dv = gradient(v);
        q = quiver(u,v,du,dv,3,'LineWidth',2);
        
        q.HandleVisibility = 'off';
        plot(u,v,'DisplayName','Simulation','color',q.Color,'DisplayName',label,'LineWidth',2);
    end

end
end
function initialise(thetaU,thetaV,a,b,c,d,alpha,beta,tau1,tau2)
%{
    Clear previous simulation from memory.
%}
%close all;
%clf;

%{
    Prepare figure.
%}
figure(1);
hold on;
legend;
param_string = sprintf('\\tau_1=%g,  \\tau_2=%g',tau1,tau2);
subtitle('Wilson & Cowan Model (with delays)');
subtitle({param_string});
xlabel('Time, t');
ylabel('Rate of change of neuronal population activity ');
end

function dudt = getDifferential(t,y,Z,thetaU,thetaV,a,b,c,d,alpha,beta)
% Differential equations function for DDEX1.
u = y(1);
v = y(2);
uLag = Z(:,1);
vLag = Z(:,2);
uLagTau1 = uLag(1);
uLagTau2 = uLag(2);
vLagTau1 = vLag(1);
vLagTau2 = vLag(2);

dudt = [
    % duDt;
    (-u) + f( (thetaU) + a*(uLagTau1) + b*(vLagTau2), beta);

    % dvDt;
    alpha * ( (-v) + f( ((thetaV) + c*(vLagTau1) + d*(uLagTau2)), beta));
    ];

end

function dydt = getDifferentialWithoutLags(t,y,thetaU,thetaV,a,b,c,d,alpha,beta)
u = y(1);
v = y(2);
uLagTau1 = 1;
vLagTau1 = 1;
vLagTau2 = 1;
uLagTau2 = 1;

dydt = [
    % duDt;
    (-u) + f( (thetaU) + a*(uLagTau1) + b*(vLagTau2), beta);

    % dvDt;
    alpha * ( (-v) + f( ((thetaV) + c*(vLagTau1) + d*(uLagTau2)), beta));
    ];
end

function s = ddex1hist(t)
% Constant history function for DDEX1.

%s = ones(2,1);
s = zeros(2,1);
end

function output = f(input,beta)
output = 1 ./ ( 1 + exp((-beta) * input) );
end
