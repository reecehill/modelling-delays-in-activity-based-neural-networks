%% Code to simulate the effect of delays on neural network activity
% Using parameters from: https://doi.org/10.1098/rsta.2008.0256

function delayed_differentials_heaviside_f( ...
    thetaU, ...
    thetaV, ...
    a, ...
    b, ...
    c, ...
    d, ...
    tau1, ...
    tau2, ...
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
    thetaU (1,1) {mustBeNumeric} = 0.7;
    thetaV (1,1) {mustBeNumeric} = 0.7;

    %% WEIGHTS
    a (1,1) {mustBeNumeric} = -1; % u self-feedback
    b (1,1) {mustBeNumeric} = -0.4; % v on u
    c (1,1) {mustBeNumeric} = -0.4; % u on v
    d (1,1) {mustBeNumeric} = -1; % v self-feedback

    %% TIME DELAY
    tau1 (1,1) {mustBePositive} = [1]; % delay for self-feedback
    tau2 (1,1) {mustBePositive} = 1.4; % delay between populations

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
    time (2,1) {mustBeNumeric} = [0,100]; % time to simulate for
end

close all;
clf;

for tau1 = tau1
initialise(thetaU,thetaV,a,b,c,d,alpha,beta,tau1,tau2);
solutions = dde23(@(t,y,Z) getDifferential(t,y,Z,thetaU,thetaV,a,b,c,d,alpha,beta),[tau1, tau2],@ddex1hist,time);
plot(solutions.x,solutions.y(1,:),'DisplayName','u population');
plot(solutions.x,solutions.y(2,:),'DisplayName','v population');
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
figure();
hold on;
legend;
param_string = sprintf('\\theta_u=%g,  \\theta_v=%g \n a=%g,  b=%g,  c=%g,  d=%g \n \\tau_1=%g,  \\tau_2=%g \n \\alpha=%g,  \\beta=%g',thetaU,thetaV,a,b,c,d,tau1,tau2,alpha,beta);
title('Wilson & Cowan Model (with delays)');
subtitle({param_string});
xlabel('Time, t');
ylabel('Rate of change of neuronal population activity ');
end

function dydt = getDifferential(t,y,Z,thetaU,thetaV,a,b,c,d,alpha,beta)
% Differential equations function for DDEX1.
u = y(1);
v = y(2);
uLag = Z(:,1);
vLag = Z(:,2);
uLagTau1 = uLag(1);
uLagTau2 = uLag(2);
vLagTau1 = vLag(1);
vLagTau2 = vLag(2);

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
%output = 1 ./ ( 1 + exp((-beta) * input) );
output = heaviside(input);

end
