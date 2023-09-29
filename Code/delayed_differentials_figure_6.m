%% Code to simulate the effect of time delays on neural network activity
%{
Using parameters from: https://doi.org/10.1098/rsta.2008.0256
This script is for exploring many (tau1,tau2) values and initial conditions
at once. It will produce nullclines, time courses, and figures that 
contain subplots for each tau pair, where relevant.

To begin, edit parameters within the arguments block.
%} 

%% Optional parameters
% Optional parameters, such as how colours are generated, are found within
% the code.

%% Code by Reece Hill.
% Code has *not* been optimised. 

function code( ...
    thetaU, ...
    thetaV, ...
    a, ...
    b, ...
    c, ...
    d, ...
    tau1Range, ...
    tau2Range, ...
    alpha, ...
    beta, ...
    time, ...
    initialUSpan, ...
    initialVSpan, ...
    runSimulation)

arguments
    %% BACKGROUND DRIVE/BIAS
    thetaU (1,1) {mustBeNumeric} = 0.2;
    thetaV (1,1) {mustBeNumeric} = 0.2;

    %% WEIGHTS
    a (1,1) {mustBeNumeric} = -6; % u self-feedback
    b (1,1) {mustBeNumeric} = 2; % v on u
    c (1,1) {mustBeNumeric} = 2; % u on v
    d (1,1) {mustBeNumeric} = -6; % v self-feedback

    %% TIME DELAY
    % Unlike previous scripts, here Tau1 and Tau2 are taken as pairs.
    % Therefore lengths must be the same.
    tau1Range (1,:) {mustBeVector} = [0.1]; % delay for self-feedback
    tau2Range (1,:) {mustBeVector} = [0.1]; % delay between populations

    %% PARAMETERS
    alpha (1,1) {mustBeNumeric} = 1;
    beta (1,1) {mustBeNumeric} = 60;

    %% SIMULATION OPTIONS
    % For machines with less RAM, consider lowering the number of
    % timepoints (Z).
    time (1,:) {mustBeNumeric} = linspace(0,25,1*10^9); % simulate from X, to Y, and generate Z timepoints between

    % Unlike previous scripts, here initialUSpan and initialVSpan are taken as pairs.
    % Therefore lengths must be the same.
    initialUSpan (1,:) {mustBeNumeric} = [0.065,0.065];
    initialVSpan (1,:) {mustBeNumeric} = [0.09,0.0901];
    runSimulation (1,1) {mustBeNumeric} = 1; % 1 = run; otherwise = don't run (for speed)
end

%% Close anything open, and reset figures.
clf;
close all;

% Time courses
figure(1);
hold on;

% Phase planes with quivers
figure(2);
hold on;

% Nullclines
figure(3);
hold on;

% Nullclines on phase planes
figure(4);
hold on;

%% Additional figures, used for fine-tuning appearance without affecting base figure.
figure(10);
hold on;
figure(11);
hold on;



% Pre-determine the number of iterations about to be performed
% (case-by-case basis, uncomment as needed)
%N = length(tau1Range) * length(tau2Range);
N =  length(initialUSpan);

% Pre-determine the colours used so that theyre consistent.
%(case-by-case basis, uncomment as needed)
%colours = rand(N,3);
colours = hsv(N);
%colours = hsv(2);

% Set initial position of loops to zero.
position = 0;

%% Begin looping through Tau1 and Tau2 pairs.

tauIndex = 0;
for tau2 = tau2Range
    tauIndex = tauIndex + 1;
    position = position + 1; %TODO (performance): is position==tauIndex??
    tau1 = tau1Range(tauIndex);

    %% Set as needed
    %rowsOfSubPlots = length(tau2Range);
    rowsOfSubPlots = 1;
    columnsOfSubPlots = length(tau1Range);

    % TODO: Consider pre-allocating AllOfU for speed.
    allOfU = 0;

    if(runSimulation == 1)
        index = 0; %TODO (performance): is index==initialIndex??
        initialIndex = 0; % index for current initialCondition pair.

        for initialU = initialUSpan
            initialIndex = initialIndex + 1;
            index = index + 1;
            initialV = initialVSpan(initialIndex);

            % Simulate delayed differential eq.
            solutions = dde23(@(t,y,Z) getDifferential(t,y,Z,thetaU,thetaV,a,b,c,d,alpha,beta),[tau1, tau2],@(t) ddex1hist(t,initialU,initialV),time);
            u = solutions.y(1,:);
            v = solutions.y(2,:);

            %% Time courses
            figure(1);
            %sgtitle('Simulations: Neuronal Activity for Populations U and V');
            subplot(rowsOfSubPlots,columnsOfSubPlots,position);
            hold on;
            axis([0 inf 0 1]);
            ax = gca;
            ax.FontSize = 16;
            leg = legend('show');
            title(leg,sprintf('population_{u(t_0),v(t_0)}'));
            xlabel('Time, t');
            ylabel('Rate of change of neuronal population activity ');
            param_string = sprintf('\\tau_1=%g,  \\tau_2=%g',tau1,tau2);
            subtitle({param_string});
            uLegendLabel = sprintf('u_{(%g,%g)}',initialU,initialV);
            vLegendLabel = sprintf('v_{(%g,%g)}',initialU,initialV);
            colour = colours(index,:);
            plot(solutions.x,u,'linewidth',2,'DisplayName',uLegendLabel,'Color',[colour 0.5]);
            plot(solutions.x,v,'linewidth',2,'DisplayName',vLegendLabel,'LineStyle',':','Color',[colour 0.5]);

            %% Additional graph, used for fine-tuning appearance for specific figures.
            figure(10);
            hold on;
            legend;
            xlabel('u population');
            ylabel('v population');
            zlabel('time');
            axis([0.1 0.2 0.1 0.2 0 inf]);
            plot(u,v,'LineWidth',2,'DisplayName',uLegendLabel);

            % Add u values to bank of preceding us for use outside of loop.
            if(allOfU == 0)
                allOfU = zeros(1,length(u));
            end
            if(length(allOfU) > length(u))
                allOfU = allOfU(:,1:length(u));
            end
            allOfU = [allOfU(:,:); u(1:length(allOfU))];

            %% Phase planes
            figure(2);
            sgtitle('Phase Planes');
            subplot(rowsOfSubPlots,columnsOfSubPlots,position);
            axis([-0.05 1.05 -0.05 1.05]);
            leg2 = legend('show');
            title(leg2,sprintf('Simulation_{u(t_0),v(t_0)}'));
            hold on;
            param_string = sprintf('\\tau_1=%g,  \\tau_2=%g',tau1,tau2);
            %subtitle({param_string});
            xlabel('u population');
            ylabel('v population');
            legendLabel = sprintf('Simulation_{(%g,%g)}',initialU,initialV);

            plot(u,v,'linewidth',1,'DisplayName',legendLabel,'Color',[colour 1]);

            du = gradient(u);
            dv = gradient(v);
            q = quiver(u,v,du,dv,2,'Color',colour);
            q.HandleVisibility = 'off';

        end
    end

    %% Prepare for Figure 11, proving SIDC.
    currentU = allOfU(2,:);
    nextU = allOfU(3,:);
    y = currentU - nextU;
    x = solutions.x;

    figure(11);

    % U1 activity over time
    subplot(2,1,1);
    axis([0 inf 0.025 0.125]);
    ax= gca;
    ax.FontSize=14;
    ylabel('u_1(t)');
    hold on;
    plot(x(1:length(currentU)),currentU);

    % Difference in u1 and u2 activity over time
    subplot(2,1,2);
    axis([0 inf -0.06 0.05]);
    ax= gca;
    ax.FontSize=14;
    hold on;
    plot(x(1:length(y)),y);
    ylabel('u_1(t) - u_2(t)');
    xlabel('time (t)');

    %% Nullclines
    plotNullClines(thetaU, thetaV, a, b, c, d, tau1,tau2, alpha, beta, rowsOfSubPlots,columnsOfSubPlots,position);

    %% Superimpose nullclines onto phase planes
    nullclinesOnSimulation(rowsOfSubPlots,columnsOfSubPlots,position);

    %% Explore quiver maps for delays and no-delays (REDUNDANT)
    %plotQuiverMap(rowsOfSubPlots,columnsOfSubPlots,position);

end

    function plotNullClines(thetaU, thetaV, a, b, c, d, tau1,tau2, alpha, beta,rowsOfSubPlots,columnsOfSubPlots,position)
        figure(3);
        sgtitle('Nullclines');

        subplot(rowsOfSubPlots,columnsOfSubPlots,position);
        param_string = sprintf('\\tau_1=%g,  \\tau_2=%g',tau1,tau2);
        subtitle({param_string});
        legend;
        ylabel('v population');
        xlabel('u population');
        hold on;
        % Ranges to plot for u and v.
        uSpan = [0:0.0000001:1];
        vSpan = uSpan;

        u = (getInverseF(vSpan,beta) - thetaV - (d .* vSpan)) ./ c;
        v = (getInverseF(uSpan,beta) - thetaU - (a .* uSpan)) ./ b;


        plot(vSpan, v,'o-','Color',[0 0 0 0.4],'LineWidth',2,'MarkerSize',4,'MarkerFaceColor',[0 0 0],'MarkerIndices',[1:800000:length(vSpan)],'DisplayName','v-nullcline');
        plot(u, uSpan,'s-','Color',[0 0 0 0.4],'linewidth',2,'MarkerSize',4,'MarkerFaceColor',[0 0 0],'MarkerIndices',[1:1200000:length(u)],'DisplayName','u-nullcline');
        axis([0 1 0 1]);
        hold on;

        function inverseF = getInverseF(z,beta)
            inverseF = log( z ./ (1-z) ) ./ beta;
        end
    end

    function nullclinesOnSimulation(rowsOfSubPlots,columnsOfSubPlots,position)
        figure(2);
        f2subplot = subplot(rowsOfSubPlots,columnsOfSubPlots,position);
        figure(3);
        f3subplot = subplot(rowsOfSubPlots,columnsOfSubPlots,position);

        figure(4);
        hold on;
        subplot(rowsOfSubPlots,columnsOfSubPlots,position);
        hold on;
        xlabel('u population');
        ylabel('v population');
        param_string = sprintf('\\tau_1=%g,  \\tau_2=%g',tau1,tau2);
        subtitle({param_string});
        sgtitle('U-V Phase Planes')
        hold on;
        axis('padded');
        axis([0 1 0 1]);
        ax = gca;
        ax.FontSize = 16

        copyobj(findobj(f2subplot, 'Type', 'line'), gca(figure(4)));
        %copyobj(findobj(f2subplot, 'Type', 'quiver'), gca(figure(4)));
        copyobj(findobj(f3subplot, 'Type', 'line'), gca(figure(4)));
        leg3 = legend('show');
        title(leg3,sprintf('Simulation_{u(t_0),v(t_0)}'));
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

    function s = ddex1hist(t,initialU,initialV)
        % Constant history function for DDEX1.

        %s = ones(2,1);
        s = zeros(2,1) + [initialU;initialV];
    end

    function output = f(input,beta)
        output = 1 ./ ( 1 + exp((-beta) * input) );
    end

    function plotQuiverMap(rowsOfSubPlots,columnsOfSubPlots,position)
        index2 = 0;
        initialQUSpan = [0,0.5,1];
        initialQVSpan = [0,0.5,1];
        du1 = zeros(length(initialUSpan),length(initialVSpan));
        dv1 = zeros(length(initialUSpan),length(initialVSpan));
        timeQ = [0,20];
        for initialQU = initialQUSpan
            index2 = index2 + 1;
            index1 = 0;
            for initialQV = initialQVSpan
                index1 = index1 + 1;

                solutionsQ = dde23(@(t,y,Z) getDifferential(t,y,Z,thetaU,thetaV,a,b,c,d,alpha,beta),[tau1, tau2],@(t) ddex1hist(t,initialQU,initialQV),timeQ);
                uQ = solutionsQ.y(1,:);
                vQ = solutionsQ.y(2,:);
                z = solutionsQ.x(:);

                firstQuarterOfuQ = uQ(1:end-(ceil(length(uQ)*0.05)));
                first5OfuQ = uQ(5:10);
                first5OfvQ = vQ(5:10);
                firstQuarterOfvQ = vQ(1:end-(ceil(length(vQ)*0.05)));
                lastQuarterOfuQ = uQ(end-(ceil(length(uQ)*0.01)):end);
                lastQuarterOfvQ = vQ(end-(ceil(length(vQ)*0.01)):end);
                middleHalfOfuQ = uQ(ceil(length(uQ)*0.2):end-ceil(length(uQ)*0.70));
                middleHalfOfvQ = vQ(ceil(length(vQ)*0.2):end-ceil(length(vQ)*0.7));
                averageFirstuQ = mean(firstQuarterOfuQ);
                averageFirstvQ = mean(firstQuarterOfvQ);
                averageLastuQ = mean(lastQuarterOfuQ);
                averageLastvQ = mean(lastQuarterOfvQ);

                %duQTemp= gradient(firstQuarterOfuQ);
                %dvQTemp= gradient(firstQuarterOfvQ);

                duQTemp= averageLastuQ - u(1);
                dvQTemp= averageLastvQ - v(1);

                du1(index2,index1) = duQTemp(1);
                dv1(index2,index1) = dvQTemp(1);

                figure(5);
                subplot(rowsOfSubPlots,columnsOfSubPlots,position);
                hold on;
                plot(firstQuarterOfuQ,firstQuarterOfvQ,'Color','r');

                %                 figure(7);
                %                 subplot(rowsOfSubPlots,columnsOfSubPlots,position);
                %                 hold on;
                %
                %                 plot3(fliplr(u),fliplr(v),fliplr(z));
                %                 xlabel('u population');
                %                 ylabel('v population');
                %                 zlabel('time');
                %                 set(gca, 'ZDir', 'reverse');
            end
        end



        figure(5);

        f5subplot = subplot(rowsOfSubPlots,columnsOfSubPlots,position);
        hold on;
        param_string = sprintf('\\tau_1=%g,  \\tau_2=%g',tau1,tau2);
        subtitle({param_string});
        ylabel('v population');
        xlabel('u population');
        axis([0 1 0 1]);

        quiver(initialQUSpan,initialQVSpan,du1,dv1,1,'Color','k');

        figure(4);
        hold on;
        subplot(rowsOfSubPlots,columnsOfSubPlots,position);
        axis([0.06 0.12 0.055 0.105]);
        %copyobj(findobj(f5subplot, 'Type', 'quiver'), gca(figure(4)));
        legend('show');
    end
end