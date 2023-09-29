function nullclines(alpha,beta,a,b,c,d)
%% Code to plot the u- and v-nullclines
%% as defined by Equation 2.1 of https://doi.org/10.1098/rsta.2008.0256

%{
    -----
    %% PARAMETER DEFINITIONS START
    -----
%}

arguments
    alpha (1,1) double = 1;
    beta (1,1) double = 1;
    a (1,1) double = 10;
    b (1,1) double = -10;
    c (1,1) double = 10;
    d (1,1) double = 2;
end

% Ranges to plot for u and v.
uSpan = linspace(0,1,9999999);
vSpan = linspace(0,1,9999999);

% Ranges to loop through theta_u and theta_v.
thetaURange = [-4.8,-3.8,-2.8];
thetaVRange = [-8];
%thetaVRange = -8; %consider only theta_v=-8;

%{
    -----
    %% PARAMETER DEFINITIONS END
    -----
%}

clf;
close all;
fh = figure(1);

for thetaU=thetaURange
    for thetaV=thetaVRange
        position = find(thetaV==thetaVRange) + find(thetaU==thetaURange) -1;
        letter = char('A'+position-1)
        subplot(length(thetaVRange),length(thetaURange), position);
        % Open figures maximised by default
        % fh.WindowState = 'maximized';
        axis([-0.25 1.05 -0.05 1.05]);
        %axis padded;
        ylabel('v');
        xlabel('u');
        param_string = sprintf(' \\theta_u=%g, \\theta_v=%g',thetaU,thetaV);
        figstr = {strcat(letter,')',param_string)};
        subtitle(figstr);
        legend;
        hold on;

        %u =         -uSpan + getF(thetaU + (a .* uSpan) + (b .* vSpan));
        %v = alpha.*(-vSpan + getF(thetaV + (c .* uSpan) + (d .* vSpan)));
        [u] = (getInverseF(vSpan) - thetaV - (d .* vSpan)) ./ c;
        [v] = (getInverseF(uSpan) - thetaU - (a .* uSpan)) ./ b;
       
        plot(v, vSpan,'Color',[0 0.514 0.792],'linewidth',2,'DisplayName','v-nullcline','LineStyle',':');
        plot(uSpan, u,'Color',[0 0.514 0.792],'linewidth',2,'DisplayName','u-nullcline','LineStyle','-');
        
    end
end

    function f = getF(z)
        f = 1 ./ (1 + exp((-beta) .* z));
    end
    
    
        function inverseF = getInverseF(z)
            inverseF = log( z ./ (1-z) ) ./ beta;
        end 
    
end
