function equation2point6(alpha,beta,a,b,c,d)
arguments
    alpha (1,1) double = 1;
    beta (1,1) double = 1;
    a (1,1) double = 10;
    b (1,1) double = -10;
    c (1,1) double = 10;
    d (1,1) double = 2;
end

% Ranges to plot for u and v.
u = linspace(0.1,1,5);
vSpan = linspace(0,1,5);


[sqroot] = sqrt(1 - ((4.*((1+alpha) ./ ( beta-a.*u.*(1-u) ))) ./ (alpha*d)));

[v] = [0.5 + (0.5.*sqroot), 0.5 - (0.5.*sqroot)];
figure();
plot( v, [u, -1.*u]);

end