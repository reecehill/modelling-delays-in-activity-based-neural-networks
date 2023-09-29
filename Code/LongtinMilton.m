%The Longtin-Milton model of the pupil light reflex
%An example of a delayed differential equation model
%Choose tau (the delay) to see delay-induced
%oscillations

function LongtinMilton = LongtinMilton(tau)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f = f( x )
%The Hill function for the model
L=1;
Ld=0.1;
theta=0.5;
n=4;
f=L*theta^n./(theta^n+x.^n) +Ld;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ddeLM = ddeLM(t,y,z,P)
%the rhs of the vector field
x=y(1);
%x delayed by tau
xdelay=z(1);
ddeLM = -P.alpha*x+P.gamma*log(f(xdelay)).*(log(f(xdelay))>0);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%parameters
P.alpha=1.0;
P.gamma=3.0;      

options=[];

%initial data
x0=0.1;

tint = linspace(0,300,50000);
y = dde23(@ddeLM,tau,x0,tint,options,P);
yint = deval(y,tint);

figure(1)
plot(tint,yint(1,:));
ylabel('x')
xlabel('t')

end
