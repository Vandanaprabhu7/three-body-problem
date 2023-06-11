close all;
clear;
clc;

% defining intial conditions
x=-1;
y=0;
vx=0;
vy=1.20000;
R0=[x;y;vx;vy];

% time duration and time step
h=0.1;
t=0:0.1:500;
l=length(t);
% defining a vector
R=zeros(4,l);
% storing initial conditions in first row
R(:,1)=R0;

% rk-4 algorithm
for i=1:l-1
    k1=h*eqn(t(i),R(:,i));
    k2=h*eqn(t(i)+h/2,R(:,i)+k1/2);
    k3=h*eqn(t(i)+h/2,R(:,i)+k2/2);
    k4=h*eqn(t(i)+h,R(:,i)+k3);
    R(:,i+1)=R(:,i)+(k1+2*k2+2*k3+k4)/6;
end

% plotting the sun and orbit of earth
plot(0,0,'or','MarkerSize',10,'MarkerFaceColor','#f80');
hold on;
axis equal;
xlabel('x');
ylabel('y');
plot(R(1,:),R(2,:),'or','MarkerSize',1,'MarkerFaceColor','red');

% the function containing the odes
function Rr=eqn(~,R)
    x=R(1);
    y=R(2);
    dxdt=R(3);
    dydt=R(4);
    d2xdt2=-x/((x^2 + y^2)^(3/2));
    d2ydt2=-y/((y^2 + x^2)^(3/2));
    Rr=[dxdt; dydt; d2xdt2; d2ydt2];
end
