clc;
clear;
close all;
% Defining initial conditions
x1 = -1;
y1 = 0;
vx1 = 0.3471128135672417;
vy1 = 0.532726851767674;
x2 = 1;
y2 = 0;
vx2 = vx1;
vy2 = vy1;
x3 = 0;
y3 = 0;
vx3 =-2*vx1;
vy3 =-2*vy1; 
% plotting the initial positions of three bodies
plot(x1,y1,'or','MarkerSize',5,'MarkerFaceColor','red');
hold on;
plot(x2,y2,'og','MarkerSize',5,'MarkerFaceColor','green');
hold on;
plot(x3,y3,'ob','MarkerSize',5,'MarkerFaceColor','blue');
hold on;
% Defining time steps and duration
t = 0:0.001:5;
h = 0.001; 
n=length(t);
% Initializing arrays to store position and velocity
pos = zeros(n, 6);
pos(1,:) = [x1, y1, x2, y2, x3, y3];
vel = zeros(n, 6);
vel(1,:) = [vx1, vy1, vx2, vy2, vx3, vy3];

% Defining constants
G = 1; 
m=1;
f=zeros(n,6);
a=zeros(n,6);
% RK-4 numerical method
for i = 2:n
    % Calculating the distance between the masses
    r12 = norm(pos(i-1,1:2) - pos(i-1,3:4));
    r13 = norm(pos(i-1,1:2) - pos(i-1,5:6));
    r23 = norm(pos(i-1,3:4) - pos(i-1,5:6));
    
    % Calculating the gravitational force on each mass
   f(i-1,1) = G*m^2*(pos(i-1,1)-pos(i-1,3))/r12^3;
   f(i-1,2) = G*m^2*(pos(i-1,2)-pos(i-1,4))/r12^3;

   f(i-1,3) = G*m^2*(pos(i-1,1)-pos(i-1,5))/r13^3;
   f(i-1,4) = G*m^2*(pos(i-1,2)-pos(i-1,6))/r13^3;

   f(i-1,5) = G*m^2*(pos(i-1,3)-pos(i-1,5))/r23^3;
   f(i-1,6) = G*m^2*(pos(i-1,4)-pos(i-1,6))/r23^3;

    % Calculating the acceleration on each mass
   
    a(i-1,1)=(-f(i-1,1)-f(i-1,3))/m;
    a(i-1,2)=(-f(i-1,2)-f(i-1,4))/m;

    a(i-1,3)=(f(i-1,1)-f(i-1,5))/m;
    a(i-1,4)=(f(i-1,2)-f(i-1,6))/m;

    a(i-1,5)=(f(i-1,3)+f(i-1,5))/m;
    a(i-1,6)=(f(i-1,4)+f(i-1,6))/m;

% Updating velocities and positions using RK-4 method
k1v = a(i-1,:);
k1r = vel(i-1,:);
k2v = a(i-1,:);
k2r = vel(i-1,:) + k1v*h/2;
k3v = a(i-1,:);
k3r = vel(i-1,:) + k2v*h/2;
k4v = a(i-1,:);
k4r = vel(i-1,:) + k3v*h;
vel(i,:) = vel(i-1,:) + (k1v + 2*k2v + 2*k3v + k4v)*h/6;
pos(i,:) = pos(i-1,:) + (k1r + 2*k2r + 2*k3r + k4r)*h/6; 
end
for i=1:n
plot(pos(i,1), pos(i,2), 'or', 'MarkerSize', 2, 'MarkerFaceColor', 'red');
hold on;
plot(pos(i,3), pos(i,4), 'og', 'MarkerSize', 2, 'MarkerFaceColor', 'green');
hold on;
plot(pos(i,5), pos(i,6), 'ob', 'MarkerSize', 2, 'MarkerFaceColor', 'blue');
drawnow;
xlabel('x');
ylabel('y');
axis([-1.5 1.5 -0.4 0.4]);
axis equal;
end
