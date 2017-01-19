close all

W = 100;
L = 200;

%1000 electrons in the box
np = 500;
dt = 3;

%generating the random positions of the electrons
xp = rand(np, 1)*L;
yp = rand(np, 1)*W;

plot(xp, yp, '*r')

%make them move; generating random velocities
vth = 100;

vx = vth*(rand(100, 1) - 0.5);
vy = vth*(rand(100, 1) - 0.5);

for t = 1:100
    dx = vx*dt;
    dy = vy*dt;
    
    x = xp + dx;
    y = yp + dy;
    
    plot (x, y, '*r')
    
    xp = x;
    yp = y;
    
end 
