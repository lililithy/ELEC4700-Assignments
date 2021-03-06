close all
clear all
clc

%dimensions of the box (nm)
W = 100e-9;
L = 200e-9;

%1000 electrons in the box (sometimes set to much lower value for testing)
np = 10;

tau = 0.2e-12;
K = 1.38e-23;
m = 0.26*9.11e-31;
T = 300;

%generating the random positions of the electrons
xp = rand(np, 1)*L; 
yp = rand(np, 1)*W;

plot(xp, yp, '.r')
hold on

line([0.7999e-7 0.80001e-7],[0 0.4e-7])
line([1.1999e-7 1.2001e-7],[0 0.4e-7])
line([0.8e-7 1.2e-7], [0.3999e-7 0.4001e-7])

line([0.7999e-7 0.80001e-7],[0.6e-7 1e-7])
line([1.1999e-7 1.2001e-7],[0.6e-7 1e-7])
line([0.8e-7 1.2e-7], [0.5999e-7 0.60001e-7])

%make them move; generating random velocities
vth = sqrt(2*K*T/m);

vx = randn(np,1)*vth/sqrt(2); 
vy = randn(np,1)*vth/sqrt(2);

dt = L/100/vth;

for t = 1:100
    
    %change in distance, ie how much an electron moves in x & y direction
    dx = dt*vx;
    dy = dt*vy;
    
    x = xp + dx; 
    y = yp + dy; 
    
    % if a particle goes too far right
    gx = x > L;
    x(gx) = x(gx) - L;
    % if partcle goes too far left
    lx = x < 0;
    x(lx) = x(lx) + L;
    
    %if particle touches the top or bottom, it must reflect (change the
    %velocity sign!)
    gy = y > W; 
    vy(gy) = -vy(gy);
        
    ly = y < 0;
    vy(ly) = -vy(ly);    
    
    % boxes
    % bottom box, left & right sides
    bbls = (x > 0.8e-7 & x < 1.2e-7 & y < 0.4e-7 & xp < x & xp < 0.8e-7);
    vx(bbls) = - vx(bbls);
    x(bbls) = 0.8e-7;
    
    bbrs = (x > 0.8e-7 & x < 1.2e-7 & y < 0.4e-7 & xp > x & xp > 1.2e-7);
    vx(bbrs) = - vx(bbrs);
    x(bbrs) = 1.2e-7;
    
    %top box, left & right sides
    tbls = (x > 0.8e-7 & x < 1.2e-7 & y > 0.6e-7 & xp < x & xp < 0.8e-7);
    vx(tbls) = - vx(tbls);
    x(tbls) = 0.8e-7;     
    
    tbrs = (x > 0.8e-7 & x < 1.2e-7 & y > 0.6e-7 & xp > x & xp > 1.2e-7);
    vx(tbrs) = - vx(tbrs);
    x(tbrs) = 1.2e-7;    
        
    % bottom box, top
    bbt = (x > 0.8e-7 & x < 1.2e-7 & y < 0.4e-7 & yp > 0.4e-7);
    vy(bbt) = -vy(bbt);
    y(bbt) = 0.4e-7;

    % top box, bottom
    tbb = (x > 0.8e-7 & x < 1.2e-7 & y > 0.6e-7 & yp < 0.6e-7);
    vy(tbb) = - vy(tbb);
    y(tbb) = 0.6e-7;

    % Scattering Code begins
    Pscat = 1 - exp(-dt/tau);

    if (Pscat > rand())
        %assign new v for a particle where Pscat is greater than random value
        vx = randn(np,1)*vth/sqrt(2);
        vy = randn(np,1)*vth/sqrt(2);
    end    
       
    %plotting trajectories 
    plot (x, y, '.r')
    
    xp = x;
    yp = y;
    
    pause(0.1)
end 

%Histogram plot of the speeds
vavg = sqrt(vx.^2 + vy.^2);
figure;
hist(vavg, 50);
