clc;clear;close all;

% simulation x length and time length
max_space_x = 401;
max_space_y = 401;
max_time = 1000;

% consts
mu = pi*4e-7;
elp = 8.85e-12;
c = sqrt(1/mu/elp); 
disp(c)


% params
wl = 523e-9;   %500nm
freq = c/wl;

% units
dx = wl/20;
dy = dx;
dt = 1/c*((1/dx)^2+(1/dy)^2)^-0.5; %Stability Conditions

% E H Conditions
Ez=zeros(max_space_x+1,max_space_y+1);
Hy=zeros(max_space_x+1,max_space_y+1);
Hx=zeros(max_space_x+1,max_space_y+1);
 
Eeta=dt/dx/elp;     %updating variable   
Heta=dt/dx/mu;      %updating variable

% Gaussian Source
t0=200;                                                   %Gaussian  start point x
ty0=200;                                                  %Gaussian  start point y
spread=15;                                                %Gaussian  width
for y = 1:max_space_y+1
    r_2 = (t0-(2:max_space_x)).^2 + (ty0-y)^2;
    Ez(2:max_space_x,y)= exp(-r_2/spread^2);  %Gaussian  function
end
Z = ((elp/mu)^0.5); % Impedance
dir_mult = 1;
travel_right = true;
if travel_right
    t0 = t0-1;
    dir_mult=-1;
end
%for y = 1:max_space_y+1
%    Hy(1:max_space_x-1,y)= dir_mult*Z*exp(-(t0-(1:max_space_x-1)).^2/spread^2);  %Gaussian  function
%end

% Plot initial frame 
%plot(Ez);
%hold on;
%plot(Hy/Z);
%hold off;
%axis([1 max_space_x -1 1])
imagesc(Ez)
axis equal tight;
title(['1D-FDTD Right & Left   ','t = 0'],'FontSize',18);
pause
% Run simulation

for n=1:max_time
    %t=0.5
    Hy(1:max_space_x-1,:) = Hy(1:max_space_x-1,:) + (dt/(dx*mu))*(Ez(2:max_space_x,:)-Ez(1:max_space_x-1,:)); %Update Hy by change in Ez across X
    Hx(:,1:max_space_y-1) = Hx(:,1:max_space_y-1) - (dt/(dx*mu))*(Ez(:,2:max_space_y)-Ez(:,1:max_space_y-1)); %Update Hx by change in Ez across Y

    %t=1
    Ez(2:max_space_x,:) = Ez(2:max_space_x,:) + (dt/(dx*elp))*(Hy(2:max_space_x,:)-Hy(1:max_space_x-1,:)); %Update Ez by change in Hy across X
    Ez(:,2:max_space_y) = Ez(:,2:max_space_y) - (dt/(dx*elp))*(Hx(:,2:max_space_x)-Hx(:,1:max_space_x-1)); %Update Ez by change in Hx across Y
      
     
    imagesc(Ez)
    axis equal tight;

    
    %plot(Ez(:,1));
    %hold on;
    %plot(Hy(:,1)/Z);
    %hold off;
    %axis([1 max_space_x -1 1])
    title(['1D-FDTD Right & Left   ','t = ',num2str(n)],'FontSize',18);
    pause(0.001)
      
    
end 