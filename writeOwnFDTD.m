% Modeled after 1D Example FDTD algorithm given by Chien Chao, National Taiwan University, 3/13/2025 Numerical Methods class

clc;clear;close all;

% simulation x length and time length
max_space = 401;
max_time = 150;

% consts
mu = pi*4e-7;
elp = 8.85e-12;
c = sqrt(1/mu/elp); 

% params
wl = 523e-9;   %500nm
freq = c/wl;

% units
dx = wl/20;
dy = dx;
dt = 1/c*((1/dx)^2+(1/dy)^2)^-0.5; %Stability Conditions

% E H Conditions
Ez=zeros(max_space+1,1);
Hy=zeros(max_space+1,1);
 
Eeta=dt/dx/elp;     %updating variable   
Heta=dt/dx/mu;      %updating variable

% Gaussian Source
t0=200;                                                   %Gaussian  start point
spread=15;                                                %Gaussian  width
Ez(2:max_space,1)= exp(-(t0-(2:max_space)).^2/spread^2);  %Gaussian  function


% Run simulation

for n=1:max_time
    %t=0.5
    Hy(1:max_space-1) = Hy(1:max_space-1) + (dt/(dx*mu))*(Ez(2:max_space)-Ez(1:max_space-1));

    %t=1
    Ez(2:max_space) = Ez(2:max_space) + (dt/(dx*elp))*(Hy(2:max_space)-Hy(1:max_space-1));
      
      
      plot(Ez);
      axis([1 max_space -1 1])
      title(['1D-FDTD Right & Left   ','t = ',num2str(n)],'FontSize',18);
      pause(0.001)
      
    
end 
