clc;clear;close all;

%map size&time
MT=500; % Max time
MX=401; % Max X

%Universal parameters
mu0 = pi*4e-7;
ep0 = 8.85e-12;
c = sqrt(1/mu0/ep0);

%Simulation parameters
wavelength = 523e-9;  
freq = c/wavelength;
dx = wavelength/20;
dt = (1/c)*(2*dx^-2)^-0.5; %to meet stability condition

%Update factors and secondary variables
Eeta=dt/dx/ep0;
Heta=dt/dx/mu0;
Z = (ep0/mu0)^0.5; %Impedance of free space

%EH Fields
Ez=zeros(MX+1,1);
Hy=zeros(MX+1,1);

%Initial starting fields, gaussian source
direction=1; %1 for right
x0_E=100;
x0_H=x0_E-1;
spread=10;

Ez(2:MX,1)= exp(-(x0_E-(2:MX)).^2/spread^2); 
Hy(1:MX-1,1)= -direction*Z*exp(-(x0_H-(1:MX-1)).^2/spread^2); 


%Boundary Condition
x_conductor = 350;

%Main time loop
for n=1:MT
    %Update Equations
    Hy(1:MX-1)=Hy(1:MX-1)+Heta.*(-Ez(1:MX-1)+Ez(2:MX));%t offset by 0.5
    Ez(2:MX-1)=Ez(2:MX-1)+Eeta.*(Hy(2:MX-1)-Hy(1:MX-2)); %x offset by 0.5

    %Boundary Condition
    Ez(x_conductor:MX)=0;

    plot(Ez);
    axis([1 MX -1 1])
    title(['1D-FDTD impinging upon conductor ','t = ',num2str(n)],'FontSize',18);
    pause(0.001)
end