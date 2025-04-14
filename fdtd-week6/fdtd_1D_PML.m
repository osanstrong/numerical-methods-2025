%Copied from class example
clc;close all;clear;
%map size&time
max_time=500;
max_space=401;

%parameters
mu = pi*4e-7;
elp = 8.85e-12;
c = sqrt(1/mu/elp); 
wavelength = 523e-9;   %500nm
freq = c/wavelength;

%Unit 
dx = wavelength/20;
dy = dx;
dt = 1/c*((1/dx)^2+(1/dy)^2)^-0.5; %Stability Conditions
resolution_factor=1;
dx = dx * resolution_factor;
max_space = floor(max_space/resolution_factor);

%E H Conditions
Ez=zeros(max_space+1,1);
Hy=zeros(max_space+1,1);
 
Eeta=dt/dx/elp;     %updating variable   
Heta=dt/dx/mu;      %updating variable

%Gaussian Source
t0=100 / resolution_factor;                                                   %Gaussian  start point
spread=15 / resolution_factor;                                                %Gaussian  width
%Ez starting field
Ez(2:max_space,1)= exp(-(t0-(2:max_space)).^2/spread^2);  %Gaussian  function

Z = ((elp/mu)^0.5); % Impedance
dir_mult = -1;
travel_right = true;
if travel_right
    t0 = t0-1;
    dir_mult=1;
end

%Hy starting field
Hy(1:max_space-1,1)= dir_mult*Z*exp(-(t0-(1:max_space-1)).^2/spread^2);  %Gaussian  function


%Boundary Condition
x_conductor = 350;
x_pml = 315; %First x of the pml boundary, then the rest are PML absorption

% PML layers
M = max_space;
ep(1:M)=elp; % permitivity array
mu_ar(1:M)=mu; % pemeability array
%---- PML absorbing boundary condition----
sigma(1:M)=0; % initialize conductivity array
d=M-x_pml; % width of PML layer
m=3; % polynomial order for grading sigma array (pp 292, Taflove)
neta=sqrt(mu/elp); 
R=1e-8; % required reflectivity
sigma_max=-(m+1)*log(R)/(2*neta*d*dx);
Pright=((1:d+1)./d).^m*sigma_max;
sigma(M-d:M)=Pright; % lossy conductivity profile
%sigma(1:d+1)=fliplr(Pright); %Left side of the boundary, skip for now
sigma_star(1:M)=sigma.*mu./elp; % Eq 7.8 Taflove, pp 275
%------------- PML constants ----------------------------------------%
A=((mu_ar-0.5*dt*sigma_star)./(mu_ar+0.5*dt*sigma_star)); 
B=(dt/dx)./(mu_ar+0.5*dt*sigma_star);                          
C=((ep-0.5*dt*sigma)./(ep+0.5*dt*sigma)); 
D=(dt/dx)./(ep+0.5*dt*sigma); 
A=transpose(A);
B=transpose(B);
C=transpose(C);
D=transpose(D);




for n=1:max_time
   
    Hy(1:M-1)=A(1:M-1).*Hy(1:M-1)-B(1:M-1).*(-Ez(1:M-1)+Ez(2:M)); %t=0.5          %equation
   
    %Hy(x_pml:x_pml+n_pml)=Hy(x_pml-1:x_pml+n_pml-1);

    Ez(2:M-1)=C(2:M-1).*Ez(2:M-1)-D(2:M-1).*(Hy(2:M-1)-Hy(1:M-2));     %t=1              %equation 
    
    %Ez(x_pml+1:x_pml+n_pml+1)=Ez(x_pml:x_pml+n_pml);

    %Ez(x_conductor:max_space) = 0;
    
    plot(Ez);
    axis([1 max_space -1 1])
    title(['1D-FDTD Right & Left   ','t = ',num2str(n)],'FontSize',18);
    pause(0.001)

end 



%}
