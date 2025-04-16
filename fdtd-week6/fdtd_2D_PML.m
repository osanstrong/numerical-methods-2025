% Rewritten following Sathya Ganta's example 
clc;clear;close all;

%Grid Dimensions
MX=100; %Max x
MY=100; %Max y
%Timestep count
MT=500; %Max time, number of timesteps

%Source locations, center and overlapping
center_x=MX/2;
center_y=MY/2;
left_x=MX/3;
right_x=MX-left_x;

%Physical parameters of free space
ep0=(1/(36*pi))*1e-9; %Vacuum permutivity, F/m
mu0=4*pi*1e-7; %Vacuum permeability, N/A²
c=3e+8; %Speed of light, m/s

%Set arbitrary dx/dy/etc of uniform 1 micron
dxy=1e-6;

%Courant Stability Condition to determine timestep
S=2^-0.5;
dt=S*dxy/c;

%Field matrices (I actually don't know why Ez is split up into x and y we
%should find out
Ez=zeros(MX,MY);
Ezx=zeros(MX,MY);
Ezy=zeros(MX,MY);
Hx=zeros(MX,MY);
Hy=zeros(MX,MY);

%Matrices of permittivity and permeability constants
ep=ep0*ones(MX,MY);
mu=mu0*ones(MX,MY);

%Electric conductivity matrices for x and y directions
sigx=zeros(MX,MY);
sigy=zeros(MX,MY);

%Perfectly matched layer boundary design
%Reference: IIT Madras team referencing http://dougneubauer.com/wp-content/uploads/wdata/yee2dpml1/yee2d_c.txt
%(An adaptation of 2-D FDTD TE code of Dr. Susan Hagness)

%Boundary width of PML in all directions
bound_width=25;

%Required reflection coefficient
max_refl=1e-6;

%Polynomial order for the σ model
sigma_order=6;
%Polynomial σ model
sigma_max=(-log10(max_refl)*(sigma_order+1)*ep0*c)/(2*bound_width*dxy);
bound_fact_1=((ep(MX/2,bound_width)/ep0)*sigma_max)/((bound_width^sigma_order)*(sigma_order+1));
bound_fact_2=((ep(MX/2,MY-bound_width)/ep0)*sigma_max)/((bound_width^sigma_order)*(sigma_order+1));
bound_fact_3=((ep(bound_width,MY/2)/ep0)*sigma_max)/((bound_width^sigma_order)*(sigma_order+1));
bound_fact_4=((ep(MX-bound_width,MY/2)/ep0)*sigma_max)/((bound_width^sigma_order)*(sigma_order+1));
x=0:1:bound_width;
for i=1:1:MX
    sigx(i,bound_width+1:-1:1)=bound_fact_1*((x+0.5*ones(1,bound_width+1)).^(sigma_order+1)-(x-0.5*[0 ones(1,bound_width)]).^(sigma_order+1));
    sigx(i,MY-bound_width:1:MY)=bound_fact_2*((x+0.5*ones(1,bound_width+1)).^(sigma_order+1)-(x-0.5*[0 ones(1,bound_width)]).^(sigma_order+1));
end
for i=1:1:MY
    sigy(bound_width+1:-1:1,i)=bound_fact_3*((x+0.5*ones(1,bound_width+1)).^(sigma_order+1)-(x-0.5*[0 ones(1,bound_width)]).^(sigma_order+1))';
    sigy(MX-bound_width:1:MX,i)=bound_fact_4*((x+0.5*ones(1,bound_width+1)).^(sigma_order+1)-(x-0.5*[0 ones(1,bound_width)]).^(sigma_order+1))';
end

%Magnetic conductivity matrix using PML condition, also split into x and y
%directions
sig_starx=(sigx.*mu)./ep;
sig_stary=(sigy.*mu)./ep;

%Multiplicartion factor matrices for updating the H matrices instead of
%doing so every frame
G=((mu-0.5*dt*sig_starx)./(mu+0.5*dt*sig_starx));
H=(dt/dxy)./(mu+0.5*dt*sig_starx);
A=((mu-0.5*dt*sig_stary)./(mu+0.5*dt*sig_stary));
B=(dt/dxy)./(mu+0.5*dt*sig_stary);

%Multiplication factor matrices for E matrix
C=((ep-0.5*dt*sigx)./(ep+0.5*dt*sigx));
D=(dt/dxy)./(ep+0.5*dt*sigx);
E=((ep-0.5*dt*sigy)./(ep+0.5*dt*sigy));
F=(dt/dxy)./(ep+0.5*dt*sigy);


%And now, the moment you've all been waiting for: The time loop begins
for n=1:1:MT

    %Update equation for Hx and Hy
    %Vectorized update in one line instead of for loo
    Hy(1:MX-1,1:MY-1)=A(1:MX-1,1:MY-1).*Hy(1:MX-1,1:MY-1)+B(1:MX-1,1:MY-1).*(Ezx(2:MX,1:MY-1)-Ezx(1:MX-1,1:MY-1)+Ezy(2:MX,1:MY-1)-Ezy(1:MX-1,1:MY-1));
    Hx(1:MX-1,1:MY-1)=G(1:MX-1,1:MY-1).*Hx(1:MX-1,1:MY-1)-H(1:MX-1,1:MY-1).*(Ezx(1:MX-1,2:MY)-Ezx(1:MX-1,1:MY-1)+Ezy(1:MX-1,2:MY)-Ezy(1:MX-1,1:MY-1));

    %Likewise update equation for Ez field
    Ezx(2:MX,2:MY)=C(2:MX,2:MY).*Ezx(2:MX,2:MY)+D(2:MX,2:MY).*(-Hx(2:MX,2:MY)+Hx(2:MX,1:MY-1));
    Ezy(2:MX,2:MY)=E(2:MX,2:MY).*Ezy(2:MX,2:MY)+F(2:MX,2:MY).*(Hy(2:MX,2:MY)-Hy(1:MX-1,2:MY));

    

    %Update with gaussian source at start
    gauss=true;
    amp_multiplier=5;
    if gauss
        if n<=42 %Why 42?
            Ezx(center_x,center_y)=amp_multiplier*(10-15*cos(n*pi/20)+6*cos(2*n*pi/20)-cos(3*n*pi/20))/64;
            Ezy(center_x,center_y)=amp_multiplier*(10-15*cos(n*pi/20)+6*cos(2*n*pi/20)-cos(3*n*pi/20))/64;
        else
            Ezx(center_x,center_y)=0;
            Ezy(center_x,center_y)=0;
        end
    else
        Ezx(center_x,center_y)=0.5;
        Ezy(center_x,center_y)=0.5;
    end


    Ez=Ezx+Ezy;

    %Display color plot of Ez
    imagesc(dxy*1e+6*(1:1:MX),(dxy*1e+6*(1:1:MY))',Ez',[-1,1]);colorbar;
    title(['\fontsize{20}Colour-scaled image plot of Ez in a spatial domain with PML boundary and at time = ',num2str(round(n*dt*1e+15)),' fs']); 
    xlabel('x (in um)','FontSize',20);
    ylabel('y (in um)','FontSize',20);
    set(gca,'FontSize',20);
    getframe;
end
