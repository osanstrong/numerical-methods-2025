%Copied from class example
clc;close all;clear;
%map size&time
max_time=500;
original_max_space=401;
original_spread=15;
original_t0=100;
infinite_space=false;%Whether to keep on shifting it to the right/left so the wave never actually stops at edges
offset=0;%How far light should have traveled in that time, compared to the current grid. Whenever this is >1, subtract 1 and move the grid to the right 

net_discrete_offset=0; %How many discrete steps the window has been offset in total

plot_each_frame=false; %Whether to plot each frame of the simulation;
plot_end=false; %Whether to plot the end


resolutions = [1,2,3,4,5,6,7,8,9,10]; % For hw for now
resolutions = [1,2,4,8,16,32,64,128];
resolutions = [1,2,5,7];
%resolutions = [1,2,4];
%resolutions = [1,0.5,0.25,0.125,0.06125];
%resolution_factor = resolutions(1);
resolution_errors = zeros(size(resolutions));



%Function to write a Gaussian source
function gauss = gaussian_source(x_points, peak, spread)
    gauss = exp(-(peak-(x_points)).^2/spread^2);
end

%Function to test RMS error
function error = RMS_error(array1, array2)
    if (size(array2,2) == size(array1,1))
        array2 = transpose(array2);
    end
    %disp(array1)
    %disp(array2)
    %disp("\n\n\n\n\n ")
    N = size(array1, 1);
    diff = array1-array2;
    %disp(diff)
    diff2 = diff.*diff;
    %disp(diff2)
    sum_err = sum(diff2,1);
    %disp(sum_err)
    error = (sum_err/N)^0.5;
end

%Main loop to test different cases for resolutions or timescales etc
%for i = 1:size(resolutions,1)
%    resolution_factor = resolutions(i);
i=0;
for resolution_factor = resolutions
    i=i+1;
%for i = 1:1
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
    dx = dx * resolution_factor;
    max_space = floor(original_max_space/resolution_factor);
    disp("Max space: "+max_space)
    
    
    
    
    %E H Conditions
    Ez=zeros(max_space+1,1);
    Hy=zeros(max_space+1,1);
     
    Eeta=dt/dx/elp;     %updating variable   
    Heta=dt/dx/mu;      %updating variable
    
    %Gaussian Source
    t0=original_t0 / resolution_factor;                                                   %Gaussian  start point
    spread=original_spread / resolution_factor;                                                %Gaussian  width
    
    
    
    %Ez starting field
    Ez(2:max_space,1)= gaussian_source(2:max_space, t0, spread);  %Gaussian  function
    
    Z = ((elp/mu)^0.5); % Impedance
    dir_mult = 1;
    travel_right = true;
    if travel_right
        t0 = t0-1;
        dir_mult=-1;
    end
    
    %Hy starting field
    Hy(1:max_space-1,1)= dir_mult*Z*exp(-(t0-(1:max_space-1)).^2/spread^2);  %Gaussian  function
    
    
    %Boundary Condition
    x_conductor = floor(350/resolution_factor);
    
    %Time loop
    current_timestamp_idx = 1;
    for n=1:max_time
        
        
        Hy(1:max_space-1)=Hy(1:max_space-1)+Heta*(-Ez(1:max_space-1)+Ez(2:max_space)); %t=0.5          %equation
        
        
        Ez(2:max_space)=Ez(2:max_space)+Eeta*(Hy(2:max_space)-Hy(1:max_space-1));     %t=1              %equation 
        
        Ez(x_conductor:max_space) = 0;
        
        

        %Update offset (if infinite) and compare to theoretical
        if infinite_space
            light_change = c*dt / dx; %How far right (in x steps) light will have traveled in a time step
            offset = offset + light_change;
            if offset > 1 %Shift everything to the right
                offset = offset-1;
                Ez(1:max_space-1) = Ez(2:max_space);
                Ez(max_space) = 0;
                Hy(1:max_space-1) = Hy(2:max_space);
                Hy(max_space) = 0;
                net_discrete_offset=net_discrete_offset+1;
            end
        end
        
      
        
        if plot_each_frame
            %Plot
            plot(Ez);
            axis([1 max_space -1 1])
            title(['1D-FDTD Right & Left   ','t = ',num2str(n)],'FontSize',18);
            pause(0.001)
        end
    end 
    
    final_light_offset = c*((max_time-1)*dt)/dx;
    rel_offset = final_light_offset-net_discrete_offset;

    theoretical = zeros(size(Ez,1),1);
    final_center = t0+rel_offset;
    final_intensity = 1;
    if final_center > x_conductor
        final_center = x_conductor - (final_center-x_conductor);
        final_intensity = final_intensity*-1;
    end
    disp("Final light center: "+(final_center))

  
    theoretical(1:max_space) = final_intensity*gaussian_source(1:max_space, final_center, spread);
    error = RMS_error(Ez, theoretical);
    disp("Error: "+error)
    resolution_errors(i)=error;
    
    if plot_end
        %Plot last frame
        plot(Ez);
        hold on;
        plot(theoretical);
        hold off;
        axis([1 max_space -1 1])
        title(['1D-FDTD Right & Left   ','t = ',num2str(n)],'FontSize',18);
        pause
    end
end

for e = resolution_errors
    disp(e)
end

plot(resolutions, resolution_errors)
title(['dx of increasing scale vs error'],'FontSize',18);

%}
