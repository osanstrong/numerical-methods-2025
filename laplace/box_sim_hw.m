clc;clear;close all;
load("box_sim_results.mat")
% Group 1 in class May 15 2025

% Grid setup
GRID_W = 280;
GRID_H = 200;

TRI_LEFT = 40;
TRI_RIGHT = TRI_LEFT + 80;
CIRC_RAD = 40;
CIRC_X = GRID_W - 40 - CIRC_RAD;
MID_Y = GRID_H / 2;

%V = zeros(GRID_W,GRID_H);
%V_locked = false(GRID_W, GRID_H);

function grid =  set_vert_line(grid, x, val)
    grid(x,:) = val;
end
function grid =  set_horz_line(grid, y, val)
    grid(:,y) = val;
end
function grid = set_circle(grid, x_0, y_0, r, val)
    sz = size(grid);
    for x = 1:1:sz(1)
        for y = 1:1:sz(2)
            if ((x-x_0)^2 + (y-y_0)^2) < r^2
                grid(x,y) = val;
            end
        end
    end
end
% Assumes theta is angle from top to bottom, not the angle relative to
% horizontal. (I.e. twice the angle from horizontal)
function grid = set_triangle(grid, x_0, x_1, y, theta_deg, val)
    theta_horz = theta_deg/2;
    slope = tand(theta_horz);
    for x = x_0:sign(x_1-x_0):x_1
        dy = floor(slope*abs(x-x_0));
        grid(x,y-dy:y+dy) = val;
    end
end

angles = [10 20 30 40 50 60];
num_angles = size(angles);
num_angles = num_angles(2);
results(1:num_angles) = {zeros(GRID_W, GRID_H)};

for i = 1:1:num_angles
    wedge_angle_deg = angles(i);
    
    %Mark which pixels shouldn't be changed (i.e. the boundary conditions)
    V_locked = false(GRID_W, GRID_H);
    V_locked = set_vert_line(V_locked, 1, true);
    V_locked = set_vert_line(V_locked, GRID_W, true);
    V_locked = set_horz_line(V_locked, 1, true);
    V_locked = set_horz_line(V_locked, GRID_H, true);
    V_locked = set_circle(V_locked, CIRC_X, MID_Y, CIRC_RAD, true);
    V_locked = set_triangle(V_locked, TRI_RIGHT, TRI_LEFT, MID_Y, wedge_angle_deg, true);
    
    %Set initial conditions, in this case at boundaries
    V = results{i};
    V = set_vert_line(V, 1, 20);
    V = set_vert_line(V, GRID_W, 5);
    V = set_horz_line(V, 1, 30);
    V = set_horz_line(V, GRID_H, 10);
    V = set_circle(V, CIRC_X, MID_Y, CIRC_RAD, 0);
    V = set_triangle(V, TRI_RIGHT, TRI_LEFT, MID_Y, wedge_angle_deg, 100);
    
    
    %Display color plot of V
    imagesc(V_locked');colorbar;
    title(['Whether each pixel is locked']); 
    xlabel('x','FontSize',20);
    ylabel('y','FontSize',20);
    set(gca,'FontSize',20);
    pbaspect([GRID_W GRID_H 1]);
    getframe;
    
    if i == 1
        pause();
    end
    
    max_iter = 3000;
    for n=1:1:max_iter
    
        V_before = V(:,:);
        for x = 1:1:GRID_W
            for y = 1:1:GRID_H
                if (~V_locked(x,y))
                    V(x,y) = (V(x-1,y) + V(x+1,y) + V(x,y-1) + V(x, y+1))/4;
                    
                end
            end
        end
    
        
    
        
    
        if (max(abs(V-V_before)) < 0.01)
            break
        end
        Display color plot of V
        imagesc(V');colorbar;
        title(['Voltage (V), theta ',num2str(wedge_angle_deg), ' (deg) ', int2str(n),' iterations']); 
        xlabel('x','FontSize',20);
        ylabel('y','FontSize',20);
        set(gca,'FontSize',20);
        pbaspect([GRID_W GRID_H 1]);
        getframe;
    end

    %results{i} = V;
    Display color plot of V
    imagesc(V');colorbar;
    title(['Voltage (V), theta ',num2str(wedge_angle_deg), ' (deg) ', int2str(n),' iterations']); 
    xlabel('x','FontSize',20);
    ylabel('y','FontSize',20);
    set(gca,'FontSize',20);
    pbaspect([GRID_W GRID_H 1]);
    getframe;
    pause();
end

for i = 1:1:num_angles
    ang_deg = angles(i);
    V = results{i}';
    x = 1:1:GRID_W;
    y = 1:1:GRID_H;
    %[X, Y] = meshgrid(x,y);
    %[E_x, E_y] = gradient(V);
    
    quiver(X, Y, E_x,E_y, 0)
    axis equal
    Display color plot of V
    imagesc(results{i}');colorbar;
    title(['E field with theta (degrees) of ',num2str(ang_deg)]); 
    xlabel('x','FontSize',20);
    ylabel('y','FontSize',20);
    set(gca,'FontSize',20);
    pbaspect([GRID_W GRID_H 1]);
    getframe;
    pause();
end

%save("box_sim_results.mat")
%HW: Write a report, vary angle of wedge, plot with quiver the electric
%field, a report and a video
%Full report of simulation due this sunday