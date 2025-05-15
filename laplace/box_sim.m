clc;clear;close all;

% Group 1 in class May 15 2025

% Grid setup
GRID_W = 100;
GRID_H = 60;

V = zeros(GRID_W,GRID_H);
V_locked = false(GRID_W, GRID_H);

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

V_locked = set_vert_line(V_locked, 1, true);
V_locked = set_vert_line(V_locked, GRID_W, true);
V_locked = set_horz_line(V_locked, 1, true);
V_locked = set_horz_line(V_locked, GRID_H, true);
V_locked = set_circle(V_locked, 0.75*GRID_W, 0.5*GRID_H, 0.33*GRID_H, true);
V_locked = set_triangle(V_locked, 0.4*GRID_W, 0.2*GRID_W, 0.5*GRID_H, 30, true);

V = set_vert_line(V, 1, 20);
V = set_vert_line(V, GRID_W, 5);
V = set_horz_line(V, 1, 30);
V = set_horz_line(V, GRID_H, 10);
V = set_circle(V, 0.75*GRID_W, 0.5*GRID_H, 0.33*GRID_H, 0);
V = set_triangle(V, 0.4*GRID_W, 0.2*GRID_W, 0.5*GRID_H, 30, 100);

%Display color plot of V
imagesc(V_locked');colorbar;
title(['Whether each pixel is locked']); 
xlabel('x','FontSize',20);
ylabel('y','FontSize',20);
set(gca,'FontSize',20);
pbaspect([GRID_W GRID_H 1]);
getframe;

pause();

max_iter = 20000000000;
for n=1:1:max_iter

    V_before = V(:,:);
    for x = 1:1:GRID_W
        for y = 1:1:GRID_H
            if (~V_locked(x,y))
                V(x,y) = (V(x-1,y) + V(x+1,y) + V(x,y-1) + V(x, y+1))/4;
                
            end
        end
    end

    

    

    if (max(abs(V-V_before)) < 0.05)
        break
    end
    %Display color plot of V
    imagesc(V');colorbar;
    title(['Voltage (V)']); 
    xlabel('x','FontSize',20);
    ylabel('y','FontSize',20);
    set(gca,'FontSize',20);
    pbaspect([GRID_W GRID_H 1]);
    getframe;
end

%Display color plot of V
imagesc(V');colorbar;
title(['Voltage (V)']); 
xlabel('x','FontSize',20);
ylabel('y','FontSize',20);
set(gca,'FontSize',20);
pbaspect([GRID_W GRID_H 1]);
getframe;