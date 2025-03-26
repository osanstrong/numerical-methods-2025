clear; clc; close all;

boulderPath = 'boulder.tmj';
bldJSON = jsondecode(fileread(boulderPath));
bldPointsRaw = bldJSON.layers.objects.polygon;

numPoints = size(bldPointsRaw,1);
bldMat = zeros(numPoints,2);
for i = 1:numPoints
    bldMat(i,1) = bldPointsRaw(i).x;
    bldMat(i,2) = bldPointsRaw(i).y;
end

rho = 2.5;
g = 1500;
dt = 0.01;
tMax = 10;

vx = 30;
vy = 0; 
ax = 0;
omega = 4;
alpha = 0;

restitution = 0.5; 
airFriction = 0.9;

groundFriction = 0.8; % Both static and kinetic

torqueFactor = 15000;

numTriangles = numPoints - 2;
triMass = zeros(numTriangles,1);
triCM = zeros(numTriangles,2);

for i = 1:numTriangles
    idxs = [1, i+1, i+2];
    tri = bldMat(idxs, 1:2);
    CM = sum(tri) / 3;
    triCM(i, :) = CM;
    triMass(i) = polyarea(tri(:,1), tri(:,2)) * rho;
end

sumMass = sum(triMass);
totalCM = sum(triCM .* triMass, 1) / sumMass;

% Utility function: Moment of inertia of a triangle around its first point, where tri is a 3x2
% matrix and rho is its density
% Using from https://physics.stackexchange.com/questions/493736/moment-of-inertia-for-an-arbitrary-polygon
function I = triangleMomentIntertia(tri,rho)
    % Subtracting the first point, so it's in coordinates around the first
    % point as (0,0)
    xa=tri(2,1)-tri(1,1);
    ya=tri(2,2)-tri(1,2);
    xb=tri(3,1)-tri(1,1);
    yb=tri(3,2)-tri(1,2);
    I = (xa*yb - xb*ya)*((xa-xb)^2 + (ya-yb)^2 + xa*xb + ya*yb) * (rho/12); 
end
% Assumes counterclockwise wound, if clockwise it will just return -1*I
function I = polygonMomentInertia(poly,rho,CM)
    I = 0;
    numSides = size(poly,1);
    for i = 1:1:numSides
        i2 = i+1;
        if i2 > numSides
            i2=1;
        end
        tri = [CM(1),poly(i,1),poly(i2,1);
               CM(2),poly(i,2),poly(i2,2)];
        tri = transpose(tri); 
        I=I+triangleMomentIntertia(tri,rho);
    end
end

sumradii = 0;
for i = 1:numPoints 
    sumradii = sumradii + length(bldMat(i,:)-totalCM);
end

totalI = abs(polygonMomentInertia(bldMat, rho, totalCM));
disp("I: "+totalI)

x0 = 500;
y0 = 1000;
bldMat = bldMat - totalCM + [x0, y0];
totalCM = [x0, y0];

fig = figure;
fig.Position = [100, 100, 800, 600]; 
hold on;
axis equal;
grid on;
xlim([-500, 2000]);
ylim([0, 1200]);

bc = [69, 78, 97] / 256;
red = [1, 0, 0];

polyHandle = fill(bldMat(:,1), bldMat(:,2), bc, "LineStyle", "none");
cmHandle = rectangle("FaceColor", red, "Position", [totalCM(1)-5, totalCM(2)-5, 10, 10], "Curvature", [1,1]);

t = 0;
collision = false;
while t < tMax  

    vx = (vx + ax * dt)*airFriction;
    vy = vy - g * dt;


    bldMatBefore = bldMat(:,:);



    omega = omega + alpha * dt;

    yNew = totalCM(2) + vy * dt;
    xNew = totalCM(1) + vx * dt;
    dy = yNew - totalCM(2);
    dx = xNew - totalCM(1);

    theta = omega * dt;
    R = [cos(theta), -sin(theta); sin(theta), cos(theta)];
    bldMat = (bldMat - totalCM) * R' + totalCM;

    bldMat(:,1) = bldMat(:,1) + dx;
    bldMat(:,2) = bldMat(:,2) + dy;
    totalCM(1) = xNew;
    totalCM(2) = yNew;

    min_point = bldMat(1,:);
    min_point_before = bldMatBefore(1,:);
    for i = 2:numPoints 
        if bldMat(i,2) < min_point(1,2)
            min_point = bldMat(i,:);
            min_point_before = bldMatBefore(i,:);
        end
    end
    
    

    if min_point(1,2) <= 0
        ax=0;
        collision = true;
        deltaY = min(bldMat(:,2));
        bldMat(:,2) = bldMat(:,2) - deltaY;
        totalCM(2) = totalCM(2) - deltaY;

        vy_before = vy;
        vy = -vy * restitution;
        
        impulse_y = abs(vy - vy_before)*sumMass;

        %max_friction_x_impulse = impulse_y*groundFriction; %TODO: Use this to check if the full slip should be negated. For now, we assume it just can't slip at all
        %arm_before = min_point_before(1,:) - totalCM;
        %arm_rot = arm_before * [0,-1;1,0];  % Radius rotated 90⁰ counterclockwise, i.e. direction of rotational velocity
        %impact_vel = [vx,vy] + arm_rot * omega;
        %new_impact_vel = [0,0];
        %d_impact_vel = new_impact_vel - impact_vel;

        slip_x = min_point(1,1)-min_point_before(1,1);
        disp(slip_x)
        bldMat(:,1) = bldMat(:,1) - slip_x;
        totalCM(1) = totalCM(1) - slip_x;

        
        %impactX = totalCM(1);

        %impactX = min(bldMat(bldMat(:,2) <= 0, 1));
        impactX = min_point(1,1);
        
        armx = impactX - totalCM(1);
        

        %disp("Arm: "+arm)
        rot_impulse = impulse_y * armx * 0.35;
        %disp("rot_impulse: " + rot_impulse)
        disp("omega before: " +omega)
        omega = omega * 0.75 + (rot_impulse/totalI);
        disp("omega after: " +omega)

        %rotateFactor = -( impactX - totalCM(1) ) / torqueFactor * omega;
        %if impactX < totalCM(1)
        %    omega = omega - abs(vy) * rotateFactor; 
        %else
        %    omega = omega + abs(vy) * rotateFactor; 
        %end
    end

    set(polyHandle, 'XData', bldMat(:,1), 'YData', bldMat(:,2));
    set(cmHandle, 'Position', [totalCM(1)-5, totalCM(2)-5, 10, 10]);
    drawnow;
    pause(dt);

    t = t + dt;
end