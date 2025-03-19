clear; clc; close all;
boulderPath = 'boulder.tmj';
bldJSON = jsondecode(fileread(boulderPath));
bldPointsRaw = bldJSON.layers.objects.polygon;

% Assumes clockwise bound, convex polygon
numPoints = size(bldPointsRaw,1);
bldMat = zeros(numPoints,2);
for i = 1:1:numPoints
    bldMat(i,1)=bldPointsRaw(i).x;
    bldMat(i,2)=bldPointsRaw(i).y;
end
disp(bldMat)

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

% Calculate physical properties of the boulder with triangulation, assuming:
% Length units: cm
% Density: 2.5 g / cm^3
rho = 2.5;
numTriangles = numPoints-2;
triMass = zeros(numTriangles,1); %Mass of each triangle (1-2-3, then 1-3-4, etc)
triCM = zeros(numTriangles,2); %Center of mass of each triangle
for i = 1:1:numTriangles
    idxs = [1,i+1,i+2];
    disp(idxs)
    tri = bldMat(idxs,1:2);
    disp(tri)
    CM = sum(tri)/3;
    disp(CM)
    triCM(i,1:2)=CM;
    
    triMass(i)=polyarea(tri(1:3,1),tri(1:3,2), 1)*rho;
end
disp(triCM)
disp(triMass)
sumMass=sum(triMass);
totalCM=sum(triCM.*triMass)/sumMass;
disp(totalCM)
bldI = -1*polygonMomentInertia(bldMat, rho, totalCM); %Negative bc our boulder is clockwise wound
disp(bldI)
disp(sumMass)

% Display Boulder
%fig=gcf;
%fig.Position(3:4)=[400,400];
%hold on;

bc=[69, 78, 97];
bc=bc./256;
red=[1,0,0];

fill(bldMat(:,1),bldMat(:,2),bc,"LineStyle","none");
rad=10;
rectangle("FaceColor",red,"Position",[totalCM(1)-0.5*rad,totalCM(2)-0.5*rad,rad,rad],"Curvature",[1,1])

%disp(size(bldPoints,1));

% Simulate the boulder falling

% 2d cross product, by extending vectors to 3d and returning the z
function cross = cross2d(a, b)
    cross = a(1,1)*b(1,2) - a(1,2)*b(1,1);
end

% Calculate velocities of a rigidbody against a horizontal surface (e.g.
% y=0, given:
% v0: Initial velocity (2x1 matrix x, y)
% w0: Initial angular velocity, counterclockwise, in radians
% m: Mass of the body
% I: Moment of inertia of the body
% collisionPoint: Where the body collided with the ground
% CM: The center of mass of the body
% restitution: coefficient of restitution, i.e. elasticity of the
% collision. Currently assumed to be 0 (inelastic)
% Cuttently assumed to be non-slip.
function vels = velsAfterHittingGround(v0, w0, m, I, collisionPoint, CM, restitution)
    %Linear impulse
    dv = v0(:,:)*-1; %This is where we would add effects of it being elastic
    v=v0+dv;
    impulse = dv*m;
    %Rotational impulse
    r=collisionPoint-CM;
    rotImp = cross2d(r,impulse);
    dw = rotImp/I;
    w=w0+dw;
    vels.v = v;
    vels.w = w;
end
% Input point as [x y;], theta in radians
function pos1 = rotatePoint(origin, pos0, theta)
    r0 = pos0 - origin;
    r1 = zeros(1,2);
    r1(1,1) = r0(1,1)*cos(theta) - r0(1,2)*sin(theta);
    r1(1,2) = r0(1,2)*cos(theta) + r0(1,1)*sin(theta);
    pos1 = r1 + origin;
end

% Simulation parameters
pos0 = [400,400]; % starting location of center of mass, in cm, same units as polygon
v0 = [-50,-400]; % cm/s
w0 = 6.28; % rad/s counterclockwise
g = -981; % gravitational acceleration, cm/s^2
groundY = 0; % y of the ground level, currently assumed to be flat
window = [-500 500 -250 750]; % Display window, x to x, then y to y
wx=[window(1), window(1), window(2), window(2)];
wy=[window(3), window(3), window(4), window(4)];
wc = [0,1,1];

sim_time = 20; % seconds
dt = 0.001; % timestep, in seconds
slowdown = 10; %What factor to slow down the simulation by


v = v0;
w = w0;
fallingBoulder = bldMat(:,:);
fallingCM = totalCM(:,:);
toPos0 = pos0-fallingCM;
fallingCM = fallingCM+toPos0;
numPoints = size(fallingBoulder,1);
m = sumMass;
I = bldI;
for i = 1:1:numPoints
    fallingBoulder(i,:) = fallingBoulder(i,:)+toPos0;
end

for t = 0:dt:sim_time
    % Update vel
    v(1,2) = v(1,2) + dt*g;
    % Update positions
    dp = v * dt;
    fallingCM = fallingCM + dp;
    dTheta = w*dt;
    for i = 1:1:numPoints
        fallingBoulder(i,:) = fallingBoulder(i,:)+dp;
        fallingBoulder(i,:) = rotatePoint(fallingCM, fallingBoulder(i,:), dTheta);
    end
    % Check for collisions
    lowestPoint = fallingBoulder(1,:);
    for i = 2:1:numPoints
        pi = fallingBoulder(i,:);
        if pi(1,2) < lowestPoint(1,2)
            lowestPoint = pi;
        end
    end
    if lowestPoint(1,2) < groundY % A collission occurs if the lowest Y is below the ground
        % Bring the boulder back up to the surface, no clipping
        target = lowestPoint(:,:);
        target(1,2) = groundY;
        dp = target - lowestPoint(1,2);
        fallingCM = fallingCM + dp;
        for i = 1:1:numPoints
            fallingBoulder(i,:) = fallingBoulder(i,:)+dp;
        end
        %Apply changes in velocity
        newVel = velsAfterHittingGround(v, w, m, I, target, fallingCM, 0);
        v = newVel.v;
        w = newVel.w;
        disp("Velocity: "+mat2str(v))
        disp("Rot Vel:"+num2str(w))
    end

    % Display
    %axis(window)
    %pbaspect([1.1,1,1])
    clf;
    fig=gcf;
    fig.Position(3:4)=[400,400];
    hold on;
    

    fill(wx, wy, wc,"LineStyle","none")
    fill(fallingBoulder(:,1),fallingBoulder(:,2),bc,"LineStyle","none");
    rad=10;
    rectangle("FaceColor",red,"Position",[fallingCM(1)-0.5*rad,fallingCM(2)-0.5*rad,rad,rad],"Curvature",[1,1])

    title(['Boulder falling','t = ',num2str(t)],'FontSize',18);
    pause(dt*slowdown)
end
