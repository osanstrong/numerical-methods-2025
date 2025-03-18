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

% Utility function: Moment of inertia of a triangle, where tri is a 3x2
% matrix and rho is its density
function I = triangleMomentIntertia(tri,rho)
    I = 0; %TODO: implement this
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

% Display Boulder
fig=gcf;
fig.Position(3:4)=[400,400];
hold on;

bc=[69, 78, 97];
bc=bc./256;
red=[1,0,0];

fill(bldMat(:,1),bldMat(:,2),bc,"LineStyle","none");
rad=10;
rectangle("FaceColor",red,"Position",[totalCM(1)-0.5*rad,totalCM(2)-0.5*rad,rad,rad],"Curvature",[1,1])

%disp(size(bldPoints,1));
