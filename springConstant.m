%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File:     springConstant.m
% Creator:  Kory Melton and Ian Besse
% Date:     11/14/17
% Purpose:  To create the spring constants for the tetrahedral
%           mesh (created with CreateTetrahedralMesh.m)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function springConstant()

tic % startTime
% Load in correct data
clear % clear previous data
load('meshInit.mat');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Info Stored in meshInit.mat (from CreateTetrahedralMesh.m)
%
% DT:           The delaunayTriangulation used to create the 
%               tetrahedral mesh
% fixedPoints:  the points in the eye that will be fixed
%               (i.e., they don't move)
% frontPoints:  the points in the eye used for initial forces
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 3D Young's Modulus
% we used 42 kPa from Hung Paper converted to N/mm^2
E3 = .000042; 
% List of points in DT (x,y,z)
Points = DT.Points;
[numPoints, ~] = size(Points);
% List of tetrahedra (p1,p2,p3,p4)
Cl = DT.ConnectivityList;
% List of Edges in DT (p1,p2)
E = edges(DT); 
% List of tetrahedron attached to each edge
A = edgeAttachments(DT, E); 
% M is the number of edges
[numEdges,~] = size(E);
% Initialize the Matrix of Natural Resting Lengths of Springs
N = zeros(numPoints, numPoints);
% Initialize the Matrix of Spring Constants
K = zeros(numPoints, numPoints);

% step through all of the edges
for i = 1:numEdges
    % find the tetrahedron for each edge and save the size
    % T is a vector with the numbers of all the tetrahedra
    % which share edge i
    T = A{i};
    numTetra = length(T);
    
    % use the two points to get length of the edge
    % E(i,1) - E(i,2) gives the node number at the endpoint
    c1 = Points(E(i,1),:);
    c2 = Points(E(i,2),:);
    c = norm(c1 - c2);
    % save the initial length of the edge
    N(E(i,1), E(i,2)) = c;
    N(E(i,2), E(i,1)) = c;
%     D(E(i,1), E(i,2),:) = c2 - c1;
%     D(E(i,2), E(i,1),:) = c1 - c2;
    
    % initialize the spring constant
    Ki = 0;
    
    % step through the tetrahedron attached to the edge
    for j = 1:numTetra
        % find/save the 4 points
        PointsVec = Cl(T(j),:);
        p1 = Points(PointsVec(1),:);
        p2 = Points(PointsVec(2),:);
        p3 = Points(PointsVec(3),:);
        p4 = Points(PointsVec(4),:);
        
        % get the three vectors for the points
        v1 = p1 - p4;
        v2 = p2 - p4;
        v3 = p3 - p4;
        
        % use the vectors to get the area of the tetrahedron
        % by means of 1/6 volume of parallelpiped
        Vol = abs((1/6)*v1*cross(v2,v3)');
        
        % Find the spring constant value for this tetrahedron
        %(E3*Vol)/c^2;
        Knew = (E3*Vol)/c^2;
        % add the spring constant
        Ki = Ki + Knew;
    end
    
    % Add Ki to the list of spring constants in the K Matrix
    K(E(i,1), E(i,2)) = Ki;
    K(E(i,2), E(i,1)) = Ki;
end

% save the values to tetrahedralMesh2.mat
save('EdgeInit.mat', 'K', 'N', 'Points', 'fixedPoints', 'frontPoints');
toc % record endTime