%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File:     createTetrahedralSpringMesh.h
% Creator:  Kory Melton and Ian Besse
% Date:     11/14/17
% Purpose:  Creates a tetrahedral spring mesh using the
%           delaunayTriangulation method already developed in
%           in matlab.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CreateTetrahedralMesh

tic % record startTime
m = 1; % steps through P
b = 12; % used to create circle for eye
step = 4; % used to control how many points there are
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Number of points for each step value
% Step 1: 7,153 points
% Step 2: 925 points
% Step 3: 257 points
% Step 4: 123 points
% Step 5: 56 points (not as useful)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fixedPoints = [];
frontPoints = [];
cutoff = 11.8;
frontCutoff = -7;

for i = -b:step:b
    for j = -b:step:b
        for k = -b:step:b
            P (m,:) = [i,j,k];
            m = m + 1;
        end
    end
end

% Grab only the points less than b away from origin
Ind = find (sum ((P.*P)') <= b^2);
P = P (Ind,:);

% Create delaunayTriangulation get numPoints
DT = delaunayTriangulation (P (:,1),P (:,2),P (:,3));
[numPoints, ~] = size (DT.Points);

% Get the distances for each point from center (origin)
for i = 1:numPoints
    distances (i) = sqrt (DT.Points (i,1)^2 + DT.Points (i,2)^2 + ...
                       DT.Points (i,3)^2);
end

% determine which points are fixed
% keep points who are farther away than cutoff and whose x coord is
% greater than frontCutoff
k = boundary(DT.Points);
bIndex = unique(k);

for i = 1:length (bIndex)
      if (DT.Points(bIndex (i),1) > frontCutoff)
          fixedPoints = [fixedPoints, bIndex(i)];  
      else
          frontPoints = [frontPoints, bIndex(i)];
      end
end

save ('MeshInit.mat', 'DT', 'fixedPoints', 'frontPoints');
toc % record endtime