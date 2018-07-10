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

fixedPoints = []; % vector for fixed points which will not move
frontPoints = []; % vector for front points for initial impact
frontCutoff = -7; % cutoff used for points on the front of the eye

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% comment out all but one of these lines to choose a method of
% descritization

P = uniformPoints; % designates descritization using uniform distribution
% P = spheremaker; % designates descritization using spheremaker
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create delaunayTriangulation get numPoints
DT = delaunayTriangulation (P (:,1),P (:,2),P (:,3));

% determine which points are fixed
% keep points who are farther away than cutoff and whose x coord is
% greater than frontCutoff
k = boundary(DT.Points,0);
bIndex = unique(k);

retinaPoints = findRetina (DT.Points);

for i = 1:length (bIndex)
      if (DT.Points(bIndex (i),1) > frontCutoff)
          fixedPoints = [fixedPoints, bIndex(i)];  
      else
          frontPoints = [frontPoints, bIndex(i)];
      end
end

caIndex = findCornealApex (DT.Points);

% save info for future use in 'Data/MeshInit.mat' file
save ('Data/MeshInit.mat', 'DT', 'fixedPoints', 'frontPoints', ...
      'retinaPoints', 'caIndex');
toc % record endtime

function retinaPoints = findRetina (points)

% length of retina from the outside of the eye
retinaOut = 1.436;
% length of retina from the inside of the eye
retinaIn = 1.716;

distCenter = sqrt(sum(points(:,1:3).^2,2)); % distance from center
maxDist = max(distCenter); % maximum distance
distOut = maxDist - distCenter;

retinaPoints = points(distOut<= (retinaIn) & distOut>=retinaOut);

function caIndex = findCornealApex (points)

% coordinates for the front point
frontCoord = [-12 0 0];

% find distance from front of eye
dist = sqrt(sum((frontCoord(:,1:3) - points(:,1:3)).^2,2));
minDist = min(dist); % get the minimum distance from front of the eye
caIndex = find(dist ==minDist);




