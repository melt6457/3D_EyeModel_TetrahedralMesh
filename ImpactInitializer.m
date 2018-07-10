function [Ind,v] = ImpactInitializer(theta)
%This program takes as input the angle (from straight on) of the impact
%force vector. It outputs the indices (from DT.Points) of the points that
%should be given an initial velocity in order to model the effects of an
%impact from this angle.

%Convert to radians (this assumes the input is in degrees)
theta = deg2rad(theta);

%Load the MeshInit data so that the set of front (movable points on the
%boundar) can be extracted.
load('Data/MeshInit.mat')

%Define M to be a p-by-4 matrix in which each of the p rows gives the
%coordinates of one of the front points in our mesh concatenated with the
%index (from DT.Points) of that point.
M = [DT.Points(frontPoints,:,:),frontPoints'];

% %Define PlanePoints to be a q-by-3 matrix in which each of the q rows gives
% %the coordinates of a front point that is in the xy-plane. Due to symmetry,
% %we need only consider impact angles within one plane that contains the
% %pupil and the center of the eye.
% PlanePoints = M(M(:,3)==0,:);

%ImpPoint is the point of impact (on the external surface of the eye) of
%the projectile that is approaching the eye at an angle of theta.
ImpPoint = -12*[cos(theta),sin(theta),0];

%k is the index (i.e. row) of the point represented in M that is nearest to
%the impact point ImpPoint.
k = dsearchn(M(:,1:3),ImpPoint);

%P is the row of M that represents the point in DT.Points that is nearest
%to the impact point ImpPoint.
P = M(k,:);

%Dist is a vector of distances between all the points in M and the point of
%M that is nearest the impact.
Dist = sqrt(sum((M(:,1:3) - P(:,1:3)).^2,2));

%ImpactZone is an r-by-4 matrix in which each of the r rows gives the
%coordinates of the points in the vicinity of the impact concatenated with
%the index (from DT.Points) of that point.
ImpactZone = M(Dist<=12,:);

%Ind is a vector of the indices of the points in DT.Points that should
%be given an initial velocity due to an impact along this angle theta.
Ind = [P(4);ImpactZone(:,4)];

%v is the vector that gives the direction of the initial velocity of these
%points.
v = -P(1:3);
