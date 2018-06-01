%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File:     pointsMovie.m
% Creator:  Kory Melton and Ian Besse
% Date:     3/1/17
% Purpose:  To create a movie
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function movietrial(P)
% load P
load('Positions.mat');

Fmax = max(max(Fmag));
Fmagp = Fmag/Fmax; % magnitude proportion
[numPoints, ~] = size (Fmag); % get number points

figure
for i = 1:numSteps
%     Cx = Fmagp(:,i).^(1/5);
%     Cy = zeros(numPoints, 1);
%     Cz = (1 - Fmagp(:,i));
%     %Cz = zeros(numPoints, 1);
%     C = [Cx Cy Cz];
    
    scatter3(P(:,i,1),P(:,i,2),P(:,i,3), 7, 'filled')
    axis([-13 13 -13 13 -11 11])
    % view(-10,20) % side view
    drawnow
    F(i) = getframe(gcf);
end

newVid = VideoWriter('movietrial.avi', 'Uncompressed AVI');
open(newVid)

for i = 1:numSteps
    writeVideo(newVid,F(i))%within the for loop saving one frame at a time
end
close(newVid)