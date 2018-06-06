%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File:     validation.m
% Creators: Kory Melton and Ian Besse
% Date:     
% Purpose:  To validate the model using corneal apex displacement
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;

dt = 0.1;
time = 30;

disp(1,:) = validate(5, time, dt);
disp(2,:) = validate(10, time, dt);
disp(3,:) = validate(15, time, dt);
disp(4,:) = validate(20, time, dt);
disp(5,:) = validate(25, time, dt);

times(:) = 0:dt:time;

plot(times, disp(1,:), times, disp(2,:), times, disp(3,:), ...
     times, disp(4,:), times, disp(5,:))

function disp = validate(Force, time, dt)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loading Data Into Workspace...
% 
% File: meshInit.mat
% Data: fixedPoints - points on the outside of the eye that don't move 
%       frontPoints - points in the front of the eye that receive impact
%       DT          - the delaunayTriangulation
%
% File: EdgeInit.mat
% Data: K      - 
%       N      - 
%       Points - 
%       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('Data/MeshInit.mat'); % Contains DT values
load('Data/EdgeInit.mat'); % Contains K, N, Points

tic % start timer

NZK = sum(sum(K~=0));
TK = sum(sum(K));
Kavg = TK/NZK;

numSteps = time / dt; % number of steps
print = numSteps / 10; % print out mod

Edges = edges(DT); % get edges from DT
[numEdges, ~] = size (Edges); % get number of edges
[numPoints, ~] = size (Points); % get number points
numFixed = length (fixedPoints);
numFront = length (frontPoints);
eyeMass = .0075;
m = eyeMass / numPoints; % set mass
c = 0.2 * sqrt(Kavg * m); % Estimate for damping coefficient

ff = [ones(numFront, 1), zeros(numFront, 2)]; % front force

Fmag = zeros(numPoints, numSteps); % used to store forces

% initialize P
P = zeros(numPoints, numSteps, 3);
P(:,1,:) = Points;
% initialize D
D = zeros(numPoints, numPoints, 3);

PS(:, :) = P(:, 1, :);
VS = 0*PS;

% set initial force
F_0 = Force / 1000;
impactTime = 0.2;
V_0 = (F_0*impactTime)/(numFront*m);

VS(frontPoints, :) = V_0*ff;
%VS (frontPoints, :) = 5*ff;

% initialize runge-kutta matrices
rungeKutta_Pos = zeros (numPoints, 3, 4);
rungeKutta_Vel = zeros (numPoints, 3, 4);
rungeKutta_Acc = zeros (numPoints, 3, 4);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Iterate for the length of the simulation (1 to final time index)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:(numSteps) 
    % get the current Position-State
    PS(:, :) = P(:, i, :);
    
    % get the previous Positions and velocities
    rungeKutta_Pos(:,:,1) = PS;
    rungeKutta_Vel(:,:,1) = VS;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Iterate for runge-kutta steps
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for j = 1:4      
        % get the vectors of the current edges
        for n = 1:numEdges
            nodes(:) = Edges(n,:);
            D(nodes(1), nodes(2), :) = rungeKutta_Pos(nodes(2), :, j) - ...
                                       rungeKutta_Pos(nodes(1), :, j);
            D(nodes(2), nodes(1), :) = rungeKutta_Pos(nodes(1), :, j) - ...
                                       rungeKutta_Pos(nodes(2), :, j);
        end
        
        % get the norm of the D vector
        D_n(:,:) = sqrt(D(:,:,1).^2 + D(:,:,2).^2 + D(:,:,3).^2);
        
        F_m = K.*(N./D_n - 1); % magnitude of forces
        F_m(isnan(F_m)) = 0; % make NaN forces 0
        F_m(isinf(F_m)) = 0; % make infinite forces 0
        
        F1 = -F_m.*D(:,:,1); % get the force vector in the x
        F2 = -F_m.*D(:,:,2); % get the force vector in the y
        F3 = -F_m.*D(:,:,3); % get the force vector in the z
        
        % should this be absolute value? -- I don't think so
        Fx = (sum(F1'))'; % get the sum of the forces in x
        Fy = (sum(F2'))'; % get the sum of the forces in y
        Fz = (sum(F3'))'; % get the sum of the forces in z
        
        F = [Fx, Fy, Fz]; % put the forces in one matrix
        % set velocity of fixed points
        F(fixedPoints, :) = zeros(numFixed, 3);
        
        rungeKutta_Acc (:, :, j) = (1/m)*(F - c*rungeKutta_Vel (:, :, j));
        
        % get the change (d) in Position/Velocity (all)
        if j < 3
            rungeKutta_Vel (:, :, (j + 1)) = rungeKutta_Vel (:, :, 1) + ...
                                         (dt/2)*rungeKutta_Acc (:, :, j);
            rungeKutta_Pos (:, :, (j + 1)) = rungeKutta_Pos (:, :, 1) + ...
                                         (dt/2)*rungeKutta_Vel (:, :, j); 
        end
        
        if j == 3
          rungeKutta_Vel (:, :, (j + 1)) = rungeKutta_Vel (:, :, 1) + ...
                                         (dt)*rungeKutta_Acc (:, :, j);
          rungeKutta_Pos (:, :, (j + 1)) = rungeKutta_Pos (:, :, 1) + ...
                                         (dt)*rungeKutta_Vel (:, :, j);   
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Complete Iteration for runge-kutta steps
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Calculate New Position-State
    PS = PS + (dt/6)*(rungeKutta_Vel (:, :, 1) + ...
         2*rungeKutta_Vel (:, :, 2) + 2*rungeKutta_Vel (:, :, 3) + ...
         rungeKutta_Vel (:, :, 4));
    Vold = VS;
    % Calculate New Velocity-State
    VS = VS + (dt/6)*(rungeKutta_Acc (:, :, 1) + ...
         2*rungeKutta_Acc (:, :, 2) + 2*rungeKutta_Acc (:, :, 3) + ...
         rungeKutta_Acc (:, :, 4));
    Fv = ((VS - Vold) / dt)*m;
    Fmag(:,i) = sqrt(sum((Fv.^2)')');
    
    % Update New Positions
    P(:, i+1, :) = PS;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Complete iteration for length of simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get vectors to determine the distance from the initial position to the
% new position at each time step for point 1
distX (:) = P(1,1,1) - P(1,:,1);
distY (:) = P(1,1,2) - P(1,:,2);
distZ (:) = P(1,1,3) - P(1,:,3);

disp (:) = sqrt(distX(:).^2 + distY(:).^2 + distZ(:).^2);

toc % end timer

end