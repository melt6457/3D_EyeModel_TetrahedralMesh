%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File:     run.m
% Creators: Kory Melton and Ian Besse
% Date:     
% Purpose:  To move the system of nodes using a 4th-order runge-kutta
%           method.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function run
tic % start timer

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
clear % Clear previous data
load('Data/MeshInit.mat'); % Contains DT values
load('Data/EdgeInit.mat'); % Contains K, N, Points

dt = .1; % set time-step to .1 millisecond
time = 20; % set time to milliseconds

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

Fmag = zeros(numPoints, numSteps); % used to store forces

% initialize P
P = zeros(numPoints, numSteps, 3);
P(:,1,:) = Points;
% initialize D
D = zeros(numPoints, numPoints, 3);

PS(:, :) = P(:, 1, :);
VS = 0*PS;

% set initial force
theta = 0;
[InitVeloPoints, v_direction] = ImpactInitializer(theta);
Force = 7.5; % in Newtons
F_0 = Force / 1000;
impactTime = 0.2;
[numInitPoints, ~] = size (InitVeloPoints);

V_0 = (F_0*impactTime)/(numInitPoints*m);
forces = [ones(numInitPoints, 3)]; % init forces
forces = forces.*v_direction;

VS(InitVeloPoints, :) = V_0*forces;

% initialize runge-kutta matrices
rungeKutta_Pos = zeros (numPoints, 3, 4);
rungeKutta_Vel = zeros (numPoints, 3, 4);
rungeKutta_Acc = zeros (numPoints, 3, 4);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Iterate for the length of the simulation (1 to (final time index - 1))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:(numSteps - 1) 
    % get the current Position-State
    PS(:, :) = P(:, i, :);
    
    % get the previous Positions and velocities
    rungeKutta_Pos(:,:,1) = PS;
    rungeKutta_Vel(:,:,1) = VS;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Iterate for runge-kutta steps
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for j = 1:4      
        % get the vectors of the current nodes
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
    
    if 0 == mod(i,print)
        %print(i)
        toc
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Complete iteration for length of simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Save Positions to a File
save('Data/Positions.mat', 'P', 'numSteps', 'Fmag');
toc % end timer

end
