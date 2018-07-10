function P = uniformPoints

m = 1; % steps through P
b = 12; % used to create circle for eye
step = 1; % used to control how many points there are
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Number of points for each step value
% Step 1: 7,153 points
% Step 2: 925 points
% Step 3: 257 points
% Step 4: 123 points
% Step 5: 56 points (not as useful)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = -b:step:b
    for j = -b:step:b
        for k = -b:step:b
            P (m,:) = [i,j,k];
            m = m + 1;
        end
    end
end

% Grab only the points less than b away from origin
Ind = find (sum((P.*P)') <= b^2);
P = P (Ind,:);