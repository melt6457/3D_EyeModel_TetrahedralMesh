%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File:     validation.m
% Creators: Kory Melton and Ian Besse
% Date:     
% Purpose:  To validate the model using corneal apex displacement
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;

dt = 0.1;
time = 5;

CreateTetrahedralMesh;
springConstant;

disp(1,:) = run(5, time, dt);
disp(2,:) = run(10, time, dt);
disp(3,:) = run(15, time, dt);
disp(4,:) = run(20, time, dt);
disp(5,:) = run(25, time, dt);

times(:) = 0:dt:time;

plot(times, disp(1,:), times, disp(2,:), times, disp(3,:), ...
     times, disp(4,:), times, disp(5,:))

