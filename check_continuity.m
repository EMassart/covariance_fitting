% This code checks whether the patchwise sectional approach is continuous.

clear all; close all; clc;

% generates the data points
for i = 1:4
    Y(:,:,i) = randn(5,2);
end

% computes the geodesic between A_2 and A_4, when computations are done in
% the section based at A_1
Y_projected = Y;
for i = 2:4
    Q = orth_pol(Y(:,:,1)'*Y(:,:,i));
    Y_projected(:,:,i) = Y(:,:,i)*Q';
end
weights = 0.5;
YS_patch1 = (1-weights)*Y_projected(:,:,2) + weights*Y_projected(:,:,4);


% second approach: interpolation in the next patch (section based at Y2).
% We change the root of the section.
Y_projected = Y;
for i = 2:4
    Q = orth_pol(Y(:,:,2)'*Y(:,:,i));
    Y_projected(:,:,i) = Y(:,:,i)*Q';
end
weights = 0.5;  % interpolation between Y2 and Y4
YS_patch2 = (1-weights)*Y_projected(:,:,2) + weights*Y_projected(:,:,4);


err = norm(YS_patch1*YS_patch1' - YS_patch2*YS_patch2','fro');
fprintf('Discontinuity gap is equal to %4.2e \n',err);