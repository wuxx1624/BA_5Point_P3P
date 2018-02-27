%% Ultimate test for selecting features to estimate fundamental matrix
clear all; close all;

% Camera matrix
fx = 2210; fy = 2196;
cx = 1292; cy = 961;
K    = [fx 0 cx; ...
        0 fy cy; ...
        0 0 1];
Kinv = [1/fx 0 -cx/fx; ...
        0 1/fy -cy/fy; ...
        0 0 1];

% Original images
I1o = imread('I10.jpg');
I2o = imread('I11.jpg');
dir = './I10I11/';
load([dir 'matches.mat']);

% Draw matches 
figure(1) ; clf ;
imagesc(cat(2, I1o, I2o)) ;
for k = 1:size(matches,1)
    xa = matches(k,1) ;
    xb = matches(k,5) + size(I1o,2) ;
    ya = matches(k,2) ;
    yb = matches(k,6) ;

    hold on ;
    h = line([xa ; xb], [ya ; yb]) ;
    set(h,'linewidth', 1, 'color', 'b') ;

    vl_plotframe(matches(k, 1:4)') ;
    corr = matches(k, 5:8) + [size(I1o,2) 0 0 0];
    vl_plotframe(corr') ;    
end
axis image off ;    

% Click and choose matches that separating 
repick = 1;
if repick == 1
    figure(1); hold on;
    P = pickpoints(8, 'y');
    uvmatrix = [];
    idx_selected_point = [];
    for i = 1:8
        res_vec = matches(:,1:2) - repmat(P(i,:),size(matches,1),1);
        res_norm = res_vec .* res_vec;
        resnorm = res_norm(:,1) + res_norm(:,2);
        [~, j] = min(resnorm);

        uvmatrix = [uvmatrix; [matches(j,1:2) matches(j,5:6)]];
        idx_selected_point = [idx_selected_point j];
    end
    save([dir 'P.mat'], 'P');
else
    load([dir 'P.mat']);
    uvmatrix = [];
    idx_selected_point = [];
    for i = 1:8
        res_vec = matches(:,1:2) - repmat(P(i,:),size(matches,1),1);
        res_norm = res_vec .* res_vec;
        resnorm = res_norm(:,1) + res_norm(:,2);
        [~, j] = min(resnorm);

        uvmatrix = [uvmatrix; [matches(j,1:2) matches(j,5:6)]];
        idx_selected_point = [idx_selected_point j];
    end
end
% Compute fundamental matrix
F = ComputeFundamentalMatrix(uvmatrix);

% Compute inlier
threshold = 60;
[inlierIdx,R,t] = ComputeInlierReprojUltimateTest(I1o, I2o, idx_selected_point, F, K, matches, threshold);
%Fls = ComputeFundamentalMatrixInlier([matches(inlierIdx,1:2) matches(inlierIdx,3:4)]);
CameraPoseCalculation(F,K,1);

