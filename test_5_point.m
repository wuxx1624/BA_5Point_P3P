% This script using 2D matches to estimate position and orientation of new
% images

% TODO: another script to grasp the index automatically

clear all; close all;
run('/home/jahid/Fei/summer/vlfeat-0.9.20/toolbox/vl_setup.m');
addpath('./JacobianModule/');
addpath('./BundleAdjustment/');
addpath('./Drawings/');
addpath('./robotics3D/');

% Camera parameters:
% Nokia Lumia
% fx = 2210; fy = 2196;
% cx = 1292; cy = 961;

% 1222x1630
% fx = 1.25833*10^3; fy = 1.2426*10^3;
% cx = 815; cy = 611;
fx = 1229; fy = 1153;
cx = 360; cy = 640;

K    = [fx 0 cx; ...
        0 fy cy; ...
        0 0 1];
Kinv = [1/fx 0 -cx/fx; ...
        0 1/fy -cy/fy; ...
        0 0 1];

CameraParams = struct;
CameraParams.fx = fx; CameraParams.fy = fy; 
CameraParams.cx = cx; CameraParams.cy = cy; 
CameraParams.K = K; CameraParams.Kinv = Kinv;

% Control variables
usePoses = [4;3];
NumberOfImages = length(usePoses);
reloadImages = 1;
workingDirectory = './Set/';
reComputeMatches = 1;
pairWiseTriangulation = 1;
constructPoseGraphMatrix = 1;



matchesVisualization = 0;
reprojVisualization = 0;
multiPosesVisualization = 0;





%% Standard operations for obtaining pairwise matches and map builder
% Load images
I_rgb = cell(NumberOfImages, 1);
for i = 1:NumberOfImages
    I_rgb{i} = imread([workingDirectory 'A' num2str(usePoses(i)) '.jpg']);
end

% Save or load feature extraction + desciptor
if reloadImages == 1
    featureExtracted = cell(NumberOfImages, 1);
    featureDescriptor = cell(NumberOfImages, 1);
    for i = 1:NumberOfImages
        I_single = single(rgb2gray(I_rgb{i}));
        [featureExtracted{i},featureDescriptor{i}] = vl_sift(I_single, 'PeakThresh', 2);    
    end
    save([workingDirectory 'featureExtracted.mat'], 'featureExtracted');
    save([workingDirectory 'featureDescriptor.mat'], 'featureDescriptor');
else
    load([workingDirectory 'featureExtracted.mat']);
    load([workingDirectory 'featureDescriptor.mat']);
end

% Extract matches pairwise
if reComputeMatches == 1
PairWiseMatches = cell(nchoosek(NumberOfImages,2),1);
for i = 1:NumberOfImages-1
    for j = i+1:NumberOfImages       
        % Make sure that i < j, otherwise it will break !!!
        PairWiseMatches{SUTrilInd(NumberOfImages,i,j)} = BFNN2directions(featureExtracted{i},featureDescriptor{i},...
                                                                         featureExtracted{j},featureDescriptor{j});
        if matchesVisualization == 1
%             VisualizeMatches(I_rgb{i}, I_rgb{j}, PairWiseMatches{SUTrilInd(NumberOfImages,i,j)}, 'b' ,'g');
%             pause;
        end
    end
end
    save([workingDirectory 'PairWiseMatches.mat'], 'PairWiseMatches');
else
    load([workingDirectory 'PairWiseMatches.mat']);
end

%% Outlier rejection + pose estimation
% RANSAC fundamental matrix computation
InlierPairWiseMatches = cell(nchoosek(NumberOfImages,2),1);
FeatureTriangulated = cell(nchoosek(NumberOfImages,2),1);
RelativePoses = cell(nchoosek(NumberOfImages,2),1);
if pairWiseTriangulation == 1
    fprintf('Extracting inlier matches using Fundamental matrix ...\n');
    for i = 1:NumberOfImages-1
        for j = i+1:NumberOfImages
            % Make sure that i < j, otherwise it will break !!!      
            if size(PairWiseMatches{SUTrilInd(NumberOfImages,i,j)},1) > 24 % Make sure sufficient matches to be inlier
%                 [Fbest, bestInlier, bestGeneratorSet, iter] = FMatrixRANSAC( PairWiseMatches{SUTrilInd(NumberOfImages,i,j)}, 300, 0.7, 1);                
                [Ci_R_Cj,Ci_t_Cj,E_est,bestInlier,bestGeneratorSet] = RANSAC_5Point(PairWiseMatches{SUTrilInd(NumberOfImages,i,j)},K,i,j);
                InlierPairWiseMatches{SUTrilInd(NumberOfImages,i,j)} = PairWiseMatches{SUTrilInd(NumberOfImages,i,j)}(bestInlier, :);
%                             [Ci_P, Ci_R_Cj, Ci_t_Cj] = LinearTriangulationF(InlierPairWiseMatches{SUTrilInd(NumberOfImages,i,j)}, ...
%                 CameraParams, Fbest, 0);   
%                 [Ci_P, Ci_R_Cj, Ci_t_Cj] = LinearTriangulationE(InlierPairWiseMatches{SUTrilInd(NumberOfImages,i,j)},CameraParams, R_est,t_est, 0);  
                Ci_P = LinearTriangulation(InlierPairWiseMatches{SUTrilInd(NumberOfImages,i,j)}, CameraParams.K*[eye(3) zeros(3,1)], CameraParams.K*[Ci_R_Cj', -Ci_R_Cj'*Ci_t_Cj]);
                FeatureTriangulated{SUTrilInd(NumberOfImages,i,j)} = Ci_P;
                RelativePoses{SUTrilInd(NumberOfImages,i,j)} = [Ci_R_Cj Ci_t_Cj];
            else
                InlierPairWiseMatches{SUTrilInd(NumberOfImages,i,j)} = zeros(1,10);  % invalid indices
            end


            if matchesVisualization == 1
                VisualizeMatches(I_rgb{i}, I_rgb{j}, PairWiseMatches{SUTrilInd(NumberOfImages,i,j)}(bestGeneratorSet, :), 'm' ,'m');
                VisualizeMatches(I_rgb{i}, I_rgb{j}, PairWiseMatches{SUTrilInd(NumberOfImages,i,j)}(bestInlier, :), 'b' ,'g');
                % VisualizeEpipolarLines(I_rgb{i}, I_rgb{j}, PairWiseMatches{SUTrilInd(NumberOfImages,i,j)}(bestInlier, :), Fbest);
%                 pause;
            end

            if reprojVisualization == 1
                rmatches1 = ReprojectionToMatches(Ci_P, eye(3), zeros(3,1),CameraParams);
                rmatches2 = ReprojectionToMatches(Ci_P, Ci_R_Cj', -Ci_R_Cj'*Ci_t_Cj,CameraParams);
                ReprojectionVisualization(I_rgb{i}, I_rgb{j}, InlierPairWiseMatches{SUTrilInd(NumberOfImages,i,j)}, [rmatches1 rmatches2] ); 
%                 pause;
            end

        end
    end

    save([workingDirectory 'InlierPairWiseMatches.mat'], 'InlierPairWiseMatches');
    save([workingDirectory 'RelativePoses.mat'], 'RelativePoses');
    save([workingDirectory 'FeatureTriangulated.mat'], 'FeatureTriangulated');
else
    load([workingDirectory 'InlierPairWiseMatches.mat']);
    load([workingDirectory 'RelativePoses.mat']);
    load([workingDirectory 'FeatureTriangulated.mat']);
end