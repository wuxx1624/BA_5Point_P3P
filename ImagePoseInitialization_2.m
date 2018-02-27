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
% iphone 6

% fx = 1229; fy = 1153;
% cx = 360; cy = 640;
% kc = [];
% fc = [fx fy];
% cc = [cx cy];
% K    = [fx 0 cx; ...
%         0 fy cy; ...
%         0 0 1];
% Kinv = [1/fx 0 -cx/fx; ...
%         0 1/fy -cy/fy; ...
%         0 0 1];

% vicon
fc = [254.997, 254.933];
cc = [326.69, 230.118];
kc = [0.925175, 0.0, 0.0, 0.0, 0.0];

fx = fc(1); fy = fc(2);
cx = cc(1); cy = cc(2);
% K    = [fx 0 cx; ...
%         0 fy cy; ...
%         0 0 1];
% Kinv = [1/fx 0 -cx/fx; ...
%         0 1/fy -cy/fy; ...
%         0 0 1];
K = eye(3);
Kinv = eye(3);

CameraParams = struct;
CameraParams.fx = fx; CameraParams.fy = fy; 
CameraParams.cx = cx; CameraParams.cy = cy; 
CameraParams.K = K; CameraParams.Kinv = Kinv;

% Control variables
% NumberOfImages = 5;
reloadImages = 0;
workingDirectory = './Set/';
reComputeMatches = 0;
pairWiseTriangulation = 0;
constructPoseGraphMatrix = 0;



matchesVisualization = 0;
reprojVisualization = 0;
multiPosesVisualization = 0;





%% Standard operations for obtaining pairwise matches and map builder
% Load images
fid = fopen('/home/jahid/documents/datasets/vicon_test/filter_result_1/image_names_1.txt','r');
bb = textscan(fid,'%s');
fclose(fid);

meaturePoses = [1;10;15];
NumberOfImages = length(meaturePoses);
I_rgb = cell(NumberOfImages, 1);
for i = 1:NumberOfImages
%     I_rgb{i} = imread([workingDirectory 'I' num2str(i) '.JPG']);
    I_rgb{i} = imread(['/home/jahid/documents',bb{1,1}{meaturePoses(i),1}]);
end

% meaturePoses = [1;2;3];
% NumberOfImages = length(meaturePoses);
% I_rgb = cell(NumberOfImages, 1);
% for i = 1:NumberOfImages
%     I_rgb{i} = imread([workingDirectory 'F' num2str(meaturePoses(i)) '.jpg']);
% end


% Save or load feature extraction + desciptor
if reloadImages == 1
    featureExtracted = cell(NumberOfImages, 1);
    featureDescriptor = cell(NumberOfImages, 1);
    for i = 1:NumberOfImages
        if size(I_rgb{i},3)~=1
            I_rgb{i} = rgb2gray(I_rgb{i});
        end
        I_single = single(I_rgb{i});
        [featureExtracted{i},featureDescriptor{i}] = vl_sift(I_single, 'PeakThresh', 2, 'edgethresh', 5);    
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
PairWiseMatches_homo = cell(nchoosek(NumberOfImages,2),1);
for i = 1:NumberOfImages-1
    for j = i+1:NumberOfImages       
        % Make sure that i < j, otherwise it will break !!!
        PairWiseMatches{SUTrilInd(NumberOfImages,i,j)} = BFNN2directions(featureExtracted{i},featureDescriptor{i},...
                                                                         featureExtracted{j},featureDescriptor{j});
        PairWiseMatches_homo{SUTrilInd(NumberOfImages,i,j)} =  PairWiseMatches{SUTrilInd(NumberOfImages,i,j)};                                                            
        homo_1 = StaticUndistort(fc,cc,kc,PairWiseMatches{SUTrilInd(NumberOfImages,i,j)}(:,1:2)');
        homo_2 = StaticUndistort(fc,cc,kc,PairWiseMatches{SUTrilInd(NumberOfImages,i,j)}(:,6:7)');
        PairWiseMatches_homo{SUTrilInd(NumberOfImages,i,j)}(:,1:2) = homo_1';
        PairWiseMatches_homo{SUTrilInd(NumberOfImages,i,j)}(:,6:7) = homo_2';
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

% change pixel value to homo 
featureExtracted_homo = cell(NumberOfImages,1);
for i = 1:NumberOfImages
    featureExtracted_homo{i}(1:2,:) = StaticUndistort(fc,cc,kc,featureExtracted{i}(1:2,:));
    featureExtracted_homo{i}(3:4,:) = featureExtracted{i}(3:4,:);
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
                [Ci_R_Cj,Ci_t_Cj,E_est,bestInlier,bestGeneratorSet] = RANSAC_5Point_pixel(PairWiseMatches_homo{SUTrilInd(NumberOfImages,i,j)},eye(3),i,j,0.01);
                InlierPairWiseMatches{SUTrilInd(NumberOfImages,i,j)} = PairWiseMatches_homo{SUTrilInd(NumberOfImages,i,j)}(bestInlier, :);                        
                Ci_P = LinearTriangulation_5pt(InlierPairWiseMatches{SUTrilInd(NumberOfImages,i,j)}, CameraParams, Ci_R_Cj, Ci_t_Cj, 0);                
%               Ci_P = LinearTriangulation(InlierPairWiseMatches{SUTrilInd(NumberOfImages,i,j)}, CameraParams.K*[eye(3) zeros(3,1)], CameraParams.K*[Ci_R_Cj', -Ci_R_Cj'*Ci_t_Cj]);
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


%% Create Pose Graph matrix
if constructPoseGraphMatrix == 1
    fprintf('Constructing PoseGraph matrix ...\n');
    [PoseGraphMatrix, pointerToImages] = ConstructPoseGraphMatrix2(featureExtracted_homo, InlierPairWiseMatches);
    save([workingDirectory 'PoseGraphMatrix.mat'], 'PoseGraphMatrix');
else
    load([workingDirectory 'PoseGraphMatrix.mat']);
end

%% Initial Feature Triangulation from Pose Graph
% 0. Absolute Poses
AbsolutePoses = zeros(3,4,NumberOfImages);

% 1. From PoseGraph, find the one having the most features used
totalFeatureCount = sum(PoseGraphMatrix~=0, 1);
[~, refPose] = max(totalFeatureCount);
remainingPoses = setdiff(1:NumberOfImages, refPose );

% 2. Find the next one which has the most overlap with this to create scale
[~, scalePoseidx] = max(sum(PoseGraphMatrix(sum(PoseGraphMatrix(:,refPose),2)~=0, ...
                                             remainingPoses)~=0, 1));
scalePose = remainingPoses(scalePoseidx);

meaturePoses = [refPose, scalePose]
remainingPoses = setdiff(1:NumberOfImages, meaturePoses);

% 3. Feature triangulated
% This data structure is the 3D feature position w.r.t to the
% reference Pose that sees it, the index is based on the nonzero index in
% the posegraph matrix, the last element is indicated if the feature is
% triangulated yet?


% 3a. Triangulate scale Pose first


[FeaturesBag, AbsolutePoses(:,:,refPose), AbsolutePoses(:,:,scalePose)] = LinearTriangulationRt(...
                                RelativePoses{SUTrilInd(NumberOfImages,refPose,scalePose)}, ...                                             
                                             PoseGraphMatrix,       ...
                                             refPose, scalePose,    ...
                                             featureExtracted_homo{refPose}, ...
                                             featureExtracted_homo{scalePose}, ...
                                             CameraParams, ...
                                             0);                                         

fprintf('Number of triangulated features: %d\n', length(find(FeaturesBag(4,:) ~= 0)));


close all;



% 3b. Triangulate for other poses using PnP 
for k = 1:NumberOfImages-2
    % Get new pose logic
    [~, newPoseidx] = max(sum(PoseGraphMatrix(FeaturesBag(4,:)~=0, ...
                                             remainingPoses)~=0, 1));
    newPose = remainingPoses(newPoseidx);
    
    % New image registration
%     [Cn_R_Cr, Cn_t_Cr, InlierIdx, ~, Cr_P, Cn_z] = PnP_RANSAC2(PoseGraphMatrix, AbsolutePoses,  ...
%                                                          featureExtracted,                      ...
%                                                          usedPoses, newPose,                      ...
%                                                          FeaturesBag, CameraParams,          ...
%                                                          300, 0.7, 30, I_rgb);
    
    [Cn_R_Cr, Cn_t_Cr, InlierIdx, ~, Cr_P, Cn_z] = P3P_RANSAC2(PoseGraphMatrix, AbsolutePoses,  ...
                                                         featureExtracted_homo,                      ...
                                                         meaturePoses, newPose,                      ...
                                                         FeaturesBag, eye(3),          ...
                                                         500, 0.7, 0.01);
%     length(InlierIdx)
    AbsolutePoses(:,:,newPose) = [Cn_R_Cr, Cn_t_Cr];
    
    fprintf('Inlier new pose I%d: %d\n', newPose, length(InlierIdx));
    
    % Visualization
    if reprojVisualization == 1
%         VisualizeAdditionalRegistration(Cn_R_Cr', -Cn_R_Cr'*Cn_t_Cr, ...
%                                         RelativePoses{SUTrilInd(NumberOfImages,1,2)}, ...
%                                         Cr_P, InlierIdx);

        reprojection_matches = ReprojectionToMatches(Cr_P, Cn_R_Cr, Cn_t_Cr, CameraParams);
        ReprojectionVisualization2(I_rgb{newPose}, Cn_z, reprojection_matches); 
    end
    
        
    % Pose refinement: TODO: fix the transpose in the last input        
    [Cn_R_Cr, Cn_t_Cr, normHist] = PnP_NL(Cn_R_Cr, Cn_t_Cr, Cr_P(:, InlierIdx), Cn_z(InlierIdx,:)',eye(3));        
%     figure;
%     semilogy(normHist, 'b*-');    
    
    % Visualization
    if reprojVisualization == 1
%         VisualizeAdditionalRegistration(Cn_R_Cr', -Cn_R_Cr'*Cn_t_Cr, ...
%                                         RelativePoses{SUTrilInd(NumberOfImages,1,2)}, ...
%                                         Cr_P, InlierIdx);        
        reprojection_matches = ReprojectionToMatches(Cr_P, Cn_R_Cr, Cn_t_Cr, CameraParams);
        ReprojectionVisualization2(I_rgb{newPose}, Cn_z, reprojection_matches); 
    end
    
    % Additional Multiposes triangulation
    FeaturesBag = MultiPoseTriangulation(FeaturesBag,           ...
                                         AbsolutePoses,         ...
                                         PoseGraphMatrix,       ...
                                         newPose, meaturePoses,    ...
                                         featureExtracted_homo,      ...                                         
                                         CameraParams,0,I_rgb);  

    fprintf('Number of triangulated features: %d\n', length(find(FeaturesBag(4,:) ~= 0)));
    
    % Visualization
    if multiPosesVisualization == 1
        VisualizeMultiPoses(AbsolutePoses, FeaturesBag, meaturePoses, newPose);
        VisualizeAllPosesReprojection(PoseGraphMatrix, AbsolutePoses, FeaturesBag, featureExtracted, meaturePoses, newPose, I_rgb, CameraParams);
    end
    
%     % Nonlinear triangulation refinement
%     FeaturesBag = NonlinearTriangulation(FeaturesBag,           ...
%                                          AbsolutePoses,         ...
%                                          PoseGraphMatrix,       ...
%                                          newPose, usedPoses,    ...
%                                          featureExtracted,      ...                                         
%                                          CameraParams,0,I_rgb);  
%     if multiPosesVisualization == 1
%         VisualizeMultiPoses(AbsolutePoses, FeaturesBag, usedPoses, newPose);
%         VisualizeAllPosesReprojection(PoseGraphMatrix, AbsolutePoses, FeaturesBag, featureExtracted, usedPoses, newPose, I_rgb, CameraParams);
%     end
    
    % Update poses been used
    meaturePoses = [meaturePoses, newPose];
    remainingPoses = setdiff(remainingPoses, newPose);    
        
    % Bundle Adjustment
    [FeaturesBag, AbsolutePoses,~,~] = BundleAdjustmentSparse(FeaturesBag, AbsolutePoses, PoseGraphMatrix, featureExtracted_homo, meaturePoses, CameraParams);
    
    if multiPosesVisualization == 1
        VisualizeMultiPoses(AbsolutePoses, FeaturesBag, meaturePoses(1:end-1), newPose);
        VisualizeAllPosesReprojection(PoseGraphMatrix, AbsolutePoses, FeaturesBag, featureExtracted, meaturePoses(1:end-1), newPose, I_rgb, CameraParams);
    end
    
%     pause;
    close all;
end

K    = [fx 0 cx; ...
        0 fy cy; ...
        0 0 1];
Kinv = [1/fx 0 -cx/fx; ...
        0 1/fy -cy/fy; ...
        0 0 1];
CameraParams.K = K; CameraParams.Kinv = Kinv;

VisualizeMultiPoses(AbsolutePoses, FeaturesBag, meaturePoses(1:end-1), newPose);
saveas(gcf,['/home/jahid/Fei/summer/BA_5point_p3p/picture/iter 12 MuliPoses.jpg']);
VisualizeAllPosesReprojection(PoseGraphMatrix, AbsolutePoses, FeaturesBag, featureExtracted, meaturePoses(1:end-1), newPose, I_rgb, CameraParams);
saveas(gcf,['/home/jahid/Fei/summer/BA_5point_p3p/picture/iter 12 AllPosesReprojection.jpg']);

