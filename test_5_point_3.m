clear all; close all;
addpath('./Drawings/');
addpath('../robotics3D');
% addpath('../SelectIn2Regions');
% addpath('../SelectIn2Regions/Regions Partition');

base_path = '/home/jahid/documents/datasets/';
dataset_path = fullfile(base_path, '/vicon_test/');
working_directory = './ComputedPose_3';

Max_iter = 1;

PnP_iter = 5;

% control variables
saveData = 0;






    
%% load the 3D features idx and position
FeatureBag_total = dlmread(fullfile(dataset_path, 'ba_result_3/landmark_id_position.txt'));
NumOfFeatures = max(FeatureBag_total(:,1)) + 1;
FeatureBag = zeros(3,NumOfFeatures);
for k = 1:size(FeatureBag_total,1)
    FeatureBag(:,FeatureBag_total(k,1)+1) = FeatureBag_total(k,2:4)';
end

%% load the 2D features idx and position of each clone
featureExtracted_total = dlmread(fullfile(dataset_path, 'ba_result_3/landmark_homogeneous_measurements.txt'));


NumOfPoses = length(unique(featureExtracted_total(:,2)));

if saveData == 1
featureExtracted = cell(NumOfPoses,1);
featureExtracted_z = cell(NumOfPoses,1);
PoseGraphMatrix = zeros(NumOfFeatures, NumOfPoses);
for l = 1:size(featureExtracted_total,1) 
    clone_idx = featureExtracted_total(l,2);
    feature_idx = featureExtracted_total(l,1);
    PoseGraphMatrix(feature_idx+1,clone_idx+1) = 1;
    % homo value
    featureExtracted{clone_idx+1} = [featureExtracted{clone_idx+1} [featureExtracted_total(l,3:4)';1]];
    % z value
    featureExtracted_z{clone_idx+1} = [featureExtracted_z{clone_idx+1} [featureExtracted_total(l,5:6)'+[1;1];1]];
end


% %% load xkk
% xkk_ba = dlmread(fullfile(dataset_path, 'ba_result_3/xkk.txt'));
% AbsolutePoses_true = zeros(3,4,size(xkk_ba,2),1);
% for i = 1:size(xkk_ba,2)
%     L_R_G = rotz(pi)*quat2rot(xkk_ba(1:4,i));
%     G_t = xkk_ba(14:16,i);
%     L_t = -L_R_G*G_t;
%     AbsolutePoses_true(:,:,i) = [L_R_G L_t];
% end

%% find common features
both_seen_feat_3D = cell(NumOfPoses,NumOfPoses);
PairWiseMatches = cell(NumOfPoses,NumOfPoses);
for i = 1:NumOfPoses
    for j = 1:NumOfPoses
        if i == j
            both_seen_feat_3D{i,j} = 0;
        else
                both_seen_feat_3D{i,j} = find(PoseGraphMatrix(:,i)+PoseGraphMatrix(:,j)==2);
        end
    end
end

%% set PoseGraphMatrix and mark the features seen by each clone
for k = 1:NumOfPoses    
    numFeatPosek = find(PoseGraphMatrix(:,k));
    for l = 1:length(numFeatPosek)
        % assign correct feature index for pose graph matrix
        PoseGraphMatrix(numFeatPosek(l),k) = l; 
    end
end

%% collect the features seen by each clone
fc = 254.997;

seen_feat = cell(NumOfPoses,1);
for k = 1:NumOfPoses
    seen_feat{k} = find(PoseGraphMatrix(:,k));
end

for i = 1:NumOfPoses
    for j = 1:NumOfPoses
        if i == j
            PairWiseMatches{i,j} = 0;
        else
                PairWiseMatches{i,j} = zeros(length(both_seen_feat_3D{i,j}),7);
                PairWiseMatches{i,j}(:,1:2) = featureExtracted{i}(1:2,PoseGraphMatrix(both_seen_feat_3D{i,j},i))';
                PairWiseMatches{i,j}(:,6:7) = featureExtracted{j}(1:2,PoseGraphMatrix(both_seen_feat_3D{i,j},j))';
        end
    end
end


    save ('/home/jahid/Fei/summer/BA_5point_p3p/database/PairWiseMatches_test.mat','PairWiseMatches');
else
    load ('/home/jahid/Fei/summer/BA_5point_p3p/database/PairWiseMatches_test.mat','PairWiseMatches');
end
%% camera figures

fc = [254.997, 254.933];
cc = [326.69, 230.118];
kc = [0.925175, 0.0, 0.0, 0.0, 0.0];

%% initial values
% initial value of camera pose
load(['/home/jahid/Fei/summer/Feature_selection_simulation/PNP_testing_FS/ComputedPose_3/AbsolutePoses_iter 5.mat'], 'AbsolutePoses');
AllPoses = AbsolutePoses;
AbsolutePoses = zeros(3,4,NumOfPoses);
AbsolutePoses_true = zeros(3,4,NumOfPoses);

% referencePose = randi(NumOfPoses,1);
% measturePose = randi(NumOfPoses,[1,1]);
% i = referencePose;
i = 1;
for j = 10:20
    % R and t is j w.r.t i
    [Ci_R_Cj,Ci_t_Cj,E_est,bestInlier,bestGeneratorSet] = RANSAC_5Point(PairWiseMatches{i,j},eye(3),i,j);
    AbsolutePoses(:,:,j) = [Ci_R_Cj,Ci_t_Cj];
    G_AllPoses_j = [AllPoses(:,1:3,j)' -AllPoses(:,1:3,j)'*AllPoses(:,4,j)];
    AllPoses_relative = [AllPoses(:,:,i);[0 0 0 1]]*[G_AllPoses_j;[0 0 0 1]];
    AbsolutePoses_true(:,:,j) = AllPoses_relative(1:3,:);
    AbsolutePoses_true(:,4,j) = AbsolutePoses_true(:,4,j)/norm(AbsolutePoses_true(:,4,j));
end


