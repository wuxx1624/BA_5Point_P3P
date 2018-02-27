function     [Rbest, tbest, bestInlier, bestGeneratorSet, featureUsed, measurementUsed] = P3P_RANSAC2(PoseGraphMatrix, ...
                                                                           AbsolutePoses,    ...
                                                                           featureExtracted, ...
                                                                           usedPoses, newPose,         ...
                                                                           FeaturesBag, K, ...
                                                                           MaxIter, inlier_ratio, threshold)  
% 2D-3D correspondences
% extract the indices of features seen in refPose in posegraph matrix that
% is already triangulated

P = []; 
w = [];
for posei = usedPoses
    % available_features = FeaturesBag(4,FeaturesBag(4,:) == posei);
    % [~, available_features_idx] = ismember(available_features, PoseGraphMatrix(:,refPose));
    available_features_idx = find(FeaturesBag(4,:) == posei);
    [w_full, common_feature_indices] = GetMatchFromPoseGraph(PoseGraphMatrix, available_features_idx, posei, newPose, ...
                               featureExtracted{posei}, featureExtracted{newPose});
%     P = [P AbsolutePoses(:,1:3,posei)'*(FeaturesBag(1:3,common_feature_indices) - ...
%                     repmat(AbsolutePoses(:,4,posei), [1 length(common_feature_indices)]))];
    P = [P FeaturesBag(1:3,common_feature_indices)];
    
    w = [w;w_full(:,6:7)];
end

                           
% % refPose_nz_elements = PoseGraphMatrix(PoseGraphMatrix(:,refPose)~=0,refPose);
% [~, available_features_idx] = ismember(w_full(:,5)', FeaturesBag(4,:));
% w = w_full(:,6:7);
% P = FeaturesBag(1:3, available_features_idx);
featureUsed = P;
measurementUsed = w;

% VisualizeMatches(I_rgb{refPose}, I_rgb{newPose}, w_full, 'b' ,'g');
% reprojection_matches = ReprojectionToMatches(P, eye(3), zeros(3,1), CameraParams);
% ReprojectionVisualization2(I_rgb{1}, w_full(:,1:2), reprojection_matches); 


% RANSAC parameter
Rbest = []; tbest = [];
MaxInlier = inlier_ratio*size(w,1);
bestInlier = [];
bestGeneratorSet = [];
iter = 1;

while iter < MaxIter
    
    % Check uniqueness
    % Picking 6 distinct points
    picked_points_idx = randi(size(w,1),1,4);
    picked_points_idx = unique(picked_points_idx);
    while length(picked_points_idx) < 4
        picked_points_idx = [picked_points_idx randi(size(w,1),1,4 - length(picked_points_idx))];
        picked_points_idx = unique(picked_points_idx);
    end
    
    % Compute pose
    [C,p] = P3P_Solver(P(1:3,picked_points_idx(1:3)), w(picked_points_idx(1:3), :)');
    if isempty(C) || isempty(p)
        continue;
    end
    [R,t,~]=TruePose_P3P(C,p,P(1:3,picked_points_idx(4)),w(picked_points_idx(4), :)',K);
    
    % Compute inlier based on reprojection
    inlierIdx = ComputeInlier3D2DReprojection(P, w, K, R, t, threshold);
    
    if length(inlierIdx) > MaxInlier
        Rbest = R;
        tbest = t;
        bestInlier = inlierIdx;
        bestGeneratorSet = picked_points_idx;
        break;
    else
        if length(inlierIdx) > length(bestInlier)
            Rbest = R;
            tbest = t;
            bestInlier = inlierIdx;
            bestGeneratorSet = picked_points_idx;
        end
    end 
    
    iter = iter + 1;
end