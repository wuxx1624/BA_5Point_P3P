function [R_est,t_est] = RANSAC_P3P(both_see_feat,PoseGraphMatrix,FeatureBag,featureExtracted,idx,K)

max_inliners = -1;
for iter = 1:1000
    % randomly select 6 common points, 5 use for computing essential matix,
    % 1 for folllowing test
    feat_select = both_see_feat(randperm(length(both_see_feat),6),:);
%     feat_select_1 = PoseGraphMatrix(feat_select,1);
    feat_select_2 = PoseGraphMatrix(feat_select,idx);
    
    [C,p] = P3P_Solver_3(FeatureBag(1:3,feat_select(1:3)),featureExtracted(:,feat_select_2(1:3)));
    [R_true,t_true,~]=TruePose_P3P(C,p,FeatureBag(1:3,feat_select(4)),featureExtracted(:,feat_select_2(4)),K);
    [num_inliners,error] = count_inliner_P3P(R_true,t_true,FeatureBag,featureExtracted,both_see_feat,PoseGraphMatrix(both_see_feat,idx),0.5,K);
     if num_inliners > max_inliners
        max_inliners = num_inliners;
        R_est = R_true;
        t_est = t_true;
     end
%      disp(['Iter: ',num2str(iter)]);  
end