% function [R_est,t_est] = RANSAC_5Point(both_see_feat,PoseGraphMatrix,featureExtracted,K,idx1,idx2)
function [R_est,t_est,E_est,bestInlier,bestGeneratorSet] = RANSAC_5Point(both_see_feat,K,idx1,idx2,thresh)
max_inliners = -1;
% thresh = 0.01;
for iter = 1:500
    % randomly select 6 common points, 5 use for computing essential matix,
    % 1 for folllowing test
    feat_select = randperm(size(both_see_feat,1),6);
    feat_select_1 = both_see_feat(feat_select,1:2);
    feat_select_2 = both_see_feat(feat_select,6:7);
    feat_select_1 = feat_select_1';
    feat_select_2 = feat_select_2';
    % compute essential matrix candidates (upsto 10)
    Evec = calibrated_fivepoint([feat_select_1(:,1:5);ones(1,5)],[feat_select_2(:,1:5);ones(1,5)]);
%     Evec = calibrated_fivepoint_GrLex(featureExtracted{1}(:,feat_select_1(1:5)),featureExtracted{2}(:,feat_select_2(1:5)));
    % check which one among all candidates is correct one
    % by q2'*E*q1 = 0 
    Q1 = [feat_select_1(:,6);1];
    Q2 = [feat_select_2(:,6);1];
    Q1 = Q1';
    Q2 = Q2';
    Q = [Q2(:,1).*Q1(:,1) , ...
         Q2(:,2).*Q1(:,1) , ...
         Q2(:,3).*Q1(:,1) , ... 
         Q2(:,1).*Q1(:,2) , ...
         Q2(:,2).*Q1(:,2) , ...
         Q2(:,3).*Q1(:,2) , ...
         Q2(:,1).*Q1(:,3) , ...
         Q2(:,2).*Q1(:,3) , ...
         Q2(:,3).*Q1(:,3) ] ; 
     Evec_error = zeros(size(Evec,2),1);
     for i=1:size(Evec,2)
         Evec_error(i) = abs(Q*Evec(:,i));
     end
     [~,Evec_idx] = min(Evec_error);
     
      % reshape essential matrix to be 3*3
      if ~isempty(Evec)
          E = reshape(Evec(:,Evec_idx),[3,3]);
      else
          E = eye(3);
      end
      % compute 4 candidates for R and t
      [R(:,:,1), t(:,:,1), R(:,:,2), t(:,:,2), R(:,:,3), t(:,:,3), R(:,:,4), t(:,:,4)] = CameraPose(E);
      [R_true,t_true,n]=TruePose(R,t,[feat_select_1(:,6);1],[feat_select_2(:,6);1],both_see_feat,thresh,K);
      [num_inliners,error,inlierIdx] = count_inliner(R_true,t_true,both_see_feat,thresh,K,idx1,idx2);
      if num_inliners > max_inliners
          max_inliners = num_inliners;
          R_est = R_true';
          t_est = -R_true'*t_true;
          E_est = E;
          bestInlier = inlierIdx;
          bestGeneratorSet = feat_select;
      end
%       disp(['Iter: ',num2str(iter)]);      
end