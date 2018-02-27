function [FeatureBag, AbsolutePoses, normHist, gradnormHist] = BundleAdjustment(FeatureBag, AbsolutePoses, PoseGraphMatrix, featureExtracted, usedPoses, CameraParams)

PoseGraphMatrix_r = PoseGraphMatrix(:, usedPoses);
NumOfPoses = length(usedPoses);
NumOfFeatures = length(FeatureBag(FeatureBag(4,:)~=0));

% figure(3); clf; hold on;
% Draw3DWorld(AbsolutePoses, FeatureBag, [eye(3) zeros(3,1)], usedPoses);
% axis([-2 2 -2 2 -1 5]);
%% Bundle Adjustment
regulator = 1;
norm_iter = 100;
grad_norm = 100;
normHist = [];
gradnormHist = [];
maxIter = 100;
iter = 1;

while iter < maxIter && norm_iter > 1e-1 && grad_norm > 10e-4
    
% 1. Build Jacobian Jp Jf corresponds to each feature
Jpose = zeros(2*nnz(PoseGraphMatrix_r), 6*NumOfPoses);
Jfeat = zeros(2*nnz(PoseGraphMatrix_r), 3*NumOfFeatures);
residual = zeros(2*nnz(PoseGraphMatrix_r), 1);
row_count = 0;
% last_row_count = 0;

featureAvailableList = find(FeatureBag(4,:));
featureAvailableList = featureAvailableList-1;

for featk = featureAvailableList %0:size(PoseGraphMatrix,1)-1
      
   poses_see_featk = find(PoseGraphMatrix_r(featk+1, :));
   
   
   for l = 0:length(poses_see_featk)-1
       
      %  Pose Jacobian
      Cr_p_f = FeatureBag(1:3 ,featk+1);      
      Cl_z_f = CameraParams.Kinv(1:2,:)*[featureExtracted{usedPoses(poses_see_featk(l+1))}(1:2,PoseGraphMatrix_r(featk+1,poses_see_featk(l+1)));1];                  
      
      
      % Pose Transformation
      Cl_T_W = AbsolutePoses(:,:,usedPoses(poses_see_featk(l+1)));
      Cr_T_W = AbsolutePoses(:,:,FeatureBag(4, featk+1));      
      Cl_T_Cr = Cl_T_W*[InversePose(Cr_T_W);zeros(1,3) 1];            
      Cl_p_f = Cl_T_Cr*[Cr_p_f;1];                   
       
       
      [Jpose(row_count+1:row_count+2, 6*(poses_see_featk(l+1)-1)+1:6*(poses_see_featk(l+1)-1)+6),...
          residual(row_count+1:row_count+2)] = ...
           dprojection_dp_dtt2(Cr_p_f, Cl_T_Cr(:,1:3), Cl_T_Cr(:,4), Cl_z_f);      
      
      Jfeat(row_count+1:row_count+2, 3*featk+1:3*featk+3) = ...
        1/Cl_p_f(3)*[1 0 -Cl_p_f(1)/Cl_p_f(3);0 1 -Cl_p_f(2)/Cl_p_f(3)] * Cl_T_Cr(:,1:3);
    
      row_count = row_count + 2;
   end   
   
   % residual(last_row_count+1:row_count) = Jpose_Jfeat_res{featk+1,3};
   
   % last_row_count = row_count;
   
   % figure(6);
   % spy([Jpose Jfeat]);
   % pause;
end

norm_iter = norm(residual);
normHist = [normHist norm_iter];
if iter > 20
    norm_iter = 1e2*abs(normHist(end-1) - normHist(end));
end

A11 = Jpose'*Jpose;
A12 = Jpose'*Jfeat;
A22 = Jfeat'*Jfeat;
res_ = [Jpose';Jfeat']*residual;

% A = [A11 A12;A12' A22];
% A = A + regulator*eye(size(A));
% x_tilde = A \ res_;

A22_inv = zeros(size(A22));
for featk = featureAvailableList
    A22_inv(3*featk+1:3*featk+3, 3*featk+1:3*featk+3) = inv(A22(3*featk+1:3*featk+3,3*featk+1:3*featk+3) + regulator*eye(3));
end

% solve:
p_tilde = (A11 +regulator*eye(size(A11))- A12*A22_inv*A12') \ (res_(1:6*NumOfPoses) - A12*A22_inv*res_(1+6*NumOfPoses:end));

% p_tilde = x_tilde(1:6*NumOfPoses);

for l = 1:NumOfPoses-1
    dq = [0.5*p_tilde(6*l+4:6*l+6);1];
    dq = dq / norm(dq);
    AbsolutePoses(:,1:3,usedPoses(l+1)) = quat2rot( quat_mul(dq, rot2quat(AbsolutePoses(:,1:3,usedPoses(l+1)))) );       
    AbsolutePoses(:,4,usedPoses(l+1)) = AbsolutePoses(:,4,usedPoses(l+1)) + p_tilde(6*l+1:6*l+3);
end

grad_norm = norm(p_tilde)^2;
 
% 4. Efficient solve for each feature f
bf_p = res_(1+6*NumOfPoses:end) - A12'*p_tilde;
for featk = featureAvailableList
   D22 = A22_inv(3*(featk)+1:3*(featk)+3, 3*(featk)+1:3*(featk)+3);
   FeatureBag(1:3,featk+1) = FeatureBag(1:3,featk+1) +  D22 * bf_p(3*featk+1:3*featk+3);
   grad_norm = grad_norm + norm(D22 * bf_p(3*featk+1:3*featk+3))^2;
   % FeatureBag(1:3,featk+1) = FeatureBag(1:3,featk+1) + x_tilde(3*featk+6*NumOfPoses+1:3*featk+6*NumOfPoses+3);
end

gradnormHist = [gradnormHist sqrt(grad_norm)];

% 6. Visualization
% figure(3); clf; hold on;
% Draw3DWorld(AbsolutePoses, FeatureBag, [eye(3) zeros(3,1)], usedPoses);
% grid on;



    fprintf('Iter :%d, norm: %f, gradnorm: %f \n', iter, norm_iter, grad_norm);
    iter = iter + 1;
end
