function [FeatureBag, AbsolutePoses, normHist, gradnormHist] = BundleAdjustmentSparse(FeatureBag, AbsolutePoses, PoseGraphMatrix, featureExtracted, usedPoses, CameraParams)

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
maxIter = 150;
iter = 1;

while iter < maxIter && norm_iter > 1e-1 && grad_norm > 1e-5
    
% 1. Build Jacobian Jp Jf corresponds to each feature
Jpose = zeros(2*nnz(PoseGraphMatrix_r), 6*NumOfPoses);
Jfeat = zeros(2*nnz(PoseGraphMatrix_r), 3*NumOfFeatures);
residual = zeros(2*nnz(PoseGraphMatrix_r), 1);
meascount = 0;
last_meas_count = 0;

featcount = 0;
last_feat_count = 0;
% Hessian build up
A22_inv = zeros(3*NumOfFeatures);

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
       
      [Jpose(2*meascount+1:2*meascount+2, 6*(poses_see_featk(l+1)-1)+1:6*(poses_see_featk(l+1)-1)+6),...
          residual(2*meascount+1:2*meascount+2)] = ...
           dprojection_dp_dtt2(Cr_p_f, Cl_T_Cr(:,1:3), Cl_T_Cr(:,4), Cl_z_f);       
       
      Jfeat_s = 1/Cl_p_f(3)*[1 0 -Cl_p_f(1)/Cl_p_f(3);0 1 -Cl_p_f(2)/Cl_p_f(3)] * Cl_T_Cr(:,1:3);
      Jfeat(2*meascount+1:2*meascount+2, 3*featcount+1:3*featcount+3) = Jfeat_s;              
      
      meascount = meascount + 1;
   end   
   
   Jfeatk = Jfeat(last_meas_count+1:2*meascount, 3*featcount+1:3*featcount+3);
   A22_inv(3*featcount+1:3*featcount+3, 3*featcount+1:3*featcount+3) = (Jfeatk'*Jfeatk + regulator*eye(3)) \ eye(3);
      
   featcount = featcount + 1;
   last_meas_count = meascount;
end

norm_iter = norm(residual);
normHist = [normHist norm_iter];
if iter > 20
    norm_iter = 1e2*abs(normHist(end-1) - normHist(end));
end

Jpose_sp = sparse(Jpose);
Jfeat_sp = sparse(Jfeat);
A11 = Jpose_sp'*Jpose_sp;
A12 = Jpose_sp'*Jfeat_sp;
A22_inv_sp = sparse(A22_inv);
res_ = [Jpose_sp'*residual;Jfeat_sp'*residual];

% solve:
p_tilde = (A11 +regulator*speye(size(A11))- A12*A22_inv_sp*A12') \ (res_(1:6*NumOfPoses) - A12*A22_inv_sp*res_(1+6*NumOfPoses:end));

for l = 1:NumOfPoses-1
    dq = [0.5*p_tilde(6*l+4:6*l+6);1];
    dq = dq / norm(dq);
    AbsolutePoses(:,1:3,usedPoses(l+1)) = quat2rot( quat_mul(dq, rot2quat(AbsolutePoses(:,1:3,usedPoses(l+1)))) );       
    AbsolutePoses(:,4,usedPoses(l+1)) = AbsolutePoses(:,4,usedPoses(l+1)) + p_tilde(6*l+1:6*l+3);
end

grad_norm = norm(p_tilde)^2;
 
% 4. Efficient solve for each feature f
bf_p = res_(1+6*NumOfPoses:end) - A12'*p_tilde;
f_tilde = A22_inv_sp * bf_p;
FeatureBag(1:3, featureAvailableList+1) = FeatureBag(1:3, featureAvailableList+1) + reshape(f_tilde, [3, NumOfFeatures]);

grad_norm = norm(f_tilde)^2 + norm(p_tilde)^2;
gradnormHist = [gradnormHist sqrt(grad_norm)];

% 6. Visualization
% figure(3); clf; hold on;
% Draw3DWorld(AbsolutePoses, FeatureBag, [eye(3) zeros(3,1)], usedPoses);
% grid on;



    fprintf('Iter :%d, norm: %f, gradnorm: %f \n', iter, norm_iter, grad_norm);
    iter = iter + 1;
end
