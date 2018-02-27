%% This script will create visual BA just for testing
clear all; close all;
addpath('../');
addpath('../Drawings/');
addpath('../JacobianModule/');

NumOfFeatures = 500;
NumOfPoses = 15;

PoseGraphMatrix = randi([0 1], NumOfFeatures, NumOfPoses);


%% 3D-WORLD creation

% Make sure every feature is observed by more than 2 poses
addMorePoseIdx = find(sum(PoseGraphMatrix,2) <= 5);
for k = addMorePoseIdx'
    listIdx = find(PoseGraphMatrix(k,:) == 0);
    PoseGraphMatrix(k, listIdx(1:2)) = 1;
end
figure(1);
spy(PoseGraphMatrix);


% To ensure f.o.v = 60 deg.
fov = 60*pi/180;
Rfeat = 1; Rcam = Rfeat/sin(fov/2);


% Random pose generator
W_T_C1 = [];
AbsolutePoses_true = zeros(3,4,NumOfPoses);
%figure(2); hold on;
for k = 1:NumOfPoses
    tt = 2*pi*rand();
    WRC = [-sin(tt) 0 -cos(tt); cos(tt) 0 -sin(tt); 0 -1 0];
    WpC = Rcam*[cos(tt);sin(tt);0];
    if k == 1
        W_T_C1 = [WRC WpC];
        AbsolutePoses_true(:,:,k) = [eye(3) zeros(3,1)];
    else
        C_T_W = [WRC' -WRC'*WpC];
        AbsolutePoses_true(:,:,k) = C_T_W*[W_T_C1;zeros(1,3) 1];
    end
%     DrawAxis(WpC, WRC, 0.2);
%     DrawCamera(WpC, WRC, 0.2, 'k');
end
% tt = 0:0.1:2*pi;
% plot3(Rcam*cos(tt),Rcam*sin(tt),0*tt);
% grid on;


% Random feature generator
FeatureBag_true = zeros(4,NumOfFeatures);
%figure(2); hold on;
for k = 1:NumOfFeatures
    tt = 2*pi*rand();
    phi = pi*rand() - pi/2;
    FeatureBag_true(:,k) = [Rfeat*[cos(phi)*cos(tt);cos(phi)*sin(tt);sin(phi)];1]; % <--- w.r.t the world
    % plot3(FeatureBag_true(1,k),FeatureBag_true(2,k),FeatureBag_true(3,k), 'bo', 'markerSize', 3, 'markerfacecolor', 'b');
    FeatureBag_true(1:3,k) = InversePose(W_T_C1)*[FeatureBag_true(1:3,k);1];
end

% Create feature projections as measurement
featureExtracted_true = cell(NumOfPoses, 1);
for k = 1:NumOfPoses    
    numFeatPosek = find(PoseGraphMatrix(:,k));
    featureExtracted_true{k} = zeros(2, length(numFeatPosek));
    for l = 1:length(numFeatPosek)
        % assign correct feature index for pose graph matrix
        PoseGraphMatrix(numFeatPosek(l),k) = l; 
        
        % Pose Transformation
        Ck_T_W = AbsolutePoses_true(:,:,k);
        Cl_T_W = AbsolutePoses_true(:,:,FeatureBag_true(4,numFeatPosek(l)));
        Ck_T_Cl = Ck_T_W*[InversePose(Cl_T_W);zeros(1,3) 1];
        
        % Image Projection
        Ck_p_fl = Ck_T_Cl*[FeatureBag_true(1:3,numFeatPosek(l));1];
        featureExtracted_true{k}(:,l) = Ck_p_fl(1:2)./Ck_p_fl(3);
    end
    
%     figure(4); clf;
%     plot(featureExtracted{k}(1,:), featureExtracted{k}(2,:), 'b*');
%     pause;
end

%% Poses and features positions pertubation
sigma_f = 0.2;
sigma_p = 0.4;
sigma_tt = 0.01;
sigma_z = 1e-4;

FeatureBag = zeros(4,NumOfFeatures);
for k = 1:NumOfFeatures
    df_k = randn(3,1)*sigma_f;
    FeatureBag(1:3,k) = FeatureBag_true(1:3,k) + df_k;   
    FeatureBag(4,k) = FeatureBag_true(4,k);
end

AbsolutePoses = AbsolutePoses_true;
for k = 2:NumOfPoses
    dp = randn(3,1)*sigma_p;
    dtt = randn(3,1)*sigma_tt;
    AbsolutePoses(:,4,k) = AbsolutePoses_true(:,4,k) + dp;
    dq = [0.5*dtt;1]; dq = dq/norm(dq);
    AbsolutePoses(:,1:3,k) = quat2rot(dq)*AbsolutePoses_true(:,1:3,k);
end
figure(2); hold on;
Draw3DWorld(AbsolutePoses_true, FeatureBag_true, W_T_C1)
tt = 0:0.1:2*pi;
plot3(Rcam*cos(tt),Rcam*sin(tt),0*tt);
grid on;

featureExtracted = cell(NumOfPoses, 1);
for k = 1:NumOfPoses
   featureExtracted{k} = featureExtracted_true{k} + sigma_z*randn(size(featureExtracted_true{k}));   
%     figure(4); clf;
%     plot(featureExtracted_true{k}(1,:), featureExtracted_true{k}(2,:), 'b*'); hold on;
%     plot(featureExtracted{k}(1,:), featureExtracted{k}(2,:), 'ro'); hold on;
%     pause;
end

%% Bundle Adjustment
regulator = 1e-1;
norm_iter = 100;
normHist = [];
maxIter = 100;
iter = 1;

while iter < maxIter && norm_iter > 1e-5

fprintf('----------------------------------------\n');
fprintf('Iter :%d, norm: %f \n', iter, norm_iter);
% 1. Build Jacobian Jp Jf corresponds to each feature
Jpose = zeros(2*nnz(PoseGraphMatrix), 6*NumOfPoses);
Jfeat = zeros(2*nnz(PoseGraphMatrix), 3*NumOfFeatures);
residual = zeros(2*nnz(PoseGraphMatrix), 1);

% Hessian build up
A22_inv = zeros(3*NumOfFeatures);


% fprintf('Building all Jacobian \n');
% tic();
% Build sparse structure
row_count = 0;
last_row_count = 0;
for featk = 0:size(PoseGraphMatrix,1)-1
   
   poses_see_featk = find(PoseGraphMatrix(featk+1, :));   
   
   for l = 0:length(poses_see_featk)-1

        %  Pose Jacobian
        Cr_p_f = FeatureBag(1:3 ,featk+1);      
        Cl_z_f = featureExtracted{poses_see_featk(l+1)}(:,PoseGraphMatrix(featk+1,poses_see_featk(l+1)));                  


        % Pose Transformation
        Cl_T_W = AbsolutePoses(:,:,poses_see_featk(l+1));
        Cr_T_W = AbsolutePoses(:,:,FeatureBag(4, featk+1));      
        Cl_T_Cr = Cl_T_W*[InversePose(Cr_T_W);zeros(1,3) 1];            
        Cl_p_f = Cl_T_Cr*[Cr_p_f;1];                   


        [Jpose(row_count+1:row_count+2, 6*(poses_see_featk(l+1)-1)+1:6*(poses_see_featk(l+1)-1)+6),...
          residual(row_count+1:row_count+2)] = ...
           dprojection_dp_dtt2(Cr_p_f, Cl_T_Cr(:,1:3), Cl_T_Cr(:,4), Cl_z_f);      

        Jfeat_s = 1/Cl_p_f(3)*[1 0 -Cl_p_f(1)/Cl_p_f(3);0 1 -Cl_p_f(2)/Cl_p_f(3)] * Cl_T_Cr(:,1:3);
        Jfeat(row_count+1:row_count+2, 3*featk+1:3*featk+3) = Jfeat_s;
                
        A22_inv(3*featk+1:3*featk+3, 3*featk+1:3*featk+3) = A22_inv(3*featk+1:3*featk+3, 3*featk+1:3*featk+3) + (Jfeat_s'*Jfeat_s);
        row_count = row_count + 2;
   end   
   
   % residual(last_row_count+1:row_count) = Jpose_Jfeat_res{featk+1,3};
   A22_inv(3*featk+1:3*featk+3, 3*featk+1:3*featk+3) = (A22_inv(3*featk+1:3*featk+3, 3*featk+1:3*featk+3) + regulator*eye(3)) \ eye(3);
   last_row_count = row_count;
   
   % figure(6);
   % spy([Jpose Jfeat]);
   % pause;
end
% toc();

norm_iter = norm(residual);
normHist = [normHist norm_iter];

% fprintf('Building Hessian \n');
% tic();
Jpose_sp = sparse(Jpose);
Jfeat_sp = sparse(Jfeat);
A11 = Jpose_sp'*Jpose_sp;
A12 = Jpose_sp'*Jfeat_sp;
res_ = [Jpose_sp'*residual;Jfeat_sp'*residual];
% toc();


% A = [A11 A12;A12' A22];
% A = A + regulator*eye(size(A));
% x_tilde = A \ res_;
% fprintf('Create sparse A22 inverse \n');
% tic();
A22_inv_sp = sparse(A22_inv);
% toc();

% solve:
p_tilde = (A11+regulator*speye(6*NumOfPoses) - A12*A22_inv_sp*A12') \ (res_(1:6*NumOfPoses) - A12*A22_inv*res_(1+6*NumOfPoses:end));

% p_tilde = x_tilde(1:6*NumOfPoses);
% fprintf('State update \n');
% tic();
 for l = 1:NumOfPoses-1
    dq = [0.5*p_tilde(6*l+4:6*l+6);1];
    dq = dq / norm(dq);
    AbsolutePoses(:,1:3,l+1) = quat2rot( quat_mul(dq, rot2quat(AbsolutePoses(:,1:3,l+1))) );       
    AbsolutePoses(:,4,l+1) = AbsolutePoses(:,4,l+1) + p_tilde(6*l+1:6*l+3);
end
% toc();

% 4. Efficient solve for each feature f
bf_p = res_(1+6*NumOfPoses:end) - A12'*p_tilde;
FeatureBag(1:3, :) = FeatureBag(1:3, :) + reshape(A22_inv_sp * bf_p, [3, NumOfFeatures]);


% 6. Visualization
figure(3); clf; hold on;
Draw3DWorld(AbsolutePoses, FeatureBag, W_T_C1)
tt = 0:0.1:2*pi;
plot3(Rcam*cos(tt),Rcam*sin(tt),0*tt);
grid on;
% pause(0.5);

iter = iter + 1;
end

figure(3); clf; hold on;
Draw3DWorld(AbsolutePoses, FeatureBag, W_T_C1)
tt = 0:0.1:2*pi;
plot3(Rcam*cos(tt),Rcam*sin(tt),0*tt);
grid on;
pause(0.5);

figure(4);
semilogy(normHist, 'b-*');

% % 2. Build A11 = Jp'*Jp, A12 = Jp'*Jf, and A22^-1 = inv(Jf'*Jf) for every
% % single block
% App = zeros(6*NumOfPoses, 6*NumOfPoses);
% Afp = zeros(3*NumOfFeatures, 6*NumOfPoses);
% Aff_inv = zeros(3*NumOfFeatures, 3); % Block diagonal structure
% bp = zeros(6*NumOfPoses,1);
% bf = zeros(3*NumOfFeatures,1);
% 
% for featk = 1:size(PoseGraphMatrix,1)
%     % a. Extract poses that see this feature
%     poses_see_featk = find(PoseGraphMatrix(featk, :));
%     
%     % b. Build pose Hessian A11
%     Jp_compact = Jpose_Jfeat_res{featk, 1};
%     App_compact = Jp_compact' * Jp_compact;
%     App = AccumulateCompactToFull(App, App_compact, poses_see_featk, poses_see_featk);
% 
%     % c. Build Schur complement part and A22_inv
%     Jfk = Jpose_Jfeat_res{featk, 2};
%     Aff_inv(3*(featk-1)+1:3*(featk-1)+3, :) = (Jfk'*Jfk + regulator*eye(size(3))) \ eye(3);
%     Jpfk = zeros(6*length(poses_see_featk), 3);
%     for l = 0:length(poses_see_featk)-1
%         Jpfk(6*l+1:6*l+6,:) = Jp_compact(:, 6*l+1:6*l+6)'*Jfk(2*l+1:2*l+2, :);
%     end
%     App_schur_compact = Jpfk * Aff_inv(3*(featk-1)+1:3*(featk-1)+3, :) * Jpfk';    
%     App = AccumulateCompactToFull(App, -App_schur_compact, poses_see_featk, poses_see_featk);
% 
%     % d. Keep Jpfk+zero extension structure
%     Afp = AccumulateCompactToFull(Afp, Jpfk', featk, poses_see_featk);     
% 
%     % e. build Residual J'*r
%     bf(3*(featk-1)+1:3*(featk-1)+3) = Jpose_Jfeat_res{featk,2}'*Jpose_Jfeat_res{featk,3};
%     for l = 0:length(poses_see_featk)-1
%         bp(6*l+1:6*l+6) = Jp_compact(:, 6*l+1:6*l+6)'*Jpose_Jfeat_res{featk,3}(2*l+1:2*l+2);
%     end
% end
% 
% 
% % 3. Solve for pose and update
% bf_ = bf;
% for featk = 0:size(PoseGraphMatrix,1)-1
%    bf_(3*featk+1:3*featk+3) = Aff_inv(3*(featk)+1:3*(featk)+3, :)*bf_(3*featk+1:3*featk+3);
% end
% bp_f = bp - Afp'*bf_;
% p_tilde = (App + regulator*eye(size(App))) \ bp_f;
% for l = 0:length(poses_see_featk)-1
%     dq = [1/2*p_tilde(6*l+1:6*l+3);1];
%     dq = dq / norm(dq);
%     AbsolutePoses(:,1:3,l+1) = quat2rot( quat_mul(dq, rot2quat(AbsolutePoses(:,1:3,l+1))) );
%     AbsolutePoses(:,4,l+1) = AbsolutePoses(:,4,l+1) + p_tilde(6*l+4:6*l+6);
% end
% 
% % 4. Efficient solve for each feature f
% bf_p = bf - Afp*p_tilde;
% for featk = 0:size(PoseGraphMatrix,1)-1
%    FeatureBag(1:3,featk+1) = FeatureBag(1:3,featk+1) + Aff_inv(3*(featk)+1:3*(featk)+3, :) \ bf_p(3*featk+1:3*featk+3);
% end
% 
% % 6. Visualization
% pause;
% figure(3); hold on;
% Draw3DWorld(AbsolutePoses, FeatureBag)
% tt = 0:0.1:2*pi;
% plot3(Rcam*cos(tt),Rcam*sin(tt),0*tt);
% grid on;
% 
