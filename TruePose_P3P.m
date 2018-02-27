function [R_true,t_true,idx]=TruePose_P3P(R,t,p,b,K)

%%%%%%%% input %%%%%%%%%%%%%%%%%%%%%%
% p: 3D features, w.r.t global frame
% b: 2D features, 3X1 vector, w.r.t camera frame
% K: intrinsic matrix

if size(b,1) ~= 3
    b = [b;ones(1,size(b,2))];
end

error = zeros(size(R,3),1);
% b = b/norm(b);
for i=1:size(R,3)
    b_est = K*R(:,:,i)'*p - K*R(:,:,i)'*t(:,i);
%     b_est = b_est/norm(b_est);
    b_est = b_est./b_est(3,:);
    error(i) = norm(b_est-b);
end
[~,idx] = min(error);
R_true = R(:,:,idx)';
t_true = -R(:,:,idx)'*t(:,idx);