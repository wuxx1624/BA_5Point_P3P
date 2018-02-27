function [inlier_idx] = ComputeInlier3D2DReprojection(P, w, K, R, t, threshold)

w_ = K*[R t]*[P(1:3,:);ones(1,size(P,2))];
w_(1,:) = w_(1,:)./w_(3,:);
w_(2,:) = w_(2,:)./w_(3,:);

dw = w_(1:2,:)' - w(:,1:2);
error_mat = sqrt(dw(:,1).*dw(:,1) + dw(:,2).*dw(:,2));

inlier_idx = find(error_mat < threshold);

