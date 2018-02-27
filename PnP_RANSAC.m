function [Rbest,tbest,bestInlier,bestGeneratorSet] = PnP_RANSAC(P, w, CameraParams, MaxIter, inlier_ratio, threshold)

% RANSAC parameter
Rbest = []; tbest = [];
MaxInlier = inlier_ratio*size(w,1);
bestInlier = [];
bestGeneratorSet = [];
iter = 1;

while iter < MaxIter
    
    % Check uniqueness
    % Picking 6 distinct points
    picked_points_idx = randi(size(w,1),1,6);
    picked_points_idx = unique(picked_points_idx);
    while length(picked_points_idx) < 6
        picked_points_idx = [picked_points_idx randi(size(w,1),1,6 - length(picked_points_idx))];
        picked_points_idx = unique(picked_points_idx);
    end
    
    % Compute pose
    [R,t] = LinearPnP(P(:,picked_points_idx), w(picked_points_idx, :), CameraParams);
    
    % Compute inlier based on reprojection
    inlierIdx = ComputeInlier3D2DReprojection(P, w, CameraParams, R, t, threshold);

    
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


% Least square PnP
% [Rbest,tbest] = LinearPnP(P(:,bestInlier), w(bestInlier, :), CameraParams);
