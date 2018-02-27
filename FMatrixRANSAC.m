function [Fbest, bestInlier, bestGeneratorSet, iter] = FMatrixRANSAC( matches, MaxIter, inlier_ratio, threshold )

% RANSAC parameter
Fbest = []; 
% MaxIter = 300;
MaxInlier = inlier_ratio*size(matches,1);
% threshold = 3;
% threshold = 5e-2;
bestInlier = [];
bestGeneratorSet = [];
iter = 1;
while iter < MaxIter
    
    % Check uniqueness
    % Picking 8 distinct points
    picked_points_idx = randi(size(matches,1),1,8);
    picked_points_idx = unique(picked_points_idx);
    while length(picked_points_idx) < 8
        picked_points_idx = [picked_points_idx randi(size(matches,1),1,8 - length(picked_points_idx))];
        picked_points_idx = unique(picked_points_idx);
    end
    
    % Check disparsity
    

    % Compute fundamental matrix
    F = ComputeFundamentalMatrix([matches(picked_points_idx,1:2) matches(picked_points_idx,6:7)]);
    if(size(F,1) == 0)
        continue;
    end
    
    % Compute inlier
    inlierIdx = ComputeInlierEpipolarDistance(F, matches, threshold);
    % find inlier Sampson Epipolar error
%     inlierIdx = [];
%     for i = 1:size(matches,1)
%         if(abs([matches(i,5:6) 1]*F*[matches(i,1:2) 1]') < threshold)            
%             inlierIdx = [inlierIdx i];
%         end
%     end
    
    if length(inlierIdx) > MaxInlier
        Fbest = F;
        %ComputeFundamentalMatrixInlier...
        % ([matches(inlierIdx,1:2) matches(inlierIdx,5:6)]);
        % Rbest = R;
        % tbest = t;
        bestInlier = inlierIdx;
        bestGeneratorSet = picked_points_idx;
        break;
    else
        if length(inlierIdx) > length(bestInlier)
            Fbest = F;
            % ComputeFundamentalMatrixInlier...
            %([matches(inlierIdx,1:2) matches(inlierIdx,5:6)]);
            % Rbest = R;
            % tbest = t;
            bestInlier = inlierIdx;
            bestGeneratorSet = picked_points_idx;
        end
    end 
    
    iter = iter + 1;
end