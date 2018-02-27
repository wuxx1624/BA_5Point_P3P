function [PoseGraphMatrix, pointerToImages] = ConstructPoseGraphMatrix(featureExtracted, PairWiseMatches)

NumberOfImages = length(featureExtracted);
threshold_matches = 30;

% Incremental map builder
total_features_extracted = 0;
for i = 1:NumberOfImages
    total_features_extracted = total_features_extracted + size(featureExtracted{i},2);
end

PoseGraphMatrix = zeros(total_features_extracted, NumberOfImages, 'uint16');
pointerToImages = zeros(NumberOfImages+1,1);
for i = 1:NumberOfImages
    if i == 1
        new_features_Imagei = 1:size(featureExtracted{i},2);
        PoseGraphMatrix(1+pointerToImages(i):size(featureExtracted{i},2)+pointerToImages(i),i) = ...
            (1:size(featureExtracted{i},2))';
    else
        % Taking new features in images i that haven't seen before
        new_features_Imagei = setdiff(1:size(featureExtracted{i},2), unique(PoseGraphMatrix(1:pointerToImages(i),i)));
        PoseGraphMatrix(1+pointerToImages(i):length(new_features_Imagei)+pointerToImages(i),i) = new_features_Imagei;
    end
    
    
    for j = i+1:NumberOfImages
        matches_ij = PairWiseMatches{SUTrilInd(NumberOfImages,i,j)};     
%         if size(matches_ij,1) < threshold_matches
%             continue;
%         end
        [idx_ij_of_new_features_a, idx_ij_of_new_features_b] = ismember(matches_ij(:,5), new_features_Imagei);
        PoseGraphMatrix(idx_ij_of_new_features_b(idx_ij_of_new_features_b~=0) + pointerToImages(i), j) = ...
            matches_ij(idx_ij_of_new_features_a,10);
        
        for k = j+1:NumberOfImages            
            matches_jk = PairWiseMatches{SUTrilInd(NumberOfImages,j,k)};
            
            [match_k_in_j_first_seen_in_ia, match_k_in_j_first_seen_in_ib] = ...
                ismember(matches_jk(:,5), PoseGraphMatrix(matches_ij(:,5) + pointerToImages(i), j));   
            % No further/longer match found
            if length(match_k_in_j_first_seen_in_ia) ~= length(match_k_in_j_first_seen_in_ib) ...
                    || isempty(match_k_in_j_first_seen_in_ib)
                continue;
            end
            
            PoseGraphMatrix(match_k_in_j_first_seen_in_ib(match_k_in_j_first_seen_in_ib~=0) + pointerToImages(i), k) = ...
                matches_jk(match_k_in_j_first_seen_in_ia,10);
        end
%         spy(PoseGraphMatrix);
%         pause;
    end
    % clean up all zero rows
    matches_row = find( ...
        sum(PoseGraphMatrix(pointerToImages(i)+1:pointerToImages(i) + length(new_features_Imagei), i+1:end), 2) ~= 0);
    new_indices_imagei = (1:length(new_features_Imagei))';
    
    % clear the rest of the mess just created
    clean_up_indices = setdiff(new_indices_imagei,matches_row);
    PoseGraphMatrix(clean_up_indices + pointerToImages(i), :) = ...
        zeros(length(clean_up_indices), NumberOfImages);

    % extract only rows with available matches
    PoseGraphMatrix(1+pointerToImages(i):length(matches_row)+pointerToImages(i), :) = ...
        PoseGraphMatrix(matches_row + pointerToImages(i), :);
    
    pointerToImages(i+1) = pointerToImages(i) + length(matches_row);
end

PoseGraphMatrix(pointerToImages(NumberOfImages+1)+1:end, :) = [];

% Last clean up -> upper traingular form:
for i = 1:NumberOfImages
    PoseGraphMatrix(1+pointerToImages(i+1):end,i) = zeros(pointerToImages(NumberOfImages+1)-pointerToImages(i+1),1);
end

%% TODO: refine multiple matches that potentially contain outliers
