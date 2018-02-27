function [PoseGraphMatrix, pointerToImages] = ConstructPoseGraphMatrix2(featureExtracted, PairWiseMatches)

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
%         idx_features_imagei_seen_b4_based_preimgs = find(PoseGraphMatrix(1:pointerToImages(i),i) ~= 0);
        idx_features_imagei_seen_b4_based_preimgs = PoseGraphMatrix(1:pointerToImages(i),i) ~= 0;
        features_imagei_seen_b4 = PoseGraphMatrix(idx_features_imagei_seen_b4_based_preimgs,i);
        
        new_features_Imagei = setdiff(1:size(featureExtracted{i},2), features_imagei_seen_b4);
        PoseGraphMatrix(1+pointerToImages(i):length(new_features_Imagei)+pointerToImages(i),i) = new_features_Imagei';
    end
    
    
    for j = i+1:NumberOfImages        
        matches_ij = PairWiseMatches{SUTrilInd(NumberOfImages,i,j)};     
%         if size(matches_ij,1) < threshold_matches
%             continue;
%         end

        % new features in image i to other images that haven't considered!
        [idx_ij_of_new_features_a, idx_ij_of_new_features_b] = ismember(matches_ij(:,5), new_features_Imagei);
        PoseGraphMatrix(idx_ij_of_new_features_b(idx_ij_of_new_features_b~=0) + pointerToImages(i), j) = ...
            matches_ij(idx_ij_of_new_features_a,10);
        
        
        
        % old features in image i that are seen before
        % clean up and look for agreement, otherwise it's a potential outlier and
        % should be removed!
        if i > 1                    
            for l = unique(features_imagei_seen_b4)'
                vp = find(PoseGraphMatrix(:,i) == l);
                % unify this common features that are seen in different
                % places but not pairwise common
                if(length(vp) > 1)
                % 1. Find the strongest support
                    [~, most_poses] = max(sum(PoseGraphMatrix(vp,:)~=0,2));
                % 2. Extract non-zero indices
                    nz_indices = find(PoseGraphMatrix(vp(most_poses),:) ~= 0);
                    nz_indices = setdiff(nz_indices, i);
                % 3. Scan through rest, if disagree, remove, else unify
                    % Extract elements in rows of other vp, but only take
                    % cols corresponds to nonzeros vp(most_pose)
                    potential_unifying_rows = setdiff(1:length(vp),most_poses);
                    coeff_matrix = PoseGraphMatrix(vp(potential_unifying_rows), nz_indices);
                    % Either 0
                    status_matrix1 = (coeff_matrix == 0);
                    % Or equal to most_pose row
                    status_matrix2 = (coeff_matrix == repmat(PoseGraphMatrix(vp(most_poses),nz_indices), [length(vp)-1,1]));
                    % 
                    status_matrix = ~(status_matrix1 | status_matrix2);                    
                    
                    % TODO: This row+col checking logic can definitely
                    % be reduced much more
                    % checking if they're all agree to update
                    for kkk = potential_unifying_rows(sum(status_matrix,2) == 0)
                        nz_indices_kkk = setdiff(find(PoseGraphMatrix(vp(kkk), :) ~= 0),i);
                        for kkkk = nz_indices_kkk
                            ss_row1 = (PoseGraphMatrix(vp(setdiff(1:length(vp),most_poses)),kkkk) == 0);
                            ss_row2 = (PoseGraphMatrix(vp(setdiff(1:length(vp),most_poses)),kkkk) == PoseGraphMatrix(vp(kkk),kkkk));
                            ss_row = ~(ss_row1 | ss_row2);
                            if sum(ss_row) == 0
                               PoseGraphMatrix(vp(most_poses), kkkk) =  PoseGraphMatrix(vp(kkk), kkkk);
                            end
                        end
                    end
                    
                    % Extract the one after unifying
                    idx_reduced_l = vp(most_poses);
                    % Clear up after unifying
                    PoseGraphMatrix(vp(setdiff(1:length(vp),most_poses)),:) = 0;
                else
                    idx_reduced_l = vp;
                end

                vm = find(matches_ij(:,5)== l);
                if isempty(vm) 
                    continue;
                end
                
                if PoseGraphMatrix(idx_reduced_l,j) == 0
                    PoseGraphMatrix(idx_reduced_l,j) = matches_ij(vm,10);     
                elseif matches_ij(vm,10) ~= PoseGraphMatrix(idx_reduced_l,j)
                    % Remove/modify?
                    if sum(PoseGraphMatrix(idx_reduced_l,:)~=0) >= 4
                        % do nothing since the matching is already across a
                        % lot of poses
                        PoseGraphMatrix(idx_reduced_l,i) = 0;
                    else
                        % remove this row from previous and replace this by
                        % the new feature from i                        
                        PoseGraphMatrix(idx_reduced_l,:) = zeros(1,NumberOfImages);
                        PoseGraphMatrix(length(new_features_Imagei)+pointerToImages(i)+1,i) = matches_ij(vm,5);
                        PoseGraphMatrix(length(new_features_Imagei)+pointerToImages(i)+1,j) = matches_ij(vm,10);
                        new_features_Imagei = [new_features_Imagei matches_ij(vm,5)];
                    end
                end
            end                        
        end        
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
    
    % Update pointers
    pointerToImages(i+1) = pointerToImages(i) + length(matches_row);

    
    % upper triangular update
    PoseGraphMatrix(1+pointerToImages(i+1):end,:) = 0;
end

PoseGraphMatrix(pointerToImages(NumberOfImages+1)+1:end, :) = [];
