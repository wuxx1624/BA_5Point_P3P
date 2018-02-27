function [matches] = BFNN2directions(f1,d1,f2,d2)
% Matches convention: 
% match(rowi) = [f1_pixel_position 2df1_idx f2_pixel_position 2df2_idx]
[idx1_2, dist1_2] = knnsearch(d1', d2', 'k', 5);
[idx2_1, dist2_1] = knnsearch(d2', d1', 'k', 5);

matches = [];
M = 5;
for k = 1:size(idx1_2,1)
    if min(dist1_2(k,:)) > 300 || dist1_2(k,1)/dist1_2(k,2) > 0.7
        continue;
    end
    
    % mindist = inf;
    getmatch = 0;
    for j = 1:M
        [~, locb] = ismember(idx2_1(idx1_2(k,j),:), k);
        locb_nz = idx2_1(idx1_2(k,j),locb~=0); %(locb~=0);                
        
        
        if length(locb_nz) >= 1
            idx = find(locb~=0);
            if dist2_1(idx1_2(k,j),1)/dist2_1(idx1_2(k,j),2) < 0.7
                matches = [matches; [f1(:,idx1_2(k,j))' idx1_2(k,j) f2(:,k)' k]];
                getmatch = 1;
            end
        end
        
        if getmatch == 1
            break;
        end
    end        
end

[~,uindex,~] = unique(matches(:,5),'stable');
matches = matches(uindex,:);