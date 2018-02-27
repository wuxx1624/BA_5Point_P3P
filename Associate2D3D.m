function [P,w] = Associate2D3D(FeatureTriangulated, PairWiseMatches, currentIdx, previousIdxs, N)

TriangulatedPairs = sort(previousIdxs);
P = []; w = [];
for i = 1:length(TriangulatedPairs)-1
    for j = i+1:length(TriangulatedPairs)
        current_i = SUTrilInd(N,currentIdx,i);
        current_j = SUTrilInd(N,currentIdx,j);
        i_j = SUTrilInd(N,i,j);
        for k = 1:size(FeatureTriangulated{i_j},2)
            li_ = find(PairWiseMatches{current_i}(:,5) == FeatureTriangulated{i_j}(4,k));
            lj_ = find(PairWiseMatches{current_j}(:,5) == FeatureTriangulated{i_j}(5,k));
            if length(li_) == 0 || length(lj_) == 0
                continue;
            end
            li = li_(1); lj = lj_(1);
            if PairWiseMatches{current_i}(li,10) == PairWiseMatches{current_j}(lj,10)
                P = [P FeatureTriangulated{i_j}(:,k)];
                w = [w; PairWiseMatches{current_i}(li,6:10) PairWiseMatches{current_i}(li,1:5) PairWiseMatches{current_j}(lj,1:5)];
            end
        end        
    end    
end


