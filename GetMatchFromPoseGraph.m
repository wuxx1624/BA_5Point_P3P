function [w, idx] =  GetMatchFromPoseGraph(PoseGraphMatrix, indices_C, A, B, fA, fB)

% extract matches in A which based on indices_C and correponding to B
indices_AB = intersect(find(PoseGraphMatrix(:,A) ~= 0), find(PoseGraphMatrix(:,B) ~= 0));
indices_ABC = intersect(indices_AB, indices_C);

w = [fA(:,PoseGraphMatrix(indices_ABC,A))' double(PoseGraphMatrix(indices_ABC,A)) ...
     fB(:,PoseGraphMatrix(indices_ABC,B))' double(PoseGraphMatrix(indices_ABC,B))];
w = double(w);
idx = indices_ABC;