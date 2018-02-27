function   [FeaturesBag] = MultiPoseTriangulation(FeaturesBag,           ...
                                                 AbsolutePoses,         ...
                                                 PoseGraphMatrix,       ...
                                                 newPose, usedPoses,    ...
                                                 featureExtracted,      ...                                         
                                                 CameraParams, visualization, I_rgb)
newposeIdx = find(PoseGraphMatrix(:, newPose));
usedposeIdx = find(sum(PoseGraphMatrix(:, usedPoses),2));
triangulateIdx = intersect(newposeIdx, usedposeIdx);
triangulateIdx = setdiff(triangulateIdx,find(FeaturesBag(4,:)~=0));


% MultiPoses linear triangulation
allPoses = [usedPoses, newPose];
for k = triangulateIdx'
    % 1. Find how many poses see this feature
    posesViewkfeatIdx = find(PoseGraphMatrix(k,allPoses));    
    
    A = zeros(3*length(posesViewkfeatIdx), 4);
    count = 0;
    for l = posesViewkfeatIdx
        featureIdxImagel = PoseGraphMatrix(k, allPoses(l));
        Cl_T_Cn = AbsolutePoses(:,:,allPoses(l))*[AbsolutePoses(:,1:3,newPose)' -AbsolutePoses(:,1:3,newPose)'*AbsolutePoses(:,4,newPose);zeros(1,3) 1];
        A(3*count+1:3*count+3, :) = skewsymm(CameraParams.Kinv*[featureExtracted{allPoses(l)}(1:2,featureIdxImagel);1])*Cl_T_Cn;        
        count = count + 1;
    end
    
%     VisualizeMatches(I_rgb{allPoses(posesViewkfeatIdx(1))},I_rgb{allPoses(posesViewkfeatIdx(2))},[featureExtracted{allPoses(posesViewkfeatIdx(1))}(:,PoseGraphMatrix(k, allPoses(posesViewkfeatIdx(1))))' ...
%         0 featureExtracted{allPoses(posesViewkfeatIdx(2))}(:,PoseGraphMatrix(k, allPoses(posesViewkfeatIdx(2))))' 0], 'b', 'g');
%     
    [~,~,V] = svd(A);
    temp_feat = [V(1:3,end)/V(4,end);newPose];
    if( ChieralityCheckMultiPoses(temp_feat, AbsolutePoses, allPoses(posesViewkfeatIdx)) )
        % Check reprojection error:
        triangulate_decision = 1;
        for l = posesViewkfeatIdx
            featureIdxImagel = PoseGraphMatrix(k, allPoses(l));
            Cl_T_Cn = AbsolutePoses(:,:,allPoses(l))*[AbsolutePoses(:,1:3,newPose)' -AbsolutePoses(:,1:3,newPose)'*AbsolutePoses(:,4,newPose);zeros(1,3) 1];
            Cl_p_f = Cl_T_Cn*[temp_feat(1:3);1];
            Cl_p_f = Cl_p_f(1:3)/Cl_p_f(3);
            reproj_error = norm(CameraParams.K(1:2,:)*Cl_p_f - ...
                featureExtracted{allPoses(l)}(1:2,featureIdxImagel));
            if reproj_error > 30
                triangulate_decision = 0;
                break;
            end
        end
        if triangulate_decision == 1
            FeaturesBag(:,k) = temp_feat;
        end
    end
end


if visualization == 1
    Ci_P = FeaturesBag(:,FeaturesBag(4,:)==newPose);
    Ci_z = featureExtracted{newPose}(1:2,PoseGraphMatrix(FeaturesBag(4,:)==newPose, newPose))';
    reprojection_matches = ReprojectionToMatches(Ci_P, eye(3), zeros(3,1), CameraParams);
	ReprojectionVisualization2(I_rgb{newPose}, Ci_z, reprojection_matches); 
end
