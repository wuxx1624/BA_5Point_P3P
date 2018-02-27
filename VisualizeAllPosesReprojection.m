function VisualizeAllPosesReprojection(PoseGraphMatrix, AbsolutePoses, FeaturesBag, featureExtracted, usedPoses, newPose, I_rgb, CameraParams)

allposes = [usedPoses, newPose];
I_full = [];
for k = allposes
    I_full = [I_full I_rgb{k}];
end

figure; imshow(I_full); hold on;
featureIdx_triangulated = find(FeaturesBag(4,:)~=0);

countImage = 0;
for k = allposes
    % 1. Extract feature index visible in this pose
    visiblefeaturesIdx = find(PoseGraphMatrix(:,k));
    visiblefeaturesIdx = intersect(visiblefeaturesIdx, featureIdx_triangulated);
        
    % 2. Find corresponding matches and adding x-shift
    Ck_z = featureExtracted{k}(1:2, PoseGraphMatrix(visiblefeaturesIdx, k))';
    Ck_z(:,1) = Ck_z(:,1) + countImage*size(I_rgb{1},2);
    
    % 3. Find feature position and reprojection and adding x-shift
    Ck_reproj = zeros(length(visiblefeaturesIdx), 2);
    countFeature = 1;
    for featIdx = visiblefeaturesIdx'
        poseViewed = FeaturesBag(4,featIdx);
        
        Ck_T_W = AbsolutePoses(:,:,k);        
        W_T_Cv = [AbsolutePoses(:,1:3,poseViewed)' -AbsolutePoses(:,1:3,poseViewed)'*AbsolutePoses(:,4,poseViewed);zeros(1,3) 1];        
        Ck_T_Cv = Ck_T_W * W_T_Cv;
        
        Ck_P = Ck_T_Cv*[FeaturesBag(1:3,featIdx);1];
        
        Ck_reproj(countFeature,:) = CameraParams.K(1:2,:)*Ck_P/Ck_P(3);         
        countFeature = countFeature + 1;
    end
    Ck_reproj(:,1) = Ck_reproj(:,1) + countImage*size(I_rgb{1},2);
    
    % 4. Display
    ReprojectionVisualization3_imageNOTshown(Ck_z, Ck_reproj);
    
    % 5. Incremental
    countImage = countImage + 1;
        
end
