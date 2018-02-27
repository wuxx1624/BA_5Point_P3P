function [FeaturesBag] = NonlinearTriangulation(FeaturesBag,           ...
                                AbsolutePoses,         ...
                                PoseGraphMatrix,       ...
                                newPose, usedPoses,    ...
                                featureExtracted,      ...                                         
                                CameraParams,visualization,I_rgb)

updateFeatIdx = find(FeaturesBag(4,:) == newPose);

for k = updateFeatIdx    
    posesViewkFeat = find(PoseGraphMatrix(k, usedPoses));
    Cn_P = FeaturesBag(1:3,k);
    
    norm_iter = zeros(3,1);
    
    if visualization == 1
        figure(1);
        subplot(1,2,1); imshow(I_rgb{newPose}); hold on;
        % reprojection_matches = ReprojectionToMatches(Cn_P, eye(3), zeros(3,1), CameraParams);
        % ReprojectionVisualization2(I_rgb{newPose}, Cn_z', reprojection_matches); 

        subplot(1,2,2); imshow(I_rgb{usedPoses(1)}); hold on;
        % reprojection_matches = ReprojectionToMatches(Cl_P, eye(3), zeros(3,1), CameraParams);
        % ReprojectionVisualization2(I_rgb{usedPoses(1)}, featureExtracted{l}(1:2, PoseGraphMatrix(k,l))', reprojection_matches); 
    end
    
    % Solving iteratively
    for iter = 1:3
        % Building simple Jacobian and res
        Jf = zeros(2*length(posesViewkFeat)+2 , 3);
        resf = zeros(2*length(posesViewkFeat)+2 , 1);

        Jf(1:2,:) = 1/Cn_P(3)*[1 0 -Cn_P(1)/Cn_P(3);0 1 -Cn_P(2)/Cn_P(3)];
        Cn_z = featureExtracted{newPose}(1:2, PoseGraphMatrix(k,newPose));
        resf(1:2) = CameraParams.Kinv(1:2,:)*[Cn_z;1] - [Cn_P(1)/Cn_P(3);Cn_P(2)/Cn_P(3)];
        count = 1;
        for l = usedPoses(posesViewkFeat)
            poseViewed = FeaturesBag(4,k);
            Cl_T_W = AbsolutePoses(:,:,l);        
            W_T_Cv = [AbsolutePoses(:,1:3,poseViewed)' -AbsolutePoses(:,1:3,poseViewed)'*AbsolutePoses(:,4,poseViewed);zeros(1,3) 1];        
            Cl_T_Cv = Cl_T_W * W_T_Cv;
            Cl_P = Cl_T_Cv * [Cn_P;1];
            Jf(2*count+1:2*count+2, :) = 1/Cl_P(3)*[1 0 -Cl_P(1)/Cl_P(3);0 1 -Cl_P(2)/Cl_P(3)] * Cl_T_Cv(:, 1:3);
            resf(2*count+1:2*count+2) = ...
            CameraParams.Kinv(1:2,:)*[featureExtracted{l}(1:2, PoseGraphMatrix(k,l));1] - ...
            [Cl_P(1)/Cl_P(3);Cl_P(2)/Cl_P(3)];
        end

        % Each iteration visualization
        if visualization == 1
            figure(1);
            subplot(1,2,1); hold on;
            reprojection_matches = ReprojectionToMatches(Cn_P, eye(3), zeros(3,1), CameraParams);
            ReprojectionVisualization3_imageNOTshown(Cn_z', reprojection_matches); 
            
            subplot(1,2,2); hold on;
            reprojection_matches = ReprojectionToMatches(Cl_P, eye(3), zeros(3,1), CameraParams);
            ReprojectionVisualization3_imageNOTshown(featureExtracted{l}(1:2, PoseGraphMatrix(k,l))', reprojection_matches); 
        end
        
        
        % Solve and update for each feature
        df = (Jf'*Jf) \ Jf'*resf;
        Cn_P = Cn_P + df;
        norm_iter(iter) = norm(resf);
        if(norm_iter(iter) < 1e-4)
            break;
        end
    end        
    
    if visualization == 1
        pause;
        close all;
    end
    % update to feature bag
    FeaturesBag(1:3,k) = Cn_P;
end


