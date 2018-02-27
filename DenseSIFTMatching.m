% Dense SIFT Matching

function [im1, im2] = DenseSIFTMatching(I1, I2, matches, C1_T_C2, CameraParams)

% rotating C'_H_C C_t = [1;0;0] => C'_R_C = R(txe3, tt) 
axis1 = skewsymm([1;0;0])*C1_T_C2(:,4);
tt1 = -atan2(norm(axis1), [1;0;0]'*C1_T_C2(:,4));
axis1 = axis1/norm(axis1);
R1 = aa2rot(axis1, tt1);
H1 = CameraParams.K*R1*CameraParams.Kinv;
    
tform1 = projective2d(H1');
ov = imref2d([size(I1) 3]);
im1 = imwarp(I1, tform1, 'OutputView', ov);

% image warping 2
axis2 = skewsymm([-1;0;0])*(-C1_T_C2(:,1:3)'*C1_T_C2(:,4));
tt2 = -atan2(norm(axis2), [-1;0;0]'*(-C1_T_C2(:,1:3)'*C1_T_C2(:,4)));
axis2 = axis2/norm(axis2);
R2 = aa2rot(axis2, tt2);
H2 = CameraParams.K*R2*CameraParams.Kinv;

tform2 = projective2d(H2');
im2 = imwarp(I2, tform2, 'OutputView', ov);

% visualization 1
figure; 
subplot(2,1,1);
imshow([I1 I2]);
subplot(2,1,2);
imshow([im1 im2]);

% visualization 2
transform_matches_im1 = [matches(:, 1:2) ones(size(matches,1), 1)];
transform_matches_im2 = [matches(:, 6:7) ones(size(matches,1), 1)];

transform_matches_im1 = H1*transform_matches_im1';
transform_matches_im1(1,:) = transform_matches_im1(1,:)./transform_matches_im1(3,:);
transform_matches_im1(2,:) = transform_matches_im1(2,:)./transform_matches_im1(3,:);

transform_matches_im2 = H2*transform_matches_im2';
transform_matches_im2(1,:) = transform_matches_im2(1,:)./transform_matches_im2(3,:);
transform_matches_im2(2,:) = transform_matches_im2(2,:)./transform_matches_im2(3,:);

transform_matches = [transform_matches_im1(1:2,:)' matches(:, 3:5) transform_matches_im2(1:2,:)' matches(:, 8:10)];
VisualizeMatches(im1, im2, transform_matches, 'b', 'g');
norm(transform_matches_im2(2,:) - transform_matches_im1(2,:))



% visualization 3: epipolar line
cmap = hsv(size(matches,1));
F = CameraParams.Kinv'* C1_T_C2(:,1:3)'* skewsymm(C1_T_C2(:,4)) *CameraParams.Kinv;

figure(5);
subplot(1,2,1); imshow(im1);
subplot(1,2,2); imshow(im2);
for k = 1:length(matches)
    u = [matches(k,1:2) 1]';
    v = [matches(k,6:7) 1]';
    Fu = inv(H2)'*F*u;
    Fpv = inv(H1)'*F'*v;
    
    x = [1;size(im2,2)];    
    C2_lu = (-Fu(1).*x - Fu(3))/Fu(2);
    figure(5);
    subplot(1,2,2); hold on;
    plot(transform_matches(k,6), transform_matches(k,7), 'go', 'markerfacecolor', 'g', 'markerSize', 5);
    plot(x, C2_lu, 'c', 'lineWidth', 1, 'Color',cmap(k,:));
    
    x = [1;size(im1,2)];
    C1_lv = (-Fpv(1).*x - Fpv(3))/Fpv(2);
    figure(5);
    subplot(1,2,1); hold on;
    plot(transform_matches(k,1), transform_matches(k,2), 'go', 'markerfacecolor', 'g', 'markerSize', 5);
    plot(x, C1_lv, 'c', 'lineWidth', 1, 'Color',cmap(k,:));
end

% %% Stereo matching refinement
% min_norm = 10000;
% best_fixed_angle = [];
% for ii = 0.35:0.01:0.45
%     for jj = 0.2:0.01:0.3
%         R1 = aa2rot(axis1, tt1 + ii);
%         H1 = CameraParams.K*R1*CameraParams.Kinv;
% 
%         tform1 = projective2d(H1');
%         ov = imref2d([size(I1) 3]);
%         % im1 = imwarp(I1, tform1, 'OutputView', ov);
% 
%         R2 = aa2rot(axis2, tt2 + jj);
%         H2 = CameraParams.K*R2*CameraParams.Kinv;
% 
%         tform2 = projective2d(H2');
%         % im2 = imwarp(I2, tform2, 'OutputView', ov);
% 
%         transform_matches_im1 = [matches(:, 1:2) ones(size(matches,1), 1)];
%         transform_matches_im2 = [matches(:, 6:7) ones(size(matches,1), 1)];
% 
%         transform_matches_im1 = H1*transform_matches_im1';
%         transform_matches_im1(1,:) = transform_matches_im1(1,:)./transform_matches_im1(3,:);
%         transform_matches_im1(2,:) = transform_matches_im1(2,:)./transform_matches_im1(3,:);
% 
%         transform_matches_im2 = H2*transform_matches_im2';
%         transform_matches_im2(1,:) = transform_matches_im2(1,:)./transform_matches_im2(3,:);
%         transform_matches_im2(2,:) = transform_matches_im2(2,:)./transform_matches_im2(3,:);
% 
%         transform_matches = [transform_matches_im1(1:2,:)' matches(:, 3:5) transform_matches_im2(1:2,:)' matches(:, 8:10)];
%         % VisualizeMatches(im1, im2, transform_matches, 'b', 'g');
% 
%         error = norm(transform_matches_im1(2,:) - transform_matches_im2(2,:))
%         if (error < min_norm)
%             min_norm = error;
%             best_fixed_angle = [ii jj];
%         end
%     end
% end
% 
% 
% R1 = aa2rot(axis1, tt1 + best_fixed_angle(1));
% H1 = CameraParams.K*R1*CameraParams.Kinv;
% 
% tform1 = projective2d(H1');
% ov = imref2d([size(I1) 3]);
% im1 = imwarp(I1, tform1, 'OutputView', ov);
% 
% R2 = aa2rot(axis2, tt2 + best_fixed_angle(2));
% H2 = CameraParams.K*R2*CameraParams.Kinv;
% 
% tform2 = projective2d(H2');
% im2 = imwarp(I2, tform2, 'OutputView', ov);
% VisualizeMatches(im1, im2, transform_matches, 'b', 'g');

%[H1, H2] = StereoMatchingRefinement(H1, H2, matches, I1, I2);

% dense sift extraction

cellsize=3;
gridspacing=1;

descriptor1 = mexDenseSIFT(im2double(im1),cellsize,gridspacing);
descriptor2 = mexDenseSIFT(im2double(im2),cellsize,gridspacing);

% load('./Set4/descriptor1.mat');
% load('./Set4/descriptor2.mat');


descriptor1 = permute(descriptor1, [3 2 1]);
descriptor2 = permute(descriptor2, [3 2 1]);
diff_matrix = zeros(128, size(im2,2));


disparity_map = zeros(size(im1,1), size(im1,2));
F = CameraParams.Kinv'*C1_T_C2(:,1:3)'*skewsymm(C1_T_C2(:,4))*CameraParams.Kinv;
for i = 10:size(disparity_map,1)-10
    
    fprintf('Progess i = %d\n', i);
    
    for j = 10:size(disparity_map,2)-200                
        
        u = [j;i;1];
        Fu = F*u;
        
        %x = 1:size(im2,2);
        x = 1000;
        C2_lu = (-Fu(1).*x - Fu(3))/Fu(2);
        C2_lu = ceil(C2_lu);
        
        if C2_lu > size(disparity_map,1) - 3
            continue;
        end
        
        min_error = 1e12;
        best_match = [];
        for w = -3:3
            diff_matrix = int16(descriptor2(:,:,C2_lu+w)) - int16(repmat(descriptor1(:,j,i), [1 size(im1,2)]));
            diff_matrix = diff_matrix .* diff_matrix;
            disparity_vect = sum(diff_matrix,1);
            [min_value, idx_min] = min(disparity_vect);
            if min_value < min_error
                min_error = min_value;
                best_match = [C2_lu+w idx_min];
            end
        end        
        disparity_map(i,j) = abs(best_match(2) - j);
    end
end

save('./Set4/disparity_map', 'disparity_map');

figure;
subplot(1,2,1);
image(disparity_map,'CDataMapping','scaled');
colorbar;
subplot(1,2,2);
imshow(im1);

end