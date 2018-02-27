function [Cr_P, RTref, RTscale] = LinearTriangulationRt(RelativePoses, PoseGraphMatrix, refPose, newPose, fr, fn, CameraParams, visualize)


if refPose < newPose
    C1_R_C2 = RelativePoses(:, 1:3);
    C1_t_C2 = RelativePoses(:, 4);    
else
    C1_R_C2 = RelativePoses(:, 1:3)';
    C1_t_C2 = -RelativePoses(:, 1:3)'*RelativePoses(:, 4);
end
RTref = [eye(3), zeros(3,1)];
RTscale = [C1_R_C2', -C1_R_C2'*C1_t_C2];


matches = GetMatchFromPoseGraph(PoseGraphMatrix, (1:size(PoseGraphMatrix))', refPose, newPose, fr, fn);

Y = LinearTriangulation(matches, CameraParams.K*[eye(3) zeros(3,1)], CameraParams.K*[C1_R_C2(:,:,1)', -C1_R_C2(:,:,1)'*C1_t_C2]);
ind = ChieralityCheck(Y, eye(3), zeros(3,1), C1_R_C2(:,:,1), C1_t_C2);
ond = setdiff(1:size(Y,2), ind);


if visualize
R(:,:,1) = C1_R_C2(:,:,1);
t = C1_t_C2;

figure;
DrawAxis(zeros(3,1), eye(3), 0.2);
DrawCamera(zeros(3,1), eye(3), 0.2);
DrawAxis(t, R(:,:,1), 0.2);
DrawCamera(t, R(:,:,1),0.2);
plot3(Y(1,ond), Y(2,ond), Y(3,ond), 'ro', 'markerFacecolor', 'r', 'markerSize', 3);
plot3(Y(1,ind), Y(2,ind), Y(3,ind), 'bo', 'markerFacecolor', 'b', 'markerSize', 3);
grid on;
end

% Cr_P = [Y(:, ind); matches(ind,5)'; matches(ind,10)']; 
% [~, update_indices] = ismember(matches(ind,5)', PoseGraphMatrix(PoseGraphMatrix(:,refPose)~=0,refPose));
Cr_P = zeros(4, size(PoseGraphMatrix,1));
[~, update_indices] = ismember(matches(ind,5)', PoseGraphMatrix(:,refPose));
Cr_P(:,update_indices) = [Y(:, ind); refPose*ones(1,length(ind))];

% Cr_P_full = zeros(4, sum(PoseGraphMatrix(:,refPose)~=0));
% Cr_P_full(:,update_indices) = Cr_P;

% REQUIREMENT: indicies of nonzero elements in Posegraph is the same as
% Cr_P_full columns
% idx = find(Cr_P_full(4,:)~=0);
% sm_vec = PoseGraphMatrix(PoseGraphMatrix(:,refPose)~=0,refPose);
% norm(double(sm_vec(idx)) - double(Cr_P_full(4,idx)'));

end

function DrawCamera(c, R, s)
p1 = c + s*R*[sqrt(2);sqrt(2);2]/2;
p2 = c + s*R*[-sqrt(2);sqrt(2);2]/2;
p3 = c + s*R*[-sqrt(2);-sqrt(2);2]/2;
p4 = c + s*R*[sqrt(2);-sqrt(2);2]/2;

hold on
plot3([c(1) p1(1)], [c(2) p1(2)], [c(3) p1(3)], 'k');
hold on
plot3([c(1) p2(1)], [c(2) p2(2)], [c(3) p2(3)], 'k');
hold on
plot3([c(1) p3(1)], [c(2) p3(2)], [c(3) p3(3)], 'k');
hold on
plot3([c(1) p4(1)], [c(2) p4(2)], [c(3) p4(3)], 'k');
hold on
plot3([p1(1) p2(1)], [p1(2) p2(2)], [p1(3) p2(3)], 'k');
hold on
plot3([p3(1) p2(1)], [p3(2) p2(2)], [p3(3) p2(3)], 'k');
hold on
plot3([p3(1) p4(1)], [p3(2) p4(2)], [p3(3) p4(3)], 'k');
hold on
plot3([p1(1) p4(1)], [p1(2) p4(2)], [p1(3) p4(3)], 'k');
hold on;
end

function DrawAxis(c, R, s)
rx = s*R(:,1); ry = s*R(:,2); rz = s*R(:,3);
hold on
plot3([c(1) c(1)+rx(1)], [c(2) c(2)+rx(2)], [c(3) c(3)+rx(3)], 'r-', 'LineWidth', 2);
hold on
plot3([c(1) c(1)+ry(1)], [c(2) c(2)+ry(2)], [c(3) c(3)+ry(3)], 'g-', 'LineWidth', 2);
hold on
plot3([c(1) c(1)+rz(1)], [c(2) c(2)+rz(2)], [c(3) c(3)+rz(3)], 'b-', 'LineWidth', 2);
axis equal
end
