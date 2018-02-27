function [P, C1_R_C2_best, C1_t_C2_best] = LinearTriangulationE(matches, CameraParams, C1_R_C2,C1_t_C2, visualize)

% [C1_R_C2,C1_t_C2] = CameraPoseCalculation(F, CameraParams.K, visualize);

Y1 = LinearTriangulation(matches, CameraParams.K*[eye(3) zeros(3,1)], CameraParams.K*[C1_R_C2(:,:,1)', -C1_R_C2(:,:,1)'*C1_t_C2]);
ind1 = ChieralityCheck(Y1, eye(3), zeros(3,1), C1_R_C2(:,:,1), C1_t_C2);
ond1 = setdiff(1:size(Y1,2), ind1);

Y2 = LinearTriangulation(matches, CameraParams.K*[eye(3) zeros(3,1)], CameraParams.K*[C1_R_C2(:,:,1)', -C1_R_C2(:,:,1)'*(-C1_t_C2)]);
ind2 = ChieralityCheck(Y2, eye(3), zeros(3,1), C1_R_C2(:,:,1), -C1_t_C2);
ond2 = setdiff(1:size(Y2,2), ind2);

Y3 = LinearTriangulation(matches, CameraParams.K*[eye(3) zeros(3,1)], CameraParams.K*[C1_R_C2(:,:,2)', -C1_R_C2(:,:,2)'*C1_t_C2]);
ind3 = ChieralityCheck(Y3, eye(3), zeros(3,1), C1_R_C2(:,:,2), C1_t_C2);
ond3 = setdiff(1:size(Y3,2), ind3);

Y4 = LinearTriangulation(matches, CameraParams.K*[eye(3) zeros(3,1)], CameraParams.K*[C1_R_C2(:,:,2)', -C1_R_C2(:,:,2)'*(-C1_t_C2)]);
ind4 = ChieralityCheck(Y4, eye(3), zeros(3,1), C1_R_C2(:,:,2), -C1_t_C2);
ond4 = setdiff(1:size(Y4,2), ind4);


if visualize
R(:,:,1) = C1_R_C2(:,:,1);
R(:,:,2) = C1_R_C2(:,:,2);
t = C1_t_C2;

figure;
subplot(2,2,1);
DrawAxis(zeros(3,1), eye(3), 0.5);
DrawCamera(zeros(3,1), eye(3), 0.5);
DrawAxis(t, R(:,:,1), 0.5);
DrawCamera(t, R(:,:,1),0.5);
plot3(Y1(1,ond1), Y1(2,ond1), Y1(3,ond1), 'ro', 'markerFacecolor', 'r', 'markerSize', 3);
plot3(Y1(1,ind1), Y1(2,ind1), Y1(3,ind1), 'bo', 'markerFacecolor', 'b', 'markerSize', 3);
grid on;


subplot(2,2,2);
DrawAxis(zeros(3,1), eye(3), 0.5);
DrawCamera(zeros(3,1), eye(3),0.5);
DrawAxis(-t, R(:,:,1), 0.5);
DrawCamera(-t, R(:,:,1),0.5);
plot3(Y2(1,ond2), Y2(2,ond2), Y2(3,ond2), 'ro', 'markerFacecolor', 'r', 'markerSize', 3);
plot3(Y2(1,ind2), Y2(2,ind2), Y2(3,ind2), 'bo', 'markerFacecolor', 'b', 'markerSize', 3);
grid on;


subplot(2,2,3);
DrawAxis(zeros(3,1), eye(3), 0.5);
DrawCamera(zeros(3,1), eye(3),0.5);
DrawAxis(t, R(:,:,2), 0.5);
DrawCamera(t, R(:,:,2),0.5);
plot3(Y3(1,ond3), Y3(2,ond3), Y3(3,ond3), 'ro', 'markerFacecolor', 'r', 'markerSize', 3);
plot3(Y3(1,ind3), Y3(2,ind3), Y3(3,ind3), 'bo', 'markerFacecolor', 'b', 'markerSize', 3);
grid on;


subplot(2,2,4);
DrawAxis(zeros(3,1), eye(3), 0.5);
DrawCamera(zeros(3,1), eye(3),0.5);
DrawAxis(-t, R(:,:,2), 0.5);
DrawCamera(-t, R(:,:,2),0.5);
plot3(Y4(1,ond4), Y4(2,ond4), Y4(3,ond4), 'ro', 'markerFacecolor', 'r', 'markerSize', 3);
plot3(Y4(1,ind4), Y4(2,ind4), Y4(3,ind4), 'bo', 'markerFacecolor', 'b', 'markerSize', 3);
grid on;
end

Y1 = [Y1(:, ind1); matches(ind1,5)'; matches(ind1,10)']; 
Y2 = [Y2(:, ind2); matches(ind2,5)'; matches(ind2,10)']; 
Y3 = [Y3(:, ind3); matches(ind3,5)'; matches(ind3,10)']; 
Y4 = [Y4(:, ind4); matches(ind4,5)'; matches(ind4,10)'];

% Choose proper relative pose and features position
[~, j] = max([length(ind1), length(ind2), length(ind3), length(ind4)]);
C1_R_C2_best = [];
C1_t_C2_best = [];

if j == 1
    P = Y1;
    C1_R_C2_best =  C1_R_C2(:,:,1); C1_t_C2_best = C1_t_C2;
elseif j == 2
    P = Y2;
    C1_R_C2_best =  C1_R_C2(:,:,1); C1_t_C2_best = -C1_t_C2;
elseif j == 3
    P = Y3;
    C1_R_C2_best =  C1_R_C2(:,:,2); C1_t_C2_best = C1_t_C2;
else
    P = Y4;
    C1_R_C2_best =  C1_R_C2(:,:,2); C1_t_C2_best = -C1_t_C2;
end




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
