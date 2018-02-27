function [R,t] = CameraPoseCalculation(F, K, visualize)

% visualize = 0;


E = K'*F*K;

[U, ~, V] = svd(E);

R = zeros(3,3,2);
W = [0 -1 0;1 0 0;0 0 1];
R(:,:,1) = V*W*U'; R(:,:,1) = R(:,:,1)/det(R(:,:,1));
R(:,:,2) = V*W'*U'; R(:,:,2) = R(:,:,2)/det(R(:,:,2));

t = V(:,3);


if visualize == 1
figure;
subplot(2,2,1);
DrawAxis(zeros(3,1), eye(3), 0.2);
DrawCamera(zeros(3,1), eye(3),0.2);
DrawAxis(t, R(:,:,1), 0.2);
DrawCamera(t, R(:,:,1),0.2);
grid on;


subplot(2,2,2);
DrawAxis(zeros(3,1), eye(3), 0.2);
DrawCamera(zeros(3,1), eye(3),0.2);
DrawAxis(-t, R(:,:,1), 0.2);
DrawCamera(-t, R(:,:,1),0.2);
grid on;


subplot(2,2,3);
DrawAxis(zeros(3,1), eye(3), 0.2);
DrawCamera(zeros(3,1), eye(3),0.2);
DrawAxis(t, R(:,:,2), 0.2);
DrawCamera(t, R(:,:,2),0.2);
grid on;


subplot(2,2,4);
DrawAxis(zeros(3,1), eye(3), 0.2);
DrawCamera(zeros(3,1), eye(3),0.2);
DrawAxis(-t, R(:,:,2), 0.2);
DrawCamera(-t, R(:,:,2),0.2);
grid on;
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