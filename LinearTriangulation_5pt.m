function P = LinearTriangulation_5pt(matches, CameraParams, C1_R_C2, C1_t_C2, visualize)

Y1 = LinearTriangulation(matches, CameraParams.K*[eye(3) zeros(3,1)], CameraParams.K*[C1_R_C2(:,:,1)', -C1_R_C2(:,:,1)'*C1_t_C2]);
ind1 = ChieralityCheck(Y1, eye(3), zeros(3,1), C1_R_C2(:,:,1), C1_t_C2);
ond1 = setdiff(1:size(Y1,2), ind1);


P = [Y1(:, ind1); matches(ind1,5)'; matches(ind1,10)']; 

if visualize
R = C1_R_C2;
t = C1_t_C2;

figure;
DrawAxis(zeros(3,1), eye(3), 0.5);
DrawCamera(zeros(3,1), eye(3), 0.5);
DrawAxis(t, R, 0.5);
DrawCamera(t, R,0.5);
plot3(Y1(1,ond1), Y1(2,ond1), Y1(3,ond1), 'ro', 'markerFacecolor', 'r', 'markerSize', 3);
plot3(Y1(1,ind1), Y1(2,ind1), Y1(3,ind1), 'bo', 'markerFacecolor', 'b', 'markerSize', 3);
grid on;

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

end