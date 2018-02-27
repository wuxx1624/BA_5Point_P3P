function VisualizeAdditionalRegistration(R,t,C1_T_C2, P, inlierIdx)

outlierIdx = setdiff(1:size(P,2), inlierIdx);
figure; hold on;
DrawAxis(zeros(3,1), eye(3), 0.2);
DrawCamera(zeros(3,1), eye(3), 0.2, 'k');
DrawAxis(C1_T_C2(:,4), C1_T_C2(:, 1:3), 0.2);
DrawCamera(C1_T_C2(:,4), C1_T_C2(:, 1:3), 0.2, 'k');
DrawAxis(t, R, 0.4);
DrawCamera(t, R, 0.5 , 'm');
plot3(P(1,outlierIdx), P(2,outlierIdx), P(3,outlierIdx), 'ro', 'markerFacecolor', 'r', 'markerSize', 3);
plot3(P(1,inlierIdx), P(2,inlierIdx), P(3,inlierIdx), 'bo', 'markerFacecolor', 'b', 'markerSize', 3);

end

function DrawCamera(c, R, s, color)
p1 = c + s*R*[sqrt(2);sqrt(2);2]/2;
p2 = c + s*R*[-sqrt(2);sqrt(2);2]/2;
p3 = c + s*R*[-sqrt(2);-sqrt(2);2]/2;
p4 = c + s*R*[sqrt(2);-sqrt(2);2]/2;

hold on
plot3([c(1) p1(1)], [c(2) p1(2)], [c(3) p1(3)], color);
hold on
plot3([c(1) p2(1)], [c(2) p2(2)], [c(3) p2(3)], color);
hold on
plot3([c(1) p3(1)], [c(2) p3(2)], [c(3) p3(3)], color);
hold on
plot3([c(1) p4(1)], [c(2) p4(2)], [c(3) p4(3)], color);
hold on
plot3([p1(1) p2(1)], [p1(2) p2(2)], [p1(3) p2(3)], color);
hold on
plot3([p3(1) p2(1)], [p3(2) p2(2)], [p3(3) p2(3)], color);
hold on
plot3([p3(1) p4(1)], [p3(2) p4(2)], [p3(3) p4(3)], color);
hold on
plot3([p1(1) p4(1)], [p1(2) p4(2)], [p1(3) p4(3)], color);
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
