function VisualizeMultiPoses(AbsolutePoses, FeaturesBag, usedPoses, newPose)

figure; hold on;
for l = usedPoses
    DrawAxis(-AbsolutePoses(:,1:3,l)'*AbsolutePoses(:,4,l), AbsolutePoses(:,1:3,l)', 0.5);
    DrawCamera(-AbsolutePoses(:,1:3,l)'*AbsolutePoses(:,4,l), AbsolutePoses(:,1:3,l)', 0.5, 'k');    
    
    P = FeaturesBag(1:3, FeaturesBag(4,:) == l);
    % Compute w.r.t the world coordinate
    P = AbsolutePoses(:,1:3,l)'*(P - repmat(AbsolutePoses(:,4,l), [1 size(P,2)]));
    plot3(P(1,:), P(2,:), P(3,:), 'bo', 'markerFacecolor', 'b', 'markerSize', 3);    
end
DrawAxis(-AbsolutePoses(:,1:3,newPose)'*AbsolutePoses(:,4,newPose), AbsolutePoses(:,1:3,newPose)', 0.5);
DrawCamera(-AbsolutePoses(:,1:3,newPose)'*AbsolutePoses(:,4,newPose), AbsolutePoses(:,1:3,newPose)', 0.5 , 'm');
P = FeaturesBag(1:3, FeaturesBag(4,:) == newPose);
% Compute w.r.t the world coordinate
P = AbsolutePoses(:,1:3,newPose)'*(P - repmat(AbsolutePoses(:,4,newPose), [1 size(P,2)]));
plot3(P(1,:), P(2,:), P(3,:), 'bo', 'markerFacecolor', 'b', 'markerSize', 3);    
axis([-10 10 -10 10 -5 30]);
grid on;

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
