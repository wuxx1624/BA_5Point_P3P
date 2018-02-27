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
