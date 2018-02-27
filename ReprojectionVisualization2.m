function ReprojectionVisualization2(I1,featureCoords,featureReproj)

figure; clf ;
imshow(I1); hold on;
plot(featureReproj(1, 1), featureReproj(1, 2), 'bo', 'markerSize', 5);   
plot(featureCoords(1, 1), featureCoords(1, 2), 'r*', 'markerSize', 3);
legend('reproj', 'measured');
for k = 2:size(featureReproj,1)
    plot(featureReproj(k, 1), featureReproj(k, 2), 'bo', 'markerSize', 5);   
end

for k = 2:size(featureCoords,1)
    plot(featureCoords(k, 1), featureCoords(k, 2), 'r*', 'markerSize', 3);
end