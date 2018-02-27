function ReprojectionVisualization(I1,I2,match,rmatch)

figure; clf ;
imagesc(cat(2, I1, I2)) ; hold on;
for k = 1:size(rmatch,1)
    plot(rmatch(k, 1), rmatch(k, 2), 'r*', 'markerSize', 3);  
    plot(match(k, 1), match(k, 2), 'bo', 'markerSize', 5);
    plot(match(k, 6)+size(I1,2), match(k, 7), 'bo', 'markerSize', 5);
    plot(rmatch(k, 3)+size(I1,2), rmatch(k, 4), 'r*', 'markerSize', 3);
end

legend('reprojection', 'original');