function VisualizeEpipolarLines(I1, I2, matches, F)

figure(10); clf;
subplot(1,2,2);
imshow(I2); hold on;
subplot(1,2,1);
imshow(I1); hold on;
cmap = hsv(size(matches,1));

for k = 1:size(matches,1)
    u = [matches((k),1:2) 1]';
    v = [matches((k),6:7) 1]';
    Fu = F*u;
    Fpv = F'*v;
    
    x = [1;size(I1,2)];    
    C2_lu = (-Fu(1).*x - Fu(3))/Fu(2);
    figure(10);
    subplot(1,2,2); hold on;
    plot(matches(k, 6), matches(k, 7), 'co', 'markerfacecolor', 'c', 'markerSize', 5);    hold on;
    plot(x, C2_lu, 'c', 'lineWidth', 2, 'Color',cmap(k,:));    
    
    
    x = [1;size(I1,2)];
    C1_lv = (-Fpv(1).*x - Fpv(3))/Fpv(2);
    figure(10);
    subplot(1,2,1); hold on;
    plot(matches(k, 1), matches(k, 2), ['co'], 'markerfacecolor', 'c' , 'markerSize', 5);    hold on;
    plot(x, C1_lv, 'c', 'lineWidth', 2, 'Color',cmap(k,:));
end