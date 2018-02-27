function VisualizeMatches3(I3, I1, I2, matches, linecolor1, linecolor2, featurecolor)

figure; clf ;
imagesc(cat(2, I1, I3, I2)) ;
for k = 1:size(matches,1)
    xa = matches(k,6) ;
    xb = matches(k,1) + size(I1,2);
    
    ya = matches(k,7) ;
    yb = matches(k,2) ;

    hold on ;
    h = line([xa ; xb], [ya ; yb]) ;
    set(h,'linewidth', 1, 'color', linecolor1) ;
    
    
    xa = matches(k,1) + size(I1,2);
    xb = matches(k,11) + size(I1,2) + size(I3,2);
    
    ya = matches(k,2) ;
    yb = matches(k,12) ;

    hold on ;
    h = line([xa ; xb], [ya ; yb]) ;
    set(h,'linewidth', 1, 'color', linecolor2) ;

    plot(matches(k, 6), matches(k, 7), [featurecolor 'o'], 'markerfacecolor',featurecolor, 'markerSize', 3);    hold on;
    plot(matches(k, 1)+size(I1,2), matches(k, 2), [featurecolor 'o'], 'markerfacecolor',featurecolor, 'markerSize', 3);    hold on;
    plot(matches(k, 11)+size(I1,2)+size(I3,2), matches(k, 12), [featurecolor 'o'], 'markerfacecolor',featurecolor, 'markerSize', 3);    hold on;
end