function VisualizeMatches(I1, I2, matches, linecolor, featurecolor)

figure(10); clf ;
imagesc(cat(2, I1, I2)) ;
for k = 1:size(matches,1)
    xa = matches(k,1) ;
    xb = matches(k,6) + size(I1,2);
    
    ya = matches(k,2) ;
    yb = matches(k,7) ;

    hold on ;
    h = line([xa ; xb], [ya ; yb]) ;
    set(h,'linewidth', 1, 'color', linecolor) ;

    plot(matches(k, 1), matches(k, 2), [featurecolor 'o'], 'markerfacecolor',featurecolor, 'markerSize', 3);    hold on;
    plot(matches(k, 6)+size(I1,2), matches(k, 7), [featurecolor 'o'], 'markerfacecolor',featurecolor, 'markerSize', 3);    hold on;
end