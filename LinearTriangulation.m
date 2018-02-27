function Y = LinearTriangulation(matches, P1, P2)

Y = zeros(3, size(matches,1));
for i = 1:size(matches,1)
    A = [skewsymm([matches(i, 1:2) 1]')*P1; skewsymm([matches(i, 6:7) 1]')*P2];
    [~,~,V] = svd(A);
    Y(:,i) = V(1:3,end)/V(4,end);
end

end