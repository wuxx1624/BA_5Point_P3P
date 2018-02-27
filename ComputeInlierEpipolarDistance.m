function [bestinlierIdx] = ComputeInlierEpipolarDistance(F, matches, threshold)

bestinlierIdx = [];
for i = 1:size(matches,1)
    u = [matches(i,1:2) 1]';
    v = [matches(i,6:7) 1]';
    Fu = F*u;
    Fpv = F'*v;
    reproj_error = abs(v'*Fu)*(1/norm(Fu(1:2)) + 1/norm(Fpv(1:2)));
    if reproj_error < threshold
        bestinlierIdx = [bestinlierIdx i];        
    end    
end
