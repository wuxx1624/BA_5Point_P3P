function [bestinlierIdx, Rbest, tbest] = ComputeInlierReproj(F, K, matches, threshold)

Kinv = inv(K);
[R,t] = CameraPoseCalculation(F,K,0);
R1 = R(:,:,1); R2 = R(:,:,2);
inlierIdx = cell(1,4);
for i = 1:size(matches,1)
    u = [matches(i,1:2) 1]';
    v = [matches(i,5:6) 1]';
    Fu = F*u;
    Fpv = F'*v;
    
    % Checking 4 cases:
    % Case 1 and 2
    b = t; A = [Kinv*u -R1*Kinv*v];
    d = A\b;
    if d(1) >= 0 && d(2) >= 0
        vproj = K*(R1*d(2)*Kinv*v + t); vproj = vproj/vproj(3);
        uproj = K*R1'*(d(1)*Kinv*u - t); uproj = uproj/uproj(3);
        reproj_error = abs(vproj'*Fu)/norm(Fu(1:2)) + abs(uproj'*Fpv)/norm(Fpv(1:2));
        if reproj_error < threshold
            inlierIdx{1} = [inlierIdx{1} i];
            continue;
        end
    elseif d(1) < 0 && d(2) < 0
        t1 = -t; d = -d;
        vproj = K*(R1*d(2)*Kinv*v + t1); vproj = vproj/vproj(3);
        uproj = K*R1'*(d(1)*Kinv*u - t1); uproj = uproj/uproj(3);
        reproj_error = abs(vproj'*Fu)/norm(Fu(1:2)) + abs(uproj'*Fpv)/norm(Fpv(1:2));
        if reproj_error < threshold
            inlierIdx{2} = [inlierIdx{2} i];
            continue;
        end
    end
    
    % Case 3 and 4
    b = t; A = [Kinv*u -R2*Kinv*v];
    d = A\b;
    if d(1) >= 0 && d(2) >= 0
        vproj = K*(R2*d(2)*Kinv*v + t); vproj = vproj/vproj(3);
        uproj = K*R2'*(d(1)*Kinv*u - t); uproj = uproj/uproj(3);
        reproj_error = abs(vproj'*Fu)/norm(Fu(1:2)) + abs(uproj'*Fpv)/norm(Fpv(1:2));
        if reproj_error < threshold
            inlierIdx{3} = [inlierIdx{3} i];
            continue;
        end
    elseif d(1) < 0 && d(2) < 0
        t1 = -t; d = -d;
        vproj = K*(R2*d(2)*Kinv*v + t1); vproj = vproj/vproj(3);
        uproj = K*R2'*(d(1)*Kinv*u - t1); uproj = uproj/uproj(3);
        reproj_error = abs(vproj'*Fu)/norm(Fu(1:2)) + abs(uproj'*Fpv)/norm(Fpv(1:2));
        if reproj_error < threshold
            inlierIdx{4} = [inlierIdx{4} i];
            continue;
        end
    end
    
end

inlier_vect = [length(inlierIdx{1}),length(inlierIdx{2}),length(inlierIdx{3}),length(inlierIdx{4})]
[~,j] = max(inlier_vect);
bestinlierIdx = inlierIdx{j};
if j == 1
    Rbest = R1; tbest = t;
elseif j == 2
    Rbest = R1; tbest = -t;
elseif j == 3
    Rbest = R2; tbest = t;
else
    Rbest = R2; tbest = -t;
end


