function [bestinlierIdx, Rbest, tbest] = ComputeInlierReprojUltimateTest(I1, I2, idx_selected, F, K, matches, threshold)

Kinv = inv(K);
[R,t] = CameraPoseCalculation(F,K,0);
R1 = R(:,:,1); R2 = R(:,:,2);
inlierIdx = cell(1,4);

figure(3) ;
imagesc(cat(2, I1, I2)) ; hold on;   

for i = 1:size(matches,1)
        
    % Necessary variables
    isthismatchinlier = 0;
    u = [matches(i,1:2) 1]';
    v = [matches(i,5:6) 1]';
    Fu = F*u;
    Fpv = F'*v;
    
    
    % Draw matches, epipole, epipolar line
    figure(4) ; clf ;
    imagesc(cat(2, I1, I2)) ;
    xa = matches(i,1) ;
    xb = matches(i,5) + size(I1,2) ;
    ya = matches(i,2) ;
    yb = matches(i,6) ;

    hold on ;
    h = line([xa ; xb], [ya ; yb]) ;
    set(h,'linewidth', 1, 'color', 'b') ;
    
    plot(matches(i, 1), matches(i, 2), 'go', 'markerfacecolor','g');    
    plot(matches(i, 5) + size(I1,2), matches(i, 6), 'go', 'markerfacecolor','g');
    
    C1_pix_C2 = null(F); C1_pix_C2 = C1_pix_C2/C1_pix_C2(3);
    C2_pix_C1 = null(F'); C2_pix_C1 = C2_pix_C1/C2_pix_C1(3);

    plot(C1_pix_C2(1), C1_pix_C2(2), 'yo', 'markerfacecolor','y');    
    plot(C2_pix_C1(1) + size(I1,2), C2_pix_C1(2), 'yo', 'markerfacecolor','y');
    
    x = [1;size(I1,2)];
    C2_lu = (-Fu(1).*x - Fu(3))/Fu(2);
    plot(x + size(I1,2), C2_lu, 'c', 'lineWidth', 2);
    
    x = [1;size(I1,2)];
    C1_lv = (-Fpv(1).*x - Fpv(3))/Fpv(2);
    plot(x, C1_lv, 'c', 'lineWidth', 2);

    
    % Checking 4 cases:
    % Case 1 and 2
    b = t; A = [Kinv*u -R1*Kinv*v];
    d = A\b;
    if d(1) >= 0 && d(2) >= 0
        fprintf('test case 1: R1 and t \n');
        fprintf('d = \n'); d
        C1_P_v = (R1*d(2)*Kinv*v + t);
        C2_P_u = R1'*(d(1)*Kinv*u - t);
        vproj = K*C1_P_v; vproj = vproj/vproj(3);
        uproj = K*C2_P_u; uproj = uproj/uproj(3);
        fprintf('error: \n');
        reproj_error = abs(v'*Fu)/norm(Fu(1:2)) + abs(u'*Fpv)/norm(Fpv(1:2))
        % reproj_error = norm(uproj - v) + norm(vproj - u)
        if reproj_error < threshold
            inlierIdx{1} = [inlierIdx{1} i];
            isthismatchinlier = 1;
        end
    elseif d(1) < 0 && d(2) < 0
        t1 = -t; d = -d;
        fprintf('test case 2: R1 and -t \n');
        fprintf('d = \n'); d
        C1_P_v = (R1*d(2)*Kinv*v + t1);
        C2_P_u = R1'*(d(1)*Kinv*u - t1);
        vproj = K*C1_P_v; vproj = vproj/vproj(3);
        uproj = K*C2_P_u; uproj = uproj/uproj(3);
        fprintf('error: \n');
        reproj_error = abs(v'*Fu)/norm(Fu(1:2)) + abs(u'*Fpv)/norm(Fpv(1:2))
        % reproj_error = norm(uproj - v) + norm(vproj - u)
        if reproj_error < threshold
            inlierIdx{2} = [inlierIdx{2} i];
            isthismatchinlier = 1;
        end
    end
    
    % Case 3 and 4
    b = t; A = [Kinv*u -R2*Kinv*v];
    d = A\b;
    if d(1) >= 0 && d(2) >= 0
        fprintf('test case 3: R2 and t \n');
        fprintf('d = \n'); d
        C1_P_v = (R2*d(2)*Kinv*v + t);
        C2_P_u = R2'*(d(1)*Kinv*u - t);
        vproj = K*C1_P_v; vproj = vproj/vproj(3);
        uproj = K*C2_P_u; uproj = uproj/uproj(3);
        fprintf('error: \n');
        reproj_error = abs(v'*Fu)/norm(Fu(1:2)) + abs(u'*Fpv)/norm(Fpv(1:2))
        % reproj_error = norm(uproj - v) + norm(vproj - u)
        if reproj_error < threshold
            inlierIdx{3} = [inlierIdx{3} i];
            isthismatchinlier = 1;
        end
    elseif d(1) < 0 && d(2) < 0
        t1 = -t; d = -d;
        fprintf('test case 4: R2 and t \n');
        fprintf('d = \n'); d
        C1_P_v = (R2*d(2)*Kinv*v + t1);
        C2_P_u = R2'*(d(1)*Kinv*u - t1);
        vproj = K*C1_P_v; vproj = vproj/vproj(3);
        uproj = K*C2_P_u; uproj = uproj/uproj(3);
        fprintf('error: \n');
        reproj_error = abs(v'*Fu)/norm(Fu(1:2)) + abs(u'*Fpv)/norm(Fpv(1:2))
        % reproj_error = norm(uproj - v) + norm(vproj - u)
        if reproj_error < threshold
            inlierIdx{4} = [inlierIdx{4} i];
            isthismatchinlier = 1;
        end
    end
    
    
    if( isthismatchinlier == 1)
        % projected point and check
        figure(4); hold on;    
        plot(vproj(1), vproj(2), 'ro', 'markerface', 'r');
        plot(uproj(1) + size(I1,2), uproj(2), 'ro', 'markerface', 'r');    
        % pause(0.2);
        
        
        figure(3); hold on;
        xa = matches(i,1) ;
        xb = matches(i,5) + size(I1,2) ;
        ya = matches(i,2) ;
        yb = matches(i,6) ;

        hold on ;
        h = line([xa ; xb], [ya ; yb]) ;
        set(h,'linewidth', 1, 'color', 'b') ;

        plot(matches(i, 1), matches(i, 2), 'go', 'markerfacecolor','g');    
        plot(matches(i, 5) + size(I1,2), matches(i, 6), 'go', 'markerfacecolor','g'); 
        inlier_vect = [length(inlierIdx{1}),length(inlierIdx{2}),length(inlierIdx{3}),length(inlierIdx{4})]    
    else
        fprintf('-------------');
    end    
end

figure(3); hold on;
for i = idx_selected
    xa = matches(i,1) ;
    xb = matches(i,5) + size(I1,2) ;
    ya = matches(i,2) ;
    yb = matches(i,6) ;

    hold on ;
    h = line([xa ; xb], [ya ; yb]) ;
    set(h,'linewidth', 1, 'color', 'm') ; hold on;
    
    plot(matches(i, 1), matches(i, 2), 'mo', 'markerfacecolor','m');    hold on;
    plot(matches(i, 5) + size(I1,2), matches(i, 6), 'mo', 'markerfacecolor','m'); hold on;
end
axis image off ; 

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

