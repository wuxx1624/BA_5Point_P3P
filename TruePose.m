function [R_true,t_true,idx]=TruePose(R,t,F1,F2,both_see_feat,thresh,K)

%%%%%%%% input %%%%%%%%%%%%%%%%%%%%%%
% F1&F2: 2D features, 3X1 vector
% K: intrinsic matrix

incorrect = 0;
best_inliners = -1;
for i=1:size(R,3)
    b1 = F1;
    b2 = F2;
%     b1 = b1/norm(b1);
%     b2 = b2/norm(b2);
    b1 = b1./b1(3);
    b2 = b2./b2(3);
%     A = [R(:,:,i)*b1 -b2 t(:,:,i)];
%     [U,S,V] = svd(A);
%     d=V(:,end)/V(end,end);
    A = [K*R(:,:,i)*b1 -b2];
    d = A\(-K*t(:,:,i));
    if d(1)>0 && d(2)>0
        [num_inliners,error,inlierIdx] = count_inliner(R(:,:,i),t(:,:,i),both_see_feat,thresh,K,0,0);
        if best_inliners < num_inliners
            R_true = R(:,:,i);
            t_true = t(:,:,i);
            idx = i;
            best_inliners = num_inliners;
        end
    else
        incorrect = incorrect + 1;
    end
end

if incorrect == 4
    R_true = zeros(3);
    t_true = zeros(3,1);
    idx = 0;
end
