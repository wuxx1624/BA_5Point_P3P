function [num_inliners,error,bestInliners] = count_inliner(R,t,both_see_feat,thresh,K,idx1,idx2)

num_inliners = 0;
error = zeros(size(both_see_feat,1),1);
bestInliners = [];
for i = 1:size(both_see_feat,1)
%     b1 = featureExtracted_true{idx1}(:,both_see_feat(i,1));
%     b2 = featureExtracted_true{idx2}(:,both_see_feat(i,2));
    b1 = both_see_feat(i,1:2);
    b2 = both_see_feat(i,6:7);
    b1 = [b1';ones(1,size(b1,1))];
    b2 = [b2';ones(1,size(b2,1))];
    b1_ = b1./b1(3,:);
    b2_ = b2./b2(3,:);
    A = [K*R*b1_ -b2_ K*t];
    [U,S,V] = svd(A);
    d=V(:,end)/V(end,end);
%     A = [K*R*b1_ -b2_];
%     d = A\(-K*t);
    Pf2 = d(1)*K*R*b1_ + K*t;
    Pf2 = Pf2./Pf2(3,:);
    Pf1 = d(2)*K*R'*b2_ - K*R'*t;
    Pf1 = Pf1./Pf1(3,:);
%     Pf2 = Pf2/norm(Pf2);
%     Pf1 = Pf1/norm(Pf1);
    error(i) = norm(Pf2-b2)+norm(Pf1-b1);
    if error(i) < thresh
        num_inliners = num_inliners + 1;
        bestInliners = [bestInliners i];
    end
end
