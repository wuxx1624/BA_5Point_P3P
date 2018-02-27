function Draw3DWorld(AbsolutePoses, FeatureBag, W_T_C1)

% for k = 1
%     W_T_C1 = [AbsolutePoses(:,1:3,k)' -AbsolutePoses(:,1:3,k)'*AbsolutePoses(:,4,k)];
%     DrawAxis(W_T_C1(:,4), W_T_C1(:,1:3), 0.2);
%     DrawCamera(W_T_C1(:,4), W_T_C1(:,1:3), 0.2, 'k');
% end

for k = 1:size(AbsolutePoses,3)
    W_T_Ck = W_T_C1*[InversePose(AbsolutePoses(:,:,k));zeros(1,3) 1];
%     W_T_Ck = W_T_C1*[C1_T_Ck;zeros(1,3) 1];
    DrawAxis(W_T_Ck(:,4), W_T_Ck(:,1:3), 0.2);
    DrawCamera(W_T_Ck(:,4), W_T_Ck(:,1:3), 0.2, 'k');
end

for k = 1:size(FeatureBag,2)
    C1_T_Ck = [AbsolutePoses(:,1:3,FeatureBag(4,k))' -AbsolutePoses(:,1:3,FeatureBag(4,k))'*AbsolutePoses(:,4,FeatureBag(4,k))];
    W_T_Ck = W_T_C1*[C1_T_Ck;zeros(1,3) 1];
    W_fk = W_T_Ck*[FeatureBag(1:3,k);1];
    plot3(W_fk(1),W_fk(2),W_fk(3), 'ro', 'markerSize', 3, 'markerfacecolor', 'r');    
end

end