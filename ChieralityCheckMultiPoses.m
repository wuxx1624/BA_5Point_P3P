function infrontof = ChieralityCheckMultiPoses(feat, AbsolutePoses, poses)

infrontof = 1;

if feat(3) <= 1e-3
    infrontof = 0;
    return;
end
for k = poses
    Rti = AbsolutePoses(:,:,k)*[AbsolutePoses(:,1:3,feat(4))' -AbsolutePoses(:,1:3,feat(4))'*AbsolutePoses(:,4,feat(4));zeros(1,3) 1];
    Ci_feat = Rti*[feat(1:3);1];
    if(Ci_feat(3) <= 1e-3)
        infrontof = 0;
        break;
    end
end
