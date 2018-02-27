function T2 = InversePose(T1)

T2 = [T1(:,1:3)' -T1(:,1:3)'*T1(:,4)];

end