function ind = ChieralityCheck(P, R1, t1, R2, t2)

C1_P = R1'*(P - repmat(t1, [1 size(P,2)]));
C2_P = R2'*(P - repmat(t2, [1 size(P,2)]));

ind = find(C1_P(3,:) > 0 & C2_P(3,:) > 0);

end
