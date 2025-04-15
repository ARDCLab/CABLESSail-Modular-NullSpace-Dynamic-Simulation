function DCM = normalizeDCM(C)

c1 = C(:,1);
c2 = C(:,2);
c1 = c1/norm(c1);
c2 = c2/norm(c2);
c3 = cross(c1,c2);
DCM = [c1, c2, c3];

end