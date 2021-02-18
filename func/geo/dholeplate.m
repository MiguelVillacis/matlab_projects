function g = dholeplate(W,L,cx1,r1,cx2,r2,d)

rc = [3 4 0 L L 0 W/2 W/2 -W/2 -W/2]';
c1 = [1 cx1 ((d/2)+r1) r1]';
c2 = [1 cx2 -((d/2)+r2) r2]';

c1 = [c1;zeros(length(rc) - length(c1),1)];
c2 = [c2;zeros(length(rc) - length(c2),1)];

gd = [rc,c1,c2];

ns = char('r1','c1','c2');

ns = ns';

sf = '(r1-c1)-c2';

[dl,bt] = decsg(gd,sf,ns);

g = csgdel(dl,bt); 

end

