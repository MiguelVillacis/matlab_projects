function g = tplate(W,L,cx,r,d)

rc1 = [3 4 0 cx cx 0 W/2 W/2 -W/2 -W/2]';
rc2 = [3 4 0 L L 0 d/2 d/2 -d/2 -d/2]';
c1 = [1 cx ((d/2)+r) r]';
c2 = [1 cx -((d/2)+r) r]';

c1 = [c1;zeros(length(rc1) - length(c1),1)];
c2 = [c2;zeros(length(rc1) - length(c2),1)];

gd = [rc1,rc2,c1,c2];

ns = char('rc1','rc2','c1','c2');

ns = ns';

sf = '((rc1+rc2)-c1)-c2';

[dl,bt] = decsg(gd,sf,ns);

g = csgdel(dl,bt); 

end

