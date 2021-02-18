function g = rhplate(W, L, r)

rc1 = [3 4 0 (L-W/2) (L-W/2) 0 W/2 W/2 -W/2 -W/2]';
c1 = [1 (L-W/2) 0  W/2]';
c2 = [1 (L-W/2) 0  r]';

c1 = [c1;zeros(length(rc1) - length(c1),1)];
c2 = [c2;zeros(length(rc1) - length(c2),1)];

gd = [rc1,c1,c2];

ns = char('rc1','c1','c2');

ns = ns';

sf = '(rc1+c1)-c2';

[dl,bt] = decsg(gd,sf,ns);

g = csgdel(dl,bt); 

end

