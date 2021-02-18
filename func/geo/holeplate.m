function g = holeplate(W,L, cx,cy,r)


r1 = [3 4 0 L L 0 W/2 W/2 -W/2 -W/2]';
c1 = [1 cx (-0.5*W+cy) r]';

c1 = [c1;zeros(length(r1) - length(c1),1)];

gd = [r1,c1];

ns = char('r1','c1');

ns = ns';

sf = 'r1-c1';

[dl,bt] = decsg(gd,sf,ns);

g = csgdel(dl,bt); % removes face boundaries


% figure
% pdegplot(dl,'EdgeLabels','on','FaceLabels','on')
% xlim([0,L*1.1])
% ylim([-(W/2)*1.2,(W/2)*1.2])
% 
% figure
% pdegplot(dl2,'EdgeLabels','on','FaceLabels','on')
% xlim([0,L*1.1])
% ylim([-(W/2)*1.2,(W/2)*1.2])


end

