function plotrimesh(ax, mesh, flag)

edges = [mesh.Elements([1,2],:)'; ...
         mesh.Elements([2,3],:)'; ...
         mesh.Elements([3,1],:)'];
     
[u,~,n] = unique(sort(edges,2),'rows'); 
counts = accumarray(n(:), 1);
uedges = u(counts==1,:);
iedges = u(counts==2,:);

hold(ax, 'on');
 axis(ax,'equal');

if flag(1) == 1
    plot(ax, mesh.Nodes(1,:),mesh.Nodes(2,:), 'k.', 'MarkerSize',1)
end

if flag(2) == 1
    for i=1:length(uedges)

        x = [mesh.Nodes(1,uedges(i,1)) ...
             mesh.Nodes(1,uedges(i,2))];
        y = [mesh.Nodes(2,uedges(i,1)) ...
             mesh.Nodes(2,uedges(i,2))];

        plot(ax, x, y, 'LineWidth', 0.5, 'Color', 'r');

    end
end

if flag(3) == 1
    for i=1:length(iedges)

        x = [mesh.Nodes(1,iedges(i,1)) ...
             mesh.Nodes(1,iedges(i,2))];
        y = [mesh.Nodes(2,iedges(i,1)) ...
             mesh.Nodes(2,iedges(i,2))];

        plot(ax, x, y, 'LineWidth', 0.5, 'Color', 'k');

    end
end
    
end



