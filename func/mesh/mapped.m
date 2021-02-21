function mpnodes = mapped(mmesh,props,load,flag)

switch flag
    case 1
        ef=1;
        efx=3;
        v=3;
    case 2
        ef=1;
        efx=2;
        v=2;
    case 3
        ef=1;
        efx = [2,3,4];
        v=2;
    case 4
        ef=[4,5];
        efx = 3;
        v=3;
end

structuralProperties(mmesh,'YoungsModulus',props(1),'PoissonsRatio',props(2));
structuralBC(mmesh,'Edge',efx,'XDisplacement',0);
structuralBC(mmesh,'Vertex',v,'YDisplacement',0);
structuralBoundaryLoad(mmesh,'Edge',ef,'SurfaceTraction',[load;0]);
generateMesh(mmesh);
mpnodes = solve(mmesh);

end

