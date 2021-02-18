% g = holeplate(40,150,100,20,7);
% g = dholeplate(40,100,40,10,70,15,20);
% g = tplate(50,250,120,45,30);
g = rhplate(40,150,15);

model = createpde;

geometryFromEdges(model,g);

phmesh = generateMesh(model, 'Hmax', 3, 'GeometricOrder','linear');
pdeplot(model)