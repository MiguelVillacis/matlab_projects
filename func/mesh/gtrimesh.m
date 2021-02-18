function gmesh = gtrimesh(g, hmax)

model = createpde;
geometryFromEdges(model,g);
gmesh = generateMesh(model, 'Hmax', hmax, 'GeometricOrder','quadratic');

end

