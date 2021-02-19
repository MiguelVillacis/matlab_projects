function model = gtrimesh(g, hmax)

model = createpde('Structural', 'static-planestress');
geometryFromEdges(model,g);
generateMesh(model, 'Hmax', hmax, 'GeometricOrder','quadratic');

end

