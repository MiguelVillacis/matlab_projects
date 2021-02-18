clear, clc;
addpath('func/');
addpath('func/geo/');
addpath('func/mesh/');
addpath('func/classes/');
addpath('src/');
addpath('ui/');

% uipt;

% g = rhplate(50,200,20);
g = holeplate(40,200,100,20,10);
props = [210000, 0.3];
boundary= [];

model = createpde;
geometryFromEdges(model, g);

mesh = TriMeshObj(g,3);
mesh.generateMesh;

plot(mesh.Nodes(1,:), mesh.Nodes(2,:), 'k.', 'MarkerSize', 2);
axis equal;



femtrisolver(mesh,boundary,props);


