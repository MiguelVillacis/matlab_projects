function femtrisolver(mesh,boundary, props)

% Constants
nel = 3;
dof = 2;
gqp = [0.5 0 1/6;0 0.5 1/6; 0.5 0.5 1/6];


E = props(1);       % Young modulus in MPa
v = props(2);       % Poisson coefficient

nnot = length(mesh.Nodes);      % Total number of nodes
nelt = length(mesh.Elements);   % Total number of elements

% Stiffness matrix initialization
k = zeros(2*nnot,2*nnot);


% Elasticity matrix initialization [b]
c = E/(1-(v^2));
b = c*[1 v 0. ;v 1 0. ;0. 0. .5*(1-v)];

% Structure of node coordinates for element
for i=1:nelt
    x = [mesh.Nodes(1, mesh.Elements(1, i)); ...
         mesh.Nodes(1, mesh.Elements(2, i)); ...
         mesh.Nodes(1, mesh.Elements(3, i))];
    y = [mesh.Nodes(2, mesh.Elements(1, i)); ...
         mesh.Nodes(2, mesh.Elements(2, i)); ...
         mesh.Nodes(2, mesh.Elements(3, i))];
     eCoord(i).e = [x,y];
end

% Degree of freedom matrix initialization
nf = [1:2:2*nnot;2:2:2*nnot]';

% Degree of freedom connectivity matrix inicialization
vdir = zeros(nelt, 6);

for i=1:nelt
    ci = 0;
    for j=1:3
        for k=1:2
            ci = ci+1;
            vdir(i,ci)=nf(mesh.Elements(j,i),k);
        end
    end
end

% General stiffness matrix integration (Gauss Quadrature) 
for i=1:nelt
    ke = zeros(6,6);        % Local stiffness matrix - Element
    
    
        
end

end
