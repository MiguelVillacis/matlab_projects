function femtrisolver(mesh,boundary, props)

% Constants
nxe = 6;
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
    
    for ig=1:3
        
        xi = gqp(ig,1);
        eta = gqp(ig,2);
        w = gqp(ig,3);
        
        %Jacobian Matrix
        x1 = eCoord(i).e(1,1);
        x2 = eCoord(i).e(3,1);
        x3 = eCoord(i).e(3,1);
        y1 = eCoord(i).e(1,2);
        y2 = eCoord(i).e(2,2);
        y3 = eCoord(i).e(3,2);

        % Jacobian matrix
        jm = [(x2-x1) (y2-y1); (x3-x1) (y3-y1)];

        % Jacobian matrix determinant
        djm = det(jm);

        % Jacobian matrix inversed
        ijm = inv(jm);
        
        
     
    end
    
        
end

end
