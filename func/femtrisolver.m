function results = femtrisolver(mesh,load, props, mmesh, flag)

% Constants
nxe = 6;
dof = 2;
% Gauss quadrature points
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

% Derivative shape matrix
ds = zeros(2,nxe);
dc = zeros(2,nxe);

% Element area matrix
ae = zeros(nelt,1);

se=zeros(nelt,nxe,nxe);
fe=zeros(nelt,nxe);
mp = mapped(mmesh,props,load,flag);

% Structure of node coordinates for element
for i=1:nelt
    x = zeros(nxe,1);
    y = zeros(nxe,1);
    for j=1:nxe
        x(j) = mesh.Nodes(1, mesh.Elements(j, i));
        y(j) = mesh.Nodes(2, mesh.Elements(j, i));
    end
     eCoord(i).e = [x,y];
end

% Degree of freedom matrix initialization
nf = [1:2:2*nnot;2:2:2*nnot]';

% Degree of freedom connectivity matrix inicialization
vdir = zeros(nelt, 12);


for i=1:nelt
    ci = 0;
    for j=1:6
        for k=1:2
            ci = ci+1;
            vdir(i,ci)=nf(mesh.Elements(j,i),k);
        end
    end
end

constrained=zeros(nnot,1);
constrained(11:15)=1;
n_const=sum(constrained);

var_map=zeros(nnot,1);
n_var=0;
for i=1:nnot
    if(constrained(i)==0)
        n_var=n_var+1;
        var_map(i)=n_var;
    end
end


% defined current densities
Je=zeros(nelt,1);
Je(1)=10;

% General stiffness matrix integration (Gauss Quadrature) 
for i=1:nelt
    ke = zeros(6,6);        % Local stiffness matrix - Element
    
    for ig=1:3
        
        xi = gqp(ig,1);
        eta = gqp(ig,2);
        w = gqp(ig,3);
        
        %Derivative Shape matrix
        ds(1,1)=4*xi-1;
        ds(1,2)=4*eta;
        ds(1,3)=4-8*xi-4*eta;
        ds(1,4)=0;
        ds(1,5)=-4*eta;
        ds(1,6)=4*eta+4*xi-3;
        ds(2,1)=0;
        ds(2,2)=4*xi;
        ds(2,3)=-4*xi;
        ds(2,4)=4*eta-1;
        ds(2,5)=4-4*xi-8*eta;
        ds(2,6)=4*eta+4*xi-3;
        

        % Jacobian matrix
        jm = ds*eCoord(i).e;

        % Jacobian matrix determinant
        djm = det(jm);

        % Jacobian matrix inversed
        ijm = inv(jm);
        
        % Area of element
        ae(i) = djm/2;
        
        dc=jm\ds;
        
        for k=1:nxe
            for j=1:nxe
                se(i,k,j)=se(i,k,j)+djm/6*(dc(1,k)*dc(1,j)+dc(2,k)*dc(2,j));
            end
        end
        
        fe(i,1)=fe(i,1)+Je(ne)*djm/6*xi*(2*xi-1);
        fe(i,2)=fe(i,2)+Je(ne)*djm/6*xi*eta*4;
        fe(i,3)=fe(i,3)+Je(ne)*djm/6*xi*(1-xi-eta)*4;
        fe(i,4)=fe(i,4)+Je(ne)*djm/6*eta*(2*eta-1);
        fe(i,5)=fe(i,5)+Je(ne)*djm/6*eta*(1-xi-eta)*4;
        fe(i,6)=fe(i,6)+Je(ne)*djm/6*(1-xi-eta)*(2*(1-xi-eta)-1);
        
    end
       
end

results = mp;

end
