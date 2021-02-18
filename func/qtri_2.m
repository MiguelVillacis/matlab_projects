clear;

% define problem

nnodes=15;
nelms=4;
% 3 nodes per element
npe=6;

% initialise
xy=zeros(nnodes,2);
node_map=zeros(nelms,npe);

% set co-ordinates
xy(1,1:2)=[0 0];
xy(2,1:2)=[0.5 0];
xy(3,1:2)=[0.5 0.5];
xy(4,1:2)=[1 0];
xy(5,1:2)=[1 0.5];
xy(6,1:2)=[1 1];
xy(7,1:2)=[1.5 0];
xy(8,1:2)=[1.5 0.5];
xy(9,1:2)=[1.5 1];
xy(10,1:2)=[1.5 1.5];
xy(11,1:2)=[2 0];
xy(12,1:2)=[2 0.5];
xy(13,1:2)=[2 1];
xy(14,1:2)=[2 1.5];
xy(15,1:2)=[2 2];

node_map(1,:)=[1 2 3 4 5 6];
node_map(2,:)=[6 5 8 4 7 11];
node_map(3,:)=[11 12 8 13 9 6];
node_map(4,:)=[6 9 10 13 14 15];

% plot mesh
xyne=zeros(npe+1,2);
figure(1)
hold on
for ne=1:nelms
    xyplot(1,:)=xy(node_map(ne,1),:);
    xyplot(2,:)=xy(node_map(ne,2),:);
    xyplot(3,:)=xy(node_map(ne,4),:);
    xyplot(4,:)=xy(node_map(ne,5),:);
    xyplot(5,:)=xy(node_map(ne,6),:);
    xyplot(6,:)=xy(node_map(ne,3),:);
    xyplot(7,:)=xy(node_map(ne,1),:);
    
    plot(xyplot(:,1),xyplot(:,2),'Color',[0 0 0],'Marker','o','MarkerFaceColor',[0 0 0]);
end
hold off
axis equal;
title('Input mesh geometry')

% define constraints on solution
% constrained values have A=0
constrained=zeros(nnodes,1);
constrained(11:15)=1;
n_const=sum(constrained);

var_map=zeros(nnodes,1);
n_var=0;
for i=1:nnodes
    if(constrained(i)==0)
        n_var=n_var+1;
        var_map(i)=n_var;
    end
end
soln=zeros(n_var,1);

% defined current densities
Je=zeros(nelms,1);
Je(1)=10;

% Construct Local Stiffness Matrices

b=zeros(3,1);  % column matrices!
c=zeros(3,1);
se=zeros(nelms,npe,npe);
area=zeros(nelms,1);
n=zeros(npe,1);
diffmat_local=zeros(2,npe);
jack=zeros(2,2);
jinv=zeros(2,2);
diffmat_cart=zeros(2,npe);
dadx=zeros(1,npe);
dady=zeros(1,npe);
dummy=zeros(npe,npe);
xy_local=zeros(npe,2);
fe=zeros(ne,npe);

for ne=1:nelms
    n(:)=node_map(ne,:);  % find global node number for local index
    % Getting coordinates of any 6 points of triangle
    for i=1:npe
        xy_local(i,:)=xy(n(i),:);
    end
    
    dummy(:,:)=0;
    
    for pnt=1:3
        if (pnt==1)
            epsilon=1/2;
            eta=0;
        elseif (pnt==2)
            epsilon=0;
            eta=1/2;
        elseif (pnt==3)
            epsilon=1/2;
            eta=1/2;
        end
        diffmat_local(1,1)=4*epsilon-1;
        diffmat_local(1,2)=4*eta;
        diffmat_local(1,3)=4-8*epsilon-4*eta;
        diffmat_local(1,4)=0;
        diffmat_local(1,5)=-4*eta;
        diffmat_local(1,6)=4*eta+4*epsilon-3;
        diffmat_local(2,1)=0;
        diffmat_local(2,2)=4*epsilon;
        diffmat_local(2,3)=-4*epsilon;
        diffmat_local(2,4)=4*eta-1;
        diffmat_local(2,5)=4-4*epsilon-8*eta;
        diffmat_local(2,6)=4*eta+4*epsilon-3;
        
        jack=diffmat_local*xy_local;
        jackdet=jack(1,1)*jack(2,2)-jack(1,2)*jack(2,1);
        jinv(1,1)=jack(2,2);
        jinv(2,1)=-jack(2,1);
        jinv(1,2)=-jack(1,2);
        jinv(2,2)=jack(1,1);
        jinv=jinv/jackdet;
        area(ne)=jackdet/2;
        diffmat_cart=jinv*diffmat_local;
        for k=1:npe
            for j=1:npe
                se(ne,k,j)=se(ne,k,j)+jackdet/6*(diffmat_cart(1,k)*diffmat_cart(1,j)+diffmat_cart(2,k)*diffmat_cart(2,j));
            end
        end
        fe(ne,1)=fe(ne,1)+Je(ne)*jackdet/6*epsilon*(2*epsilon-1);
        fe(ne,2)=fe(ne,2)+Je(ne)*jackdet/6*epsilon*eta*4;
        fe(ne,3)=fe(ne,3)+Je(ne)*jackdet/6*epsilon*(1-epsilon-eta)*4;
        fe(ne,4)=fe(ne,4)+Je(ne)*jackdet/6*eta*(2*eta-1);
        fe(ne,5)=fe(ne,5)+Je(ne)*jackdet/6*eta*(1-epsilon-eta)*4;
        fe(ne,6)=fe(ne,6)+Je(ne)*jackdet/6*(1-epsilon-eta)*(2*(1-epsilon-eta)-1);
 
    end

end
 
% global stiffness matrix
sg=zeros(n_var,n_var);

for ne=1:nelms
    n(:)=node_map(ne,:);  % find global node number for local index
    for i=1:npe
        if (constrained(n(i))==0)
            nv=var_map(n(i));
            for j=1:npe
                if (constrained(n(j))==0)
                    nvj=var_map(n(j));
                    sg(nvi,nvj)=sg(nvi,nvj)+se(ne,i,j);
                end
            end
        end
    end
end

% global forcing matrix
fg=zeros(n_var,1);
for ne=1:nelms
    n(:)=node_map(ne,:);  % find global node number for local index
    for i=1:npe
        if (constrained(n(i))==0) 
            fg(n(i))=fg(n(i))+fe(ne,i);
        end
    end
end


soln=sg\fg;

a=zeros(nnodes,1);
for i=1:nnodes
    if (var_map(i)~=0)
        a(i)=soln(var_map(i));
    end
end
a;

    

