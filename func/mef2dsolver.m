function [dpl,S,Ep, res] = mef2dsolver(E,Lx,vu,Load_X,Ly,Nelx,Nely,thick)


%
%% %%%%%%%%%%%%%%%% P R E P R O C E S A D O %%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%*************************************************************************%
%                        Entrada de datos                                 %
%*************************************************************************%


%**************************************************************************
%-----------------------------------------
%  Inicializacion de variables
%-----------------------------------------
%
res = Result();

Nodof=2;                % Número de grados de libertad por nodo.
NNdx=Nelx+1;            % Número de nodos en X
NNdy=Nely+1;            % Número de nodos en Y
Dof=Nodof*4;            % Número de grado de libertad por elemento.

Loadx=Load_X*thick*Ly;  % Esfuerzo en la cara lateral F = P*A = P*thick*Ly [N]
LoadNdx=Loadx/NNdy;     % Fuerza aplicada por nodo [N]

Nel=Nelx*Nely;          % Número de elementos en la malla.
Node_Total=Nel*4;       % Número total de nodos por el número de elementos.

NNtotal=NNdy*NNdx;      % Número total de nodos en la malla.
scalex=Lx*0.3;          % Siguiente cuadro de imágenes x
scaley=Ly*0.3;          % Siguiente marco de imagen y

KK=zeros(Nodof*NNtotal,Nodof*NNtotal); %Inicialización de la matriz de rigidez total


res.Nodof = Nodof;
res.NNdx = NNdx;
res.NNdy = NNdy;
res.Dof = Dof;
res.Loadx = Loadx;
res.LoadNdx = LoadNdx;
res.Nel = Nel;
res.Node_Total = Node_Total;
res.NNtotal = NNtotal;
res.scalex = scalex;
res.scaley = scaley;


%*************************************************************************%
%           La inicialización de la matriz de elasticidad [D]             %
%*************************************************************************%
%
   c=E/(1.-vu*vu);
   D=c*[1 vu 0. ;vu 1 0. ;0. 0. .5*(1.-vu)];

%
%*************************************************************************%
%                      Discretización del dominio                      %
%*************************************************************************%
%
%--------------------------------------------------------------
% Tabla de coordenadas global por nodos                
%--------------------------------------------------------------

% En este punto se obtienen las cordenadas de cada nodo.
% La matriz horizontal resultado tien la forma  [nodo cx cy]';

Cx=0;
Cy=0;
NumNd=0;
% T_Coord = zeros([NNtotal 3]);

for Ix=1:NNdx
    Cx=(Lx/Nelx)*(Ix-1);
    for Iy=1:NNdy
        Cy=(Ly/Nely)*(Iy-1);
        NumNd=NumNd+1;
        T_Coord(1,NumNd)=NumNd;
        T_Coord(2,NumNd)=Cx;
        T_Coord(3,NumNd)=Cy;
    end
end

% Posicion de nodos.
res.t_coord = T_Coord;

x=T_Coord(2,:);
y=T_Coord(3,:);
Coord=[x;y]';                   %% NOTA: No se entiende donde utiliza esta matriz

%*************************************************************************%
%                  El cálculo de la matriz de rigidez total [K]           %
%*************************************************************************%
%
%--------------------------------------------------------------
%  Matriz de conexión global de elementos [T_Conex]                     
%--------------------------------------------------------------
%
% NOTA: SE CREA SIN DEFINICION PREVIO T_Conex, y se le asigna 

count=0;  % Contador de columnas
count1=0; % Contador de elementos por columna
count2=0; % Contador de elementos
T_Conex = zeros(Nel, 5);

for elementc = 1:Nel
    count2=count2+1;
    T_Conex(elementc,1) = count2;                      %Elemento de referencia
    T_Conex(elementc,2) = elementc + count;            %El mismo nodo de referencia
    T_Conex(elementc,3) = elementc + Nely + 1 + count; %Elemento adyacente derecho
    T_Conex(elementc,4) = elementc + Nely + 2 + count; %Elemento Diagonal derecho
    T_Conex(elementc,5) = elementc + 1 + count;        %Elemento Superior
    count1 = count1+1;
    
    if count1 == Nely % Reinicia el contador de columna
       count = count + 1;
       count1 = 0;
    end

end

% Nota: esta seccion crea una matriz vertical con la siguiente estructura [nref, n1, n2, n3, n4]
% Que representan el nodo de referencia y los ccuatro nodos adyacentes proximos.

%--------------------------------------------------------------
% Tabla de coordenadas por elementos [elemnts]                     
%--------------------------------------------------------------

countn=0;
for i=1:Nel 
    for j=2:5
        countn=countn+1;
        nd(1,countn) = T_Conex(i,j);
        Coonec(i,j-1)=  T_Conex(i,j);
    end
end

% nd es un ventor fila que lista secuencialmente los nodos que conforman un elemento 
% nd = [n1r1 n2r1 n3r1 n4r1 n1r2 n2r2 n3r2 n4r2 ...]; Secuencia lineal por elemento.
% Coonec es un vector que elimina la primera column de T_Conex vertical.
% Coonec = T_Conex(:,2:5)

T_Coord=T_Coord';               % Convierte en una matriz vertical T_Coord
Node_line=nd';                  % Convierte en un vector horizontal el vector nd
Node_Coord=zeros(NNtotal,3);    % Inicializa un vector de tres columnas y NNtotal filas.


% Obtener secuencia de 4 x 2, con las cordenadas de los nodos por elemento 
for i=1:Node_Total
     Node_Coord(i,:)=T_Coord(Node_line(i),:);
end

% Crear una estructura que contiene los elementos y las cordenadas de 
% los nodos que lo conforman.
    % e[x1,y1
    %   x2,y2
    %   x3,y3
    %   x4,y4]

count1=1;
count2=4;
for i=1:Nel
    element(i).e = Node_Coord(count1:count2,2:3);
    count1=count1+4;
    count2=count2+4;
end


%----------------------------------------
% La tabla de conectividad [vdir]
%----------------------------------------
% Crean una matriz de 2 columnas y NTotal filas y numera los grados de libertad
nf=zeros(NNtotal,Nodof); %Inicialización del vector Grados de libertad
countn=0;
 for i=1:NNtotal
    for j=1:Nodof
        countn=countn+1;
        nf(i,j)=countn;
    end
 end 

vdir=zeros(1,Dof); %Inicializando la tabla de conectividad

%Se crea una tabla con el ID de los grados de libertad an base a T_Conex
 for i=1:Nel
     countn=0;
    for j=2:5
        for k=1:2
            countn=countn+1;
            vdir(i,countn)=nf(T_Conex(i,j),k);
        end
    end
 end
%

 for i=1:Nel
     Ke=zeros(Dof,Dof); %La matriz de rigidez local inicializada por ceros
   
%-----------------------------------------------------------------------------
% Cuadratura gaussiana para dos puntos  (Gqp) xi et eta, wi et wj 
%-----------------------------------------------------------------------------

     Gqp=[-0.5774 1; %Dos puntos con sus pesos en cuadratura gaussiana para
          0.5774 1]; %Integración numerica
  
    for ig=1:2
     wi=Gqp(ig,2);
    for jg=1:2
     wj=Gqp(jg,2);  
     
          xi=Gqp(ig,1);
          eta=Gqp(jg,1);
     
%------------------------------------
% La función de forma (fun)
%------------------------------------

    fun = 0.25*[(1.- xi - eta + xi*eta);(1.+ xi - eta - xi*eta);...
    (1.+ xi + eta + xi*eta);...
    (1.- xi + eta - xi*eta)];
%
%----------------------------------------------------------------
% Derivada de la forma función (de) respeto por xi y eta
%----------------------------------------------------------------
%
    der = 0.25*[-(1-eta) (1-eta) (1+eta) -(1+eta);...
                -(1-xi) -(1+xi) (1+xi) (1-xi)];
%
%---------------------------------------------
% Cálculo de la matriz jacobiana (jac)
%---------------------------------------------
%
    jac=der*element(i).e; 
%
%------------------------------------------
% Cálculo del determinante de (jac)
%------------------------------------------
%
    d=det(jac); 
%
%------------------------------------
% Cálculo de la inversa de (jac)
%------------------------------------
%
    jac_inv=inv(jac); 
%
%-----------------------------------------------------
% Cálculo de la matriz base [Base]
%-----------------------------------------------------
%     
      Base = zeros(4,Dof);
          for m=1:4
          k=2*m;
          l=k-1;
          dx=der(1,m);
          Base(1,l)=dx;
          Base(3,k)=dx;
          dy=der(2,m);
          Base(2,l)=dy;
          Base(4,k)=dy;
          end
%
%-----------------------------------------------------------------
% Cálculo de la matriz de transformación geométrica [Tgeo]
%-----------------------------------------------------------------
%       
          Tgeo = zeros(3,4);
          for j=1:2
              k=j*2;
              l=k-1;
              r=2+j;
              dx=jac_inv(1,j);
              dy=jac_inv(2,j);
              Tgeo(1,j)=dx;
              Tgeo(3,r)=dx;
              Tgeo(2,r)=dy;
              Tgeo(3,j)=dy; 
          end
%
%-----------------------------------------------------
% Cálculo de la matriz elemental local [Ke]]
%-----------------------------------------------------
%
      Ke=Ke + thick*Base'*Tgeo'*D*Tgeo*Base*d*wi*wj;
      
    end
    end
%
%-----------------------------------------------------
% Ensamblaje de la matriz de rigidez total [K]
%-----------------------------------------------------
% 
        for k=1:Dof
          for j=1:Dof
           KK(vdir(i,k),vdir(i,j))= KK(vdir(i,k),vdir(i,j))+ Ke(k,j);
          end
        end
 end
 format shortE

%*************************************************************************%
%                    Montaje del sistema [K] [q] = [F]                    %
%*************************************************************************%

% Condición de frontera
%------------------------------------
% Desplazamiento: Incrustación en Lx = 0 
%------------------------------------

for i=1:NNtotal
    if T_Coord(i,2) == 0
       nf(i,:) = [0 0];
    end
end

nfm=nf;
nf=nf';  
nf=nf(:);  %El vector de grado reducido de libertad
n=length(nf);
m=0;

for i=1:n
    if nf(i)== 0
        m = m+1;
    end
end

% Condición de frontera
%---------------------------------------
% Vector de fuerza impuesta [f]
%---------------------------------------

ff=zeros(NNtotal,Nodof); %Inicialización vectorial
for i=1:NNtotal
    if T_Coord(i,2) == Lx
       ff(i,:) = [LoadNdx 0];
    end
end

ff=ff';
ff=ff(:); %Impulsor general de los esfuerzos

Nx=n-m;
f=zeros(1,Nx)';
Ny=0;

for i=m+1:n
    Ny=Ny+1;
    f(Ny,1)=ff(i,1); %El vector de fuerza reducida
end

%-----------------------------------------------------------------
% La matriz global de rigidez reducida
%-----------------------------------------------------------------

K=KK(m+1:n,m+1:n);
%
format short g
%
%% %%%%%%%%%%%%%%%%%%%% P R O C E S A D O %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%*************************************************************************%
%                    Resolución del sistema [K] [q] = [F]                %
%*************************************************************************%
%
%--------------------------------------------------------------
% El cálculo de desplazamientos desconocidos
%--------------------------------------------------------------

q = K\f; %Resolución por el método de eliminación de Gauss 

n=length(nfm);
Node=0;
count=0;

for i=1:n
    Node=Node+1;
    if nfm(i,1)==0
        deltax = 0;
        deltay = 0;
    else
        count=count+1;
        deltax=q(count,1);
        deltay=q(count+1,1);
        count=count+1;
    end    
        displ(i,1)= Node;
        Noeud(i,1)= Node;
        displ(i,2)= deltax;
        Ux(i,1)= deltax;
        displ(i,3)= deltay;
        Uy(i,1)=deltay;
end

disp('------------------------------------ ')
disp('             Los resultados ')
disp('------------------------------------ ')
disp(' ')
disp('El desplazamiento [q] por nodo');
disp(' ')
tableu =table(Noeud,Ux,Uy);
res.TablaQ = tableu;
dpl=[Ux Uy];


%--------------------------------------------------------------
% El cálculo de la reacción en el nodo [F]
%--------------------------------------------------------------

displ=displ(1:n,2:3);
displ=displ';
U=displ(:);

F=KK*U;
Force=reshape(F,2,n);
Force=Force';
for i=1:n
    
Fnodal(i,1)=i;
     if abs(Force(i,1))< 0.01
         Fnodal(i,2)=0;
     else
Fnodal(i,2)=Force(i,1);
     end
     if abs(Force(i,2))< 0.01
          Fnodal(i,3)=0;
     else
Fnodal(i,3)=Force(i,2);
     end
end

res.Fnodal = Fnodal;
%format short g
%fprintf('Reaction [F] by node\n');
%fprintf ('------------------------------------------\n');
%fprintf ('           Node        Fx        Fy\n');
%fprintf ('------------------------------------------\n');
%disp(Fnodal)   


%% %%%%%%%%%%%%%% P O S T - P R O C E S A D O %%%%%%%%%%%%%%%%%%%%%%%%%%
%
%*************************************************************************%
%       El cálculo de la deformación [Ep] y la tensión [Sigma]            %
%*************************************************************************%

count=0;
for i=1:Nel
    Uel=zeros(Dof,1); 
    count=count+1;
 for j=1:Dof
     if vdir(count,j)==0 
        Uel(j)=0.;
     else %
        Uel(j) =U(vdir(count,j),1); %Reorganización de desplazamiento vector por elemento
     end
 end


Gqp=[0 2]; %Un punto de Gauss y sus pesos.

%for i=1:Nel;
for ig=1:1
    wi=Gqp(ig,2);
    for jg=1:1
        wj=Gqp(jg,2); 

        xi=Gqp(ig,1);
        eta=Gqp(jg,1);
        %     
        %------------------------------------
        % La función de forma (fun) 
        %------------------------------------
        %
            fun = 0.25*[(1.- xi - eta + xi*eta);(1.+ xi - eta - xi*eta);...
            (1.+ xi + eta + xi*eta);...
            (1.- xi + eta - xi*eta)];
        %
        %----------------------------------------------------------------
        % La derivada de las formas focntion (der) respeto por xi y eta
        %----------------------------------------------------------------

            der = 0.25*[-(1-eta) (1-eta) (1+eta) -(1+eta);...
                        -(1-xi) -(1+xi) (1+xi) (1-xi)];

        %---------------------------------------------
        % Cálculo de la matriz jacobiana (jac)
        %---------------------------------------------

            jac=der*element(i).e; 

        %------------------------------------
        % Cálculo de la inversa de (jac)
        %------------------------------------

            jac_inv=inv(jac); 

        %-----------------------------------------------------
        % El cálculo de la matriz base [Base]
        %-----------------------------------------------------
        %     
        Base = zeros(4,Dof);
        for m=1:4
            k=2*m;
            l=k-1;
            dx=der(1,m);
            Base(1,l)=dx;
            Base(3,k)=dx;
            dy=der(2,m);
            Base(2,l)=dy;
            Base(4,k)=dy;
        end       
        %
        %-----------------------------------------------------------------
        % Cálculo de la matriz de transformación geométrica [Tgeo]
        %-----------------------------------------------------------------
        %       
        Tgeo = zeros(3,4);
        for j=1:2
          k=j*2;
          l=k-1;
          r=2+j;
          dx=jac_inv(1,j);
          dy=jac_inv(2,j);
          Tgeo(1,j)=dx;
          Tgeo(3,r)=dx;
          Tgeo(2,r)=dy;
          Tgeo(3,j)=dy; 
        end
        %          
        %-----------------------------------------------------
        %   Cálculo de la deformación [Epsilon]
        %-----------------------------------------------------
        % 
            Epsilon = Tgeo*Base*Uel;  

        %-----------------------------------------------------
        %   Cálculo de la tensión [Sigma]
        %-----------------------------------------------------

            sigma=D*Epsilon;
       
    end
end
       Ep(i,:)=Epsilon;
       S(i,:)=sigma; %
end


%-----------------------------------------------------
%    La tabla de resultados
%-----------------------------------------------------

for i=1:Nel
       Element(i,1)=i;
       Sx(i,1)= S(i,1);
       Ex(i,1)= Ep(i,1);
       Sy(i,1)= S(i,2);
       Ey(i,1)= Ep(i,2);
       Sxy(i,1)= S(i,3);
       Exy(i,1)= Ep(i,3);
end
% %
% format short
% %
% disp(' ')
% disp('Deformación en el centroide del elemento de referencia.')
% disp(' ')
tableuDef =table(Element,Ex,Ey,Exy);
% disp(tableu)
res.TablaDef = tableuDef;

% disp(' ')
% disp('Restricción al elemento central de referencia.')
% disp(' ')
tableuRes =table(Element,Sx,Sy,Sxy);
% disp(tableu)
res.TableRes = tableuRes;

% %*************************************************************************%
% %                            Plot displacement                            %
% %*************************************************************************%

%figure
%hold on
pini=[x ; y];
pini=pini';
dispx=displ(1,:);
dispy=displ(2,:);
pfin=[dispx ; dispy];
pfin=pfin';
DISP=pini+pfin;
x=DISP(:,1);
y=DISP(:,2);
Coordv=[x y];
U2=DISP(:,1);
x=x';
y=y';
% for i=NNtotal
%     plot(x,y,'.');
% end
cmin=min(x);
cmax=max(x);
%caxis([cmin cmax]);
sx=x(1,end)+scalex;
sy=y(1,end)+scaley;
%axis ([0 sx -scaley sy]);
sx=x(1,end)+0.1;
sy=y(1,end)+0.1;
%axis ([0 sx -0.1 sy]);

%patch('Faces',Coonec,'Vertices',Coordv,'FaceVertexCData',U2,'FaceColor','interp','Marker','.');
%colorbar;

%title ('Deplacement Ux');

% count=0;
% for i=1:Nelx+1
%     count=count+1;
%     for j=1:Nely 
%         px1=x(1,count);
%         py1=y(1,count);
%         px2=x(1,count+1);
%         py2=y(1,count+1);
%         plot([px1 px2],[py1 py2],'r','LineWidth',2)
%         count=count+1;
%     end
% end
% 
% count=0;
% 
% for i=1:Nely+1
%     count=i;
%     for j=1:Nelx
%         px1=x(1,count);
%         px2=x(1,count+Nely+1);
%         py1=y(1,count);
%         py2=y(1,count+Nely+1);
%         plot([px1 px2],[py1 py2],'r','LineWidth',2)
%         count=count+Nely+1;
%     end
% end
% 
% for ti=1:NNtotal
%     
%     text(Coordv(ti,1),Coordv(ti,2),num2str(ti),'Color','k','FontSize',14)
%     
% end

end

