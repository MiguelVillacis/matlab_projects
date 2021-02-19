clear, clc;
addpath('func/');
addpath('func/geo/');
addpath('func/mesh/');
addpath('func/classes/');
% addpath('src/');
addpath('ui/');

% uipt;

% g = holeplate(40,200,100,20,15);
% g = dholeplate(40,200,100,20,100,20,15);
g = tplate(40,200,100,20,15);
% g = rhplate(40,200,15);

props = [200E3 0.25];
load = 100;
flag = 3;
mmesh = gtrimesh(g,3);


mp = mapped(mmesh,props,load,flag);

plot(mp.Displacement.x);


