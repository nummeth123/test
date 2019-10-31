clear

freq = 1.0e6;
nCyc = 3;
c0 = 1500;
rho0 = 1000;


radCyl = 1e-3;%0.1e-3*1.0;
cCyl = 1500;
rhoCyl = 10;


radArray = 15e-3;
nArray = 120;


lambda = c0/freq;

dx = 0.05e-3*1;


arrayAlpha = (0:nArray-1)/nArray*2*pi;

%transducer location
tx = radArray*cos(arrayAlpha);
ty = radArray*sin(arrayAlpha);

transPoints = radArray*[cos(arrayAlpha.') sin(arrayAlpha.')].';

boundPoints = [-20 -20
    20 -20
    20 20
    -20 20].'*1e-3;
nBoundPoints = size(boundPoints,2);


circCyl = 2*pi*radCyl;
nCyl = ceil(circCyl/dx);
% return
cylAlpha = (0:nCyl-1)/nCyl*2*pi;

circPoints = radCyl*[cos(cylAlpha.') sin(cylAlpha.')].';

allPoints = [transPoints boundPoints circPoints];
plot(allPoints)
%%

boundSegs = [1 2; 
    2 3;
    3 4;
    4 1].'+nArray;

cylSegs = [1:nCyl; [2:nCyl 1]]+nArray+nBoundPoints;

allSegs = [boundSegs cylSegs];

%holes = [0.0 0.0].';   
holes = [];

savePoly( 'temp.poly', allPoints, allSegs, holes );



system(sprintf('pogoMesh temp.poly -s %f',dx))

delete('temp.poly')

m = loadVtuMesh('temp.vtu');

% figure
% plotMeshTri(m)
% axis equal
% drawnow
% return
%%

delete('temp.vtu')



courant = 0.3; % mix of sizes so need to be conservative


m.dt = dx/c0*courant; 
m.nt = round((radArray*3)/c0/m.dt);

nEls = size(m.elNodes,2);
m.matTypeRefs = ones(nEls,1);

%do cyl here
[ex, ey] = getElCents(m);
cylEls = find(((ex).^2 + (ey).^2) < radCyl.^2);



m.matTypeRefs(cylEls) = 2;

m.matTypes = cell(2,1);
m.matTypes{1}.paramsType = 4;
m.matTypes{1}.paramValues = [c0, rho0];
m.matTypes{2}.paramsType = 4;
m.matTypes{2}.paramValues = [cCyl, rhoCyl];

m.elTypes{1}.name = 'AC2D3';

%put on an absorbing boundary
xLims = [-20 -17 17 20]*1e-3;
yLims = [-20 -17 17 20]*1e-3;

m = addAbsBound(m,xLims,yLims,[],[],[],c0,freq);    

m = genPogoHannSignal(m,nCyc,freq,1,1);

m.frames{1}.sigs{1}.nodeSpec = 1;
m.frames{1}.sigs{1}.dofSpec = 1;

m.frames{1}.sigs{1}.sigAmps = 1;
m.frames{1}.sigs{1}.sigType = 0;

m.measFreq = 2;
m.measStart = 1;

m.measSets{1}.name = 'main';
m.measSets{1}.measNodes = 1:nArray;
m.measSets{1}.measDof = ones(1,nArray);

gap = round(m.nt/80);
if gap < 1
    gap = 1;
end

% m.fieldStoreIncs - which increments to output the field at
m.fieldStoreIncs = 1:gap:m.nt;

mblank = m;

% m = deleteEls(m,cylEls);
savePogoInp(sprintf('free2d.pogo-inp'),m);

system(sprintf('pogoBlock free2d'));

system(sprintf('pogoSolve free2d -o'));


mblank.matTypeRefs(((ex).^2 + (ey).^2) < radCyl.^2) = 1;

savePogoInp(sprintf('free2dblank.pogo-inp'),mblank);

system(sprintf('pogoBlock free2dblank'));
system(sprintf('pogoSolve free2dblank -o'));

