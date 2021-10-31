%%
close all
clear
clc


%%
FDTD = InitFDTD('NrTS',100, 'EndCriteria',0, 'OverSampling',50);
FDTD = SetSinusExcite(FDTD,10e6);
FDTD = SetBoundaryCond(FDTD,{'PMC' 'PMC' 'PEC' 'PEC' 'MUR' 'MUR'});

%%
CSX = InitCSX();

%%
mesh.x = -10:10;
mesh.y = -10:10;
mesh.z = -10:30;
CSX = DefineRectGrid(CSX, 1, mesh);

%%
CSX = AddExcitation(CSX,'excitation',0,[0 1 0]);
CSX = AddBox(CSX,'excitation',0,[-10 -10 0],[10 10 0]);

%%

CSX = AddDump(CSX,'Et');
CSX = AddBox(CSX,'Et',0,[-10 0 -10],[10 0 30]);

%%
mkdir('tmp')
WriteOpenEMS('tmp/tmp.xml',FDTD,CSX);

%%
CSXGeomPlot( 'tmp/tmp.xml' );


%%
RunOpenEMS('tmp','tmp.xml','');