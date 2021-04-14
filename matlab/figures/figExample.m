clearvars
home
addpath('../')

% this script generates a true vector map
% it requires an input map as pdf and geo-referenced png
% procedure is as follows:
%   * read in vector map
%   * do some plotting on top of it
%   * hide the map
%   * save to temporary transparent pdf
%   * use external tool (pdftk) to place temporary pdf on top of existing
%     vector map

% define file locations
mapVec  = './mapArrows.pdf';
mapRast = './mapArrows.png';
mapWrl  = './mapArrows.pgw';
saveTo  = './figExample.pdf';

% initialize the figure
fig             = figure(1);
fig             = clf(fig);
fig.Units       = "centimeters";
fig.PaperUnits  = "centimeters";
fig.PaperSize   = [15,12];
fig.Position    = [0,0,15,12];

% plot the rasterized map in background (will not be saved later!)
[baseImage1,cmap1]  = imread(mapRast);
R1                  = worldfileread(mapWrl,'planar',size(baseImage1));
mp                  = mapshow(baseImage1,cmap1,R1);

% store x and y limits, to reset later on
xlims = xlim();
ylims = ylim();

%     ___                                 
%    /   |  ____ ___  ____ ___  ___  _____
%   / /| | / __ `__ \/ __ `__ \/ _ \/ ___/
%  / ___ |/ / / / / / / / / / /  __/ /    
% /_/  |_/_/ /_/ /_/_/ /_/ /_/\___/_/     
%                                         

% model definition and calculations
exA         = fpAna;
exA.L       = 3000;
exA.wMax    = 600;
exA.wMin    = 175;
exA.h2      = 349; % 9
exA.h1      = 341;
exA.Tx      = 5e-5;
exA.Ty      = 5e-5;
exA.qNorth  = -2.5e-8;
exA.por     = 2*0.1;
exA.shape   = "bump";
[~]         = fprintf('_______AMMER________\n');
[~]         = modelSummary(exA);

% plot the model
hold on
plot(exA,'outline',true,'divide',false)
exA.stor.go(1).LineColor  = [0.6 0.61 0.63];
exA.stor.go(1).LineWidth  = 1.5;

% redefine contours to range within [0,1]
exA.stor.go(2).ZData      = (exA.stor.go(2).ZData-exA.h1)/(exA.h2-exA.h1);
exA.stor.go(2).LevelList  = (exA.stor.go(2).LevelList-exA.h1)/(exA.h2-exA.h1);

% plot translucent areas
pS        = isolines(exA,exA.psi(exA.L,0),'type',"psi");
pgon      = polyshape(pS.x,pS.y);
pg        = plot(pgon,'FaceColor',[105 186 201]/255,...
                      'FaceAlpha',0.7,...
                      'EdgeAlpha',0);
pS1       = isolines(exA,exA.psi(exA.L,exA.wMin),'type',"psi");
pS2x      = flip(linspace(0,exA.L,100)');
pS2y      = exA.fNorth(pS2x);
pgon2     = polyshape([0;pS1.x;pS2x],[0;pS1.y;pS2y],'simplify',false);
pg2       = plot(pgon2,'FaceColor',[201 174 105]/255,...
                       'FaceAlpha',0.7,...
                       'EdgeAlpha',0);

% group graphics objects to move and rotate them as one
h         = [flip(exA.stor.go) pg pg2];
ax        = gca;
t1        = hgtransform('Parent',ax);
set(h,'Parent',t1);

% center, rotate and translate the groupt to coincide with the map
matCen      = makehgtform('translate',[-0.5*exA.L -0.5*exA.wMax 0]);
matRot      = makehgtform('zrotate',0.925*pi);
matTra      = makehgtform('translate',[R1.firstCornerX+4800 R1.firstCornerY-1600 0]);
t1.Matrix   = matTra * matRot * matCen;
t1.Children = h;

%     _   __          __             
%    / | / /__  _____/ /______ ______
%   /  |/ / _ \/ ___/ //_/ __ `/ ___/
%  / /|  /  __/ /__/ ,< / /_/ / /    
% /_/ |_/\___/\___/_/|_|\__,_/_/     
%                            

% model definition and calculations       
exN         = fpAna();
exN.L       = 6500;
exN.wMax    = 1750;
exN.wMin    = 500;
exN.h1      = 345;
exN.h2      = 324;
exN.Tx      = 1.25e-2;
exN.Ty      = 1.25e-2;
Q0          = (exN.h1-exN.h2)/exN.L*exN.Tx*(exN.wMax-exN.wMin);
exN.qNorth  = -7.5e-7;
exN.por     = 5*0.15;
exN.shape   = "cosinusoidal";
[~]         = fprintf('______NECKAR________\n');
[~]         = modelSummary(exN);

% plot the model
hold on
plot(exN,'outline',true,'divide',false)
exN.stor.go(1).LineColor = [0.6 0.61 0.63];
exN.stor.go(1).LineWidth = 1.5;

% redefine contours to range within [0,1]
exN.stor.go(2).ZData      = (exN.stor.go(2).ZData-exN.h2)/(exN.h1-exN.h2);
exN.stor.go(2).LevelList  = (exN.stor.go(2).LevelList-exN.h2)/(exN.h1-exN.h2);

% plot translucent areas
pS        = isolines(exN,exN.psi(0,0),'type',"psi");
pgon      = polyshape(pS.x,pS.y);
pg        = plot(pgon,'FaceColor',[105 186 201]/255,...
                      'FaceAlpha',0.7,...
                      'EdgeAlpha',0);
pS1       = isolines(exN,exN.psi(0,exN.wMin),'type',"psi");
pS2x      = flip(linspace(0,exN.L,100)');
pS2y      = exN.fNorth(pS2x);
pgon2     = polyshape([pS1.x;pS2x],[pS1.y;pS2y],'simplify',false);
pg2       = plot(pgon2,'FaceColor',[201 174 105]/255,...
                       'FaceAlpha',0.7,...
                       'EdgeAlpha',0);

% group graphics objects to move and rotate them as one
h         = [flip(exN.stor.go) pg pg2];
t1        = hgtransform('Parent',ax);
set(h,'Parent',t1)

% group graphics objects to move and rotate them as one
matCen      = makehgtform('translate',[-0.5*exN.L -0.5*exN.wMax 0]);
matRot      = makehgtform('zrotate',0.15*pi);
matTra      = makehgtform('translate',[R1.firstCornerX+4200 R1.firstCornerY-4550 0]);
t1.Matrix   = matTra * matRot * matCen;
t1.Children = h;


%     _______       _      __    _            
%    / ____(_)___  (_)____/ /_  (_)___  ____ _
%   / /_  / / __ \/ / ___/ __ \/ / __ \/ __ `/
%  / __/ / / / / / (__  ) / / / / / / / /_/ / 
% /_/   /_/_/ /_/_/____/_/ /_/_/_/ /_/\__, /  
%                                    /____/   

% reset the axis limits to the map extent
xlim(xlims)
ylim(ylims)
set(gca,'Position',[0 0 1 1]);

% hide the map (because we want to use the shiny vector map!)
mp.Visible = 'off';

% set color range from 0 to 1
caxis([0 1])
axis off

% add new empty axis spanning all the figure --> use full paper for print
ax              = axes;
ax.Position     = [0,0,1,1];
ax.Color        = 'None';
ax.XTick        = [0 1];
ax.XColor       = 'white';
ax.YTick        = [0 1];
ax.YColor       = 'white';

% export to temporary pdf
tempFig = fullfile(tempdir,'tempFig.pdf');
exportgraphics(fig,tempFig,'BackgroundColor','none','ContentType','vector')

% use systems pdftk to overlap with existing vector map
system(sprintf('pdftk %s background %s output %s',tempFig,mapVec,saveTo));
% exportgraphics(fig,"../tex/figExample.pdf",'BackgroundColor','none','ContentType','vector')