% avoid mess
clearvars
home
addpath('../');

% define parameters for example models
h1      = 10;
h2      = 7;
qNorth  = -0.0/seconds(years(1));
wMin    = [50 150 250 350];
L       = 1000;
K       = 1e-4;
wDiff   = 0.6*150;


%% initialize the figure
pltW            = 17;
pltH            = 7;
fig             = figure(1);
fig             = clf(fig);
fig.Units       = "centimeters";
fig.PaperUnits  = "centimeters";
fig.PaperSize   = [pltW,pltH];
fig.Position    = [0,0,pltW,pltH];
t               = tiledlayout('flow');
t.TileSpacing   = 'tight';
t.Padding       = 'tight';
alphabet        = 'abcdefghijklmnopqrstuvwxyz';

for i = 1:numel(wMin)
  % define the model
  mdl = fpAna('h1',h1,'h2',h2,'L',L,'wMin',wMin(i),...
                'wMax',wMin(i)+wDiff,'Tx',K,'Ty',K,...
                'qNorth',qNorth,'shape',"bump",'autoSolve',true);
  
  % plot the model
  nexttile()
  mdl.plot('outline',true,'divide',false)
  hold on
  
  % adjust line widths
  mdl.stor.go(1).LineWidth = 1.5;
  mdl.stor.go(3).LineWidth = 1.5;
  
  % plot exchange zone/hillslope zone as area
  pS    = isolines(mdl,mdl.psi(0,0),'type',"psi");
  pgon  = polyshape(pS.x,pS.y);
  plot(pgon,'FaceColor',[105 186 201]/255,'FaceAlpha',0.4,'EdgeAlpha',0);
  
  % put area not on top
  h = get(gca,'Children');
    % boundary
    % contourH
    % polygons
    % contourPsi
  h = [h(2) h(3) h(4) h(1) h(5)];
  set(gca,'Children',h);
  
  % add a title
  title(sprintf('{\\bf%s}: %.2f',alphabet(i),mdl.area/mdl.L/mdl.L),...
      'FontName','Helvetica','FontWeight','Normal','FontSize',10)

  ylim([0 440])
  axis off
end

% add dummy axis to make exportgraphics use the full size
ax              = axes;
ax.Position     = [0,0,1,1];
ax.Color        = 'None';
ax.XTick        = [0 1];
ax.XColor       = 'white';
ax.YTick        = [0 1];
ax.YColor       = 'white';

% export the figure
exportgraphics(fig,"./figSeparation.pdf",'ContentType','vector')
