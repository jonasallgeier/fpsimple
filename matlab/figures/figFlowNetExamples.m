% avoid mess
clearvars
home
addpath('../');

% define parameters for example models
h1      = 10;
h2      = 7;
qNorth  = -0.0/seconds(years(1));
wMin    = 150;
L       = 1000;
K       = 1e-4;
wMax    = wMin.*1.6;

% define arrow positions
x = linspace(0,L,21);
y = 0*x;

% initialize the figure
pltW = 17;
pltH = 7;
fig = figure(1);
fig = clf(fig);
fig.Units       = "centimeters";
fig.PaperUnits  = "centimeters";
fig.PaperSize   = [pltW,pltH];
fig.Position    = [0,0,pltW,pltH];
t = tiledlayout('flow');
t.TileSpacing = 'tight';
t.Padding = 'tight';
alphabet = 'abcdefghijklmnopqrstuvwxyz';

% iterate through shape types
shapes = ["bump", "composite", "cosinusoidal", "cosinusoidal"];
for i = 1:numel(shapes)
  % define the model
  mdl = fpAna('h1',h1,'h2',h2,'L',L,'wMin',wMin,...
                'wMax',wMax,'Tx',K,'Ty',K,...
                'qNorth',qNorth,'shape',"bump",'autoSolve',true);
  mdl.shape = shapes(i);
  
  % in the last case: add a northern influx
  if i == numel(shapes)
    mdl.qNorth = -1e-8;
  end
  
  % plot the model
  nexttile()
  mdl.plot('outline',true,'divide',false)
  hold on
  
  % adjust line widths
  mdl.stor.go(1).LineWidth = 1.5;
  mdl.stor.go(3).LineWidth = 1.5;
  
  % plot exchange zone/hillslope zone as area
  pS = isolines(mdl,mdl.psi(0,0),'type',"psi");
  pgon = polyshape(pS.x,pS.y);
  plot(pgon,'FaceColor',[105 186 201]/255,'FaceAlpha',0.4,'EdgeAlpha',0);
  
  % if there is northern influx: show the area
  if mdl.qNorth~=0
    pS1 = isolines(mdl,mdl.psi(0,mdl.wMin),'type',"psi");
    pS2x = flip(linspace(0,mdl.L,100)');
    pS2y = mdl.fNorth(pS2x);
    pgon2 = polyshape([pS1.x;pS2x],[pS1.y;pS2y],'simplify',false);
    plot(pgon2,'FaceColor',[201 174 105]/255,'FaceAlpha',0.4,'EdgeAlpha',0);
  else
    % add dummy plot to comply with graphics object order
    plot(NaN,NaN);
  end
  
  % set ylim to be able to show arrows
  ylim([-0.33*mdl.wMax mdl.wMax])
  axis off
  
  % add arrows
  sc = 0.33*mdl.wMax*mdl.qy(x,y)./max(abs(mdl.qy(x,y))); % <-- scaling
  sz = 5*abs(sc)./max(abs(sc));
  for ii = 1:numel(x)
    xAr = [x(ii) x(ii)];
    yAr = [y(ii) y(ii)+sc(ii)];
    if mdl.psi(x(ii),0) < min(mdl.psi([0 mdl.L],[0 0]))
      clr = 'k';
    else
      clr = [0.75 0.75 0.75];
    end
    ah = annotation('arrow',[0 1],[0 1],'Color',clr,'headStyle','cback2',...
                    'HeadLength',sz(ii),'HeadWidth',sz(ii),'LineWidth',1,'LineStyle','none');
    set(ah,'parent',gca);
    set(ah,'position',[x(ii) y(ii) 0 sc(ii)]);
    % we draw arrow stems by ourselves to avoid ugly overlapping of stem
    % and head
    plot([x(ii) x(ii)],[0 sign(sc(ii))*(abs(sc(ii))-4*sz(ii))],'Color',clr,'LineWidth',(sz(ii)+1)/3);
  end

  % change order such that arrows are not on top
  h = get(gca,'Children');
    % boundary
    % arrows
    % contourH
    % polygons
    % contourPsi
  h = [h(end-3); h(1:end-6); h(end-1); h(end-4); h(end-5); h(end); h(end-2);];
  set(gca,'Children',h);

  % add a title
  title(sprintf('{\\bf%s}: %s',alphabet(i),shapes(i)),...
      'FontName','Helvetica','FontWeight','Normal','FontSize',10)
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
exportgraphics(fig,"./figFlowNetExamples.pdf",'ContentType','vector')
