% avoid mess
clearvars
home
addpath('../');

% parameters
h1      = 10;
h2      = 7;
qNorth  = -0.0/seconds(years(1));
wMin    = 150;
L       = 1000;
K       = 1e-4;
wMax    = wMin.*1.6;

x = linspace(0,L,21);
y = 0*x;

pltW = 17;
pltH = 7;

fig = figure(1);
clf
fig.Units       = "centimeters";
fig.PaperUnits  = "centimeters";
fig.PaperSize   = [pltW,pltH];
fig.Position    = [0,0,pltW,pltH];
t = tiledlayout('flow');
t.TileSpacing = 'tight';
t.Padding = 'tight';

shapes = ["bump", "composite", "cosinusoidal", "cosinusoidal"];
alphabet = 'abcdefghijklmnopqrstuvwxyz';


for i = 1:numel(shapes)
  mdl = fpAna('h1',h1,'h2',h2,'L',L,'wMin',wMin,...
                'wMax',wMax,'Tx',K,'Ty',K,...
                'qNorth',qNorth,'shape',"bump",'autoSolve',true);
  nexttile()
  mdl.shape = shapes(i);
  if i == numel(shapes)
    mdl.qNorth = -1e-8;
  end
  plot(mdl)
  mdl.stor.go(1).LineWidth = 1.5;
  mdl.stor.go(3).LineWidth = 1.5;
  
  hold on
  % nexttile()
  sc = 0.33*mdl.wMax*mdl.qy(x,y)./max(abs(mdl.qy(x,y)));

  mdl.plot('outline',true,'divide',false)
  
  % plot exchange zone/hillslope zone as area
  pS = isolines(mdl,mdl.psi(0,0),'type',"psi");
  pgon = polyshape(pS.x,pS.y);
  plot(pgon,'FaceColor',[105 186 201]/255,'FaceAlpha',0.4,'EdgeAlpha',0);
  
  if mdl.qNorth~=0
    pS1 = isolines(mdl,mdl.psi(0,mdl.wMin),'type',"psi");
    pS2x = flip(linspace(0,mdl.L,100)');
    pS2y = mdl.fNorth(pS2x);
    pgon2 = polyshape([pS1.x;pS2x],[pS1.y;pS2y],'simplify',false);
    plot(pgon2,'FaceColor',[201 174 105]/255,'FaceAlpha',0.4,'EdgeAlpha',0);
  else
    pgon2 = polyshape([NaN;NaN;NaN],[NaN;NaN;NaN],'simplify',false);
    plot(pgon2,'FaceColor',[201 174 105]/255,'FaceAlpha',0.4,'EdgeAlpha',0);
  end
  
  ylim([-0.33*mdl.wMax mdl.wMax])
  axis off
  
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

%   title(alphabet(i),'Position',[0 mdl.wMin+0.5*(mdl.wMax-mdl.wMin)])
  
  title(sprintf('{\\bf%s}: %s',alphabet(i),shapes(i)),...
      'FontName','Helvetica','FontWeight','Normal','FontSize',10)
% title('this text is {\bfbold} this text is {\ititalizized}','FontWeight','Normal')
end
% 
ax              = axes;
ax.Position     = [0,0,1,1];
ax.Color        = 'None';
ax.XTick        = [0 1];
ax.XColor       = 'white';
ax.YTick        = [0 1];
ax.YColor       = 'white';

exportgraphics(fig,"./figFlowNetExamples.pdf",'ContentType','vector')
