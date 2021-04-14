% avoid mess
clearvars
home
addpath('../');

% constant parameters
h1      = 10;
h2      = 7;
Tx      = 1e-3;
Ty      = 1e-3;
qNorth  = -0.0/seconds(years(1));

% changing parameters
wMinVec     = 50:50:350;
LVec        = 500:250:3500;
rWVec       = 1.0:0.25:2.5;
wMin        = wMinVec;
L           = LVec;
rW          = rWVec;
[wMin,L,rW] = meshgrid(wMin,L,rW);
wMax        = wMin.*rW;

mdl     = fpAna('h1',h1,'h2',h2,'L',L(1),'wMin',wMin(1),...
                'wMax',wMax(1),'Tx',Tx,'Ty',Ty,...
                'qNorth',qNorth,'shape',"cosinusoidal",'autoSolve',false);

Qex     = NaN(size(L));
Atot    = NaN(size(L));
AL2     = NaN(size(L));
for idx = 1:numel(L)
  % convert linear index to subscripts
  [i,j,k] = ind2sub(size(L),idx);
  % modify geometry
  [~]   = set(mdl,'wMin',wMin(i,j,k),'wMax',wMax(i,j,k),'L',L(i,j,k));
  % solve model and store results
  [~]           = mdl.solve();
  Qex(i,j,k)    = mdl.Qex;
  Atot(i,j,k)   = mdl.area();
  AL2(i,j,k)    = L(i,j,k).^2;
end

% create a figure with properties
pltW                    = 17;
pltH                    = 10;
fig                     = figure(1);
fig                     = clf(fig);
fig.Units               = "centimeters";
fig.PaperUnits          = "centimeters";
fig.PaperSize           = [pltW,pltH];
fig.Position            = [0,0,pltW,pltH];
t                       = tiledlayout('flow');
t.TileSpacing           = 'tight';
t.Padding               = 'compact';
fig.Color               = 'White';
set(fig,'DefaultAxesFontSize',12,'DefaultAxesFontName','Crimson Text');

% nexttile([2 1])
nexttile()
scatter3(wMin(:),L(:),rW(:),25,Qex(:),'filled')
hold on
slice(wMin,L,rW,Qex,[],[],rWVec)
pbaspect([0.3 0.5 0.8])
view([-45 20])
shading interp
lighting gouraud
material default
camlight('right','infinite')
camlight('headlight','local')

colormap(colorcet('L17','n',20))
box on
ax          = gca;
ax.BoxStyle = 'full';
grid on
% grid minor
set(gca,'Color',[0.85 0.85 0.85])

xlabel('$w_\mathrm{min}\;\mathrm{in\;m}$','Interpreter','LaTeX');
ylabel('$L\;\mathrm{in\;m}$','Interpreter','LaTeX');
zlabel('$\frac{w_\mathrm{max}}{w_\mathrm{min}}$','Interpreter','LaTeX','Rotation',0,'FontSize',14);
ztickformat('%.1f')

% ylabel(cb,'$\frac{Q_\mathrm{ex}}{m(h_1-h_2)K}$','Interpreter','LaTeX','FontSize',14);
xlim([min(wMin(:)) max(wMin(:))]);
ylim([min(L(:)) max(L(:))]);
zlim([min(rW(:)) max(rW(:))]);
title('{\bfa}: three-dimensional','FontName','Helvetica','FontWeight','normal',...
      'FontSize',10)

nexttile()
% x     = wMin./L.*wMax./L;
x     = (Atot./L.^2).*sqrt(Tx./Ty);
y     = Qex.*L./(wMax-wMin)./(h1-h2)./Tx;

randOrd = randperm(numel(x));
scatter(x(randOrd),y(randOrd),25,Qex(randOrd),'filled')
set(gca,'Color',[0.85 0.85 0.85])

[~]   = xlabel('$\kappa\frac{w_\mathrm{mean}}{L}$','Interpreter','latex','FontSize',14);
[~]   = ylabel('$\frac{Q_\mathrm{ex}\cdot L}{\Delta w \cdot \Delta h \cdot  T_x}$','Interpreter','latex','FontSize',14);

xT = x(:);
yT = y(:);
good = ~isnan(yT) & ~isinf(yT);
modelFun      =  @(p,x) sech(p(1).*x);
startingVals  = 1;
nlModel       = fitnlm(xT(good),yT(good),modelFun,startingVals);
xgrid         = linspace(0,1.1*max(x(:)),75)';
hold on
plt = plot(xgrid,predict(nlModel,xgrid),'-','LineWidth',2,'Color',[0.0 0.0 0.0 0.7]);

title('{\bfb}: transformed','FontName','Helvetica','FontWeight','normal',...
      'FontSize',10)

ytickformat('%.1f')
xtickformat('%.1f')

cb = colorbar('FontSize',12);
cb.Layout.Tile = 'east';
cb.Ruler.TickLabelFormat = '%.1f';
ylabel(cb,'$Q_\mathrm{ex}\;\mathrm{in\;m^3/s}$','Interpreter','LaTeX','FontSize',12);

box on
grid on

ax = axes;
ax.Position = [0,0,1,1];
ax.Color = 'None';
ax.XTick = [0 1];
ax.XColor = 'white';
ax.YTick = [0 1];
ax.YColor = 'white';

exportgraphics(fig,"./figGeometry.pdf",'ContentType','vector')

