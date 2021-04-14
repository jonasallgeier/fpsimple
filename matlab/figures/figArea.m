clearvars
home
addpath('../');

% number of model runs
nM = 100;

cmap = colorcet('L17','n',25);
colormap(cmap(2:end-1,:))
shapes = ["cosinusoidal","bump","composite"];
% shapes = "cosinusoidal";
rmses            = NaN(numel(shapes),1);

alphabet  = 'abcdefghijklmnopqrstuvwxyz';
f = waitbar(0,'Please be patient.','Name','Running simulations...');

% create the halton set
QNmax     = 3;
hSet      = getHaltonSet(nM,'qNstar',[-QNmax 0]); % create the halton set
L         = repmat(hSet.L,2,1);
h1        = repmat(hSet.h1,2,1);
h2        = repmat(hSet.h2,2,1);
Tx        = repmat(hSet.Tx,2,1);
Ty        = repmat(hSet.Ty,2,1);
wMax      = repmat(hSet.wMax,2,1);
wMin      = repmat(hSet.wMin,2,1);
Q0        = repmat(hSet.Q0,2,1);
qNorth    = [0*hSet.qNorth; hSet.qNorth];

pltW                    = 17;
pltH                    = 5;
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

for iS = 1:numel(shapes)
  % pre-allocate some outputs
  Qex     = NaN(numel(L),1);
  Aex     = NaN(numel(L),1);
  Atot    = NaN(numel(L),1);
  AL2     = NaN(numel(L),1);

  for i = 1:numel(L)
    shape = shapes(iS);

    mdl = fpAna('h1',h1(i),'h2',h2(i),'L',L(i),'wMin',wMin(i),...
      'wMax',wMax(i),'Tx',Tx(i),'Ty',Ty(i),...
      'qNorth',qNorth(i),'shape',shape,'autoSolve',false);

    % solve and read results
    [~]       = mdl.solve();
    Atot(i)   = mdl.area();
    AL2(i)    = L(i).^2;
    try
      Qex(i)    = mdl.Qex;
      if Qex(i) > 1e-11
        Aex(i)    = mdl.areaEx();
      else
        Aex(i) = 0;
      end
    catch
    end
    waitbar(((iS-1).*nM+i)/(numel(shapes).*nM),f)
  end

  Abump   = Atot-wMin.*L;
  QexNorm = Qex./Q0;
  QNstar  = abs(qNorth./Q0.*L);
  
  nexttile()
%   x = sqrt(Tx./Ty).*Atot./AL2;
  %   y = sqrt(Kx./Ky).*Atot./AL2;
  %   z = Aex./(wMax-wMin)./L;
  y = Aex./Abump;
  x = QexNorm./(sqrt(1+QNstar));
  
  scatter(x,y,35,QNstar,'filled')
  
  rmses(iS) = sqrt(mean((y-x).^2));
  
  xlabel('$\tilde{Q}_\mathrm{ex}/\sqrt{1+|\tilde{Q}_\mathrm{north}|}$','Interpreter','LaTeX');
  ytickformat('%.1f')
  xtickformat('%.1f')
  
  if iS > 1
    set(gca,'Yticklabel',[])
  end
  
  grid on
  box on
  axis equal
  xlim([0 1]);
  ylim([0 1]);
%   pbaspect([1 1 1]);
  
  title(sprintf('{\\bf%s}: %s',alphabet(iS),shapes(iS)),...
        'FontName','Helvetica','FontWeight','Normal','FontSize',10)
  caxis([0 QNmax])
  set(gca,'Color',[0.85 0.85 0.85])
  

end
% the following line uses an ugly trick (empty second line) to get spacing
% right; cause of problem is "axis equal" with "tiledlayout"
ylabel(t,{'$\tilde{A}$' ''},'Interpreter','LaTeX','FontSize',12);

close(f)

cb = colorbar('FontSize',12);
cb.Layout.Tile = 'east';
cb.Ticks = 0:1:QNmax;
cb.Ruler.TickLabelFormat = '%.1f';
ylabel(cb,'$|\tilde{Q}_\mathrm{north}|$','Interpreter','LaTeX','FontSize',12);

ax = axes;
ax.Position = [0,0,1,1];
ax.Color = 'None';
ax.XTick = [0 1];
ax.XColor = 'white';
ax.YTick = [0 1];
ax.YColor = 'white';
t.Padding = 'loose';
t.Padding = 'tight';
t.TileSpacing = 'compact';
exportgraphics(fig,"./figArea.pdf",'ContentType','vector')