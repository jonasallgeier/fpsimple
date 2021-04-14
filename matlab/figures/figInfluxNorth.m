clearvars
home
addpath('../');

nM              = 1500; % number of model runs

fig             = figure(1);
fig             = clf(fig);
pltW            = 17;
pltH            = 6;
fig.Units       = "centimeters";
fig.PaperUnits  = "centimeters";
fig.PaperSize   = [pltW,pltH];
fig.Position    = [0,0,pltW,pltH];
t               = tiledlayout(1,3);
fig.Color       = 'White';
set(gcf,'DefaultAxesFontSize',12,'DefaultAxesFontName', 'Crimson Text')
       
colormap(colorcet('L17','n',20))

shapes          = ["cosinusoidal","bump","composite"];
colNames        = ["a1", "a1u", "a2", "a2u", "a3",  "a3u", "aRsq","bRsq"];
data            = NaN(numel(shapes),numel(colNames));

alphabet  = 'abcdefghijklmnopqrstuvwxyz';
f         = waitbar(0,'Please be patient.','Name','Running simulations...');
QNmax     = 3;
hSet      = getHaltonSet(nM,'qNstar',[-QNmax 0]); % create the halton set
L         = hSet.L;
h1        = hSet.h1;
h2        = hSet.h2;
Tx        = hSet.Tx;
Ty        = hSet.Ty;
wMax      = hSet.wMax;
wMin      = hSet.wMin;
Q0        = hSet.Q0;
qNorth    = hSet.qNorth;

for iS = 1:numel(shapes)
  % pre-allocate some outputs
  Qex  = NaN(numel(L),2);
  Atot = NaN(numel(L),1);
  AL2  = NaN(numel(L),1);
  
  for iR = 1:2
    if iR == 1
      rVal = 0.*qNorth;
    else
      rVal = qNorth;
    end
    
    for i = 1:nM
      shape = shapes(iS);
      mdl = fpAna('h1',h1(i),'h2',h2(i),'L',L(i),'wMin',wMin(i),...
        'wMax',wMax(i),'Tx',Tx(i),'Ty',Ty(i),...
        'qNorth',rVal(i),'shape',shape,'autoSolve',false);

      % solve and read results
      [~]       = mdl.solve();
      Qex(i,iR) = mdl.Qex;
      Atot(i)   = mdl.area();
      AL2(i)    = L(i).^2;

      waitbar(((iS-1).*nM+i)/(numel(shapes).*nM),f)
    end
  end
  
  y = abs(qNorth./Q0.*L);
  z = Qex./Q0;
  x = sqrt(Tx./Ty).*Atot./AL2;

  % fit the first run (r=0)
  zM            = z(:,1);
  good          = ~isnan(zM) & ~isinf(zM);
  modelA        =  @(p,x) sech(p(1).*x);
  startingVals  = 1;
  nlModA        = fitnlm(x(good),zM(good),modelA,startingVals);
  pA            = nlModA.Coefficients.Estimate';
  pAu           = errorProp(x(good),zM(good),modelA,pA);
  fitA          = nlModA.RMSE;
  
  % now fit the second run (r~=0)
  zM            = 1-z(:,2)./z(:,1);
  zM(z(:,2)==0 | z(:,1)==0) = NaN;
  zM            = zM./y;
  
  modelB        =  @(p,x) x(:,3).*max(0,1- p(1)*x(:,2).*cosh(p(2)*x(:,1)) );
  startingVals  = [0.37 3.9];
  mBxIn         = [x(:),y(:),modelA(pA,x(:))];
  nlModB        = fitnlm(mBxIn,z(:,2),modelB,startingVals); %,'Weights',w);
  pB            = nlModB.Coefficients.Estimate';
  pBu           = errorProp(mBxIn,z(:,2),modelB,pB);
  fitB          = nlModB.RMSE;
 
  % plot the second run (first plot all with qNorth == 0; then the rest)
  [~,v] = sort(y,'ascend');
  v1    = flip(v(ismember(v,find(z(:,2)==0))));
  v2    = v(ismember(v,find(z(:,2)~=0)));
  nexttile()
  
  scatter(x(v1),z(v1,2),25,y(v1),'filled')
  hold on
  scatter(x(v1),z(v1,1),25,0*y(v1),'filled')
  scatter(x(v2),z(v2,2),25,y(v2),'filled')
  scatter(x(v2),z(v2,1),25,0*y(v2),'filled')
  caxis([0 QNmax])
  set(gca,'Color',[0.85 0.85 0.85])
  
  title(sprintf('{\\bf%s}: %s',alphabet(iS),shapes(iS)),...
        'FontName','Helvetica','FontWeight','Normal','FontSize',10)
  xlabel('$\tilde{x}$','Interpreter','LaTeX');
  ytickformat('%.2f')
  xtickformat('%.1f')

  grid on
  box on
  
  xgrid   = linspace(0,1.5,75)';
  QnorthIso = linspace(0,QNmax,7);
  Q1 = NaN(numel(xgrid),numel(QnorthIso));
  for i = 1:numel(QnorthIso)
    Q1(:,i)      = modelB(pB,[xgrid,QnorthIso(i).*ones(size(xgrid)), modelA(pA,xgrid)]);
  end
  plt     = plot(xgrid,Q1,'-','LineWidth',2,'Color',[0.0 0.0 0.0 0.7]);
  
  xlim([0 1.5])
  ylim([0 1])
  
  if iS > 1
    set(gca,'Yticklabel',[])
  end
  % determine total Rsquare
  
%   QtilMod  = modelB(pB,[x,y, modelA(pA,x)]);
%   Qtil     = z(:,2);
%   fitTotal = sqrt(sum((Qtil-QtilMod).^2)/(numel(Qtil)-numel([pB,pA])));
  
  % store the fit results
  data(iS,:)    = [pA, pAu, pB(1), pBu(1), pB(2), pBu(2), fitA, fitB];
end
ylabel(t,'$\tilde{Q}_\mathrm{ex}$','Interpreter','LaTeX');
close(f)

cb = colorbar('FontSize',12);
cb.Layout.Tile = 'east';
% cb = colorbar('Location','layout','FontSize',12);
% cb.Layout.Tile = 7;
% cb.Layout.Span = [1 3];
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
exportgraphics(fig,"./figInfluxNorth.pdf",'ContentType','vector')
% results = array2table(round(data,2),'RowNames',shapes,'VariableNames',colNames);
% disp(results)
% writetable(results,'fittedCoeffs.txt','Delimiter',' ','WriteRowNames',true);
