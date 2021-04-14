% avoid mess
clearvars
home
addpath('../');

% define model parameters
h1                = 10;
h2                = 7;
Tx                = 1e-3;
Ty                = 1e-3;
wMin              = 50:50:350;
L                 = 500:500:3500;
RstarO            = -linspace(0,3,15);
rW                = 1.0:0.25:2.5;
[wMin,L,rW,Rstar] = ndgrid(wMin,L,rW,RstarO);
wMax              = wMin.*rW;
shape             = "cosinusoidal";

% determine Q0 for all realizations
Q0              = (h1-h2)./L.*Tx.*(wMax-wMin);

% create initial model
mdl = fpAna('h1',h1,'h2',h2,'L',L(1),'wMin',wMin(1),...
          'wMax',wMax(1),'Tx',Tx,'Ty',Ty,...
          'qNorth',Rstar(1).*Q0(1)./L(1),'shape',shape,'autoSolve',false);

% initialize output
Qex  = NaN(size(L));
Atot = NaN(size(L));
AL2  = NaN(size(L));
% iterate through all model realizations
for i1 = 1:size(L,1)
  for i2 = 1:size(wMin,2)
    for i3 = 1:size(wMax,3)
      % modify geometry
      [~]   = set(mdl,'wMin',wMin(i1,i2,i3,1),'wMax',wMax(i1,i2,i3,1),'L',L(i1,i2,i3,1));
    
      for i4 = 1:size(Rstar,4)
        % modify Ky
        [~]                 = set(mdl,'qNorth',Rstar(i1,i2,i3,i4).*Q0(i1,i2,i3,i4)./L(i1,i2,i3,i4));
        % solve and read results
        [~]                 = mdl.solve();
        Qex(i1,i2,i3,i4)    = mdl.Qex;
        Atot(i1,i2,i3,i4)   = mdl.area();
        AL2(i1,i2,i3,i4)    = L(i1,i2,i3,i4).^2;
      end
    end
  end
end

% prepare parameters for plotting
clrs  = colorcet('L17','n',20);
clrs  = clrs(2:end-1,:);
colormap(clrs)

% initialize figure
fig             = figure(1);
fig             = clf(fig);
pltW            = 8;
pltH            = 15;
fig.Units       = "centimeters";
fig.PaperUnits  = "centimeters";
fig.PaperSize   = [pltW,pltH];
fig.Position    = [0,0,pltW,pltH];
t               = tiledlayout('flow');
fig.Color       = 'White';
set(gcf,'DefaultAxesFontSize',12,'DefaultAxesFontName', 'Crimson Text')
nexttile()
alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ';

% create 3D scatter variables
x   = Atot./L.^2.*sqrt(Tx./Ty);
y   = abs(Rstar);
z   = Qex./Q0;

% create 3D scatter plot & rotate to perpendicular view
sctr3           = scatter3(x(:),y(:),z(:),15,y(:),'filled');
xlabel('x')
ylabel('y')
zlabel('z')
view([0 0])
set(gca,'Color',[0.85 0.85 0.85])
xlim([0 1.5])
zlim([0 1])

zlabel('$\tilde{Q}_\mathrm{ex}$','Interpreter','LaTeX');
xlabel('$\tilde{x}$','Interpreter','LaTeX');

title(alphabet(1),'FontName','Helvetica')

z(z==0) = NaN;
z = z./z(:,:,:,1);
z = 1-z;
z = z./y; %.^0.925;

% new subfigure/tile for the fitting of transformed values
nexttile()
sc = scatter(x(:),z(:),36,y(:),'filled');
colormap(clrs)
set(gca,'Color',[0.85 0.85 0.85])
grid on

ylabel('$\frac{1}{|\tilde{Q}_\mathrm{north}|}\left(1-\frac{\tilde{Q}_\mathrm{ex}}{\tilde{Q}_\mathrm{ex}^0}\right)$','Interpreter','LaTeX');
xlabel('$\tilde{x}$','Interpreter','LaTeX');
  
% fit a model through the points
modelB        =  @(p,x) p(1)*cosh(p(2)*x);
startingVals  = [1.0 1.0];
nlModel       = fitnlm(x(:),z(:),modelB,startingVals);
pB            = nlModel.Coefficients.Estimate';
fitB          = nlModel.Rsquared.Ordinary;

% add a fitted line
xgrid = linspace(0,1,50);
ygrid = predict(nlModel,xgrid');
hold on
plot(xgrid,ygrid,'-','LineWidth',2,'Color',[0.0 0.0 0.0 0.7])
ylim([0 5])
xlim([0 1])

% add colorbar
cb = colorbar('southoutside','FontSize',12);
cb.Ticks = 0:1:max(abs(RstarO));
cb.Ruler.TickLabelFormat = '%.1f';
ylabel(cb,'$|\tilde{Q}_\mathrm{north}|$','Interpreter','LaTeX','FontSize',12);

% add title
title(alphabet(2),'FontName','Helvetica')

% dummy axis, to make exportgraphics export everything
ax = axes;
ax.Position = [0,0,1,1];
ax.Color = 'None';
ax.XTick = [0 1];
ax.XColor = 'white';
ax.YTick = [0 1];
ax.YColor = 'white';

% export the result
t.Padding = 'loose';
t.Padding = 'compact';
t.TileSpacing = 'loose';
exportgraphics(fig,"./figSuppQnorth.pdf",'ContentType','vector')
