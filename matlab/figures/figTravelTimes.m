clearvars
home
addpath('../');

% number of model runs
nM = 1500;

% number of points on travel time distribution
nTT = 100;

% get a color map
cmap = colorcet('L17','n',25);
colormap(cmap(2:end-1,:))
shapes = ["cosinusoidal","bump","composite"];

% define subfigure labels
alphabet = 'abcdefghijklmnopqrstuvwxyz';

% create the halton set
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

% do everything twice: first time without northern influx; second time with it
for iFig = 1:2
  if iFig ==2
    qNorth    = 0*hSet.qNorth;
  end
  
  % initialize the figure
  fig             = figure(1);
  fig             = clf(fig);
  pltW            = 17.15;
  pltH            = 20;
  fig.Units       = "centimeters";
  fig.PaperUnits  = "centimeters";
  fig.PaperSize   = [pltW,pltH];
  fig.Position    = [0,0,pltW,pltH];
  t               = tiledlayout('flow');
  fig.Color       = 'White';
  set(fig,'DefaultAxesFontSize',12,'DefaultAxesFontName', 'Crimson Text');

  % request travel time distribution values within the interval (0,2)
  TT = linspace(0,2,nTT);
  
  % iterate through all shapes
  f = waitbar(0,'Please be patient.','Name','Running simulations...');
  for iS = 1:numel(shapes)
    % pre-allocate some outputs
    Qex     = NaN(numel(L),1);
    Atot    = NaN(numel(L),1);
    AL2     = NaN(numel(L),1);
    FF      = NaN(numel(L),nTT);
    
    % iterate through all realizations
    for i = 1:nM
      % create the model
      shape = shapes(iS);
      mdl   = fpAna('h1',h1(i),'h2',h2(i),'L',L(i),'wMin',wMin(i),...
                    'wMax',wMax(i),'Tx',Tx(i),'Ty',Ty(i),...
                    'qNorth',qNorth(i),'shape',shape,'autoSolve',false);
      
      % solve and read results
      [~]       = mdl.solve();
      Qex(i)    = mdl.Qex;
      Atot(i)   = mdl.area();
      AL2(i)    = L(i).^2;
      
      % retrieve travel time distribution of exchange zone (if it exists)
      try
        if abs(mdl.Qex) > 1e-12 && mdl.areaEx/mdl.area > 1e-3
          [f2,x2] = mdl.ttd("k",100,"normalize","mean",'xq',TT);
          FF(i,:) = f2;
        end
      catch
      end
      
      % show progress with waitbar
      waitbar(((iS-1).*nM+i)/(numel(shapes).*nM),f)
    end
    
    % first make a plot of all cdfs
    nexttile()
    plot(TT,FF','Color',[0.6 0.6 0.6 max(1/nM,0.1)])
    hold on
    % add median and 5th/95th percentiles
    pctls   = prctile(FF,[5 95]);
    mdn     = prctile(FF,50);
    p1 = plot(TT,mdn,'Color',[0.86 0.13 0.32],'LineWidth',2);
    p2 = plot(TT,pctls,'k','LineWidth',2);
    
    % add legend for first shape type
    if iS ==1
      legend([p1,p2(1)],{'$50^\mathrm{th}$','$5^\mathrm{th}/95^\mathrm{th}$'},'Location','Northwest','Interpreter','LaTeX');
    end
    
    % now fit beta distributions, and store fitted coefficients
    pFit    = NaN(nM,3);
    for i = 1:nM
      xS = TT;
      yS = FF(i,:);
      if any(~isnan(yS))
        model = @(p,x) betacdf(min(max((x-0)./(p(3)-0),0),1),p(1),p(2));
        startingVals  = [2.0 2.0 2.0];
        nlModel       = fitnlm(xS(:),yS(:),model,startingVals);
        pFit(i,:)     = nlModel.Coefficients.Estimate';
        
        % yFit = predict(nlModel,xS');
        % Qvec = repmat(QexNorm(i),1,numel(xS));
        % plot3(xS,Qvec,yFit,'k-','LineWidth',2.5)
      end
    end
    
    % add title/ticks etc.
    title(sprintf('{\\bf%s}: %s',alphabet(2*(iS-1)+1),shapes(iS)),...
        'FontName','Helvetica','FontWeight','Normal','FontSize',10)
    ylabel('$F(\tilde{t})$','Interpreter','LaTeX');
    ytickformat('%.1f')
    xtickformat('%.1f')
    if iS == numel(shapes)
      xlabel('$\tilde{t}$','Interpreter','LaTeX');
    else
      set(gca,'Xticklabel',[])
    end
    grid on
    box on
    
    % secondly, plot how the fitted parameters relate to tildeQex
    nexttile();
    tildeQex = Qex./Q0;
    [~,ord] = sort(tildeQex);
    s1 = scatter(tildeQex(ord),pFit(ord,1),15,[102,194,165]/255,'filled');
    hold on
    s2 = scatter(tildeQex(ord),pFit(ord,2),15,[252,141,98]/255,'filled');
    s3 = scatter(tildeQex(ord),pFit(ord,3),15,[141,160,203]/255,'filled');
    ylims = ylim();
    ylim([0 6])
    xlim([0 1])
    
    % add legend on first subfigure
    if iS ==1
      legend([s1,s2,s3],{'$\alpha$','$\beta$','$\tilde{t}_\mathrm{max}$'},'Location','Northwest','Interpreter','LaTeX');
    end
    
    % add title, ticks etc.
    title(sprintf('{\\bf%s}: %s',alphabet(2*iS),shapes(iS)),...
        'FontName','Helvetica','FontWeight','Normal','FontSize',10)
    xtickformat('%.1f')
    if iS == numel(shapes)
      xlabel('$\tilde{Q}_\mathrm{ex}$','Interpreter','LaTeX');
    else
      set(gca,'Xticklabel',[])
    end
    ylabel('$\alpha,\beta,\tilde{t}_\mathrm{max}$','Interpreter','LaTeX');
    grid on
    box on
  end
  close(f)
  
  % adjust tile spacing
  t.Padding = 'compact';
  t.TileSpacing = 'compact';
  
  % disable legend updates % hide legends
  set(findall(gcf,'type','Legend'),'AutoUpdate','off')
  set(findall(gcf,'type','Legend'),'Visible','off')
  
  % to reduce the file size of the output, we replace the vectorized data
  % points of all subfigures, with raster images of the data at the exact
  % same place. this way, we can achieve real vector output (e.g. for the
  % axes, labels etc) while having small files.
  
  % iterate through tiles
  noTiles = numel(findobj( get(t,'Children'), '-depth', 1, 'type', 'axes'));
  for i = 1:noTiles
    % create a second, temporary figure with another axis to create the raster image
    tempFig = figure();
    tempFig = clf(tempFig);
    ax2 = axes(tempFig);
    
    % in the original figure, select the current tile
    figure(1)
    ax = nexttile(i);

    % change and store various settings, to prevent auto-changes of axis
    clX = ax.XColor;
    clY = ax.YColor;
    clT = ax.Title.Color;
    ax.GridColor = [0.15 0.15 0.15];
    ax.XAxis.LimitsMode = 'manual';
    ax.XAxis.MinorTickValuesMode = 'manual';
    ax.XAxis.TickDirectionMode = 'manual';
    ax.XAxis.TickLabelsMode = 'manual';
    ax.XAxis.TickLabelRotationMode= 'manual';
    ax.XAxis.TickValuesMode = 'manual';
    
    ax.YAxis.LimitsMode = 'manual';
    ax.YAxis.MinorTickValuesMode = 'manual';
    ax.YAxis.TickDirectionMode = 'manual';
    ax.YAxis.TickLabelsMode = 'manual';
    ax.YAxis.TickLabelRotationMode= 'manual';
    ax.YAxis.TickValuesMode = 'manual';
    
    axCh = ax.Children;
    % copy all the plot data from original axis, to secondary figure
    copyobj(axCh,ax2);
    
    % store limits, aspect ratio, position
    ar = pbaspect(ax);
    ax.Units = 'centimeters';
    pos = ax.Position;
    xlims = xlim();
    ylims = ylim();
    
    % set temporary figure to same size as original tile
    tempFig.Units = 'centimeters';
    tempFig.Position = [pos(1) pos(2) pos(3) pos(3)./ar(1).*ar(2)];
    drawnow
    
    % set temporary axis properties
    ax2.Units = 'normalized';
    ax2.Position = [0,0,1,1];
    ax2.Color = 'None';
    ax2.XColor = 'black';
    ax2.YColor = 'black';
    ax2.XTick = ax.XTick;
    ax2.YTick = ax.YTick;
    xlim(ax2,xlims);
    ylim(ax2,ylims);
    grid(ax2,'on');
    box(ax2,'on');
    ax2.GridColor = [0.15 0.15 0.15];
    ax2.Units = 'centimeters';
    
    % create a temporary file with the raster image
    filename = fullfile(tempdir,sprintf('temp_%d_%d.png',iFig,i));
    exportgraphics(ax2,filename,'Resolution',300);
    
    % clear original axis data
    cla(ax);

    % read the image and place it where the previous data was
    I = imread(filename);
    DX = xlims(2)-xlims(1);
    DY = ylims(2)-ylims(1);
    nx = size(I,1);
    ny = size(I,2);
    
    xP = [xlims(1)+DX/(nx*2) xlims(2)-DX/(nx*2)];
    yP = [ylims(2)-DY/(ny*2) ylims(1)+DY/(ny*2)];
    
    img = image(ax,'XData',xP,'YData',yP,'CData',I);
    xlim(xlims)
    ylim(ylims)
    
    % now adjust/recreate some details
    grid(ax,'off');
    box(ax,'off');
    ax.YColor       = 'white';
    ax.YLabel.Color = 'black';
    ax.XColor       = 'white';
    ax.XLabel.Color = 'black';
    
    % recreate ticklabels for y axis
    ticklabels      = get(ax,'YTickLabel');
    ticklabels_new  = cell(size(ticklabels));
    for ii = 1:length(ticklabels)
      ticklabels_new{ii} = ['\color{black} ' ticklabels{ii}];
    end
    set(ax, 'YTickLabel', ticklabels_new);
    
    % recreate ticklabels for x axis
    ticklabels = get(ax,'XTickLabel');
    ticklabels_new = cell(size(ticklabels));
    for ii = 1:length(ticklabels)
      ticklabels_new{ii} = ['\color{black} ' ticklabels{ii}];
    end
    set(ax, 'XTickLabel', ticklabels_new);
    
    % delete the temporary figure
    close(tempFig)
  end
  
  % re-enable the legends
  set(findall(gcf,'type','Legend'),'Visible','on')
  drawnow
  
  % add dummy axis, to ensure that exportgraphics exports full figure
  ax = axes;
  ax.Position = [0,0,1,1];
  ax.Color = 'None';
  ax.XTick = [0 1];
  ax.XColor = 'white';
  ax.YTick = [0 1];
  ax.YColor = 'white';
  
  % export the figure
  if iFig==1
    exportgraphics(fig,"./figTravelTimesMean.pdf",'ContentType','vector')
  else
    exportgraphics(fig,"./figTravelTimesQ0Mean.pdf",'ContentType','vector')
  end
end