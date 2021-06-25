clearvars
clc
clf
tic

% do the analysis for all 3 shapes
shapes    = ["cosinusoidal" "bump" "composite"];
sobl      = NaN(7,3,2,2); % nPar, nShapes, nTargets, nResults
nM        = 1500;

% gather all the data we need for the plot
for iS = 1:numel(shapes)
  shape = shapes(iS);
  
  filename = sprintf('datasetSobol_%s.mat',shape);
    
  % load the ensemble data if present; if not: create & load it
  if ~exist(filename,'file')
    createDataSet(nM,shape);
  end
  load(filename,'A','B','AB',...
                'fAflux','fBflux','fABflux',...
                'fAarea','fBarea','fABarea');

  % determine sensitivities for exchange flux and exchange area
  for j = 1:2
    if j == 1
      fB  = fBflux;
      fA  = fAflux;
      fAB = fABflux;
    elseif j== 2
      fB  = fBarea;
      fA  = fAarea;
      fAB = fABarea;      
    else
      error('You should not end up here. Check code!')
    end
    
    % number of parameters / dimension
    d = width(A);
    % number of realizations
    N = height(A);

    for i = 1:d
      varT  = 1/N * sum(  fB.* ( fAB(:,i)  - fA ));
      EX(i)    = 1/(2*N) * sum((fA-fAB(:,i)).^2); %#ok<SAGROW>
      V(i) = varT; %#ok<SAGROW>
    end

    vTot = var([fA; fB]);

    % first order Sobol indices
    S = V./vTot;
    % total effect indices
    ST = EX/vTot;
    
    sobl(:,iS,j,1) = S;
    sobl(:,iS,j,2) = ST;
  end
end

% create the plot

% set up the figure
pltW                    = 17;
pltH                    = 17;
fig                     = figure(1);
fig                     = clf(fig);
fig.Units               = "centimeters";
fig.PaperUnits          = "centimeters";
fig.PaperSize           = [pltW,pltH];
fig.Position            = [0,0,pltW,pltH];
t                       = tiledlayout(2,2);
t.TileSpacing           = 'compact';
t.Padding               = 'compact';
fig.Color               = 'White';
alphabet                = 'abcdefghijklmnopqrstuvwxyz';
set(fig,'DefaultAxesFontSize',12,'DefaultAxesFontName','Crimson Text');

% redefine labels
labels = string(A.Properties.VariableNames);
labels(labels=="I") = "$\Delta{}h/L$";
labels(labels=="Tmean") = "$\sqrt{T_x T_Y}$";
labels(labels=="aniso") = "$T_x/T_Y$";
labels(labels=="L") = "$L$";
labels(labels=="ratMin") = "$w_\mathrm{min}/w_\mathrm{max}$";
labels(labels=="ratMax") = "$w_\mathrm{max}/L$";
labels(labels=="QNstar") = "$\tilde{Q}_\mathrm{north}$";

% define bar colors
cmap        = [128,177,211;
              255,237,111;
              179,222,105]/255;
cmap2       = brighten(cmap,-0.7);

% iterate through all four subplot and create them
ctr = 1;
for i = 1:2
  for j = 1:2
    nexttile()
    % create bar chart
    X = categorical(labels);
    X = reordercats(X,labels);
    Y = sobl(:,:,j,i);
    b = bar(X,Y,0.85, 'EdgeColor','Flat',...
                      'FaceColor','Flat',...
                      'EdgeAlpha',1.0,...
                      'LineWidth',1);
    
    % adjust colors
    colormap(cmap)
    for k = 1:size(Y,2)
      b(k).CData = k;
      b(k).EdgeColor = cmap2(k,:);
    end
    
    % adjust axis limits
    if i==1
      ylim([0 0.3])
    else
      ylim([0 0.75])
    end
    
    % print a legend on the first subplot
    if (i == 1) && (j==1)
      legend(shapes,'FontName','Helvetica','FontWeight','Normal','FontSize',10)
    end
    
    % add titles
    if (i==1) && (j==1)
      ylabel('1^{st} Order Sobol Indices',...
              'FontName','Helvetica',...
              'FontWeight','Normal',...
              'FontSize',10);
      title(sprintf('{\\bf%s}: exchange flux',alphabet(ctr)),...
              'FontName','Helvetica',...
              'FontWeight','Normal',...
              'FontSize',10)
    elseif (i==2) && (j==1)
      ylabel('Total Effect Sobol Indices',...
              'FontName','Helvetica',...
              'FontWeight','Normal',...
              'FontSize',10);
      title(sprintf('{\\bf%s}: exchange flux',alphabet(ctr)),...
              'FontName','Helvetica',...
              'FontWeight','Normal',...
              'FontSize',10)
    elseif (j==2)
      title(sprintf('{\\bf%s}: exchange zone area',alphabet(ctr)),...
              'FontName','Helvetica',...
              'FontWeight','Normal',...
              'FontSize',10)  
    end
    
    % create nice ticks and tick labels
    ax = gca;
    ax.TickLabelInterpreter = "LaTeX";
    ytickformat('%.2f')
    if j ~= 1
      set(gca,'Yticklabel',[])
    end
    
    % add a grid and a box
    ax.YGrid = 'on';
    box on
    
    % step the counter
    ctr = ctr+1;
  end
end

% add an invisible axis on top, to export all of the figure without crop
ax = axes;
ax.Position = [0,0,1,1];
ax.Color = 'None';
ax.XTick = [0 1];
ax.XColor = 'white';
ax.YTick = [0 1];
ax.YColor = 'white';

% export the result
exportgraphics(fig,"./sobol.pdf",'ContentType','vector')
























