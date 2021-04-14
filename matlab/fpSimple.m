classdef fpSimple < matlab.mixin.SetGet & matlab.mixin.CustomDisplay
  %fpSimple Abstract simplified floodplain model.
  %   This class summarizes common properties and methods of analytical and
  %   numerical simplified floodplain models. It cannot be instantiated but
  %   serves as the parent class for fpAna.
  %
  %   See also FPANA.
  properties (AbortSet,SetObservable)
    % Domain length in x direction (dim L).
    L           (1,1) double {mustBeFinite,mustBePositive}          = 1000
    % Maximum domain length in y direction (dim L).
    wMax        (1,1) double {mustBeFinite,mustBePositive}          = 360
    % Minimum domain length in y direction (dim L).
    wMin        (1,1) double {mustBeFinite,mustBePositive}          = 250
    % Hydraulic head at inlet x = 0 (dim L).
    h1          (1,1) double {mustBeFinite}                         = 10
    % Hydraulic head at inlet x = L (dim L).
    h2          (1,1) double {mustBeFinite}                         = 7
    % Hydraulic conductivity in x-direction (dim L^2/T).
    Tx          (1,1) double {mustBeFinite,mustBePositive}          = 1e-3
    % Hydraulic conductivity in y-direction (dim L^2/T).
    Ty          (1,1) double {mustBeFinite,mustBePositive}          = 1e-3
    % Influx at y = fNorth in y-direction (dim L^2/T).
    qNorth      (1,1) double {mustBeFinite}                         = 0
    % Depth-integrated porosity (m*poros) (dim L).
    por         (1,1) double {mustBeFinite}                         = 1
    % Boundary shape at y = fNorth. Possible values are "parabolic", 
    % "sinusoidal", "cosinusoidal", "ellipsoidal", "linear", and 
    % "asymmetrical".
    shape       (1,1) string {mustBeMember(shape,["parabolic";
                                                  "sinusoidal";
                                                  "cosinusoidal";
                                                  "ellipsoidal";
                                                  "linear";
                                                  "asymmetrical";
                                                  "bump";
                                                  "bump2";
                                                  "composite";
                                                  "composite2"])}  = ...
                                                                "cosinusoidal"
    % Control on lateral resolution. This is the number of points used to
    % approximate the northern boundary.
    res         (1,1) double {mustBeInteger,mustBeGreaterThan(res,5)} = 40
  end
  
  properties (Dependent)
    % Square root of anisotropy ratio sqrt(Kx/Ky) (dimless).
    kappa
  end
  
  properties (Dependent,Abstract)
    % Total absolute unique exchange discharge across y = 0 (dim L^3/T).
    Qex
    % Net discharge across x = L (dim L^3/T).
    Qeast
    % Net discharge across x = 0 (dim L^3/T).
    Qwest
  end
  
  properties (SetAccess = protected, Hidden)
    % Structure containing all current contours and patches of an object.
    % There are three fields:
    %   1. go: stores the actual graphics objects.
    %   2. axis: stores the axis object to decide whether a plot is current
    %      or not.
    %   3. id: stores strings to differentiate between graphics objects.
    stor = struct('go',gobjects(0),'axis',[],'id',string.empty);
    % Logical that indicates whether the analytical or numerical model has
    % been solved for the current properties.
    isReady       = false;
  end
  
  properties (Hidden)
    % Listener that triggers automatic re-solution of models when
    % properties change.
    changeListener
    % Logical switch to allow or prohibit auto-updating of plots.
    autoUpdate  (1,1) logical                                       = true
    % Logical switch to allow or prohibit auto-solving of models.
    autoSolve     (1,1) logical                                     = true
  end
  
  methods (Abstract)
    % Hydraulic head query function h(x,y).
    h
    % Stream function value query function psi(x,y).
    psi
    % Darcy-velocity query function for x-direction qx(x,y).
    qx
    % Darcy-velocity query function for y-direction qy(x,y).
    qy
    % Solve/run the model.
    solve
  end
  
  methods (Access = protected, Abstract)
    % This method is called whenever a settable property is changed.
    propChange
  end
  
  methods
    function obj = fpSimple()
      % This is the constructor of the fpSimple class. It can only be
      % called implicitly by subclasses.
      
      % Make sure that a recent version of Matlab is used.
      if verLessThan('matlab','9.10.0.1538726')
        warning(['This code was developed with Matlab R2021a.'...
                 'Consider updating to ensure full compatibility.']);
      end
    end
    
    function set.h2(obj,val)
      % make sure that h1 > h2
      if val==obj.h1 %#ok<MCSUP>
        warning('Setting h1==h2 does not yield accurate results.');
      end
      obj.h2 = val;
    end
%     
    function value = get.kappa(obj)
      value = sqrt(obj.Tx/obj.Ty);
    end
    
    function set.autoSolve(obj,input)
      % If autoSolve is turned on, solve the model directly.
      if input == true
        obj.autoSolve = true;
        obj.solve();
      else
        obj.autoSolve = false;
      end
    end
    
    function plot(obj,opts)
      % This function can plot the geometry and solution of a fpSimple
      % model. It accepts name value arguments to hide or show specific
      % aspects of model (e.g, 'psi', 'h' or 'outline'). Each of these
      % names must be followed by a logical value indicating whether it
      % should be plotted or not. Subclasses of fpSimple might extend the
      % list of allowable names.
      arguments
        obj
        opts.psi      (1,1) logical = true
        opts.h        (1,1) logical = true
        opts.outline  (1,1) logical = true
        opts.divide   (1,1) logical = true
      end
      
      wasHold = ishold();
      hold on
      
      [X,Y] = obj.meshgrid(); % get meshgrid of points for contours
      
      for mode = ["psi","h","divide"]
        switch mode
          case "psi"
          	Z   = obj.psi(X,Y);
            lvl = linspace(min(Z(:)),max(Z(:)),11);
            lvl = lvl(2:end-1);
            clr = [0.9 0.91 0.93];
          case "h"
            Z   = obj.h(X,Y);
            lvl = linspace(min(Z(:)),max(Z(:)),11);
            clr = 'flat';
          case "divide"
            Z   = obj.psi(X,Y);
            lvl = obj.psi([0;0;obj.L;obj.L],[obj.wMin;0;0;obj.wMin]);
            clr = [0.7 0.2 0.6];
        end
        
        [flag,pos] = onlyRefresh(obj,mode);
        if ~obj.isReady
          % model is still unsolved; delete old plots
          delete(obj.stor.go(pos));
          obj.stor.go(pos) = [];
          obj.stor.id(pos) = [];
        elseif flag && obj.isReady
          % an updatable plot exists --> update
          set(obj.stor.go(pos),'XData',X,'YData',Y,'ZData',Z,'LevelList',lvl);
        else
          % no updatable plot exists --> create a new one
          [~,c] = contour(X,Y,Z,lvl,'LineWidth',2,'Color',clr);
          obj.stor.go(end+1) = c;
          obj.stor.id(end+1) = mode;
        end
      end
      colormap(colorcet('D1A'));
      caxis([min(obj.h2,obj.h1) max(obj.h2,obj.h1)]);
      
      % create an outline of the geometry defined by some points
      p           = outline(obj);
      [flag,pos]  = onlyRefresh(obj,"outline");
      if flag
        % an updatable plot exists --> update
        set(obj.stor.go(pos),'XData',p(:,1),'YData',p(:,2));
      else
        % no updatable plot exists --> create a new one
        thePatch = patch('XData',p(:,1),'YData',p(:,2),'LineWidth',2,...
                          'EdgeColor',[0.4 0.4 0.4],'FaceAlpha',0);
        obj.stor.go(end+1)  = thePatch;
        obj.stor.id(end+1)  = "outline";
      end

      obj.togglePlots(opts); % switch contours/patches on or off
      obj.stor.axis = gca; % save the current axis for later identification
      
      % set basic axis properties
      xlim([0 obj.L])
      ylim([0 obj.wMax])
      daspect([1 1 1])
      if ~wasHold
        hold off
      end
    end

    function A = area(obj)
      % Surface area of the domain.
      p   = outline(obj);
      A   = polyarea(p(:,1),p(:,2));
    end
    
    function l = arcLength(obj)
      p   = outline(obj);
      lTot = sum(sqrt(sum((p(1:end-1,:)-p(2:end,:)).^2,2)));
      l = lTot - obj.L -2*obj.wMin;
    end
    
    function area = areaEx(obj)
      % Surface area of the exchange zone.
      Pini    = [0, 0; obj.L, 0];
      if obj.h1 > obj.h2
        lvl     = min(obj.psi(Pini(:,1),Pini(:,2)));
      else
        lvl     = max(obj.psi(Pini(:,1),Pini(:,2)));
      end
      s       = obj.isolines(lvl,'type',"psi");
      points  = [];
      for i = 1:numel(s)
        % we do not want northern oscillation fragments
        if min(s(i).y./obj.wMin) < 0.5
          points = [ points; s(i).x, s(i).y]; %#ok<AGROW>
        end
      end
      area = polyarea(points(:,1),points(:,2));
    end
    
    function [f,t,pd] = ttd(obj,opt)
      % Travel time distribution (cdf) of water within the exchange zone.
      arguments
        obj
        opt.k     = 50
        opt.normalize (1,1) string {mustBeMember(opt.normalize,[ "none";
                                                  "median";
                                                  "mean";
                                                  "max"])}  = "median"
        opt.xq    = [];
      end
      k       = opt.k;
      norm    = opt.normalize;
      
      % get psi-equally-spaced contour lines between min and max
      psiH    = min(obj.psi([0 obj.L],[0 0]));
      xL      = linspace(0,obj.L,obj.res*3);
      psiL    = min(obj.psi(xL,0*xL));
      lvl     = linspace(psiL,psiH,k);
      s       = obj.isolines(lvl(2:end),'type',"psi");
      
      % determine travel times by numerical integration
      t       = NaN(k,1);
      flux    = NaN(k,1);
      t(1)    = 0;
      flux(1) = 0;
      for i = 1:numel(s)
        x         = s(i).x;
        y         = s(i).y;
        v         = sqrt((obj.qx(x,y)/obj.por).^2+(obj.qy(x,y)/obj.por).^2);
        ds        = sqrt(diff(x).^2+diff(y).^2);
        vM        = mean([v(1:end-1), v(2:end)],2);
%         vM        = (v(2:end)-v(1:end-1))./log(v(2:end)./v(1:end-1));
%         vM(isnan(vM)) = v(isnan(vM));
        dt        = ds./vM;
        t(i+1)    = sum(dt);
        flux(i+1) = s(i).level-psiL;
      end
      % determine flux fractions
      f = flux./flux(end);
      
      % normalize if requested
      switch norm
        case "median"
          t = t./median(t);
        case "mean"
          t = t./mean(t);
        case "max"
          t = t./max(t);
      end
      
      pd = makedist('PiecewiseLinear','x',t,'Fx',f);
      if ~isempty(opt.xq)
        t = opt.xq;
        f = pd.cdf(t);
      end
    end
    
    function [x,y] = fNorth(obj,x)
      % Function defining the northern boundary of the model geometry. If
      % requested with a single output argument, it provides y at requested
      % x locations. If two output arguments are requested, they take the
      % format [x,y].
      arguments
        obj
        x = linspace(0,obj.L,obj.res)
      end
      
      yMi = obj.wMin;
      yMa = obj.wMax;
      s   = obj.L;
      
      switch obj.shape
        case "parabolic"
          y = yMi+(yMa-yMi)*4/s^2*x.*(s-x);
        case "sinusoidal"
          y = yMi+(yMa-yMi)*sin(pi*x/s);
        case "cosinusoidal"
          y = yMi+0.5*(yMa-yMi)+0.5*(yMi-yMa)*cos(2*pi*x/s);
        case "ellipsoidal"
          y = yMi+(yMa-yMi)*2/s*sqrt((s/2)^2-(x-s/2).^2);
        case "linear"
          y = -abs(x-s/2)/s*2*(yMa-yMi)+yMa;
        case "asymmetrical"
          % we solve a system of equations to get a polynomial of order 4
          % that goes through the points (0|wMin), (L|wMin), (p|wMax) and
          % has a slope of q at x = 0
          p   = s*6/10;
          q   = 0.00;

          M   = [ 0       0       0     0   1;
                  p^4     p^3     p^2   p   1;
                  s^4     s^3     s^2   s   1;
                  4*p^3   3*p^2   2*p   1   0;
                  0       0       0     1   0];

          rhs = [yMi; yMa; yMi; 0; q];
          asy = M\rhs;
          str = sprintf('@(x) %g*x.^4+%g*x.^3 +%g*x.^2 +%g*x+%g',asy);
          f   = str2func(str);
          y   = f(x);
        case "bump"
          y = yMi+(yMa-yMi)*exp(-1./(1-((2*x/s-1).^2)))./exp(-1);
        case "bump2"
          y = yMi+(yMa-yMi)*exp(-1./(1-((2*x/s-1).^4)))./exp(-1);
        case {"composite","composite2"}
          switch obj.shape
            case "composite"
              LS = 0.35;
              LB = 0.6;
            case "composite2"
              LS = 0.2;
              LB = 0.7;
          end

          A = max(0.5-0.5*LB-0.5*LS,0);
          B = min(0.5-0.5*LB+0.5*LS,0.5);
          C = max(0.5+0.5*LB-0.5*LS,0.5);
          D = min(0.5+0.5*LB+0.5*LS,1);

          f = @(x) 0.5*( (1-cos(pi/(B-A)*(x-A)) ).*(x>=A & x < B) + ...
            (1+cos(pi/(D-C)*(x-C)) ).*(x>=C & x <= D) + 2*double(x >=B & x < C) );
          y = yMi+(yMa-yMi)*f(x/s);
        otherwise
          error('Northern boundary outline shape is unknown.');
      end
      
      % if single output request: give y only
      if nargout == 1
        x = y;
      end
    end
  end
  
  methods (Hidden)
    % We hide the following methods to de-clutter Matlab's autocomplete
    % suggestions for the class.    
    function getdisp(obj,varargin)
      obj.disp();
    end
    
    function setdisp(obj,varargin)
      obj.disp();
    end
    
    function s = isolines(obj,lvl,opt)
      % Extract isolines from a field.
      arguments
        obj
        lvl (:,1)
        opt.type = "h"
      end
      
      [X,Y,xq,yq] = obj.meshgrid();
      
      switch opt.type
        case "h"
          Z = obj.h(X,Y);
        case "psi"
          Z = obj.psi(X,Y);
      end
      
      % contourc does not understand scalar levels, so expand if needed
      if isscalar(lvl)
        lvl = lvl.*ones(1,2);
      end
      c     = contourc(xq,yq,Z,lvl);
      
      k     = 0;
      col   = 1;
      while col<size(c,2)
        col = col+c(2,col)+1;
        k   = k+1;
      end
      s(k)  = struct();
      col   = 1;
      for i = 1:k
         s(i).level = c(1,col);
         num        = c(2,col);
         idx        = col+1:col+num;
         s(i).x     = c(1,idx).';
         s(i).y     = obj.fNorth(c(1,idx).').* c(2,idx).';
         col        = col+num+1;
      end
    end
  end
 
  methods (Static)
    function help(~)
      % A help reference function enhancing user experience.
      %
      % See also HELP.
      help fpSimple;
      methods fpSimple;
    end
    
    function doc(~)
      % A doc reference function enhancing user experience.
      %
      % See also DOC.
      doc fpSimple;
    end
  end
  
  methods (Access = protected)
    function togglePlots(obj,opts)
      % Toggle visibility of different plot elements.
      
      names = string(fieldnames(opts));
      for i = 1:numel(names)
        pos = find(obj.stor.id==names(i),1);
        if isempty(pos)
          continue
        end
        if opts.(names(i))
          obj.stor.go(pos).Visible = 'on';
        else
          obj.stor.go(pos).Visible = 'off';
        end
      end
    end
    
    function p = outline(obj)
      % Get points defining the domain geometry outline.
      [x,y]   = obj.meshgrid('y',1);
      p       = [x' y'; obj.L 0; 0 0];      

      % arrange points in clockwise manner
      p       = unique(p,'rows');
      ctr     = [obj.L/2 0];
      angles  = atan2d(p(:,2)-ctr(:,2),p(:,1)-ctr(:,1));
      [~,ord] = sort(angles,'descend');
      p       = p(ord,:);
      p       = [p; p(1,:)];
    end
    
    function [X,Y,x,y] = meshgrid(obj,opt)
      % Get a meshgrid of points within the domain.
      arguments
        obj
        opt.x (:,1) {mustBeLessThanOrEqual(opt.x,1)} = linspace(0,1,obj.res)
        opt.y (:,1) {mustBeLessThanOrEqual(opt.y,1)} = linspace(0,1,...
                                        max(ceil(obj.res*obj.wMin/obj.L),3))
      end
      x = opt.x;
      y = opt.y;
      
      x         = x.*obj.L;
      [X,Y]     = meshgrid(x,y);
      [~,fac]   = obj.fNorth(X(:));
      Y(:)      = Y(:).*fac; 
    end
    
    function obj = check(obj)
      % An internal helper function to verify that model has been solved.
      if ~obj.isReady
        str = input('Model is not solved yet. Do it now? (y/n)','s');
        if isempty(str) || strcmp(str,'y')
          obj.solve();
        else
          error('Model needs to be solved first.');
        end
      end
    end
    
    function updatePlot(obj)
      % A function to update an existing plot if it is still the active
      % one. It replaces existing contours/patches of the current axes if
      % they are still in the storage and requested to be updated.
      if obj.autoUpdate...
            && ~isempty(obj.stor.go)...
            && isequal(gca,obj.stor.axis)
        % in case we have some old contours hanging around, delete them
        if ~obj.isReady
          toDelete = false(1,numel(obj.stor.go));
          for i = 1:numel(obj.stor.go)
            if isa(obj.stor.go(i),'matlab.graphics.chart.primitive.Contour')
              delete(obj.stor.go(i));
              toDelete(i) = true;
            end
          end
          obj.stor.go(toDelete) = [];
          obj.stor.id(toDelete) = [];
        end
        plot(obj);
      end
    end
    
    function [flag,pos] = onlyRefresh(obj,id)
      % A helper function to decide whether a particular plot element can
      % be updated or if it has to be constructed from scratch.
      toDelete = ~isvalid(obj.stor.go);
      obj.stor.go(toDelete) = [];
      obj.stor.id(toDelete) = [];
      
      % A graphics object can be updated if the autoUpdate switch is 'on',
      % the object's axis is the current axis and if the object is still
      % valid (i.e., not yet deleted).
      flag  = true;
      flag  = flag && obj.autoUpdate;
      flag  = flag && isequal(gca,obj.stor.axis);
      pos   =  find(obj.stor.id==id,1);
      if ~isempty(pos)
        flag  = flag && isvalid(obj.stor.go(pos));
      else
        flag  = false;
      end
    end
      
    function header = getHeader(obj)
      % A Matlab helper function to display fpSimple objects nicely.
      hStr = matlab.mixin.CustomDisplay.getClassNameForHeader(obj);
      header = sprintf('%s describing a simplified floodplain.\n\t',hStr);
    end
    
    function propgrp = getPropertyGroups(~)
      % A Matlab helper function to display fpSimple objects nicely.
      proplist = {'L','wMax','wMin','por'};
      propgrp(1) = matlab.mixin.util.PropertyGroup(proplist);
      
      proplist = {'h1','h2','Tx','Ty','qNorth'};
      propgrp(end+1) = matlab.mixin.util.PropertyGroup(proplist);
    end
    
    function footer = getFooter(~)
      % A Matlab helper function to display fpSimple objects nicely.
      footer = sprintf('\nUse plot to visualize.\n');
    end
  end
end


