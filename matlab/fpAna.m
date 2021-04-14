classdef fpAna < fpSimple
  %fpAna Simplified floodplain model with semi-analytical solution method.
  %   Objects of this class describe simplified floodplain models that are
  %   solved by a semi-analytical method.
  %
  %   See also FPSIMPLE.
  properties (AbortSet,SetObservable)
    % Number of members that should be considered in the infinite series.
    % If too small, the solution is inaccurate. If too large, the problem
    % might become instable.
    N (1,1) double {isinteger} = 10
    % Switch between equidistant or chebyshev node placement for
    % constructing the system of equations. Chebyshev nodes are supposed to
    % help with oscillations. Equidistant nodes are cleaner.
    placement (1,1) string {mustBeMember(placement,["chebyshev";
                                                  "equidistant"])}  = ...
                                                              "equidistant"
  end
  
  properties (Dependent)
    Qex
    Qeast
    Qwest
  end
  
  properties (SetAccess = private)
    % Coefficients of the infinite series. A(1) is A_0, the rest follows
    % incrementally. --> A(n+1) = A_n.
    A
  end
  
  methods
    function obj = fpAna(opt)
      % This is the constructor of the fpAna class. It can be called
      % without arguments, which leads to the usage of default values. It
      % can also be called with name/value syntax, where all settable
      % properties are allowed as names. Input argument validation ensures
      % that all inputs are valid.
      arguments
        opt.?fpAna
      end
      
      % pass properties forward to object
      fn = fieldnames(opt);
      for i = 1:numel(fn)
        obj.(fn{i}) = opt.(fn{i});
      end
      
      obj.isReady         = false;
      
      % add listener for dynamic resolving/replotting
      m                   = ?fpAna;
      props               = m.PropertyList([m.PropertyList.SetObservable]);
      obj.changeListener  = addlistener(obj,props,'PostSet',@obj.propChange);
      
      if obj.autoSolve
        obj = obj.solve();
      end
    end
    
    function Qw = get.Qwest(obj)
      % define shorthands
      f0    = obj.fNorth(0);
      yMa   = obj.wMax;
      kap   = obj.kappa;
      
      % determine Q(n) and sum up
      Qw    = obj.Tx*(obj.h1-obj.h2)/obj.L*f0;
      for n = 1:obj.N
        c   = n*pi/obj.L;
        Qw  = Qw - obj.Tx*obj.A(n+1)*(cosh(c*kap*f0)-1)/...
                                                    (kap*cosh(c*kap*yMa));
      end
    end

    function Qe = get.Qeast(obj)
      % define shorthands
      fL    = obj.fNorth(obj.L);
      yMa   = obj.wMax;
      kap   = obj.kappa;
      
      % determine Q(n) and sum up
      Qe    = obj.Tx*(obj.h1-obj.h2)/obj.L*fL;
      for n = 1:obj.N
        c   = n*pi/obj.L;
        Qe  = Qe - obj.Tx*obj.A(n+1)*cos(n*pi)*(cosh(c*kap*fL)-1)/...
                                                    (kap*cosh(c*kap*yMa));
      end
    end

    function Qs = get.Qex(obj)
      % evaluate stream function at some positions
      xT = linspace(0,obj.L,obj.res*3);
      yT = 0*xT;
      ps = obj.psi(xT,yT);
      
      % determine total unique exchange discharge
      if obj.h2 < obj.h1
        Qs = min([ps(1) ps(end)]) - min(ps);
      else
        Qs = max(ps)-max([ps(1) ps(end)]);
      end
    end
      
    function obj = solve(obj)
      % This method solves the given floodplain model with a
      % semi-analytical method. The results are stored in the object itself
      % (property "A"). Solving is necessary ahead of plotting or 
      % requesting results. This method does not take any arguments.
      M     = (0.5*obj.N)^2;
      
      % determine node positions
      switch obj.placement
        case "chebyshev"
          x     = 0.5*obj.L+0.5*obj.L*cos((2*(1:M)'-1)/(2*M)*pi);
        case "equidistant"
          x     = linspace(0,obj.L,M)';
      end
      y     = obj.fNorth(x);
      
      w = 1/numel(x);
      w = w/sum(w);

      % evaluate U and F for all points and all n
      UM  = obj.U(0:obj.N,x,y);
      Fvec  = obj.F(x,y);

      % assemble and solve system of equations
      rhs         = (UM.*w')*Fvec;
      lhs         = (UM.*w')*UM';
      obj.A       = lhs\rhs;
      
      % mark the model as solved
      obj.isReady = true;
    end
    
    function h = h(obj,x,y)
      % Hydraulic head query function h(x,y).
      arguments
        obj
        x (:,:) double
        y (:,:) double {mustBeEqualSize(y,x)}
      end
      
      % validate that model has been solved and is ready
      obj.check();
      
      % evaluate n terms and sum up
      h     = obj.h1 + (obj.h2-obj.h1)/obj.L*x;
      for n = 1:obj.N
        h   = h + obj.A(n+1)*obj.V(n,x,y);
      end
    end

    function ps = psi(obj,x,y)
      % Stream function value query function psi(x,y).
      arguments
        obj
        x (:,:) double
        y (:,:) double {mustBeEqualSize(y,x)}
      end
      
      % validate that model has been solved and is ready
      obj.check();
      
      % evaluate n terms and sum up
      ps    = obj.A(1) + (obj.h2-obj.h1)/obj.L*y;
      An    = obj.A(2:obj.N+1)';
      U     = obj.U(1:obj.N,x,y);
      ps    = ps + shiftdim(pagemtimes(An,U),1);
      
      % multiply with constant factor
      ps    = -ps*obj.Tx;
    end
        
    function qx = qx(obj,x,y)
      % Darcy-velocity query function for x-direction qx(x,y). in L^2/T!
    
      % validate that model has been solved and is ready
      obj.check();

      % define shorthands
      yMa   = obj.wMax;
      kap   = obj.kappa;
      
      % evaluate n terms and sum up
      qx    = (obj.h2-obj.h1)/obj.L;
      for n = 1:obj.N
        c   = n*pi/obj.L;
        qx  = qx + obj.A(n+1)*cos(c*x).*sinh(c*kap*y)*c/cosh(c*kap*yMa);
      end
      
      % multiply with constant factor
      qx    = -obj.Tx * qx;
    end

    function qy = qy(obj,x,y)
      % Darcy-velocity query function for y-direction qy(x,y). in L^2/T!

      % validate that model has been solved and is ready
      obj.check();

      % define shorthands
      yMa   = obj.wMax;
      kap   = obj.kappa;
      
      % evaluate n terms and sum up
      qy    = 0;
      for n = 1:obj.N
        c   = n*pi/obj.L;
        qy  = qy + obj.A(n+1)*sin(c*x).*cosh(c*kap*y)*c*kap/cosh(c*kap*yMa);
      end
      
      % multiply with constant factor
      qy    = -obj.Ty * qy;
    end
  end
  
  methods (Access = private)
    function value = R(obj,x)
      % This is a helper function for the R term.
      value = -obj.qNorth*x;
    end
    
    function F = F(obj,x,y)
      % This is a helper function for right-hand-side term.
      F = -obj.R(x)/obj.Tx+(obj.h1-obj.h2)/obj.L*y;
    end
    
    function v = U(obj,n,x,y)
      % This is a helper function for the left-handside term.
      arguments
        obj
        n (:,1)
        x
        y
      end
      c     = n*pi/obj.L;
      kap   = sqrt(obj.Tx/obj.Ty);
      yMa   = obj.wMax;
      
      % add another dimension on the left (size 1)
      x     = shiftdim(x,-1);
      y     = shiftdim(y,-1);
      b1    = pagemtimes(c,kap*(y-yMa));
      b2    = pagemtimes(c,kap*(-y-yMa));
      b3    = pagemtimes(c,x);
      
      v     = ( exp(b1) + exp(b2) )./( 1 + exp(-2*kap*c*yMa) ) .*cos(b3);
      v     = v./kap;
    end
    
    function v = V(obj,n,x,y)
      % This is a helper function to reconstruct the heads in the end.
      c     = n*pi/obj.L;
      yMa   = obj.wMax;
      kap   = obj.kappa;
      v     = ( exp(c*kap*(y-yMa)) - exp(c*kap*(-y-yMa)) )./...
                                    ( 1 + exp(-2*kap*c*yMa) ) .*sin(c*x);
    end
  end
  
  methods (Access = protected)
    function propChange(obj,~,~)
      obj.isReady      = false;
      
      if obj.autoSolve
        obj.solve();
      end
      
      updatePlot(obj);
    end
  end
end

function mustBeEqualSize(a,b)
  % Test for equal size
  if ~isequal(size(a),size(b))
    eid = 'Size:notEqual';
    msg = 'Size of first input must equal size of second input.';
    throwAsCaller(MException(eid,msg))
  end
end
