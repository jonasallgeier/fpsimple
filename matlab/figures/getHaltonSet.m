function [out,p] = getHaltonSet(nM,opts)
  % get a halton set of requested size and within requested parameter
  % bounds
  arguments
    nM          (1,1) double {mustBePositive}     = 1000;
    opts.I      (1,2) double {mustBeNonnegative}  = [0 3e-2];
    opts.T      (1,2) double {mustBeNonnegative}  = [1e-6 5e-3];
    opts.aniso  (1,2) double {mustBeNonnegative}  = [0.1 10];
    opts.L      (1,2) double {mustBePositive}     = [100 3000];
    opts.rMin   (1,2) double {mustBeNonnegative}  = [2/5 1];
    opts.rMax   (1,2) double {mustBeNonnegative}  = [1/10 1/2];
    opts.qNstar (1,2) double                      = [-4 0];
  end
  
  I = opts.I;           % <-- bounds of ambient hydraulic gradient
  L = opts.L;           % <-- bounds of domain length
  T = log(opts.T);      % <-- bounds of average transmissivity
  a = log(opts.aniso);  % <-- bounds of anisotropy ratio
  r = opts.rMin;        % <-- bounds of minimum width ratio (wMin/wMax)
  R = opts.rMax;        % <-- bounds of maximum width ratio (wMax/L)
  q = opts.qNstar;      % <-- bounds of normalized northern discharge
  
  % create the scrambled halton sequence
  nPar    = 7;
  P       = haltonset(nPar,'Leap',nthprime(nPar+3)-1,'Skip',10*nPar);
  P       = scramble(P,'RR2');
  
  % rescale the sequence to obtain primary parameters
  % ambI | Kmean | anisotropy ratio | L | wMin/wMax | wMax/L | R (maybe in terms of Kmean?)
  p.I       =      I(1)+ (I(2)-I(1)).*P(1:nM,1);
  p.Tmean   = exp( T(1)+ (T(2)-T(1)).*P(1:nM,2) );
  p.aniso   = exp( a(1)+ (a(2)-a(1)).*P(1:nM,3) );
  p.L       =      L(1)+ (L(2)-L(1)).*P(1:nM,4);
  p.ratMin  =      r(1)+ (r(2)-r(1)).*P(1:nM,5);
  p.ratMax  =      R(1)+ (R(2)-R(1)).*P(1:nM,6);
  p.QNstar  =      q(1)+ (q(2)-q(1)).*P(1:nM,7);
  
  % transform to get real input parameters
  out.L       = p.L;                        % domain length
  out.h1      = p.I.*p.L;                     % first fixed head
  out.h2      = 0.*out.h1;                % second fixed head
  out.Tx      = sqrt(p.Tmean.^2.*p.aniso);    % transmissivity in x
  out.Ty      = sqrt(p.Tmean.^2./p.aniso);    % transmissivity in y
  out.wMax    = p.ratMax.*p.L;                % maximum width
  out.wMin    = p.ratMin.*out.wMax;         % minimum width
  out.Q0      = (out.h1-out.h2)./p.L.*out.Tx.*(out.wMax-out.wMin); % normalization discharge
  out.qNorth  = p.QNstar.*out.Q0./p.L;        % northern influx (L³/L/T)
end