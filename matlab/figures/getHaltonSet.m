function out = getHaltonSet(nM,opts)
  arguments
    nM          = 1000;
    opts.I      = [0 3e-2];
    opts.T      = [1e-6 5e-3];
    opts.aniso  = [0.1 10];
    opts.L      = [100 3000];
    opts.rMin   = [2/5 1];
    opts.rMax   = [1/10 1/2];
    opts.qNstar = [-4 0];
  end
  
  I = opts.I;
  L = opts.L;
  T = log(opts.T);
  a = log(opts.aniso);
  r = opts.rMin;
  R = opts.rMax;
  q = opts.qNstar;
  
  nPar    = 7;
  P       = haltonset(nPar,'Leap',nthprime(nPar+3)-1,'Skip',10*nPar);
  P       = scramble(P,'RR2');
  
  % primary parameters
  % ambI | Kmean | anisotropy ratio | L | wMin/wMax | wMax/L | R (maybe in terms of Kmean?)
  I       =      I(1)+ (I(2)-I(1)).*P(1:nM,1);
  Tmean   = exp( T(1)+ (T(2)-T(1)).*P(1:nM,2) );
  aniso   = exp( a(1)+ (a(2)-a(1)).*P(1:nM,3) );
  L       =      L(1)+ (L(2)-L(1)).*P(1:nM,4);
  ratMin  =      r(1)+ (r(2)-r(1)).*P(1:nM,5);
  ratMax  =      R(1)+ (R(2)-R(1)).*P(1:nM,6);
  QNstar  =      q(1)+ (q(2)-q(1)).*P(1:nM,7);
  
  % transform to get real input parameters
  out.L       = L;
  out.h1      = I.*L;
  out.h2      = 0.*out.h1;
  out.Tx      = sqrt(Tmean.^2.*aniso);
  out.Ty      = sqrt(Tmean.^2./aniso);
  out.wMax    = ratMax.*L;
  out.wMin    = ratMin.*out.wMax;
  out.Q0      = (out.h1-out.h2)./L.*out.Tx.*(out.wMax-out.wMin);
  out.qNorth  = QNstar.*out.Q0./L;
end