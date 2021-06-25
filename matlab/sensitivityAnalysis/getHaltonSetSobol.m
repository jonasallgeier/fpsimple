function [A,B,AB] = getHaltonSetSobol(nM,opts)
  % get a halton set of requested size and within requested parameter
  % bounds for a Sobol sensitivity analysis
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
  
  % create the scrambled halton sequence
  nPar    = 7;
  P       = haltonset(2*nPar,'Leap',nthprime(2*nPar+3)-1,'Skip',10*2*nPar);
  P       = scramble(P,'RR2');
  
  P = P(1:nM,:);
  
  A = P(:,1:nPar);
  B = P(:,nPar+1:end);
  
  A = transformToTable(A,opts);
  B = transformToTable(B,opts);
  
  % unfortunately 3D tables do not exist, so we have to fall back to
  % regular 3D arrays for the AB matrices
  for i = 1:nPar
    AB(:,:,i) = A{:,:}; %#ok<AGROW>
    AB(:,i,i) = B{:,i}; %#ok<AGROW>
  end
end

function out = transformToTable(P,opts)
  I = opts.I;           % <-- bounds of ambient hydraulic gradient
  L = opts.L;           % <-- bounds of domain length
  T = log(opts.T);      % <-- bounds of average transmissivity
  a = log(opts.aniso);  % <-- bounds of anisotropy ratio
  r = opts.rMin;        % <-- bounds of minimum width ratio (wMin/wMax)
  R = opts.rMax;        % <-- bounds of maximum width ratio (wMax/L)
  q = opts.qNstar;      % <-- bounds of normalized northern discharge
  
  I       =      I(1)+ (I(2)-I(1)).*P(:,1);
  Tmean   = exp( T(1)+ (T(2)-T(1)).*P(:,2) );
  aniso   = exp( a(1)+ (a(2)-a(1)).*P(:,3) );
  L       =      L(1)+ (L(2)-L(1)).*P(:,4);
  ratMin  =      r(1)+ (r(2)-r(1)).*P(:,5);
  ratMax  =      R(1)+ (R(2)-R(1)).*P(:,6);
  QNstar  =      q(1)+ (q(2)-q(1)).*P(:,7);
  
  out = table(I,Tmean,aniso,L,ratMin,ratMax,QNstar);
end