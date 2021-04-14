function [sigma,corr] = errorProp(xObs,yObs,mdl,pars)
  % obtain error estimates with linearized uncertainty estimation
  
  % x-Values
  xPr = xObs;
  
  % predicted y-Values with fitted parameters
  yPr = mdl(pars,xPr);
  
  % determine jacobian through numerical perturbation analysis
  J = NaN(numel(yObs),numel(pars));
  for j = 1:size(J,2)
    P1      = pars;
    P2      = P1;
    perturb = 1e-3 * P2(j);
    P2(j)   = P2(j)+perturb;
    J(:,j)  = ( mdl(P2,xPr) - mdl(P1,xPr) )/perturb;
  end
  
  % approximate variance of "measurements"
  varY = sum((yObs-yPr).^2)./(numel(yObs)-numel(pars));
  
  % covariance matrix of "measurement" errors
  Cyy = eye(numel(yObs))*varY;
  
  % the Hessian
  H = J'*(Cyy\J); % <-- J' * inv(Cyy) * J
  
  % covariance matrix of estimated coefficients
  Cpp = inv(H);
  
  % back out the standard deviations of parameters
  sigma = sqrt(diag(Cpp));
  
  % back out the correlations between parameter standard deviations
  corr = Cpp./(sigma * sigma');
end