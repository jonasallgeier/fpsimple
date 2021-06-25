function createDataSet(nM,shape)
  arguments
    nM (1,1) double = 1500 % number of model runs
    shape (1,1) string = "cosinusoidal" % shape type
  end
  oldPath = addpath('../');
  
  % create the halton set
  [A,B,AB]    = getHaltonSetSobol(nM,'qNstar',[-3 0]);
  
  % pre-allocate some outputs
  fAflux      = NaN(height(A),1);
  fAarea      = NaN(height(A),1);
  fBflux      = NaN(height(B),1);
  fBarea      = NaN(height(B),1);
  fABflux     = NaN(size(AB,1),size(AB,3));
  fABarea     = NaN(size(AB,1),size(AB,3));
    
  % iterate through all realizations of Halton set
  f           = waitbar(0,'Please be patient.','Name','Running simulations...');
  
  nM  = height(A)+height(B)+size(AB,1)*size(AB,3);
  for i = 1:height(A)
    parSet      = A(i,:);
    [Qex, Aex]  = runModel(shape,parSet);
    fAflux(i)   = Qex;
    fAarea(i)   = Aex;
    
    % show progress with waitbar
    waitbar(i/nM,f)
  end
  
  for i = 1:height(B)
    parSet      = B(i,:);    
    [Qex, Aex]  = runModel(shape,parSet);
    fBflux(i)   = Qex;
    fBarea(i)   = Aex;
    
    % show progress with waitbar
    waitbar((i+height(A))/nM,f)
  end
  
  for i = 1:size(AB,1)
    for j = 1:size(AB,3)
      parSet = array2table(AB(i,:,j),'VariableNames',A.Properties.VariableNames);
      
      [Qex, Aex]    = runModel(shape,parSet);
      fABflux(i,j)  = Qex;
      fABarea(i,j)  = Aex;
      
      % show progress with waitbar
      waitbar(((i-1)*size(AB,3)+j+height(A)+height(B))/nM,f)
    end
  end
  close(f)
  
  filename = sprintf('datasetSobol_%s.mat',shape);
  save(filename,'A','B','AB',...
                'fAflux','fBflux','fABflux',...
                'fAarea','fBarea','fABarea');
  
  path(oldPath)
end

function [Qex, Aex] = runModel(shape,p)
  % transform to get real input parameters
  L       = p.L;                        % domain length
  h1      = p.I.*p.L;                     % first fixed head
  h2      = 0.*h1;                % second fixed head
  Tx      = sqrt(p.Tmean.^2.*p.aniso);    % transmissivity in x
  Ty      = sqrt(p.Tmean.^2./p.aniso);    % transmissivity in y
  wMax    = p.ratMax.*p.L;                % maximum width
  wMin    = p.ratMin.*wMax;         % minimum width
  Q0      = (h1-h2)./p.L.*Tx.*(wMax-wMin); % normalization discharge
  qNorth  = p.QNstar.*Q0./p.L;        % northern influx (L³/L/T)
  
  mdl = fpAna('h1',h1,'h2',h2,'L',L,'wMin',wMin,...
    'wMax',wMax,'Tx',Tx,'Ty',Ty,...
    'qNorth',qNorth,'shape',shape,'autoSolve',false);
  
  % solve and read results
  [~]         = mdl.solve();
  
  try
    Qex       = mdl.Qex;
    if Qex > 1e-11
      Aex     = mdl.areaEx();
    else
      Aex     = 0;
    end
  catch
    Aex = 0;
  end
end

% NUMERICAL GRADIENTS
%   for j = 1:width(delQex)
%     temp        = parSet;
%     switch inT.Properties.VariableNames{j}
%       case "Tmean"
%         delta       = 1e-3*abs(log(temp{1,j}));
%         temp{1,j}   = exp(log(temp{1,j})+delta);
%       otherwise
%         delta       = 1e-3*temp{1,j};
%         temp{1,j}   = temp{1,j}+delta;
%     end
%
%     [QexN, AexN]  = runModel(shape,temp);
%
%     delQex{i,j} = (QexN-Qex)/delta;
%     delAex{i,j} = (AexN-Aex)/delta;
%   end




