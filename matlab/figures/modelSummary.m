function [out] = modelSummary(mdl)

  switch mdl.shape
    case "bump"
      a1 = 5.852;
      a2 = 0.355;
      a3 = 4.607;
    case "cosinusoidal"
      a1 = 6.242;
      a2 = 0.434;
      a3 = 4.121;
    case "composite"
      a1 = 5.515;
      a2 = 0.331;
      a3 = 4.755;
  end
  
Anorth  = mdl.area-mdl.L.*mdl.wMin;
wMean   = mdl.area./mdl.L;
Q0      = abs(mdl.h1-mdl.h2)/mdl.L*mdl.Tx*(mdl.wMax-mdl.wMin);
xT      = sqrt(mdl.Tx./mdl.Ty)*wMean./mdl.L;
QnT     = mdl.L*mdl.qNorth./Q0;
QexT    = sech(a1*xT)*max(0,1-a2*abs(QnT)*cosh(a3*xT));
AexT    = QexT/sqrt(1+abs(QnT));
Aex     = AexT*Anorth;
Qex     = QexT*Q0;

fprintf('wMmean       = %.3gm\n',wMean);
fprintf('A_north      = %.3gm²\n',Anorth);
fprintf('Q_0          = %.3em^³/s\n',Q0);
fprintf('tildeX       = %.3g\n',xT);
fprintf('tildeQ_north = %.3g\n',QnT);
fprintf('tildeQ_ex    = %.3g\n',QexT);
fprintf('Q_ex         = %.3em³/s\n',Qex);
fprintf('tildeA_ex    = %.3g \n',AexT);
fprintf('A_ex         = %.3gm²\n',Aex);
fprintf('tMean        = %.3es approx %.1fa\n',Aex*mdl.por/Qex,years(seconds(Aex*mdl.por/Qex)));
fprintf('realQ_ex     = %.3em³/s\n',mdl.Qex);
fprintf('realA_ex     = %.3gm²\n',mdl.areaEx);
fprintf('realtMean    = %.1fa\n',years(seconds(mdl.areaEx*mdl.por/mdl.Qex)));

out = 1;

end