% 
% [AUX,Wflag] = chkaux(AUXin,i,u(1,1),y(1,1),flag,W,sil)
% 
% Description:
%          Compatibility of AUXin check
%          AUX = [] when nothing is useful
%          Wflag = 1 when R information is OK, but weight information not
%          
% References:     
%          None
%
% Copyright: 
%          Peter Van Overschee, December 1995
%          peter.vanoverschee@esat.kuleuven.ac.be
%
%

function [AUX,Wflag] = chkaux(AUXin,i,u,y,flag,W,sil);

AUX = AUXin;
Wflag = [];

if ~isempty(AUXin);
  info = AUXin(1,:);Wflag = 0;
  if (info(1) ~= flag);
    if (flag == 1)
      mydisp(sil,'      Warning: AUXin is neglected: only valid for stochastic models');
    end
    if (flag == 2)  
      mydisp(sil,'      Warning: AUXin is neglected: only valid for deterministic models');
    end   
    AUX = [];
  end
  if (info(2) ~= i);
    mydisp(sil,'      Warning: AUXin is neglected: Incompatible i');
    AUX = [];
  end
  if ~isempty(u)
    if (info(3) ~= u);
      mydisp(sil,'      Warning: AUXin is neglected: Incompatible input');
      AUX = [];
    end  
  end
  if (info(4) ~= y);
    mydisp(sil,'      Warning: AUXin is neglected: Incompatible output');
    AUX = [];
  end  
  if ~isempty(W) & (W ~= 0) 
    if (info(5) ~= W);
      mydisp(sil,'      Warning: Weighting part in AUXin is neglected: Incompatible weight');
      Wflag = 1;
    end  
  end
end


