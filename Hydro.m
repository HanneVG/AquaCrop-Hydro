% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script performs calculations of catchment hydrology
% based on inputs that are calculated by AquaCrop, and parameters that are
% specified in an input txt file
%
%
% Author: Hanne Van Gaelen
% Last update: 1/08/2015
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Q_MBF,Q_MIF,Q_MOF,Q_MTF,area,f,Wrmin,Wrmax,pbf]=Hydro(Par,SoilPar,nAC,ROCatch,DPCatch,Wr2Catch)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 1. DEFINE ALL INPUT PARAMETERS AND VARIABLES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    n=nAC; % number of simulated timesteps by Aquacrop
    
    kbf=Par(1,1); % recession constant for baseflow
    kif=Par(2,1);% recession constant for interflow
    kof=Par(3,1);% recession constant for overland flow
    
    qinibf=Par(1,2); % initial discharge for baseflow
    qiniif=Par(2,2);% initial discharge for interflow
    qiniof=Par(3,2);% initial discharge for overland flow
    
    cinibf=Par(1,3);% initial contribution of baseflow
    ciniif=Par(2,3);% initial contribution of interflow
    ciniof=Par(3,3);% initial contribution of overland flow
    
    ResOF=Par(1,5); % number of parallell reservoirs to rout overland flow (max 3)
    
    WrFC=2*SoilPar(1,3); % Soil water content at field capacity for 2 m soil depth
    
    pbfmin=Par(1,4); % min fraction of deep percolation that goes to baseflow
    pbfmax=Par(2,4); % max fraction of deep percolation that goes to baseflow
    Wrmax=max(Wr2Catch);% upper threshold is highest SWC value found in sim period
    Wrmin=WrFC;% lower threshold is field capacity
        
    area=Par(1,6); % area of catchment (km2)
    f =Par(1,7); % number of seconds *10^3 within one timestep
     
    DP=DPCatch; % Deep percolation for these calculations is the one of the whole catchment
    RO=ROCatch; % Runoff for these calculations is the one of the whole catchment
    
 
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 2. SPLIT DEEP PERCOLATION INTO INTERFLOW AND BASEFLOW FRACTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Predefine variables 
pbf=NaN(n,1);
DPif=NaN(n,1);
DPbf=NaN(n,1);

% Linair function for split DP into IF and BF
for i = 1:n
    if Wr2Catch(i)<Wrmin
        pbf(i)=pbfmin;
    elseif Wr2Catch(i)>Wrmax
        pbf(i)=pbfmax;
    else
    pbf(i)=pbfmin + ((pbfmax-pbfmin)/(Wrmax-Wrmin))*(Wr2Catch(i)-Wrmin);% derive pbf on that day depending on soil water content with lineair relation
    end
    
    DPif(i)=DP(i)*(1-pbf(i));%derive contributions of DP to baseflow and interflow
    DPbf(i)=DP(i)*pbf(i);
end
   

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 3. ROUTING MODEL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% Predefine variables 
    MBF=zeros(n,1);
    MIF=zeros(n,1);
    MOF=zeros(n,1);
    MOF1=zeros(n,1);
    MOF2=zeros(n,1);
    MOF3=zeros(n,1);   
    
% Initial values for each subflow
    MBF(1)=((exp(-1/kbf))*qinibf)+((1-exp(-1/kbf))*0.5*(cinibf+DPbf(1)));   
    MIF(1)=((exp(-1/kif))*qiniif)+((1-exp(-1/kif))*0.5*(ciniif+DPif(1)));     
    MOF1(1)=((exp(-1/kof))*qiniof)+((1-exp(-1/kof))*0.5*(ciniof+RO(1)));
    MOF2(1)=0;
    MOF3(1)=0;     
        
% Other timesteps for each subflow using lineair reservoir functions
for i=2:n
     MBF(i)=(exp(-1/kbf)*MBF(i-1))+((1-exp(-1/kbf))*0.5*(DPbf(i-1)+DPbf(i))); 
     MIF(i)=(exp(-1/kif)*MIF(i-1))+((1-exp(-1/kif))*0.5*(DPif(i-1)+DPif(i))); 
        
      if ResOF==1
       MOF1(i)=(exp(-1/kof)*MOF1(i-1))+((1-exp(-1/kof))*0.5*(RO(i-1)+RO(i)));
       clear MOF2 MOF3
       MOF(i)=MOF1(i);
      elseif ResOF==2
       MOF1(i)=(exp(-1/kof)*MOF1(i-1))+((1-exp(-1/kof))*0.5*(RO(i-1)+RO(i))); 
       MOF2(i)=(exp(-1/kof)*MOF(i-1))+((1-exp(-1/kof))*0.5*(MOF1(i-1)+MOF1(i)));
       clear MOF3
       MOF(i)=MOF2(i);
      else
        MOF1(i)=(exp(-1/kof)*MOF1(i-1))+((1-exp(-1/kof))*0.5*(RO(i-1)+RO(i))); 
        MOF2(i)=(exp(-1/kof)*MOF(i-1))+((1-exp(-1/kof))*0.5*(MOF1(i-1)+MOF1(i)));
        MOF3(i)=(exp(-1/kof)*MOF(i-1))+((1-exp(-1/kof))*0.5*(MOF2(i-1)+MOF2(i)));
        MOF(i)=MOF3(i);
      end
            
end
                 
% Conversion from mm/timestep to m³/s 
    Q_MBF=MBF*(area/f);
    Q_MIF=MIF*(area/f);
    Q_MOF=MOF*(area/f);
    Q_MTF=Q_MBF+Q_MIF+Q_MOF;
    
   
end