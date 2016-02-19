%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DESCRIPTION: This script containes small pieces of script that can be
% used in the Main_ACHydro_Scenario.m script. These are just additional
% functionalities to analyze or calculate things, that are normally not
% requireed when running Main_ACHydro_Scenario.m
%
% Author: Hanne Van Gaelen
% Last updated: 15/02/2016
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% -----------------------------------------------------------------------
% EXTRA 1: DETERMINATION OF EXTENSION OF GDD 
%------------------------------------------------------------------------

% This script is used to determine how much the growing period could be
% extended in future climatic conditions

% 11.1 Median maturity and sowing date per group
% ----------------------------------------------------------------------- 

MATactstatsg(1,1:ncrop) =MATactstats(1,1:ncrop); 
MATpotstatsg(1,1:ncrop) =MATpotstats(1,1:ncrop); 
sowingstatsg(1,1:ncrop) =MATactstats(1,1:ncrop); 

for c=1:ncrop
    for g=1:ngroup2
        n=groupnames2(g);
        index=find(strcmp(groupmat,n)==1);
        MATactstatsg{2,c}(1,g)=median(MATactstats{2,c}(2,index));
        MATpotstatsg{2,c}(1,g)=median(MATpotstats{2,c}(2,index));
        sowingstatsg{2,c}(1,g)=sowing{2,c}(1,index(1)); % no need to calculate median because always the same
    end
end

clear c g index n

% 11.2 Calculate maturity target
% ----------------------------------------------------------------------- 
% TARGET = MEDIAN ACTUAL MATURITY DATE OF HISTORICAL TIMESERIES

maturitytarget(1,1:ncrop)=MATactstats(1,1:ncrop); 

indexhist=find(strcmp(groupnames2,'Hist')==1);
for c=1:ncrop
    targetdate=MATactstatsg{2,c}(1,indexhist);
    m=month(targetdate(1,1));
    d=day(targetdate(1,1));
    maturitytarget{2,c}=datetime(year(maturityact{2,c}(:,2)),m,d);
end

clear indexhist m d targetdate c

% 11.2 Calculate GDD between sowing and target
% ----------------------------------------------------------------------- 
indexfut=find(strcmp(groupmat,'Adapted')==1);
targetGDD=maturitytarget(1,1:ncrop); 

for sc=1:length(indexfut)
    for c=1:ncrop
        cropindex=find(strcmp(Crop(:,1),targetGDD{1,c})==1,1,'first'); %first landunit with this crop
        for y=1:length(maturitytarget{2,c})
            start=sowing{2,c}(y,indexfut(sc));
            startT=find(Date(:,indexfut(sc))==start);
            stop=maturitytarget{2,c}(y,1);
            stopT=find(Date(:,indexfut(sc))==stop);
            targetGDD{2,c}(y,indexfut(sc))=sum(GDD{2,indexfut(sc)}(startT:stopT,cropindex));
        end
    end
end

clear c sc y indexfut start startT stop stopT

%calculate stats
targetGDDstats(1,1:ncrop) =targetGDD(1,1:ncrop);
for c=1:ncrop
       targetGDDstats{2,c}(1,1:nsc)=nanmean(targetGDD{2,c}(:,1:nsc));
       targetGDDstats{2,c}(2,1:nsc)=nanmedian(targetGDD{2,c}(:,1:nsc));
       targetGDDstats{2,c}(3,1:nsc)=nanstd(targetGDD{2,c}(:,1:nsc));
       targetGDDstats{2,c}(4,1:nsc)=min(targetGDD{2,c}(:,1:nsc));
       targetGDDstats{2,c}(5,1:nsc)=max(targetGDD{2,c}(:,1:nsc));
        
end

clear c 

%calculate stats per group
targetGDDstatsg(1,1:ncrop) =targetGDD(1,1:ncrop); 

for c=1:ncrop
    for g=1:ngroup2
        n=groupnames2(g);
        index=find(strcmp(groupmat,n)==1);
        targetGDDstatsg{2,c}(1,g)=median(targetGDDstats{2,c}(2,index));
    end
end

clear c n g index

%% -----------------------------------------------------------------------
% EXTRA 2. FIND YEAR OF MEDIAN CROP YIELD 
%------------------------------------------------------------------------
% This script can be used to find out in which year the median crop yield
% was obtained


name={'Maize','WinterWheat','Potato','Sugarbeet','Pea'};

for sc=1:nsc
year2(:,sc)=unique(Year(:,sc));
end

for c=1:length(name)
    cindex=find(strcmp(Cropnames(1,1:ncrop),name{1,c})==1);
    medianyear{1,c}=name{1,c};
    
    med=Yactstats{2,cindex}(2,1:nsc);
    [a,~]=size(Yact{2,cindex}(:,1:nsc));
    med=repmat(med,a,1);
    dif=abs(Yact{2,cindex}(:,1:nsc)-med);
    
    for sc=1:nsc
    y=find(dif(:,sc)==min(dif(:,sc)));
    
    medianyear{2,c}(:,sc)=min(year2(y,sc));
    end
    
    clear med dif a cindex sc y 
end

clear c year2 name 