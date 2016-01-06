%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DESCRIPTION: This script compares different scenarios that are simulated with
% AquaCrop-Hydro. It runs each scenario and compares crop yield, soil water balance and river
% discharge (cumulative volumes) for every scenario
%
% TO DO before running the script: 
%   1. Run AquaCrop simulations for all land units in the catchment for
%      each scenario
%   2. Ensure that all AquaCrop output files are numbered with format "01", "02",
%      "03", "10","20" according to the landunit number
%   3. Ensure that all AquaCrop ouput files are organized in subfolders
%      with the name of the scenario
%   4. Prepare all required input files for this script 
% 
%   
%
% ASSUMPTIONS: 
%   1. Affecting results of final yield:
%       - Even simulation runs are main season, odd simulation runs are after
%       season crops (see line X)
%   2. the first scnenario is always the baseline scenario to which all scenarios are compared 
%
%
%
% Author: Hanne Van Gaelen
% Last updated: 22/12/2015
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ------------------------------------------------------------------------
% 1. DEFINE DATA PATHS (BY USER)                                                    
%--------------------------------------------------------------------------      

% specify in which path all AquaCrop simulation output files are stored for
% all scenarios (each scenario in different subfolder)

     DatapathAC=uigetdir('C:\','Select directory of scenario subfloders with AquaCrop output files for all landunits');     
     
% specify in which path all inputs for AquaCrop-Hydro are stored including 
%    a) additional information on AquaCrop simulations of each landunits (SimInfo.txt)
%    b) the parameters for the hydrological model (Parameters.txt)
%    c) the maximum root depth for every landunit and sim run (Zrx.txt)
%    d) the soil parameters of each soil type present (SoilPar.txt)

     DatapathACHinput = uigetdir('C:\','Select directory with all input files for AquaCrop-Hydro ');
     
% specify in which path all output should be stored including 
%   a) Flow values (baseflow, interflow, overland flow, total flow) as
%   simulated by AquaCrop-Hydro
%   b) Yield values as simulated for each landunit during the main cropping season of each simulated year

     DatapathACHOutput=uigetdir('C:\','Select directory to store AquaCrop-Hydro output (flows and yield values)');        

% specify in which path all information on the scenario comparison is stored including  
 %    a) information on the different scenarios (Scenario.txt)    
     
     DatapathScenIn = uigetdir('C:\','Select directory with all information on the scenario comparison ');     

% specify in which path all output of scenario comparison (figure and workspace variables) should be stored   
     
     DatapathScenOut = uigetdir('C:\','Select directory to store scenario comparison output (figures and variables) ');     

% specify in which path all inputs for date calculations are stored (each scenario in different subfolder) including *Temp.Txt and *PrChar.txt) 

     DatapathDate= uigetdir('C:\','Select directory of scenario subfloders with all input files for date calculations');

% specify the AquaCrop mode that was used (1= normal AquaCrop, 2= plugin version)
     ACMode=inputdlg('Did you use AquaCrop normal (1) or stand-alone plugin version (2)?','AquaCrop Mode');
     ACMode=cell2mat(ACMode);
     ACMode=str2double(ACMode);

     if ACMode==1 || ACMode==2
         %continue
     else
         error('invalid AquaCrop mode selected');
     end
 
     
%% ------------------------------------------------------------------------
% 2. LOAD EXTRA INFORMATION ON SIIMULATIONS AND LANDUNITS                                                  
%--------------------------------------------------------------------------        
        
% 2.1 Load information on different landunits
%-------------------------------------------------------------------------    
 
 % Select file with info on simulations
    name='SimInfo.txt';
    file = fullfile(DatapathACHinput, name);        
    A= importdata(file); 
    clear name file
    
 % Read all data of each simulation unit
    SimName=A.textdata(:,1);    % Name 
    SimNr=A.data(:,2);          % Numeric code of this simulation unit 
    SimType=A.data(:,3);        % Type of this simulation unit (1= other landuse (e.g forest, water, urban), 2= agricultural landuse, 999= no simulation)
    Clim=A.textdata(:,2);       % Climate 
    Soil=A.textdata(:,3);       % Soil type
    Crop=A.textdata(:,4);       % Crop type grown in main season
    CropAfter=A.textdata(:,5);  % Crop type grown after main season crop
    CropRot=A.textdata(:,6);    % Crop rotation (= main crop + after crop)
    SimArea=A.data(:,1);        % Relative area of this simulation unit in the catchment
    [nlu,~]=size(SimType(SimType(:,1)<900));   % real land units
    clear A

% 2.2 Load information on different scenarios
%-------------------------------------------------------------------------   
    name='Scenario.txt';
    file = fullfile(DatapathScenIn, name);        
    A= importdata(file); 
    clear name file   
    
    scnumb=A.data(:,1).'; 
    nsc=max(scnumb);
    ScenarioName=A.textdata(:,1).';
    ClimModel=A.textdata(:,2).';
    RCP=A.textdata(:,3).';
    StartDate=datetime((A.data(:,2).'),'ConvertFrom','excel');
    EndDate=datetime((A.data(:,3).'),'ConvertFrom','excel');
    nTime=daysact(datenum(StartDate(1,1)), datenum(EndDate(1,1)))+1;
    clear A 
    
%% -----------------------------------------------------------------------
% 3. RUN AQUACROP-HYDRO FOR ALL SCENARIOS & SAVE OUTPUT OF EACH SCENARIO
%-------------------------------------------------------------------------

% % 3.1 Check if AquaCrop results for all scenarios can be found
% % -----------------------------------------------------------------------
% Folders_DatapathAC=dir(DatapathAC);    % Folder structure of scenarios
% 
% %clean up
%  for k = length(Folders_DatapathAC):-1:1 
%     % remove non-folders
%     if ~Folders_DatapathAC(k).isdir
%         Folders_DatapathAC(k) = [ ];
%         continue
%     end
% 
%     % remove folders starting with .
%     fname = Folders_DatapathAC(k).name;
%     if fname(1) == '.'
%         Folders_DatapathAC(k) = [ ];
%     end
% end
% 
% % order scenarios
%     %convert to matrix      
%     FolderNames = fieldnames(Folders_DatapathAC);
% 
%     FoldersMat=struct2cell(Folders_DatapathAC);
%     FoldersMat=FoldersMat';
% 
%     %sort
%     FoldersMatSort=sortrows(FoldersMat,2);   % use date timestamp to sort files (=proxy for order of scenarios) 
% 
%     %convert back to structure
%     Folders_DatapathACSort=cell2struct(FoldersMatSort,FolderNames,2);
% 
%     %replace old datastructure with sorted
%     Folders_DatapathAC=Folders_DatapathACSort;      
% 
% 
% % check relevant information
% [nsc2,~]=size(Folders_DatapathAC);      % Number of scenarios
% if nsc2==nsc
%     % continue as the number of simulated scenarios is correct
% else 
%     error('the number of simulated scenarios does not match the number of scenarios in information file')
% end
% 
% 
% ScenarioName2= {Folders_DatapathAC.name};
% for sc=1:nsc
%     if strcmp(ScenarioName2(1,sc),ScenarioName(1,sc))
%         %continue as the names of the simulated scenarios are correct 
%     else
%         error(['The scenarioname of simulated scenario ',num2str(sc),' does not match the scenario name in the information file'])
%     end
% end
% 
% clear FoldersMatSortt FoldersMat Folders_DatapathACSort FolderNames fname k Folders_DatapathAC nsc2 ScenarioName2

% 3.2 Run all scenarios and save results
% -----------------------------------------------------------------------             

%initialize variables
    %time variables
    Day=NaN(nTime,nsc);    % Day number
    Month=NaN(nTime,nsc);  % Month number
    Year=NaN(nTime,nsc);   % Year number

    %catchment-scale results
    TrCatch=NaN(nTime,nsc);     % Crop transpiration (actual) (mm)
    TrxCatch=NaN(nTime,nsc);    % Potential (maximum) crop transpiration (mm)
    ECatch=NaN(nTime,nsc);      % Evaporation(actual)(mm)
    ExCatch=NaN(nTime,nsc);     % Potential (maximum) evaporation(mm)
    ETaCatch=NaN(nTime,nsc);    % Actual evapotranspiration(mm)
    ETxCatch=NaN(nTime,nsc);    % Potential (maximum) evapotranspiration(mm)
    ROCatch=NaN(nTime,nsc);     % Runoff(mm)
    DPCatch=NaN(nTime,nsc);     % Deep percolation(mm)
    CRCatch=NaN(nTime,nsc);     % Capilary rise(mm)
    BundWatCatch=NaN(nTime,nsc); % Water between bunds (mm)
    Wr2Catch=NaN(nTime,nsc);   % Soil water content in 2 m soil depth (mm)
    CCCatch=NaN(nTime,nsc);    % Canopy Cover (%)

   % AquaCrop results per simulation unit 
    Tr=cell(2,nsc) ;    % Crop transpiration(actual)(mm)
    Trx=cell(2,nsc) ;   % Potential (maximum) crop transpiration(mm)
    E=cell(2,nsc) ;     % Evaporation (actual)(mm)
    Ex=cell(2,nsc) ;    % Potential (maximum) evaporation(mm)
    ETa=cell(2,nsc) ;   % Actual evapotranspiration(mm)
    ETx=cell(2,nsc) ;   % Potential (maximum) evapotranspiration(mm)
    RO=cell(2,nsc) ;    % Runoff(mm)
    DP=cell(2,nsc) ;  % Deep percolation(mm)
    CR=cell(2,nsc) ;    % Capilary rise(mm)
    BundWat=cell(2,nsc) ; % Water between bunds (mm)
    Wr2=cell(2,nsc) ; % Soil water content in 2 m soil depth (mm)
    CC=cell(2,nsc) ;   % Canopy Cover (%)
    B=cell(2,nsc) ;    % Dry aboveground biomass during growing season (ton/ha)
    Bfin=cell(2,nsc) ; % Final dry aboveground biomass at maturity (ton/ha)
    Y=cell(2,nsc) ;    % Dry yield at maturity for each sim run season(ton/ha)
    Ymain=cell(2,nsc) ;     % Dry yield at maturity for only main season(ton/ha)(even sim runs correspond to the main season, odd run numbers are after-season crops)   

   % the catchment hydrology results   
    Q_MBF=NaN(nTime,nsc);
    Q_MIF=NaN(nTime,nsc);
    Q_MOF=NaN(nTime,nsc);      
    Q_MTF=NaN(nTime,nsc);

  % Add headers to matrices 
    Tr(1,1:nsc)= ScenarioName(1,1:nsc);  
    Trx(1,1:nsc)= ScenarioName(1,1:nsc);  
    E(1,1:nsc)= ScenarioName(1,1:nsc);  
    Ex(1,1:nsc)= ScenarioName(1,1:nsc);  
    ETa(1,1:nsc)=  ScenarioName(1,1:nsc);  
    ETx(1,1:nsc)=ScenarioName(1,1:nsc);   
    RO(1,1:nsc)= ScenarioName(1,1:nsc);  
    DP(1,1:nsc)=ScenarioName(1,1:nsc);         
    CR(1,1:nsc)=  ScenarioName(1,1:nsc);  
    BundWat(1,1:nsc)=ScenarioName(1,1:nsc);   
    Wr2(1,1:nsc)=   ScenarioName(1,1:nsc);   
    CC(1,1:nsc)=  ScenarioName(1,1:nsc);    
    B(1,1:nsc)= ScenarioName(1,1:nsc);  
    Bfin(1,1:nsc)= ScenarioName(1,1:nsc);  
    Y(1,1:nsc)=ScenarioName(1,1:nsc);  
    Ymain(1,1:nsc)=ScenarioName(1,1:nsc);        
        
for sc=1:nsc %loop trough all scenarios
% show progress
disp(['Now running scenario ',num2str(sc),' of the ',num2str(nsc),' scenarios']);
 
% Extract scenario name
Name=ScenarioName{1,sc};

% Datapaths for this scenario
  DatapathACSC=fullfile(DatapathAC,Name);
  DatapathInputSC=DatapathACHinput;                % input is the same for all scenarios   
  DatapathOutputSC=fullfile(DatapathACHOutput,Name);  

% Check if AquaCrop results for this scenario can be found  
    
    if exist(DatapathACSC) == 7
        %continue as the AquaCrop results for a scneario with this name
        can be found
    else
        error(['The AquaCrop results for scenario ',num2str(sc),' with scenarioname ',Name,' could not be found'])
    end

% Run AquaCrop-Hydro for this scenario
[Q_MBFsc,Q_MIFsc,Q_MOFsc,Q_MTFsc,area,f,Wrmin,Wrmax,pbf,SoilPar,SimACOutput,CatchACOutput,Par]=AquaCropHydro(DatapathACSC, DatapathInputSC,ACMode);

% extract output for this scenario       

        % Save time variables        
        Day(:,sc)=CatchACOutput(:,13);    % Day number
        Month(:,sc)=CatchACOutput(:,14);  % Month number
        Year(:,sc)=CatchACOutput(:,15);   % Year number
        Date(:,sc)=datetime(Year(:,sc),Month(:,sc),Day(:,sc)); % Date 
        
        %Check number of timesteps
        nt=length(CatchACOutput(:,1));   % Number of timesteps
        if nt==nTime 
         %continue and use nTime as timemarker
        else
         error('number of simulated days does not match number of days in specified baseline simulation period');
        end
               
        % Save catchment-scale AquaCrop results
        TrCatch(1:nTime,sc)=CatchACOutput(1:nTime,1);     % Crop transpiration (actual) (mm)
        TrxCatch(1:nTime,sc)=CatchACOutput(1:nTime,2);    % Potential (maximum) crop transpiration (mm)
        ECatch(1:nTime,sc)=CatchACOutput(1:nTime,3);      % Evaporation(actual)(mm)
        ExCatch(1:nTime,sc)=CatchACOutput(1:nTime,4);     % Potential (maximum) evaporation(mm)
        ETaCatch(1:nTime,sc)=CatchACOutput(1:nTime,5);    % Actual evapotranspiration(mm)
        ETxCatch(1:nTime,sc)=CatchACOutput(1:nTime,6);    % Potential (maximum) evapotranspiration(mm)
        ROCatch(1:nTime,sc)=CatchACOutput(1:nTime,7);     % Runoff(mm)
        DPCatch(1:nTime,sc)=CatchACOutput(1:nTime,8);     % Deep percolation(mm)
        CRCatch(1:nTime,sc)=CatchACOutput(1:nTime,9);     % Capilary rise(mm)
        BundWatCatch(1:nTime,sc)=CatchACOutput(1:nTime,10); % Water between bunds (mm)
        Wr2Catch(1:nTime,sc)=CatchACOutput(1:nTime,11);   % Soil water content in 2 m soil depth (mm)
        CCCatch(1:nTime,sc)=CatchACOutput(1:nTime,12);    % Canopy Cover (%)
                
       % Save original AquaCrop results per simulation unit 
        Tr{2,sc}=SimACOutput{1,1};    % Crop transpiration(actual)(mm)
        Trx{2,sc}=SimACOutput{1,2};   % Potential (maximum) crop transpiration(mm)
        E{2,sc}=SimACOutput{1,3};     % Evaporation (actual)(mm)
        Ex{2,sc}=SimACOutput{1,4};    % Potential (maximum) evaporation(mm)
        ETa{2,sc}=SimACOutput{1,5};   % Actual evapotranspiration(mm)
        ETx{2,sc}=SimACOutput{1,6};   % Potential (maximum) evapotranspiration(mm)
        RO{2,sc}=SimACOutput{1,7};    % Runoff(mm)
        DP{2,sc}=SimACOutput{1,8};    % Deep percolation(mm)
        CR{2,sc}=SimACOutput{1,9};    % Capilary rise(mm)
        BundWat{2,sc}=SimACOutput{1,10}; % Water between bunds (mm)
        Wr2{2,sc}=SimACOutput{1,11};  % Soil water content in 2 m soil depth (mm)
        CC{2,sc}=SimACOutput{1,12};   % Canopy Cover (%)
        B{2,sc}=SimACOutput{1,13};    % Dry aboveground biomass during growing season (ton/ha)
        Bfin{2,sc}=SimACOutput{1,14}; % Final dry aboveground biomass at maturity (ton/ha)
        Y{2,sc}=SimACOutput{1,15};    % Dry yield at maturity for each sim run season(ton/ha)
        Ymain{2,sc}=Y{2,sc}(2:2:end,:);     % Dry yield at maturity for only main season(ton/ha)(even sim runs correspond to the main season, odd run numbers are after-season crops)                  
        
      % Save the catchment hydrology results   
        Q_MBF(1:nTime,sc)=Q_MBFsc(1:nTime,1);
        Q_MIF(1:nTime,sc)=Q_MIFsc(1:nTime,1);
        Q_MOF(1:nTime,sc)=Q_MOFsc(1:nTime,1);      
        Q_MTF(1:nTime,sc)=Q_MTFsc(1:nTime,1);
 clear Q_MTFsc Q_MOFsc Q_MIFsc QMBFsc
 
% write output for this scenario to excel      
      % Combine output in matrix if necessary
        HeadersFlow={'Date','Baseflow','Interflow','Overland flow', 'Total flow'};
        FlowOutput=[exceltime(Date(1:nTime,1)),Q_MBFsc(1:nTime,1),Q_MIFsc(1:nTime,1),Q_MOFsc(1:nTime,1),Q_MTFsc(1:nTime,1)];

        HeadersWabalCatch={'Date','Tr','Trx','E','Ex','ETa','ETx','RO','DP','CR','BundWat','Wr2','CC'};

       % Write output to one excel tabsheet
        xlname='FlowSimResults.xlsx';
        filename = fullfile(DatapathOutputSC,xlname);
        xlswrite(filename,HeadersFlow,'SimFlow','A1');
        xlswrite(filename,FlowOutput,'SimFlow','A2');
    
        xlname='WabalSimResults.xlsx';
        filename = fullfile(DatapathOutputSC,xlname);
        xlswrite(filename,HeadersWabalCatch,'SimWabal','A1');
        xlswrite(filename,exceltime(Date(1:nTime,1)),'SimWabal','A2');
        xlswrite(filename,CatchACOutput(1:nTime,1:12),'SimWabal','B2'); 
        
        xlname='YieldResults.xlsx';
        filename = fullfile(DatapathOutputSC,xlname);
        xlswrite(filename,Ymain{2,sc},'SimYield','A1'); 
        
clear SimACOutput CatchOutput FlowOutput WabalCatchOutput YieldOutput   DatapathOutputSC DatapathInputSC   
end

clear sc

% save workspace variables so that you can skip this the first part of the
% code next time
filename=['workspace ',datestr(date)];
filename=fullfile(DatapathScenOutput,filename);
save(filename)
clear filename 
 

%% -----------------------------------------------------------------------
% 4. YIELD IMPACT 
%------------------------------------------------------------------------

% 4.1 Compose yield matrix per crop type
%---------------------------------------
Crop2 = unique(Crop(:,1)); % number of unique crops in catchment
Crop2=Crop2(strcmp(Crop2(:,1),'Unknown')==0); % remove the landunits with unknown crop 
Crop2=Crop2(strcmp(Crop2(:,1),'Impervious')==0); % remove the landunits with unknown crop 
Crop2=Crop2(strcmp(Crop2(:,1),'Water')==0); % remove the landunits with unknown crop 
[ncrop,~]=size(Crop2); % number of real crops

Yall=cell(2,ncrop);% initialize

for c=1:ncrop% loop trough each crop and make subset of crop yields 
    % write away crop name
    Yall(1,c)=Crop2(c); 
    
    %search all projects with this crop
    index=find(strcmp(Crop(:,1),Crop2(c))==1);
   
    for sc=1:nsc %loop trough al scenarios
        ysub=Ymain{2,sc}(:,index); % write data away
        Yall{2,c}(:,sc)=mean(ysub.'); % take average of all data of same crop
    end        
end
clear sc c 

% 4.2 Calculate statistics
%---------------------------------------
% stats for each crop over different years 
    Ystats(1,1:ncrop)=Yall(1,1:ncrop); 
    
    for c=1:ncrop
        Ystats{2,c}(1,1:nsc)=mean(Yall{2,c}(:,1:nsc));
        Ystats{2,c}(2,1:nsc)=median(Yall{2,c}(:,1:nsc));
        Ystats{2,c}(3,1:nsc)=std(Yall{2,c}(:,1:nsc));
        Ystats{2,c}(4,1:nsc)=min(Yall{2,c}(:,1:nsc));
        Ystats{2,c}(5,1:nsc)=max(Yall{2,c}(:,1:nsc));
    end
clear c

% 4.3 Mean changes of statistics
%---------------------------------------
% Calculate changes of stats (change of avg, change of median)
    YDeltastats(1,1:ncrop)=Yall(1,1:ncrop); 

    for c=1:ncrop
        for stat=1:2
        YDeltastats{2,c}(stat,1:nsc)=(Ystats{2,c}(stat,1:nsc)-Ystats{2,c}(stat,1))./Ystats{2,c}(stat,1);
        end
    end

clear c
   
% 4.4 Vizualize yield impact
%--------------------------------------- 
NameMat= {'RCP 4.5','RCP 8.5'};

%search indices of crops you want to show
maize=find(strcmp(Yall(1,:),'Maize')==1);
wwheat=find(strcmp(Yall(1,:),'WinterWheat')==1);
sugarbeet=find(strcmp(Yall(1,:),'Sugarbeet')==1);
potato=find(strcmp(Yall(1,:),'Potato')==1);
pea=find(strcmp(Yall(1,:),'Pea')==1);

figure('name','Median yield changes')
        sub(1)=subplot(2,5,1,'fontsize',10);
        boxplot(YDeltastats{2,maize}(2,2:16)*100,RCP(1,2:16),'grouporder',{'RCP4.5','RCP8.5'},'labels',NameMat,'colors',['r','b']);
        line(xlim,[0,0],'Color','k','LineStyle','--')
        ylabel('Median yield change (%)')
        axis([xlim, -30,30])
        set(gca,'box','off')
        title('Maize')
        
        sub(2)=subplot(2,5,2,'fontsize',10);
        boxplot(YDeltastats{2,wwheat}(2,2:16)*100,RCP(1,2:16),'grouporder',{'RCP4.5','RCP8.5'},'labels',NameMat);
        line(xlim,[0,0],'Color','k','LineStyle','--')
        set(gca,'box','off')
        title('Winter Wheat')
        
        sub(3)=subplot(2,5,3,'fontsize',10);
        boxplot(YDeltastats{2,potato}(2,2:16)*100,RCP(1,2:16),'grouporder',{'RCP4.5','RCP8.5'},'labels',NameMat);
        line(xlim,[0,0],'Color','k','LineStyle','--')
        set(gca,'box','off')
        title('Potato')
        
        sub(4)=subplot(2,5,4,'fontsize',10);
        boxplot(YDeltastats{2,sugarbeet}(2,2:16)*100,RCP(1,2:16),'grouporder',{'RCP4.5','RCP8.5'},'labels',NameMat);
        line(xlim,[0,0],'Color','k','LineStyle','--')
        set(gca,'box','off')
        title('Sugarbeet')
        
        sub(5)=subplot(2,5,5,'fontsize',10);
        boxplot(YDeltastats{2,pea}(2,2:16)*100,RCP(1,2:16),'grouporder',{'RCP4.5','RCP8.5'},'labels',NameMat);
        line(xlim,[0,0],'Color','k','LineStyle','--')
        set(gca,'box','off')
        title('Peas')    

        linkaxes(sub,'y')% link y axis of different plots (so that they change simultaneously
        
        h(1)=subplot(2,5,6,'fontsize',10);
        bar(Ystats{2,maize}(2,1))
        ylabel('Median historical yield (ton/ha)')
        axis([xlim , 0 ,15])
        title('Maize')
        set(gca,'XTickLabel',{' '})
        
        h(2)=subplot(2,5,7,'fontsize',10);
        bar(Ystats{2,wwheat}(2,1))
        title('Winter Wheat')
        set(gca,'XTickLabel',{' '})
        
        h(3)=subplot(2,5,8,'fontsize',10);
        bar(Ystats{2,potato}(2,1))
        title('Potato')
        set(gca,'XTickLabel',{' '})
         
        h(4)=subplot(2,5,9,'fontsize',10);
        bar(Ystats{2,sugarbeet}(2,1))
        title('Sugarbeet')
        set(gca,'XTickLabel',{' '})
                
        h(5)=subplot(2,5,10,'fontsize',10);
        bar(Ystats{2,pea}(2,1))
        title('Peas')
        set(gca,'XTickLabel',{' '})
        
        linkaxes(h,'y')% link y axis of different plots (so that they change simultaneously
 
clear sub h  

%% -----------------------------------------------------------------------
% 5. IMPACT ON LENGTH OF THE GROWING CYCLE 
%-------------------------------------------------------------------------

% 5.1 Calculate cycle lengths
%-------------------------------------------------------------------------
CycleLength(1,1:nsc)=ScenarioName(1,1:nsc);

for sc=1:nsc
    DatapathDateSC=fullfile(DatapathDate,ScenarioName{1,sc});
    DatapathOutputSC=fullfile(DatapathACHOutput,ScenarioName{1,sc});
    RotationDateSC=CalcDate(DatapathDateSC,DatapathOutputSC,StartDate(1,sc),EndDate(1,sc)); 
    
    for lu=1:nlu
        if SimType(lu,1)==2
            nrun=length(RotationDateSC{1,lu}(:,3));
            r=(nrun-rem(nrun,2))/2;     
            SowingDateSC(1:r,lu)= RotationDateSC{1,lu}(2:2:nrun,3);   
            MaturityDateSC(1:r,lu)= RotationDateSC{1,lu}(2:2:nrun,4);
        else 
            % skip this landunit (not agriculture)
        end
    end
    
    LengthSC=(datenum(MaturityDateSC)-datenum(SowingDateSC))+1;
    
    CycleLength{2,sc}=LengthSC; % per scenario one matrix with cycle lengths per landunit
    
end

clear sc SowingDateSC MaturityDateSC LengthSC

% 5.1 Compose cycle length matrix per crop type
%--------------------------------------------------------------------------
CycleLengthCrop=cell(2,ncrop);% initialize

for c=1:ncrop% loop trough each crop 
    % write away crop name
    CycleLengthCrop(1,c)=Crop2(c); 
    
    %search all projects with this crop
    index=find(strcmp(Crop(:,1),Crop2(c))==1);
   
    for sc=1:nsc %loop trough al scenarios
        subset= CycleLength{2,sc}(:,index);
        CycleLengthCrop{2,c}(:,sc)=mean(subset.'); % take average of all data of same crop
    end        
end

clear sc c subset 

% 5.3 Calculate statistics (over different year)
%-------------------------------------------------------------------------
 CLengthCropstats(1,1:ncrop)=CycleLengthCrop(1,1:ncrop); 
    
    for c=1:ncrop
        CLengthCropstats{2,c}(1,1:nsc)=nanmean(CycleLengthCrop{2,c}(:,1:nsc));
        CLengthCropstats{2,c}(2,1:nsc)=nanmedian(CycleLengthCrop{2,c}(:,1:nsc));
        CLengthCropstats{2,c}(3,1:nsc)=nanstd(CycleLengthCrop{2,c}(:,1:nsc));
        CLengthCropstats{2,c}(4,1:nsc)=min(CycleLengthCrop{2,c}(:,1:nsc));
        CLengthCropstats{2,c}(5,1:nsc)=max(CycleLengthCrop{2,c}(:,1:nsc));
    end
clear c

% 5.4 Mean changes of statistics
%-------------------------------------------------------------------------
% Calculate changes of stats (change of avg, change of median)
     CLengthCropDeltastats(1,1:ncrop)=CycleLengthCrop(1,1:ncrop);

    for c=1:ncrop
        for stat=1:2
        CLengthCropDeltastats{2,c}(stat,1:nsc)=(CLengthCropstats{2,c}(stat,1:nsc)-CLengthCropstats{2,c}(stat,1));
        end
    end

clear c

% 4.5 Vizualize growing cycle length impact
%------------------------------------------------------------------------- 
NameMat= {'RCP 4.5','RCP 8.5'};

%search indices of crops you want to show
maize=find(strcmp(CycleLengthCrop(1,:),'Maize')==1);
wwheat=find(strcmp(CycleLengthCrop(1,:),'WinterWheat')==1);
sugarbeet=find(strcmp(CycleLengthCrop(1,:),'Sugarbeet')==1);
potato=find(strcmp(CycleLengthCrop(1,:),'Potato')==1);
pea=find(strcmp(CycleLengthCrop(1,:),'Pea')==1);

figure('name','Median growing cycle length changes')
        sub(1)=subplot(2,5,1,'fontsize',10);
        boxplot(CLengthCropDeltastats{2,maize}(2,2:16),RCP(1,2:16),'grouporder',{'RCP4.5','RCP8.5'},'labels',NameMat,'colors',['r','b']);
        line(xlim,[0,0],'Color','k','LineStyle','--')
        ylabel('Change of growing cycle length (days)')
        axis([xlim, -50,10])
        set(gca,'box','off');
        title('Maize')
        
        sub(2)=subplot(2,5,2,'fontsize',10);
        boxplot(CLengthCropDeltastats{2,wwheat}(2,2:16),RCP(1,2:16),'grouporder',{'RCP4.5','RCP8.5'},'labels',NameMat);
        line(xlim,[0,0],'Color','k','LineStyle','--')
        set(gca,'box','off');
        title('Winter Wheat')
        
        sub(3)=subplot(2,5,3,'fontsize',10);
        boxplot(CLengthCropDeltastats{2,potato}(2,2:16),RCP(1,2:16),'grouporder',{'RCP4.5','RCP8.5'},'labels',NameMat);
        line(xlim,[0,0],'Color','k','LineStyle','--')
        set(gca,'box','off');
        title('Potato')
        
        sub(4)=subplot(2,5,4,'fontsize',10);
        boxplot(CLengthCropDeltastats{2,sugarbeet}(2,2:16),RCP(1,2:16),'grouporder',{'RCP4.5','RCP8.5'},'labels',NameMat);
        line(xlim,[0,0],'Color','k','LineStyle','--')
        set(gca,'box','off');
        title('Sugarbeet')
        
        sub(5)=subplot(2,5,5,'fontsize',10);
        boxplot(CLengthCropDeltastats{2,pea}(2,2:16),RCP(1,2:16),'grouporder',{'RCP4.5','RCP8.5'},'labels',NameMat);
        line(xlim,[0,0],'Color','k','LineStyle','--')
        set(gca,'box','off');
        title('Peas')    

        linkaxes(sub,'y')% link y axis of different plots (so that they change simultaneously
        
        h(1)=subplot(2,5,6,'fontsize',10);
        bar(CycleLengthCrop{2,maize}(2,1))
        ylabel('Medain historical growing cycle length (days)')
        axis([xlim , 0 ,300])
        title('Maize')
        set(gca,'XTickLabel',{' '});
        
        h(2)=subplot(2,5,7,'fontsize',10);
        bar(CycleLengthCrop{2,wwheat}(2,1))
        title('Winter Wheat')
        set(gca,'XTickLabel',{' '});
        
        h(3)=subplot(2,5,8,'fontsize',10);
        bar(CycleLengthCrop{2,potato}(2,1))
        title('Potato')
        set(gca,'XTickLabel',{' '});
         
        h(4)=subplot(2,5,9,'fontsize',10);
        bar(CycleLengthCrop{2,sugarbeet}(2,1))
        title('Sugarbeet')
        set(gca,'XTickLabel',{' '});
                
        h(5)=subplot(2,5,10,'fontsize',10);
        bar(CycleLengthCrop{2,pea}(2,1))
        title('Peas')
        set(gca,'XTickLabel',{' '});
        
        linkaxes(h,'y')% link y axis of different plots (so that they change simultaneously
 
clear sub h  


%% -----------------------------------------------------------------------
% 6. HYDROLOGICAL IMPACT 
%------------------------------------------------------------------------

% 6.1 Subtotal calculation 
% ----------------------------------------------------------------------- 
Q_MTFcum=cumsum(Q_MTF); % daily cumulative

[Q_MTFyear,Q_MTFmonth,~,Q_MTFseason,Q_MTFymonth,~,Q_MTFyseason]=ClimSubtotal(Date(:,1),Q_MTF*f/area,'sum');


% 6.2 Yearly Subtotal analysis and vizualization
% ----------------------------------------------------------------------- 

% put data in matrix (not structure)
    nyear=year(EndDate(1,1))-year(StartDate(1,1))+1;
    Q_MTFyear2=NaN(nyear,nsc);

    for sc=1:nsc
    Q_MTFyear2(:,sc)=Q_MTFyear{1,sc}(:,2);
    end

% stats for each crop over different years 
    Q_MTFyearstats=NaN(5,nsc);
     
    Q_MTFyearstats(1,1:nsc)=mean(Q_MTFyear2(:,1:nsc));
    Q_MTFyearstats(2,1:nsc)=median(Q_MTFyear2(:,1:nsc));
    Q_MTFyearstats(3,1:nsc)=std(Q_MTFyear2(:,1:nsc));
    Q_MTFyearstats(4,1:nsc)=min(Q_MTFyear2(:,1:nsc));
    Q_MTFyearstats(5,1:nsc)=max(Q_MTFyear2(:,1:nsc));
       
%  Mean changes of statistics (mean and median)
    Q_MTFyearDeltastats=NaN(2,nsc);
    
    for stat=1:2
        Q_MTFyearDeltastats(stat,1:nsc)=(Q_MTFyearstats(stat,1:nsc)-Q_MTFyearstats(stat,1))./Q_MTFyearstats(stat,1);
    end

% vizualization    
    NameMat= {'RCP 4.5','RCP 8.5'};
    figure('name','Median yearly discharge changes')
            sub(1)=subplot(2,1,2,'fontsize',10);
            boxplot(Q_MTFyearDeltastats(1,2:16)*100,RCP(1,2:16),'grouporder',{'RCP4.5','RCP8.5'},'labels',NameMat,'colors',['r','b']);
            line(xlim,[0,0],'Color','k','LineStyle','--')
            ylabel('Median annual flow change (%)')
            axis([xlim, -10,30])
            set(gca,'box','off')
            
            sub(2)=subplot(2,1,1,'fontsize',10);
            P=plot(Date(:,1),Q_MTFcum(1:nTime,2:nsc)*f/area,'color',[0.6 0.6 0.6]);
            hold on
            plot(Date(:,1),Q_MTFcum(1:nTime,1)*f/area,'color','k','LineWidth',2);
            ylabel('Total river flow (mm)')

% 6.3 Month Subtotal analysis and vizualization
% ----------------------------------------------------------------------- 

% reorganize data in one matrix   
    Q_MTFymonth2=cell(1,12);

    for m=1:12
        mindex=find(Q_MTFymonth{1,1}(:,2)==m);
        for sc =1:nsc
        Q_MTFymonth2{1,m}(:,sc)=Q_MTFymonth{1,sc}(mindex,3);
        end
    end

% stats for each month and each scenario
    Q_MTFmonthstats=cell(2,12);
    Q_MTFmonthstats(1,1:12)={1:12};
    
    for m=1:12 
        Q_MTFmonthstats{2,m}(1,1:nsc)=mean(Q_MTFymonth2{1,m}(:,1:nsc));
        Q_MTFmonthstats{2,m}(2,1:nsc)=median(Q_MTFymonth2{1,m}(:,1:nsc));
        Q_MTFmonthstats{2,m}(3,1:nsc)=std(Q_MTFymonth2{1,m}(:,1:nsc));
        Q_MTFmonthstats{2,m}(4,1:nsc)=min(Q_MTFymonth2{1,m}(:,1:nsc));
        Q_MTFmonthstats{2,m}(5,1:nsc)=max(Q_MTFymonth2{1,m}(:,1:nsc));
    end
    
%  Mean changes of statistics (mean and median)for each month and each
%  scenario
    Q_MTFmonthDeltastats=cell(2,12);
    Q_MTFmonthDeltastats(1,1:12)=Q_MTFmonthstats(1,1:12);
   
    for m=1:12 
        for stat=1:2
            Q_MTFmonthDeltastats{2,m}(stat,1:nsc)=(Q_MTFmonthstats{2,m}(stat,1:nsc)-Q_MTFmonthstats{2,m}(stat,1))./Q_MTFmonthstats{2,m}(stat,1);
        end
    end
    
%% vizualization    
   
    figure('name','Median monthly discharge changes')
            sub(1)=subplot('Position',[0.05, 0.4, 0.055,0.55],'fontsize',10); 
            boxplot(Q_MTFmonthDeltastats{2,1}(1,2:16)*100,RCP(1,2:16),'grouporder',{'RCP4.5','RCP8.5'},'colors',['r','b']);
            line(xlim,[0,0],'Color','k','LineStyle','--')
            ylabel('Median monthly flow change (%)')
            xlabel('January')
            axis([xlim, -60,60])
            set(gca,'XTick',[])
            set(gca,'box','off')
            
            sub(2)=subplot('Position',[0.16, 0.4, 0.04,0.55],'fontsize',10); 
            boxplot(Q_MTFmonthDeltastats{2,2}(1,2:16)*100,RCP(1,2:16),'grouporder',{'RCP4.5','RCP8.5'},'colors',['r','b']);
            line(xlim,[0,0],'Color','k','LineStyle','--')
            xlabel('February')           
            set(gca,'YTick',[],'XTick',[])
            set(gca,'box','off')
            
            sub(3)=subplot('Position',[0.2208, 0.4, 0.04,0.55],'fontsize',10); 
            boxplot(Q_MTFmonthDeltastats{2,3}(1,2:16)*100,RCP(1,2:16),'grouporder',{'RCP4.5','RCP8.5'},'colors',['r','b']);
            line(xlim,[0,0],'Color','k','LineStyle','--')
            xlabel('March')
            set(gca,'YTick',[],'XTick',[])
            set(gca,'box','off')
            
            sub(4)=subplot('Position',[0.28, 0.4, 0.04,0.55],'fontsize',10); 
            boxplot(Q_MTFmonthDeltastats{2,4}(1,2:16)*100,RCP(1,2:16),'grouporder',{'RCP4.5','RCP8.5'},'colors',['r','b']);
            line(xlim,[0,0],'Color','k','LineStyle','--')
            xlabel('April')
            set(gca,'YTick',[],'XTick',[])
            set(gca,'box','off')

            sub(5)=subplot('Position',[0.3355, 0.4, 0.04,0.55],'fontsize',10); 
            boxplot(Q_MTFmonthDeltastats{2,5}(1,2:16)*100,RCP(1,2:16),'grouporder',{'RCP4.5','RCP8.5'},'colors',['r','b']);
            line(xlim,[0,0],'Color','k','LineStyle','--')
            xlabel('May')
             set(gca,'YTick',[],'XTick',[])
            set(gca,'box','off')
            
            sub(6)=subplot('Position',[0.401, 0.4, 0.04,0.55],'fontsize',10); 
            boxplot(Q_MTFmonthDeltastats{2,6}(1,2:16)*100,RCP(1,2:16),'grouporder',{'RCP4.5','RCP8.5'},'colors',['r','b']);
            line(xlim,[0,0],'Color','k','LineStyle','--')
            xlabel('June')
            set(gca,'YTick',[],'XTick',[])
            set(gca,'box','off')
            
            title('Changes of median monthly discharge')            
            sub(7)=subplot('Position',[0.46, 0.4, 0.04,0.55],'fontsize',10); 
            boxplot(Q_MTFmonthDeltastats{2,7}(1,2:16)*100,RCP(1,2:16),'grouporder',{'RCP4.5','RCP8.5'},'colors',['r','b']);
            line(xlim,[0,0],'Color','k','LineStyle','--')
            xlabel('July')
            set(gca,'YTick',[],'XTick',[])
            set(gca,'box','off')
          
            sub(8)=subplot('Position',[0.525, 0.4, 0.04,0.55],'fontsize',10); 
            boxplot(Q_MTFmonthDeltastats{2,8}(1,2:16)*100,RCP(1,2:16),'grouporder',{'RCP4.5','RCP8.5'},'colors',['r','b']);
            line(xlim,[0,0],'Color','k','LineStyle','--')
            xlabel('August')
            set(gca,'YTick',[],'XTick',[])
            set(gca,'box','off')
          
            sub(9)=subplot('Position',[0.587, 0.4, 0.04,0.55],'fontsize',10); 
            boxplot(Q_MTFmonthDeltastats{2,9}(1,2:16)*100,RCP(1,2:16),'grouporder',{'RCP4.5','RCP8.5'},'colors',['r','b']);
            line(xlim,[0,0],'Color','k','LineStyle','--')
            xlabel('September')
            set(gca,'YTick',[],'XTick',[])
            set(gca,'box','off')
              
            sub(10)=subplot('Position',[0.64, 0.4, 0.04,0.55],'fontsize',10); 
            boxplot(Q_MTFmonthDeltastats{2,10}(1,2:16)*100,RCP(1,2:16),'grouporder',{'RCP4.5','RCP8.5'},'colors',['r','b']);
            line(xlim,[0,0],'Color','k','LineStyle','--')
            xlabel('October')  
            set(gca,'YTick',[],'XTick',[])
            set(gca,'box','off')
          
            sub(11)=subplot('Position',[0.725, 0.4, 0.04,0.55],'fontsize',10); 
            boxplot(Q_MTFmonthDeltastats{2,11}(1,2:16)*100,RCP(1,2:16),'grouporder',{'RCP4.5','RCP8.5'},'colors',['r','b']);
            line(xlim,[0,0],'Color','k','LineStyle','--')
            xlabel('November')
            set(gca,'YTick',[],'XTick',[])
            set(gca,'box','off')
          
            sub(12)=subplot('Position',[0.785, 0.4, 0.04,0.55],'fontsize',10); 
            boxplot(Q_MTFmonthDeltastats{2,12}(1,2:16)*100,RCP(1,2:16),'grouporder',{'RCP4.5','RCP8.5'},'colors',['r','b']);
            line(xlim,[0,0],'Color','k','LineStyle','--')
            xlabel('December')
            set(gca,'YTick',[],'XTick',[])
            set(gca,'box','off')     
            linkaxes(sub,'y')% link y axis of different plots (so that they change simultaneously
            
            subplot('Position',[0.06, 0.1, 0.76,0.2],'fontsize',12);
                for i=1:12
                HistValue(i,1)=Q_MTFmonthstats{2,i}(1,2);
                end
            bar(1:12,HistValue)
            ylabel('Median monthly flow (mm/month)','fontsize',10)
            axis([0.5,12.5, ylim])
            set(gca,'XTick',[1:12])
            set(gca,'box','off')
            
            
            