%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DESCRIPTION: This script compares different scenarios that are simulated with
% AquaCrop-Hydro. It runs each scenario and compares crop yield, crop water productivity, length of the growing cycle, temperature and water stress
% and river discharge (cumulative volumes) for every scenario
%
% TO DO before running the script: 
%   1. Run AquaCrop simulations for all land units in the catchment for
%      each scenario (with and without heat stress)
%   2. Ensure that all AquaCrop output files are numbered with format "01", "02",
%      "03", "10","20" according to the landunit number
%   3. Ensure that all AquaCrop ouput files are organized in subfolders
%      with the name of the scenario
%   4. Prepare all required input files for this script
%   5. Make sure that the first scenario (specified in Scenario.txt) is
%   the baseline scenario 
%   
%
% Author: Hanne Van Gaelen
% Last updated: 31/03/2016
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ------------------------------------------------------------------------
% 1. DEFINE DATA PATHS (BY USER)                                                    
%--------------------------------------------------------------------------      

% Add the paths where the DateCalc.m en ClimSubtotal.m file are
%   located to the searchpath of matlab
    newpath1='C:\DATA_HanneV\~Onderzoek\DEEL 2 - Ecohydro model\DEEL IIB  -VHM-AC\Github\AquaCrop_Extra';
    newpath2='C:\DATA_HanneV\~Onderzoek\DEEL 2 - Ecohydro model\DEEL IIB  -VHM-AC\Github\Climate_Processing';
    path(path,newpath1) %adds the newpath folder to the bottom of the search path. 
    path(path,newpath2)

    clear newpath1 newpath2
    
% specify in which path all AquaCrop simulation output files are stored for
% all scenarios (each scenario in different subfolder). Make sure also to
% store the AquaCrop climate files (Temp, ET0, rainfall)

     DatapathAC=uigetdir('C:\','Select directory of scenario subfloders with AquaCrop output files for all landunits');     

% specify in which path all AquaCrop simulation output files are stored for
% all scenarios (each scenario in different subfolder) FOR THE AQUACROP SIMULATION WITHOUT HEAT STRESS. Make sure also to
% store the AquaCrop climate files (Temp, ET0, rainfall)

     DatapathACHeat=uigetdir('C:\','Select directory of scenario subfloders with AquaCrop output files for all landunits NO HEAT STRESS');     
     
% specify in which path all inputs for AquaCrop-Hydro are stored including 
%    a) additional information on AquaCrop simulations of each landunits (SimInfo.txt)
%    b) the parameters for the hydrological model (Parameters.txt)
%    c) the maximum root depth for every landunit and sim run (Zrx.txt)
%    d) the soil parameters of each soil type present (SoilPar.txt)

     DatapathACHinput = uigetdir('C:\','Select directory with all input files for AquaCrop-Hydro ');
     
% specify in which path all output should be stored including 
%   a) Flow values (baseflow, interflow, overland flow, total flow) as
%   simulated by AquaCrop-Hydro
%   b) Catchment water balance values as simulated by AquaCrop-Hydro

     DatapathACHOutput=uigetdir('C:\','Select directory with scenario subfloders to store AquaCrop-Hydro output (flows and water balance values)');        

% specify in which path all information on the scenario comparison is stored including  
 %    a) information on the different scenarios (Scenario.txt)    
     
     DatapathScenIn = uigetdir('C:\','Select directory with all input information on the scenario comparison ');     

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
     
% specify the category that you want to group scenarios in
     GroupCat=inputdlg('Do you want to group scenario analysis per RCP scenario (1), climate model (2) or management category (3)?','Grouping category');
     GroupCat=cell2mat(GroupCat);
     GroupCat=str2double(GroupCat);

     if GroupCat==1 || GroupCat==2 || GroupCat==3
         %continue
     else
         error('invalid grouping category selected - it can only be 1-2-3');
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
    nsc=length(scnumb);
    ScenarioName=A.textdata(:,1).';
    ClimModel=A.textdata(:,2).';
    RCP=A.textdata(:,3).';
    Manag=A.textdata(:,4).';
    StartDate=datetime((A.data(:,2).'),'ConvertFrom','excel');
    EndDate=datetime((A.data(:,3).'),'ConvertFrom','excel');
    nTime=daysact(datenum(StartDate(1,1)), datenum(EndDate(1,1)))+1;
    clear A 
    
% 2.3 Create groups and labels for vizualization 
%-------------------------------------------------------------------------     

    if GroupCat==1 % grouping per RCP
        groupnames=unique(RCP(1,2:nsc),'stable'); %baseline excluded
        groupnames2=unique(RCP,'stable'); % baseline not excluded
        ngroup=length(groupnames);
        ngroup2=length(groupnames2);   
        groupmat=RCP;
    elseif GroupCat==2 % grouping per climate model
        groupnames=unique(ClimModel(1,2:nsc),'stable');
        groupnames2=unique(ClimModel,'stable');
        ngroup=length(groupnames) ;
        ngroup2=length(groupnames2);   
        groupmat=ClimModel;
    elseif GroupCat==3 % grouping per management scenario
        groupnames=unique(Manag(1,2:nsc),'stable');
        groupnames2=unique(Manag,'stable');
        ngroup=length(groupnames);  
        ngroup2=length(groupnames2);   
        groupmat=Manag;
    else
        error('grouping category is not well defined');
    end
    
    % format of potential groups
    linesstructall={'-','-','--',':','-.','-'}; 
    linewstructall={1.5,0.5,0.5,0.5,0.5,0.5};
    
    if  GroupCat==2
        colorstructall={'[0 0 0]','[0.6 0.6 0.6]','[1 0 0]','[1 0 0]','[1 0 0]','[1 0 0]'};
    else
        colorstructall={'[0 0 0]','[0.6 0.6 0.6]','[0.6 0.6 0.6]','[0.3 0.3 0.3]','[0.3 0.3 0.3]','[0.3 0.3 0.3]'};
    end

    % format of actual groups
    linesstruct=cell(nsc,1); 
    colorstruct=cell(nsc,1);
    linewstruct=cell(nsc,1);
    
    linesstructg=cell(ngroup2,1); 
    colorstructg=cell(ngroup2,1);
    linewstructg=cell(ngroup2,1);

    for g=1:ngroup2 
        index=strcmp(groupmat,groupnames2(1,g));
        linesstruct(index==1,1)=linesstructall(1,g);
        colorstruct(index==1,1)=colorstructall(1,g);
        linewstruct(index==1,1)=linewstructall(1,g);
    end
    
    for g=1:ngroup2 
        linesstructg(g,1)=linesstructall(1,g);
        colorstructg(g,1)=colorstructall(1,g);
        linewstructg(g,1)=linewstructall(1,g);
    end
    
    
    boxplotcolor=[0 0 0;0.6 0.6 0.6;0.4 0.4 0.4];
    
%% -----------------------------------------------------------------------
% 3. RUN AQUACROP-HYDRO FOR ALL SCENARIOS & SAVE OUTPUT OF EACH SCENARIO
%-------------------------------------------------------------------------            

% 3.1 initialize variables
%-------------------------------------------------------------------------
    %time variables
    Day=NaN(nTime,nsc);    % Day number
    Month=NaN(nTime,nsc);  % Month number
    Year=NaN(nTime,nsc);   % Year number
    %Date=NaT(nTime,nsc); % function works not in Matlab 2015a
    
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
    Wr2Catch=NaN(nTime,nsc);    % Soil water content in 2 m soil depth (mm)
    CCCatch=NaN(nTime,nsc);     % Canopy Cover (%)

   % AquaCrop results per simulation unit 
    Tr=cell(2,nsc) ;        % Crop transpiration(actual)(mm)
    Trx=cell(2,nsc) ;       % Potential (maximum) crop transpiration(mm)
    E=cell(2,nsc) ;         % Evaporation (actual)(mm)
    Ex=cell(2,nsc) ;        % Potential (maximum) evaporation(mm)
    ETa=cell(2,nsc) ;       % Actual evapotranspiration(mm)
    ETx=cell(2,nsc) ;       % Potential (maximum) evapotranspiration(mm)
    RO=cell(2,nsc) ;        % Runoff(mm)
    DP=cell(2,nsc) ;        % Deep percolation(mm)
    CR=cell(2,nsc) ;        % Capilary rise(mm)
    BundWat=cell(2,nsc) ;   % Water between bunds (mm)
    Wr2=cell(2,nsc) ;       % Soil water content in 2 m soil depth (mm)
    CC=cell(2,nsc) ;        % Canopy Cover (%)
    B=cell(2,nsc) ;         % Dry aboveground biomass during growing season (ton/ha)
    GDD=cell(2,nsc);        % Growing degree days accumulated on that day (°C-day)
   
   % Crop results for main season (average over catchment)
     Prod=cell(2,nsc);  % crop production variabiles for main season crops
                        % for each scenario and each crop there is a matrix with 
                            % Bfinact=Actual simulated dry aboveground biomass at maturity (ton/ha)
                            % Bfinpot= Potential simulated biomass at maturity if no stresses(ton/ha)
                            % Bfinrel= Bfinact/Bfinpot (%)
                            % Yact= Actual simulated final yield at maturity (ton/ha)  
                            % HIact= Actual simulated Harvest index at maturity (%) as affected by stresses   
                            % LGPact=  Actual simulated length of growing period (days) as affected by early senescence
    Prodh=cell(2,nsc);
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
    Prod(1,1:nsc)= ScenarioName(1,1:nsc);
    Prodh(1,1:nsc)= ScenarioName(1,1:nsc);
    GDD(1,1:nsc)=ScenarioName(1,1:nsc);

% 3.2 Run ACHydro and save output
%-------------------------------------------------------------------------
    
for sc=9:nsc %loop trough all scenarios

% show progress
disp(['Now running scenario ',num2str(sc),' of  ',num2str(nsc)]);
 
% Extract scenario name
Name=ScenarioName{1,sc};

% Datapaths for this scenario
  DatapathACSC=fullfile(DatapathAC,Name);
  DatapathACSCh=fullfile(DatapathACHeat,Name);
  DatapathInputSC=DatapathACHinput;                % input is the same for all scenarios   
  DatapathOutputSC=fullfile(DatapathACHOutput,Name);  

% Check if AquaCrop results for this scenario can be found  
    
    if exist(DatapathACSC) == 7 %#ok<EXIST>
        %continue as the AquaCrop results for a scneario with this name
        %can be found
    else
        error(['The AquaCrop results for scenario ',num2str(sc),' with scenarioname ',Name,' could not be found'])
    end   
    
    if exist(DatapathACSCh) == 7 %#ok<EXIST>
        %continue as the AquaCrop results for a scneario with this name
        %can be found
    else
        error(['The AquaCrop results for scenario ',num2str(sc),' with scenarioname ',Name,' (NO HEAT STRESS) could not be found'])
    end  
    
    
% Run AquaCrop-Hydro for this scenario
[Q_MBFsc,Q_MIFsc,Q_MOFsc,Q_MTFsc,area,f,Wrmin,Wrmax,pbf,SoilPar,SimACOutput,CatchACOutput,CropCatchACOutput,Par]=AquaCropHydro(DatapathACSC, DatapathInputSC,ACMode); % with heat stress
[~,~,~,~,~,~,~,~,~,~,~,~,CropCatchACOutputh,~]=AquaCropHydro(DatapathACSCh, DatapathInputSC,ACMode); % without heat stress for pollination

% extract output for this scenario       
        
        % Save time variables        
        Day(:,sc)=CatchACOutput(:,13);    % Day number
        Month(:,sc)=CatchACOutput(:,14);  % Month number
        Year(:,sc)=CatchACOutput(:,15);   % Year number
        Date(:,sc)=datetime(Year(:,sc),Month(:,sc),Day(:,sc)); %#ok<SAGROW> % Date 
                
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
        GDD{2,sc}=SimACOutput{1,19};    % Dry aboveground biomass during growing season (ton/ha)
        
      % Save crop production results for main season   
        Prod{2,sc}=CropCatchACOutput(1:2,:); 
        Prodh{2,sc}=CropCatchACOutputh(1:2,:); 
        
      % Save the catchment hydrology results   
        Q_MBF(1:nTime,sc)=Q_MBFsc(1:nTime,1);
        Q_MIF(1:nTime,sc)=Q_MIFsc(1:nTime,1);
        Q_MOF(1:nTime,sc)=Q_MOFsc(1:nTime,1);      
        Q_MTF(1:nTime,sc)=Q_MTFsc(1:nTime,1);

 
% write output for this scenario (with heat stress) to excel      
      % Combine output in matrix if necessary
        HeadersFlow={'Date','Baseflow','Interflow','Overland flow', 'Total flow'};
        FlowOutput=[exceltime(Date(1:nTime,1)),Q_MBFsc(1:nTime,1),Q_MIFsc(1:nTime,1),Q_MOFsc(1:nTime,1),Q_MTFsc(1:nTime,1)];

        HeadersWabalCatch={'Date','Tr','Trx','E','Ex','ETa','ETx','RO','DP','CR','BundWat','Wr2'};

       % Write output to one excel tabsheet
        xlname='FlowSimResults.xlsx';
        filename = fullfile(DatapathOutputSC,xlname);
        xlswrite(filename,HeadersFlow,'SimFlow','A1');
        xlswrite(filename,FlowOutput,'SimFlow','A2');
    
        xlname='WabalSimResults.xlsx';
        filename = fullfile(DatapathOutputSC,xlname);
        xlswrite(filename,HeadersWabalCatch,'SimWabal','A1');
        xlswrite(filename,exceltime(Date(1:nTime,1)),'SimWabal','A2');
        xlswrite(filename,CatchACOutput(1:nTime,1:11),'SimWabal','B2'); 
        
 
% save workspace variables (to be sure)
filename=['workspace ',datestr(date),' untilscenario',num2str(sc)];
filename=fullfile(DatapathScenOut,filename);
save(filename)
clear filename 
 
        
clear SimACOutput CatchOutput FlowOutput WabalCatchOutput DatapathOutputSC DatapathInputSC   
clear Q_MTFsc Q_MOFsc Q_MIFsc Q_MBFsc
end

clear sc

% 3.3 Reorganize production variables per crop
%-------------------------------------------------------------------------
% Put yield, WPET, DSI, WSI,BWSI TSI en LGP in one structure 

Cropnames= Prod{2,1}(1,:);
[~,ncrop]=size(Cropnames);

Yact(1,1:ncrop)=Cropnames(1,1:ncrop);
DSI(1,1:ncrop)=Cropnames(1,1:ncrop);  % drought stress index (based on yield loss)
WSI(1,1:ncrop)=Cropnames(1,1:ncrop);  % water stress index (based on transpiration)
BWSI(1,1:ncrop)=Cropnames(1,1:ncrop); % water stress index (based on biomass loss)
WP(1,1:ncrop)=Cropnames(1,1:ncrop); % water productivity
TSI(1,1:ncrop)=Cropnames(1,1:ncrop); % temperature (heat) stress index
HSI(1,1:ncrop)=Cropnames(1,1:ncrop);  % cold stres index
LGPact(1,1:ncrop)=Cropnames(1,1:ncrop); % actual length of growing period

% initialize subsets
nyear=NaN(ncrop,1);
for c=1:ncrop
[nyear(c,1),~]=size(Prod{2,1}{2,c}(:,1));
end
subsetY=NaN(max(nyear(:)),nsc);
subsetYh=NaN(max(nyear(:)),nsc);% yield for conditions where heat stress is not considered
subsetDSI=NaN(max(nyear(:)),nsc);
subsetWSI=NaN(max(nyear(:)),nsc);
subsetBWSI=NaN(max(nyear(:)),nsc);
subsetWP=NaN(max(nyear(:)),nsc);
subsetTS=NaN(max(nyear(:)),nsc);
subsetLGP=NaN(max(nyear(:)),nsc);

for c=1:ncrop % loop trough each crop
    for sc=1:nsc %loop trough each scenario
        subsetY(1:nyear(c,1),sc)=Prod{2,sc}{2,c}(:,4); 
        subsetYh(1:nyear(c,1),sc)=Prodh{2,sc}{2,c}(:,4); 
        subsetDSI(1:nyear(c,1),sc)=Prod{2,sc}{2,c}(:,8);
        subsetWSI(1:nyear(c,1),sc)=Prod{2,sc}{2,c}(:,11);
        subsetBWSI(1:nyear(c,1),sc)=100-Prod{2,sc}{2,c}(:,3); % 100- Brel
        subsetWP(1:nyear(c,1),sc)=Prod{2,sc}{2,c}(:,9);
        subsetTS(1:nyear(c,1),sc)=Prod{2,sc}{2,c}(:,10);
        subsetLGP(1:nyear(c,1),sc)=Prod{2,sc}{2,c}(:,7);    
    end
    
    Yact{2,c}=subsetY;
    DSI{2,c}=subsetDSI;
    WSI{2,c}=subsetWSI;
    BWSI{2,c}=subsetBWSI;
    WP{2,c}=subsetWP;
    TSI{2,c}=subsetTS;
    HSI{2,c}=(subsetYh-subsetY)./subsetYh;
    LGPact{2,c}=subsetLGP;  
    
    clear subsetY;
    clear subsetYh;
    clear subsetDSI;
    clear subsetWSI;
    clear subsetBWSI;
    clear subsetWP;
    clear subsetTS;
    clear subsetLGP;
end

clear c sc subsetY subsetDSI subsetWSI subsetBWSI subsetWP subsetTS subsetLGP subset Yh nyear

% 3.4 Define index of crops you want to show
%------------------------------------------------------------------------- 
maize=find(strcmp(Cropnames(1,1:ncrop),'Maize')==1);
wwheat=find(strcmp(Cropnames(1,1:ncrop),'WinterWheat')==1);
sugarbeet=find(strcmp(Cropnames(1,1:ncrop),'Sugarbeet')==1);
potato=find(strcmp(Cropnames(1,1:ncrop),'Potato')==1);
pea=find(strcmp(Cropnames(1,1:ncrop),'Pea')==1);

% 3.5 Save workspace
%-------------------------------------------------------------------------
% save workspace variables so that you can skip this the first part of the
% code next time

filename=['workspace ',datestr(date),' all scenarios'];
filename=fullfile(DatapathScenOut,filename);
save(filename)
clear filename 
 

%% -----------------------------------------------------------------------
% 4. YIELD IMPACT 
%------------------------------------------------------------------------

% 4.1 Calculate statistics
%-------------------------------------------------------------------------
% stats for each crop over different years 
    Yactstats(1,1:ncrop)=Yact(1,1:ncrop); 
    
    for c=1:ncrop
        Yactstats{2,c}(1,1:nsc)=mean(Yact{2,c}(:,1:nsc));
        Yactstats{2,c}(2,1:nsc)=median(Yact{2,c}(:,1:nsc));
        Yactstats{2,c}(3,1:nsc)=std(Yact{2,c}(:,1:nsc));
        Yactstats{2,c}(4,1:nsc)=min(Yact{2,c}(:,1:nsc));
        Yactstats{2,c}(5,1:nsc)=max(Yact{2,c}(:,1:nsc));
        Yactstats{2,c}(6,1:nsc)=Yactstats{2,c}(3,1:nsc)./Yactstats{2,c}(1,1:nsc);
        Yactstats{2,c}(7,1:nsc)=Yactstats{2,c}(5,1:nsc)-Yactstats{2,c}(4,1:nsc);
    end
clear c

% Calculate changes of stats (change of avg, change of median)
    YactDeltastats(1,1:ncrop)=Yact(1,1:ncrop); 
    YactDeltastats2(1,1:ncrop)=Yact(1,1:ncrop); 
    
    for c=1:ncrop
        for stat=1:2
        YactDeltastats{2,c}(stat,1:nsc)=(Yactstats{2,c}(stat,1:nsc)-Yactstats{2,c}(stat,1))./Yactstats{2,c}(stat,1);
        YactDeltastats2{2,c}(stat,1:nsc)=(Yactstats{2,c}(stat,1:nsc)-Yactstats{2,c}(stat,1));
        end
    end

clear c
   
% 4.2 Vizualize yield impact with boxplots
%-------------------------------------------------------------------------
           
f1=figure('name','Median yield changes');%(boxplot= variation over different GCMs) 
        sub(1)=subplot(2,5,1,'fontsize',10);
        boxplot(YactDeltastats{2,maize}(2,2:nsc)*100,groupmat(1,2:nsc),'grouporder',groupnames,'labels',groupnames);
        line(xlim,[0,0],'Color','k','LineStyle','--')
        ylabel('Median yield change (%)')
        axis([xlim, -40,50])
        set(gca,'box','off')
        title('Maize')
        
        sub(2)=subplot(2,5,2,'fontsize',10);
        boxplot(YactDeltastats{2,wwheat}(2,2:nsc)*100,groupmat(1,2:nsc),'grouporder',groupnames,'labels',groupnames);
        line(xlim,[0,0],'Color','k','LineStyle','--')
        set(gca,'box','off')
        title('Winter Wheat')
        
        sub(3)=subplot(2,5,3,'fontsize',10);
        boxplot(YactDeltastats{2,potato}(2,2:nsc)*100,groupmat(1,2:nsc),'grouporder',groupnames,'labels',groupnames);
        line(xlim,[0,0],'Color','k','LineStyle','--')
        set(gca,'box','off')
        title('Potato')
        
        sub(4)=subplot(2,5,4,'fontsize',10);
        boxplot(YactDeltastats{2,sugarbeet}(2,2:nsc)*100,groupmat(1,2:nsc),'grouporder',groupnames,'labels',groupnames);
        line(xlim,[0,0],'Color','k','LineStyle','--')
        set(gca,'box','off')
        title('Sugarbeet')
        
        sub(5)=subplot(2,5,5,'fontsize',10);
        boxplot(YactDeltastats{2,pea}(2,2:nsc)*100,groupmat(1,2:nsc),'grouporder',groupnames,'labels',groupnames);
        line(xlim,[0,0],'Color','k','LineStyle','--')
        set(gca,'box','off')
        title('Peas')    

        linkaxes(sub,'y')% link y axis of different plots (so that they change simultaneously
        
        h(1)=subplot(2,5,6,'fontsize',10);
        bar(Yactstats{2,maize}(2,1))
        ylabel('Median historical yield (ton/ha)')
        axis([xlim , 0 ,15])
        title('Maize')
        set(gca,'XTickLabel',{' '})
        
        h(2)=subplot(2,5,7,'fontsize',10);
        bar(Yactstats{2,wwheat}(2,1))
        title('Winter Wheat')
        set(gca,'XTickLabel',{' '})
        
        h(3)=subplot(2,5,8,'fontsize',10);
        bar(Yactstats{2,potato}(2,1))
        title('Potato')
        set(gca,'XTickLabel',{' '})
         
        h(4)=subplot(2,5,9,'fontsize',10);
        bar(Yactstats{2,sugarbeet}(2,1))
        title('Sugarbeet')
        set(gca,'XTickLabel',{' '})
                
        h(5)=subplot(2,5,10,'fontsize',10);
        bar(Yactstats{2,pea}(2,1))
        title('Peas')
        set(gca,'XTickLabel',{' '})
        
        linkaxes(h,'y')% link y axis of different plots (so that they change simultaneously

        clear sub h 
        
% 4.3 Analyse with focus on interannual variability as well
% -----------------------------------------------------------------------        

% Change of yield as compared to historical median for all years
    YallDelta(1,1:ncrop)=Yact(1,1:ncrop); 
    
    YallDelta{2,maize}=(Yact{2,maize}-Yactstats{2,maize}(2,1))./Yactstats{2,maize}(2,1);   
    YallDelta{2,wwheat}=(Yact{2,wwheat}-Yactstats{2,wwheat}(2,1))./Yactstats{2,wwheat}(2,1);   
    YallDelta{2,sugarbeet}=(Yact{2,sugarbeet}-Yactstats{2,sugarbeet}(2,1))./Yactstats{2,sugarbeet}(2,1);   
    YallDelta{2,potato}=(Yact{2,potato}-Yactstats{2,potato}(2,1))./Yactstats{2,potato}(2,1);   
    YallDelta{2,pea}=(Yact{2,pea}-Yactstats{2,pea}(2,1))./Yactstats{2,pea}(2,1);  
    
f2=figure('name','Median yield changes- GCM&year variation');% boxplot = variation over different GCMs & over 30 different year)   
    sub(1)=subplot(1,5,1,'fontsize',10);
    boxplot(YallDelta{2,maize}(:,1:nsc)*100,groupmat(1,1:nsc),'grouporder',groupnames2,'labels',groupnames2)
    line(xlim,[0,0],'Color','k','LineStyle','--')
    ylabel('Annual yield from historical median (%)')
    title('maize')
    axis([xlim, -100,100])
    set(gca,'box','off')

    sub(2)=subplot(1,5,2,'fontsize',10);
    boxplot(YallDelta{2,wwheat}(:,1:nsc)*100,groupmat(1,1:nsc),'grouporder',groupnames2,'labels',groupnames2)
    line(xlim,[0,0],'Color','k','LineStyle','--')
    title('winter wheat')
    set(gca,'box','off','YTick',[])

    sub(3)=subplot(1,5,3,'fontsize',10);
    boxplot(YallDelta{2,sugarbeet}(:,1:nsc)*100,groupmat(1,1:nsc),'grouporder',groupnames2,'labels',groupnames2)
    line(xlim,[0,0],'Color','k','LineStyle','--')
    title('sugarbeet')
    set(gca,'box','off','YTick',[])

    sub(4)=subplot(1,5,4,'fontsize',10);
    boxplot(YallDelta{2,potato}(:,1:nsc)*100,groupmat(1,1:nsc),'grouporder',groupnames2,'labels',groupnames2)
    line(xlim,[0,0],'Color','k','LineStyle','--')
    title('potato')
    set(gca,'box','off','YTick',[])

    sub(5)=subplot(1,5,5,'fontsize',10);
    boxplot(YallDelta{2,pea}(:,1:nsc)*100,groupmat(1,1:nsc),'grouporder',groupnames2,'labels',groupnames2)
    line(xlim,[0,0],'Color','k','LineStyle','--')
    title('pea')
    set(gca,'box','off','YTick',[])
    
    linkaxes(sub,'y')   
    clear sub
      
% 4.4 Vizualize yield impact with cumulative distribution function
%-------------------------------------------------------------------------

% normality check 
[notnormalmaize,~]=NormalityCheck(Yact{2,maize},'lillie',0.05);
[notnormalwwheat,~]=NormalityCheck(Yact{2,wwheat},'lillie',0.05);
[notnormalsbeet,~]=NormalityCheck(Yact{2,sugarbeet},'lillie',0.05);
[notnormalpotato,~]=NormalityCheck(Yact{2,potato},'lillie',0.05);
[notnormalpea,~]=NormalityCheck(Yact{2,pea},'lillie',0.05);

if isempty(notnormalmaize)==1 && isempty(notnormalwwheat)==1 && isempty(notnormalsbeet)==1 && isempty(notnormalpotato)==1 && isempty(notnormalpea)==1
    disp('Yield values for all crops and all scenarios are normally distributed')
else
    if isempty(notnormalmaize)==0
    warning(['Maize yield is not normally distributed for scenarios: ',num2str(notnormalmaize.')])
    end
    
    if isempty(notnormalwwheat)==0
    warning(['Winter wheat yield is not normally distributed for scenarios: ',num2str(notnormalwwheat.')])
    end
    
    if isempty(notnormalsbeet)==0
    warning(['Sugar beet yield is not normally distributed for scenarios: ',num2str(notnormalsbeet.')])
    end
    
    if isempty(notnormalpotato)==0
    warning(['Potato yield is not normally distributed for scenarios: ',num2str(notnormalpotato.')])
    end
    
    if isempty(notnormalpea)==0
    warning(['Pea yield is not normally distributed for scenarios: ',num2str(notnormalpea.')])
    end
end

clear notnormalwwheat notnormalmaize notnormalpotato notnormalsbeet notnormalpea

% fit theoretical normal distributions
xrangemaize=0:0.5:max(Yact{2,maize}(:))+1.5;
xrangewwheat=0:0.5:max(Yact{2,wwheat}(:))+1.5;
xrangesbeet=0:0.5:max(Yact{2,sugarbeet}(:))+1.5;
xrangepotato=0:0.5:max(Yact{2,potato}(:))+1.5;
xrangepea=0:0.1:max(Yact{2,pea}(:))+1;

probabilitiesmaize=NaN(length(xrangemaize),nsc);
probabilitieswwheat=NaN(length(xrangewwheat),nsc);
probabilitiessbeet=NaN(length(xrangesbeet),nsc);
probabilitiespotato=NaN(length(xrangepotato),nsc);
probabilitiespea=NaN(length(xrangepea),nsc);

for sc=1:nsc
pdsc=fitdist(Yact{2,maize}(:,sc),'Normal');
probabilitiesmaize(:,sc)=cdf(pdsc,xrangemaize);

pdsc=fitdist(Yact{2,wwheat}(:,sc),'Normal');
probabilitieswwheat(:,sc)=cdf(pdsc,xrangewwheat);

pdsc=fitdist(Yact{2,sugarbeet}(:,sc),'Normal');
probabilitiessbeet(:,sc)=cdf(pdsc,xrangesbeet);

pdsc=fitdist(Yact{2,potato}(:,sc),'Normal');
probabilitiespotato(:,sc)=cdf(pdsc,xrangepotato);

pdsc=fitdist(Yact{2,pea}(:,sc),'Normal');
probabilitiespea(:,sc)=cdf(pdsc,xrangepea);
end

clear pdsc sc

% vizualize
f3=figure('name','Seasonal yield theoretical CDF');
    subplot(3,2,1,'fontsize',10);
    P=plot(xrangemaize,probabilitiesmaize(:,1:nsc));   
    set(P,{'Color'},colorstruct,{'LineStyle'},linesstruct,{'LineWidth'},linewstruct)
    ylabel('Cumulative probability','fontsize',8);
    xlabel('Annual yield (ton/ha)','fontsize',8);
    title('Maize')
    set(gca,'box','off')

    subplot(3,2,2,'fontsize',10);
    P=plot(xrangewwheat,probabilitieswwheat(:,1:nsc));   
    set(P,{'Color'},colorstruct,{'LineStyle'},linesstruct,{'LineWidth'},linewstruct)
    ylabel('Cumulative probability','fontsize',8);
    xlabel('Annual yield (ton/ha)','fontsize',8);
    title('Winter Wheat')
    set(gca,'box','off') 
    
    subplot(3,2,3,'fontsize',10);
    P=plot(xrangesbeet,probabilitiessbeet(:,1:nsc));   
    set(P,{'Color'},colorstruct,{'LineStyle'},linesstruct,{'LineWidth'},linewstruct)
    ylabel('Cumulative probability','fontsize',8);
    xlabel('Annual yield (ton/ha)','fontsize',8);
    title('Sugarbeet')
    set(gca,'box','off')
    
    subplot(3,2,4,'fontsize',10);
    P=plot(xrangepotato,probabilitiespotato(:,1:nsc));   
    set(P,{'Color'},colorstruct,{'LineStyle'},linesstruct,{'LineWidth'},linewstruct)
    ylabel('Cumulative probability','fontsize',8);
    xlabel('Annual yield (ton/ha)','fontsize',8);
    title('Potato')
    set(gca,'box','off')
    
    subplot(3,2,5,'fontsize',10);
    P=plot(xrangepea,probabilitiespea(:,1:nsc));   
    set(P,{'Color'},colorstruct,{'LineStyle'},linesstruct,{'LineWidth'},linewstruct)
    ylabel('Cumulative probability','fontsize',8);
    xlabel('Annual yield (ton/ha)','fontsize',8);
    title('pea')
    set(gca,'box','off')
    clear P
    
f4=figure('name','Seasonal yield emperical CDF');
    subplot(3,2,1,'fontsize',10);
    P=NaN(nsc,1);
    for i=1:nsc
        P(i)=cdfplot(Yact{2,maize}(:,i));
        hold on 
    end
    set(P,{'Color'},colorstruct,{'LineStyle'},linesstruct,{'LineWidth'},linewstruct)
    ylabel('Cumulative probability','fontsize',8);
    xlabel('Seasonal yield (ton/ha)','fontsize',8);
    title('Maize')
    set(gca,'box','off')
    grid off

    subplot(3,2,2,'fontsize',10);
    for i=1:nsc
        P(i)=cdfplot(Yact{2,wwheat}(:,i));
        hold on 
    end
    set(P,{'Color'},colorstruct,{'LineStyle'},linesstruct,{'LineWidth'},linewstruct)
    ylabel('Cumulative probability','fontsize',8);
    xlabel('Seasonal yield (ton/ha)','fontsize',8);
    title('Winter Wheat')
    set(gca,'box','off') 
    grid off
    
    subplot(3,2,3,'fontsize',10);
    for i=1:nsc
        P(i)=cdfplot(Yact{2,sugarbeet}(:,i));
        hold on 
    end  
    set(P,{'Color'},colorstruct,{'LineStyle'},linesstruct,{'LineWidth'},linewstruct)
    ylabel('Cumulative probability','fontsize',8);
    xlabel('Seasonal yield (ton/ha)','fontsize',8);
    title('Sugar beet')
    set(gca,'box','off')
    grid off
    
    subplot(3,2,4,'fontsize',10);
    for i=1:nsc
        P(i)=cdfplot(Yact{2,potato}(:,i));
        hold on 
    end    
    set(P,{'Color'},colorstruct,{'LineStyle'},linesstruct,{'LineWidth'},linewstruct)
    ylabel('Cumulative probability','fontsize',8);
    xlabel('Seasonal yield (ton/ha)','fontsize',8);
    title('Potato')
    set(gca,'box','off')
    grid off

    subplot(3,2,5,'fontsize',10);
    for i=1:nsc
        P(i)=cdfplot(Yact{2,pea}(:,i));
        hold on 
    end    
    set(P,{'Color'},colorstruct,{'LineStyle'},linesstruct,{'LineWidth'},linewstruct)
    ylabel('Cumulative probability','fontsize',8);
    xlabel('Seasonal yield (ton/ha)','fontsize',8);
    title('Pea')
    set(gca,'box','off')
    grid off
    clear P
    
f5=figure('name','Seasonal yield emperical CDF -Simple');
    subplot(3,2,1,'fontsize',10);
    P=NaN(ngroup2,1);
    for g=1:ngroup2
        gindex=find(strcmp(groupmat(1,1:nsc),groupnames2{1,g})==1);
        P(g)=cdfplot(reshape(Yact{2,maize}(:,gindex),[],1));
        hold on 
    end
    set(P,{'Color'},colorstructg,{'LineStyle'},linesstructg,{'LineWidth'},linewstructg)
    ylabel('Cumulative probability','fontsize',8);
    xlabel('Seasonal yield (ton/ha)','fontsize',8);
    title('Maize')
    set(gca,'box','off')
    grid off

    subplot(3,2,2,'fontsize',10);
    for g=1:ngroup2
        gindex=find(strcmp(groupmat(1,1:nsc),groupnames2{1,g})==1);
        P(g)=cdfplot(reshape(Yact{2,wwheat}(:,gindex),[],1));
        hold on 
    end
    set(P,{'Color'},colorstructg,{'LineStyle'},linesstructg,{'LineWidth'},linewstructg)
    ylabel('Cumulative probability','fontsize',8);
    xlabel('Seasonal yield (ton/ha)','fontsize',8);
    title('Winter Wheat')
    set(gca,'box','off') 
    grid off
    
    subplot(3,2,3,'fontsize',10);
    for g=1:ngroup2
        gindex=find(strcmp(groupmat(1,1:nsc),groupnames2{1,g})==1);
        P(g)=cdfplot(reshape(Yact{2,sugarbeet}(:,gindex),[],1));
        hold on 
    end  
    set(P,{'Color'},colorstructg,{'LineStyle'},linesstructg,{'LineWidth'},linewstructg)
    ylabel('Cumulative probability','fontsize',8);
    xlabel('Seasonal yield (ton/ha)','fontsize',8);
    title('Sugar beet')
    set(gca,'box','off')
    grid off
    
    subplot(3,2,4,'fontsize',10);
    for g=1:ngroup2
        gindex=find(strcmp(groupmat(1,1:nsc),groupnames2{1,g})==1);
        P(g)=cdfplot(reshape(Yact{2,potato}(:,gindex),[],1));
        hold on 
    end    
    set(P,{'Color'},colorstructg,{'LineStyle'},linesstructg,{'LineWidth'},linewstructg)
    ylabel('Cumulative probability','fontsize',8);
    xlabel('Seasonal yield (ton/ha)','fontsize',8);
    title('Potato')
    set(gca,'box','off')
    grid off

    subplot(3,2,5,'fontsize',10);
    for g=1:ngroup2
        gindex=find(strcmp(groupmat(1,1:nsc),groupnames2{1,g})==1);
        P(g)=cdfplot(reshape(Yact{2,pea}(:,gindex),[],1));
        hold on 
    end    
    set(P,{'Color'},colorstructg,{'LineStyle'},linesstructg,{'LineWidth'},linewstructg)
    ylabel('Cumulative probability','fontsize',8);
    xlabel('Seasonal yield (ton/ha)','fontsize',8);
    title('Pea')
    set(gca,'box','off')
    grid off
       
    
clear xrangepotato xrangemaize xrangewwheat xrangesbeet xrangepea 
clear i P gindex

% 4.5 Save vizualization
%-------------------------------------------------------------------------
filename='Yield median changes - GCMboxplots';
filename=fullfile(DatapathScenOut,filename);
savefig(f1,filename)

filename='Yield empirical CDF ';
filename=fullfile(DatapathScenOut,filename);
savefig(f4,filename)

filename='Yield empirical CDF - simple';
filename=fullfile(DatapathScenOut,filename);
savefig(f5,filename)

clear filename f1 f2 f3 f4 

%% -----------------------------------------------------------------------
% 5. WPET IMPACT 
%------------------------------------------------------------------------

% 5.1 Calculate WPET statistics
%-------------------------------------------------------------------------
% stats for each crop over different years 

    WPstats(1,1:ncrop)=WP(1,1:ncrop); 
    
    for c=1:ncrop
        WPstats{2,c}(1,1:nsc)=mean(WP{2,c}(:,1:nsc));
        WPstats{2,c}(2,1:nsc)=median(WP{2,c}(:,1:nsc));
        WPstats{2,c}(3,1:nsc)=std(WP{2,c}(:,1:nsc));
        WPstats{2,c}(4,1:nsc)=min(WP{2,c}(:,1:nsc));
        WPstats{2,c}(5,1:nsc)=max(WP{2,c}(:,1:nsc));
        WPstats{2,c}(6,1:nsc)=WPstats{2,c}(3,1:nsc)./WPstats{2,c}(1,1:nsc);
        WPstats{2,c}(7,1:nsc)=WPstats{2,c}(5,1:nsc)-WPstats{2,c}(4,1:nsc);
    end
clear c

% Calculate changes of stats (change of avg, change of median)
    WPDeltastats(1,1:ncrop)=WP(1,1:ncrop); 
    WPDeltastats2(1,1:ncrop)=WP(1,1:ncrop); 
    
    for c=1:ncrop
        for stat=1:2
        WPDeltastats{2,c}(stat,1:nsc)=(WPstats{2,c}(stat,1:nsc)-WPstats{2,c}(stat,1))./WPstats{2,c}(stat,1);
        WPDeltastats2{2,c}(stat,1:nsc)=(WPstats{2,c}(stat,1:nsc)-WPstats{2,c}(stat,1));
        end
    end

clear c stat

% 5.2 Vizualize WPET impact with boxplots
%-------------------------------------------------------------------------

f1=figure('name','Median WPET changes');%(boxplot= variation over different GCMs) 
        sub(1)=subplot(2,5,1,'fontsize',10);
        boxplot(WPDeltastats{2,maize}(2,2:nsc)*100,groupmat(1,2:nsc),'grouporder',groupnames,'labels',groupnames);
        line(xlim,[0,0],'Color','k','LineStyle','--')
        ylabel('Median WPET change (%)')
        axis([xlim, -30,150])
        set(gca,'box','off')
        title('Maize')
        
        sub(2)=subplot(2,5,2,'fontsize',10);
        boxplot(WPDeltastats{2,wwheat}(2,2:nsc)*100,groupmat(1,2:nsc),'grouporder',groupnames,'labels',groupnames);
        line(xlim,[0,0],'Color','k','LineStyle','--')
        set(gca,'box','off')
        title('Winter Wheat')
        
        sub(3)=subplot(2,5,3,'fontsize',10);
        boxplot(WPDeltastats{2,potato}(2,2:nsc)*100,groupmat(1,2:nsc),'grouporder',groupnames,'labels',groupnames);
        line(xlim,[0,0],'Color','k','LineStyle','--')
        set(gca,'box','off')
        title('Potato')
        
        sub(4)=subplot(2,5,4,'fontsize',10);
        boxplot(WPDeltastats{2,sugarbeet}(2,2:nsc)*100,groupmat(1,2:nsc),'grouporder',groupnames,'labels',groupnames);
        line(xlim,[0,0],'Color','k','LineStyle','--')
        set(gca,'box','off')
        title('Sugarbeet')
        
        sub(5)=subplot(2,5,5,'fontsize',10);
        boxplot(WPDeltastats{2,pea}(2,2:nsc)*100,groupmat(1,2:nsc),'grouporder',groupnames,'labels',groupnames);
        line(xlim,[0,0],'Color','k','LineStyle','--')
        set(gca,'box','off')
        title('Peas')    

        linkaxes(sub,'y')% link y axis of different plots (so that they change simultaneously
        
        h(1)=subplot(2,5,6,'fontsize',10);
        bar(WPstats{2,maize}(2,1))
        ylabel('Median historical WP(%)')
        axis([xlim , 0 ,4])
        title('Maize')
        set(gca,'XTickLabel',{' '})
        
        h(2)=subplot(2,5,7,'fontsize',10);
        bar(WPstats{2,wwheat}(2,1))
        title('Winter Wheat')
        set(gca,'XTickLabel',{' '})
        
        h(3)=subplot(2,5,8,'fontsize',10);
        bar(WPstats{2,potato}(2,1))
        title('Potato')
        set(gca,'XTickLabel',{' '})
         
        h(4)=subplot(2,5,9,'fontsize',10);
        bar(WPstats{2,sugarbeet}(2,1))
        title('Sugarbeet')
        set(gca,'XTickLabel',{' '})
                
        h(5)=subplot(2,5,10,'fontsize',10);
        bar(WPstats{2,pea}(2,1))
        title('Peas')
        set(gca,'XTickLabel',{' '})
        
        linkaxes(h,'y')% link y axis of different plots (so that they change simultaneously

        clear sub h 
        
% 5.3 Vizualize WPET with cumulative distribution function
%-------------------------------------------------------------------------
f2=figure('name','WP empirical CDF');
    subplot(3,2,1,'fontsize',10);
    P=NaN(nsc,1);
    for i=1:nsc
        P(i)=cdfplot(WP{2,maize}(:,i));
        hold on 
    end
    set(P,{'Color'},colorstruct,{'LineStyle'},linesstruct,{'LineWidth'},linewstruct)
    ylabel('Cumulative probability','fontsize',8);
    xlabel('WPET(kg/m³)','fontsize',8);
    title('Maize')
    set(gca,'box','off')
    grid off

    subplot(3,2,2,'fontsize',10);
    for i=1:nsc
        P(i)=cdfplot(WP{2,wwheat}(:,i));
        hold on 
    end
    set(P,{'Color'},colorstruct,{'LineStyle'},linesstruct,{'LineWidth'},linewstruct)
    ylabel('Cumulative probability','fontsize',8);
    xlabel('WPET (kg/m³)','fontsize',8);
    title('Winter Wheat')
    set(gca,'box','off') 
    grid off
    
    subplot(3,2,3,'fontsize',10);
    for i=1:nsc
        P(i)=cdfplot(WP{2,sugarbeet}(:,i));
        hold on 
    end  
    set(P,{'Color'},colorstruct,{'LineStyle'},linesstruct,{'LineWidth'},linewstruct)
    ylabel('Cumulative probability','fontsize',8);
    xlabel('WPET (kg/m³)','fontsize',8);
    title('Sugar beet')
    set(gca,'box','off')
    grid off
    
    subplot(3,2,4,'fontsize',10);
    for i=1:nsc
        P(i)=cdfplot(WP{2,potato}(:,i));
        hold on 
    end    
    set(P,{'Color'},colorstruct,{'LineStyle'},linesstruct,{'LineWidth'},linewstruct)
    ylabel('Cumulative probability','fontsize',8);
    xlabel('WPET (kg/m³)','fontsize',8);
    title('Potato')
    set(gca,'box','off')
    grid off

    subplot(3,2,5,'fontsize',10);
    for i=1:nsc
        P(i)=cdfplot(WP{2,pea}(:,i));
        hold on 
    end    
    set(P,{'Color'},colorstruct,{'LineStyle'},linesstruct,{'LineWidth'},linewstruct)
    ylabel('Cumulative probability','fontsize',8);
    xlabel('WPET (kg/m³)','fontsize',8);
    title('Pea')
    set(gca,'box','off')
    grid off
    
    clear P i 
    
    
% 5.4 save vizualization
%-------------------------------------------------------------------------
filename='WPET median changes - GCMboxplots';
filename=fullfile(DatapathScenOut,filename);
savefig(f1,filename)

filename='WPET empirical CDF ';
filename=fullfile(DatapathScenOut,filename);
savefig(f2,filename)

clear filename 
clear f1 f2 
%% -----------------------------------------------------------------------
% 6. DROUGHT STRESS INDEX (DSI), WATER STRESS INDEX (WSI) and TEMPERATURE STRESS(TSI) IMPACT 
%------------------------------------------------------------------------

% 6.1 Calculate DSI, WSI en TSI statistics
%-------------------------------------------------------------------------
% stats for each crop over different years 
    DSIstats(1,1:ncrop)=DSI(1,1:ncrop); 
    TSIstats(1,1:ncrop)=TSI(1,1:ncrop);
    HSIstats(1,1:ncrop)=HSI(1,1:ncrop);
    WSIstats(1,1:ncrop)=WSI(1,1:ncrop);
    BWSIstats(1,1:ncrop)=BWSI(1,1:ncrop);
    
    for c=1:ncrop
        DSIstats{2,c}(1,1:nsc)=mean(DSI{2,c}(:,1:nsc));
        DSIstats{2,c}(2,1:nsc)=median(DSI{2,c}(:,1:nsc));
        DSIstats{2,c}(3,1:nsc)=std(DSI{2,c}(:,1:nsc));
        DSIstats{2,c}(4,1:nsc)=min(DSI{2,c}(:,1:nsc));
        DSIstats{2,c}(5,1:nsc)=max(DSI{2,c}(:,1:nsc));
        DSIstats{2,c}(6,1:nsc)=prctile(DSI{2,c}(:,1:nsc),95);
        
        TSIstats{2,c}(1,1:nsc)=mean(TSI{2,c}(:,1:nsc));
        TSIstats{2,c}(2,1:nsc)=median(TSI{2,c}(:,1:nsc));
        TSIstats{2,c}(3,1:nsc)=std(TSI{2,c}(:,1:nsc));
        TSIstats{2,c}(4,1:nsc)=min(TSI{2,c}(:,1:nsc));
        TSIstats{2,c}(5,1:nsc)=max(TSI{2,c}(:,1:nsc));
        TSIstats{2,c}(6,1:nsc)=prctile(TSI{2,c}(:,1:nsc),95);
  
        HSIstats{2,c}(1,1:nsc)=mean(HSI{2,c}(:,1:nsc));
        HSIstats{2,c}(2,1:nsc)=median(HSI{2,c}(:,1:nsc));
        HSIstats{2,c}(3,1:nsc)=std(HSI{2,c}(:,1:nsc));
        HSIstats{2,c}(4,1:nsc)=min(HSI{2,c}(:,1:nsc));
        HSIstats{2,c}(5,1:nsc)=max(HSI{2,c}(:,1:nsc));
        HSIstats{2,c}(6,1:nsc)=prctile(HSI{2,c}(:,1:nsc),95);
        
        WSIstats{2,c}(1,1:nsc)=mean(WSI{2,c}(:,1:nsc));
        WSIstats{2,c}(2,1:nsc)=median(WSI{2,c}(:,1:nsc));
        WSIstats{2,c}(3,1:nsc)=std(WSI{2,c}(:,1:nsc));
        WSIstats{2,c}(4,1:nsc)=min(WSI{2,c}(:,1:nsc));
        WSIstats{2,c}(5,1:nsc)=max(WSI{2,c}(:,1:nsc));
        WSIstats{2,c}(6,1:nsc)=prctile(WSI{2,c}(:,1:nsc),95);
        
        BWSIstats{2,c}(1,1:nsc)=mean(BWSI{2,c}(:,1:nsc));
        BWSIstats{2,c}(2,1:nsc)=median(BWSI{2,c}(:,1:nsc));
        BWSIstats{2,c}(3,1:nsc)=std(BWSI{2,c}(:,1:nsc));
        BWSIstats{2,c}(4,1:nsc)=min(BWSI{2,c}(:,1:nsc));
        BWSIstats{2,c}(5,1:nsc)=max(BWSI{2,c}(:,1:nsc));
        BWSIstats{2,c}(6,1:nsc)=prctile(BWSI{2,c}(:,1:nsc),95);
    end
clear c

% Calculate changes of stats (change of avg, change of median)
    DSIDeltastats(1,1:ncrop)=DSI(1,1:ncrop); 
    TSIDeltastats(1,1:ncrop)=TSI(1,1:ncrop); 
    HSIDeltastats(1,1:ncrop)=HSI(1,1:ncrop); 
    WSIDeltastats(1,1:ncrop)=WSI(1,1:ncrop); 
    BWSIDeltastats(1,1:ncrop)=BWSI(1,1:ncrop); 
    
    for c=1:ncrop
        for stat=1:2
        DSIDeltastats{2,c}(stat,1:nsc)=(DSIstats{2,c}(stat,1:nsc)-DSIstats{2,c}(stat,1));
        HSIDeltastats{2,c}(stat,1:nsc)=(HSIstats{2,c}(stat,1:nsc)-HSIstats{2,c}(stat,1));       
        TSIDeltastats{2,c}(stat,1:nsc)=(TSIstats{2,c}(stat,1:nsc)-TSIstats{2,c}(stat,1));
        WSIDeltastats{2,c}(stat,1:nsc)=(WSIstats{2,c}(stat,1:nsc)-WSIstats{2,c}(stat,1));      
        BWSIDeltastats{2,c}(stat,1:nsc)=(BWSIstats{2,c}(stat,1:nsc)-BWSIstats{2,c}(stat,1));  
        end
    end

clear c stat


% 6.2 Vizualize DSI, WSI, TSI & HSI impact with boxplots
%-------------------------------------------------------------------------

f1=figure('name','Median DSI changes');%(boxplot= variation over different GCMs) 
        sub(1)=subplot(2,5,1,'fontsize',10);
        boxplot(DSIDeltastats{2,maize}(2,2:nsc),groupmat(1,2:nsc),'grouporder',groupnames,'labels',groupnames);
        line(xlim,[0,0],'Color','k','LineStyle','--')
        ylabel('Median DSI change (%)')
        axis([xlim, -5,20])
        set(gca,'box','off')
        title('Maize')
        
        sub(2)=subplot(2,5,2,'fontsize',10);
        boxplot(DSIDeltastats{2,wwheat}(2,2:nsc),groupmat(1,2:nsc),'grouporder',groupnames,'labels',groupnames);
        line(xlim,[0,0],'Color','k','LineStyle','--')
        set(gca,'box','off')
        title('Winter Wheat')
        
        sub(3)=subplot(2,5,3,'fontsize',10);
        boxplot(DSIDeltastats{2,potato}(2,2:nsc),groupmat(1,2:nsc),'grouporder',groupnames,'labels',groupnames);
        line(xlim,[0,0],'Color','k','LineStyle','--')
        set(gca,'box','off')
        title('Potato')
        
        sub(4)=subplot(2,5,4,'fontsize',10);
        boxplot(DSIDeltastats{2,sugarbeet}(2,2:nsc),groupmat(1,2:nsc),'grouporder',groupnames,'labels',groupnames);
        line(xlim,[0,0],'Color','k','LineStyle','--')
        set(gca,'box','off')
        title('Sugarbeet')
        
        sub(5)=subplot(2,5,5,'fontsize',10);
        boxplot(DSIDeltastats{2,pea}(2,2:nsc),groupmat(1,2:nsc),'grouporder',groupnames,'labels',groupnames);
        line(xlim,[0,0],'Color','k','LineStyle','--')
        set(gca,'box','off')
        title('Peas')    

        linkaxes(sub,'y')% link y axis of different plots (so that they change simultaneously
        
        h(1)=subplot(2,5,6,'fontsize',10);
        bar(DSIstats{2,maize}(2,1))
        ylabel('Median historical DSI(%)')
        axis([xlim , 0 ,20])
        title('Maize')
        set(gca,'XTickLabel',{' '})
        
        h(2)=subplot(2,5,7,'fontsize',10);
        bar(DSIstats{2,wwheat}(2,1))
        title('Winter Wheat')
        set(gca,'XTickLabel',{' '})
        
        h(3)=subplot(2,5,8,'fontsize',10);
        bar(DSIstats{2,potato}(2,1))
        title('Potato')
        set(gca,'XTickLabel',{' '})
         
        h(4)=subplot(2,5,9,'fontsize',10);
        bar(DSIstats{2,sugarbeet}(2,1))
        title('Sugarbeet')
        set(gca,'XTickLabel',{' '})
                
        h(5)=subplot(2,5,10,'fontsize',10);
        bar(DSIstats{2,pea}(2,1))
        title('Peas')
        set(gca,'XTickLabel',{' '})
        
        linkaxes(h,'y')% link y axis of different plots (so that they change simultaneously

        clear sub h       

f1b=figure('name','Median WSI changes');%(boxplot= variation over different GCMs) 
        sub(1)=subplot(2,5,1,'fontsize',10);
        boxplot(WSIDeltastats{2,maize}(2,2:nsc),groupmat(1,2:nsc),'grouporder',groupnames,'labels',groupnames);
        line(xlim,[0,0],'Color','k','LineStyle','--')
        ylabel('Median WSI change (%)')
        axis([xlim, -20,10])
        set(gca,'box','off')
        title('Maize')
        
        sub(2)=subplot(2,5,2,'fontsize',10);
        boxplot(WSIDeltastats{2,wwheat}(2,2:nsc),groupmat(1,2:nsc),'grouporder',groupnames,'labels',groupnames);
        line(xlim,[0,0],'Color','k','LineStyle','--')
        set(gca,'box','off')
        title('Winter Wheat')
        
        sub(3)=subplot(2,5,3,'fontsize',10);
        boxplot(WSIDeltastats{2,potato}(2,2:nsc),groupmat(1,2:nsc),'grouporder',groupnames,'labels',groupnames);
        line(xlim,[0,0],'Color','k','LineStyle','--')
        set(gca,'box','off')
        title('Potato')
        
        sub(4)=subplot(2,5,4,'fontsize',10);
        boxplot(WSIDeltastats{2,sugarbeet}(2,2:nsc),groupmat(1,2:nsc),'grouporder',groupnames,'labels',groupnames);
        line(xlim,[0,0],'Color','k','LineStyle','--')
        set(gca,'box','off')
        title('Sugarbeet')
        
        sub(5)=subplot(2,5,5,'fontsize',10);
        boxplot(WSIDeltastats{2,pea}(2,2:nsc),groupmat(1,2:nsc),'grouporder',groupnames,'labels',groupnames);
        line(xlim,[0,0],'Color','k','LineStyle','--')
        set(gca,'box','off')
        title('Peas')    

        linkaxes(sub,'y')% link y axis of different plots (so that they change simultaneously
        
        h(1)=subplot(2,5,6,'fontsize',10);
        bar(WSIstats{2,maize}(2,1))
        ylabel('Median historical WSI(%)')
        axis([xlim , 0 ,100])
        title('Maize')
        set(gca,'XTickLabel',{' '})
        
        h(2)=subplot(2,5,7,'fontsize',10);
        bar(WSIstats{2,wwheat}(2,1))
        title('Winter Wheat')
        set(gca,'XTickLabel',{' '})
        
        h(3)=subplot(2,5,8,'fontsize',10);
        bar(WSIstats{2,potato}(2,1))
        title('Potato')
        set(gca,'XTickLabel',{' '})
         
        h(4)=subplot(2,5,9,'fontsize',10);
        bar(WSIstats{2,sugarbeet}(2,1))
        title('Sugarbeet')
        set(gca,'XTickLabel',{' '})
                
        h(5)=subplot(2,5,10,'fontsize',10);
        bar(WSIstats{2,pea}(2,1))
        title('Peas')
        set(gca,'XTickLabel',{' '})
        
        linkaxes(h,'y')% link y axis of different plots (so that they change simultaneously

        clear sub h         
f1c=figure('name','Median BWSI changes');%(boxplot= variation over different GCMs) 
        sub(1)=subplot(2,5,1,'fontsize',10);
        boxplot(BWSIDeltastats{2,maize}(2,2:nsc),groupmat(1,2:nsc),'grouporder',groupnames,'labels',groupnames);
        line(xlim,[0,0],'Color','k','LineStyle','--')
        ylabel('Median BWSI change (%)')
        axis([xlim, -20,20])
        set(gca,'box','off')
        title('Maize')
        
        sub(2)=subplot(2,5,2,'fontsize',10);
        boxplot(BWSIDeltastats{2,wwheat}(2,2:nsc),groupmat(1,2:nsc),'grouporder',groupnames,'labels',groupnames);
        line(xlim,[0,0],'Color','k','LineStyle','--')
        set(gca,'box','off')
        title('Winter Wheat')
        
        sub(3)=subplot(2,5,3,'fontsize',10);
        boxplot(BWSIDeltastats{2,potato}(2,2:nsc),groupmat(1,2:nsc),'grouporder',groupnames,'labels',groupnames);
        line(xlim,[0,0],'Color','k','LineStyle','--')
        set(gca,'box','off')
        title('Potato')
        
        sub(4)=subplot(2,5,4,'fontsize',10);
        boxplot(BWSIDeltastats{2,sugarbeet}(2,2:nsc),groupmat(1,2:nsc),'grouporder',groupnames,'labels',groupnames);
        line(xlim,[0,0],'Color','k','LineStyle','--')
        set(gca,'box','off')
        title('Sugarbeet')
        
        sub(5)=subplot(2,5,5,'fontsize',10);
        boxplot(BWSIDeltastats{2,pea}(2,2:nsc),groupmat(1,2:nsc),'grouporder',groupnames,'labels',groupnames);
        line(xlim,[0,0],'Color','k','LineStyle','--')
        set(gca,'box','off')
        title('Peas')    

        linkaxes(sub,'y')% link y axis of different plots (so that they change simultaneously
        
        h(1)=subplot(2,5,6,'fontsize',10);
        bar(BWSIstats{2,maize}(2,1))
        ylabel('Median historical BWSI(%)')
        axis([xlim , 0 ,100])
        title('Maize')
        set(gca,'XTickLabel',{' '})
        
        h(2)=subplot(2,5,7,'fontsize',10);
        bar(BWSIstats{2,wwheat}(2,1))
        title('Winter Wheat')
        set(gca,'XTickLabel',{' '})
        
        h(3)=subplot(2,5,8,'fontsize',10);
        bar(BWSIstats{2,potato}(2,1))
        title('Potato')
        set(gca,'XTickLabel',{' '})
         
        h(4)=subplot(2,5,9,'fontsize',10);
        bar(BWSIstats{2,sugarbeet}(2,1))
        title('Sugarbeet')
        set(gca,'XTickLabel',{' '})
                
        h(5)=subplot(2,5,10,'fontsize',10);
        bar(BWSIstats{2,pea}(2,1))
        title('Peas')
        set(gca,'XTickLabel',{' '})
        
        linkaxes(h,'y')% link y axis of different plots (so that they change simultaneously

        clear sub h         
        
f2=figure('name','Median Temperature stress changes');%(boxplot= variation over different GCMs) 
        sub(1)=subplot(2,5,1,'fontsize',10);
        boxplot(TSIDeltastats{2,maize}(2,2:nsc),groupmat(1,2:nsc),'grouporder',groupnames,'labels',groupnames);
        line(xlim,[0,0],'Color','k','LineStyle','--')
        ylabel('Median TSI change (%)')
        axis([xlim, -20,5])
        set(gca,'box','off')
        title('Maize')
        
        sub(2)=subplot(2,5,2,'fontsize',10);
        boxplot(TSIDeltastats{2,wwheat}(2,2:nsc),groupmat(1,2:nsc),'grouporder',groupnames,'labels',groupnames);
        line(xlim,[0,0],'Color','k','LineStyle','--')
        set(gca,'box','off')
        title('Winter Wheat')
        
        sub(3)=subplot(2,5,3,'fontsize',10);
        boxplot(TSIDeltastats{2,potato}(2,2:nsc),groupmat(1,2:nsc),'grouporder',groupnames,'labels',groupnames);
        line(xlim,[0,0],'Color','k','LineStyle','--')
        set(gca,'box','off')
        title('Potato')
        
        sub(4)=subplot(2,5,4,'fontsize',10);
        boxplot(TSIDeltastats{2,sugarbeet}(2,2:nsc),groupmat(1,2:nsc),'grouporder',groupnames,'labels',groupnames);
        line(xlim,[0,0],'Color','k','LineStyle','--')
        set(gca,'box','off')
        title('Sugarbeet')
        
        sub(5)=subplot(2,5,5,'fontsize',10);
        boxplot(TSIDeltastats{2,pea}(2,2:nsc),groupmat(1,2:nsc),'grouporder',groupnames,'labels',groupnames);
        line(xlim,[0,0],'Color','k','LineStyle','--')
        set(gca,'box','off')
        title('Peas')    

        linkaxes(sub,'y')% link y axis of different plots (so that they change simultaneously
        
        h(1)=subplot(2,5,6,'fontsize',10);
        bar(TSIstats{2,maize}(2,1))
        ylabel('Median historical TSI(%)')
        axis([xlim , 0 ,40])
        title('Maize')
        set(gca,'XTickLabel',{' '})
        
        h(2)=subplot(2,5,7,'fontsize',10);
        bar(TSIstats{2,wwheat}(2,1))
        title('Winter Wheat')
        set(gca,'XTickLabel',{' '})
        
        h(3)=subplot(2,5,8,'fontsize',10);
        bar(TSIstats{2,potato}(2,1))
        title('Potato')
        set(gca,'XTickLabel',{' '})
         
        h(4)=subplot(2,5,9,'fontsize',10);
        bar(TSIstats{2,sugarbeet}(2,1))
        title('Sugarbeet')
        set(gca,'XTickLabel',{' '})
                
        h(5)=subplot(2,5,10,'fontsize',10);
        bar(TSIstats{2,pea}(2,1))
        title('Peas')
        set(gca,'XTickLabel',{' '})
        
        linkaxes(h,'y')% link y axis of different plots (so that they change simultaneously
        clear sub h         

f2b=figure('name','Median Heat stress changes');%(boxplot= variation over different GCMs) 
        sub(1)=subplot(2,5,1,'fontsize',10);
        boxplot(HSIDeltastats{2,maize}(2,2:nsc),groupmat(1,2:nsc),'grouporder',groupnames,'labels',groupnames);
        line(xlim,[0,0],'Color','k','LineStyle','--')
        ylabel('Median HSI change (%)')
        axis([xlim, -20,5])
        set(gca,'box','off')
        title('Maize')
        
        sub(2)=subplot(2,5,2,'fontsize',10);
        boxplot(HSIDeltastats{2,wwheat}(2,2:nsc),groupmat(1,2:nsc),'grouporder',groupnames,'labels',groupnames);
        line(xlim,[0,0],'Color','k','LineStyle','--')
        set(gca,'box','off')
        title('Winter Wheat')
        
        sub(3)=subplot(2,5,3,'fontsize',10);
        boxplot(HSIDeltastats{2,potato}(2,2:nsc),groupmat(1,2:nsc),'grouporder',groupnames,'labels',groupnames);
        line(xlim,[0,0],'Color','k','LineStyle','--')
        set(gca,'box','off')
        title('Potato')
        
        sub(4)=subplot(2,5,4,'fontsize',10);
        boxplot(HSIDeltastats{2,sugarbeet}(2,2:nsc),groupmat(1,2:nsc),'grouporder',groupnames,'labels',groupnames);
        line(xlim,[0,0],'Color','k','LineStyle','--')
        set(gca,'box','off')
        title('Sugarbeet')
        
        sub(5)=subplot(2,5,5,'fontsize',10);
        boxplot(HSIDeltastats{2,pea}(2,2:nsc),groupmat(1,2:nsc),'grouporder',groupnames,'labels',groupnames);
        line(xlim,[0,0],'Color','k','LineStyle','--')
        set(gca,'box','off')
        title('Peas')    

        linkaxes(sub,'y')% link y axis of different plots (so that they change simultaneously
        
        h(1)=subplot(2,5,6,'fontsize',10);
        bar(HSIstats{2,maize}(2,1))
        ylabel('Median historical HSI(%)')
        axis([xlim , 0 ,40])
        title('Maize')
        set(gca,'XTickLabel',{' '})
        
        h(2)=subplot(2,5,7,'fontsize',10);
        bar(HSIstats{2,wwheat}(2,1))
        title('Winter Wheat')
        set(gca,'XTickLabel',{' '})
        
        h(3)=subplot(2,5,8,'fontsize',10);
        bar(HSIstats{2,potato}(2,1))
        title('Potato')
        set(gca,'XTickLabel',{' '})
         
        h(4)=subplot(2,5,9,'fontsize',10);
        bar(HSIstats{2,sugarbeet}(2,1))
        title('Sugarbeet')
        set(gca,'XTickLabel',{' '})
                
        h(5)=subplot(2,5,10,'fontsize',10);
        bar(HSIstats{2,pea}(2,1))
        title('Peas')
        set(gca,'XTickLabel',{' '})
        
        linkaxes(h,'y')% link y axis of different plots (so that they change simultaneously
        clear sub h        
      
% 6.3 Vizualize DSI & TSI with cumulative distribution function
%-------------------------------------------------------------------------
f3=figure('name','DSI empirical CDF');
    P=NaN(nsc,1);
    subplot(5,2,1,'fontsize',10);
    for i=1:nsc
        P(i)=cdfplot(DSI{2,maize}(:,i));
        hold on 
    end
    set(P,{'Color'},colorstruct,{'LineStyle'},linesstruct,{'LineWidth'},linewstruct)
    ylabel('Cumulative probability','fontsize',8);
    xlabel('Drought stress index (%)','fontsize',8);
    title('Maize')
    set(gca,'box','off')
    grid off

    subplot(5,2,3,'fontsize',10);
    for i=1:nsc
        P(i)=cdfplot(DSI{2,wwheat}(:,i));
        hold on 
    end
    set(P,{'Color'},colorstruct,{'LineStyle'},linesstruct,{'LineWidth'},linewstruct)
    ylabel('Cumulative probability','fontsize',8);
    xlabel('Drought stress index (%)','fontsize',8);
    title('Winter Wheat')
    set(gca,'box','off') 
    grid off
    
    subplot(5,2,5,'fontsize',10);
    for i=1:nsc
        P(i)=cdfplot(DSI{2,sugarbeet}(:,i));
        hold on 
    end  
    set(P,{'Color'},colorstruct,{'LineStyle'},linesstruct,{'LineWidth'},linewstruct)
    ylabel('Cumulative probability','fontsize',8);
    xlabel('Drought stress index (%)','fontsize',8);
    title('Sugar beet')
    set(gca,'box','off')
    grid off
    
    subplot(5,2,7,'fontsize',10);
    for i=1:nsc
        P(i)=cdfplot(DSI{2,potato}(:,i));
        hold on 
    end    
    set(P,{'Color'},colorstruct,{'LineStyle'},linesstruct,{'LineWidth'},linewstruct)
    ylabel('Cumulative probability','fontsize',8);
    xlabel('Drought stress index (%)','fontsize',8);
    title('Potato')
    set(gca,'box','off')
    grid off

    subplot(5,2,9,'fontsize',10);
    for i=1:nsc
        P(i)=cdfplot(DSI{2,pea}(:,i));
        hold on 
    end    
    set(P,{'Color'},colorstruct,{'LineStyle'},linesstruct,{'LineWidth'},linewstruct)
    ylabel('Cumulative probability','fontsize',8);
    xlabel('Drought stress index (%)','fontsize',8);
    title('Pea')
    set(gca,'box','off')
    grid off
    
    subplot(5,2,2,'fontsize',10);
    for i=1:nsc
        P(i)=cdfplot(TSI{2,maize}(:,i));
        hold on 
    end
    set(P,{'Color'},colorstruct,{'LineStyle'},linesstruct,{'LineWidth'},linewstruct)
    ylabel('Cumulative probability','fontsize',8);
    xlabel('Temperature stress index (%)','fontsize',8);
    title('Maize')
    set(gca,'box','off')
    grid off

    subplot(5,2,4,'fontsize',10);
    for i=1:nsc
        P(i)=cdfplot(TSI{2,wwheat}(:,i));
        hold on 
    end
    set(P,{'Color'},colorstruct,{'LineStyle'},linesstruct,{'LineWidth'},linewstruct)
    ylabel('Cumulative probability','fontsize',8);
    xlabel('Temperature stress index (%)','fontsize',8);
    title('Winter Wheat')
    set(gca,'box','off') 
    grid off
    
    subplot(5,2,6,'fontsize',10);
    for i=1:nsc
        P(i)=cdfplot(TSI{2,sugarbeet}(:,i));
        hold on 
    end  
    set(P,{'Color'},colorstruct,{'LineStyle'},linesstruct,{'LineWidth'},linewstruct)
    ylabel('Cumulative probability','fontsize',8);
    xlabel('Temperature stress index (%)','fontsize',8);
    title('Sugar beet')
    set(gca,'box','off')
    grid off
    
    subplot(5,2,8,'fontsize',10);
    for i=1:nsc
        P(i)=cdfplot(TSI{2,potato}(:,i));
        hold on 
    end    
    set(P,{'Color'},colorstruct,{'LineStyle'},linesstruct,{'LineWidth'},linewstruct)
    ylabel('Cumulative probability','fontsize',8);
    xlabel('Temperature stress index (%)','fontsize',8);
    title('Potato')
    set(gca,'box','off')
    grid off

    subplot(5,2,10,'fontsize',10);
    for i=1:nsc
        P(i)=cdfplot(TSI{2,pea}(:,i));
        hold on 
    end    
    set(P,{'Color'},colorstruct,{'LineStyle'},linesstruct,{'LineWidth'},linewstruct)
    ylabel('Cumulative probability','fontsize',8);
    xlabel('Temperature stress index (%)','fontsize',8);
    title('Pea')
    set(gca,'box','off')
    grid off    
    
    clear P i
    
% 6.4 save vizualization
%-------------------------------------------------------------------------
filename='DSI median changes - GCMboxplots';
filename=fullfile(DatapathScenOut,filename);
savefig(f1,filename)

filename='WSI median changes - GCMboxplots';
filename=fullfile(DatapathScenOut,filename);
savefig(f1b,filename)

filename='BWSI median changes - GCMboxplots';
filename=fullfile(DatapathScenOut,filename);
savefig(f1c,filename)

filename='TSI median changes - GCMboxplots';
filename=fullfile(DatapathScenOut,filename);
savefig(f2,filename)

filename='DSI_TSI empirical CDF ';
filename=fullfile(DatapathScenOut,filename);
savefig(f3,filename)

clear filename 
clear f1 f1b f2 f3 
%% -----------------------------------------------------------------------
% 7. IMPACT ON LENGTH OF THE POTENTIAL LENGTH OF GROWING CYCLE 
%-------------------------------------------------------------------------
% Potential length of growing cycle = the length based purely on
% temperature and crop cycle requirements. If a crop dies of early due to
% water stress this is not taken into account

% 7.1 Calculate potential length of growing period & sowing/maturity dates
%-------------------------------------------------------------------------
LGPpotAll(1,1:nsc)=ScenarioName(1,1:nsc);
sowingAll(1,1:nsc)=ScenarioName(1,1:nsc);
maturitypotAll(1,1:nsc)=ScenarioName(1,1:nsc);

for sc=1:nsc
    DatapathDateSC=fullfile(DatapathDate,ScenarioName{1,sc});
    DatapathOutputSC=fullfile(DatapathACHOutput,ScenarioName{1,sc});
    RotationDateSC=CalcDate(DatapathDateSC,DatapathOutputSC,StartDate(1,sc),EndDate(1,sc)); 
    
    for lu=1:nlu
        if SimType(lu,1)==2
            nrun=length(RotationDateSC{1,lu}(:,3));
            r=(nrun-rem(nrun,2))/2;     
            SowingDateSC(1:r,lu)= RotationDateSC{1,lu}(2:2:nrun,3);    %#ok<SAGROW>
            MaturitypotDateSC(1:r,lu)= RotationDateSC{1,lu}(2:2:nrun,4); %#ok<SAGROW>
        else 
            % skip this landunit (not agriculture)
        end
    end
    
    LengthSC=(datenum(MaturitypotDateSC)-datenum(SowingDateSC))+1;
    
    LGPpotAll{2,sc}=LengthSC; % per scenario one matrix with cycle lengths per landunit
    sowingAll{2,sc}=SowingDateSC;
    maturitypotAll{2,sc}=MaturitypotDateSC;
end

clear sc lu SowingDateSC MaturitypotDateSC LengthSC RotationDateSC

% 7.2 Compose LGPpotential/sowing date  matrix per crop type
%--------------------------------------------------------------------------
Cropnames= Prod{2,1}(1,:);
[~,ncrop]=size(Cropnames);

LGPpot=cell(2,ncrop);% initialize
sowing=cell(2,ncrop);
LGPpot(1,1:ncrop)=Cropnames(1:1:ncrop); % write away crop name
sowing(1,1:ncrop)=Cropnames(1:1:ncrop); 
maturitypot(1,1:ncrop)=Cropnames(1:1:ncrop);

for c=1:ncrop% loop trough each crop 
     
    %search all projects with this crop
    index=find(strcmp(Crop(:,1),Cropnames(c))==1);
   
    for sc=1:nsc %loop trough al scenarios
        subset= LGPpotAll{2,sc}(:,index(1));% all landunits with same crop will have same LGPpot want same temp 
        subset=subset(all(~isnan(subset),2),:); 
        LGPpot{2,c}(:,sc)=subset;
        clear subset
        
        subset=sowingAll{2,sc}(:,index(1));
        subset=subset(all(~isnat(subset),2),:); 
        sowing{2,c}(:,sc)= subset;
        clear subset
            
        subset=maturitypotAll{2,sc}(:,index(1));
        subset=subset(all(~isnat(subset),2),:); 
        maturitypot{2,c}(:,sc)= subset;
        clear subset
    end        
end

clear sc c subset

% 7.3 Calculate LGP statistics (over different year)
%-------------------------------------------------------------------------
 LGPpotstats(1,1:ncrop)=LGPpot(1,1:ncrop); 
    
    for c=1:ncrop
        LGPpotstats{2,c}(1,1:nsc)=nanmean(LGPpot{2,c}(:,1:nsc));
        LGPpotstats{2,c}(2,1:nsc)=nanmedian(LGPpot{2,c}(:,1:nsc));
        LGPpotstats{2,c}(3,1:nsc)=nanstd(LGPpot{2,c}(:,1:nsc));
        LGPpotstats{2,c}(4,1:nsc)=min(LGPpot{2,c}(:,1:nsc));
        LGPpotstats{2,c}(5,1:nsc)=max(LGPpot{2,c}(:,1:nsc));
    end
clear c

% Calculate changes of stats (change of avg, change of median)
     LGPpotDeltastats(1,1:ncrop)=LGPpot(1,1:ncrop);

    for c=1:ncrop
        for stat=1:2
        LGPpotDeltastats{2,c}(stat,1:nsc)=(LGPpotstats{2,c}(stat,1:nsc)-LGPpotstats{2,c}(stat,1));
        end
    end

clear c stat

% Calculate average/median maturity dates
MATpotstats(1,1:ncrop)=LGPpot(1,1:ncrop); 

for c=1:ncrop
    
MATpotstats{2,c}(1,1:nsc)=sowing{2,c}(1,1:nsc)+(LGPpotstats{2,c}(1,1:nsc)-1);
MATpotstats{2,c}(2,1:nsc)=sowing{2,c}(1,1:nsc)+(LGPpotstats{2,c}(2,1:nsc)-1);
end

clear c

% 7.4 Vizualize potential length of growing period 
%------------------------------------------------------------------------- 

f1=figure('name','Median LGPpot changes'); %(boxplot= variation over different GCMs) 
        sub(1)=subplot(2,5,1,'fontsize',10);
        boxplot(LGPpotDeltastats{2,maize}(2,2:nsc),groupmat(1,2:nsc),'grouporder',groupnames,'labels',groupnames);
        line(xlim,[0,0],'Color','k','LineStyle','--')
        ylabel('Change of growing cycle length (days)')
        axis([xlim, -50,10])
        set(gca,'box','off');
        title('Maize')
        
        sub(2)=subplot(2,5,2,'fontsize',10);
        boxplot(LGPpotDeltastats{2,wwheat}(2,2:nsc),groupmat(1,2:nsc),'grouporder',groupnames,'labels',groupnames);
        line(xlim,[0,0],'Color','k','LineStyle','--')
        set(gca,'box','off');
        title('Winter Wheat')
        
        sub(3)=subplot(2,5,3,'fontsize',10);
        boxplot(LGPpotDeltastats{2,potato}(2,2:nsc),groupmat(1,2:nsc),'grouporder',groupnames,'labels',groupnames);
        line(xlim,[0,0],'Color','k','LineStyle','--')
        set(gca,'box','off');
        title('Potato')
        
        sub(4)=subplot(2,5,4,'fontsize',10);
        boxplot(LGPpotDeltastats{2,sugarbeet}(2,2:nsc),groupmat(1,2:nsc),'grouporder',groupnames,'labels',groupnames);
        line(xlim,[0,0],'Color','k','LineStyle','--')
        set(gca,'box','off');
        title('Sugarbeet')
        
        sub(5)=subplot(2,5,5,'fontsize',10);
        boxplot(LGPpotDeltastats{2,pea}(2,2:nsc),groupmat(1,2:nsc),'grouporder',groupnames,'labels',groupnames);
        line(xlim,[0,0],'Color','k','LineStyle','--')
        set(gca,'box','off');
        title('Peas')    

        linkaxes(sub,'y')% link y axis of different plots (so that they change simultaneously
        
        h(1)=subplot(2,5,6,'fontsize',10);
        bar(LGPpot{2,maize}(2,1))
        ylabel('Medain historical growing cycle length (days)')
        axis([xlim , 0 ,300])
        title('Maize')
        set(gca,'XTickLabel',{' '});
        
        h(2)=subplot(2,5,7,'fontsize',10);
        bar(LGPpot{2,wwheat}(2,1))
        title('Winter Wheat')
        set(gca,'XTickLabel',{' '});
        
        h(3)=subplot(2,5,8,'fontsize',10);
        bar(LGPpot{2,potato}(2,1))
        title('Potato')
        set(gca,'XTickLabel',{' '});
         
        h(4)=subplot(2,5,9,'fontsize',10);
        bar(LGPpot{2,sugarbeet}(2,1))
        title('Sugarbeet')
        set(gca,'XTickLabel',{' '});
                
        h(5)=subplot(2,5,10,'fontsize',10);
        bar(LGPpot{2,pea}(2,1))
        title('Peas')
        set(gca,'XTickLabel',{' '});
        
        linkaxes(h,'y')% link y axis of different plots (so that they change simultaneously

        clear sub h 

       
f2=figure('name','LGPpot emperical CDF');
    subplot(3,2,1,'fontsize',10);
    P=NaN(nsc,1);
    for i=1:nsc
        P(i)=cdfplot(LGPpot{2,maize}(:,i));
        hold on 
    end
    set(P,{'Color'},colorstruct,{'LineStyle'},linesstruct,{'LineWidth'},linewstruct)
    ylabel('Cumulative probability','fontsize',8);
    xlabel('Potential LGP (days)','fontsize',8);
    title('Maize')
    set(gca,'box','off')
    grid off

    subplot(3,2,2,'fontsize',10);
    for i=1:nsc
        P(i)=cdfplot(LGPpot{2,wwheat}(:,i));
        hold on 
    end
    set(P,{'Color'},colorstruct,{'LineStyle'},linesstruct,{'LineWidth'},linewstruct)
    ylabel('Cumulative probability','fontsize',8);
    xlabel('Potential LGP (days)','fontsize',8);
    title('Winter Wheat')
    set(gca,'box','off') 
    grid off
    
    subplot(3,2,3,'fontsize',10);
    for i=1:nsc
        P(i)=cdfplot(LGPpot{2,sugarbeet}(:,i));
        hold on 
    end  
    set(P,{'Color'},colorstruct,{'LineStyle'},linesstruct,{'LineWidth'},linewstruct)
    ylabel('Cumulative probability','fontsize',8);
    xlabel('Potential LGP (days)','fontsize',8);
    title('Sugar beet')
    set(gca,'box','off')
    grid off
    
    subplot(3,2,4,'fontsize',10);
    for i=1:nsc
        P(i)=cdfplot(LGPpot{2,potato}(:,i));
        hold on 
    end    
    set(P,{'Color'},colorstruct,{'LineStyle'},linesstruct,{'LineWidth'},linewstruct)
    ylabel('Cumulative probability','fontsize',8);
    xlabel('Potential LGP (days)','fontsize',8);
    title('Potato')
    set(gca,'box','off')
    grid off

    subplot(3,2,5,'fontsize',10);
    for i=1:nsc
        P(i)=cdfplot(LGPpot{2,pea}(:,i));
        hold on 
    end    
    set(P,{'Color'},colorstruct,{'LineStyle'},linesstruct,{'LineWidth'},linewstruct)
    ylabel('Cumulative probability','fontsize',8);
    xlabel('Potential LGP (days)','fontsize',8);
    title('Pea')
    set(gca,'box','off')
    grid off     
    
    clear P i

%save figure
filename='LGPpot median changes - GCMboxplot';
filename=fullfile(DatapathScenOut,filename);
savefig(f1,filename)
 
filename='LGPpot empirical CDF';
filename=fullfile(DatapathScenOut,filename);
savefig(f2,filename)

clear f1 f2 
%% -----------------------------------------------------------------------
% 8. IMPACT ON ACTUAL LENGTH OF GROWING CYCLE 
%-------------------------------------------------------------------------

% 8.1 Calculate LGPact and Maturity act 
%-------------------------------------------------------------------------
% Actual length of growing cycle = the length based  on
% temperature and crop cycle requirements, as well as water stress which causes early senescence. 

% LGPact= output AquaCrop

% Maturity act: calculate from sowing date and LGPact
% ATTENTION/ THIS IS NOT COOMPLETELY CORRECT. LGPact is period between
% emergence and maturity, not sowing and maturity 
maturityact(1,1:ncrop)=LGPact(1,1:ncrop); 

for c=1:ncrop
    maturityact{2,c}=(sowing{2,c}+(LGPact{2,c}-1));
    
end

% LGP gap
LGPgap(1,1:ncrop)=LGPpot(1,1:ncrop);

for c=1:ncrop
    LGPgap{2,c}=LGPpot{2,c}-LGPact{2,c};
end


clear c

% 8.2 Calculate statistics (over different year)
%-------------------------------------------------------------------------

LGPactstats(1,1:ncrop)=LGPact(1,1:ncrop); 
    
    for c=1:ncrop
        LGPactstats{2,c}(1,1:nsc)=nanmean(LGPact{2,c}(:,1:nsc));
        LGPactstats{2,c}(2,1:nsc)=nanmedian(LGPact{2,c}(:,1:nsc));
        LGPactstats{2,c}(3,1:nsc)=nanstd(LGPact{2,c}(:,1:nsc));
        LGPactstats{2,c}(4,1:nsc)=min(LGPact{2,c}(:,1:nsc));
        LGPactstats{2,c}(5,1:nsc)=max(LGPact{2,c}(:,1:nsc));
        
    end
clear c

LGPgapstats(1,1:ncrop)=LGPact(1,1:ncrop); 
    for c=1:ncrop
        LGPgapstats{2,c}(1,1:nsc)=nanmean(LGPgap{2,c}(:,1:nsc));
        LGPgapstats{2,c}(2,1:nsc)=nanmedian(LGPgap{2,c}(:,1:nsc));
        LGPgapstats{2,c}(3,1:nsc)=nanstd(LGPgap{2,c}(:,1:nsc));
        LGPgapstats{2,c}(4,1:nsc)=min(LGPgap{2,c}(:,1:nsc));
        LGPgapstats{2,c}(5,1:nsc)=max(LGPgap{2,c}(:,1:nsc));
        
    end
    
    
% Calculate changes of stats (change of avg, change of median)
     LGPactDeltastats(1,1:ncrop)=LGPact(1,1:ncrop);

    for c=1:ncrop
        for stat=1:2
        LGPactDeltastats{2,c}(stat,1:nsc)=(LGPactstats{2,c}(stat,1:nsc)-LGPactstats{2,c}(stat,1));    
        end
    end

clear c

% Calculate average/median maturity dates
MATactstats(1,1:ncrop)=LGPact(1,1:ncrop); 

for c=1:ncrop
    
MATactstats{2,c}(1,1:nsc)=sowing{2,c}(1,1:nsc)+(LGPactstats{2,c}(1,1:nsc)-1);
MATactstats{2,c}(2,1:nsc)=sowing{2,c}(1,1:nsc)+(LGPactstats{2,c}(2,1:nsc)-1);
end

clear c

% 8.3 Vizualize actual & potential  length of growing period 
%------------------------------------------------------------------------- 

f1=figure('name','Median LGPact changes'); %(boxplot= variation over different GCMs) 
        sub(1)=subplot(2,5,1,'fontsize',10);
        boxplot(LGPactDeltastats{2,maize}(2,2:nsc),groupmat(1,2:nsc),'grouporder',groupnames,'labels',groupnames);
        line(xlim,[0,0],'Color','k','LineStyle','--')
        ylabel('Change of growing cycle length (days)')
        axis([xlim, -50,10])
        set(gca,'box','off');
        title('Maize')
        
        sub(2)=subplot(2,5,2,'fontsize',10);
        boxplot(LGPactDeltastats{2,wwheat}(2,2:nsc),groupmat(1,2:nsc),'grouporder',groupnames,'labels',groupnames);
        line(xlim,[0,0],'Color','k','LineStyle','--')
        set(gca,'box','off');
        title('Winter Wheat')
        
        sub(3)=subplot(2,5,3,'fontsize',10);
        boxplot(LGPactDeltastats{2,potato}(2,2:nsc),groupmat(1,2:nsc),'grouporder',groupnames,'labels',groupnames);
        line(xlim,[0,0],'Color','k','LineStyle','--')
        set(gca,'box','off');
        title('Potato')
        
        sub(4)=subplot(2,5,4,'fontsize',10);
        boxplot(LGPactDeltastats{2,sugarbeet}(2,2:nsc),groupmat(1,2:nsc),'grouporder',groupnames,'labels',groupnames);
        line(xlim,[0,0],'Color','k','LineStyle','--')
        set(gca,'box','off');
        title('Sugarbeet')
        
        sub(5)=subplot(2,5,5,'fontsize',10);
        boxplot(LGPactDeltastats{2,pea}(2,2:nsc),groupmat(1,2:nsc),'grouporder',groupnames,'labels',groupnames);
        line(xlim,[0,0],'Color','k','LineStyle','--')
        set(gca,'box','off');
        title('Peas')    

        linkaxes(sub,'y')% link y axis of different plots (so that they change simultaneously
        
        h(1)=subplot(2,5,6,'fontsize',10);
        bar(LGPact{2,maize}(2,1))
        ylabel('Medain historical growing cycle length (days)')
        axis([xlim , 0 ,300])
        title('Maize')
        set(gca,'XTickLabel',{' '});
        
        h(2)=subplot(2,5,7,'fontsize',10);
        bar(LGPact{2,wwheat}(2,1))
        title('Winter Wheat')
        set(gca,'XTickLabel',{' '});
        
        h(3)=subplot(2,5,8,'fontsize',10);
        bar(LGPact{2,potato}(2,1))
        title('Potato')
        set(gca,'XTickLabel',{' '});
         
        h(4)=subplot(2,5,9,'fontsize',10);
        bar(LGPact{2,sugarbeet}(2,1))
        title('Sugarbeet')
        set(gca,'XTickLabel',{' '});
                
        h(5)=subplot(2,5,10,'fontsize',10);
        bar(LGPact{2,pea}(2,1))
        title('Peas')
        set(gca,'XTickLabel',{' '});
        
        linkaxes(h,'y')% link y axis of different plots (so that they change simultaneously

        clear sub h  

f11=figure('name','Median LGPact en LGP pot'); %(boxplot= variation over different GCMs) % PUBLICATION FIGURE

        sub(1)=subplot(2,5,1,'fontsize',10);
        boxplot(LGPpotDeltastats{2,maize}(2,2:nsc),groupmat(1,2:nsc),'grouporder',groupnames,'labels',groupnames,'colorgroup', groupnames,'colors',boxplotcolor);
        line(xlim,[0,0],'Color','k','LineStyle','--')
        a=xlim;
        text(a(1,1),1,[num2str(LGPpot{2,maize}(2,1)),' days'],'HorizontalAlignment','left')
        ylabel('Change length of growing period (days)')
        axis([xlim, -50,10])
        set(gca,'box','off');
        title('Maize')
        
        sub(2)=subplot(2,5,2,'fontsize',10);
        boxplot(LGPpotDeltastats{2,wwheat}(2,2:nsc),groupmat(1,2:nsc),'grouporder',groupnames,'labels',groupnames,'colorgroup', groupnames,'colors',boxplotcolor);
        line(xlim,[0,0],'Color','k','LineStyle','--')
        a=xlim;
        text(a(1,1),1,[num2str(LGPpot{2,wwheat}(2,1)),' days'],'HorizontalAlignment','left')
        set(gca,'box','off');
        title('Winter wheat')
        
        sub(3)=subplot(2,5,3,'fontsize',10);
        boxplot(LGPpotDeltastats{2,potato}(2,2:nsc),groupmat(1,2:nsc),'grouporder',groupnames,'labels',groupnames,'colorgroup', groupnames,'colors',boxplotcolor);
        line(xlim,[0,0],'Color','k','LineStyle','--')
        a=xlim;
        text(a(1,1),1,[num2str(LGPpot{2,potato}(2,1)),' days'],'HorizontalAlignment','left')
        set(gca,'box','off');
        title('Potato')
        
        sub(4)=subplot(2,5,4,'fontsize',10);
        boxplot(LGPpotDeltastats{2,sugarbeet}(2,2:nsc),groupmat(1,2:nsc),'grouporder',groupnames,'labels',groupnames,'colorgroup', groupnames,'colors',boxplotcolor);
        line(xlim,[0,0],'Color','k','LineStyle','--')
        a=xlim;
        text(a(1,1),1,[num2str(LGPpot{2,sugarbeet}(2,1)),' days'],'HorizontalAlignment','left')
        set(gca,'box','off');
        title('Sugar beet')
        
        sub(5)=subplot(2,5,5,'fontsize',10);
        boxplot(LGPpotDeltastats{2,pea}(2,2:nsc),groupmat(1,2:nsc),'grouporder',groupnames,'labels',groupnames,'colorgroup', groupnames,'colors',boxplotcolor);
        line(xlim,[0,0],'Color','k','LineStyle','--')
        a=xlim;
        text(a(1,1),1,[num2str(LGPpot{2,pea}(2,1)),' days'],'HorizontalAlignment','left')
        set(gca,'box','off');
        title('Peas')    
        
        linkaxes(sub,'y')% link y axis of different plots (so that they change simultaneously
        
        sub(1)=subplot(2,5,6,'fontsize',10);
        boxplot(LGPactDeltastats{2,maize}(2,2:nsc),groupmat(1,2:nsc),'grouporder',groupnames,'labels',groupnames,'colorgroup', groupnames,'colors',boxplotcolor);
        line(xlim,[0,0],'Color','k','LineStyle','--')
        a=xlim;
        text(a(1,1),1,[num2str(LGPact{2,maize}(2,1)),' days'],'HorizontalAlignment','left','fontsize',7)
        ylabel('Change of length of growing period (days)')
        axis([xlim, -50,10])
        set(gca,'box','off');
        title('Maize')
        
        sub(2)=subplot(2,5,7,'fontsize',10);
        boxplot(LGPactDeltastats{2,wwheat}(2,2:nsc),groupmat(1,2:nsc),'grouporder',groupnames,'labels',groupnames,'colorgroup', groupnames,'colors',boxplotcolor);
        line(xlim,[0,0],'Color','k','LineStyle','--')
        a=xlim;
        text(a(1,1),1,[num2str(LGPact{2,wwheat}(2,1)),' days'],'HorizontalAlignment','left','fontsize',7)
        set(gca,'box','off');
        title('Winter wheat')
        
        sub(3)=subplot(2,5,8,'fontsize',10);
        boxplot(LGPactDeltastats{2,potato}(2,2:nsc),groupmat(1,2:nsc),'grouporder',groupnames,'labels',groupnames,'colorgroup', groupnames,'colors',boxplotcolor);
        line(xlim,[0,0],'Color','k','LineStyle','--')
        a=xlim;
        text(a(1,1),1,[num2str(LGPact{2,potato}(2,1),3),' days'],'HorizontalAlignment','left','fontsize',7)
        set(gca,'box','off');
        title('Potato')
        
        sub(4)=subplot(2,5,9,'fontsize',10);
        boxplot(LGPactDeltastats{2,sugarbeet}(2,2:nsc),groupmat(1,2:nsc),'grouporder',groupnames,'labels',groupnames,'colorgroup', groupnames,'colors',boxplotcolor);
        line(xlim,[0,0],'Color','k','LineStyle','--')
        a=xlim;
        text(a(1,1),1,[num2str(LGPact{2,sugarbeet}(2,1)),' days'],'HorizontalAlignment','left','fontsize',7)
        set(gca,'box','off');
        title('Sugar beet')
        
        sub(5)=subplot(2,5,10,'fontsize',10);
        boxplot(LGPactDeltastats{2,pea}(2,2:nsc),groupmat(1,2:nsc),'grouporder',groupnames,'labels',groupnames,'colorgroup', groupnames,'colors',boxplotcolor);
        line(xlim,[0,0],'Color','k','LineStyle','--')
        a=xlim;
        text(a(1,1),1,[num2str(LGPact{2,pea}(2,1)),' days'],'HorizontalAlignment','left','fontsize',7)
        set(gca,'box','off');
        title('Peas')    

        linkaxes(sub,'y')% link y axis of different plots (so that they change simultaneously

        clear sub 
              
        
% 8.4 Vizualize LGPact with cumulative distribution function
%-------------------------------------------------------------------------

f2=figure('name','LGPact emperical CDF');
    subplot(3,2,1,'fontsize',10);
    P=NaN(nsc,1);
    for i=1:nsc
        P(i)=cdfplot(LGPact{2,maize}(:,i));
        hold on 
    end
    PP=cdfplot(LGPpot{2,maize}(:,1));
    set(P,{'Color'},colorstruct,{'LineStyle'},linesstruct,{'LineWidth'},linewstruct)
    set(PP,{'Color'},{'k'},{'LineStyle'},{'--'},{'LineWidth'},{1.5})
    ylabel('Cumulative probability','fontsize',8);
    xlabel('Actual LGP (days)','fontsize',8);
    title('Maize')
    set(gca,'box','off')
    grid off

    subplot(3,2,2,'fontsize',10);
    for i=1:nsc
        P(i)=cdfplot(LGPact{2,wwheat}(:,i));
        hold on 
    end
    PP=cdfplot(LGPpot{2,wwheat}(:,1));
    set(P,{'Color'},colorstruct,{'LineStyle'},linesstruct,{'LineWidth'},linewstruct)
    set(PP,{'Color'},{'k'},{'LineStyle'},{'--'},{'LineWidth'},{1.5})
    ylabel('Cumulative probability','fontsize',8);
    xlabel('Actual LGP (days)','fontsize',8);
    title('Winter Wheat')
    set(gca,'box','off') 
    grid off
    
    subplot(3,2,3,'fontsize',10);
    for i=1:nsc
        P(i)=cdfplot(LGPact{2,sugarbeet}(:,i));
        hold on 
    end
    PP=cdfplot(LGPpot{2,sugarbeet}(:,1));
    set(P,{'Color'},colorstruct,{'LineStyle'},linesstruct,{'LineWidth'},linewstruct)
    set(PP,{'Color'},{'k'},{'LineStyle'},{'--'},{'LineWidth'},{1.5})
    ylabel('Cumulative probability','fontsize',8);
    xlabel('Actual LGP (days)','fontsize',8);
    title('Sugar beet')
    set(gca,'box','off')
    grid off
    
    subplot(3,2,4,'fontsize',10);
    for i=1:nsc
        P(i)=cdfplot(LGPact{2,potato}(:,i));
        hold on 
    end    
    PP=cdfplot(LGPpot{2,potato}(:,1));
    set(P,{'Color'},colorstruct,{'LineStyle'},linesstruct,{'LineWidth'},linewstruct)
    set(PP,{'Color'},{'k'},{'LineStyle'},{'--'},{'LineWidth'},{1.5})
    ylabel('Cumulative probability','fontsize',8);
    xlabel('Actual LGP (days)','fontsize',8);
    title('Potato')
    set(gca,'box','off')
    grid off

    subplot(3,2,5,'fontsize',10);
    for i=1:nsc
        P(i)=cdfplot(LGPact{2,pea}(:,i));
        hold on 
    end    
    PP=cdfplot(LGPpot{2,pea}(:,1));
    set(P,{'Color'},colorstruct,{'LineStyle'},linesstruct,{'LineWidth'},linewstruct)
    set(PP,{'Color'},{'k'},{'LineStyle'},{'--'},{'LineWidth'},{1.5})
    ylabel('Cumulative probability','fontsize',8);
    xlabel('Actual LGP (days)','fontsize',8);
    title('Pea')
    set(gca,'box','off')
    grid off   
    
    clear P PP i 
    
f3=figure('name','LGPgap emperical CDF'); %#ok<*NASGU>
    subplot(3,2,1,'fontsize',10);
    P=NaN(nsc,1);
    for i=1:nsc
        P(i)=cdfplot(LGPgap{2,maize}(:,i));
        hold on 
    end
    set(P,{'Color'},colorstruct,{'LineStyle'},linesstruct,{'LineWidth'},linewstruct)
    ylabel('Cumulative probability','fontsize',8);
    xlabel('LGPgap (days)','fontsize',8);
    title('Maize')
    set(gca,'box','off')
    grid off

    subplot(3,2,2,'fontsize',10);
    for i=1:nsc
        P(i)=cdfplot(LGPgap{2,wwheat}(:,i));
        hold on 
    end
    set(P,{'Color'},colorstruct,{'LineStyle'},linesstruct,{'LineWidth'},linewstruct)
    ylabel('Cumulative probability','fontsize',8);
    xlabel('LGPgap (days)','fontsize',8);
    title('Winter Wheat')
    set(gca,'box','off') 
    grid off
    
    subplot(3,2,3,'fontsize',10);
    for i=1:nsc
        P(i)=cdfplot(LGPgap{2,sugarbeet}(:,i));
        hold on 
    end  
    set(P,{'Color'},colorstruct,{'LineStyle'},linesstruct,{'LineWidth'},linewstruct)
    ylabel('Cumulative probability','fontsize',8);
    xlabel('LGPgap (days)','fontsize',8);
    title('Sugar beet')
    set(gca,'box','off')
    grid off
    
    subplot(3,2,4,'fontsize',10);
    for i=1:nsc
        P(i)=cdfplot(LGPgap{2,potato}(:,i));
        hold on 
    end    
    set(P,{'Color'},colorstruct,{'LineStyle'},linesstruct,{'LineWidth'},linewstruct)
    ylabel('Cumulative probability','fontsize',8);
    xlabel('LGPgap (days)','fontsize',8);
    title('Potato')
    set(gca,'box','off')
    grid off

    subplot(3,2,5,'fontsize',10);
    for i=1:nsc
        P(i)=cdfplot(LGPgap{2,pea}(:,i));
        hold on 
    end    
    set(P,{'Color'},colorstruct,{'LineStyle'},linesstruct,{'LineWidth'},linewstruct)
    ylabel('Cumulative probability','fontsize',8);
     xlabel('LGPgap (days)','fontsize',8);
    title('Pea')
    set(gca,'box','off')
    grid off     
    
    clear P i            
              
% 8.5 Save vizualizations
%-------------------------------------------------------------------------
filename='LGPact median changes - gcmboxplot';
filename=fullfile(DatapathScenOut,filename);
savefig(f1,filename)
 
filename='LGPact empirical CDF';
filename=fullfile(DatapathScenOut,filename);
savefig(f2,filename)

clear f1 f2

%% -----------------------------------------------------------------------
% 9. TOTAL DISCHARGE IMPACT 
%------------------------------------------------------------------------

% 9.1 aggregation of results (subtotals and cumulative) 
% ----------------------------------------------------------------------- 

% cumulative over days
Q_MTFcum=cumsum(Q_MTF); 
Q_MBFcum=cumsum(Q_MBF); 
Q_MIFcum=cumsum(Q_MIF); 
Q_MOFcum=cumsum(Q_MOF); 

% subtotals
[Q_MTFyear,Q_MTFmonth,~,Q_MTFseason,Q_MTFymonth,~,Q_MTFyseason]=ClimSubtotal(Date(:,1),Q_MTF*f/area,'sum');
[Q_MBFyear,Q_MBFmonth,~,Q_MBFseason,Q_MBFymonth,~,Q_MBFyseason]=ClimSubtotal(Date(:,1),Q_MBF*f/area,'sum');
[Q_MIFyear,Q_MIFmonth,~,Q_MIFseason,Q_MIFymonth,~,Q_MIFyseason]=ClimSubtotal(Date(:,1),Q_MIF*f/area,'sum');
[Q_MOFyear,Q_MOFmonth,~,Q_MOFseason,Q_MOFymonth,~,Q_MOFyseason]=ClimSubtotal(Date(:,1),Q_MOF*f/area,'sum');


% 9.2 put data in matrix (not structure)
% ----------------------------------------------------------------------- 

    nyear=year(EndDate(1,1))-year(StartDate(1,1))+1;
    Q_MTFyear2=NaN(nyear,nsc);
    Q_MBFyear2=NaN(nyear,nsc);
    Q_MIFyear2=NaN(nyear,nsc);
    Q_MOFyear2=NaN(nyear,nsc);
    
    Q_MTFyseason1=NaN(nyear,nsc);
    Q_MTFyseason2=NaN(nyear,nsc);
    Q_MTFyseason3=NaN(nyear,nsc);
    Q_MTFyseason4=NaN(nyear,nsc);
    
    for sc=1:nsc
    Q_MTFyear2(:,sc)=Q_MTFyear{1,sc}(:,2);
    Q_MBFyear2(:,sc)=Q_MBFyear{1,sc}(:,2);
    Q_MIFyear2(:,sc)=Q_MIFyear{1,sc}(:,2);
    Q_MOFyear2(:,sc)=Q_MOFyear{1,sc}(:,2);
    
    Q_MTFyseason1(:,sc)=Q_MTFyseason{1,sc}(Q_MTFyseason{1,1}(:,2)==1,3);
    Q_MTFyseason2(:,sc)=Q_MTFyseason{1,sc}(Q_MTFyseason{1,1}(:,2)==2,3);
    Q_MTFyseason3(:,sc)=Q_MTFyseason{1,sc}(Q_MTFyseason{1,1}(:,2)==3,3);
    Q_MTFyseason4(:,sc)=Q_MTFyseason{1,sc}(Q_MTFyseason{1,1}(:,2)==4,3);
    end
        
    clear sc 
    
% 9.3 Analyse with focus on variation between GCMs
% ----------------------------------------------------------------------- 
% stats per GCM over different years 
    Q_MTFyearstats=NaN(7,nsc);
    Q_MBFyearstats=NaN(7,nsc);
    Q_MIFyearstats=NaN(7,nsc);        
    Q_MOFyearstats=NaN(7,nsc);
    Q_MTFyseason1stats=NaN(7,nsc);
    Q_MTFyseason3stats=NaN(7,nsc);
        
    Q_MTFyearstats(1,1:nsc)=mean(Q_MTFyear2(:,1:nsc));
    Q_MTFyearstats(2,1:nsc)=median(Q_MTFyear2(:,1:nsc));
    Q_MTFyearstats(3,1:nsc)=std(Q_MTFyear2(:,1:nsc));
    Q_MTFyearstats(4,1:nsc)=min(Q_MTFyear2(:,1:nsc));
    Q_MTFyearstats(5,1:nsc)=max(Q_MTFyear2(:,1:nsc));
    Q_MTFyearstats(6,1:nsc)=Q_MTFyearstats(3,1:nsc)./Q_MTFyearstats(1,1:nsc);
    Q_MTFyearstats(7,1:nsc)=Q_MTFyearstats(5,1:nsc)-Q_MTFyearstats(4,1:nsc);
    
    Q_MBFyearstats(1,1:nsc)=mean(Q_MBFyear2(:,1:nsc));
    Q_MBFyearstats(2,1:nsc)=median(Q_MBFyear2(:,1:nsc));
    Q_MBFyearstats(3,1:nsc)=std(Q_MBFyear2(:,1:nsc));
    Q_MBFyearstats(4,1:nsc)=min(Q_MBFyear2(:,1:nsc));
    Q_MBFyearstats(5,1:nsc)=max(Q_MBFyear2(:,1:nsc));
    Q_MBFyearstats(6,1:nsc)=Q_MBFyearstats(3,1:nsc)./Q_MBFyearstats(1,1:nsc);
    Q_MBFyearstats(7,1:nsc)=Q_MBFyearstats(5,1:nsc)-Q_MBFyearstats(4,1:nsc);
    
    Q_MIFyearstats(1,1:nsc)=mean(Q_MIFyear2(:,1:nsc));
    Q_MIFyearstats(2,1:nsc)=median(Q_MIFyear2(:,1:nsc));
    Q_MIFyearstats(3,1:nsc)=std(Q_MIFyear2(:,1:nsc));
    Q_MIFyearstats(4,1:nsc)=min(Q_MIFyear2(:,1:nsc));
    Q_MIFyearstats(5,1:nsc)=max(Q_MIFyear2(:,1:nsc));       
    Q_MIFyearstats(6,1:nsc)=Q_MIFyearstats(3,1:nsc)./Q_MIFyearstats(1,1:nsc);
    Q_MIFyearstats(7,1:nsc)=Q_MIFyearstats(5,1:nsc)-Q_MIFyearstats(4,1:nsc);
    
    Q_MOFyearstats(1,1:nsc)=mean(Q_MOFyear2(:,1:nsc));
    Q_MOFyearstats(2,1:nsc)=median(Q_MOFyear2(:,1:nsc));
    Q_MOFyearstats(3,1:nsc)=std(Q_MOFyear2(:,1:nsc));
    Q_MOFyearstats(4,1:nsc)=min(Q_MOFyear2(:,1:nsc));
    Q_MOFyearstats(5,1:nsc)=max(Q_MOFyear2(:,1:nsc));      
    Q_MOFyearstats(6,1:nsc)=Q_MOFyearstats(3,1:nsc)./Q_MOFyearstats(1,1:nsc);
    Q_MOFyearstats(7,1:nsc)=Q_MOFyearstats(5,1:nsc)-Q_MOFyearstats(4,1:nsc);
    
    Q_MTFyseason1stats(1,1:nsc)=mean(Q_MTFyseason1(:,1:nsc));
    Q_MTFyseason1stats(2,1:nsc)=median(Q_MTFyseason1(:,1:nsc));
    Q_MTFyseason1stats(3,1:nsc)=std(Q_MTFyseason1(:,1:nsc));
    Q_MTFyseason1stats(4,1:nsc)=min(Q_MTFyseason1(:,1:nsc));
    Q_MTFyseason1stats(5,1:nsc)=max(Q_MTFyseason1(:,1:nsc));
    Q_MTFyseason1stats(6,1:nsc)=Q_MTFyseason1stats(3,1:nsc)./Q_MTFyseason1stats(1,1:nsc);
    Q_MTFyseason1stats(7,1:nsc)=Q_MTFyseason1stats(5,1:nsc)-Q_MTFyseason1stats(4,1:nsc);
    
    Q_MTFyseason3stats(1,1:nsc)=mean(Q_MTFyseason3(:,1:nsc));
    Q_MTFyseason3stats(2,1:nsc)=median(Q_MTFyseason3(:,1:nsc));
    Q_MTFyseason3stats(3,1:nsc)=std(Q_MTFyseason3(:,1:nsc));
    Q_MTFyseason3stats(4,1:nsc)=min(Q_MTFyseason3(:,1:nsc));
    Q_MTFyseason3stats(5,1:nsc)=max(Q_MTFyseason3(:,1:nsc));
    Q_MTFyseason3stats(6,1:nsc)=Q_MTFyseason3stats(3,1:nsc)./Q_MTFyseason3stats(1,1:nsc);
    Q_MTFyseason3stats(7,1:nsc)=Q_MTFyseason3stats(5,1:nsc)-Q_MTFyseason3stats(4,1:nsc);
    
%  Changes of mean and median (as compared to historical mean and median) for each GCM
    Q_MTFyearDeltastats=NaN(2,nsc);
    Q_MBFyearDeltastats=NaN(2,nsc);
    Q_MIFyearDeltastats=NaN(2,nsc);
    Q_MOFyearDeltastats=NaN(2,nsc);    
    
    Q_MTFyseason1Deltastats=NaN(2,nsc);
    Q_MTFyseason3Deltastats=NaN(2,nsc);
        
    for stat=1:2
        Q_MTFyearDeltastats(stat,1:nsc)=(Q_MTFyearstats(stat,1:nsc)-Q_MTFyearstats(stat,1))./Q_MTFyearstats(stat,1);
        Q_MBFyearDeltastats(stat,1:nsc)=(Q_MBFyearstats(stat,1:nsc)-Q_MBFyearstats(stat,1))./Q_MBFyearstats(stat,1);
        Q_MIFyearDeltastats(stat,1:nsc)=(Q_MIFyearstats(stat,1:nsc)-Q_MIFyearstats(stat,1))./Q_MIFyearstats(stat,1);
        Q_MOFyearDeltastats(stat,1:nsc)=(Q_MOFyearstats(stat,1:nsc)-Q_MOFyearstats(stat,1))./Q_MOFyearstats(stat,1);

        Q_MTFyseason1Deltastats(stat,1:nsc)=(Q_MTFyseason1stats(stat,1:nsc)-Q_MTFyseason1stats(stat,1))./Q_MTFyseason1stats(stat,1);
        Q_MTFyseason3Deltastats(stat,1:nsc)=(Q_MTFyseason3stats(stat,1:nsc)-Q_MTFyseason3stats(stat,1))./Q_MTFyseason3stats(stat,1);       
    end

    clear stat
% vizualization 
    f1=figure('name','Cumulative discharge');
            subplot(2,2,1,'fontsize',10);
            plot(Date(:,1),Q_MTFcum(1:nTime,2:nsc)*f/area,'color',[0.6 0.6 0.6]);
            hold on
            plot(Date(:,1),Q_MTFcum(1:nTime,1)*f/area,'color','k','LineWidth',2);
            ylabel('Total flow (mm)')
            set(gca,'box','off');
            
            subplot(2,2,2,'fontsize',10);
            plot(Date(:,1),Q_MBFcum(1:nTime,2:nsc)*f/area,'color',[0.6 0.6 0.6]);
            hold on
            plot(Date(:,1),Q_MBFcum(1:nTime,1)*f/area,'color','k','LineWidth',2);
            ylabel('Baseflow (mm)')
            set(gca,'box','off');
            
            subplot(2,2,3,'fontsize',10);
            plot(Date(:,1),Q_MIFcum(1:nTime,2:nsc)*f/area,'color',[0.6 0.6 0.6]);
            hold on
            plot(Date(:,1),Q_MIFcum(1:nTime,1)*f/area,'color','k','LineWidth',2);
            ylabel('Interflow (mm)')
            set(gca,'box','off');
            
            subplot(2,2,4,'fontsize',10);
            plot(Date(:,1),Q_MOFcum(1:nTime,2:nsc)*f/area,'color',[0.6 0.6 0.6]);
            hold on
            plot(Date(:,1),Q_MOFcum(1:nTime,1)*f/area,'color','k','LineWidth',2);
            ylabel('Overland flow (mm)')
            set(gca,'box','off');    
    
     f2=figure('name','Median yearly discharge changes - GCM variation'); % (PUBLICATION FIGURE)
            boxplot(Q_MTFyearDeltastats(2,2:nsc)*100,groupmat(1,2:nsc),'grouporder',groupnames,'labels',groupnames,'colorgroup', groupnames,'colors',boxplotcolor,'symbol','+');
            line(xlim,[0,0],'Color','k','LineStyle','--')
            ylabel('Median annual flow change (%)')
            title('Total flow')
            axis([xlim, -15,40])
            set(gca,'box','off')
            
            sub(2)=subplot(1,4,2,'fontsize',10);
            boxplot(Q_MBFyearDeltastats(2,2:nsc)*100,groupmat(1,2:nsc),'grouporder',groupnames,'labels',groupnames,'colorgroup', groupnames,'colors',boxplotcolor,'symbol','+');
            line(xlim,[0,0],'Color','k','LineStyle','--')
            title('Baseflow')
            set(gca,'box','off','YTick',[])
            
            sub(3)=subplot(1,4,3,'fontsize',10);
            boxplot(Q_MIFyearDeltastats(2,2:nsc)*100,groupmat(1,2:nsc),'grouporder',groupnames,'labels',groupnames,'colorgroup', groupnames,'colors',boxplotcolor,'symbol','+');
            line(xlim,[0,0],'Color','k','LineStyle','--')
            title('Interflow')
            set(gca,'box','off','YTick',[])
            
            sub(4)=subplot(1,4,4,'fontsize',10);
            boxplot(Q_MOFyearDeltastats(2,2:nsc)*100,groupmat(1,2:nsc),'grouporder',groupnames,'labels',groupnames,'colorgroup', groupnames,'colors',boxplotcolor,'symbol','+');
            line(xlim,[0,0],'Color','k','LineStyle','--')
            title('Overland flow')
            set(gca,'box','off','YTick',[])
            
            linkaxes(sub,'y')      
            
            clear sub
            
% 9.4 Analyse with focus on GCM & interannual variability 
% -----------------------------------------------------------------------     
% stats for all GCMs & different years together (historical set excluded)
    Q_MTBIOFyearstats=NaN(5,4);
    
    Q_MTFyear2fut=Q_MTFyear2(:,2:nsc);%select all columns of future
    Q_MBFyear2fut=Q_MBFyear2(:,2:nsc);
    Q_MIFyear2fut=Q_MIFyear2(:,2:nsc);
    Q_MOFyear2fut=Q_MOFyear2(:,2:nsc);
    
    Q_MTBIOFyearstats(1,1)=mean(Q_MTFyear2fut(:));
    Q_MTBIOFyearstats(2,1)=median(Q_MTFyear2fut(:));
    Q_MTBIOFyearstats(3,1)=std(Q_MTFyear2fut(:));
    Q_MTBIOFyearstats(4,1)=min(Q_MTFyear2fut(:));
    Q_MTBIOFyearstats(5,1)=max(Q_MTFyear2fut(:));
    
    Q_MTBIOFyearstats(1,2)=mean(Q_MBFyear2fut(:));
    Q_MTBIOFyearstats(2,2)=median(Q_MBFyear2fut(:));
    Q_MTBIOFyearstats(3,2)=std(Q_MBFyear2fut(:));
    Q_MTBIOFyearstats(4,2)=min(Q_MBFyear2fut(:));
    Q_MTBIOFyearstats(5,2)=max(Q_MBFyear2fut(:));
    
    Q_MTBIOFyearstats(1,3)=mean(Q_MIFyear2fut(:));
    Q_MTBIOFyearstats(2,3)=median(Q_MIFyear2fut(:));
    Q_MTBIOFyearstats(3,3)=std(Q_MIFyear2fut(:));
    Q_MTBIOFyearstats(4,3)=min(Q_MIFyear2fut(:));
    Q_MTBIOFyearstats(5,3)=max(Q_MIFyear2fut(:));       
    
    Q_MTBIOFyearstats(1,4)=mean(Q_MOFyear2fut(:));
    Q_MTBIOFyearstats(2,4)=median(Q_MOFyear2fut(:));
    Q_MTBIOFyearstats(3,4)=std(Q_MOFyear2fut(:));
    Q_MTBIOFyearstats(4,4)=min(Q_MOFyear2fut(:));
    Q_MTBIOFyearstats(5,4)=max(Q_MOFyear2fut(:)); 
    
    clear Q_MOFyear2fut Q_MBFyear2fut Q_MTFyear2fut Q_MIFyear2fut

% Change of yearly discharge as compared to historical median for all years

    Q_MTFyearDelta=(Q_MTFyear2-Q_MTFyearstats(2,1))./Q_MTFyearstats(2,1);
    Q_MBFyearDelta=(Q_MBFyear2-Q_MBFyearstats(2,1))./Q_MBFyearstats(2,1);
    Q_MIFyearDelta=(Q_MIFyear2-Q_MIFyearstats(2,1))./Q_MIFyearstats(2,1);
    Q_MOFyearDelta=(Q_MOFyear2-Q_MOFyearstats(2,1))./Q_MOFyearstats(2,1);

% vizualization                  
f3=figure('name','Median yearly discharge absolute- GCM&year variation');% boxplot = variation over different GCMs & over 30 different year)   

        sub(1)=subplot(1,6,1:2,'fontsize',10);
        boxplot(Q_MTFyear2(:,1:nsc),groupmat(1,1:nsc),'grouporder',groupnames2,'labels',groupnames2)
        line(xlim,[Q_MTFyearstats(2,1),Q_MTFyearstats(2,1)],'Color','k','LineStyle','--')
        ylabel('Annual flow(mm/year)')
        title('total flow')
        axis([xlim, 0,550])
        set(gca,'box','off')

        sub(2)=subplot(1,5,3,'fontsize',10);
        boxplot(Q_MBFyear2(:,1:nsc),groupmat(1,1:nsc),'grouporder',groupnames2,'labels',groupnames2)
        line(xlim,[Q_MBFyearstats(2,1),Q_MBFyearstats(2,1)],'Color','k','LineStyle','--')
        title('baseflow')
        set(gca,'box','off','YTick',[])

        sub(3)=subplot(1,5,4,'fontsize',10);
        boxplot(Q_MIFyear2(:,1:nsc),groupmat(1,1:nsc),'grouporder',groupnames2,'labels',groupnames2)
        line(xlim,[Q_MIFyearstats(2,1),Q_MIFyearstats(2,1)],'Color','k','LineStyle','--')
        title('interflow')
        set(gca,'box','off','YTick',[])

        sub(4)=subplot(1,5,5,'fontsize',10);
        boxplot(Q_MOFyear2(:,1:nsc),groupmat(1,1:nsc),'grouporder',groupnames2,'labels',groupnames2)
        line(xlim,[Q_MOFyearstats(2,1),Q_MOFyearstats(2,1)],'Color','k','LineStyle','--')
        title('overland flow')
        set(gca,'box','off','YTick',[])

        linkaxes(sub,'y')
        clear sub

f4=figure('name','Median yearly discharge changes- GCM&year variation');% boxplot = variation over different GCMs & over 30 different year)   

        sub(1)=subplot(1,6,1:2,'fontsize',10);
        boxplot(Q_MTFyearDelta(:,1:nsc)*100,groupmat(1,1:nsc),'grouporder',groupnames2,'labels',groupnames2)
        line(xlim,[0,0],'Color','k','LineStyle','--')
        ylabel('Annual flow change from historical median (%)')
        title('total flow')
        axis([xlim, -100,100])
        set(gca,'box','off')

        sub(2)=subplot(1,5,3,'fontsize',10);
        boxplot(Q_MBFyearDelta(:,1:nsc)*100,groupmat(1,1:nsc),'grouporder',groupnames2,'labels',groupnames2)
        line(xlim,[0,0],'Color','k','LineStyle','--')
        title('baseflow')
        set(gca,'box','off','YTick',[])

        sub(3)=subplot(1,5,4,'fontsize',10);
        boxplot(Q_MIFyearDelta(:,1:nsc)*100,groupmat(1,1:nsc),'grouporder',groupnames2,'labels',groupnames2)
        line(xlim,[0,0],'Color','k','LineStyle','--')
        title('interflow')
        set(gca,'box','off','YTick',[])

        sub(4)=subplot(1,5,5,'fontsize',10);
        boxplot(Q_MOFyearDelta(:,1:nsc)*100,groupmat(1,1:nsc),'grouporder',groupnames2,'labels',groupnames2)
        line(xlim,[0,0],'Color','k','LineStyle','--')
        title('overland flow')
        set(gca,'box','off','YTick',[])

        linkaxes(sub,'y')  
        clear sub
        
% 9.5 Analyse probabilities with cumulative distribution function  
% -----------------------------------------------------------------------
% normality check 
[notnormalTF,~]=NormalityCheck(Q_MTFyear2,'lillie',0.05);
[notnormalBF,~]=NormalityCheck(Q_MBFyear2,'lillie',0.05);
[notnormalIF,~]=NormalityCheck(Q_MIFyear2,'lillie',0.05);
[notnormalOF,~]=NormalityCheck(Q_MOFyear2,'lillie',0.05);

if isempty(notnormalTF)==1 && isempty(notnormalBF)==1 && isempty(notnormalIF)==1 && isempty(notnormalOF)==1
    disp('All flow values for all scenarios are normally distributed')
else
    if isempty(notnormalTF)==0
    warning(['Total flow is not normally distributed for scenarios: ',num2str(notnormalTF.')])
    end
    
    if isempty(notnormalBF)==0
    warning(['Baseflow is not normally distributed for scenarios: ',num2str(notnormalBF.')])
    end
    
    if isempty(notnormalIF)==0
    warning(['Interflow is not normally distributed for scenarios: ',num2str(notnormalIF.')])
    end
    
    if isempty(notnormalOF)==0
    warning(['Overland flow is not normally distributed for scenarios: ',num2str(notnormalOF.')])
    end
end       
    
clear notnormalTF notnormalBF notnormalIF notnormalOF 

% fit theoretical normal distributions
xrangeTF=0:5:max(Q_MTFyear2(:));
xrangeBF=0:5:max(Q_MBFyear2(:));
xrangeIF=0:5:max(Q_MIFyear2(:));
xrangeOF=0:5:max(Q_MOFyear2(:));

probabilitiesTF=NaN(length(xrangeTF),nsc);
probabilitiesBF=NaN(length(xrangeBF),nsc);
probabilitiesIF=NaN(length(xrangeIF),nsc);
probabilitiesOF=NaN(length(xrangeOF),nsc);

for sc=1:nsc
pdsc=fitdist(Q_MTFyear2(:,sc),'Normal');
probabilitiesTF(:,sc)=cdf(pdsc,xrangeTF);

pdsc=fitdist(Q_MBFyear2(:,sc),'Normal');
probabilitiesBF(:,sc)=cdf(pdsc,xrangeBF);

pdsc=fitdist(Q_MIFyear2(:,sc),'Normal');
probabilitiesIF(:,sc)=cdf(pdsc,xrangeIF);

pdsc=fitdist(Q_MOFyear2(:,sc),'Normal');
probabilitiesOF(:,sc)=cdf(pdsc,xrangeOF);

end
clear pdsc sc

% vizualize
f5=figure('name','Yearly discharge theoretical CDF');
    subplot(2,2,1,'fontsize',10);
    P=plot(xrangeTF,probabilitiesTF(:,1:nsc));   
    set(P,{'Color'},colorstruct,{'LineStyle'},linesstruct,{'LineWidth'},linewstruct)
    ylabel('Cumulative probability','fontsize',8);
    xlabel('Annual total flow (mm/year)','fontsize',8);
    title('Total flow')
    set(gca,'box','off')

    subplot(2,2,2,'fontsize',10);
    P=plot(xrangeBF,probabilitiesBF(:,1:nsc));   
    set(P,{'Color'},colorstruct,{'LineStyle'},linesstruct,{'LineWidth'},linewstruct)
    ylabel('Cumulative probability','fontsize',8);
    xlabel('Annual baseflow (mm/year)','fontsize',8);
    title('Baseflow')
    set(gca,'box','off') 
    
    subplot(2,2,3,'fontsize',10);
    P=plot(xrangeIF,probabilitiesIF(:,1:nsc));   
    set(P,{'Color'},colorstruct,{'LineStyle'},linesstruct,{'LineWidth'},linewstruct)
    ylabel('Cumulative probability','fontsize',8);
    xlabel('Annual interflow (mm/year)','fontsize',8);
    title('Interflow')
    set(gca,'box','off')
    
    subplot(2,2,4,'fontsize',10);
    P=plot(xrangeOF,probabilitiesOF(:,1:nsc));   
    set(P,{'Color'},colorstruct,{'LineStyle'},linesstruct,{'LineWidth'},linewstruct)
    ylabel('Cumulative probability','fontsize',8);
    xlabel('Annual overland flow (mm/year)','fontsize',8);
    title('Overland flow')
    set(gca,'box','off')

    clear P
    
f6=figure('name','Yearly discharge emperical CDF');
    subplot(2,2,1,'fontsize',10);
    P=NaN(nsc,1);
    for i=1:nsc
        P(i)=cdfplot(Q_MTFyear2(:,i));
        hold on 
    end
    set(P,{'Color'},colorstruct,{'LineStyle'},linesstruct,{'LineWidth'},linewstruct)
    ylabel('Cumulative probability','fontsize',8);
    xlabel('Annual total flow (mm/year)','fontsize',8);
    title('Total flow')
    set(gca,'box','off')
    grid off

    subplot(2,2,2,'fontsize',10);
    for i=1:nsc
        P(i)=cdfplot(Q_MBFyear2(:,i));
        hold on 
    end
    set(P,{'Color'},colorstruct,{'LineStyle'},linesstruct,{'LineWidth'},linewstruct)
    ylabel('Cumulative probability','fontsize',8);
    xlabel('Annual baseflow (mm/year)','fontsize',8);
    title('Baseflow')
    set(gca,'box','off') 
    grid off
    
    subplot(2,2,3,'fontsize',10);
    for i=1:nsc
        P(i)=cdfplot(Q_MIFyear2(:,i));
        hold on 
    end  
    set(P,{'Color'},colorstruct,{'LineStyle'},linesstruct,{'LineWidth'},linewstruct)
    ylabel('Cumulative probability','fontsize',8);
    xlabel('Annual interflow (mm/year)','fontsize',8);
    title('Interflow')
    set(gca,'box','off')
    grid off
    
    subplot(2,2,4,'fontsize',10);
    for i=1:nsc
        P(i)=cdfplot(Q_MOFyear2(:,i));
        hold on 
    end    
    set(P,{'Color'},colorstruct,{'LineStyle'},linesstruct,{'LineWidth'},linewstruct)
    ylabel('Cumulative probability','fontsize',8);
    xlabel('Annual overland flow (mm/year)','fontsize',8);
    title('Overland flow')
    set(gca,'box','off')
    grid off
    clear P 
    
 f7=figure('name','Yearly discharge emperical CDF - simple');
    subplot(2,2,1,'fontsize',10);
    P=NaN(ngroup2,1);
    for g=1:ngroup2
        gindex=find(strcmp(groupmat(1,1:nsc),groupnames2{1,g})==1);
        P(g)=cdfplot(reshape(Q_MTFyear2(:,gindex),[],1)); 
        hold on 
    end
    set(P,{'Color'},colorstructg,{'LineStyle'},linesstructg,{'LineWidth'},linewstructg)
    ylabel('Cumulative probability','fontsize',8);
    xlabel('Annual total flow (mm/year)','fontsize',8);
    title('')
    set(gca,'box','off')
    grid off

    subplot(2,2,2,'fontsize',10);
    for g=1:ngroup2
        gindex=find(strcmp(groupmat(1,1:nsc),groupnames2{1,g})==1);
        P(g)=cdfplot(reshape(Q_MBFyear2(:,gindex),[],1)); 
        hold on 
    end
    set(P,{'Color'},colorstructg,{'LineStyle'},linesstructg,{'LineWidth'},linewstructg)
    ylabel('Cumulative probability','fontsize',8);
    xlabel('Annual baseflow (mm/year)','fontsize',8);
    title('')
    set(gca,'box','off') 
    grid off
    
    subplot(2,2,3,'fontsize',10);
    for g=1:ngroup2
        gindex=find(strcmp(groupmat(1,1:nsc),groupnames2{1,g})==1);
        P(g)=cdfplot(reshape(Q_MIFyear2(:,gindex),[],1)); %#ok<*FNDSB>
        hold on 
    end
    set(P,{'Color'},colorstructg,{'LineStyle'},linesstructg,{'LineWidth'},linewstructg)
    ylabel('Cumulative probability','fontsize',8);
    xlabel('Annual interflow (mm/year)','fontsize',8);
    title('')
    set(gca,'box','off')
    grid off
    
    subplot(2,2,4,'fontsize',10);
    for g=1:ngroup2
        gindex=find(strcmp(groupmat(1,1:nsc),groupnames2{1,g})==1);
        P(g)=cdfplot(reshape(Q_MOFyear2(:,gindex),[],1));
        hold on 
    end
    set(P,{'Color'},colorstructg,{'LineStyle'},linesstructg,{'LineWidth'},linewstructg)
    ylabel('Cumulative probability','fontsize',8);
    xlabel('Annual overland flow (mm/year)','fontsize',8);
    title('')
    set(gca,'box','off')
    grid off
  
 clear xrangeTF xrangeBF xrangeIF xrangeOF   P i
 
 f8=figure('name','Yearly & seasonal discharge emperical CDF - simple');
    subplot(2,2,1,'fontsize',10);
    P=NaN(ngroup2,1);
    for g=1:ngroup2
        gindex=find(strcmp(groupmat(1,1:nsc),groupnames2{1,g})==1);
        P(g)=cdfplot(reshape(Q_MTFyear2(:,gindex),[],1));
        hold on 
    end
    set(P,{'Color'},colorstructg,{'LineStyle'},linesstructg,{'LineWidth'},linewstructg)
    ylabel('Cumulative probability','fontsize',8);
    xlabel('Annual total flow (mm/year)','fontsize',8);
    title('')
    set(gca,'box','off')
    grid off
            
    subplot(2,2,3,'fontsize',10);
    for g=1:ngroup2
        gindex=find(strcmp(groupmat(1,1:nsc),groupnames2{1,g})==1);
        P(g)=cdfplot(reshape(Q_MTFyseason1(:,gindex),[],1));
        hold on 
    end
    set(P,{'Color'},colorstructg,{'LineStyle'},linesstructg,{'LineWidth'},linewstructg)
    ylabel('Cumulative probability','fontsize',8);
    xlabel('Winter total flow (mm/season)','fontsize',8);
    title('')
    set(gca,'box','off')
    grid off
                
    subplot(2,2,4,'fontsize',10);
    for g=1:ngroup2
        gindex=find(strcmp(groupmat(1,1:nsc),groupnames2{1,g})==1);
        P(g)=cdfplot(reshape(Q_MTFyseason3(:,gindex),[],1));
        hold on 
    end
    set(P,{'Color'},colorstructg,{'LineStyle'},linesstructg,{'LineWidth'},linewstructg)
    ylabel('Cumulative probability','fontsize',8);
    xlabel('Summer total flow (mm/season)','fontsize',8);
    title('')
    set(gca,'box','off')
    grid off   
    
 f9=figure('name','Yearly and seasonal discharge emperical CDF'); %PUBLICATION FIGURE
    subplot(2,2,1,'fontsize',10);
    for i=1:nsc
        P(i)=cdfplot(Q_MTFyear2(:,i));
        hold on 
    end
    set(P,{'Color'},colorstruct,{'LineStyle'},linesstruct,{'LineWidth'},linewstruct)
    ylabel('Cumulative probability','fontsize',8);
    xlabel('Annual total flow (mm/year)','fontsize',8);
    title(' ')
    set(gca,'box','off')
    grid off
    text(5,0.95,'(a)','HorizontalAlignment','left')
    
    subplot(2,2,3,'fontsize',10);
    for i=1:nsc
        P(i)=cdfplot(Q_MTFyseason1(:,i));
        hold on 
    end
    set(P,{'Color'},colorstruct,{'LineStyle'},linesstruct,{'LineWidth'},linewstruct)
    ylabel('Cumulative probability','fontsize',8);
    xlabel('Winter total flow (mm/season)','fontsize',8);
    set(gca,'box','off')
    title(' ')
    grid off   
    text(5,0.95,'(b)','HorizontalAlignment','left')
    
    subplot(2,2,4,'fontsize',10);
    for i=1:nsc
        P(i)=cdfplot(Q_MTFyseason3(:,i));
        hold on 
    end
    set(P,{'Color'},colorstruct,{'LineStyle'},linesstruct,{'LineWidth'},linewstruct)
    ylabel('Cumulative probability','fontsize',8);
    xlabel('Summer total flow (mm/season)','fontsize',8);
    set(gca,'box','off')
    title(' ')
    grid off   
    text(5,0.95,'(c)','HorizontalAlignment','left')
    
    clear P i 
    
% 9.6 save results (vizualizations)
% ----------------------------------------------------------------------- 
filename='Discharge - cumulative';
filename=fullfile(DatapathScenOut,filename);
savefig(f1,filename)

filename='Yearly discharge median changes-gcmboxplot';
filename=fullfile(DatapathScenOut,filename);
savefig(f2,filename)

filename='Yearly discharge median changes-yeargcmboxplot';
filename=fullfile(DatapathScenOut,filename);
savefig(f4,filename)

filename='Yearly discharge emperical CDF';
filename=fullfile(DatapathScenOut,filename);
savefig(f6,filename)

filename='Yearly discharge emperical CDF - simple';
filename=fullfile(DatapathScenOut,filename);
savefig(f7,filename)

clear filename 
clear f1 f2 f3 f4 f5 f6 f7 f8 f9

%% -----------------------------------------------------------------------
% 10. MONTHLY DISCHARGE IMPACT (Only QTF) 
%------------------------------------------------------------------------

% 10.1 reorganize data in one matrix  
% ----------------------------------------------------------------------- 
Q_MTFymonth2=cell(1,12);

    for m=1:12
        mindex=find(Q_MTFymonth{1,1}(:,2)==m);
        for sc =1:nsc
        Q_MTFymonth2{1,m}(:,sc)=Q_MTFymonth{1,sc}(mindex,3);
        end
    end
    
    clear m sc

% 10.2 Calculate stats
% -----------------------------------------------------------------------

% stats for each month and each scenario
    Q_MTFmonthstats=cell(2,12);
    Q_MTFmonthstats(1,1:12)={'jan','febr','march', 'april', 'may', 'june', 'july', 'august', 'september', 'october', 'november', 'december'};
    
    for m=1:12 
        Q_MTFmonthstats{2,m}(1,1:nsc)=mean(Q_MTFymonth2{1,m}(:,1:nsc));
        Q_MTFmonthstats{2,m}(2,1:nsc)=median(Q_MTFymonth2{1,m}(:,1:nsc));
        Q_MTFmonthstats{2,m}(3,1:nsc)=std(Q_MTFymonth2{1,m}(:,1:nsc));
        Q_MTFmonthstats{2,m}(4,1:nsc)=min(Q_MTFymonth2{1,m}(:,1:nsc));
        Q_MTFmonthstats{2,m}(5,1:nsc)=max(Q_MTFymonth2{1,m}(:,1:nsc));
    end
    clear m
    
%  Mean changes of statistics (mean and median)for each month and each
%  scenario
    Q_MTFmonthDeltastats=cell(2,12);
    Q_MTFmonthDeltastats(1,1:12)=Q_MTFmonthstats(1,1:12);
   
    for m=1:12 
        for stat=1:2
            Q_MTFmonthDeltastats{2,m}(stat,1:nsc)=(Q_MTFmonthstats{2,m}(stat,1:nsc)-Q_MTFmonthstats{2,m}(stat,1))./Q_MTFmonthstats{2,m}(stat,1);
        end
    end
    clear m 
    
% 10.3 show in graphs
% -----------------------------------------------------------------------    
f1= figure('name','Median monthly discharge changes');%(boxplot= variation over different GCMs)    
            sub(1)=subplot('Position',[0.05, 0.4, 0.055,0.55],'fontsize',10); 
            boxplot(Q_MTFmonthDeltastats{2,1}(2,2:nsc)*100,groupmat(1,2:nsc),'grouporder',groupnames,'labels',groupnames);
            line(xlim,[0,0],'Color','k','LineStyle','--')
            ylabel('Median monthly flow change (%)')
            xlabel('January')
            axis([xlim, -60,60])
            set(gca,'XTick',[])
            set(gca,'box','off')
            
            sub(2)=subplot('Position',[0.16, 0.4, 0.04,0.55],'fontsize',10); 
            boxplot(Q_MTFmonthDeltastats{2,2}(2,2:nsc)*100,groupmat(1,2:nsc),'grouporder',groupnames,'labels',groupnames);
            line(xlim,[0,0],'Color','k','LineStyle','--')
            xlabel('February')           
            set(gca,'YTick',[],'XTick',[])
            set(gca,'box','off')
            
            sub(3)=subplot('Position',[0.2208, 0.4, 0.04,0.55],'fontsize',10); 
            boxplot(Q_MTFmonthDeltastats{2,3}(2,2:nsc)*100,groupmat(1,2:nsc),'grouporder',groupnames,'labels',groupnames);
            line(xlim,[0,0],'Color','k','LineStyle','--')
            xlabel('March')
            set(gca,'YTick',[],'XTick',[])
            set(gca,'box','off')
            
            sub(4)=subplot('Position',[0.28, 0.4, 0.04,0.55],'fontsize',10); 
            boxplot(Q_MTFmonthDeltastats{2,4}(2,2:nsc)*100,groupmat(1,2:nsc),'grouporder',groupnames,'labels',groupnames);
            line(xlim,[0,0],'Color','k','LineStyle','--')
            xlabel('April')
            set(gca,'YTick',[],'XTick',[])
            set(gca,'box','off')

            sub(5)=subplot('Position',[0.3355, 0.4, 0.04,0.55],'fontsize',10); 
            boxplot(Q_MTFmonthDeltastats{2,5}(2,2:nsc)*100,groupmat(1,2:nsc),'grouporder',groupnames,'labels',groupnames);
            line(xlim,[0,0],'Color','k','LineStyle','--')
            xlabel('May')
             set(gca,'YTick',[],'XTick',[])
            set(gca,'box','off')
            
            sub(6)=subplot('Position',[0.401, 0.4, 0.04,0.55],'fontsize',10); 
            boxplot(Q_MTFmonthDeltastats{2,6}(2,2:nsc)*100,groupmat(1,2:nsc),'grouporder',groupnames,'labels',groupnames);
            line(xlim,[0,0],'Color','k','LineStyle','--')
            xlabel('June')
            set(gca,'YTick',[],'XTick',[])
            set(gca,'box','off')
            
            title('Changes of median monthly discharge')            
            sub(7)=subplot('Position',[0.46, 0.4, 0.04,0.55],'fontsize',10); 
            boxplot(Q_MTFmonthDeltastats{2,7}(2,2:nsc)*100,groupmat(1,2:nsc),'grouporder',groupnames,'labels',groupnames);
            line(xlim,[0,0],'Color','k','LineStyle','--')
            xlabel('July')
            set(gca,'YTick',[],'XTick',[])
            set(gca,'box','off')
          
            sub(8)=subplot('Position',[0.525, 0.4, 0.04,0.55],'fontsize',10); 
            boxplot(Q_MTFmonthDeltastats{2,8}(1,2:nsc)*100,groupmat(1,2:nsc),'grouporder',groupnames,'labels',groupnames);
            line(xlim,[0,0],'Color','k','LineStyle','--')
            xlabel('August')
            set(gca,'YTick',[],'XTick',[])
            set(gca,'box','off')
          
            sub(9)=subplot('Position',[0.587, 0.4, 0.04,0.55],'fontsize',10); 
            boxplot(Q_MTFmonthDeltastats{2,9}(2,2:nsc)*100,groupmat(1,2:nsc),'grouporder',groupnames,'labels',groupnames);
            line(xlim,[0,0],'Color','k','LineStyle','--')
            xlabel('September')
            set(gca,'YTick',[],'XTick',[])
            set(gca,'box','off')
              
            sub(10)=subplot('Position',[0.64, 0.4, 0.04,0.55],'fontsize',10); 
            boxplot(Q_MTFmonthDeltastats{2,10}(2,2:nsc)*100,groupmat(1,2:nsc),'grouporder',groupnames,'labels',groupnames);
            line(xlim,[0,0],'Color','k','LineStyle','--')
            xlabel('October')  
            set(gca,'YTick',[],'XTick',[])
            set(gca,'box','off')
          
            sub(11)=subplot('Position',[0.725, 0.4, 0.04,0.55],'fontsize',10); 
            boxplot(Q_MTFmonthDeltastats{2,11}(2,2:nsc)*100,groupmat(1,2:nsc),'grouporder',groupnames,'labels',groupnames);
            line(xlim,[0,0],'Color','k','LineStyle','--')
            xlabel('November')
            set(gca,'YTick',[],'XTick',[])
            set(gca,'box','off')
          
            sub(12)=subplot('Position',[0.785, 0.4, 0.04,0.55],'fontsize',10); 
            boxplot(Q_MTFmonthDeltastats{2,12}(2,2:nsc)*100,groupmat(1,2:nsc),'grouporder',groupnames,'labels',groupnames);
            line(xlim,[0,0],'Color','k','LineStyle','--')
            xlabel('December')
            set(gca,'YTick',[],'XTick',[])
            set(gca,'box','off')     
            linkaxes(sub,'y')% link y axis of different plots (so that they change simultaneously
            
            subplot('Position',[0.06, 0.1, 0.76,0.2],'fontsize',12);
                HistValue=NaN(12,1);
                for i=1:12
                HistValue(i,1)=Q_MTFmonthstats{2,i}(1,2);
                end
            bar(1:12,HistValue)
            ylabel('Median monthly flow (mm/month)','fontsize',10)
            axis([0.5,12.5, ylim])
            set(gca,'XTick',1:12)
            set(gca,'box','off')
            
            clear sub i HistValue
               
% 10.4 save vizualization
% ----------------------------------------------------------------------- 
filename='Monthly Discharge - median changes GCM boxlplots';
filename=fullfile(DatapathScenOut,filename);
savefig(f1,filename)            
            
clear filename f1

%% -----------------------------------------------------------------------
% 11. IMPACT ON ET coefficient
%-------------------------------------------------------------------------

% 11.1 Read climate input 
%-------------------------------------------------------------------------
%initialize climate variables
Tmin=NaN(nTime,nsc); % minimum temperature (°C)
Tmax=NaN(nTime,nsc); % maximum temperature (°C)
Rain=NaN(nTime,nsc); % rainfall (mm)
ETo=NaN(nTime,nsc);  % reference evapotranspiration (mm)

%read climate variables
for sc=1:nsc
Name=ScenarioName{1,sc};% Extract scenario name
DatapathACSC=fullfile(DatapathAC,Name);% Datapaths for this scenario

Tempstr=ReadACTempInput(DatapathACSC);
Tmin(:,sc)=Tempstr{2,1}(:,1);
Tmax(:,sc)=Tempstr{2,1}(:,2);
Rainstr=ReadACPluInput(DatapathACSC);
Rain(:,sc)=Rainstr{2,1};
ETostr=ReadACEToInput(DatapathACSC);
ETo(:,sc)=ETostr{2,1};
end

% 11.2 Calculate ET coefficient
%-------------------------------------------------------------------------
kET=ETaCatch./ETo;
%Wr2CatchRel=Wr2Catch./Wrmax;
%kET=(ETaCatch./ETo)./Wr2CatchRel;
kET(isnan(kET))=1; % all days with ET0 and ETa is zero kET is 1

% 11.3 Calculate subtotals
%-------------------------------------------------------------------------
[~,~,~,~,kETymonth,~,~]=ClimSubtotal(Date(:,1),kET,'mean');

% 11.4 Calculate stats
%------------------------------------------------------------------------- 
%reorganize data
kETymonth2=cell(1,12);

    for m=1:12
        mindex=find(kETymonth{1,1}(:,2)==m);
        for sc =1:nsc
        kETymonth2{1,m}(:,sc)=kETymonth{1,sc}(mindex,3);
        end
    end
    
    clear m sc

% stats for each month and each scenario
    kETmonthstats=cell(2,12);
    kETmonthstats(1,1:12)={'jan','febr','march', 'april', 'may', 'june', 'july', 'august', 'september', 'october', 'november', 'december'};
    
    for m=1:12 
        kETmonthstats{2,m}(1,1:nsc)=mean(kETymonth2{1,m}(:,1:nsc));
        kETmonthstats{2,m}(2,1:nsc)=median(kETymonth2{1,m}(:,1:nsc));
        kETmonthstats{2,m}(3,1:nsc)=std(kETymonth2{1,m}(:,1:nsc));
        kETmonthstats{2,m}(4,1:nsc)=min(kETymonth2{1,m}(:,1:nsc));
        kETmonthstats{2,m}(5,1:nsc)=max(kETymonth2{1,m}(:,1:nsc));
    end
    clear m
    
%  Mean changes of statistics (mean and median)for each month and each
%  scenario
    kETmonthDeltastats=cell(2,12);
    kETmonthDeltastats(1,1:12)=kETmonthstats(1,1:12);
   
    for m=1:12 
        for stat=1:2
            kETmonthDeltastats{2,m}(stat,1:nsc)=(kETmonthstats{2,m}(stat,1:nsc)-kETmonthstats{2,m}(stat,1));
            %kETmonthDeltastats{2,m}(stat,1:nsc)=(kETmonthstats{2,m}(stat,1:nsc)-kETmonthstats{2,m}(stat,1))./kETmonthstats{2,m}(stat,1);
        end
    end
    clear m 
    
% 11.5 show in graphs
% -----------------------------------------------------------------------    
f1= figure('name','Median monthly discharge changes');%(boxplot= variation over different GCMs)  
            sub(1)=subplot('Position',[0.05, 0.4, 0.055,0.55],'fontsize',10); 
            boxplot(kETmonthDeltastats{2,1}(2,2:nsc)*100,groupmat(1,2:nsc),'grouporder',groupnames,'labels',groupnames);
            line(xlim,[0,0],'Color','k','LineStyle','--')
            ylabel('Median monthly flow change (%)')
            xlabel('January')
            axis([xlim, -60,60])
            set(gca,'XTick',[])
            set(gca,'box','off')
            
            sub(2)=subplot('Position',[0.16, 0.4, 0.04,0.55],'fontsize',10); 
            boxplot(kETmonthDeltastats{2,2}(2,2:nsc)*100,groupmat(1,2:nsc),'grouporder',groupnames,'labels',groupnames);
            line(xlim,[0,0],'Color','k','LineStyle','--')
            xlabel('February')           
            set(gca,'YTick',[],'XTick',[])
            set(gca,'box','off')
            
            sub(3)=subplot('Position',[0.2208, 0.4, 0.04,0.55],'fontsize',10); 
            boxplot(kETmonthDeltastats{2,3}(2,2:nsc)*100,groupmat(1,2:nsc),'grouporder',groupnames,'labels',groupnames);
            line(xlim,[0,0],'Color','k','LineStyle','--')
            xlabel('March')
            set(gca,'YTick',[],'XTick',[])
            set(gca,'box','off')
            
            sub(4)=subplot('Position',[0.28, 0.4, 0.04,0.55],'fontsize',10); 
            boxplot(kETmonthDeltastats{2,4}(2,2:nsc)*100,groupmat(1,2:nsc),'grouporder',groupnames,'labels',groupnames);
            line(xlim,[0,0],'Color','k','LineStyle','--')
            xlabel('April')
            set(gca,'YTick',[],'XTick',[])
            set(gca,'box','off')

            sub(5)=subplot('Position',[0.3355, 0.4, 0.04,0.55],'fontsize',10); 
            boxplot(kETmonthDeltastats{2,5}(2,2:nsc)*100,groupmat(1,2:nsc),'grouporder',groupnames,'labels',groupnames);
            line(xlim,[0,0],'Color','k','LineStyle','--')
            xlabel('May')
             set(gca,'YTick',[],'XTick',[])
            set(gca,'box','off')
            
            sub(6)=subplot('Position',[0.401, 0.4, 0.04,0.55],'fontsize',10); 
            boxplot(kETmonthDeltastats{2,6}(2,2:nsc)*100,groupmat(1,2:nsc),'grouporder',groupnames,'labels',groupnames);
            line(xlim,[0,0],'Color','k','LineStyle','--')
            xlabel('June')
            set(gca,'YTick',[],'XTick',[])
            set(gca,'box','off')
            
            title('Changes of median monthly discharge')            
            sub(7)=subplot('Position',[0.46, 0.4, 0.04,0.55],'fontsize',10); 
            boxplot(kETmonthDeltastats{2,7}(2,2:nsc)*100,groupmat(1,2:nsc),'grouporder',groupnames,'labels',groupnames);
            line(xlim,[0,0],'Color','k','LineStyle','--')
            xlabel('July')
            set(gca,'YTick',[],'XTick',[])
            set(gca,'box','off')
          
            sub(8)=subplot('Position',[0.525, 0.4, 0.04,0.55],'fontsize',10); 
            boxplot(kETmonthDeltastats{2,8}(1,2:nsc)*100,groupmat(1,2:nsc),'grouporder',groupnames,'labels',groupnames);
            line(xlim,[0,0],'Color','k','LineStyle','--')
            xlabel('August')
            set(gca,'YTick',[],'XTick',[])
            set(gca,'box','off')
          
            sub(9)=subplot('Position',[0.587, 0.4, 0.04,0.55],'fontsize',10); 
            boxplot(kETmonthDeltastats{2,9}(2,2:nsc)*100,groupmat(1,2:nsc),'grouporder',groupnames,'labels',groupnames);
            line(xlim,[0,0],'Color','k','LineStyle','--')
            xlabel('September')
            set(gca,'YTick',[],'XTick',[])
            set(gca,'box','off')
              
            sub(10)=subplot('Position',[0.64, 0.4, 0.04,0.55],'fontsize',10); 
            boxplot(kETmonthDeltastats{2,10}(2,2:nsc)*100,groupmat(1,2:nsc),'grouporder',groupnames,'labels',groupnames);
            line(xlim,[0,0],'Color','k','LineStyle','--')
            xlabel('October')  
            set(gca,'YTick',[],'XTick',[])
            set(gca,'box','off')
          
            sub(11)=subplot('Position',[0.725, 0.4, 0.04,0.55],'fontsize',10); 
            boxplot(kETmonthDeltastats{2,11}(2,2:nsc)*100,groupmat(1,2:nsc),'grouporder',groupnames,'labels',groupnames);
            line(xlim,[0,0],'Color','k','LineStyle','--')
            xlabel('November')
            set(gca,'YTick',[],'XTick',[])
            set(gca,'box','off')
          
            sub(12)=subplot('Position',[0.785, 0.4, 0.04,0.55],'fontsize',10); 
            boxplot(kETmonthDeltastats{2,12}(2,2:nsc)*100,groupmat(1,2:nsc),'grouporder',groupnames,'labels',groupnames);
            line(xlim,[0,0],'Color','k','LineStyle','--')
            xlabel('December')
            set(gca,'YTick',[],'XTick',[])
            set(gca,'box','off')     
            linkaxes(sub,'y')% link y axis of different plots (so that they change simultaneously
            
            subplot('Position',[0.06, 0.1, 0.76,0.2],'fontsize',12);
                HistValue=NaN(12,1);
                for i=1:12
                HistValue(i,1)=kETmonthstats{2,i}(1,2);
                end
            bar(1:12,HistValue)
            ylabel('Median monthly flow (mm/month)','fontsize',10)
            axis([0.5,12.5, ylim])
            set(gca,'XTick',1:12)
            set(gca,'box','off')
            
            clear sub i HistValue
               
% 11.6 save vizualization
% ----------------------------------------------------------------------- 
filename='Monthly ET coefficient - median changes GCM boxlplots';
filename=fullfile(DatapathScenOut,filename);
savefig(f1,filename)            
            
clear filename f1


%% -----------------------------------------------------------------------
% 12. SAVE SUMMARY TABLES
%-------------------------------------------------------------------------

xlname='SummaryTables.xlsx';
filename = fullfile(DatapathScenOut,xlname);

% 12.1 Summary of annual discharge 
% -----------------------------------------------------------------------  
%  Baseline median and relative changes to median

letters = 'AGMSY';
letters2='BHNTZ';

for g=1:ngroup
    
    gindex=find(strcmp(groupmat(1,1:nsc),groupnames{1,g})==1);   
    xlswrite(filename,groupnames(g),'Discharge',[letters2(g) num2str(1)]);

    HeadersColumns={' ','Historical median (mm/year)','Future median change-min (%)','Future median change-med(%)','Future median change-max (%)'};
    HeadersRows={'Annual Qtot';'Annual QBF';'Annual QIF';'Annual QOF'};
    Row1=[Q_MTFyearstats(2,1),min(Q_MTFyearDeltastats(2,gindex))*100,median(Q_MTFyearDeltastats(2,gindex))*100,max(Q_MTFyearDeltastats(2,gindex))*100];
    Row2=[Q_MBFyearstats(2,1),min(Q_MBFyearDeltastats(2,gindex))*100,median(Q_MBFyearDeltastats(2,gindex))*100,max(Q_MBFyearDeltastats(2,gindex))*100];
    Row3=[Q_MIFyearstats(2,1),min(Q_MIFyearDeltastats(2,gindex))*100,median(Q_MIFyearDeltastats(2,gindex))*100,max(Q_MIFyearDeltastats(2,gindex))*100];
    Row4=[Q_MOFyearstats(2,1),min(Q_MOFyearDeltastats(2,gindex))*100,median(Q_MOFyearDeltastats(2,gindex))*100,max(Q_MOFyearDeltastats(2,gindex))*100];
    DataMatrix=[Row1;Row2;Row3;Row4];

    xlswrite(filename,HeadersColumns,'Discharge',[letters(g) num2str(2)]);
    xlswrite(filename,HeadersRows,'Discharge',[letters(g) num2str(3)]);
    xlswrite(filename,DataMatrix,'Discharge',[letters2(g) num2str(3)]);

    clear gindex

end

clear HeadersColumns HeadersRows
clear Row1 Row2 Row3 Row4 DataMatrix
clear g gindex letters letters2

% 12.2 Summary of monthly total flow 
% ----------------------------------------------------------------------- 
%  Relative changes to baseline median

xlswrite(filename,groupnames2,'QTFmonth','B1');

for m=1:12
    xlswrite(filename,Q_MTFmonthDeltastats(1,m),'QTFmonth',['A' num2str(m+1)]);
    DataMatrix=NaN(1,ngroup2);    
    for g=1:ngroup2
       gindex=find(strcmp(groupmat(1,1:nsc),groupnames2{1,g})==1);         
       DataMatrix(1,g)=median(Q_MTFmonthDeltastats{2,m}(2,gindex));
    end     
       xlswrite(filename,DataMatrix,'QTFmonth',['B' num2str(m+1)]);

end

clear m g gindex DataMatrix 

% 12.3 Summary of Crop development and production parameters
% -----------------------------------------------------------------------     
% Baseline median and absolute changes to median

name={'Maize','WinterWheat','Potato','Sugarbeet','Pea'};
HeadersColumns={'Historical median','Future median change - min','Future median change - median','Future median change - max'};

letters = 'AGMSY';
letters2='BHNTZ';

for g=1:ngroup
    
    gindex=find(strcmp(groupmat(1,1:nsc),groupnames{1,g})==1);
    
    xlswrite(filename,groupnames(g),'CropVar',[letters(g) num2str(1)]);
    xlswrite(filename,HeadersColumns,'CropVar',[letters2(g) num2str(2)]);
    row=3;

    for c=1:length(name)
        cindex=find(strcmp(Cropnames(1,1:ncrop),name{1,c})==1);
        HeadersRows={name{1,c};'Seasonal crop yield (ton/ha)';'WPET (kg/m³)';'LGPact (days)';'BWSI (%)';'TSI (%)'};
        row1=[Yactstats{2,cindex}(2,1),min(YactDeltastats2{2,cindex}(2,gindex)),median(YactDeltastats2{2,cindex}(2,gindex)),max(YactDeltastats2{2,cindex}(2,gindex))];
        row2=[WPstats{2,cindex}(2,1),min(WPDeltastats2{2,cindex}(2,gindex)),median(WPDeltastats2{2,cindex}(2,gindex)),max(WPDeltastats2{2,cindex}(2,gindex))];
        row3=[LGPactstats{2,cindex}(2,1),min(LGPactDeltastats{2,cindex}(2,gindex)),median(LGPactDeltastats{2,cindex}(2,gindex)),max(LGPactDeltastats{2,cindex}(2,gindex))];
        row4=[BWSIstats{2,cindex}(2,1),min(BWSIDeltastats{2,cindex}(2,gindex)),median(BWSIDeltastats{2,cindex}(2,gindex)),max(BWSIDeltastats{2,cindex}(2,gindex))];
        row5=[TSIstats{2,cindex}(2,1),min(TSIDeltastats{2,cindex}(2,gindex)),median(TSIDeltastats{2,cindex}(2,gindex)),max(TSIDeltastats{2,cindex}(2,gindex))];
        DataMatrix=[row1;row2;row3;row4;row5];
         
        xlswrite(filename,HeadersRows,'CropVar',[letters(g) num2str(row)]);
        xlswrite(filename,DataMatrix,'CropVar',[letters2(g) num2str(row+1)]);
        row=row+6  ;
    end
    clear gindex
end

clear HeadersColumns HeadersRows
clear row1 row2 row3 row4 row 5 DataMatrix row
clear c g letters letters2 name gindex cindex

% 12.4 Summary of LGP gap
% ----------------------------------------------------------------------- 
% Baseline median and future median

name={'Maize','WinterWheat','Potato','Sugarbeet','Pea'};

xlswrite(filename,groupnames2,'LGPgap','B1');


for c=1:length(name)
    cindex=find(strcmp(Cropnames(1,1:ncrop),name{1,c})==1);
    xlswrite(filename,name(1,c),'LGPgap',['A' num2str(c+1)]);
    DataMatrix=NaN(1,ngroup2);   
    for g=1:ngroup2
            gindex=find(strcmp(groupmat(1,1:nsc),groupnames2{1,g})==1);
                  
            DataMatrix(1,g)=median(LGPgapstats{2,cindex}(2,gindex));
    end     
            xlswrite(filename,DataMatrix,'LGPgap',['B' num2str(c+1)]);

end

clear c cindex g gindex DataMatrix 
clear name 

% 12.5 Summary of CV yield, WPET and flow
% ----------------------------------------------------------------------- 
% Baseline and future median

%YIELD
name={'Maize','WinterWheat','Potato','Sugarbeet','Pea'};

xlswrite(filename,{'Yield'},'CV','A1');
xlswrite(filename,groupnames2,'CV','B1');


for c=1:length(name)
    cindex=find(strcmp(Cropnames(1,1:ncrop),name{1,c})==1);
    xlswrite(filename,name(1,c),'CV',['A' num2str(c+1)]);
    DataMatrix=NaN(1,ngroup2);    
    for g=1:ngroup2
            gindex=find(strcmp(groupmat(1,1:nsc),groupnames2{1,g})==1);
                  
            DataMatrix(1,g)=median(Yactstats{2,cindex}(6,gindex));
    end     
            xlswrite(filename,DataMatrix,'CV',['B' num2str(c+1)]);

end

clear DataMatrix c g gindex cindex 

%WPET
xlswrite(filename,{'WPET'},'CV','A8');
xlswrite(filename,groupnames2,'CV','B8');
for c=1:length(name)
    cindex=find(strcmp(Cropnames(1,1:ncrop),name{1,c})==1);
    xlswrite(filename,name(1,c),'CV',['A' num2str(c+7+1)]);
    DataMatrix=NaN(1,ngroup2);   
    for g=1:ngroup2
            gindex=find(strcmp(groupmat(1,1:nsc),groupnames2{1,g})==1);
                  
            DataMatrix(1,g)=median(WPstats{2,cindex}(6,gindex));
    end     
            xlswrite(filename,DataMatrix,'CV',['B' num2str(c+7+1)]);

end

clear DataMatrix c g gindex cindex 

%FLOW
xlswrite(filename,{'FLOW'},'CV','A15');
xlswrite(filename,groupnames2,'CV','B15');
xlswrite(filename,{'QTF';'QBF';'QIF';'QOF'},'CV','A16');

row1=NaN(1,ngroup2);
row2=NaN(1,ngroup2);
row3=NaN(1,ngroup2);
row4=NaN(1,ngroup2);

    for g=1:ngroup2
       gindex=find(strcmp(groupmat(1,1:nsc),groupnames2{1,g})==1);
       row1(1,g)=median(Q_MTFyearstats(6,gindex));
       row2(1,g)=median(Q_MBFyearstats(6,gindex));
       row3(1,g)=median(Q_MIFyearstats(6,gindex));
       row4(1,g)=median(Q_MOFyearstats(6,gindex));
    end 
    
    DataMatrix=[row1;row2;row3;row4];
    
    xlswrite(filename,DataMatrix,'CV','B16');

clear c cindex g gindex DataMatrix 
clear row1 row2 row3 row4
clear name

% 12.6 Summary of variability (range) of yield, WPET and flow
% ----------------------------------------------------------------------- 
% Baseline and future median
name={'Maize','WinterWheat','Potato','Sugarbeet','Pea'};

%YIELD
xlswrite(filename,{'Yield'},'range','A1');
xlswrite(filename,groupnames2,'range','B1');

for c=1:length(name)
    cindex=find(strcmp(Cropnames(1,1:ncrop),name{1,c})==1);
    xlswrite(filename,name(1,c),'range',['A' num2str(c+1)]);
    DataMatrix=NaN(1,ngroup2);    
    for g=1:ngroup2
            gindex=find(strcmp(groupmat(1,1:nsc),groupnames2{1,g})==1);
                  
            DataMatrix(1,g)=median(Yactstats{2,cindex}(7,gindex));
    end     
            xlswrite(filename,DataMatrix,'range',['B' num2str(c+1)]);

end

clear DataMatrix c g gindex cindex 

%WPET
xlswrite(filename,{'WPET'},'range','A8');
xlswrite(filename,groupnames2,'range','B8');
for c=1:length(name)
    cindex=find(strcmp(Cropnames(1,1:ncrop),name{1,c})==1);
    xlswrite(filename,name(1,c),'range',['A' num2str(c+7+1)]);
    DataMatrix=NaN(1,ngroup2);     
    for g=1:ngroup2
            gindex=find(strcmp(groupmat(1,1:nsc),groupnames2{1,g})==1);
                  
            DataMatrix(1,g)=median(WPstats{2,cindex}(7,gindex));
    end     
            xlswrite(filename,DataMatrix,'range',['B' num2str(c+7+1)]);

end
clear DataMatrix c g gindex cindex 

%FLOW
xlswrite(filename,{'FLOW'},'range','A15');
xlswrite(filename,groupnames2,'range','B15');
xlswrite(filename,{'QTF';'QBF';'QIF';'QOF'},'range','A16');

row1=NaN(1,ngroup2);
row2=NaN(1,ngroup2);
row3=NaN(1,ngroup2);
row4=NaN(1,ngroup2);

    for g=1:ngroup2
       gindex=find(strcmp(groupmat(1,1:nsc),groupnames2{1,g})==1);
       row1(1,g)=median(Q_MTFyearstats(7,gindex));
       row2(1,g)=median(Q_MBFyearstats(7,gindex));
       row3(1,g)=median(Q_MIFyearstats(7,gindex));
       row4(1,g)=median(Q_MOFyearstats(7,gindex));
    end 
    
    DataMatrix=[row1;row2;row3;row4];
    
    xlswrite(filename,DataMatrix,'range','B16');

clear  c cindex g gindex DataMatrix 
clear row1 row2 row3 row4

% 12.7 Soil water balance
% ----------------------------------------------------------------------- 
%calculate water balances
WabalAll=NaN(8,nsc);     
WabalGroup=NaN(8,ngroup2); % each column is one group, for each group the sum of P, ETo, E, Tr, DP, RO, Storage of soil and bunds water, total

for sc=1:nsc
WabalAll(1,sc)=sum(Rain(:,sc));
WabalAll(2,sc)=sum(ETo(:,sc));
WabalAll(3,sc)=-sum(ECatch(:,sc));
WabalAll(4,sc)=-sum(TrCatch(:,sc));
WabalAll(5,sc)=-sum(DPCatch(:,sc));
WabalAll(6,sc)=-sum(ROCatch(:,sc));
WabalAll(7,sc)=(BundWatCatch(end,sc)-BundWatCatch(1,sc))+(Wr2Catch(end,sc)-Wr2Catch(1,sc));
WabalAll(8,sc)=WabalAll(1,sc)+sum(WabalAll(3:7,sc));
end

for g=1:ngroup2
    gindex=find(strcmp(groupmat(1,1:nsc),groupnames2{1,g})==1);
    for comp=1:8
    WabalGroup(comp,g)=median(WabalAll(comp,gindex));
    end
end

% write away table
xlswrite(filename,groupnames2,'Wabal','B1');
xlswrite(filename,{'Rain';'ETo';'E';'Tr';'DP';'RO'; 'Storage';'Total'},'Wabal','A2:A9');
xlswrite(filename,WabalGroup,'Wabal','B2');

% 12.7 KET
% ----------------------------------------------------------------------- 
kETmat=NaN(12,nsc);
for m=1:12
kETmat(m,:)=kETmonthDeltastats{2,m}(2,:);
end

xlswrite(filename,kETmat,'kET','A1');

clear filename kETmat
%% -----------------------------------------------------------------------
% 13. PUBLICATION FIGURES
%-------------------------------------------------------------------------
close all 

% FIGURE 1: BOXPLOT ANNUAL FLOW & SUBFLOW CHANGES
%-------------------------------------------------------------------------
f1=figure('name','Figure 1 - Annual Flow boxplot');
    sub(1)=subplot(1,4,1,'fontsize',8);
    boxplot(Q_MTFyearDeltastats(2,2:nsc)*100,groupmat(1,2:nsc),'grouporder',groupnames,'labels',groupnames,'colorgroup', groupnames,'colors',boxplotcolor,'symbol','+');
    line(xlim,[0,0],'Color','k','LineStyle','--')
    ylabel('Change of median annual flow (%)','fontsize',8);
    title('Total flow','fontsize',8);
    axis([xlim, -15,40])
    set(gca,'box','off')

    sub(2)=subplot(1,4,2,'fontsize',10);
    boxplot(Q_MBFyearDeltastats(2,2:nsc)*100,groupmat(1,2:nsc),'grouporder',groupnames,'labels',groupnames,'colorgroup', groupnames,'colors',boxplotcolor,'symbol','+');
    line(xlim,[0,0],'Color','k','LineStyle','--')
    title('Baseflow','fontsize',8);
    set(gca,'box','off','YTick',[])

    sub(3)=subplot(1,4,3,'fontsize',10);
    boxplot(Q_MIFyearDeltastats(2,2:nsc)*100,groupmat(1,2:nsc),'grouporder',groupnames,'labels',groupnames,'colorgroup', groupnames,'colors',boxplotcolor,'symbol','+');
    line(xlim,[0,0],'Color','k','LineStyle','--')
    title('Interflow','fontsize',8);
    set(gca,'box','off','YTick',[])

    sub(4)=subplot(1,4,4,'fontsize',10);
    boxplot(Q_MOFyearDeltastats(2,2:nsc)*100,groupmat(1,2:nsc),'grouporder',groupnames,'labels',groupnames,'colorgroup', groupnames,'colors',boxplotcolor,'symbol','+');
    line(xlim,[0,0],'Color','k','LineStyle','--')
    title('Overland flow','fontsize',8);
    set(gca,'box','off','YTick',[])

    linkaxes(sub,'y')      

    clear sub
    fig=gcf;
    fig.PaperUnits='centimeters';
    fig.PaperPosition=[0 0 18 8];
    fig.PaperSize=[18 8];
    filename=fullfile(DatapathScenOut,'FlowBox_600dpi');
    filename2=fullfile(DatapathScenOut,'FlowBox_150dpi');
    print(filename,'-dpdf','-r600')
    print(filename2,'-dpdf','-r150') 

% FIGURE 2: CDF PLOTS ANNUAL FLOW & SUBFLOWS
%-------------------------------------------------------------------------
f2=figure('name','Figure 2 - Annual flow CDF');
    subplot(2,2,1,'fontsize',10);
    P=NaN(nsc,1);
    for i=1:nsc
        P(i)=cdfplot(Q_MTFyear2(:,i));
        hold on 
    end
    set(P,{'Color'},colorstruct,{'LineStyle'},linesstruct,{'LineWidth'},linewstruct)
    ylabel('Cumulative probability','fontsize',8);
    xlabel('Annual total flow (mm/y)','fontsize',8);
    title('','fontsize',8);
    set(gca,'box','off')
    grid off

    subplot(2,2,2,'fontsize',10);
    for i=1:nsc
        P(i)=cdfplot(Q_MBFyear2(:,i));
        hold on 
    end
    set(P,{'Color'},colorstruct,{'LineStyle'},linesstruct,{'LineWidth'},linewstruct)
    ylabel('Cumulative probability','fontsize',8);
    xlabel('Annual baseflow (mm/y)','fontsize',8);
    title('','fontsize',8);
    set(gca,'box','off') 
    grid off
    
    subplot(2,2,3,'fontsize',10);
    for i=1:nsc
        P(i)=cdfplot(Q_MIFyear2(:,i));
        hold on 
    end  
    set(P,{'Color'},colorstruct,{'LineStyle'},linesstruct,{'LineWidth'},linewstruct)
    ylabel('Cumulative probability','fontsize',8);
    xlabel('Annual interflow (mm/y)','fontsize',8);
    title('','fontsize',8);
    set(gca,'box','off')
    grid off
    
    subplot(2,2,4,'fontsize',10);
    for i=1:nsc
        P(i)=cdfplot(Q_MOFyear2(:,i));
        hold on 
    end    
    set(P,{'Color'},colorstruct,{'LineStyle'},linesstruct,{'LineWidth'},linewstruct)
    ylabel('Cumulative probability','fontsize',8);
    xlabel('Annual overland flow (mm/y)','fontsize',8);
    title('','fontsize',8);
    set(gca,'box','off')
    grid off
    
    for i=1:4
        subplot(2,2,i)
        text(1.2,0.97,['(',char(i+96),')'],'color','k','fontsize',8)
    end
    clear P i
    
    fig=gcf;
    fig.PaperUnits='centimeters';
    fig.PaperPosition=[0 0 14 11];
    fig.PaperSize=[14 11];
    filename=fullfile(DatapathScenOut,'FlowCDF_600dpi');
    filename2=fullfile(DatapathScenOut,'FlowCDF_150dpi');
    print(filename,'-dpdf','-r600')
    print(filename2,'-dpdf','-r150') 
%
% FIGURE 3: BOXPLOT MONTHLY TOTAL FLOW CHANGES
%-------------------------------------------------------------------------
f3=figure('name','Figure 3 - Montlhly flow boxplot');
    sub(1)=subplot('Position',[0.064, 0.15, 0.06,0.80],'fontsize',8); 
    boxplot(Q_MTFmonthDeltastats{2,1}(2,2:nsc)*100,groupmat(1,2:nsc),'grouporder',groupnames,'colorgroup',groupnames,'colors',boxplotcolor );
    line(xlim,[0,0],'Color','k','LineStyle','--')
    ylabel('Change of median monthly flow (%)','fontsize',10)
    xlabel('1')
    axis([xlim, -60,80])
    set(gca,'XTick',[])
    set(gca,'box','off')

    sub(2)=subplot('Position',[0.16, 0.15, 0.03,0.8],'fontsize',8); 
    boxplot(Q_MTFmonthDeltastats{2,2}(2,2:nsc)*100,groupmat(1,2:nsc),'grouporder',groupnames,'colorgroup',groupnames,'colors',boxplotcolor );
    line(xlim,[0,0],'Color','k','LineStyle','--')
    xlabel('2')           
    set(gca,'YTick',[],'XTick',[])
    set(gca,'box','off')

    sub(3)=subplot('Position',[0.23, 0.15, 0.03,0.8],'fontsize',8); 
    boxplot(Q_MTFmonthDeltastats{2,3}(2,2:nsc)*100,groupmat(1,2:nsc),'grouporder',groupnames,'colorgroup',groupnames,'colors',boxplotcolor );
    line(xlim,[0,0],'Color','k','LineStyle','--')
    xlabel('3')
    set(gca,'YTick',[],'XTick',[])
    set(gca,'box','off')

    sub(4)=subplot('Position',[0.29, 0.15, 0.03,0.8],'fontsize',10); 
    boxplot(Q_MTFmonthDeltastats{2,4}(2,2:nsc)*100,groupmat(1,2:nsc),'grouporder',groupnames,'colorgroup',groupnames,'colors',boxplotcolor );
    line(xlim,[0,0],'Color','k','LineStyle','--')
    xlabel('4')
    set(gca,'YTick',[],'XTick',[])
    set(gca,'box','off')

    sub(5)=subplot('Position',[0.36, 0.15, 0.03,0.8],'fontsize',10); 
    boxplot(Q_MTFmonthDeltastats{2,5}(2,2:nsc)*100,groupmat(1,2:nsc),'grouporder',groupnames,'colorgroup',groupnames,'colors',boxplotcolor );
    line(xlim,[0,0],'Color','k','LineStyle','--')
    xlabel('5')
    set(gca,'YTick',[],'XTick',[])
    set(gca,'box','off')

    sub(6)=subplot('Position',[0.43, 0.15, 0.03,0.8],'fontsize',10); 
    boxplot(Q_MTFmonthDeltastats{2,6}(2,2:nsc)*100,groupmat(1,2:nsc),'grouporder',groupnames,'colorgroup',groupnames,'colors',boxplotcolor );
    line(xlim,[0,0],'Color','k','LineStyle','--')
    xlabel('6')
    set(gca,'YTick',[],'XTick',[])
    set(gca,'box','off')

    sub(7)=subplot('Position',[0.50, 0.15, 0.03,0.8],'fontsize',10); 
    boxplot(Q_MTFmonthDeltastats{2,7}(2,2:nsc)*100,groupmat(1,2:nsc),'grouporder',groupnames,'colorgroup',groupnames,'colors',boxplotcolor );
    line(xlim,[0,0],'Color','k','LineStyle','--')
    xlabel('7')
    set(gca,'YTick',[],'XTick',[])
    set(gca,'box','off')

    sub(8)=subplot('Position',[0.57, 0.15, 0.03,0.8],'fontsize',10); 
    boxplot(Q_MTFmonthDeltastats{2,8}(2,2:nsc)*100,groupmat(1,2:nsc),'grouporder',groupnames,'colorgroup',groupnames,'colors',boxplotcolor );
    line(xlim,[0,0],'Color','k','LineStyle','--')
    xlabel('8')
    set(gca,'YTick',[],'XTick',[])
    set(gca,'box','off')

    sub(9)=subplot('Position',[0.64, 0.15, 0.03,0.8],'fontsize',10); 
    boxplot(Q_MTFmonthDeltastats{2,9}(2,2:nsc)*100,groupmat(1,2:nsc),'grouporder',groupnames,'colorgroup',groupnames,'colors',boxplotcolor );
    line(xlim,[0,0],'Color','k','LineStyle','--')
    xlabel('9')
    set(gca,'YTick',[],'XTick',[])
    set(gca,'box','off')

    sub(10)=subplot('Position',[0.71, 0.15, 0.03,0.8],'fontsize',10); 
    boxplot(Q_MTFmonthDeltastats{2,10}(2,2:nsc)*100,groupmat(1,2:nsc),'grouporder',groupnames,'colorgroup',groupnames,'colors',boxplotcolor );
    line(xlim,[0,0],'Color','k','LineStyle','--')
    xlabel('10')  
    set(gca,'YTick',[],'XTick',[])
    set(gca,'box','off')

    sub(11)=subplot('Position',[0.78, 0.15, 0.03,0.8],'fontsize',10); 
    boxplot(Q_MTFmonthDeltastats{2,11}(2,2:nsc)*100,groupmat(1,2:nsc),'grouporder',groupnames,'colorgroup',groupnames,'colors',boxplotcolor );
    line(xlim,[0,0],'Color','k','LineStyle','--')
    xlabel('11')
    set(gca,'YTick',[],'XTick',[])
    set(gca,'box','off')

    sub(12)=subplot('Position',[0.85, 0.15, 0.03,0.8],'fontsize',10); 
    boxplot(Q_MTFmonthDeltastats{2,12}(2,2:nsc)*100,groupmat(1,2:nsc),'grouporder',groupnames,'colorgroup',groupnames,'colors',boxplotcolor );
    line(xlim,[0,0],'Color','k','LineStyle','--')
    xlabel('12')
    set(gca,'YTick',[],'XTick',[])
    set(gca,'box','off')  
    
    text(-10,-78,'Month','VerticalAlignment', 'bottom','HorizontalAlignment', 'center')

    linkaxes(sub,'y')
    clear sub
    
    fig=gcf;
    fig.PaperUnits='centimeters';
    fig.PaperPosition=[0 0 16 9];
    fig.PaperSize=[16 9];
    filename=fullfile(DatapathScenOut,'FlowMonth_600dpi');
    filename2=fullfile(DatapathScenOut,'FlowMonth_150dpi');
    print(filename,'-dpdf','-r600')
    print(filename2,'-dpdf','-r150') 

% FIGURE 4: BOXPLOT/SCATTERPLOT YIELD AND WPET CHANGES (per crop)
%-------------------------------------------------------------------------
% f4=figure('name','Figure 4 - Productivity boxplot');
%         sub(1)=subplot(2,5,1,'fontsize',10);
%         boxplot(YactDeltastats{2,maize}(2,2:nsc)*100,groupmat(1,2:nsc),'grouporder',groupnames,'labels',groupnames,'colorgroup',groupnames,'colors',boxplotcolor,'symbol','+' );
%         line(xlim,[0,0],'Color','k','LineStyle','--')
%         ylabel('Median yield change (%)')
%         axis([xlim, -20,40])
%         set(gca,'box','off')
%         title('Maize','fontsize',8)
%         
%         sub(2)=subplot(2,5,2,'fontsize',10);
%         boxplot(YactDeltastats{2,wwheat}(2,2:nsc)*100,groupmat(1,2:nsc),'grouporder',groupnames,'labels',groupnames,'colorgroup',groupnames,'colors',boxplotcolor,'symbol','+' );
%         line(xlim,[0,0],'Color','k','LineStyle','--')
%         set(gca,'box','off')
%         title('Winter wheat','fontsize',8)
%         
%         sub(3)=subplot(2,5,3,'fontsize',10);
%         boxplot(YactDeltastats{2,potato}(2,2:nsc)*100,groupmat(1,2:nsc),'grouporder',groupnames,'labels',groupnames,'colorgroup',groupnames,'colors',boxplotcolor,'symbol','+' );
%         line(xlim,[0,0],'Color','k','LineStyle','--')
%         set(gca,'box','off')
%         title('Potato','fontsize',8)
%         
%         sub(4)=subplot(2,5,4,'fontsize',10);
%         boxplot(YactDeltastats{2,sugarbeet}(2,2:nsc)*100,groupmat(1,2:nsc),'grouporder',groupnames,'labels',groupnames,'colorgroup',groupnames,'colors',boxplotcolor,'symbol','+' );
%         line(xlim,[0,0],'Color','k','LineStyle','--')
%         set(gca,'box','off')
%         title('Sugar beet','fontsize',8)
%         
%         sub(5)=subplot(2,5,5,'fontsize',10);
%         boxplot(YactDeltastats{2,pea}(2,2:nsc)*100,groupmat(1,2:nsc),'grouporder',groupnames,'labels',groupnames,'colorgroup',groupnames,'colors',boxplotcolor,'symbol','+');
%         line(xlim,[0,0],'Color','k','LineStyle','--')
%         set(gca,'box','off')
%         title('Peas','fontsize',8)    
% 
%         linkaxes(sub,'y')% link y axis of different plots (so that they change simultaneously
%         clear sub
%         
%         sub(1)=subplot(2,5,6,'fontsize',10);
%         boxplot(WPDeltastats{2,maize}(2,2:nsc)*100,groupmat(1,2:nsc),'grouporder',groupnames,'labels',groupnames,'colorgroup',groupnames,'colors',boxplotcolor,'symbol','+' );
%         line(xlim,[0,0],'Color','k','LineStyle','--')
%         h1=ylabel('Median WP_{ET} change (%)');
%         set(h1,'interpreter','tex')
%         axis([xlim, -10,70])
%         set(gca,'box','off')
%         
%         sub(2)=subplot(2,5,7,'fontsize',10);
%         boxplot(WPDeltastats{2,wwheat}(2,2:nsc)*100,groupmat(1,2:nsc),'grouporder',groupnames,'labels',groupnames,'colorgroup',groupnames,'colors',boxplotcolor,'symbol','+' );
%         line(xlim,[0,0],'Color','k','LineStyle','--')
%         set(gca,'box','off')
%         
%         sub(3)=subplot(2,5,8,'fontsize',10);
%         boxplot(WPDeltastats{2,potato}(2,2:nsc)*100,groupmat(1,2:nsc),'grouporder',groupnames,'labels',groupnames,'colorgroup',groupnames,'colors',boxplotcolor,'symbol','+' );
%         line(xlim,[0,0],'Color','k','LineStyle','--')
%         set(gca,'box','off')
%         
%         sub(4)=subplot(2,5,9,'fontsize',10);
%         boxplot(WPDeltastats{2,sugarbeet}(2,2:nsc)*100,groupmat(1,2:nsc),'grouporder',groupnames,'labels',groupnames,'colorgroup',groupnames,'colors',boxplotcolor,'symbol','+' );
%         line(xlim,[0,0],'Color','k','LineStyle','--')
%         set(gca,'box','off')
%         
%         sub(5)=subplot(2,5,10,'fontsize',10);
%         boxplot(WPDeltastats{2,pea}(2,2:nsc)*100,groupmat(1,2:nsc),'grouporder',groupnames,'labels',groupnames,'colorgroup',groupnames,'colors',boxplotcolor,'symbol','+' );
%         line(xlim,[0,0],'Color','k','LineStyle','--')
%         set(gca,'box','off')
% 
%         linkaxes(sub,'y')% link y axis of different plots (so that they change simultaneously
%         clear sub 

f4=figure('name','Figure 4 - Productivity response scatter');
        col={[0,0,0],[0.6,0.6,0.6]};
        P=NaN(ngroup,1);
        sub(1)=subplot(3,2,1,'fontsize',8);
        line([-20,40],[-20,40],'Color','k','LineStyle','--')
        line([-20,40],[0,0],'Color','k')
        line([0,0],[-20,60],'Color','k')
        axis([-20,40,-20,60])
        hold on
        for g=1:ngroup % to display the median
            gindex=find(strcmp(groupmat(1,1:nsc),groupnames{1,g})==1);
            c=col{g};
            P(g)=scatter(YactDeltastats{2,maize}(2,gindex)*100,WPDeltastats{2,maize}(2,gindex)*100,20,c,'o');
            hold on
            P(g)=scatter(median(YactDeltastats{2,maize}(2,gindex)*100),median(WPDeltastats{2,maize}(2,gindex)*100),20,c,'o','filled');
            hold on
        end
        %xlabel('Median yield response(%)','fontsize',8)
        ylabel('Median WP_{ET} response (%)','fontsize',8)
        set(gca,'box','off')
        title('Maize','fontsize',7)
                 
        sub(2)=subplot(3,2,2,'fontsize',8);
        line([-20,40],[-20,40],'Color','k','LineStyle','--')
        line([-20,40],[0,0],'Color','k')
        line([0,0],[-20,60],'Color','k')
        hold on
        for g=1:ngroup % to display the median
            gindex=find(strcmp(groupmat(1,1:nsc),groupnames{1,g})==1);
            c=col{g};
            P(g)=scatter(YactDeltastats{2,wwheat}(2,gindex)*100,WPDeltastats{2,wwheat}(2,gindex)*100,20,c,'o');
            hold on
            P(g)=scatter(median(YactDeltastats{2,wwheat}(2,gindex)*100),median(WPDeltastats{2,wwheat}(2,gindex)*100),20,c,'o','filled');
            hold on
        end
       % xlabel('Median yield response(%)','fontsize',8)
       % ylabel('Median WP_{ET} response (%)','fontsize',8)
        set(gca,'box','off')
        title('Winter Wheat','fontsize',7)
        
        sub(3)=subplot(3,2,3,'fontsize',8);
        line([-20,40],[-20,40],'Color','k','LineStyle','--')
        line([-20,40],[0,0],'Color','k')
        line([0,0],[-20,60],'Color','k')
        hold on
        for g=1:ngroup % to display the median
            gindex=find(strcmp(groupmat(1,1:nsc),groupnames{1,g})==1);
            c=col{g};
            P(g)=scatter(YactDeltastats{2,potato}(2,gindex)*100,WPDeltastats{2,potato}(2,gindex)*100,20,c,'o');
            hold on
            P(g)=scatter(median(YactDeltastats{2,potato}(2,gindex)*100),median(WPDeltastats{2,potato}(2,gindex)*100),20,c,'o','filled');
            hold on
        end
        %xlabel('Median yield response(%)','fontsize',8)
        ylabel('Median WP_{ET} response (%)','fontsize',8)
        set(gca,'box','off')
        title('Potato','fontsize',7)
        
        sub(4)=subplot(3,2,4,'fontsize',8);
        line([-20,40],[-20,40],'Color','k','LineStyle','--')
        line([-20,40],[0,0],'Color','k')
        line([0,0],[-20,60],'Color','k')
        hold on
        for g=1:ngroup % to display the median
            gindex=find(strcmp(groupmat(1,1:nsc),groupnames{1,g})==1);
            c=col{g};
            P(g)=scatter(YactDeltastats{2,sugarbeet}(2,gindex)*100,WPDeltastats{2,sugarbeet}(2,gindex)*100,20,c,'o');
            hold on
            P(g)=scatter(median(YactDeltastats{2,sugarbeet}(2,gindex)*100),median(WPDeltastats{2,sugarbeet}(2,gindex)*100),20,c,'o','filled');
            hold on
        end
        xlabel('Median yield response (%)','fontsize',8)
        %ylabel('Median WP_{ET} response (%)','fontsize',8)
        set(gca,'box','off')
        title('Sugar beet','fontsize',7)
        
        sub(5)=subplot(3,2,5,'fontsize',8);
        line([-20,40],[-20,40],'Color','k','LineStyle','--')
        line([-20,40],[0,0],'Color','k')
        line([0,0],[-20,60],'Color','k')
        hold on
        for g=1:ngroup % to display the median
            gindex=find(strcmp(groupmat(1,1:nsc),groupnames{1,g})==1);
            c=col{g};
            P(g)=scatter(YactDeltastats{2,pea}(2,gindex)*100,WPDeltastats{2,pea}(2,gindex)*100,20,c,'o');
            hold on
            P(g)=scatter(median(YactDeltastats{2,pea}(2,gindex)*100),median(WPDeltastats{2,pea}(2,gindex)*100),20,c,'o','filled');
            hold on
        end
        xlabel('Median yield response (%)','fontsize',8)
        ylabel('Median WP_{ET} response (%)','fontsize',8)
        set(gca,'box','off')
        title('Peas','fontsize',7)
        
        linkaxes(sub,'xy')% 
        clear sub  
        
        fig=gcf;
        fig.PaperUnits='centimeters';
        fig.PaperPosition=[0 0 12 15];
        fig.PaperSize=[12 15];
        filename=fullfile(DatapathScenOut,'ProdResponse_600dpi');
        filename2=fullfile(DatapathScenOut,'ProdResponse_150dpi');
        print(filename,'-dpdf','-r600')
        print(filename2,'-dpdf','-r150')

% FIGURE 5: CDF PLOTS YIELD AND WPET (traditional man)
%-------------------------------------------------------------------------
f5=figure('name','Figure 5 - Productivity CDF');
    subplot(5,2,1,'fontsize',8);
    P=NaN(nsc,1);
    for i=1:nsc
        P(i)=cdfplot(Yact{2,maize}(:,i));
        hold on 
    end
    set(P,{'Color'},colorstruct,{'LineStyle'},linesstruct,{'LineWidth'},linewstruct)
    axis([0 15 0 1])
    ylabel('Cum. probability','fontsize',8);
    xlabel('','fontsize',8);
    title('Maize','fontsize',8)
    set(gca,'box','off')
    grid off

    subplot(5,2,3,'fontsize',8);
    for i=1:nsc
        P(i)=cdfplot(Yact{2,wwheat}(:,i));
        hold on 
    end
    set(P,{'Color'},colorstruct,{'LineStyle'},linesstruct,{'LineWidth'},linewstruct)
    ax = gca;
    ax.XTick = [0 5 10 15 20];
    axis([0 20 0 1])
    ylabel('Cum. probability','fontsize',8);
    xlabel('','fontsize',8);
    title('Winter wheat','fontsize',8)
    set(gca,'box','off') 
    grid off
    
    subplot(5,2,5,'fontsize',8);
    for i=1:nsc
        P(i)=cdfplot(Yact{2,sugarbeet}(:,i));
        hold on 
    end  
    set(P,{'Color'},colorstruct,{'LineStyle'},linesstruct,{'LineWidth'},linewstruct)
    axis([0 25 0 1])
    ylabel('Cum. probability','fontsize',8);
    xlabel('','fontsize',8);
    title('Sugar beet','fontsize',8)
    set(gca,'box','off')
    axis([0,25,0,1])
    grid off
    
    subplot(5,2,7,'fontsize',8);
    for i=1:nsc
        P(i)=cdfplot(Yact{2,potato}(:,i));
        hold on 
    end    
    set(P,{'Color'},colorstruct,{'LineStyle'},linesstruct,{'LineWidth'},linewstruct)
    axis([0 20 0 1])
    ylabel('Cum. probability','fontsize',8);
    xlabel('','fontsize',8);
    title('Potato','fontsize',8)
    set(gca,'box','off')
    grid off

    subplot(5,2,9,'fontsize',8);
    for i=1:nsc
        P(i)=cdfplot(Yact{2,pea}(:,i));
        hold on 
    end    
    set(P,{'Color'},colorstruct,{'LineStyle'},linesstruct,{'LineWidth'},linewstruct)
    axis([0 4 0 1])
    ylabel('Cum. probability','fontsize',8);
    xlabel('Seasonal yield (ton/ha)','fontsize',8);
    title('Peas','fontsize',8)
    set(gca,'box','off')
    grid off
    
    subplot(5,2,2,'fontsize',8);
    for i=1:nsc
        P(i)=cdfplot(WP{2,maize}(:,i));
        hold on 
    end
    set(P,{'Color'},colorstruct,{'LineStyle'},linesstruct,{'LineWidth'},linewstruct)
    axis([0 4 0 1])
    ylabel('','fontsize',8);
    xlabel('','fontsize',8);
    title('Maize','fontsize',8)
    set(gca,'box','off')
    grid off

    subplot(5,2,4,'fontsize',8);
    for i=1:nsc
        P(i)=cdfplot(WP{2,wwheat}(:,i));
        hold on 
    end
    set(P,{'Color'},colorstruct,{'LineStyle'},linesstruct,{'LineWidth'},linewstruct)
    ylabel('','fontsize',8);
    xlabel('','fontsize',8);
    title('Winter wheat','fontsize',8)
    set(gca,'box','off') 
    grid off
    
    subplot(5,2,6,'fontsize',8);
    for i=1:nsc
        P(i)=cdfplot(WP{2,sugarbeet}(:,i));
        hold on 
    end  
    set(P,{'Color'},colorstruct,{'LineStyle'},linesstruct,{'LineWidth'},linewstruct)
    axis([0 6 0 1])
    ylabel('','fontsize',8);
    xlabel('','fontsize',8);
    title('Sugar beet','fontsize',8)
    set(gca,'box','off')
    grid off
    
    subplot(5,2,8,'fontsize',8);
    for i=1:nsc
        P(i)=cdfplot(WP{2,potato}(:,i));
        hold on 
    end    
    set(P,{'Color'},colorstruct,{'LineStyle'},linesstruct,{'LineWidth'},linewstruct)
    ylabel('','fontsize',8);
    xlabel('','fontsize',8);
    title('Potato','fontsize',8)
    set(gca,'box','off')
    grid off

    subplot(5,2,10,'fontsize',6);
    for i=1:nsc
        P(i)=cdfplot(WP{2,pea}(:,i));
        hold on 
    end    
    set(P,{'Color'},colorstruct,{'LineStyle'},linesstruct,{'LineWidth'},linewstruct)
    axis([0 2 0 1])
    ylabel('','fontsize',8);
    xlabel('WP_{ET} (kg/m³)','fontsize',8);
    title('Peas','fontsize',8)
    set(gca,'box','off')
    grid off
     
    clear P i
    
    fig=gcf;
    fig.PaperUnits='centimeters';
    fig.PaperPosition=[0 0 14 18];
    fig.PaperSize=[14 18];
    filename=fullfile(DatapathScenOut,'ProdCDF_600dpi');
    filename2=fullfile(DatapathScenOut,'ProdCDF_150dpi');
    print(filename,'-dpdf','-r600')
    print(filename2,'-dpdf','-r150') 

% FIGURE 6: SIMPLIFIED CDF PLOTS YIELD AND WPET (2 management trategies compared)
%-------------------------------------------------------------------------
f6=figure('name','Figure 6 - Productivity CDF simple');
    subplot(5,2,1,'fontsize',10);
    P=NaN(ngroup2,1);
    for g=1:ngroup2
        gindex=find(strcmp(groupmat(1,1:nsc),groupnames2{1,g})==1);
        P(g)=cdfplot(reshape(Yact{2,maize}(:,gindex),[],1));
        hold on 
    end
    set(P,{'Color'},colorstructg,{'LineStyle'},linesstructg,{'LineWidth'},linewstructg)
    axis([0 15 0 1])
    ylabel('Cum. probability','fontsize',8);
    xlabel('','fontsize',8);
    title('Maize','fontsize',8)
    set(gca,'box','off')
    grid off

    subplot(5,2,3,'fontsize',10);
    for g=1:ngroup2
        gindex=find(strcmp(groupmat(1,1:nsc),groupnames2{1,g})==1);
        P(g)=cdfplot(reshape(Yact{2,wwheat}(:,gindex),[],1));
        hold on 
    end
    set(P,{'Color'},colorstructg,{'LineStyle'},linesstructg,{'LineWidth'},linewstructg)
    axis([0 20 0 1])
    ylabel('Cum. probability','fontsize',8);
    xlabel('','fontsize',8);
    title('Winter wheat','fontsize',8)
    set(gca,'box','off') 
    grid off
    
    subplot(5,2,5,'fontsize',10);
    for g=1:ngroup2
        gindex=find(strcmp(groupmat(1,1:nsc),groupnames2{1,g})==1);
        P(g)=cdfplot(reshape(Yact{2,sugarbeet}(:,gindex),[],1));
        hold on 
    end  
    set(P,{'Color'},colorstructg,{'LineStyle'},linesstructg,{'LineWidth'},linewstructg)
    axis([0 25 0 1])
    ylabel('Cum. probability','fontsize',8);
    xlabel('','fontsize',8);
    title('Sugar beet','fontsize',8)
    set(gca,'box','off')
    grid off
    
    subplot(5,2,7,'fontsize',10);
    for g=1:ngroup2
        gindex=find(strcmp(groupmat(1,1:nsc),groupnames2{1,g})==1);
        P(g)=cdfplot(reshape(Yact{2,potato}(:,gindex),[],1));
        hold on 
    end    
    set(P,{'Color'},colorstructg,{'LineStyle'},linesstructg,{'LineWidth'},linewstructg)
    ylabel('Cum. probability','fontsize',8);
    xlabel('','fontsize',8);
    title('Potato','fontsize',8)
    set(gca,'box','off')
    grid off

    subplot(5,2,9,'fontsize',10);
    for g=1:ngroup2
        gindex=find(strcmp(groupmat(1,1:nsc),groupnames2{1,g})==1);
        P(g)=cdfplot(reshape(Yact{2,pea}(:,gindex),[],1));
        hold on 
    end    
    set(P,{'Color'},colorstructg,{'LineStyle'},linesstructg,{'LineWidth'},linewstructg)
    ylabel('Cum. probability','fontsize',8);
    xlabel('Seasonal yield (ton/ha)','fontsize',8);
    title('Peas','fontsize',8)
    set(gca,'box','off')
    grid off
       
    subplot(5,2,2,'fontsize',10);
    for g=1:ngroup2
        gindex=find(strcmp(groupmat(1,1:nsc),groupnames2{1,g})==1);
        P(g)=cdfplot(reshape(WP{2,maize}(:,gindex),[],1));
        hold on 
    end
    set(P,{'Color'},colorstructg,{'LineStyle'},linesstructg,{'LineWidth'},linewstructg)
    axis([0 4 0 1])
    ylabel('','fontsize',8);
    xlabel('','fontsize',8);
    title('Maize','fontsize',8)
    set(gca,'box','off')
    grid off

    subplot(5,2,4,'fontsize',10);
    for g=1:ngroup2
        gindex=find(strcmp(groupmat(1,1:nsc),groupnames2{1,g})==1);
        P(g)=cdfplot(reshape(WP{2,wwheat}(:,gindex),[],1));
        hold on 
    end
    set(P,{'Color'},colorstructg,{'LineStyle'},linesstructg,{'LineWidth'},linewstructg)
    axis([0 6 0 1])
    ylabel('','fontsize',8);
    xlabel('','fontsize',8);
    title('Winter wheat','fontsize',8)
    set(gca,'box','off') 
    grid off
    
    subplot(5,2,6,'fontsize',10);
    for g=1:ngroup2
        gindex=find(strcmp(groupmat(1,1:nsc),groupnames2{1,g})==1);
        P(g)=cdfplot(reshape(WP{2,sugarbeet}(:,gindex),[],1));
        hold on 
    end  
    set(P,{'Color'},colorstructg,{'LineStyle'},linesstructg,{'LineWidth'},linewstructg)
    axis([0 6 0 1])
    ylabel('','fontsize',8);
    xlabel('','fontsize',8);
    title('Sugar beet','fontsize',8)
    set(gca,'box','off')
    grid off
    
    subplot(5,2,8,'fontsize',10);
    for g=1:ngroup2
        gindex=find(strcmp(groupmat(1,1:nsc),groupnames2{1,g})==1);
        P(g)=cdfplot(reshape(WP{2,potato}(:,gindex),[],1));
        hold on 
    end    
    set(P,{'Color'},colorstructg,{'LineStyle'},linesstructg,{'LineWidth'},linewstructg)
    axis([0 6 0 1])
    ylabel('','fontsize',8);
    xlabel('','fontsize',8);
    title('Potato','fontsize',8)
    set(gca,'box','off')
    grid off

    subplot(5,2,10,'fontsize',10);
    for g=1:ngroup2
        gindex=find(strcmp(groupmat(1,1:nsc),groupnames2{1,g})==1);
        P(g)=cdfplot(reshape(WP{2,pea}(:,gindex),[],1));
        hold on 
    end    
    set(P,{'Color'},colorstructg,{'LineStyle'},linesstructg,{'LineWidth'},linewstructg)
    axis([0 2 0 1])
    ylabel('','fontsize',8);
    xlabel('WP_{ET} (kg/m³)','fontsize',8);
    title('Peas','fontsize',8)
    set(gca,'box','off')
    grid off

    clear P g gindex
    
    fig=gcf;
    fig.PaperUnits='centimeters';
    fig.PaperPosition=[0 0 14 18];
    fig.PaperSize=[14 18];
    filename=fullfile(DatapathScenOut,'ProdCDFMan_600dpi');
    filename2=fullfile(DatapathScenOut,'ProdCDFMan_150dpi');
    print(filename,'-dpdf','-r600')
    print(filename2,'-dpdf','-r150') 
%
% FIGURE 7: BOXPLOT LGP (actual) CHANGES
%-------------------------------------------------------------------------
f7=figure('name','Figure 7 - LGP boxplot'); 
        sub(1)=subplot(1,5,1,'fontsize',8);
        boxplot(LGPactDeltastats{2,maize}(2,2:nsc),groupmat(1,2:nsc),'grouporder',groupnames,'labels',groupnames,'colorgroup', groupnames,'colors',boxplotcolor);
        line(xlim,[0,0],'Color','k','LineStyle','--')
        a=xlim;
        text(a(1,1),1,[num2str(LGPact{2,maize}(2,1)),' days'],'HorizontalAlignment','left','fontsize',8)
        ylabel('Change of median LGP (days)','fontsize',8)
        axis([xlim, -50,10])
        set(gca,'box','off');
        title('Maize','fontsize',8);
        
        sub(2)=subplot(1,5,2,'fontsize',8);
        boxplot(LGPactDeltastats{2,wwheat}(2,2:nsc),groupmat(1,2:nsc),'grouporder',groupnames,'labels',groupnames,'colorgroup', groupnames,'colors',boxplotcolor);
        line(xlim,[0,0],'Color','k','LineStyle','--')
        a=xlim;
        text(a(1,1),1,[num2str(LGPact{2,wwheat}(2,1)),' days'],'HorizontalAlignment','left','fontsize',8)
        set(gca,'box','off');
        title('Winter wheat','fontsize',8);
        
        sub(3)=subplot(1,5,3,'fontsize',8);
        boxplot(LGPactDeltastats{2,potato}(2,2:nsc),groupmat(1,2:nsc),'grouporder',groupnames,'labels',groupnames,'colorgroup', groupnames,'colors',boxplotcolor);
        line(xlim,[0,0],'Color','k','LineStyle','--')
        a=xlim;
        text(a(1,1),1,[num2str(LGPact{2,potato}(2,1),3),' days'],'HorizontalAlignment','left','fontsize',8)
        set(gca,'box','off');
        title('Potato','fontsize',8);
        
        sub(4)=subplot(1,5,4,'fontsize',8);
        boxplot(LGPactDeltastats{2,sugarbeet}(2,2:nsc),groupmat(1,2:nsc),'grouporder',groupnames,'labels',groupnames,'colorgroup', groupnames,'colors',boxplotcolor);
        line(xlim,[0,0],'Color','k','LineStyle','--')
        a=xlim;
        text(a(1,1),1,[num2str(LGPact{2,sugarbeet}(2,1)),' days'],'HorizontalAlignment','left','fontsize',8)
        set(gca,'box','off');
        title('Sugar beet','fontsize',8);
        
        sub(5)=subplot(1,5,5,'fontsize',8);
        boxplot(LGPactDeltastats{2,pea}(2,2:nsc),groupmat(1,2:nsc),'grouporder',groupnames,'labels',groupnames,'colorgroup', groupnames,'colors',boxplotcolor);
        line(xlim,[0,0],'Color','k','LineStyle','--')
        a=xlim;
        text(a(1,1),1,[num2str(LGPact{2,pea}(2,1)),' days'],'HorizontalAlignment','left','fontsize',8)
        set(gca,'box','off');
        title('Peas','fontsize',8);  

        linkaxes(sub,'y')
        clear sub a
        
        fig=gcf;
        fig.PaperUnits='centimeters';
        fig.PaperPosition=[0 0 16 8];
        fig.PaperSize=[16 8];
        filename=fullfile(DatapathScenOut,'LGP_600dpi');
        filename2=fullfile(DatapathScenOut,'LGP_150dpi');
        print(filename,'-dpdf','-r600')
        print(filename2,'-dpdf','-r150') 
%
% FIGURE 8: BOXPLOT WSI CHANGES 
%-------------------------------------------------------------------------
f8=figure('name','Figure 8 - WSI boxplot');
        sub(1)=subplot(1,5,1,'fontsize',8);
        boxplot(BWSIDeltastats{2,maize}(2,2:nsc),groupmat(1,2:nsc),'grouporder',groupnames,'labels',groupnames,'colorgroup',groupnames,'colors',boxplotcolor,'symbol','+');
        line(xlim,[0,0],'Color','k','LineStyle','--')
        ylabel('Change of median WSI (%)','fontsize',8)
        axis([xlim, -5,20])
        set(gca,'box','off')
        title('Maize','fontsize',8)
        
        sub(2)=subplot(1,5,2,'fontsize',8);
        boxplot(BWSIDeltastats{2,wwheat}(2,2:nsc),groupmat(1,2:nsc),'grouporder',groupnames,'labels',groupnames,'colorgroup',groupnames,'colors',boxplotcolor,'symbol','+' );
        line(xlim,[0,0],'Color','k','LineStyle','--')
        set(gca,'box','off')
        title('Winter wheat','fontsize',8)
        
        sub(3)=subplot(1,5,3,'fontsize',8);
        boxplot(BWSIDeltastats{2,potato}(2,2:nsc),groupmat(1,2:nsc),'grouporder',groupnames,'labels',groupnames,'colorgroup',groupnames,'colors',boxplotcolor,'symbol','+' );
        line(xlim,[0,0],'Color','k','LineStyle','--')
        set(gca,'box','off')
        title('Potato','fontsize',8)
        
        sub(4)=subplot(1,5,4,'fontsize',8);
        boxplot(BWSIDeltastats{2,sugarbeet}(2,2:nsc),groupmat(1,2:nsc),'grouporder',groupnames,'labels',groupnames,'colorgroup',groupnames,'colors',boxplotcolor,'symbol','+' );
        line(xlim,[0,0],'Color','k','LineStyle','--')
        set(gca,'box','off')
        title('Sugar beet','fontsize',8)
        
        sub(5)=subplot(1,5,5,'fontsize',8);
        boxplot(BWSIDeltastats{2,pea}(2,2:nsc),groupmat(1,2:nsc),'grouporder',groupnames,'labels',groupnames,'colorgroup',groupnames,'colors',boxplotcolor,'symbol','+' );
        line(xlim,[0,0],'Color','k','LineStyle','--')
        set(gca,'box','off')
        title('Peas','fontsize',8)    

        linkaxes(sub,'y')
        clear sub
        
        fig=gcf;
        fig.PaperUnits='centimeters';
        fig.PaperPosition=[0 0 16 8];
        fig.PaperSize=[16 8];
        filename=fullfile(DatapathScenOut,'WSI_600dpi');
        filename2=fullfile(DatapathScenOut,'WSI_150dpi');
        print(filename,'-dpdf','-r600')
        print(filename2,'-dpdf','-r150') 

% FIGURE 9: BOXPLOT TSI (CSI & HSI) CHANGES
%-------------------------------------------------------------------------        
f9=figure('name','Figure 9 - TSI boxplot');
        sub(1)=subplot(1,5,1,'fontsize',8);
        boxplot(TSIDeltastats{2,maize}(2,2:nsc),groupmat(1,2:nsc),'grouporder',groupnames,'labels',groupnames,'colorgroup',groupnames,'colors',boxplotcolor,'symbol','+' );
        line(xlim,[0,0],'Color','k','LineStyle','--')
        ylabel('Change of median CSI (%)','fontsize',8)
        axis([xlim, -20,5])
        set(gca,'box','off')
        title('Maize','fontsize',8)    
        
        sub(2)=subplot(1,5,2,'fontsize',8);
        boxplot(TSIDeltastats{2,wwheat}(2,2:nsc),groupmat(1,2:nsc),'grouporder',groupnames,'labels',groupnames,'colorgroup',groupnames,'colors',boxplotcolor,'symbol','+' );
        line(xlim,[0,0],'Color','k','LineStyle','--')
        set(gca,'box','off')
        title('Winter wheat','fontsize',8)    
        
        sub(3)=subplot(1,5,3,'fontsize',8);
        boxplot(TSIDeltastats{2,potato}(2,2:nsc),groupmat(1,2:nsc),'grouporder',groupnames,'labels',groupnames,'colorgroup',groupnames,'colors',boxplotcolor,'symbol','+' );
        line(xlim,[0,0],'Color','k','LineStyle','--')
        set(gca,'box','off')
        title('Potato','fontsize',8)    
        
        sub(4)=subplot(1,5,4,'fontsize',8);
        boxplot(TSIDeltastats{2,sugarbeet}(2,2:nsc),groupmat(1,2:nsc),'grouporder',groupnames,'labels',groupnames,'colorgroup',groupnames,'colors',boxplotcolor,'symbol','+' );
        line(xlim,[0,0],'Color','k','LineStyle','--')
        set(gca,'box','off')
        title('Sugar beet','fontsize',8)    
        
        sub(5)=subplot(1,5,5,'fontsize',8);
        boxplot(TSIDeltastats{2,pea}(2,2:nsc),groupmat(1,2:nsc),'grouporder',groupnames,'labels',groupnames,'colorgroup',groupnames,'colors',boxplotcolor,'symbol','+' );
        line(xlim,[0,0],'Color','k','LineStyle','--')
        set(gca,'box','off')
        title('Peas','fontsize',8)        

        linkaxes(sub,'y')
        clear sub
 
%         sub(1)=subplot(2,5,6,'fontsize',10);
%         boxplot(HSIDeltastats{2,maize}(2,2:nsc),groupmat(1,2:nsc),'grouporder',groupnames,'labels',groupnames,'colorgroup',groupnames,'colors',boxplotcolor,'symbol','+' );
%         line(xlim,[0,0],'Color','k','LineStyle','--')
%         ylabel('Median HSI change (%)','fontsize',8)
%         axis([xlim, -20,5])
%         set(gca,'box','off')
%         
%         sub(2)=subplot(2,5,7,'fontsize',10);
%         boxplot(HSIDeltastats{2,wwheat}(2,2:nsc),groupmat(1,2:nsc),'grouporder',groupnames,'labels',groupnames,'colorgroup',groupnames,'colors',boxplotcolor,'symbol','+' );
%         line(xlim,[0,0],'Color','k','LineStyle','--')
%         set(gca,'box','off')
%         
%         sub(3)=subplot(2,5,8,'fontsize',10);
%         boxplot(HSIDeltastats{2,potato}(2,2:nsc),groupmat(1,2:nsc),'grouporder',groupnames,'labels',groupnames,'colorgroup',groupnames,'colors',boxplotcolor,'symbol','+' );
%         line(xlim,[0,0],'Color','k','LineStyle','--')
%         set(gca,'box','off')
%         
%         sub(4)=subplot(2,5,9,'fontsize',10);
%         boxplot(HSIDeltastats{2,sugarbeet}(2,2:nsc),groupmat(1,2:nsc),'grouporder',groupnames,'labels',groupnames,'colorgroup',groupnames,'colors',boxplotcolor,'symbol','+' );
%         line(xlim,[0,0],'Color','k','LineStyle','--')
%         set(gca,'box','off')
%         
%         sub(5)=subplot(2,5,10,'fontsize',10);
%         boxplot(HSIDeltastats{2,pea}(2,2:nsc),groupmat(1,2:nsc),'grouporder',groupnames,'labels',groupnames,'colorgroup',groupnames,'colors',boxplotcolor,'symbol','+' );
%         line(xlim,[0,0],'Color','k','LineStyle','--')
%         set(gca,'box','off')
%         
%         linkaxes(sub,'y')
%         clear sub    

        fig=gcf;
        fig.PaperUnits='centimeters';
        fig.PaperPosition=[0 0 16 8];
        fig.PaperSize=[16 8];
        filename=fullfile(DatapathScenOut,'CSI_600dpi');
        filename2=fullfile(DatapathScenOut,'CSI_150dpi');
        print(filename,'-dpdf','-r600')
        print(filename2,'-dpdf','-r150') 
   
% FIGURE 10: BOXPLOT kET CHANGES 
%-------------------------------------------------------------------------           
 f10=figure('name','Figure 10 - Montlhly kET boxplot');
    sub(1)=subplot('Position',[0.063, 0.15, 0.06,0.80],'fontsize',8); 
    boxplot(kETmonthDeltastats{2,1}(2,2:nsc)*100,groupmat(1,2:nsc),'grouporder',groupnames,'colorgroup',groupnames,'colors',boxplotcolor,'symbol','+');
    line(xlim,[0,0],'Color','k','LineStyle','--')
    ylabel('Change of median kET (%)','fontsize',10)
    xlabel('1')
    axis([xlim, -30,5])
    set(gca,'XTick',[])
    set(gca,'box','off')

    sub(2)=subplot('Position',[0.17, 0.15, 0.03,0.8],'fontsize',8); 
    boxplot(kETmonthDeltastats{2,2}(2,2:nsc)*100,groupmat(1,2:nsc),'grouporder',groupnames,'colorgroup',groupnames,'colors',boxplotcolor,'symbol','+' );
    line(xlim,[0,0],'Color','k','LineStyle','--')
    xlabel('2')           
    set(gca,'YTick',[],'XTick',[])
    set(gca,'box','off')

    sub(3)=subplot('Position',[0.24, 0.15, 0.028,0.8],'fontsize',8); 
    boxplot(kETmonthDeltastats{2,3}(2,2:nsc)*100,groupmat(1,2:nsc),'grouporder',groupnames,'colorgroup',groupnames,'colors',boxplotcolor,'symbol','+' );
    line(xlim,[0,0],'Color','k','LineStyle','--')
    xlabel('3')
    set(gca,'YTick',[],'XTick',[])
    set(gca,'box','off')

    sub(4)=subplot('Position',[0.295, 0.15, 0.03,0.8],'fontsize',10); 
    boxplot(kETmonthDeltastats{2,4}(2,2:nsc)*100,groupmat(1,2:nsc),'grouporder',groupnames,'colorgroup',groupnames,'colors',boxplotcolor,'symbol','+' );
    line(xlim,[0,0],'Color','k','LineStyle','--')
    xlabel('4')
    set(gca,'YTick',[],'XTick',[])
    set(gca,'box','off')

    sub(5)=subplot('Position',[0.36, 0.15, 0.03,0.8],'fontsize',10); 
    boxplot(kETmonthDeltastats{2,5}(2,2:nsc)*100,groupmat(1,2:nsc),'grouporder',groupnames,'colorgroup',groupnames,'colors',boxplotcolor,'symbol','+' );
    line(xlim,[0,0],'Color','k','LineStyle','--')
    xlabel('5')
    set(gca,'YTick',[],'XTick',[])
    set(gca,'box','off')

    sub(6)=subplot('Position',[0.43, 0.15, 0.03,0.8],'fontsize',10); 
    boxplot(kETmonthDeltastats{2,6}(2,2:nsc)*100,groupmat(1,2:nsc),'grouporder',groupnames,'colorgroup',groupnames,'colors',boxplotcolor,'symbol','+' );
    line(xlim,[0,0],'Color','k','LineStyle','--')
    xlabel('6')
    set(gca,'YTick',[],'XTick',[])
    set(gca,'box','off')

    sub(7)=subplot('Position',[0.50, 0.15, 0.03,0.8],'fontsize',10); 
    boxplot(kETmonthDeltastats{2,7}(2,2:nsc)*100,groupmat(1,2:nsc),'grouporder',groupnames,'colorgroup',groupnames,'colors',boxplotcolor,'symbol','+' );
    line(xlim,[0,0],'Color','k','LineStyle','--')
    xlabel('7')
    set(gca,'YTick',[],'XTick',[])
    set(gca,'box','off')

    sub(8)=subplot('Position',[0.57, 0.15, 0.03,0.8],'fontsize',10); 
    boxplot(kETmonthDeltastats{2,8}(2,2:nsc)*100,groupmat(1,2:nsc),'grouporder',groupnames,'colorgroup',groupnames,'colors',boxplotcolor,'symbol','+' );
    line(xlim,[0,0],'Color','k','LineStyle','--')
    xlabel('8')
    set(gca,'YTick',[],'XTick',[])
    set(gca,'box','off')

    sub(9)=subplot('Position',[0.64, 0.15, 0.03,0.8],'fontsize',10); 
    boxplot(kETmonthDeltastats{2,9}(2,2:nsc)*100,groupmat(1,2:nsc),'grouporder',groupnames,'colorgroup',groupnames,'colors',boxplotcolor,'symbol','+' );
    line(xlim,[0,0],'Color','k','LineStyle','--')
    xlabel('9')
    set(gca,'YTick',[],'XTick',[])
    set(gca,'box','off')

    sub(10)=subplot('Position',[0.70, 0.15, 0.03,0.8],'fontsize',10); 
    boxplot(kETmonthDeltastats{2,10}(2,2:nsc)*100,groupmat(1,2:nsc),'grouporder',groupnames,'colorgroup',groupnames,'colors',boxplotcolor,'symbol','+' );
    line(xlim,[0,0],'Color','k','LineStyle','--')
    xlabel('10')  
    set(gca,'YTick',[],'XTick',[])
    set(gca,'box','off')

    sub(11)=subplot('Position',[0.78, 0.15, 0.03,0.8],'fontsize',10); 
    boxplot(kETmonthDeltastats{2,11}(2,2:nsc)*100,groupmat(1,2:nsc),'grouporder',groupnames,'colorgroup',groupnames,'colors',boxplotcolor,'symbol','+' );
    line(xlim,[0,0],'Color','k','LineStyle','--')
    xlabel('11')
    set(gca,'YTick',[],'XTick',[])
    set(gca,'box','off')

    sub(12)=subplot('Position',[0.85, 0.15, 0.03,0.8],'fontsize',10); 
    boxplot(kETmonthDeltastats{2,12}(2,2:nsc)*100,groupmat(1,2:nsc),'grouporder',groupnames,'colorgroup',groupnames,'colors',boxplotcolor,'symbol','+' );
    line(xlim,[0,0],'Color','k','LineStyle','--')
    xlabel('12')
    set(gca,'YTick',[],'XTick',[])
    set(gca,'box','off')  
    
    text(-8,-36,'Month','VerticalAlignment', 'bottom','HorizontalAlignment', 'center')

    linkaxes(sub,'y')
    clear sub   
    
    fig=gcf;
    fig.PaperUnits='centimeters';
    fig.PaperPosition=[0 0 16 9];
    fig.PaperSize=[16 9];
    filename=fullfile(DatapathScenOut,'KET_600dpi');
    filename2=fullfile(DatapathScenOut,'KET_150dpi');
    print(filename,'-dpdf','-r600')
    print(filename2,'-dpdf','-r150') 
%
%% FIGURE 11: CDF PLOTS ANNUAL and SEASONAL TOTAL FLOW (with synthetic scenarios)
%-------------------------------------------------------------------------
f11=figure('name','Figure 11 - Flow synth CDF'); %PUBLICATION FIGURE
    subplot(2,2,1,'fontsize',10);
    P=NaN(nsc,1);
    for i=1:nsc
        P(i)=cdfplot(Q_MTFyear2(:,i));
        hold on 
    end
    set(P,{'Color'},colorstruct,{'LineStyle'},linesstruct,{'LineWidth'},linewstruct)
    ylabel('Cumulative probability','fontsize',8);
    xlabel('Annual total flow (mm/year)','fontsize',8);
    title(' ')
    set(gca,'box','off')
    grid off
    text(5,0.95,'(a)','HorizontalAlignment','left')
    
    subplot(2,2,3,'fontsize',10);
    for i=1:nsc
        P(i)=cdfplot(Q_MTFyseason1(:,i));
        hold on 
    end
    set(P,{'Color'},colorstruct,{'LineStyle'},linesstruct,{'LineWidth'},linewstruct)
    ylabel('Cumulative probability','fontsize',8);
    xlabel('Winter total flow (mm/season)','fontsize',8);
    set(gca,'box','off')
    title(' ')
    grid off   
    text(5,0.95,'(b)','HorizontalAlignment','left')
    
    subplot(2,2,4,'fontsize',10);
    for i=1:nsc
        P(i)=cdfplot(Q_MTFyseason3(:,i));
        hold on 
    end
    set(P,{'Color'},colorstruct,{'LineStyle'},linesstruct,{'LineWidth'},linewstruct)
    ylabel('Cumulative probability','fontsize',8);
    xlabel('Summer total flow (mm/season)','fontsize',8);
    set(gca,'box','off')
    title(' ')
    grid off   
    text(295,0.95,'(c)','HorizontalAlignment','right')
    
    clear P i 
    
    fig=gcf;
    fig.PaperUnits='centimeters';
    fig.PaperPosition=[0 0 14 11];
    fig.PaperSize=[14 11];
    filename=fullfile(DatapathScenOut,'FlowCDFSynth_600dpi');
    filename2=fullfile(DatapathScenOut,'FlowCDFSynth_150dpi');
    print(filename,'-dpdf','-r600')
    print(filename2,'-dpdf','-r150') 

% FIGURE 12: BOXPLOT YIELD CHANGES (with synthetic scenarios)
%-------------------------------------------------------------------------
f12=figure('name','Figure 12 - Yield synth boxplot');
        sub(1)=subplot(1,5,1,'fontsize',8);
            gindex=find(strcmp(groupmat(1,1:nsc),'Ens')==1);
            boxplot(YactDeltastats{2,maize}(2,gindex)*100,'colors',[0.6 0.6 0.6],'symbol','+');
            ylabel('Change of median yield (%)','fontsize',8);
            axis([xlim, -40,45])
            set(gca,'box','off','XTickLabel',{' '})
            title('Maize','fontsize',8);
            hold on
            %line(xlim,[0,0],'Color','k','LineStyle','-')
            hold on
            gindex=find(strcmp(groupmat(1,1:nsc),'Hs')==1);
            line(xlim,[YactDeltastats{2,maize}(2,gindex)*100,YactDeltastats{2,maize}(2,gindex)*100],'Color','r','LineStyle','--')
            hold on
            gindex=find(strcmp(groupmat(1,1:nsc),'Hw')==1);
            line(xlim,[YactDeltastats{2,maize}(2,gindex)*100,YactDeltastats{2,maize}(2,gindex)*100],'Color','r','LineStyle',':')
            hold on
            gindex=find(strcmp(groupmat(1,1:nsc),'L')==1);
            line(xlim,[YactDeltastats{2,maize}(2,gindex)*100,YactDeltastats{2,maize}(2,gindex)*100],'Color','r','LineStyle','-.')        
            hold on
            gindex=find(strcmp(groupmat(1,1:nsc),'M')==1);
            line(xlim,[YactDeltastats{2,maize}(2,gindex)*100,YactDeltastats{2,maize}(2,gindex)*100],'Color','r','LineStyle','-')

        
        sub(2)=subplot(1,5,2,'fontsize',8);
            gindex=find(strcmp(groupmat(1,1:nsc),'Ens')==1);
            boxplot(YactDeltastats{2,wwheat}(2,gindex)*100,'colors',[0.6 0.6 0.6],'symbol','+');
            axis([xlim, min(YactDeltastats{2,wwheat}(2,:)*100)-5,max(YactDeltastats{2,wwheat}(2,:)*100)+5])
            set(gca,'box','off','XTickLabel',{' '},'YTickLabel',{' '})
            title('Winter wheat','fontsize',8);
            hold on
            line(xlim,[0,0],'Color','k','LineStyle','-')
            hold on
            gindex=find(strcmp(groupmat(1,1:nsc),'Hs')==1);
            line(xlim,[YactDeltastats{2,wwheat}(2,gindex)*100,YactDeltastats{2,wwheat}(2,gindex)*100],'Color','r','LineStyle','--')
            hold on
            gindex=find(strcmp(groupmat(1,1:nsc),'Hw')==1);
            line(xlim,[YactDeltastats{2,wwheat}(2,gindex)*100,YactDeltastats{2,wwheat}(2,gindex)*100],'Color','r','LineStyle',':')
            hold on
            gindex=find(strcmp(groupmat(1,1:nsc),'L')==1);
            line(xlim,[YactDeltastats{2,wwheat}(2,gindex)*100,YactDeltastats{2,wwheat}(2,gindex)*100],'Color','r','LineStyle','-.')        
            hold on
            gindex=find(strcmp(groupmat(1,1:nsc),'M')==1);
            line(xlim,[YactDeltastats{2,wwheat}(2,gindex)*100,YactDeltastats{2,wwheat}(2,gindex)*100],'Color','r','LineStyle','-')
        
        sub(3)=subplot(1,5,3,'fontsize',8);
            gindex=find(strcmp(groupmat(1,1:nsc),'Ens')==1);
            boxplot(YactDeltastats{2,potato}(2,gindex)*100,'colors',[0.6 0.6 0.6],'symbol','+');
            axis([xlim, min(YactDeltastats{2,potato}(2,:)*100)-5,max(YactDeltastats{2,potato}(2,:)*100)+5])
            set(gca,'box','off','XTickLabel',{' '},'YTickLabel',{' '})
            title('Potato','fontsize',8);
            hold on
            line(xlim,[0,0],'Color','k','LineStyle','-')
            hold on
            gindex=find(strcmp(groupmat(1,1:nsc),'Hs')==1);
            line(xlim,[YactDeltastats{2,potato}(2,gindex)*100,YactDeltastats{2,potato}(2,gindex)*100],'Color','r','LineStyle','--')
            hold on
            gindex=find(strcmp(groupmat(1,1:nsc),'Hw')==1);
            line(xlim,[YactDeltastats{2,potato}(2,gindex)*100,YactDeltastats{2,potato}(2,gindex)*100],'Color','r','LineStyle',':')
            hold on
            gindex=find(strcmp(groupmat(1,1:nsc),'L')==1);
            line(xlim,[YactDeltastats{2,potato}(2,gindex)*100,YactDeltastats{2,potato}(2,gindex)*100],'Color','r','LineStyle','-.')        
            hold on
            gindex=find(strcmp(groupmat(1,1:nsc),'M')==1);
            line(xlim,[YactDeltastats{2,potato}(2,gindex)*100,YactDeltastats{2,potato}(2,gindex)*100],'Color','r','LineStyle','-')
        
        sub(4)=subplot(1,5,4,'fontsize',8);
            gindex=find(strcmp(groupmat(1,1:nsc),'Ens')==1);
            boxplot(YactDeltastats{2,sugarbeet}(2,gindex)*100,'colors',[0.6 0.6 0.6],'symbol','+');
            axis([xlim, min(YactDeltastats{2,sugarbeet}(2,:)*100)-5,max(YactDeltastats{2,sugarbeet}(2,:)*100)+5])
            set(gca,'box','off','XTickLabel',{' '},'YTickLabel',{' '})
            title('Sugar beet','fontsize',8);
            hold on
            line(xlim,[0,0],'Color','k','LineStyle','-')
            hold on
            gindex=find(strcmp(groupmat(1,1:nsc),'Hs')==1);
            line(xlim,[YactDeltastats{2,sugarbeet}(2,gindex)*100,YactDeltastats{2,sugarbeet}(2,gindex)*100],'Color','r','LineStyle','--')
            hold on
            gindex=find(strcmp(groupmat(1,1:nsc),'Hw')==1);
            line(xlim,[YactDeltastats{2,sugarbeet}(2,gindex)*100,YactDeltastats{2,sugarbeet}(2,gindex)*100],'Color','r','LineStyle',':')
            hold on
            gindex=find(strcmp(groupmat(1,1:nsc),'L')==1);
            line(xlim,[YactDeltastats{2,sugarbeet}(2,gindex)*100,YactDeltastats{2,sugarbeet}(2,gindex)*100],'Color','r','LineStyle','-.')        
            hold on
            gindex=find(strcmp(groupmat(1,1:nsc),'M')==1);
            line(xlim,[YactDeltastats{2,sugarbeet}(2,gindex)*100,YactDeltastats{2,sugarbeet}(2,gindex)*100],'Color','r','LineStyle','-')
        
        sub(5)=subplot(1,5,5,'fontsize',8);
            gindex=find(strcmp(groupmat(1,1:nsc),'Ens')==1);
            boxplot(YactDeltastats{2,pea}(2,gindex)*100,'colors',[0.6 0.6 0.6],'symbol','+');
            axis([xlim, min(YactDeltastats{2,pea}(2,:)*100)-5,max(YactDeltastats{2,pea}(2,:)*100)+5])
            set(gca,'box','off','XTickLabel',{' '},'YTickLabel',{' '})
            title('Peas','fontsize',8);
            hold on
            line(xlim,[0,0],'Color','k','LineStyle','-')
            hold on
            gindex=find(strcmp(groupmat(1,1:nsc),'Hs')==1);
            line(xlim,[YactDeltastats{2,pea}(2,gindex)*100,YactDeltastats{2,pea}(2,gindex)*100],'Color','r','LineStyle','--')
            hold on
            gindex=find(strcmp(groupmat(1,1:nsc),'Hw')==1);
            line(xlim,[YactDeltastats{2,pea}(2,gindex)*100,YactDeltastats{2,pea}(2,gindex)*100],'Color','r','LineStyle',':')
            hold on
            gindex=find(strcmp(groupmat(1,1:nsc),'L')==1);
            line(xlim,[YactDeltastats{2,pea}(2,gindex)*100,YactDeltastats{2,pea}(2,gindex)*100],'Color','r','LineStyle','-.')        
            hold on
            gindex=find(strcmp(groupmat(1,1:nsc),'M')==1);
            line(xlim,[YactDeltastats{2,pea}(2,gindex)*100,YactDeltastats{2,pea}(2,gindex)*100],'Color','r','LineStyle','-')   
            
            linkaxes(sub,'y')
            clear sub gindex
 
            fig=f12;
            fig.PaperUnits='centimeters';
            fig.PaperPosition=[0 0 16 8];
            fig.PaperSize=[16 8];
            filename=fullfile(DatapathScenOut,'ProdBoxSynth_600dpi');
            filename2=fullfile(DatapathScenOut,'ProdBoxSynth_150dpi');
            print(filename,'-dpdf','-r600')
            print(filename2,'-dpdf','-r150') 
%%
% SAVE FIGURES
%-------------------------------------------------------------------------

filename=fullfile(DatapathScenOut,'Fig -  Annual Flow boxplot');
savefig(f1,filename)                      
filename=fullfile(DatapathScenOut,'Fig -   Annual flow CDF');
savefig(f2,filename)      
filename=fullfile(DatapathScenOut,'Fig -  Montlhly flow boxplot');
savefig(f3,filename)      
filename=fullfile(DatapathScenOut,'Fig -   Productivity response scatter');
savefig(f4,filename)      
filename=fullfile(DatapathScenOut,'Fig -   Productivity CDF');
savefig(f5,filename)      
filename=fullfile(DatapathScenOut,'Fig -   Productivity CDF simple');
savefig(f6,filename)      
filename=fullfile(DatapathScenOut,'Fig -   LGP boxplot');
savefig(f7,filename)      
filename=fullfile(DatapathScenOut,'Fig -   WSI boxplot');
savefig(f8,filename)      
filename=fullfile(DatapathScenOut,'Fig -   TSI boxplot');
savefig(f9,filename)    
filename=fullfile(DatapathScenOut,'Fig -   Monthly kET boxplot');
savefig(f10,filename)    
filename=fullfile(DatapathScenOut,'Fig -   Flow synth CDF');
savefig(f11,filename)      
filename=fullfile(DatapathScenOut,'Fig -   Yield synth boxplot');
savefig(f12,filename)    

clear f1 f2 f3 f4 f5 f6 f7 f8 f9 f10 f12 