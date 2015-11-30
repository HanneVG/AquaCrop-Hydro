% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  This is the script for the AquaCrop-Hydro model and its data processing
%  This scripts brings together AquaCrop simulation output and simulates catchment hydrology.
%  Moreover, this script vizualizes all input & results, and does a performance analysis by comparing sim results with observations. 
% 
%  Theoretical background as well as a case-study example can be consulted
%  in the article:
%  Van Gaelen, H, Willems, P., Diels J. & Raes D. 2016. Modeling water
%  availability from field to catchemnt scale with a simple
%  agro-hydrological model. Environmental Modelling & Software [Submitted]
%
%  Documentation on the script can be found in the directory "Information",
%  Examples of the required input textfiles can be found in the directory
%  "ExampleInput"
%
%  TO DO before running the script: 
%   1. Run AquaCrop simulations for all land units in the catchment
%   2. Ensure that all AquaCrop output files are numbered with format "01", "02",
%        "03", "10","20" according to the landunit number
%   3. Prepare all required input files for this script 
%       (see example in 'ExampleInput' directory)
%   4. Specify all required datapaths (see script section 0)
%   5. Adapt the script for the soil types that occur in the catchment
%        (see CatchmentOutput.m, section 5)
%   6. Specify the AquaCrop model you uses (norma versus plugin)
%        (see AquaCropHydro.m, section 1)
%
%  Author: Hanne Van Gaelen
%  Last update: 12/11/2015
%  Version: AquaCrop-Hydro version 2015/11
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
clear % clears all variables from a previous run
close all % close all figures of a previous run

%% ------------------------------------------------------------------------
% 0. DEFINE DATA PATHS                                                    %
%--------------------------------------------------------------------------      

% specify in which path all inputs for this script are stored including 
%    a) the observations (Climate.txt & Flow.txt)
%    b) additional information on AquaCrop simulations (SimInfo.txt)
%    c) the parameters for the hydrological model (Parameters.txt)
%    d) the maximum root depth for every landunit and sim run (Zrx.txt)
%    e) the soil parameters of each soil type present (SoilPar.txt)
    
     DatapathInput='C:\DATA_HanneV\~Onderzoek\DEEL 2 - Ecohydro model\DEEL IIB  -VHM-AC\Testing\1. Kruishoutem\Models\AquaCrop-Hydro\Input\';

% specify in which path all AquaCrop simulation output files are stored

     DatapathAC='C:\DATA_HanneV\~Onderzoek\DEEL 2 - Ecohydro model\DEEL IIB  -VHM-AC\Testing\1. Kruishoutem\Simulation results\AquaCrop\Catchment sim files\Normal Ksat - Normal runoff';      

% specify in which path all outputs of this script should be stored
% including
%    a) performance statistics
%    b) yield data for every AquaCrop simulation unit and run

     DatapathOutput='C:\DATA_HanneV\~Onderzoek\DEEL 2 - Ecohydro model\DEEL IIB  -VHM-AC\Testing\1. Kruishoutem\Models\AquaCrop-Hydro\Output\';

%% ------------------------------------------------------------------------ 
% 1. LOAD ALL OBSERVATIONS                                                %
%--------------------------------------------------------------------------      
   
%%%% 1.1 Load AquaCrop climate observations (those are used as AquaCrop model input)
       name='Climate.txt';
       file = [DatapathInput name];
       A = importdata(file);
        
       Time1=A(:,1);    % Timesteps (days)
       Date1=A(:,2);    % Datenumbers (excel format)
       Prec=A(:,3);     % Rainfall (mm)
       ETo=A(:,4);      % Reference evapotranspiration (mm)
       Tmin=A(:,5);     % Daily minimum temperature (°C)
       Tmax=A(:,6);     % Daily maximum temperature (°C)
       
       a=size(Date1);
       nAC=a(1,1);
       StartDateAC=Date1(1,1);
       StopDateAC=Date1(end,1);
       
       clear a A name file

%%%% 1.2 Load de flow observation 
       name='Flow.txt';
       file = [DatapathInput name];
       F = importdata(file);
  
       Time2=F(:,1);    % timesteps (days)
       Date2=F(:,2);    % datenumbers (excel format)
       ObsTot=F(:,3);   % observed total discharge (m³/s)
       FiltBF=F(:,4);   % filtered baseflow (m³/s)
       FiltIF=F(:,5);   % filtered interflow (m³/s)
       FiltOF=F(:,6);   % filtered overland flow (m³/s)
       FiltTot=F(:,7);  % filtered total flow (m³/s) (= sum FiltBF, FiltOF, FiltIF)
       TiniPOTQF=F(:,8);% timestep that quickflow peak begins
       TPOTQF=F(:,9);   % timesteps that quickflow peak takes place
       TiniPOTSF=F(:,10);%timestep that slowflow peak begins
       TPOTSF=F(:,11);  % timesteps that slowflow peak takes place
       
       A=[Date2 ObsTot];
       A=A(~isnan(A(:,2)),:);% Remove lines with missing value for observed flow 
       a=size(Date2);
       nQobs=a(1,1);
       StartDateQobs=min(A(:,1));
       StopDateQobs=max(A(:,1));
       
       clear A a F name file
            
%% ------------------------------------------------------------------------- 
% 2. RUN ALL SUBMODELS & EXTRACT OUTPUT                                                  %
%--------------------------------------------------------------------------

% Run AquaCrop-Hydro with all its submodels

[Q_MBF,Q_MIF,Q_MOF,Q_MTF,area,f,Wrmin,Wrmax,pbf,SoilPar,SimACOutput,CatchACOutput,Par]=AquaCropHydro(DatapathAC, DatapathInput,nAC);

% Save important variables       
       % Save time variables        
        Day=CatchACOutput(:,13);    % Day number
        Month=CatchACOutput(:,14);  % Month number
        Year=CatchACOutput(:,15);   % Year number
        DateCombi=[Day,Month,Year]; % Date 
        
       % Save catchment-scale results
        TrCatch=CatchACOutput(:,1);     % Crop transpiration (actual) (mm)
        TrxCatch=CatchACOutput(:,2);    % Potential (maximum) crop transpiration (mm)
        ECatch=CatchACOutput(:,3);      % Evaporation(actual)(mm)
        ExCatch=CatchACOutput(:,4);     % Potential (maximum) evaporation(mm)
        ETaCatch=CatchACOutput(:,5);    % Actual evapotranspiration(mm)
        ETxCatch=CatchACOutput(:,6);    % Potential (maximum) evapotranspiration(mm)
        ROCatch=CatchACOutput(:,7);     % Runoff(mm)
        DPCatch=CatchACOutput(:,8);     % Deep percolation(mm)
        CRCatch=CatchACOutput(:,9);     % Capilary rise(mm)
        BundWatCatch=CatchACOutput(:,10); % Water between bunds (mm)
        Wr2Catch=CatchACOutput(:,11);   % Soil water content in 2 m soil depth (mm)
        CCCatch=CatchACOutput(:,12);    % Canopy Cover (%)
                
       % Save original AquaCrop results per simulation unit 
        Tr=SimACOutput{1,1};    % Crop transpiration(actual)(mm)
        Trx=SimACOutput{1,2};   % Potential (maximum) crop transpiration(mm)
        E=SimACOutput{1,3};     % Evaporation (actual)(mm)
        Ex=SimACOutput{1,4};    % Potential (maximum) evaporation(mm)
        ETa=SimACOutput{1,5};   % Actual evapotranspiration(mm)
        ETx=SimACOutput{1,6};   % Potential (maximum) evapotranspiration(mm)
        RO=SimACOutput{1,7};    % Runoff(mm)
        DP=SimACOutput{1,8};    % Deep percolation(mm)
        CR=SimACOutput{1,9};    % Capilary rise(mm)
        BundWat=SimACOutput{1,10}; % Water between bunds (mm)
        Wr2=SimACOutput{1,11};  % Soil water content in 2 m soil depth (mm)
        CC=SimACOutput{1,12};   % Canopy Cover (%)
        B=SimACOutput{1,13};    % Dry aboveground biomass during growing season (ton/ha)
        Bfin=SimACOutput{1,14}; % Final dry aboveground biomass at maturity (ton/ha)
        Y=SimACOutput{1,15};    % Dry yield at maturity (ton/ha)

        
%% ------------------------------------------------------------------------
% 3. SOME DATA PROCESSING                                                %
%--------------------------------------------------------------------------      

%%%% 3.1 TIME PERIOD DEFINITION 
   
       % maximal timeperiod
        StartDate=min(StartDateAC,StartDateQobs);
        StopDate=max(StopDateAC,StopDateQobs);
        n=(StopDate-StartDate)+1;
        Time=1:n;
        Time=Time.';
       
       % define calibration / validation period
        StartDateCalib=StartDateQobs;
        StartTimeCalib=(StartDateQobs-StartDate)+1;
        StopTimeCalib=Par(1,8)-1;
        StopDateCalib=(StartDate+StopTimeCalib)-1;
        nCalib=(StopDateCalib-StartDateCalib)+1;

        StartDateValid=StopDateCalib+1;
        StartTimeValid=StopTimeCalib+1;
        StopDateValid=StopDateQobs;
        StopTimeValid=(StopDateQobs-StartDate)+1;
        nValid=(StopDateValid-StartDateValid)+1;
    
%%%% 3.2 YIELD DATA PROCESSING       

% Extract yield data for the main growing season (main crop).
% For some landunits (non-agricultural units) this yield data do not make sense

Ymain=Y(2:2:end,:); % Only even sim runs correspond to the main season, odd run numbers are after-season crops

%%%% 3.3 POT PROCESSING        

% Delete empty values in loaded POT matrix & define number of  POT values
          TPOTQF = TPOTQF(~isnan(TPOTQF));
          a=size(TPOTQF);
          nPOTQF=a(1,1);
    
          TPOTSF=TPOTSF(~isnan(TPOTSF));
          a=size(TPOTSF);
          nPOTSF=a(1,1);
    
% Calculate observed discharge for  POT 
            QPOTQF=zeros(nPOTQF,1); % predefining result matrix
            QPOTSF=zeros(nPOTSF,1);
            for k=1:nPOTQF;
                QPOTQF(k)=ObsTot(TPOTQF(k));
            end

            for k=1:nPOTSF;
                QPOTSF(k)=ObsTot(TPOTSF(k));
            end
            
 % Calculate modelled POT values & do box cox transformation
     
    lambda=0.25; % suggested lambda for Box cox transformation
    
    % Predefine matrixes
    ModPOTQF=zeros(nPOTQF,1);
    TModPOTQF=zeros(nPOTQF,1);
    resQF=zeros(nPOTQF,1);
    SEQF=zeros(nPOTQF,1);
    QPOTQF_BC=zeros(nPOTQF,1);
    ModPOTQF_BC=zeros(nPOTQF,1);
    resQF_BC=zeros(nPOTQF,1);
    SEQF_BC=zeros(nPOTQF,1);
    
    ModPOTSF=zeros(nPOTSF,1);
    TModPOTSF=zeros(nPOTSF,1);
    resSF=zeros(nPOTSF,1);
    SESF=zeros(nPOTSF,1);
    QPOTSF_BC=zeros(nPOTSF,1);
    ModPOTSF_BC=zeros(nPOTSF,1);
    resSF_BC=zeros(nPOTSF,1);
    SESF_BC=zeros(nPOTSF,1);
            
    % QF POT values 
     for k=1:nPOTQF-1
         [M,I]=max(Q_MTF(TiniPOTQF(k):TiniPOTQF(k+1)));% find maximum between peak k and begin of next peak (M is value of maximum, I is the index of the maximum)
         ModPOTQF (k)= M; 
         TModPOTQF(k)= TiniPOTQF(k)+(I-1);
         resQF(k)=(ModPOTQF(k)-QPOTQF(k));%calculate residual
         SEQF(k)=(resQF(k))^2;%squared residual 

         QPOTQF_BC(k)=(QPOTQF(k)^lambda-1)/lambda;%box cox transformation
         ModPOTQF_BC(k)=(ModPOTQF(k)^lambda-1)/lambda;
         resQF_BC(k)=(ModPOTQF_BC(k)-QPOTQF_BC(k));%calculate residual
         SEQF_BC(k)=(resQF_BC(k))^2;
     end

     for k=nPOTQF 
         [M,I]=max(Q_MTF(TiniPOTQF(k):n));% last timestep you need to look untill the end of the timeseries
         ModPOTQF (k)=M;
         TModPOTQF(k)= TiniPOTQF(k)+(I-1);
         resQF(k)=(ModPOTQF(k)-QPOTQF(k));%squared residual 
         SEQF(k)=(resQF(k))^2;%squared residual 

         QPOTQF_BC(k)=(QPOTQF(k)^lambda-1)/lambda;
         ModPOTQF_BC(k)=(ModPOTQF(k)^lambda-1)/lambda;
         resQF_BC(k)=(ModPOTQF_BC(k)-QPOTQF_BC(k));%squared residual 
         SEQF_BC(k)=(resQF_BC(k))^2;
     end
 
  % SF POT values 
     for k=1:nPOTSF-1
         [M,I]=max(Q_MTF(TiniPOTSF(k):TiniPOTSF(k+1)));% find maximum between peak k and begin of next peak (M is value of maximum, I is the index of the maximum)
         ModPOTSF (k)=M;
         TModPOTSF(k)=TiniPOTSF(k)+(I-1);
         resSF(k)=(ModPOTSF(k)-QPOTSF(k));%calculate residual
         SESF(k)=(resSF(k))^2;%squared residual 

         QPOTSF_BC(k)=(QPOTSF(k)^lambda-1)/lambda;%box cox transformation
         ModPOTSF_BC(k)=(ModPOTSF(k)^lambda-1)/lambda;
         resSF_BC(k)=(ModPOTSF_BC(k)-QPOTSF_BC(k));%calculate residual
         SESF_BC(k)=(resSF_BC(k))^2;
     end

     for k=nPOTSF 
         [M,I]=max(Q_MTF(TiniPOTSF(k):n));
         ModPOTSF (k)=M;% last timestep you need to look untill the end of the timeseries
         TModPOTSF(k)=TiniPOTSF(k)+(I-1);
         resSF(k)=(ModPOTSF(k)-QPOTSF(k));%calculate residual
         SESF(k)=(resSF(k))^2;%squared residual

         QPOTSF_BC(k)=(QPOTSF(k)^lambda-1)/lambda;
         ModPOTSF_BC(k)=(ModPOTSF(k)^lambda-1)/lambda;
         resSF_BC(k)=(ModPOTSF_BC(k)-QPOTSF_BC(k));%calculate residual
         SESF_BC(k)=(resSF_BC(k))^2;
     end
  
 % Calculate summary statistics
    DEVQF_BC=mean(resQF_BC); % mean deviation "DEV" = average of all residuals
    DEVSF_BC=mean(resSF_BC);

    sigmaQF_BC=std(resQF_BC); %standard deviation of all residuals
    sigmaSF_BC=std(resSF_BC);
                         
%%%% 3.4 CUMULATIVE TIMESERIES
 
 % Remark: For this cumulative series we only start at the moment observations are available
 % and for dates that do not have observations we add 0 to the cumulative sum of both observed and simulated
 % values.
        
 indices = find(isnan(ObsTot) == 0);% alle NaN in the observed series
 start=min(indices);% The first observation (not NaN) is the timestep the calculation should start 
        
% Predefining cum matrixes
    ObsTotcum=zeros(n,1);
    FiltBFcum=zeros(n,1);
    FiltOFcum=zeros(n,1);
    FiltIFcum=zeros(n,1);
    FiltTotcum=zeros(n,1);
    
    Q_MBFcum=zeros(n,1);
    Q_MIFcum=zeros(n,1);
    Q_MOFcum=zeros(n,1);
    Q_MTFcum=zeros(n,1);
    
    Preccum=zeros(n,1);
    ETocum=zeros(n,1);
    ETaCatchcum=zeros(n,1);
    TrCatchcum=zeros(n,1);
    ECatchcum=zeros(n,1);
    DPCatchcum=zeros(n,1);
    ROCatchcum=zeros(n,1);
    Wr2Catchcum=zeros(n,1);
        
% Initial values               
    ObsTotcum(1)=ObsTot(1);
    FiltBFcum(1)=FiltBF(1);
    FiltOFcum(1)=FiltOF(1);
    FiltIFcum(1)=FiltIF(1);
    FiltTotcum(1)=FiltTot(1);
            
    Q_MBFcum(1)=Q_MBF(1);
    Q_MIFcum(1)=Q_MIF(1);
    Q_MOFcum(1)=Q_MOF(1);
    Q_MTFcum(1)=Q_MTF(1);
        
    Preccum(1)=Prec(1);
    ETaCatchcum(1)=ETaCatch(1);
    ETocum(1)=ETo(1);
    TrCatchcum(1)=TrCatch(1);
    ECatchcum(1)=ECatch(1);
    DPCatchcum(1)=DPCatch(1);
    ROCatchcum(1)=ROCatch(1);
    Wr2Catchcum(1)=Wr2Catch(1);
        
% Values for other timesteps
   Novergeslagen=0;% Counts skipped timesteps
   overgeslagen=[]; % Makes list of skipped timesteps
            
   for i=1:n
       if i<start % For this cumulative series we only start at the moment observations are available (before start cum value = 0)
          ObsTotcum(i)=0;
          FiltBFcum(i)=0;
          FiltOFcum(i)=0;
          FiltIFcum(i)=0;
          FiltTotcum(i)=0;

          Q_MBFcum(i)=0;
          Q_MIFcum(i)=0;
          Q_MOFcum(i)=0;
          Q_MTFcum(i)=0;
                
          Preccum(i)=0;
          ETaCatchcum(i)=0;
          ETocum(i)=0;
          TrCatchcum(i)=0;
          ECatchcum(i)=0;
          DPCatchcum(i)=0;
          ROCatchcum(i)=0;
          Wr2Catchcum(i)=0;
       else    
           if isnan(ObsTot(i))==1 % We do not want to include days for which Qtot obs are missing (in stead we want to sum up o then)      
            Novergeslagen=Novergeslagen+1;
            overgeslagen=[overgeslagen;i]; 

            ObsTotcum(i)=ObsTotcum(i-1)+0;
            FiltBFcum(i)=FiltBFcum(i-1)+0;
            FiltOFcum(i)=FiltOFcum(i-1)+0;
            FiltIFcum(i)=FiltIFcum(i-1)+0;
            FiltTotcum(i)=FiltTotcum(i-1)+0;

            Q_MBFcum(i)=Q_MBFcum(i-1)+0;
            Q_MIFcum(i)=Q_MIFcum(i-1)+0;
            Q_MOFcum(i)=Q_MOFcum(i-1)+0;
            Q_MTFcum(i)=Q_MTFcum(i-1)+0;

            Preccum(i)=Preccum(i-1)+0;
            ETaCatchcum(i)=ETaCatchcum(i-1)+0;
            ETocum(i)=ETocum(i-1)+0;
            TrCatchcum(i)=TrCatchcum(i-1)+0;
            ECatchcum(i)=ECatchcum(i-1)+0;
            DPCatchcum(i)=DPCatchcum(i-1)+0;
            ROCatchcum(i)=ROCatchcum(i-1)+0;
            Wr2Catchcum(i)=Wr2Catchcum(i-1)+0;
           else    
            ObsTotcum(i)=ObsTotcum(i-1)+ObsTot(i);
            FiltBFcum(i)=FiltBFcum(i-1)+FiltBF(i);
            FiltOFcum(i)=FiltOFcum(i-1)+FiltOF(i);
            FiltIFcum(i)=FiltIFcum(i-1)+FiltIF(i);
            FiltTotcum(i)=FiltTotcum(i-1)+FiltTot(i);

            Q_MBFcum(i)=Q_MBFcum(i-1)+Q_MBF(i);
            Q_MIFcum(i)=Q_MIFcum(i-1)+Q_MIF(i);
            Q_MOFcum(i)=Q_MOFcum(i-1)+Q_MOF(i);
            Q_MTFcum(i)=Q_MTFcum(i-1)+Q_MTF(i);

            Preccum(i)=Preccum(i-1)+Prec(i);
            ETaCatchcum(i)=ETaCatchcum(i-1)+ETaCatch(i);
            ETocum(i)=ETocum(i-1)+ETo(i);
            TrCatchcum(i)=TrCatchcum(i-1)+TrCatch(i);
            ECatchcum(i)=ECatchcum(i-1)+ECatch(i);
            DPCatchcum(i)=DPCatchcum(i-1)+DPCatch(i);
            ROCatchcum(i)=ROCatchcum(i-1)+ROCatch(i);
            Wr2Catchcum(i)=Wr2Catchcum(i-1)+Wr2Catch(i);
          end
       end
   end     
        
%%%% 3.5 CLIPPED TIMESERIES

% we want to create a dataseries which only contains the dates on which all
% observations values are present. Obviously simulation values are present for all
% timesteps of simulation period, so they don't need to be checked for
% missing values
 
% Put all data together in matrix
    AquaCropMatrix=[Day,Month,Year,Prec,ETo,TrCatch,TrxCatch,ECatch,ExCatch,ETaCatch,ETxCatch,ROCatch,DPCatch,CRCatch,BundWatCatch,Wr2Catch,CCCatch];
    HydroMatrix=[ObsTot,FiltTot,FiltBF,FiltIF,FiltOF,Q_MTF,Q_MBF,Q_MIF,Q_MOF,Time,Date2,pbf];

    a=size(AquaCropMatrix);
    b=size(HydroMatrix);
    k=a(1,2)+b(1,2);

    AllMatrix=[AquaCropMatrix,HydroMatrix];        

% Remove lines with missing values
    AllMatrix_cl=AllMatrix; 
    for i=1:k
    AllMatrix_cl=AllMatrix_cl(~isnan(AllMatrix_cl(:,i)),:);% remove lines with missing value for first column
    end

% Create clipped matrixes for each variable seperately
    Day_cl=AllMatrix_cl(:,1);
    Month_cl=AllMatrix_cl(:,2);
    Year_cl=AllMatrix_cl(:,3);
    DateCombi_cl=[AllMatrix_cl(:,1),AllMatrix_cl(:,2),AllMatrix_cl(:,3)];
    Prec_cl=AllMatrix_cl(:,4);
    ETo_cl=AllMatrix_cl(:,5);
    TrCatch_cl=AllMatrix_cl(:,6);
    TrxCatch_cl=AllMatrix_cl(:,7);
    ECatch_cl=AllMatrix_cl(:,8);
    ExCatch_cl=AllMatrix_cl(:,9);
    ETaCatch_cl=AllMatrix_cl(:,10);
    ETxCatch_cl=AllMatrix_cl(:,11);
    ROCatch_cl=AllMatrix_cl(:,12);
    DPCatch_cl=AllMatrix_cl(:,13);
    CRCatch_cl=AllMatrix_cl(:,14);
    BundWatCatch_cl=AllMatrix_cl(:,15);
    Wr2Catch_cl=AllMatrix_cl(:,16);
    CCCatch_cl=AllMatrix_cl(:,17);
    ObsTot_cl=AllMatrix_cl(:,18);
    FiltTot_cl=AllMatrix_cl(:,19);
    FiltBF_cl=AllMatrix_cl(:,20);
    FiltIF_cl=AllMatrix_cl(:,21);
    FiltOF_cl=AllMatrix_cl(:,22);
    Q_MTF_cl=AllMatrix_cl(:,23);
    Q_MBF_cl=AllMatrix_cl(:,24);
    Q_MIF_cl=AllMatrix_cl(:,25);
    Q_MOF_cl=AllMatrix_cl(:,26);
    Time_cl=AllMatrix_cl(:,27);
    Date2_cl=AllMatrix_cl(:,28);
    pbf_cl=AllMatrix_cl(:,29);
    a=size(Time_cl);
    n_cl=a(1,1);
      
clear a b r k i AquaCropMatrix HydroMatrix AllMatrix ;      

%%%% 3.6 CLIPPED TIMESERIES FOR PERIODE BEFORE MISSING DATA

AllMatrix_clCompl=AllMatrix_cl;
cond1= AllMatrix_clCompl(:,27)<min(overgeslagen);
AllMatrix_clCompl(cond1,:)=[];
clear cond1

% Seperate variables for complete period 
    Day_clCompl=AllMatrix_clCompl(:,1);
    Month_clCompl=AllMatrix_clCompl(:,2);
    Year_clCompl=AllMatrix_clCompl(:,3);
    DateCombi_clCompl=[AllMatrix_clCompl(:,1),AllMatrix_clCompl(:,2),AllMatrix_clCompl(:,3)];
    Prec_clCompl=AllMatrix_clCompl(:,4);
    ETo_clCompl=AllMatrix_clCompl(:,5);
    TrCatch_clCompl=AllMatrix_clCompl(:,6);
    TrxCatch_clCompl=AllMatrix_clCompl(:,7);
    ECatch_clCompl=AllMatrix_clCompl(:,8);
    ExCatch_clCompl=AllMatrix_clCompl(:,9);
    ETaCatch_clCompl=AllMatrix_clCompl(:,10);
    ETxCatch_clCompl=AllMatrix_clCompl(:,11);
    ROCatch_clCompl=AllMatrix_clCompl(:,12);
    DPCatch_clCompl=AllMatrix_clCompl(:,13);
    CRCatch_clCompl=AllMatrix_clCompl(:,14);
    BundWatCatch_clCompl=AllMatrix_clCompl(:,15);
    Wr2Catch_clCompl=AllMatrix_clCompl(:,16);
    CCCatch_clCompl=AllMatrix_clCompl(:,17);
    ObsTot_clCompl=AllMatrix_clCompl(:,18);
    FiltTot_clCompl=AllMatrix_clCompl(:,19);
    FiltBF_clCompl=AllMatrix_clCompl(:,20);
    FiltIF_clCompl=AllMatrix_clCompl(:,21);
    FiltOF_clCompl=AllMatrix_clCompl(:,22);
    Q_MTF_clCompl=AllMatrix_clCompl(:,23);
    Q_MBF_clCompl=AllMatrix_clCompl(:,24);
    Q_MIF_clCompl=AllMatrix_clCompl(:,25);
    Q_MOF_clCompl=AllMatrix_clCompl(:,26);
    Time_clCompl=AllMatrix_clCompl(:,27);
    Date2_clCompl=AllMatrix_clCompl(:,28);
    pbf_clCompl=AllMatrix_clCompl(:,29);
    a=size(Time_clCompl);
    n_clCompl=a(1,1);
clear a

%%%% 3.7 CLIPPED TIMESERIES FOR CALIBRATION AND VALIDATION PERIOD
%%%% SEPERATELY

AllMatrix_clCalib=AllMatrix_cl; % remove rows that do not belong to calib or valid period
cond1= AllMatrix_clCalib(:,27)>=StartTimeValid;
AllMatrix_clCalib(cond1,:)=[];

AllMatrix_clValid=AllMatrix_cl;
cond1= AllMatrix_clValid(:,27)<StartTimeValid;
AllMatrix_clValid(cond1,:)=[];

% Seperate variables for calibration 
    Day_clCalib=AllMatrix_clCalib(:,1);
    Month_clCalib=AllMatrix_clCalib(:,2);
    Year_clCalib=AllMatrix_clCalib(:,3);
    DateCombi_clCalib=[AllMatrix_clCalib(:,1),AllMatrix_clCalib(:,2),AllMatrix_clCalib(:,3)];
    Prec_clCalib=AllMatrix_clCalib(:,4);
    ETo_clCalib=AllMatrix_clCalib(:,5);
    TrCatch_clCalib=AllMatrix_clCalib(:,6);
    TrxCatch_clCalib=AllMatrix_clCalib(:,7);
    ECatch_clCalib=AllMatrix_clCalib(:,8);
    ExCatch_clCalib=AllMatrix_clCalib(:,9);
    ETaCatch_clCalib=AllMatrix_clCalib(:,10);
    ETxCatch_clCalib=AllMatrix_clCalib(:,11);
    ROCatch_clCalib=AllMatrix_clCalib(:,12);
    DPCatch_clCalib=AllMatrix_clCalib(:,13);
    CRCatch_clCalib=AllMatrix_clCalib(:,14);
    BundWatCatch_clCalib=AllMatrix_clCalib(:,15);
    Wr2Catch_clCalib=AllMatrix_clCalib(:,16);
    CCCatch_clCalib=AllMatrix_clCalib(:,17);
    ObsTot_clCalib=AllMatrix_clCalib(:,18);
    FiltTot_clCalib=AllMatrix_clCalib(:,19);
    FiltBF_clCalib=AllMatrix_clCalib(:,20);
    FiltIF_clCalib=AllMatrix_clCalib(:,21);
    FiltOF_clCalib=AllMatrix_clCalib(:,22);
    Q_MTF_clCalib=AllMatrix_clCalib(:,23);
    Q_MBF_clCalib=AllMatrix_clCalib(:,24);
    Q_MIF_clCalib=AllMatrix_clCalib(:,25);
    Q_MOF_clCalib=AllMatrix_clCalib(:,26);
    Time_clCalib=AllMatrix_clCalib(:,27);
    Date2_clCalib=AllMatrix_clCalib(:,28);
    pbf_clCalib=AllMatrix_clCalib(:,29);
    a=size(Time_clCalib);
    n_clCalib=a(1,1);
clear a

% Seperate variables for validation
    Day_clValid=AllMatrix_clValid(:,1);
    Month_clValid=AllMatrix_clValid(:,2);
    Year_clValid=AllMatrix_clValid(:,3);
    DateCombi_clValid=[AllMatrix_clValid(:,1),AllMatrix_clValid(:,2),AllMatrix_clValid(:,3)];
    Prec_clValid=AllMatrix_clValid(:,4);
    ETo_clValid=AllMatrix_clValid(:,5);
    TrCatch_clValid=AllMatrix_clValid(:,6);
    TrxCatch_clValid=AllMatrix_clValid(:,7);
    ECatch_clValid=AllMatrix_clValid(:,8);
    ExCatch_clValid=AllMatrix_clValid(:,9);
    ETaCatch_clValid=AllMatrix_clValid(:,10);
    ETxCatch_clValid=AllMatrix_clValid(:,11);
    ROCatch_clValid=AllMatrix_clValid(:,12);
    DPCatch_clValid=AllMatrix_clValid(:,13);
    CRCatch_clValid=AllMatrix_clValid(:,14);
    BundWatCatch_clValid=AllMatrix_clValid(:,15);
    Wr2Catch_clValid=AllMatrix_clValid(:,16);
    CCCatch_clValid=AllMatrix_clValid(:,17);
    ObsTot_clValid=AllMatrix_clValid(:,18);
    FiltTot_clValid=AllMatrix_clValid(:,19);
    FiltBF_clValid=AllMatrix_clValid(:,20);
    FiltIF_clValid=AllMatrix_clValid(:,21);
    FiltOF_clValid=AllMatrix_clValid(:,22);
    Q_MTF_clValid=AllMatrix_clValid(:,23);
    Q_MBF_clValid=AllMatrix_clValid(:,24);
    Q_MIF_clValid=AllMatrix_clValid(:,25);
    Q_MOF_clValid=AllMatrix_clValid(:,26);
    Time_clValid=AllMatrix_clValid(:,27);
    Date2_clValid=AllMatrix_clValid(:,28);
    pbf_clValid=AllMatrix_clValid(:,29);
    a=size(Time_clValid);
    n_clValid=a(1,1);
clear a

%%%% 3.8 CLIPPED SERIES FILLED WITH NaN
% This series contains a NaN value for dates on which observations are missing. That way the clipped series can be plotted in a linegraph

    ObsTot_fill=NaN(n,1);
    FiltTot_fill=NaN(n,1);
    FiltBF_fill=NaN(n,1);
    FiltIF_fill=NaN(n,1);
    FiltOF_fill=NaN(n,1);

    for i=1:n_cl
        ObsTot_fill(Time_cl(i))=ObsTot_cl(i);
        FiltTot_fill(Time_cl(i))=FiltTot_cl(i);
        FiltBF_fill(Time_cl(i))=FiltBF_cl(i);
        FiltIF_fill(Time_cl(i))=FiltIF_cl(i);
        FiltOF_fill(Time_cl(i))=FiltOF_cl(i);

    end

clear a

%%%% 3.9 SUBTOTALS FOR AGGREGATED TIMESTEPS (in mm)

% Determine decades (10-daily period)
    Dec_cl=NaN(n_cl,1);

    for i=1:n_cl;
        if Day_cl(i,1)/10<=1
            Dec_cl(i,1)=1;
        elseif Day_cl(i,1)/10<=2
            Dec_cl(i,1)=2;
        else
            Dec_cl(i,1)=3;
        end
    end   

% Calculate subtotals
    Hydro_cl=[Day_cl, Month_cl, Year_cl,Prec_cl, ETo_cl,ObsTot_cl,FiltTot_cl,FiltBF_cl,FiltIF_cl,FiltOF_cl,Q_MTF_cl,Q_MBF_cl,Q_MIF_cl,Q_MOF_cl,Dec_cl];


 for i=4:14%loop trough all data columns
     % Yearly subtotals
         HydroCl_YearSub{i-3}=accumarray(Year_cl,Hydro_cl(:,i)*f/area,[],[],NaN);

     % Monthly subtotals
         HydroCl_MonthSub{i-3}=accumarray([Year_cl,Month_cl],Hydro_cl(:,i)*f/area,[],[],NaN);

     % 10-daily subtotals
         HydroCl_DecadeSub{i-3}=accumarray([Year_cl,Month_cl, Dec_cl],Hydro_cl(:,i)*f/area,[],[],NaN);
        
 end
 
 % Determine position of result data
     minYear=DateCombi_cl(1,3); 
     minMonth=DateCombi_cl(1,2);
     minDay=DateCombi_cl(1,1);
            if minDay/10<=1
                minDec=1;
            elseif minDay/10<=2
                 minDec=2;
            else
                 minDec=3;
            end
     maxYear=DateCombi_cl(end,3);
     nYear=(maxYear-minYear)+1;
 
 % Clean result matrix
       for i=1:11 % loop trough all variables and remove empty years
          HydroCl_YearSub{1,i}=HydroCl_YearSub{1,i}(minYear:maxYear,:);
          HydroCl_MonthSub{1,i}=HydroCl_MonthSub{1,i}(minYear:maxYear,:);
          HydroCl_DecadeSub{1,i}=HydroCl_DecadeSub{1,i}(minYear:maxYear,:,:);
       end
      
       for i=1:11 % loop trough all variables and add year numbers 
         HydroCl_YearSub{1,i}(:,2)= HydroCl_YearSub{1,i}(:,1);
         HydroCl_YearSub{1,i}(:,1)= minYear:maxYear; 
       end
  
  % Put data in one column for monthly totals
  for i=1:11 % loop trough all 11 variables
      
      j=1;% first year (maybe not all months complete)
      startindex=1;
      HydroCl_MonthSub2{1,i}(1:startindex+((12-minMonth)),1)=minYear;% Write year
      HydroCl_MonthSub2{1,i}(1:startindex+((12-minMonth)),2)=minMonth:12;% Write month
      HydroCl_MonthSub2{1,i}(1:startindex+((12-minMonth)),3)=HydroCl_MonthSub{1,i}(j,minMonth:12);% Read all months of first year and writes them away
      startindex=startindex+((12-minMonth)+1);
      
      for j=2:nYear % loop trough all other years
            HydroCl_MonthSub2{1,i}(startindex:startindex+11,1)=(minYear+j)-1;% Write year
            HydroCl_MonthSub2{1,i}(startindex:startindex+11,2)=1:12;% Write month
            HydroCl_MonthSub2{1,i}(startindex:startindex+11,3)=HydroCl_MonthSub{1,i}(j,1:12);% Read all 12 months of year j and writes them away
            startindex=startindex+12;
      end
      
  end
  HydroCl_MonthSub=HydroCl_MonthSub2;
  clear HydroCl_MonthSub2;
  

  % Put data in one column for 10-daily data
    for i=1:11 % loop trough the 11 variabelen
      
      j=1;% first year (maybe not all months complete)
      startindex=1;
      HydroCl_DecadeSub2{1,i}(1:startindex+((12-minMonth)),1,1:3)=minYear;% Write year
      HydroCl_DecadeSub2{1,i}(1:startindex+((12-minMonth)),2,1:3)=[minMonth:12;minMonth:12;minMonth:12].';% Write month
      HydroCl_DecadeSub2{1,i}(1:startindex+((12-minMonth)),3,1:3)=HydroCl_DecadeSub{1,i}(j,minMonth:12,1:3);% Read all months of first year and writes 3 corresponding decades
      startindex=startindex+((12-minMonth)+1);
      
      for j=2:nYear % loop trough all other years
            HydroCl_DecadeSub2{1,i}(startindex:startindex+11,1,1:3)=(minYear+j)-1;% Write year
            HydroCl_DecadeSub2{1,i}(startindex:startindex+11,2,1:3)=[1:12;1:12;1:12].';% Write month
            HydroCl_DecadeSub2{1,i}(startindex:startindex+11,3,1:3)=HydroCl_DecadeSub{1,i}(j,1:12,1:3);% Read all 12 months of year j and writes 3 corresponding decades
            startindex=startindex+12;
      end

      %loop trough all year-month combinations
      a=size(HydroCl_DecadeSub2{1,1});%Determine number of months
      nm=a(1,1);
          m=1;% first month (maybe not all three decades complete)
             startindex=1;
             HydroCl_DecadeSub3{1,i}(1:startindex+2,1)=HydroCl_DecadeSub2{1,i}(m,1);% Write year
             HydroCl_DecadeSub3{1,i}(1:startindex+2,2)=HydroCl_DecadeSub2{1,i}(m,2);% Write month n
             HydroCl_DecadeSub3{1,i}(1:startindex+2,3)=minDec:3;% Write decade numbers
             HydroCl_DecadeSub3{1,i}(1:startindex+2,4)=HydroCl_DecadeSub2{1,i}(m,3,1:3); % Write data of the decades
             startindex=startindex+((3-minDec)+1);
          for m=2:nm % for all other months            
             HydroCl_DecadeSub3{1,i}(startindex:startindex+2,1)=HydroCl_DecadeSub2{1,i}(m,1);% Write year
             HydroCl_DecadeSub3{1,i}(startindex:startindex+2,2)=HydroCl_DecadeSub2{1,i}(m,2);% Write month
             HydroCl_DecadeSub3{1,i}(startindex:startindex+2,3)=1:3;% Write decade numbers
             HydroCl_DecadeSub3{1,i}(startindex:startindex+2,4)=HydroCl_DecadeSub2{1,i}(m,3,1:3);% Write data of the 3 decades 
             startindex=startindex+3;% go to next three decades
          end      
    end
  HydroCl_DecadeSub=HydroCl_DecadeSub3;
 
   clear i k nYear minYear maxYear minMonth HydroCl_MonthSub2 Hydro_cl HydroCl_DecadeSub2 HydroCl_DecadeSub3 

%% ------------------------------------------------------------------------
% 4. PERFORMANCE ANALYSIS                                                 %
%--------------------------------------------------------------------------          
 
%%%% 4.1 WATER BALANCES AND VOLUMES  

 %%% 4.1.1 Water balance AquaCrop (soil water balance)
 % data before routing (so in mm and not m³/s)

    % CALIB+VALID PERIOD 
    % Only use data for which also observation are available         
        CumWabalmm=zeros(8,3); % 3 columns (input, output, difference), 8 rows (P, E, Tr, DP, RO, Soil water,Bunds water, Total)

        CumWabalmm(1,1)=sum(Prec_cl);
        CumWabalmm(6,1)=Wr2Catch_cl(1,1);
        CumWabalmm(6,2)=Wr2Catch_cl(end,1);
        CumWabalmm(7,1)=BundWatCatch_cl(1,1);
        CumWabalmm(7,2)=BundWatCatch_cl(end,1);
        CumWabalmm(2,2)=sum(ECatch_cl);
        CumWabalmm(3,2)=sum(TrCatch_cl);
        CumWabalmm(4,2)=sum(DPCatch_cl);
        CumWabalmm(5,2)=sum(ROCatch_cl);
        CumWabalmm(8,1)=sum(CumWabalmm(1:7,1));
        CumWabalmm(8,2)=sum(CumWabalmm(2:7,2));
        CumWabalmm(:,3)=NaN;
        CumWabalmm(8,3)=CumWabalmm(8,1)-CumWabalmm(8,2);

    % WHOLE SIM PERIODE 
    % Use all data, also those for which no ObsTot observations are available        
        CumWabalmm2=zeros(8,3);% 3 columns (input, output, difference), 8 rows (P, E, Tr, DP, RO, Soil water,Bunds water, Total)

        CumWabalmm2(1,1)=sum(Prec);
        CumWabalmm2(6,1)=Wr2Catch(1,1);
        CumWabalmm2(6,2)=Wr2Catch(end,1);
        CumWabalmm2(7,1)=BundWatCatch(1,1);
        CumWabalmm2(7,2)=BundWatCatch(end,1);
        CumWabalmm2(2,2)=sum(ECatch);
        CumWabalmm2(3,2)=sum(TrCatch);
        CumWabalmm2(4,2)=sum(DPCatch);
        CumWabalmm2(5,2)=sum(ROCatch);
        CumWabalmm2(8,1)=sum(CumWabalmm2(1:7,1));
        CumWabalmm2(8,2)=sum(CumWabalmm2(2:7,2));
        CumWabalmm2(:,3)=NaN;
        CumWabalmm2(8,3)=CumWabalmm2(8,1)-CumWabalmm2(8,2);

 %%% 4.1.2 Water balance after routing (catchment water balance)
 % data after routing are converted from m³/s to mm/period using the
 % catchment area
 
    % CALIBRATION +VALIDATION PERIOD
    % Only use data for which also observation are available     
        CumWabalflow=NaN(9,7);% 7 columns (Observed, Observed filtered, Modelled, Obs-Mod,Obs-Mod error %, ObsFilt-Mod,ObsFilt-Mod error%)
                              % 9 rows (n days,P,ETa,TF,OF,IF+BF=DP,BF,IF,Total of sublfows)  

        CumWabalflow(1,1)=sum(Prec_cl); 
        CumWabalflow(1,3)=sum(Prec_cl);
        CumWabalflow(2,3)=sum(ETaCatch_cl);

        CumWabalflow(3,1)=sum(f*ObsTot_cl)/(area); 
        CumWabalflow(3,2)=sum(f*FiltTot_cl)/(area);
        CumWabalflow(3,3)=sum(f*Q_MTF_cl)/(area);

        CumWabalflow(4,2)=sum(f*FiltOF_cl)/(area);
        CumWabalflow(6,2)=sum(f*FiltBF_cl)/(area);
        CumWabalflow(7,2)=sum(f*FiltIF_cl)/(area);

        CumWabalflow(4,3)=sum(f*Q_MOF_cl)/(area);
        CumWabalflow(6,3)=sum(f*Q_MBF_cl)/(area);
        CumWabalflow(7,3)=sum(f*Q_MIF_cl)/(area);

        CumWabalflow(5,:)=sum(CumWabalflow(6:7,:));
        CumWabalflow(8,:)=sum(CumWabalflow(4:5,:));

        a=size(ObsTot_cl);
        nn=a(1,1);
        CumWabalflow(9,1:3)=nn;
        clear a nn 

        CumWabalflow(:,4)=CumWabalflow(:,1)-CumWabalflow(:,3);% absolute error
            for i=1:8
                if isnan(CumWabalflow(i,4))==0;
                CumWabalflow(i,5)=CumWabalflow(i,4)*100/CumWabalflow(i,1);% percent error
                end
            end
        CumWabalflow(:,6)=CumWabalflow(:,2)-CumWabalflow(:,3); % absolute error
            for i=1:8
                if isnan(CumWabalflow(i,6))==0;
                CumWabalflow(i,7)=CumWabalflow(i,6)*100/CumWabalflow(i,2);% percent error
                end
            end

    % CALIBRATION+VALIDATION PERIOD (ALL DATA)
    % Use all data, also those for which no ObsTot observations are available  
        CumWabalflow2=NaN(8,1);% 8 rows (P,ETa,E,Tr TF,OF,BF,IF)  

        CumWabalflow2(1,1)=sum(Prec);
        CumWabalflow2(2,1)=sum(ETaCatch);
        CumWabalflow2(3,1)=sum(ECatch);
        CumWabalflow2(4,1)=sum(TrCatch);
        CumWabalflow2(5,1)=sum(f*Q_MTF)/(area);
        CumWabalflow2(6,1)=sum(f*Q_MOF)/(area);
        CumWabalflow2(7,1)=sum(f*Q_MBF)/(area);
        CumWabalflow2(8,1)=sum(f*Q_MIF)/(area);

    % COMPLETE PERIOD 
    % period before the first missing observation
        CumWabalflowCompl=NaN(9,7);% 7 columns (Observed, Observed filtered, Modelled, Obs-Mod,Obs-Mod error %, ObsFilt-Mod,ObsFilt-Mod error%)
                                   % 9 rows (n days,P,ETa,TF,OF,IF+BF=DP,BF,IF,Total of sublfows) 

        CumWabalflowCompl(1,1)=sum(Prec_clCompl); 
        CumWabalflowCompl(1,3)=sum(Prec_clCompl);
        CumWabalflowCompl(2,3)=sum(ETaCatch_clCompl);

        CumWabalflowCompl(3,1)=sum(f*ObsTot_clCompl)/(area); 
        CumWabalflowCompl(3,2)=sum(f*FiltTot_clCompl)/(area);
        CumWabalflowCompl(3,3)=sum(f*Q_MTF_clCompl)/(area);

        CumWabalflowCompl(4,2)=sum(f*FiltOF_clCompl)/(area);
        CumWabalflowCompl(6,2)=sum(f*FiltBF_clCompl)/(area);
        CumWabalflowCompl(7,2)=sum(f*FiltIF_clCompl)/(area);

        CumWabalflowCompl(4,3)=sum(f*Q_MOF_clCompl)/(area);
        CumWabalflowCompl(6,3)=sum(f*Q_MBF_clCompl)/(area);
        CumWabalflowCompl(7,3)=sum(f*Q_MIF_clCompl)/(area);

        CumWabalflowCompl(5,:)=sum(CumWabalflowCompl(6:7,:));
        CumWabalflowCompl(8,:)=sum(CumWabalflowCompl(4:5,:));

        a=size(ObsTot_clCompl);
        nn=a(1,1);
        CumWabalflowCompl(9,1:3)=nn;
        clear a nn 

        CumWabalflowCompl(:,4)=CumWabalflowCompl(:,1)-CumWabalflowCompl(:,3);% absolute error
            for i=1:8
                if isnan(CumWabalflowCompl(i,4))==0;
                CumWabalflowCompl(i,5)=CumWabalflowCompl(i,4)*100/CumWabalflowCompl(i,1);% percent error
                end
            end
        CumWabalflowCompl(:,6)=CumWabalflowCompl(:,2)-CumWabalflowCompl(:,3); % absolute error 
            for i=1:8
                if isnan(CumWabalflowCompl(i,6))==0;
                CumWabalflowCompl(i,7)=CumWabalflowCompl(i,6)*100/CumWabalflowCompl(i,2);% percent error
                end
            end

    % CALIBRATION PERIOD
    % Only use data for which also observation are available                    
        CumWabalflowCalib=NaN(9,7);% 7 columns (Observed, Observed filtered, Modelled, Obs-Mod,Obs-Mod error %, ObsFilt-Mod,ObsFilt-Mod error%)
                                   % 9 rows (n days,P,ETa,TF,OF,IF+BF=DP,BF,IF,Total of sublfows) 


        CumWabalflowCalib(1,1)=sum(Prec_clCalib); 
        CumWabalflowCalib(1,3)=sum(Prec_clCalib);
        CumWabalflowCalib(2,3)=sum(ETaCatch_clCalib);

        CumWabalflowCalib(3,1)=sum(f*ObsTot_clCalib)/(area); 
        CumWabalflowCalib(3,2)=sum(f*FiltTot_clCalib)/(area);
        CumWabalflowCalib(3,3)=sum(f*Q_MTF_clCalib)/(area);

        CumWabalflowCalib(4,2)=sum(f*FiltOF_clCalib)/(area);
        CumWabalflowCalib(6,2)=sum(f*FiltBF_clCalib)/(area);
        CumWabalflowCalib(7,2)=sum(f*FiltIF_clCalib)/(area);

        CumWabalflowCalib(4,3)=sum(f*Q_MOF_clCalib)/(area);
        CumWabalflowCalib(6,3)=sum(f*Q_MBF_clCalib)/(area);
        CumWabalflowCalib(7,3)=sum(f*Q_MIF_clCalib)/(area);

        CumWabalflowCalib(5,:)=sum(CumWabalflowCalib(6:7,:));
        CumWabalflowCalib(8,:)=sum(CumWabalflowCalib(4:5,:));

        a=size(ObsTot_clCalib);
        nn=a(1,1);
        CumWabalflowCalib(9,1:3)=nn;
        clear a nn 

        CumWabalflowCalib(:,4)=CumWabalflowCalib(:,1)-CumWabalflowCalib(:,3);% absolute error
            for i=1:8
                if isnan(CumWabalflowCalib(i,4))==0;
                CumWabalflowCalib(i,5)=CumWabalflowCalib(i,4)*100/CumWabalflowCalib(i,1);% percent error
                end
            end
        CumWabalflowCalib(:,6)=CumWabalflowCalib(:,2)-CumWabalflowCalib(:,3); % absolute error
            for i=1:8
                if isnan(CumWabalflowCalib(i,6))==0;
                CumWabalflowCalib(i,7)=CumWabalflowCalib(i,6)*100/CumWabalflowCalib(i,2);% percent error
                end
            end
            
    % VALIDATION PERIOD
    % Only use data for which also observation are available 
        CumWabalflowValid=NaN(9,7);% 7 columns (Observed, Observed filtered, Modelled, Obs-Mod,Obs-Mod error %, ObsFilt-Mod,ObsFilt-Mod error%)
                                   % 9 rows (n days,P,ETa,TF,OF,IF+BF=DP,BF,IF,Total of sublfows) 
        
                                   CumWabalflowValid(1,1)=sum(Prec_clValid); 
        CumWabalflowValid(1,3)=sum(Prec_clValid);
        CumWabalflowValid(2,3)=sum(ETaCatch_clValid);

        CumWabalflowValid(3,1)=sum(f*ObsTot_clValid)/(area); 
        CumWabalflowValid(3,2)=sum(f*FiltTot_clValid)/(area);
        CumWabalflowValid(3,3)=sum(f*Q_MTF_clValid)/(area);

        CumWabalflowValid(4,2)=sum(f*FiltOF_clValid)/(area);
        CumWabalflowValid(6,2)=sum(f*FiltBF_clValid)/(area);
        CumWabalflowValid(7,2)=sum(f*FiltIF_clValid)/(area);

        CumWabalflowValid(4,3)=sum(f*Q_MOF_clValid)/(area);
        CumWabalflowValid(6,3)=sum(f*Q_MBF_clValid)/(area);
        CumWabalflowValid(7,3)=sum(f*Q_MIF_clValid)/(area);

        CumWabalflowValid(5,:)=sum(CumWabalflowValid(6:7,:));
        CumWabalflowValid(8,:)=sum(CumWabalflowValid(4:5,:));

        a=size(ObsTot_clValid);
        nn=a(1,1);
        CumWabalflowValid(9,1:3)=nn;
        clear a nn 

        CumWabalflowValid(:,4)=CumWabalflowValid(:,1)-CumWabalflowValid(:,3);% absolute error
            for i=1:8
                if isnan(CumWabalflowValid(i,4))==0;
                CumWabalflowValid(i,5)=CumWabalflowValid(i,4)*100/CumWabalflowValid(i,1);% percent error
                end
            end
        CumWabalflowValid(:,6)=CumWabalflowValid(:,2)-CumWabalflowValid(:,3); % absolute error
            for i=1:8
                if isnan(CumWabalflowValid(i,6))==0;
                CumWabalflowValid(i,7)=CumWabalflowValid(i,6)*100/CumWabalflowValid(i,2);% percent error
                end
            end

%%%% 4.2 PERFORMANCE STATISTICS 

%%% 4.2.1 Yearly subtotals
    Sets = 1; % Stats only need to be calculated for 1 set (=calibration and validatio period together) 
    StartTimeSet2=0; % value does not matter, will be ignored if Sets <2

    StatMatrixAllYear=NaN(6,4);% 4 columns : {'n', NSE, CVRMSE, R2}
                               % 6 rows : {QTot', 'QTotFilt', 'QBF+QIF','QBF', 'QIF' ,'QOF'}

    [Result]=PerformanceStat(HydroCl_YearSub{1,3}(:,1),HydroCl_YearSub{1,3}(:,2),HydroCl_YearSub{1,8}(:,2),Sets,StartTimeSet2);% ObsTot versus QMTF   
    StatMatrixAllYear(1,:)=[Result(:,1), Result(:,3), Result(:,6), Result(:,2)];
    [Result]=PerformanceStat(HydroCl_YearSub{1,3}(:,1),HydroCl_YearSub{1,4}(:,2),HydroCl_YearSub{1,8}(:,2),Sets,StartTimeSet2);%FiltTot_cl versus Q_MTF_cl
    StatMatrixAllYear(2,:)=[Result(:,1), Result(:,3), Result(:,6), Result(:,2)];
    [Result]=PerformanceStat(HydroCl_YearSub{1,3}(:,1),HydroCl_YearSub{1,5}(:,2)+HydroCl_YearSub{1,6}(:,2),HydroCl_YearSub{1,9}(:,2)+HydroCl_YearSub{1,10}(:,2),Sets,StartTimeSet2);%FiltBF_cl+FiltIF_cl versus Q_MBF_cl+Q_MIF_cl
    StatMatrixAllYear(3,:)=[Result(:,1), Result(:,3), Result(:,6), Result(:,2)];
    [Result]=PerformanceStat(HydroCl_YearSub{1,3}(:,1),HydroCl_YearSub{1,5}(:,2),HydroCl_YearSub{1,9}(:,2),Sets,StartTimeSet2);%FiltBF_cl versus Q_MBF_cl
    StatMatrixAllYear(4,:)=[Result(:,1), Result(:,3), Result(:,6), Result(:,2)];
    [Result]=PerformanceStat(HydroCl_YearSub{1,3}(:,1),HydroCl_YearSub{1,6}(:,2),HydroCl_YearSub{1,10}(:,2),Sets,StartTimeSet2);%FiltIF_cl versus Q_MIF_cl
    StatMatrixAllYear(5,:)=[Result(:,1), Result(:,3), Result(:,6), Result(:,2)];
    [Result]=PerformanceStat(HydroCl_YearSub{1,3}(:,1),HydroCl_YearSub{1,7}(:,2),HydroCl_YearSub{1,11}(:,2),Sets,StartTimeSet2);%FiltOF_cl versus Q_MOF_cl
    StatMatrixAllYear(6,:)=[Result(:,1), Result(:,3), Result(:,6), Result(:,2)];

%%% 4.2.2 Monthly subtotals
    Sets = 1; % Stats only need to be calculated for 1 set (=calibration and validatio period together) 
    StartTimeSet2=0; % value does not matter, will be ignored if Sets <2

    StatMatrixAllMonth=NaN(6,4);% 4 columns : {'n', NSE, CVRMSE, R2}
                                % 6 rows : {QTot', 'QTotFilt', 'QBF+QIF','QBF', 'QIF' ,'QOF'}

    [Result]=PerformanceStat(HydroCl_MonthSub{1,3}(:,1),HydroCl_MonthSub{1,3}(:,3),HydroCl_MonthSub{1,8}(:,3),Sets,StartTimeSet2);%ObsTot versus QMTF  
    StatMatrixAllMonth(1,:)=[Result(:,1), Result(:,3), Result(:,6), Result(:,2)];
    [Result]=PerformanceStat(HydroCl_MonthSub{1,3}(:,1),HydroCl_MonthSub{1,4}(:,3),HydroCl_MonthSub{1,8}(:,3),Sets,StartTimeSet2);%FiltTot_cl versus Q_MTF_cl
    StatMatrixAllMonth(2,:)=[Result(:,1), Result(:,3), Result(:,6), Result(:,2)];
    [Result]=PerformanceStat(HydroCl_MonthSub{1,3}(:,1),HydroCl_MonthSub{1,5}(:,3)+HydroCl_MonthSub{1,6}(:,3),HydroCl_MonthSub{1,9}(:,3)+HydroCl_MonthSub{1,10}(:,3),Sets,StartTimeSet2);%FiltBF_cl+FiltIF_cl versus Q_MBF_cl+Q_MIF_cl
    StatMatrixAllMonth(3,:)=[Result(:,1), Result(:,3), Result(:,6), Result(:,2)];
    [Result]=PerformanceStat(HydroCl_MonthSub{1,3}(:,1),HydroCl_MonthSub{1,5}(:,3),HydroCl_MonthSub{1,9}(:,3),Sets,StartTimeSet2);%FiltBF_cl versus Q_MBF_cl
    StatMatrixAllMonth(4,:)=[Result(:,1), Result(:,3), Result(:,6), Result(:,2)];
    [Result]=PerformanceStat(HydroCl_MonthSub{1,3}(:,1),HydroCl_MonthSub{1,6}(:,3),HydroCl_MonthSub{1,10}(:,3),Sets,StartTimeSet2);%FiltIF_cl versus Q_MIF_cl
    StatMatrixAllMonth(5,:)=[Result(:,1), Result(:,3), Result(:,6), Result(:,2)];
    [Result]=PerformanceStat(HydroCl_MonthSub{1,3}(:,1),HydroCl_MonthSub{1,7}(:,3),HydroCl_MonthSub{1,11}(:,3),Sets,StartTimeSet2);%FiltOF_cl versus Q_MOF_cl
    StatMatrixAllMonth(6,:)=[Result(:,1), Result(:,3), Result(:,6), Result(:,2)];

%%% 4.2.3 10-daily subtotals
    Sets = 1; % Stats only need to be calculated for 1 set (=calibration and validatio period together) 
    StartTimeSet2=0; % value does not matter, will be ignored if Sets <2

    StatMatrixAllDecade=NaN(6,4);% 4 columns : {'n', NSE, CVRMSE, R2}
                                 % 6 rows : {QTot', 'QTotFilt', 'QBF+QIF','QBF', 'QIF' ,'QOF'}

    [Result]=PerformanceStat(HydroCl_DecadeSub{1,3}(:,1),HydroCl_DecadeSub{1,3}(:,4),HydroCl_DecadeSub{1,8}(:,4),Sets,StartTimeSet2);%ObsTot versus QMTF  
    StatMatrixAllDecade(1,:)=[Result(:,1), Result(:,3), Result(:,6), Result(:,2)];
    [Result]=PerformanceStat(HydroCl_DecadeSub{1,3}(:,1),HydroCl_DecadeSub{1,4}(:,4),HydroCl_DecadeSub{1,8}(:,4),Sets,StartTimeSet2);%FiltTot_cl versus Q_MTF_cl
    StatMatrixAllDecade(2,:)=[Result(:,1), Result(:,3), Result(:,6), Result(:,2)];
    [Result]=PerformanceStat(HydroCl_DecadeSub{1,3}(:,1),HydroCl_DecadeSub{1,5}(:,4)+HydroCl_DecadeSub{1,6}(:,4),HydroCl_DecadeSub{1,9}(:,3)+HydroCl_DecadeSub{1,10}(:,3),Sets,StartTimeSet2);%FiltBF_cl+FiltIF_cl versus Q_MBF_cl+Q_MIF_cl
    StatMatrixAllDecade(3,:)=[Result(:,1), Result(:,3), Result(:,6), Result(:,2)];
    [Result]=PerformanceStat(HydroCl_DecadeSub{1,3}(:,1),HydroCl_DecadeSub{1,5}(:,4),HydroCl_DecadeSub{1,9}(:,4),Sets,StartTimeSet2);%FiltBF_cl versus Q_MBF_cl
    StatMatrixAllDecade(4,:)=[Result(:,1), Result(:,3), Result(:,6), Result(:,2)];
    [Result]=PerformanceStat(HydroCl_DecadeSub{1,3}(:,1),HydroCl_DecadeSub{1,6}(:,4),HydroCl_DecadeSub{1,10}(:,4),Sets,StartTimeSet2);%FiltIF_cl versus Q_MIF_cl
    StatMatrixAllDecade(5,:)=[Result(:,1), Result(:,3), Result(:,6), Result(:,2)];
    [Result]=PerformanceStat(HydroCl_DecadeSub{1,3}(:,1),HydroCl_DecadeSub{1,7}(:,4),HydroCl_DecadeSub{1,11}(:,4),Sets,StartTimeSet2);%FiltOF_cl versus Q_MOF_cl
    StatMatrixAllDecade(6,:)=[Result(:,1), Result(:,3), Result(:,6), Result(:,2)];

%%% 4.2.4 Daily values
    StatMatrixAllDay=NaN(8,4);  % 4 columns : {'n', NSE, CVRMSE, R2}
                                % 8 rows : {QTot', 'QTotFilt', 'QBF+QIF','QBF', 'QIF' ,'QOF','QFPOT','SFPOT'}
    StatMatrixCalibDay=NaN(6,4);% 4 columns : {'n', NSE, CVRMSE, R2}
                                % 6 rows : {QTot', 'QTotFilt', 'QBF+QIF','QBF', 'QIF' ,'QOF')
    StatMatrixValidDay=NaN(6,4);% 4 columns : {'n', NSE, CVRMSE, R2}
                                % 6 rows : {QTot', 'QTotFilt', 'QBF+QIF','QBF', 'QIF' ,'QOF')
    
    Sets = 2; %stats need to be calculated for calibration and validation period seperately
    StartTimeSet2=StartTimeValid;
    [Result]=PerformanceStat(Time_cl,ObsTot_cl,Q_MTF_cl,Sets,StartTimeSet2);
    StatMatrixAllDay(1,:)=[Result(1,1), Result(1,3), Result(1,6), Result(1,2)];
    StatMatrixCalibDay(1,:)=[Result(2,1), Result(2,3), Result(2,6), Result(2,2)];
    StatMatrixValidDay(1,:)=[Result(3,1), Result(3,3), Result(3,6), Result(3,2)];
    [Result]=PerformanceStat(Time_cl,FiltTot_cl,Q_MTF_cl,Sets,StartTimeSet2);
    StatMatrixAllDay(2,:)=[Result(1,1), Result(1,3), Result(1,6), Result(1,2)];
    StatMatrixCalibDay(2,:)=[Result(2,1), Result(2,3), Result(2,6), Result(2,2)];
    StatMatrixValidDay(2,:)=[Result(3,1), Result(3,3), Result(3,6), Result(3,2)];
    [Result]=PerformanceStat(Time_cl,FiltBF_cl+FiltIF_cl,Q_MBF_cl+Q_MIF_cl,Sets,StartTimeSet2);
    StatMatrixAllDay(3,:)=[Result(1,1), Result(1,3), Result(1,6), Result(1,2)];
    StatMatrixCalibDay(3,:)=[Result(2,1), Result(2,3), Result(2,6), Result(2,2)];
    StatMatrixValidDay(3,:)=[Result(3,1), Result(3,3), Result(3,6), Result(3,2)];
    [Result]=PerformanceStat(Time_cl,FiltBF_cl,Q_MBF_cl,Sets,StartTimeSet2);
    StatMatrixAllDay(4,:)=[Result(1,1), Result(1,3), Result(1,6), Result(1,2)];
    StatMatrixCalibDay(4,:)=[Result(2,1), Result(2,3), Result(2,6), Result(2,2)];
    StatMatrixValidDay(4,:)=[Result(3,1), Result(3,3), Result(3,6), Result(3,2)];
    [Result]=PerformanceStat(Time_cl,FiltIF_cl,Q_MIF_cl,Sets,StartTimeSet2);
    StatMatrixAllDay(5,:)=[Result(1,1), Result(1,3), Result(1,6), Result(1,2)];
    StatMatrixCalibDay(5,:)=[Result(2,1), Result(2,3), Result(2,6), Result(2,2)];
    StatMatrixValidDay(5,:)=[Result(3,1), Result(3,3), Result(3,6), Result(3,2)];
    [Result]=PerformanceStat(Time_cl,FiltOF_cl,Q_MOF_cl,Sets,StartTimeSet2);
    StatMatrixAllDay(6,:)=[Result(1,1), Result(1,3), Result(1,6), Result(1,2)];
    StatMatrixCalibDay(6,:)=[Result(2,1), Result(2,3), Result(2,6), Result(2,2)];
    StatMatrixValidDay(6,:)=[Result(3,1), Result(3,3), Result(3,6), Result(3,2)];

    Sets = 1; % Stats only need to be calculated for 1 set (=calibration and validatio period together) 
    StartTimeSet2=0; % value does not matter, will be ignored if Sets <2
    [Result]=PerformanceStat(TModPOTQF,QPOTQF, ModPOTQF,Sets,StartTimeSet2);
    StatMatrixAllDay(7,:)=[Result(:,1), Result(:,3), Result(:,6), Result(:,2)];
    [Result]=PerformanceStat(TModPOTSF,QPOTSF, ModPOTSF,Sets,StartTimeSet2);
    StatMatrixAllDay(8,:)=[Result(:,1), Result(:,3), Result(:,6), Result(:,2)];

%% ------------------------------------------------------------------------
% 5. VIZUALIZATION OF MODEL INPUT & RESULTS                               %
%--------------------------------------------------------------------------       
 
%%%% 5.1. WETSPRO FILTERING 
    figure('name','WETSPRO filter and POT');
        sub(1)=subplot (3,1,1,'fontsize',7);
        plot(Time,ObsTot, Time,FiltBF, Time,(FiltBF+FiltIF), Time,(FiltTot));
        xlabel('Time (days)','fontsize',7);
        ylabel('Flow (m3/s)','fontsize',7);
        legend ('Total discharge','Filtered baseflow','Filtered baseflow+interflow', 'Filtered total flow','Location','northwest')
        title('\bf\fontsize{12} discharge and subflows')
        a=axis;%asks for axislimits of y axis
        maxy=a(1,4);% neemt maximum y limit
        clear a
        line([StartTimeCalib,StartTimeCalib],[0,maxy],'Color','k','LineStyle','--'); % lines for calibration and validation period af te schermen
        line([StartTimeValid,StartTimeValid],[0,maxy],'Color','k','LineStyle','--');
        text(StartTimeCalib+nCalib/2,maxy-(maxy/20),'Calibration');%text om aan te duiden wat welke periode is
        text(StartTimeValid+nValid/2,maxy-(maxy/20),'Validation');

        sub(2)=subplot (3,1,2,'fontsize',7);
        plot (Time, ObsTot);
        hold
        scatter(TPOTSF, QPOTSF,'red');
        xlabel('Time (days)','fontsize',7);
        ylabel('Flow(m³/s)','fontsize',7);
        legend('Observed total discharge','Slow flow POT','Location','northwest');
        title('\bf\fontsize{12} POT selection')
        a=axis;%asks for axislimits of y axis
        maxy=a(1,4);% neemt maximum y limit
        clear a
        line([StartTimeCalib,StartTimeCalib],[0,maxy],'Color','k','LineStyle','--'); % lines for calibration and validation period af te schermen
        line([StartTimeValid,StartTimeValid],[0,maxy],'Color','k','LineStyle','--');
        text(StartTimeCalib+nCalib/2,maxy-(maxy/20),'Calibration');%text om aan te duiden wat welke periode is
        text(StartTimeValid+nValid/2,maxy-(maxy/20),'Validation');

        sub(3)=subplot (3,1,3,'fontsize',7);
        plot (Time, ObsTot);
        hold
        scatter(TPOTQF, QPOTQF,'magenta');
        xlabel('Time (days)','fontsize',7);
        ylabel('Flow(m³/s)','fontsize',7);
        legend('Observed total discharge', 'Quick flow POT','Location','northwest');
        a=axis;%asks for axislimits of y axis
        maxy=a(1,4);% neemt maximum y limit
        clear a
        line([StartTimeCalib,StartTimeCalib],[0,maxy],'Color','k','LineStyle','--'); % lines for calibration and validation period af te schermen
        line([StartTimeValid,StartTimeValid],[0,maxy],'Color','k','LineStyle','--');
        text(StartTimeCalib+nCalib/2,maxy-(maxy/20),'Calibration');%text om aan te duiden wat welke periode is
        text(StartTimeValid+nValid/2,maxy-(maxy/20),'Validation');

        linkaxes(sub,'xy')

%%%% 5.2. AQUACROP CATCHMENT SOIL WATER BALANCE 
     
    figure('name','Simulated catchment soil water balance');
        sub(1)=subplot (4,1,1,'fontsize',7);
        plot(Time,Prec, Time,ETo);
        xlabel('Time','fontsize',7);
        ylabel('Precipitation or ET0 (mm)','fontsize',7);
        legend ('Precipitation','Reference evapotranspiration')
        title('\bf\fontsize{12} climate input')
        
        sub(2)=subplot (4,1,2, 'fontsize',7);
        P=plot(Time,Wr2Catch,Time,repmat(2*SoilPar(1,end),nAC,1), Time, repmat(2*SoilPar(2,3),nAC,1));% graph of soil water content
        xlabel('Time (days)','fontsize',7);
        ylabel('Soil water content in 2 m soil (mm)','fontsize',7);
        legend ('soil water content', 'Field capacity','Permanent wilting point')
        title('\bf\fontsize{12} simulated soil water content')
        NameArray = {'LineStyle'};
        ValueArray = {'-','--','-.'}';
        set(P,NameArray,ValueArray)
        NameArray = {'Color'};
        ValueArray = {'b','k','k'}';
        set(P,NameArray,ValueArray);
    
        sub(3)=subplot (4,1,3,'fontsize',7);
        plot(Time,-DPCatch, Time,ROCatch);
        xlabel('Time (days)','fontsize',7);
        ylabel('mm','fontsize',7);
        legend ('Deep percolation','Runoff')
        title('\bf\fontsize{12} simulated water losses')
    
        sub(4)=subplot (4,1,4,'fontsize',7);
        plot(Time,[TrCatch,ECatch,ETaCatch]);
        xlabel('Time (days)','fontsize',7);
        ylabel('mm','fontsize',7);
        legend ('Transpiration','Evaporation','Evapotranspiration')
        title('\bf\fontsize{12} simulated evapotranspiration')
        
        linkaxes(sub,'x')

%%%% 5.3. DEEP PERCOLATION SPLIT 
 
 % Create variables for "observed" split
    FiltDP_cl=FiltIF_cl+FiltBF_cl;
    pbfObs_cl=FiltBF_cl./FiltDP_cl;
    pbfObs_cl(isnan(pbfObs_cl(:,1)),:)=1;% replace Nan (because of 0/0) by 1 

    FiltDP_clCalib=FiltIF_clCalib+FiltBF_clCalib;
    pbfObs_clCalib=FiltBF_clCalib./FiltDP_clCalib;
    pbfObs_clCalib(isnan(pbfObs_clCalib(:,1)),:)=1;% replace Nan (because of 0/0) by 1 

    FiltDP_clValid=FiltIF_clValid+FiltBF_clValid;
    pbfObs_clValid=FiltBF_clValid./FiltDP_clValid;
    pbfObs_clValid(isnan(pbfObs_clValid(:,1)),:)=1;% replace Nan (because of 0/0) by 1 
 
 % split series in Wr> Wrmin en Wr< Wrmin
    XAll=[Wr2Catch_cl,pbf_cl, pbfObs_cl];
    XCalib=[Wr2Catch_clCalib,pbf_clCalib,pbfObs_clCalib];
    XValid=[Wr2Catch_clValid,pbf_clValid,pbfObs_clValid];

    RAll=XAll;% matrix for high Wr (Right side graph, Wr< Wrmin)
    LAll=XAll;% matrix for small Wr (Left side graph, Wr> Wrmin)
    condleft=XAll(:,1)<Wrmin;
    condright=XAll(:,1)>Wrmin;
    RAll(condleft,:)=[];
    LAll(condright,:)=[];

    RCalib=XCalib;
    LCalib=XCalib;
    condleft=XCalib(:,1)<Wrmin;
    condright=XCalib(:,1)>Wrmin;
    RCalib(condleft,:)=[];
    LCalib(condright,:)=[];

    RValid=XValid;
    LValid=XValid;
    condleft=XValid(:,1)<Wrmin;
    condright=XValid(:,1)>Wrmin;
    RValid(condleft,:)=[];
    LValid(condright,:)=[];    

    Wr2Catch_clR=RAll(:,1); 
    pbf_clR=RAll(:,2); 
    pbfObs_clR=RAll(:,3);

    Wr2Catch_clCalibR=RCalib(:,1); 
    pbf_clCalibR=RCalib(:,2); 
    pbfObs_clCalibR=RCalib(:,3); 

    Wr2Catch_clValidR=RValid(:,1); 
    pbf_clValidR=RValid(:,2); 
    pbfObs_clValidR=RValid(:,3); 

    Wr2Catch_clL=LAll(:,1); 
    pbf_clL=LAll(:,2); 
    pbfObs_clL=LAll(:,3);

    Wr2Catch_clCalibL=LCalib(:,1); 
    pbf_clCalibL=LCalib(:,2); 
    pbfObs_clCalibL=LCalib(:,3); 

    Wr2Catch_clValidL=LValid(:,1); 
    pbf_clValidL=LValid(:,2); 
    pbfObs_clValidL=LValid(:,3);  

    clear condleft condright
    clear XAll XValid XCalib
    clear LCalib LValid LAll
    clear RCalib RValid RAll

% select only days on which RO>0
     pbf_clDP=pbf_cl;
     pbfObs_clDP=pbfObs_cl;
     Time_clDP=Time_cl;
     cond=DPCatch_cl(:,1)<=0.0072;
     pbf_clDP(cond,:)=[];
     pbfObs_clDP(cond,:)=[];
     Time_clDP(cond,:)=[];
     
 %figure result split
 figure('name','Split of DP in IF and BF');
    sub(1)=subplot (3,1,1,'fontsize',7);
    P=plot(Time,Wr2Catch, Time, repmat(Wrmin,nAC,1), Time, repmat(Wrmax,nAC,1),Time,repmat(2*SoilPar(1,end),nAC,1));% graph of soil water content
    NameArray = {'LineStyle'};
    ValueArray = {'-','--','--',':'}';
    set(P,NameArray,ValueArray)
    NameArray = {'Color'};
    ValueArray = {'b','k','k','r'}';
    set(P,NameArray,ValueArray);
    xlabel('Time','fontsize',7);
    ylabel('Soil water content in 2 m depth (mm)','fontsize',7);
    title('\bf\fontsize{12} soil water content')
    
    sub(2)=subplot (3,1,2, 'fontsize',7);
    plot(Time,pbf, Time_cl,pbfObs_cl)
    xlabel('Time (days)','fontsize',7);
    ylabel('Pbf(-)','fontsize',7);
    legend ('Modeled','Observed')
    title('\bf\fontsize{12} fraction of DP (IF+BF) going to baseflow')
    a=axis;%asks for axislimits of y axis
    maxy=a(1,4);% neemt maximum y limit
    line([StartTimeCalib,StartTimeCalib],[0,maxy],'Color','k','LineStyle','--'); % lines for calibration and validation period af te schermen
    line([StartTimeValid,StartTimeValid],[0,maxy],'Color','k','LineStyle','--');
       
    sub(3)=subplot(3,1,3, 'fontsize',7);
    plot(Time,Q_MBF+Q_MIF,Time_cl,FiltBF_cl+FiltIF_cl)
    xlabel('Time (days)','fontsize',7);
    ylabel('BF+IF(m³/s)','fontsize',7);
    legend ('Modeled','Observed')
    title('\bf\fontsize{12} Baseflow+Interflow')
    a=axis;%asks for axislimits of y axis
    maxy=a(1,4);% neemt maximum y limit
    clear a 
    line([StartTimeCalib,StartTimeCalib],[0,maxy],'Color','k','LineStyle','--'); % lines for calibration and validation period af te schermen
    line([StartTimeValid,StartTimeValid],[0,maxy],'Color','k','LineStyle','--');
    
    linkaxes(sub,'x')    
     
%figure of calibration pbf function
figure('name','pBF SWC relatie');
    sub(1)=subplot(4,1,1, 'fontsize',7);
    plot(Wr2Catch,DPCatch,'+');    
    xlabel('Soil water content in 2 m depth(m)','fontsize',7);
    ylabel('Deep percolation (mm)','fontsize',7);
    line([Wrmin Wrmin],[0 20],'Color','k');
    line([Wrmax Wrmax],[0 20],'Color','k');
    title('\bf\fontsize{12} DP for each water content')
    
    sub(2)=subplot(4,1,2, 'fontsize',7);
    plot(Wr2Catch_clL,pbf_clL,'+',Wr2Catch_clL,pbfObs_clL,'o',Wr2Catch_clR,pbf_clR,'+',Wr2Catch_clR,pbfObs_clR,'o');    
    line([Wrmin Wrmin],[0 1],'Color','k');
    line([Wrmax Wrmax],[0 1],'Color','k');
    hold on
        fit=polyfit(Wr2Catch_clR,pbfObs_clR,1);% creer de beste fit functie (1e graad)
        fit2=polyval(fit,Wr2Catch_clR);% gebruik de coefficienten om waarden te berekenen voor alle Wr waarden
        fit3=polyval(fit,Wrmin); % gebruik de coefficienten om waarden voor Wrmin en Wrmax te berekenen
        fit4=polyval(fit,Wrmax);
        plot(Wr2Catch_clR,fit2); % plot de fit functie
        polyfit_str = ['y = ' num2str(fit(1,1)) ' *x + ' num2str(fit(1,2))]; % creer tekst voor op functie te zetten
        Wr_str=['pbf(Wrmin) = ' num2str(fit3) ' & pbf(Wrmax) = ' num2str(fit4)];
        text(Wrmin,0.2,polyfit_str); % zet text op de grafiek
        text(Wrmin,0.05,Wr_str);
     hold on
        fit=polyfit(Wr2Catch_clL,pbfObs_clL,1);
        fit2=polyval(fit,Wr2Catch_clL);
        fit3=polyval(fit,Wrmin);
        plot(Wr2Catch_clL,fit2);
        polyfit_str = ['y = ' num2str(fit(1,1)) ' *x + ' num2str(fit(1,2))];
        Wr_str=['pbf(Wrmin) = ' num2str(fit3)];
        text(600,0.2,polyfit_str);
        text(600,0.05,Wr_str);
    xlabel('Soil water content in 2 m depth(m)','fontsize',7);
    ylabel('Fraction of BF in DP(-)','fontsize',7);
    title('\bf\fontsize{12} pbf-Wr relation for whole sim period')
    
    sub(3)=subplot(4,1,3, 'fontsize',7);
    plot(Wr2Catch_clCalibL,pbf_clCalibL,'+',Wr2Catch_clCalibL,pbfObs_clCalibL,'o',Wr2Catch_clCalibR,pbf_clCalibR,'+',Wr2Catch_clCalibR,pbfObs_clCalibR,'o'); 
    line([Wrmin Wrmin],[0 1],'Color','k');
    line([Wrmax Wrmax],[0 1],'Color','k');
    hold on
        fit=polyfit(Wr2Catch_clCalibR,pbfObs_clCalibR,1);
        fit2=polyval(fit,Wr2Catch_clCalibR);
        fit3=polyval(fit,Wrmin); 
        fit4=polyval(fit,Wrmax);
        plot(Wr2Catch_clCalibR,fit2);
        polyfit_str = ['y = ' num2str(fit(1,1)) ' *x + ' num2str(fit(1,2))];
         Wr_str=['pbf(Wrmin) = ' num2str(fit3) ' & pbf(Wrmax) = ' num2str(fit4)];
        text(Wrmin,0.2,polyfit_str); % zet text op de grafiek
        text(Wrmin,0.05,Wr_str);
     hold on
        fit=polyfit(Wr2Catch_clCalibL,pbfObs_clCalibL,1);
        fit2=polyval(fit,Wr2Catch_clCalibL);
        fit3=polyval(fit,Wrmin);
        plot(Wr2Catch_clCalibL,fit2);
        polyfit_str = ['y = ' num2str(fit(1,1)) ' *x + ' num2str(fit(1,2))];
        Wr_str=['pbf(Wrmin) = ' num2str(fit3)];
        text(600,0.2,polyfit_str);
        text(600,0.05,Wr_str);
    xlabel('Soil water content in 2 m depth(m)','fontsize',7);
    ylabel('Fraction of BF in DP(-)','fontsize',7);
    title('\bf\fontsize{12} pbf-Wr relation for calibration period')
    
    sub(4)=subplot(4,1,4, 'fontsize',7);
    plot(Wr2Catch_clValidL,pbf_clValidL,'+',Wr2Catch_clValidL,pbfObs_clValidL,'o',Wr2Catch_clValidR,pbf_clValidR,'+',Wr2Catch_clValidR,pbfObs_clValidR,'o'); 
    line([Wrmin Wrmin],[0 1],'Color','k');
    line([Wrmax Wrmax],[0 1],'Color','k');
    hold on
        fit=polyfit(Wr2Catch_clValidR,pbfObs_clValidR,1);
        fit2=polyval(fit,Wr2Catch_clValidR);
        fit3=polyval(fit,Wrmin); 
        fit4=polyval(fit,Wrmax);
        plot(Wr2Catch_clValidR,fit2);
        polyfit_str = ['y = ' num2str(fit(1,1)) ' *x + ' num2str(fit(1,2))];
        Wr_str=['pbf(Wrmin) = ' num2str(fit3) ' & pbf(Wrmax) = ' num2str(fit4)];
        text(Wrmin,0.2,polyfit_str); % zet text op de grafiek
        text(Wrmin,0.05,Wr_str);
     hold on
        fit=polyfit(Wr2Catch_clValidL,pbfObs_clValidL,1);
        fit2=polyval(fit,Wr2Catch_clValidL);
        fit3=polyval(fit,Wrmin);
        plot(Wr2Catch_clValidL,fit2);
        polyfit_str = ['y = ' num2str(fit(1,1)) ' *x + ' num2str(fit(1,2))];
        Wr_str=['pbf(Wrmin) = ' num2str(fit3)];
        text(600,0.2,polyfit_str);
        text(600,0.05,Wr_str);
    xlabel('Soil water content in 2 m depth(m)','fontsize',7);
    ylabel('Fraction of BF in DP(-)','fontsize',7);
    title('\bf\fontsize{12} pbf-Wr relation for validation period')
    
    linkaxes(sub,'x');
    
%figure of matching pbf function    
figure('name','pBF fit between obs and simulated'); 
    scatter(pbfObs_clDP,pbf_clDP) 
    xlabel('Observed baseflow fraction','fontsize',7);
    ylabel('Simulated baseflow fraction','fontsize',7);
    hold on
        [fit,s]=polyfit(pbfObs_clDP,pbf_clDP,1);
        fit2=polyval(fit,pbfObs_clDP);
        R2=1 - s.normr^2 / norm(pbfObs_clDP-mean(pbfObs_clDP))^2;
        plot(pbfObs_clDP,fit2);
        tekst1 = ['y = ' num2str(fit(1,1)) ' *x + ' num2str(fit(1,2))];
        tekst2=['R2= ' num2str(R2)];
        text(0.5,0.5,tekst1); % zet text op de grafiek
        text(0.5,0.45,tekst2); % zet text op de grafiek
     
 
clear R Wr2CatchR  Wr2Catch_clCalibR Wr2Catch_clValidR pbfR pbf_clCalibR pbf_clValidR pbfObsR  pbfObs_clCalibR pbfObs_clValidR
clear L Wr2CatchL   Wr2Catch_clCalibL Wr2Catch_clValidL pbfL  pbf_clCalibL pbf_clValidL pbfObsL  pbfObs_clCalibL pbfObs_clValidL
clear FiltDP FiltDP_clCalib FiltDP_clValid
   

%%%% 5.4 TOTAL FLOW 
figure('name','Total flow');
    plot(Time,[Q_MTF,ObsTot_fill]);% graph of total flow
    xlabel('Time','fontsize',7);
    ylabel('Total flow (m³/s)','fontsize',7);
    legend ('Modeled', 'Observed','Location','northwest','Orientation','vertical')
    title('\bf\fontsize{12} Total flow')
    a=axis;%asks for axislimits of y axis
    maxy=a(1,4);% neemt maximum y limit
    clear a
    line([StartTimeCalib,StartTimeCalib],[0,maxy],'Color','k','LineStyle','--'); % lines for calibration and validation period af te schermen
    line([StartTimeValid,StartTimeValid],[0,maxy],'Color','k','LineStyle','--');
    text(StartTimeCalib+nCalib/2,maxy-(maxy/8),'Calibration');%text om aan te duiden wat welke periode is
    text(StartTimeValid+nValid/2,maxy-(maxy/8),'Validation');
    text(StartTimeCalib/2,maxy-(maxy/8),'Warm-up');
 
%%%% 5.5 SUBFLOWS  
 figure('name','Flows comparison');
    sub(1)=subplot (5,1,1,'fontsize',7);
    plot(Time,Q_MBF,Time,FiltBF_fill);
    xlabel('Time','fontsize',7);
    ylabel('baseflow (m³/s)','fontsize',7);
    legend ('Modeled','Filtered', 'Location','northwest','Orientation','vertical')
    title('\bf\fontsize{12} Baseflow')
    a=axis;%asks for axislimits of y axis
    maxy=a(1,4);% neemt maximum y limit
    line([StartTimeCalib,StartTimeCalib],[0,maxy],'Color','k','LineStyle','--'); % lines for calibration and validation period af te schermen
    line([StartTimeValid,StartTimeValid],[0,maxy],'Color','k','LineStyle','--');
    text(StartTimeCalib+nCalib/20,maxy-(maxy/8),'Calibration');%text om aan te duiden wat welke periode is
    text(StartTimeValid+nValid/20,maxy-(maxy/8),'Validation');
    
    sub(2)=subplot (5,1,2,'fontsize',7);
    plot(Time,[Q_MIF,FiltIF_fill]);% graph of interflow
    xlabel('Time','fontsize',7);
    ylabel('inter flow (m³/s)','fontsize',7);
    legend ('Modeled','Filtered','Location','northwest','Orientation','vertical');
    title('\bf\fontsize{12} Interflow');
    a=axis;%asks for axislimits of y axis
    maxy=a(1,4);% neemt maximum y limit
    line([StartTimeCalib,StartTimeCalib],[0,maxy],'Color','k','LineStyle','--'); % lines for calibration and validation period af te schermen
    line([StartTimeValid,StartTimeValid],[0,maxy],'Color','k','LineStyle','--');
    text(StartTimeCalib+nCalib/20,maxy-(maxy/8),'Calibration');%text om aan te duiden wat welke periode is
    text(StartTimeValid+nValid/20,maxy-(maxy/8),'Validation');
    
    sub(3)=subplot (5,1,3,'fontsize',7);
    plot(Time,[Q_MIF+Q_MBF,FiltIF_fill+FiltBF_fill]);% graph of interflow
    xlabel('Time','fontsize',7);
    ylabel('inter flow (m³/s)','fontsize',7);
    legend ('Modeled','Filtered','Location','northwest','Orientation','vertical');
    title('\bf\fontsize{12} Baseflow+Interflow=DP');
    a=axis;%asks for axislimits of y axis
    maxy=a(1,4);% neemt maximum y limit
    line([StartTimeCalib,StartTimeCalib],[0,maxy],'Color','k','LineStyle','--'); % lines for calibration and validation period af te schermen
    line([StartTimeValid,StartTimeValid],[0,maxy],'Color','k','LineStyle','--');
    text(StartTimeCalib+nCalib/20,maxy-(maxy/8),'Calibration');%text om aan te duiden wat welke periode is
    text(StartTimeValid+nValid/20,maxy-(maxy/8),'Validation');
       
    sub(4)=subplot (5,1,4,'fontsize',7);
    plot(Time,[Q_MOF,FiltOF_fill]);% graph of overland flow
    xlabel('Time','fontsize',7);
    ylabel('overland flow (m³/s)','fontsize',7);
    legend ('Modeled','Filtered','Location','northwest','Orientation','vertical');
    title('\bf\fontsize{12} Overlandflow');
    a=axis;%asks for axislimits of y axis
    maxy=a(1,4);% neemt maximum y limit
    line([StartTimeCalib,StartTimeCalib],[0,maxy],'Color','k','LineStyle','--'); % lines for calibration and validation period af te schermen
    line([StartTimeValid,StartTimeValid],[0,maxy],'Color','k','LineStyle','--');
    text(StartTimeCalib+nCalib/20,maxy-(maxy/8),'Calibration');%text om aan te duiden wat welke periode is
    text(StartTimeValid+nValid/20,maxy-(maxy/8),'Validation');
    
    sub(5)=subplot (5,1,5,'fontsize',7);
    plot(Time,[Q_MTF,FiltTot_fill, ObsTot_fill]);% graph of total flow
    xlabel('Time','fontsize',7);
    ylabel('Total flow (m³/s)','fontsize',7);
    legend ('Modeled','Filtered', 'Observed','Location','northwest','Orientation','vertical')
    title('\bf\fontsize{12} Total flow')
    a=axis;%asks for axislimits of y axis
    maxy=a(1,4);% neemt maximum y limit
    line([StartTimeCalib,StartTimeCalib],[0,maxy],'Color','k','LineStyle','--'); % lines for calibration and validation period af te schermen
    line([StartTimeValid,StartTimeValid],[0,maxy],'Color','k','LineStyle','--');
    text(StartTimeCalib+nCalib/20,maxy-(maxy/8),'Calibration');%text om aan te duiden wat welke periode is
    text(StartTimeValid+nValid/20,maxy-(maxy/8),'Validation');
    clear a
    
    linkaxes(sub,'x')% link x axis of different plots (so that they change simultaneously

 %%%% 5.6 SUBFLOW PROBLEMS 
  figure('name','Baseflow problem');
  
        sub(1)=subplot (4,1,1,'fontsize',7);
        plot(Time,Prec, Time,ETo);
        xlabel('Time','fontsize',7);
        ylabel('Precipitation or ET0 (mm)','fontsize',7);
        legend ('Precipitation','Reference evapotranspiration')
        title('\bf\fontsize{12} climate input')
        
        sub(2)=subplot (4,1,2,'fontsize',7);
        P=plot(Time,Wr2Catch, Time, repmat(Wrmin,nAC,1), Time, repmat(Wrmax,nAC,1),Time,repmat(2*SoilPar(1,end),nAC,1));% graph of soil water content
        NameArray = {'LineStyle'};
        ValueArray = {'-','--','--',':'}';
        set(P,NameArray,ValueArray)
        NameArray = {'Color'};
        ValueArray = {'b','k','k','r'}';
        set(P,NameArray,ValueArray);
        xlabel('Time','fontsize',7);
        ylabel('Soil water content in 2 m depth (mm)','fontsize',7);
        title('\bf\fontsize{12} soil water content')

        sub(3)=subplot (4,1,3, 'fontsize',7);
        plot(Time,DPCatch)
        xlabel('Time','fontsize',7);
        ylabel('Deep percolation (mm)','fontsize',7);
        title('\bf\fontsize{12} deep percolation')
        
        sub(4)=subplot(4,1,4, 'fontsize',7);
        plot(Time,Q_MBF+Q_MIF,Time_cl,FiltBF_cl+FiltIF_cl)
        xlabel('Time (days)','fontsize',7);
        ylabel('BF+IF(m³/s)','fontsize',7);
        legend ('Modeled','Observed')
        title('\bf\fontsize{12} Baseflow+Interflow')
        a=axis;%asks for axislimits of y axis
        maxy=a(1,4);% neemt maximum y limit
        clear a 
        line([StartTimeCalib,StartTimeCalib],[0,maxy],'Color','k','LineStyle','--'); % lines for calibration and validation period af te schermen
        line([StartTimeValid,StartTimeValid],[0,maxy],'Color','k','LineStyle','--');

        linkaxes(sub,'x')    
        
   figure ('name', 'baseflow fit')
        plot(Q_MBF_cl,FiltBF_cl,'+',[0,0.15],[0,0.15],'-k');
        [fit,s]=polyfit(Q_MBF_cl,FiltBF_cl,2);
        fit2=polyval(fit,Q_MBF_cl);
        R2=1 - s.normr^2 / norm(Q_MBF_cl-mean(Q_MBF_cl))^2;
        hold on;
        plot(Q_MBF_cl,fit2,'--b');
        tekst = ['y = ' num2str(fit(1,1)) ' *x + ' num2str(fit(1,2))];
        tekst2=['R2= ' num2str(R2)];
        text(0.05,0.1,tekst);
        text(0.05,0.08,tekst2);       
        
   figure('name','Overlandflow problem');
        sub(1)=subplot (3,1,1,'fontsize',7);
        plot(Time,Prec, Time,ETo);% graph of climate data
        xlabel('Time','fontsize',7);
        ylabel('Precipitation or ET0 (mm)','fontsize',7);
        legend ('Precipitation','Reference evapotranspiration')
        title('\bf\fontsize{12} climate input')
        a=axis;%asks for axislimits of y axis
        maxy=a(1,4);% neemt maximum y limit
        line([StartTimeCalib,StartTimeCalib],[0,maxy],'Color','k','LineStyle','--'); % lines for calibration and validation period af te schermen
        line([StartTimeValid,StartTimeValid],[0,maxy],'Color','k','LineStyle','--');

        sub(2)=subplot (3,1,2, 'fontsize',7);
        P=plot(Time,Wr2Catch,Time,repmat(2*SoilPar(1,end),nAC,1), Time, repmat(2*SoilPar(2,end),nAC,1),Time, repmat(2*(SoilPar(2,end)+0.5*(SoilPar(1,end)-SoilPar(2,end))),nAC,1));% graph of soil water content
        xlabel('Time (days)','fontsize',7);
        ylabel('Soil water content in 2 m soil (mm)','fontsize',7);
        legend ('soil water content', 'Field capacity=CNIII','Permanent wilting point=CNII','CNII')
        title('\bf\fontsize{12} simulated soil water content')
        NameArray = {'LineStyle'};
        ValueArray = {'-','--','-.',':'}';
        set(P,NameArray,ValueArray)
        NameArray = {'Color'};
        ValueArray = {'b','k','k','k'}';
        set(P,NameArray,ValueArray);
        a=axis;%asks for axislimits of y axis
        maxy=a(1,4);% neemt maximum y limit
        line([StartTimeCalib,StartTimeCalib],[0,maxy],'Color','k','LineStyle','--'); % lines for calibration and validation period af te schermen
        line([StartTimeValid,StartTimeValid],[0,maxy],'Color','k','LineStyle','--');

        sub(3)=subplot (3,1,3,'fontsize',7);
        plot(Time,[Q_MOF,FiltOF]);% graph of overland flow
        xlabel('Time','fontsize',7);
        ylabel('overland flow (m³/s)','fontsize',7);
        legend ('Modeled','Filtered','Location','northwest','Orientation','vertical');
        title('\bf\fontsize{12} Overlandflow');
        a=axis;%asks for axislimits of y axis
        maxy=a(1,4);% neemt maximum y limit
        line([StartTimeCalib,StartTimeCalib],[0,maxy],'Color','k','LineStyle','--'); % lines for calibration and validation period af te schermen
        line([StartTimeValid,StartTimeValid],[0,maxy],'Color','k','LineStyle','--');
        linkaxes(sub,'x')% link x axis of different plots (so that they change simultaneously)
        clear a    
        
%%%% 5.7 CUMULATIVE FLOW IN TIME 

    % processing to visualize missing days at bottom of graph
        M=0;
        Z=NaN(n,1); 
        Z2=NaN((n-StartTimeCalib)+1,1);
        for i=1:Novergeslagen
            Z(overgeslagen(i,1),1)=M;
            Z2((overgeslagen(i,1)-StartTimeCalib)+1,1)=M;
         end

     figure('name','Cumulative flows original');
        P=plot(Time,[Q_MBFcum+Q_MIFcum,FiltBFcum+FiltIFcum,Q_MBFcum,FiltBFcum, Q_MOFcum,FiltOFcum, Q_MIFcum,FiltIFcum, Q_MTFcum,FiltTotcum,ObsTotcum]);% graph of cumulatives
        xlabel('Time','fontsize',7);
        ylabel('cumulative volume (m³/s)','fontsize',10);
        legend ('BF+IF Modeled','BF+IF Filtered','BF Modeled','BF Filtered','OF Modeled','OF Filtered','IF Modeled','IF Filtered', 'TOT Modeled','TOT Filtered','TOT Observed','Excluded timesteps','Location','northwest','Orientation','vertical')
        title('\bf\fontsize{12} Cumulative flow')
        NameArray = {'LineStyle'};
        ValueArray = {'-','--','-','--','-','--','-','--','-','--','-.'}';
        set(P,NameArray,ValueArray)
        NameArray = {'Color'};
        ValueArray = {'b','b','m','m','g','g','c','c','k','k','k'}';
        set(P,NameArray,ValueArray);
        a=axis;%asks for axislimits of y axis
        maxy=a(1,4);% neemt maximum y limit
        clear a
        hold on
        scatter(Time, Z,5,'x','k')
        line([StartTimeCalib,StartTimeCalib],[0,maxy],'Color','k','LineStyle','--'); % lines for calibration and validation period af te schermen
        line([StartTimeValid,StartTimeValid],[0,maxy],'Color','k','LineStyle','--');

    %  stop=min(overgeslagen);% stoppen met cum maken voor 1 grafiek op moment van eerste overgeslagen getal    
    %  figure('name','Cumulative flows period before missing data');

    %      P=plot(Time(1:stop,:),[Q_MBFcum(1:stop,:)+Q_MIFcum(1:stop,:),FiltBFcum(1:stop,:)+FiltIFcum(1:stop,:),Q_MBFcum(1:stop,:),FiltBFcum(1:stop,:), Q_MOFcum(1:stop,:),FiltOFcum(1:stop,:), Q_MIFcum(1:stop,:),FiltIFcum(1:stop,:), Q_MTFcum(1:stop,:),FiltTotcum(1:stop,:),ObsTotcum(1:stop,:)], Time,Z);% graph of cumulatives
    %     xlabel('Time','fontsize',7);
    %     ylabel('cumulative volume (m³/s)','fontsize',10);
    %     legend ('BF+IF Modeled','BF+IF Filtered','BF Modeled','BF Filtered','OF Modeled','OF Filtered','IF Modeled','IF Filtered', 'TOT Modeled','TOT Filtered','TOT Observed','Excluded timesteps','Location','northwest','Orientation','vertical')
    %     title('\bf\fontsize{12} Cumulative flow')
    %     NameArray = {'LineStyle'};
    %     ValueArray = {'-','--','-','--','-','--','-','--','-','--','-.',':'}';
    %     set(P,NameArray,ValueArray)
    %     NameArray = {'Color'};
    %     ValueArray = {'b','b','m','m','g','g','c','c','k','k','k','r'}';
    %     set(P,NameArray,ValueArray);
    %     a=axis;%asks for axislimits of y axis
    %     maxy=a(1,4);% neemt maximum y limit
    %     clear a
    %     line([StartTimeCalib,StartTimeCalib],[0,maxy],'Color','k','LineStyle','--'); % lines for calibration and validation period af te schermen
    %     line([StartTimeValid,StartTimeValid],[0,maxy],'Color','k','LineStyle','--');
    %     clear stop

%%% 5.8 VOLUME COMPARISON BAR PLOT 

% Create variables for stacked plots   
    plotof=[sum(f*FiltOF_cl)/(area),sum(f*Q_MOF_cl)/(area),sum(f*FiltOF_clCompl)/(area),sum(f*Q_MOF_clCompl)/(area),sum(f*FiltOF_clCalib)/(area),sum(f*Q_MOF_clCalib)/(area),sum(f*FiltOF_clValid)/(area),sum(f*Q_MOF_clValid)/(area)].';
    plotbf=[sum(f*FiltBF_cl)/(area),sum(f*Q_MBF_cl)/(area),sum(f*FiltBF_clCompl)/(area),sum(f*Q_MBF_clCompl)/(area),sum(f*FiltBF_clCalib)/(area),sum(f*Q_MBF_clCalib)/(area),sum(f*FiltBF_clValid)/(area),sum(f*Q_MBF_clValid)/(area)].';
    plotif=[sum(f*FiltIF_cl)/(area),sum(f*Q_MIF_cl)/(area),sum(f*Q_MIF_clCompl)/(area),sum(f*FiltIF_clCompl)/(area),sum(f*FiltIF_clCalib)/(area),sum(f*Q_MIF_clCalib)/(area),sum(f*FiltIF_clValid)/(area),sum(f*Q_MIF_clValid)/(area)].';
  
figure('name','Subflow stacked plot'); 
    P=bar([plotbf plotif plotof], 0.5,'stack');
    NameArray = {'FaceColor'};
    ValueArray = {'b','c','g'}';
    set(P,NameArray,ValueArray);
    set(gca,'XTickLabel',{'Obs all', 'Modall', 'Obs first part','Mod first part', 'Obs calib', 'Mod Calib','Obs Valid', 'Mod Valid'}); % zorgt dat op de x-as text staat en geen nummertjes
    ylabel('Cumulative Flow (mm/n days)','fontsize',7);
    legend ('baseflow','inter flow','overland flow','Location','eastoutside','Orientation','vertical');
 
%%% 5.9 SUBTOTALS PER YEAR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%     a=size(HydroCl_YearSub{1,3}(:,1));% alle years have stacked group
%     NumGroupsPerAxis =  a(1,1);
%     NumStacksPerGroup = 2;%observed versus modelled
%     NumStackElements = 3; % baseflow, interflow, overlandflow
%     
%     groupLabels = HydroCl_YearSub{1,3}(:,1);% labels to use on tick marks for groups
% 
%             Baseflow=[HydroCl_YearSub{1,5}(:,2),HydroCl_YearSub{1,9}(:,2)];
%             Interflow=[HydroCl_YearSub{1,6}(:,2),HydroCl_YearSub{1,10}(:,2)];
%             Overlandflow=[HydroCl_YearSub{1,7}(:,2),HydroCl_YearSub{1,11}(:,2)];
%             All=NaN(NumGroupsPerAxis,NumStacksPerGroup,NumStackElements); 
%             All(:,:,1)=Baseflow;
%             All(:,:,2)=Interflow;
%             All(:,:,3)=Overlandflow;
%             plotBarStackGroups(All, groupLabels);
%     
%  
%     figure('name','Yearly subtotals');
%     subplot (5,1,1, 'fontsize',7);
%     P=bar(HydroCl_YearSub{1,3}(:,1),[HydroCl_YearSub{1,3}(:,2) HydroCl_YearSub{1,8}(:,2)]);
%     title('\bf\fontsize{12} Total flow');
%     NameArray = {'FaceColor'};
%     ValueArray = {'b','c'}';
%     set(P,NameArray,ValueArray);
%     xlabel('Year','fontsize',7);
%     ylabel('Yearly flow volume (mm/year)','fontsize',7);
%     legend ('Observed','Modeled','Location','west');   
%     
%     subplot (5,1,2, 'fontsize',7);
%     P=bar(HydroCl_YearSub{1,3}(:,1),[HydroCl_YearSub{1,5}(:,2)+HydroCl_YearSub{1,6}(:,2),HydroCl_YearSub{1,9}(:,2)+HydroCl_YearSub{1,10}(:,2)]);
%     title('\bf\fontsize{12} Baseflow+Interflow=DP');
%     NameArray = {'FaceColor'};
%     ValueArray = {'b','c'}';
%     set(P,NameArray,ValueArray);
%     xlabel('Year','fontsize',7);
%     ylabel('Yearly flow volume (mm/year)','fontsize',7);
%     legend ('Observed','Modeled','Location','west');   
%     
%     subplot (5,1,3, 'fontsize',7);
%     P=bar(HydroCl_YearSub{1,3}(:,1),[HydroCl_YearSub{1,5}(:,2),HydroCl_YearSub{1,9}(:,2)]);
%     title('\bf\fontsize{12} Baseflow');
%     NameArray = {'FaceColor'};
%     ValueArray = {'b','c'}';
%     set(P,NameArray,ValueArray);
%     xlabel('Year','fontsize',7);
%     ylabel('Yearly flow volume (mm/year)','fontsize',7);
%     legend ('Observed','Modeled','Location','west');   
%     
%     subplot (5,1,4, 'fontsize',7);
%     P=bar(HydroCl_YearSub{1,3}(:,1),[HydroCl_YearSub{1,6}(:,2),HydroCl_YearSub{1,10}(:,2)]);
%     title('\bf\fontsize{12} Interflow');
%     NameArray = {'FaceColor'};
%     ValueArray = {'b','c'}';
%     set(P,NameArray,ValueArray);
%     xlabel('Year','fontsize',7);
%     ylabel('Yearly flow volume (mm/year)','fontsize',7);
%     legend ('Observed','Modeled','Location','west');   
%     
%     subplot (5,1,5, 'fontsize',7);
%     P=bar(HydroCl_YearSub{1,3}(:,1),[HydroCl_YearSub{1,7}(:,2),HydroCl_YearSub{1,11}(:,2)]);
%     title('\bf\fontsize{12} Overland flow');
%     NameArray = {'FaceColor'};
%     ValueArray = {'b','c'}';
%     set(P,NameArray,ValueArray);
%     xlabel('Year','fontsize',7);
%     ylabel('Yearly flow volume (mm/year)','fontsize',7);
%     legend ('Observed','Modeled','Location','west');   
%   

% %%% 5.10 SUBTOTALS PER MONTH (missing days excluded)%%%%%%%%%%%%%%%%%%%
%    
%     ZZ=num2str(HydroCl_MonthSub{1,3}(:,2)); % xas labels maken
%     YY=num2str(HydroCl_MonthSub{1,3}(:,1));
%     XAS=strcat(ZZ,'/',YY);
% 
%   %total stacked figure 
%     
%     a=size(HydroCl_MonthSub{1,3}(:,3));% alle years have stacked group
%     NumGroupsPerAxis =  a(1,1);
%     NumStacksPerGroup = 2;%observed versus modelled
%     NumStackElements = 3; % baseflow, interflow, overlandflow
%     
%     groupLabels = XAS;% labels to use on tick marks for groups
% 
%             Baseflow=[HydroCl_MonthSub{1,5}(:,3),HydroCl_MonthSub{1,9}(:,3)];
%             Interflow=[HydroCl_MonthSub{1,6}(:,3),HydroCl_MonthSub{1,10}(:,3)];
%             Overlandflow=[HydroCl_MonthSub{1,7}(:,3),HydroCl_MonthSub{1,11}(:,3)];
%             All=NaN(NumGroupsPerAxis,NumStacksPerGroup,NumStackElements); 
%             All(:,:,1)=Baseflow;
%             All(:,:,2)=Interflow;
%             All(:,:,3)=Overlandflow;
%             figure('name','Monthly subtotals - stacked components');
%             plotBarStackGroups(All, groupLabels);
%             
%     %figuur met subplots      
%     figure('name','Monthly subtotals');
%     subplot (5,1,1, 'fontsize',7);
%     P=bar([HydroCl_MonthSub{1,3}(:,3) HydroCl_MonthSub{1,8}(:,3)]);
%     title('\bf\fontsize{12} Total flow');
%     NameArray = {'FaceColor'};
%     ValueArray = {'b','c'}';
%     set(P,NameArray,ValueArray);
%     set(gca,'XTickLabel',XAS); % zorgt dat op de x-as text staat en geen nummertjes
%     xlabel('Month','fontsize',7);
%     ylabel('Monthly flow volume (mm/Month)','fontsize',7);
%     legend ('Observed','Modeled','Location','west');   
%     
%     subplot (5,1,2, 'fontsize',7);
%     P=bar([HydroCl_MonthSub{1,5}(:,3)+HydroCl_MonthSub{1,6}(:,3),HydroCl_MonthSub{1,9}(:,3)+HydroCl_MonthSub{1,10}(:,3)]);
%     title('\bf\fontsize{12} Baseflow+Interflow=DP');
%     NameArray = {'FaceColor'};
%     ValueArray = {'b','c'}';
%     set(P,NameArray,ValueArray);
%     set(gca,'XTickLabel',XAS); % zorgt dat op de x-as text staat en geen nummertjes
%     xlabel('Month','fontsize',7);
%     ylabel('Monthly flow volume (mm/Month)','fontsize',7);
%     legend ('Observed','Modeled','Location','west');      
%     
%         
%     subplot (5,1,3, 'fontsize',7);
%     P=bar([HydroCl_MonthSub{1,5}(:,3),HydroCl_MonthSub{1,9}(:,3)]);
%     title('\bf\fontsize{12} Baseflow');
%     NameArray = {'FaceColor'};
%     ValueArray = {'b','c'}';
%     set(P,NameArray,ValueArray);
%     set(gca,'XTickLabel',XAS); % zorgt dat op de x-as text staat en geen nummertjes
%     xlabel('Month','fontsize',7);
%     ylabel('Monthly flow volume (mm/Month)','fontsize',7);
%     legend ('Observed','Modeled','Location','west');   
%     
%     subplot (5,1,4, 'fontsize',7);
%     P=bar([HydroCl_MonthSub{1,6}(:,3),HydroCl_MonthSub{1,10}(:,3)]);
%     title('\bf\fontsize{12} Interflow');
%     NameArray = {'FaceColor'};
%     ValueArray = {'b','c'}';
%     set(P,NameArray,ValueArray);
%     set(gca,'XTickLabel',XAS); % zorgt dat op de x-as text staat en geen nummertjes
%     xlabel('Month','fontsize',7);
%     ylabel('Monthly flow volume (mm/Month)','fontsize',7);
%     legend ('Observed','Modeled','Location','west');   
%     
%     subplot (5,1,5, 'fontsize',7);
%     P=bar([HydroCl_MonthSub{1,7}(:,3),HydroCl_MonthSub{1,11}(:,3)]);
%     title('\bf\fontsize{12} Overland flow');
%     NameArray = {'FaceColor'};
%     ValueArray = {'b','c'}';
%     set(P,NameArray,ValueArray);
%     set(gca,'XTickLabel',XAS); % zorgt dat op de x-as text staat en geen nummertjes
%     xlabel('Month','fontsize',7);
%     ylabel('Monthly flow volume (mm/Month)','fontsize',7);
%     legend ('Observed','Modeled','Location','west');   
%   
%     clear XX YY XAS
% 

% %% 5.11 PEAK COMPARISON 
%     
% % present graph of peak values in time
%     
%     figure('name','POT values');  
%     subplot (2,1,1,'fontsize',7);
%     P=plot (Time, [ObsTot,Q_MTF]);
%     NameArray = {'Color'};
%     ValueArray = {'b','r'}';
%     set(P,NameArray,ValueArray);
%     hold
%     scatter(TPOTSF,QPOTSF,'o','b')
%     scatter(TModPOTSF,ModPOTSF,'o','r');
%     xlabel('Time (days)','fontsize',7);
%     ylabel('Flow(m³/s)','fontsize',7);
%     legend('Observed total discharge','Modeled total discharge','Observed SF POT','Modeled SF POT');
%     title('\bf\fontsize{12} Slow flow peaks')
%     
%     subplot (2,1,2,'fontsize',7);
%     P=plot (Time, [ObsTot,Q_MTF]);
%     NameArray = {'Color'};
%     ValueArray = {'b','r'}';
%     set(P,NameArray,ValueArray);
%     hold
%     scatter(TPOTQF, QPOTQF,'o','b')
%     scatter(TModPOTQF,ModPOTQF,'o','r');
%     xlabel('Time (days)','fontsize',7);
%     ylabel('Flow(m³/s)','fontsize',7);
%     legend('Observed total discharge','Modeled total discharge','Observed QF POT','Modeled QF POT');  
%     title('\bf\fontsize{12} Quick flow peaks')
%   
%  %present graph of observed versus modelled peak values
%   
%   figure('name','Box Cox Peaks'); 
%  
%   rank(1)=-50; % random waarden om zeker ver genoeg rechte lijn te tekenen
%   rank(2)=50;
%   DEVLine=rank+DEVQF_BC;
%   SigLine=DEVLine+sigmaQF_BC;
%   SigLine2=DEVLine-sigmaQF_BC;
%   
%   subplot(1,2,1,'fontsize',7);%graph for quickflow
%   plot (QPOTQF_BC, ModPOTQF_BC,'+',rank, rank, '-k',rank, DEVLine,'-b', rank, SigLine,'--b',rank, SigLine2,'--b');
%   xlabel('BC of Observed QF peak (m³/s)','fontsize',7);
%   ylabel('BC of Modelled QF peak (m³/s)','fontsize',7);
%   title('\bf\fontsize{12} Quickflow events');
%   minx=min(QPOTQF_BC-sigmaQF_BC);
%   maxx=max(QPOTQF_BC+sigmaQF_BC);
%   miny=min(ModPOTQF_BC-sigmaQF_BC);
%   maxy=max(ModPOTQF_BC+sigmaQF_BC);
%   axis([minx maxx miny maxy]);
%   legend ('Peaks','1:1 line','mean of residuals','mean residual +- standard deviation','Location','southeast','Orientation','vertical'); 
% 
%   
%   rank(1)=-50; % random waarden om zeker ver genoeg rechte lijn te tekenen
%   rank(2)=50;
%   DEVLine=rank+DEVSF_BC;
%   SigLine=DEVLine+sigmaSF_BC;
%   SigLine2=DEVLine-sigmaSF_BC;
%   
%   subplot(1,2,2,'fontsize',7);%graph for quickflow
%   plot (QPOTQF_BC, ModPOTQF_BC,'+',rank, rank, '-k', rank, DEVLine,'-b', rank, SigLine,'--b',rank, SigLine2,'--b');
%   xlabel('BC of Observed SF peak (m³/s)','fontsize',7);
%   ylabel('BC of Modelled SF peak (m³/s)','fontsize',7);
%   title('\bf\fontsize{12} Slowflow events');
%   minx=min(QPOTQF_BC-sigmaSF_BC);
%   maxx=max(QPOTQF_BC+sigmaSF_BC);
%   miny=min(ModPOTQF_BC-sigmaSF_BC);
%   maxy=max(ModPOTQF_BC+sigmaSF_BC);
%   axis([minx maxx miny maxy]);
%   legend ('Peaks','1:1 line','mean of residuals','mean residual +- standard deviation','Location','southeast','Orientation','vertical'); 


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 6. WRITE SOME OUTPUT TO EXCEL                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

name='Performance analysis.xls';
filename = [DatapathOutput name];

xlswrite(filename,StatMatrixCalibDay,'Statistics','C16:F21');
xlswrite(filename,StatMatrixValidDay,'Statistics','C24:F29');
xlswrite(filename,StatMatrixAllDay,'Statistics','C6:F13');

xlswrite(filename,StatMatrixAllMonth,'Statistics','I6:L11');
xlswrite(filename,StatMatrixAllYear,'Statistics','O6:R11');
xlswrite(filename,StatMatrixAllDecade,'Statistics','U6:X11');


xlswrite(filename,CumWabalmm,'Balances','B6:D13');
xlswrite(filename,CumWabalmm2,'Balances','B28:D35');
xlswrite(filename,CumWabalflow,'Balances','J6:P14');
xlswrite(filename,CumWabalflow2,'Balances','T26:T35');
xlswrite(filename,CumWabalflowCompl,'Balances','J21:P29');
xlswrite(filename,CumWabalflowCalib,'Balances','J35:P43');
xlswrite(filename,CumWabalflowValid,'Balances','J49:P57');

name2='Yield data.xls';
filename2 = [DatapathOutput name2];
xlswrite(filename2, Ymain,'allYdata','B5:AF19');

clear name name2 filename filename2

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 7. PUBLICATION FIGURES                                                  %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
date=x2mdate(Date2);
date_cl=x2mdate(Date2_cl);
%%% FIGURE 3 - CUMULATIVE COMPARISON

%Black and white figure
figure('name','Cumulative flows publication');
        P=plot(date(StartTimeCalib):date(end),[Q_MBFcum(StartTimeCalib:end,1)*f/area,FiltBFcum(StartTimeCalib:end,1)*f/area, Q_MOFcum(StartTimeCalib:end,1)*f/area,FiltOFcum(StartTimeCalib:end,1)*f/area, Q_MIFcum(StartTimeCalib:end,1)*f/area,FiltIFcum(StartTimeCalib:end,1)*f/area, Q_MTFcum(StartTimeCalib:end,1)*f/area,FiltTotcum(StartTimeCalib:end,1)*f/area,ObsTotcum(StartTimeCalib:end,1)*f/area]);% graph of cumulatives
        xlabel('Time','fontsize',10);
        ylabel('Cumulative volume (mm)','fontsize',10);
        NameArray = {'LineStyle'};
        ValueArray = {'-','--','-','--','-','--','-','--',':'}';
        set(P,NameArray,ValueArray)
        NameArray = {'Color'};
        ValueArray = {'[0.75, 0.75,0.75]','[0.75, 0.75 ,0.75]','[0.3, 0.3 ,0.3]','[0.3, 0.3 ,0.3]','[0.6, 0.6 ,0.6]','[0.6, 0.6 ,0.6]','k','k','k'}';
        set(P,NameArray,ValueArray)
        a=axis;%asks for axislimits of y axis
        maxy=a(1,4);% takes maximum y limit
        clear a
        hold on
        scatter(date(StartTimeCalib):date(end), Z2,10,'x','k')
        line([date(StartTimeValid),date(StartTimeValid)],[0,maxy],'Color','k','LineStyle','--');
        text(date(1700),maxy-(maxy/20),'\bf Calibration');%text to show which period
        text(date(4200),maxy-(maxy/20),'\bf Validation');
        text(date(end)+60,3500,'TOT')
        text(date(end)+60,2100,'BF')
        text(date(end)+60,600,'IF')
        text(date(end)+60,950,'OF')
        axis([date(1),date(end),0,maxy]);
        datetick('x','mm/yyyy')
        set(gca,'box','off');%remove ticks in upper right side
        
%Color figure
figure('name','Cumulative flows publication');
        P=plot(date(StartTimeCalib):date(end),[Q_MBFcum(StartTimeCalib:end,1)*f/area,FiltBFcum(StartTimeCalib:end,1)*f/area, Q_MOFcum(StartTimeCalib:end,1)*f/area,FiltOFcum(StartTimeCalib:end,1)*f/area, Q_MIFcum(StartTimeCalib:end,1)*f/area,FiltIFcum(StartTimeCalib:end,1)*f/area, Q_MTFcum(StartTimeCalib:end,1)*f/area,FiltTotcum(StartTimeCalib:end,1)*f/area,ObsTotcum(StartTimeCalib:end,1)*f/area]);% graph of cumulatives
        xlabel('Time','fontsize',10);
        ylabel('Cumulative volume (mm)','fontsize',10);
        NameArray = {'LineStyle'};
        ValueArray = {'-','-.','-','-.','-','-.','-','-.','-'}';
        set(P,NameArray,ValueArray)
        NameArray = {'Color'};
        ValueArray = {'b','r','b','r','b','r','b','r','r'}';
        set(P,NameArray,ValueArray)
        a=axis;%asks for axislimits of y axis
        maxy=a(1,4);% takes maximum y limit
        clear a
        hold on
        scatter(date(StartTimeCalib):date(end), Z2,10,'x','k')
        line([date(StartTimeValid),date(StartTimeValid)],[0,maxy],'Color','k','LineStyle','--');
        text(date(1700),maxy-(maxy/20),'\bf Calibration');%text to show which period
        text(date(4200),maxy-(maxy/20),'\bf Validation');
        text(date(end)+60,3500,'TOT')
        text(date(end)+60,2100,'BF')
        text(date(end)+60,600,'IF')
        text(date(end)+60,950,'OF')
        axis([date(1),date(end),0,maxy]);
        datetick('x','mm/yyyy')
        set(gca,'box','off');%remove ticks in upper right side

%FIG file
figure('name','Cumulative flows publication');
        P=plot(date(StartTimeCalib):date(end),[Q_MBFcum(StartTimeCalib:end,1)*f/area,FiltBFcum(StartTimeCalib:end,1)*f/area, Q_MOFcum(StartTimeCalib:end,1)*f/area,FiltOFcum(StartTimeCalib:end,1)*f/area, Q_MIFcum(StartTimeCalib:end,1)*f/area,FiltIFcum(StartTimeCalib:end,1)*f/area, Q_MTFcum(StartTimeCalib:end,1)*f/area,FiltTotcum(StartTimeCalib:end,1)*f/area,ObsTotcum(StartTimeCalib:end,1)*f/area]);% graph of cumulatives
        hold on
        p2=plot(date(StartTimeCalib):date(end),[Q_MBFcum(StartTimeCalib:end,1)*f/area,FiltBFcum(StartTimeCalib:end,1)*f/area,ObsTotcum(StartTimeCalib:end,1)*f/area]);
        xlabel('Time','fontsize',10);
        ylabel('Cumulative volume (mm)','fontsize',10);
        NameArray = {'LineStyle'};
        ValueArray = {'-','-.','-','-.','-','-.','-','-.','-'}';
        set(P,NameArray,ValueArray)
        NameArray = {'Color'};
        ValueArray = {'b','r','b','r','b','r','b','r','r'}';
        set(P,NameArray,ValueArray)
        NameArray = {'LineStyle'};
        ValueArray = {'-','-.','-'}';
        set(p2,NameArray,ValueArray)
        NameArray = {'Color'};
        ValueArray = {'b','r','r'}';
        set(p2,NameArray,ValueArray)       
        a=axis;%asks for axislimits of y axis
        maxy=a(1,4);% takes maximum y limit
        clear a
        hold on
        scatter(date(StartTimeCalib):date(end), Z2,10,'x','k')
        line([date(StartTimeValid),date(StartTimeValid)],[0,maxy],'Color','k','LineStyle','--');
        text(date(1700),maxy-(maxy/20),'\bf Calibration');%text to show which period
        text(date(4200),maxy-(maxy/20),'\bf Validation');
        text(date(end)+60,3500,'TOT')
        text(date(end)+60,2100,'BF')
        text(date(end)+60,600,'IF')
        text(date(end)+60,950,'OF')
        axis([date(1),date(end),0,maxy]);
        datetick('x','mm/yyyy')
        set(gca,'box','off');%remove ticks in upper right side
        legend (p2,'Modelled', 'Observed filtered','Observed','Location','west');

%%% FIGURE 4 -  SUBFLOW COMPARISON

%Black and white figure
figure('name','Flow comparison');
        sub(1)=subplot (3,1,1:2,'fontsize',10);
        p=plot(date,[ObsTot_fill,Q_MTF]);% graph of total flow
        xlabel('Time','fontsize',10);
        ylabel('Total discharge (m³/s)','fontsize',10);
        axis tight;
        NameArray = {'Color'};
        ValueArray = {'k',[0.6 0.6 0.6]}';
        set(p,NameArray,ValueArray);
        line([date(StartTimeCalib),date(StartTimeCalib)],[0,1.35],'Color','k','LineStyle','--'); % lines for calibration and validation period af te schermen
        line([date(StartTimeValid),date(StartTimeValid)],[0,1.35],'Color','k','LineStyle','--');
        text(date(2400),1.2,'Calibration');%text om aan te duiden wat welke periode is
        text(date(4300),1.2,'Validation');
        text(date(400),1.2,'Warm-up');
        set(gca,'box','off');%remove ticks in upper right side
        tickDates = datenum(2000:1:2016,1,1) ; %// creates a vector of tick positions
        set(gca, 'XTick' , tickDates, 'XTickLabel' , datestr(tickDates,'yyyy') ) %//

        sub(2)=subplot (3,1,3);
        p=plot(date,[FiltBF_fill,Q_MBF]);
        xlabel('Time','fontsize',10);
        ylabel('Baseflow (m³/s)','fontsize',10);
        axis tight;
        NameArray = {'Color'};
        ValueArray = {'k',[0.6 0.6 0.6]}';
        set(p,NameArray,ValueArray);
        a=axis;%asks for axislimits of y axis
        maxy=a(1,4);% neemt maximum y limit
        line([date(StartTimeCalib),date(StartTimeCalib)],[0,maxy],'Color','k','LineStyle','--'); % lines for calibration and validation period af te schermen
        line([date(StartTimeValid),date(StartTimeValid)],[0,maxy],'Color','k','LineStyle','--');
        set(gca,'box','off');%remove ticks in upper right side
        tickDates = datenum(2000:1:2016,1,1) ; %// creates a vector of tick positions
        set(gca, 'XTick' , tickDates, 'XTickLabel' , datestr(tickDates,'yyyy') ) %//
        
        clear a p maxy
        linkaxes(sub,'x')% link x axis of different plots (so that they change simultaneously)
%Color figure
figure('name','Flow comparison');
        sub(1)=subplot (3,1,1:2,'fontsize',10);
        p=plot(date,[ObsTot_fill,Q_MTF]);% graph of total flow
        xlabel('Time','fontsize',10);
        ylabel('Total discharge (m³/s)','fontsize',10);
        axis tight;
        NameArray = {'Color'};
        ValueArray = {'r','b'}';
        set(p,NameArray,ValueArray);
        line([date(StartTimeCalib),date(StartTimeCalib)],[0,1.35],'Color','k','LineStyle','--'); % lines for calibration and validation period af te schermen
        line([date(StartTimeValid),date(StartTimeValid)],[0,1.35],'Color','k','LineStyle','--');
        text(date(2400),1.2,'Calibration');%text om aan te duiden wat welke periode is
        text(date(4300),1.2,'Validation');
        text(date(400),1.2,'Warm-up');
        set(gca,'box','off');%remove ticks in upper right side
        tickDates = datenum(2000:1:2016,1,1) ; %// creates a vector of tick positions
        set(gca, 'XTick' , tickDates, 'XTickLabel' , datestr(tickDates,'yyyy') ) %//

        sub(2)=subplot (3,1,3);
        p=plot(date,[FiltBF_fill,Q_MBF]);
        xlabel('Time','fontsize',10);
        ylabel('Baseflow (m³/s)','fontsize',10);
        axis tight;
        NameArray = {'Color'};
        ValueArray = {'r','b'}';
        set(p,NameArray,ValueArray);
        a=axis;%asks for axislimits of y axis
        maxy=a(1,4);% neemt maximum y limit
        line([date(StartTimeCalib),date(StartTimeCalib)],[0,0.13],'Color','k','LineStyle','--'); % lines for calibration and validation period af te schermen
        line([date(StartTimeValid),date(StartTimeValid)],[0,0.13],'Color','k','LineStyle','--');
        set(gca,'box','off');%remove ticks in upper right side
        tickDates = datenum(2000:1:2016,1,1) ; %// creates a vector of tick positions
        set(gca, 'XTick' , tickDates, 'XTickLabel' , datestr(tickDates,'yyyy') ) %//
        
        clear a p maxy
        linkaxes(sub,'x')% link x axis of different plots (so that they change simultaneously)

%FIG file
figure('name','Flow comparison');
        sub(1)=subplot (3,1,1:2,'fontsize',10);
        p=plot(date,[ObsTot_fill,Q_MTF]);% graph of total flow
        xlabel('Time','fontsize',10);
        ylabel('Total discharge (m³/s)','fontsize',10);
        axis tight;
        NameArray = {'Color'};
        ValueArray = {'r','b'}';
        set(p,NameArray,ValueArray);
        line([date(StartTimeCalib),date(StartTimeCalib)],[0,1.35],'Color','k','LineStyle','--'); % lines for calibration and validation period af te schermen
        line([date(StartTimeValid),date(StartTimeValid)],[0,1.35],'Color','k','LineStyle','--');
        text(date(2400),1.2,'Calibration');%text om aan te duiden wat welke periode is
        text(date(4300),1.2,'Validation');
        text(date(400),1.2,'Warm-up');
        set(gca,'box','off');%remove ticks in upper right side
        tickDates = datenum(2000:1:2016,1,1) ; %// creates a vector of tick positions
        set(gca, 'XTick' , tickDates, 'XTickLabel' , datestr(tickDates,'yyyy') ) %//
        legend ('Observed total discharge', 'Modelled total discharge', 'Orientation', 'vertical', 'Location', 'northeast');

        sub(2)=subplot (3,1,3);
        p=plot(date,[FiltBF_fill,Q_MBF]);
        xlabel('Time','fontsize',10);
        ylabel('Baseflow (m³/s)','fontsize',10);
        axis tight;
        NameArray = {'Color'};
        ValueArray = {'r','b'}';
        set(p,NameArray,ValueArray);
        a=axis;%asks for axislimits of y axis
        maxy=a(1,4);% neemt maximum y limit
        line([date(StartTimeCalib),date(StartTimeCalib)],[0,maxy],'Color','k','LineStyle','--'); % lines for calibration and validation period af te schermen
        line([date(StartTimeValid),date(StartTimeValid)],[0,maxy],'Color','k','LineStyle','--');
        set(gca,'box','off');%remove ticks in upper right side
        tickDates = datenum(2000:1:2016,1,1) ; %// creates a vector of tick positions
        set(gca, 'XTick' , tickDates, 'XTickLabel' , datestr(tickDates,'yyyy') ) %//
        legend ('Observed baseflow', 'Modelled baseflow','Orientation', 'vertical', 'Location', 'northeast');
        
        clear a p maxy
        linkaxes(sub,'x')% link x axis of different plots (so that they change simultaneously)
        
%%% FIGURE 5 -  OVERLAND FLOW PROBLEM

%black and white figure
    figure('name','Overlandflow problem');
        subplot (2,1,1);
        p=plot(date,Prec);% graph of climate data
        xlabel('Time','fontsize',10);
        ylabel('Rainfall (mm)','fontsize',10);
        axis tight;
        NameArray = {'Color'};
        ValueArray = {'k'}';
        set(p,NameArray,ValueArray);
        xmin=date(3300);
        xmax=date(3340);
        axis([xmin,xmax,0,25]);
        set(gca,'box','off');%remove ticks in upper right side
        tickDates = xmin:10:xmax ; %// creates a vector of tick positions
        set(gca, 'XTick' , tickDates, 'XTickLabel' , datestr(tickDates,'dd/mm/yyyy') ) %//

        subplot (2,1,2);
        p=plot(date,[FiltOF,Q_MOF]);% graph of overland flow
        xlabel('Time','fontsize',10);
        ylabel('Overland flow (m³/s)','fontsize',10);
        axis tight;
        NameArray = {'Color'};
        ValueArray = {'k',[0.6 0.6 0.6]}';
        set(p,NameArray,ValueArray);
        xmin=date(3300);
        xmax=date(3340);
        axis([xmin,xmax,0,0.7]);
        set(gca,'box','off');%remove ticks in upper right side
        tickDates = xmin:10:xmax ; %// creates a vector of tick positions
        set(gca, 'XTick' , tickDates, 'XTickLabel' , datestr(tickDates,'dd/mm/yyyy') ) %//
        clear a p sub   
        
%color figure
    figure('name','Overlandflow problem');
        subplot (2,1,1);
        p=plot(date,Prec);% graph of climate data
        xlabel('Time','fontsize',10);
        ylabel('Rainfall (mm)','fontsize',10);
        axis tight;
        NameArray = {'Color'};
        ValueArray = {'k'}';
        set(p,NameArray,ValueArray);
        xmin=date(3300);
        xmax=date(3340);
        axis([xmin,xmax,0,25]);
        set(gca,'box','off');%remove ticks in upper right side
        tickDates = xmin:10:xmax ; %// creates a vector of tick positions
        set(gca, 'XTick' , tickDates, 'XTickLabel' , datestr(tickDates,'dd/mm/yyyy') ) %//

        subplot (2,1,2);
        p=plot(date,[FiltOF,Q_MOF]);% graph of overland flow
        xlabel('Time','fontsize',10);
        ylabel('Overland flow (m³/s)','fontsize',10);
        axis tight;
        NameArray = {'Color'};
        ValueArray = {'r','b'}';
        set(p,NameArray,ValueArray);
        xmin=date(3300);
        xmax=date(3340);
        axis([xmin,xmax,0,0.7]);
        set(gca,'box','off');%remove ticks in upper right side
        tickDates = xmin:10:xmax ; %// creates a vector of tick positions
        set(gca, 'XTick' , tickDates, 'XTickLabel' , datestr(tickDates,'dd/mm/yyyy') ) %//
        clear a p sub 
        
%FIG file
     figure('name','Overlandflow problem');
        subplot (2,1,1);
        p=plot(date,Prec);% graph of climate data
        xlabel('Time','fontsize',10);
        ylabel('Rainfall (mm)','fontsize',10);
        axis tight;
        NameArray = {'Color'};
        ValueArray = {'k'}';
        set(p,NameArray,ValueArray);
        xmin=date(3300);
        xmax=date(3340);
        axis([xmin,xmax,0,25]);
        set(gca,'box','off');%remove ticks in upper right side
        tickDates = xmin:10:xmax ; %// creates a vector of tick positions
        set(gca, 'XTick' , tickDates, 'XTickLabel' , datestr(tickDates,'dd/mm/yyyy') ) %//
        legend ('Rainfall');

        subplot (2,1,2);
        p=plot(date,[FiltOF,Q_MOF]);% graph of overland flow
        xlabel('Time','fontsize',10);
        ylabel('Overland flow (m³/s)','fontsize',10);
        axis tight;
        NameArray = {'Color'};
        ValueArray = {'r','b'}';
        set(p,NameArray,ValueArray);
        xmin=date(3300);
        xmax=date(3340);
        axis([xmin,xmax,0,0.7]);
        set(gca,'box','off');%remove ticks in upper right side
        tickDates = xmin:10:xmax ; %// creates a vector of tick positions
        set(gca, 'XTick' , tickDates, 'XTickLabel' , datestr(tickDates,'dd/mm/yyyy') ) %//
        legend ('Observed overland flow', 'Modelled overland flow');
        clear a p sub 
            
 
%%% FIGURE 6 -  BASEFLOW PROBLEM           

% black and white figure
 figure('name','Baseflow problem');
        sub(1)=subplot (3,1,1,'fontsize',10);
        [ax,h1,h2]=plotyy(date,[Prec,ETo],date, DPCatch);             
        NameArray = {'Color'};
        ValueArray = {'w','k'}'; 
        set(h1,NameArray,ValueArray);
        NameArray = {'LineStyle'};
        ValueArray = {'-','-'}'; 
        set(h1,NameArray,ValueArray);
        NameArray = {'Color'};
        ValueArray = {'[0.6 0.6 0.6]'}'; 
        set(h2,NameArray,ValueArray);        
        ylabel(ax(1),'Rainfall or ET0 (mm)','fontsize',10);
        ylabel(ax(2),'Deep percolation (mm)','fontsize',10);       
        xlabel(ax(1),'Time','fontsize',10);        
        xmin=date(1497);
        xmax=date(1800);
        axis(ax(1),[xmin,xmax,0,50]);
        axis(ax(2),[xmin,xmax,0,20]);        
        tickDates = xmin:60:xmax ; %// creates a vector of tick positions
        set(gca, 'XTick' , tickDates, 'XTickLabel' , datestr(tickDates,'mm/yyyy') ) %//
        set(gca,'box','off');%remove ticks in upper right side
        set(ax(2),'YDir','reverse','YColor','[0.6 0.6 0.6]')
        text(xmax-155,47,'(a)','fontsize',10);
        hold on
        b=bar(date, Prec);
        set(b,'FaceColor','k','edgecolor','k')

        
        sub(2)=subplot (3,1,2,'fontsize',10);
        P=plot(date,Wr2Catch,date,repmat(2*SoilPar(1,end),nAC,1));% graph of soil water content
        NameArray = {'LineStyle'};
        ValueArray = {'-','--'}';
        set(P,NameArray,ValueArray)
        NameArray = {'Color'};
        ValueArray = {'k','k'}';
        set(P,NameArray,ValueArray);
        xlabel('Time','fontsize',10);
        ylabel('Soil water content (mm/2 m soil depth)','fontsize',10);
        xmin=date(1497);
        xmax=date(1800);
        a=axis;%asks for axislimits of y axis
        maxy=a(1,4);% neemt maximum y limit
        axis([xmin,xmax,500,maxy]);
        set(gca,'box','off');%remove ticks in upper right side
        tickDates = xmin:60:xmax ; %// creates a vector of tick positions
        set(gca, 'XTick' , tickDates, 'XTickLabel' , datestr(tickDates,'mm/yyyy') ) %//
        text(xmax-155,maxy-10,'(b)','fontsize',10);
        
        sub(3)=subplot(3,1,3, 'fontsize',10);
        p=plot(date,Q_MBF,date_cl,FiltBF_cl);
        xlabel('Time','fontsize',10);
        ylabel('Baseflow (m³/s)','fontsize',10);
        NameArray = {'Color'};
        ValueArray = {'[0.6 0.6 0.6]','k'}';
        set(p,NameArray,ValueArray);
        xmin=date(1497);
        xmax=date(1800);
        axis([xmin,xmax,0,0.2]);
        set(gca,'box','off');%remove ticks in upper right side
        tickDates = xmin:60:xmax ; %// creates a vector of tick positions
        set(gca, 'XTick' , tickDates, 'XTickLabel' , datestr(tickDates,'mm/yyyy') ) %//
        text(xmax-155,0.19,'(c)','fontsize',10);
 
 
 % color figure
 figure('name','Baseflow problem');
        sub(1)=subplot (3,1,1,'fontsize',10);
        [ax,h1,h2]=plotyy(date,[Prec,ETo],date, DPCatch);             
        NameArray = {'Color'};
        ValueArray = {'w','k'}'; 
        set(h1,NameArray,ValueArray);
        NameArray = {'LineStyle'};
        ValueArray = {'-','-'}'; 
        set(h1,NameArray,ValueArray);
        NameArray = {'Color'};
        ValueArray = {'r'}'; 
        set(h2,NameArray,ValueArray);        
        ylabel(ax(1),'Rainfall or ET0 (mm)','fontsize',10);
        ylabel(ax(2),'Deep percolation (mm)','fontsize',10);       
        xlabel(ax(1),'Time','fontsize',10);        
        xmin=date(1497);
        xmax=date(1800);
        axis(ax(1),[xmin,xmax,0,50]);
        axis(ax(2),[xmin,xmax,0,20]);        
        tickDates = xmin:60:xmax ; %// creates a vector of tick positions
        set(gca, 'XTick' , tickDates, 'XTickLabel' , datestr(tickDates,'mm/yyyy') ) %//
        set(gca,'box','off');%remove ticks in upper right side
        set(ax(2),'YDir','reverse','YColor','r')
        text(xmax-155,47,'(a)','fontsize',10);
        hold on
        b=bar(date, Prec);
        set(b,'FaceColor','k','edgecolor','k')

        
        sub(2)=subplot (3,1,2,'fontsize',10);
        P=plot(date,Wr2Catch,date,repmat(2*SoilPar(1,end),nAC,1));% graph of soil water content
        NameArray = {'LineStyle'};
        ValueArray = {'-','--'}';
        set(P,NameArray,ValueArray)
        NameArray = {'Color'};
        ValueArray = {'k','k'}';
        set(P,NameArray,ValueArray);
        xlabel('Time','fontsize',10);
        ylabel('Soil water content (mm/2 m soil depth)','fontsize',10);
        xmin=date(1497);
        xmax=date(1800);
        a=axis;%asks for axislimits of y axis
        maxy=a(1,4);% neemt maximum y limit
        axis([xmin,xmax,500,maxy]);
        set(gca,'box','off');%remove ticks in upper right side
        tickDates = xmin:60:xmax ; %// creates a vector of tick positions
        set(gca, 'XTick' , tickDates, 'XTickLabel' , datestr(tickDates,'mm/yyyy') ) %//
        text(xmax-155,maxy-10,'(b)','fontsize',10);
        
        sub(3)=subplot(3,1,3, 'fontsize',10);
        p=plot(date,Q_MBF,date_cl,FiltBF_cl);
        xlabel('Time','fontsize',10);
        ylabel('Baseflow (m³/s)','fontsize',10);
        NameArray = {'Color'};
        ValueArray = {'b','r'}';
        set(p,NameArray,ValueArray);
        xmin=date(1497);
        xmax=date(1800);
        axis([xmin,xmax,0,0.2]);
        set(gca,'box','off');%remove ticks in upper right side
        tickDates = xmin:60:xmax ; %// creates a vector of tick positions
        set(gca, 'XTick' , tickDates, 'XTickLabel' , datestr(tickDates,'mm/yyyy') ) %//
        text(xmax-155,0.19,'(c)','fontsize',10);
 %FIG file
 figure('name','Baseflow problem');
        sub(1)=subplot (3,1,1,'fontsize',10);
        [ax,h1,h2]=plotyy(date,[Prec,ETo],date, DPCatch);             
        NameArray = {'Color'};
        ValueArray = {'w','k'}'; 
        set(h1,NameArray,ValueArray);
        NameArray = {'LineStyle'};
        ValueArray = {'-','-'}'; 
        set(h1,NameArray,ValueArray);
        NameArray = {'Color'};
        ValueArray = {'r'}'; 
        set(h2,NameArray,ValueArray);        
        ylabel(ax(1),'Rainfall or ET0 (mm)','fontsize',10);
        ylabel(ax(2),'Deep percolation (mm)','fontsize',10);       
        xlabel(ax(1),'Time','fontsize',10);        
        xmin=date(1497);
        xmax=date(1800);
        axis(ax(1),[xmin,xmax,0,50]);
        axis(ax(2),[xmin,xmax,0,20]);        
        tickDates = xmin:60:xmax ; %// creates a vector of tick positions
        set(gca, 'XTick' , tickDates, 'XTickLabel' , datestr(tickDates,'mm/yyyy') ) %//
        set(gca,'box','off');%remove ticks in upper right side
        set(ax(2),'YDir','reverse','YColor','r')
        text(xmax-155,47,'(a)','fontsize',10);
        hold on
        b=bar(date, Prec);
        set(b,'FaceColor','k','edgecolor','k')
        legend ([h1(2),h2,b],'Reference evapotranspiration','Deep percolation','Precipitation','Location','west')
        
        sub(2)=subplot (3,1,2,'fontsize',10);
        P=plot(date,Wr2Catch,date,repmat(2*SoilPar(1,end),nAC,1));% graph of soil water content
        NameArray = {'LineStyle'};
        ValueArray = {'-','--'}';
        set(P,NameArray,ValueArray)
        NameArray = {'Color'};
        ValueArray = {'k','k'}';
        set(P,NameArray,ValueArray);
        xlabel('Time','fontsize',10);
        ylabel('Soil water content (mm/2 m soil depth)','fontsize',10);
        xmin=date(1497);
        xmax=date(1800);
        a=axis;%asks for axislimits of y axis
        maxy=a(1,4);% neemt maximum y limit
        axis([xmin,xmax,500,maxy]);
        set(gca,'box','off');%remove ticks in upper right side
        tickDates = xmin:60:xmax ; %// creates a vector of tick positions
        set(gca, 'XTick' , tickDates, 'XTickLabel' , datestr(tickDates,'mm/yyyy') ) %//
        legend ('Soil water content','Field capacity','Location','west')
        text(xmax-155,maxy-10,'(b)','fontsize',10);
        
        sub(3)=subplot(3,1,3, 'fontsize',10);
        p=plot(date,Q_MBF,date_cl,FiltBF_cl);
        xlabel('Time','fontsize',10);
        ylabel('Baseflow (m³/s)','fontsize',10);
        legend ('Modeled baseflow','Observed baseflow', 'Location','west')
        NameArray = {'Color'};
        ValueArray = {'b','r'}';
        set(p,NameArray,ValueArray);
        xmin=date(1497);
        xmax=date(1800);
        axis([xmin,xmax,0,0.2]);
        set(gca,'box','off');%remove ticks in upper right side
        tickDates = xmin:60:xmax ; %// creates a vector of tick positions
        set(gca, 'XTick' , tickDates, 'XTickLabel' , datestr(tickDates,'mm/yyyy') ) %//
        text(xmax-155,0.19,'(c)','fontsize',10);
       
 %%% GRAPH ABSTRACT - TOTAL FLOW
  figure('name','Flow comparison');
        p=plot(date(StartTimeCalib):1:date(end),[ObsTot_fill(StartTimeCalib:end),Q_MTF(StartTimeCalib:end)]);% graph of total flow
        xlabel('Time','fontsize',10);
        ylabel('Total flow (m³/s)','fontsize',10);
        datetick('x','mm/yyyy')
        axis tight;
        NameArray = {'Color'};
        ValueArray = {'r','b'}';
        set(p,NameArray,ValueArray);
        set(gca,'box','off');%remove ticks in upper right side
        legend ('Observed discharge', 'Modelled discharge', 'Orientation', 'vertical', 'Location', 'northeast');
        
        
  %%% GRAPH PRESENTATION - AGGREGATED VOLUMES    
        
figure('name','DAILY');
sub(1)=subplot(3,1,1, 'fontsize',10);%DAILY
        p=plot(date(StartTimeCalib):1:date(end),[ObsTot_fill(StartTimeCalib:end)*(f/area),Q_MTF(StartTimeCalib:end)*(f/area)]);% graph of total flow
        xlabel('Time','fontsize',10);
        ylabel('Total flow (mm/day)','fontsize',10);
        datetick('x','mm/yyyy')
        axis tight;
        NameArray = {'Color'};
        ValueArray = {'r','b'}';
        set(p,NameArray,ValueArray);
        set(gca,'box','off');%remove ticks in upper right side
sub(2)= subplot(3,1,2, 'fontsize',10);%DECADE
        p=plot(1:1:453,[HydroCl_DecadeSub{1,3}(:,4),HydroCl_DecadeSub{1,8}(:,4)]);% graph of total flow
        xlabel('Time','fontsize',10);
        ylabel('Total flow (mm/10-days)','fontsize',10);
        axis tight;
        NameArray = {'Color'};
        ValueArray = {'r','b'}';
        set(p,NameArray,ValueArray);
        set(gca,'box','off');%remove ticks in upper right side
sub(3)=subplot(3,1,3, 'fontsize',10);%MONTHLY
        p=plot(1:1:151,[HydroCl_MonthSub{1,3}(:,3),HydroCl_MonthSub{1,8}(:,3)]);% graph of total flow
        xlabel('Time','fontsize',10);
        ylabel('Total flow (mm/month)','fontsize',10);
        axis tight;
        NameArray = {'Color'};
        ValueArray = {'r','b'}';
        set(p,NameArray,ValueArray);
        set(gca,'box','off');%remove ticks in upper right side   


test=HydroCl_DecadeSub{1,3};
test(:,5)=test(:,4)*(f/area);
test(:,6)=test(:,4)*(f/(11*area));