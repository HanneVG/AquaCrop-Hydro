% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script runs the AquaCrop-Hydro model
%
%
% Theoretical background as well as a case-study example can be consulted
%  in the article:
%  Van Gaelen, H, Willems, P., Diels J. & Raes D. 2016. Modeling water
%  availability from field to catchemnt scale with a simple
%  agro-hydrological model. Environmental Modelling & Software [Submitted]
%
%  Documentation on the script can be found in the directory "Information",
%  Examples of the required input textfiles can be found in the directory
%  "Input"
%
%  TO DO before running the script: 
%   1. Run AquaCrop simulations for all land units in the catchment
%   2. Ensure that all AquaCrop output files are numbered with format "01", "02",
%        "03", "10","20" according to the landunit number
%   3. Prepare all required input files for this script:
%      SoilPar.txt, SimInfo.txt, Zrx.txt and Parameters.txt     
%   4. Adapt the script for the soil types that occur in the catchment
%        (see CatchmentOutput.m, section 5)
%   5. Specify the AquaCrop model you use (norma versus plugin)
%        (see AquaCropHydro.m, section 1)
%
%  Author: Hanne Van Gaelen
%  Last update: 30/11/2015
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
clear % clears all variables from a previous run
clc % clear command window

%% ------------------------------------------------------------------------
% 0. DEFINE DATA PATHS (BY USER)                                                    
%--------------------------------------------------------------------------      

% specify in which path all inputs for this script are stored including 
%    a) additional information on AquaCrop simulations (SimInfo.txt)
%    b) the parameters for the hydrological model (Parameters.txt)
%    d) the maximum root depth for every landunit and sim run (Zrx.txt)
%    d) the soil parameters of each soil type present (SoilPar.txt)
     DatapathInput = uigetdir('C:\','Select directory with all input files for AquaCrop-Hydro');
 
% specify in which path all AquaCrop simulation output files are stored

     DatapathAC=uigetdir('C:\','Select directory with AquaCrop output files for all landunits');

% specify in which path all AquaCrop-Hydro output (flow values) should be stored

     DatapathOutput=uigetdir('C:\','Select directory to store AquaCrop-Hydro output (flows)');
     
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
% 1. RUN AQUACROP-HYDRO                                                    %
%--------------------------------------------------------------------------  

[Q_MBF,Q_MIF,Q_MOF,Q_MTF,~,~,~,~,~,~,SimACOutput,CatchACOutput,~]=AquaCropHydro(DatapathAC, DatapathInput,ACMode);


%% ------------------------------------------------------------------------
% 2. WRITE OUTPUT TO EXCEL                                                   %
%-------------------------------------------------------------------------- 

% Combine output
OutputHeaders={'Baseflow','Interflow','Overland flow', 'Total flow'};
FlowOutput=[Q_MBF,Q_MIF,Q_MOF,Q_MTF];

% write output to one excel tabsheet
xlname='FlowSimResults.xlsx';
filename = fullfile(DatapathOutput,xlname);
xlswrite(filename,OutputHeaders,'SimFlow','A1');
xlswrite(filename,FlowOutput,'SimFlow','A2');

clear xlname filename FlowOutput OutputHeaders DatapathOutput DatapathInput DatapathAC 
