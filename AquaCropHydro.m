% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This scripts runs the submodels of the AquaCrop-Hydro model
%
% Author: Hanne Van Gaelen
% Last update: 30/11/2015
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function[Q_MBF,Q_MIF,Q_MOF,Q_MTF,area,f,Wrmin,Wrmax,pbf,SoilPar,SimACOutput,CatchACOutput,Par]=AquaCropHydro(DatapathAC, DatapathInput)


%% ------------------------------------------------------------------------
% 1. AquaCrop model
%--------------------------------------------------------------------------
       
       % AquaCrop is ran outside matlab (every project simulation is one landunit)
       % AquaCrop simulations need to be ran before running this script
       % output of AquaCrop should be stored in DataPathAC
       
       ACMode=2; % 1= normal AquaCrop, 2= plugin version
                    
%% ------------------------------------------------------------------------
% 2. Upscaling : fielscale to catchment scale
%------------------------------------------------------------------------
        
       % Load AquaCrop output for all simulation runs & calculate the summarized results for the catchment
        [SimACOutput,CatchACOutput, SoilPar,nTime]=CatchmentOutput(DatapathAC, DatapathInput,ACMode);  
        
       % Extract relevant soil water balance components for hydro model
        Wr2Catch=CatchACOutput(:,11);   % Soil water content in 2 m soil depth (mm)
        ROCatch=CatchACOutput(:,7);     % Runoff(mm)
        DPCatch=CatchACOutput(:,8);     % Deep percolation(mm)
      
%% ------------------------------------------------------------------------
% 3. Water routing model (hydrological model) 
%--------------------------------------------------------------------------
                
       % Load the parameters that are needed for hydrological model
        name='Parameters.txt';
        file = [DatapathInput name];
        Par = importdata(file); 
        clear name file
                           
       % Run the hydrological model
        [Q_MBF,Q_MIF,Q_MOF,Q_MTF,area,f,Wrmin,Wrmax,pbf]=Hydro(Par,SoilPar,nTime,ROCatch,DPCatch,Wr2Catch);
end