% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script performs calculations of performance statistics
% based on  series of observed versus modelled values
% 
% Warning: if either an observation or modelled value is missing, the
% record is removed from the series to calculate the statistics
%
%
% Author: Hanne Van Gaelen
% Last update: 1/08/2015
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Result]=PerformanceStat(Time,Obs,Mod,Sets,StartTimeSet2)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. PREPARE TO CALCULATE STATISTICS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
            Var=[Time, Obs, Mod]; % Put Time, Observed and Modeled in one matrix together
            
            Var=Var(~isnan(Var(:,2)),:);% Remove lines with missing values (=NaN value) for observed
            Var=Var(~isnan(Var(:,3)),:);% Remove lines with missing values (=NaN value) for modeled
            
            % split series in whole evaluation period, calibration period
            % and validation period
            Var_all=Var;
            
            if Sets==2 
            Var_calib=Var; 
            cond1= Time>=StartTimeSet2;
            Var_calib(cond1,:)=[];

            Var_valid=Var;
            cond2= Time<StartTimeSet2;
            Var_valid(cond2,:)=[];
            end
            
            clear cond1 cond2 Var
 %%           
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 % 2. CALCULATE STATS FOR WHOLE OBSERVATION PERIOD
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       Var_calc=Var_all; % define for which set statistics are calculated here
       Obs=Var_calc(:,2);% split matrix back into observed and modeled values
       Mod=Var_calc(:,3);
 
 %2.1 Calculate basics components for statistics
          E = Obs(:,1) - Mod(:,1);   % O-P 
          e = Mod(:,1) - Obs(:,1); % P-O
          SumE=sum(E);%Sum(O-P)
          SE = E.^2; % (O-P)²
          Se = e.^2; % (P-O)²
          AbsE = abs (E); % |O-P|
          SSE = sum(E.^2); %Sum(O-P)²
          SSe = sum(e.^2); %Sum(P-O)²
          SAbsE = sum (abs(E)); % Sum |O-P|
          u = mean(Obs(:,1)); %Oav
          v = mean(Mod(:,1)); %Pav
          SumObs=sum(Obs(:,1)); %Sum O
          U = (Obs(:,1) - u); %O-Oav
          P = (Mod (:,1) - u); %P -Oav 
          AbsP = abs (P);
          AbsU= abs (U);
          SU = (Obs(:,1) - u).^2; %(O-Oav)²
          SSU = sum((Obs(:,1) - u).^2);%Sum(O-Oav)²
          J = (Mod(:,1) - v); %P-Pav
          SJ = (Mod(:,1) - v).^2; %(P-Pav)²
          SSJ = sum((Mod(:,1) - v).^2);%Sum(P-Pav)²
          
          % define number of records for stats calculations
          a=size(Obs);
          b=size(Mod);
          n=a(1,1);
          n2=b(1,1);
          clear a b

            if n==n2% check match number of observed and modeled values
              clear n2
            else
              disp(' the number of observations is not equal to the number of simulations')
            end

%2.2 Calculate statistics      
        %Calculate NSE
            if n >= 2 
               NSE = 1 - SSE/SSU;

            else % cannot compute statistics
                error('Intesecting data resulted in too few elements to compute. \n Function has been terminated. If this is unexpected, \n check your index vectors of the two arrays.');
            end
        %Calculate RMSE
            if n >= 2 
              RMSE = (SSe/n).^0.5;

            else % cannot compute statistics
                error('Intesecting data resulted in too few elements to compute. \n Function has been terminated. If this is unexpected, \n check your index vectors of the two arrays.');
            end
        %Calculate (CV)RMSE 
            if n >= 2 
              CVRMSE = RMSE/u;

            else % cannot compute statistics
                error('Intesecting data resulted in too few elements to compute. \n Function has been terminated. If this is unexpected, \n check your index vectors of the two arrays.');
            end
        %Calculate d
            if n >= 2 

                teller = SSe;
                noemer = sum ((AbsP + AbsU).^2);
                d= 1-(teller./noemer);        
            else % cannot compute statistics
                error('Intesecting data resulted in too few elements to compute. \n Function has been terminated. If this is unexpected, \n check your index vectors of the two arrays.');
            end

        %Calculate ABSERR
            if n >= 2 
               ABSERR= SAbsE/n;       
            else % cannot compute statistics
                error('Intesecting data resulted in too few elements to compute. \n Function has been terminated. If this is unexpected, \n check your index vectors of the two arrays.');
            end
        %Calculate R²
           if n >= 2 

                    teller = sum(U.*J);
                    noemer = (SSU.*SSJ).^0.5;
                    R2= (teller./noemer).^2;        
                else % cannot compute statistics
                    error('Intesecting data resulted in too few elements to compute. \n Function has been terminated. If this is unexpected, \n check your index vectors of the two arrays.');
            end
        %Calculate MSE
            if n >= 2
               MSE=SSE/n; 
            else
            error('Intesecting data resulted in too few elements to compute. \n Function has been terminated. If this is unexpected, \n check your index vectors of the two arrays.');
            end
         %calculate WABAL error   
         if n >= 2
               WBERR=SumE/SumObs; 
            else
            error('Intesecting data resulted in too few elements to compute. \n Function has been terminated. If this is unexpected, \n check your index vectors of the two arrays.');
         end
         
% 2.3 Write stats in a results matrix
         Result=[n, R2, NSE, d, RMSE, CVRMSE, ABSERR, MSE, WBERR];
 %%        
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
 % 3. CALCULATE STATS FOR CALIBRATION AND VALIDATION PERIOD SEPERATELY
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 if Sets==2    
        
 for i=1:Sets % loop trough the different subsets
        
        if i==1; % define which set to use for every loop
            Var_calc=Var_calib;
        elseif i==2;
            Var_calc=Var_valid;
        end
            
       Obs=Var_calc(:,2);% split matrix back into observed and modeled values
       Mod=Var_calc(:,3);
                        
%3.1 Calculate basics components for statistics
          E = Obs(:,1) - Mod(:,1);   % O-P 
          e = Mod(:,1) - Obs(:,1); % P-O
          SumE=sum(E);%Sum(O-P)
          SE = E.^2; % (O-P)²
          Se = e.^2; % (P-O)²
          AbsE = abs (E); % |O-P|
          SSE = sum(E.^2); %Sum(O-P)²
          SSe = sum(e.^2); %Sum(P-O)²
          SAbsE = sum (abs(E)); % Sum |O-P|
          u = mean(Obs(:,1)); %Oav
          v = mean(Mod(:,1)); %Pav
          SumObs=sum(Obs(:,1)); %Sum O
          U = (Obs(:,1) - u); %O-Oav
          P = (Mod (:,1) - u); %P -Oav 
          AbsP = abs (P);
          AbsU= abs (U);
          SU = (Obs(:,1) - u).^2; %(O-Oav)²
          SSU = sum((Obs(:,1) - u).^2);%Sum(O-Oav)²
          J = (Mod(:,1) - v); %P-Pav
          SJ = (Mod(:,1) - v).^2; %(P-Pav)²
          SSJ = sum((Mod(:,1) - v).^2);%Sum(P-Pav)²
          
          % define number of records for stats calculations
          a=size(Obs);
          b=size(Mod);
          n=a(1,1);
          n2=b(1,1);
          clear a b

            if n==n2 % check match number of observed and modeled values
              clear n2
            else
              error('the number of observations is not equal to the number of simulations')
            end
    
%3.2 Calculate statistics      
        %Calculate NSE
            if n >= 2 
               NSE = 1 - SSE/SSU;

            else % cannot compute statistics
                error('Intesecting data resulted in too few elements to compute. \n Function has been terminated. If this is unexpected, \n check your index vectors of the two arrays.');
            end
        %Calculate RMSE
            if n >= 2 
              RMSE = (SSe/n).^0.5;

            else % cannot compute statistics
                error('Intesecting data resulted in too few elements to compute. \n Function has been terminated. If this is unexpected, \n check your index vectors of the two arrays.');
            end
        %Calculate (CV)RMSE 
            if n >= 2 
              CVRMSE = RMSE/u;

            else % cannot compute statistics
                error('Intesecting data resulted in too few elements to compute. \n Function has been terminated. If this is unexpected, \n check your index vectors of the two arrays.');
            end
        %Calculate d
            if n >= 2 % if only one row you can not calculate RMSE

                teller = SSe;
                noemer = sum ((AbsP + AbsU).^2);
                d= 1-(teller./noemer);        
            else % cannot compute statistics
                error('Intesecting data resulted in too few elements to compute. \n Function has been terminated. If this is unexpected, \n check your index vectors of the two arrays.');
            end

        %Calculate ABSERR
            if n >= 2 
               ABSERR= SAbsE/n;       
            else % cannot compute statistics
                error('Intesecting data resulted in too few elements to compute. \n Function has been terminated. If this is unexpected, \n check your index vectors of the two arrays.');
            end
        %Calculate R²
            if n >= 2 

                    teller = sum(U.*J);
                    noemer = (SSU.*SSJ).^0.5;
                    R2= (teller./noemer).^2;        
                else % cannot compute statistics
                    error('Intesecting data resulted in too few elements to compute. \n Function has been terminated. If this is unexpected, \n check your index vectors of the two arrays.');
            end
        %Calculate MSE
            if n >= 2
               MSE=SSE/n; 
            else
            error('Intesecting data resulted in too few elements to compute. \n Function has been terminated. If this is unexpected, \n check your index vectors of the two arrays.');
            end
       %Calculate WABAL error   
         if n >= 2
               WBERR=SumE/SumObs; 
            else
            error('Intesecting data resulted in too few elements to compute. \n Function has been terminated. If this is unexpected, \n check your index vectors of the two arrays.');
         end
         
%3.3 Add stats to the results matrix
         Result=[Result;n, R2, NSE, d, RMSE, CVRMSE, ABSERR, MSE, WBERR]; 
 end
 end
 
end
 


