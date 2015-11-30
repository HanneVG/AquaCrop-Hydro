% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script reads all the AquaCrop output from ....Prof.OUT files. 
%
% Warning: This script is built to read output files of simulations ran
% with AquaCrop version 5,normal interface (NOT Plugin version)
%
% Author: Hanne Van Gaelen
% Last update: 1/08/2015
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [output]= ReadACProfOutput(Datapath)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. DEFINE THE TYPE OF FILE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Define which type of files should be read (* is wild character)
Datafile=dir(fullfile(Datapath,'*Prof.out'));
%Define the format of data in this files
Readingformat = '%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f' ;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. READ ALL FILES & WRITE DATA IN A MATRIX
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
for filenumb=1:length(Datafile)%loop over all files with extension *Crop.OUT
   
    filename=Datafile(filenumb).name; %retrieve filename
    filenamefull=fullfile(Datapath, filename); % create exact reference to file (with folders)
    output{1,filenumb}= filename; %save filename into output structure
    
    %2.1 open file for reading    
        fid = fopen(filenamefull); 
            if fid==-1 % check if file was really opened
                disp ('File could not be opened')
            else    
                %carry on, file can now be read
            end

    %2.2 read this file    
       %2.2.1 read headerlines
         for i=1:6 
            TextLine=fgetl(fid);% read the first 6 lines but skip them
            clear TextLine
         end

       %2.2.1 read real data line by line
            linecounter = 1; %linecounter is an indicator of on which row data should be save
           
            for numberofruns=1:999 %999 is used as a the maximum number of runs (code will break out if loop is there is no more data)
                
                TextLine=fgetl(fid); %get a textline from the file 
                 
                %Loop over a single AquaCrop simulation run
                    while isempty(TextLine)==0 %as long there is no blank space in TextLine (blank space indicating the end of a run)
                        try
                        LineData=sscanf(TextLine,Readingformat); % read number from the textline
                        output{2,filenumb} (linecounter,1)=numberofruns; % write the runnumber for this line
                        output{2,filenumb} (linecounter,2:16) = LineData'; % write the data of this line
                        TextLine=fgetl(fid); %get the next line 
                        linecounter = linecounter+1;% go to the next line
                        catch
                            break %if an error occurs in the block between try and catch (should only be at end of file, but can also be caused by invalid data)
                        end
                    end

                %If a blank space has been detected the simulation run has ended
                    for i=1:3
                    TextLine=fgetl(fid); % just read the next three lines but do not save them
                    clear TextLine
                    end
            end % go to the next run
                   
%2.3 close the file again before next file is read
        fclose (fid);
      
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3. CLEAN UP RESULTING MATRIX
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% remove the last line of the outputfiles (this line is zeros & is an unnecessary artefact of the way the file is read

a=size(output);
nsim=a(1,2);
clear a;
    
for i=1:nsim %loop over all the simulaties
 output{2,i} (end,:)=[]; % remove last line
end

  
       

end










