% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script reads all the AquaCrop output from ....Run.OUT files. 
%
% Warning: This script is built to read output files of simulations ran
% with AquaCrop version 5,normal interface (NOT Plugin version)
%
% Author: Hanne Van Gaelen
% Last update: 12/01/2016
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [output]= ReadACPlugSeasonOutput(Datapath)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. DEFINE THE TYPE OF FILE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Define which type of files should be read (* is wild character)
Datafile=dir(fullfile(Datapath,'*season.out'));

%Define the format of data in this files
Readingformat = '%*s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %*s' ;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. READ ALL FILES & WRITE DATA IN A MATRIX
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for filenumb=1:length(Datafile)%loop over all files with extension *Run.OUT

% Save filename
filename=Datafile(filenumb).name; %retrieve filename
filenamefull=fullfile(Datapath, filename); % create exact reference to file (with folders)
output{1,filenumb}= filename; %save filename into output structure

% Read file   
    %2.1 open file for reading    
        fid = fopen(filenamefull); 
            if fid==-1 % check if file was really opened
                disp ('File could not be opened')
            else    
                %carry on, file can now be read
            end

    %2.2 read this file    
       % read headerlines
         for i=1:4
            TextLine=fgetl(fid);% read the first 4 lines but skip them
            clear TextLine
         end

       % read real data line by line
            linecounter = 1; %linecounter is an indicator of on which row data should be saved
                                     
            TextLine=fgetl(fid); %get a textline from the file 

               while isempty(TextLine)==0 %as long there is no blank space in TextLine (blank space indicating the end of a project)
                    try
                    LineData=sscanf(TextLine,Readingformat); % read number from the textline
                    output{2,filenumb} (linecounter,1) = linecounter;
                    output{2,filenumb} (linecounter,2:35) = LineData'; % write the data of this line
                    TextLine=fgetl(fid); %get the next line 
                    linecounter = linecounter+1;% go to the next line
                    catch
                        break %if an error occurs in the block between try and catch (should only be at end of file, but can also be caused by invalid data)
                    end
               end
                   
    %2.3 close the file again before next file is read
            fclose (fid);
    
end
end

