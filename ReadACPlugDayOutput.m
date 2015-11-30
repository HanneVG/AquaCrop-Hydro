% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script reads all the AquaCrop output from ....day.OUT files. 
%
% Warning: This script is built to read output files of simulations ran
% with AquaCrop version 5, PLUGIN VERSION (requesting daily output1-2-3)
%
%
%
% Author: Hanne Van Gaelen
% Last update: 30/11/2015
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [output]= ReadACPlugDayOutput(Datapath)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. DEFINE THE TYPE OF FILE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Define which type of files should be read (* is wild character)
Datafile=dir(fullfile(Datapath,'*day.out'));

%Define the format of data in this files
Readingformat = '%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f' ; 

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. READ ALL FILES & WRITE DATA IN A MATRIX
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for filenumb=1:length(Datafile)%loop over all files with extension *day.OUT

% Save filename
filename=Datafile(filenumb).name; %retrieve filename
filenamefull=fullfile(Datapath, filename); % create exact reference to file (with folders)
output{1,filenumb}= filename; %save filename into output structure

% Read file
[~,name,~]=fileparts(filenamefull); 
k=strfind(name,'PRMday');
if isempty(k)==1 %file is PROday.OUT (There are no runnumbers inside (all is one run)
     
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
            linecounter = 1; %linecounter is an indicator of on which row data should be save
                                     
            TextLine=fgetl(fid); %get a textline from the file 

               while isempty(TextLine)==0 %as long there is no blank space in TextLine (blank space indicating the end of a project)
                    try
                    LineData=sscanf(TextLine,Readingformat); % read number from the textline
                    output{2,filenumb} (linecounter,1)=1; % there is only one run, so runnr is 1
                    output{2,filenumb} (linecounter,2:52) = LineData'; % write the data of this line
                    TextLine=fgetl(fid); %get the next line 
                    linecounter = linecounter+1;% go to the next line
                    catch
                        break %if an error occurs in the block between try and catch (should only be at end of file, but can also be caused by invalid data)
                    end
               end

                   
    %2.3 close the file again before next file is read
            fclose (fid);

elseif isempty(k)==0 %file is PRMday.OUT
    
    %2.1 open file for reading    
        fid = fopen(filenamefull); 
            if fid==-1 % check if file was really opened
                disp ('File could not be opened')
            else    
                %carry on, file can now be read
            end

    %2.2 read this file    
       % read headerlines
         for i=1:5 
            TextLine=fgetl(fid);% read the first 4 lines but skip them
            clear TextLine
         end
    
     %2.3 read real data line by line
        linecounter = 1; %linecounter is an indicator of on which row data should be save

        for numberofruns=1:999 %999 is used as a the maximum number of runs (code will break out if loop is there is no more data)
            TextLine=fgetl(fid); %get a textline from the file 

            %Loop over a single AquaCrop simulation run
                while isempty(TextLine)==0 %as long there is no blank space in TextLine (blank space indicating the end of a run)
                    try
                    LineData=sscanf(TextLine,Readingformat); % read number from the textline
                    output{2,filenumb} (linecounter,1)=numberofruns; % write the runnumber for this line
                    output{2,filenumb} (linecounter,2:52) = LineData'; % write the data of this line
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
    
    
end
end


end