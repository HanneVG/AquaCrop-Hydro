%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test for normality for each column in input vector
%
% Testname = Lillietest or Jarque-Bera test
% Ho= data are normal distributed
% test=1: Ho rejected -> data not normally distributed
% test=0: Ho accepted -> data normally distributed
%
% sign=significance level 
%
% Author: Hanne Van Gaelen
% Last updated: 08/01/2016
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[notnormal,normal]=NormalityCheck(Input,TestName,sign)

% number of columns
[~,ncol]=size(Input);
 normal=NaN; 
 notnormal=NaN;
 
for c=1:ncol
% do the requested statistical test
    if strcmp(TestName,'lillie')==1
        test=lillietest(Input(:,c),'Alpha',sign);
    elseif strcmp(TestName,'jb')==1
        test=jbtest(Input(:,c),sign);
    else
        error('Type of normality check test could not be recognized. Please try "lillie" or "jb"')
    end
        
% provide information on result in command window            

    
    if test==0 
        normal=[normal;c];
    elseif  test==1       
        notnormal=[notnormal;c];
        warning(['Data in column ' num2str(c) ' are not normally distributed'])        
    end   

end
    normal=normal(2:end,1);  
    notnormal=notnormal(2:end,1);


end
    
 