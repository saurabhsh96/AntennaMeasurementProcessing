function [S,Nport,Frequencies]=readsp2(FileName,Nport)

% Function to read Zeland *.sp files (Touchstone format). S-parameter files of
% Touchstone can also be read.
%
% Format:      [S,Nport,Frequencies]=readsp2(FileName,Nport)
%
% Input variables:
%
%       FileName:    Name of the SP-file (must include full path and .sp extension) 
%       Nport:       Number of ports. Must be specified for reading
%                    Touchstone files. Can be 0 for Zeland files
%                    (Nport will be read from a comment line in that case).
%
% Output variables:
%
%       S:           Multidimensional matrix containing the S-parameters
%                    Format: S(Freq_index,ReflPort,IncPort).
%       Nport:       Number of ports. is the same as the input variable Nport
%                    if that variable is not equal to 0.
%       Frequencies: Frequencies for which S-parameters are available.
%
A=[];
fid=fopen(FileName);
Line=fgets(fid);
while Line~=-1
   if Line(1)~='!' && Line(1)~='#' && length(Line)~=10
      Temp=sscanf(Line,'%e');
      if ~strcmp(Temp,'')
          A=[A ; Temp];
      end
   end
   if Line(1)=='#'
      FirstLine=Line;
   end
   Line=fgets(fid);
   if strncmp(Line,'! Nport =',9)
      if Nport==0
         Nport=sscanf(Line(10:length(Line)),'%e');
      end
   end
end
fclose(fid);
%
% Rearrange data into S-matrix
%
CountA=1;
N_freq=length(A)./(1+2.*Nport.^2);    % Determine number of frequencies
S=zeros(N_freq,Nport,Nport);
for Count_freq=1:N_freq
   Frequencies(Count_freq)=A(CountA);
   CountA=CountA+1;
   for ReflPort=1:Nport
      for IncPort=1:Nport
         if ~isempty(findstr(lower(FirstLine),'ma'))
            AbsValue=A(CountA);
         	Phase=A(CountA+1)./180.*pi;
         	CountA=CountA+2;
            S(Count_freq,IncPort,ReflPort)=AbsValue.*exp(j.*Phase);
         end
         if ~isempty(findstr(lower(FirstLine),'ri'))
            RealValue=A(CountA);
         	ImagValue=A(CountA+1);
         	CountA=CountA+2;
            S(Count_freq,IncPort,ReflPort)=RealValue+j.*ImagValue;
         end
         if ~isempty(findstr(lower(FirstLine),'db'))
            AbsValue=10^(A(CountA)./20);
         	Phase=A(CountA+1)./180.*pi;
         	CountA=CountA+2;
            S(Count_freq,IncPort,ReflPort)=AbsValue.*exp(j.*Phase);
         end
      end
   end
end
clear A;
end