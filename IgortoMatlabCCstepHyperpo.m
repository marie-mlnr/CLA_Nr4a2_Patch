%addpath '/Volumes/Marie/scripts'
%{ %This hidden part is used to run the function across different cells folder if needed
parentpath=pwd;

FilesToSave = dir('2*');
for ko = 1:length(FilesToSave) 
dirname=FilesToSave(ko).folder;   
dirnamebis=FilesToSave(ko).name;
Cellname=convertCharsToStrings(dirnamebis);
subdir=fullfile(dirname,dirnamebis);
cd(subdir)
clearvars -except parentpath FilesToSave ko

load('workspaceACSF.mat')
  if exist('HyperpoMinus20pA')
 %} 
 
 if exist('HyperpoMinus20pA') %about the hyperpolarizing step of -20pA
     
     NS= size(HyperpoMinus20pA);
    NbrofSteps=NS(1);

   if NbrofSteps>1 
   MeanHyperpoMinus20pA=mean(HyperpoMinus20pA(:,:));
   else 
   MeanHyperpoMinus20pA=(HyperpoMinus20pA(:,:));
   end 
   plot(MeanHyperpoMinus20pA)
   
   RestingMembranePot=mean(MeanHyperpoMinus20pA(1:9000))*1000;    % carefull this value isn't the real RMP as the cell is supposed to be around -70mV
   SteadyStateMembranePot=mean(MeanHyperpoMinus20pA(25000:55000))*1000; %cell potential when current is given
   DeltaMP= abs(SteadyStateMembranePot-RestingMembranePot); 
   InputResistance=(DeltaMP/20)*1000; %in Megaohm
   
   Width1stPeak=mean(nonzeros(ParametersofPeaks(:,3))); 
   MeanWidthAllPeaks=mean(nonzeros(ParametersofPeaks(:,4)));
   Amp1stPeak=mean(nonzeros(ParametersofPeaks(:,5)));
   MeanAmpAllPeaks=mean(nonzeros(ParametersofPeaks(:,6)));
   
  PercentofMax=(1-exp(-1))*DeltaMP ;
  PotentialOf1Tau=(RestingMembranePot-PercentofMax)/1000;
  
   MinMatrix(1,:)=abs(MeanHyperpoMinus20pA(10000:15000)-PotentialOf1Tau);
  [M,I] = min(MinMatrix); 
  IndexOfTau=10000+I; %since the index is looked for from 10000 so it nedds to be added here
  
  %plot(IndexOfTau,PotentialOf1Tau,'*r')
  
  Tauinms=(I*200)/10000;  %in ms (10000 point are 200ms)
  Capacitance=Tauinms/InputResistance;   %in nF
        
  Curvetofit=MeanHyperpoMinus20pA(10000:20000)'; %take only the decay part of the curve
  Curvetofitat0=Curvetofit-min(Curvetofit);
  
  TimeMatrix(1,1)=0.02;
  for k=2:length(Curvetofitat0)
      TimeMatrix(1,k)=TimeMatrix(1,(k-1))+0.02;
  end
   
  TimeMatrix=TimeMatrix';
  
  [ExpFit,gof3]=fit(TimeMatrix,Curvetofitat0,'exp1'); %matlab build in function
  if gof3.rsquare<0.8
      disp('CAREFULL rsquare is below 80%')
      carefulsentence=('CAREFULL rsquare is below 80%');
      save('BADEXPFIT','carefulsentence') 
  end 
  
  
  Tauinms=-1/(ExpFit.b);
  Capacitance=Tauinms/InputResistance; %in nF
  
  figure(2)
   plot(ExpFit,TimeMatrix,Curvetofitat0)
   RS=num2str(gof3.rsquare);
   legend('rsquare=',RS)
   title('Exp fit ACSF protocol')
  
 
 %save('CellIntrinsicProperties','RestingMembranePot', 'DeltaMP','InputResistance','Tauinms','Capacitance','RheoBase','Width1stPeak','MeanWidthAllPeaks','Amp1stPeak','MeanAmpAllPeaks')
 
 end 
 
% cd(parentpath) %needed to go through several cell folder
%end 
