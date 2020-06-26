function [] = waverror(ficollconc,data1,data2,data3,col11,col12,col21,col22,col31,col32,erfitcolor,ercolor,titration)
d1 = wavdataplot2(data1,col11,col12,erfitcolor,titration); % calls wavdataplot2 function, refer wavdataplot2.m
d2 = wavdataplot2(data2,col21,col22,erfitcolor,titration);
d3 = wavdataplot2(data3,col31,col32,erfitcolor,titration);
titr = xlsread(titration); % reads excel file
lig = titr(:,2); % extracts second column which should contain ligand concentrations
pro = titr(:,1);
serr = [];
avg = [];
for i = 1:length(d1)
sd = std([d1(i),d2(i),d3(i)]);
m = mean([d1(i),d2(i),d3(i)]);
se = sd/sqrt(3-1);
serr = [serr;se];
avg = [avg;m];
end
% semilogx(lig,avg,'o',lig,(avg-serr),'+',lig,(avg+serr),'+')
if (ishold == 0)
hold
end
errorbar(lig,avg,serr,ercolor)
errorbarlogx(0.01)
p = [pro;pro;pro];
x = [lig;lig;lig];
y = [d1;d2;d3] ; % fdata is the lambda max values obtained from wavdataplot2 function
options = fitoptions('Method','NonlinearLeastSquares',...
'Lower',[0.00000000,100,100],...
'Upper',[Inf,Inf,Inf],...
'Startpoint',[0 200 200],'algorithm','Trust-Region',...
'DiffMinChange',(1.0E-8),'DiffMaxChange',(0.1),'MaxFunEvals',(600),'MaxIter',(400),...
'TolFun',(1.0E-6),'TolX',(1.0E-6));
% Edit fitting options such as Startpoint, Lower and Upper limits etc.
f = fittype('y0 + (ym - y0)*((p+x+k)-sqrt((p+x+k).^2-4*p.*x))./(2*p)',...
'dependent',{'y'},'independent',{'x','p'},...
'coefficients',{'k', 'y0', 'ym'});
% equation maybe customized. Current equation is for single site quadratic
% binding model with double independent parameters, i.e. total ligand and protein
% concentrations.
[c,gof,output] = fit([x,p],y,f,options); % fits the data with x and p as independents and f as the model.
% Returns parameters in c, goodness of fit in gof and other statistical
% parameters in output
k = c.k;
y0 = c.y0;
ym = c.ym;
yfit = y0 + (ym - y0)*((p+x+k)-sqrt((p+x+k).^2-4*p.*x))./(2*p); % plotting experimental values using newly obtained parameters
% fbfit = (yfit-y0)/(ym-y0); % use by uncommenting if calculating fraction bound values
% fb = (y-y0)/(ym-y0);
fbfit = yfit(1:length(yfit)/3); % comment out if using fraction bound
fb = y;
%'RSquared is',gof.rsquare
%resid = output.residuals;
% 'Number of iterations:',output.iterations
% 'THE DISSOCIATION CONSTANT IS and F0 and Fmax:',[k y0 ym]
conf = confint(c,0.66)
kerr = [num2str(k),'+/-',num2str(((k-conf(1,1))+(conf(2,1)-k))/2)]
y0err = [num2str(y0),'+/-',num2str(((y0-conf(1,2))+(conf(2,2)-y0))/2)]
ymerr = [num2str(ym),'+/-',num2str(((ym-conf(1,3))+(conf(2,3)-ym))/2)]
semilogx(lig,fbfit,erfitcolor)
xlabel('Ligand concentration in uM')
ylabel('Wavelength in nm') % type in fraction bound if using it
line([k k],[min(fbfit) max(fbfit)],'Color','black','Linestyle','--') % plots a line where the value of k is
line([conf(1,1) conf(2,1)],[max(fbfit) max(fbfit)],'Color','black','Linestyle','-')
% if (conf(1,1)>0)
% line([conf(1,1) conf(2,1)],[(min(fbfit)+max(fbfit))/2 (min(fbfit)+max(fbfit))/2],'Color','black','Linestyle','--')
% else
% line([0.00000001 conf(2,1)],[(min(fbfit)+max(fbfit))/2 (min(fbfit)+max(fbfit))/2],'Color','black','Linestyle','--')
% end
title(['Kd = ',kerr,' Fmin = ',y0err,' Fmax = ',ymerr])
hold
ficollconcuM = [];
89
for i = 1:length(lig)
ficollconcuM = [ficollconcuM;ficollconc*10^6/70000];
end
col_header={'[protein]','[ligand]','[crowder]','wavpeak','fitwavpeak','Kd = '}; %Row cell array (for column labels)
xlswrite([num2str(ficollconc),'gperLf-','wavdata.xls'],col_header,'Sheet1','A1');
xlswrite([num2str(ficollconc),'gperLf-','wavdata.xls'],[pro lig ficollconcuM avg],'Sheet1','A2'); %Write data
% col_header={'fitwavpeak','Kd = '}; %Row cell array (for column labels)
xlswrite([num2str(ficollconc),'gperLf-','wavdata.xls'],fbfit,'Sheet1','E2'); %Write data
xlswrite([num2str(ficollconc),'gperLf-','wavdata.xls'],k,'Sheet1','F2');
% xlswrite([num2str(ficollconc),'gperLf-','wavdata.xls'],col_header,'Sheet1','E1'); %Write column header
end
