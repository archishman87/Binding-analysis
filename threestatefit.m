function [kdc,kdm,lp,lc,lm,b] = threestatefit(data1,data2,data3)
n=17;
dat1 = xlsread(data1);
dat2 = xlsread(data2);
dat3 = xlsread(data3);
mt = dat1(:,1);
ct = dat1(:,2);
d1 = dat1(:,3);
d2 = dat2(:,3);
d3 = dat3(:,3);
ptot = 0.2;
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
mtot = [mt;mt;mt];
ctot = [ct;ct;ct];
lobs = [d1;d2;d3];
options = fitoptions('Method','NonlinearLeastSquares',...
'Lower',[0.0001,0.0001,200,200,-Inf,-Inf],...
'Upper',[Inf,Inf,400,400,Inf,Inf],...
'Startpoint',[1 1 300 300 300 1],'algorithm','Trust-Region',...
'DiffMinChange',(1.0E-8),'DiffMaxChange',(0.1),'MaxFunEvals',(600),'MaxIter',(400),...
'TolFun',(1.0E-6),'TolX',(1.0E-6));
% Edit fitting options such as Startpoint, Lower and Upper limits etc.
f = fittype('(lp./(1+ctot/kdc)+(lc*ctot/kdc)./(1+ctot/kdc))+((lm+b*ctot)-(lp./(1+ctot/kdc)+(lc*ctot/kdc)./(1+ctot/kdc))).*((0.2+mtot+(kdm*(1+ctot/kdc)))-sqrt((0.2+mtot+(kdm*(1+ctot/kdc))).^2-4*0.2.*mtot))./(2*0.2)',...
'dependent',{'lobs'},'independent',{'mtot','ctot'},...
'coefficients',{'kdc', 'kdm', 'lp','lc','lm','b'});
% equation maybe customized. Current equation is for single site quadratic
% binding model with double independent parameters, i.e. total ligand and protein
% concentrations.
[c,gof,output] = fit([mtot,ctot],lobs,f,options); % fits the data with x and p as independents and f as the model.
% Returns parameters in c, goodness of fit in gof and other statistical
% parameters in output
kdc = c.kdc;
kdm = c.kdm;
lp = c.lp;
lc = c.lc;
lm = c.lm;
b = c.b;
conf = confint(c,0.66);
kdcerr = [num2str(kdc),'+/-',num2str(((kdc-conf(1,1))+(conf(2,1)-kdc))/2)];
kdmerr = [num2str(kdm),'+/-',num2str(((kdm-conf(1,2))+(conf(2,2)-kdm))/2)];
lperr = [num2str(lp),'+/-',num2str(((lp-conf(1,3))+(conf(2,3)-lp))/2)];
lcerr = [num2str(lc),'+/-',num2str(((lc-conf(1,4))+(conf(2,4)-lc))/2)];
lmerr = [num2str(lm),'+/-',num2str(((lm-conf(1,5))+(conf(2,5)-lm))/2)];
berr = [num2str(b),'+/-',num2str(((lm-conf(1,6))+(conf(2,6)-lm))/2)];
lobsfit = (lp./(1+ctot/kdc)+(lc*ctot/kdc)./(1+ctot/kdc))+((lm+b*ctot)-(lp./(1+ctot/kdc)+(lc*ctot/kdc)./(1+ctot/kdc))).*((0.2+mtot+(kdm*(1+ctot/kdc)))-
sqrt((0.2+mtot+(kdm*(1+ctot/kdc))).^2-4*0.2.*mtot))./(2*0.2); % plotting experimental values using newly obtained parameters
%lobsfit
ctot = ct; mtot = mt;
% all data
plot3d_errorbars(ctot(1:n),mtot(1:n),avg(1:n),serr(1:n),'-b','bo') % 0 g/L
hold on
plot3d_errorbars(ctot(n+1:2*n),mtot(n+1:2*n),avg(n+1:2*n),serr(n+1:2*n),'-c','co') % 50 g/L
plot3d_errorbars(ctot(2*n+1:3*n),mtot(2*n+1:3*n),avg(2*n+1:3*n),serr(2*n+1:3*n),'-g','go') % 100 g/L
plot3d_errorbars(ctot(3*n+1:4*n),mtot(3*n+1:4*n),avg(3*n+1:4*n),serr(3*n+1:4*n),'-m','mo') % 200 g/L
plot3d_errorbars(ctot(4*n+1:5*n),mtot(4*n+1:5*n),avg(4*n+1:5*n),serr(4*n+1:5*n),'-r','ro') % 300 g/L
plot3d_errorbars(ctot(5*n+1:5*n+20),mtot(5*n+1:5*n+20),avg(5*n+1:5*n+20),serr(5*n+1:5*n+20),'-k','ko') % apo ficoll titration
%plot3(ctot(5*n+21:end),mtot(5*n+21:end),avg(5*n+21:end),'ko') % holo ficoll titration optional
% all fits
plot3(ctot(1:n),mtot(1:n),lobsfit(1:n),'k')
plot3(ctot(n+1:2*n),mtot(n+1:2*n),lobsfit(n+1:2*n),'k')
plot3(ctot(2*n+1:3*n),mtot(2*n+1:3*n),lobsfit(2*n+1:3*n),'k')
plot3(ctot(3*n+1:4*n),mtot(3*n+1:4*n),lobsfit(3*n+1:4*n),'k')
plot3(ctot(4*n+1:5*n),mtot(4*n+1:5*n),lobsfit(4*n+1:5*n),'k')
plot3(ctot(5*n+1:5*n+20),mtot(5*n+1:5*n+20),lobsfit(5*n+1:5*n+20),'k')
%plot3(ctot(5*n+21:end),mtot(5*n+21:end),lobsfit(5*n+21:end),'k')
mtot(1)=0.001;
mtot(n+1)=0.001;
mtot(2*n+1)=0.001;
mtot(3*n+1)=0.001;
mtot(4*n+1)=0.001;
for i = 2:n
x = linspace(ctot(i),ctot(i+4*n),100); x = x';
% xi = interp(x,100);
y = linspace(mtot(i),mtot(i+4*n),100); y = y';
% yi = interp(y,100);
z = (lp./(1+x/kdc)+(lc*x/kdc)./(1+x/kdc))+((lm+b*x)-(lp./(1+x/kdc)+(lc*x/kdc)./(1+x/kdc))).*((0.2+y+(kdm*(1+x/kdc)))-sqrt((0.2+y+(kdm*(1+x/kdc))).^2-4*0.2.*y))./(2*0.2); % plotting experimental values using newly obtained parameters
%lobsfit
% z=[lobsfit(i),lobsfit(i+12),lobsfit(i+2*12),lobsfit(i+3*12)];
% zi = interp(z,100);
plot3(x,y,z,'k')
end
%title(['Kdc = ',kdcerr,' Kdm = ',kdmerr,' Fcmin = ',lcerr,' Fmmin = ',lperr,' Fmax = ',lmerr,' berr = ',berr])
set(gca, 'YScale', 'log')
% axis([0 4500 0.0001 10 343 351])
rotate3d on
end
