% FOR CALLING WITHIN wavdataplot2.m
% if hyper and hypochromic shifts is the measured parameter instead of
% batho or hysochromic shifts then replace xmax with ymax
% function [ymax,ymmax,y0max] = wavdataimport2(data,col1,col2,col3,col4,col5,col6) % input csv file and mention columns
function xmax = wavdataimport2(data,col1,col2)
d = importdata(data);
dd = d.data;
x = dd(:,col1);
y = dd(:,col2);
p = polyfit(x,y,5); % fits to 5th order polynomial
% Numerical maximum finding
xi = interp(x,100);% interpolation refines data
q = polyval(p,xi);
% Comment out the following two lines if NOT using wavelength shifts
ymax = max(q);
xmax = xi(q==ymax);
% Comment out the following 3 lines if NOT using Intensity shifts and
% uncomment the previous one for wavelength
%[~,i] = min(abs(xi-335));
%xmax = xi(i);
%ymax = q(i);
%plot(xi,q,color)
end
