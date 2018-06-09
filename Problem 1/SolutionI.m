%Empirical Macro I - First Homework
clear
path(path,'U:\EmpiricalMacroI')
%
%Question I
%Uploading the Data
X=xlsread('U:\EmpiricalMacroI\Winter2017\Homeworks\DatasetHomeworkI.xls','Dataset', 'A16:D173')
[T,N]=size(X);
%Column 1
Time=X(2:T,1);
%Column 2
logRealGDP=log(X(:,2));
RealGDPGrowth=diff(logRealGDP)*400;
%colum 3
logGDPDeflator=log(X(:,3));
Inflation=diff(logGDPDeflator)*400;
%Column 4
logM2=log(X(:,4));
M2Growth=diff(logM2)*400;
%The Data Matrix
X=[RealGDPGrowth Inflation M2Growth];
[T,N]=size(X);
%
%1)
%selecting the lag order with Akaike and Schwartz criteria
MaxLagOrder=8
AIC=varplagselect(X,'Y',MaxLagOrder,'AIC');
SIC=varplagselect(X,'Y',MaxLagOrder,'SIC');
LagOrder=max(AIC,SIC)
%
%Estimation of VAR for the log-differences of the three series
[B,VARB,U,V,vechV,VARvechV]=varp(X,LagOrder,'Y');
U=U';
%Checking for stationarity
[f,F,VAR]=Companion(B,V,N,LagOrder);
[VV,DD] = eig(F);
MaxAbsEigenvalue=max(abs(diag(DD)));
if MaxAbsEigenvalue<1
    disp('Estimated VAR is stationary')
else
    disp('Estimated VAR is non-stationary')
end
%
%2)
randn('state',123);
NN=10000;
BB=zeros(size(B,1),size(B,2),NN);
VV=zeros(size(V,1),size(V,2),NN);
%
xx=1;
while xx<=NN
    xx
    % (i) Bootstrapping X:
    BootstrappedX=zeros(N,LagOrder+T+100);
    BootstrappedX(:,1:LagOrder)=X(1:LagOrder,:)';
    for tt=LagOrder+1:LagOrder+T+100
        BootstrappedX(:,tt)=B*[1; vec(fliplr(BootstrappedX(:,tt-LagOrder:tt-1)))]+bootstrap(U)';
    end
    BootstrappedX=BootstrappedX(:,LagOrder+100+1:LagOrder+100+T)';
    % (ii) Estimating the VAR:
    [b,varb,u,v,vechv,varvechv]=varp(BootstrappedX,LagOrder,'Y');
    % (iii) Checking for stationarity:
    [f,ff,var]=Companion(b,v,N,LagOrder);
    [vv,dd] = eig(ff);
    maxabseigenvalue=max(abs(diag(dd)));
    if maxabseigenvalue<1
        BB(:,:,xx)=b;
        VV(:,:,xx)=v;
        xx=xx+1;
    end
end
%Centering
MedianBB=median(BB,3);
MedianVV=median(VV,3);
%
for hh=1:size(B,1)
    for kk=1:size(B,2)
        BB(hh,kk,:)=BB(hh,kk,:)-MedianBB(hh,kk)+B(hh,kk);
    end
end
for hh=1:size(V,1)
    for kk=1:size(V,2)
        VV(hh,kk,:)=VV(hh,kk,:)-MedianVV(hh,kk)+V(hh,kk);
    end
end
%
% (i) Comparing the innovation variances of GDP growth and inflation:
BootstrappedV11=squeeze(VV(1,1,:));
BootstrappedV22=squeeze(VV(2,2,:));
% Then we check:
sum(BootstrappedV11>BootstrappedV22)/NN
disp('Its not statistically significant at the conventional levels ')
%
% (ii) Comparing the innovation variances of the M2Growth and Inflation
BootstrappedV33=squeeze(VV(3,3,:));
sum(BootstrappedV33>BootstrappedV22)/NN
disp('Its not statistically significant at conventional levels')
%
disp('End of the first Question')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Question 2
%Uploading the Data
clear
path(path,'U:\EmpiricalMacroI')
X=xlsread('U:\EmpiricalMacroI\Winter2017\Homeworks\DatasetHomeworkI.xls','Dataset', 'G16:K283')
[T,N]=size(X);
%
%Column 1
Time=X(2:T,1);
%Column 2
logRealGDP=log(X(:,2));
RealGDPGrowth=diff(logRealGDP)*400;
%colum 3
logGDPDeflator=log(X(:,3));
Inflation=diff(logGDPDeflator)*400;
%colum 4
UnemploymentRate=X(2:T,4);
%Column 5
TreeMonthTreasuryBill=X(2:T,5);
%The Data Matrix
X=[RealGDPGrowth Inflation UnemploymentRate TreeMonthTreasuryBill];
[T,N]=size(X);
%
%1)
%selecting the lag order with Akaike and Schwartz criteria
MaxLagOrder=8
AIC=varplagselect(X,'Y',MaxLagOrder,'AIC');
SIC=varplagselect(X,'Y',MaxLagOrder,'SIC');
LagOrder=max(AIC,SIC)
%
%Estimation of VAR for the log-differences of the three series
[B,VARB,U,V,vechV,VARvechV]=varp(X,LagOrder,'Y');
U=U'
%
%Checking for stationarity
[f,F,VAR]=Companion(B,V,N,LagOrder);
[VV,DD] = eig(F);
MaxAbsEigenvalue=max(abs(diag(DD)));
if MaxAbsEigenvalue<1
    disp('Estimated VAR is stationary')
else
    disp('Estimated VAR is non-stationary')
end
%
% 2) Generating Forecasts
ForecastHorizon=10*4;
NumberOfBootstrapReplications=10000;
ForecastedX=zeros(1+ForecastHorizon,N,NumberOfBootstrapReplications);
MU=B(:,1);
BB=B(:,2:size(B,2));
B1=BB(:,1:N);
%
xx=1;
while xx<=NumberOfBootstrapReplications
    forecastedx(:,1:LagOrder)=X(T-LagOrder+1:T,:)';
    for tt=LagOrder+1:LagOrder+1+ForecastHorizon
        shocks=bootstrap(U)';
        forecastedx(:,tt)=MU+BB*vec(fliplr(forecastedx(:,tt-LagOrder:tt-1)))+shocks;
    end
    ForecastedX(:,:,xx)=forecastedx(:,LagOrder+1:LagOrder+1+ForecastHorizon)';
    xx=xx+1;
end
%
% Extracting various projections:
% The order of the variables:
%X=[RealGDPGrowth Inflation UnemploymentRate TreeMonthTreasuryBill];
ForecastedRealGDPGrowth=squeeze(ForecastedX(:,1,:));
ForecastedInflation=squeeze(ForecastedX(:,2,:));
ForecastedUnemploymentRate=squeeze(ForecastedX(:,3,:));
ForecastedTreeMonthTreasuryBill=squeeze(ForecastedX(:,4,:));
%
% Ordering the projections:
ForecastedRealGDPGrowth=sort(ForecastedRealGDPGrowth,2);
ForecastedInflation=sort(ForecastedInflation,2);
ForecastedUnemploymentRate=sort(ForecastedUnemploymentRate,2);
ForecastedTreeMonthTreasuryBill=sort(ForecastedTreeMonthTreasuryBill,2);
%
% Getting the percentiles:
Percentiles=NumberOfBootstrapReplications*[0.16 0.84 0.05 0.95]';
ForecastedRealGDPGrowth=ForecastedRealGDPGrowth(:,Percentiles);
ForecastedInflation=ForecastedInflation(:,Percentiles);
ForecastedUnemploymentRate=ForecastedUnemploymentRate(:,Percentiles);
ForecastedTreeMonthTreasuryBill=ForecastedTreeMonthTreasuryBill(:,Percentiles);
%
TimeForecast=(Time(T)):0.25:(Time(T)+0.25*ForecastHorizon);
%
% Plotting the 16th, 84th, 5th and the 95th percentile of the bootstrapped
%distributions for the forecasts:
%
Zero=zeros(size(Time));
ZeroForecast=zeros(size(TimeForecast));
%
figure(1)
subplot(2,2,1)
plot(Time,X(:,1),'k',TimeForecast,ForecastedRealGDPGrowth(:,1:2),'r:',TimeForecast,ForecastedRealGDPGrowth(:,3:4),'b',Time,Zero,'k:',TimeForecast,ZeroForecast,'k:','LineWidth',2)
xlim([Time(1) TimeForecast(length(TimeForecast))])
title('Real GDP growth')
subplot(2,2,2)
plot(Time,X(:,2),'k',TimeForecast,ForecastedInflation(:,1:2),'r:',TimeForecast,ForecastedInflation(:,3:4),'b',Time,Zero,'k:',TimeForecast,ZeroForecast,'k:','LineWidth',2)
xlim([Time(1) TimeForecast(length(TimeForecast))])
title('Inflation')
subplot(2,2,3)
plot(Time,X(:,3),'k',TimeForecast,ForecastedUnemploymentRate(:,1:2),'r:',TimeForecast,ForecastedUnemploymentRate(:,3:4),'b',Time,Zero,'k:',TimeForecast,ZeroForecast,'k:','LineWidth',2)
xlim([Time(1) TimeForecast(length(TimeForecast))])
title('UnemploymentRate')
subplot(2,2,4)
plot(Time,X(:,4),'k',TimeForecast,ForecastedTreeMonthTreasuryBill(:,1:2),'r:',TimeForecast,ForecastedTreeMonthTreasuryBill(:,3:4),'b',Time,Zero,'k:',TimeForecast,ZeroForecast,'k:','LineWidth',2)
xlim([Time(1) TimeForecast(length(TimeForecast))])
title('TreeMonthTreasuryBill')
%
%Probability that, in 2020Q1, Inflation is going to be negative:
ForecastedInflation2020Q1=squeeze(ForecastedX(20,2,:));
ForecastedInfaltion2020Q1=sort(ForecastedInflation2020Q1)
figure(2)
hist(ForecastedInflation2020Q1,1000)
sum(ForecastedInflation2020Q1<0)/NumberOfBootstrapReplications
disp('The probability that, in 2020Q1, inflation is negative is 13.82%')
%
disp('End of the second Question')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

