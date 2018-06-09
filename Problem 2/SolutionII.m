%Homework 2 - Empirical Macro I
clear
path(path,'U:\EmpMJS');
X=xlsread('U:\EmpMJS\DatasetHomeworkII.xls','Dataset','A17:K552');
[T,N]=size(X);
Time=X(2:T,1);

%log-differences
IndustrialProductionIndex=X(:,2);
LogDiffIndustrialProductionIndex=((diff(log(IndustrialProductionIndex))+1).^12-1)*100;
LogStockPriceIndex=((diff(log(X(:,6)))+1).^12-1)*100;
PCEDeflator=X(:,3)
LogPCEDeflator=((diff(log(X(:,3)))+1).^12-1)*100;
MedianHousePrices=X(:,10);
RealHousePrices=MedianHousePrices./PCEDeflator;
LogRealHousePrices=((diff(log(RealHousePrices))+1).^12-1)*100;

%Levels
UnemploymentRate=X(2:T,4);
VacacyRate=X(2:T,9);
CivilianNonindustrialPopulation=X(2:T,7);
AggregateWeeklyHours=X(2:T,8);
LogAggregateWeeklyHoursPerCapita=log(AggregateWeeklyHours./CivilianNonindustrialPopulation);
LogHousingStarts=log(X(2:T,11));
FederalFundsRate=X(2:T,5);

%The Dataset 
X=[LogDiffIndustrialProductionIndex LogStockPriceIndex LogPCEDeflator LogRealHousePrices UnemploymentRate VacacyRate LogAggregateWeeklyHoursPerCapita LogHousingStarts FederalFundsRate];
[T,N]=size(X);

%Identify a monetary policy shock --> reordering the variables because 
%the structure of Cholesky is lower triangular. The variables which do not
%react comtemporaneously to the shock must be ordered before the Fed Funds
%Rate (used to identify the shock) and the variables which react
%contemporaneously the the shock must be ordered after the interest rate. 

%The ordered Dataset
X=[LogDiffIndustrialProductionIndex LogPCEDeflator LogRealHousePrices UnemploymentRate VacacyRate LogAggregateWeeklyHoursPerCapita LogHousingStarts FederalFundsRate LogStockPriceIndex];
[T,N]=size(X);

%Identifying the the monetary shock on the restriction that , within a
%month it only impacts upon
%i) the federal funds rate
%ii) stock prices

%Estimating the Var with 6 lags
LagOrder=6;
[B,VARB,U,V,vechV,VARvechV]=varp(X,LagOrder,'Y');

%Checking for Stationarity
[f,F,VAR]=Companion(B,V,N,LagOrder);
% Here the eigenvalue-eigenvector decomposition of F:
[VV,DD] = eig(F);
MaxAbsEigenvalue=max(abs(diag(DD)));
if MaxAbsEigenvalue<1
    disp('Estimated VAR is stationary')
else
    disp('Estimated VAR is non-stationary')
end

% Extracting the relevant objects
MU=B(:,1);
B=B(:,2:size(B,2));

%A0 is the structural impact matrix
A0_SimpleEstimate=chol(V)';

%Checking if it is a matrix of zeros 
V;
A0_SimpleEstimate*A0_SimpleEstimate'-V 

%Identifying the shocks
A0_SimpleEstimate;
% The monetary policy shock has a positive effect on the FedFundsRate and a
% negative effect on the Stock Price Growth.

% the structural shocks
E=A0_SimpleEstimate\U;

%Now the impulse-response function
Horizon=10*12; % 10 years ahead
IRFs=GetIRFs(B,A0_SimpleEstimate,N,LagOrder,Horizon,N-1);

%Now we have to identify highly volatile variables --> plotting
figure(1)
for xx=1:N
    subplot(2,5,xx)
    plot(Time,X(:,xx))
    if xx==1
        title('LogDiffIndustrialProductionIndex')
    elseif xx==2
        title('LogPCEDeflator')
    elseif xx==3
        title('LogRealHousePrices')
    elseif xx==4
        title('UnemploymentRate')
    elseif xx==5
        title('VacacyRate')
    elseif xx==6
        title('LogAggregateWeeklyHoursPerCapita')  
    elseif xx==7
        title('LogHousingStarts')
    elseif xx==8
        title('FederalFundsRate')
    elseif xx==9
        title('LogStockPriceIndex')
    end
end
%The variables on the positions 1,2,3 and 9 seem to be highly volatile.
XX=X;
Index=[1 2 3 N]'; %the highly volatile variables
for xx=1:size(Index)
    XX(:,Index(xx))=RollingMean([zeros(11,1); XX(:,Index(xx))],12);
end

for xx=1:size(Index)
    IRFs(:,Index(xx))=RollingMean([zeros(11,1); IRFs(:,Index(xx))],12);
end
IRFs_SimpleEstimate=IRFs;

%Variance Decomposition
FractionsOfVariance_SimpleEstimate=GetFractionsOfVariance([MU B],A0_SimpleEstimate,N,LagOrder,Horizon,N-1);

%Bootstrapping
% 1)Characterizing uncertainty about the estimated impulse-response
% functions and fraction of forecast error variance by bootstrapping the
% VAR
NumberBootstraps=2000;
IRFs_Bootstrapped=zeros(size(IRFs_SimpleEstimate,1),size(IRFs_SimpleEstimate,2),NumberBootstraps);
FEVs_Bootstrapped=zeros(size(FractionsOfVariance_SimpleEstimate,1),size(FractionsOfVariance_SimpleEstimate,2),NumberBootstraps);

kk=1;
while kk<=NumberBootstraps
    kk
    YY=BootstrapSVARP([MU B],A0_SimpleEstimate,E',LagOrder,'Y',T,X);
    % Estimating the VAR:
    [b,varb,u,cov,vechv,varvechv]=varp(YY,LagOrder,'Y');
    % Checking that the VAR is stationary
    maxabsrootsvar=max(abs(varroots(LagOrder,N,b)));
    if maxabsrootsvar<1 & DetectImag(b)==0
        mu=b(:,1);
        b=b(:,2:size(b,2));
        a0=chol(cov)';
        e=a0\u;
        % impulse-response functions:
        irfs=GetIRFs(b,a0,N,LagOrder,Horizon,N-1);
        for xx=1:size(Index)
            irfs(:,Index(xx))=RollingMean([zeros(11,1); irfs(:,Index(xx))],12);
        end
        FEVs_Bootstrapped(:,:,kk)=GetFractionsOfVariance([mu b],a0,N,LagOrder,Horizon,N-1);
        IRFs_Bootstrapped(:,:,kk)=irfs;
        kk=kk+1;
    end
end

%Centering the bootstrapped IRFs
MedianIRFsBootstrapped=median(IRFs_Bootstrapped,3);

%Plotting the IFRs and confidence bands
kk=1;
while kk<=NumberBootstraps
    IRFs_Bootstrapped(:,:,kk)=IRFs_Bootstrapped(:,:,kk)-MedianIRFsBootstrapped+IRFs_SimpleEstimate;
    kk=kk+1;
end

HOR=(0:1:Horizon)';

for xx=1:N
    IRFSimpleEstimate=IRFs_SimpleEstimate(:,xx);
    PercentilesIRFs=ExtractPercentiles(sort(squeeze(IRFs_Bootstrapped(:,xx,:))'),[0.5 0.16 0.84 0.05 0.95]')';
    figure(2)
    subplot(2,5,xx)
    plot(HOR,IRFSimpleEstimate,'b',HOR,PercentilesIRFs(:,1),'k',HOR,PercentilesIRFs(:,2:3),'r:',HOR,PercentilesIRFs(:,4:5),'r',HOR,zeros(size(HOR)),'k:','LineWidth',2)
    xlim([0 Horizon])
    if xx==1
        title('LogDiffIndustrialProductionIndex')
    elseif xx==2
        title('LogPCEDeflator')
    elseif xx==3
        title('LogRealHousePrices')
    elseif xx==4
        title('UnemploymentRate ')
    elseif xx==5
        title('VacacyRate')
    elseif xx==6
        title('LogAggregateWeeklyHoursPerCapita')
    elseif xx==7
        title('LogHousingStarts')
    elseif xx==8
        title('FederalFundsRate')
    elseif xx==9
        title('LogStockPriceIndex')    
    end
end

% 2) The probability that the IFR of the unemployment rate is positive 4
%years after a monetary policy shock
%squeezing the Unemployment rate of the IRFS Bootstrapped (position 4)
ForYearsInTheFuture=4*12;
IRFsBootstrappedUnemploymentRate=squeeze(IRFs_Bootstrapped(ForYearsInTheFuture,4,:));
figure(3);
hist(IRFsBootstrappedUnemploymentRate, 200);
%Calculating the result in percent
(sum(IRFsBootstrappedUnemploymentRate>0)/NumberBootstraps)*100;
%The probability is 58.8%.

% 3) The probability that the IRF Housing starts is negative 1.5 years after
%a monetary shock (position 7)
OneAndAHalfYearsInTheFuture=1.5*12;
IRFsBootstrappedHousingStarts=squeeze(IRFs_Bootstrapped(OneAndAHalfYearsInTheFuture,7,:));
figure(4);
hist(IRFsBootstrappedHousingStarts, 200);
%Calculating the result in percent
(sum(IRFsBootstrappedHousingStarts<0)/NumberBootstraps)*100;
%The probability is 98%.

% 4) Performing a counterfactual simulation in which we re-run history by killing of
%monetary shocks
%Re-running history and comparing it with the original series
ReRunningHistory=X';
for tt=LagOrder+1:T
    Shocks=E(:,tt-LagOrder);
    ReRunningHistory(:,tt)=MU+B*vec(fliplr(ReRunningHistory(:,tt-LagOrder:tt-1)))+A0_SimpleEstimate*Shocks;
end
ReRunningHistory=ReRunningHistory';

for xx=1:size(Index)
    ReRunningHistory(:,Index(xx))=RollingMean([zeros(11,1);ReRunningHistory(:,Index(xx))],12);
end

Counterfactual_SimpleEstimate=ReRunningHistory;
CounterfactualMinusActual_SimpleEstimate=Counterfactual_SimpleEstimate-XX;
%seems to be very close to zero --> everything is fine

%Let's kill of the monetary shocks
KillingMPShock=X';
for tt=LagOrder+1:T
    Shocks=E(:,tt-LagOrder);
    Shocks(N-1)=0;
    KillingMPShock(:,tt)=MU+B*vec(fliplr(KillingMPShock(:,tt-LagOrder:tt-1)))+A0_SimpleEstimate*Shocks;
end

KillingMPShock=KillingMPShock';

for xx=1:size(Index)
    KillingMPShock(:,Index(xx))=RollingMean([zeros(11,1); KillingMPShock(:,Index(xx))],12);
end
CounterfactualMPShockKilledOffSimpleEstimate=KillingMPShock;
CounterfactualMPMinusActualSimpleEstimate=CounterfactualMPShockKilledOffSimpleEstimate-XX;

%Now let's characterize uncertainty about the re-runned history without
%the monetary policy shocks (using the results from the bootstrapped above)
Counterfactual_Bootstrapped=zeros(size(Counterfactual_SimpleEstimate,1),size(Counterfactual_SimpleEstimate,2),NumberBootstraps);
CounterfactualMinusActual_Bootstrapped=zeros(size(CounterfactualMinusActual_SimpleEstimate,1),size(CounterfactualMinusActual_SimpleEstimate,2),NumberBootstraps);

kk=1;
while kk<=NumberBootstraps
    kk
    YY=BootstrapSVARP([MU B],A0_SimpleEstimate,E',LagOrder,'Y',T,X);
    % Estimating the VAR:
    [b,varb,u,cov,vechv,varvechv]=varp(YY,LagOrder,'Y');
    % Here we check that the VAR is stationary, based on the roots of the
    % companion form of the VAR:
    maxabsrootsvar=max(abs(varroots(LagOrder,N,b)));
    if maxabsrootsvar<1 & DetectImag(b)==0
        mu=b(:,1);
        b=b(:,2:size(b,2));
        a0=chol(cov)';
        e=a0\u;
    
    %killing of the identified monetary policy shock
    yy=YY';
    for tt=LagOrder+1:T
        shocks=e(:,tt-LagOrder);
        shocks(N-1)=0;
        yy(:,tt)=mu+b*vec(fliplr(yy(:,tt-LagOrder:tt-1)))+a0*shocks;
    end
    yy=yy';
    for xx=1:size(Index)
        yy(:,Index(xx))=RollingMean([zeros(11,1); yy(:,Index(xx))],12);
    end
    counterfactual=yy;
    counterfactualminusactual=counterfactual-YY;
        
    Counterfactual_Bootstrapped(:,:,kk)=counterfactual;
    CounterfactualMinusActual_Bootstrapped(:,:,kk)=counterfactualminusactual;
        
    kk=kk+1;
    end
end

% Plotting the results
% Centering the bootstrapped Counterfactuals:
MedianCounterfactualBootstrapped=median(Counterfactual_Bootstrapped,3);

kk=1;
while kk<=NumberBootstraps
    Counterfactual_Bootstrapped(:,:,kk)=Counterfactual_Bootstrapped(:,:,kk)-MedianCounterfactualBootstrapped+CounterfactualMPShockKilledOffSimpleEstimate;
    kk=kk+1;
end

hor=(1:1:T)';

for xx=1:N
    Counterfact_SimpleEstimate=CounterfactualMPShockKilledOffSimpleEstimate(:,xx);
    PercentilesCounterfactual=ExtractPercentiles(sort(squeeze(Counterfactual_Bootstrapped(:,xx,:))'),[0.5 0.16 0.84 0.05 0.95]')';
    figure(5)
    subplot(2,5,xx)
    plot(hor,Counterfact_SimpleEstimate,'b',hor,PercentilesCounterfactual(:,1),'k',hor,PercentilesCounterfactual(:,2:3),'r:',hor,PercentilesCounterfactual(:,4:5),'r',hor,zeros(size(hor)),'k:','LineWidth',2)
    xlim([1 T])
    % X=[IndustrialProductionGrowth Inflation RealHousePricesGrowth UnemploymentRate VacancyRate LogAggWeeklyHoursPerC LogHousingStarts FedFundsRate StockPricesGrowth];
    if xx==1
        title('LogDiffIndustrialProductionIndex')
    elseif xx==2
        title('LogPCEDeflator')
    elseif xx==3
        title('LogRealHousePrices')
    elseif xx==4
        title('UnemploymentRate')
    elseif xx==5
        title('VacancyRate')
    elseif xx==6
        title('LogAggregateWeeklyHoursPerCapita')
    elseif xx==7
        title('LogHousingStarts')
    elseif xx==8
        title('FederalFundsRate')
    elseif xx==9
        title('LogStockPriceIndex')
    end
end


