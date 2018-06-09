%%Homework 2 - Jonas Schwery
clear
path(path,'C:\Users\Jonas Schwery\Desktop\Studium\EmpiricalMacro');
[X,B]=xlsread('C:\Users\Jonas Schwery\Desktop\Studium\EmpiricalMacro\DatasetHomeworkIII.xls','MonthlyData','A18:J343');
[T,N]=size(X);
Time=X(2:T,1);

%log differences
LogIndustrialProduction=log(X(:,2));
IndustrialProductionGrowth=((diff(LogIndustrialProduction)+1).^12-1)*100;
LogSEP=log(X(:,8));
StockPricesGrowth=((diff(LogSEP)+1).^12-1)*100;
LogPCEDeflator=log(X(:,3));
Inflation=((diff(LogPCEDeflator)+1).^12-1)*100;
CommercialandIndustrialLoans=log(X(:,7));
DiffcommercialandIndustrialLoans=((diff(CommercialandIndustrialLoans)+1).^12-1)*100;

%levels
UnemploymentRate=X(2:T,4);
VacancyRate=X(2:T,9)*100;
FEDFUNDS=X(2:T,5);
ExcessBondPremium=X(2:T,10);

%The Data Set
X=[IndustrialProductionGrowth StockPricesGrowth Inflation DiffcommercialandIndustrialLoans UnemploymentRate VacancyRate FEDFUNDS ExcessBondPremium];
[T,N]=size(X);

%First, we need to reorder the variables
X=[IndustrialProductionGrowth Inflation UnemploymentRate VacancyRate ExcessBondPremium FEDFUNDS StockPricesGrowth DiffcommercialandIndustrialLoans];
[T,N]=size(X);
PositionExcessBondPremium=5;

%Identifying highly volatile variables by plotting them
figure(1)
for xx=1:N
    subplot(2,4,xx)
    plot(Time,X(:,xx))
    if xx==1
        title('IndustrialProductionGrowth')
    elseif xx==2
        title('Inflation')
    elseif xx==3
        title('UnemploymentRate')
    elseif xx==4
        title('VacancyRate')
    elseif xx==5
        title('ExcessBondPremium')
    elseif xx==6
        title('FEDFUNDS')  
    elseif xx==7
        title('LogHousingStarts')
    elseif xx==8
        title('DiffcommercialandIndustrialLoans')
    end
end

%The variables on the positions 1, 2, 7 and 8 seem to be highly volatile
Index=[1 2 7 8]';

%Then we estimate a ver with 6 lags
LagOrder=6;
[B,VARB,U,V,vechV,VARvechV]=varp(X,LagOrder,'Y');

%and check for starionarity
[f,F,VAR]=Companion(B,V,N,LagOrder);
% Eigenvalue-eigenvector decomposition of F:
[VV,DD] = eig(F);
MaxAbsEigenvalue=max(abs(diag(DD)));
if MaxAbsEigenvalue<1
    disp('Estimated VAR is stationary')
else
    disp('Estimated VAR is non-stationary')
end

NumberOfDraws=2000; 
Stationary='Y'; 
Intercept='Y'; 
Trend=0; 
LagOrder=6;

% Estimating a Bayesian reduced-form VAR and the identifying a credit shock based on the restriction that,
% within the month, it only impacts upon the EBP, the Federal Funds
% rate, stock prices, and loans, whereas it does not have any impact o any other
% variable; and the impact at t=0 on the EBP is positive, whereas the impacts on the Federal
% Funds rate, stock prices, and loans are all negative.
[B,Sigma,U,TimeTrend]=BayesianVARp(X,LagOrder,Intercept,Trend,NumberOfDraws,Stationary);
if Trend>0
    CoefficientsOnTimeTrend=B(:,size(B,2),:);
    B=B(:,1:size(B,2)-1,:);
end
if Intercept=='Y'
    MU=B(:,1,:);
    B=B(:,2:size(B,2),:);
end

% Smoothing the original series
XX=X;
for xx=1:size(Index)
    XX(:,Index(xx))=RollingMean([zeros(11,1); XX(:,Index(xx))],12);
end

% The expected signs:
SIGNS=[1 -1 -1 -1]';

% The maximum number of randomly drawn matrices
MaxNumberTry=2000; 

%The structural impact matrices
AA00=zeros(N,N,NumberOfDraws); 

% This is an indicator which is equal to 1 if we draw a matrix which
% satisfies the signs we want to impose (=1), and it is equal to zero otherwise
A0IN=zeros(NumberOfDraws,1);

% We will use this to check if A0*A0'-var is close to zero (up to machine
% precision)
test=zeros(NumberOfDraws,1);  % to check if A0*A0'-var is close to zero 
VerySmallValue = 0.0000001;

ss=1;
while ss<=NumberOfDraws
    ss
    % Here we get, for each single draw from 
    % the posterior, the elements of the VAR:
    mu=MU(:,:,ss);
    bb=B(:,:,ss);
    var=Sigma(:,:,ss);
    % Computing the initial estimate of the A0 matrix via Cholesky, so that we impose the
    % zeros:
    a0start=chol(var)'; 
    
    IND=0; 
      
    xx=1;
    while IND==0 & xx<MaxNumberTry 
        % First, we get a randomly draw (N-PositionExcessBondPremium+1)x(N-PositionExcessBondPremium+1) matrix from the Normal(0,1) distribution:
        K=NormalRandomVector(0,1,N-PositionExcessBondPremium+1,N-PositionExcessBondPremium+1);
        % Then, from this we construct a random rotation matrix:
        [Q,R]=qr(K);
        for i=1:N-PositionExcessBondPremium+1
            if R(i,i)<0
                Q(:,i)=-Q(:,i);
            end
        end
        RotationMatrix=eye(N);  
        RotationMatrix(PositionExcessBondPremium:N,PositionExcessBondPremium:N)=Q';
        % Candidate impact matrix A0:
        A0=a0start*RotationMatrix;
        % Check out that: has to be 0 up to machine precision A0*A0'-var
        % is negligible. 
        if A0*A0'-var < VerySmallValue
            test(ss,1)=1;
        end 
        % Checking the signs:
        SignImpactAt0=sign(A0(PositionExcessBondPremium:N,PositionExcessBondPremium));
        % If the signs are fine, we terminate the search:
        if any(SignImpactAt0-SIGNS)==0
            IND=1;
            break
        elseif any(-SignImpactAt0-SIGNS)==0
            A0=-A0;
            IND=1;
            break
        else
        end
        xx=xx+1;
    end
    if IND==1 % Storing the results if we are successfull
        AA00(:,:,ss)=A0;
        A0IN(ss)=1;
    end
    ss=ss+1;
end
BBBB=B;
SSSS=Sigma;
UUUU=U;
if Intercept=='Y'
    MUMU=MU;
else
    MUMU=-9999;
end
Indices=ExtractIndices(A0IN);
NN=length(Indices);
NumberOfSuccessfulDraws=NN;
Percentiles=fix(NN*[0.5 0.16 0.84 0.05 0.95]');
N=size(AA00,1);

BBBB=BBBB(:,:,Indices);
SSSS=SSSS(:,:,Indices);
UUUU=UUUU(:,:,Indices);
MUMU=MUMU(:,:,Indices);
AA00=AA00(:,:,Indices);

% Test if A0*t(A0)-var is very very close to zero (up to machine percision)
sum(test,1)== NumberOfDraws;
% Yes, seems to be fine --> logical =1 (true)
  
% 1) Now we comptue the IRFs
NN=size(MUMU,3);
Percentiles=fix(NN*[0.5 0.16 0.84 0.05 0.95]');
% The horizon for the IRFs, in months (10 years)
Horizon=10*12; 

IRFsCreditShock=zeros(Horizon+1,N,NN);

ss=1;
while ss<=NN
    ss
    % The matrices B and A0 for draw ss:
    bb=BBBB(:,:,ss);
    a0=AA00(:,:,ss);
    % Here we get the IRFs:
    irfs=GetIRFs(bb,a0,N,LagOrder,Horizon,PositionExcessBondPremium);
    % Smoothing the variables which are especially volatile (1,2,7 and 8):
    for xx=1:size(Index)
        irfs(:,Index(xx))=RollingMean([zeros(11,1); irfs(:,Index(xx))],12);
    end
    IRFsCreditShock(:,:,ss)=irfs;
    %
    ss=ss+1;
end

HOR=(0:1:Horizon)';

Variable=1;
while Variable<=N
    IRF=sort(squeeze(IRFsCreditShock(:,Variable,:)),2);
    IRF=IRF(:,Percentiles);
    figure(2)
    subplot(2,4,Variable)
    plot(HOR,IRF(:,2:3),'r:',HOR,IRF(:,4:5),'r',HOR,IRF(:,1),'k',HOR,zeros(size(HOR)),'b:','LineWidth',2)
    xlim([0 Horizon])
    if Variable==1
        title('IndustrialProductionGrowth')
    elseif Variable==2
        title('Inflation')
    elseif Variable==3
        title('UnemploymentRate')
    elseif Variable==4
        title('VacancyRate')
    elseif Variable==5
        title('ExcessBondPremium')
    elseif Variable==6
        title('FEDFUNDS')
    elseif Variable==7
        title('StockPricesGrowth')
    elseif Variable==8
        title('DiffcommercialandIndustrialLoans')    
    else
    end
    Variable=Variable+1;
end

% And now the variance decomposition
FractionsOfVarianceCreditShock=zeros(Horizon+1,N,NN);

ss=1;
while ss<=NN
    ss
    mu=MUMU(:,:,ss);
    bb=BBBB(:,:,ss);
    a0=AA00(:,:,ss);
    
    FractionsOfVarianceCreditShock(:,:,ss)=GetFractionsOfVariance([mu bb],a0,N,LagOrder,Horizon,PositionExcessBondPremium);
    
    ss=ss+1;
end

Variable=1;
while Variable<=N
    FRAC=sort(squeeze(FractionsOfVarianceCreditShock(:,Variable,:)),2);
    FRAC=FRAC(:,Percentiles);
    figure(3)
    subplot(2,4,Variable)
    plot(HOR,FRAC(:,2:3),'r:',HOR,FRAC(:,4:5),'r',HOR,FRAC(:,1),'k','LineWidth',2)
    axis([0 Horizon 0 1])
    if Variable==1
        title('IndustrialProductionGrowth')
    elseif Variable==2
        title('Inflation')
    elseif Variable==3
        title('UnemploymentRate')
    elseif Variable==4
        title('VacancyRate')
    elseif Variable==5
        title('ExcessBondPremium')
    elseif Variable==6
        title('FEDFUNDS')
    elseif Variable==7
        title('StockPricesGrowth')
    elseif Variable==8
        title('DiffcommercialandIndustrialLoans') 
    else
    end
    Variable=Variable+1;
end

% 2) Performing and plotting a counterfactual simulation in which we re-run
%history by ‘killing off’ the identified credit shocks.

% The position of the credit shocks (5:ExcessBondPremium)
ShocksEBP=5;
CounterfactualXX=zeros(T,N,NN);
CounterfactualMinusActualXX=zeros(T,N,NN);
check=zeros(NN,1);

ss=1;
while ss<=NN
    A0=AA00(:,:,ss);
    StructuralShocks(:,:,ss)=A0\UUUU(:,:,ss);
    ss=ss+1;
end

ss=1;
while ss<=NN
    ss
    mu=MUMU(:,:,ss);
    bb=BBBB(:,:,ss);
    a0=AA00(:,:,ss);
    
    CounterfactualO=X';
    
    tt=LagOrder+1;
    while tt<=T
        Shocks=StructuralShocks(:,tt-LagOrder,ss);
        % trying to recover the original series by not killing of any shock
        CounterfactualO(:,tt)=mu+bb*vec(fliplr(CounterfactualO(:,tt-LagOrder:tt-1)))+a0*Shocks;
        tt=tt+1;
    end
    % checking if the difference between counterfactual and original series is
    % very close to zero
    Diff=CounterfactualO'-X;
    hh=1;
    ww=1;
    for hh=1:size(CounterfactualO',1)
        for ww=1:size(CounterfactualO',2)
            if Diff(hh,ww)< VerySmallValue
                 check(ss,1)=1;
            end 
        end 
    end 
    % check is 0, if the condition is not satisfied; otherwise it is 1
    % (--> see the test sum(check)==NN after the loop)
    CounterfactualX=X';
    tt=LagOrder+1;
    while tt<=T
        % These are the identified shocks for month tt
        % and draw ss from the posterior dostribution:
        Shocks=StructuralShocks(:,tt-LagOrder,ss);
        % Here we kill off the shock we are not interested in:
        Shocks(ShocksEBP)=0;
        % Here we re-run history after having killed off the relevant shock:
        CounterfactualX(:,tt)=mu+bb*vec(fliplr(CounterfactualX(:,tt-LagOrder:tt-1)))+a0*Shocks;
        tt=tt+1;
    end
    %Smoothing the variables which seem to be very volatile
    CounterfactualX=CounterfactualX';
    for xx=1:size(Index)
       CounterfactualX(:,Index(xx))=RollingMean([zeros(11,1); CounterfactualX(:,Index(xx))],12);
    end
    
    % Here we store it:
    CounterfactualXX(:,:,ss)=CounterfactualX;
    CounterfactualXXMinusActualX(:,:,ss)=CounterfactualXX(:,:,ss)-XX;
    ss=ss+1;
end

%Test
sum(check)== NN;
%logical of 1 --> everything is fine

% Getting the percentiles of the counterfactual series:
CounterfactualSeries=sort(CounterfactualXX,3);
PercentilesCounterfactualSeries=CounterfactualSeries(:,:,Percentiles);

% Getting the percentiles of the difference between counterfactual and actual series:
CounterfactualMinusActualSeries=sort(CounterfactualXXMinusActualX,3);
PercentilesCounterfactualMinusActualSeries=CounterfactualMinusActualSeries(:,:,Percentiles);

% Plotting the results
xx=1;
while xx<=N
    Counterfactual=squeeze(PercentilesCounterfactualSeries(:,xx,:));
    figure(4)
    subplot(2,4,xx)
    plot(Time,Counterfactual(:,1),'k',Time,Counterfactual(:,4:5),'r',Time,X(:,xx),'b',Time,zeros(size(Time)),'b:')
    xlim([Time(1) Time(T)])
    if xx==1
        title('IndustrialProductionGrowth')
    elseif xx==2
        title('Inflation')
    elseif xx==3
        title('UnemploymentRate')
    elseif xx==4
        title('VacancyRate')
    elseif xx==5
        title('ExcessBondPremium')
    elseif xx==6
        title('FEDFUNDS')
    elseif xx==7
        title('StockPricesGrowth')
    elseif xx==8
        title('DiffcommercialandIndustrialLoans')
    else
    end
    
    CounterfactualMinusActual=squeeze(PercentilesCounterfactualMinusActualSeries(:,xx,:));
    
    figure(5)
    subplot(2,4,xx)
    plot(Time,CounterfactualMinusActual(:,1),'k',Time,CounterfactualMinusActual(:,4:5),'r',Time,zeros(size(Time)),'b:')
    xlim([Time(1) Time(T)])
   if xx==1
        title('IndustrialProductionGrowth')
    elseif xx==2
        title('Inflation')
    elseif xx==3
        title('UnemploymentRate')
    elseif xx==4
        title('VacancyRate')
    elseif xx==5
        title('ExcessBondPremium')
    elseif xx==6
        title('FEDFUNDS')
    elseif xx==7
        title('StockPricesGrowth')
    elseif xx==8
        title('DiffcommercialandIndustrialLoans') 
    else
    end
    %
    xx=xx+1;
end

% 3) Comptuing the probability that, 12 months after a credit shock, the IRFs of
% industrial production and the vacancy rate are both negative, and the IRF of
% the unemployment rate is positive

% Defining the Horizon of 12 Months and squeezing the corresponding
% variables
Horizon1=12;
IrfsCreditShockOfIndustrialProductionGrowthHorizonOfMonths=squeeze(IRFsCreditShock(Horizon1,1,:));
IrfsCreditShockOfVacancyRateHorizonOfMonths=squeeze(IRFsCreditShock(Horizon1,4,:));
IrfsCreditShockOfUnemploymentRateHorizonOfMonths=squeeze(IRFsCreditShock(Horizon1,3,:));

% Lets compute the probability, that Industrial Production Growth and
% Vacancy Rate are positive and the unemployment Rate is negative
(sum(IrfsCreditShockOfIndustrialProductionGrowthHorizonOfMonths<0 & IrfsCreditShockOfVacancyRateHorizonOfMonths<0 & IrfsCreditShockOfUnemploymentRateHorizonOfMonths>0)/NN)*100;
% The probability is 89.9%. 

% 4)Computing the probability that, 10 years ahead, credit shocks explain at least 10
% per cent of the forecast error variance of both the unemployment rate, and the
% vacancy rate

% Defining the Horizon 10 years ahead (12*10=120)
Horizon2=10*12;
FractionsOfVarianceCreditShockUnemploymentRateHor2=squeeze(FractionsOfVarianceCreditShock(Horizon2,3,:));
FractionsOfVarianceCreditShockVacancyRateHor2=squeeze(FractionsOfVarianceCreditShock(Horizon2,4,:));

% And now computing the probability that 10 years ahead, credit shocks explain at least 10
% per cent of the forecast error variance of both the unemployment rate, and the
% vacancy rate
(sum(FractionsOfVarianceCreditShockUnemploymentRateHor2>0.1 & FractionsOfVarianceCreditShockVacancyRateHor2>0.1)/NN)*100;
% The probability is 26.5%.


