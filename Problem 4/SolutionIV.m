%% Empirical Macro - Homework IV Jonas Schwery
clear
path(path,'C:\Users\Jonas Schwery\Desktop\Studium\EmpiricalMacro');
[X,B]=xlsread('C:\Users\Jonas Schwery\Desktop\Studium\EmpiricalMacro\DatasetHomeworkIV.xls','Dataset','A17:F238');
[T,N]=size(X);
Time=X(2:T,1);
% The log difference of the real gross domestic product
LogRealGDP=log(X(:,2));
RealGDPGrowth=diff(LogRealGDP);
%Variables that enter in levels
UnemploymentRate=X(2:T,3)/100;
ThreeMonthTreasuryBill=X(2:T,4)/400;
logGDPDeflator=log(X(:,5));
GDPDeflatorInflation=diff(logGDPDeflator);
TenYearRate=X(2:T,6)/400;
Spread=TenYearRate-ThreeMonthTreasuryBill;

X=[RealGDPGrowth UnemploymentRate ThreeMonthTreasuryBill GDPDeflatorInflation Spread];
[T,N]=size(X);

%Defining some things to estimate the Bayesian reduced-form Var
NumberOfDraws=2000; 
MaxNumberTry=2000;
Stationary='Y'; 
Intercept='Y'; 
Trend=0; 
LagOrder=4;

% Estimating a Bayesian reduced-form VAR based on 2000 draws as defined
% before
[B,Sigma,U,TimeTrend]=BayesianVARp(X,LagOrder,Intercept,Trend,NumberOfDraws,Stationary);
if Trend>0
    CoefficientsOnTimeTrend=B(:,size(B,2),:);
    B=B(:,1:size(B,2)-1,:);
end
if Intercept=='Y'
    MU=B(:,1,:);
    B=B(:,2:size(B,2),:);
end

%Now we define the signs we want to impose from the Excel Sheet
%X=[RealGDPGrowth UnemploymentRate ThreeMonthTreasuryBill GDPDeflatorInflation Spread]
%Shocks:			
%		                     Technology	Monetary	Taste	Markup
%					
%Real Gross Domestic Product      1	      -1	      1	      -1
%Unemployment Rate		          1	       1	     -1	       1					
%3-Month Treasury Bill		     -1	       1	      1	       1
%GDP deflator inflation		     -1	      -1	      1	       1
%The vector signs is going to be:
SIGNS=vec([ 1 -1 1 -1; 1 1 -1 1; -1 1 1 1; -1 -1 1 1]);
%Notice: 1th shock: technology shock, 2th shock: monetary shock, 3th shock:
%taste shock and 4th shock: mark up shock..
AA00=zeros(N,N,NumberOfDraws); 
% This is an indicator which is equal to 1 if we draw a matrix which
% satisfies the signs we want to impose, and it is equal to zero otherwise: 
A0IN=zeros(NumberOfDraws,1);
%
ss=1;
while ss<=NumberOfDraws
    ss
    % Getting the elements of the Var, for each single draw from the
    % posterior
    mu=MU(:,:,ss);
    bb=B(:,:,ss);
    var=Sigma(:,:,ss);
    % Getting the matrix C = MyInverse[I - B(1) - B(2) - ... - B(p)]
    % where B(1), B(2), ... , B(p) are the coefficient matrices of the VAR
    if LagOrder>1
        c=MyInverse(eye(N)-sum(reshape(bb,N,N,LagOrder),3));
    else
        c=MyInverse(eye(N)-bb);
    end
    % the starting estimate of the a0 matrix:
    [vv,dd]=eig(var);
    a0start=vv*dd.^0.5;
    % The corresponding starting estimate of the long-run impact matrix:
    longrunimpactstart=c*a0start;
    % Notice that longrunimpactstart has no zero entries, that is: all
    % shocks have a long-run impact on the relevant variables.
    %
    % Now we want to impose the restriction that there is only 1 shock that 
    % has a long-run impact on GDP  (which is the first variable).
    % This means that, in the first row of the matrix of long-run impacts, it has to be the case that
    % (1) the entries in the first column is non-zero, and
    % (2) the entries in the last 4 columns are equal to zero.
    % 
    % Let's partition longrunimpactstart into four blocks as:
    %
    % longrunimpactstart = [a b]
    %                      [c d]
    %
    % with a being (1x1) and b being (1X4):
    a=longrunimpactstart(1:1,1:1);
    b=longrunimpactstart(1:1,2:N);
    % Then, we want to set b=zeros(1,4), and we do it in this way:
    %
    BB=(a\b)';
    AA=chol(MyInverse(eye(1)+BB'*BB));
    CC=chol(MyInverse(eye(N-1)+BB*BB'));
    RR=[AA AA*BB'; -CC*BB CC]';
    %
    % Checking that RR truly is a rotation matrix, that is, both RR*RR' and
    % RR'*RR are equal to the identity matrix (I did the check after the
    % loop...
    %
    longrunimpact=longrunimpactstart*RR;
    % Check that longrunimpact has precisely the pattern of zeros we want
    % to impose ...
    %
    % Then we get the implied structural impact matrix a0, by exploiting the fact that:
    % longrunimpact = c * a0
    % which automatically implies that
    % a0 = c\longrunimpact
    a0=c\longrunimpact;
    % Given the structure of the matrix longrunimpact, the last 4 shocks
    % have no long-run impact on GDP.
    % Therefore, the shock we are looking for is the remaining one, that
    % is: the first shock. In order to separate the shock, now we randomly
    % rotate the shock in order to impose the sign restrictions.
    a0start=a0;
    %
    % This is an index which keeps track of the fact that we have not found the
    % matrix we are looking for yet: when we find it, we set the index to
    % 1, and the loop automaically stops ...
    IND=0; 
    %   
    xx=1;
    while IND==0 & xx<MaxNumberTry
        % Drawing a random rotation matrix:
        %
        % First, we get a randomly draw (4x4) matrix from the Normal(0,1) distribution:
        K=NormalRandomVector(0,1,4,4);
        % Constructing a a random rotation matrix for this
        [Q,R]=qr(K);
        for i=1:4
            if R(i,i)<0
                Q(:,i)=-Q(:,i);
            end
        end
        RotationMatrix=eye(N);
        RotationMatrix(2:5,2:5)=Q';
        % Candidate impact matrix A0:
        A0=a0start*RotationMatrix;
        % Checking out that:
        % A0*A0'-var
        % is negligible. This is indeed the case ... (once again i made the
        % check after the loop...)
        % Checking the signs:
        SignImpactAt0=vec(sign(A0(1:4,2:5)));
        % Terminating the search if the signs are fine
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
    if IND==1 %Storing the results if we are successfull
        AA00(:,:,ss)=A0;
        A0IN(ss)=1;
    end
    ss=ss+1;
end

% Checking that both RR*RR' and RR'*RR are equal to the identity matrix
RR*RR';
RR'*RR;

%Checking out that: A0*A0'-var is negligible.
A0*A0'-var;
%This is indeed the case.

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

% Calculating the fraction of successfull draws
NumberOfSuccessfulDraws/NumberOfDraws;

%Extracting the successfull draws
BBBB=BBBB(:,:,Indices);
SSSS=SSSS(:,:,Indices);
UUUU=UUUU(:,:,Indices);
MUMU=MUMU(:,:,Indices);
AA00=AA00(:,:,Indices);

%% Now we compute the IRFS
% The horizon for the IRFs, in quarters (20 years)
Horizon=20*4; 

% Defining some empty matrices which we fill in later on
IRFsPermanentShock=zeros(Horizon+1,N,NN);
IRFsTechnologyShock=zeros(Horizon+1,N,NN);
IRFsMonetaryShock=zeros(Horizon+1,N,NN);
IRFsTasteShock=zeros(Horizon+1,N,NN);
IRFsMarkupShock=zeros(Horizon+1,N,NN);
%
ss=1;
while ss<=NN
    ss
    % The matrices B and A0 for draw ss:
    bb=BBBB(:,:,ss);
    a0=AA00(:,:,ss);
    irfs=GetIRFs(bb,a0,N,LagOrder,Horizon,1); %permanent
    % Here we cumulate the IRFs for RealGDPGrowth
    irfs(:,1)=cumsum(irfs(:,1));
    irfs=irfs*sign(irfs(1+Horizon,1));
    IRFsPermanentShock(:,:,ss)=irfs;
    % Here we cumulate the IRFs for RealGDPGrowth
    irfs=GetIRFs(bb,a0,N,LagOrder,Horizon,2);
    irfs(:,1)=cumsum(irfs(:,1));
    IRFsTechnologyShock(:,:,ss)=irfs;
    % Here we cumulate the IRFs for RealGDPGrowth
    irfs=GetIRFs(bb,a0,N,LagOrder,Horizon,3);
    irfs(:,1)=cumsum(irfs(:,1));
    IRFsMonetaryShock(:,:,ss)=irfs;
    % Here we cumulate the IRFs for RealGDPGrowth
    irfs=GetIRFs(bb,a0,N,LagOrder,Horizon,4);
    irfs(:,1)=cumsum(irfs(:,1));
    IRFsTasteShock(:,:,ss)=irfs;
    % Here we cumulate the IRFs for RealGDPGrowth
    irfs=GetIRFs(bb,a0,N,LagOrder,Horizon,5);
    irfs(:,1)=cumsum(irfs(:,1));
    IRFsMarkupShock(:,:,ss)=irfs;
    ss=ss+1;
end

HOR=(0:1:Horizon)';

Variable=1;
while Variable<=N
    IRF=sort(squeeze(IRFsPermanentShock(:,Variable,:)),2);
    IRF=IRF(:,Percentiles);
    figure(2)
    subplot(5,N,Variable)
    plot(HOR,IRF(:,2:3),'r:',HOR,IRF(:,4:5),'r',HOR,IRF(:,1),'k',HOR,zeros(size(HOR)),'b:','LineWidth',2)
    xlim([0 Horizon])
    if Variable==1
        title('RealGDPGrowth')
        ylabel('Permanent shock')
    elseif Variable==2
        title('UnemploymentRate')
    elseif Variable==3
        title('ThreeMonthTreasuryBill')
    elseif Variable==4
        title('GDPDeflatorInflation')
    elseif Variable==5
        title('Spread')
    else
    end
    %X=[RealGDPGrowth UnemploymentRate ThreeMonthTreasuryBill GDPDeflatorInflation Spread];
    
    IRF=sort(squeeze(IRFsTechnologyShock(:,Variable,:)),2);
    IRF=IRF(:,Percentiles);
    figure(2)
    subplot(5,N,N+Variable)
    plot(HOR,IRF(:,2:3),'r:',HOR,IRF(:,4:5),'r',HOR,IRF(:,1),'k',HOR,zeros(size(HOR)),'b:','LineWidth',2)
    xlim([0 Horizon])
    if Variable==1
        ylabel('Technology Shock')
    else
    end
    
    IRF=sort(squeeze(IRFsMonetaryShock(:,Variable,:)),2);
    IRF=IRF(:,Percentiles);
    figure(2)
    subplot(5,N,N*2+Variable)
    plot(HOR,IRF(:,2:3),'r:',HOR,IRF(:,4:5),'r',HOR,IRF(:,1),'k',HOR,zeros(size(HOR)),'b:','LineWidth',2)
    xlim([0 Horizon])
    if Variable==1
        ylabel('Monetary shock')
    else
    end
    
    IRF=sort(squeeze(IRFsTasteShock(:,Variable,:)),2);
    IRF=IRF(:,Percentiles);
    figure(2)
    subplot(5,N,N*3+Variable)
    plot(HOR,IRF(:,2:3),'r:',HOR,IRF(:,4:5),'r',HOR,IRF(:,1),'k',HOR,zeros(size(HOR)),'b:','LineWidth',2)
    xlim([0 Horizon])
    if Variable==1
        ylabel('Taste shock')
    else
    end
    
    IRF=sort(squeeze(IRFsMarkupShock(:,Variable,:)),2);
    IRF=IRF(:,Percentiles);
    figure(2)
    subplot(5,N,N*4+Variable)
    plot(HOR,IRF(:,2:3),'r:',HOR,IRF(:,4:5),'r',HOR,IRF(:,1),'k',HOR,zeros(size(HOR)),'b:','LineWidth',2)
    xlim([0 Horizon])
    if Variable==1
        ylabel('Markup shock')
    else
    end
    
    Variable=Variable+1;
end

%% Computing the variance decomposition
FractionsOfVariancePermanentShock=zeros(Horizon+1,N,NN);
FractionsOfVarianceTechnologyShock=zeros(Horizon+1,N,NN);
FractionsOfVarianceMonetaryShock=zeros(Horizon+1,N,NN);
FractionsOfVarianceTasteShock=zeros(Horizon+1,N,NN);
FractionsOfVarianceMarkupShock=zeros(Horizon+1,N,NN);
%
ss=1;
while ss<=NN
    ss
    mu=MUMU(:,:,ss);
    bb=BBBB(:,:,ss);
    a0=AA00(:,:,ss);
    
    FractionsOfVariancePermanentShock(:,:,ss)=GetFractionsOfVariance([mu bb],a0,N,LagOrder,Horizon,1);
    FractionsOfVarianceTechnologyShock(:,:,ss)=GetFractionsOfVariance([mu bb],a0,N,LagOrder,Horizon,2);
    FractionsOfVarianceMonetaryShock(:,:,ss)=GetFractionsOfVariance([mu bb],a0,N,LagOrder,Horizon,3);
    FractionsOfVarianceTasteShock(:,:,ss)=GetFractionsOfVariance([mu bb],a0,N,LagOrder,Horizon,4);
    FractionsOfVarianceMarkupShock(:,:,ss)=GetFractionsOfVariance([mu bb],a0,N,LagOrder,Horizon,5);
    
    ss=ss+1;
end

Variable=1;
while Variable<=N
    
    FRAC=sort(squeeze(FractionsOfVariancePermanentShock(:,Variable,:)),2);
    FRAC=FRAC(:,Percentiles);
    figure(3)
    subplot(5,N,Variable)
    plot(HOR,FRAC(:,2:3),'r:',HOR,FRAC(:,4:5),'r',HOR,FRAC(:,1),'k','LineWidth',2)
    axis([0 Horizon 0 1])
    if Variable==1
        title('RealGDPGrowth')
        ylabel('Permanent shock')
    elseif Variable==2
        title('UnemploymentRate')
    elseif Variable==3
        title('ThreeMonthTreasuryBill')
    elseif Variable==4
        title('GDPDeflatorInflation')
    elseif Variable==5
        title('Spread')
    else
    end
    %X=[RealGDPGrowth UnemploymentRate ThreeMonthTreasuryBill GDPDeflatorInflation Spread];
    FRAC=sort(squeeze(FractionsOfVarianceTechnologyShock(:,Variable,:)),2);
    FRAC=FRAC(:,Percentiles);
    figure(3)
    subplot(5,N,N+Variable)
    plot(HOR,FRAC(:,2:3),'r:',HOR,FRAC(:,4:5),'r',HOR,FRAC(:,1),'k','LineWidth',2)
    axis([0 Horizon 0 1])
    if Variable==1
        ylabel('Technology shock')
    else
    end
    FRAC=sort(squeeze(FractionsOfVarianceMonetaryShock(:,Variable,:)),2);
    FRAC=FRAC(:,Percentiles);
    figure(3)
    subplot(5,N,N*2+Variable)
    plot(HOR,FRAC(:,2:3),'r:',HOR,FRAC(:,4:5),'r',HOR,FRAC(:,1),'k','LineWidth',2)
    axis([0 Horizon 0 1])
    if Variable==1
        ylabel('Monetary shock')
    else
    end
    
    FRAC=sort(squeeze(FractionsOfVarianceTasteShock(:,Variable,:)),2);
    FRAC=FRAC(:,Percentiles);
    figure(3)
    subplot(5,N,N*3+Variable)
    plot(HOR,FRAC(:,2:3),'r:',HOR,FRAC(:,4:5),'r',HOR,FRAC(:,1),'k','LineWidth',2)
    axis([0 Horizon 0 1])
    if Variable==1
        ylabel('Taste shock')
    else
    end
   
    FRAC=sort(squeeze(FractionsOfVarianceMarkupShock(:,Variable,:)),2);
    FRAC=FRAC(:,Percentiles);
    figure(3)
    subplot(5,N,N*4+Variable)
    plot(HOR,FRAC(:,2:3),'r:',HOR,FRAC(:,4:5),'r',HOR,FRAC(:,1),'k','LineWidth',2)
    axis([0 Horizon 0 1])
    if Variable==1
        ylabel('Markup shock')
    else
    end
    
    Variable=Variable+1;
end

%% Performing and plotting counterfactual simulations in which we re-run history by ‘killing off’ the identified shocks, ONE AT A TIME!
% We have to kill off the shocks ONE AT A TIME, so let's begin with the
% first one..
StructuralShocks=zeros(N,T-LagOrder,NN);
IdentifiedShocks=1;
CounterfactualX1=zeros(T,N,NN);
CounterfactualMinusActualX1=zeros(T,N,NN);
check=zeros(NN,1);  

% to check if A0*A0'-var is close to zero we define a small value
VerySmallValue = 0.0000001;

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
        Shocks=StructuralShocks(:,tt-LagOrder,ss);
        % Here we kill off the shock we are not interested in:
        Shocks(IdentifiedShocks)=0;
        % Here we re-run history after having killed off the relevant shock:
        CounterfactualX(:,tt)=mu+bb*vec(fliplr(CounterfactualX(:,tt-LagOrder:tt-1)))+a0*Shocks;
        tt=tt+1;
    end
    CounterfactualX1(:,:,ss)=CounterfactualX';
    CounterfactualMinusActualX1(:,:,ss)=CounterfactualX1(:,:,ss)-X;
    ss=ss+1;
end

%Test
sum(check)== NN;
%logical of 1 --> everything is fine

% Getting the percentiles of the counterfactual series:
CounterfactualSeries=sort(CounterfactualX1,3);
PercentilesCounterfactualSeries=CounterfactualSeries(:,:,Percentiles);

% Getting the percentiles of the difference between counterfactual and actual series:
CounterfactualMinusActualSeries=sort(CounterfactualMinusActualX1,3);
PercentilesCounterfactualMinusActualSeries=CounterfactualMinusActualSeries(:,:,Percentiles);

% Plotting the results
xx=1;
while xx<=N
    Counterfactual=squeeze(PercentilesCounterfactualSeries(:,xx,:));
    figure(4)
    subplot(2,5,xx)
    plot(Time,Counterfactual(:,1),'k',Time,Counterfactual(:,4:5),'r',Time,X(:,xx),'b',Time,zeros(size(Time)),'b:')
    xlim([Time(1) Time(T)])
    if xx==1
        title('RealGDPGrowth')
    elseif xx==2
        title('UnemploymentRate')
    elseif xx==3
        title('ThreeMonthTreasuryBill')
    elseif xx==4
        title('GDPDeflatorInflation')
    elseif xx==5
        title('Spread')
     %X=[RealGDPGrowth UnemploymentRate ThreeMonthTreasuryBill GDPDeflatorInflation Spread];    
    else
    end
    
    CounterfactualMinusActual=squeeze(PercentilesCounterfactualMinusActualSeries(:,xx,:));
    
    figure(4)
    subplot(2,5,xx+5)
    plot(Time,CounterfactualMinusActual(:,1),'k',Time,CounterfactualMinusActual(:,4:5), 'r',Time,CounterfactualMinusActual(:,4:5),'r:',Time,zeros(size(Time)),'b:')
    xlim([Time(1) Time(T)])
   if xx==1
        title('RealGDPGrowth')
    elseif xx==2
        title('UnemploymentRate')
    elseif xx==3
        title('ThreeMonthTreasuryBill')
    elseif xx==4
        title('GDPDeflatorInflation')
    elseif xx==5
        title('Spread') 
    else
    end
    %
    xx=xx+1;
end

%% Now the same for Shock 2
IdentifiedShocks=2;
CounterfactualX2=zeros(T,N,NN);
CounterfactualMinusActualX2=zeros(T,N,NN);
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
        Shocks=StructuralShocks(:,tt-LagOrder,ss);
        % Here we kill off the shock we are not interested in:
        Shocks(IdentifiedShocks)=0;
        % Here we re-run history after having killed off the relevant shock:
        CounterfactualX(:,tt)=mu+bb*vec(fliplr(CounterfactualX(:,tt-LagOrder:tt-1)))+a0*Shocks;
        tt=tt+1;
    end
    CounterfactualX2(:,:,ss)=CounterfactualX';
    CounterfactualMinusActualX2(:,:,ss)=CounterfactualX2(:,:,ss)-X;
  
    ss=ss+1;
end

%Test
sum(check)== NN;
%logical of 1 --> everything is fine

% Getting the percentiles of the counterfactual series:
CounterfactualSeries=sort(CounterfactualX2,3);
PercentilesCounterfactualSeries=CounterfactualSeries(:,:,Percentiles);

% Getting the percentiles of the difference between counterfactual and actual series:
CounterfactualMinusActualSeries=sort(CounterfactualMinusActualX2,3);
PercentilesCounterfactualMinusActualSeries=CounterfactualMinusActualSeries(:,:,Percentiles);

% Plotting the results
xx=1;
while xx<=N
    Counterfactual=squeeze(PercentilesCounterfactualSeries(:,xx,:));
    figure(5)
    subplot(2,5,xx)
    plot(Time,Counterfactual(:,1),'k',Time,Counterfactual(:,4:5),'r',Time,X(:,xx),'b',Time,zeros(size(Time)),'b:')
    xlim([Time(1) Time(T)])
    if xx==1
        title('RealGDPGrowth')
    elseif xx==2
        title('UnemploymentRate')
    elseif xx==3
        title('ThreeMonthTreasuryBill')
    elseif xx==4
        title('GDPDeflatorInflation')
    elseif xx==5
        title('Spread')
     %X=[RealGDPGrowth UnemploymentRate ThreeMonthTreasuryBill GDPDeflatorInflation Spread];    
    else
    end
    
    CounterfactualMinusActual=squeeze(PercentilesCounterfactualMinusActualSeries(:,xx,:));
    
    figure(5)
    subplot(2,5,xx+5)
    plot(Time,CounterfactualMinusActual(:,1),'k',Time,CounterfactualMinusActual(:,4:5),'r',Time,zeros(size(Time)),'b:')
    xlim([Time(1) Time(T)])
   if xx==1
        title('RealGDPGrowth')
    elseif xx==2
        title('UnemploymentRate')
    elseif xx==3
        title('ThreeMonthTreasuryBill')
    elseif xx==4
        title('GDPDeflatorInflation')
    elseif xx==5
        title('Spread') 
    else
    end
    %
    xx=xx+1;
end


%% Now for the 3th shock
IdentifiedShocks=3;
CounterfactualX3=zeros(T,N,NN);
CounterfactualMinusActualX3=zeros(T,N,NN);
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
        Shocks=StructuralShocks(:,tt-LagOrder,ss);
        % Here we kill off the shock we are not interested in:
        Shocks(IdentifiedShocks)=0;
        % Here we re-run history after having killed off the relevant shock:
        CounterfactualX(:,tt)=mu+bb*vec(fliplr(CounterfactualX(:,tt-LagOrder:tt-1)))+a0*Shocks;
        tt=tt+1;
    end
  
     CounterfactualX3(:,:,ss)=CounterfactualX';
     CounterfactualMinusActualX3(:,:,ss)=CounterfactualX3(:,:,ss)-X;
   
    ss=ss+1;
end

%Test
sum(check)== NN;
%logical of 1 --> everything is fine

% Getting the percentiles of the counterfactual series:
CounterfactualSeries=sort(CounterfactualX3,3);
PercentilesCounterfactualSeries=CounterfactualSeries(:,:,Percentiles);

% Getting the percentiles of the difference between counterfactual and actual series:
CounterfactualMinusActualSeries=sort(CounterfactualMinusActualX3,3);
PercentilesCounterfactualMinusActualSeries=CounterfactualMinusActualSeries(:,:,Percentiles);

% Plotting the results
xx=1;
while xx<=N
    Counterfactual=squeeze(PercentilesCounterfactualSeries(:,xx,:));
    figure(6)
    subplot(2,5,xx)
    plot(Time,Counterfactual(:,1),'k',Time,Counterfactual(:,4:5),'r',Time,X(:,xx),'b',Time,zeros(size(Time)),'b:')
    xlim([Time(1) Time(T)])
    if xx==1
        title('RealGDPGrowth')
    elseif xx==2
        title('UnemploymentRate')
    elseif xx==3
        title('ThreeMonthTreasuryBill')
    elseif xx==4
        title('GDPDeflatorInflation')
    elseif xx==5
        title('Spread')
     %X=[RealGDPGrowth UnemploymentRate ThreeMonthTreasuryBill GDPDeflatorInflation Spread];    
    else
    end
    
    CounterfactualMinusActual=squeeze(PercentilesCounterfactualMinusActualSeries(:,xx,:));
    
    figure(6)
    subplot(2,5,xx+5)
    plot(Time,CounterfactualMinusActual(:,1),'k',Time,CounterfactualMinusActual(:,4:5),'r',Time,zeros(size(Time)),'b:')
    xlim([Time(1) Time(T)])
   if xx==1
        title('RealGDPGrowth')
    elseif xx==2
        title('UnemploymentRate')
    elseif xx==3
        title('ThreeMonthTreasuryBill')
    elseif xx==4
        title('GDPDeflatorInflation')
    elseif xx==5
        title('Spread') 
    else
    end
    %
    xx=xx+1;
end

%% Now the 4th shock
IdentifiedShocks=4;
CounterfactualX4=zeros(T,N,NN);
CounterfactualMinusActualX4=zeros(T,N,NN);
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
        Shocks=StructuralShocks(:,tt-LagOrder,ss);
        % Here we kill off the shock we are not interested in:
        Shocks(IdentifiedShocks)=0;
        % Here we re-run history after having killed off the relevant shock:
        CounterfactualX(:,tt)=mu+bb*vec(fliplr(CounterfactualX(:,tt-LagOrder:tt-1)))+a0*Shocks;
        tt=tt+1;
    end
 
    CounterfactualX4(:,:,ss)=CounterfactualX';
    CounterfactualMinusActualX4(:,:,ss)=CounterfactualX4(:,:,ss)-X;
    
    ss=ss+1;
end

%Test
sum(check)== NN;
%logical of 1 --> everything is fine

% Getting the percentiles of the counterfactual series:
CounterfactualSeries=sort(CounterfactualX4,3);
PercentilesCounterfactualSeries=CounterfactualSeries(:,:,Percentiles);

% Getting the percentiles of the difference between counterfactual and actual series:
CounterfactualMinusActualSeries=sort(CounterfactualMinusActualX4,3);
PercentilesCounterfactualMinusActualSeries=CounterfactualMinusActualSeries(:,:,Percentiles);

% Plotting the results
xx=1;
while xx<=N
    Counterfactual=squeeze(PercentilesCounterfactualSeries(:,xx,:));
    figure(7)
    subplot(2,5,xx)
    plot(Time,Counterfactual(:,1),'k',Time,Counterfactual(:,4:5),'r',Time,X(:,xx),'b',Time,zeros(size(Time)),'b:')
    xlim([Time(1) Time(T)])
    if xx==1
        title('RealGDPGrowth')
    elseif xx==2
        title('UnemploymentRate')
    elseif xx==3
        title('ThreeMonthTreasuryBill')
    elseif xx==4
        title('GDPDeflatorInflation')
    elseif xx==5
        title('Spread')
     %X=[RealGDPGrowth UnemploymentRate ThreeMonthTreasuryBill GDPDeflatorInflation Spread];    
    else
    end
    
    CounterfactualMinusActual=squeeze(PercentilesCounterfactualMinusActualSeries(:,xx,:));
    
    figure(7)
    subplot(2,5,xx+5)
    plot(Time,CounterfactualMinusActual(:,1),'k',Time,CounterfactualMinusActual(:,4:5),'r',Time,zeros(size(Time)),'b:')
    xlim([Time(1) Time(T)])
   if xx==1
        title('RealGDPGrowth')
    elseif xx==2
        title('UnemploymentRate')
    elseif xx==3
        title('ThreeMonthTreasuryBill')
    elseif xx==4
        title('GDPDeflatorInflation')
    elseif xx==5
        title('Spread') 
    else
    end
    %
    xx=xx+1;
end


%% Now the 5th shock

IdentifiedShocks=5;
CounterfactualX5=zeros(T,N,NN);
CounterfactualMinusActualX5=zeros(T,N,NN);
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
        Shocks=StructuralShocks(:,tt-LagOrder,ss);
        % Here we kill off the shock we are not interested in:
        Shocks(IdentifiedShocks)=0;
        % Here we re-run history after having killed off the relevant shock:
        CounterfactualX(:,tt)=mu+bb*vec(fliplr(CounterfactualX(:,tt-LagOrder:tt-1)))+a0*Shocks;
        tt=tt+1;
    end

    CounterfactualX5(:,:,ss)=CounterfactualX';
    CounterfactualMinusActualX5(:,:,ss)=CounterfactualX5(:,:,ss)-X;
    
    ss=ss+1;
end

%Test
sum(check)== NN;
%logical of 1 --> everything is fine

% Getting the percentiles of the counterfactual series:
CounterfactualSeries=sort(CounterfactualX5,3);
PercentilesCounterfactualSeries=CounterfactualSeries(:,:,Percentiles);

% Getting the percentiles of the difference between counterfactual and actual series:
CounterfactualMinusActualSeries=sort(CounterfactualMinusActualX5,3);
PercentilesCounterfactualMinusActualSeries=CounterfactualMinusActualSeries(:,:,Percentiles);

% Plotting the results
xx=1;
while xx<=N
    Counterfactual=squeeze(PercentilesCounterfactualSeries(:,xx,:));
    figure(8)
    subplot(2,5,xx)
    plot(Time,Counterfactual(:,1),'k',Time,Counterfactual(:,4:5),'r',Time,X(:,xx),'b',Time,zeros(size(Time)),'b:')
    xlim([Time(1) Time(T)])
    if xx==1
        title('RealGDPGrowth')
    elseif xx==2
        title('UnemploymentRate')
    elseif xx==3
        title('ThreeMonthTreasuryBill')
    elseif xx==4
        title('GDPDeflatorInflation')
    elseif xx==5
        title('Spread')
     %X=[RealGDPGrowth UnemploymentRate ThreeMonthTreasuryBill GDPDeflatorInflation Spread];    
    else
    end
    
    CounterfactualMinusActual=squeeze(PercentilesCounterfactualMinusActualSeries(:,xx,:));
    
    figure(8)
    subplot(2,5,xx+5)
    plot(Time,CounterfactualMinusActual(:,1),'k',Time,CounterfactualMinusActual(:,4:5),'r',Time,zeros(size(Time)),'b:')
    xlim([Time(1) Time(T)])
   if xx==1
        title('RealGDPGrowth')
    elseif xx==2
        title('UnemploymentRate')
    elseif xx==3
        title('ThreeMonthTreasuryBill')
    elseif xx==4
        title('GDPDeflatorInflation')
    elseif xx==5
        title('Spread') 
    else
    end
    %
    xx=xx+1;
end