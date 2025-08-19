function []=pcaplots(model,A1,A2,plottype,spec)
% []=pcaplots(model,A1,A2,plottype)
%
%
% Plots:
% scores for A1 and A2,
% loadings for A1 and A2
% explained variance for all components
% normal probability plot of residuals
%
% INPUT
% model       model struct from mypls
% A1,A2       PCs to plot for scores and loadings
% plottype    Optional.   =1: plots one figure with subplots (default)
%                         =2: plots each plot in a separate figure
% spec        Optional.   =1: plots loadings as curves

if ~exist('A1','var')
    A1=1;
end

if ~exist('A2','var')
    A2=2;
    if size(model.T,2)==1; A2=1; end
end

if ~exist('plottype','var')
    plottype=1;
end

if ~exist('spec','var')
    spec=0;
end



%Scores
if plottype==1
    figure('Name','PCA overview');
subplot(2,2,1)
else
    figure('Name','Scores');
end

if size(model.T,1)<300
H=plot(model.T(:,A1),model.T(:,A2),'.');
hold on
T1range=max(model.T(:,A1))-min(model.T(:,A1));
text(model.T(:,A1)+T1range/100,model.T(:,A2),model.X.i)

else
        dscatter(model.T(:,A1),model.T(:,A2))
    
end
title('Scores')
xlabel(['PC ' num2str(A1) ', ' num2str(round(model.ExpVar(A1+1)-model.ExpVar(A1))) '%'])
ylabel(['PC ' num2str(A2) ', ' num2str(round(model.ExpVar(A2+1)-model.ExpVar(A2))) '%'])
grid on
h=gca;
if get(h,'Xlim')==0
    set(h,'Xtick',0)
else
set(h,'Xtick',sort(unique([get(h,'Xlim') 0])))
end
if get(h,'Ylim')==0
    set(h,'Ytick',0)
else
set(h,'Ytick',sort([get(h,'Ylim') 0]))
end

%Correlation Loadings
if plottype==1
    subplot(2,2,2)
else
    figure('Name','Loadings');
end

if spec~=1
hold off
plot(corr(model.T(:,A1),model.X.d),corr(model.T(:,A2),model.X.d),'.')
hold on
H=text(corr(model.T(:,A1),model.X.d)+2/100,corr(model.T(:,A2),model.X.d),model.X.v);
set(H,'Color','b')
title('Correlation Loadings')
xlabel(['PC ' num2str(A1)])
ylabel(['PC ' num2str(A2)])
%circle 100%
[X,Y] = pol2cart(linspace(0,2*pi,100),ones(1,100));
plot(X,Y,'k');
%circle 50%
[X,Y] = pol2cart(linspace(0,2*pi,100),ones(1,100)*sqrt(0.5));
plot(X,Y,'k');
H=line([0 0],[-1 1]);
set(H,'linestyle',':','Color','k')
H=line([-1 1],[0 0]);
set(H,'linestyle',':','Color','k')
else
    plot(str2num(model.X.v),model.P(:,1:model.Aopt)); axis tight
    legend(strcat('PC',num2str((1:model.Aopt)')),'AutoUpdate','off')
    title('Loadings')
end

%Explained variance
if plottype==1
    subplot(2,2,3)
else
    figure('Name','Explained variance');
end


plot(0:size(model.T,2),model.ExpVar)
hold on
plot(model.Aopt,model.ExpVar(model.Aopt+1),'o')
title('Explained variance, Total')
xlabel('Components')
H=line([0 size(model.T,2)],[0 0]);
set(H,'color','k')
set(gca,'ylim',[0 100])

% leverage plot
if plottype==1
    subplot(2,2,4)
else
    figure('Name','Leverage');
end

if size(model.T,1)<1000
plot(model.Hi(:,model.Aopt),model.SampRes(:,model.Aopt+1),'.')
Hirange=max(model.Hi(:,model.Aopt))-min(model.Hi(:,model.Aopt));
text(model.Hi(:,model.Aopt)+Hirange/100,model.SampRes(:,model.Aopt+1),model.X.i)
else
     dscatter(model.Hi(:,model.Aopt),model.SampRes(:,model.Aopt+1))
end
title(['Influence at ' num2str(model.Aopt) ' PCs'])
xlabel('Leverage')
ylabel('Residual')