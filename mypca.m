function model=mypca(X,A,preproc)
% model=mypca(X,A,preproc)
% 
% 
% INPUT
% X         Descriptor data matrix, saisir structure
% A         Maximum number of components
% preproc  optional. cell array containing preprocessing methods for X
% 
% OUTPUT
% model
% 
% Created: 10/11 2008
% Modified: 
% Status:  works 10/11 2008
%
% Ingrid Måge


[n,p]=size(X.d);


%set labels if they are not given in saisir
if  ~isfield(X,'v')
    X.v=num2str((1:px)');
elseif isempty(X.v)
    X.v=num2str((1:px)');
end

if ~isfield(X,'i')
    X.i=num2str((1:nx)');
elseif isempty(X.i) 
    X.i=num2str((1:nx)');
end

model.type = 'pca';
model.X=X; %save raw data

if nargin==3
    X=mypreprocess(X,preproc);
    model.Xpp=X.pp;
    X=rmfield(X,'pp');
end


[U,S,P]=svd(X.d,0);
T=U*S;
T=T(:,1:A);
P=P(:,1:A);
totvar=sum((diag(S).^2));
S=S(1:A,1:A);

ExpVar=[0; cumsum(diag(S.^2)./repmat(totvar,A,1)*100)]';

E=cell(A+1,1);
ResVar=zeros(A+1,1);
ResVarVal=zeros(A+1,1);
SampRes=zeros(n,A+1);
VarRes=zeros(p,A+1);
Hi=zeros(n,A+1);
CorrLoads = nan(size(P));
for i=0:A
    E{i+1}=X.d-T(:,1:i)*P(:,1:i)';
    ResVar(i+1)=mean(mean(E{i+1}.^2));
    SampRes(:,i+1)=mean((E{i+1}.^2)')';
    VarRes(:,i+1)=mean((E{i+1}.^2))';
    if i>0
        Hi(:,i+1)=Hi(:,i)+T(:,i).^2/(T(:,i)'*T(:,i));
    end
    ResVarVal(i+1)=mean(mean((E{i+1}).^2./repmat((1-Hi(:,i+1)).^2,1,p)));

    if i>0
    CorrLoads(:,i) = corr(T(:,i),X.d);
    end
end
   
Hi(:,1)=[];
Hi=Hi+1/n;

model.T=T;
model.P=P;
model.CorrLoads = CorrLoads;
%model.E=E;
model.ResVar=ResVar;
model.ResVarVal=ResVarVal;
model.ExpVar=ExpVar;
model.SampRes=SampRes;
model.VarRes=VarRes;
model.VarExp=100-VarRes'./repmat(VarRes(:,1)',A+1,1)*100;
model.Hi=Hi;

crit = (1:A)'*0.02*mean(ResVar(1)) + mean(ResVar(2:end),2);
[dummy, Aopt] = min(crit);
model.Aopt = Aopt;
model.LevLim=3;
model.p=0.95;

