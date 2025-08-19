function model=projectPCA(model,X)
% model=projectPCA(model,X)
% 
% Modified 17/12-08

[n,p]=size(X.d);
model.projected.X=X;

if p~=size(model.X.d,2)
    error('Number of variables does not match model');
end
    
%preprocess data
if isfield(model,'Xpp')
X.d=(X.d-repmat(model.Xpp.mean,n,1))./repmat(model.Xpp.std,n,1);
end

model.projected.T=X.d*model.P(:,1:model.Aopt);
model.projected.Res=X.d-model.projected.T*model.P(:,1:model.Aopt)';

Hi=zeros(n,1);
for i=1:model.Aopt
    Hi=Hi+model.projected.T(:,i).^2/(model.T(:,i)'*model.T(:,i));
end
Hi=Hi+1/size(model.X.d,1);

model.projected.Hi=Hi;
model.projected.Si=sqrt(mean(model.projected.Res.^2,2));
model.projected.SampRes = mean(model.projected.Res.^2,2);
model.projected.HiLimit=model.LevLim*(model.Aopt+1)/size(model.X.d,1);
model.projected.S0=sqrt(model.ResVarVal(model.Aopt+1));
model.projected.SiLimit=model.projected.S0*sqrt(finv(model.p,1,size(model.X.d,1)-model.Aopt-1));


