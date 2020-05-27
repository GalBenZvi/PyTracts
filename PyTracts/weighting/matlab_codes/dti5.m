function [FA, MD, DT, Eigval, Eigvec]=dti5(bval, bvec, data, bv, mask)

for i=1:length(bval)
    bvec2(i,:)=bvec(i,:)./norm(bvec(i,:));
end

Bvalue=bval*1000;
b0locs=find(Bvalue==0);
bvlocs=find(Bvalue==bv*1000);

S0=data(:,:,:,b0locs);
S=data(:,:,:,bvlocs);
Bvalue=Bvalue(bvlocs);
H=bvec2(bvlocs,:);

S0=mean(S0, 4);

b=zeros([3 3 size(H,1)]);
for i=1:size(H,1),
    b(:,:,i)=Bvalue(i)*H(i,:)'*H(i,:);
end

Slog=zeros(size(S),'double');
for i=1:size(H,1),
    Slog(:,:,:,i)=log((S(:,:,:,i)./S0)+eps);
end

Bv=squeeze([b(1,1,:),2*b(1,2,:),2*b(1,3,:),b(2,2,:),2*b(2,3,:),b(3,3,:)])';

DT=zeros([size(S0) 6],'single');
Eigval=zeros([size(S0) 3],'single');
FA=zeros(size(S0),'single');
MD=zeros(size(S0),'single');
Eigvec=zeros([size(S0) 3],'single');

[xlocs ylocs zlocs]=ind2sub(size(mask), find(mask==1));

for i=1:length(xlocs)
    Z=-squeeze(Slog(xlocs(i),ylocs(i),zlocs(i),:));
    M=Bv\Z;
    DiffusionTensor=[M(1) M(2) M(3); M(2) M(4) M(5); M(3) M(5) M(6)];
    
    [EigenVectors,D]=eig(DiffusionTensor); EigenValues=diag(D);
    [t,index]=sort(EigenValues);
    EigenValues=EigenValues(index)*1000; EigenVectors=EigenVectors(:,index);
    EigenValues_old=EigenValues;
    if((EigenValues(1)<0)&&(EigenValues(2)<0)&&(EigenValues(3)<0)), EigenValues=abs(EigenValues);end
    if(EigenValues(1)<=0), EigenValues(1)=eps; end
    if(EigenValues(2)<=0), EigenValues(2)=eps; end
    
    MDv=(EigenValues(1)+EigenValues(2)+EigenValues(3))/3;
    
    FAv=sqrt(1.5)*( sqrt((EigenValues(1)-MDv).^2+(EigenValues(2)-MDv).^2+(EigenValues(3)-MDv).^2)./sqrt(EigenValues(1).^2+EigenValues(2).^2+EigenValues(3).^2) );
    
    MD(xlocs(i), ylocs(i), zlocs(i))=MDv;
    Eigval(xlocs(i), ylocs(i), zlocs(i),:)=EigenValues;
    DT(xlocs(i), ylocs(i), zlocs(i),:)=[DiffusionTensor(1:3) DiffusionTensor(5:6) DiffusionTensor(9)];
    FA(xlocs(i), ylocs(i), zlocs(i))=FAv;
    Eigvec(xlocs(i), ylocs(i), zlocs(i),:)=EigenVectors(:,end)*EigenValues_old(end);
end

