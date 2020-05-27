function [decay, MD0]=predictH(DTmaps1, mask, bval, grad_dirs, MD1000, MD3000, MD5000)
[xlocs, ylocs, zlocs]=ind2sub(size(mask), find(mask==1));
MD0=zeros(size(mask));
DTmaps1=1000.*DTmaps1;
decay=zeros([size(mask) length(bval)]);
DTmaps(:,:,:,1,1)=DTmaps1(:,:,:,1);
DTmaps(:,:,:,1,2)=DTmaps1(:,:,:,2);
DTmaps(:,:,:,1,3)=DTmaps1(:,:,:,3);
DTmaps(:,:,:,2,1)=DTmaps1(:,:,:,2);
DTmaps(:,:,:,2,2)=DTmaps1(:,:,:,4);
DTmaps(:,:,:,2,3)=DTmaps1(:,:,:,5);
DTmaps(:,:,:,3,1)=DTmaps1(:,:,:,3);
DTmaps(:,:,:,3,2)=DTmaps1(:,:,:,5);
DTmaps(:,:,:,3,3)=DTmaps1(:,:,:,6);    
xm(1)=1000; xm(2)=3000; xm(3)=5000;
midi=double(cat(4,MD1000, MD3000, MD5000));
for i=1:length(zlocs)
    if midi(xlocs(i),ylocs(i),zlocs(i),1)*midi(xlocs(i),ylocs(i),zlocs(i),2)*midi(xlocs(i),ylocs(i),zlocs(i),3)~=0
        f = fit(double(xm'),squeeze(midi(xlocs(i), ylocs(i), zlocs(i),:)),'exp1');
        fac=f.a./midi(xlocs(i), ylocs(i), zlocs(i),1);
        MD0(xlocs(i), ylocs(i), zlocs(i))=1.5.*f.a;
        D_mat=squeeze(1.5*fac.*DTmaps(xlocs(i), ylocs(i), zlocs(i), : ,:));
        for j=1:length(bval)
            decay(xlocs(i), ylocs(i), zlocs(i), j)=exp(-4.*(grad_dirs(j,:)*D_mat*grad_dirs(j,:)'));
        end
    else
        continue
    end
end










