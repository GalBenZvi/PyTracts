function CalcAxV6(fname)
grad_dirs=[]; vec=[]; rDWI=[];
load('AxCal3Ddata','bval', 'DT1000','mask', 'grad_dirs', 'MD1000', 'phi_q','theta_q', 'R_q', 'rDWI','Param','vec', 'MD250');
load([fname,'_alldata.mat'], 'MD3000', 'MD5000')
pfr=zeros(size(mask));
ph=zeros(size(mask));
pcsf=zeros(size(mask));
pasi=zeros(size(mask));
paxsi=zeros([size(mask) 160]);
CMDfr=zeros(size(mask));
CMDfh=zeros(size(mask));
pixpredictCSF=zeros(1,88);


mask2=mask;
B0s=squeeze(mean(rDWI(:,:,:,find(bval==0)),4));
[decayH, MD0]=predictH(DT1000, mask2, bval, grad_dirs, MD1000, MD3000, MD5000);

D_mat=[4 0 0; 0 4 0; 0 0 4];
for k=1:length(bval)
    pixpredictCSF(k)=exp(-4.*(grad_dirs(k,:)*D_mat*grad_dirs(k,:)'));
end
prcsf=pbrain(MD250);
ax=0.1:0.2:32;
[xlocs, ylocs, zlocs]=ind2sub(size(mask), find(mask2==1));
maskmap=zeros(size(mask));
A=ones(1,162);
b=1;
N = 160;
L = eye(N)+[[zeros(1,N-1);-eye(N-1)],zeros(N,1)];
L = [L,zeros(N,2)];
Aeq=ones(1,162);
beq=1;

h1 = waitbar(0,'Charmed/AxCaliber Analysis...');
for i=1:length(zlocs)
    waitbar(i / length(zlocs));
    maskmap(xlocs(i), ylocs(i), zlocs(i))=i;
    fvec=squeeze(vec(xlocs(i),ylocs(i),zlocs(i),:));
    B0sp=squeeze(B0s(xlocs(i), ylocs(i), zlocs(i)));
    ydata=squeeze(rDWI(xlocs(i), ylocs(i), zlocs(i),:));
    pixpredictH=squeeze(decayH(xlocs(i), ylocs(i), zlocs(i),:));
    vH=pixpredictH(:);
    vCSF=pixpredictCSF(:);
    decayR=predictR_singleR(theta_q, phi_q, fvec, R_q, Param.bigdel, Param.smalldel, grad_dirs, B0sp, ax/2, MD0(xlocs(i), ylocs(i), zlocs(i)));
    vR=double(decayR);
    vR=vR./max(vR(:));
    
    yd=gampdf(ax,2.5, 2.5);
    yd=yd./sum(yd);
    vRes=vR*yd';
    x0=double([0.5  1000]);
    min_val=[ 0  0];
    max_val=[ 1  5000];
    h=optimset('DiffMaxChange',1e-1,'DiffMinChange',1e-3,'MaxIter',20000,...
        'MaxFunEvals',20000,'TolX',1e-6,...
        'TolFun',1e-6, 'Display', 'off');
    [parameter_hat,~,~,~,~]=lsqnonlin('regression_function',...
        x0,min_val,max_val,h, ydata, vH, vRes, vCSF, prcsf(xlocs(i), ylocs(i), zlocs(i))) ;
    CMDfh(xlocs(i), ylocs(i), zlocs(i))=parameter_hat(1);
    CMDfr(xlocs(i), ylocs(i), zlocs(i))=1-parameter_hat(1)-prcsf(xlocs(i), ylocs(i), zlocs(i));
    vdata=ydata./parameter_hat(2);
    preds=[vR vH vCSF];
    lb=zeros(size(yd'));
    ub=ones(size(yd'));
    lb(161)=parameter_hat(1)-0.02;
    ub(161)=parameter_hat(1)+0.02;
    lb(162)=prcsf(xlocs(i), ylocs(i), zlocs(i))-0.02;
    ub(162)=prcsf(xlocs(i), ylocs(i), zlocs(i))+0.02;
    Lambda=1;
    Xprim=[preds; sqrt(Lambda)*L];
    yprim=[vdata;zeros(160,1)];
    options=optimoptions('lsqlin', 'Display',  'off');
    x = lsqlin(Xprim,yprim,A,b,Aeq,beq,lb,ub, yd, options);
    if isempty(x)==1
        x=zeros(160,1);
    end
    x(find(x<0))=0;
    a_h=parameter_hat(1);
    a_csf=prcsf(xlocs(i), ylocs(i), zlocs(i));
    a_fr=1-a_csf-a_h;
    if a_fr<0
        a_fr=0;
    end
    ph(xlocs(i), ylocs(i), zlocs(i))=a_h;
    pcsf(xlocs(i), ylocs(i), zlocs(i))=a_csf;
    pfr(xlocs(i), ylocs(i), zlocs(i))=sum(a_fr);
    pasi(xlocs(i), ylocs(i), zlocs(i))=sum(x(1:130).*ax(1:130)');
    paxsi(xlocs(i), ylocs(i), zlocs(i), :)=x(1:160);   
end
close(h1)

i1=load_untouch_nii([fname,'.nii']);
i1.img=pfr.*100;
i1.hdr.dime.dim(5)=1;
i1.hdr.dime.dim(1)=3;
i1.hdr.dime.datatype=64;
i1.hdr.dime.bitpix=64;
i1.hdr.dime.scl_slopre=1/100;
save_untouch_nii(i1, [fname,'_pfrS3.nii']);
i1.img=ph.*100;
save_untouch_nii(i1, [fname,'_pfhS3.nii']);
i1.img=pcsf.*100;
save_untouch_nii(i1, [fname,'_pcsfS3.nii']);
i1.img=pasi;
save_untouch_nii(i1, [fname,'_pasiS3.nii']);

save([fname,'_AxCalS3.mat']);

