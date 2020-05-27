function AxCaliberV6(fname,bval_fname,bvec_fname,inpsmldel, inpbigdel, inpgmax)
% in TheBase: (folder, 15.5, 60, 7.2)
a= strsplit(fname,filesep);
main_path = [filesep,fullfile(a{1:end-1})];
% load data

i1=load_untouch_nii(fname);
data=double(i1.img);

% A=fopen([bvalsname,'.bval']);
A = fopen(bval_fname);
bval=fscanf(A,'%f');
% A=fopen([bvalsname,'.bvec']);
A = fopen(bvec_fname);
bvec=fscanf(A,'%f')';
bvec=reshape(bvec,[length(bval) 3]);
bvec=bvec(:,[1 2 3]);
[sx sy sz sd]=size(data);
bval=bval./1000;



% DTI1000 FA calc
b1000locs=find(bval==1);
data1000=data(:,:,:,[1,b1000locs']);
rdata1000=zeros(size(data1000));
rdata1000(:,:,:,1)=data1000(:,:,:,1);
bvec1000=bvec([1 b1000locs'],:);
bval1000=bval([1 b1000locs']);
nbvec1000=zeros(size(bvec1000));
sd1000=length(bval1000);
[optimizer,metric] = imregconfig('multimodal');
optimizer.InitialRadius = optimizer.InitialRadius/300;
optimizer.MaximumIterations = 200;
h = waitbar(0,'Motion Correction for DTI...');
for i=1:(length(b1000locs)+1)
    waitbar(i/(length(b1000locs)+1));
    [rdata1000(:,:,:,i)] = imregister(data1000(:,:,:,i), rdata1000(:,:,:,1), 'affine', optimizer, metric);
    rdata1000(:,:,:,i)=smooth3(rdata1000(:,:,:,i), 'gaussian');
    tform = imregtform(data1000(:,:,:,i), rdata1000(:,:,:,1), 'affine', optimizer, metric);
    nbvec1000(i,:)=bvec1000(i,:)*tform.T(1:3,1:3);
end
close(h);

tmpim=mean(rdata1000(:,:,:,2:(length(b1000locs)+1)),4);
mask=zeros(size(tmpim));
mask(find(tmpim>40))=1;
for i=1:sz
    A1=mask(:,:,i);
    A1=imfill(A1, 'holes');
    mask(:,:,i)=A1;
end

[FA1000, MD1000, DT1000, Eigval1000, Eigvec1000]=dti5(bval1000, bvec1000, rdata1000, 1, mask);

i1.hdr.dime.dim(5)=1;
i1.hdr.dime.dim(1)=3;
i1.hdr.dime.datatype=64;
i1.hdr.dime.bitpix=64;
i1.hdr.dime.scl_slopre=1/100;
i1.img=FA1000.*100;
save_untouch_nii(i1, [fname,'_FA.nii']);
i1.img=MD1000.*100;
save_untouch_nii(i1, [fname,'_MD.nii']);


% Experimental Parameters and calculations
Param.nb0 = 1;
bvfac=sqrt(bval./max(bval));
grad_dirs=bvec.*repmat(bvfac, [1 3]);
gamma=4257 ;
Param.smalldel = inpsmldel; %in ms
Param.bigdel = inpbigdel; % in ms
Param.TE = 125; %in ms
%Gmax calc: sqrt(bval*100/(7.178e8*0.03^2*(0.06-0.01)))
Param.maxG = inpgmax; % in G/cm = 1/10 mT/m
Param.ndir = length(grad_dirs);
Param.maxq = gamma.*Param.smalldel.*Param.maxG./10e6;
Param.qdirs = grad_dirs.*Param.maxq;
[phi_q, theta_q R_q]=cart2sph(grad_dirs(:,1), grad_dirs(:,2), -grad_dirs(:,3));
R_q=R_q.*Param.maxq;
Param.bval=4.*pi^2.*R_q.^2.*(Param.bigdel-Param.smalldel/3);
Param.phi=phi_q;
Param.theta=theta_q;
Param.R=R_q;
ax=0.1:0.2:32;
gw3=gampdf(ax,2,2);
gw3=gw3./sum(gw3);

% Simulate CHARMED for motion correction
FixedB0=data(:,:,:,1,1);
[optimizer,metric] = imregconfig('multimodal');
optimizer.InitialRadius = optimizer.InitialRadius/1000;
optimizer.MaximumIterations = 300;
FA=FA1000;
vec=Eigvec1000;
mag=FA;
DTmaps=DT1000; %(:,:,:,[1 2 3 5 6 9]);

B0s=data(:,:,:,1:Param.nb0,1);
for i=1:Param.nb0
    [rB0(:,:,:,i)] = imregister(FixedB0, B0s(:,:,:,i), 'affine', optimizer, metric);
    rB0(:,:,:,i)=smooth3(rB0(:,:,:,i), 'box', [3 3 3]);
end
B0map=mean(rB0(:,:,:,1:Param.nb0),4);
simdwis = simchm(B0map, FA, DTmaps, mag, vec, grad_dirs, mask, Param.maxq, Param.bigdel, Param.smalldel, theta_q, phi_q, ax/2, gw3, Param.R, bval, MD1000);
osimdwis=simdwis;
for i=1:length(bval)
    A1=simdwis(:,:,:,i);
    A2=A1(:);
    A2=A2(find(A2>0));
    mA2=sort(A2);
    locl=round(0.995*length(mA2));
    A1(find(A1>(mA2(locl))))=mA2(locl);
    simdwis(:,:,:,i)=A1;
end
% UNDISTORT - Registration and gradient reorientation

[optimizer,metric] = imregconfig('multimodal');
optimizer.InitialRadius = optimizer.InitialRadius/100;
optimizer.MaximumIterations = 200;

h = waitbar(0,'UNDISTORT...');
parfor i=1:length(grad_dirs)
    A1=simdwis(:,:,:,i);
    A1=smooth3(A1,'box',[3 3 3]);
    maxi=max(max(max(A1)));
    J = imnoise(double(A1)./maxi,'gaussian', 0.05, 0.005);
    KK=smooth3(data(:,:,:,i), 'box', [3 3 3]);
    skk=sort(KK(:));
    locl=round(0.8*length(skk));
    KK(find(KK>(skk(locl))))=skk(locl);
    KK=KK./skk(locl);
     tform = imregtform(KK, J, 'affine', optimizer, metric);
    newgrad(i,:)=[grad_dirs(i,:) sqrt(sum(grad_dirs(i,:).^2))]*tform.T;
    newgrad2(i,:)=grad_dirs(i,:)*tform.T(1:3,1:3);
    rDWI(:,:,:,i)=imwarp(data(:,:,:,i), tform, 'OutputView', imref3d(size(KK)));
    waitbar(i/length(grad_dirs));
end
close(h)
for i=1:length(newgrad);
    grad_dirs(i,:)=newgrad(i,4).*newgrad(i,1:3)./sqrt(sum(newgrad(i,1:3).^2));
end
grad_dirs=newgrad2;
grad_dirs(1:Param.nb0,:)=zeros(Param.nb0,3);

for i=1:length(bval)
    rDWI(:,:,:,i)=smooth3(rDWI(:,:,:,i), 'gaussian');
end

[FA250, MD250, DT250, Eigval250, Eigvec250]=dti5(bval, grad_dirs, rDWI, 0.25, mask);
[FA3000, MD3000, DT3000, Eigval3000, Eigvec3000]=dti5(bval, grad_dirs, rDWI, 3, mask);
[FA5000, MD5000, DT5000, Eigval5000, Eigvec5000]=dti5(bval, grad_dirs, rDWI, 5, mask);
save([fname,'_alldata.mat']);
% AxCalibering!
cd(main_path);
save('AxCal3Ddata.mat','mask', 'rDWI', 'Param', 'grad_dirs', 'theta_q', 'phi_q', 'vec', 'R_q', 'bval', 'DTmaps', 'FA', 'DT1000', 'MD1000', 'DT250', 'MD250');
CalcCHMV6(fname);

