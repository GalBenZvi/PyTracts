function correlation_map = rspipe(func_name)
% func_name = '/home/gal/derivatives/sub-12/func/func.feat/filtered_func_data.nii.gz';
i1=load_untouch_nii(func_name);
data = i1.img;

[sx sy sz sr]=size(data);
%% step 2 registration

% [optimizer,metric] = imregconfig('multimodal');
% optimizer.InitialRadius = optimizer.InitialRadius/300;
% optimizer.MaximumIterations = 200;
% %h = waitbar(0,'Please wait...');
% 
% [sx sy sz sr]=size(data);
% 
% fixed=data(:,:,:,1);
% tic
% for i=1:sr
%     disp(i/sr)
%     [data(:,:,:,i)] = imregister(data(:,:,:,i), fixed, 'affine', optimizer, metric);
% end
% toc
disp('step 1 - division to small regions ~2000')
mask = data(:,:,:,1);
mask(mask~=0)=1;
mean_img=mean(data,4);
mean_img=mean_img.*mask;
[divided_img_template, numl]=superpixels3(mean_img,10000);
divided_img_template(find(mask==0))=0;
divided_img=zeros(size(divided_img_template));
ind=1;
for i=1:10000
    disp(i/10000)
    temp_group=divided_img_template;
    temp_group(find(temp_group~=i))=0;
    if sum(temp_group(:))>0
        divided_img(find(divided_img_template==i))=ind;
        ind=ind+1;
    end
end


disp('step 2 - ICA on mean of 2000 regions')

grouped_timecourse=repmat(divided_img,[1 1 1 sr]);
for i=1:max(divided_img(:))
    disp(i/max(divided_img(:)))
    temp_group=data;
    temp_group(find(grouped_timecourse~=i))=0;
    npix=length(divided_img(find(divided_img==i)));
    for j=1:sr
        compo(j,i)=mean(mean(mean(temp_group(:,:,:,j))));
    end
    compo(:,i)=compo(:,i)./npix;
end
zcompo=zscore(compo);
icompo=rica(zcompo',15);
correlation_map=zeros([size(mask) 15]);
[xlocs ylocs zlocs]=ind2sub(size(mask), find(mask==1));

disp('step 3 - voxel correlation with 15 ICA componentns')

for i=1:15
    disp(i)
    vec1=icompo.TransformWeights(:,i);
    for j=1:length(xlocs)
        [r p]=corrcoef(vec1,zscore(double(squeeze(data(xlocs(j), ylocs(j), zlocs(j),:)))));
        correlation_map(xlocs(j), ylocs(j), zlocs(j),i)=r(2);
    end
end

% %% display images. 
% cmap=gray(100);
% cmap(101:200,:)=hot(100);
% 
% for j=1:15
%     bsl=load_untouch_nii('/home/gal/Trials/sub-12/func/func.feat/mean_func.nii.gz');
%     tim=bsl.img;
%     tim=tim./max(tim(:));
%     net1=cormap(:,:,:,j);
%     net1(find(net1<0))=0;
%     net1=net1./max(net1(:));
%     net1(find(net1>0.6))=0.6;
%     net1=net1./max(net1(:));
%     net1(find(net1<0.2))=0;
%     tim(find(net1>0))=1;
%     ctim=tim+net1;
%     ctim(1,1,:)=2;
%     ctim(1,:,1)=2;
%     ctim(:,1,:)=2;
%     figure;
%     ind=1;
%     for i=16:6:82
%         A1=squeeze(ctim(:,i,:));
%         subplot(3,4,ind); imagesc(imrotate(A1,-90)); axis equal; axis off
%         colormap(cmap);
%         ind=ind+1;
%     end
% end