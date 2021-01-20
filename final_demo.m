%------------------------------------------------------------------ 
% Demo of the TSR algorithm for the following paper: 
% Gifani, Parisa, et al. "Temporal Super Resolution Enhancement of Echocardiographic Images 
% based on Sparse Representations." (2015).

% this Demo is created for an IVTC signal.  
%------------------------------------------------------------------ 
 %%
clc
clear all
close all

warning OFF

root_image='data';
root=[root_image ];% dictionary folder
folderName='Dictionary'; % the name of the Dictionary Folder
mkdir([root '\' folderName '\rgb' ])

aa=dir(root_image);
T1=length(aa)-3;%Length of the Original IVTC signal (number of original frames)
T2=65; %Length of Temporal Super Enhanced signal (TSR signal)

L1=T1;
L2=T2;



%**********Cerating Dictionaries**************************
%%
create_dict_LR(root,folderName,T1)
create_dict_HR(root,folderName,T1,T2)

%*********************************************************
%%

%**********Loading Dictionaries***************************
disp('Loading Dictionaries...')
load([root '\' folderName '\dict_db4_' num2str(T1) '.mat'])
load([root '\' folderName '\dict_dmey_' num2str(T1) '.mat'])
load([root '\' folderName '\dict_sym2_' num2str(T1) '.mat'])
load([root '\' folderName '\dict_sym4_' num2str(T1) '.mat'])
load([root '\' folderName '\dict_sin_' num2str(T1) '.mat'])
load([root '\' folderName '\dict_cos_' num2str(T1) '.mat'])

LR_dict_Wavelet_Part=[dict_db4,dict_sym2,dict_sym4,dict_dmey ];
LR_dict=[dict_db4,dict_sym2,dict_sym4,dict_dmey ,dict_sin,dict_cos];

load([root '\' folderName '\dict_db4_' num2str(T2) '.mat'])
load([root '\' folderName '\dict_dmey_' num2str(T2) '.mat'])
load([root '\' folderName '\dict_sym2_' num2str(T2) '.mat'])
load([root '\' folderName '\dict_sym4_' num2str(T2) '.mat'])
load([root '\' folderName '\dict_sin_' num2str(T2) '.mat'])
load([root '\' folderName '\dict_cos_' num2str(T2) '.mat'])

HR_dict_Wavelet_Part=[dict_db4,dict_sym2,dict_sym4,dict_dmey];
HR_dict=[dict_db4,dict_sym2,dict_sym4,dict_dmey ,dict_sin,dict_cos];
%*************************************************************
TSR_signal=zeros(1,T2);

mydict1=LR_dict;
mydict2=HR_dict;
mydict = sparse(mydict1);


im1=imread([ root_image '\im (' num2str(1) ').jpg']);
[h,w,z]=size(im1);

all_im=uint8(zeros(h,w,T1));
all_im_post=uint8(zeros(h,w,T2));
cnt=1;
for i=1:L1
    im=(rgb2gray(imread([ root '\im (' num2str(i) ').jpg'])));
    all_im(:,:,cnt)=im(:,:);
    cnt=cnt+1;
end

wave_dic=HR_dict_Wavelet_Part;
mydict1_wave=LR_dict_Wavelet_Part;

for ii=1:h
    ii
    %     tic
    for jj=1:w
%         jj
        lineH=double(all_im(ii,jj,:));
        xh=reshape(lineH,1,L1);
        
        xd=zeros(1,L1);
        for dd=1:L1-1
            xd(dd)=xh(dd+1)-xh(dd);
        end
        if ~any(xh)  ||   ~any(xd) %sum(xh)==L1 
            
            all_im_post(ii,jj,:)=0;
            
        else
            
            [ x_BCS,parameter_error, dynamics_error,used]=bayesian_oscillation3(xh', mydict1, 1); % baysian sparse learning
            
            
            
            sig1=zeros(L1,1);
            for im=1:size(mydict1,2)
                
                sig=x_BCS(im).*mydict1(:,im);
                sig1=sig+sig1;
                
            end
           
            
            
            ratio=L2/L1;
           
            [mm,nn]=find(x_BCS);
            new_mm1=(mm*ratio);
            new_mm=round(mm*ratio);
           
            value_mm=x_BCS(mm);
            BCS_new=zeros(floor(size(x_BCS,1)*ratio),1);
           
            BCS_new(new_mm)=value_mm;
 
            sig2=zeros(L2,1);
            for j=1:size(wave_dic,2)
                if abs(BCS_new(j))>0
                    sig=BCS_new(j)*mydict2(:,j);
                    sig2=sig+sig2;
                    
                end
            end
            
            sdw1=size(mydict1_wave,2);
            [m,n]=find(x_BCS);
            [m1,n1]=find(m>sdw1);
            
            sparse_coefs=x_BCS(m);
            sparse_coefs_sin_cos=x_BCS(m(m1));
            active_atoms=zeros(L2,size(m1,1));
            x = 1:1:L2;
            cnt=1;
            for j=1:size(m)
                if m(j)<sdw1+126 && m(j)>sdw1
                    atom='sin';
                    atom_sin=sin((m(j)-sdw1)*x/L2);
                    active_atoms(:,cnt)=atom_sin;
                    cnt=cnt+1;
                end
                
                
                if m(j)<sdw1+251 && m(j)>sdw1+125
                    atom='cos';
                    atom_cos=cos((m(j)-(sdw1+125))*x/L2);
                    active_atoms(:,cnt)=atom_cos;
                    cnt=cnt+1;
                    
                end
                
                
                
            end
            %
            S1 = sum(active_atoms.*active_atoms,1);
            coefs_active_atoms= active_atoms./repmat(S1.^0.5,L2,1);
            
            sig3=zeros(L2,1);
            for q=1:size(m1)
                sig=sparse_coefs_sin_cos(q)*coefs_active_atoms(:,q);
                sig3=sig+sig3;
                               
            end
            
            
            
            final=sig3+sig2;
            
            ratio2=max(sig1)/max(final);
            final2=final*ratio2;
            
            all_im_post(ii,jj,:)=uint8(final2);
            %
%                     figure;
%                     subplot(3,1,1);plot(xh);xlim([0,L1]);ylim([0,200])
%                     subplot(3,1,2);plot(sig1);xlim([0,L1]);ylim([0,200])
%                     subplot(3,1,3);plot(final2);xlim([0,L2]);ylim([0,200])
                    clear sparse_coefs_sin_cos coefs_active_atoms active_atoms cnt m m1 n n1
                    clear final sig1 sig2 sig3 sdw1 BCS_new value_mm x_BCS mm nn new_mm ratio2
%                     pause
        end
        %     toc
    end
end  
    
%save([root '\' folderName '\data' num2str(L1) '_' num2str(L2) '.mat'],'all_im_post')
% % % 
% for k=1:360
%     figure(1)
% imshow(all_im_post(:,:,k))
% end
% % % 
for k=1:L2
    imwrite(all_im_post(:,:,k),[[root '\' folderName] '\im' num2str(k) '.bmp'])
    
end

rgb_im=uint8(zeros(h,w,3,T2));

for k=1:L2
    rgb_im(:,:,1,k)=all_im_post(:,:,k);
    rgb_im(:,:,2,k)=all_im_post(:,:,k);
    rgb_im(:,:,3,k)=all_im_post(:,:,k);
    imwrite(rgb_im(:,:,:,k),[[root '/' folderName '/rgb'] '\im' num2str(k) '.jpg'])
    
end

    