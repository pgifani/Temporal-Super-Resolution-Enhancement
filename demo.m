clc
clear all
close all

T1=81;
T2=405;

% root1='E:\';%images database folder

root2=;% dictionary folder


load([root2 '\dict_db4_' num2str(L1) '.mat'])
load([root2 '\dict_dmey_' num2str(L1) '.mat'])
load([root2 '\dict_sym1_' num2str(L1) '.mat'])
load([root2 '\dict_sym4_' num2str(L1) '.mat'])
load([root2 '\dict_sin_' num2str(L1) '.mat'])
load([root2 '\dict_cos_' num2str(L1) '.mat'])

LR_dict_Wavelet_Part=[dict_db4,dict_sym1,dict_sym4,dict_dmey ];
LR_dict=[dict_db4,dict_sym1,dict_sym4,dict_dmey ,dict_sin,dict_cos];

load([root2 '\dict_db4_' num2str(L2) '.mat'])
load([root2 '\dict_dmey_' num2str(L2) '.mat'])
load([root2 '\dict_sym1_' num2str(L2) '.mat'])
load([root2 '\dict_sym4_' num2str(L2) '.mat'])
load([root2 '\dict_sin_' num2str(L2) '.mat'])
load([root2 '\dict_cos_' num2str(L2) '.mat'])

HR_dict_Wavelet_Part=[dict_db4,dict_sym1,dict_sym4,dict_dmey];
HR_dict=[dict_db4,dict_sym1,dict_sym4,dict_dmey ,dict_sin,dict_cos];

TSR_signal=zeros(1,T2);

xd=zeros(1,T1);
for dd=1:T1-1
    xd(dd)=xh(dd+1)-xh(dd);
end

if ~any(xh)  ||   ~any(xd) 
    
     TSR_signal(1:T2)=0;
    
else
    
     Sparse_Sig=bayesian_sparse_coding(xh', LR_dict);
    
     TSR_signal(1:T2)=TSR_Gifani(HR_dict_Wavelet_Part,HR_dict,LR_dict_Wavelet_Part,LR_dict,Sparse_Sig,T1,T2);
    
  
end





