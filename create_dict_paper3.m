clc
clear all
close all
folderName='46-230';
% load(['E:\compress sensing\BCS\Oscillator\Oscillator\paper4\' folderName '.mat'])
mkdir(['E:\gifani\Dr Ranjbar\Dict\' folderName])
root=['E:\gifani\Dr Ranjbar\Dict\' folderName];
L1=46;
L2=46;
a_SIG=2;

ratio=L1/L2;
% ratio=4;
a=2*ratio;
for i=1:L1
    i
ySIG=zeros(1,L1);
ySIG(i)=1;
wname='db4';
wtype = wavemngr('type',wname);
iter=10;
switch wtype
  case {1,3} , [phi,psi,xval] = wavefun(wname,iter);
  case 2     , [phi1,psi1,phi2,psi2,xval] = wavefun(wname,iter);
  case {4,5} , [psi,xval] = wavefun(wname,iter);
end

step = xval(2)-xval(1);
switch wtype
  case {1,3,4,5} , out1 = cumsum(psi)*step;
  case 2         , out1 = cumsum(psi1)*step; out3 = cumsum(psi2)*step;
end
 

xval = xval-xval(1);
l_sig=size(xval,2);
xMaxWAV = xval(end);
stepWAV = xval(2)-xval(1);
step1=stepWAV;
R_step=1/stepWAV;
j=1+floor(0:R_step/(2*ratio):l_sig);
% f = fliplr(out1(j));
f = (psi(j));
f_phi = (phi(j));

coefs(i,:) = sqrt(a)*wkeep1((wconv1(ySIG,f)),L1);
phi_coefs(i,:) = wkeep1((wconv1(ySIG,f_phi)),L1);

end

% figure;plot(coefs(37,:))
% figure;plot(phi_coefs(37,:))

% figure;plot(psi)
% figure;plot(phi)

S1 = sum(coefs.*coefs,1);
coefs_db4 = coefs./repmat(S1.^0.5,L1,1);


S2 = sum(phi_coefs.*phi_coefs,1);
phi_coefs_db4 = phi_coefs./repmat(S2.^0.5,L1,1);

dict_db4=[phi_coefs_db4,coefs_db4];
save( [root '\dict_db4_' num2str(L1) '.mat'] ,'dict_db4')
clear coefs phi_coefs phi psi out1

for i=1:L1
    i
ySIG=zeros(1,L1);
ySIG(i)=1;

wname='sym2';
wtype = wavemngr('type',wname);
iter=10;
switch wtype
  case {1,3} , [phi,psi,xval] = wavefun(wname,iter);
  case 2     , [phi1,psi1,phi2,psi2,xval] = wavefun(wname,iter);
  case {4,5} , [psi,xval] = wavefun(wname,iter);
end

step = xval(2)-xval(1);
switch wtype
  case {1,3,4,5} , out1 = cumsum(psi)*step;
  case 2         , out1 = cumsum(psi1)*step; out3 = cumsum(psi2)*step;
end
 

xval = xval-xval(1);
l_sig=size(xval,2);
xMaxWAV = xval(end);
stepWAV = xval(2)-xval(1);
step1=stepWAV;
R_step=1/stepWAV;

j=1+floor(0:R_step/(2*ratio):l_sig);
f = (psi(j));
f_phi = (phi(j));

coefs(i,:) = sqrt(a)*wkeep1((wconv1(ySIG,f)),L1);
phi_coefs(i,:) = wkeep1((wconv1(ySIG,f_phi)),L1);

end

S1 = sum(coefs.*coefs,1);
coefs_sym1 = coefs./repmat(S1.^0.5,L1,1);

% save('coefs_sym1.mat','coefs_sym1')
% clear coefs_sym1 coefs S1

S2 = sum(phi_coefs.*phi_coefs,1);
phi_coefs_sym1 = phi_coefs./repmat(S2.^0.5,L1,1);

% save('phi_coefs_sym1.mat','phi_coefs_sym1')
% clear phi_coefs_sym1 phi_coefs S2

dict_sym1=[phi_coefs_sym1,coefs_sym1];
save([root '\dict_sym1_' num2str(L1) '.mat'],'dict_sym1')
clear coefs phi_coefs phi psi out1

for i=1:L1
    i
ySIG=zeros(1,L1);
ySIG(i)=1;

wname='sym4';
wtype = wavemngr('type',wname);
iter=10;
switch wtype
  case {1,3} , [phi,psi,xval] = wavefun(wname,iter);
  case 2     , [phi1,psi1,phi2,psi2,xval] = wavefun(wname,iter);
  case {4,5} , [psi,xval] = wavefun(wname,iter);
end

step = xval(2)-xval(1);
switch wtype
  case {1,3,4,5} , out1 = cumsum(psi)*step;
  case 2         , out1 = cumsum(psi1)*step; out3 = cumsum(psi2)*step;
end
 

xval = xval-xval(1);
l_sig=size(xval,2);
xMaxWAV = xval(end);
stepWAV = xval(2)-xval(1);
step1=stepWAV;
R_step=1/stepWAV;
j=1+floor(0:R_step/(2*ratio):l_sig);
f = (psi(j));
f_phi = (phi(j));

coefs(i,:) = sqrt(a)*wkeep1((wconv1(ySIG,f)),L1);
phi_coefs(i,:) = wkeep1((wconv1(ySIG,f_phi)),L1);

end

S1 = sum(coefs.*coefs,1);
coefs_sym4 = coefs./repmat(S1.^0.5,L1,1);


S2 = sum(phi_coefs.*phi_coefs,1);
phi_coefs_sym4 = phi_coefs./repmat(S2.^0.5,L1,1);

dict_sym4=[phi_coefs_sym4,coefs_sym4];
save([root '\dict_sym4_' num2str(L1) '.mat'],'dict_sym4')
clear coefs phi_coefs phi psi out1



for i=1:L1
    i
ySIG=zeros(1,L1);
ySIG(i)=1;

wname='dmey';
wtype = wavemngr('type',wname);
iter=10;
switch wtype
  case {1,3} , [phi,psi,xval] = wavefun(wname,iter);
  case 2     , [phi1,psi1,phi2,psi2,xval] = wavefun(wname,iter);
  case {4,5} , [psi,xval] = wavefun(wname,iter);
end

step = xval(2)-xval(1);
switch wtype
  case {1,3,4,5} , out1 = cumsum(psi)*step;
  case 2         , out1 = cumsum(psi1)*step; out3 = cumsum(psi2)*step;
end
 

xval = xval-xval(1);
l_sig=size(xval,2);
xMaxWAV = xval(end);
stepWAV = xval(2)-xval(1);
step1=stepWAV;
R_step=1/stepWAV;
j=1+floor(0:R_step/(2*ratio):l_sig);
f = (psi(j));
f_phi = (phi(j));

coefs(i,:) = sqrt(a)*wkeep1((wconv1(ySIG,f)),L1);
phi_coefs(i,:) = wkeep1((wconv1(ySIG,f_phi)),L1);

end

S1 = sum(coefs.*coefs,1);
coefs_dmey = coefs./repmat(S1.^0.5,L1,1);


S2 = sum(phi_coefs.*phi_coefs,1);
phi_coefs_dmey = phi_coefs./repmat(S2.^0.5,L1,1);

dict_dmey=[phi_coefs_dmey,coefs_dmey];
save([root '\dict_dmey_' num2str(L1) '.mat'],'dict_dmey')
clear coefs phi_coefs phi psi out1




x = 1:1:L1;
for k=1:125

dict_cos(:,k)=cos(k*x/L1);

end

S2 = sum(dict_cos.*dict_cos,1);
dict_cos = dict_cos./repmat(S2.^0.5,L1,1);

save([root '\dict_cos_' num2str(L1) '.mat'],'dict_cos')

x = 1:1:L1;
for k=1:125

dict_sin(:,k)=sin(k*x/L1);

end

S1 = sum(dict_sin.*dict_sin,1);
dict_sin = dict_sin./repmat(S1.^0.5,L1,1);

save([root '\dict_sin_' num2str(L1) '.mat'],'dict_sin')





















% for i=1:L1
% figure(1);plot(phi_coefs_db4(i,:)); 
% pause
% end
% 
% for i=1:L1
% figure(2);plot(coefs_db4(i,:)); 
% pause
% end

