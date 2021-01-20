function TSR_signal=TSR_Gifani(wave_dic,mydict2,mydict1_wave,mydict1,Sparse_Sig,T1,T2)

L1=T1;
L2=T2;
x_BCS=Sparse_Sig;
ratio=L2/L1;

sig1=zeros(L1,1);
for im=1:size(mydict1,2)
    
    sig=x_BCS(im).*mydict1(:,im);
    sig1=sig+sig1;
    
end

[mm,nn]=find(x_BCS);

new_mm=round(mm*ratio);

value_mm=x_BCS(mm);
BCS_new=zeros(size(x_BCS,1)*ratio,1);
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

% sparse_coefs=x_BCS(m);
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
TSR_signal=final*ratio2;

