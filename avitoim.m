clc
clear all
close all

folderName='1';

root='E:\gifani\CS';
mkdir([root '\' folderName])


xyloObj = VideoReader([root '\' folderName '.avi']);
vidFrames = read(xyloObj);
nFrames = xyloObj.NumberOfFrames;


for k = 1 : nFrames
    mov(k).cdata = read(xyloObj, k);
end

for kk=1:nFrames
    u=(rgb2gray(mov(kk).cdata));
%     u1=u(71:300,126:260);
    imwrite(u,[root '\' folderName '\im' num2str(kk) '.bmp'])
    
end

