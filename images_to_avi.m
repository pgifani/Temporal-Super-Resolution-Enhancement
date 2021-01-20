clc
% clear all

froot='\81-162';
list = dir(sprintf('%s\\*.bmp', froot));
num_im=length(list);
[xi ,yi, zi]=size((imread([ froot '\im' num2str(1) '.bmp'])));


movie=uint8(zeros(xi, yi,zi,num_im));
for frame=1:num_im
[movie(:,:,frame),map] = (imread([ froot '\im' num2str(frame) '.bmp']));
end
mov = immovie(movie,map);
%
movie2avi(mov, [froot '\81-162.avi'], 'compression', 'None','fps',60);
