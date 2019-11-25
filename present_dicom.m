function present_dicom(folder, date)

%clc; clear all;

% info =  dicominfo('test1/test1.20191119.1.dcm')
% [x, cmap] = dicomread(info)
%  montage(x, cmap, 'Size', [64 64])
% % 
% fprintf('go')
% folder = 'test1';
% files = dir(folder);
% names = strings(1,length(files)-2);
% for k=3:length(files)
%     names(1,k-2) = files(k).name;
%     k = k+1;
% end
% names

%% Print all dicom files in succession to view the sideview of an artery
%folder = 'test1';
%date = 20191120;
files = dir(folder);
for k=1:(length(files)-3)
    info = dicominfo(sprintf('ba_normal/ba_normal.%0.f.%0.f.dcm', date, k));
    Y = dicomread(info);
    figure(1)
    imshow(Y,[], 'InitialMagnification', 'fit');
    truesize([500 500]);
    pause(0.5)
end


% %% Print all dicom images in a singular figure
% folder = 'test1';
% files = dir(folder);
% figure(2)
% for k=1:(length(files)-3)
%     info = dicominfo(sprintf('test1/test1.20191120.%0.f.dcm',k));
%     Y = dicomread(info);
%     %subplot(4,(length(files)-2)/4,k)
%     imshow(Y,[], 'InitialMagnification', 'fit');
%     truesize([500 500]);
%     pause(0.5)
% end


% info = dicominfo('test1/test1.20191119.1.dcm');
% Y = dicomread(info);
% figure(1)
% imshow(Y,[]);

end