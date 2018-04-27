clear all
image_files = dir(fullfile('condor\','*.png'));
n = length(image_files); % number of samples
[h,w,~] = size(imread(fullfile('condor\',image_files(1).name)));
J = zeros(h,w,n);
for i = 1:n
    X = imread(fullfile('condor\',image_files(i).name));
    J(:,:,i) = rgb2gray(X);
    %J(:,:,i) = poissrnd(double(X)); % original image with poisson noise
end
set = J;