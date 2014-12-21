%creates a size*size matrix with gaussian function centered at size/2,size/2                    
function blur = gauss2dconvolve(matrix,sigma2)

sigma=sigma2^0.5;

G=fspecial('gaussian',[101 101],sigma);

% G=fspecial('gaussian',size(matrix),sigma);
blur = imfilter(matrix,G,'same');


% size=40;
% 
% gauss2d     = zeros(size);
% mid = round(size/2);
% 
% for i=1:size
%     for j=1:size
%         x=(i^2+j^2)^0.5;
%         gauss2d(i,j) = exp(-(  (j-mid)^2/(2*sigma)  +  (i-mid)^2/(2*sigma)  ));
%     end
% end
% 
% 
% blur = conv2(matrix,gauss2d);