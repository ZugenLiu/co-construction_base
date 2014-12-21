%creates a size*size matrix with gaussian function centered at size/2,size/2                    
function mat = localpert(matrix,sigma2,alpha,xloc,yloc)

if ~(sigma2==0) && ~(alpha==0)

    size=length(matrix);
    sigma=sigma2^0.5;
    gauss = zeros(size);

    for i=1:size
        for j=1:size
            x=(i^2+j^2)^0.5;
            gauss(i,j) = alpha*exp(-(  (j-xloc)^2/(2*sigma)  +  (i-yloc)^2/(2*sigma)  )) + 1;
        end
    end

    mat=matrix.*gauss;

else
    mat=matrix;
end