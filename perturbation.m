clear all
close all

addpath('/Users/eagmon/Desktop/multicore')


% localized gaussian perturbation
savdir = '/Users/eagmon/Dropbox/Alife special edition/autopoiesis_9-26/local_gauss/localgauss1.mat';
resultCell=[0 0 0 0 0 0 0];
save(savdir, 'resultCell');

sigma2 =[0:1:5];
alpha  =[-1:1:4];
xloc = 20;
yloc = 15;
nSigma=length(sigma2);
nAlpha=length(alpha);


parameterCell = cell(1, nSigma*nAlpha);
for k = 1:nSigma
    for j = 1:nAlpha
        n=(k-1)*nSigma+j;
      parameterCell{1,n} = {sigma2(k), alpha(j), xloc, yloc, savdir};
    end
end

resultCell = startmulticoremaster(@localgaussian, parameterCell);




addpath('/Users/eagmon/Dropbox/Alife special edition/autopoiesis_9-26')
addpath('/Users/eagmon/Desktop/multicore')
startmulticoreslave;


% % gaussian blur perturbation
% savdir = '/Users/eagmon/Desktop/AlifeModel_9-24/blur/blur4.mat';
% resultCell=[0 0 0 0];
% save(savdir, 'resultCell');
% 
% sigma2=[0.05:0.05:3];
% nrOfEvals=length(sigma2);
% 
% 
% parameterCell = cell(1, nrOfEvals);
% for k = 1:nrOfEvals
%   parameterCell{1,k} = {sigma2(k), savdir};
% end
% 
% resultCell = startmulticoremaster(@gaussianblur, parameterCell);


% for k = 1:numel(parameterCell) 
%     parameterCell(k)
% %     resultCell(k,:) = gaussianblur(parameterCell(k)); 
%   [sigma2,timetostable,change,viable] = gaussianblur(parameterCell(k)); 
%   resultCell(k,:)=[sigma2 timetostable change viable];
%   
%   %save to file
%   save(fullfile(savdir,str),'resultCell','-append')
% end
%     
    


