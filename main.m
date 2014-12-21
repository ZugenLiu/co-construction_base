% Run simulation
%
% Author: Eran Agmon, agmon.eran@gmail.com
% Affilitation: Cognitive Science Program and Informatics Department,
% Indiana University
% Last updated: 10/15/2014

function main() 

clear all
close all

global KEY_IS_PRESSED
global nx Diff DiffF ff1 ff2 ka km flowF saturF Frefill 
global long short neighbor_rep samesite_rep self_rep_A self_rep_M attraction

KEY_IS_PRESSED = 0;
nx=40;          %set spatial size

Diff = 0.25;   DiffF = Diff;    %diffusion
ff1  = 0.25;%0.25 
ff2  = 0.3;    %forward reaction rates
ka   = 0.005; 
km   = 0.005;  %decay

flowF=0.8; saturF=0.18; %food parameters

long=20; short=1;
neighbor_rep=3;
samesite_rep=3;     %repulsion at same-site (a fraction of neighbor rep)
self_rep_A=0.1;     %repulsion between same molecular type
self_rep_M=0.1;    %repulsion between same molecular type

attraction=0.01;%0.1;     %for alignment

%% initialization

fD=15;      %food refill distance

%food refill
Frefill=ones(nx);
F_circ=zeros(nx);
for r=0:0.01:fD   
    F_circ=drawcircle(F_circ,ceil(nx/2),ceil(nx/2),r);
end
Frefill=Frefill-F_circ;
Frefill = Frefill > 0.5;

load('SC')
m0=MM;a0=AA;f0=FF;w0=WW;d0=DIR; 


%% more set up

%vectorize
MM=reshape(m0,nx^2,1); AA=reshape(a0,nx^2,1); FF=reshape(f0,nx^2,1); WW=reshape(w0,nx^2,1); DIR=reshape(d0,nx^2,1); 
z0=[MM;AA;FF;WW;DIR];

Mvec=[sum(sum(MM))/(nx*nx)];Avec=[sum(sum(AA))/(nx*nx)];Wvec=[sum(sum(WW))/(nx*nx)];Fvec=[sum(sum(FF))/(nx*nx)]; Tvec=[0];
dMvec=[0];dAvec=[0];dWvec=[0];dFvec=[0]; 


qX=zeros(nx*nx,1); qY=zeros(nx*nx,1);
ind=1;
for xx=1:nx
    for yy=1:nx
        qX(ind)=xx;
        qY(ind)=yy;        
        ind=ind+1;
    end
end

    dir=reshape(DIR,nx^2,1); mc=reshape(MM,nx^2,1);
    qU=mc.*cos(dir); qV=mc.*sin(dir);

%% integration and graphics

perturb=20.0;   %set perturbation time
sigma2=0.5;
alpha =2;
xloc  =16;
yloc  =20;


zz=z0; 
time=0; tt=10; %start time, timevisualize
fig = figure; set(fig,'Position', [200, 200, 1600, 800]); set(gcf, 'KeyPressFcn', @myKeyPressFcn)
while ~KEY_IS_PRESSED
    
    if time==perturb
        m_soln = z(size(z,1),1:nx^2);
        a_soln = z(size(z,1),(nx^2+1):(2*nx^2));
        f_soln = z(size(z,1),(2*nx^2+1):(3*nx^2));
        w_soln = z(size(z,1),(3*nx^2+1):(4*nx^2));
        d_soln = z(size(z,1),(4*nx^2+1):(5*nx^2));

        MM =reshape(m_soln,nx,nx); 
        AA =reshape(a_soln,nx,nx); 
        FF =reshape(f_soln,nx,nx);
        WW =reshape(w_soln,nx,nx); 
        DIR=reshape(d_soln,nx,nx); 
    
        
%         MM = 0.5*MM;

          MM = localpert(MM,sigma2,alpha,xloc,yloc);
          AA = localpert(AA,sigma2,alpha,xloc,yloc);
          WW = localpert(WW,sigma2,alpha,xloc,yloc);
          FF = localpert(FF,sigma2,alpha,xloc,yloc);

        MMP=reshape(MM,nx^2,1); AAP=reshape(AA,nx^2,1); FFP=reshape(FF,nx^2,1); WWP=reshape(WW,nx^2,1); DIRP=reshape(DIR,nx^2,1); 
        zz=[MMP;AAP;FFP;WWP;DIRP];         
    end
    
    
    [t,z] = ode45(@autop, [0 tt], zz); %Runge-Kutta one-step solver
    
    time=time+tt;
    zz=z(size(z,1),:);
        
    m_soln = z(size(z,1),1:nx^2);
    a_soln = z(size(z,1),(nx^2+1):(2*nx^2));
    f_soln = z(size(z,1),(2*nx^2+1):(3*nx^2));
    w_soln = z(size(z,1),(3*nx^2+1):(4*nx^2));
    d_soln = z(size(z,1),(4*nx^2+1):(5*nx^2));
    
    Mvec=[Mvec sum(sum(MM))/(nx*nx)];Avec=[Avec sum(sum(AA))/(nx*nx)];
    Wvec=[Wvec sum(sum(WW))/(nx*nx)];Fvec=[Fvec sum(sum(FF))/(nx*nx)];
    Tvec=[Tvec time];
    
    MM =reshape(m_soln,nx,nx); 
    AA =reshape(a_soln,nx,nx); 
    FF =reshape(f_soln,nx,nx);
    WW =reshape(w_soln,nx,nx); 
    DIR=reshape(d_soln,nx,nx); 
    
    
    qU=cos(d_soln)'; qV=sin(d_soln)';
    
    subplot(2,4,1); surf(AA); view(2); caxis([0 1]); xlabel('autocatalyst','FontSize',20); shading flat; %colorbar
    set(gca,'XTick',[],'YTick',[])
    str = sprintf('time = %f', time);
    title(str,'FontSize',40)

    subplot(2,4,2); surf(MM); view(2); caxis([0 1]); xlabel('membrane','FontSize',20); shading flat; %colorbar
    set(gca,'XTick',[],'YTick',[])

    subplot(2,4,3); surf(WW); view(2); caxis([0 1]); xlabel('water','FontSize',20); shading flat; %colorbar
    set(gca,'XTick',[],'YTick',[])

    subplot(2,4,4); surf(FF); view(2); caxis([0 1]); xlabel('food','FontSize',20); shading flat; %colorbar
    set(gca,'XTick',[],'YTick',[])

    subplot(2,4,5); quiver(qX,qY,qU,qV)
    xlim([0 nx]); ylim([0 nx]); xlabel('orientation','FontSize',20);
    set(gca, 'CLim', [0 1],'xlim',[1 nx],'ylim',[1 nx],'xtick',[],'ytick',[]); 
    

    subplot(2,4,6); plot(Tvec,Mvec,'k','LineWidth',2); hold on; 
    plot(Tvec,Avec,'r','LineWidth',2); plot(Tvec,Wvec,'b','LineWidth',2);plot(Tvec,Fvec,'g','LineWidth',2); hold off
    xlabel('time','FontSize',20); ylabel('total concentration','FontSize',20);legend('Membrane','Autocatalyst','Water','Food');
    ylim([0 0.15]); 
    
    drawnow

end
disp('loop ended')


function myKeyPressFcn(hObject, event)
global KEY_IS_PRESSED
KEY_IS_PRESSED  = 1;
disp('key is pressed')
