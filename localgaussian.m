function timetostable = localgaussian(parameter)
close all

param=cell2mat(parameter);
sigma2=param(1);
alpha =param(2);
xloc  =param(3);
yloc  =param(4);

global savdir



global nx Diff DiffF ff1 ff2 ka km flowF saturF Frefill 
global long short neighbor_rep samesite_rep self_rep_A self_rep_M attraction

nx=40;          %set spatial size
Diff = 0.25;   DiffF = Diff;    %diffusion
ff1  = 0.25; ff2  = 0.3; ka   = 0.005; km   = 0.005;  %decay
flowF=0.8; saturF=0.18; %food parameters
long=20; short=1;
neighbor_rep=3; samesite_rep=3; %repulsion at same-site (a fraction of neighbor rep)
self_rep_A=0.1; self_rep_M=0.1; %repulsion between same molecular type

attraction=0.01;%0.1;     %for alignment

%% initialization

avgM=0.9;%0.8;  %initial concentration of M
avgA=0.9;  %initial concentration of A

mR=3;%3;      %membraneRadius
mT=4;%4;      %membraneThickness
fD=8;%6;%4;      %food refill distance

%food refill
Frefill=ones(nx);
F_circ=zeros(nx);
for r=0:0.01:mR+mT+fD   
    F_circ=drawcircle(F_circ,ceil(nx/2),ceil(nx/2),r);
end
Frefill=Frefill-F_circ;
Frefill = Frefill > 0.5;


load('f9')
m0=MM;a0=AA;a0=AA;f0=FF;w0=WW;d0=DIR; 


%% more set up

%vectorize
MM=reshape(m0,nx^2,1); AA=reshape(a0,nx^2,1); FF=reshape(f0,nx^2,1); WW=reshape(w0,nx^2,1); DIR=reshape(d0,nx^2,1); 
z0=[MM;AA;FF;WW;DIR];
z0=z0';

%% integration and graphics

tt=10;          %integration length
epsilon=0.05;    %epsilon constant, scaled to the sample timestep
tk=50;%5000;        %how long it must remain stable
perturb=20.0;   %set perturbation time

zz=z0; zinitial=z0;
difference=0;
time=0;         %initialize time
tstable=0;      %initialize stability length counter
timetostable=0; %initialize time to stabilization counter
viable=1;       %it starts of viable;


epsilon=epsilon*2/tt;
% savefig=[100:100:200000];%[10:10:1000];
% fig = figure; set(fig,'Position', [200, 200, 500, 500]);
% figure
while tstable<=tk && viable==1
    
    [t,z] = ode45(@autop, [0 tt], zz);  %Runge-Kutta one-step solver
    
    time=time+tt;
    zz=z(size(z,1),:);
    

    difference=sum(abs(zz(1:(2*nx^2))-z0(1:(2*nx^2)))); %difference between membrane and autocatalyst 
    z0=zz;	%reset
    
    if difference<epsilon && tstable>0         %if the system is still stable
        tstable=tstable+tt;
    elseif difference<epsilon && tstable==0    %if the system is stable for the first time
        tstable=tstable+tt; timetostable=time;
    elseif difference>epsilon && tstable>0     %if the system destabilized
        tstable=0;
    end
        
    m_soln = z(size(z,1),1:nx^2); 
    a_soln = z(size(z,1),(nx^2+1):(2*nx^2));
    f_soln = z(size(z,1),(2*nx^2+1):(3*nx^2)); 
    w_soln = z(size(z,1),(3*nx^2+1):(4*nx^2));
    d_soln = z(size(z,1),(4*nx^2+1):(5*nx^2));
    
    
    if sum(m_soln+a_soln)<10
        viable=0;       %cell has died
        tstable=0;
    end
    
    MM =reshape(m_soln,nx,nx); 
    AA =reshape(a_soln,nx,nx); 
    FF =reshape(f_soln,nx,nx);
    WW =reshape(w_soln,nx,nx); 
    DIR=reshape(d_soln,nx,nx); 
    
    if time==perturb

          MM = localpert(MM,sigma2,alpha,xloc,yloc);
          AA = localpert(AA,sigma2,alpha,xloc,yloc);
          WW = localpert(WW,sigma2,alpha,xloc,yloc);
          FF = localpert(FF,sigma2,alpha,xloc,yloc);

        MMP=reshape(MM,nx^2,1); AAP=reshape(AA,nx^2,1); FFP=reshape(FF,nx^2,1); WWP=reshape(WW,nx^2,1); DIRP=reshape(DIR,nx^2,1); 
        zz=[MMP;AAP;FFP;WWP;DIRP];         
    end
    
    
%     surf(MM); view(2); caxis([0 1]); xlabel('membrane','FontSize',20); shading flat; %colorbar
%     set(gca,'XTick',[],'YTick',[])
%     strt = sprintf('t = %f, sig2 = %f, tstable = %f, dsum = %f', time, sigma2, tstable, difference);
%     title(strt,'FontSize',15)
%     drawnow
    
end

change=sum(abs(zz-zinitial));%z-zinitial;

% results=[timetostable change];
% DateString=datestr(clock);
% save(['var%f.mat', 'sigma2'],'MM','AA','WW','FF','DIR');

results=[sigma2 alpha xloc yloc timetostable change viable];

load(savdir);
resultCell=[resultCell;results];
save(savdir, 'resultCell', '-append');




