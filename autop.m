% find derivative
%
% Author: Eran Agmon, agmon.eran@gmail.com
% Affilitation: Cognitive Science Program and Informatics Department,
% Indiana University
% Last updated: 10/15/2014

function zp = autop(t,z)

global Diff DiffF ff1 ff2 ka km nx attraction Frefill flowF saturF

%reformat 
zp=z; 
MM   =reshape(z(1:nx^2,1)            ,nx,nx);
AA   =reshape(z((nx^2+1):(2*nx^2))   ,nx,nx);
FF   =reshape(z((2*nx^2+1):(3*nx^2)) ,nx,nx);
WW   =reshape(z((3*nx^2+1):(4*nx^2)) ,nx,nx);
DIR  =reshape(z((4*nx^2+1):(5*nx^2)) ,nx,nx);

%repulsion
[M_d,A_d,W_d]=potential_ellipse(MM,AA,WW,DIR); %no self-repulsion

%diffusion
M_d2=Diff*diffusion(MM,M_d);
A_d2=Diff*diffusion(AA,A_d);
W_d2=Diff*diffusion(WW,W_d);

%diffusion without potential
F_d2=DiffF*diffusion_F(FF);

%alignment
dd=attraction*align(MM,DIR);

%reactions
dm = M_d2 + ff2*FF.*AA.*AA - km*MM;
da = A_d2 + ff1*FF.*AA.*AA - ka*AA;
df = F_d2 - ff1*FF.*AA.*AA - ff2*FF.*AA.*AA;
dw = W_d2;

% food refill
df(Frefill) = df(Frefill) + flowF*(saturF-FF(Frefill));

%vectorize output
zp(1:nx^2)              =reshape(dm,nx^2,1);
zp((nx^2+1):(2*nx^2))   =reshape(da,nx^2,1);
zp((2*nx^2+1):(3*nx^2)) =reshape(df,nx^2,1);
zp((3*nx^2+1):(4*nx^2)) =reshape(dw,nx^2,1);
zp((4*nx^2+1):(5*nx^2)) =reshape(dd,nx^2,1);


