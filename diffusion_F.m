% finds diffusion for food molecule, which has no repulsion
%
% Author: Eran Agmon, agmon.eran@gmail.com
% Affilitation: Cognitive Science Program and Informatics Department,
% Indiana University
% Last updated: 10/15/2014

function [p_move]=diffusion_F(conc) %concentration and potential
    global nx
    
    %preallocate - make all potential differences one
    ge_out=ones(nx,nx,8);
    ge_in=ones(nx,nx,8);

    % perform circular shift on potentials and concentration
    N_conc  = [conc(2:nx,:) ;      conc(1,:)];
    W_conc  = [conc(:,nx)     conc(:,1:nx-1)];
    S_conc  = [conc(nx,:) ;   conc(1:nx-1,:)];
    E_conc  = [conc(:,2:nx)        conc(:,1)];   
    NE_conc = [N_conc(:,2:nx)    N_conc(:,1)];
    NW_conc = [N_conc(:,nx) N_conc(:,1:nx-1)];
    SE_conc = [S_conc(:,2:nx)    S_conc(:,1)];
    SW_conc = [S_conc(:,nx) S_conc(:,1:nx-1)];

    
    % multiply by concentrations to get total outflow
    outflow = conc.*sum(ge_out,3);

    % multiply by concentrations to get total inflow
    ge_in(:,:,1) = ge_in(:,:,1).* NE_conc;  
    ge_in(:,:,2) = ge_in(:,:,2).*  N_conc;   
    ge_in(:,:,3) = ge_in(:,:,3).* NW_conc;    
    ge_in(:,:,4) = ge_in(:,:,4).*  W_conc;   
    ge_in(:,:,5) = ge_in(:,:,5).* SW_conc;   
    ge_in(:,:,6) = ge_in(:,:,6).*  S_conc;   
    ge_in(:,:,7) = ge_in(:,:,7).* SE_conc;   
    ge_in(:,:,8) = ge_in(:,:,8).*  E_conc;    
    
    inflow = sum(ge_in,3);
    
    % total flow = inflow - outflow
    p_move = inflow - outflow;
    

