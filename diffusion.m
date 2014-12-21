% transforms concentrations and repulsion potentials to find diffusion
%
% Author: Eran Agmon, agmon.eran@gmail.com
% Affilitation: Cognitive Science Program and Informatics Department,
% Indiana University
% Last updated: 10/15/2014

function [p_move]=diffusion(conc,pot) %concentration and potential

    global nx
    
    %preallocate
    ge_out=zeros(nx,nx,8);
    ge_in=zeros(nx,nx,8);
    deltP=zeros(nx,nx,8);
    
    % perform circular shift on potentials and concentration
    N_pot  = [pot(2:nx,:) ;      pot(1,:)];
    S_pot  = [pot(nx,:) ;   pot(1:nx-1,:)];
    W_pot  = [pot(:,nx)     pot(:,1:nx-1)];
    E_pot  = [pot(:,2:nx)        pot(:,1)];
    NE_pot = [N_pot(:,2:nx)    N_pot(:,1)];
    NW_pot = [N_pot(:,nx) N_pot(:,1:nx-1)];
    SE_pot = [S_pot(:,2:nx)    S_pot(:,1)];
    SW_pot = [S_pot(:,nx) S_pot(:,1:nx-1)];

    N_conc  = [conc(2:nx,:) ;      conc(1,:)];
    W_conc  = [conc(:,nx)     conc(:,1:nx-1)];
    S_conc  = [conc(nx,:) ;   conc(1:nx-1,:)];
    E_conc  = [conc(:,2:nx)        conc(:,1)];   
    NE_conc = [N_conc(:,2:nx)    N_conc(:,1)];
    NW_conc = [N_conc(:,nx) N_conc(:,1:nx-1)];
    SE_conc = [S_conc(:,2:nx)    S_conc(:,1)];
    SW_conc = [S_conc(:,nx) S_conc(:,1:nx-1)];
    
    % outflow
    deltP(:,:,1) = pot - NE_pot;  
    deltP(:,:,2)  = pot - N_pot;
    deltP(:,:,3) = pot - NW_pot; 
    deltP(:,:,4)  = pot - W_pot;
    deltP(:,:,5) = pot - SW_pot;
    deltP(:,:,6)  = pot - S_pot;
    deltP(:,:,7) = pot - SE_pot;
    deltP(:,:,8)  = pot - E_pot;
    
    ge_out = (deltP)./(1-exp(-deltP)); 
    ge_out(isnan(ge_out)) = 1;  %remove NANs here, from when deltaP=0
    ge_out(isinf(ge_out)) = 1;
    
    % multiply by concentrations to get total outflow
    outflow = conc.*sum(ge_out,3);

    
    % inflow use -delta to reverse direction of flow
    ge_in = (-deltP)./(1-exp(deltP)); 
    ge_in(isnan(ge_in)) = 1;
    ge_in(isinf(ge_in)) = 1;
    
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
    

