
%   transforms orientation based on neighbor's orientation and
%   concentration

function [d_theta]=align(conc,dir) %concentration and direction
    global nx
    
    % perform circular shift 
    N_dir  = [dir(2:nx,:) ;      dir(1,:)];
    W_dir  = [dir(:,nx)     dir(:,1:nx-1)];
    S_dir  = [dir(nx,:) ;   dir(1:nx-1,:)];
    E_dir  = [dir(:,2:nx)        dir(:,1)];
    NE_dir = [N_dir(:,2:nx)    N_dir(:,1)];
    NW_dir = [N_dir(:,nx) N_dir(:,1:nx-1)];
    SE_dir = [S_dir(:,2:nx)    S_dir(:,1)];
    SW_dir = [S_dir(:,nx) S_dir(:,1:nx-1)]; 
    
    N_conc  = [conc(2:nx,:) ;      conc(1,:)];
    W_conc  = [conc(:,nx)     conc(:,1:nx-1)];
    S_conc  = [conc(nx,:) ;   conc(1:nx-1,:)];
    E_conc  = [conc(:,2:nx)        conc(:,1)];
    NE_conc = [N_conc(:,2:nx)    N_conc(:,1)];
    NW_conc = [N_conc(:,nx) N_conc(:,1:nx-1)];
    SE_conc = [S_conc(:,2:nx)    S_conc(:,1)];
    SW_conc = [S_conc(:,nx) S_conc(:,1:nx-1)]; 
    
    N_psi  = sin(2*(N_dir - dir));
    W_psi  = sin(2*(W_dir - dir));
    S_psi  = sin(2*(S_dir - dir));
    E_psi  = sin(2*(E_dir - dir));
    NE_psi = sin(2*(NE_dir - dir));
    NW_psi = sin(2*(NW_dir - dir));
    SE_psi = sin(2*(SE_dir - dir));
    SW_psi = sin(2*(SW_dir - dir));  
           

    d_theta = N_psi.*N_conc   + W_psi.*W_conc   + S_psi.*S_conc   + E_psi.*E_conc ...
            + NE_psi.*NE_conc + NW_psi.*NW_conc + SE_psi.*SE_conc + SW_psi.*SW_conc;
  
end




