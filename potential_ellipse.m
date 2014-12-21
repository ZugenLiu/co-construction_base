% calculate total repulsion based on repulsion ellipse
%
% Author: Eran Agmon, agmon.eran@gmail.com
% Affilitation: Cognitive Science Program and Informatics Department,
% Indiana University
% Last updated: 10/15/2014

function [Mp,Ap,Wp]=potential_ellipse(MM,AA,WW,DIR)

    global nx long short neighbor_rep samesite_rep self_rep_A self_rep_M
    
    %preallocate
    Mp=zeros(nx);    %membrane potential    
    Wp=zeros(nx);    %water potential
    Ap=zeros(nx);    %water potential
    
    Ar=zeros(nx,nx,4);  %area
    
    %find area
    a1 = long*short/2 * (  pi/8-DIR - atan((short-long)*sin(2*(  pi/8-DIR)) ./ (short+long + (short-long)*cos(2*(  pi/8-DIR)))));
    a2 = long*short/2 * (3*pi/8-DIR - atan((short-long)*sin(2*(3*pi/8-DIR)) ./ (short+long + (short-long)*cos(2*(3*pi/8-DIR)))));
    a3 = long*short/2 * (5*pi/8-DIR - atan((short-long)*sin(2*(5*pi/8-DIR)) ./ (short+long + (short-long)*cos(2*(5*pi/8-DIR)))));
    a4 = long*short/2 * (7*pi/8-DIR - atan((short-long)*sin(2*(7*pi/8-DIR)) ./ (short+long + (short-long)*cos(2*(7*pi/8-DIR)))));
    a5 = long*short/2 * (9*pi/8-DIR - atan((short-long)*sin(2*(9*pi/8-DIR)) ./ (short+long + (short-long)*cos(2*(9*pi/8-DIR)))));   
    
    Ar(:,:,1)=(a2-a1)/(long*short*pi);
    Ar(:,:,2)=(a3-a2)/(long*short*pi);
    Ar(:,:,3)=(a4-a3)/(long*short*pi);
    Ar(:,:,4)=(a5-a4)/(long*short*pi);
    
    Ar_N  = [Ar(2:nx,:,:) ;    Ar(1,:,:)];
    Ar_W  = [Ar(:,nx,:)   Ar(:,1:nx-1,:)];
    Ar_S  = [Ar(nx,:,:) ; Ar(1:nx-1,:,:)];
    Ar_E  = [Ar(:,2:nx,:)      Ar(:,1,:)];
    Ar_NE = [Ar_N(:,2:nx,:)    Ar_N(:,1,:)];
    Ar_NW = [Ar_N(:,nx,:) Ar_N(:,1:nx-1,:)];
    Ar_SE = [Ar_S(:,2:nx,:)    Ar_S(:,1,:)];
    Ar_SW = [Ar_S(:,nx,:) Ar_S(:,1:nx-1,:)];
    
    M_N  = [MM(2:nx,:) ;    MM(1,:)];
    M_W  = [MM(:,nx)   MM(:,1:nx-1)];
    M_S  = [MM(nx,:) ; MM(1:nx-1,:)];
    M_E  = [MM(:,2:nx)      MM(:,1)];
    M_NE = [M_N(:,2:nx)    M_N(:,1)];
    M_NW = [M_N(:,nx) M_N(:,1:nx-1)];
    M_SE = [M_S(:,2:nx)    M_S(:,1)];
    M_SW = [M_S(:,nx) M_S(:,1:nx-1)];
    
    A_N  = [AA(2:nx,:) ;    AA(1,:)];
    A_W  = [AA(:,nx)   AA(:,1:nx-1)];
    A_S  = [AA(nx,:) ; AA(1:nx-1,:)];
    A_E  = [AA(:,2:nx)      AA(:,1)];
    A_NE = [A_N(:,2:nx)    A_N(:,1)];
    A_NW = [A_N(:,nx) A_N(:,1:nx-1)];
    A_SE = [A_S(:,2:nx)    A_S(:,1)];
    A_SW = [A_S(:,nx) A_S(:,1:nx-1)];
    
    W_N  = [WW(2:nx,:) ;    WW(1,:)];
    W_W  = [WW(:,nx)   WW(:,1:nx-1)];
    W_S  = [WW(nx,:) ; WW(1:nx-1,:)];
    W_E  = [WW(:,2:nx)      WW(:,1)];
    W_NE = [W_N(:,2:nx)    W_N(:,1)];
    W_NW = [W_N(:,nx) W_N(:,1:nx-1)];
    W_SE = [W_S(:,2:nx)    W_S(:,1)];
    W_SW = [W_S(:,nx) W_S(:,1:nx-1)];
    
    
    Ap = neighbor_rep*(M_NE.*Ar_NE(:,:,1) + M_N.*Ar_N(:,:,2) + M_NW.*Ar_NW(:,:,3) + M_W.*Ar_W(:,:,4)  ...
                     + M_SW.*Ar_SW(:,:,1) + M_S.*Ar_S(:,:,2) + M_SE.*Ar_SE(:,:,3) + M_E.*Ar_E(:,:,4) + samesite_rep*MM) ...
        + self_rep_A*((A_NE + A_N + A_NW + A_W + A_SW + A_S + A_SE + A_E)/8 + samesite_rep*AA);
            
    
    Wp = neighbor_rep*(M_NE.*Ar_NE(:,:,1) + M_N.*Ar_N(:,:,2) + M_NW.*Ar_NW(:,:,3) + M_W.*Ar_W(:,:,4)  ...
                     + M_SW.*Ar_SW(:,:,1) + M_S.*Ar_S(:,:,2) + M_SE.*Ar_SE(:,:,3) + M_E.*Ar_E(:,:,4) + samesite_rep*MM);
            
                 
%     Mp = neighbor_rep*((W_NE+A_NE).*Ar(:,:,1) + (W_N+A_N).*Ar(:,:,2) + (W_NW+A_NW).*Ar(:,:,3) + (W_W+A_W).*Ar(:,:,4)  ...
%                      + (W_SW+A_SW).*Ar(:,:,1) + (W_S+A_S).*Ar(:,:,2) + (W_SE+A_SE).*Ar(:,:,3) + (W_E+A_E).*Ar(:,:,4) + samesite_rep*(WW+AA)) ...
%         + self_rep_M*(M_NE.*(Ar_NE(:,:,1)+Ar(:,:,1))/2 + M_N.*(Ar_N(:,:,2)+Ar(:,:,2))/2 + M_NW.*(Ar_NW(:,:,3)+Ar(:,:,3))/2 + M_W.*(Ar_W(:,:,4)+Ar(:,:,4))/2 ...
%                   + M_SW.*(Ar_SW(:,:,1)+Ar(:,:,1))/2 + M_S.*(Ar_S(:,:,2)+Ar(:,:,2))/2 + M_SE.*(Ar_SE(:,:,3)+Ar(:,:,3))/2 + M_E.*(Ar_E(:,:,4)+Ar(:,:,4))/2 + samesite_rep*MM);
    
    Mp = neighbor_rep*((W_NE+A_NE).*Ar(:,:,1) + (W_N+A_N).*Ar(:,:,2) + (W_NW+A_NW).*Ar(:,:,3) + (W_W+A_W).*Ar(:,:,4)  ...
        + (W_SW+A_SW).*Ar(:,:,1) + (W_S+A_S).*Ar(:,:,2) + (W_SE+A_SE).*Ar(:,:,3) + (W_E+A_E).*Ar(:,:,4) + samesite_rep*(WW+AA)) ...
        + self_rep_M*((M_NE + M_N + M_NW + M_W + M_SW + M_S + M_SE + M_E)/8 + samesite_rep*MM);
            
      

end




