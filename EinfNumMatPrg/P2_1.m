%% Aufgabe P2.1 - Gruppe 7 - Alexander Glock, Jannis Röder
clearvars
%A = [1 3 0; 3 2 6; 0 6 5]%[1 2 0; 4 5 6; 0 8 9]


%[L,R,P]=lr_zerlegung(A)
%C = cholesky_zerlegung(A)
Laplace_Mat(6)

%--------------------------------------------------------------------------
%---------------------  Laplace-Matrix  -----------------------------------

function A = Laplace_Mat(m)
    B = zeros(m)+diag(4*ones(1,m))+diag((-1)*ones(1,m-1),-1)+diag((-1)*ones(1,m-1),1);
    Id = eye(m);
    Z = zeros(m,m);
    
    % Null Matrix + Hauptdiagonale
    A = zeros(m^2,m^2)+diag(4*ones(1,m^2));
    % Nebendiagonalen mit Unterbrechungen bauen und addieren:
    n_diag=[0,ones(1,m-1)];
    dg=n_diag;
    for ind = 2:1:m
        dg=[dg,n_diag];
    end 
    %n_diag%diag(n_diag(2:end),-1)
    A=A+diag(dg(2:end),-1)+diag(dg(2:end),1);
    % Diagonalen der Einheitsmatrizen addieren:
    A=A+diag(ones(1,m^2-m),m)+diag(ones(1,m^2-m),-m);
end

%--------------------------------------------------------------------------
%---------------------  LR-Zerlegung    -----------------------------------

function [L,R,P] = lr_zerlegung(A)
    % Thomas Algorithmus
    n = size(A,1);
    a = diag(A)';
    b = [0, diag(A,-1)'];
    c = diag(A, 1)';
    
    alpha=zeros(1,n)';
    beta=zeros(1,n)';
    
    alpha(1)=a(1);
    jj=(1:1:n);
    for j=jj(2:end)
        beta(j)=(b(j)/alpha(j-1));
        alpha(j)=a(j)-beta(j)*c(j-1);
    end
    
    L = eye(n)+diag(b(2:end),-1)
    R = zeros(n)+diag(alpha)+diag(c,1)
    P = zeros(n)
end

%--------------------------------------------------------------------------
%---------------------  cholesky-Zerlegung    -----------------------------

function C = cholesky_zerlegung(A)
n= size(A,1);

C = zeros(n,n);
    for i = 1:1:n
        for j = 1:1:i
            if (i == j)
                C1 = 0;
                for k = 1:1:(i-1)
                    C1 = C1 + (C(i,k))^2;
                end
                C(i,i) = sqrt( A(i,i) -  C1 );
             else
                C2 = 0;
                for r = 1:1:(i-1)
                    C2 = C2 + C(i,r)*C(j,r);
                end
                C(i,j) = (A(i,j) - C2)*1/C(j,j);
             end
        end
    end 
end