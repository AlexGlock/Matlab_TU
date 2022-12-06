%% Aufgabe P3.1 - Gruppe 7 - Alexander Glock, Jannis Röder
clearvars

ausgleich_stabilitaet()

%---------------------  Funktionen-Tests  ---------------------------------

%R = [1 2 3 4; 0 5 6 7; 0 0 8 9; 0 0 0 10];
%L = R';
%b = [11 12 13 14]';
%x=rueckwaerts(R,b)
%x=vorwaerts(L,b)

%A = [1 0;1 1;1 4];
%b = [4 -2 6]';
%x = normalengleichung(A,b)
%x=qr_ausgleich(A,b)

%--------------------------------------------------------------------------
%--------------------- a) rueckwaertseinsetzen  ---------------------------

function x = rueckwaerts(R,b)
    n = length(R);
    x = zeros(n, 1);
    for i = n: -1 : 1
        if i < n
           b(i) = b(i) - sum( R(i, i+1: n) * x(i+1 : n, 1) ); 
        end
        x(i) = b(i) / R(i, i);
    end    
end

%--------------------------------------------------------------------------
%--------------------- b) vorwaertseinsetzen  -----------------------------

function x = vorwaerts(L,b)
    n = length(L);
    x = zeros(n, 1);
    for i = 1: 1 : n
        if i > 1
           b(i) = b(i) - sum( L(i, 1: i) * x(1 : i, 1) ); 
        end
        x(i) = b(i) / L(i, i);
    end    
end

%--------------------------------------------------------------------------
%--------------------- c) normalengleichung  ------------------------------

function x = normalengleichung(A,b)
    % A*x=b
    % A'*A*x=A'*b
    % B*x=c
    B = A'*A;
    c = A'*b;
    % R*R'*x=c
    R=chol(B);

    %R'*y=c
    y=vorwaerts(R',c);
    %R*x=y
    x=rueckwaerts(R,y);
end

%--------------------------------------------------------------------------
%--------------------- d) qr-ausgleichsproblem  ---------------------------

function x = qr_ausgleich(A,b)
    % A'*A*x = A'*b
    % AN*x=bN
    n=rank(A);
    [Q, R]=qr(A);
    c=Q'*b;

    % R_t*x=c_t
    R_t=R(1:n,:);
    c_t=c(1:n);
    x=rueckwaerts(R_t,c_t);
end

%--------------------------------------------------------------------------
%--------------------- e) stabilitaet  ------------------------------------

function ausgleich_stabilitaet()
    err1=0;
    err2=0;
    ii=(-7:1:7);
    for i=ii(1:end)
        a=10^i;
        A=[1 1;1 1+a;1 1+a;];
        b=[2 2+a 2+a]';
        
        % analytisch:
        % I.    x1+x2 = 2
        % II.   x1+(1+a)*x2 = 2+a
        % aus I: x1 = 2-x2
        % in II: x2*a = a
        % --> x2 = 1
        % --> x1 = 1
        x = [1; 1]
        
        x_ng = normalengleichung(A,b)
        x_qr = qr_ausgleich(A,b)
    
        err1=[err1,norm(x-x_ng,'inf')];
        err2=[err2,norm(x-x_qr,'inf')];
    end

    plot(ii,err1(2:end),'-b',ii,err2(2:end),'--r')
    ylim([-0.01 0.08])
    title('Abweichung beim Lösen von Ax=b')
    legend('Normalengleichung','QR-Zerlegung')
    xlabel('Potenz i')
    ylabel('Fehler in \infty-Norm ')
end

