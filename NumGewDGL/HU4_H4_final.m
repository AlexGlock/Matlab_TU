%% Hausübung HU4 - Aufgabe H4
% Code von Alexander Glock

% Das explizite Verfahren von Heun konvertgiert als einziges Verfahren mit
% der Ordnung 2 (h^2). Die anderen beiden konvergieren lediglich linear,
% also mir der Ordnung 1.
% Beim Verfahren von Heun führt ein Erhöhen der Systemgröße n bei sonst
% gleichen Parametern zum früheren erreichen der kritischen Schrittweitengröße 
% (= Instabilität des num. Verfahrens). 
% Die anderen beiden Verfahren sind implizit und deshalb nicht von diesem
% Effekt betroffen.
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

% Wahl numerisches Verfahren

% Verfahrensauswahl mit make_RK(sel):
% sel: 1=Heun a)   2=Gauss-2 b)   3=imp. Trapez c) 
[A,b,g]=make_RK(1);


%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% Dimesion A-Matrix zur DGL
n = 20;

% Anfangswertproblem auswählen: 
% --> make_task1(n) = System aus Aufgabenteil (I)
% --> make_task2(n) = System aus Aufgabenteil (II)
[Amat,y0,t0,T]=make_task2(n); %make_task2(n)

% Parameter für Zeitintervall und Zeitschritte
n_t   = 800;

% Erstellen der A-Matrix zur DGL
Amat = (-1)*Amat;
% Funktion zur DGL
f = @(t,y) (Amat*y);

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

% Kalkulation der analytischen homogenen Lösung
[yhom,~] = yhomogen(Amat,y0,n,t0,T,n_t);

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

% Experimentelle Konvergenzuntersuchung 

narray = (100:10:200);
h      = T ./ narray;
y_sup  = []; 

for n_t = narray
    
    % Numerische Lösung
    [t,y] = myRK(y0,f,t0,T,n_t,A,b,g);
    
    % Exakte Lösung zum Vergleich
    [yhom,~] = yhomogen(Amat,y0,n,t0,T,n_t);
    
    %y_sup = [y_sup,max(abs(yhom(:,row)-y(:,row)))]; %#ok<AGROW>
    y_sup = [y_sup,norm(yhom-y)]; %#ok<AGROW>
end

% Konvergenz plot
h1=@(h) h;
h2=@(h) h.^2;
loglog(h,y_sup,'-r',h,h1(h),'--b',h,h2(h),'--g')
title('Numerical error vs Stepwidth in log scale')
xlabel('size of h')
ylabel('Numerical error || y-y_h ||')
legend('Num error','lin. konv (h)','quad. konv (h^2)')
grid on


%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

% Funktion zur Berechnung der analytischen, homogenen Lösung

%--------------------------------------------------------------------------

function [yhom,t] = yhomogen(Amat,y0,n,t0,T,n_t)

    % Berechnung der Schrittweitte h
    h = (T-t0)/n_t;

    % Berechnung der Eigenwerte und -vektoren zu Amat
    [eigvec,eigval] = eig(Amat);
    
    % Erstelle yhom leer und t 
    t       = t0:h:T;
    yhom    = zeros(n_t,n+1);
    
    % Führe Berechnung fort, solange geometrische Vielfalt 
    % der Eigenvektoren passt
    if det(eigvec) ~= 0
        
        % Bestimmung der Konstanten zur homogenen Lösung
        cons = ones(1,n+1);
        F = @(cons) (cons .*eigvec * ones(n+1,1) - y0);
        cons = fsolve(F,cons,optimoptions('fsolve','Display','none'));
    
        % Multiplikation der Konstanten mit Eigenvektoren
        vec_hom = cons.*eigvec;

       % Berechne homogene Lösung
        for i=1:1:n_t+1 

            ti = t(1,i);
            solvevec = zeros(n+1,1);
            for val=1:1:n+1
                solvevec = solvevec + vec_hom(:,val) * exp(eigval(val,val)*ti);
            end
            yhom(i,:) = solvevec';   

        end
    
    end

end

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

% Funktion zum Lösen von DGL mit RKV

%--------------------------------------------------------------------------

function [t,y] = myRK(y0,f,t0,T,n,A,b,g) 

    % Berechnung von Schrittweite h und Vektor t
    h       = T/n;
    t       = t0:h:T;

    % Erstellung von y
    dimy    = max(size(y0));
    y       = zeros(n+1,dimy);
    y(1,:)  = y0';

    % Ermittle Anzahl der Stufen des Verfahrens
    stuf = max(size(g));

    % Prüfe ob explizit, implizit oder diagonalimplizit 
    % und führe dann entsprechendes RKV-Verfahren durch
    testvec = ones(stuf,1);
    diagsum = diag(A)' * testvec;
    trilsum = ((A .* triu(ones(stuf,stuf),1)) * testvec)' * testvec;

    if diagsum == 0 && trilsum == 0

        % Explizites RKV
        [y,t] = expRKV(y,f,t,h,n,A,b,g,dimy,stuf);

    elseif trilsum == 0

        % Implizites RKV (diagonalimpliziter Fall)
        [y,t] = diagimpRKV(y,f,t,h,n,A,b,g,dimy,stuf);

    else

        % Implizites RKV
        [y,t] = impRKV(y,f,t,h,n,A,b,g,dimy,stuf);

    end

end

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

% Explizites RKV

%--------------------------------------------------------------------------

function [y,t] = expRKV(y,f,t,h,n,A,b,g,dimy,stuf)

    for i = 1:1:n

        % Grundlegende Parameter für Zeitschritt
        ti  = t(i);
        yi  = y(i,:)';

        % Erstelle kvec
        kvec = zeros(dimy,stuf);
        for j=1:1:stuf
            kvec(:,j) = f(ti + g(j,1) * h, yi + h * kvec * A(j,:)');
        end

        % Berechne y an nächster Stützstelle
        y(i+1,:) = (yi + h * kvec * b)';

    end
    
end

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

% RKV für diagonalimpliziten Fall

%--------------------------------------------------------------------------

function [y,t] = diagimpRKV(y,f,t,h,n,A,b,g,dimy,stuf)

    % Anwendung implizites RKV (diagonal-impliziter Fall)
    for i = 1:1:n

        % Grundlegende Parameter für Zeitschritt
        ti  = t(i);
        yi  = y(i,:)';

        % Startwert für  kvec für aktuellen Zeitschritt
        kvec = zeros(dimy,stuf);

        % Erstelle kvec
        for j=1:1:stuf

            % Definition des zu lösenden Systems für den Iterationsschritt
            f_opt = @(v) Fdiag(v,kvec,j,ti,h,yi,A,g,f); 

            % Startwert von v für Iteration
            v_start = ones(dimy,1);

            % Lösen des Systems -> kvec wird um v erweitert
            v = fsolve(f_opt,v_start,optimoptions('fsolve','Display','none')); 
            kvec(:,j) = v;

        end

        % Berechne y an nächster Stützstelle
        y(i+1,:) = (yi + h * kvec * b)';

    end
    
    
    % Funktion zur Aufstellung des Gleichungssystems zum Verfahren
    function mindiag = Fdiag(v,kvec,j,ti,h,yi,A,g,f)

        % Erstelle kvec für Iterationsschritt anhand von Lösungsvektor
        kvec(:,j) = v;

        % Zu lösendes System
        mindiag = v - f(ti + g(j,1) * h, yi + h * kvec * A(j,:)');

    end

end

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

% Implizites RKV

%--------------------------------------------------------------------------

function [y,t] = impRKV(y,f,t,h,n,A,b,g,dimy,stuf)

    % Startwert kvec für ersten Zeitschritt
    kstart = ones(dimy,stuf);

    % Anwendung implizites RKV
    for i = 1:1:n

        % Grundlegende Parameter für Zeitschritt
        ti  = t(i);
        yi  = y(i,:)';

        % Definition des zu lösenden Systems für den Zeitschritt
        f_opt = @(kvec) F(kvec,dimy,stuf,ti,h,yi,A,g,f); 

        % Lösen des Systems -> Lösung wird als Startwert für nächsten
        % Zeitschritt verwendet
        kvec = fsolve(f_opt,kstart,optimoptions('fsolve','Display','none')); 
        kstart = kvec;

        % Berechne y an nächster Stützstelle
        y(i+1,:) = (yi + h * kvec * b)';

    end

    
    % Funktion zur Aufstellung des Gleichungssystems zum Verfahren
    function minF = F(kvec,dimy,stuf,ti,h,yi,A,g,f)

        minF = zeros(dimy,stuf);
        for j=1:1:stuf
            minF(:,j) = kvec(:,j) - f(ti + g(j,1) * h, yi + h * kvec * A(j,:)');
        end

    end

end

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

function [A_awp,y0,t0,T] = make_task1(n)
    % A Diagonale erstellen - Dim n+1
    d0 = 2 + zeros(1, n+1);
    d0(1)=0;d0(n+1)=0;
    % A Nebendiagonalen erstellen - Dim n
    d1 = -1 + zeros(1, n);
    d1(1)=0;d1(n)=0;
    % Systemmatrix A generieren - Dim (n+1)x(n+1)
    A_awp = n^2*(diag(d0)+diag(d1,-1)+diag(d1,1));
    
    % Startvektor y0 generieren - Dim (n+1)
    y0=zeros(1,n+1)';
    jj=(1:1:n+1);
    for j=jj(1:end)
        if (n+1)/4 < j && j < 3*(n+1)/4 
            y0(j)=sin((j-(n+1)/4)*2*pi/(n+1));
        end
    end
    t0=0;
    T=0.1;
end

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

function [A_awp,y0,t0,T] = make_task2(n)
    % A Diagonale erstellen - Dim n+1
    d0 = ones(1, n+1);
    % A Nebendiagonalen erstellen - Dim n+1
    d1 = -1 + zeros(1, n);
    % Systemmatrix A generieren - Dim (n+1)x(n+1)
    A_awp = (n)*(diag(d0)+diag(d1,-1));
    A_awp(1,n+1)=n*(-1);
    
    % Startvektor y0 generieren - Dim (n+1)
    y0=zeros(1,n+1)';
    jj=(1:1:n+1);
    for j=jj(1:end)
        if (n+1)/4 < j && j < 3*(n+1)/4 
            y0(j)=sin((j-(n+1)/4)*2*pi/(n+1));
        end
    end
    t0=0;
    T=2;
end

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------


function [A,b,g] = make_RK(sel)
    switch sel
        case 3 % imp. Trapez
            A =[0 0; 1/2 2/2];  
            b =[1/2; 1/2];
            g =[0; 1];
        case 2 % Gauss-2
            A =[1/4 1/4-sqrt(3)/6; 1/4+sqrt(3)/6 1/4];  
            b =[1/2; 1/2];
            g =[1/2-sqrt(3)/6; 1/2+sqrt(3)/6];
        otherwise % Heun
            A =[0 0; 1 0];  
            b =[1/2; 1/2];
            g =[0; 1];
    end
end