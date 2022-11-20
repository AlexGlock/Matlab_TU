%% Aufgabe H4 aus Hausübung HU4
% Code von Alexander Glock
clearvars

% Vorgabe der AWP Dimension A_awp dim := (n_awp+1)x(n_awp+1)
n_awp=50;

% Anfangswertproblem Paramter: 
% --> make_task1(n) = System aus Aufgabenteil (I)
% --> make_task2(n) = System aus Aufgabenteil (II)
[A_awp,y0,t0,T]=make_task2(n_awp); %make_task2(n)
f = @(t,y) (-1)*A_awp*y;

% Verfahrensauswahl mit make_RK(sel):
% sel: 1=Heun (a)   2=Gauss-2 (b)   3=imp. Trapez (c) 
[A,b,g]=make_RK(3);

% experimentelle Ordnungsermittlung:
test = convPlot(y0,f,t0,T,A,b,g,A_awp);


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
function [A_awp,y0,t0,T] = make_task2(n)
    % A Diagonale erstellen - Dim n+1
    d0 = ones(1, n+1);
    % A Nebendiagonalen erstellen - Dim n+1
    d1 = -1 + zeros(1, n);
    % Systemmatrix A generieren - Dim (n+1)x(n+1)
    A_awp = n*(diag(d0)+diag(d1,-1));
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
function [A,b,g] = make_RK(sel)
    switch sel
        case 3 % imp. Trapez
            A =[0 0; 1/2 2/2];  
            b =[1/2 1/2]';
            g =[0 1]';
        case 2 % Gauss-2
            A =[1/4 1/4-sqrt(3)/6; 1/4+sqrt(3)/6 1/4];  
            b =[1/2 1/2]';
            g =[1/2-sqrt(3)/6 1/2+sqrt(3)/6]';
        otherwise % Heun
            A =[0 0; 1 0];  
            b =[1/2 1/2]';
            g =[0 1]';
    end
end
function dY = myODE(Y,A_awp)
    dY= (-1)*A_awp*Y;
end
function test = convPlot(y0,f,t0,T,A,b,g,A_awp)
        h1=@(h) h;
        h2=@(h) h.^2;
      
        hmax = T/10;                    % maximale Zeitschrittweite
        hmin = T/1000;                  % minimale Zeitschrittweite
        hh = flip(hmin:15*hmin:hmax);   % Zeitschwrittgrößen von hmax bis hmin
        cc = 0;                         % norm fehlervektor initialisieren

    for h=hh(2:end)
        n=round(T/h);
        % Approximation mit RKV
        [tt, yy] = myRK(y0,f,t0,T,n,A,b,g);
        % "exakte" Lösung des DGL Systems mit ODE berechnen;
        sol=ode45(@(t,Y) myODE(Y,A_awp),tt,y0);%[~, yyex]
        [yyex,~]=deval(sol,tt,1);

        % norm. Fehler des RKV
        c = norm((yyex-yy(1,:)),'Inf');
        cc = [cc,c];
    end
    
    % Konvergenz plot
    loglog(hh,cc,'-r',hh,h1(hh),'--b',hh,h2(hh),'--g')
    title('Numerical error vs Stepwidth in double log scale')
    xlabel('size of h')
    ylabel('Numerical error || y-y_h ||')
    legend('Num error','lin. konv','quad. konv.')
    grid on

    test='fertig';
end