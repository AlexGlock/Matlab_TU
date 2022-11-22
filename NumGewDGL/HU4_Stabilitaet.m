%% Aufgabe H4 aus Hausübung HU4
% Code von Alexander Glock
clearvars

% Vorgabe der Matrix Dimension des AWP := (n_awp+1)x(n_awp+1)
n_awp=30;
% Anfangswertproblem generieren: 
% --> make_task1(n) = System aus Aufgabenteil (I)
% --> make_task2(n) = System aus Aufgabenteil (II)
[A_awp,y0,t0,T]=make_task2(n_awp); %make_task2(n)

% ODE der Form:
f = @(t,y) (-1)*A_awp*y;

% Verfahrensauswahl mit make_RK(sel):
% sel: 1=Heun a)   2=Gauss-2 b)   3=imp. Trapez c) 
[A,b,g]=make_RK(2);

% analytische LSG vorbereiten:
[EV, D] = eig(-A_awp);  % Eigenvektoren              
E=diag(D);              % Eigenwerte
C= linsolve(EV,y0);     % Konstanten mittels y(0) Bed.

% experimentelle Ordnungsermittlung mit Plot:
test = convPlot(y0,f,t0,T,A,b,g,EV,E,C);

% schriftliche Ergebnisse ganz unten 
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
function Y = Yex(C,E,EV,t)
    Y=zeros(1,length(E))';
    for k=(1:1:length(E))
        Y_k= C(k)*EV(k)*exp(E(k)*t);
        Y=Y+Y_k;
    end
end
function test = convPlot(y0,f,t0,T,A,b,g,EV,E,C)
        h0=@(h) 1;
        h1=@(h) h;
      
        hmax = T/100;                    % maximale Zeitschrittweite
        hmin = T/10000;                  % minimale Zeitschrittweite
        hh = flip(hmin:10*hmin:hmax);    % Zeitschwrittgrößen von hmax bis hmin
        cc = 0;                          % norm fehlervektor initialisieren

    for h=hh(2:end)
        n=round(T/h);
        % Approximation mit RKV
        [tt, yy] = myRK(y0,f,t0,T,n,A,b,g);
        % "exakte" Lösung des DGL Systems zu tt berechnen:
        yyex=Yex(C,E,EV,tt);

        % norm. Fehler des RKV
        c = norm((yyex(1,:)-yy(1,:)),'Inf');
        cc = [cc,c];
    end
    
    % Konvergenz plot
    loglog(hh,cc,'-r',hh,h0(hh),'--b',hh,h1(hh),'--g')
    title('Numerical error vs Stepwidth in double log scale')
    xlabel('size of h')
    ylabel('Numerical error || y-y_h ||')
    legend('Num error','konst 1','lin. konv')
    grid on

    test='fertig';
end

%% Ergebnis:
% Mit den gegebenen Anfangswertproblemen besitzen alle drei Verfahren die
% Konvergenzordnung 0, da sie nicht konvergieren und damit einen konstanten
% Fehler zur analytischen Lösung vorweisen der niemals verschwindet.

% Die impliziten Verfahren (b & c) sind für alle Schrittweiten h stabil.
% Das äußert sich im Konvergenzplot durch einen gleichmäßigen Verlauf des numerischen
% Fehlers mit endlicher Steigung. Beim expliziten Verfahren (a) entwickelt
% sich ab einem h von ungefähr 10^(-2) eine Instabilität. Zu sehen ist das am
% numerische Fehler der an diesem Punkt im Konvergenzplot rasant ansteigt und sehr groß wird.



