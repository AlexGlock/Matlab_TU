
%% Aufgabe A3
hmax = 0.1; % maaximale Zeitschrittweite
hmin = 0.00001; % minimale Zeitschrittweite

hh = flip(hmin:hmin*15:hmax); % Zeitschwrittgrößen von hmax bis 0.0001
cc = 0; % konstante initialisieren
for h = hh(2:end)

    T = 1; % Endzeit
    tt = 0:h:T; % Zeitgitter
    % Problemdefinition
    f = @(t,y) [y(2); -pi^2*y(1)]; y0=[1; 0];
    yex = @(t) [cos(pi*t); -pi*sin(pi*t)];
    yy = y0; y=y0;
    
    for t=tt(2:end)
        y = y + h*f(t,y); % explizites Eulerverfahren
        yy = [yy,y];
    end
    % exakte Lösung
    % 1
    yyex=yex(tt);
    % erste Komponente der Lösung
    %plot(tt,yyex(1,:),'b-',tt,yy(1,:),'r.')
    c = norm((yyex(1,:)-yy(1,:)),Inf)/h;
    cc = [cc,c];
end
% zweite Komponente der Lösung
% konvergenz
plot(hh(2:end),cc(2:end),'r-')
title('Konvergenzverhalten des Eulerverfahrens')
xlabel('Zeitschrittgröße h')
ylabel('Fehler-Konstante C')
set ( gca, 'xdir', 'reverse' )

%% Aufgabe A4
clearvars

%Parameter
a = 0; b = 2;
m = 2; c = 3; k = 1;

%lsg des char. polynoms:
k1=(sqrt(c^2-4*m*k)-c)/2*m
k2=(-sqrt(c^2-4*m*k)-c)/2*m 
%konstanten nach RB:
c2 = (b-a*k1)*(1/(k2-k1));
c1 = a-c2;

h = 0.01; % Schrittweite
T = 1; % Endzeit
tt = 0:h:T; % Zeitgitter
% Problemdefinition
f = @(t,y) [y(2); -(k/m)*y(1)-(c/m)*y(2)]; y0=[a; b];
yex = @(t) [c1*exp(k1*t)+c2*exp(k2*t); c1*k1*exp(k1*t)+c2*k2*exp(k2*t)];
yy = y0; y=y0;
for t=tt(2:end)
    y = y + h*f(t,y); % explizites Eulerverfahren
    yy = [yy,y];
end

yyex=yex(tt);
C = norm((yyex(1,:)-yy(1,:)),Inf)/h
% erste Komponente der Lösung
plot(tt,yyex(1,:),'b-',tt,yy(1,:),'r.')
title('Näherung von Aufgabe A1')
xlabel('Zeit t in sec')
ylabel('Lösung y(t)')
