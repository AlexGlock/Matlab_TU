%% Vorlesung 1 FEM Beispiel - Alexander Glock
% In diesem Skript wird ein einfacher FEM solver am Beispiel der
% eindimensionalen (Zeit) Wärmeverteilung implemetiert.
clearvars

%% Part 1 Diskretisierung der Domäne (Zeit) und Input
h = 0.1;               % Schrittweite in sekunden
t_end = 1;              % Endzeitpunkt der Berechnung
T_0 = 100;              % Temperatur Startwert
rhs=@(t) 100*exp(-t);   % rechte Seite der DGL

d_t = 0:h:t_end;
x_i = rhs(d_t);

%% Gewichtsfunktionen erster Ordnung für Index (ind) ausgewertet an t
%syms w(ind, t)
%w(ind, t) = piecewise((ind-1)*h < t) & (t < ind*h),1,(ind*h < t) & (t < (ind+1)*h),-1,0);
syms w(x)
w(x) = piecewise((x >= 0) & (x < 1),x,(x < 2) & (x >= 1),2-x,0) 

plot((-1:h:3), w((-1:h:3)))

%% Aufstellen der Systemmatrix nach Diskretisierungsvorschrift




%% Lösen des Systems