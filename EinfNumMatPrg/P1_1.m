%% Aufgabe P1.1 - Gruppe 7 - Alexander Glock, Jannis Röder
% Die exakte Lösung der Reihe wird symbolisch berechnet (Z.49)
% Zu diesem Zweck wird die Matlab symbolic Toolbox benutzt.
clearvars

% Parametervorgabe
N = 2000000;
plot_partial_sums(N)

% d) Versuchserkenntnisse:
%   Die Summationsrichtung spielt zunächst keine Rolle (kommutativ).
%   Ab ~n=210000 konvergiert der Vorwärtssummenfehler nicht weiter.  
%   - Die Vorwärtssumme beginnt mit der größten Zahl(=1). Inkrement wird zunehmend kleiner
%   => Maschinengenauigkeit wird erreicht, es wird nichts mehr addiert
%   - Die Rückwärtssumme beginnt mit der kleinsten Zahl. Inkrement wächst nurnoch
%   => günstiger für begrenzten Stellengenauigkeit von double

% a) Vorwärtssumme - Fkt.
function S = forward_sum(N)

    kk=(1:1:N);
    S=0;
    for k = kk(1:end)
        S=S+(1/k^3);
    end
end

% b) Rückwärtssumme - Fkt.
function S = backward_sum(N)

    kk=flip(1:1:N);
    S=0;
    for k = kk(1:end)
        S=S+(1/k^3);
    end
end

% c) Abweichungsplot - Fkt.
function plot_partial_sums(N)
    % check ob N definiert, sonst 2000000
    if not(exist('N','var'))N=2000000; end
    % erzeuge logspace vektor mit 200 GANZEN Zahlen zwischen 1 und N 
    n_log=round(logspace(log10(1),log10(N),200));

    % exakte Lösung der Reihenkonvergenz ist die zeta 3 funktion
    %S_inf = zeta(3);
    syms k
    S_inf=symsum(1/k^3,k,1,Inf);
    % Initialisierung der Ergebnisvektoren auf 0
    En_up=0;
    En_down=0;

    for n = n_log
        % berechnen der Partialsummen
        S_up = forward_sum(n);
        S_down = backward_sum(n);
    
        % Berechnen der Abweichungen & anhängen an Ergebnisvektor
        En_up = [En_up,abs(S_inf-S_up)];
        En_down = [En_down,abs(S_inf-S_down)];
    end

    % loglog Plot, ab zweitem Wert weil Wert1 = 0
    loglog(n_log,En_up(2:end),n_log,En_down(2:end))
    title('Abweichung der Partialsummen')
    legend('Vorwärtssummenfehler E_n \uparrow','Rückwärtssummenfehler E_n \downarrow')
    xlabel('N')
    ylabel('|S_\infty-S_n|')
    grid on
end



