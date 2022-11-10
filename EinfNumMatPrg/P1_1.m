%% Aufgabe P1.1 - Gruppe 7 - Alexander Glock, Jannis Röder
clearvars
% Parameter
N=1000;

%t=logspace(0,3,3);
%t
%plot_partial_sums(N)



% vorwärtssumme
function S = forward_sum(N)

    kk=(1:1:N);
    S=0;
    for k = kk(1:end)
        S=S+(1/k^3);
    end
end

% rückwärtssumme
function S = backward_sum(N)
    kk=flip(1:1:N);
    S=0;
    for k = kk(1:end)
        S=S+(1/k^3);
    end
end

% plot fkt.
function plot_partial_sums(N)
    % check ob N definiert, sonst 2000000
    if not(exist('N','var'))N=2000000; end

    nn = (1:1:N);
    % analytische Lösung ist die zeta 3 funktion
    S_inf = zeta(3);
    En_up=0;
    En_down=0;

    for n = nn(1:end)
        % berechnen der Partialsummen
        S_up = forward_sum(n);
        S_down = backward_sum(n);
    
        % Berechnen der Fehler & anhängen an vektor
        En_up = [En_up,abs(S_inf-S_up)];
        En_down = [En_down,abs(S_inf-S_down)];
    end

    % Konvergenz plot

    %En_up_log=En_up(logspace(1,N,200))
    %En_down_log=En_down(logspace(1,N,200))
    %nn_log=nn(logspace(1,N,200))
    loglog(nn,En_up(2:end),nn,En_down(2:end))
    title('Konvergenzverhalten der Partialsummen')
    legend('Vorwärtssumme','Rückwärtssumme')
    xlabel('N')
    ylabel('|S-S_n|')
    grid on
end



