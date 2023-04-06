%% Aufgabe P5.1 - Gruppe 7 - Alexander Glock, Jannis Röder
% Die Aufgabe 1 benötigt das Paket chebfun: https://www.chebfun.org/download/
clearvars

% Beispiel funktion chebfun, interpo. ordnung n:
% f = chebfun('sqrt(abs(x-0.5))', [-1 1],'splitting','on');
% N = 10;

plot_bestapproximationen()

% Glattheit der Funktion führt zu geringerem Interpolationsfehler, bei
% vergleichbarem Polynomgrad => Bessere Abbildung durch
% Interpolationspolynom möglich

%--------------------------------------------------------------------------
%---------------------  b) L2 best-appr.  ---------------------------------

function pn = L2_best_approximation(f,n)

    % min | f - u |_L2 := tscheby interpol.
    % N = n+1 = Datapoints = polyn of deg n
    % intervall [-1, 1] = default
    pn = chebfun(f, n+1); 
    
    % plot(f,'b',pn,'r'), grid on

end

%--------------------------------------------------------------------------
%---------------------  c) Linf best-appr.  -------------------------------

function qn = Linf_best_approximation(f,n)

    % one line
    qn = minimax(f,n);

    % clf, plot(f,'b',qn,'r'), grid on

end

%--------------------------------------------------------------------------
%---------------------  d) bestapprox  ------------------------------------

function bestapproximationen(f,N)
    % loop over all n
    nn = (1:1:N);
    err_qn_L2 = zeros(N, 1);
    err_pn_L2 = zeros(N, 1);
    err_qn_inf = zeros(N, 1);
    err_pn_inf = zeros(N, 1);
    for n = nn(1:N)
        
        % calc interpolations
        pn = L2_best_approximation(f, n);
        qn = Linf_best_approximation(f, n);
    
        % calc errors
        err_pn_L2(n) = norm(pn-f);
        err_qn_L2(n) = norm(qn-f);
        err_pn_inf(n) = norm(pn-f, 'inf');
        err_qn_inf(n) = norm(qn-f, 'inf');
    
        if n == N
            % plot max interpolations and errors 1. - 4.
            figure
            subplot(2,2,1)
            plot(f,'b',pn,'r'), grid on
            legend('f(x)','p_n(x)')
            title('L2 Bestapproximation für max n')
            subplot(2,2,2)
            plot(f,'b',qn,'r'), grid on
            legend('f(x)','q_n(x)')
            title('Linf Bestapproximation für max n')
            subplot(2,2,3)
            plot(abs(f-pn),'g'), grid on
            subplot(2,2,4)
            plot(abs(f-qn),'g'), grid on
            set(gcf,'Position',[100 100 1000 500])
 
        end
    end
    
    % plot errors over rising n
    figure
    subplot(1,2,1)
    plot(nn, err_pn_L2,'b',nn, err_qn_L2,'r'), grid on
    legend('p_n(x)','q_n(x)')
    title('Fehler in L2 norm über n')
    subplot(1,2,2)
    plot(nn, err_pn_inf,'b',nn, err_qn_inf,'r'), grid on
    legend('p_n(x)','q_n(x)')
    title('Fehler in Linf norm über n')
    set(gcf,'Position',[100 100 1000 500])

end

%--------------------------------------------------------------------------
%---------------------  e) bestapprox multi  ------------------------------

function plot_bestapproximationen()

    N = 20;
    
    % tanh()
    f = chebfun('tanh(10.0*x-5.0)', [-1 1],'splitting','on');
    bestapproximationen(f, N);
    
    % 1/2*|x|
    f = chebfun('(1.0/2.0)*abs(x)', [-1 1],'splitting','on');
    bestapproximationen(f, N);
    
    % sin(2*pi*x)
    f = chebfun('sin(2.0*pi*x)', [-1 1],'splitting','on');
    bestapproximationen(f, N);

end
