%% Aufgabe P4.1 - Gruppe 7 - Alexander Glock, Jannis Röder
clearvars


interpolation_time()


% Testfunktionswerte von tangens
%xi = [-1.5 , -0.75, 0, 0.75, 1.5]';
%fi = [-14.101420, -0.931596, 0, 0.931596, 14.101420]';

%--------------------------------------------------------------------------
%--------------------- a) Lagrange-Basispolynome  -------------------------

function y = Lagrange_basis(xi,i,x)
    % init vars
    x_ind_list = (1:1:length(xi));
    y=1;

    % loop over index
    for x_ind=x_ind_list
        if xi(x_ind)~=xi(i)
            y=y*(x-xi(x_ind))/(xi(i)-xi(x_ind));
        end
    end
end

%--------------------------------------------------------------------------
%--------------------- b) Lagrange-Interpolation an x  --------------------

function y = Lagrange(xi,fi,x)
    nn = (1:1:length(xi));
    y=0;
    for n=nn
        L = fi(n)*Lagrange_basis(xi, n, x);
        y = y+L;
    end
end

%--------------------------------------------------------------------------
%--------------------- c) Poly. Auswertung HORNER  ------------------------

function y = Horner_eval(c,x)
    nn=flip(2:1:length(c));
    y=0;
    for n=nn
        y= (y + c(n))*x;
    end
    y = y + c(1);
end

%--------------------------------------------------------------------------
%--------------------- d) dividierte Differenzen  -------------------------

function c = dividierte_differenzen(xi,fi)
n = length(xi);
nn=(1:1:n);
F = zeros(n,n);
F(:,1)=fi;

% Compute the rest of the table % k=2 i=1
for k=nn(2:end)
    for i=nn(1:end-k+1)
      F(i,k) = F(i + 1,k - 1) - F(i,k - 1) / (xi(i + 1) - xi(i));
    end
end
c = F(1,:)';
end

%--------------------------------------------------------------------------
%--------------------- e) Newton mit Horner Ausw.  ------------------------

function y = Newton_horner(c,xi,x)
% x differenzen bestimmen
x_d = x - xi(1:end-1);
% horner auswertung
nn=flip(2:1:length(c));
y=0;
for n=nn
    y= (y + c(n))*x_d(n-1); 
end
y = y + c(1);
end

%--------------------------------------------------------------------------
%--------------------- f) Zeitmessungen  ----------------------------------

function interpolation_time()
% Ergbenisvektoren initialisieren:
Lagrange_ges = zeros(1, 90);
LGS_prep = Lagrange_ges;
horner_eval = Lagrange_ges;
newton_prep = Lagrange_ges;
newton_ges = Lagrange_ges;

% 1000 äquidistante Auswertungspunkte im Intervall:
nn=(0:1/1000:1);
n_pl=(10:1:100);
for n=n_pl

% i) Stützstellen u. Funktionswerte auf [0,1]:
f =@(x) sin(x);
xi=(0:1/n:1);
fi = f(xi)';
% ii) Zeit für 1000 Lagrange Auswertungen:
tic
for ni=nn
y1 = Lagrange(xi, fi, ni);
end
Lagrange_ges(n-9) = toc;
% iii) Zeit für LGS aufstellen und lösen:
tic
v_mat = vander(xi);
X = v_mat \ fi;
LGS_prep(n-9) = toc;
% iv) Zeit für 1000 Polynom Auswertung mit Horner:
tic
for ni=nn
y2 = Horner_eval(X, ni);
end
horner_eval(n-9) = toc;
% v) Zeit für dividierte Differenzen:
tic
c = dividierte_differenzen(xi, fi);
newton_prep(n-9) = toc;
% vi) Zeit für 1000 Newton(mit Horner) Auswertung:
tic
for ni=nn
y3 = Newton_horner(c,xi,ni);
end
newton_ges(n-9) = toc;

end

% Auswertung erfolgt bei allen "horner artig"
newton_eval = horner_eval;
Lagrange_eval = horner_eval;
LGS_eval = horner_eval;
% Gesamtdauer wird bei Newton und Lagrange direkt gemessen
LGS_ges = LGS_prep+LGS_eval;
% Vorbereitung wird bei LGS und Newton (divDiff) direkt gemessen:
Lagrange_prep = Lagrange_ges-horner_eval;

figure(1)
semilogx(n_pl, Lagrange_prep, n_pl, LGS_prep, n_pl, newton_prep)
title('Dauer der Vorbereitungsrechnungen')
legend('Lagrange Interpolation','LGS - Ansatz', 'Newton Interpolation')
xlabel('Stützstellenanzahl n')
ylabel('Zeit in Sekunden')
grid

figure(2)
semilogx(n_pl, horner_eval)
title('Dauer der Funktionsauswertungen mit Horner')
legend('für alle Ansätze gleich')
xlabel('Stützstellenanzahl n')
ylabel('Zeit in Sekunden')
grid

figure(3)
semilogx(n_pl, Lagrange_ges, n_pl, LGS_ges, n_pl, newton_ges)
title('Gesamtdauer der Ausgleichsrechnungen')
legend('Lagrange Interpolation','LGS - Ansatz','Newton Interpolation')
xlabel('Stützstellenanzahl n')
ylabel('Zeit in Sekunden')
grid

end

