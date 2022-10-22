%% Aufgabe H1.3 b)
% bestimmen der Maschinengenauigkeit von Matlab
clearvars

b = 10;                  % Basis des dualen Zahlensystems
p = 0;                  % start exponent 0
m = 1;                  % mantisse ist statisch 1 

z=2;                    % startwert mit exponent 0 = 1

while z>1
    p=p-1;
    z=single(1+b^(p));             % geb Zahl mit aktuallem Exp. an
end

b^(p+1)
eps("single")