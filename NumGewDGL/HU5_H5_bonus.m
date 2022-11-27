%% Aufgabe H2 aus Hausübung 5
% Code von Alexander Glock
clearvars

% Paramter
t0=0;
T=4*pi;                 % 2 Perioden = 4 pi
y0=[1;0];               % 1 cos Lsg
n_max=1000;             % ganzzahl größer 100

% Verfahrensauswahl order: 3 oder 4 mit impl=true/false
order = 4;
impl= false;

% AWP und Verfahren initialisieren
[A,beta,gamma,msg]=select_RK(order,impl)
[f, S] = create_DES();

% erzeuge logspace vektor mit 200 GANZEN Zahlen zwischen 100 und n_max 
n_log=round(logspace(log10(100),log10(n_max),50));

% fehlerbestimmung + plot
hh=0;
error=0;
n_end=n_log(end);
for n=n_log(1:end)
    if n==n_end
        tic
        [tt, yy]=myRK(y0,f,t0,T,n,A,beta,gamma);
        time = toc
    else
        [tt, yy]=myRK(y0,f,t0,T,n,A,beta,gamma);
    end
    error=[error,norm(yy(:,1)-yy(:,end))];
    hh=[hh,(T-t0)/n];
end
min_error = error(end)
loglog(hh(2:end),error(2:end))
title('Fehler in Abhänigkeit zur Gitterweite')
xlabel('Gitterweite h')
ylabel('Verfahrensfehler nach 2 Perioden')
grid on


% Ergebnis:
% Bei allen gegebenen Verfahren wird der kleinste Fehler bei der kleinsten
% Gitterweite h erreicht. Bei der Verwendung von n_max=1000 ergeben sich
% diese minimale Gitterweite zu ~ 10^(-2) und die dortigen Fehler:
%   2.6*10^(-9)   RK4 - explizit, ordnung 4             0.0153 sec
%   0.0387        Simpson - explizit, ordnung 3         0.0117 sec
%   4.34*10^(-10) Gauss2 - implizit, ordnung 4          0.7091 sec
%   3.4*10^(-7)   Radau-IA - implizit, ordnung 3        0.7056 sec
%
% Die Zeiten wurden immer beim Durchlauf mit der kleinsten Gitterweite
% und ohne Änderungen am myRK-solver gemessen (TolX,TolFun=1e-14).
% Anschließend wurden die Paramter auf Tolx, TolFun = 1e-10 geändert.
% Dadurch konnten für die impliziten verfahren folgende Ergebnisse erreicht
% werden:
%   4.6*10^(-10)  Gauss2 - implizit, ordnung 4          0.6698 sec
%   3.5*10^(-7)   Radau-IA - implizit, ordnung 3        0.6222 sec
% Bei einer Verschärfung der Toleranzen auf Tolx, TolFun 1e-18 hingegen:
%   4.35*10^(-10)  Gauss2 - implizit, ordnung 4         0.7641 sec
%   3.46*10^(-7)   Radau-IA - implizit, ordnung 3       0.7262 sec
%
% Fsolve ist ebenfalls ein iteratives Verfahren. Dieses Verfahren liefert
% kein absolutes Ergebnis sondern eine qualitative Näherung als Lösung für
% die nichtlinearen Gleichungssysteme bei impliziten RKV. Dieser iterative
% Gleichungslöser besitzt Abbruchbedingungen (Tolx, TolFun) mit denen das
% Näherungsverfahren ab einem bestimmten Punkt als ausreichend genau
% hingenommen wird. Bei der Verschärfung der Abbruchbedingungen steigt die
% Qualität der fsolve approximation aber eben auch die Zahl der Iterationen
% und damit die Dauer des RK-Verfahrens


%--------------------------------------------------------------------------
%------------------ Differentialgleichungssystem  y'' = -y  ---------------

function [f, S] = create_DES()
    
    syms y(t)
    de = diff(y,2) == -y;
    [V, S]=odeToVectorField(de);
    % Struktur von sym Y = [y; diff(y)]
    f = matlabFunction(V, 'vars', {'t','Y'});

end

%--------------------------------------------------------------------------
%-------------------------  Verfahren definieren --------------------------

function [A,beta,gamma,msg]=select_RK(order,impl_true)
    if order==4 && impl_true
        %Gauss Two-Step - implizites Verf.
        A = [[1/4,1/4 - sqrt(3)/6];[1/4 + sqrt(3)/6,1/4]];
        beta = [1/2,1/2]';
        gamma = [1/2-sqrt(3)/6,1/2 + sqrt(3)/6]';
        msg = "selected gaussII order 4";
    
    elseif order==3 && impl_true
        %Radau-IA
        A = [[1/4,-1/4];[1/4, 5/12]];
        beta = [1/4,3/4]';
        gamma = [0, 2/3]';
        msg = "selected Radau-IA order 3";
    
    elseif order==4 && not(impl_true)
        %klassisches RK - explizites Verf.
        A = [0 0 0 0; 1/2 0 0 0; 0 1/2 0 0; 0 0 1 0];
        beta = [1/6,1/3,1/3,1/6]';
        gamma = [0, 1/2, 1/2, 1]';
        msg = "selected Runge-Kutta4 order 4";
    
    elseif order==3 && not(impl_true)
        %Simpsonregel- explizites Verf.
        A = [0 0 0 ; 1/2 0 0 ; -1 2 0];
        beta = [1/4,0,3/4]';
        gamma = [0, 1/2, 1]';
        msg = "selected Simpson order 3";
    else
        error("unsupported order input (3 or 4)")
    end
end

%--------------------------------------------------------------------------
%-------------------------  myRK solver -----------------------------------

function [t,y] = myRK(y0,f,t0,T,n,A,b,g)
%%Input: 
%y0 initial value as colum vector dx1 with d the Dimension (vertical)
%f right hand side. goes from (1) X (d) -> (d) (vertical)
%t0 the initial time 
%T the desired stopping time
%n the number of steps to be done
%A the RK matrix
%b,g column vectors
%
%% Output:
%Array t with the list of times (horizontal)
%Array y with the list of approximated values of y(t_i) (d)x(n+1)-Shaped



%Preallocate the arrays for t,y
t = zeros(1,n+1);
t(1) = t0;

y = zeros(length(y0),n+1);
y(:,1) = y0;


%Calculate step width h 
h = (T-t0)/n;

   %Einstellungen für fsolve (implizit) definieren:
        %Für eine "gute" Implementation erlaubt man die Wahl der toleranzen
        %als Variable statt es im Code festzulegen. 
 options = optimoptions('fsolve','Display','none','TolX',1e-14,'TolFun',1e-14);

%Main loop for repeatedly evaluating a single step of the RKV
%Repeat a single step of the RKV as many times as required



%%
for curStep = 2:n+1
    
 y(:,curStep) = generalRKStep(t(curStep-1),y(:,curStep-1),A,g,b,h,f,options);

 t(curStep) = t(curStep-1)+h;
end

end
%% Auxillary Function below for single-step RKV
function yi = generalRKStep(ti,yi,A,gamma,beta,h,f,options)
    %matlabfunktion für einen einzelnen Schritt eines beliebigen Runge-Kutta
    %verfahren 
    %Input: Butchertableau, die Schrittweite h, die
    %rechten Seite der DGL sowie das aktuelle Wertepaar (ti,yi)
    
    %Rückgabe: y_(i+1) die neue Approximation bei ti+h

    %Konventionen für input:
    %f ist Spaltenvektor [ ; ; ...]
    %beta ist Spaltenvektor
    %gamma ist Spaltenvektor
    %y ist (wie f) Spaltenvektor

    beta = transpose(beta);
    %Intern ein mal transponieren weil ich Zeilenvektoren bevorzuge
    
   
    szRHS = size(yi,1); %Dimension der rechten Seite der DGL
    
    numStages = size(A,1); %Stufenzahl des Verfahrens
                           %Bemerkung: Hier wird nicht noch geprüft dass
                           %die eingaben des Nutzers "Sinnvoll" sind, also
                           %z.B. gamma und A kompatible Dimensionen haben.
                           %Das wäre in einer "richtigen" Implementation
                           %sicherlich hilfreich
    
    yi = yi';% Ich bevorzuge wieder Zeilenvektoren
   
    if istril(A) && (trace(abs(A)) == 0) %istril: Prüft auf (schwache) untere Dreiecksform
                                         % Summe über betrag der Spur um
                                         % Diagonale zu prüfen
         mode = 0; %explicit      
    else 
        if istril(A)
            mode = 2;%diagonal implicit
        else
            mode = 1; % implicit     
        end
    end
                     
    if mode == 0 %Explizit lösen wir in g-form
         fbar = @(t,y) f(t',y')';%Ist nur ein kosmetischer Trick weil ich lieber
                                 %Mit y, f als Zeilenvektoren
                                 %arbeite
                                 
        g = zeros(numStages,szRHS);%g's als Matrix der richtigen Dimension anlegen
            
                                %g nach Formel explizit berechnen
        for i = 1:numStages
            g(i,:) = yi;
            for j = 1:(i-1)
                g(i,:) = g(i,:) + h.*A(i,j).*fbar(ti + h*gamma(j),g(j,:));
            end
        end
        for i = 1:numStages
           yi = yi + h*beta(i)*fbar(ti + h*gamma(i),g(i,:)); %y_(i+1) nach Formel
        end
    elseif mode == 2
        %Implizit bietet sich die k-Form des RKV an
        kvec = zeros(numStages,szRHS);%Initialisieren der k. Die k sind ebenfalls
                                      %Genauso wie y Zeilenvektoren, die
                                      %Stufen in die Spaltenrichtung
                                      %angeordnet

        fbar = @(t,y) f(t',y')';      %Wieder wie oben: fbar ist aus kosmetischen 
                                      %Gründen mit diesem doppelten
                                      %Transponiert versehen damit man y,f,k,..
                                      %alles als Zeilenvektoren schreiben kann
        %Aufgrund der diagonal-imliziten Struktur erhalen wir eine Schleife
        %über die Stufen
        for i=1:numStages
            %Erst bekannte Stufenwerte aufaddieren wegen diagonal-imp
            kSum = yi;
            for j=1:(i-1)
                kSum = kSum + h*A(i,j)*kvec(j,:);
            end
            %Definieren der Fixpunktfunktion für diese Stufe
            Psi = @(k) fbar(ti + gamma(i)*h,kSum + h*A(i,i)*k);
            %Umschreiben in Nullstellenform für fsolve
            G = @(k) k - Psi(k); 
            %Lösen der nichtlinearen Gleichung:
            kvec(i,:) = fsolve(G,Psi(kvec(i,:)),options); 
            %Diese Rechnung kann man noch beschleunigen z.B. durch das
            %explizite Ausrechnen der Ableitungsmatrix von G und dem Übergeben
            %als Argument an fsolve. Hier wird stattdessen die Ableitung numerisch für die
            %internen Algorithmen nur geschätzt   
        end
        for i = 1:numStages
           yi = yi + h*beta(i)*kvec(i,:); %y_(i+1) nach Formel (in k-form)
        end
    else
        %Allgemeiner impliziter Löser
        %Implizit bietet sich die k-Form des RKV an
        kvec = zeros(numStages,szRHS);%Initialisieren der k. Die k sind ebenfalls
                                      %Genauso wie y Zeilenvektoren, die
                                      %Stufen in die Spaltenrichtung
                                      %angeordnet
       
         fbar = @(t,y) f(t',y')'; %Wieder wie oben: fbar ist aus kosmetischen 
                                  %Gründen mit diesem doppelten
                                  %Transponiert versehen damit man y,f,k,..
                                  %alles als Zeilenvektoren schreiben kann
         
        %Definition der Fixpunkt-funktion k = Psi(k)
        %Geschieht weiter unten da nested function definitions nicht
        %innerhalb von if/while/for etc auftauchen dürfen
        
        %Umschreiben in Nullstellenproblem für Newton/Fsolve
        G = @(kvec) kvec - Psi_imp(kvec);
        
        %Lösen der nichtlinearen Gleichung:
        kvec_new = fsolve(G,Psi_imp(kvec),options); 
        %Diese Rechnung kann man noch beschleunigen z.B. durch das
        %explizite Ausrechnen der Ableitungsmatrix von G und dem Übergeben
        %als Argument an fsolve. Hier wird stattdessen die Ableitung numerisch für die
        %internen Algorithmen nur geschätzt.
        
        %Klassischer Update-Schritt in k-Form:
        yi = yi + (h.*beta*kvec_new);
    end
    yi = yi'; %Zurücktransponieren der yi for Ende der Funktion damit das Ausgabe
              %Argument wieder Spaltenvektor genauso wie die Eingabe ist.
    function psiK = Psi_imp(kVec)
        psiK = zeros(size(kVec));
        %Wir nennen die variable i_loc um auf den lokalen Scope der
        %Variable in der "nested function" hinzuweisen.
        for i_loc = 1:numStages
            kSumLoc = yi;
            for j_loc = 1:numStages
                kSumLoc = kSumLoc + h*A(i_loc,j_loc)*kVec(j_loc,:);
            end
            psiK(i_loc,:) = fbar(ti + gamma(i_loc)*h,kSumLoc);
        end
    end
end