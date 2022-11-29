%% Aufgabe H2 aus Hausübung 5
% Code von Alexander Glock
clearvars

% AWP und Verfahren initialisieren
[AWP, RK4, gaussII] = createSetting();

% DGL System aus Beispiel 6.1 definieren
[M, mu] = create_DES();

% ============ ACHTUNG: ===================================================
% Das implizite Verfahren braucht sehr lang zur Berechung der Lösung mit der
% gewünschten Genauigkeit. In den Verfahrens-Structs ist jeweils die
% benötigte Schrittzahl zum geforderten Fehler als .n Attribut hinterlegt
% =========================================================================
% RK4.n, gaussII.n oder beliebiger Ganzzahlenwert
n=RK4.n  %n=gaussII.n  %n=10000  

% numerische Lösung mit Runge Kutta und dabei Zeit stoppen
%  ----> entweder gaussII oder RK4 auswählen:
tic
%[tt,yy]=myRK(AWP.Y0,M,AWP.t0,AWP.Tper,n,gaussII.A,gaussII.beta,gaussII.gamma);
[tt,yy]=myRK(AWP.Y0,M,AWP.t0,AWP.Tper,n,RK4.A,RK4.beta,RK4.gamma);
time = toc

% Plot der Lösung 
plot(yy(3,:),yy(1,:))
title('Satellitenlaufbahn innerhalb einer Periode')
xlabel('X-Koordinate')
ylabel('Y-Koordinate')

% Fehler und Schrittweite der berechneten Lösung auf Konsole ausgeben
h=(AWP.Tper-AWP.t0)/n
error=norm(AWP.Y0-yy(:,end),'inf')

% Energiebilanz der Lösung errechnen: E_end-E_anf
% yy: [y, y', x, x']
E_end = 1/2*(yy(4,end)^2+yy(2,end)^2)-1/2*(yy(3,end)^2+yy(1,end)^2)-(1-mu)/sqrt((yy(3,end)+mu)^2+(yy(1,end)))-mu/sqrt((yy(3,end)-(1-mu))^2+yy(1,end)^2);
E_anf = 1/2*(yy(4,1)^2+yy(2,1)^2)-1/2*(yy(3,1)^2+yy(1,1)^2)-(1-mu)/sqrt((yy(3,1)+mu)^2+(yy(1,1)))-mu/sqrt((yy(3,1)-(1-mu))^2+yy(1,1)^2);
E_diff = E_end-E_anf

% Ergebnisse:
% =========================================================================
% Die Lösung des DGL-Systems mit dem expliziten RK4-Verfahren benötigt
% ungefähr eine Schrittzahl (n) von 845200 und damit eine ungefähre
% Schrittgröße von h ~ 2.02*10^(-5) um die geforderte Genauigkeit zu erreichen.
% Die Lösung mittels des impliziten Gauss2-Verfahrens benötigt dagegen nur
% rund 605000 Schritte um ein vergleichbares Ergebnis zu liefern. Damit
% ergibt sich für dieses Verfahren eine konstante Schrittgröße von h ~
% 2.82x10^(-5)
% Das implizite Verfahren (Gauss2) benötigt trotz der geringeren
% Schrittzahl erheblich länger (~550 sec vs ~18 sec) um die Lösung
% zu berechnen, was auf die Lösung von impliziten Gleichungssytemen in jedem 
% Zeitschritt zurück zu führen ist. 
%
% Bei beiden Verfahren ist die Energiedeifferenz über den berechneten Zeitraum negativ.
% Die Enegie des Systems wird also durch die numerische Näherung (mit diesen Verf.)
% systematisch reduziert.

%--------------------------------------------------------------------------
%-- Differentialgleichungssystem  aus Beispiel 6.1 aufbauen ---------------

function [M, mu] = create_DES()
    
    syms x(t) y(t)
    mu=0.012277471;
    D_E = @(x, y) ((x+mu)^2+y^2)^(3/2);
    D_M = @(x, y) ((x-(1-mu))^2+y^2)^(3/2);
    de1 = diff(x,2) == x+2*diff(y)-(1-mu)*((x+mu)/D_E)-mu*((x-1+mu)/D_M);
    de2 = diff(y,2) == y-2*diff(x)-(1-mu)*(y/D_E)-mu*(y/D_M);
    des = [de1 de2];
    [V, S]=odeToVectorField(des);
    % Struktur von sym Y = [y; diff(y); x; diff(x)]
    M = matlabFunction(V, 'vars', {'t','Y'});

end

%--------------------------------------------------------------------------
%-------------------------  Verfahren und AWP definieren ------------------

function [AWP, RK4, gaussII] = createSetting()

    % AWP Definition
    x0 = [0.994; 0];
    y0 = [0; -2.00158510637908252240];
    AWP.Y0 = [y0;x0];
    AWP.t0 = 0;
    AWP.Tper = 17.065216560157962;
    
    %Gauss Two-Step - implizites Verf.
    gaussII.A = [[1/4,1/4 - sqrt(3)/6];[1/4 + sqrt(3)/6,1/4]];
    gaussII.beta = [1/2,1/2]';
    gaussII.gamma = [1/2-sqrt(3)/6,1/2 + sqrt(3)/6]';
    gaussII.Order = 4;
    gaussII.n = 605000;
    gaussII.Name = "gaussII";
    
    %klassisches RK - explizites Verf.
    RK4.A = [0 0 0 0; 1/2 0 0 0; 0 1/2 0 0; 0 0 1 0];
    RK4.beta = [1/6,1/3,1/3,1/6]';
    RK4.gamma = [0, 1/2, 1/2, 1]';
    RK4.Order = 4;
    RK4.n = 845200;
    RK4.Name = "klassisches Runge-Kutta Verfahren";

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