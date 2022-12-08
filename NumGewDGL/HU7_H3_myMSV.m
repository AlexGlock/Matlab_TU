%% numdgl Hausübung 7 Aufgabe H3 - myMSV solver
% Code von Alexander Glock und Luisa Emrich
clearvars


% Problemdefinition
f = @(t,y) [y(2); -pi^2*y(1)];
y0=[1; 0];
t0=0;
T=pi;

% Adams-Bashford 3 explizit
%beta=[23/12 -16/12 5/12 0];
%alpha=[0 0 -1 1];

% Adams-Bashford 2 explizit
%beta=[-1/2 3/2 0];
%alpha=[0 -1 1];

% Adams-Moulton 2 implizit
beta=[-1/12 8/12 5/12];
alpha=[0 -1 1];

%-----------------------------  solver  -----------------------------------
n=10000;
[t,y]=myMSV(alpha,beta,y0,f,t0,T,n);



% Plot der Lösung 
plot(t,y)
title('Lösung des AWP mit myMSV solver')
xlabel('Zeit t in Sekunden')
ylabel('y(t)')

%--------------------------------------------------------------------------
%-------------------------  myMSV solver ----------------------------------

function [t,y]=myMSV(alpha,beta,y0,f,t0,T,n)
    h=T/n;
    % Fallunterscheidung impl/ expl BEM beta(1) entspricht ana. b(0)
    if beta(end)==0 && alpha(end)==1
        is_imp=false;
    elseif alpha(end)==1
        is_imp=true;
    else
        msgbox("MSV ist nicht linear! Bitte alpha und beta prüfen")
    end

    l = length(beta)-1;
    % l-1 Zeitschritte mit RK4 berechnen = Anlaufrechnung
    % lösen mit RK4
    order = 4;
    impl= false;
    [A,b,g, ~]=select_RK(order,impl);
    T_esv = (l-1)*h;
    n_esv = l-1;
    [t_esv,y_esv, ~]=myRK(y0,f,t0,T_esv,n_esv,A,b,g);
        
    % Ergebnisverktoren an Anlaufrechnung anhängen
    n_msv=n-n_esv-1;
    y=cat(2,y_esv,zeros(length(y0),n_msv));
    t=cat(2,t_esv,zeros(1,n_msv));
    %--------------------------------------------------------------------------    
    % Löser für expliziten Fall -> b0 == 0 
    if not(is_imp) 
        % init f storage
        f_stor=zeros(length(y0),l-1);
        for ll=(1:1:l)
            ti_s=t0+ll*h;
            f_stor(:,ll)=f(ti_s,y(:,ll));
        end
        % calculate 
        b_flip=beta(1:end-1);
        a_flip=-alpha(1:end-1);
        for i = (l:1:n)
            ti=t0+i*h;
            fi=f(ti,y(:,i));
            % add new to f storage
            f_stor(:,l) = fi;
            % calc step
            y_sum = sum(a_flip.*y(:,(i-l+1:i)),2);
            y(:,i+1)=y_sum+h*sum(b_flip.*f_stor(:,(1:end)),2);
            t(i+1)=ti;
            % remove oldest storage
            f_stor = f_stor(:,(2:end)); 
        end
    %--------------------------------------------------------------------------
    % Löser für impliziten Fall -> b0 != 0
    else
        % init f storage
        f_stor=zeros(length(y0),l-1);
        for ll=(1:1:l)
            ti_s=t0+ll*h;
            f_stor(:,ll)=f(ti_s,y(:,ll));
        end
        % calculate 
        b_flip=beta(1:end-1);
        a_flip=-alpha(1:end-1);
        for i = (l:1:n)
            ti=t0+i*h;
            fi=f(ti,y(:,i));
            % add new to f storage
            f_stor(:,l) = fi;
            % calc step
            y_sum =sum(a_flip.*y(:,(i-l+1:i)),2);
            options = optimset('Display','off');
            fun=@(ys) y_sum+h*(sum(b_flip.*f_stor(:,(1:end)),2)+ beta(end)*f(ti,ys))-ys;
            y(:,i+1) = fsolve(fun,y(:,i),options);
            t(i+1)=ti;
            % remove oldest storage
            f_stor = f_stor(:,(2:end)); 
        end
    end
end

%--------------------------------------------------------------------------
%-------------------------  ESV definieren --------------------------------

function [A,beta,gamma,n_max]=select_RK(order,impl_true)
    if order==4 && impl_true
        %Gauss Two-Step - implizites Verf.
        A = [[1/4,1/4 - sqrt(3)/6];[1/4 + sqrt(3)/6,1/4]];
        beta = [1/2,1/2]';
        gamma = [1/2-sqrt(3)/6,1/2 + sqrt(3)/6]';
        n_max=200000;
    
    elseif order==3 && impl_true
        %Radau-IA
        A = [1/4 -1/4; 1/4 5/12];
        beta = [1/4 3/4]';
        gamma = [0 2/3]';
        n_max=300000;
    
    elseif order==4 && not(impl_true)
        %klassisches RK - explizites Verf.
        A = [0 0 0 0; 1/2 0 0 0; 0 1/2 0 0; 0 0 1 0];
        beta = [1/6,1/3,1/3,1/6]';
        gamma = [0, 1/2, 1/2, 1]';
        n_max=1000000;
    
    elseif order==3 && not(impl_true)
        %Simpsonregel- explizites Verf.
        A = [0 0 0 ; 1/2 0 0 ; -1 2 0];
        beta = [1/6,4/6,1/6]';
        gamma = [0, 1/2, 1]';
        n_max=6000000;
    else
        error("unsupported order-impl combination")
    end
    
end

%--------------------------------------------------------------------------
%-- Differentialgleichungssystem  aus Beispiel 6.1 aufbauen ---------------

function [f, AWP] = create_DES()
    
    syms x(t) y(t)
    mu=0.012277471;
    D_E = @(x, y) ((x+mu)^2+y^2)^(3/2);
    D_M = @(x, y) ((x-(1-mu))^2+y^2)^(3/2);
    de1 = diff(x,2) == x+2*diff(y)-(1-mu)*((x+mu)/D_E)-mu*((x-1+mu)/D_M);
    de2 = diff(y,2) == y-2*diff(x)-(1-mu)*(y/D_E)-mu*(y/D_M);
    des = [de1 de2];
    [V, ~]=odeToVectorField(des);
    % Struktur von sym Y = [y; diff(y); x; diff(x)]
    f = matlabFunction(V, 'vars', {'t','Y'});
   
    % AWP Definition
    x0 = [0.994; 0];
    y0 = [0; -2.00158510637908252240];
    AWP.y0 = [y0;x0];
    AWP.t0 = 0;
    AWP.T = 17.065216560157962;
    
end

%--------------------------------------------------------------------------
%------------------------------ myRK-solver -------------------------------

function [t,y,f_calls] = myRK(y0,f,t0,T,n,A,b,g)
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


f_calls=0;
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




    for curStep = 2:n+1
        
      [y(:,curStep), f_calls] = generalRKStep(t(curStep-1),y(:,curStep-1),A,g,b,h,f,options,f_calls);
    
      t(curStep) = t(curStep-1)+h;
    end

end

%% Auxillary Function below for single-step RKV

function [yi, f_calls] = generalRKStep(ti,yi,A,gamma,beta,h,f,options, f_calls)
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
                f_calls=f_calls+1;
                g(i,:) = g(i,:) + h.*A(i,j).*fbar(ti + h*gamma(j),g(j,:));
            end
        end
        for i = 1:numStages
            f_calls=f_calls+1;
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

