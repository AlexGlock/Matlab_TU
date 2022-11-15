%% HU3 Aufgabe H4 - RKV
clearvars

% Eingabeparam
A = [1/2 0; 1/2 0] %[1/2 -1/2; 1/2 1/2] % [0 0 0 0; 1/3 0 0 0; -1/3 1 0 0; 1 -1 1 0] 
b = [1/2 1/2] % [1/2 1/2] % [1/8 3/8 3/8 1/8] 
g = [0 1] % [0 1] % [0 1/3 2/3 1] 

% Problemdefinition
f = @(t,y) [y(2); -pi^2*y(1)];
yex = @(t) [cos(pi*t); -pi*sin(pi*t)];
y0=[1; 0];
t0=0;
T=1;
%n=30; n Vorgabe durch for loop

hmax = 0.1; % maaximale Zeitschrittweite
hmin = 0.001; % minimale Zeitschrittweite
hh = flip(hmin:hmin*15:hmax); % Zeitschwrittgrößen von hmax bis 0.0001
cc = 0; % konstante initialisieren
% Iteration über kleiner werdendes h
for h = hh(2:end)
    % Reset der Startwerte/ Vektoren
    yy=y0; 
    y=y0;
    n=T/h;
    
    % num. Approximation mit RKV
    [tt,yy]=myRK(y0,f,t0,T,n,A,b,g);
    % exakte Lösung ausgewertet am Zeitgitter
    yyex=yex(tt);
    % norm. Fehler der num. Approximation
    c = norm((yyex(1,:)-yy(1,:)),'Inf');
    cc = [cc,c];
end

% Konvergenz plot
loglog(hh,cc)
title('Konvergenzverhalten - RKV')
xlabel('h')
ylabel('|| y-y_h ||')
grid on

% exakte Lösung
%yyex=yex(t);
% Plot der ersten Komponente der Lösung
%plot(t,yyex(2,:),'b-',t,yy(2,:),'r.')
%title('numerische Lösung durch gegebenes RKV')
%legend('exakte LSG y(t)','RKV LSG y_h(t)')
%xlabel('Zeit t')
%ylabel('y(t)')


function [t,yy] = myRK(y0,f,t0,T,n,A,b,g)
    y=y0;
    yy=y0;
    % Stufenzahl + vektor
    s = abs(size(A,1));
    ss = (1:1:s);
    % äqidistante Zeitschritte
    h=(T-t0)/n;
    t =(t0:h:T);

    % strikte Dreiecksmatrix?
    if istril(A) & (norm(diag(A))==0)
        % -- bei explizitem Verfahren -- FUNKTIONIERT %
        for ti=t(2:end)
            S_bk=zeros(1,size(y0,1)).';
            KK = f(ti,y);                  % Liste mit allen k [k_1, k_2, ...]
            % k's höherer Stufe berechnen:
            for j=ss(2:end)
                S_ak=zeros(1,size(y0,1)).';
                for k=ss(1:j-1)   % Summe über alpha_k * k_k = S_ak
                    S_ak = S_ak+(A(j,k)*KK(:,k));
                end    
                KK=[KK,f(ti+g(j)*h,y+h*S_ak)]
            end
            for j=ss    % Summe über beta_j * k_j
                S_bk = S_bk+b(j)*KK(:,j);
            end
            %S_bk = symsum(b(l)*K(:,l),l,1,2)
            y = y+h*S_bk;
            yy = [yy,y];
        end

    else % implizite RKV

        % -- Matrix ist diagonal implizit -- FUNKTIONIERT%
        if istril(A) 
            for ti=t(2:end)
                S_bk=zeros(1,size(y0,1)).';
                % K = Liste mit allen k [k_1, k_2, ...]
                fun= @(sk) [sk-f(ti+g(1)*h,y+h*A(1,1)*sk)];
                [kval, ~]=fsolve(fun,y0);
                KK = kval;
                % k's höherer Stufe berechnen:
                for j=ss(2:end)
                    S_ak=zeros(1,size(y0,1)).';
                    for k=ss(1:j)   % Summe über alpha_k * k_k = S_ak
                        if k==j
                            S= @(sk) [y+h*(S_ak+A(j,k)*sk)];
                        else
                            S_ak = S_ak+(A(j,k)*KK(:,k));
                        end
                    end
                    fun= @(sk) [sk-f(ti+g(j)*h,S(sk))];
                    [kval, ~]=fsolve(fun,KK(:,j-1));
                    KK=[KK,kval]
                end
                for j=ss    % Summe über beta_j * k_j
                    S_bk = S_bk+b(j)*KK(:,j);
                end
                %S_bk = symsum(b(l)*K(:,l),l,1,2)
                y = y+h*S_bk;
                yy = [yy,y];
            end 

        % -- Matrix hat keine besondere symmetrie/aufbau -- FUNKTIONIERT NICHT %  
        else 
            % s dimensionales vektorsymbol k = [k1 k2 ... ks]
            sym('Ks',[1 s])
            sym('F',[1 s])
            
            for ti=t(2:end)
                S_bk=zeros(1,size(y0,1)).';
                
                for j=ss(2:end)
                    S_aK=zeros(1,size(y0,1)).';
                    for k=ss(1:j)   % Summe über alpha_k * k_k = S_ak
                       S_aK= @(Ks) [y + h*(S_ak+A(j,k)*Ks(k))];
                       if k==1
                        F= @(Ks) [Ks(k)-f(ti+g(j)*h,y+h*S_aK(Ks))];
                       else
                        Fnew = @(Ks) [Ks(k)-f(ti+g(j)*h,y+h*S_aK(Ks))];
                        F = @(Ks) {F,Fnew}
                       end
                    end
                end
                Ks0=zeros(1:s)
                [KK, ~]=fsolve(F,Ks0)
                for j=ss    % Summe über beta_j * k_j
                    S_bk = S_bk+b(j)*KK(:,j);
                end
                %S_bk = symsum(b(l)*K(:,l),l,1,2)
                y = y+h*S_bk;
                yy = [yy,y];
            end 
       
        end
    end
end

