%% HU3 Aufgabe H4 - RKV
clearvars
% Eingabeparam
%A = [0 0; 1/2 0]
%b = [0 1]
%g = [0 1/2]
A = [0 0 0; 1/3 0 0; 0 2/3 0]
b = [1/4 0 3/4]
g = [0 1/3 2/3]
g(2)
% Problemdefinition
f = @(t,y) [y(2); -pi^2*y(1)];
yex = @(t) [cos(pi*t); -pi*sin(pi*t)];

y0=[1; 0];
yy=y0; 
y=y0;

t0=0;
T=5;
n=200;


% Stufenzahl + vektor
s = size(A,1);
ss = (1:1:s);
% äqidistante Zeitschritte
h=abs((T-t0)/n);
t =(t0:h:T);

% strikte Dreiecksmatrix? -> explizit
if istril(A) & (norm(diag(A))==0)
    for ti=t(2:end)
        S_bk=zeros(1,size(y0,1)).';
        G = y;                  % Liste mit allen g. g1 = y
        % g's höherer Stufe berechnen:
        for j=ss(2:end)
            S_ak=zeros(1,size(y0,1)).';
            % bestimmen von Summe über alpha * k = S_ak
            for k=ss(1:ss)
                S_ak = S_ak+A(j,k)*f(ti+g(k)*h,G(:,k));
            end    
            G = [G,y+h*S_ak]
            
        end
        for j=ss    % Summe über beta_j * k_j
            S_bk = S_bk+b(j)*f(ti+g(j)*h,G(:,j));
        end
        %S_bk = symsum(b(l)*K(:,l),l,1,2)
        y = y+h*S_bk;
        yy = [yy,y];
    end

else % implizite RKV
    syms 
    for ti=t(2:end)
        S_bk=zeros(1,size(y0,1)).';
        G = y;                  % Liste mit allen g. g1 = y
        % g's höherer Stufe berechnen:
        for j=ss(2:end)
            S_ak=zeros(1,size(y0,1)).';
            % bestimmen von Summe über alpha * k = S_ak
            for k=ss(1:ss)
                S_ak = A(j,k)*f(ti+g(k)*h,G(:,k));
            end    
            G = [G,y+h*S_ak]
            
        end
        for j=ss    % Summe über beta_j * k_j
            S_bk = S_bk+b(j)*f(ti+g(j)*h,G(:,j));
        end
        %S_bk = symsum(b(l)*K(:,l),l,1,2)
        y = y+h*S_bk;
        yy = [yy,y];
    end
    
end

% exakte Lösung
% 1
yyex=yex(t);
%c = norm((yyex(1,:)-yy(1,:)),Inf)/h
% erste Komponente der Lösung
plot(t,yyex(1,:),'b-',t,yy(1,:),'r.')
title('numerische Lösung durch gegebenes RKV')
legend('exakte LSG y(t)','RKV LSG y_h(t)')
xlabel('Zeit t')
ylabel('y(t)')