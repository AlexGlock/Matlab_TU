 h  = 0.01;  % Schrittweite 
 T  = 1;     % Endzeit 
 tt = 0:h:T; % Zeitgitter

 % Problemdefinition 
 f   = @(t,y) [y(2); -pi^2*y(1)]; y0=[1; 0]; 
 yex = @(t) [cos(pi*t); -pi*sin(pi*t)];

 yy = y0; y=y0; 
% num. Approximation mit implizitem Euler Verfahren
for t=tt(2:end)
    fun = @(y1) [y-y1-h*(f(t+h,y1))]; % impliziter Euler
    [y,fval]=fsolve(fun,y0);
    yy = [yy,y];
end

 % exakte Lösung
 yyex=yex(tt);

 % erste Komponente der Lösung
 plot(tt,yyex(1,:),'b-',tt,yy(1,:),'r.')

 % maximaler Fehler in der ersten Komponente
 ee = yy(1,:)-yyex(1,:);
 fprintf('h=%f, err=%f\n',h,norm(ee,inf));