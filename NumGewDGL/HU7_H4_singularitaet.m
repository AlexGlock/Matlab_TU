%% numdgl Hausübung 7 Aufgabe H4 - Singularitäten
% Code von Alexander Glock und Luisa Emrich
clearvars
% In den DGLs wurde an Stelle von der Vairablen x die Variable y genutzt,
% deshalb ist im Plot auch die Lösung mit y(t) beschrieben...

% Ergebnisse:
% Die DGL aus Aufgabenteil (i) lässt sich mit dem DOPRI Verfahren bis zu
% einer Dauer von T~2 gut lösen. Dabei wird jedoch in den Plots schon 
% ersichtlich, dass größere Lösungszeiträume poblematisch werden da die
% Lösungsfunktion y(t) stark gegen null abfällt und die akzeptierten Schrittweiten
% deshalb zunehmend kleiner werden bis kein Fortschritt in der Zeit mehr gemacht wird und der Solver stecken bleibt.
%
% Die DGL aus Aufgabenteil (ii) besitzt bereits um T~0.7 ihren kritischen
% Punkt, ab dem die Lösungsfunktion so stark absinkt das die akzeptierten Schrittweiten
% schnell abfallen und der Solver nicht mehr weiter kommt.
%
% Lediglich die DGL (iii) lässt sich mit dem Solver über einen sehr langen Zeitrum lösen (T>50),
% da ihre (vektorwertige) Lösungsfunktion stetig und mit einer endlichen Steigung wächst und
% die Solver-Schrittweiten deshalb groß genug bleiben um von der
% Maschinengenauigkeit aufgelöst werden zu können, woduch noch effektive Fortschritte
% im Lösungszeitraum für große T gemacht werden.
%
% ZUSAMMENFASSUNG:
% Differentialgleichungen mit starkem Dämpfungscharakter in ihrer
% Lösungsfunktion (Vgl. Teil i und ii) lassen sich mit einem variablen
% Schrittweitensolver nur bedingt lösen, da die Schrittweite ab einem gewissen Punkt
% nicht mehr aufgelöst werden kann und es zum Stillstand kommt.
% DGL Systeme mit (stetig) wachsenden Lösungen sind dagegen gut lösbar
% und unproblematisch für den Schrittweitenalgorithmus.


%------------------- DGL-System Auswahl ----------------------------------

[f, ~, y0, t0, T]=create_DES_i()
%[f, ~, y0, t0, T]=create_DES_ii()
%[f, S, y0, t0, T]=create_DES_iii()

%------------------- DOPRIsolver settings ---------------------------------
h=1e-5;
facMin = 0.8;
facMax = 1.2;
fac = 0.9; %Gives 4 Rejects, 603 Steps
tol = 1e-8;
%------------------- lösen und plotten ------------------------------------

[t,y, hh]=myDOPRI(y0,f,t0,T,tol,facMin,facMax,fac,h);

plot(t,y)
title('Lösung des AWP durch DOPRI')
xlabel('Zeit t in Sekunden')
ylabel('y(t)')

figure
plot(t,hh)
title('Schrittweiten über Lösungszeitraum')
xlabel('Zeit t in Sekunden')
ylabel('Schrittweite h')

%--------------------------------------------------------------------------
%--------------- i) Differentialgleichungssystem  y' = -y^(1/2)  ----------

function [f, S, y0, t0, T] = create_DES_i()
    
    syms y(t)
    de = diff(y) == -y^(1/2);
    [V, S]=odeToVectorField(de);
    % Struktur von sym Y = [y; diff(y)]
    f = matlabFunction(V, 'vars', {'t','Y'});
    y0=1;
    t0=0;
    T=2;

end

%--------------------------------------------------------------------------
%-------------- ii) Differentialgleichungssystem  y' = sin(1/y)-2  --------

function [f, S, y0, t0, T] = create_DES_ii()
    
    syms y(t)
    de = diff(y) == sin(1/y)-2;
    [V, S]=odeToVectorField(de);
    % Struktur von sym Y = [y; diff(y)]
    f = matlabFunction(V, 'vars', {'t','Y'});
    y0=1;
    t0=0;

    T=0.75;

end

%--------------------------------------------------------------------------
%------------- iii) Differentialgleichungssystem  y'' = 5(1 − y^2)y' + y --

function [f, S, y0, t0, T] = create_DES_iii()
    
    syms y(t)
    de = diff(y, 2) == 5*(1-y^2)*diff(y)+y;
    [V, S]=odeToVectorField(de);
    % Struktur von sym Y = [y; diff(y)]
    f = matlabFunction(V, 'vars', {'t','Y'});
    y0 = [0.1; 0];
    t0=0;
    T = 15;

end

%--------------------------------------------------------------------------
%-------------------------  myDOPRI solver --------------------------------

function [t,y, hh]=myDOPRI(y0,f,t0,T,tol,facMin,facMax,fac,h)

    %% Input
        %y0: Vertical Vector containing initial values
        % f function handle of the form f= @(t,y)
        % t0: initial Time
        % T: stopping Time
        % tol: desired accuracy
        % facMin <1, facMax>1, fac <1: Parameters to adapt step size
        % h: initial h to start with

%%
   % A Simple way of handling dynamic memory allocation
   % This is not the best way of doing things, but a simple improvement
   % For our purposes we will still try and avoid the matlab syntax of yy =
   % [yy,y_new] as this is very slow. We will instead allocate some initial
   % amount of points and once these run out we can then append new points
   % (get more memory).
   initSize = 100;
   
   
%Preallocate the arrays for t,y
t = zeros(1,initSize);
hh = zeros(1,initSize);
t(1) = t0;

y = zeros(length(y0),initSize);
y(:,1) = y0; 
hh(1)=h;
szRHS = length(y0);
%% Dopri Tableau definition
dopri.A = [[0,0,0,0,0,0,0];
           [1/5,0,0,0,0,0,0];
           [3/40,9/40,0,0,0,0,0];
           [44/45,-56/15,32/9,0,0,0,0];
           [19372/6561,-25360/2187,64448/6561,-212/729,0,0,0];
           [9017/3168,-355/33,46732/5247,49/176,-5103/18656,0,0]];
%Note that we actually do not need to keep the 7th Row in the DOPRI butcher
%tableau. This is because g6 is actually the 4th order embedded method
%already and as such there is no need to calculate it twice as one would do
%classically for a RK method
dopri.beta5 = [35/384,0,500/1113,125/192,-2187/6784,11/84,0]'; %5th Order accurate Solution
dopri.beta4 = [5179/57600,0,7571/16695,393/640,-92097/339200,187/2100,1/40]'; %4th Order accurate Solution
dopri.gamma = [0,1/5,3/10,4/5,8/9,1,1]';


%% Some internal housekeeping for variables we will need


tCur = t0;             %Current Time
yCur = y0;             %Current point in Solution
hCur = h;              %Current step width
totalRejects = 0;      %Stats on how much is rejected. Remove here and below in main 
                       %Loop if not needed.
totalSteps = 1;        %Stats on total steps taken
sizeLeft = initSize-1; %We already filled the first entry in y
                       %Keeps track of how much preallocated Memory we have
                       %left
numReallocs = 1;       %Helper Variable to allocate more memory per allocation 
                       %if many allocations are needed
                       
                       
    endFlag = 1;      %signalling variable for when to exit the main program loop
   
%% Main Program Loop
    while endFlag
       
        if sizeLeft == 0
            %If we have no preallocated memory left we are forced to
            %reallocate some memory
            y = [y,zeros(length(y0),initSize*2^(numReallocs-1))];
            t = [t,zeros(1,initSize*2^(numReallocs-1))];
            hh = [hh,zeros(1,initSize*2^(numReallocs-1))];
            %
            %NumReallocs serves as a way to allocate progressively bigger
            %amounts of Memory as we need more reallocations. This is a
            % simple heuristic because a problem needing many Reallocations is more
            %likely to need even bigger amounts of memory than we have
            %supplied in the previous reallocation of memory 
            %This way, we only need logarithmically many reallocations.
            numReallocs = numReallocs +1; 
            sizeLeft = initSize*2^(numReallocs-1);
        end
        
        
         if tCur +hCur  >T%This really simple check serves to give the last Step
                         %of the DOPRI method at EXACTLY the T specified,
                         %not only somewhere *around* the specified T
                         
            [isFine,y5,hCurInt] = verifyEnd(yCur,f,tCur,T,tol,facMin,facMax,fac,szRHS,dopri);
                         
                         %Then do the Regular keep/don't keep process
            if isFine
                totalSteps = totalSteps +1;
                
                t(totalSteps) = tCur + hCurInt;

                hh(totalSteps) = hCurInt;
                
                y(:,totalSteps) = y5; %We will in practice want to use the
                % 5th order accurate approximator for our
                % solution even though we control the
                % time stepping using the 4th Order
                % estimator. This is because BOTH
                % solutions are guaranteed to be within
                % tol, and y5 is on average a better
                % solution
                break
                
            else
                totalRejects = totalRejects +1; 
            end
         end
        
         %singleDopriStep calculates the 4th and 5th Order solutions for a
         %given h
        [y5,y4] = singleDopriStep(yCur,f,tCur,hCur,szRHS,dopri);
        
        
        %Check if step is okay and also update hCur with the new optimal h
        [keepStep, hNew] = acceptReject(y5,y4,hCur,tol,facMin,facMax,fac);        
              
        if keepStep
            %Increment all the variables that are just there for
            %housekeeping
            totalSteps = totalSteps +1;
            sizeLeft = sizeLeft-1;
            t(totalSteps) = tCur + hCur;
            tCur = t(totalSteps);
            hh(totalSteps) = hCur;
            y(:,totalSteps) = y4; %We will in practice want to use the
                                  % 5th order accurate approximator for our
                                  % solution even though we control the
                                  % time stepping using the 4th Order
                                  % estimator. This is because BOTH
                                  % solutions are guaranteed to be within
                                  % tol, and y5 is on average a better
                                  % solution
          
          yCur = y4;
          hCur = hNew;
        else
            %Only propose a new h, but do not accept the given step due to
            %error bounds.
            hCur = hNew;
            totalRejects = totalRejects +1; %Interesting diagnostic
        end
        
      
    end
    %% Cleanup overallocated arrays and we are done
    
    %We only need to keep entries that contain actual steps. Due to the way
    %we handle dynamic memory allocation, our arrays for t and y will have
    %some empty entries that we can discard now that we are finished
    t = t(1:totalSteps);
    hh = hh(1:totalSteps);
    y = y(:,1:totalSteps);
    disp("Total Rejects:");
    disp(totalRejects);

end

function [y5,y4] = singleDopriStep(yCur,f,tCur,hCur,szRHS,meth)

%Method for calculating a single set of 4th and 5th order solutions from a
%given current solution and step size

 kVec = zeros(szRHS,7); %Again, note that we only really need 6 k and not the
                        %full 7 k from the DOPRI tableau because g7 = 5th
                        %order accurate solution already.
                        %We still allocate 7 rows because we will need one
                        % extra later.
                        
                        %Also, we choose the k-Form because it helps save
                        %on function evaluations (minimum: 7 Evaluations)
    kVec(:,1) = f(tCur,yCur);
    
    %Calculate all the Stage Values using k-Form
    for stage = 2:6
        stageVal = 0;
       for colInd = 1:(stage-1)
           %Note that we only sum the alpha_ij*kj part first and then
           %multiply by h once. In theory this should help save on roundoff
           %error and a few floating point evaluations but nothing major
           stageVal = stageVal + meth.A(stage,colInd).*kVec(:,colInd);
       end
       kVec(:,stage) = f(tCur + meth.gamma(stage)*hCur,yCur + hCur.*stageVal);
    end
    
    %Calculate the 5th order accurate solution
    %Again, we only multiply by h once to reduce roundoff errors, which
    %makes the order of calculations somewhat unfamiliar.
    y5 = 0;
    for stage = 1:6
       y5 = y5 + meth.beta5(stage).*kVec(:,stage); 
    end
    y5 = hCur.*y5 + yCur; %This is now essentially "g7" and also the 5th order accurate solution
    
    %This is where we need the extra row in kVec:
    kVec(:,7) = f(tCur + meth.gamma(7)*hCur,y5);
    
    %Calculate 4th Order accurate solution
     y4 = 0;
    for stage = 1:7
       y4 = y4 + meth.beta4(stage).*kVec(:,stage); 
    end
    y4 = hCur.*y4 + yCur;
end

function [isOkay, hNew] = acceptReject(y5,y4,hCur,tol,facMin,facMax,fac)
    %This function checks a given proposed set of 4th and 5th order
    %Solutions against the error bounds and proposes a new step size given
    %the current error estimate

        tauEstimate = norm(y5-y4,"inf")/hCur; % Estimator for current truncation error
                                              % We Choose the infinity norm
                                              % here, but the choice of
                                              % norm is arbitrary and you
                                              % can use any other ones
                                              % you'd like as well.
        
        if tauEstimate < tol 
            %Accept Step
            isOkay = 1;
        else
            %Reject Step
            isOkay = 0;
        end
        
        %Propose new h, regardless if accepted or rejected step size
        hNew = hCur * min(facMax,max(facMin,fac *(tol/tauEstimate)^(1/4)));
        %We will be keeping the 5th order accurate solution but our error
        %estimator is indeed actually only given for the 4th order
        %accurate solution so we have p = 4 and NOT p = 5
        
        %Using the estimate for the truncation error instead of the
        %estimate for the step error will on average lead to more frequent
        %step rejection, but it improves accuracy as well.
     

end

function [isFine,y,hCur] = verifyEnd(yCur,f,tCur,T,tol,facMin,facMax,fac,szRHS,dopri)
            %Method to check wether the last step to end at exactly T is in
            %fact admissible or not.
 
            hCur = T-tCur;
            
            %Do another step with the exact h needed to reach the exact
            %endpoint
            [y5,y4] = singleDopriStep(yCur,f,tCur,hCur,szRHS,dopri);
            %Do a verification it is acceptably accurate
            [isFine, ~] = acceptReject(y5,y4,hCur,tol,facMin,facMax,fac);
            
           
            y = y5;                 
         
end
