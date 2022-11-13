
plot_costs(4);

function x = Thomas_solve(A,d)

    [a,b,c] = get_Vectors(A);
    n= size(b);
    
    for j=2 :n
        w = a(j)/b(j-1);
        b(j) = b(j) -w*c(j-1);
        d(j) = d(j) -w*d(j-1);
    end
    
    x = d./b;
    
    for j=(n-1):-1:1
        x(j)=(d(j)-c(j)*x(j+1))/b(j);
    end
end

function [a,b,c] = get_Vectors(A)
        n = size(A); % get dimension of Matrix
        n = n(1);
        a = zeros(n,1);
        b = zeros(n,1);
        c = zeros(n,1);

        for i=1 :n
            b(i) = A(i,i); 
            if i>= 2
                a(i) = A(i,i-1);
            end
            if i<= (n-1)
                c(i) = A(i,i+1);
            end
        end
end

function [A, Asparse, d] = setup_Data(k)

    n = 10^k;
    A = zeros(n);
    d = zeros(n,1);
    for j=1 :n
        if j<n
            A(j+1,j)=1; % setze a
            A(j,j+1)=1; % setze c
        end
        A(j,j)=2;       %setze b
        d(j)=2;         %setze d
    end
    Asparse = sparse(A);
end

function R = compare_Mehods(A, Asparse, d)
R=zeros(3,2); n=length(d);

x = zeros(n,1);
tStart_Backslash_Full = tic;
    x=A\d;
tEnd_Backslash_Full = toc(tStart_Backslash_Full);
R(1,1) = tEnd_Backslash_Full;

x = zeros(n,1);
tStart_Backslash_Sparse = tic;
    x=Asparse\d;
tEnd_Backslash_Sparse = toc(tStart_Backslash_Sparse);
R(1,2) = tEnd_Backslash_Sparse;

x = zeros(n,1);
tStart_Inv_Full = tic;
    A_inv = inv(A);
    x = A_inv*d;
tEnd_Inv_Full = toc(tStart_Inv_Full);
R(2,1) = tEnd_Inv_Full;

x = zeros(n,1);
tStart_Inv_Sparse = tic;
    A_inv_sparse = inv(Asparse);
    x = A_inv_sparse*d;
tEnd_Inv_Sparse = toc(tStart_Inv_Sparse);
R(2,2) = tEnd_Inv_Sparse;
    
x = zeros(n,1);
tStart_Thomas_Full = tic;
    x = Thomas_solve(A,d);
tEnd_Thomas_Full = toc(tStart_Thomas_Full);
R(3,1) = tEnd_Thomas_Full;

x = zeros(n,1);
tStart_Thomas_Sparse = tic;
    x = Thomas_solve(Asparse,d);
tEnd_Thomas_Sparse = toc(tStart_Thomas_Sparse);
R(3,2) = tEnd_Thomas_Sparse;
end

function plot_costs(k)
    t = zeros(k,6);
    for i=1 :k
        [A, Asparse, d] = setup_Data(i);
        R = compare_Mehods(A, Asparse, d);
        t(i,1)=R(1,1);
        t(i,2)=R(1,2);
        t(i,3)=R(2,1);
        t(i,4)=R(2,2);
        t(i,5)=R(3,1);
        t(i,6)=R(3,2);
    end
    x = logspace(1,4,4);
    loglog(x,t(:,1),x,t(:,2),x,t(:,3),x,t(:,4),x,t(:,5),x,t(:,6))
    grid on
    legend('Back-Full','Back-Sp','InvFull','InvSp','ThFull','ThSp')
end