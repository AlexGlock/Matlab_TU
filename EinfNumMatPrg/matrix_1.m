function M1 = matrix_1(n)

assert(n<=1000);

T1 = -5 + 10*rand(n);
T2 = T1 - diag(diag(T1));
M1 = T2 + (-1).^(1:n).*diag(1.1*sum(abs(T2),2));

end