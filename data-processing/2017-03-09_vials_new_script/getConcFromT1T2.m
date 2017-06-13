function [conc, n, c_min] = getConcFromT1T2(sol, T1, T2)
    R1 = 1000./T1;
    R2 = 1000./T2;
    A = [sol(1:3)' ; sol(4:6)'];
    b = [R1 - sol(7); R2 - sol(8)];
    
    % Initial solution
    conc = A \ b;
    n = null(A);
    c = -conc./n;
    c_min = min(c);
    
    %LB = zeros(3,1);
    %w = [1, 2, 1000]';
    %x_cost = linprog(w, [], [], A, b, LB, [])
end