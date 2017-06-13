function t = getT1T2FromConc(sol, nicl, mncl, agar)
    A = [sol(1:3)' ; sol(4:6)'];
    k = [nicl ; mncl; agar];
    r = A*k + [sol(7);sol(8)];
    t = 1000./r;
end