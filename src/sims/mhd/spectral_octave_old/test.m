    n = 32;
    
    k = 0:n-1;
    k(k>n/2) = k(k>n/2) - n;

    mask = abs(k) <= n/3;
    mask = mask & mask.';

    kx = k;
    ky = k.';

    kx