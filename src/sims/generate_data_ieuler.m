addpath('mhd')
addpath('navier_stokes')

butcher = struct();
butcher.A = [1];
butcher.b = [1];

n = 64;
timesteps = 100;
nu = 2^-6;
eta = 2^-6;
skip = 10;

truncate = 1;

# actual data
rng(1);
omega = @(x, y) rand(n,n);

rng(2);
A = @(x, y) rand(n,n);
forcing = @(x, y) sin(x) + sin(y);

dts = linspace(0.001, 0.005, 9);

for i = 1:length(dts)
    [x, y, t, u, v, p, Bx, By] = mhd_rk_implicit(omega, A, forcing, n, dts(i), timesteps, nu, eta, 1, skip, butcher, true);
    result = struct("x",x,"y",y,"t",t(truncate:end),"u",u(:,:,truncate:end),"v",v(:,:,truncate:end),"p",p(:,:,truncate:end),"Bx",Bx(:,:,truncate:end),"By",By(:,:,truncate:end));
    save("-hdf5", ["data-" num2str(dts(i)) ".hdf5"], "result");
end