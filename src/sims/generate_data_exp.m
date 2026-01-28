clear;
pkg load symbolic;

addpath('mhd')
addpath('navier_stokes')


n = 64;
dt = 0.01;
timesteps = 100;
nu = 0.002;
eta = 0.01;
k1 = 1.0;
k2 = 1.0
skip = 10;

truncate = 1;

# random noise
rng(1)
omega = @(x, y) rand(n,n) #1/k1 * cos(k1*x) + 1/k2 * cos(k2*y)

rng(2)
A = @(x, y) rand(n,n);
forcing = @(x, y) 0;

butcher_rk2 = [
    0 0 0;
    1/2 1/2 0;
    0 0 1]
butcher_rk3 = [
    0 0 0 0;
    1/2 1/2 0 0;
    3/4 0 3/4 0;
    0 2/9 1/3 4/9
]
butcher_rk4 = [
    0 0 0 0 0;
    1/2 1/2 0 0 0;
    1/2 0 1/2 0 0;
    1 0 0 1 0;
    0 1/6 1/3 1/3 1/6
]
butcher_rk5 = [
    0 0 0 0 0 0 0;
    1/3 1/3 0 0 0 0 0;
    2/5 4/25 6/25 0 0 0 0;
    1 1/4 -3 15/4 0 0 0;
    2/3 2/27 10/9 -50/81 8/81 0 0;
    4/5 2/25 12/25 2/15 8/75 0 0;
    0 23/192 0 125/192 0 -27/62 125/192
]

schemes={'RK2E', 'RK3E', 'RK4E', 'RK5E'};

butchers={butcher_rk2, butcher_rk3, butcher_rk4, butcher_rk5};

for i = 1:length(schemes)
    [x, y, t, u, v, p, Bx, By] = mhd_rk_explicit(omega, A, forcing, n, dt, timesteps, nu, eta, 1, skip, butchers{i}, true);
    result = struct("x",x,"y",y,"t",t(truncate:end),"u",u(:,:,truncate:end),"v",v(:,:,truncate:end),"p",p(:,:,truncate:end),"Bx",Bx(:,:,truncate:end),"By",By(:,:,truncate:end));
    save("-hdf5", ["data-" schemes{i} ".hdf5"], "result");
end