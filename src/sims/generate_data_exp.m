clear;
pkg load symbolic;

addpath('mhd')
addpath('navier_stokes')


n = 64;
dt = 2^-6;
timesteps = 200;
nu = 2^-6;
eta = 2^-6;
skip = 10;

truncate = 101;

# random noise
rng(1)
omega = @(x, y) rand(n,n) #1/k1 * cos(k1*x) + 1/k2 * cos(k2*y)

rng(2)
A = @(x, y) rand(n,n);
forcing = @(x, y) 0;

# https://www.ams.org/journals/mcom/1962-16-080/S0025-5718-1962-0150954-0/S0025-5718-1962-0150954-0.pdf
butcher_rk2 = [
    0 0 0;
    2/3 2/3 0;
    0 1/4 3/4]
butcher_rk3 = [
    0 0 0 0;
    1/2 1/2 0 0;
    3/4 0 3/4 0;
    0 2/9 1/3 4/9
]
rt5 = sqrt(5)
butcher_rk4 = [
    0 0 0 0 0;
    2/5 2/5 0 0 0;
    (14-3*rt5)/16 (-2889+1428*rt5)/1024 (3785-1620*rt5)/1024 0 0;
    1 (-3365+2094*rt5)/6040 (-975-3046*rt5)/2552 (467040+203968*rt5)/(240845) 0;
    0 (263+24*rt5)/1812 (125-1000*rt5)/3828 (3426304+1661952*rt5)/(5924787) (30-4*rt5)/123
]

schemes={'RK2', 'RK3', 'RK4'};

butchers={butcher_rk2, butcher_rk3, butcher_rk4};

for i = 1:length(schemes)
    [x, y, t, u, v, p, Bx, By] = mhd_rk_explicit(omega, A, forcing, n, dt, timesteps, nu, eta, 1, skip, butchers{i}, true);
    result = struct("x",x,"y",y,"t",t(truncate:end),"u",u(:,:,truncate:end),"v",v(:,:,truncate:end),"p",p(:,:,truncate:end),"Bx",Bx(:,:,truncate:end),"By",By(:,:,truncate:end));
    save("-hdf5", ["data-" schemes{i} ".hdf5"], "result");
end