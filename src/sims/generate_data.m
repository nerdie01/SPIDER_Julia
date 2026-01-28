clear;
pkg load symbolic;

addpath('mhd')
addpath('navier_stokes')


n = 64;
dt = 0.02;
timesteps = 100;
nu = 2^-6;
eta = 2^-6;
skip = 10;

truncate = 1;

# random noise
rng(1)
omega = @(x, y) rand(n,n)

rng(2)
A = @(x, y) rand(n,n);
forcing = @(x, y) 0;

butcher_rk1 = struct();
butcher_rk1.A = [1];
butcher_rk1.b = [1]

butcher_rk2 = struct();
butcher_rk2.A = [1/2]; # implicit midpoint / gauss-legendre order 2
butcher_rk2.b = [1];

butcher_rk3 = struct();
butcher_rk3.A = [
    1/4 -1/4;
    1/4 5/12]; # radau order 3
butcher_rk3.b = [1/4 3/4];

butcher_rk4 = struct();
butcher_rk4.A = [
    1/4 (1/4-sqrt(3)/6);
    (1/4+sqrt(3)/6) 1/4]; # gauss-legendre order 4
butcher_rk4.b = [1/2 1/2];

butcher_rk5 = struct();
butcher_rk5.A = [
    1/9 ((-1-sqrt(6))/18) ((-1+sqrt(6))/18);
    1/9 (11/45+7*sqrt(6)/360) (11/45-43*sqrt(6)/360);
    1/9 (11/45+43*sqrt(6)/360) (11/45-7*sqrt(6)/360)];
butcher_rk5.b = [1/9 (4/9+sqrt(6)/36) (4/9-sqrt(6)/36)];

butcher_rk6 = struct();
butcher_rk6.A = [
    5/36 (2/9-sqrt(15)/15) (5/36-sqrt(15)/30);
    (5/36+sqrt(15)/24) 2/9 (5/36-sqrt(15)/24);
    (5/36+sqrt(15)/30) (2/9+sqrt(15)/15) 5/36
];
butcher_rk6.b = [5/18 4/9 5/18];

schemes={'RK1_bs'};

butchers={butcher_rk1};

for i = 1:length(schemes)
    [x, y, t, u, v, p, Bx, By] = mhd_rk_implicit(omega, A, forcing, n, dt, timesteps, nu, eta, 1, skip, butchers{i}, true);
    result = struct("x",x,"y",y,"t",t(truncate:end),"u",u(:,:,truncate:end),"v",v(:,:,truncate:end),"p",p(:,:,truncate:end),"Bx",Bx(:,:,truncate:end),"By",By(:,:,truncate:end));
    save("-hdf5", ["data-" schemes{i} ".hdf5"], "result");
end