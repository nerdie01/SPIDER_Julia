clear;
pkg load symbolic;

addpath('mhd')
addpath('navier_stokes')


n = 128;
dt = 0.01;
timesteps = 1000;
nu = 0.002;
eta = 0.01;
k1 = 1.0;
k2 = 1.0
skip = 10;

# random noise
rng = 42
omega = @(x, y) rand(n,n) #1/k1 * cos(k1*x) + 1/k2 * cos(k2*y)

rng = 21
A = @(x, y) rand(n,n);
forcing = @(x, y) 0;

butcher_rk2 = struct();
butcher_rk2.A = [1/2] # implicit midpoint / gauss-legendre order 2
butcher_rk2.b = [1]

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

schemes={'RK2', 'RK3', 'RK4', 'RK5', 'RK6'};

butchers={butcher_rk6};

for i = 1:length(schemes)
    [x, y, t, u, v, p, Bx, By] = mhd_rk_implicit(omega, A, forcing, n, dt, timesteps, nu, eta, 1, skip, butchers{i}, true);
    #result = struct("x",x,"y",y,"t",t,"u",u,"v",v,"p",p,"Bx",Bx,"By",By);
    #save("-hdf5", ["orszag-tang-" schemes{i} ".hdf5"], "result");
end