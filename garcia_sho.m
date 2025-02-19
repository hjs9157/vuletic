%**************************************************
%**************
% Program 4: Find several lowest eigenmodes V(x) and
% eigenenergies E of 1D Schrodinger equation
%         -1/2*hbar^2/m(d2/dx2)V(x) + U(x)V(x) = EV(x)
% for arbitrary potentials U(x)
%**************************************************
%**************
% Parameters for solving problem in the interval -L < x < L
% PARAMETERS:
L = 5;                      % Interval Length
N = 10000;                   % No of points
x = linspace(-L,L,N)';      % Coordinate vector
dx = x(2) - x(1);           % Coordinate step
% POTENTIAL, choose one or make your own
U = 1/2*100*x.^(2);    % quadratic harmonic oscillator potential
%U = 1/2*x.^(4);       % quartic potential
% Finite square well of width 2w and depth given
%w = L/50;
%U = -500*(heaviside(x+w)-heaviside(x-w));
% Two finite square wells of width 2w and distance 2a apart
%w = L/50; a=3*w;
%U = -200*(heaviside(x+w-a) - heaviside(x-w-a) ...
%              + heaviside(x+w+a) - heaviside(x-w+a));
% Three-point finite-difference representation of Laplacian
% using sparse matrices, where you save memory by only
% storing non-zero matrix elements
e = ones(N,1); Lap = spdiags([e -2*e e],[-1 0 1],N,N)/dx^2;
% Total Hamiltonian
hbar = 1; m = 1;      % constants for Hamiltonian
H = -1/2*(hbar^2/m)*Lap + spdiags(U,0,N,N);
% Find lowest nmodes eigenvectors and eigenvalues of sparse matrix
nmodes = 3; options.disp = 0;options.maxit = 10000;options.tol=1e-10;
[V,E] = eigs(H,nmodes,'sa',options);   % find eigs
[E,ind] = sort(diag(E));% convert E to vector and sort low to high
V = V(:,ind);           % rearrange corresponding eigenvectors
% Generate plot of lowest energy eigenvectors V(x) and U(x)
Usc = U*max(abs(V(:)))/max(abs(U));       % rescale U for plotting
plot(x,V,x,Usc,'--k');          % plot V(x) and rescaled U(x)
% Add legend showing Energy of plotted V(x)
lgnd_str = [repmat('E = ',nmodes,1),num2str(E)];
legend(lgnd_str)                % place lengend string on plot
shg