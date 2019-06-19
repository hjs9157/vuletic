function [p,se,S,mu] = wpolyfit(x,y,W,k)
%%% reconstruction of MATLAB (R2018b) built-in function polyfit to:
%%% 1) take weight into account,
%%% 2) specify polynomial orders: {k} for y = \sum_{k} \beta_k x^k, and
%%% 3) output residual standard error of polynomial coefficient
%%% Last modified: 2019/04/30

%POLYFIT Fit polynomial to data.
%   P = POLYFIT(X,Y,N) finds the coefficients of a polynomial P(X) of
%   degree N that fits the data Y best in a least-squares sense. P is a
%   row vector of length N+1 containing the polynomial coefficients in
%   descending powers, P(1)*X^N + P(2)*X^(N-1) +...+ P(N)*X + P(N+1).
%
%   [P,S] = POLYFIT(X,Y,N) returns the polynomial coefficients P and a
%   structure S for use with POLYVAL to obtain error estimates for
%   predictions.  S contains fields for the triangular factor (R) from a QR
%   decomposition of the Vandermonde matrix of X, the degrees of freedom
%   (df), and the norm of the residuals (normr).  If the data Y are random,
%   an estimate of the covariance matrix of P is (Rinv*Rinv')*normr^2/df,
%   where Rinv is the inverse of R.
%
%   [P,S,MU] = POLYFIT(X,Y,N) finds the coefficients of a polynomial in
%   XHAT = (X-MU(1))/MU(2) where MU(1) = MEAN(X) and MU(2) = STD(X). This
%   centering and scaling transformation improves the numerical properties
%   of both the polynomial and the fitting algorithm.
%
%   Warning messages result if N is >= length(X), if X has repeated, or
%   nearly repeated, points, or if X might need centering and scaling.
%
%   Example: simple linear regression with polyfit
%
%     % Fit a polynomial p of degree 1 to the (x,y) data:
%       x = 1:50;
%       y = -0.3*x + 2*randn(1,50);
%       p = polyfit(x,y,1);
%
%     % Evaluate the fitted polynomial p and plot:
%       f = polyval(p,x);
%       plot(x,y,'o',x,f,'-')
%       legend('data','linear fit')
%
%   Class support for inputs X,Y:
%      float: double, single
%
%   See also POLY, POLYVAL, ROOTS, LSCOV.

%   Copyright 1984-2017 The MathWorks, Inc.

if ~isequal(size(x),size(y),size(W))
    error(message('MATLAB:polyfit:XYSizeMismatch'))
end

x = x(:);
y = y(:);
W = W(:);

if any(W < 0)
    error('w: Negative weights are not allowed.')
end

if any(k < 0)
    error('k: Negative orders are not allowed.')
end

% if length(k) > length(x)% if length(k) > length(x)
%     error('k: the number of parameters must not be bigger than the number of data.')
% end

%     error('k: the number of parameters must not be bigger than the number of data.')
% end

k = k(:).';

if nargout > 3
    mu = [mean(x); std(x)];
    x = (x - mu(1))/mu(2);
end

% get sqrt weight and scale it
w = sqrt(W);
w = w./mean(w);

% Construct the Vandermonde matrix V = [x.^k(1) ... x.^k(end)]
V = x.^k;

% weight variables
wy = w.*y;
wV = V;
for i = 1:length(y)
    wV(i,:) = w(i)*wV(i,:);
end

% Solve least squares problem p = V\y to get polynomial coefficients p.
[wQ,wR] = qr(wV,0);
oldws = warning('off','all');   % Turn all warnings off before solving
try
    p = wR\(wQ.'*wy);               % Same as p = wV\wy
catch ME
    warning(oldws);             % Restore initial warning state
    throw(ME);
end
warning(oldws);                 % Restore initial warning state

% Issue warnings.
PolyNotUniqueFlag = false;
if size(wR,2) > size(wR,1)
    PolyNotUniqueFlag = true;
    warning(message('MATLAB:polyfit:PolyNotUnique'))
elseif warnIfLargeConditionNumber(wR)
    if nargout > 2
        warning(message('MATLAB:polyfit:RepeatedPoints'));
    else
        warning(message('MATLAB:polyfit:RepeatedPointsOrRescale'));
    end
end

% get standard error se of polynomial coefficients p.
if PolyNotUniqueFlag
    se = Inf(size(k));
else
    wR1 = wR(1:size(wR,2),:);
    IwR1 = inv(wR1);
    wQ = IwR1*IwR1.'; % variance-covariance matrix
    
    wr = wy - wV*p;
    df = max(0,length(y) - length(k));
    ssq = sum(wr.^2)/df; % s^2; residual standard estimate of variance of weighted errors
    
    se = sqrt(ssq*diag(wQ).');  
    
end

% if nargout > 2
    r = y - V*p; % sqrt-weighted residual
    % S is a structure containing four elements: the triangular factor
    % from a QR decomposition of the Vandermonde matrix, the degrees of
    % freedom, the norm of the residuals, and average of weights.
    S.R = wR;
    S.df = max(0,length(y) - length(k));
    S.normr = norm(r);
    S.meanW = mean(w.^2);
% end

p = p.'; % Polynomial coefficients are row vectors by convention.

function flag = warnIfLargeConditionNumber(R)
if isa(R, 'double')
    flag = (condest(R) > 1e+10);
else
    flag = (condest(R) > 1e+05);
end
