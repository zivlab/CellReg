function [sigma, mu] = gaussfit( x, y, sigma0, mu0 )
% [sigma, mu] = gaussfit( x, y, sigma0, mu0 )
% Fits a guassian probability density function into (x,y) points using iterative 
% LMS method. Gaussian p.d.f is given by: 
% y = 1/(sqrt(2*pi)*sigma)*exp( -(x - mu)^2 / (2*sigma^2))
% The results are much better than minimazing logarithmic residuals
%
% INPUT:
% sigma0 - initial value of sigma (optional)
% mu0 - initial value of mean (optional)
%
% OUTPUT:
% sigma - optimal value of standard deviation
% mu - optimal value of mean
%
% REMARKS:
% The function does not always converge in which case try to use initial
% values sigma0, mu0. Check also if the data is properly scaled, i.e. p.d.f
% should approx. sum up to 1
% 
% VERSION: 23.02.2012
% 
% EXAMPLE USAGE:
% x = -10:1:10;
% s = 2;
% m = 3;
% y = 1/(sqrt(2*pi)* s ) * exp( - (x-m).^2 / (2*s^2)) + 0.02*randn( 1, 21 );
% [sigma,mu] = gaussfit( x, y )
% xp = -10:0.1:10;
% yp = 1/(sqrt(2*pi)* sigma ) * exp( - (xp-mu).^2 / (2*sigma^2));
% plot( x, y, 'o', xp, yp, '-' );


% Maximum number of iterations
Nmax = 50;

if( length( x ) ~= length( y ))
    fprintf( 'x and y should be of equal length\n\r' );
    exit;
end

n = length( x );
x = reshape( x, n, 1 );
y = reshape( y, n, 1 );

%sort according to x
X = [x,y];
X = sortrows( X );
x = X(:,1);
y = X(:,2);

%Checking if the data is normalized
dx = diff( x );
dy = 0.5*(y(1:length(y)-1) + y(2:length(y)));
s = sum( dx .* dy );
if( s > 1.5 | s < 0.5 )
    fprintf( 'Data is not normalized! The pdf sums to: %f. Normalizing...\n\r', s );
    y = y ./ s;
end

% try to estimate mean mu from the location of the maximum
[ymax,index]=max(y);
mu = x(index);

% estimate sigma
sigma = 1/(sqrt(2*pi)*ymax);

if( nargin == 3 )
    sigma = sigma0;
end

if( nargin == 4 )
    mu = mu0;
end

%xp = linspace( min(x), max(x) );

% iterations
for i=1:Nmax
%    yp = 1/(sqrt(2*pi)*sigma) * exp( -(xp - mu).^2 / (2*sigma^2));
%    plot( x, y, 'o', xp, yp, '-' );

    dfdsigma = -1/(sqrt(2*pi)*sigma^2)*exp(-((x-mu).^2) / (2*sigma^2));
    dfdsigma = dfdsigma + 1/(sqrt(2*pi)*sigma).*exp(-((x-mu).^2) / (2*sigma^2)).*((x-mu).^2/sigma^3);

    dfdmu = 1/(sqrt(2*pi)*sigma)*exp(-((x-mu).^2)/(2*sigma^2)).*(x-mu)/(sigma^2);

    F = [ dfdsigma dfdmu ];
    a0 = [sigma;mu];
    f0 = 1/(sqrt(2*pi)*sigma).*exp( -(x-mu).^2 /(2*sigma^2));
    a = (F'*F)^(-1)*F'*(y-f0) + a0;
    sigma = a(1);
    mu = a(2);
    
    if( sigma < 0 )
        fprintf( 'Instability detected! Rerun with initial values sigma0 and mu0! \n\r' );
        warning( 'Check if your data is properly scaled! p.d.f should approx. sum up to \n\r' );        
    end
end

