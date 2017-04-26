function [ u] = fschelter2009_vap1( nPoints, nDiscard)

% Schelter, Winterhalder, Eichler, Peifer,Hellwig, Guschlbauer, L?cking,
% Dahlhaus & Timmer. Testing for directed influences among neural signals 
% using partial directed coherence. J. Neurosci Methods 152:210-9, 2005.
%         [http://dx.doi.org/10.1016/j.jneumeth.2005.09.001]
% 
% Example Eq. (5) Five-dimensional VAR[4]-process 
%

disp('======================================================================');
disp('               Schelter et al. J Neurosci Methods (2009).')
disp('                    Five-dimensional VAR[3]-process')
disp('           x1-->x3  x1-->x4 x2-->x1 x2-->x3 x4==>x5 x5-->x4 ');
disp('======================================================================');


if (nargin == 0),
   nPoints = 1000;
   nDiscard = 5000;
   disp(['Adopting default ' int2str(nDiscard) ' discarding points, and'])
   disp(['generating ' int2str(nPoints) ' simulation data point.'])
   
elseif   (nargin < 2),
   nDiscard = 5000;
   disp(['Adopting default ' int2str(nDiscard) ' discarding points.'])
end;

if (nDiscard < 1),
   nDiscard = 5000;
   disp(['Adopting default ' int2str(nDiscard) ' discarding points.'])
end;

if nPoints < 10,
   nPoints = 100;
   disp(['Adopting default ' int2str(nPoints) ' simulation data points.'])
end;


N = nDiscard + nPoints; % number of simulated points

randn('state', sum(100*clock))

randn('state', sum(100*clock))
ei=randn(5,N);
x1=zeros(1,N);
x2=zeros(1,N);
x3=zeros(1,N);
x4=zeros(1,N);
x5=zeros(1,N);
% Variables initialization
for t=1:4,
   x1(t)=randn(1); x2(t)=randn(1); x3(t)=randn(1); x4(t)=randn(1);
   x5(t)=randn(1);
end;

for t=5:N,
   x1(t) = 0.9*x1(t-1) + 0.3*x2(t-2) + ei(1,t);
   x2(t) = 1.3*x2(t-1) - 0.8*x2(t-2) + ei(2,t);
   x3(t) = 0.3*x1(t-2) + 0.6*x2(t-1) + ei(3,t);
   x4(t) = -.7*x4(t-3) - 0.7*x1(t-3) + 0.3*x5(t-3) + ei(4,t);
   x5(t) = 1.0*x5(t-1) - 0.4*x5(t-2) + 0.3*x4(t-2) + ei(5,t);
end;


y=[x1' x2' x3' x4' x5']; % data must be organized column-wise
u=y(nDiscard+1:N,:);
end

