%% BACCALA & SAMESHIMA 2001A, EXAMPLE 3, PAGE 468
% DESCRIPTION:
%
% Example *baccala2001_ex3.m*:
%
% Linear five-dimensional  VAR(3) model with feedback
%
% Baccala & Sameshima. Partial directed coherence: a new concept in neural
% structure determination. _Biol. Cybern._ *84*:463-474, 2001.
%
% <http://dx.doi.org/10.1007/PL00007990>
%                
% Example 3 (pag.468)
%                       VAR(3) with feedback between x4 and x5
%
%     x1-->x2  x1-->x3 x1-->x4 x4-->x5 x5-->x4


%% Data sample generation

clear; clc; format compact; format short

nDiscard = 5000;    % number of points discarded at beginning of simulation
nPoints  = 1000;   % number of analyzed samples points

u = fbaccala2001a_ex3(nPoints, nDiscard);

chLabels = []; % or = {'x_1';'x_2';'x_3';'x_4';'x_5'};

fs = 1;

%% Interaction diagram
%
% <<fig_baccala2001a_ex3_graph.png>>
%
% Figure 2a from Baccala & Sameshima. _Biol. Cybern._ *84*:463-474, 2001.

%% Equation
%
% <<fig_baccala2001a_ex3_eq.png>>
%

%%
% Data pre-processing: detrending and normalization options

flgDetrend = 1;     % Detrending the data set
flgStandardize = 0; % No standardization
[nChannels,nSegLength] =size(u);
if nChannels > nSegLength, 
   u = u.'; 
   [nChannels,nSegLength]=size(u);
end;
if flgDetrend,
   for i=1:nChannels, u(i,:)=detrend(u(i,:)); end;
   disp('Time series were detrended.');
end;
if flgStandardize,
   for i=1:nChannels, u(i,:)=u(i,:)/std(u(i,:)); end;
   disp('Time series were scale-standardized.');
end;

%% 
% MVAR model estimation


maxIP = 30;         % maximum model order to consider.
alg = 1;            % 1: Nutall-Strand MVAR estimation algorithm
criterion = 1;      % 1: AIC, Akaike Information Criteria

disp('Running MVAR estimation routine.')

[IP,pf,A,pb,B,ef,eb,vaic,Vaicv] = mvar(u,maxIP,alg,criterion);

disp(['Number of channels = ' int2str(nChannels) ' with ' ...
    int2str(nSegLength) ' data points; MAR model order = ' int2str(IP) '.']);

%%
% Testing for adequacy of MAR model fitting through Portmanteau test

h = 20; % testing lag
MVARadequacy_signif = 0.05; % VAR model estimation adequacy significance
                            % level
aValueMVAR = 1 - MVARadequacy_signif; % Confidence value for the testing

flgPrintResults = 1;

[Pass,Portmanteau,st,ths] = mvarresidue(ef,nSegLength,IP,aValueMVAR,h,...
                                           flgPrintResults);

%%
% Granger causality test (GCT) and instantaneous GCT

gct_signif  = 0.01;  % Granger causality test significance level
igct_signif = 0.01;  % Instantaneous GCT significance level

metric = 'diag'; % euc  = original PDC or DTF;
                 % diag = generalized PDC (gPDC) or DC;
                 % info = information PDC (iPDC) or iDTF.
flgPrintResults = 1;

[Tr_gct, pValue_gct, Tr_igct, pValue_igct] = gct_alg(u,A,pf,gct_signif, ...
                                              igct_signif,flgPrintResults);
                                                       
%% Original PDC estimation
%
% PDC analysis results are saved in *c* structure.
% See asymp_dtf.m or issue 
%   >> help asymp_pdc 
% command for more detail.

nFreqs = 128;
metric = 'euc';
alpha = 0.01;
c = asymp_pdc(u,A,pf,nFreqs,metric,alpha); % Estimate PDC and asymptotic statistics

%%
% $|PDC(\lambda)|^2 Matrix Layout Plotting

flgPrinting = [1 1 1 3 0 0 3]; % overriding default setting
flgColor = 0;
w_max=fs/2;

strTitle1 = ['5-dimensional linear VAR[3] Model: '];
strTitle2 = ['[N=' int2str(nSegLength) 'pts; IP=' int2str(c.p) '; ' ...
   datestr(now) ']'];
strTitle =[strTitle1 strTitle2];
h=figure;
set(h,'NumberTitle','off','MenuBar','none', ...
   'Name', 'Baccala & Sameshima (2001) Example 3')
[hxlabel hylabel] = xplot(c,...
                          flgPrinting,fs,w_max,chLabels,flgColor);
xplot_title(alpha,metric);
[ax,hT]=suplabel( strTitle, 't' );
set(hT,'FontSize',10)


%% information PDC estimation
%
% PDC analysis results are saved in *c* structure.
% See asymp_dtf.m or issue 
%   >> help asymp_pdc 
% command for more detail.
nFreqs = 128;
metric = 'info';
alpha = 0.01;
d = asymp_pdc(u,A,pf,nFreqs,metric,alpha); % Estimate PDC and asymptotic statistics

%%
% $|_iPDC(\lambda)|^2 Matrix Layout Plotting

flgPrinting = [1 1 1 3 0 0 3]; % overriding default setting
flgColor = 0;
w_max=fs/2;

strTitle1 = ['5-dimensional linear VAR[3] Model: '];
strTitle2 = ['[N=' int2str(nSegLength) 'pts; IP=' int2str(c.p) '; ' ...
   datestr(now) ']'];
strTitle =[strTitle1 strTitle2];
h=figure;
set(h,'NumberTitle','off','MenuBar','none', ...
   'Name', 'Baccala & Sameshima (2001) Example 3')
[hxlabel hylabel] = xplot(d,...
                          flgPrinting,fs,w_max,chLabels,flgColor);
xplot_title(alpha,metric);
[ax,hT]=suplabel( strTitle, 't' );
set(hT,'FontSize',10)

%% Original DTF estimation
%
% PDC analysis results are saved in *d* structure.
% See asymp_dtf.m or issue 
%   >> help asymp_dtf 
% on command line for more detail.
metric = 'info';
d = asymp_dtf(u,A,pf,nFreqs,metric,alpha);


%%
% $|DTF(\lambda)|^2$ Matrix Layout Plotting

flgColor = 0;
flgScale = 1;
flgMax = 'dtf';
flgSignifColor = 3;

strTitle1 = ['5-dimensional linear VAR[3] Model: '];
strTitle2 = ['[N=' int2str(nSegLength) 'pts; IP=' int2str(c.p) '; ' ...
   datestr(now) ']'];
strTitle =[strTitle1 strTitle2];

h=figure;
set(h,'NumberTitle','off','MenuBar','none', ...
   'Name', 'Baccala & Sameshima (2001) Example 3')

[hxlabel,hylabel] = xplot(d,flgPrinting,fs,w_max,chLabels, ...
                                 flgColor,flgScale,flgMax,flgSignifColor);
xplot_title(alpha,metric, 'dtf');
[ax,hT]=suplabel( strTitle, 't' );
set(hT,'FontSize',10)


%% Result from the original article, Baccala & Sameshima (2001) 
% Figure 2, page 469 from article.
%
% <<fig_baccala2001a_ex3abc.png>>
% 

%% Some remarks:
% 
% # Check & compare Fig. 2b, page 468, Baccala & Sameshima (2001).
% # Note that in the original article the amplitude PDC has been plotted.
%   Here we preferred to graph squared-PDC.

%%
% This completes the Example 3 (Baccala & Sameshima, 2001)'
