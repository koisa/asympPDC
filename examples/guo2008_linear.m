%% GUO ET AL. (2008) 5-dimensional VAR[3] with large common exogenous input
%
% Guo, Wu, Ding & Feng. Uncovering interactions in frequency domains.
%         PLoS Computational Biology, 4(5):1-10, February 8, 2008.
%           [http://dx.plos.org/10.1371/journal.pcbi.1000087] 
% 
% Page 2 Toy Model Example 5-dimensional VAR[3] with large common exogenous input

%% Data sample generation
%
clear; clc; format compact; format short

nDiscard = 5000;    % number of points discarded at beginning of simulation
nPoints  = 5000;   % number of analyzed samples points

u = fguo2008_linear(nPoints, nDiscard);

chLabels = []; % or = {'x_1';'x_2';'x_3';'x_4';'x_5'};
fs = 1;

%% Interaction diagram
%
%
% <<fig_guo2008_graph.png>>
%

%% Equation system  
%
% <<fig_guo2008_eq.png>>
%

%%
% Data pre-processing: detrending and standardization options

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

flgPrintResults = 1;

[Tr_gct, pValue_gct, Tr_igct, pValue_igct] = gct_alg(u,A,pf,gct_signif, ...
                                              igct_signif,flgPrintResults);

%% Original PDC estimation
%
% PDC analysis results are saved in *c* structure.
% See asymp_dtf.m or issue 
%   >> help asymp_pdc 
% command for more detail.
%
nFreqs = 128;
metric = 'euc'; % euc  = original PDC or DTF;
                 % diag = generalized PDC (gPDC) or DC;
                 % info = information PDC (iPDC) or iDTF.
alpha = 0.01;
c = asymp_pdc(u,A,pf,nFreqs,metric,alpha); % Estimate PDC and asymptotic statistics

%%
% $|PDC(\lambda)|^2 Matrix Layout Plotting

flgPrinting = [1 1 1 3 0 0 1]; % overriding default setting
flgColor = 1;
w_max = fs/2;

for kflgColor = flgColor,
   h=figure;
   set(h,'NumberTitle','off','MenuBar','none', ...
      'Name', 'Guo et al.(2008) Linear model with large common exogenous inputs')
   [hxlabel hylabel] = xplot(c,...
      flgPrinting,fs,w_max,chLabels,kflgColor);
   xplot_title(alpha,metric,'pdc');;
end;

%% Information PDC estimation
%
% PDC analysis results are saved in *c* structure.
% See asymp_dtf.m or issue 
%   >> help asymp_pdc 
% command for more detail.
%

metric = 'info';
c = asymp_pdc(u,A,pf,nFreqs,metric,alpha); % Estimate PDC and asymptotic statistics

%%
% iPDC Matrix Layout Plotting=======================

flgScale = 2;
flgMax= 'TCI';
flgSignifColor = 2;
for kflgColor = flgColor,
   h=figure;
   set(h,'NumberTitle','off','MenuBar','none', ...
      'Name', 'Guo et al.(2008) Linear model with large common exogenous inputs')

   [hxlabel,hylabel] = xplot(c,flgPrinting,fs,w_max,chLabels, ...
                                    kflgColor,flgScale,flgMax,flgSignifColor);   
   xplot_title(alpha,metric,'pdc');
end;

%% Result from Figure 1 (bottom) Guo et al.(2006) 
% Figure 1, page 2.
%
% <<fig_guo2008_pdc_wrong_results.png>>
% 
% Guo et al. (2006) argued that PDC was not able to correctly estimate the 
% right connectivity pattern as their partial Granger causality (PGC) was (see top half of figure above).
% 
% However the PDC estimates shown on bottom of Figure 1 are completely wrong,
% as easily one can perceive by simple one-to-one comparison. It is easy to see that the PDC estimates are very similar to PGC results. 
%

%% Original DTF estimation
%
% DTF analysis results are saved in *e* structure.
% See asymp_dtf.m or issue 
%   >> help asymp_pdc 
% command for more detail.
%

nFreqs = 128;
metric = 'euc';
alpha = 0.01;
e = asymp_dtf(u,A,pf,nFreqs,metric,alpha); % Estimate DTF and asymptotic statistics

%%
% DTF Matrix Layout Plotting
%
w_max=fs/2;
for kflgColor = flgColor,
   h=figure;
   set(h,'NumberTitle','off','MenuBar','none', ...
      'Name', 'Guo et al.(2008) Linear model with large common exogenous inputs')
   [hxlabel hylabel] = xplot(e,...
      flgPrinting,fs,w_max,chLabels,kflgColor);
   xplot_title(alpha,metric,'dtf');;
end;

%% Concluding remarkds;
%
% * Guo et al. (2008)''s five-dimensional VAR[3] process with large 
% common exogenous inputs, a modified version of a example  
% from Baccala & Sameshima (2001). The exogenous inputs 
% introduce large common variance that overpower the magnitude of 
% directed interactions. As you may have noticed, there are 
% significant instantaneous Granger causality between all pair of 
% variables. 
% Due to large common exogenous white noise probably one may see  
% false-positive and false-negative connectivity in some  
% simulations, which will depend on your choice of alpha. 
% 
% * The PDC and gPDC figures do not resemble PDC plot in Fig.1 
% in Guo et al. (2008). Our best guess is that Guo & colleagues 
% used an incorrect PDC estimator. Actually you may notice 
% that PDC and gPDC estimates are very similar to PGC  
% in Guo et al (2008). 
%
% * In all three PDC formulations, the significant PDC frequency   
% range is similar. iPDC gives a measure of size effect, which is 
% very small in this case. Even so iPDC is significant in the same 
% frequency range as PDC and gPDC. 

% * In conclusion: the statement by Guo and collaborators that PDC can 
% not uncover the connectivity pattern in large common noise does 
% not hold. For nonlinear systems, in some cases PDC and other linear 
% methods can uncover correct connectivity pattern, but in other 
% models PDC and GCT will simply fail.  
%
%
%This completes Guo et al. (2008) example. ;
