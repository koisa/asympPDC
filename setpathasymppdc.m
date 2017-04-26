% Directory under which the AsympPDC package is installed
% We used this for using Octave.
% In MATLAB set path using menu option.

vroot = '/Users/koichi/Dropbox/octave'; % base directory where your Octave packages and AsympPDC are located

% Edit this to reflect your Octave package installation.
%addpath([vroot '/signal-1.3.0'])
%addpath([vroot '/signal-1.3.0/x86_64-apple-darwin13.0.0-api-v49+'])
%addpath([vroot '/statistics-1.2.3'])
%addpath([vroot '/io-2.2.6'])

% Control package is not necessary.
%addpath([vroot '/control-2.6.6'])
%addpath([vroot '/control-2.6.6/x86_64-apple-darwin13.0.0-api-v49+'])

addpath([vroot '/asymp_package_v3'])
addpath([vroot '/asymp_package_v3/supporting'])
addpath([vroot '/asymp_package_v3/supporting/arfit'])
addpath([vroot '/asymp_package_v3/supporting/extras'])
%addpath([vroot '/asymp_package_v3/supporting/shadedErrorBar'])
addpath([vroot '/asymp_package_v3/supporting/shadedplot'])
addpath([vroot '/asymp_package_v3/supporting/suplabel'])
addpath([vroot '/asymp_package_v3/supporting/boundedline'])
addpath([vroot '/asymp_package_v3/routines'])
addpath([vroot '/asymp_package_v3/examples'])

format short
warning('off'); more off
