% Calculate SCoTMI and DISCOH measures (un-thresholded) for estimating the
% topology of time series networks. 
%
% Add the directory for this file to the matlab path.
% 
% 
% [SCoTMI_raw,DISCOH_raw, NWvarTopo] = SCoTMI_DISCOH_calc(Y, varargin)
%
% Inputs: 
%
% Y   : 2D matrix (double), 
%   size== T x N ( time samples x # network nodes )
% opt : (optional) struct
%   fields: (default values shown) (any fields may be empty)
%       DO_PARALLEL             % # computing cores, if parallel toolbox available
%       EBIC                    % tuning parameter, 0 < EBIC < 1
%       vb                      % verbose
%       ar_o_max                % maximum VAR model order to search
%       NFFT                    % number of FFT points
%       cvg_tol                 % convergence tolerance
%       maxIter                 % maximum number of gL0LS-SC iterations
%       h_numIter               % gL0LS-SC tuning parameter search resolution
%       h_initSearchScale       % gL0LS-SC tuning parameter search scale
%       ignoreSelfNodePen       % do not penalise network self-connections
%
% Outputs:
% 
% SCoTMI_raw : lower-triangular 2D matrix (double), size == N x N
%       Network link strengths (SCoTMI is a measure of similarity)
% DISCOH_raw : lower-triangular 2D matrix (double), size == N x N
%       Network link strengths (DISCOH is a measure of dissimilarity)
% NWvarTopo : 2D matrix, size == N x N
%       Estimated network VAR model orders at each node, and location of
%       non-zero VAR parameters. This gives an indication if the EBIC
%       tuning parameter is too low, if any rows of this matrix are all
%       non-zero. Useful for sanity checking results.
%
% Example : unconnected white noise network, 100 time points, 5 nodes
%     Y=randn(100,5); 
%     opt.vb = true; 
%     opt.ar_o_max= 2;
%     [SCoTMI_raw,DISCOH_raw, NWvarTopo] = SCoTMI_DISCOH_calc(Y, opt);
%
% v0.1
% Ben Cassidy, (UNSW, NeuRA) Dec 2014



function [SCoTMI_raw,DISCOH_raw, NWvarTopo] = SCoTMI_DISCOH_calc(Y, varargin)
%% path spec.
% fdir = fileparts(which('SCoTMI_DISCOH_calc'));
% path(fullfile(fdir, 'auxFiles' ),path)

%% input checking
narginchk(1, 2)
if nargin ==2, 
    opt = varargin{1};
    if  ~isstruct(opt), opt = struct; end
else opt = struct;
end

[Y, opt]= SCoTMI_DISCOH_validateInputs(Y, opt);
opt.L1_init = 0; % bug, remove later


%% network system identification

[T,nNodes] = size(Y);
Tm = T-opt.ar_o_max;

if opt.vb, fprintf('Starting gL0LS VAR network ID....'), end



[nodeResid, varPara, BIC, NWvarTopo] ...
    = gL0LS_SC_VAR...
        (Y, opt.ar_o_max, opt.h_initSearchScale, ...
        opt.h_numIter, opt.L1_init, opt.ignoreSelfNodePen, opt.EBIC,...
        opt.maxIter, ...
        opt.cvg_tol, opt.vb, opt.DO_PARALLEL);

    if opt.vb, fprintf('Finished gL0LS VAR network ID....'), end
%% partial coherence 


Aomega = cell(nNodes,1);
for nodeLP=1:nNodes
    % important:: fft explicitly along columns.
    Aomega{nodeLP} =    -   (1./sqrt(opt.NFFT))*fft(full(varPara{nodeLP}), opt.NFFT, 1);
    Aomega{nodeLP}(:,nodeLP) = Aomega{nodeLP}(:,nodeLP) + 1; %already did minus, now add 1.
end

[partCoh, CHOLFLAG] = partCohCalc(nodeResid, Aomega, opt.NFFT);


%% SCoTMI & DISCOH calculation
SCoTMI_raw = tMIcalc(partCoh, opt.NFFT, Tm);

DISCOH_raw = DISCOHprojected_calc(partCoh, opt.NFFT);
% DISCOH_raw = DISCOH_calc(partCoh, opt.NFFT);

end

