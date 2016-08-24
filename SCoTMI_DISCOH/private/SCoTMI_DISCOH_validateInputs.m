% Ben Cassidy (UNSW, NeuRA) Dec 2014

function [Y_o, opt_o]= SCoTMI_DISCOH_validateInputs(Y, varargin)
%%
narginchk(1, 2)
if nargin ==2, 
    opt = varargin{1};
    if  ~isstruct(opt), opt = struct; end
else opt = struct;
end

%% Need to do 2 stages since matlab can't handle extra conditions on struct inputs
% (Ridiculous!) 
%% Stage 1 - make sure all struct fields exist

YParser = inputParser;
p = inputParser; % separate because opt is a struct

%                % Defaults %
DEF.DO_PARALLEL =0;         % # computing cores, if parallel toolbox available
DEF.EBIC = 0.5;             % tuning parameter, 0 < EBIC < 1
DEF.vb = 0;                 % verbose
DEF.ar_o_max = 5;           % maximum VAR model order to search
DEF.NFFT = 2^nextpow2(1e3); % number of FFT points
DEF.cvg_tol = 1e-3;         % convergence tolerance
DEF.maxIter = 1e4;          % maximum number of gL0LS-SC iterations
DEF.h_numIter = 50;         % gL0LS-SC tuning parameter search resolution
DEF.h_initSearchScale = 10; % gL0LS-SC tuning parameter search scale
DEF.ignoreSelfNodePen = 1;  % do not penalise network self-connections

addRequired(YParser,'Y',@isnumeric);
addOptional(p,'DO_PARALLEL', DEF.DO_PARALLEL);
addOptional(p,'EBIC', DEF.EBIC);
addOptional(p,'vb', DEF.vb);
addOptional(p,'ar_o_max', DEF.ar_o_max);
addOptional(p,'NFFT', DEF.NFFT);
addOptional(p,'cvg_tol', DEF.cvg_tol);
addOptional(p,'maxIter', DEF.maxIter);
addOptional(p,'h_numIter', DEF.h_numIter);
addOptional(p,'h_initSearchScale', DEF.h_initSearchScale);
addOptional(p,'ignoreSelfNodePen', DEF.ignoreSelfNodePen);
% optParser.St`ructExpand = false;
parse(p, opt);
parse(YParser, Y);

Y_o = YParser.Results.Y;
opt_o_buf = p.Results;
%% Stage 2 - check validity for struct fields
optParser2 = inputParser;
DO_PARALLEL = opt_o_buf.DO_PARALLEL;
EBIC = opt_o_buf.EBIC;
vb = opt_o_buf.vb;
ar_o_max = opt_o_buf.ar_o_max;
NFFT = opt_o_buf.NFFT;
cvg_tol = opt_o_buf.cvg_tol;
maxIter = opt_o_buf.maxIter;
h_numIter = opt_o_buf.h_numIter;
h_initSearchScale = opt_o_buf.h_initSearchScale;
ignoreSelfNodePen = opt_o_buf.ignoreSelfNodePen;

isPosInteger = @(x) ...
    ~isempty(x) ...
    && isnumeric(x) ...
    && isreal(x) ...
    && isfinite(x) ...
    && (x == fix(x)) ...
    && x > 0;

isPosNumber = @(x) ...
    ~isempty(x) ...
    && isnumeric(x) ...
    && isreal(x) ...
    && isfinite(x) ...
    && x > 0;

isGTzeroInteger = @(x) ...
    ~isempty(x) ...
    && isnumeric(x) ...
    && isreal(x) ...
    && isfinite(x) ...
    && (x == fix(x));
        
is_0to1 = @(x)(0<=x && x<=1) ;

addOptional(optParser2, 'DO_PARALLEL', DEF.DO_PARALLEL, isGTzeroInteger)
addOptional(optParser2, 'EBIC', DEF.EBIC, is_0to1)
addOptional(optParser2, 'vb', DEF.vb, @(x)(x==true || x==false || x==1 || x==0))
addOptional(optParser2, 'ar_o_max', DEF.ar_o_max, isPosInteger)
addOptional(optParser2, 'NFFT', DEF.NFFT, isPosInteger)
addOptional(optParser2, 'cvg_tol', DEF.cvg_tol, isPosNumber)
addOptional(optParser2, 'maxIter', DEF.maxIter, isPosInteger)
addOptional(optParser2, 'h_numIter', DEF.h_numIter, isPosInteger)
addOptional(optParser2, 'h_initSearchScale', DEF.h_initSearchScale, ...
    isPosInteger)
addOptional(optParser2, 'ignoreSelfNodePen', DEF.ignoreSelfNodePen, ...
    @(x)(x==true || x==false))

parse(optParser2, DO_PARALLEL, EBIC, vb, ar_o_max, NFFT, cvg_tol, ...
    maxIter, h_numIter, h_initSearchScale, ignoreSelfNodePen)
opt_o = optParser2.Results;

end

