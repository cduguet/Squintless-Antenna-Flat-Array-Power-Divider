function [x, resnorm] = fit_nl_ex(model, x0, xdata, ydata, options)
% Uses fmincon or fminsearch to perform MLE parameter estimation with the
% additonal possibility to penalize the likelihood.
%
% This function uses the power of Matlab functions fmincon or fminsearch to
% minimize a constructed negative log. likelihood function which calculates
% for a given noise model (Gaussian/Poisson) the negative log. likelihood
% of parameters x0 and a model function to give a measurement ydata. It is
% therefore equivalent to fit parameters x0 of a function to data, but
% unlike lsqcurvefit() it is not restricted to least square problems.
% Additionally the possibility to use built-in functions (see fitfunc.m)
% and additional penalties which depend on the parameters in the
% likelihood. With this a-priori information about distributions of parameters
% can be incorporated.
%
% Syntax
%   function [x, resnorm] = fit_nl_ex(model, x0, xdata, ydata, options)
%
% Input parameters
%   model   Defines the function which is used for the fit. This is either
%           a keyword for built-in functions (see fit_func.m) or a function
%           pointer to a function equivalent to fun used by fmincon/fminsearch,
%           i.d. function F = myfun(x, xdata)
%   x0      The initial guess for the model parameter to be estimated.
%           For the nubmer of meaning of parameter for the built-in
%           functions, see fit_func.m
%   xdata   Input data for the model function/indepedent variables.
%           For the xdata which is required by the built-in functions see
%           fit_func.m
%   ydata   Output data to be matched by the model function
%
%   Optional parameter
%   options A struct with optional fields
%       fixed   Array of 0 or 1 determining which parameters should stay fixed
%               in the fitting process.
%               default: set to zeros(size(x0) - no fixed parameter
%               Note: Even if some parameter are set fixed the model function
%               still should take all parameter (including the fixed).
%       lb      Lower bound for x
%               default: set to -Infinity
%       ub      Upper bound for x
%               default: set to +Infinity
%       opt     the options which can be forwarded to fmincon/fminsearch
%               default: all warnings and printing will be suppressed
%       use_fmincon     ={'yes'(default),'no'} determines which function
%                       from Matlab to use
%       penalty_fun     penalty function pointer for a function of type
%                       L = fun(x) which computes an additional factor for
%                       the neg.log.likelihood depending on the actual x
%                       default: not set
%       likelihood      ={'Gaussian'(default), 'Poisson'} determines which
%                       probability noise model is used to compute the
%                       likelihood
%       A, b, Aeq, beq, nonlcon
%           Like in the call to fmincon, are ignored if fminsearch is used
%           instead.
%
% Output parameters
%   x       Resulting paramters that minimize the distance between model
%           and ydata
%   resnorm Quadratic norm of (model(x, xdata) - ydata)
%
% Comment
%   When using fmincon this function can do everything that fit_nl() can do
%   and extends this even a little bit.

% parameter checking and setting to default values ------------------------

% need at least 4 parameter
if nargin < 4
    error('Not enough arguments given!');
end

x0    = double(x0); % fmincon/fminsearch need double input
ydata = double(ydata);

s = warning('query', 'all'); % save state of warnings

if nargin < 5
    options = [];
end

% developer-comment: we cannot use
% ~isfield(options, 'fixed') | isempty(options.fixed) because Matlab will
% always evaluate also the second condition even if the first already gave
% true and therefore we need more code which we pack into a function
if isnotfieldorisempty(options, 'fixed')
    options.fixed = zeros(size(x0));
end
if isnotfieldorisempty(options, 'lb')
    options.lb = zeros(size(x0)) - Inf;
end
if isnotfieldorisempty(options, 'ub')
    options.ub = zeros(size(x0)) + Inf;
end
if isnotfieldorisempty(options, 'opt')
    options.opt = optimset('Display','off');
    warning('off', 'all');
end
if isnotfieldorisempty(options, 'use_fmincon')
    options.use_fmincon = 'yes';
end
if isnotfieldorisempty(options, 'penalty_fun')
    options.penalty_fun = 0;
    % a check with isa(.., 'function_handle') will then result in false
end
if isnotfieldorisempty(options, 'likelihood')
    options.likelihood = 'Gaussian';
end
% for options.A/b/Aeq/beq/nonlcon the default is []
if ~isfield(options, 'A')
    options.A = [];
end
if ~isfield(options, 'b')
    options.b = [];
end
if ~isfield(options, 'Aeq')
    options.Aeq = [];
end
if ~isfield(options, 'beq')
    options.beq = [];
end
if ~isfield(options, 'nonlcon')
    options.nonlcon = [];
end

if size(options.fixed) ~= size(x0) | size(options.lb) ~= size(x0) | size(options.ub) ~= size(x0)
    error('Parameter x0 and options.fixed/lb/ub have not same size!');
end

% we do not check for integers, although Poissonian noise only gives
% integers, but sometimes some devices deliver significant decimal places
% and it does not hurt the algorithm, negative values however will surely
% result in havoc
if strcmp(options.likelihood, 'Poisson')
    if any(ydata < 0)
        error('Likelihood type is Poisson but image data is negative!');
    end
end

% assemble 'global' struct accessible from subfunction F

p = [];
p.fix = options.fixed;              % fixed parameters
p.x   = x0;                         % initial parameters (needed for knowing which are fixed)
p.fun = model;                      % the model function (either keyword or function pointer)
p.lik = options.likelihood;         % the likelihood type
p.ydata = ydata;                    % the measured data
p.xdata = xdata;                    % the independent variable
p.pen_fun = options.penalty_fun;    % penalty function pointer

% reduce x0 to parameters not fixed (also lb and ub)
x0 = x0(options.fixed == 0);
options.lb = options.lb(options.fixed == 0);
options.ub = options.ub(options.fixed == 0);

% the call for the internal Matlab function
switch options.use_fmincon
    case 'yes'
        % all options are given
        [x, resnorm] = fmincon(@L, x0, options.A, options.b, options.Aeq, options.beq, options.lb, options.ub, options.nonlcon, options.opt);
    case 'no'
        % with fminsearch, some options must be ignored
        [x, resnorm] = fminsearch(@L, x0, options.opt);
    otherwise
        error('Unknown content of variable options.use_fmincon!');
end

% mix fixed parameters in again
p.x(options.fixed == 0) = x;
x = p.x;

% restore warning settings
warning(s);

% functions ---------------------------------------------------------------

    function f = L(x)
        % Computes the neg. log. likelihood for given parameters x where
        % fixed parameters can be included. The model function, the noise
        % likelihood model and the likelihood penalty are specified in
        % variable p. Built-in functions can be used (see fit_func.m)
       
        % computes the likelihood of an n-D function with fixed and variable
        % parameters and a dataset (both given in the global variable p)
        
        % first mix fixed and variable parameters
        h = x;
        x = p.x;
        x(p.fix == 0) = h;

        % calculate the model function (external/built-in)
        y = fit_func(p.fun, x, p.xdata);        
        
        % calculate the likelihood between model and data and add penalty
        f = Likelihood(p.lik, p.ydata, y);
        if isa(p.pen_fun, 'function_handle')
            f = f + p.pen_fun(x);
        end
    end

    function L = Likelihood(which, img, mod)
        % Computes neg. log. likelihood between image (img) and model (mod)
        % for the given noise model (which)
        switch which
            case 'Gaussian'
                % Gaussian negative logarithmic likelihood is proportional
                % to sum((img-mod).^2)
                L = (img - mod).^2;
            case 'Poisson'
                % Poisson negative logarithmic likelihood is defined as:
                % L = sum(mod - img * ln(mod) + ln(img!)) with the Poisson
                % probability P(v, n) = exp(-v)*v^n/n! of obtaining n(>=0)
                % events if the mean number of events is v
                % img must be a nonnegative integer array and mod must be
                % strictly positive
                % we set mod to greater than 1E-6 to avoid zeros.
                mod = max(1E-6, mod);
                L = mod - img .* log(mod);
                % however the additional summand ln(img!) is not dependent
                % on the parameters and therefore does not need to be
                % included in the minimization
            otherwise
                error('Likelihood type unsupported!');
        end
        L = sum(L(:));
        % finally sum up
    end

    function b = isnotfieldorisempty(s, fieldname)
        % Gives true if fieldname is not a field of struct s or
        % if s.fieldname is empty

        b = 1;
        if isfield(s, fieldname)
            if ~isempty(getfield(s, fieldname))
                b = 0;
            end
        end
    end

end
