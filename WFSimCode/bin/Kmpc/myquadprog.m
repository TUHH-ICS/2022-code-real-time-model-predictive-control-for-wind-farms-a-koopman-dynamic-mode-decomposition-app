function [X,fval] = myquadprog(H,f,A,B)
%QUADPROG Quadratic programming. 
%   X = QUADPROG(H,f,A,b) attempts to solve the quadratic programming 
%   problem:
%
%            min 0.5*x'*H*x + f'*x   subject to:  A*x <= b 
%             x    
%
%   X = QUADPROG(H,f,A,b,Aeq,beq) solves the problem above while
%   additionally satisfying the equality constraints Aeq*x = beq. (Set A=[]
%   and B=[] if no inequalities exist.)
%
%   X = QUADPROG(H,f,A,b,Aeq,beq,LB,UB) defines a set of lower and upper
%   bounds on the design variables, X, so that the solution is in the 
%   range LB <= X <= UB. Use empty matrices for LB and UB if no bounds 
%   exist. Set LB(i) = -Inf if X(i) is unbounded below; set UB(i) = Inf if 
%   X(i) is unbounded above.
%
%   X = QUADPROG(H,f,A,b,Aeq,beq,LB,UB,X0) sets the starting point to X0.
%
%   X = QUADPROG(H,f,A,b,Aeq,beq,LB,UB,X0,OPTIONS) minimizes with the
%   default optimization parameters replaced by values in OPTIONS, an
%   argument created with the OPTIMOPTIONS function. See OPTIMOPTIONS for
%   details.
%
%   X = QUADPROG(PROBLEM) finds the minimum for PROBLEM. PROBLEM is a
%   structure with matrix 'H' in PROBLEM.H, the vector 'f' in PROBLEM.f,
%   the linear inequality constraints in PROBLEM.Aineq and PROBLEM.bineq,
%   the linear equality constraints in PROBLEM.Aeq and PROBLEM.beq, the
%   lower bounds in PROBLEM.lb, the upper bounds in PROBLEM.ub, the start
%   point in PROBLEM.x0, the options structure in PROBLEM.options, and 
%   solver name 'quadprog' in PROBLEM.solver. Use this syntax to solve at 
%   the command line a problem exported from OPTIMTOOL. 
%
%   [X,FVAL] = QUADPROG(H,f,A,b) returns the value of the objective 
%   function at X: FVAL = 0.5*X'*H*X + f'*X.
%
%   [X,FVAL,EXITFLAG] = QUADPROG(H,f,A,b) returns an EXITFLAG that
%   describes the exit condition. Possible values of EXITFLAG and the
%   corresponding exit conditions are
%
%   All algorithms:
%     1  First order optimality conditions satisfied.
%     0  Maximum number of iterations exceeded.
%    -2  No feasible point found.
%    -3  Problem is unbounded.
%   Interior-point-convex only:
%     2  Solver stalled at feasible point.
%    -6  Non-convex problem detected.
%    -8  Unable to compute step direction; no further progress can be made.
%   Trust-region-reflective only:
%     3  Change in objective function too small.
%    -4  Current search direction is not a descent direction; no further 
%        progress can be made.
%
%   [X,FVAL,EXITFLAG,OUTPUT] = QUADPROG(H,f,A,b) returns a structure
%   OUTPUT with the number of iterations taken in OUTPUT.iterations,
%   maximum of constraint violations in OUTPUT.constrviolation, the 
%   type of algorithm used in OUTPUT.algorithm, the number of conjugate
%   gradient iterations (if used) in OUTPUT.cgiterations, a measure of 
%   first order optimality (large-scale algorithm only) in 
%   OUTPUT.firstorderopt, and the exit message in OUTPUT.message.
%
%   [X,FVAL,EXITFLAG,OUTPUT,LAMBDA] = QUADPROG(H,f,A,b) returns the set of 
%   Lagrangian multipliers LAMBDA, at the solution: LAMBDA.ineqlin for the 
%   linear inequalities A, LAMBDA.eqlin for the linear equalities Aeq, 
%   LAMBDA.lower for LB, and LAMBDA.upper for UB.
%
%   See also LINPROG, LSQLIN.

%   Copyright 1990-2018 The MathWorks, Inc.

defaultopt = struct( ...
    'Algorithm','interior-point-convex', ...
    'Diagnostics','off', ...
    'Display','off', ... %final
    'HessMult',[], ... 
    'MaxIter',[], ...    
    'MaxPCGIter','max(1,floor(numberOfVariables/2))', ...   
    'PrecondBandWidth',0, ... 
    'TolCon',1e-8, ...
    'TolFun',[], ...
    'TolFunValue', [], ...
    'TolPCG',0.1, ...    
    'TolX',100*eps, ...
    'TypicalX','ones(numberOfVariables,1)', ...  
    'LinearSolver', 'auto' ...
    );

% If just 'defaults' passed in, return the default options in X
% if nargin == 1 && nargout <= 1 && strcmpi(H,'defaults')
%    X = defaultopt;
%    return
% end

% Handle missing arguments
%if nargin < 10
options = [];
% if nargin < 9
X0 = [];
% if nargin < 8
ub = [];
% if nargin < 7
lb = [];
%  if nargin < 6
Beq = [];
% if nargin < 5
Aeq = [];
% if nargin < 4
%     B = [];
%     if nargin < 3
%         A = [];
%     end
% end

% Detect problem structure input
% if nargin == 1
%    if isa(H,'struct')
%        [H,f,A,B,Aeq,Beq,lb,ub,X0,options] = separateOptimStruct(H);
%    else % Single input and non-structure.
%        error(message('optim:quadprog:InputArg'));
%    end
% end

% No options passed. Set options directly to defaultopt after
allDefaultOpts = isempty(options);

% Prepare the options for the solver
options = prepareOptionsForSolver(options, 'quadprog');

% if nargin == 0 
%    error(message('optim:quadprog:NotEnoughInputs'))
% end

% Set options to default if no options were passed.
% if allDefaultOpts
    % Options are all default
    options = defaultopt;
% end

% Check for non-double inputs
% msg = isoptimargdbl('QUADPROG', {'H','f','A','b','Aeq','beq','LB','UB','X0'}, ...
%                                   H,  f,  A,  B,  Aeq,  Beq,  lb,  ub,  X0);
% if ~isempty(msg)
%     error('optim:quadprog:NonDoubleInput',msg);
% end
                     
% Set up constant strings
trustRegReflect = 'trust-region-reflective';
interiorPointConvex = 'interior-point-convex';

% if nargout > 4
%    computeLambda = true;
% else 
   computeLambda = false;
% end
% if nargout > 3
%    computeConstrViolation = true;
%    computeFirstOrderOpt = true;
% else 
   computeConstrViolation = false;
   computeFirstOrderOpt = false;
%end

% Options setup
Algorithm = optimget(options,'Algorithm',defaultopt,'fast',allDefaultOpts); 

diagnostics = strcmpi(optimget(options,'Diagnostics',defaultopt,'fast',allDefaultOpts),'on');
display = optimget(options,'Display',defaultopt,'fast',allDefaultOpts);
detailedExitMsg = contains(display,'detailed');
% switch display
% case {'off', 'none'}
   verbosity = 0;
% case {'iter','iter-detailed'}
%    verbosity = 2;
% case {'final','final-detailed'}
%    verbosity = 1;
% case 'testing'
%    verbosity = 3;
% otherwise
%    verbosity = 1;
% end

% Determine algorithm user chose via options. (We need this now to set
% OUTPUT.algorithm in case of early termination due to inconsistent
% bounds.) 
% if strcmpi(Algorithm,'active-set')
%     % Error: active-set algorithm is removed
%    [linkTag,endLinkTag] = linkToAlgDefaultChangeCsh('quadprog_warn_will_error'); % links to context sensitive help
%    error(message('optim:quadprog:ActSetRemoved','active-set','quadprog', ...
%             linkTag,endLinkTag,'interior-point-convex','trust-region-reflective'));
% elseif strcmpi(Algorithm,'interior-point-convex')
    output.algorithm = interiorPointConvex;
% elseif strcmpi(Algorithm,'trust-region-reflective')
%     output.algorithm = trustRegReflect;
% else
%     error(message('optim:quadprog:InvalidAlgorithm'));
% end 

mtxmpy = optimget(options,'HessMult',defaultopt,'fast',allDefaultOpts);
% Check for name clash
functionNameClashCheck('HessMult',mtxmpy,'hessMult_optimInternal','optim:quadprog:HessMultNameClash');
% if isempty(mtxmpy)
    % Internal Hessian-multiply function
    mtxmpy = @hessMult_optimInternal;
    usrSuppliedHessMult = false;     
% else
%     usrSuppliedHessMult = true;
% end

% Set the constraints up: defaults and check size
[nineqcstr,numberOfVariablesineq] = size(A);
[neqcstr,numberOfVariableseq] = size(Aeq);
% if isa(H,'double') && ~usrSuppliedHessMult
%    % H must be square and have the correct size 
nColsH = size(H,2);
%    if nColsH ~= size(H,1)
%       error(message('optim:quadprog:NonSquareHessian'));
%    end
% else % HessMult in effect, so H can be anything
%    nColsH = 0;
% end

% Check the number of variables. The check must account for any combination of these cases:
% * User provides HessMult
% * The problem is linear (H = zeros, or H = [])
% * The objective has no linear component (f = [])
% * There are no linear constraints (A,Aeq = [])
% * There are no, or partially specified, bounds 
% * There is no X0
numberOfVariables = ...
    max([length(f),nColsH,numberOfVariablesineq,numberOfVariableseq]);

% if numberOfVariables == 0
%     % If none of the problem quantities indicate the number of variables,
%     % check X0, even though some algorithms do not use it.
%     if isempty(X0)
%         error(message('optim:quadprog:EmptyProblem'));
%     else
%         % With all other data empty, use the X0 input to determine
%         % the number of variables.
%         numberOfVariables = length(X0);
%     end
% end

% if isempty(f)
%     f = zeros(numberOfVariables,1);
% else 
%     % Make sure that the number of rows/columns in H matches the length of
%     % f under the following conditions:
%     % * The Hessian is passed in explicitly (no HessMult)
%     % * There is a non-empty Hessian
%     if ~usrSuppliedHessMult && ~isempty(H)
%         if numel(f) ~= nColsH
%             error(message('optim:quadprog:MismatchObjCoefSize'));
%         end
%     end
% end
if isempty(H)
    H = sparse(numberOfVariables,numberOfVariables);
end
if isempty(A)
    A = zeros(0,numberOfVariables);
end
if isempty(B)
    B = zeros(0,1);
end
%if isempty(Aeq)
    Aeq = zeros(0,numberOfVariables); 
%end
%if isempty(Beq)
    Beq = zeros(0,1);
%end

% Expect vectors
f = f(:);
B = B(:);
Beq = Beq(:);

% if ~isequal(length(B),nineqcstr)
%     error(message('optim:quadprog:InvalidSizesOfAAndB'))
% elseif ~isequal(length(Beq),neqcstr)
%     error(message('optim:quadprog:InvalidSizesOfAeqAndBeq'))
% elseif ~isequal(length(f),numberOfVariablesineq) && ~isempty(A)
%     error(message('optim:quadprog:InvalidSizesOfAAndF'))
% elseif ~isequal(length(f),numberOfVariableseq) && ~isempty(Aeq)
%     error(message('optim:quadprog:InvalidSizesOfAeqAndf'))
% end

[X0,lb,ub,msg] = checkbounds(X0,lb,ub,numberOfVariables);
if ~isempty(msg)
   exitflag = -2;
   X=X0; fval = []; lambda = [];
   output.iterations = 0;
   output.constrviolation = [];
   output.firstorderopt = [];
   output.cgiterations = []; 
   output.linearsolver = [];
   output.message = msg;
   if verbosity > 0
      disp(msg)
   end
   return
end

% Check that all data is real
% if ~(isreal(H) && isreal(A) && isreal(Aeq) && isreal(f) && ...
%      isreal(B) && isreal(Beq) && isreal(lb) && isreal(ub) && isreal(X0))
%     error(message('optim:quadprog:ComplexData'))
% end

isLPProblem = false;
% Perform checks on H
if isa(H,'double') && ~usrSuppliedHessMult
   % check if H is all zeros
   if ~any(H(:))
      % Really an LP problem. Warn, but continue to solve using quadprog. 
      warning(message('optim:quadprog:NullHessian'))
      isLPProblem = true;
   else
      % Make sure Hessian matrix is symmetric
      if ~issymmetric(H)
         if verbosity > -1
            warning(message('optim:quadprog:HessianNotSym'))
         end
         H = (H+H')*0.5;
      end
   end
end

% If user passed HessMult (and no Hessian matrix) and chose an algorithm
% other than trust-region-reflectve, error out.
if ~isa(H,'double') || usrSuppliedHessMult &&  ...
        ~strcmpi(output.algorithm,trustRegReflect)
    error(message('optim:quadprog:NoHessMult', Algorithm))
end

% if diagnostics 
%    % Do diagnostics on information so far
%    gradflag = []; hessflag = []; constflag = false; gradconstflag = false; 
%    non_eq=0;non_ineq=0; lin_eq=size(Aeq,1); lin_ineq=size(A,1); 
%    XOUT=ones(numberOfVariables,1); funfcn{1} = []; confcn{1}=[];
%    diagnose('quadprog',output,gradflag,hessflag,constflag,gradconstflag,...
%       XOUT,non_eq,non_ineq,lin_eq,lin_ineq,lb,ub,funfcn,confcn);
% end

if strcmpi(output.algorithm,interiorPointConvex)
  
    defaultopt.MaxIter = 200;
    defaultopt.TolFun = 1e-8;
    defaultopt.TolX = 1e-12;
    % Set ConvexCheck to notify solver that the problem should be monitored
    % for non-convexity. However, if a user has passed an LP, then the
    % problem is convex, so no need to check.
    if isLPProblem
        defaultopt.ConvexCheck = 0;
    else
        defaultopt.ConvexCheck = 1;
    end
    
    % If the output structure is requested, we must reconstruct the
    % Lagrange multipliers in the postsolve. Therefore, set computeLambda
    % to true if the output structure is requested.
    flags.computeLambda = computeFirstOrderOpt; 
    flags.detailedExitMsg = detailedExitMsg;
    flags.verbosity = verbosity;
    flags.caller = 'quadprog';
    
    % Check which solver the user requested
    linearSolver = optimget(options,'LinearSolver',defaultopt,'fast',allDefaultOpts); 
%     autoSelect = strcmp(linearSolver, 'auto');
    
    % Full QP
    % if (autoSelect && ~issparse(H)) || strcmp(linearSolver, 'dense')
    [X, fval, exitflag, output, lambda] = ...
        ipqpdense(full(H), f, A, B, Aeq, Beq, lb, ub, X0, flags, ...
        options, defaultopt);
    output.linearsolver = 'dense';

    % Sparse QP    
%     else
%         [X, fval, exitflag, output, lambda] = ...
%             ipqpcommon(sparse(H), f, A, B, Aeq, Beq, lb, ub, X0, flags, ...
%             options, defaultopt);
%         output.linearsolver = 'sparse';
%     end
    
    % Presolve may have removed variables and constraints from the problem.
    % Postsolve will re-insert the primal and dual solutions after the main
    % algorithm has run. Therefore, constraint violation and first-order
    % optimality must be re-computed.
    %  
    % If no initial point was provided by the user and the presolve has
    % declared the problem infeasible or unbounded, X will be empty. The
    % lambda structure will also be empty, so do not compute constraint
    % violation or first-order optimality if lambda is missing.
    
    % Compute constraint violation if the output structure is requested
    if computeFirstOrderOpt && ~isempty(lambda)
        output.constrviolation = norm([Aeq*X-Beq;max([A*X - B;X - ub;lb - X],0)],Inf);        
    end
elseif strcmpi(output.algorithm,trustRegReflect)
    % Determine whether trust-region-reflective can handle problem
    hasIneqs = (nineqcstr > 0);  % Does the problem have any inequalities?
    hasEqsAndBnds = (neqcstr > 0) && (any(isfinite(ub)) || any(isfinite(lb))); % Does the problem have both equalities and bounds?
    hasMoreEqsThanVars = (neqcstr > numberOfVariables); % Does the problem have more equalities than variables?
    hasNoConstrs = (neqcstr == 0) && (nineqcstr == 0) && ...
        all(eq(ub, inf)) && all(eq(lb, -inf)); % Does the problem not have equalities, bounds, or inequalities?
    
    linkToDoc = addLink('Choosing the Algorithm', 'optim', 'helptargets.map', ...
        'choose_algorithm', false);
    if (hasIneqs || hasEqsAndBnds || hasMoreEqsThanVars || hasNoConstrs)
        error(message('optim:quadprog:ConstrTRR', linkToDoc)) 
    else
        if ~usrSuppliedHessMult
            H = sparse(H);
        end
        A = sparse(A); Aeq = sparse(Aeq);
        % Call sqpmin when just bounds or just equalities
        [X,fval,output,exitflag,lambda] = sqpmin(f,H,mtxmpy,X0,Aeq,Beq,lb,ub,verbosity, ...
            options,defaultopt,computeLambda,computeConstrViolation,varargin{:});
        
        if exitflag == -10  % Problem not handled by sqpmin at this time: dependent rows
            error(message('optim:quadprog:DepConstrTRR', linkToDoc))
        end
    end
    output.linearsolver = [];
end
        
% Compute fval and first-order optimality if the interior-point-convex
% algorithm was run (not stopped in presolve)
if (strcmpi(output.algorithm,interiorPointConvex) && ~isempty(lambda))
    % Compute objective function value
    fval = 0.5*X'*(H*X)+f'*X;
   
   % Compute first order optimality if needed
   if computeFirstOrderOpt && ~isempty(lambda)
      output.firstorderopt = computeKKTErrorForQPLP(H,f,A,B,Aeq,Beq,lb,ub,lambda,X); 
   else
      output.firstorderopt = []; 
   end
   output.cgiterations = [];  
end
