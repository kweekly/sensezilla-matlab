function [Xest, West, Xtraj] = SIRFilter( X0, W0, Y, T, transition_model, observation_model )
%[Xest, West] = SIRFILTER( X0, Y, transition_model, observation_model )
%   Implements the Particle Filtering algorithm via
%   Sequential-Importance-Resampling ala. (Gordon et.al. 1993)
%
%   Numerical arguments:
%   X0 is a matrix of the initial state estimates. 
%           size(X,1) is the dimension of the state.
%           size(X,2) is the number of particles to be used.
%   W0 is a vector of the initial estimate weights.
%           size(W0,1) must equal size(X,2)
%   Y is a matrix of observations. 
%           size(Y,1) is the dimension of the observations. 
%           size(Y,2) is the number of observations.
%   T is a vector of times that the observations arrived.
%           size(T,1) must equal size(Y,2)
%
%   Functional arguments:
%   [Xnext] = transition_model( X, dt ) : Takes a matrix of states, X (each
%   column is a state vector), and returns a matrix estimated next states.
%   dt is the time since the last time it was called (or -1 if unknown).
%
%   [Pobs] = observation_model( Y, X )  : Takes an observation vector, Y,
%   and matrix of states, X, and returns a vector of probabilities, Pobs,
%   one for each state provided.
%
%   Return values:
%   Xest : The latest state estimates of the particles
%   West : The weights of the particles returned
%   Xtraj: Trajectory of the highest-weighted particle

if (~isvector(T))
    error('T must be a vector')
end
if (size(T,1) ~= size(Y,2))
    error('size(T,1) must equal size(Y,2)')
end
if (~isvector(W0))
    error('W0 must be a vector')
end
if( size(W0,1) ~= size(X0,2))
   error('size(W0,1) must equal size(X0,2)') 
end
if (~isa(transition_model,'function_handle') || nargin(transition_model) ~= 2 || nargout(transition_model) ~= 1 && nargout(transition_model) ~= -1)
   error('transition_model must be a function handle with two inputs and one output')
end
if (~isa(observation_model,'function_handle') || nargin(observation_model) ~= 2 || nargout(observation_model) ~= 1  && nargout(observation_model) ~= -1)
   error('observation_model must be a function handle with two inputs and one output')
end

% initialization (nothing to do, already done for us)
Xest = X0;
West = W0;
dt = -1;
Xtraj = zeros(size(X0,1),size(T,1));
for tidx=1:size(T)
    % normalize the importance weights
    if ( sum(West) == 0)
        West = 1/size(West,1);
    else
        West = West / sum(West);
    end
    % sample according to belief distribution
    Xidx = discretesample(West, size(Xest,2));
    Xest = Xest(:,Xidx);
    % predict the next state
    Xest = transition_model(Xest, dt);
    % update the importance weights
    West = observation_model(Y(:,tidx),Xest);
    
    [w,i] = max(West);
    Xtraj(:,tidx) = Xest(:,i);
    
    if (tidx < size(T))
        dt = T(tidx+1) - T(tidx);
    end
end