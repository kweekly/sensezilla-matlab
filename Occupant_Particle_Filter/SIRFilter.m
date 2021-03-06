function [Xest, West, Xtraj] = SIRFilter( X0, W0, Y, T, transition_model, observation_model, params )
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
%   [Xnext] = transition_model( X, dt, params ) : Takes a matrix of states, X (each
%   column is a state vector), and returns a matrix estimated next states.
%   dt is the time since the last time it was called (or -1 if unknown).
%   params is a parameter structure.
%
%   [Pobs] = observation_model( Y, X, params )  : Takes an observation vector, Y,
%   and matrix of states, X, and returns a vector of probabilities, Pobs,
%   one for each state provided. params is a parameter structure.
%
%   Structural arguments:
%   params : parameter structure to pass to transition_model and
%   observation_model
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
if (~isa(params,'struct'))
   error('params must be a structure') 
end
if (~isa(transition_model,'function_handle') || nargin(transition_model) ~= 3 || nargout(transition_model) ~= 1 && nargout(transition_model) ~= -1)
   error('transition_model must be a function handle with two inputs and one output')
end
if (~isa(observation_model,'function_handle') || nargin(observation_model) ~= 3 || nargout(observation_model) ~= 1  && nargout(observation_model) ~= -1)
   error('observation_model must be a function handle with two inputs and one output')
end

doframes = 0;
if isfield(params,'framedir')
    doframes = 1;
   
    framefig = figure;
    set(framefig,'Renderer','zbuffer');
    hold on;
    [obsj,obsi] = find(params.obsmap_hex > 0.2);
    hplot = hexPlot(obsi,obsj,'k',params.grid);
    view([0, 90]);
    axis equal;
    frameplt = plot(X0(1,:),X0(2,:),'g.','MarkerSize',1);
    gtplt = plot(params.gtdata_x(1,1),params.gtdata_x(1,2),'bo','MarkerSize',10);
    maxplt = plot(X0(1,1),X0(1,1),'rx','MarkerSize',10);
    drawnow;
    if 0==exist(params.framedir,'dir')
       mkdir(params.framedir);
    end
    [ffname,junk] = sprintf('%s/%05d.png',params.framedir,0);
    saveas(framefig,ffname,'png');
end

% initialization (nothing to do, already done for us)
Xest = X0;
West = W0;
dt = -1;
Xtraj = zeros(size(X0,1),size(T,1));
for tidx=1:size(T,1)
    params.sim_t = T(tidx);
    if (mod(tidx,1) == 0)
       fprintf([num2str(tidx) ' of ' num2str(size(T,1)) ' samples processed.\n']); 
    end
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
    Xest = transition_model(Xest, dt, params);  
    % update the importance weights
    West = observation_model(Y(:,tidx),Xest, params);
    % fprintf(['\t' num2str(sum(West)) '\n']);
    
    % finding maximum weight particle
    %[w,i] = max(West);
    %Xtraj(:,tidx) = Xest(:,i);
    
    % finding mean of particles
    %Xtraj(:,tidx) = mean(Xest(:,~isnan(Xest(1,:)) & ~isnan(Xest(2,:))),2);
    
    %
    % weighted mean of particles
    Xwonan = Xest(:,~isnan(Xest(1,:)) & ~isnan(Xest(2,:)) & ~isnan(West'));
    Wwonan = West(~isnan(Xest(1,:)) & ~isnan(Xest(2,:)) & ~isnan(West'));
    Xtraj(:,tidx) = sum(Xwonan .* repmat(Wwonan',2,1),2)./sum(Wwonan);
    
    if doframes==1
        set(frameplt, 'XData',Xest(1,:), 'YData',Xest(2,:));
        t = T(tidx);
        set(gtplt, 'XData', interp1(params.gtdata_t,params.gtdata_x(:,1),t,'linear','extrap'), ...
                   'YData', interp1(params.gtdata_t,params.gtdata_x(:,2),t,'linear','extrap'));
        set(maxplt, 'XData', Xtraj(1,tidx), 'YData', Xtraj(2,tidx));
        drawnow;
        [ffname,junk] = sprintf('%s/%05d.png',params.framedir,tidx);
        saveas(framefig,ffname,'png');
    end
    
    if (tidx < size(T,1))
        dt = T(tidx+1) - T(tidx);
    end
end