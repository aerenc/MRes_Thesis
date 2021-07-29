function f = BLP_GMM_demand(x0)

global W tol_inner Kbeta iterativedelta v Indicator share IV A price constant g indexxx error


theta1              = x0(1:Kbeta,1);                                       % Mean tastes
theta2              = x0(Kbeta+1:end,1);                                   % Deviation from mean tastes

% theta1 = rand(Kbeta,1);
% theta2 = rand(Ktheta,1);


ii                  = 0;
norm_max            = 1;
deltavalue          = iterativedelta;    %   ~999*1   matrix              % With this variable, we will be able to iterate and 
                                                                           % update delta values alongside the inner loop
                                                                                              
muvalue             = -price*(theta2 .* v);                      % This represents our mu value in RC model. Here v enters
                                                                           % and allows to produce simulations for different consumers
                        %999*2   *      1*1  * 1*100         ~999*100 matrix
     
while norm_max > tol_inner && ii < 10000                                                      % "While" loop until convergence

    simulatedshare_market            = Indicator*(exp(deltavalue + muvalue)) + 1;             % 200*1 vector containing shares of markets + 1 in simulation.
                                                                                              % The variable with exp(.) contains the mean utilities 
                                                                                              % of products for consumers
                                                                                        
    simulatedshare_product           = repelem(simulatedshare_market,indexxx,1);              % Here we are basically stretching the row dimension of 
                                                                                              % above matrix by repelem() function: Each market share is
                                                                                              % expanded in column by the number of its products. 
                                                                                              
    simulatedshare                   = (exp(deltavalue + muvalue)) ./ simulatedshare_product; % This contains estimated shares of each products in every 
                                                                                              % market by simulation  
                                                                                              
    deltavaluenew                    = deltavalue + log(share) - log(mean(simulatedshare,2)); % Our main equation in inner loop: We are basically trying 
                                                                                              % to obtain a new delta such that the distance between the 
                                                                                              % old and new one could be considered negligible. The latter 
                                                                                              % variable represents estimated mean shares across consumers
                                                                                              
    norm_max                         = max(abs(deltavaluenew - deltavalue));                  % When this is below tol_inner, we are done  
    
    deltavalue                       = deltavaluenew;                                         % Here, we arrive to the end of the loop and a new delta
                                                                                              % value is produced by the reference of the old one
                                                                                              
    ii                               = ii + 1;                                                % Iterating for the next period
    
end

equilibriumdelta  = deltavaluenew;                                         % Gives us equilibrium delta value reached at the end of inner loop
error             = equilibriumdelta - [constant A price]*theta1;          % Gives us econometric error term (unobserved product specifications)
g                 = IV' * error;                                           % Moment function to be minimized
f                 = g'*W*g;                                                % Objective function


end

