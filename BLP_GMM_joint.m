function ff = BLP_GMM_joint(x0)

global W tol_inner Kbeta iterativedelta Indicator share IV v A price constant indexxx error_J Ktheta z gg zero marginalcost zeta Total M T unobs_cost ns alpha_i simulatedshare RCderivative


theta1              = x0(1:Kbeta,1);                                                       % Mean tastes
theta2              = x0(Kbeta+1:Kbeta+Ktheta,1);                                          % Deviation from mean tastes
gamma               = x0(Kbeta+Ktheta+1:end,1);                                            % MC cost parameters (constant term & cost shifters)

% theta1= rand(Kbeta,1);
% theta2= rand(Ktheta,1);
% gamma  = rand(Kmc,1);

ii                  = 0;
norm_max            = 1;
deltavalue          = iterativedelta;    %   ~999*1   matrix              % With this variable, we will be able to iterate and 
                                                                           % update delta values alongside the inner loop
                                                                                              
muvalue             = -price*(theta2 .* v);                      % This represents our mu value in RC model. Here v enters
                                                                           % and allows to produce simulations for different consumers
                        %999*1   *      1*1  * 1*100         ~999*100 matrix
     
while norm_max > tol_inner && ii < 10000                                                      % "While" loop until convergence

    simulatedshare_market            = Indicator*(exp(deltavalue + muvalue)) + 1;             % 300*1 vector containing shares of markets + 1 in simulation.
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
error_J           = equilibriumdelta - [constant A price]*theta1;          % Gives us econometric error term (unobserved product specifications)
gg                = IV' * error_J;                                         % Moment function to be minimized


%% Supply Side
                                
alpha_i   = theta1(5,1).*ones(1,ns) - theta2(1,1).*v(1,:);         % Individual specific alpha term 

RCderivative = zeros(Total-M*T-1,1);  % YOKSA MINUS MI KOYMAM GEREKIYOR "mean" IFADESININ BASINA??
for        j = 1:Total-M*T-1
           RCderivative(j,:) = mean(alpha_i.*simulatedshare(j,:).*(1 - simulatedshare(j,:))); % gathering ds/dp from RC elasticity definition
end                           

unobs_cost   = price + share./RCderivative - [constant z]*gamma;  % Getting unonbserved costs

zeta         = IV' * unobs_cost;                                  % Moment function to be minimized
jointmoment  = [gg;zeta];                                         % Getting moments together 
ff           = jointmoment'*[W zero;zero W]*jointmoment;          % Objective Function

marginalcost = unobs_cost + [constant z]*gamma;                   % Resulting MC that we will use in 3rd question

end