% Define the number of agents and tasks
num_agents = 5;
num_tasks = 5;

% Create random locations for agents and tasks
agent_locations = rand(num_agents, 2);
task_locations = rand(num_tasks, 2);

% Initialize agent and task matrices
agents = (1:num_agents)'; % Assign each agent to its index-matched task
tasks = (1:num_tasks)';

% Initialize the cost matrix
cost_matrix = zeros(num_agents, num_tasks);

% Calculate the cost matrix
for i = 1:num_agents
    for j = 1:num_tasks
        % Calculate the Euclidean distance (Pythagorean distance)
        diff = agent_locations(i, :) - task_locations(j, :);
        cost_matrix(i, j) = norm(diff);
    end
end

% Display the initial cost matrix
disp("Cost Matrix:");
disp(cost_matrix);

% Initialize a flag to check if all tasks are assigned
all_tasks_assigned = false;
task_allocations = zeros(num_tasks, 2);

for auctioneer = 1:num_agents
        % Announce the task to other robots
        announced_task = tasks(auctioneer);
        disp(['Robot r' num2str(auctioneer) ' announces Task t' num2str(announced_task)]);
        
        % Bidding Phase
        bids = zeros(1, num_agents);
        for bidder = 1:num_agents
            if bidder == auctioneer
                continue; % Skip the auctioneer
            end
            
            % Calculate the bid (e.g., you can use a function of cost and distance)
            bid = 1 / cost_matrix(bidder, announced_task);
            bids(bidder) = bid;
            
            % Display bids
            disp(['  Robot r' num2str(bidder) ' bids ' num2str(bid)]);
        end
        
        % Determine the winner
        [max_bid, winner] = max(bids);
        
        % Assign the task to the winning robot
        assigned_task = announced_task;
        task_allocations(auctioneer, :) = [assigned_task, winner];
        disp(['  Robot r' num2str(winner) ' wins the auction and is assigned Task t' num2str(assigned_task)]);
        
        % Mark the task as assigned
        tasks(auctioneer) = "";
        
end

disp('task_allocations:');
disp(task_allocations);
