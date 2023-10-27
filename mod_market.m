clear all; close all; clc;
% Define the number of agents and tasks
num_agents = 5;
num_tasks = 5;

% Create random locations for agents and tasks
agent_locations = rand(num_agents, 2);
task_locations = rand(num_tasks, 2);

% Initialize agent and task matrices
agents = (1:num_agents)'; % Assign each agent to its index-matched task
tasks = (1:num_tasks)';

% Initialize the cost matrix %rows are task cost for one agent
cost_matrix = zeros(num_agents, num_tasks);

% Initialize the assigned task array
assigned_agents_storage = [];

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
disp("Note: rows are task cost for one agent")
disp(' ');

% Initialize a flag to check if all tasks are assigned
all_tasks_assigned = false;
task_allocations = zeros(num_tasks, 2);
task_allocations_exact_assigment = zeros(num_tasks, 2);
while all_tasks_assigned == false
    for auctioneer = 1:num_agents
            % Announce the task to other robots
            announced_task = tasks(auctioneer);
            disp(['Robot r' num2str(auctioneer) ' announces Task t' num2str(announced_task)]);
            
            % Bidding Phase
            bids = zeros(1, num_agents);
            for bidder = 1:num_agents
                %if bidder == auctioneer
                 %   continue; % Skip the auctioneer
                %end
                
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
            
            %Remove Assigned agents from agent list
            if ~ismember(winner, assigned_agents_storage)
                if numel(assigned_agents_storage) <= num_agents
                    assigned_agents_storage(end+1) = winner;
                end
            end
    
            %[unique_ones,~,idx] = unique(task_allocations(:,2));
    
            task_allocations_table = array2table(task_allocations, 'VariableNames', {'Task', 'Allocation'});
           
    
            % Mark the task as assigned
            %tasks(auctioneer) = "";
            
    end
    
    
    disp('task_allocations:');
    disp(task_allocations_table);
    
    disp('assigned tasks:');
    disp(assigned_agents_storage)
    
    tasks_to_agents = [];
    correct_assignment_agents = [];
    correct_assignment_tasks = [];
    unassigned_tasks = [];
    unsasinged_agents = [];
    
    for agent = 1:num_agents
        colNames = {'Task', 'Allocation'};
        results = array2table(zeros(0, length(colNames)), 'VariableNames', colNames);
        rows = task_allocations_table.Allocation == agent;
        result = task_allocations_table(rows,:);
        results = [results;result];
        %result_array = table2array(result);
        disp(results);
        results_array = table2array(result);
        if size(results_array, 1)>1
            for i = 1:size(results_array, 1)
                tasks_to_agents = [tasks_to_agents; results_array(i,:)];
            end
        % Add Already One-Matched Agents To Them Task Array
        elseif size(results_array, 1)==1
            % Already has 1 result so nothing wrong with using (1,1) constant
            % index parameter
            correct_assignment_tasks(end + 1) = results_array(1,1);
            correct_assignment_agents(end + 1) = agent;
        else
            unsasinged_agents(end + 1) = agent;
        end
    
        %tasks_to_agents(end + 1) = table2array(result);
        
    end
    
    disp('tasks to agents:');
    disp(tasks_to_agents)
    
    disp('Correct Agent Assigments:');
    disp(correct_assignment_agents)
    disp('Correct Task Assigments:');
    disp(correct_assignment_tasks)
    
    % Get Task With Lowest Cost If There Is An Agent Assigned to Multiple Tasks
    for i = 1:num_agents
        results_table = array2table(tasks_to_agents, 'VariableNames', {'Task', 'Allocation'});
        selected_rows = results_table.Task(results_table.Allocation == i);
        disp("Multi-Tasks For R" + i);
        disp(selected_rows);
        min_cost = inf;
        if size(selected_rows, 1)>0
            for index = 1:size(selected_rows, 1)
                iter_cost = cost_matrix(i,selected_rows(index));
                if iter_cost<min_cost 
                    min_cost = iter_cost;
                    min_cost_task = selected_rows(index); 
                end
            end
            disp("Min Cost is " + num2str(min_cost) + " and task is " + num2str(min_cost_task));
            % Add Lowest Task For Current Agent To Them Task Array
            correct_assignment_tasks(end + 1) = min_cost_task;
            correct_assignment_agents(end + 1) = i;
            % Now Remove Assigned Tasks row from table and add remanined rows
            % to unassigned
            %assigned_row = results_table.Task == min_cost_task;
            %results_table(results_table.Task == min_cost_task,:) = [];
            unassigned_tasks = setdiff(tasks, correct_assignment_tasks);
            unassigned_agents = setdiff(agents, correct_assignment_agents);
            %disp("Multi-Tasks For R" + i)
        end
    end
    if isempty(unassigned_agents) && isempty(unassigned_tasks)
        all_tasks_assigned = true;
    else
        num_agents = length(unassigned_agents);
        num_tasks = length(unassigned_tasks);
        tasks = unassigned_tasks;
        agents = unassigned_agents;
    end
end

disp(' ');
disp("After Simplfication Of Multi-Assigned Agents");
disp('Correct Agent Assigments:');
disp(correct_assignment_agents)
disp('Correct Task Assigments:');
disp(correct_assignment_tasks)
disp('Unassigned Agents:');
disp(unassigned_agents)
disp('Unassigned Tasks:');
disp(unassigned_tasks)


