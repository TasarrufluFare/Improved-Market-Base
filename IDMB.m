% Creator: Tasarruflu Fare
% Date: 28/10/2023 - 11.37 PM

clear all; close all; clc;
% Define the number of agents and tasks
num_agents = 50;
num_agents_start = 50;
num_tasks = 50;
num_tasks_start = 50;

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
cost_matrix_start = cost_matrix;

% Initialize a flag to check if all tasks are assigned
all_tasks_assigned = false;

task_allocations_exact_assigment = zeros(num_tasks, 2);
iter = 0;

% Arrays To Save Correct Matches
correct_assignment_agents = [];
correct_assignment_tasks = [];

% Arrays To Calculate Non-Correct Matches
unassigned_tasks = 1:num_tasks_start;
unsasinged_agents = 1:num_agents_start;

while all_tasks_assigned == false
    iter = iter + 1;
    disp(["While Loop Iteration: "  num2str(iter) " ----------------------------------- "])
    disp(" ")
    disp(" ")
    disp(" ")

    if num_agents < 1
        continue
    end
    
    task_allocations = zeros(num_tasks, 2); % Array To Save Matches In Current Iteration. Later I will create Tables using this array 
    for auctioneer = 1:num_agents
            %auctioneer = agents(auctioneer_index);
            % Announce the task to other robots
            announced_task = tasks(auctioneer);
            task_index = auctioneer;
            if announced_task == 0
                continue
            end
            
            disp(['Robot r' num2str(agents(auctioneer)) ' announces Task t' num2str(announced_task)]);
            
            % Bidding Phase
            bids = zeros(1, num_agents);
            for bidder_index = 1:num_agents
                bidder = agents(bidder_index);
                
                % Calculate the bid (e.g., you can use a function of cost and distance)
                bid = 1 / cost_matrix(bidder_index, task_index);
                bids(bidder_index) = bid;
                
                % Display bids
                disp(['  Robot r' num2str(bidder) ' bids ' num2str(bid)]);
            end
            
            % Determine the winner
            [max_bid, winner_index] = max(bids);

            
            % Assign the task to the winning robot
            assigned_task = announced_task;
            if tasks(auctioneer) ~= 0
                task_allocations(tasks(auctioneer), :) = [assigned_task, agents(winner_index)];
            end
            
            disp(['  Robot r' num2str(agents(winner_index)) ' wins the auction and is assigned Task t' num2str(assigned_task)]);
            
            %Remove Assigned agents from agent list
            if ~ismember(agents(winner_index), assigned_agents_storage)
                if numel(assigned_agents_storage) < num_agents_start
                    % To Store Assigned Agents
                    assigned_agents_storage(end+1) = agents(winner_index);
                end
            end
    
    
            task_allocations_table = array2table(task_allocations, 'VariableNames', {'Task', 'Allocation'});    
    end
    disp('task_allocations:');
    disp(task_allocations_table);
    
    disp('assigned agents:');
    disp(assigned_agents_storage)

    tasks_to_agents = [];
    for agent = 1:num_agents
        colNames = {'Task', 'Allocation'};
        results = array2table(zeros(0, length(colNames)), 'VariableNames', colNames);
        rows = task_allocations_table.Allocation == agents(agent);
        result = task_allocations_table(rows,:);
        results = [results;result];
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
         correct_assignment_agents(end + 1) = agents(agent);
       end
        
        
    end
    
    disp('tasks to agents:');
    disp(tasks_to_agents)
    
    disp('Correct Agent Assigments - 1:');
    disp(correct_assignment_agents)
    disp('Correct Task Assigments - 1:');
    disp(correct_assignment_tasks)

    % Get Task With Lowest Cost If There Is An Agent Assigned to Multiple Tasks
    for agent_index = 1:num_agents
        i = agents(agent_index);

        if ~isempty(tasks_to_agents)
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
                unassigned_tasks = setdiff((1:num_tasks_start), correct_assignment_tasks);
                unassigned_agents = setdiff((1:num_agents_start), correct_assignment_agents);

            end
        end

        
        
    end

    if (numel(correct_assignment_agents) == num_agents_start || iter >= 15)
        all_tasks_assigned = true;
        disp("Results:")
        disp("Assigned Agents In Order:")
        disp(correct_assignment_agents)
        disp(" ")
        disp("Assigned Tasks In Order:")
        disp(correct_assignment_tasks)

    else
        disp('Unassigned Agents - 1:');
        disp(unassigned_agents)
        disp('Unassigned Tasks - 1:');
        disp(unassigned_tasks)

        agents = unassigned_agents;
        tasks = unassigned_tasks;
        num_agents = length(agents);
        num_tasks = length(tasks);
        new_cost_matrix = zeros(num_agents, num_tasks);
        for agent_index = 1:num_agents
            for task_index = 1:num_tasks
                new_cost_matrix(agent_index,task_index) = cost_matrix(agents(agent_index), tasks(task_index));
            end
        end
        disp("New Cost Matrix - Assigned")
        disp(new_cost_matrix)
    end
end

% Improved Market Based
swap_count = 0;
for agent_ref_index = 1:numel(correct_assignment_agents)
    for agent_comp_index = 1:numel(correct_assignment_agents)
        if agent_ref_index == agent_comp_index
            continue
        end
        index_ref_task = agent_ref_index;
        index_comp_task = agent_comp_index;
        cost_ref = cost_matrix_start(correct_assignment_agents(agent_ref_index), correct_assignment_tasks(index_ref_task));
        cost_ref_swap = cost_matrix_start(correct_assignment_agents(agent_ref_index), correct_assignment_tasks(index_comp_task));
        cost_comp = cost_matrix_start(correct_assignment_agents(agent_comp_index), correct_assignment_tasks(index_comp_task));
        cost_comp_swap = cost_matrix_start(correct_assignment_agents(agent_comp_index), correct_assignment_tasks(index_ref_task));
        cost_sum = cost_ref + cost_comp;
        cost_sum_swap = cost_ref_swap + cost_comp_swap;
        if cost_sum_swap < cost_sum
            swap_count = swap_count+1;
            % Swap Elements In Tasks Array If Necessary
            temp = correct_assignment_tasks(index_ref_task);
            correct_assignment_tasks(index_ref_task) = correct_assignment_tasks(index_comp_task);
            correct_assignment_tasks(index_comp_task) = temp;
        end
    end
end

disp("Swap Results:")
disp("Assigned And Swapped Agents In Order:")
disp(correct_assignment_agents)
disp(" ")
disp("Assigned And Swapped Tasks In Order:")
disp(correct_assignment_tasks)
disp("Swap Count:")
disp(swap_count)

% Calculate Total Cost
total_cost = 0;
for i = 1:length(correct_assignment_agents)
    cost_to_add = cost_matrix_start(correct_assignment_agents(i), correct_assignment_tasks(i));
    total_cost = total_cost + cost_to_add;
end

disp("Total Cost")
disp(total_cost)
disp(" ")

% Draw Graphic And Visualize
% Create figure
figure;
hold on;

% Draw Agent Locations using blue circles
scatter(agent_locations(:, 1), agent_locations(:, 2), 100, 'filled', 'MarkerFaceColor', 'b');

% Draw Task Locations using red triangles
scatter(task_locations(:, 1), task_locations(:, 2), 100, 'filled', 'MarkerFaceColor', 'r', 'Marker', '^');


% Draw Green Lines To Visualize Matches
for i = 1:length(correct_assignment_agents)
    agent_idx = correct_assignment_agents(i);
    task_idx = correct_assignment_tasks(i);
    agent_loc = agent_locations(agent_idx, :);
    task_loc = task_locations(task_idx, :);
    plot([agent_loc(1), task_loc(1)], [agent_loc(2), task_loc(2)], 'g');
end

for i = 1:length(correct_assignment_tasks)
    task_number = correct_assignment_tasks(i);
    pos = task_locations(task_number, :);
    text(pos(1)+0.01, pos(2)+0.01, ['T ' num2str(task_number)], 'Color', 'b', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
end

for i = 1:length(correct_assignment_agents)
    agent_number = correct_assignment_agents(i);
    pos = agent_locations(agent_number, :);
    text(pos(1)+0.01, pos(2)+0.01, ['A ' num2str(agent_number)], 'Color', 'r', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
end

% Set Axis And Titles
xlabel('X Coordinates');
ylabel('Y Coordinates');
title('Agent And Task Locations');

% Meanings
legend('Agents', 'Tasks', 'Matches');

% Set Axis
axis([0 1.02 0 1.02])

% Show On Screen
hold off;
