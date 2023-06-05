clear;
clc;
close all;

%% Generating weight vector w
%heat, mass, cost
w_array = [0.9 0.05 0.05;
    0.05 0.9 0.05;
    0.05 0.05 0.9;
    1/3 1/3 1/3;
    0.6 0.2 0.2;
    0.2 0.6 0.2;
    0.2 0.2 0.6;];

%adding random points
for randpoints = 1:5
    r = rand(1,3);
    r = r/sum(r);
    w_array = [w_array;r];
end


%% Initial uniform sample generation
m_init = 1e5;
m = 400;
m_elite = 200;

%Initial bounds
min_MLI = 0.0001;
max_MLI = 1;

min_rego = 0.0001;
max_rego = 8;

min_aero = 0.0001;
max_aero = 2;

% Running through all weights for all Pareto optimal points
for pareto = 1:size(w_array,1)
    w = w_array(pareto,:);
    for metaiter = 1:20
%         metaiter

        %Generating uniform samples for first iteration
        r_MLI = unifrnd(min_MLI, max_MLI, 1, m_init);
        r_rego = unifrnd(min_rego, max_rego, 1, m_init);
        r_aero= unifrnd(min_aero, max_aero, 1, m_init);
        
        samples = horzcat(r_MLI', r_rego', r_aero');
        
        for i = 1:m
            x = samples(i,:);
            y(i) = obj(x,w);
        end
        
        %Iterating to find minumum
        iter = 100;
    %     figure;
        
        for j = 1:iter
    %%Plotting of sample point cloud
    %         if j > 1
    %             scatter3(samples(:,1), samples(:,2), samples(:,3))          
    %             hold on;
    %             scatter3(samples_elite(:,1), samples_elite(:,2), samples_elite(:,3))
    %             hold off;
    % 
    %             title("Sample Points")
    %             xlim([min_MLI max_MLI]);
    %             xlabel('MLI thickness')
    %             ylim([min_rego max_rego]);
    %             ylabel('Regolith thickness')
    %             zlim([min_aero, max_aero]);
    %             zlabel('Aerogel thickness');
    %             legend("Samples", "Elite Samples")
    %             
    %             drawnow;
    %             pause(0.1);
    %         end
    
        
        %Sorting samples by objective value
            [y_sort ind_sort] = sort(y);
            samples_sort = samples(ind_sort,:);
        %Taking elite samples
            y_elite = y_sort(1:m_elite);
            samples_elite = samples_sort(1:m_elite,:);
        % Generating new distribution from elite samples
            sample_elite_cov = cov(samples_elite);
            sample_elite_mean = mean(samples_elite);
            newsamples = mvnrnd(sample_elite_mean, sample_elite_cov, m);
        
            samples = newsamples;
            for i = 1:m
                x = samples(i,:);
                y(i) = obj(x,w);
            end
        end
        x_best(metaiter,:) = sample_elite_mean;
        y_best(metaiter) = obj(sample_elite_mean,w);
        q_best(metaiter) = heat_loss(sample_elite_mean);
        m_best(metaiter) = mass_calc(sample_elite_mean);
        c_best(metaiter) = cost_calc(sample_elite_mean);
    end
%%Plotting of best points per metaiteration, for debugging
%     figure;
%     scatter3(x_best(:,1),x_best(:,2),x_best(:,3))
%     figure;
%     plot(y_best,'ro')
%     title('y')
%     figure;
%     plot(q_best,'ro')
%     title('q')
%     figure;
%     plot(m_best,'ro')
%     title('m')
%     figure;
%     plot(c_best,'ro')
%     title('c')
    
% Saving all Pareto optimal designs
    [y_final(pareto), ind_final] = min(y_best);
    x_final(pareto,:) = x_best(ind_final,:);
    q_final(pareto) = q_best(ind_final);
    m_final(pareto) = m_best(ind_final);
    c_final(pareto) = c_best(ind_final);
    


end
%Compiling results into a table
y_final = y_final';
q_final = q_final';
m_final = m_final';
c_final = c_final';
T = table(w_array, x_final, y_final, q_final, m_final, c_final)
%% Results Plotting
figure;
scatter3(x_final(:,1),x_final(:,2),x_final(:,3));
title('Design Space, Optimal Points')
xlabel('MLI Thickness [m]');
ylabel('Regolith Thickness [m]');
zlabel('Aerogel Thickness [m]');

figure;
s = scatter3(q_final, m_final, c_final);
s.LineWidth = 0.6;
s.MarkerEdgeColor = 'b';
s.MarkerFaceColor = [0 0.5 0.5];
title('Objective Space, Optimal Points')
xlabel('Heat Loss [W]');
ylabel('Mass [kg]');
zlabel('Cost [$]');




% %% Debugging
% mlipoints = ones(100,1)*0.1;
% regopoints = linspace(min_rego, max_rego, 100);
% aeropoints = linspace(min_aero, max_aero, 100);
% 
% [regomesh, aeromesh] = meshgrid(regopoints, aeropoints);
% 
% % xtest = [mlipoints regopoints' aeropoints'];
% 
% for i = 1:100
%     for j = 1:100
%         x = [0.001 regomesh(i,j) aeromesh(i,j)];
%         z(i,j) = obj(x,w);
%     end
% end
% figure;
% surf(regomesh, aeromesh, z)






function func = func(x, w)

    q_max = 10000;
    m_max = 2.7e8;
    c_max = 2e9;

    q = (heat_loss(x)/q_max)*100;
    m = (mass_calc(x)/m_max)*100;
    c = (cost_calc(x)/c_max)*100;

    func = dot(w,[q m c]);
end
function p = p(x)
    p = 0;
    q_max = 10000;
    m_max = 2.7e8;
    c_max = 2e9;

    if heat_loss(x)>q_max
        p = p+1;
    end    
    if mass_calc(x)>m_max
        p = p+1;
    end
    if cost_calc(x)>c_max
        p = p+1;
    end    
    if sum(x<0)
        p = p+1;
    end

end

function p = pquad(x)
    q_max = 10000;
    m_max = 2.7e8;
    c_max = 2e9;

    q = (heat_loss(x)/q_max)*100;
    m = (mass_calc(x)/m_max)*100;
    c = (cost_calc(x)/c_max)*100;
    
    p = max(q-100, 0)^2 + max(m-100, 0)^2 + max(c-100,0)^2;

end
function obj = obj(x,w)
    obj = func(x,w) + 1e4*p(x) + 10*pquad(x);
end






