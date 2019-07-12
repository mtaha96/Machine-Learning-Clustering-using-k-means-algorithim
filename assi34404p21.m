clc;
close all;
clear all;

% Scanning input from the files
fid1=fopen("C:\\Users\\Mohammed\\Desktop\\York Final Semester\\EECS 4404\\twodpoints.txt");
s1=textscan(fid1,'%f,%f');
fclose(fid1);
input = s1{1};
output = s1{2};
inputs = [input output];
figure();




%Random Points
index = randi([1,200],1,1);
center1 = inputs(index, :);
index = randi([1,200],1,1);
center2 = inputs(index, :);
index = randi([1,200],1,1);
center3 = inputs(index, :);
center = [center1;center2;center3];
cost(1) = findClusters(s1, center, 3, "Random Points");
figure();


% Second chose centers by random points
n = size(input,1);
index = randi([1,200],1,1);
center1 = [input(index,1), output(index,1)];
center(1,:) = center1;
cost(1) = findClusters(s1, center, 1, "Maximum distance of Centers");

% Second chose centers by random points
n = 120;
index = randi([1,120],1,1);
center1 = inputs(index, :);
center(1,:) = center1;
cost(1) = findClusters(s1, center, 1, "Maximum distance of Centers")
for (k = 2:10)
    for m = 2:k
        max = -10000;
        for q=1:n
            dist = 0;
            newPoint = inputs(q,:);
            already = 0;
            for (w = 1:m-1)
                if (newPoint == center(w,:))
                    already = 1;
                end
            end
            if (already == 0)
                for (w = 1:m-1)
                    dist = dist +  pdist([center(w,:);newPoint]);
                end
                if (dist > max)
                    center2 = newPoint;
                    max = dist;
                end
            end
        end
        center(m,:) = center2;
    end
    cost(k) = findClusters(s1, center, k, "Maximum distance of Centers");
end
plot(center(:,1), center(:,2), 'o', 'MarkerFaceColor', 'g')
hold on 
hold off
figure();
plot([1:10], cost);
xlabel('K - Number of clusters'); 
ylabel('Clustering Cost'); 
title('Cost vs number of clusters');
grid on
grid minor

%plot(input, output, '.');
    
function [result] = findClusters(data, centers, k, label)
initialcenters = centers;
input = data{1};
output = data{2};
% Second chose centers by random points
n = size(input,1);

% Number of iteration used for stopping criteria
N = 10;
%Cluster Assigning
c = zeros(n,1);
for i = 1:N
    %Assign Clusters
    for q = 1:n
        minD = 10000000;
        for m = 1:k
            dist = pdist([input(q),output(q);centers(m,1),centers(m,2)]);
            if (dist < minD)
                c(q) = m;
                minD = dist;
            end
        end    
    end   
    % Assign new centers
    
    % Now for each cluster find new center
    for m = 1:k
        % for each cluster
        sum = [0,0];
        count = 0;
       for q = 1:n
           % find points belonging to that cluster
           if c(q)  == m
               sum = sum + [input(q), output(q)];
               count = count +1;
           end
       end
       centers(m,:) = sum/count;
    end
    
    % Assign clusters after new centers
     for q = 1:n
        minD = 10000000;
        for m = 1:k
            dist3 = pdist([input(q),output(q);centers(m,1),centers(m,2)]);
            if (dist3 < minD)
                c(q) = m;
                minD = dist3;
            end
        end    
     end 
end

    % Finding the cost
    cost = 0;
     for (q=1:n)
         clus = c(q);
         cost = cost + pdist([input(q),output(q);centers(clus,:)]);
     end
    result = cost;
    for (p = 1:k)
        clus1 = [0,0];
        c1 = 1;
        for j=1:n
            if (c(j) == p)
                clus1(c1,:) = [input(j), output(j)];
                c1 = c1 + 1;
            end
        end
        
        % Plot the different clusters
        plot(clus1(:,1), clus1(:,2), 'x');
        hold on
    end
    title(label);
end
