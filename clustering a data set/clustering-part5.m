clc;
close all;
clear all;

% Here are all the variables that we need to start
fid1=fopen('seeds_dataset.txt');
s1=textscan(fid1,'%f %f %f %f %f %f %f %d');
inputs = [s1{1} s1{2} s1{3} s1{4} s1{5} s1{6} s1{7}];
classes = s1{8};

idx = kmeans(inputs, 3);
% Second chose centers by random points
n = 210;
index = randi([1,120],1,1);
center1 = inputs(index, :);
center(1,:) = center1;
[cost(1),c] = findClusters(inputs, center, 1, "Maximum distance of Centers")
for (k = 2:3)
    for m = 2:k
        max = -10000;
        for q=1:n
            dist = 0;
            newPoint = inputs(q,:);
            already = 0;
            for (w = 1:m-1)
                dist = dist +  pdist([center(w,:);newPoint]);
            end
            if (dist > max)
                center2 = newPoint;
                max = dist;
            end
        end
        center(m,:) = center2;
    end
    [cost(k), c] = findClusters(inputs, center, k, "Maximum distance of Centers");
end
plot3(center(:,1), center(:,2), center(:,3), 'o', 'MarkerFaceColor', 'g');
hold on 
hold off
figure();
%plot([1:10], cost);
xlabel('K - Number of clusters'); 
ylabel('Clustering Cost'); 
title('Cost vs number of clusters');
grid on
grid minor
count = 0;
for (i=1:n)
    if (classes(i) ~= c(i))
        count = count + 1;
    end
end
count
count/210
%plot(input, output, '.');
    
function [result, c] = findClusters(inputs, centers, k, label)
initialcenters = centers;
% Second chose centers by random points
n = 210;

% Number of iteration used for stopping criteria
N = 50;
%Cluster Assigning
c = zeros(n,1);
for i = 1:N
    %Assign Clusters
    for q = 1:n
        minD = 10000000;
        for m = 1:k
            dist = pdist([inputs(q,:);centers(m,:)]);
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
        sum = [0,0,0,0,0,0,0];
        count = 0;
       for q = 1:n
           % find points belonging to that cluster
           if c(q)  == m
               sum = sum + inputs(q,:);
               count = count +1;
           end
       end
       centers(m,:) = sum/count;
    end
     for q = 1:n
        minD = 10000000;
        for m = 1:k
            dist3 = pdist([inputs(q,:);centers(m,:)]);
            if (dist3 < minD)
                c(q) = m;
                minD = dist3;
            end
        end    
     end 
end
    cost = 0;
     for (q=1:n)
         clus = c(q);
         cost = cost + pdist([inputs(q,:);centers(clus,:)]);
     end
    result = cost;
    for (p = 1:k)
        clus1 = [0,0,0,0,0,0,0];
        c1 = 1;
        for j=1:n
            if (c(j) == p)
                clus1(c1,:) = inputs(j,:);
                c1 = c1 + 1;
            end
        end
        plot3(clus1(:,1), clus1(:,2), clus1(:,3), 'x');
        hold on
    end
end
