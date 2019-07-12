clc;
close all;
clear all;

% Scanning input from the files
fid1=fopen("C:\\Users\\Mohammed\\Desktop\\York Final Semester\\EECS 4404\\twodpoints.txt");
s1=textscan(fid1,'%f,%f');
fclose(fid1);
input = s1{1}*100;
output = s1{2}*100;
inputs = [input output];
plot(inputs(:,1), inputs(:,2), 'x');
figure();
% k-means algorithim
k = 3;
% The centers that we have
centers = zeros(k, 2);
centers(1,:) = [-4,-4];
centers(2,:) = [3,-2];
centers(3,:) = [2,5];
findClusters(s1, centers, k, "By hand Centers selection 1");

% k-means algorithim
k = 3;
% The centers that we have
centers = zeros(k, 2);
centers(1,:) = [4,4];
centers(2,:) = [4,2];
centers(3,:) = [4,7];
findClusters(s1, centers, k, "By hand Centers selection 2");


% Second chose centers by random points
n = size(input,1);
% The centers that we have
centers = zeros(k, 2);
index = randi([1,200],1,1);
center1 = [input(index,1), output(index,1)];
centers(1,:) = center1;
center2 = [0,0];
max = -10000;
% Second chose centers by random points
index = randi([1,200],1,1);
center1 = inputs(index, :);
center(1,:) = center1;
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
    if k == 3
    findClusters(s1, center, k, "Maximum distance of Centers");
    end
end
    
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
    c1 = 1;
    c2 = 1;
    c3 = 1;
    for j=1:n
        if (c(j) == 1)
            clus1(c1,:) = [input(j), output(j)];
            c1 = c1 + 1;
        elseif (c(j) == 2)
            clus2(c2,:) = [input(j), output(j)];
            c2 = c2 + 1;
        elseif (c(j) == 3)
            clus3(c3,:) = [input(j), output(j)];
            c3 = c3 + 1;
        end
    end
        figure();
        
    plot(clus1(:,1), clus1(:,2), ['.', 'r'], clus2(:,1), clus2(:,2), ['.', 'b'], clus3(:,1), clus3(:,2), ['.', 'g']);
    hold on
        plot(centers(1,1),centers(1,2),['*','r'], centers(2,1),centers(2,2),['*','b'], centers(3,1),centers(3,2),['*','g']);
end
    figure();
    plot(clus1(:,1), clus1(:,2), ['.', 'r'], clus2(:,1), clus2(:,2), ['.', 'b'], clus3(:,1), clus3(:,2), ['.', 'g']);
    hold on
    plot(centers(1,1),centers(1,2),['*','r'], centers(2,1),centers(2,2),['*','b'], centers(3,1),centers(3,2),['*','g']);
    hold on
    cent1 = plot(initialcenters(1,1),initialcenters(1,2),'d', 'MarkerFaceColor', 'r', 'DisplayName','initial center for Cluster 1');
    hold on
    cent2 = plot(initialcenters(2,1),initialcenters(2,2),'d', 'MarkerFaceColor', 'b', 'DisplayName','initial center for Cluster 2');
    hold on
    cent3 = plot(initialcenters(3,1),initialcenters(3,2),'d', 'MarkerFaceColor', 'g', 'DisplayName','initial center for Cluster 3');
    legend([cent1 cent2 cent3])
    %, initialcenters(2,1),initialcenters(2,2),['o', 'b'], initialcenters(3,1),initialcenters(3,2),['o', 'k']);
    title(label);
    result = c;
end
