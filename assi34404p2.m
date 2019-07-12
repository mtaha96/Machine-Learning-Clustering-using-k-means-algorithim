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
plot(input, output, 'x');
title("data points scatter diagram");
figure();

k=3;
% Second chose centers by random points
n = size(input,1);
index = randi([1,200],1,1);
center1 = [input(index,1), output(index,1)];
center(1,:) = center1;
for (k = 2:3)
    for m = 2:k
        newCenter = [0,0];
        for (s = 1:(m-1))
            newCenter = newCenter + center(s,:);
        end
        newCenter = newCenter / (m-1)
        max = -10000;
        for q=1:n
            newPoint = [input(q),output(q)];
            already = 0;
            for (w = 1:m-1)
                if (newPoint == center(w,:))
                    already = 1;
                end
            end
            if (already == 0)
                dist3 = pdist([newCenter;input(q),output(q)]);
                if (dist3 > max)
                    center2 = [input(q), output(q)];
                    max = dist3;
                end
            end
        end
        center(m,:) = center2;
    end
    cost(k) = findClusters(s1, center, k, "Maximum distance of Centers");
end
plot(center(:,1), center(:,2), 'o', 'MarkerFaceColor', 'g')
hold on 
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
    cost = 0;
     for (q=1:n)
         clus = c(q);
         cost = cost + pdist([input(q),output(q);centers(clus,:)]);
     end
    result = cost;
end
