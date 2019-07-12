clc;
close all;
clear all;

% Scanning input from the files
% Here are all the variables that we need to start
fid1=fopen('bg.txt');
s1=textscan(fid1,'%f,%f,%f,%f,%f');
fclose(fid1);
label1 = zeros(100, 5);
label0 = zeros(100, 5);
allFeatures = zeros(200, 5);
for k=1:5
input = s1{k};
label1(:,k) = input(1:100);
label0(:,k) = input(101:200);

allFeatures(:,k) = input;
end
allFeatures(:,6) = 1;
% Now here we start the soft svm code
% Set number of operations
T = 100;
u = (1:1:T);

% Set gamma here
gamma = 0.001;

% The weight vector is with 5 different variables
w = zeros(T,6);

% This is the subgradient
v = zeros(100, 6);

% Theta
theta = zeros(100, 6);

% Binay loss
binaryLoss = zeros(T, 1);
hingeLoss = zeros(T, 1);


% Now set Xi and Ti 
Xi = allFeatures(1,:);
Ti = 1; 
for j = 1:(T-1)
    % Current wait is w
    w(j,:) = (1/(gamma*j))*theta(j, :);
    wj = w(j,:);
    
    % Choosing a random variable
    index = randi([1 200]);
   
    % The random variable for X
    Xj = allFeatures(index,:);
    
    % Find if t is 1 or -1
    if (index > 100)
        t = -1;
    else 
        t = 1;
    end
    
    % Dot product is t multiplied dot product of x and w
    dotProduct = t * dot(Xj, wj);    
    % Lets set theta first
    if ((1-dotProduct) > 0)
        v(j,:) = -t * Xj;
        summation = 0;
        for p=1:j
            summation = summation + v(p,:);
        end
        theta(j,:) = -1 * summation;
    else 
        v(j,:) = 0;
        summation = 0;
        for p=1:j
            summation = summation + v(p,:);
        end
        theta(j,:) = -1 * summation;    end 
    
    
    
    % If X is misclassified then set theta accordingy
    if dotProduct < 0
        theta(j+1, :) = theta(j,:) + t * (Xj);
    else 
        theta(j+1, :) = theta(j,:);
    end
        
    % Binary Loss finding
    for (m=1:200)
        X2 = allFeatures(m,:);
        if (m > 100)
            t2 = -1;
        else 
            t2 = 1;
        end
        m
        dotProduct2 = t2 * dot(X2, wj)
        if (dotProduct2 <= 0)
            binaryLoss(j) = binaryLoss(j) + 1;
        end
    end
    % Hinge Loss finding
    for (m=1:200)
        X2 = allFeatures(m,:);
        if (m > 100)
            t2 = -1;
        else 
            t2 = 1;
        end
        dotProduct2 = t2 * dot(X2, wj);
        if (dotProduct2 > 0)
            hingeLoss(j) = hingeLoss(j) + 0;
        else
            hingeLoss(j) = hingeLoss(j) + (1 - (dotProduct2));
        end
    end
end
hingeLoss(T) = hingeLoss(T-1);
hingeLoss = hingeLoss/200;
figure();
plot(u, hingeLoss);
xlabel('iterations'); 
ylabel('Hinge Loss'); 
title('Hinge Loss vs iteration'); 
binaryLoss(T) = binaryLoss(T-1);

figure();
plot(u, binaryLoss);
xlabel('iterations'); 
ylabel('Binary Loss'); 
title('Binary Loss vs iteration');
