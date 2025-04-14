
%% Visualizing Microstate Prototypes with Directed Connectivity (Transition Probabilities)
% This script creates visualizations of EEG microstate prototypes arranged in a circular layout.
% Each figure highlights one microstate (MS) at a time as the "seed" and draws directional arrows
% pointing from this seed to all other microstates. The arrow thickness represents a scaled Bayes
% Factor (BF) indicating the strength of the directed effect from the seed to the other microstates.

%% 
% Figure for MS-B
% Define image file paths
MSA = 'path_to_microstate_prototypes/MS4.jpg';
MSB = 'path_to_microstate_prototypes/MS3.jpg';
MSC = 'path_to_microstate_prototypes/MS2.jpg';
MSD = 'path_to_microstate_prototypes/MS1.jpg';
MSE = 'path_to_microstate_prototypes/MS5.jpg';

image_paths = {MSB, MSA, MSE, MSD, MSC};
image_size = [75, 100]; % Define a uniform image size [width, height]

% Load and resize images
num_images = length(image_paths);
images = cell(1, length(image_paths));
for i = 1:length(image_paths)
    images{i} = imresize(imread(image_paths{i}), image_size);
end

% Define the circle parameters
circle_radius = 2 * max(image_size); % Adjust the circle radius
theta = linspace(0, 2*pi, num_images + 1); % Divide the circle into equal parts
theta(end) = []; % Remove the last element to avoid a duplicate of the first image

% Define the figure size dynamically based on the circle radius and image size
padding = max(image_size); % Space around the circle
figure_size = 2 * (circle_radius + padding);
fig = figure('Color', 'white', 'Position', [100, 100, figure_size, figure_size]);
hold on;
axis off

% Calculate positions and display images
for i = 1:num_images

    % Calculate image (x,y) positions
    x = circle_radius * cos(theta(i));
    y = circle_radius * sin(theta(i));

    % Calculate the normalized position
    x_norm = (x + figure_size/2 - image_size(1)/2) / figure_size;
    y_norm = (y + figure_size/2 - image_size(2)/2) / figure_size;

    % Create an axes for the image
    ax = axes('Position', [x_norm, y_norm, image_size(1)/figure_size, image_size(2)/figure_size]);

    % Display the image
    imshow(images{i}, 'Parent', ax);

    % Remove axis ticks
    set(ax, 'XTick', [], 'YTick', []);
end

axis equal;
axis off;


% Draw lines from one seed image to all other images

% Calculate the normalized start position MS-1: MS-B
ms = 1;
xdisplace = [image_size(1)/3 image_size(1)/2 image_size(1)/2 image_size(1)/3];
ydisplace = [-image_size(2)/3 -10 10 image_size(2)/3];
start_x = ((circle_radius * cos(theta(ms)) + figure_size/2 - image_size(1)/3) / figure_size);
start_y = ((circle_radius * sin(theta(ms)) + figure_size/2) / figure_size);

BFvals = [2.583050e+08 1.714031e+05 5.670481e+10 1.656750e+07];
BFvals = 1 + BFvals/sum(BFvals)*3;

for j = 2:num_images
    if j ~= ms
        % Calculate the normalized end position for the arrows (center of the target image)
        end_x = (circle_radius * cos(theta(j)) + figure_size/2 + xdisplace(j-1)) / figure_size;
        end_y = (circle_radius * sin(theta(j)) + figure_size/2 + ydisplace(j-1)) / figure_size;

        % Draw the arrow
        annotation(fig, 'arrow', [start_x, end_x], [start_y, end_y], 'HeadStyle', 'cback1','Color','black', 'LineWidth', BFvals(j-1));
    end
end

%% Figure for second MS-A

% Define image file paths
MSA = 'path_to_microstate_prototypes/MS4.jpg';
MSB = 'path_to_microstate_prototypes/MS3.jpg';
MSC = 'path_to_microstate_prototypes/MS2.jpg';
MSD = 'path_to_microstate_prototypes/MS1.jpg';
MSE = 'path_to_microstate_prototypes/MS5.jpg';

image_paths = {MSB, MSA, MSE, MSD, MSC};
image_size = [75, 100]; % Define a uniform image size [width, height]

% Load and resize images
num_images = length(image_paths);
images = cell(1, length(image_paths));
for i = 1:length(image_paths)
    images{i} = imresize(imread(image_paths{i}), image_size);
end

% Define the circle parameters
circle_radius = 2 * max(image_size); % Adjust the circle radius
theta = linspace(0, 2*pi, num_images + 1); % Divide the circle into equal parts
theta(end) = []; % Remove the last element to avoid a duplicate of the first image

% Define the figure size dynamically based on the circle radius and image size
padding = max(image_size); % Space around the circle
figure_size = 2 * (circle_radius + padding);
fig = figure('Color', 'white', 'Position', [100, 100, figure_size, figure_size]);
hold on;
axis off

% Calculate positions and display images
for i = 1:num_images

    % Calculate image (x,y) positions
    x = circle_radius * cos(theta(i));
    y = circle_radius * sin(theta(i));

    % Calculate the normalized position
    x_norm = (x + figure_size/2 - image_size(1)/2) / figure_size;
    y_norm = (y + figure_size/2 - image_size(2)/2) / figure_size;

    % Create an axes for the image
    ax = axes('Position', [x_norm, y_norm, image_size(1)/figure_size, image_size(2)/figure_size]);

    % Display the image
    imshow(images{i}, 'Parent', ax);

    % Remove axis ticks
    set(ax, 'XTick', [], 'YTick', []);
end

axis equal;
axis off;

% Calculate the normalized start position MS-2: MS-A
ms = 2;
xdisplace = [-image_size(1)/3 0 image_size(1)/2 image_size(1)/3 0];
ydisplace = [image_size(2)/3 0 10 20 image_size(2)/3];
start_x = ((circle_radius * cos(theta(ms)) + figure_size/2) / figure_size);
start_y = ((circle_radius * sin(theta(ms)) + figure_size/2 - 25) / figure_size);

BFvals = [0.59 0 0.09 0.13 4.84];
BFvals = 1 + BFvals/sum(BFvals)*3;

for j = 1:num_images
    if j ~= ms
        % Calculate the normalized end position for the arrows (center of the target image)
        end_x = (circle_radius * cos(theta(j)) + figure_size/2 + xdisplace(j)) / figure_size;
        end_y = (circle_radius * sin(theta(j)) + figure_size/2 + ydisplace(j)) / figure_size;

        % Draw the arrow
        annotation(fig, 'arrow', [start_x, end_x], [start_y, end_y], 'HeadStyle', 'cback1', 'LineWidth', BFvals(j));
    end
end


%% Figure for third MS-E

% Define image file paths
MSA = 'path_to_microstate_prototypes/MS4.jpg';
MSB = 'path_to_microstate_prototypes/MS3.jpg';
MSC = 'path_to_microstate_prototypes/MS2.jpg';
MSD = 'path_to_microstate_prototypes/MS1.jpg';
MSE = 'path_to_microstate_prototypes/MS5.jpg';

image_paths = {MSB, MSA, MSE, MSD, MSC};
image_size = [75, 100]; % Define a uniform image size [width, height]

% Load and resize images
num_images = length(image_paths);
images = cell(1, length(image_paths));
for i = 1:length(image_paths)
    images{i} = imresize(imread(image_paths{i}), image_size);
end

% Define the circle parameters
circle_radius = 2 * max(image_size); % Adjust the circle radius
theta = linspace(0, 2*pi, num_images + 1); % Divide the circle into equal parts
theta(end) = []; % Remove the last element to avoid a duplicate of the first image

% Define the figure size dynamically based on the circle radius and image size
padding = max(image_size); % Space around the circle
figure_size = 2 * (circle_radius + padding);
fig = figure('Color', 'white', 'Position', [100, 100, figure_size, figure_size]);
hold on;
axis off

% Calculate positions and display images
for i = 1:num_images

    % Calculate image (x,y) positions
    x = circle_radius * cos(theta(i));
    y = circle_radius * sin(theta(i));

    % Calculate the normalized position
    x_norm = (x + figure_size/2 - image_size(1)/2) / figure_size;
    y_norm = (y + figure_size/2 - image_size(2)/2) / figure_size;

    % Create an axes for the image
    ax = axes('Position', [x_norm, y_norm, image_size(1)/figure_size, image_size(2)/figure_size]);

    % Display the image
    imshow(images{i}, 'Parent', ax);

    % Remove axis ticks
    set(ax, 'XTick', [], 'YTick', []);
end

axis equal;
axis off;

% Calculate the normalized start position MS-3: MS-E
ms = 3;
xdisplace = [-25 -20 0 5 -10];
ydisplace = [10 -10 0 30 image_size(2)/3];
start_x = ((circle_radius * cos(theta(ms)) + figure_size/2 + 25) / figure_size);
start_y = ((circle_radius * sin(theta(ms)) + figure_size/2 -10) / figure_size);

BFvals = [.09 .16 0 76.46 .18];
BFvals = 1 + BFvals/sum(BFvals)*3;

for j = 1:num_images
    if j ~= ms
        % Calculate the normalized end position for the arrows (center of the target image)
        end_x = (circle_radius * cos(theta(j)) + figure_size/2 + xdisplace(j)) / figure_size;
        end_y = (circle_radius * sin(theta(j)) + figure_size/2 + ydisplace(j)) / figure_size;

        % Draw the arrow
        annotation(fig, 'arrow', [start_x, end_x], [start_y, end_y], 'HeadStyle', 'cback1', 'LineWidth', BFvals(j));
    end
end


%% Figure for fourth MS-D

% Define image file paths
MSA = 'path_to_microstate_prototypes/MS4.jpg';
MSB = 'path_to_microstate_prototypes/MS3.jpg';
MSC = 'path_to_microstate_prototypes/MS2.jpg';
MSD = 'path_to_microstate_prototypes/MS1.jpg';
MSE = 'path_to_microstate_prototypes/MS5.jpg';

image_paths = {MSB, MSA, MSE, MSD, MSC};
image_size = [75, 100]; % Define a uniform image size [width, height]

% Load and resize images
num_images = length(image_paths);
images = cell(1, length(image_paths));
for i = 1:length(image_paths)
    images{i} = imresize(imread(image_paths{i}), image_size);
end

% Define the circle parameters
circle_radius = 2 * max(image_size); % Adjust the circle radius
theta = linspace(0, 2*pi, num_images + 1); % Divide the circle into equal parts
theta(end) = []; % Remove the last element to avoid a duplicate of the first image

% Define the figure size dynamically based on the circle radius and image size
padding = max(image_size); % Space around the circle
figure_size = 2 * (circle_radius + padding);
fig = figure('Color', 'white', 'Position', [100, 100, figure_size, figure_size]);
hold on;
axis off

% Calculate positions and display images
for i = 1:num_images

    % Calculate image (x,y) positions
    x = circle_radius * cos(theta(i));
    y = circle_radius * sin(theta(i));

    % Calculate the normalized position
    x_norm = (x + figure_size/2 - image_size(1)/2) / figure_size;
    y_norm = (y + figure_size/2 - image_size(2)/2) / figure_size;

    % Create an axes for the image
    ax = axes('Position', [x_norm, y_norm, image_size(1)/figure_size, image_size(2)/figure_size]);

    % Display the image
    imshow(images{i}, 'Parent', ax);

    % Remove axis ticks
    set(ax, 'XTick', [], 'YTick', []);
end

axis equal;
axis off;

% Calculate the normalized start position MS-4: MS-D
ms = 4;
xdisplace = [-25 -10 5 5 -25];
ydisplace = [-5 -20 -25 30 10];
start_x = ((circle_radius * cos(theta(ms)) + figure_size/2 + 25) / figure_size);
start_y = ((circle_radius * sin(theta(ms)) + figure_size/2 + 15) / figure_size);

BFvals = [.12 .10 .26 0 .11];
BFvals = 1 + BFvals/sum(BFvals)*3;

for j = 1:num_images
    if j ~= ms
        % Calculate the normalized end position for the arrows (center of the target image)
        end_x = (circle_radius * cos(theta(j)) + figure_size/2 + xdisplace(j)) / figure_size;
        end_y = (circle_radius * sin(theta(j)) + figure_size/2 + ydisplace(j)) / figure_size;

        % Draw the arrow
        annotation(fig, 'arrow', [start_x, end_x], [start_y, end_y], 'HeadStyle', 'cback1', 'LineWidth', BFvals(j));
    end
end


%% Figure for fifth MS-C

% Define image file paths
MSA = 'path_to_microstate_prototypes/MS4.jpg';
MSB = 'path_to_microstate_prototypes/MS3.jpg';
MSC = 'path_to_microstate_prototypes/MS2.jpg';
MSD = 'path_to_microstate_prototypes/MS1.jpg';
MSE = 'path_to_microstate_prototypes/MS5.jpg';

image_paths = {MSB, MSA, MSE, MSD, MSC};
image_size = [75, 100]; % Define a uniform image size [width, height]

% Load and resize images
num_images = length(image_paths);
images = cell(1, length(image_paths));
for i = 1:length(image_paths)
    images{i} = imresize(imread(image_paths{i}), image_size);
end

% Define the circle parameters
circle_radius = 2 * max(image_size); % Adjust the circle radius
theta = linspace(0, 2*pi, num_images + 1); % Divide the circle into equal parts
theta(end) = []; % Remove the last element to avoid a duplicate of the first image

% Define the figure size dynamically based on the circle radius and image size
padding = max(image_size); % Space around the circle
figure_size = 2 * (circle_radius + padding);
fig = figure('Color', 'white', 'Position', [100, 100, figure_size, figure_size]);
hold on;
axis off

% Calculate positions and display images
for i = 1:num_images

    % Calculate image (x,y) positions
    x = circle_radius * cos(theta(i));
    y = circle_radius * sin(theta(i));

    % Calculate the normalized position
    x_norm = (x + figure_size/2 - image_size(1)/2) / figure_size;
    y_norm = (y + figure_size/2 - image_size(2)/2) / figure_size;

    % Create an axes for the image
    ax = axes('Position', [x_norm, y_norm, image_size(1)/figure_size, image_size(2)/figure_size]);

    % Display the image
    imshow(images{i}, 'Parent', ax);

    % Remove axis ticks
    set(ax, 'XTick', [], 'YTick', []);
end

axis equal;
axis off;

% Calculate the normalized start position MS-5: MS-C
ms = 5;
xdisplace = [-20 0 20 30 -25];
ydisplace = [-25 -25 -20 -5 10];
start_x = ((circle_radius * cos(theta(ms)) + figure_size/2) / figure_size);
start_y = ((circle_radius * sin(theta(ms)) + figure_size/2 + 25) / figure_size);

BFvals = [27.32 295.80 24.27 .13 0];
BFvals = 1 + BFvals/sum(BFvals)*3;

for j = 1:num_images
    if j ~= ms
        % Calculate the normalized end position for the arrows (center of the target image)
        end_x = (circle_radius * cos(theta(j)) + figure_size/2 + xdisplace(j)) / figure_size;
        end_y = (circle_radius * sin(theta(j)) + figure_size/2 + ydisplace(j)) / figure_size;

        % Draw the arrow
        annotation(fig, 'arrow', [start_x, end_x], [start_y, end_y], 'HeadStyle', 'cback1', 'LineWidth', BFvals(j));
    end
end

%% For all to all lines

% % Now draw lines
% % After placing the images, draw arrows between them
% for i = 1:num_images
%     % Calculate the normalized start position for the arrows (center of the image)
%     start_x = ((circle_radius * cos(theta(i)) + figure_size/2) / figure_size);
%     start_y = ((circle_radius * sin(theta(i)) + figure_size/2) / figure_size);
%
%     for j = 1:num_images
%         if i ~= j
%             % Calculate the normalized end position for the arrows (center of the target image)
%             end_x = (circle_radius * cos(theta(j)) + figure_size/2) / figure_size;
%             end_y = (circle_radius * sin(theta(j)) + figure_size/2) / figure_size;
%
%             % Draw the arrow
%             annotation(fig, 'arrow', [start_x, end_x], [start_y, end_y], 'HeadStyle', 'cback1', 'LineWidth', 1);
%         end
%     end
% end
