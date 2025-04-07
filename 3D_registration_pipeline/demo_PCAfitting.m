% Generate synthetic 3D data
rng(42); % For reproducibility
num_points = 100;

% Generate points lying near a tilted plane: z = 0.5x - 0.3y + 2
X = rand(num_points, 1) * 10; % Random x-coordinates
Y = rand(num_points, 1) * 10; % Random y-coordinates
Z = 0.5 * X - 0.3 * Y + 2 + randn(num_points, 1) * 0.5; % Add noise to z-coordinates

% Combine into a single matrix
centroids = [X, Y, Z];

%% Perform PCA
[coeff, ~, ~] = pca(centroids); 

% Extract the normal vector (3rd principal component)
normal_vector = coeff(:, 3);

% Find the mean of the points
mean_point = mean(centroids, 1);

% Visualize the data and the fitted plane
figure;
scatter3(centroids(:, 1), centroids(:, 2), centroids(:, 3), 'b', 'filled'); 
hold on;

% Define the plane based on the PCA normal vector
[X_plane, Y_plane] = meshgrid(linspace(min(X), max(X), 10), linspace(min(Y), max(Y), 10));
Z_plane = mean_point(3) - ...
    normal_vector(1)/normal_vector(3) * (X_plane - mean_point(1)) - ...
    normal_vector(2)/normal_vector(3) * (Y_plane - mean_point(2));

% Plot the plane
surf(X_plane, Y_plane, Z_plane, 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'FaceColor', 'r');
colormap autumn;

% Annotate the plot
xlabel('X'); ylabel('Y'); zlabel('Z');
title('3D Points and PCA Fitted Plane');
grid on; axis equal;
legend('Data points', 'Fitted plane');
view(3);
hold off;

% Display the normal vector and its angles
disp('Normal vector of the plane:');
disp(normal_vector);

% Calculate angles with coordinate planes
angle_with_xy = acosd(abs(normal_vector(3)) / norm(normal_vector));
angle_with_xz = acosd(abs(normal_vector(2)) / norm(normal_vector));
angle_with_yz = acosd(abs(normal_vector(1)) / norm(normal_vector));

fprintf('Angle with XY-plane: %.2f degrees\n', angle_with_xy);
fprintf('Angle with XZ-plane: %.2f degrees\n', angle_with_xz);
fprintf('Angle with YZ-plane: %.2f degrees\n', angle_with_yz);

%% 
% Center the data
centered_data = centroids - mean_point;

% PCA transformation matrix (columns are new axes)
T = coeff; % Transformation matrix based on PCA components

% Transform the data into the new coordinate system
transformed_data = centered_data * T;

% Extract the coordinates in the new XY' plane
XY_plane_data = transformed_data(:, 1:2); % Only X' and Y' components

% Calculate distances in the new XY' plane
% Example: Distance from the first point to all others
distances = sqrt(sum((XY_plane_data - XY_plane_data(1, :)).^2, 2));

% Visualize the points in the new coordinate system
figure;
scatter(XY_plane_data(:, 1), XY_plane_data(:, 2), 'b', 'filled');
xlabel('X'''); ylabel('Y''');
title('Points in the Transformed XY Plane');
grid on; axis equal;

