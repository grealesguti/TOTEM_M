% Generate synthetic data
rng(1); % Set random seed for reproducibility
data = randn(100, 3); % 100 samples with 3 features

% Perform PCA
[coeff, score, latent, ~, explained] = pca(data);

% Display results
disp('Principal components (coefficients):');
disp(coeff);

disp('Principal component scores:');
disp(score);

disp('Variance explained by each principal component:');
disp(explained);

% Plot the cumulative explained variance
figure;
plot(cumsum(explained), 'o-');
xlabel('Number of Principal Components');
ylabel('Cumulative Explained Variance');
title('Cumulative Explained Variance');

% Select the number of components based on a desired explained variance threshold
desiredExplainedVariance = 95; % Set your desired threshold
numComponents = find(cumsum(explained) >= desiredExplainedVariance, 1);

disp(['Number of principal components to achieve ' num2str(desiredExplainedVariance) '% explained variance: ' num2str(numComponents)]);

% Use the selected number of components for dimensionality reduction
reducedData = score(:, 1:numComponents);

% You can now use 'reducedData' for further analysis or visualization
