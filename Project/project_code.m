Eiffel = im2double(rgb2gray(imread('EiffelTower.jpg')));

[U,S,V] = svd(Eiffel);

% Extract singular values
sigmas = diag(S);

% Plot singular values on a log scale
figure;
plot(log10(sigmas));
title('Singular Values (Log10 Scale)');
xlabel('Index');
ylabel('log_{10}(\sigma_i)');

% Plot cumulative percent of singular values
figure;
plot(cumsum(sigmas)/sum(sigmas));
title('Cumulative Percent of Total Sigmas');
xlabel('Number of Singular Values');
ylabel('Cumulative Percentage');

% Determine approximately how many singular values
% retain 75% of the information
cum_percent = cumsum(sigmas)/sum(sigmas);
rank75 = find(cum_percent >= 0.75,1);

fprintf('Approximate rank for 75%% information: %d\n', rank75);

% Display original image
figure;
imshow(Eiffel);
title('Full-Rank Eiffel Tower');

% Low-rank approximations
ranks = [200,150,100,50,30,20,10];

for i = 1:length(ranks)

    r = ranks(i);

    approx_sigmas = sigmas;
    approx_sigmas(r+1:end) = 0;

    ns = length(sigmas);

    approx_S = S;
    approx_S(1:ns,1:ns) = diag(approx_sigmas);

    approx_Eiffel = U * approx_S * V';

    figure;
    imshow(approx_Eiffel);
    title(sprintf('Rank %d Eiffel Tower', r));

end