numSpins = 100;        
numSamples = 1000;      
numEpochs = 500;         
learningRate = 0.01;

data = sign(randn(numSamples, numSpins));
data(data == 0) = 1;

J = 0.01 * randn(numSpins);    
J = triu(J,1);                 
J = J + J';                  
h = zeros(1, numSpins);     

% Contrastive Divergence 
for epoch = 1:numEpochs
    % Positive phase
    pos_corr = (data' * data) / numSamples;

    % Negative phase
    model_samples = data;
    for i = 1:numSpins
        local_field = model_samples * J(:,i) + h(i);
        prob = 1 ./ (1 + exp(-2 * local_field));  % P(s=+1)
        model_samples(:,i) = sign(prob - rand(numSamples,1));
        model_samples(model_samples(:,i)==0, i) = 1;  % map 0 to +1
    end

    % Correlations
    neg_corr = (model_samples' * model_samples) / numSamples;

    dJ = learningRate * (pos_corr - neg_corr);
    dJ = triu(dJ,1);
    J = J + dJ + dJ';
    dh = learningRate * (mean(data,1) - mean(model_samples,1));
    h = h + dh;

    energy = -sum(sum((data * J) .* data)) / numSamples;
    fprintf('Epoch %d: Avg Energy = %.4f\n', epoch, energy);
end