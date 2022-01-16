
rng(0);
impulse_response = [1 1 1 1; 1 1 0 1; 1 0 0 0];
% message = randi([0 1], 1, 20);
message = [1 1 1 0 0];
encoded_data = convolutional_enc(message, impulse_response);

message2 = [1 0 0 0 1 0 1 1 0 1 0 1 1 0 0 0 1 0 1 0 1 0 0 0];
decoded_data = convolutional_dec(message2, impulse_response);

function encoded_data = convolutional_enc(binary_data, impulse_response)
    [row, col] = size(impulse_response);
    len = length(binary_data);
    for ii = 1 : row
        encoded_data(ii,:) = mod(conv(binary_data, impulse_response(ii,:)), 2);
    end
    encoded_data = reshape(encoded_data,1,[]);
end

function decoded_data = convolutional_dec(binary_data, impulse_response)
    [row, col] = size(impulse_response);
    num_state = 2^(col-1);
    states = de2bi([0:num_state-1], col-1, 'left-msb');
    Cu = zeros(num_state, row);
    Cd = zeros(num_state, row);
    for ii = 1 : row
        Cu(:,ii) = mod(sum([states zeros(num_state,1)].*impulse_response(ii,:),2),2);
        Cd(:,ii) = mod(sum([states ones(num_state,1)].*impulse_response(ii,:),2),2);
    end

    num_level = length(binary_data) / row;
    cost = inf(num_state, 1);
    cost(1) = 0;
    binary_data = reshape(binary_data,row,[])';
    idu = mod([0:num_state-1]'*2, num_state) + 1;
    idd = mod([0:num_state-1]'*2+1, num_state) + 1;
    previous = zeros(num_state, num_level);
    for ii = 1 : num_level
        costu = cost(idu) + sum(Cu ~= binary_data(ii,:), 2);
        costd = cost(idd) + sum(Cd ~= binary_data(ii,:), 2);
        [cost, idx] = min([costu costd], [], 2);
        idx = idx - 1;
        previous(:, ii) = idu.*(1 - idx) + idd.*(idx);
    end

    place = 1;
    decoded_data(num_level) = 0;
    for ii = num_level : -1 : 2
        place = previous(place, ii);
        decoded_data(ii-1) = (place > num_state/2);
    end
    decoded_data = decoded_data(1:(end-row));
end