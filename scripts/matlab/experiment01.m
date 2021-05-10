sequences = ["V1_01_easy"; "V1_02_medium"; "V1_03_difficult"; "V2_01_easy"; "V2_02_medium"; "V2_03_difficult"; "MH_01_easy"; "MH_02_easy"; "MH_03_medium"; "MH_04_difficult"; "MH_05_difficult"];

k_values = [5 10 20 50 75];

N = size(sequences, 1);
time = zeros(N, 5);
scale_errors = zeros(N, 5);
gyro_bias_errors = zeros(N, 5);
gyro_bias_errors2 = zeros(N, 5);
acc_bias_errors = zeros(N, 5);
acc_bias_errors2 = zeros(N, 5);
gravity_errors = zeros(N, 5);

disp('Iterative w/o prior, Iterative, MQH, Ours w/o prior, Ours')
for k_index = 1:length(k_values)
    k = k_values(k_index);
    for idx = 1:N
        s = sequences(idx);
        
        data_iterative_noprior = csvread("data/experiment01/" + s + "_" + int2str(k) + "_iterative_noprior.csv");
        scale_errors(idx, 1) = mean(data_iterative_noprior(:, 1));
        gyro_bias_errors(idx, 1) = mean(data_iterative_noprior(:, 2));
        gyro_bias_errors2(idx, 1) = mean(data_iterative_noprior(:, 3));
        acc_bias_errors(idx, 1) = mean(data_iterative_noprior(:, 4));
        acc_bias_errors2(idx, 1) = mean(data_iterative_noprior(:, 5));
        gravity_errors(idx, 1) = mean(data_iterative_noprior(:, 6));
        
        data_iterative = csvread("data/experiment01/" + s + "_" + int2str(k) + "_iterative.csv");
        scale_errors(idx, 2) = mean(data_iterative(:, 1));
        gyro_bias_errors(idx, 2) = mean(data_iterative(:, 2));
        gyro_bias_errors2(idx, 2) = mean(data_iterative(:, 3));
        acc_bias_errors(idx, 2) = mean(data_iterative(:, 4));
        acc_bias_errors2(idx, 2) = mean(data_iterative(:, 5));
        gravity_errors(idx, 2) = mean(data_iterative(:, 6));
        
        data_mqh = csvread("data/experiment01/" + s + "_" + int2str(k) + "_mqh.csv");
        scale_errors(idx, 3) = mean(data_mqh(:, 1));
        gyro_bias_errors(idx, 3) = mean(data_mqh(:, 2));
        gyro_bias_errors2(idx, 3) = mean(data_mqh(:, 3));
        acc_bias_errors(idx, 3) = mean(data_mqh(:, 4));
        acc_bias_errors2(idx, 3) = mean(data_mqh(:, 5));
        gravity_errors(idx, 3) = mean(data_mqh(:, 6));
        
        data_ours_noprior = csvread("data/experiment01/" + s + "_" + int2str(k) + "_ours_noprior.csv");
        scale_errors(idx, 4) = mean(data_ours_noprior(:, 1));
        gyro_bias_errors(idx, 4) = mean(data_ours_noprior(:, 2));
        gyro_bias_errors2(idx, 4) = mean(data_ours_noprior(:, 3));
        acc_bias_errors(idx, 4) = mean(data_ours_noprior(:, 4));
        acc_bias_errors2(idx, 4) = mean(data_ours_noprior(:, 5));
        gravity_errors(idx, 4) = mean(data_ours_noprior(:, 6));
        
        data_ours = csvread("data/experiment01/" + s + "_" + int2str(k) + "_ours.csv");
        scale_errors(idx, 5) = mean(data_ours(:, 1));
        gyro_bias_errors(idx, 5) = mean(data_ours(:, 2));
        gyro_bias_errors2(idx, 5) = mean(data_ours(:, 3));
        acc_bias_errors(idx, 5) = mean(data_ours(:, 4));
        acc_bias_errors2(idx, 5) = mean(data_ours(:, 5));
        gravity_errors(idx, 5) = mean(data_ours(:, 6));
    end

    disp("K:")
    disp(k)
    
    disp("Scale (%):")
    disp(scale_errors)
    disp("Mean:")
    disp(mean(scale_errors))

    disp("Gyro bias (%):")
    disp(gyro_bias_errors)
    disp("Mean:")
    disp(mean(gyro_bias_errors))
    
    disp("Gyro bias (º):")
    disp(gyro_bias_errors2)
    disp("Mean:")
    disp(mean(gyro_bias_errors2))

    disp("Acc bias (%):")
    disp(acc_bias_errors)
    disp("Mean:")
    disp(mean(acc_bias_errors))
    
    disp("Acc bias (º):")
    disp(acc_bias_errors2)
    disp("Mean:")
    disp(mean(acc_bias_errors2))

    disp("Gravity (º):")
    disp(gravity_errors)
    disp("Mean:")
    disp(mean(gravity_errors))
end
