sequences = ["V1_01_easy"; "V1_02_medium"; "V1_03_difficult"; "V2_01_easy"; "V2_02_medium"; "V2_03_difficult"; "MH_01_easy"; "MH_02_easy"; "MH_03_medium"; "MH_04_difficult"; "MH_05_difficult"];

N = size(sequences, 1);
time = zeros(N, 5);
scale_errors = zeros(N, 5);
gyro_errors = zeros(N, 5);
gyro_errors2 = zeros(N, 5);
acc_errors = zeros(N, 5);
acc_errors2 = zeros(N, 5);
gravity_errors = zeros(N, 5);

for idx = 1:N
    s = sequences(idx);
    
    data_iterative = csvread("data/experiment02/" + s + "_iterative_noprior.csv");
    time(idx, 1) = mean(data_iterative(:, 1))*1e-6;
    scale_errors(idx, 1) = mean(data_iterative(:, 2));
    gyro_errors(idx, 1) = mean(data_iterative(:, 3));
    gyro_errors2(idx, 1) = mean(data_iterative(:, 4));
    acc_errors(idx, 1) = mean(data_iterative(:, 5));
    acc_errors2(idx, 1) = mean(data_iterative(:, 6));
    gravity_errors(idx, 1) = mean(data_iterative(:, 7));
    
    data_iterative = csvread("data/experiment02/" + s + "_iterative.csv");
    time(idx, 2) = mean(data_iterative(:, 1))*1e-6;
    scale_errors(idx, 2) = mean(data_iterative(:, 2));
    gyro_errors(idx, 2) = mean(data_iterative(:, 3));
    gyro_errors2(idx, 2) = mean(data_iterative(:, 4));
    acc_errors(idx, 2) = mean(data_iterative(:, 5));
    acc_errors2(idx, 2) = mean(data_iterative(:, 6));
    gravity_errors(idx, 2) = mean(data_iterative(:, 7));
    
    data_mqh = csvread("data/experiment02/" + s + "_mqh.csv");
    time(idx, 3) = mean(data_mqh(:, 1))*1e-6;
    scale_errors(idx, 3) = mean(data_mqh(:, 2));
    gyro_errors(idx, 3) = mean(data_mqh(:, 3));
    gyro_errors2(idx, 3) = mean(data_mqh(:, 4));
    acc_errors(idx, 3) = mean(data_mqh(:, 5));
    acc_errors2(idx, 3) = mean(data_mqh(:, 6));
    gravity_errors(idx, 3) = mean(data_mqh(:, 7));
    
    data_ours = csvread("data/experiment02/" + s + "_ours_noprior.csv");
    time(idx, 4) = mean(data_ours(:, 1))*1e-6;
    scale_errors(idx, 4) = mean(data_ours(:, 2));
    gyro_errors(idx, 4) = mean(data_ours(:, 3));
    gyro_errors2(idx, 4) = mean(data_ours(:, 4));
    acc_errors(idx, 4) = mean(data_ours(:, 5));
    acc_errors2(idx, 4) = mean(data_ours(:, 6));
    gravity_errors(idx, 4) = mean(data_ours(:, 7));
    
    data_ours = csvread("data/experiment02/" + s + "_ours.csv");
    time(idx, 5) = mean(data_ours(:, 1))*1e-6;
    scale_errors(idx, 5) = mean(data_ours(:, 2));
    gyro_errors(idx, 5) = mean(data_ours(:, 3));
    gyro_errors2(idx, 5) = mean(data_ours(:, 4));
    acc_errors(idx, 5) = mean(data_ours(:, 5));
    acc_errors2(idx, 5) = mean(data_ours(:, 6));
    gravity_errors(idx, 5) = mean(data_ours(:, 7));
end

disp('Iterative w/o prior, Iterative, MQH, Ours w/o prior, Ours')

disp("Time (ms):")
disp(time)
disp("Mean:")
disp(mean(time))

disp("Scale (%):")
disp(scale_errors)
disp("Mean:")
disp(mean(scale_errors))

disp("Gyro bias (%):")
disp(gyro_errors)
disp("Mean:")
disp(mean(gyro_errors))

disp("Gyro bias (º):")
disp(gyro_errors2)
disp("Mean:")
disp(mean(gyro_errors2))

disp("Acc bias (%):")
disp(acc_errors)
disp("Mean:")
disp(mean(acc_errors))

disp("Acc bias (º):")
disp(acc_errors2)
disp("Mean:")
disp(mean(acc_errors2))

disp("Gravity (º):")
disp(gravity_errors)
disp("Mean:")
disp(mean(gravity_errors))
