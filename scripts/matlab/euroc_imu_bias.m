sequences = ["V1_01_easy"; "V1_02_medium"; "V1_03_difficult"; "V2_01_easy"; "V2_02_medium"; "V2_03_difficult"; "MH_01_easy"; "MH_02_easy"; "MH_03_medium"; "MH_04_difficult"; "MH_05_difficult"];
N = size(sequences, 1);

gyro_bias = Inf;
acc_bias = Inf;
for idx = 1:N
    s = sequences(idx);
    state_data = csvread('data/euroc/' + s + '/state_groundtruth_estimate0/data.csv', 1, 0);
    
    % timestamp, p_RS_R_x [m], p_RS_R_y [m], p_RS_R_z [m], q_RS_w [], q_RS_x [], q_RS_y [], q_RS_z [], v_RS_R_x [m s^-1], v_RS_R_y [m s^-1], v_RS_R_z [m s^-1], b_w_RS_S_x [rad s^-1], b_w_RS_S_y [rad s^-1], b_w_RS_S_z [rad s^-1], b_a_RS_S_x [m s^-2], b_a_RS_S_y [m s^-2], b_a_RS_S_z [m s^-2]
    gyro_bias = min([gyro_bias, min(vecnorm(state_data(:, 12:14), 2, 2))]);
    acc_bias = min([acc_bias, min(vecnorm(state_data(:, 15:17), 2, 2))]);
end

disp('Min gyroscope bias magnitude:')
disp(gyro_bias)

disp('Min accelerometer bias magnitude:')
disp(acc_bias)
