sequences = ["V1_01_easy"; "V1_02_medium"; "V1_03_difficult"; "V2_01_easy"; "V2_02_medium"; "V2_03_difficult"; "MH_01_easy"; "MH_02_easy"; "MH_03_medium"; "MH_04_difficult"; "MH_05_difficult"];

config = aconfig;
config.Colormap = 'gray';

afigure(1, config);
hold on

frames = [20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240]';
hz = 4;

X = frames./hz;
Y_iterative = zeros(size(frames, 1), 1);
Y_mqh = zeros(size(frames, 1), 1);
Y_ours = zeros(size(frames, 1), 1);
velocities = zeros(size(frames, 1), 1);
preintegration = zeros(size(frames, 1), 1);
for i = 1:size(frames, 1)
    k = frames(i);
    sum_times_iterative = 0;
    sum_times_mqh = 0;
    sum_times_ours = 0;
    sum_vel_times_ours = 0;
    sum_int_times = 0;
    N = 0;
    for j = 1:size(sequences, 1)
        s = sequences(j);
        
        data_iterative = csvread("data/experiment01c/" + s + "_" + int2str(k) + "_iterative.csv");
        sum_times_iterative = sum_times_iterative + 1e-6*median(data_iterative(:, 1));
        
        data_mqh = csvread("data/experiment01c/" + s + "_" + int2str(k) + "_mqh.csv");
        sum_times_mqh = sum_times_mqh + 1e-6*median(data_mqh(:, 1)) + 1e-6*median(data_mqh(:, 2));
        
        data_ours = csvread("data/experiment01c/" + s + "_" + int2str(k) + "_ours.csv");
        sum_times_ours = sum_times_ours + 1e-6*median(data_ours(:, 1)) + 1e-6*median(data_ours(:, 2));
        sum_vel_times_ours = sum_vel_times_ours + 1e-6*median(data_ours(:, 2));
        
        % Preintegration time
        sum_int_times = sum_int_times + 1e-6*median(data_ours(:, 3));
        
        N = N + 1;
    end
    Y_iterative(i) = sum_times_iterative / N;
    Y_mqh(i) = sum_times_mqh / N;
    Y_ours(i) = sum_times_ours / N;
    velocities(i) = sum_vel_times_ours / N;
    preintegration(i) = sum_int_times / N;
end
plot(X, Y_iterative, 'LineStyle', ':', 'LineWidth', 3); %, 'Color', [0 0 0])
plot(X, Y_mqh, 'LineStyle', '--', 'LineWidth', 3); %, 'Color', [1 0.25 0.15])
plot(X, Y_ours, 'LineStyle', '-', 'LineWidth', 3); %, 'Color', [0.9 0.75 0.1])

xlim([5, 60])

xlabel('Initialization time [s]')
ylabel('Solving time [ms]')
legend('Iterative inertial-only [3]', 'Fast initialization [2]', 'Analytical solution', 'Location', 'northwest')

disp('Velocity estimation for 20 frames [ms]:')
disp(velocities(1))

disp('Pre-integration for 20 frames [ms]:')
disp(preintegration(1))
