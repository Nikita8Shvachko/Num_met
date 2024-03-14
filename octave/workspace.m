
% Путь к папке с файлами (относительный путь)
folder_path = 'cmake-build-debug/';

% Формирование путей к файлам для функции 1
files_function_1 = {'function_1_data_3.csv', 'function_1_data_4.csv', ...
    'function_1_data_5.csv', 'function_1_data_6.csv', ...
    'function_1_data_7.csv', 'function_1_data_8.csv', ...
    'function_1_data_9.csv', 'function_1_data_10.csv'};

% Формирование путей к файлам для функции 2
files_function_2 = {'function_2_data_3.csv', 'function_2_data_4.csv', ...
    'function_2_data_5.csv', 'function_2_data_6.csv', ...
    'function_2_data_7.csv', 'function_2_data_8.csv', ...
    'function_2_data_9.csv', 'function_2_data_10.csv'};


% Чтение данных для функции 1
for i = 1:length(files_function_1)
    [x, y] = readCSV(fullfile(folder_path, files_function_1{i}));
    eval(['x1_' num2str(i) ' = x;']);
    eval(['y1_' num2str(i) ' = y;']);
end

% Чтение данных для функции 2
for i = 1:length(files_function_2)
    [x, y] = readCSV(fullfile(folder_path, files_function_2{i}));
    eval(['x2_' num2str(i) ' = x;']);
    eval(['y2_' num2str(i) ' = y;']);
end

% Построение графиков
figure;

% Графики для функции 1
subplot(2, 1, 1); % Создание верхнего графика
for i = 1:length(files_function_1)
    eval(['plot(x1_' num2str(i) ', y1_' num2str(i) ', ''Color'', rand(1, 3), ''LineWidth'', 1.5);']);
    hold on;
end
title('Function 1');
xlabel('x');
ylabel('y');
legend('3 nodes', '4 nodes', '5 nodes', '6 nodes', '7 nodes', '8 nodes', '9 nodes', '10 nodes');
grid on;

% Графики для функции 2
subplot(2, 1, 2); % Создание нижнего графика
for i = 1:length(files_function_2)
    eval(['plot(x2_' num2str(i) ', y2_' num2str(i) ', ''Color'', rand(1, 3), ''LineWidth'', 1.5);']);
    hold on;
end
title('Function 2');
xlabel('x');
ylabel('y');
legend('3 nodes', '4 nodes', '5 nodes', '6 nodes', '7 nodes', '8 nodes', '9 nodes', '10 nodes');
grid on;

