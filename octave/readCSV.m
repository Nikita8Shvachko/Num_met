% Функция для чтения данных из CSV-файла
function [x, y] = readCSV(filename)
    data = csvread(filename);
    x = data(:, 1);
    y = data(:, 2);
end

