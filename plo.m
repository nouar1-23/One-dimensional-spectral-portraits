function [] = plo()
% Функция вычисляет спектральные портреты для трех матриц,

    [A0, A1, A3] = matre_x();
    A = {A0, A1, A3};

    % --- Вычисление спектральных портретов ---
    for i = 1:3
        [F1{i}, F2{i}, F3{i}, dd1{i}, current_mm1{i}, dd2{i}, current_mm2{i}] = solve_lyapunov11(0.1, 3, A{i});
    end

    % --- Построение графиков спектральных портретов ---
    figure;
    colors = {'b', 'm', 'g'};
    for i = 1:3
        plot(F1{i}, log10(F2{i}), [colors{i}, '-'], 'LineWidth', 2); hold on;
        plot(F1{i}, log10(F3{i}), [colors{i}, ':'], 'LineWidth', 2);
    end

    % Линия порога точности
    x_limits = xlim; 
    plot(x_limits, [-15, -15], 'r-', 'LineWidth', 2);

    % Настройка оси Y (цветовое оформление -15)
    current_ticks = get(gca, 'ytick');
    new_ticks = unique(sort([current_ticks, -15]));
    set(gca, 'ytick', new_ticks);
    
    tick_labels = cell(size(new_ticks));
    for i = 1:length(new_ticks)
        if new_ticks(i) == -15
            tick_labels{i} = '\color{red}-15';
        else
            tick_labels{i} = num2str(new_ticks(i));
        end
    end
    set(gca, 'yticklabel', tick_labels);

    xlabel('r');
    set(gca, 'FontSize', 25);

    % Легенда
    legend({'$||H_{A_{0}}(r)||$', '$\delta_{A_{0}}(r)$', ...
            '$||H_{A_{1}}(r)||$', '$\delta_{A_{1}}(r)$', ...
            '$||H_{A_{3}}(r)||$', '$\delta_{A_{3}}(r)$'}, ...
            'Interpreter', 'latex', 'FontSize', 16, 'Location', 'northeast');

    