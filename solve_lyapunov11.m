function [D1, D2, D3, dd1, current_mm1, dd2, current_mm2] = solve_lyapunov11(a, b, A_original)
%
% Входные параметры:
%   a          - Начальное значение параметра лямбда
%   b          - Конечное значение параметра лямбда
%   A_original - Исходная матрица системы
%

    % --- Константы и параметры ---
    h      = 0.01;   % Шаг по параметру
    w_max  = 1e20;   % Ограничение нормы решения
    u_max  = 1e20;   % Ограничение обусловленности
    ip     = 1e-15;  % Точность вычислений
    
    [n, ~] = size(A_original);
    I      = eye(n);
    j      = 1;
    i      = 1;
    
    % Инициализация массивов
    dd1 = []; dd2 = []; current_mm1 = []; current_mm2 = [];
    e1 = []; m1 = []; mm = [];
    
    inside_zone = false;
    m0 = round(log2(- (1 + w_max) * log(ip / ((2 + 2 * ip) * sqrt(w_max)))));
    
    r1_val = a; 

    % --- Основной цикл ---
    while (r1_val < b)
        skip_step = false; 
        A0 = A_original' / r1_val;
        B0 = I;
        p_mat = 2 * I;
        d_count = 0;
        
        % Внутренний цикл уточнения проектора
        while (max(abs(p_mat * p_mat - p_mat), [], 'all') > ip * max(abs(p_mat), [], 'all') && d_count <= m0)
            d_count = d_count + 1;
            for i1 = 1:3
                S_temp = [-1 * B0; A0];
                [Q_mat, ~] = qr(S_temp);
                QQ = Q_mat';
                Q1_sub = QQ(n+1:2*n, 1:n);
                Q2_sub = QQ(n+1:2*n, n+1:2*n);
                A0 = Q1_sub * A0;
                B0 = Q2_sub * B0;
            end

            if (max(cond(A0 - B0), cond(A0 + B0)) > u_max)
                skip_step = true;
                break; 
            end
            p_mat = -inv(A0 - B0) * B0; 
        end
        
        % Проверка корректности проектора
        if skip_step || isempty(p_mat) || any(isnan(p_mat(:)))
            e1(i) = r1_val;
            m1(i) = NaN; 
            r1_val = r1_val + h;
            i = i + 1;
            fprintf('Ошибка: собственное значение близко к окружности r = %f\n', r1_val);
            continue 
        end
        
        % Вычисление H1 и параметров спектра
        inv_su = inv(B0 + A0);
        H1 = inv_su * inv_su';
        m1(i) = max(abs(H1), [], 'all');
        nor_A = max(abs(A_original), [], 'all');

        
       % Оценка параметра m(A), характеризующего близость к неустойчивости
        H_mu  = m1(i) * r1_val / (2 * nor_A) + sqrt((m1(i) * r1_val)^2 / nor_A^2 + 4 * m1(i));
        mm(i) = 1 / H_mu;

        % --- Алгоритм идентификации зон спектральной неустойчивости ---
        current_mm = mm(i); 
        if current_mm > 1e-12
            jump = h; 
        elseif current_mm <= 1e-12 && current_mm > ip
            if inside_zone
                dd2(j) = r1_val;
                current_mm2(j) = mm(i);
                j = j + 1;
                inside_zone = false; 
            end
            jump = 1.01 * 1e-4;
        else
            jump = 1.01 * 1e-4;
            if ~inside_zone
                dd1(j) = e1(i-1);
                current_mm1(j) = mm(i-1);
                inside_zone = true;
            end
        end

        e1(i) = r1_val;
        i = i + 1;
        r1_val = r1_val + h;
    end

    % --- Возврат результатов ---
    D1 = e1;
    D2 = m1;
    D3 = mm;
