function [] = Mal(a, b, n)
% MAL Реализация метода Малышева для вычисления одномерных спектральных портретов.
% 
% Входные параметры:
%   a - Начальное значение параметра λ
%   b - Конечное значение параметра λ
%   n - Размерность матрицы A
%


    % --- Инициализация параметров ---
    r1_val = a;         % Текущее значение радиуса (λ)
    h      = 0.001;     % Шаг по параметру λ
    i      = 1;         % Счетчик итераций
    w_max  = 1e14;      % Ограничение нормы решения
    u_max  = 1e10;      % Ограничение числа обусловленности
    ip     = 1e-15;     % Точность вычислений
    
    A      = randn(n, n); % Случайная матрица
    I      = eye(n);      % Единичная матрица
    
     % Расчет количества итераций для уточнения
    m0 = round(log2(- (1 + w_max) * log(ip / ((2 + 2 * ip) * sqrt(w_max)))));


    % --- Основной цикл по параметру r1_val---
    while (r1_val < b)
        skip_step = false; 
        A0 = A' / r1_val;
        B0 = I;
        p  = 2 * I;
        d  = 0;
        
        % Внутренний цикл уточнения проектора
        while (max(abs(p * p - p), [], 'all') > ip && d <= m0 / 3)
            d = d + 1;
            for i1 = 1:3
                S_temp = [-B0; A0];
                [Q_mat, ~] = qr(S_temp);
                QQ     = Q_mat';
                Q1_sub = QQ(n+1:2*n, 1:n);
                Q2_sub = QQ(n+1:2*n, n+1:2*n);
                A0     = Q1_sub * A0;
                B0     = Q2_sub * B0;
            end

            % Проверка обусловленности
            if (max(cond(A0 - B0), cond(A0 + B0)) > u_max)
                skip_step = true;
                break; 
            end
            
            p = -inv(A0 - B0) * B0; 
        end
        
        % Обработка ошибок сходимости
        if skip_step || isempty(p) || any(isnan(p(:)))
            e(i) = r1_val;
            m(i) =NaN;
            i    = i + 1;
            r1_val = r1_val + h;
            fprintf('Внимание: не удалось получить проектор при r = %f\n', r1_val);
            continue 
        end
        
        % Вычисление решения
        inv_su = inv(B0 + A0);
        H     = inv_su * inv_su';
        m(i)  = max(abs(H), [], 'all');

        e(i)  = r1_val;
        i      = i + 1;
        r1_val = r1_val + h;
    end
    
    % --- Визуализация результатов ---
    figure;
    plot(e, log10(m), 'k', 'LineWidth', 2);
    xlabel('r', 'FontSize', 14, 'FontWeight', 'bold');
    ylabel('log_{10}||H||', 'FontSize', 14, 'FontWeight', 'bold');
    set(gca, 'FontSize', 15);
    
   
