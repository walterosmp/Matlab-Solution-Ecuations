format long
%               Исходние матрицы
A=[ 2.0000    1.0000    0.6667    0.5000    0.4000 ;
    1.0000    0.6621    0.5000    0.4000    0.3321 ;
    0.6621    0.5000    0.4021    0.3335    0.2857 ;
    0.5000    0.4000    0.3335    0.2857    0.2523 ;
    0.4000    0.3321    0.2857    0.2523    0.2218 ];

b = [ 3.7680
    -0.1333
    2.4980
    -11.1570
    5.2990];
%   Создаем возмущенную матрицу A_perturbed
A_p = A;
A_p(4, 1) = A_p(4, 1) + 0.01;
A_p(1, 4) = A_p(1, 4) + 0.01;


%   Cначала создаем функции для базовые операции 

%   Создаем функцию для нормализации вектора

function normValue = norm(A)
    res = 0;
    % Определяем количество элементов в A
    [m, n] = size(A);  % m - количество строк, n - количество столбцов

    % Суммируем квадраты всех элементов матрицы A
    for i = 1:m
        for j = 1:n
            res = res + A(i,j)^2;
        end
    end

    % Вычисляем 2-норму
    normValue = abs(sqrt(res));
end

%   Создаем единичной матрицы

function I =myeye(N)
    I = zeros(N);
    for i=1: N
        I(i, i) = 1;
    end
end


%Создание диагональной матрицы 

function diagMat = diag1(A)
    % Длина вектора
    n = length(A);

    % Создаем матрицу нулей размером n x n
    diagMat = zeros(n);

    % Заполняем главную диагональ элементами вектора
    for i = 1:n
        diagMat(i, i) = A(i, i);
    end
end

%   Создаем транспонированную матрицу

function T = trans(A)
[row, col] = size(A);
    A_transposed = zeros(col, row);
    for i = 1:row
        for j = 1:col
            A_transposed(j, i) = A(i, j); % Перестановка индексов
        end
    end
    T=A_transposed;
end

%   Создаем функцию для умножения матриц

function M = mult(A, B)
    [rowsA, colsA] = size(A);
    [rowsB, colsB] = size(B);
    M = zeros(rowsA, colsB);
    for i = 1:rowsA
        for j = 1:colsB
            for k = 1:colsA
               M(i, j) = M(i, j) + A(i, k) * B(k, j);
            end
        end
    end
end

%   Метод гаусса

function X = methodG(A,b)
    N=length(b);
    %   Объединение матрицы A и столбца B
    Ab=[A b];
    %   Начинаем цикл начинаяя с 1 до N (5)
    for i = 1:N   
        Ab(i,:) = Ab(i,:)/Ab(i,i);   %Нормализация строки по ведущему элементу
        %каждое число разделяется на себя получая 1 
        for j = i+1:N   
            L = -Ab(j,i);   %Ведущий элемент
            Ab(j, :) = Ab(j, :) + L * Ab(i, :); 
            %Тут проходим по столбцам с j ого элемента, добавляем к нему
            %его противоположен получая нуль.
        end
    end
    X=Ab;
end

%   Обратный ход метода Гаусса 

function X = reverseG(Ab)
    N = 5;
    %   Обратный ход, снизу верх 
    X = zeros(N, 1); %Создаем матрица столбец нули
    X(N) = Ab(N,6) /Ab(N,N);
    for i = N-1:-1:1
        X(i) = Ab(i,6) - Ab(i, i+1:N)* X(i+1:N);   
        X(i) = X(i) / Ab(i, i);
        %{
        Ab(i,6)- это текущая правая чаcть уравнения
        Ab(i, i+1:N)-Это значения i ого строка, c предыдущей операцией (i+1)
        до первого (N)
        %}
    end
end

%   Функция разложения матрицы M с использованием отражений Хаусхолдера

function [Q, R] = hholder(M)

    % Выход: Q -ортогональная матрица, R - верхнетреугольная матрица
    % Размеры матрицы
    [rows, cols] = size(M); 
    % Инициализация Q как единичной матрицы
    Q = eye(rows);
    R = M;

    for i = 1: cols -1
        % Выделяем вектор x из текущего столбца
        x = R(i:end, i);
        % Создаем вектор e (единичный вектор)
        e = zeros(length(x), 1);
        e(1) = 1;
        % Норма вектора x
        norm_x = norm(x);
        % Вычисляем вектор u для отражения
        u = x + norm_x * e;
        u = u / norm(u); % Нормализация вектора u
        % Построение матрицы Хаусхолдера
        P = eye(rows);
        P(i:end, i:end) = eye(length(x))-2*mult(u,trans(u))/mult(trans(u),u);
        % Обновление Q и R

        Q = Q *trans(P);
        R = P * R;
    end
end

% Нахождения наибольшее собственное число и соответствующое ему собственное вектор

function [eigenvalue, eigenvector] = power_iteration(A, tolerance, max_iterations)
    v = randn(size(A, 1), 1);
    %{
    Случайный начальный вектор-столбец размерностю 5*1
    tolerance — задаёт критерий сходимости,
    насколько близко собственный вектор соответствует условию
    Av=lambda*v
    %} 
    for i = 1:max_iterations
        Av = A * v;     %5*5  5*1 ----> 5*1
        % Вычисление нормы

        eigenvalue = norm(Av);

        % Нормализация вектора
        v = Av/eigenvalue;

        %Проверка сходимости
        if norm(A*v-eigenvalue*v) < tolerance
            break;
        end
    end
    eigenvector = v;
end

% Нахождения наименьшее собственное число и соответствующое ему собственное вектор

function [lambda_min, v_min] = qr_method_min_vector(A, tolerance, max_iterations)
    % Инициализация
    n = size(A, 1);  % Размерность матрицы A
    A_k = A;  % Создаем копию А которую будем модифицировать на каждой итерации
    Q_accum = myeye(n);  % Единичная матрица для накопления Q
    i = 0;
    %Начинаем цикл, который будет выполняться до тех пор,
    % пока не достигнем максимального числа итераций (max_iterations),
    % или пока не сработает условие сходимости
    while i < max_iterations
        i = i+ 1;
        
        % Выполнение QR-разложения
        [Q, R] = hholder(A_k); 
        
        % Обновляем A_k как R * Q Это операция, которая постепенно приводит
        % A_k к диагональной матрице с собственными значениями на диагонали.
        A_k = R * Q;
        
        % Накопление преобразований Q
        Q_accum = Q_accum * Q;

        % Проверка сходимости
        % Вычисляем норму элементов матрицы A_k
        off_diagonal_norm = norm(A_k - diag1(A_k));
        if off_diagonal_norm < tolerance
            break;
        end
    end

    % Собственные числа находятся на диагонали матрицы A_k
    eigenvalues = diag(A_k);
    
    % Минимальное собственное число
    lambda_min = min(eigenvalues);
    
    % Индекс минимального собственного числа
    [~, idx_min] = min(eigenvalues);
    
    % Собственный вектор, соответствующий минимальному числу
    v_min = -Q_accum(:, idx_min);
end

%-------------Обработка данных---------------------------
disp("Решение методом Гаусса для исходной системы " );
disp("Оригинальная система' " );
disp(A);
disp("Имеется возмущенная система " );
disp(A_p);

disp("Решение методом Гаусса для исходной системы   " );
G=methodG(A,b);
G1=reverseG(G);
disp(G1);
disp("Решение методом Гаусса для возмущенной системы   " );
G_p=methodG(A_p,b);
G_p1=reverseG(G_p);
disp(G_p1);

%Невязки

e = mult(A,G1);
e_p = mult(A_p,G_p1);

%{
Определяем невязку исходя из того что AX=b,
мы получили вектор X, сравниваем AX с тем что 
должно было получится
%}
disp("Определяем невязку " );

res_err = norm (b-e);
res_err_p = norm (b-e_p);
disp(" для исходной системы:  " );
disp(res_err);
disp(" для возмущенной системы: " );
disp(res_err_p);


%{
 QR раздожение
A=QR---> QRx=b --->Rx=Q(t)b
    Для исходной систамы
%}


[Q, R] = hholder(A);
Qb=mult(trans(Q),b);
R_Qb=[R Qb];
disp("----- QR раздожением---------" );
disp("Решение для исходной системы:" );
%Решение будет:
XQR = reverseG(R_Qb);
disp(XQR);

% Невязка для решения методом QR-разложения
bqr=mult(A,XQR);


%   Для возмущенной системы  
[Qp, Rp] = hholder(A_p);
Qb_p=mult(trans(Qp),b);
R_Qb_p=[Rp Qb_p];

disp("Решение для возмущенной системы:" );
XQR_p = reverseG(R_Qb_p);
disp(XQR_p);
disp("Невязка для решения методом QR-разложения" );

bqr_p=mult(A_p,XQR_p);

res_err_qr = norm (b-bqr);
res_err_qr_p = norm (b-bqr_p);

disp(" для исходной системы:  " );
disp(res_err_qr );
disp(" для возмущенной системы: " );
disp(res_err_qr_p);


disp("Собственные числа и Собственные вектора " ); 
tolerance = 10e-10;
max_iterations = 100;
disp(" Максимальные собственные числа и соответствующие им собственные векторы " );
[maxeigenvalue, maxeigenvector ] = power_iteration(A, tolerance, max_iterations);

%   Для возмущенной системы  
[maxeigenvalue_p, maxeigenvector_p] = power_iteration(A_p, tolerance, max_iterations);

disp(" для исходной системы:  " );
disp(maxeigenvalue);
disp(maxeigenvector);
disp(" для возмущенной системы: " );
disp(maxeigenvalue_p);
disp(maxeigenvector_p );

disp(" Минимальное собственные числа и соответствующие им собственные векторы " );
tolerance = 10e-10;
max_iterations = 10000;

[mineigenvalue, mineigenvector ] = qr_method_min_vector(A, tolerance, max_iterations);

%   Для возмущенной системы  
[mineigenvalue_p, mineigenvector_p] = qr_method_min_vector(A_p, tolerance, max_iterations);

disp(" для исходной системы:  " );
disp(mineigenvalue); 
disp(mineigenvector);
disp(" для возмущенной системы: " );
disp(mineigenvalue_p);
disp(mineigenvector_p);

disp("Определяем невязку для максимального значения" );
% Невязка для максимального значения Av - lambdav = 0
Av = A*maxeigenvector;
Av_p = A_p*maxeigenvector_p;

lamv = maxeigenvalue*maxeigenvector;
lamv_p = maxeigenvalue_p*maxeigenvector_p;

maxerroreigen = Av-lamv ;
maxerroreigen_p = Av_p-lamv_p ;

disp(" для исходной системы:  " );
disp(norm(maxerroreigen));
disp(" для возмущенной системы: " );
disp(norm(maxerroreigen_p));

disp("Определяем невязку для минимального значения" );

Avmin = A*mineigenvector;
Av_pmin = A_p*mineigenvector_p;

lamvmin = mineigenvalue*mineigenvector;
lamv_pmin = mineigenvalue_p*mineigenvector_p;

minerroreigen = Avmin-lamvmin;
minerroreigen_p = Av_pmin-lamv_pmin;

disp(" для исходной системы:  " );
disp(norm(minerroreigen));
disp(" для возмущенной системы: " );
disp(norm(minerroreigen_p ));