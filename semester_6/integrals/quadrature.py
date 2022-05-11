import numpy as np
from math import acos, sqrt, cos, pi, log

TRUE_INTEGRAL_WITH_WEIGHT_VALUE = 20.7302711095


def simple_quadratures(function, a, b, n_dots=100, method="") -> float:
    """Вычисление интеграла с помощью простейших квадратурных формул"""
    integral_value = 0
    for i in range(n_dots):
        x_curr = a + ((b - a) / n_dots) * i
        x_next = a + ((b - a) / n_dots) * (i + 1)
        if method == "mean":
            integral_value += function((x_next + x_curr) / 2) * (x_next - x_curr)
        elif method == "left":
            integral_value += function(x_curr) * (x_next - x_curr)
        elif method == "right":
            integral_value += function(x_next) * (x_next - x_curr)
        elif method == "trapezoid":
            integral_value += (function(x_next) + function(x_curr)) * (x_next - x_curr) / 2
        elif method == "Simpson":
            integral_value += (function(x_next) + function(x_curr) + 4 * function((x_next + x_curr) / 2)) * (
                    x_next - x_curr) / 6
        else:
            assert Exception, "Неправильное название метода."
    return integral_value


def split_interval(a, b, n_dots) -> list[float]:
    return [a + ((b - a) / (n_dots - 1)) * i for i in range(n_dots)]


def moment_j(j, z1, z2, a, b) -> float:
    """Вычислление j-го момента для функции"""
    alpha = 2 / 3
    return ((z2 - a) ** (j - alpha + 1) - (z1 - a) ** (j - alpha + 1)) / (j - alpha + 1)


def newton_cotes(function, a, b, n_dots) -> float:
    """Составная квадратурная формула на базе 3-точечной формулы Ньютона-Котеса"""
    splitted_interval = split_interval(a, b, n_dots)
    integral_value = 0
    for i in range(0, n_dots - 1):
        # выбираем 3 равноотстоящих узла
        nodes = [splitted_interval[i], (splitted_interval[i] + splitted_interval[i + 1]) / 2, splitted_interval[i + 1]]
        mu_vec = [moment_j(j, nodes[0], nodes[2], a, b) for j in range(3)]
        # вычисление матрицы x
        X_matrix = [[(nodes[k] - a) ** j for k in range(3)] for j in range(3)]
        # вычисление вектора A методом гаусса решения СЛАУ
        A_vec = np.linalg.solve(X_matrix, mu_vec)
        integral_value += A_vec @ [function(nodes[j]) for j in range(3)]

    return integral_value


def kardano(a, b, c, d) -> tuple[float, float, float]:
    """Формула Кардано для нахождения корней уравнения  3й степени"""
    p = (3 * a * c - b * b) / (9 * a * a)
    q = ((2 * b * b * b) / (27 * a * a * a) - (b * c) / (3 * a * a) + d / a) / 2

    r = sqrt(abs(p))
    r *= (-1 if q < 0 else 1)

    fi = acos(q / (r * r * r))

    return (-2 * r * cos(fi / 3) - b / (3 * a), 2 * r * cos(pi / 3 - fi / 3) - b / (3 * a),
            2 * r * cos(pi / 3 + fi / 3) - b / (3 * a))


def count_gauss_nodes(x_1, x_2, a, b) -> tuple[float, float, float]:
    """Вычисление узлов для трёхточечного метода Гаусса"""
    mu_vec_right = [-moment_j(j, x_1, x_2, a, b) for j in range(3, 6)]
    mu_matrix_left = [[moment_j(j + s, x_1, x_2, a, b) for j in range(3)] for s in range(3)]
    # решаем слау для того чтобы найти коэфиценты полинома
    coeffs = [1] + list(np.linalg.solve(mu_matrix_left, mu_vec_right))
    coeff_a, coeff_d, coeff_c, coeff_b = map(float, coeffs)

    return kardano(coeff_a, coeff_b, coeff_c, coeff_d)


def gauss(function, a, b, n_dots) -> float:
    """Составная квадратурная формула на базе 3-точечной формулы Гаусса"""
    splitted_interval = split_interval(a, b, n_dots)
    integral_value = 0
    for i in range(0, n_dots - 1):
        # получение узлов
        nodes = count_gauss_nodes(splitted_interval[i], splitted_interval[i + 1], a, b)

        mu_vec = [moment_j(j, splitted_interval[i], splitted_interval[i + 1], a, b) for j in range(3)]
        # вычисление матрицы x
        X_matrix = [[(nodes[k]) ** j for k in range(3)] for j in range(3)]
        # вычисление вектора A методом гаусса решения СЛАУ
        A_vec = np.linalg.solve(X_matrix, mu_vec)
        integral_value += A_vec @ [function(nodes[j] + a) for j in range(3)]

    return integral_value


def richardson_error(steps, results_list, m) -> tuple[float, float]:
    """Оценка погрешности методом Ричардсона"""
    # Задаем сетку шагов
    # r = min(len(steps), 5)
    steps_matrix = [[i ** (m + j) for j in range(len(steps))] for i in steps]
    # Вычисляем коэфиценты С_m в методе Ричардсона путем решения СЛАУ
    C = np.linalg.solve(steps_matrix, results_list)
    return abs(results_list[-1] - C @ steps_matrix[-1]), results_list[-1]


def precise_three_dots_quadratures(function, a, b, eps=1e-6, L=2, method="Newton") -> tuple[float, float]:
    """Вычисление определенного интеграла с заданной точностью"""
    steps = [(b - a) / L, (b - a) / (L * L), (b - a) / (L * L * L)]
    num_splits_lst = [L, L * L, L ** 3]
    res_list = None
    if method == "Newton":
        res_list = [newton_cotes(function, a, b, k) for k in num_splits_lst]
    elif method == "Gauss":
        res_list = [gauss(function, a, b, k) for k in num_splits_lst]

    m = -log(abs((res_list[-1] - res_list[-2]) / (res_list[-2] - res_list[-3]))) / log(L ** len(steps))
    err, val = map(float, richardson_error(steps, res_list, m))
    print("Порядок погрешности: ", m)
    while err > eps:
        print("Порядок погрешности: ", m)
        steps.append(steps[-1] / L)
        if method == "Newton":
            res_list.append(newton_cotes(function, a, b, (L ** len(steps))))
        elif method == "Gauss":
            res_list.append(gauss(function, a, b, (L ** len(steps))))

        m = -log(abs((res_list[-1] - res_list[-2]) / (res_list[-2] - res_list[-3]))) / log(L ** len(steps))
        err, val = map(float, richardson_error(steps, res_list, m))

    return steps[-1], len(steps), abs(val - TRUE_INTEGRAL_WITH_WEIGHT_VALUE)


def get_optimal_step(function, a, b, eps=1e-6, L=2):
    steps_grid = [(b - a) / 2, (b - a) / 3, (b - a) / 4]
    sh1 = gauss(function, a, b, 2)
    sh2 = gauss(function, a, b, 3)
    sh3 = gauss(function, a, b, 4)
    m = -log(abs((sh3 - sh2) / (sh2 - sh1)))
    h_opt = steps_grid[0] * pow((eps * (1 - L ** (-m)) / abs(sh2 - sh1)), 1 / m)

    # считаем интеграл, начиная с оптимального шага
    steps = steps_grid + [h_opt]
    res_list = [-sh1, -sh2, -sh3, -gauss(function, a, b, round((b - a) / h_opt))]
    print(m)
    err, val = map(float, richardson_error(steps, res_list, m))
    iters = 1
    while err > eps:
        print(res_list)
        m = -log(abs((res_list[-1] - res_list[-2]) / (res_list[-2] - res_list[-3]))) / log(L ** len(steps))
        err, val = map(float, richardson_error(steps, res_list, m))

        steps.append(steps[-1] / L)
        res_list.append(-gauss(function, a, b, (L ** len(steps))))
        iters += 1

    return steps[-1], iters, abs(-val - TRUE_INTEGRAL_WITH_WEIGHT_VALUE)
