import numpy as np

X_0 = 0
X_K = np.pi
A = 1 / 12
B = 1 / 20
C_2 = 1 / 12
Y_10 = B * np.pi
Y_20 = A * np.pi


def get_first_step(A, B, y_1, y_2, x_0, x_k, s, eps=1e-4):
    f = [A * y_2, -B * y_1]
    delta = (1 / (max(x_0, x_k))) ** (s + 1) + np.linalg.norm(f) ** (s + 1)
    return (eps / delta) ** (1 / (1 + s))


def get_runge_err(res1, res2, s):
    """Оценка погрешности по методу Рунге"""
    res1 = np.array(res1)
    res2 = np.array(res2)
    return (res2 - res1) / (2 ** s - 1)


def get_runge_method_res(h, x_0, y_01, y02, x_k, CalculationScheme, s, log=False, to_calc=False):
    """
    Метод рунге
    Начало в x_0, конец в x_k
    """
    scheme = CalculationScheme(x_0, y_01, y02, A, B, C_2, s)
    scheme.x_0 = x_0
    x_i = x_0
    res_list = []
    while x_i < x_k:
        if x_i + h > x_k:
            h = x_k - x_i
        x_i += h
        y_1i = scheme.y_1(x_i)
        y_2i = scheme.y_2(x_i)
        # print(f"В точке x_i = {x_i} получили значения y_1(x_i) = {y_1i}, y_2(x_i) = {y_2i}")
        scheme.x_0 = x_i
        res_list.append((x_i, np.array([y_1i, y_2i]), np.array([scheme.y_10, scheme.y_20])))
        scheme.y_10 = y_1i
        scheme.y_20 = y_2i
    if log:
        return res_list
    if to_calc:
        return np.array([scheme.y_1(x_k), scheme.y_2(x_k)]), scheme.num_of_calc
    return np.array([scheme.y_1(x_k), scheme.y_2(x_k)])


def runge_method_const_step(CalculationScheme, s, eps=1e-4):
    step = get_first_step(A, B, Y_10, Y_20, X_0, X_K, s, eps=eps)
    res1 = get_runge_method_res(step, X_0, Y_10, Y_20, X_K, CalculationScheme, s)
    res2 = get_runge_method_res(step / 2, X_0, Y_10, Y_20, X_K, CalculationScheme, s)
    while np.linalg.norm(get_runge_err(res1, res2, s)) > eps:
        step /= 2
        res1 = get_runge_method_res(step, X_0, Y_10, Y_20, X_K, CalculationScheme, s)
        res2 = get_runge_method_res(step / 2, X_0, Y_10, Y_20, X_K, CalculationScheme, s)
        print(f"Результаты с шагом step: {res1}\nРезультаты с шагом step/2: {res2}")
        print(f"Шаг: {step}\nОшибка: {get_runge_err(res1, res2, s)}")
    return step / 2, res2  # без уточнения + get_runge_err(res1, res2, s)


def runge_method_auto_step(CalculationScheme, s, eps=1e-5, log=False, to_calc=False):
    h = get_first_step(A, B, Y_10, Y_20, X_0, X_K, s, eps=eps)

    x_i = X_0
    y_10 = Y_10
    y_20 = Y_20

    res_lst = []
    num_of_calc = 0

    while x_i < X_K:
        if x_i + h > X_K:
            h = X_K - x_i
        # результаты (y1(x_i+h), y2(x_i+h)) в  точке x_i + h с шагом h
        res1 = get_runge_method_res(h, x_i, y_10, y_20, x_i + h, CalculationScheme, s, to_calc=to_calc)
        # результаты (y1(x_i+h), y2(x_i+h)) в точке x_i + h с шагом h/2 (делается 2 половинных шага)
        res2 = get_runge_method_res(h / 2, x_i, y_10, y_20, x_i + h, CalculationScheme, s, to_calc=to_calc)

        if to_calc:
            num_of_calc += res1[1] + res2[1]
            res1 = res1[0]
            res2 = res2[0]
        added = get_runge_err(res1, res2, s)
        r = np.linalg.norm(added)
        y_m10, y_m20 = y_10, y_20
        if r > eps * (2 ** s):
            h /= 2  # остаемся в x_i, но шаг уменьшаем
            continue
        elif eps < r <= eps * (2 ** s):
            x_i += h
            h /= 2  # уменьшаем шаг и в качестве приблежения выбираем результаты с половинным шагом
            y_10, y_20 = map(float, res2)  # без уточнения + added
        elif eps / (2 ** (s + 1)) <= r:
            x_i += h  # шаг такой же, приближение выбирается с полным шагом
            y_10, y_20 = map(float, res1)  # без уточнения + added
        else:
            x_i += h
            h *= 2  # увеличиваем шаг в 2 раза, приближение с полным шагом
            y_10, y_20 = map(float, res1)  # без уточнения + added
        res_lst.append((x_i, h, np.array([y_10, y_20]), added, np.array([y_m10, y_m20])))
    if log:
        return res_lst
    if to_calc:
        return num_of_calc
    return h, np.array([y_10, y_20])
