import numpy as np


class CalculationScheme:
    # Класс расчётной схемы двухэтапного метода Рунге-Кутты
    def __init__(self, x_0, y_10, y_20, A=1 / 12, B=1 / 20, x_k=np.pi):
        self.s = 2
        self.c_2 = 1 / 12
        self.a_21 = 1 / 12
        self.b_2 = 6
        self.b_1 = -5

        self.A = A
        self.B = B

        self.x_0 = x_0
        self.x_k = x_k
        self.y_10 = y_10
        self.y_20 = y_20

    def y_1(self, x):
        h = abs(self.x_0 - x)
        if h < 1e-9:
            return self.y_10
        return self.y_1(self.x_0) + self.b_1 * self.k_1(1, h) + self.b_2 * self.k_1(2, h)

    def y_2(self, x):
        h = abs(self.x_0 - x)
        if h < 1e-9:
            return self.y_20
        return self.y_2(self.x_0) + self.b_1 * self.k_2(1, h) + self.b_2 * self.k_2(2, h)

    def k_1(self, j, h):
        # j - воспринимать как k_1j, для удобства
        if j == 1:
            return self.A * h * self.y_2(self.x_0)

        return self.A * h * (self.y_2(self.x_0) + self.a_21 * self.k_2(1, h))

    def k_2(self, j, h):
        # j - воспринимать как k_2j, для удобства
        if j == 1:
            return (-self.B) * h * self.y_1(self.x_0)

        return (-self.B) * h * (self.y_1(self.x_0) + self.a_21 * self.k_1(1, h))


def get_first_step(A, B, y_1, y_2, x_0, x_k, s, eps=1e-4):
    f = [A * y_2, -B * y_1]
    delta = (1 / (max(x_0, x_k))) ** (s + 1) + np.linalg.norm(f) ** (s + 1)
    return (eps / delta) ** (1 / (1 + s))


def get_runge_err(res1, res2, s):
    """Оценка погрешности по методу Рунге"""
    res1 = np.array(res1)
    res2 = np.array(res2)
    return np.linalg.norm((res2 - res1) / (2 ** s - 1))


def get_runge_method_res(A, B, h, x_0, x_k):
    """
    Двухэтапный метод рунге
    Начало в x_0, конец в x_k
    """
    scheme = CalculationScheme(x_0, B * np.pi, A * np.pi)
    scheme.x_0 = x_0
    x_i = x_0
    while x_i < x_k:
        if x_i + h > x_k:
            h = x_k - x_i
        x_i += h
        y_1i = scheme.y_1(x_i)
        y_2i = scheme.y_2(x_i)
        # print(f"В точке x_i = {x_i} получили значения y_1(x_i) = {y_1i}, y_2(x_i) = {y_2i}")
        scheme.x_0 = x_i
        scheme.y_10 = y_1i
        scheme.y_20 = y_2i
    return scheme.y_1(x_k), scheme.y_2(x_k)


def runge_method_const_step(eps=1e-4):
    A = 1 / 12
    B = 1 / 20
    s = 2
    x_0 = 0
    x_k = np.pi
    y_10 = B * np.pi
    y_20 = A * np.pi

    step = get_first_step(A, B, y_10, y_20, x_0, x_k, s)
    res1 = get_runge_method_res(A, B, step, x_0, x_k)
    res2 = get_runge_method_res(A, B, step / 2, x_0, x_k)
    while get_runge_err(res1, res2, s) > eps:
        step /= 2
        res1 = get_runge_method_res(A, B, step, x_0, x_k)
        res2 = get_runge_method_res(A, B, step / 2, x_0, x_k)
        print(f"Результаты с шагом step: {res1}\nРезультаты с шагом step/2: {res2}")
        print(f"Шаг: {step}\nОшибка: {get_runge_err(res1, res2, s)}")
    return step / 2, res2


def runge_method_auto_step(eps=1e-5):
    A = 1 / 12
    B = 1 / 20
    s = 2
    x_0 = 0
    x_k = np.pi
    y_10 = B * np.pi
    y_20 = A * np.pi
    h = get_first_step(A, B, y_10, y_20, x_0, x_k, s)

    x_i = x_0

    while x_i < x_k:
        if x_i + h > x_k:
            h = x_k - x_i
        res1 = get_runge_method_res(A, B, h, x_i, x_i + h)
        res2 = get_runge_method_res(A, B, h / 2, x_i, x_i + h)
        x_i += h

        r = get_runge_err(res1, res2, s)
        # TODO: написать условия дла автоматического выбора шага
        if r > eps * (2 ** s):
            pass
        elif eps < r <= eps * (2 ** s):
            pass
        elif eps / (2 ** (s + 1)) <= eps:
            pass
        else:
            pass


print(runge_method_const_step(eps=1e-4))
