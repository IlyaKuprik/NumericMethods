import numpy as np


class CalculationScheme:
    # TODO: Класс расчётной схемы двухэтапного метода Рунге-Кутты
    def __init__(self):
        self.c_2 = 1 / 12
        self.a_21 = 1 / 12
        self.b_2 = 6
        self.b_1 = -5
        self.A = 1 / 12
        self.B = 1 / 20
        self.x_0 = 0

    def y_1(self, x):
        h = abs(self.x_0 - x)
        if h < 1e-9:
            return self.B * np.pi
        return self.y_1(self.x_0) + self.b_1 * self.k_1(1, h) + self.b_2 * self.k_1(2, h)

    def y_2(self, x):
        h = abs(self.x_0 - x)
        if h < 1e-9:
            return self.A * np.pi
        return self.y_2(self.x_0) + self.b_1 * self.k_2(1, h) + self.b_2 * self.k_2(2, h)

    def k_1(self, j, h):
        # j - воспринимать как k_1j, для удобства
        if j == 1:
            return self.A * h * self.y_2(self.x_0)

        return self.A * h * (self.y_2(self.x_0) + self.a_21 * self.k_2(1, h))

    def k_2(self, j, h):
        # j - воспринимать как k_2j, для удобства
        if j == 1:
            return -self.B * h * self.y_1(self.x_0)

        return -self.B * h * (self.y_1(self.x_0) + self.a_21 * self.k_1(1, h))


scheme = CalculationScheme()
