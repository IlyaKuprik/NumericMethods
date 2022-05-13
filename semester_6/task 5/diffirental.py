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
        self.x = 0

    def y_1(self, x):
        if abs(x - self.x) < 1e-9:
            return self.B * np.pi
        # return self.y_1(x) + self.b_1*self.k11(h) +

    def y_2(self, x):
        pass

    def k_1(self, j, h):
        # j - воспринимать как k_1j, для удобства
        pass

    def k_2(self, j, h):
        # j - воспринимать как k_1j, для удобства
        pass
