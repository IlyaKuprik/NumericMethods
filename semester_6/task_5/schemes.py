class CalculationScheme:
    def __init__(self, x_0, y_10, y_20, A, B, c_2, s):
        self.s = s
        self.c_2 = c_2
        self.a_21 = c_2
        self.b_2 = 1 / (2 * c_2)
        self.b_1 = 1 - self.b_2

        self.A = A
        self.B = B

        self.x_0 = x_0
        self.y_10 = y_10
        self.y_20 = y_20


class TwoStageCalculationScheme(CalculationScheme):
    """Класс расчётной схемы двухэтапного метода Рунге-Кутты"""

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


class FourStageCalculationScheme(CalculationScheme):
    """Класс расчётной схемы четырёхэтапного метода Рунге-Кутты"""
    def y_1(self, x):
        h = abs(self.x_0 - x)
        if h < 1e-9:
            return self.y_10
        return self.y_1(self.x_0) + (1/6) * (
                self.k_1(1, h) + 2 * self.k_1(2, h) + 2 * self.k_1(3, h) + self.k_1(4, h))

    def y_2(self, x):
        h = abs(self.x_0 - x)
        if h < 1e-9:
            return self.y_20
        return self.y_2(self.x_0) + (1 / 6) * (
                self.k_2(1, h) + 2 * self.k_2(2, h) + 2 * self.k_2(3, h) + self.k_2(4, h))

    def k_1(self, j, h):
        # j - воспринимать как k_1j, для удобства
        if j == 1:
            return self.A * h * self.y_2(self.x_0)
        if j == 2:
            return self.A * h * (self.y_2(self.x_0) + 0.5*self.k_2(1, h))
        if j == 3:
            return self.A * h * (self.y_2(self.x_0) + 0.5 * self.k_2(2, h))
        return self.A * h * (self.y_2(self.x_0) + self.k_2(3, h))

    def k_2(self, j, h):
        # j - воспринимать как k_2j, для удобства
        if j == 1:
            return (-self.B) * h * self.y_1(self.x_0)
        if j == 2:
            return (-self.B) * h * (self.y_1(self.x_0) + 0.5 * self.k_1(1, h))
        if j == 3:
            (-self.B) * h * (self.y_1(self.x_0) + 0.5 * self.k_1(2, h))
        return (-self.B) * h * (self.y_1(self.x_0) + self.k_1(3, h))
