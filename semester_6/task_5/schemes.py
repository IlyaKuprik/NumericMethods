
class TwoStageCalculationScheme:
    # Класс расчётной схемы двухэтапного метода Рунге-Кутты
    def __init__(self, x_0, y_10, y_20, A, B, c_2):
        self.s = 2
        self.c_2 = c_2
        self.a_21 = c_2
        self.b_2 = 1/(2*c_2)
        self.b_1 = 1 - self.b_2

        self.A = A
        self.B = B

        self.x_0 = x_0
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