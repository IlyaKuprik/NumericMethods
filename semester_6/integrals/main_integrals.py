from math import cos, exp, sin
import semester_6.integrals.quadrature as quadratures
import seaborn as sns
import matplotlib.pyplot as plt

# a = 2.5, b = 3.3, α = 2/3, β = 0;
a, b = 2.5, 3.3

TRUE_INTEGRAL_WITH_WEIGHT_VALUE = 20.7302711095
TRUE_INTEGRAL_VALUE = 7.71356710860381


def my_function(x):
    return 3 * cos(1.5 * x) * exp(x / 4) + 4 * sin(3.5 * x) * exp(x * (-3)) + 4 * x


def weight(x):
    return pow(x - a, 2 / 3)


def plot_simple_quadratures_graphics(max_n=100):
    methods_list = ["mean", "left", "right", "trapezoid", "Simpson"]
    num_of_dots = [i for i in range(3, max_n + 1)]
    methods_results = []
    for method in methods_list:
        methods_results.append([quadratures.simple_quadratures(my_function,
                                                               a,
                                                               b,
                                                               n_dots=n_dots,
                                                               method=method) - TRUE_INTEGRAL_VALUE
                                for n_dots in range(3, max_n + 1)])

    for result in methods_results:
        sns.lineplot(x=num_of_dots, y=result)
    plt.title("График абсолютной погрешности простейших квадратурных методов")
    plt.legend(["Погрешность метода средних прямоугольников",
                "Погрешность метода левых прямоугольников",
                "Погрешность метода правых прямоугольников",
                "Погрешность метода трапеции",
                "Погрешность метода Симпсона"])
    plt.show()


def plot_tree_dots_quadratures_methods(max_n=100):
    num_of_dots = [i for i in range(3, max_n + 1)]
    sns.lineplot(x=num_of_dots, y=[quadratures.newton_cotes(my_function,
                                                            a,
                                                            b,
                                                            n_dots=n_dots) - TRUE_INTEGRAL_WITH_WEIGHT_VALUE
                                   for n_dots in range(3, max_n + 1)])
    sns.lineplot(x=num_of_dots, y=[quadratures.gauss(my_function,
                                                     a,
                                                     b,
                                                     n_dots=n_dots) - TRUE_INTEGRAL_WITH_WEIGHT_VALUE
                                   for n_dots in range(3, max_n + 1)])
    plt.title("График абсолютной погрешности метода Гаусса и Ньютона-Котеса")
    plt.legend(["Погрешность метода Ньютона-Котеса",
                "Погрешность метода Гаусса"])
    plt.show()


# plot_simple_quadratures_graphics(100)

# plot_tree_dots_quadratures_methods(100)

print("Шаг, с которым сошёлся метод Ньютона-Котса, количество итераций и абсолютная ошибка:",
      quadratures.precise_three_dots_quadratures(my_function, a, b, 1e-6, method="Newton"))
print("Шаг, с которым сошёлся метод Гаусса, количество итераций  и абсолютная ошибка:",
      quadratures.precise_three_dots_quadratures(my_function, a, b, 1e-6, method="Gauss"))

print("Шаг, с которым сошёлся метод Гаусса, количество итераций и абсолютная ошибка, начиная с оптимального шага:",
      quadratures.get_optimal_step(my_function, a, b, eps=1e-6))
