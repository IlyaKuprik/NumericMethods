import seaborn as sns
import matplotlib.pyplot as plt
from semester_5.utils import read_data
from math import acos


def my_func(x):
    return x * x + 1 - acos(x)


def plot_big_graph(method='lagrange', key='eq'):
    original = read_data("C:/Users/Ilya/PycharmProjects/numeric_methods/interpolation/data/original.txt")
    L_n = []
    nums = [3, 10, 20, 30, 40, 50]

    for n in nums:
        L_n.append(
            read_data(f"C:/Users/Ilya/PycharmProjects/numeric_methods/interpolation/data/{key}_{method}_{n}.txt"))

    sns.lineplot(x=original[0], y=original[1])
    for i in range(len(L_n)):
        sns.lineplot(x=L_n[i][0], y=L_n[i][1])

    legend_list = ["f(x)"]
    for n in nums:
        legend_list.append(f"{method}_{n}_{key}(x)")

    plt.legend(legend_list)
    plt.show()


def plot_deviation_graph(method='lagrange', key='eq'):
    #original = read_data("C:/Users/Ilya/PycharmProjects/numeric_methods/interpolation/data/original.txt")

    data = read_data(f"C:/Users/Ilya/PycharmProjects/numeric_methods/interpolation/data/{key}_{method}_{40}.txt")
    print(data)
    for i in range(len(data[0])):
        data[1][i] = abs(my_func(data[0][i]) - data[1][i])
    sns.lineplot(x=data[0], y=data[1])
    plt.title("график распределения абсолютной погрешности по оптимальным узлам")
    plt.legend([f"Погрешность {method}"])
    plt.show()


#plot_deviation_graph(method='lagrange', key='opt')

#plot_deviation_graph(method='spline_32', key='opt')

plot_big_graph(method="spline_21", key = "eq")
plot_big_graph(method="spline_21", key = "opt")