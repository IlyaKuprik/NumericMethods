import seaborn as sns
import matplotlib.pyplot as plt
from semester_5.utils import read_data


def plot_graph(nums=None, method="normal", save_flag=False):
    if nums is None:
        nums = [1]
    original = read_data("data/original.txt")
    graph = []
    for n in nums:
        graph.append(read_data(f"data/{method}_{n}.txt"))
    fig, ax = plt.subplots(figsize=(16, 10))
    ax.grid(which="major", linewidth=1.2)
    ax.grid(which="minor", linestyle="--", color="gray", linewidth=0.5)

    sns.lineplot(x=original[0], y=original[1])
    for i in range(len(nums)):
        sns.lineplot(x=graph[i][0], y=graph[i][1])

    legend_list = ["f(x)"]
    for n in nums:
        legend_list.append(f"{method}_{n}(x)")

    noise_dots = read_data("data/noise.txt")
    sns.scatterplot(x=noise_dots[0], y=noise_dots[1], color=".5", marker='.')
    plt.legend(legend_list)
    if save_flag:
        plt.savefig(
            f'C:/Users/Ilya/Desktop/Учеба/семестр 5/численные методы/ТЕМА №3 _ Аппроксимация функций/графики/{method}_{nums[0]}.png')
    else:
        plt.show()


plot_graph([1, 2, 3, 4, 5], "normal")
