from semester_6.task_5.schemes import TwoStageCalculationScheme, ThreeStageCalculationScheme
from semester_6.task_5.diffirental import *
import seaborn as sns
import matplotlib.pyplot as plt

TRUE_Y_VALUE = np.array([0.2219309915117, 0.2319295303698])


def y_val(x):
    y1 = (20 / (np.sqrt(12 * 20))) * (
            (np.pi / 12) * np.sin(x / (np.sqrt(12 * 20))) + (np.sqrt(12 * 20) / 400) * np.pi * np.cos(
        x / (np.sqrt(12 * 20))))
    y2 = (np.pi / 12) * np.cos(x / (np.sqrt(12 * 20))) - (np.sqrt(12 * 20) / 400) * np.pi * np.sin(
        x / (np.sqrt(12 * 20)))
    return np.array([y1, y2])


def const_step_err_graphic(h2, h3):
    res2 = get_runge_method_res(h2, X_0, Y_10, Y_20, X_K, TwoStageCalculationScheme, 2, log=True)
    res3 = get_runge_method_res(h3, X_0, Y_10, Y_20, X_K, ThreeStageCalculationScheme, 3, log=True)
    err_lst2 = []
    x_lst2 = []
    for i in range(len(res2)):
        err_lst2.append(np.linalg.norm(res2[i][1] - y_val(res2[i][0])))
        x_lst2.append(res2[i][0])

    err_lst3 = []
    x_lst3 = []
    for i in range(len(res3)):
        err_lst3.append(np.linalg.norm(res3[i][1] - y_val(res3[i][0])))
        x_lst3.append(res3[i][0])

    sns.lineplot(x=x_lst2, y=err_lst2)
    sns.lineplot(x=x_lst3, y=err_lst3)
    plt.title('Графики зависимости истинной полной погрешности от значения независимой переменной x')
    plt.xlabel('x')
    plt.ylabel('погрешность')
    plt.legend(['Двухэтапный метод Рунге', 'Трёхэтапный метод Рунге'])
    plt.show()


def auto_step_val_graphic():
    res2 = runge_method_auto_step(TwoStageCalculationScheme, 2, eps=1e-5, log=True)
    res3 = runge_method_auto_step(ThreeStageCalculationScheme, 3, eps=1e-5, log=True)
    h_lst2 = []
    x_lst2 = []
    for i in range(len(res2)):
        h_lst2.append(res2[i][1])
        x_lst2.append(res2[i][0])

    h_lst3 = []
    x_lst3 = []
    for i in range(len(res3)):
        h_lst3.append(res3[i][1])
        x_lst3.append(res3[i][0])

    sns.lineplot(x=x_lst2, y=h_lst2)
    sns.lineplot(x=x_lst3, y=h_lst3)
    plt.title('Графики зависимости величины шага интегрирования от значения независимой переменной x')
    plt.xlabel('x')
    plt.ylabel('Величина шага')
    plt.legend(['Двухэтапный метод Рунге', 'Трёхэтапный метод Рунге'])
    plt.show()


def auto_step_local_err_graphic():
    res2 = runge_method_auto_step(TwoStageCalculationScheme, 2, eps=1e-5, log=True)
    res3 = runge_method_auto_step(ThreeStageCalculationScheme, 3, eps=1e-5, log=True)
    loc_approximate_err_lst2 = []
    loc_real_err_lst2 = []
    x_lst2 = []
    for i in range(len(res2)):
        loc_approximate_err_lst2.append(res2[i][3])
        loc_real_err_lst2.append(res2[i][2] - y_val(res2[i][0]))
        x_lst2.append(res2[i][0])

    loc_approximate_err_lst3 = []
    loc_real_err_lst3 = []
    x_lst3 = []
    for i in range(len(res3)):
        loc_approximate_err_lst3.append(res3[i][3])
        loc_real_err_lst3.append(res3[i][2] - y_val(res3[i][0]))
        x_lst3.append(res3[i][0])

    loc_approximate_err_lst2 = np.array(loc_approximate_err_lst2)
    loc_approximate_err_lst3 = np.array(loc_approximate_err_lst3)
    loc_real_err_lst2 = np.array(loc_real_err_lst2)
    loc_real_err_lst3 = np.array(loc_real_err_lst3)

    sns.lineplot(x=x_lst2, y=np.linalg.norm(loc_real_err_lst2 / loc_approximate_err_lst2))
    sns.lineplot(x=x_lst3, y=np.linalg.norm(loc_real_err_lst3 / loc_approximate_err_lst3))
    plt.title('Графики зависимости отношения локальной погрешности к оценке погрешности от значения независимой '
              'переменной x')
    plt.xlabel('x')
    plt.ylabel('Велечина отношения')
    plt.legend(['Двухэтапный метод Рунге', 'Трёхэтапный метод Рунге'])
    plt.show()


def auto_step_calc_num_graphic(s, CalculationScheme):
    n_1 = runge_method_auto_step(CalculationScheme, s, eps=1e-1, log=True)[1]
    n_2 = runge_method_auto_step(CalculationScheme, s, eps=1e-2, log=True)[1]
    n_3 = runge_method_auto_step(CalculationScheme, s, eps=1e-3, log=True)[1]
    n_4 = runge_method_auto_step(CalculationScheme, s, eps=1e-4, log=True)[1]
    n_5 = runge_method_auto_step(CalculationScheme, s, eps=1e-5, log=True)[1]
    n_6 = runge_method_auto_step(CalculationScheme, s, eps=1e-6, log=True)[1]
    n_7 = runge_method_auto_step(CalculationScheme, s, eps=1e-7, log=True)[1]



def print_result():
    res_1 = runge_method_const_step(TwoStageCalculationScheme, s=2, eps=1e-4)
    print(f"Двухэтапный метод Рунге-Кутты:\n\tС постоянным шагом:\n\t\tШаг: {res_1[0]}"
          f"\n\t\tПогрешность: {np.linalg.norm(y_val(np.pi) - res_1[1])}")

    res_2 = runge_method_auto_step(TwoStageCalculationScheme, s=2, eps=1e-5)
    print(f"\tС атоматическим шагом:\n\t\tШаг: {res_2[0]}"
          f"\n\t\tПогрешность: {np.linalg.norm(y_val(np.pi) - res_2[1])}")

    res_3 = runge_method_const_step(ThreeStageCalculationScheme, s=3, eps=1e-4)
    print(f"\nТрёхэтапный метод Рунге-Кутты:\n\tС постоянным шагом:\n\t\tШаг: {res_3[0]}"
          f"\n\t\tПогрешность: {np.linalg.norm(y_val(np.pi) - res_3[1])}")

    res_4 = runge_method_auto_step(ThreeStageCalculationScheme, s=3, eps=1e-5)
    print(f"\tС атоматическим шагом:\n\t\tШаг: {res_4[0]}"
          f"\n\t\tПогрешность: {np.linalg.norm(y_val(np.pi) - res_4[1])}")


# print_result()  # пункт 1.2, 2.1, 3.2, 3.3
#
# const_step_err_graphic(h2=0.14580103006207398/2, h3=0.3141570538911852/2)  # пункт 3.2
#
# auto_step_val_graphic()  # пункт 3.3.1
#
# auto_step_local_err_graphic()  # пункт 3.3.2
