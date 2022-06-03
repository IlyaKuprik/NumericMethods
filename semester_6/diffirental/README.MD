# Численные методы решения задачи Коши для СОДУ
Цель задания: практическое освоение методов управления численным процессом с целью 
обеспечения заданной точности решения задачи Коши для СОДУ за приемлемое время, используя 
явные методы Рунге-Кутты (ЯМРК) с методами оценки полной и локальной погрешности.

## Условия и навигация по файлам
Условия задач находятся в файле _Задание №5_ЧИСЛЕННЫЕ МЕТОДЫ РЕШЕНИЯ ЗАДАЧИ КОШИ ДЛЯ СОДУ.pdf_

Вариант 8.
Значения констант:
  _ksi_ = 1/12, _A_ = 1/12, _B_ = 1/20

### Навигация по файлам

* Расчётные схемы методов Рунге-Кутты 2, 3 и 4 порядка: [schemes.py](https://github.com/IlyaKuprik/NumericMethods/blob/master/semester_6/diffirental/schemes.py) (4-этапная схема не работает)
* Реализация методов Рунге-Кутты с постоянным и автоматическим шагом на основе задаваемой расчётной схемы: [diffirental.py](https://github.com/IlyaKuprik/NumericMethods/blob/master/semester_6/diffirental/diffirental.py)
* Вывод результатов и построение графиков: [main_diffirental.py](https://github.com/IlyaKuprik/NumericMethods/blob/master/semester_6/diffirental/main_diffirental.py)

## Анализ эффективности

В этом разделе представлены все графики, которые требовались в задании.

<img src="https://github.com/IlyaKuprik/NumericMethods/blob/master/images/task5/график%20зависимости%20истинной%20полной%20погрешности.png" alt="normal_1" width="800"/>

<img src="https://github.com/IlyaKuprik/NumericMethods/blob/master/images/task5/график%20зависимости%20величины%20шага%20интегрирования.png" alt="normal_1" width="800"/>

<img src="https://github.com/IlyaKuprik/NumericMethods/blob/master/images/task5/график%20зависимости%20отношения%20локальной%20погрешности.png" alt="normal_1" width="800"/>

<img src="https://github.com/IlyaKuprik/NumericMethods/blob/master/images/task5/график%20зависимочти%20количества%20вычислений%20от%20погрешности.png" alt="normal_1" width="800"/>
