import matplotlib.pyplot
import math

import scipy.misc

import approx


def input_resolver():
    answer = input("Ввести данные с консоли или из файла?\nc/f >")
    data = []
    if answer == 'c':
        data = read_input_from_console()
    elif answer == 'f':
        data = read_input_from_file()
    else:
        input_resolver()
    func = function_dispatcher(data)
    while True:
        x = float(input("Введите значение для интерполяции:"))
        print("Интерполяция членом Лагранжа:")
        print(str(func[0](x)) + "±" + str(lagrange_accuracy(data, func[0], x)))
        print("Интерполяция методом Ньютона для равноотстоящих узлов:")
        print(str(func[1](x)) + "±" + str(newton_accuracy(data, x)))
        show_plot(data, func, x)


def show_plot(data, function, interplot):
    x = []
    y = []
    for i in range(len(data)):
        x.append(data[i][0])
        y.append(data[i][1])
    matplotlib.pyplot.figure()
    matplotlib.pyplot.subplot()
    matplotlib.pyplot.scatter(x, y)
    matplotlib.pyplot.scatter(interplot, function[0](interplot))
    func_x = []
    func_y = []
    a = data[0][0]
    b = data[len(data) - 1][0]
    if interplot > b:
        b = interplot
    if interplot < a:
        a = interplot
    step = (b - a) / 100
    i = a
    while i <= b:
        func_x.append(i)
        func_y.append(function[0](i))
        i += step
    matplotlib.pyplot.plot(func_x, func_y)
    matplotlib.pyplot.figure()
    matplotlib.pyplot.scatter(x, y)
    matplotlib.pyplot.scatter(interplot, function[1](interplot))
    func_x = []
    func_y = []
    a = data[0][0]
    b = data[len(data) - 1][0]
    if interplot > b:
        b = interplot
    if interplot < a:
        a = interplot
    step = (b - a) / 100
    i = a
    while i <= b:
        func_x.append(i)
        func_y.append(function[1](i))
        i += step
    matplotlib.pyplot.plot(func_x, func_y)
    matplotlib.pyplot.show()


# Метод для чтения данных с консоли
def read_input_from_console():
    values = []
    dot_count = int(input('Введите количество точек: '))
    for i in range(dot_count):
        x = input('Введите координаты X и Y ' + str(i + 1) + '-ой точки через пробел: ')
        values.append([float(x.split(" ")[0]), float(x.split(" ")[1])])
    return values


# Метод для чтения данных из файла
def read_input_from_file():
    values = []
    path_to_file = input('Введите путь до файла: ')
    file = open(path_to_file, 'r')
    lines = file.readlines()
    for line in lines:
        values.append([float(line.split(" ")[0]), float(line.split(" ")[1])])
    return values


def function_dispatcher(data):
    lagrange_f = lambda x: lagrange_polynomial(data, x)
    newton_f = lambda x: newton_polynomial(data, x)
    return [lagrange_f, newton_f]


def lagrange_polynomial(data, x):
    result = 0
    i = 0
    for i in range(len(data)):
        j = 0
        denominator = 1
        numerator = 1
        for j in range(len(data)):
            if i != j:
                numerator = numerator * (x - data[j][0])
                denominator = denominator * (data[i][0] - data[j][0])
        result = result + data[i][1] * (numerator / denominator)
    return result


def lagrange_accuracy(data, function, x):
    f = approx.function_dispatcher(data)[0]
    accuracy = max_on_derivative(f, len(data) + 1, data[0][0], data[len(data) - 1][0])
    accuracy = accuracy / math.factorial(len(data) + 1)
    i = 0
    for dot in data:
        accuracy = accuracy * abs(x - dot[0])
    return accuracy


def max_on_derivative(function, n, a, b):
    step = (b - a) / 10
    maximum = derivative(function, n, a)
    iterator = a
    while iterator <= b:
        der = derivative(function, n, iterator)
        if der > maximum:
            maximum = der
        iterator = iterator + step
    return maximum


def derivative(function, n, x):
    e = 1E-3
    result = 0
    i = 0
    for i in range(n):
        result = result + ((-1) ** i) * (math.factorial(n) / (math.factorial(n - i)
                                                              * math.factorial(i))) * function(x + (n/2 - i) * e)
    return result


def newton_polynomial(data, x):
    starter_i = 0
    i = 0
    for i in range(len(data)):
        if x >= data[i][0]:
            starter_i = i
    h = data[1][0] - data[0][0]
    result = 0
    i = 0
    if starter_i < len(data) / 2:
        t = (x - data[starter_i][0]) / h
        for i in range(len(data) - starter_i):
            j = 0
            nominator = 1
            for j in range(i):
                nominator = nominator * (t - j)
            result = result + (nominator / math.factorial(i)) * delta_y(data, starter_i, i)
    else:
        t = (x - data[len(data) - 1][0]) / h
        for i in range(len(data)):
            j = 0
            nominator = 1
            for j in range(i):
                nominator = nominator * (t + j)
            result = result + (nominator / math.factorial(i)) * delta_y(data, len(data) - 1 - i, i)
    return result


def delta_y(data, i, n):
    if n == 0:
        return data[i][1]
    else:
        return delta_y(data, i + 1, n - 1) - delta_y(data, i, n - 1)


def newton_accuracy(data, x):
    starter_i = 0
    i = 0
    for i in range(len(data)):
        if x >= data[i][0]:
            starter_i = i
    h = data[1][0] - data[0][0]
    if starter_i < len(data) / 2:
        result = delta_y(data, starter_i, len(data) - starter_i - 1)
        t = (x - data[starter_i][0]) / h
        i = 0
        for i in range(len(data) + 1):
            result = result * abs(t - i)
        result = result / math.factorial(len(data) + 1)
        return result
    else:
        result = delta_y(data, starter_i, len(data) - starter_i - 1)
        t = (x - data[len(data) - 1][0]) / h
        i = 0
        for i in range(len(data) + 1):
            result = result * abs(t + i)
        result = result / math.factorial(len(data) + 1)
        return result


# Точка входа в программу
input_resolver()
