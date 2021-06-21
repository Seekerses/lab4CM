import copy
import matplotlib.pyplot
import math


def input_resolver():
    answer = input("Ввести данные с консоли или из файла?\nc/f >")
    data = []
    if answer == 'c':
        data = read_input_from_console()
    elif answer == 'f':
        data = read_input_from_file()
    else:
        input_resolver()
    approximation = function_dispatcher(data)
    print("Функция с наименьшим отколнением - " + approximation[3] + " с кэффициентами " + str(approximation[2])
          + "\nСКО = " + str(math.sqrt(deviate(approximation[0], data) / len(data))))
    if approximation[3] == "linear: ax + b":
        print("Коэффициент корреляции = " + str(correlation(data)))
    show_plot(data, approximation[0])


def show_plot(data, function):
    x = []
    y = []
    for i in range(len(data)):
        x.append(data[i][0])
        y.append(data[i][1])
    matplotlib.pyplot.scatter(x, y)
    func_x = []
    func_y = []
    a = data[0][0]
    b = data[len(data) - 1][0]
    step = (b - a) / 100
    a = a - step*10
    b = b + step*10
    i = a
    while i < b:
        func_x.append(i)
        func_y.append(function(i))
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
    approximations = []
    try:
        linear_approximation_result = linear_approximation(data)
        approximations.append([linear_approximation_result[0], linear_approximation_result[1],
                               linear_approximation_result[2], "linear: ax + b"])
    except ValueError:
        pass
    try:
        polynomial_approximation_result = polynomial_approximation(data)
        approximations.append([polynomial_approximation_result[0], polynomial_approximation_result[1],
                               polynomial_approximation_result[2], "polynomial: ax^2 + bx + c"])
    except ValueError:
        pass
    try:
        exponential_approximation_result = exponential_approximation(data)
        approximations.append([exponential_approximation_result[0], exponential_approximation_result[1],
                               exponential_approximation_result[2], "exponential: a e^(bx)"])
    except ValueError:
        pass
    try:
        logarithmic_approximation_result = logarithmic_approximation(data)
        approximations.append([logarithmic_approximation_result[0], logarithmic_approximation_result[1],
                               logarithmic_approximation_result[2], "logarithmic: a ln(x) + b"])
    except ValueError:
        pass
    try:
        power_approximation_result = power_approximation(data)
        approximations.append([power_approximation_result[0], power_approximation_result[1],
                               power_approximation_result[2], "power: ax^b"])
    except ValueError:
        pass
    approximations = sorted(approximations, key=lambda x: x[1])
    return approximations[0]


def linear_approximation(data):
    coefficients = [[sum_x(data, 2), sum_x(data, 1)], [sum_x(data, 1), len(data)]]
    vector = [sum_xy(data, 1, 1), sum_y(data, 1)]
    roots = cramer_method(coefficients, vector)
    function = lambda x: roots[0] * x + roots[1]
    dev = deviate(function, data)
    return [function, dev, roots]


def polynomial_approximation(data):
    matrix_coefficients = [[len(data), sum_x(data, 1), sum_x(data, 2)],
                           [sum_x(data, 1), sum_x(data, 2), sum_x(data, 3)],
                           [sum_x(data, 2), sum_x(data, 3), sum_x(data, 4)]]
    matrix_vector = [sum_y(data, 1), sum_xy(data, 1, 1), sum_xy(data, 2, 1)]
    coefficients = cramer_method(matrix_coefficients, matrix_vector)
    coefficients.reverse()
    function = lambda x: coefficients[0] * (x ** 2) + coefficients[1] * x + coefficients[2]
    dev = deviate(function, data)
    return [function, dev, coefficients]


def exponential_approximation(data):
    coefficients = [[sum(map(lambda x: x[0] ** 2, data)), sum(map(lambda x: x[0], data))],
                    [sum(map(lambda x: x[0], data)), len(data)]]
    vector = [sum(map(lambda x: math.log(x[1]) * x[0], data)), sum(map(lambda x: math.log(x[1]), data))]
    res = cramer_method(coefficients, vector)
    res.reverse()
    function = lambda x: (math.e ** res[0]) * (math.e ** (res[1] * x))
    dev = deviate(function, data)
    return [function, dev, res]


def logarithmic_approximation(data):
    coefficients = [[sum(map(lambda x: math.log(x[0]) ** 2, data)), sum(map(lambda x: math.log(x[0]), data))],
                    [sum(map(lambda x: math.log(x[0]), data)), len(data)]]
    vector = [sum(map(lambda x: math.log(x[0]) * x[1], data)), sum(map(lambda x: x[1], data))]
    res = cramer_method(coefficients, vector)
    function = lambda x: res[0] * math.log(x) + res[1]
    dev = deviate(function, data)
    return [function, dev, res]


def power_approximation(data):
    coefficients = [[sum(map(lambda x: math.log(x[0]) ** 2, data)), sum(map(lambda x: math.log(x[0]), data))],
                    [sum(map(lambda x: math.log(x[0]), data)), len(data)]]
    vector = [sum(map(lambda x: math.log(x[1]) * math.log(x[0]), data)), sum(map(lambda x: math.log(x[1]), data))]
    res = cramer_method(coefficients, vector)
    res.reverse()
    function = lambda x: (math.e ** res[0]) * (x ** res[1])
    dev = deviate(function, data)
    return [function, dev, res]


def deviate(function, data):
    s = 0
    for dot in data:
        s += (function(dot[0]) - dot[1]) ** 2
    return s


def correlation(data):
    xy_sum = 0
    xx_sum = 0
    yy_sum = 0
    x_average = sum_x(data, 1) / len(data)
    y_average = sum_y(data, 1) / len(data)
    for dot in data:
        xy_sum += (dot[0] - x_average) * (dot[1] - y_average)
        xx_sum += (dot[0] - x_average) ** 2
        yy_sum += (dot[1] - y_average) ** 2
    return xy_sum / math.sqrt(xx_sum * yy_sum)


def sum_x(data, power):
    return sum(map(lambda x: x[0] ** power, data))


def sum_y(data, power):
    return sum(map(lambda x: x[1] ** power, data))


def sum_xy(data, power_x, power_y):
    return sum(map(lambda x: (x[0] ** power_x) * (x[1] ** power_y), data))


def cramer_method(coefficients, vector):
    roots = []
    det = determinant(coefficients)
    for i in range(len(vector)):
        det_root = determinant(swap_matrix_columns(coefficients, vector, i))
        roots.append(det_root / det)
    return roots


def swap_matrix_columns(coefficients, vector, index):
    matrix = copy.deepcopy(coefficients)
    for i in range(len(coefficients)):
        matrix[i][index] = copy.deepcopy(vector[i])
    return matrix


def matrix_minor(matrix, index_1, index_2):
    new_matrix = copy.deepcopy(matrix)
    for i in range(len(matrix)):
        new_matrix[i].pop(index_2)
    new_matrix.pop(index_1)
    return new_matrix


def determinant(matrix):
    det = 0
    if len(matrix) > 2:
        pos = 1
        for i in range(len(matrix)):
            det += pos * matrix[0][i] * determinant(matrix_minor(matrix, 0, i))
            pos *= -1
        return det
    det = matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]
    return det


# Точка входа в программу
input_resolver()
