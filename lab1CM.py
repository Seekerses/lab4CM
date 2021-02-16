import copy

# Variable section
# Порядок матрицы
#matrixRank
# Свободные коэффициенты
#absoluteTerms = []
# Коэффициенты при неизвестных
#coefficients = [[]]
# Необходимая точность расчетов
#accuracy
# Максимальное число итераций
#maxIterationsCount

# Input section

# Метод для разрешения способа ввода данных
def inputResolver():
	answer = input("Ввести данные с консоли или из файла?\nc/f >")
	if (answer == 'c'):
		data = readInputFromConsole()
	elif (answer == 'f'):
		data = readInputFromFile()
	else:
		inputResolver()
	resultVector = gaussSeidelMethod(data[1], data[2], data[3], data[4])
	if resultVector:
		print('Вектор результата: ' + str(resultVector[0]))
		print('Количество итераций: ' + str(resultVector[1]))
		print('Вектор погрешностей: ' + str(resultVector[2]))

# Метод для чтения данных с консоли
def readInputFromConsole():
	coefficients = []
	absoluteTerms = []
	matrixRank = int(input('Введите порядок матрицы: '))
	i=0
	for i in range(matrixRank):
		b = float(input('Введите свободный член ' + str(i+1) +'-го уравнения: '))
		absoluteTerms.append(b)
	i=0
	j=0
	for i in range(matrixRank):
		coefficients.append([])
		for j in range(matrixRank):
			a = float(input('Введите коэффициент с индексом ' + str(i + 1) + '' + str(j + 1) + ' : '))
			coefficients[i].append(a)
	accuracy = float(input('Введите желаемую точность: '))
	maxIterationsCount = int(input('Введите максимальное число итераций: '))
	return [matrixRank, absoluteTerms, coefficients, accuracy, maxIterationsCount]

# Метод для чтения данных из файла
def readInputFromFile():
	coefficients = []
	absoluteTerms = []
	pathToInputFile = input('Введите путь к файлу: ')
	file = open(pathToInputFile)
	matrixRank = int(file.readline())
	absoluteTerms = list(map(float, file.readline().split(',')))
	i = 0
	for i in range(matrixRank):
		coefficients.append(list(map(float, file.readline().split(','))))
	accuracy = float(file.readline())
	maxIterationsCount = int(file.readline())
	return [matrixRank, absoluteTerms, coefficients, accuracy, maxIterationsCount]

# Calculater section

# Класс репрезентирующий СЛАУ
class MatrixSLAE:

	def __init__(self, absoluteTerms, coefficients):
		self.absoluteTerms = absoluteTerms
		self.coefficients = coefficients

	def swapRows(self, firsRowIndex, secondRowIndex):
		self.coefficients[firsRowIndex], self.coefficients[secondRowIndex] = self.coefficients[secondRowIndex], self.coefficients[firsRowIndex]
		self.absoluteTerms[firsRowIndex], self.absoluteTerms[secondRowIndex] = self.absoluteTerms[secondRowIndex], self.absoluteTerms[firsRowIndex]

	def toString(self):
		result = ''
		i = 0
		for i in range(len(self.coefficients)):
			result += str(self.coefficients[i]) + " " + str(self.absoluteTerms[i]) + '\n'
		return result

	def rank(self):
		return len(self.absoluteTerms)

# Метод курирующий расчеты, основной метод
def gaussSeidelMethod(absoluteTerms, coefficients, accuracy, maxIterationsCount):
	matrix = MatrixSLAE(absoluteTerms, coefficients)
	matrix = toDiagonalPredominance(matrix)
	if (matrix == None):
		print('Не соблюдается диагональное преобладание...')
		return None
	matrix = toCalculativeForm(matrix)
	resultVector = calculateMatrix(matrix, accuracy, maxIterationsCount)
	if(resultVector == None):
		print('Не достигнут результат за заданное количество итераций...')
	return resultVector

# Метод, приводящий СЛАУ к форме с диагональным преобладанием, возвращает None в случае неудачи
def toDiagonalPredominance(matrix):
	matrix = copy.copy(matrix)
	counterOfIndexes = []
	currentPermutation = []
	maxValueCountInRow = []
	isOneStrict = False

	i = 0
	for i in range(matrix.rank()):
		counterOfIndexes.append([])
		currentPermutation.append(i)
		maxValueCountInRow.append(0)

	i = 0
	for i in range(matrix.rank()):
		maxValueInRow = max(map(abs,matrix.coefficients[i]))
		rowWithoutMaxValue = list(map(abs,matrix.coefficients[i]))
		rowWithoutMaxValue.remove(maxValueInRow)

		if(maxValueInRow == 0):
			return None
		if(maxValueInRow < sum(rowWithoutMaxValue)):
			return None
		if(maxValueInRow > sum(rowWithoutMaxValue)):
			isOneStrict = True

		indexOfIteratedElement = -1
		maxValueCountInRow[i] = matrix.coefficients[i].count(maxValueInRow)
		j = 0
		for j in range(maxValueCountInRow[i]):
			indexOfIteratedElement = matrix.coefficients[i].index(maxValueInRow, indexOfIteratedElement + 1, len(matrix.coefficients[i])+1)
			counterOfIndexes[indexOfIteratedElement].append(i)

	i = 0
	while (max(maxValueCountInRow) > 0):
		for i in range(len(maxValueCountInRow)):
			tempList = list(maxValueCountInRow)
			j = 0
			for j in range(tempList.count(0)):
				tempList.remove(0)
			if (maxValueCountInRow[i] == min(tempList)):
				j = 0
				k = 0
				for j in range(len(counterOfIndexes)):
					if counterOfIndexes[j].count(i) != 0:
						for k in range(len(counterOfIndexes[j])):
							maxValueCountInRow[counterOfIndexes[j][k]] -= 1
						counterOfIndexes[j] = [i]
						break;
			if (max(maxValueCountInRow) <= 0):
				break;

	tempList = counterOfIndexes
	counterOfIndexes = []
	i = 0
	for i in range (len(tempList)):
		counterOfIndexes.append(tempList[i][0])

	i = 0
	for i in range(len(counterOfIndexes)):
		if counterOfIndexes[i] != currentPermutation[i]:
			indexForSwap = currentPermutation.index(counterOfIndexes[i])
			matrix.swapRows(i, indexForSwap)
			currentPermutation[i], currentPermutation[indexForSwap] = currentPermutation[indexForSwap], currentPermutation[i]
	return matrix

# Метод, который преобразует матрицу к виду, пригодному для проведения расчетов
def toCalculativeForm(matrix):
	matrix = copy.copy(matrix)

	i = 0
	for i in range(matrix.rank()):
		matrix.absoluteTerms[i], matrix.coefficients[i][i] = -matrix.coefficients[i][i], -matrix.absoluteTerms[i]
		j = 0
		for j in range (matrix.rank()):
			matrix.coefficients[i][j] = matrix.coefficients[i][j]/matrix.absoluteTerms[i]
		matrix.absoluteTerms[i] = matrix.coefficients[i][i]
		matrix.coefficients[i][i] = 0
	return matrix

# Метод, который проводит итерационные вычесления
def calculateMatrix(matrix, accuracy, maxIterationsCount):
	matrix = copy.copy(matrix)
	previousResults = []
	resultVector = list(matrix.absoluteTerms)
	d = list(matrix.absoluteTerms)
	currentAccuracy = []

	i = 0
	for i in range(matrix.rank()):
		currentAccuracy.append(1)

	iterationCounter = 1
	for iterationCounter in range(maxIterationsCount+1):
		previousResults = list(resultVector)
		i = 0
		for i in range(matrix.rank()):
			j = 0
			newValue = 0
			for j in range(matrix.rank()):
				if (i != j):
					newValue += matrix.coefficients[i][j] * resultVector[j]
			newValue += d[i]
			resultVector[i] = newValue
		i = 0
		for i in range(matrix.rank()):
			currentAccuracy[i] = abs(resultVector[i] - previousResults[i])
		if (max(currentAccuracy) < accuracy):
			return [resultVector, iterationCounter, currentAccuracy]
	return None

# Точка входа в программу
inputResolver()
