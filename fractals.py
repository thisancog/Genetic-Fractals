import os, sys, getopt, datetime, math, json
import numpy as np
from scipy.signal import convolve2d
from PIL import Image, ImageDraw
from PIL.ImageChops import multiply
from random import randint, choice

# Options #
numShades = 30				# number of shades / rules per ruleset (minimum: 2)
populationSize = 2000			# number of rulesets to initially create
generationSize = 200			# number of rulesets to create per generation
targetConvergance = 100			# fitness per Pixel to reach to define a ruleset converged to the target
					# -> real convergance = targetConvergance * 3 ** targetIteration
maxGenerations = 50000			# maximum amount of generations to breed
breedingSize = 20			# number of best adopted rulesets to breed with each other
sparedIndividuals = 4			# number of top individuals to keep from the previous generation
numNewIndividuals = 5			# number of completely random rulesets to add with each generation
defaultFitnessType = 'average'		# fitness measurement definition: 'average' or 'squaredDiff'

targetIteration = 4			# which iteration the target image should be comapred to (resolution: 3^n x 3^n)
startShade = 0				# the shade to start the first iteration with
shadeSwapProb = 0.2			# probability to swap two shades in a rule
genomePenetrationRatio = 0.5		# ratio of inherited genes from parent A to parent B
mutationProb = 0.0001			# probability of a shade to randomly mutate
allowCloning = False			# states if a rule can be copied/cloned if there are always have to be two parent rules
discardClones = False			# states if identical individuals in the breeding pool should be discarded

saveFolder = 'results'			# name of the folder to store results
progressFrequency = 1			# the number of generations to pass before saving a progress picture
autosaveFrequency = 100			# the number of generations to pass before autosaving the rulesets
tintColor = False			# the color all shades should be tinted in, tuple of rgb values [0-255],
					# False for no shading
					# before: (143, 201, 255)

# Global variables #
rulesets = []
currentGeneration = 0
currentFitness = None
startTime = datetime.datetime.now()
targetFile = ''
target = []
targetConv = []


# Grayscale and posterize the target image to n shades and parse as list of numbers
def load_target_image():
	global target, targetConv

	size = 3 ** targetIteration
	img = Image.open(targetFile)

	if (img.size[0] != size or img.size[1] != size):
		print("Invalid image. You need to specify a target image that is %d x %d pixels large." % (size, size))
		print("fractals.py -target <targetfile>")
		sys.exit(2)

	gray = img.convert("L")
	gray.save(saveFolder + '/target.png')

	target = np.array(list(gray.getdata())).reshape(size, size)
	target = np.floor_divide(target, 255 / (numShades - 1))
	targetConv = conv_iteration(target)

	print("Target file %s parsed." % targetFile)


# Generate a new, random ruleset for all shades
def new_ruleset():
	ruleset = { 'fitness': 0, 'rules': [] }
	for i in range(numShades):
		rule = []
		for c in range(9):
			rule.append(randint(0, numShades - 1))
		ruleset['rules'].append(rule)
	return ruleset


# Populate list of rulesets with random rules
def init_rulesets():
	global rulesets

	for i in range(populationSize):
		rulesets.append(new_ruleset())

# Create new rules from 
def breed_rulesets(ruleseta, rulesetb):
	ruleseta = ruleseta['rules']
	rulesetb = rulesetb['rules']
	newRuleset = []

	# completely swap two shades in ruleset A
	if (randint(0, 100) < shadeSwapProb * 100):
		ruleA = randint(0, numShades - 1)
		ruleB = ruleA
		while (ruleA == ruleB):
			ruleB = randint(0, numShades - 1)

		ruleTemp = ruleseta[ruleA]
		ruleseta[ruleA] = ruleseta[ruleB]
		ruleseta[ruleB] = ruleTemp

		for i in range(numShades):
			for n in range(9):
				if (ruleseta[n] == ruleA):
					ruleseta[n] == ruleB
				elif (ruleseta[n] == ruleB):
					ruleseta[n] = ruleA
		newRuleset = ruleseta

	else:
		for i in range(numShades):
			if (randint(0, 100) < genomePenetrationRatio * 100):
				newRuleset.append(ruleseta[i])
			else:
				newRuleset.append(rulesetb[i])

	# mutate some of the numbers
	if mutationProb * numShades * 9 >= 1:
		for i in range(int(mutationProb * numShades * 9)):
			firstIndex = randint(0, numShades - 1)
			secondIndex = randint(0, 8)
			newRuleset[firstIndex][secondIndex] = randint(0, numShades - 1)

	return { 'rules': newRuleset, 'fitness': 0 }


# Iterate a rule n times
def iterate_ruleset(ruleset, iteration = targetIteration):
	dimension = 3 ** iteration
	initial = startShade * math.floor(255 / (numShades - 1))
	result = Image.new('L', (dimension, dimension), initial)

	for i in range(iteration):
		oldBlocks = 3 ** (iteration - i)
		newBlocks = 3 ** (iteration - i - 1)

		for block in range(9 ** i):
			pY = oldBlocks * math.floor(block / (3 ** (i)))
			pX = oldBlocks * (block % (3 ** (i)))

			oldShade = result.getpixel((pX + 1, pY + 1))
			rule = ruleset['rules'][oldShade]

			for n in range(9):
				newY = pY + newBlocks * math.floor(n / 3)
				newX = pX + newBlocks * (n % 3)

				draw = ImageDraw.Draw(result)
				draw.rectangle([newX, newY, newX + newBlocks, newY + newBlocks], rule[n], None)
				del draw

	return np.array(result.getdata()).reshape(dimension, dimension)


def measure_fitness():
	if (defaultFitnessType == 'average'):
		for ruleset in rulesets:
			ruleset['fitness'] = get_fitness_average(ruleset)
	elif (defaultFitnessType == 'squaredDiff'):
		for ruleset in rulesets:
			ruleset['fitness'] = get_fitness_squaredDiff(ruleset)


# Compare iteration of a rule with the target
# with smaller values indicating better fitness
def get_fitness(ruleset, fitnessType = defaultFitnessType):
	if (fitnessType == 'average'):
		fitness = get_fitness_average(ruleset)
	elif (fitnessType == 'squaredDiff'):
		fitness = get_fitness_squaredDiff(ruleset)

	ruleset['fitness'] = fitness
	return fitness


# squaredDiff approach: squared difference
def get_fitness_squaredDiff(ruleset):
	return np.sum(np.square(target - iterate_ruleset(ruleset)))


def get_fitness_average(ruleset):
	iteration = conv_iteration(iterate_ruleset(ruleset))
	return np.sum(np.square(targetConv - iteration))


def conv_iteration(iteration):
	# https://stackoverflow.com/a/30082326
	# Pad around the input array to take care of boundary conditions
	arr_pad = np.lib.pad(iteration, (1,1), 'wrap')

	R,C = np.where(iteration==0)   # Row, column indices for zero elements in input array
	N = arr_pad.shape[1]     # Number of rows in input array

	offset = np.array([-N, -1, 1, N])
	idx = np.ravel_multi_index((R+1,C+1),arr_pad.shape)[:,None] + offset

	arr_out = iteration.copy()
	arr_out[R,C] = arr_pad.ravel()[idx].sum(1)/4
	return arr_out




def store_rulesets():
	filename = saveFolder + '/save-' + str(currentGeneration) + '.json'

	data = {
		'rulesets': rulesets,
		'currentGenerations': currentGeneration,
		'currentFitness': currentFitness,

		'numShades': numShades,
		'populationSize': populationSize,
		'generationSize': generationSize,
		'targetConvergance': targetConvergance,
		'maxGenerations': maxGenerations,
		'breedingSize': breedingSize,
		'sparedIndividuals': sparedIndividuals,
		'numNewIndividuals': numNewIndividuals,
		'defaultFitnessType': defaultFitnessType,
		'targetIteration': targetIteration,
		'startShade': startShade,
		'shadeSwapProb': shadeSwapProb,
		'genomePenetrationRatio': genomePenetrationRatio,
		'mutationProb': mutationProb,
		'allowCloning': allowCloning,
		'discardClones': discardClones,
		'saveFolder': saveFolder,
		'progressFrequency': progressFrequency,
		'autosaveFrequency': autosaveFrequency,
		'tintColor': tintColor,
		'startTime': startTime.isoformat(),
		'targetFile': targetFile,
		'target': target.tolist(),
		'targetConv': targetConv.tolist()
	}

	with open(filename, 'w', encoding="utf8") as outfile:
		json.dump(data, outfile, indent = 4, separators = (',', ': '))


def new_generation():
	global rulesets

	rulesets = rulesets[:breedingSize]
	newGeneration = rulesets

	# get rid of clones
	if discardClones:
		newGeneration = [rulesets[0]]
		i = 1

		while (len(newGeneration) < breedingSize and i < len(rulesets)):
			if (rulesets[i] != newGeneration[i - 1]):
				newGeneration.append(rulesets[i])
		
	# keep best individuals from previous generation
	newGeneration = newGeneration[:sparedIndividuals]

	# add random new individuals
	randomNewIndividuals = []
	for i in range(numNewIndividuals):
		randomNewIndividuals.append(new_ruleset())
	newGeneration.extend(randomNewIndividuals)


	# breed parents and fill generation with their children
	while (len(newGeneration) < generationSize):
		parenta = choice(rulesets)
		parentb = choice(rulesets)

		if not allowCloning:
			while (parentb == parenta):
				parentb = choice(rulesets)

		newIndividual = breed_rulesets(parenta, parentb)
		newGeneration.append(newIndividual)

	rulesets = newGeneration


def breed_generation():
	global rulesets, currentFitness, currentGeneration

	measure_fitness()
	rulesets.sort(key = lambda x: x['fitness'])
	currentFitness = rulesets[0]['fitness']

	if (currentGeneration > 0 and currentGeneration % progressFrequency == 0):
		construct_image(rulesets[0], 'progress-' + str(currentGeneration), targetIteration)
		print('Breeding generation %d, current fitness: %d' % (currentGeneration, currentFitness))

	if (currentGeneration > 0 and currentGeneration % autosaveFrequency == 0):
		store_rulesets()

	new_generation()
	currentGeneration += 1


def construct_image(ruleset, filename, iteration = targetIteration):
	data = iterate_ruleset(ruleset, iteration)
	img = Image.fromarray((data * 255 / (numShades - 1)).astype(np.uint8), 'L')

	if (tintColor != False):
		img = multiply(img.convert('RGB'), Image.new('RGB', img.size, tintColor))
	img.save(saveFolder + '/' + filename + '.png')


def construct_animation():
	return


# log reuslts to the console
def report_results():
	stopTime = datetime.datetime.now()
	delta = stopTime - startTime

	if (currentGeneration == maxGenerations):
		print("Training is over, maximum generations of %d reached." % maxGenerations)
	else:
		print("Training was stopped, %d generations reached." % currentGeneration)
	print("Resulting fitness: %d" % currentFitness)
	print("Time elapsed: %s days, %s hours, %s minutes, %s seconds" % (delta.days, delta.seconds//3600, (delta.seconds//60) % 60, (delta.seconds//3600) % 60))


# run script and check for system arguments
def main():
	global targetFile, saveFolder

	try:
		opts, args = getopt.getopt(sys.argv[1:], "t:", ["target="])
	except getopt.GetoptError:
		print("Invalid call. Run script like this:")
		print("fractals.py -target <targetfile>")
		sys.exit(2)
	for opt, arg in opts:
		if opt in ("-target", "--target", "-t", "--t"):
			targetFile = arg

	if (targetFile == ''):
		print("Invalid call. You need to specify a target file.")
		print("fractals.py -target <targetfile>")
		sys.exit(2)

	saveFolder = os.path.join(os.path.dirname(__file__), saveFolder + '-' + os.path.splitext(targetFile)[0])
	if not os.path.exists(saveFolder):
		os.makedirs(saveFolder)

	load_target_image()
	init_rulesets()


	while (currentGeneration < maxGenerations or currentFitness > targetConvergance * 3 ** targetIteration):
		try:
			breed_generation()
		except (KeyboardInterrupt, SystemExit):
			store_rulesets()
			report_results()
			raise
		except:
			print(sys.exc_info()[1])

	store_rulesets()
	report_results()


def test():
	global saveFolder, rulesets

	saveFolder = os.path.join(os.path.dirname(__file__), saveFolder + '-' + os.path.splitext(targetFile)[0])
	if not os.path.exists(saveFolder):
		os.makedirs(saveFolder)

	rulesets.append({ 'fitness': 0, 'rules': ([0,0,0,0,1,0,0,0,0], [1,1,1,1,1,1,1,1,1]) })
	construct_image(rulesets[0], 'sierpinski', 6)


if __name__ == "__main__":
#	test()
	main()

