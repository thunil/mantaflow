import uniio
import numpy as np
import scipy.misc
import os
import shutil
import copy
from shutil import copyfile
from random import randint, shuffle, random

basePath = '../data/'

# arrays of tiles, excluding discarded ones
tile_inputs_all = []
tile_outputs_all = []

# arrays of all tiles
tile_inputs_all_complete = []
tile_outputs_all_complete = []

combined_inputs_all = []
combined_outputs_all = []

combined_inputs_all_complete = []
combined_outputs_all_complete = []

tile_outputs_all_croped = []
img_inputs_all = []
img_outputs_all = []

tile_data = {
	'inputs_train': [],
	'inputs_test': [],
	'inputs_val': [],

	'outputs_train': [],
	'outputs_test': [],
	'outputs_val': []
}

paths = {
	'base': '',
	'sim': '',
	'frame': '',
	'frame_low_uni': '',
	'frame_high_uni': '',
	'tiles': '',
	'tile_low_uni': '',
	'tile_high_uni': ''
}

def updatePaths(simNo=None, frameNo=None, tileNo=None, tile_size_x=0, tile_size_y=0, overlapping=0, data_type=None):
	paths['base'] = basePath
	paths['sim'] = paths['base'] + 'sim_%04d/' % simNo
	paths['frame'] = paths['sim'] + 'frame_%04d/' % frameNo
	paths['frame_low_uni'] = paths['frame'] + data_type + '_low_%04d_%04d.uni' % (simNo, frameNo)
	paths['frame_high_uni'] = paths['frame'] + data_type + '_high_%04d_%04d.uni' % (simNo, frameNo)
	if overlapping > 0:
		paths['tiles'] = paths['frame'] + 'tiles_%02dx%02d_o%02d/' % (tile_size_x, tile_size_y, overlapping)
	else:
		paths['tiles'] = paths['frame'] + 'tiles_%02dx%02d/' % (tile_size_x, tile_size_y)
	paths['tile_low_uni'] = paths['tiles'] + data_type + '_low_%04d_%04d_%04d.uni' % (simNo, frameNo, tileNo)
	paths['tile_high_uni'] = paths['tiles'] + data_type + '_high_%04d_%04d_%04d.uni' % (simNo, frameNo, tileNo)

def uniToArray(uniPath, is_vel=False):
	head, content = uniio.readUni(uniPath)

	imageHeight = head['dimX']
	imageWidth  = head['dimY']
	if not is_vel:
		fixedArray = np.zeros((imageHeight, imageWidth), dtype='f')
	else:
		fixedArray = np.zeros((imageHeight, imageWidth, 3), dtype='f')
	if imageHeight == imageWidth:
		for x in range(0, imageHeight):
			for y in range(0, imageWidth):
				fixedArray[x][y] = content[(imageHeight - 1) - x][y][0]
	else:
		if not is_vel:
			fixedArray = np.reshape(content, [imageWidth, imageHeight])
			fixedArray = fixedArray[::-1]
		else:
			fixedArray = np.reshape(content, [imageWidth, imageHeight, 3])
			fixedArray = fixedArray[::-1]

	#print("%s img size  %d %d " % (uniPath, imageWidth, imageHeight ) )
	return fixedArray

def arrayToUni(input, savePath, motherUniPath, imageHeight, imageWidth, is_vel=False):
	head, _ = uniio.readUni(motherUniPath)
	head['dimX'] = imageWidth
	head['dimY'] = imageHeight

	if not is_vel:
		fixedArray = np.zeros((imageHeight, imageWidth), dtype='f')
		for x in range(0, imageHeight):
			for y in range(0, imageWidth):
				fixedArray[x][y] = input[(imageHeight - 1) - x][y]
	else:
		fixedArray = np.zeros((imageHeight, imageWidth, 3), dtype='f')
		for x in range(0, imageHeight):
			for y in range(0, imageWidth):
				fixedArray[x][y] = input[(imageHeight - 1) - x][y]

	uniio.writeuni(savePath, head, fixedArray)

def createTiles(input, imageHeight, imageWidth, tileHeight, tileWidth, overlapping=0):
	if ((imageHeight % tileHeight) != 0 |
		imageWidth % tileWidth != 0):
		print('Error: Image and tile size do not match.')

	tilesVertical = imageHeight // tileHeight
	tilesHorizontal = imageWidth // tileWidth

	tiles = []

	for currTileVert in range(0, tilesVertical):
		for currTileHor in range(0, tilesHorizontal):

			currTile = np.zeros((tileHeight + overlapping * 2, tileWidth + overlapping * 2), dtype='f')

			for currPixelVert in range(-overlapping, tileHeight + overlapping):
				for currPixelHor in range(-overlapping, tileWidth + overlapping):
					indexX = currPixelVert + (currTileVert * tileWidth)
					indexY = currPixelHor + (currTileHor * tileHeight)

					if (indexX >= 0) and (indexY >= 0) and (indexX < len(input)) and (indexY < len(input[0])):
						currTile[currPixelVert + overlapping][currPixelHor + overlapping] = input[indexX][indexY]
					else:
						currTile[currPixelVert + overlapping][currPixelHor + overlapping] = 0

			tiles.append(currTile)

	return tiles

def combineTiles(input, imageHeight, imageWidth, tileHeight, tileWidth):
	if ((imageHeight % tileHeight) != 0 |
			imageWidth % tileWidth != 0):
		print('Error: Image and tile size do not match.')

	tilesVertical = imageHeight // tileHeight
	tilesHorizontal = imageWidth // tileWidth

	resultArray = np.zeros((imageHeight, imageWidth), dtype='f')

	for currTileVert in range(0, tilesVertical):
		for currTileHor in range(0, tilesHorizontal):

			for currPixelVert in range(0, tileHeight):
				for currPixelHor in range(0, tileWidth):
					indexX = currPixelHor + (currTileHor * tileHeight)
					indexY = currPixelVert + (currTileVert * tileWidth)

					resultArray[indexX][indexY] = (input[currTileHor * tilesHorizontal + currTileVert])[currPixelHor][currPixelVert]

	return resultArray

def createPngFromUni(uniPath, save_path=''):
	if save_path == '':
		savePath = uniPath[:-3] + 'png'
	else:
		savePath = save_path

	array = uniToArray(uniPath)
	createPngFromArray(array, savePath)

def createPngFromArray(input, savePath):
	imageHeight = len(input)
	imageWidth = len(input[0])

	fixedArray = np.zeros((imageHeight, imageWidth), dtype='f')
	for x in range(0, imageHeight):
		for y in range(0, imageWidth):
			if (input[x][y] < 0.0):
				fixedArray[x][y] = 0.0
			else:
				fixedArray[x][y] = input[x][y]

	scipy.misc.toimage(fixedArray, cmin=0.0, cmax=1.0).save(savePath)

def createTestData(simNo, tileSize, lowResSize, upScalingFactor, overlapping=0, createPngs=False, with_vel=False, with_pos=False):
	frameNo = 0
	tileNo = 0

	data_type = 'density'

	updatePaths(simNo, frameNo, tileNo, tileSize, tileSize, overlapping, data_type)

	# for each frame: create tiles + folders
	# print(paths['tiles'])
	while os.path.exists(paths['frame']):
		# print('Creating tiles for sim %04d, frame %04d' % (simNo, frameNo))
		if not os.path.exists(paths['tiles']):
			os.makedirs(paths['tiles'])

		# create tiles

		if with_vel:
			lowArray = combine_density_with_vel(uniToArray(paths['frame_low_uni']), uniToArray(paths['frame_low_uni'].replace('density', 'vel'), is_vel=True))
			lowTiles = createTiles(lowArray, lowResSize * 2, lowResSize * 2, tileSize * 2, tileSize * 2, overlapping * 2)
		elif with_pos:
			lowArray = add_position_to_density_vel(uniToArray(paths['frame_low_uni']), uniToArray(paths['frame_low_uni'].replace('density', 'vel'), is_vel=True))
			lowTiles = createTiles(lowArray, lowResSize * 3, lowResSize * 3, tileSize * 3, tileSize * 3, overlapping * 3)
		else:
			lowArray = uniToArray(paths['frame_low_uni'])
			lowTiles = createTiles(lowArray, lowResSize, lowResSize, tileSize, tileSize, overlapping)
		highArray = uniToArray(paths['frame_high_uni'])
		highTiles = createTiles(highArray, lowResSize * upScalingFactor, lowResSize * upScalingFactor, tileSize * upScalingFactor, tileSize * upScalingFactor, overlapping)

		for currTile in range(0, len(lowTiles)):
			if with_vel:
				arrayToUni(lowTiles[currTile], paths['tile_low_uni'].replace('density', 'density_with_vel'), paths['frame_low_uni'], (tileSize + overlapping * 2) * 2, (tileSize + overlapping * 2) * 2)
			elif with_pos:
				arrayToUni(lowTiles[currTile], paths['tile_low_uni'].replace('density', 'density_with_vel_pos'), paths['frame_low_uni'], (tileSize + overlapping * 2) * 3, (tileSize + overlapping * 2) * 3)
			else:
				arrayToUni(lowTiles[currTile], paths['tile_low_uni'], paths['frame_low_uni'], tileSize + overlapping * 2, tileSize + overlapping * 2)
			arrayToUni(highTiles[currTile], paths['tile_high_uni'], paths['frame_high_uni'], tileSize * upScalingFactor, tileSize * upScalingFactor)
			if createPngs:
				createPngFromUni(paths['tile_low_uni'].replace('density', 'density_with_vel'))
				createPngFromUni(paths['tile_high_uni'])

			tileNo += 1
			updatePaths(simNo, frameNo, tileNo, tileSize, tileSize, overlapping, data_type)

		frameNo += 1
		tileNo = 0
		updatePaths(simNo, frameNo, tileNo, tileSize, tileSize, overlapping, data_type)

def createTestDataFromTo(simFrom, simTo, tileSize, lowResSize, upScalingFactor, overlapping=0, createPngs=False):
	for sim in range(simFrom, simTo + 1):
		createTestData(sim, tileSize, lowResSize, upScalingFactor, overlapping, createPngs)
		print('Created test data for sim %04d' % sim)

def selectRandomTiles(selectionSize, isTraining=True, is_combined=False, croped=False):
	# returns #selectionsize elements of inputs_train/test and outputs_train/test
	allInputs = tile_data['inputs_train']
	if not croped:
		allOutputs = tile_data['outputs_train']
	else:
		allOutputs = tile_outputs_all_croped

	if not isTraining:
		allInputs = tile_data['inputs_test']
		allOutputs = tile_data['outputs_test']

	selectedInputs = []
	selectedOutputs = []
	inputSize = len(allInputs) - 1
	for currNo in range(0, selectionSize):
		randNo = randint(0, inputSize)
		selectedInputs.append(allInputs[randNo])
		selectedOutputs.append(allOutputs[randNo])

	return selectedInputs, selectedOutputs

def combine_density_with_vel(low_tile_density, low_tile_velocity):
	# Combines one tile of density with one tile of velocity
	output_tile = np.zeros((len(low_tile_density) * 2, len(low_tile_density) * 2), dtype='f')

	for curr_pos_x in range(len(low_tile_density)):
		for curr_pos_y in range(len(low_tile_density)):
			output_tile[curr_pos_x * 2	][curr_pos_y * 2	] = low_tile_density[curr_pos_x][curr_pos_y]
			output_tile[curr_pos_x * 2 + 1][curr_pos_y * 2	] = low_tile_velocity[curr_pos_x][curr_pos_y][0]
			output_tile[curr_pos_x * 2	][curr_pos_y * 2 + 1] = low_tile_velocity[curr_pos_x][curr_pos_y][1]
			output_tile[curr_pos_x * 2 + 1][curr_pos_y * 2 + 1] = low_tile_velocity[curr_pos_x][curr_pos_y][2]

	createPngFromArray(output_tile, basePath + 'balabubbel.png')
	return output_tile

def add_position_to_density_vel(low_tile_density, low_tile_velocity):
	output_tile = np.zeros((len(low_tile_density) * 3, len(low_tile_density) * 3), dtype='f')
	pos = []

	for x in range(len(low_tile_density)):
		pos.append(np.linspace(0., 1., len(low_tile_density)))

	for curr_pos_x in range(len(low_tile_density)):
		for curr_pos_y in range(len(low_tile_density)):
			output_tile[curr_pos_x * 3	][curr_pos_y * 3	] = low_tile_density[curr_pos_x][curr_pos_y]
			output_tile[curr_pos_x * 3 + 1][curr_pos_y * 3	] = low_tile_velocity[curr_pos_x][curr_pos_y][0]
			output_tile[curr_pos_x * 3	][curr_pos_y * 3 + 1] = low_tile_velocity[curr_pos_x][curr_pos_y][1]
			output_tile[curr_pos_x * 3 + 1][curr_pos_y * 3 + 1] = low_tile_velocity[curr_pos_x][curr_pos_y][2]
			output_tile[curr_pos_x * 3	][curr_pos_y * 3 + 2] = pos[curr_pos_x][curr_pos_y]
			output_tile[curr_pos_x * 3 + 1][curr_pos_y * 3 + 2] = pos[curr_pos_y][curr_pos_x]

	# createPngFromArray(output_tile, basePath + 'blablubbel.png')
	return output_tile

def loadTestData(fromSim, toSim, densityMinimum, tileSizeLow, overlapping, partTrain=3, partTest=1, partVal=1, load_img=False, load_pos=False, load_vel=False, to_frame=200, low_res_size=64, upres=2):
	total_tiles_all_sim = 0
	discarded_tiles_all_sim = 0 
	data_type = 'density'

	for simNo in range(fromSim, toSim + 1):
		print('\nLoading sim %04d' % simNo)
		frameNo = 0
		tileNo = 0

		updatePaths(simNo, frameNo, tileNo, tileSizeLow, tileSizeLow, overlapping, data_type)
		print('from ' + paths['tiles'])
		# check if right tiles are available - create them if not
		if load_vel:
			look_for_tiles_path = paths['tile_low_uni'].replace('density', 'density_with_vel')
		elif load_pos:
			look_for_tiles_path = paths['tile_low_uni'].replace('density', 'density_with_vel_pos')
		else:
			look_for_tiles_path = paths['tile_low_uni']

		if not os.path.exists(paths['tiles']) or not os.path.exists(look_for_tiles_path):
			print('Could not find tiles for sim %04d. Creating new ones.' % simNo)
			createTestData(simNo, tileSizeLow, low_res_size, upres, overlapping, with_vel=load_vel, with_pos=load_pos)
			updatePaths(simNo, frameNo, tileNo, tileSizeLow, tileSizeLow, overlapping, data_type)
		totalTiles = 0
		discardedTiles = 0

		while frameNo < to_frame:
			if load_img:
				lowArray = uniToArray(paths['frame_low_uni'])
				highArray = uniToArray(paths['frame_high_uni'])
				img_inputs_all.append(lowArray.flatten())
				img_outputs_all.append(highArray.flatten())

			#print('Loading frame %04d' % frameNo)
			while os.path.exists(paths['tile_high_uni']):
				# check if tile is empty. If so, don't add it
				accumulatedDensity = 0.0
				if load_vel:
					lowTile = uniToArray(paths['tile_low_uni'].replace('density', 'density_with_vel'))
				elif load_pos:
					lowTile = uniToArray(paths['tile_low_uni'].replace('density', 'density_with_vel_pos'))
				else:
					lowTile = uniToArray(paths['tile_low_uni'])
				highTile = uniToArray(paths['tile_high_uni'])

				for value in lowTile.flatten():
					accumulatedDensity += value
				if accumulatedDensity > (densityMinimum * tileSizeLow * tileSizeLow):
					tile_inputs_all.append(lowTile.flatten())
					tile_outputs_all.append(highTile.flatten())
				else:
					# print('Discarded empty tile.')
					discardedTiles += 1

				tile_inputs_all_complete.append(lowTile.flatten())
				tile_outputs_all_complete.append(highTile.flatten())

				totalTiles += 1
				tileNo += 1
				updatePaths(simNo, frameNo, tileNo, tileSizeLow, tileSizeLow, overlapping, data_type)

			tileNo = 0
			frameNo += 1
			updatePaths(simNo, frameNo, tileNo, tileSizeLow, tileSizeLow, overlapping, data_type)

		print('Total Tiles: %d' % totalTiles)
		print('Discarded Tiles: %d' % discardedTiles)
		print('Used Tiles: %d' % (totalTiles - discardedTiles))
		total_tiles_all_sim += totalTiles
		discarded_tiles_all_sim += discardedTiles

	print('Used Tiles All Sim: %d' % (total_tiles_all_sim - discarded_tiles_all_sim))

	# split into train, test und val data

	parts_complete = partTrain + partTest + partVal
	end_train = (len(tile_inputs_all) // parts_complete) * partTrain
	end_test = end_train + (len(tile_inputs_all) // parts_complete) * partTest
	# shuffleSeed = random()
	# shuffle(tile_inputs_all, lambda: shuffleSeed)
	# shuffle(tile_outputs_all, lambda: shuffleSeed)
	# tile_data['inputs_train'] = tile_inputs_all
	# tile_data['outputs_train'] = tile_outputs_all
	#print( "Debug out   %f   %f  %f  " % (len(tile_inputs_all), len(tile_inputs_all) // parts_complete, parts_complete ) )
	print( "Ranges: len tile_inputs_all %d, end_train %d, end_test %d " % (len(tile_inputs_all), end_train, end_test) )
	tile_data['inputs_train'],  tile_data['inputs_test'],  tile_data['inputs_val']  = np.split(tile_inputs_all,  [end_train, end_test])
	tile_data['outputs_train'], tile_data['outputs_test'], tile_data['outputs_val'] = np.split(tile_outputs_all, [end_train, end_test])
	print('Training Sets: {}'.format(len(tile_data['inputs_train'])))
	print('Testing Sets:  {}'.format(len(tile_data['inputs_test'])))

def shuffle_in_unison(a, b):
	rng_state = np.random.get_state()
	np.random.shuffle(a)
	np.random.set_state(rng_state)
	np.random.shuffle(b)
	return a, b

def load_combined_test_data(fromSim, toSim, from_bad_sim, to_bad_sim, densityMinimum, tileSizeLow, overlapping, partTrain=3, partTest=1, partVal=1, load_img=False, load_vel=False, to_frame=200):
	loadTestData(fromSim, toSim, densityMinimum, tileSizeLow, overlapping, partTrain=partTrain, partTest=partTest,
				 partVal=partVal, load_img=load_img, load_vel=load_vel, to_frame=to_frame)
	good_data_len_all = len(tile_inputs_all)
	good_data_len_all_complete = len(tile_inputs_all_complete)
	loadTestData(from_bad_sim, to_bad_sim, densityMinimum, tileSizeLow, overlapping, partTrain=partTrain, partTest=partTest,
				 partVal=partVal, load_img=load_img, load_vel=load_vel, to_frame=to_frame)
	bad_data_len_all = len(tile_inputs_all) - good_data_len_all
	bad_data_len_all_complete = len(tile_inputs_all_complete) - good_data_len_all_complete

	# create labels
	inputs_permutated_all = []
	outputs_permutated_all = []

	for x in range(good_data_len_all):
		outputs_permutated_all.append(1)
	for x in range(bad_data_len_all):
		outputs_permutated_all.append(0)

	for x in range(good_data_len_all_complete):
		combined_outputs_all_complete.append(1)
	for x in range(bad_data_len_all_complete):
		combined_outputs_all_complete.append(0)


	# create input data
	for curr_tile in range(len(outputs_permutated_all)):
		inputs_permutated_all.append(np.concatenate((tile_inputs_all[curr_tile], tile_outputs_all[curr_tile])))
	for curr_tile in range(len(combined_outputs_all_complete)):
		combined_inputs_all_complete.append(np.concatenate((tile_inputs_all_complete[curr_tile], tile_outputs_all_complete[curr_tile])))

	# permutate and then add to true arrays
	inputs_permutated_all, outputs_permutated_all = shuffle_in_unison(inputs_permutated_all, outputs_permutated_all)
	for x in range(len(outputs_permutated_all)):
		combined_outputs_all.append(outputs_permutated_all[x])
	for curr_tile in range(len(inputs_permutated_all)):
		combined_inputs_all.append(inputs_permutated_all[curr_tile])

	parts_complete = partTrain + partTest + partVal
	end_train = (len(tile_inputs_all) // parts_complete) * partTrain
	end_test = end_train + (len(tile_inputs_all) // parts_complete) * partTest

	tile_data['inputs_train'], tile_data['inputs_test'], tile_data['inputs_val'] = np.split(combined_inputs_all, [end_train, end_test])
	tile_data['outputs_train'], tile_data['outputs_test'], tile_data['outputs_val'] = np.split(combined_outputs_all, [end_train, end_test])
	print('Training Sets: {}'.format(len(tile_data['inputs_train'])))
	print('Testing Sets: {}'.format(len(tile_data['inputs_test'])))

def debugOutputPngs(input, expected, output, tileSizeLow, tileSizeHigh, imageSizeLow, imageSizeHigh, path, imageCounter=0):
	expectedArray = []
	outputArray = []
	inputArray = []

	tileCounter = 0
	# imageCounter = 0
	for currTile in range(0, len(input)):
		expectedArray.append(np.reshape(expected[currTile], (tileSizeHigh, tileSizeHigh)))
		outputArray.append(np.reshape(output[currTile], (tileSizeHigh, tileSizeHigh)))
		# inputArray.append(np.reshape(input[currTile], (tileSizeLow, tileSizeLow)))

		tileCounter += 1

		if tileCounter == ((imageSizeLow // tileSizeLow) ** 2):
			imageCounter += 1
			tileCounter = 0

			expectedImage = combineTiles(expectedArray, imageSizeHigh, imageSizeHigh, tileSizeHigh, tileSizeHigh)
			outputImage = combineTiles(outputArray, imageSizeHigh, imageSizeHigh, tileSizeHigh, tileSizeHigh)
			# inputImage = combineTiles(inputArray, imageSizeLow, imageSizeLow, tileSizeLow, tileSizeLow)

			createPngFromArray(expectedImage, path + 'full_%04d_2.png' % imageCounter)
			createPngFromArray(outputImage, path + 'full_%04d_3.png' % imageCounter)
			# createPngFromArray(inputImage, path + 'full_%04d_1.png' % imageCounter)

			expectedArray = []
			outputArray = []
			inputArray = []

# Edited for croping: Added this method. Nearly the same as the original one, didn't want to ruin the original one.
def debugOutputPngs_for_croping(input, expected, output, tileSizeLow, tileSizeHigh, imageSizeLow, imageSizeHigh, path,
					imageCounter=0, cut_output_to=-1, tiles_in_image=-1):
	expectedArray = []
	outputArray = []
	# inputArray = []

	image_counter_init = imageCounter

	tileCounter = 0
	imageCounter = image_counter_init
	for currTile in range(0, len(output)):
		if cut_output_to == -1:
			outputArray.append(np.reshape(output[currTile], (tileSizeHigh, tileSizeHigh)))
		else:
			from_value = (tileSizeHigh - cut_output_to) / 2
			to_value = tileSizeHigh - from_value
			output_tile = np.reshape(output[currTile], (tileSizeHigh, tileSizeHigh))
			outputArray.append(output_tile[from_value: to_value, from_value: to_value])

		tileCounter += 1

		if tileCounter == (tiles_in_image):
			imageCounter += 1
			tileCounter = 0

			if cut_output_to == -1:
				outputImage = combineTiles(outputArray, imageSizeHigh, imageSizeHigh, tileSizeHigh, tileSizeHigh)
			else:
				outputImage = combineTiles(outputArray, imageSizeHigh, imageSizeHigh, cut_output_to, cut_output_to)

			createPngFromArray(outputImage, path + 'full_%04d_3.png' % imageCounter)
			outputArray = []

def combine_tiles_easy(input, image_size, tile_size):
	# Combines a batch of 1D tiles to a 2D image. Only quadratic.
	outputArray = []
	for currTile in range(0, len(input)):
		outputArray.append(np.reshape(input[currTile], (tile_size, tile_size)))

	return combineTiles(outputArray, image_size, image_size, tile_size, tile_size)

def create_bad_sim_data(output, from_sim, to_sim, to_frame=-1):
	# Output has to be image, quadratic and size length in power of 2. Input is copied.
	# Creates for every input sim a "bad output sim" folder, with a simulation number >1000
	# to_frame has to be the number of frames on which was trained

	# find next free sim path (1000+)
	folder_no = 1000
	path_addition = 'sim_%04d/' % folder_no
	while os.path.exists(basePath + path_addition):
		folder_no += 1
		path_addition = 'sim_%04d/' % folder_no
	simNo = folder_no

	image_size = len(output[0])

	print("\n")
	# Copy inputs (density and velocity) and save output
	for curr_input_sim_no in range(from_sim, to_sim + 1):
		sim_path = basePath + 'sim_%04d/' % simNo
		os.makedirs(sim_path)

		# for every frame in input
		curr_input_frame_no = 0
		while (os.path.exists(basePath + 'sim_%04d/' % curr_input_sim_no
									   + 'frame_%04d/' % curr_input_frame_no) and not curr_input_frame_no == to_frame):

			new_frame_path = sim_path + 'frame_%04d/' % curr_input_frame_no
			if not (os.path.exists(new_frame_path)):
				os.makedirs(new_frame_path)

			curr_input_frame_path_density = basePath + 'sim_%04d/' % curr_input_sim_no \
													 + 'frame_%04d/' % curr_input_frame_no \
													 + 'density_low_%04d_%04d.uni' % (curr_input_sim_no, curr_input_frame_no)
			curr_input_frame_path_vel	 = basePath + 'sim_%04d/' % curr_input_sim_no \
													 + 'frame_%04d/' % curr_input_frame_no \
													 + 'vel_low_%04d_%04d.uni' % (curr_input_sim_no, curr_input_frame_no)

			new_input_frame_path_density  = basePath + 'sim_%04d/' % simNo \
													 + 'frame_%04d/' % curr_input_frame_no \
													 + 'density_low_%04d_%04d.uni' % (simNo, curr_input_frame_no)
			new_input_frame_path_vel	  = basePath + 'sim_%04d/' % simNo \
													 + 'frame_%04d/' % curr_input_frame_no \
													 + 'vel_low_%04d_%04d.uni' % (simNo, curr_input_frame_no)

			copyfile(curr_input_frame_path_density, new_input_frame_path_density)
			copyfile(curr_input_frame_path_vel, new_input_frame_path_vel)

			new_output_frame_path_density = basePath + 'sim_%04d/' % simNo \
													 + 'frame_%04d/' % curr_input_frame_no \
													 + 'density_high_%04d_%04d.uni' % (simNo, curr_input_frame_no)

			arrayToUni(output[curr_input_frame_no], new_output_frame_path_density, curr_input_frame_path_density, image_size, image_size)

			curr_input_frame_no += 1

		print("Created bad output for sim %04d." % simNo)

		simNo += 1
	print("\n")

def create_png_for_sims(from_sim, to_sim, target_folder):
	# target_folder with / at the end
	frame_global_no = 0

	for simNo in range(from_sim, to_sim + 1):
		frameNo = 0

		updatePaths(simNo, frameNo, 0, 0, 0, 0, 'density')

		while os.path.exists(paths['frame']):
			createPngFromUni(paths['frame_low_uni'], target_folder + 'frame_%04d_2.png' % frame_global_no)
			createPngFromUni(paths['frame_high_uni'], target_folder + 'frame_%04d_3.png' % frame_global_no)

			frameNo += 1
			frame_global_no += 1
			updatePaths(simNo, frameNo, 0, 0, 0, 0, 'density')

def create_bad_mixed_data(from_sim, to_sim):
	bad_data_array = []

	for simNo in range(from_sim, to_sim + 1):
		frameNo = 0

		updatePaths(simNo, frameNo, 0, 0, 0, 0, 'density')

		while os.path.exists(paths['frame']):
			bad_data_array.append(uniToArray(paths['frame_high_uni']))

			frameNo += 1
			updatePaths(simNo, frameNo, 0, 0, 0, 0, 'density')

	np.random.shuffle(bad_data_array)
	create_bad_sim_data(bad_data_array, from_sim, to_sim)


# deep copy of sim data tree (only dens & vel, no tiles)
def copySimData(fromSim, toSim, to_frame=200):
	total_tiles_all_sim = 0
	discarded_tiles_all_sim = 0 
	dt_dens = 'density'
	dt_vel  = 'vel'

	#for simNo in range(fromSim, toSim + 1):
	simNo = fromSim
	if 1:
		print('\nCopying sim %04d' % simNo)
		frameNo = 0
		tileNo = 0

		updatePaths(simNo, frameNo, tileNo, 0, 0, 0, dt_dens )
		pathsFrom = copy.deepcopy(paths)
		updatePaths(toSim, frameNo, tileNo, 0, 0, 0, dt_dens )
		pathsTo = copy.deepcopy(paths)
		print('To ' + pathsTo['sim'])

		if not os.path.isdir( pathsTo['sim'] ):
			os.makedirs(pathsTo['sim'])

		while frameNo < to_frame:
			updatePaths(simNo, frameNo, tileNo, 0, 0, 0, dt_dens )
			pathsFrom = copy.deepcopy(paths)
			updatePaths(toSim, frameNo, tileNo, 0, 0, 0, dt_dens )
			pathsTo = copy.deepcopy(paths)

			if not os.path.isdir( pathsFrom['frame'] ):
				frameNo += 1
				continue;

			if not os.path.isdir( pathsTo['frame'] ):
				os.makedirs(pathsTo['frame'])

			if os.path.isfile( pathsFrom['frame_low_uni'] ):
				print('   copy ' + pathsFrom['frame_low_uni']+ " " + pathsTo['frame_low_uni'] )
				shutil.copyfile( pathsFrom['frame_low_uni'],  pathsTo['frame_low_uni'] )
			if os.path.isfile( pathsFrom['frame_high_uni'] ):
				print('   copy ' + pathsFrom['frame_high_uni']+ " " + pathsTo['frame_high_uni'] )
				shutil.copyfile( pathsFrom['frame_high_uni'], pathsTo['frame_high_uni'] )

			updatePaths(simNo, frameNo, tileNo, 0, 0, 0, dt_vel )
			pathsFrom = copy.deepcopy(paths)
			updatePaths(toSim, frameNo, tileNo, 0, 0, 0, dt_vel )
			pathsTo = copy.deepcopy(paths)

			if os.path.isfile( pathsFrom['frame_low_uni'] ):
				print('   copy vel ' + pathsFrom['frame_low_uni']+ " " + pathsTo['frame_low_uni'] )
				shutil.copyfile( pathsFrom['frame_low_uni'] , pathsTo['frame_low_uni'] )
			if os.path.isfile( pathsFrom['frame_high_uni'] ):
				print('   copy vel ' + pathsFrom['frame_high_uni']+ " " + pathsTo['frame_high_uni'] )
				shutil.copyfile( pathsFrom['frame_high_uni'] , pathsTo['frame_high_uni'] )

			frameNo += 1

	# copy sim done




# test

# create_bad_mixed_data(65, 90)
# create_png_for_sims(554, 555, '/home/sunija/manta_output/debug/')

# simNo, tileSize, lowResSize, upScalingFactor, overlapping=0, createPngs=False, with_vel=False, with_pos=False
# createTestData(0, 16, 64, 2, 0, with_pos=True)
# loadTestData(fromSim=0, toSim=0, densityMinimum=0.01, tileSizeLow=16, overlapping=0, load_vel=True)
# print(uniToArray('/home/sunija/manta_output/sim_0000/frame_0000/vel_low_0000_0000.uni', is_vel=True))
# print(combine_density_with_vel(uniToArray('/home/sunija/manta_output/sim_0000/frame_0000/density_low_0000_0000.uni'), uniToArray('/home/sunija/manta_output/sim_0000/frame_0000/vel_low_0000_0000.uni', is_vel=True)))
# updatePaths(1, 0, 0, 8, 8, 'density')
# createTestData(0, 0, 64, 2, 4, False)
# createTestDataFromTo(54, 56, 8, 64, 2)
# createTestDataFromTo(simFrom=0, simTo=0, tileSize=8, lowResSize=64, overlapping=0, upScalingFactor=2, createPngs=True)

# testPath = '/home/sunija/manta_output/sim_0000/frame_0064/density_low_0000_0064.uni'
# testArray = uniToArray(testPath)
# tiles = createTiles(testArray, 64, 64, 8, 8, 4)
#
# for currTile in range(len(tiles)):
#	 createPngFromArray(tiles[currTile], basePath + 'testTile_%04d.png' % currTile)




