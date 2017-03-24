Here's a quick guide how to run this example. Note that the scripts here assume
that the executable is located in ".../mantaflow/build/", and most scripts were
primarily developed for Unix systems (Windows will work, but require path modifications)..

- first generate data by calling
	>>> manta manta_genSimData.py
- you can also use the tf_genManySims.py script to generate 10 data sets in one go.

- then use it to train a first model with tensorflow. This will take ca. 2 min. E.g., with: 
	>>>python tf_model_np.py out 0  fromSim 1000 toSim -1  trainingEpochs 1000 
  now you should have a trained model checkpoint for test_0001 in the data directory (../data by default).

- you can then use the model to generate an output with the command below. It assumes the trained
  model was generated in ../data/test_0001 , and that the sim data to apply it to is number 1000:
	>>>python tf_model_np.py out 1  fromSim 1000 toSim -1  loadModelTest 1 loadModelNo 0

- for an automated test run, you can call runTests.sh , which runs through a few different version
  but requires (possibly shortened) sim test data sets 2007, 2008

