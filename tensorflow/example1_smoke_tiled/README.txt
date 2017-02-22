Here's a quick guide how to run this example:

- first generate data by calling
	>>> manta manta_genSimData.py
- you can also use the gen_many_doublesim_outputs.sh to generate 10 data sets in one go (Linux)

- then use it to train a first model with tensorflow. This will take ca. 2 min. E.g., with: 
	>>>python modelNp.py out 0 trainingEpochs 1000 
  now you should have a trained model checkpoint for test_0001 in the data directory (../data by default).

- you can then use the model to generate an output with the command below. It assumes the trained
  model was generated in ../data/test_0001 , and that the sim data to apply it to is number 1000:
	>>>python modelNp.py out 1  fromSim 1000 toSim -1  load_model_test 1 load_model_no 0

- for an automated test run, you can call runTests.sh , which runs through a few different version
  but requires (possibly shortened) sim test data sets 2007, 2008

