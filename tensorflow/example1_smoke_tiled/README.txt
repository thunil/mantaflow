Here's a quick guide how to run this example. Note that the scripts here assume
that the executable is located in ".../mantaflow/build/", and most scripts were
primarily developed for Unix systems (Windows will work, but require path
modifications)..

- first generate data by calling
	>>> manta manta_genSimData.py
	# windows: ../../build/Release/manta.exe manta_genSimData.py
- you can also use the tf_genManySims.py script to generate 10 data sets in one go.

- then use it to train a first model with tensorflow. This will take ca. 2 min.
  E.g., with: 
	>>>python tf_train_pressure.py  out 0 fromSim 1000 toSim -1 trainingEpochs 10000
  now you should have a trained model checkpoint for test_0001 in the data
  directory (../data by default).

- you can then use the model to generate an output with the command below. It
  assumes the trained model was generated in ../data/test_0001 , and that the
  sim data to apply it to the same sim data set:
	>>>python tf_train_pressure.py  out 1 fromSim 1000 toSim -1 loadModelTest 1 

- for an automated test run, you can call runTests.sh , which runs through a
  few different version but requires (possibly shortened) sim test data 
  sets 2007, 2008

---

This example also contains scripts to train networks for pressure instead of
density.  To run these tests, use manta_genSimData_pressure.py and
tf_train_pressure.py instead of the corresponding scritps above. Note that the
data gen script outputs frames in intervals, and thus definitely should be run
multiple times. The examples below assume 6 times. (Also note that density and
pressure data directories are not interchangeable.)

tf_train_pressure needs to be called with additional parameters, e.g. for training:
	>>>python tf_train_pressure.py  out 0 fromSim 1000 toSim 1005 useVelocities 1 outName pressure bWidth 1 trainingEpochs 10000
and for generating outputs:
	>>>python tf_train_pressure.py  out 1 fromSim 1006 toSim -1   useVelocities 1 outName pressure bWidth 1 loadModelTest 1 brightenOutput 10

