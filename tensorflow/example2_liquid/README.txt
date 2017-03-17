Here's a quick guide how to run this example:

- first generate simulation data by calling
	>>> manta manta_flip.py
- then extract the training data for tensorflow with
	>>> manta manta_gendata.py

- then use it to train a first model with:
	>>>python tf_train.py /tmp/manta-flip/training_data/ 
  the last parameter contains the extracted data from the previous step, you can also enter multiple ones here

- once the model is trained, you can use it in a new sim with
	>>> manta manta_mlflip.py

note that all data is stored in /tmp/ by default. so better copy it somewhere else if you want to keep it...

