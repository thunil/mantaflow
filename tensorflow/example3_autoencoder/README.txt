
... OLD update

This is a brief overview and getting-started guide for the source code of tempoGAN.

Main source code directories
.../tensorflow/datagen: contains scene files for generating 2D/3D training data
.../tensorflow/tools:   contains necessary tools for inputs, outputs, neural networks operation, and etc.
.../tensorflow/GAN:     contains the tempoGAN model.
And two data directories
.../tensorflow/2ddata_sim:     contains the training and test data
.../tensorflow/2ddata_gan:     outputs will be written here

First generate simulation data with the following command. e.g.
	>>> manta ../datagen/gen_sim_data.py basePath ../2ddata_sim/ reset 1 savenpz 1
Also generate the sample plume data (gen_sim_2006.py for 2D, gen_sim_3006.py for 3D) into the 2ddata_sim directory.

Then you can start to train a GAN using
	>>> python example_run_training.py
This trains four models, for a quick test disable the later three.

After you trained a GAN model, you can use the model to generate new results, e.g.
	>>>python example_run_output.py
By default, these examples run on simulation "2006".

Note: all the commands above are just examples, please check parameters when running them (esp. paths, simulation ID ranges etc.)
