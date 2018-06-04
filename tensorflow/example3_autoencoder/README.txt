This example contains a res-net autoencoder with data augmentation.

First compile mantaflow, then generate data, e.g., by running the following call ten times:
	>>> manta ./manta_genSimData2.py basePath ../data/ reset 1 savenpz 1 gui 0  dim 2  

Then you can train a model, e.g., with
	>>> python  tf_resnet.py   out 0  trainingEpochs 1000  alwaysSave 1 fromSim 1000 toSim 1009 outputInterval 200 genTestImg 1   dataDim 2  batchSize 16  genModel gen_resnetSm   basePath  ../data/  loadPath ../data/  data_fraction 0.2   dataAugmentation 1 rot 1     saveInterval 200  simSize 64  

You can e.g. modify this test by using a larger network (genModel gen_resnet), or by changing data generation and training run to 3D.

