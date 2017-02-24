import os
from subprocess import call

scenefile = "./manta_genSimData_ra.py"
mantaexe = "../../build/manta"
my_env = os.environ.copy()
my_env["MANTA_DISABLE_UI"]="1"

for i in range(1,11):
	# call(args, *, stdin=None, stdout=None, stderr=None, shell=False, env=os.environ, ...)
	# call([mantaexe, scenefile], env=my_env) # it is not reproducible without a specific npSeed
	# call([mantaexe, scenefile, "npSeed","%d"%i], env=my_env)
	call([mantaexe, scenefile, "npSeed","%d"%i, "basePath","../data2D_a/"], env=my_env)