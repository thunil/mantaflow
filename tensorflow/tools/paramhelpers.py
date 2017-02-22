#
# Helpers for handling command line parameters and the like
# example: path = getParam("path", "path.uni")
#
import sys
import os
import shutil

# global for checking used params
paramUsed = []


# ======================================================================================================================
# read parameters

#! check for a specific parameter, note returns strings, no conversion; not case sensitive! all converted to lower case
def getParam(name, default):
	global paramUsed
	while( len(paramUsed)<len(sys.argv) ):
		paramUsed.append(0);
	for iter in range(1, len(sys.argv)):
		#if(iter <  len(sys.argv)-1): print("Param %s , used %d, val %s" %( sys.argv[iter].lower(), paramUsed[iter] , sys.argv[iter+1]) ); # debug
		if(sys.argv[iter].lower() == name.lower()) and (iter+1<len(paramUsed)):
			paramUsed[iter] = paramUsed[iter+1] = 1;
			return sys.argv[iter+1];
	return default;

def checkUnusedParams():
	global paramUsed
	err = False;
	for iter in range(1, len(sys.argv)):
		if(paramUsed[iter]==0):
			print("Error: param %d '%s' not used!" % (iter,sys.argv[iter]) );
			err = True;
	if err:
		exit(1);

# others

def backupSources(name):
	sceneFile = sys.argv[0];
	shutil.copyfile( sceneFile, '%s_scene.py' % (name) )
	# double check path, might not exist everywhere
	#srcFile = "xxx.cpp";
	#if (os.path.isfile(srcFile)):
	#	shutil.copyfile( srcFile, '%s_of.cpp' % (name) )

