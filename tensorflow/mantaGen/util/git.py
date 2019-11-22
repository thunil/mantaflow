#******************************************************************************
#
# MantaGen
# Copyright 2018 Steffen Wiewel, Moritz Becher, Nils Thuerey
#
# This program is free software, distributed under the terms of the
# Apache License, Version 2.0 
# http://www.apache.org/licenses/LICENSE-2.0
#
#******************************************************************************

from subprocess import check_output

def revision():
    return check_output(["git", "rev-parse", "--short", "HEAD"], universal_newlines=True).rstrip()

def status():
    return check_output(["git", "status", "-s"], universal_newlines=True)

def is_clean():
    return not bool(status())