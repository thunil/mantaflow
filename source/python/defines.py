################################################################################
#
# MantaFlow fluid solver framework
# Copyright 2011 Tobias Pfaff, Nils Thuerey 
#
# This program is free software, distributed under the terms of the
# GNU General Public License (GPL) 
# http://www.gnu.org/licenses
#
# Defines some constants for use in python subprograms
#
#################################################################################

# mantaflow conventions
Real = float

# grid flags
FlagFluid = 1
FlagObstacle = 2
FlagEmpty = 4
FlagStick = 128
FlagReserved = 256

# integration mode
IntEuler = 0
IntRK2 = 1
IntRK4 = 2
