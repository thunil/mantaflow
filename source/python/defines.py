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

# some defines to make C code and scripts more alike...
false = False
true  = True
Vec3  = vec3
Vec4  = vec4
Vec3Grid = VecGrid

# grid flags
FlagFluid    = 1
FlagObstacle = 2
FlagEmpty    = 4
FlagOutflow  = 16
FlagStick    = 128
FlagReserved = 256
# and same for FlagGrid::CellType enum names:
TypeFluid    = 1
TypeObstacle = 2
TypeEmpty    = 4
TypeStick    = 128
TypeReserved = 256

# integration mode
IntEuler = 0
IntRK2   = 1
IntRK4   = 2


