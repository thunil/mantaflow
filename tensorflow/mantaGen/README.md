## mantagen

# Agreed on naming scheme (180823):

Directories should be of the form: .../ID/sim_XXXXXX/field_YYYYYY.npz

Use consecutive numbered ranges for X & Y with format %06d.

ID - global identifier string, eg, "liquid_2d" , "smoke_3d_low" +  "smoke_3d_hi"
sim - fixed, dont change
field - name of quantity, e.g., "density" , "velocity"

Per field npz file:
Single grid content, 2D [Y,X,c] , 3D [Z,Y,X,c]

# sample call

.../manta create_dataset.py --name=TEST --type smoke_simple -n 3  -s 10 -w 5 --seed 10 --dimension 2 --resolution_x 40 --resolution_y 40 

will generate 2D smoke datasets in the directories:
datasets/TEST/sim_000000,
datasets/TEST/sim_000001, 
...,
in total "-n" times.

For 3D you can use:

.../manta create_dataset.py --name=TEST --type smoke_simple -n 3  -s 10 -w 5 --seed 10 --resolution_x 40 --resolution_y 40 --resolution_z 40

