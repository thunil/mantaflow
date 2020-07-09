# a selection of test runs

# path to executable
MANTA=manta
MANTA=../mbbManta/buildSecondStd/manta

# optional GUI
GUI=
#GUI=-g

# 2D
${MANTA} create_dataset.py --name=TEST_smoke_simple -t smoke_simple -n 2 -s 10 -w 15 -d 2 ${GUI}

${MANTA} create_dataset.py --name=TEST_smoke_buoyant -t smoke_buoyant -n 2 -s 10 -w 15 -d 2 ${GUI}

${MANTA} create_dataset.py --name=TESTSIM -t flip -n 5 -s 10 -w 15 -d 2 --max_obstacle_count 3 --max_drop_count 5 ${GUI}

# 3D
${MANTA} create_dataset.py --name=TEST_smoke_simple -t smoke_simple -n 2 -s 10 -w 15 -d 3 ${GUI}

${MANTA} create_dataset.py --name=TEST_smoke_buoyant -t smoke_buoyant -n 2 -s 10 -w 15 -d 3 ${GUI}

${MANTA} create_dataset.py --name=TESTSIM -t flip -n 3 -s 10 -w 15 -d 3 --max_obstacle_count 0 --max_drop_count 3 ${GUI}
