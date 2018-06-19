# ----------------------------------------------------------------------------
#
# MantaFlow fluid solver framework
# Copyright 2018 Kiwon Um, Nils Thuerey
#
# This program is free software, distributed under the terms of the
# Apache License, Version 2.0 
# http://www.apache.org/licenses/LICENSE-2.0
#
# Benchmark case for wave with fluid-implicit-particle (FLIP) simulation
#
# ----------------------------------------------------------------------------

import os
import numpy as np

measured = np.load(os.path.dirname(os.path.abspath(__file__))+'/benchmark_wave_measured.npz')
def g_tweaked(t, grav):
        r = np.interp(t, measured['t'], measured['r'])
        return vec3(-np.sin(-r)*grav, np.cos(-r)*grav, 0)

guion = True
pause = True

# default solver parameters
params               = {}
params['dim']        = 2                  # dimension
params['sres']       = 2                  # sub-resolution per cell
params['dx']         = 1.0/params['sres'] # particle spacing (= 2 x radius)
params['res']        = 10                 # reference resolution
params['len']        = 0.06               # reference length
params['bnd']        = 4                  # boundary cells
params['gref']       = -9.8               # real-world gravity
params['cgaccuracy'] = 1e-3               # cg solver's threshold
params['jitter']     = 0.5                # jittering particles
params['gfm']        = True               # 2nd order fluid-empty BC
params['fps']        = 30                 # frames per second
params['t_end']      = 5.0                # quit simulation: NOTE: don't set this greater than 8; measured data are not available
params['sdt']        = None               # fix timestep size

# scale unit in regard to the manta world
scaleToManta   = float(params['res'])/params['len']
# NOTE: the original test is slightly modified for the sake of convenience in discretization
# 900 x 508 x 62 -> 900 x 504 x 60 (= 60x(15 x 8.4 x 1))
params['gs'] = [15*params['res']+params['bnd']*2, round(8.4*params['res'])+params['bnd']*2, params['res']+params['bnd']*2 if params['dim']==3 else 1]
params['grav'] = params['gref']*scaleToManta

s             = Solver(name="FLIP", gridSize=vec3(params['gs'][0], params['gs'][1], params['gs'][2]), dim=params['dim'])
s.cfl         = 1
s.frameLength = 1.0/float(params['fps'])
s.timestepMin = 0
s.timestepMax = s.frameLength
s.timestep    = s.frameLength

# prepare grids and particles
gFlags  = s.create(FlagGrid)
gV      = s.create(MACGrid)
gVold   = s.create(MACGrid)
gP      = s.create(RealGrid)
gPhiSld = s.create(LevelsetGrid)

pp    = s.create(BasicParticleSystem)
pT    = pp.create(PdataInt)
pV    = pp.create(PdataVec3)
pVtmp = pp.create(PdataVec3)

mesh = s.create(name='mesh', type=Mesh) if (params['dim']==3 and guion) else None

paramSolvePressure = dict(flags=gFlags, vel=gV, pressure=gP, cgAccuracy=params['cgaccuracy'])
if params['gfm']:               # for the free-surface boundary condition
        gPhi    = s.create(LevelsetGrid)
        gIdxSys = s.create(ParticleIndexSystem)
        gIdx    = s.create(IntGrid)
        paramSolvePressure.update(phi=gPhi)

# boundary setup
gFlags.initDomain(params['bnd']-1)
bndBox = s.create(Box, p0=vec3(0), p1=vec3(params['gs'][0], params['gs'][1], params['gs'][2]))
inBox  = s.create(Box, p0=vec3(params['bnd'], params['bnd'], params['bnd'] if params['dim']==3 else 0), p1=vec3(params['gs'][0]-params['bnd'], params['gs'][1]-params['bnd'], (params['gs'][0]-params['bnd']) if params['dim']==3 else 1))
gPhiSld.join(bndBox.computeLevelset(notiming=True), notiming=True)
gPhiSld.subtract(inBox.computeLevelset(notiming=True), notiming=True)

# fluid setup: wave
pool_c = vec3(s.getGridSize().x*0.5, 0.093*0.5*scaleToManta+params['bnd'], s.getGridSize().z*0.5 if (params['dim']==3) else 0)
pool_s = vec3(s.getGridSize().x*0.5, 0.093*0.5*scaleToManta, s.getGridSize().z*0.5 if (params['dim']==3) else s.getGridSize().z)
fld = s.create(Box, center=pool_c, size=pool_s)
fld.applyToGrid(grid=gFlags, value=FlagFluid, respectFlags=gFlags)

begin = pp.pySize()
sampleShapeWithParticles(shape=fld, flags=gFlags, parts=pp, discretization=params['sres'], randomness=0, notiming=True)
end = pp.pySize()
pT.setConstRange(s=FlagFluid, begin=begin, end=end, notiming=True)

if guion:
        gui = Gui()
        gui.show()
        if pause: gui.pause()

while (s.timeTotal<params['t_end']): # main loop
        mapPartsToMAC(vel=gV, flags=gFlags, velOld=gVold, parts=pp, partVel=pV, ptype=pT, exclude=FlagEmpty)

        if params['sdt'] is None: s.adaptTimestep(gV.getMaxAbs())
        else: s.adaptTimestepByDt(params['sdt'])

        addGravityNoScale(flags=gFlags, vel=gV, gravity=g_tweaked(s.timeTotal, params['grav']))

        if params['gfm']:
                gridParticleIndex(parts=pp, flags=gFlags, indexSys=gIdxSys, index=gIdx)
                unionParticleLevelset(parts=pp, indexSys=gIdxSys, flags=gFlags, index=gIdx, phi=gPhi, radiusFactor=1.0)
                extrapolateLsSimple(phi=gPhi, distance=4, inside=True)

        setWallBcs(flags=gFlags, vel=gV)
        solvePressure(**paramSolvePressure)
        setWallBcs(flags=gFlags, vel=gV)
        extrapolateMACSimple(flags=gFlags, vel=gV)

        # update velocity (general update from FLIP and individual update for Lagrangian particles)
        flipVelocityUpdate(vel=gV, velOld=gVold, flags=gFlags, parts=pp, partVel=pV, flipRatio=0.97, ptype=pT, exclude=FlagEmpty)
        addForcePvel(vel=pV, a=vec3(0, params['grav'], 0), dt=s.timestep, ptype=pT, exclude=FlagFluid)

        # update position
        pp.getPosPdata(target=pVtmp)
        pp.advectInGrid(flags=gFlags, vel=gV, integrationMode=IntRK4, deleteInObstacle=False, ptype=pT, exclude=FlagEmpty)
        eulerStep(parts=pp, vel=pV, ptype=pT, exclude=FlagFluid)
        pp.projectOutOfBnd(flags=gFlags, bnd=params['bnd']+params['dx']*0.5, plane='xXyYzZ', ptype=pT)
        pushOutofObs(parts=pp, flags=gFlags, phiObs=gPhiSld, thresh=params['dx']*0.5, ptype=pT)

        # update velocity of the Lagrangian particles
        updateVelocityFromDeltaPos(parts=pp, vel=pV, x_prev=pVtmp, dt=s.timestep, ptype=pT, exclude=FlagFluid)

        # We don't need to solve the pressure for isolated cells.
        markFluidCells(parts=pp, flags=gFlags, ptype=pT)
        setPartType(parts=pp, ptype=pT, mark=FlagFluid, stype=FlagEmpty, flags=gFlags, cflag=FlagFluid)
        markIsolatedFluidCell(flags=gFlags, mark=FlagEmpty)
        setPartType(parts=pp, ptype=pT, mark=FlagEmpty, stype=FlagFluid, flags=gFlags, cflag=FlagEmpty)

        if params['dim']==3 and guion:
                gridParticleIndex(parts=pp, flags=gFlags, indexSys=gIdxSys, index=gIdx)
                unionParticleLevelset(parts=pp, indexSys=gIdxSys, flags=gFlags, index=gIdx, phi=gPhi, radiusFactor=1.0)
                extrapolateLsSimple(phi=gPhi, distance=4, inside=True)
                gPhi.createMesh(mesh)

        s.step()
