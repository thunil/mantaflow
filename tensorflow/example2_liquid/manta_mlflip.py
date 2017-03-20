# ----------------------------------------------------------------------------
#
# MantaFlow fluid solver framework
# Copyright 2017 Kiwon Um, Nils Thuerey
#
# This program is free software, distributed under the terms of the
# GNU General Public License (GPL)
# http://www.gnu.org/licenses
#
# Fluid implicit particle (FLIP) with Machine Learning
#
# ----------------------------------------------------------------------------

import os, sys, argparse, pickle

def path_to_frame(outdir, frame):
    return '{}/{:05d}/'.format(outdir, frame)

def save_frame(outdir, frame, saving_funcs):
    if (outdir is None) or (params['frame_saved']==frame): return

    path = path_to_frame(outdir, frame)
    os.path.isdir(path) or os.makedirs(path)
    for save, name in saving_funcs: save(path+name, notiming=True)

    params['frame_saved'] = frame
    print('Frame #{} was saved in {}\n'.format(frame, path))

parser = argparse.ArgumentParser(description='FLIP with ML', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument(      '--load',   default='/tmp/tf/',  help='path to the trained tensorflow model directory')
parser.add_argument('-w', '--window', default=1, type=int, help='window size for sampling features; 1 (default) means 3x3, 2 means 5x5, so on.')
pargs = parser.parse_known_args()[0]

onphi  = True
ongeom = False

if pargs.load is None: sys.exit('You have to specify the path to the trained model.')
pargs.load = os.path.normpath(pargs.load)

with open(pargs.load + '/run_args.pickle', 'rb') as f: tfopt = pickle.load(f)
with open(pargs.load + '/scale.pickle', 'rb') as f: scale = pickle.load(f)

import numpy as np
dtype_real = np.float32         # NOTE: if double precision, use float64
dtype_int  = np.int32           # NOTE: if int in C is 64bits, use int64
np.random.seed(1)

import tensorflow as tf
tf_sess = tf.InteractiveSession()

import tf_network
dlayers = list(map(int, tfopt['dnet'].split('-')))
mlayers = list(map(int, tfopt['mnet'].split('-')))
dact    = list(map(tf_network.parse_act, tfopt['dact'].split('-')))
mact    = list(map(tf_network.parse_act, tfopt['mact'].split('-')))
x       = tf.placeholder(tf.float32, shape=[None, dlayers[0]], name='x-input')
y_,  y  = tf_network.build_network(dlayers, dact, input_x_holder=x, bn=tfopt['bn'], is_training=False, scope='detector/')[1:]
y2_, y2 = tf_network.build_network(mlayers, mact, input_x_holder=x, bn=tfopt['bn'], is_training=False, scope='modifier/')[1:]
if tfopt['mve']:
    sd  = tf_network.build_network(mlayers, mact, input_x_holder=x, input_y_holder=y2_, bn=tfopt['bn'], is_training=False, scope='modifier_var/')[2]

tf_saver = tf.train.Saver()
modelfile = pargs.load + '/model.ckpt'
tf_saver.restore(tf_sess, modelfile)
print('Pre-trained model {} loaded\n'.format(modelfile))

import manta as mt
mt.tFluid    = FlagFluid
mt.tObstacle = FlagObstacle

nogui       = False
pause       = True
output      = None
#output      = '/tmp/manta-mlflip'
savingFuncs = []

# default solver parameters
params                = {}
params['dim']         = 2                  # dimension
params['sres']        = 2                  # sub-resolution per cell
params['dx']          = 1.0/params['sres'] # particle spacing (= 2 x radius)
params['res']         = 50                 # reference resolution
params['len']         = 1.0                # reference length
params['bnd']         = 2                  # boundary cells
params['grav']        = 0                  # applied gravity (mantaflow scale); recomputed later
params['gref']        = -9.8               # real-world gravity
params['jitter']      = 0.1                # jittering particles
#params['stref']       = 0.073              # surface tension (reference scale [m]; e.g., 0.073)
params['cgaccuracy']  = 1e-3               # cg solver's threshold
params['fps']         = 30
params['t_end']       = 6.0
params['sdt']         = None
params['frame_saved'] = -1
params['frame_last']  = -1

scaleToManta = float(params['res'])/params['len']
params['gs']    = [params['res']*1.5+params['bnd']*2, params['res']+params['bnd']*2, params['res']+params['bnd']*2 if params['dim']==3 else 1]
params['grav']  = params['gref']*scaleToManta
#params['stens'] = params['stref']*scaleToManta

s             = Solver(name='MLFLIP', gridSize=vec3(params['gs'][0], params['gs'][1], params['gs'][2]), dim=params['dim'])
s.cfl         = 1
s.frameLength = 1.0/float(params['fps'])
s.timestepMin = 0
s.timestepMax = s.frameLength
s.timestep    = s.frameLength

# prepare grids and particles
gFlags   = s.create(FlagGrid)
gV       = s.create(MACGrid)
gVold    = s.create(MACGrid)
gP       = s.create(RealGrid)
gPhi     = s.create(LevelsetGrid)
gFlagTmp = s.create(FlagGrid)

pp      = s.create(BasicParticleSystem)
gIdxSys = s.create(ParticleIndexSystem)
gIdx    = s.create(IntGrid)

pT    = pp.create(PdataInt)
pV    = pp.create(PdataVec3)
pVtmp = pp.create(PdataVec3)
pVtm2 = pp.create(PdataVec3)
pVtm3 = pp.create(PdataVec3)
pItmp = pp.create(PdataInt)

mesh = s.create(name='mesh', type=Mesh) if (params['dim']==3 and not nogui) else None

savingFuncs.append([pp.save, 'particles.uni'])
savingFuncs.append([pV.save, 'particlesVel.uni'])
savingFuncs.append([pT.save, 'particlesType.uni'])

fv_N_axis = 2*pargs.window + 1
fv_N_stn  = fv_N_axis*fv_N_axis*(fv_N_axis if params['dim']==3 else 1)
fv_N_row  = params['dim']*fv_N_stn + (fv_N_stn if onphi else 0) + (fv_N_stn if ongeom else 0)
fv_vscale = params['len']/float(params['res'])

# boundary
gFlags.initDomain(params['bnd']-1)

# fluid dam setup
a = vec3(params['res']*0.2+params['bnd'], params['res']*0.2+params['bnd'], params['res']*0.2+params['bnd'] if (params['dim']==3) else 0)
b = vec3(params['res']*0.2, params['res']*0.2, params['res']*0.2 if (params['dim']==3) else params['gs'][2])
fld = s.create(Box, center=a, size=b)

begin = pp.pySize()
sampleShapeWithParticles(shape=fld, flags=gFlags, parts=pp, discretization=params['sres'], randomness=params['jitter'], notiming=True)
end = pp.pySize()
pT.setConstRange(s=FlagFluid, begin=begin, end=end, notiming=True)
markFluidCells(parts=pp, flags=gFlags, ptype=pT, exclude=FlagObstacle)

gui = None
if not nogui:
    gui = Gui()
    gui.show()
    if pause: gui.pause()

if output:
    save_frame(output, s.frame, savingFuncs)
    with open(output+'/params.pickle', 'wb') as f: pickle.dump(params, f)    

stats = {'candidate': 0, 'decision': 0, 'reverted': 0, 'splashed': 0}
np_pTimer = np.zeros(pp.pySize(), dtype=dtype_real)
while (s.timeTotal<params['t_end']): # main loop

    mapPartsToMAC(vel=gV, flags=gFlags, velOld=gVold, parts=pp, partVel=pV, ptype=pT, exclude=FlagEmpty)
    if params['sdt'] is None: s.adaptTimestep(gV.getMax())
    else: s.adaptTimestepByDt(params['sdt'])

    addGravityNoScale(flags=gFlags, vel=gV, gravity=vec3(0, params['grav'], 0))

    gridParticleIndex(parts=pp, flags=gFlags, indexSys=gIdxSys, index=gIdx)
    unionParticleLevelset(parts=pp, indexSys=gIdxSys, flags=gFlags, index=gIdx, phi=gPhi, radiusFactor=1.0)

    setWallBcs(flags=gFlags, vel=gV)
    solvePressure(flags=gFlags, vel=gV, pressure=gP, cgAccuracy=params['cgaccuracy'], phi=gPhi)
    setWallBcs(flags=gFlags, vel=gV)
    extrapolateMACSimple(flags=gFlags, vel=gV)

    # BEGIN: machine learning part >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    # get candidate particles
    gFlagTmp.copyFrom(gFlags)
    extendRegion(flags=gFlagTmp, region=FlagEmpty, exclude=FlagObstacle, depth=1)
    pItmp.copyFrom(pT)
    setPartType(parts=pp, ptype=pItmp, mark=0, stype=FlagEmpty, flags=gFlagTmp, cflag=FlagEmpty|FlagFluid) # already individual? then, kill
    setPartType(parts=pp, ptype=pItmp, mark=FlagEmpty, stype=FlagFluid, flags=gFlagTmp, cflag=FlagEmpty)   # mark surface particles
    candidate = np.zeros(pp.pySize(), dtype=dtype_int)
    particleDataImplToNumpyInt(n=candidate, p=pItmp)
    candidate = (candidate==FlagEmpty)
    N_candidate = np.count_nonzero(candidate)
    stats['candidate'] += N_candidate

    # extract features -> numpy array
    inputs_c = np.zeros(pp.pySize()*fv_N_row, dtype=dtype_real)
    off_feature = 0
    extractFeatureVel(fv=inputs_c, N_row=fv_N_row, off_begin=off_feature, p=pp, vel=gV, scale=fv_vscale, ptype=pItmp, exclude=FlagObstacle|FlagFluid, window=pargs.window)
    off_feature = params['dim']*fv_N_stn
    if onphi:
        extractFeaturePhi(fv=inputs_c, N_row=fv_N_row, off_begin=off_feature, p=pp, phi=gPhi, scale=1.0, ptype=pItmp, exclude=FlagObstacle|FlagFluid, window=pargs.window)
        off_feature += fv_N_stn
    if ongeom:
        extractFeatureGeo(fv=inputs_c, N_row=fv_N_row, off_begin=params['dim']*fv_N_stn, p=pp, flag=gFlags, scale=1.0, ptype=pItmp, exclude=FlagObstacle|FlagFluid, window=pargs.window)
        off_feature += fv_N_stn

    inputs_c = inputs_c.reshape((-1, fv_N_row))[candidate]

    # approximate using tf: detection and modification
    if tfopt['mve']:    dtct_c, dv_c, appx_s_c = tf_sess.run([y, y2, sd], feed_dict={x: inputs_c})
    else:               dtct_c, dv_c           = tf_sess.run([y, y2],     feed_dict={x: inputs_c})
    if tfopt['nosmax']: dtct_c = (np.random.random(size=(N_candidate, 1))<dtct_c).reshape(-1)
    else:               dtct_c = (np.argmax(dtct_c, axis=1)==0) # NOTE: when 0th value is larger, it means splashing

    # mark new decision for the splashing particles (as FlagEmpty)
    N_splashing = np.count_nonzero(dtct_c)
    decision = np.zeros(pp.pySize(), dtype=dtype_int)
    decision[candidate] = dtct_c
    decision[(decision==1)] = FlagEmpty
    decision[(decision!=FlagEmpty)] = FlagFluid
    stats['decision'] += N_splashing

    if params['dim']==2:
        dv_c = np.append(dv_c, np.zeros((N_candidate, 1), dtype=dtype_real), axis=1)
        if tfopt['mve']: appx_s_c = np.append(appx_s_c, np.zeros((N_candidate, 1), dtype=dtype_real), axis=1)

    if tfopt['mve']: dv_c += appx_s_c*np.random.normal(size=(N_candidate, 3))

    dv = np.zeros((pp.pySize(), 3), dtype=dtype_real)
    dv[candidate] = dv_c

    # END: machine learning part <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    # update velocity; general update from FLIP
    flipVelocityUpdate(vel=gV, velOld=gVold, flags=gFlags, parts=pp, partVel=pV, flipRatio=0.97, ptype=pT, exclude=FlagObstacle|FlagEmpty)

    # update velocity; individual update for Lagrangian particles
    addForcePvel(vel=pV, a=vec3(0, params['grav'], 0), dt=s.timestep, ptype=pT, exclude=FlagObstacle|FlagFluid)

    # 0. let's predict splash in framestep, which is the same with the training data's timestep
    pp.getPosPdata(target=pVtmp)
    c_dt = s.timestep
    s.timestep = s.frameLength

    # 1. mark fluid region after moving without splash-correction
    gFlagTmp.copyFrom(gFlags)
    pp.advectInGrid(flags=gFlagTmp, vel=gV, integrationMode=IntRK4, deleteInObstacle=False, ptype=pT, exclude=FlagObstacle|FlagEmpty)
    eulerStep(parts=pp, vel=pV, ptype=pT, exclude=FlagObstacle|FlagFluid)
    markFluidCells(parts=pp, flags=gFlagTmp, ptype=pT, exclude=FlagObstacle|FlagEmpty)
    extendRegion(flags=gFlagTmp, region=FlagEmpty, exclude=FlagObstacle, depth=1)

    # 2. try to move splashing particles only, so check if it's really splashing; revert the wrong decisions
    pp.setPosPdata(source=pVtmp)
    pVtm2.copyFrom(pV)
    numpyToParticleDataImplVec3(p=pVtm3, n=dv.reshape(-1, 1))
    pVtm2.add(pVtm3)
    numpyToParticleDataImplInt(p=pItmp, n=decision)
    eulerStep(parts=pp, vel=pVtm2, ptype=pItmp, exclude=FlagObstacle|FlagFluid)
    setPartType(parts=pp, ptype=pItmp, mark=FlagFluid, stype=FlagEmpty, flags=gFlagTmp, cflag=FlagFluid|FlagObstacle) # empty -> fluid if they are not acturally splashing.
    particleDataImplToNumpyInt(n=decision, p=pItmp)

    # 3. final decision and velocity modification
    stats['splashed'] += np.count_nonzero(decision==FlagEmpty)
    dv[(decision!=FlagEmpty)] = 0
    dv *= scale['modvel']
    np_pTimer[(decision==FlagEmpty)] = s.frameLength # set judgement timer

    # 4. roll-back
    s.timestep = c_dt
    pp.setPosPdata(source=pVtmp)

    # mark splashing particles and modify the velocities so that they can flow individually
    numpyToParticleDataImplInt(p=pItmp, n=decision)
    pT.setConstIntFlag(s=FlagEmpty, t=pItmp, flag=FlagEmpty)
    numpyToParticleDataImplVec3(p=pVtm2, n=dv.reshape(-1, 1))
    pV.add(pVtm2)

    # update position
    pp.advectInGrid(flags=gFlags, vel=gV, integrationMode=IntRK4, deleteInObstacle=False, ptype=pT, exclude=FlagObstacle|FlagEmpty)
    eulerStep(parts=pp, vel=pV, ptype=pT, exclude=FlagObstacle|FlagFluid)
    pp.projectOutOfBnd(flags=gFlags, bnd=params['bnd']+params['dx']*0.5)
    markFluidCells(parts=pp, flags=gFlags, ptype=pT, exclude=FlagObstacle|FlagEmpty)
    markIsolatedFluidCell(flags=gFlags, mark=FlagEmpty)

    # update velocity of the Lagrangian particles
    updateVelocityFromDeltaPos(parts=pp, vel=pV, x_prev=pVtmp, dt=s.timestep, ptype=pT, exclude=FlagObstacle|FlagFluid)

    # NOTE: We don't need to solve the pressure for isolated cells.
    setPartType(parts=pp, ptype=pT, mark=FlagFluid, stype=FlagEmpty, flags=gFlags, cflag=FlagFluid) # empty -> fluid if they enter again.
    setPartType(parts=pp, ptype=pT, mark=FlagEmpty, stype=FlagFluid, flags=gFlags, cflag=FlagEmpty) # fluid -> empty if they escape

    # keep the valid splashing judgements; the particles may still stay inside the flow (due to a small timestep size)
    curr = np.zeros(pp.pySize(), dtype=dtype_int)
    particleDataImplToNumpyInt(n=curr, p=pT)
    np_pTimer = np_pTimer - s.timestep
    np_pTimer[(np_pTimer<=0.0)] = 0.0 # time-over of a splashing judgement
    keep = np.zeros(pp.pySize(), dtype=dtype_int)
    keep[(curr==FlagFluid)&(np_pTimer>0.0)] = FlagEmpty # keep valid judgement
    numpyToParticleDataImplInt(p=pItmp, n=keep)
    pT.setConstIntFlag(s=FlagEmpty, t=pItmp, flag=FlagEmpty) # judgement is still valid -> splashing (empty)

    s.step()

    if output: save_frame(output, s.frame, savingFuncs)
    
tf_sess.close()

stats['fraction'] = float(stats['decision'])/stats['candidate']
stats['reverted'] = stats['decision']-stats['splashed']
if output:
    with open(output+'/run_tf_stats.txt', 'w') as f:
        for key,value in sorted(stats.items()): f.write('{}: {}\n'.format(key, value))
