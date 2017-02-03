## FEM - uzawa algorithm
#modules
import numpy as np
import pylan as pn
from scipy import sparse
from scipy.sparse import linalg
import time
import matplotlib.pyplot as plt

#path and functions
path = '/Users/milan/python/FiniteElements/uzawa/'
exec(open(path+'call_functions.py').read())

## define domain
nx,ny = 35,35
param = set_param(nx,ny)
param['path'] = path

## triangulation of domain
tic = time.time()
param = triangulation(param)
print('triangulation in '+str(time.time() - tic)[:5]+'s.')

## nonlinear integral matrices
tic = time.time()
op = nlintmat(param)
print('nonlinear intmat in '+str(time.time() - tic)[:5]+'s.')

## stiffness&mass matrix
tic = time.time()
op = stiffmass(param,op)
print('stiffness&mass in '+str(time.time() - tic)[:5]+'s.')

## divergence matrix
tic = time.time()
op = divmat(param,op)
print('divmat in '+str(time.time() - tic)[:5]+'s.')

## initial conditions and forcing
tic = time.time()
vars = set_vars(param,op)
print('initial and forcing in '+str(time.time() - tic)[:5]+'s.')

## uzawa with or without optimal stepwidth, linear/nonlinear
eval = dict()

tic = time.time()
vars,eval = call_uzawa_alg(param,op,vars,eval)
print('Newton done in '+str(time.time() - tic)[:5]+'s.')

## test
print('\ntest the two equations...')
#the following values should be small in order of the equations to be fulfilled
print('mom eq. '+format(np.linalg.norm(op['A'].dot(vars['v']) \
+ op['Bt'].dot(vars['p']) - vars['f']),'.2e'))
print('con eq. '+format(np.linalg.norm(op['B'].dot(vars['v'])),'.2e'))

## projection back to normal space
vars = reproj(param,vars)

## error calculation
eval = calc_error(param,vars,eval)

## plotting
plot_uzawa(param,vars,eval)

