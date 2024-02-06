import numpy as np

def run(params):
    last_step = None
    steps = []
    for t in range(params['target_steps']):
        last_step = step(params, last_step, t)
        if not last_step:
             break
        steps.append(last_step)
    return result(params, steps)

def defaults():
    return {
        'H0': 300,
        'KH': 1000,
        'lambdaH': 0.3,
        'P0': 100,
        'su': 0.2,
        'Ebase': 20,
        'c': 40,
        'minl': 0.4,
        'a': 1,
        'target_steps': 200
    }

def F(EE, P, H, KH, minl, Ebase, c):
    if EE <= Ebase:
        l = 1
    else:
        l = minl + (1 - minl) * np.exp(-(EE - Ebase))
    if P > 0:
        eggs = min([EE, l * c * (H / KH)])
    else:
        eggs = 0
    return eggs

def step(params, last_step, t):
    if t == 0:
        return {
            'H': params['H0'],
            'P': params['P0'],
            'EE': params['Ebase']
        }
    H = last_step['H']
    P = last_step['P']
    EE = last_step['EE']
    FF = F(EE, P, H, params['KH'], params['minl'], params['Ebase'], params['c'])
    H = max(H + H * params['lambdaH'] * (1 - H / params['KH']) - FF * P, 100)
    P = np.floor(params['su'] * FF * P)
    EE = np.floor(EE - params['a'] if FF < EE else EE + params['a'])
    return {
        'H': H,
        'P': P,
        'EE': EE,
        'FF': FF
    }

def result(params, steps):
    if len(steps) != params['target_steps']:
        time = len(steps)
    else:
        time = params['target_steps']
    last_step = steps[-1]
    return {
        'final_H': last_step['H'],
        'final_P': last_step['P'],
        'final_EE': last_step['EE'],
        'final_FF': last_step['FF'],
        'steps_to_time': time
    }
