import numpy as np

# Defaults
def defaults() -> dict:
    return {"target_steps":150, "dt":0.1, "pn0":0.1, "mn0":0.1, "pi0":0, "mi0":0, 
            "rp_n":0.02, "rp_i":0.02, "mup_n":0.1, "mup_i":0.1, "mum_n":0.1, "mum_i":0.1, "d":2,
            "qhp_n":5, "qhp_i":5, "qcp_n":1, "qcp_i":1, "qcm_n":1, "qcm_i":1, "qhm_n":1, "qhm_i":1,
            "cp_ni":0.0, "cp_in":0.0, "cm_ni":0.0, "cm_in":0.0,
            "alpha_nn":0.4, "beta_nn":0.4, "alpha_ii":0.4, "beta_ii":0.4, "alpha_ni":0, "beta_in":0, "alpha_in":0, "beta_ni":0}

# A single step in the simulation
def step(params: dict, last_step: dict, t: int) -> dict:
    if last_step is None:
        return {'pn': params['pn0'], 'mn' : params['mn0'], 'pi' : params['pi0'], 'mi' : params['mi0']}

    pn = last_step['pn']
    mn = last_step['mn']
    pi = last_step['pi']
    mi = last_step['mi']

    Qpnn = params['qhp_n'] * params['alpha_nn'] / params['d'] - params['qcp_n'] * params['beta_nn']
    Qpii = params['qhp_i'] * params['alpha_ii'] / params['d'] - params['qcp_i'] * params['beta_ii']
    Qpin = params['qhp_n'] * params['alpha_in'] / params['d'] - params['qcp_n'] * params['beta_ni']
    Qpni = params['qhp_i'] * params['alpha_ni'] / params['d'] - params['qcp_i'] * params['beta_in']

    Qmnn = params['qcm_n'] * params['beta_nn'] - params['qhm_n'] * params['alpha_nn'] / params['d']
    Qmii = params['qcm_i'] * params['beta_ii'] - params['qhm_i'] * params['alpha_ii'] / params['d']
    Qmin = params['qcm_n'] * params['beta_in'] - params['qhm_n'] * params['alpha_ni'] / params['d']
    Qmni = params['qcm_i'] * params['beta_ni'] - params['qhm_i'] * params['alpha_in'] / params['d']

    for sub_step in range(int(1/params['dt'])):
        pn_next = pn + (params['rp_n'] * pn + (mn * Qpnn + mi * Qpin) * pn / (mn + mi + pn / params['d'] + pi / params['d']) - params['cp_in'] * pn * pi - params['mup_n'] * pn**2) * params['dt']
        mn_next = mn + ((pn * Qmnn + pi * Qmin) * mn / (mn + mi + pn / params['d'] + pi / params['d']) - params['cm_in'] * mn * mi - params['mum_n'] * mn**2) * params['dt']
        pi_next = pi + (params['rp_i'] * pi + (mi * Qpii + mn * Qpni) * pi / (mn + mi + pn / params['d'] + pi / params['d']) - params['cp_ni'] * pn * pi - params['mup_i'] * pi**2) * params['dt']
        mi_next = mi + ((pi * Qmii + pn * Qmni) * mi / (mn + mi + pn / params['d'] + pi / params['d']) - params['cm_ni'] * mn * mi - params['mum_i'] * mi**2) * params['dt']
        pn = pn_next
        mn = mn_next
        pi = pi_next
        mi = mi_next
    
    return {'pn': pn, 'mn': mn, 'pi': pi, 'mi': mi}

# Simulation
def run(params: dict) -> dict:
    steps = []
    for t in np.arange(params['dt'], params['target_steps'], params['dt']):
        last_step = step(params, steps[-1] if len(steps) > 0 else None, t)
        if not last_step:
            break
        steps.append(last_step)
    return result(params, steps)

# Result
def result(params: dict, steps: list) -> dict:
    # For this model, we can just return the final step
    return steps[-1]