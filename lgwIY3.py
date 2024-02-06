import numpy as np

model_outcomes = ['Collapse', 'Spillover', 'Failure (resistance)', 'Rescue (resistance)', 'Failure', 'Gene swamping', 'Oscillations', 'Differential targeting', 'Suppression']

def defaults():
	return {
		's': 0.5,
		'c': 1,
		'h': 1,
		'd': 10,
		'r0': 2,
		'm': 0.01,
		'q1': 0.8,
		'q2': 0,
		'mu': 0,
		'nhej': 0,
		's_b': 0,
		'K': 1,
		'conversion': 'Gametic',
		'N_model': 'Beverton-Holt',
		'r_s': 'Eq. 6',
		'epsilon': 0.01,
		'target_steps': 100,
		'repeats': 1,
		'output': 'Outcome',
		'outcome': 'Suppression'
	}

def migration_network(params, graph):
	n = graph.shape[0]
	M = graph * params['m']
	M[range(n), range(n)] = 1 - np.sum(M, axis=1)
	return M

def initial(params):
	edges = np.array([[0, 1], [1, 0]])
	M = migration_network(params, edges)
	q = np.zeros(M.shape[0])
	p = np.zeros(M.shape[0])
	N = np.full(M.shape[0], params['K'] * (1 - params['epsilon']))
	q[0] = params['q1']
	q[1] = params['q2']
	return {'q1': q[0], 'q2': q[1], 'p1': p[0], 'p2': p[1], 'N1': N[0], 'N2': N[1], 'M': M}

def r_s(params, r0, s):
	if params['r_s'] == 'Fixed':
		return 1 / r0
	elif params['N_model'] == 'Logistic difference equation':
		if params['r_s'] == 'Eq. S1':
			return -s * params['d']
		else:
			return r0 - (1 - np.exp(-s * params['d'])) * (r0 + s)
	elif params['N_model'] == 'Beverton-Holt':
		if params['r_s'] == 'Eq. S1':
			return 1 - s * params['d']
		else:
			return r0 - (1 - np.exp(-s * params['d'])) * (r0 + s - 1)

def step(params, _step, t):
	if t == 0:
		return initial(params)
	M = np.array(_step['M'])
	s = params['s']
	c = params['c']
	h = params['h']
	m = params['m']
	K = params['K']
	s_b = params['s_b']
	q_m = np.array([
		((1 - m) * _step['q1'] + m * (_step['N2']/_step['N1']) * _step['q2']) * 1/(1 - m + m * (_step['N2']/_step['N1'])),
		((1 - m) * _step['q2'] + m * (_step['N1']/_step['N2']) * _step['q1']) * 1/(1 - m + m * (_step['N1']/_step['N2']))
	])
	p_m = np.array([
		((1 - m) * _step['p1'] + m * (_step['N2']/_step['N1']) * _step['p2']) * 1/(1 - m + m * (_step['N2']/_step['N1'])),
		((1 - m) * _step['p2'] + m * (_step['N1']/_step['N2']) * _step['p1']) * 1/(1 - m + m * (_step['N1']/_step['N2']))
	])
	N_m = np.dot([_step['N1'], _step['N2']], M)

	if params['mu'] > 0 or params['nhej'] > 0:
		for i in range(M.shape[0]):
			rate = 	params['mu'] * 2 * N_m[i] * ((1 - q_m[i])**2 + q_m[i] * (1 - q_m[i]) * (1 - c)) + \
				params['nhej'] * 2 * N_m[i] * c * q_m[i] * (1 - q_m[i])
			resistance_alleles = 0 if rate == 0 else np.random.poisson(rate)
			p_m[i] += resistance_alleles / (2 * N_m[i])

	w_bar = 2 * q_m * (1 - q_m - p_m) * (1 - h * s) + \
		2 * p_m * (1 - q_m - p_m) * (1 - s_b / 2) + \
		2 * p_m * q_m * (1 - (s * h + s_b * (1 - h))) + \
		q_m ** 2 * (1 - s) + \
		p_m ** 2 * (1 - s_b) + \
		(1 - q_m - p_m) ** 2

	q_n = 	(q_m * (1 - q_m - p_m) * (1 + c) * (1 - h * s) + \
		p_m * q_m * (1 - (s * h + s_b * (1 - h))) + \
		q_m ** 2 * (1 - s)) / w_bar

	p_n = 	(p_m * (1 - q_m - p_m) * (1 - s_b / 2) + \
		p_m * q_m * (1 - (s * h + s_b * (1 - h))) + \
		p_m ** 2 * (1 - s_b)) / w_bar

	r0 = params['r0'] if params['N_model'] == 'Beverton-Holt' else params['r0'] / 2

	R = 	(2 * q_m * (1 - q_m - p_m) * (1 - h * s) * r_s(params, r0, h * s) + \
		2 * p_m * (1 - q_m - p_m) * (1 - s_b / 2) * r_s(params, r0, s_b / 2) + \
		2 * p_m * q_m * (1 - (s * h + s_b * (1 - h))) * r_s(params, r0, (s * h + s_b * (1 - h))) + \
		q_m ** 2 * (1 - s) * r_s(params, r0, s) + \
		p_m ** 2 * (1 - s_b) * r_s(params, r0, s_b) + \
		(1 - q_m - p_m) ** 2 * r0) / w_bar

	if params['conversion'] == 'Zygotic':
		w_bar = 2 * q_m * (1 - q_m - p_m) * (1 - c) * (1 - h * s) + \
			2 * p_m * (1 - q_m - p_m) * (1 - s_b / 2) + \
			2 * p_m * q_m * (1 - (s * h + s_b * (1 - h))) + \
			(q_m ** 2 + 2 * q_m * (1 - q_m - p_m) * c) * (1 - s) + \
			p_m ** 2 * (1 - s_b) + \
			(1 - q_m - p_m) ** 2
		q_n = 	(q_m * (1 - q_m - p_m) * (1 - c) * (1 - h * s) + \
			p_m * q_m * (1 - (s * h + s_b * (1 - h))) + \
			(q_m ** 2 + 2 * q_m * (1 - q_m - p_m) * c) * (1 - s)) / w_bar
		R = 	(2 * q_m * (1 - q_m - p_m) * (1 - c) * (1 - h * s) * r_s(params, r0, h * s) + \
			2 * p_m * (1 - q_m - p_m) * (1 - s_b / 2) * r_s(params, r0, s_b / 2) + \
			2 * p_m * q_m * (1 - (s * h + s_b * (1 - h))) * r_s(params, r0, (s * h + s_b * (1 - h))) + \
			(q_m ** 2 + 2 * q_m * (1 - q_m - p_m) * c) * (1 - s) * r_s(params, r0, s) + \
			p_m ** 2 * (1 - s_b) * r_s(params, r0, s_b) + \
			(1 - q_m - p_m) ** 2 * r0) / w_bar

	if params['N_model'] == 'Logistic difference equation':
		N_n = N_m + R * N_m * (1 - (N_m / K)) 
	else:
		N_n = (R * K * N_m) / (K + (R - 1) * N_m)

	return {'q1': q_n[0], 'q2': q_n[1], 'p1': p_n[0], 'p2': p_n[1], 'N1': N_n[0], 'N2': N_n[1], 'M': M}

def determine_outcome(params, steps):
	gen_threshold = 0.5
	pop_threshold = 0.9
	res_threshold = 0.1
	q = [steps[-1]['q1'], steps[-1]['q2']]
	p = [steps[-1]['p1'], steps[-1]['p2']]
	N = np.array([[step['N1'], step['N2']] for step in steps]) / params['K']
	N_min = np.min(N, axis=0)
	N_diff = np.diff(N[:,0])
	min_max = N[np.any([np.all([np.concatenate(([0], N_diff)) > 0, np.concatenate((-N_diff, [0])) > 0], axis=0), np.all([np.concatenate(([0], N_diff)) < 0, np.concatenate((-N_diff, [0])) < 0], axis=0)], axis=0),0]
	oscillating = np.abs(np.diff(min_max))[-1] > 0.5 if len(min_max) > 1 else False
	if N_min[0] < pop_threshold and N_min[1] < pop_threshold:
		return 'Collapse'
	elif q[0] >= gen_threshold and q[1] >= gen_threshold:
		return 'Spillover'
	elif q[0] < gen_threshold and q[1] < gen_threshold and (p[0] >= res_threshold or p[1] >= res_threshold) and N_min[0] >= pop_threshold and N_min[1] >= pop_threshold:
		return 'Failure (resistance)'
	elif q[0] < gen_threshold and q[1] < gen_threshold and (p[0] >= res_threshold or p[1] >= res_threshold) and N_min[0] < pop_threshold and N[-1][0] >= pop_threshold:
		return 'Rescue (resistance)'
	elif q[0] < gen_threshold and q[1] < gen_threshold and N_min[0] < pop_threshold and N[-1][0] >= pop_threshold:
		return 'Gene swamping'
	elif q[0] < gen_threshold and q[1] < gen_threshold and N_min[0] >= pop_threshold and N_min[1] >= pop_threshold:
		return 'Failure'
	elif q[0] > 0 and q[1] < gen_threshold and oscillating:
		return 'Oscillations'
	elif q[1] < gen_threshold and N_min[0] < pop_threshold:
		return 'Suppression'
	elif q[0] >= gen_threshold and q[1] < gen_threshold and N_min[0] >= pop_threshold and N_min[1] >= pop_threshold:
		return 'Differential targeting'
	else:
		return 'undefined'

def result(params, steps):
	outcome = determine_outcome(params, steps)
	N1 = np.array([step['N1'] for step in steps])
	pop_threshold = 0.9
	return {
		'outcome': outcome,
		'max_suppression': 1 - np.min(N1) if outcome == 'Gene swamping' else 0,
		'sum_suppression': np.sum(1 - N1) / len(N1) if outcome == 'Gene swamping' else 0,
		'dur_suppression': np.sum(N1 < pop_threshold) / len(N1) if outcome == 'Gene swamping' else 0,
		'non_target_freq': np.max([step['q2'] for step in steps])
	}

def single_run(params):
	steps = []
	last_step = False
	for t in range(0, params['target_steps'] + 1):
		last_step = step(params, last_step, t)
		if not last_step:
			break
		steps.append(last_step)
	return result(params, steps)

def outcome_prop(params, res = 101):
	s = np.linspace(0, 1, res)
	c = np.linspace(0, 1, res)
	grid_x, grid_y = np.meshgrid(s, c)
	grid = np.array([grid_x.reshape(res * res), grid_y.reshape(res * res)]).T
	outcomes = np.array([single_run({**params, **{'s': x_y[0], 'c': x_y[1]}})['outcome'] for x_y in grid])
	outcome_labels, counts = np.unique(outcomes, return_counts=True)
	return {**{outcome: 0 for outcome in model_outcomes}, **{outcome_label: counts[i] / len(grid) for i, outcome_label in enumerate(outcome_labels)}}

def mcrit(params):
	m = 0
	incr = -2
	precision = -5
	while True:
		run_params = params.copy()
		run_params.update({'m': m})
		result = single_run(run_params)
		if result['outcome'] in ['Differential targeting', 'Suppression', 'Oscillations', 'Gene swamping']:
			m += 10**incr
		elif m != 0:
			m -= 10**incr
			if incr == precision:
				return {'mcrit': np.min([1, m])}
			incr -= 1
			m += 10**incr
		else:
			return {'mcrit': np.min([1, m])}

def run(params, **kwargs):
	if params['output'] == 'm*':
		return mcrit(params)
	elif params['output'] == 'Outcome proportion':
		return outcome_prop(params)
	elif params['repeats'] == 1:
		return single_run(params)
	else:
		runs = np.array([single_run(params)['outcome'] for i in range(params['repeats'])])
		return {'outcome_proportion': len(runs[runs == params['outcome']]) / params['repeats']}

