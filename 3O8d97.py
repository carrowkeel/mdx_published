import numpy as np

def defaults():
	return {
		's': 0.6,
		'c': 0.8,
		'h': 0.5,
		'q0': 0.8,
		'm': 0.01,
		'target_steps': 100
	}

def generate_graph(params):
	M = np.array([[0, 1], [1, 0]])
	return M, np.array([])

def migration_network(params, graph):
	n = graph.shape[0]
	M = graph * params['m']
	M[range(n), range(n)] = 0
	M[range(n), range(n)] = 1 - np.sum(M, axis=1)
	return M

def initial(params):
	edges, nodes = generate_graph(params)
	M = migration_network(params, edges)
	q = np.zeros([1, M.shape[0]])
	q[:,0] = params['q0']
	nodes_q = np.insert(nodes, 2, np.maximum(q[0], 0.1), axis=1) if len(nodes) > 0 else np.ndarray([0,3])
	return {'q': q, 'q_target': params['q0'], 'q_non_target': 0, 'spillover': 0, 'M': M}

def step(params, _step, t):
	if t == 0:
		return initial(params)
	M = np.array(_step['M']) ## Migration network, generated in first step by "initial"
	s = np.array([[params['s'], params['s']]])
	c = params['c']
	h = params['h']
	s_c = (1 - s) * c
	s_n = 0.5 * (1 - h * s) * (1 - c)
	q_tilde = np.array([np.dot(np.array(_step['q'][0]), M)])
	w_bar = q_tilde**2 * (1 - s) + 2 * q_tilde * (1 - q_tilde) * (s_c + 2 * s_n) + (1 - q_tilde)**2
	q_tag = (q_tilde**2 * (1 - s) + 2 * q_tilde * (1 - q_tilde) * (s_c + s_n)) / w_bar
	return {'q': q_tag, 'q_target': np.mean(q_tag[:,0]), 'q_non_target': (np.mean(np.sum(q_tag, axis=1)) - np.mean(q_tag[:,0])) / (M.shape[0] - 1), 'spillover': np.mean([len(q_i[q_i >= 0.5]) / M.shape[0] for q_i in q_tag]), 'M': M}

def determine_outcome(params, last_step):
	target_q = last_step['q_target']
	non_target = last_step['q_non_target']
	M = np.array(last_step['M'])
	threshold = 0.5
	if target_q >= threshold and last_step['spillover'] <= 0:
		return 'Differential targeting'
	elif target_q < threshold and non_target < threshold:
		return 'Failure'
	else:
		return 'Spillover'

def result(params, last_step):
	outcome = determine_outcome(params, last_step)
	return {'outcome': outcome}

def single_run(params):
	last_step = None
	for t in range(0, params['target_steps'] + 1):
		last_step = step(params, last_step, t)
		if not last_step:
			break
	return result(params, last_step)

def run(params, **kwargs):
	return single_run(params)
