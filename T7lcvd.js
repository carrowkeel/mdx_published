/* Generated by ModelRxiv, editing could affect compatibility */

export const defaults = () => ({
	'target_steps': 100,
	's': '0.6',
	'c': '0.8',
	'h': '0.5',
	'q': '0.5'
});

export const step = (params, _step, t) => {
	const variables = Object.assign({}, params, _step);
	variables['sn'] = 0.5 * (1 - variables['c']) * (1 - variables['h'] * variables['s']);
	variables['sc'] = variables['c'] * (1 - variables['s']);
	variables['w'] = variables['q']**2 * (1 - variables['s']) + 2 * variables['q'] * (1 - variables['q']) * (2 * variables['sn'] + variables['sc']) + (1 - variables['q'])**2;
	variables['q'] = (variables['q']**2 * (1 - variables['s']) + 2 * variables['q'] * (1 - variables['q']) * (variables['sn'] + variables['sc'])) / variables['w'];
	return {'q': variables['q']};
};

export const result = (params, steps) => {
	const variables = Object.assign({}, params, steps[steps.length - 1]);
	return {'q': variables['q']};
};

export const run = (params) => {
	let prev = undefined;
	for(let t = 0;t <= params.target_steps;t++) {
		const r = step(params, prev, t);
		if (!r)
			break;
		prev = r;
	}
	return result(params, [prev]);
};