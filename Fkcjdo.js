
const mean = arr => arr.length === 0 ? 0 : arr.reduce((sum, value) => sum + value, 0) / arr.length;

export const defaults = () => ({
	'target_steps': 100,
	'q0': '0.01',
	's': '0.1',
	'h': '0.5',
	'N': '100',
	'L': '10',
	'rep': '100'
});

export const step = (params, _step, t) => { // frequency in loci in reps
	const variables = Object.assign({}, params, _step);
	if (t === 0)
		variables['q'] = Array.from(new Array(variables['rep'])).map(() => new Array(variables['L']).fill(variables['q0']));
	const selection = variables['q'].map(locus => locus.map(q => (q**2*(1-variables['s']) + q*(1-q)*(1-variables['h']*variables['s']))/(q**2*(1-variables['s']) + 2*q*(1-q)*(1-variables['h']*variables['s']) + (1-q)**2)));
	const drift = selection.map(locus => locus.map(q => Array.from(new Array(2*variables['N'])).map(v => Math.random() <= q).reduce((s,q) => s + q, 0)/(2*variables['N'])));
	return {'q': drift, 'purged': mean(drift.map(locus => locus.filter(q => q === 0).length / variables['L'])), 'fixed': mean(drift.map(locus => locus.filter(q => q === 1).length / variables['L'])), 'fitness': drift.map(locus => locus.map(q => q**2*(1-variables['s']) + 2*q*(1-q)*(1-variables['h']*variables['s']) + (1-q)**2).reduce((f,w) => f*w, 1))};
};

export const result = (params, steps) => {
	const variables = Object.assign({}, params, steps[steps.length - 1]);
	return {'purged': variables['purged']};
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
