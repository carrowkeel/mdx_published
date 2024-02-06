
const range = (start,end) => Array.from(Array(end-start)).map((v,i)=>i+start);
const sum = arr => arr.reduce((a,v) => a+v, 0);
const mean = arr => arr.length === 0 ? 0 : arr.reduce((a,v) => a+v, 0) / arr.length;
const binomial = (trials, p) => sum(Array.from(new Array(trials)).map(() => Math.random() < p));

export const defaults = () => ({
	'target_steps': 100,
	'p0': [0.1, 0.1],
	's': [0.1, 0.1],
	'h': 0.5,
	'm': 0.01,
	'N': 10,
	'rep': 100
});

export const step = (params, _step, t) => { // allele frequencies in repeats
	const variables = Object.assign({}, params, _step);
	if (t === 0) {
		variables['p'] = Array.from(new Array(params['rep'])).map(v => variables['p0'].slice());
		const hs = range(0, 2).map(i => variables['p'].map(pop => 2*pop[i]*(1-pop[i])));
		const ht = variables['p'].map(pop => 2*mean(pop)*(1-mean(pop)));
		const fst = ht.map((ht,i) => ht === 0 ? 0 : 1 - mean([hs[0][i], hs[1][i]]) / ht);
		const w = range(0, 2).map(i => variables['p'].map(p => p[i]**2 * (1-variables['s'][i]) + 2 * p[i] * (1 - p[i]) * (1 - variables['h'] * variables['s'][i]) + (1 - p[i])**2));
		return {'p': variables['p'], p1: variables['p'].map(p => p[0]), p2: variables['p'].map(p => p[1]), w1: w[0], w2: w[1], 'hs1': hs[0], 'hs2': hs[1], ht, fst};
	}
	const migration = variables['p'].map(pop => [
		(1 - variables['m']) * pop[0] + variables['m'] * pop[1],
		(1 - variables['m']) * pop[1] + variables['m'] * pop[0]
	]);
	const selection = migration.map(repeat => repeat.map((p, i) => {
		const w = p**2 * (1-variables['s'][i]) + 2 * p * (1 - p) * (1 - variables['h'] * variables['s'][i]) + (1 - p)**2;
		return (p**2 * (1-variables['s'][i]) + p * (1 - p) * (1 - variables['h'] * variables['s'][i]))/w;
	}));
	const drift = selection.map(repeat => repeat.map(pop => binomial(2*variables['N'], pop) / (2 * variables['N'])));
	const hs = range(0, 2).map(i => drift.map(pop => 2*pop[i]*(1-pop[i])));
	const ht = drift.map(pop => 2*mean(pop)*(1-mean(pop)));
	const fst = ht.map((ht,i) => ht === 0 ? 0 : 1 - mean([hs[0][i], hs[1][i]]) / ht);
	const w = range(0, 2).map(i => drift.map(p => p[i]**2 * (1-variables['s'][i]) + 2 * p[i] * (1 - p[i]) * (1 - variables['h'] * variables['s'][i]) + (1 - p[i])**2));
	return {'p': drift, p1: drift.map(p => p[0]), p2: drift.map(p => p[1]), w1: w[0], w2: w[1], 'hs1': hs[0], 'hs2': hs[1], 'hs3': hs[2], ht, fst};
};

export const result = (params, steps) => {
	const variables = Object.assign({}, params, steps[steps.length - 1]);
	return {'fst': variables['fst']};
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
