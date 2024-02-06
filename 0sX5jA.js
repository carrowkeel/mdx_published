
const range = (start,end) => Array.from(Array(end-start)).map((v,i)=>i+start);
const sum = arr => arr.reduce((a,v) => a+v, 0);
const mean = arr => arr.length === 0 ? 0 : arr.reduce((a,v) => a+v, 0) / arr.length;
const binomial = (trials, p) => sum(Array.from(new Array(trials)).map(() => Math.random() < p));

export const defaults = () => ({
	'target_steps': 100,
	'p0': [0.1, 0.1, 0.1],
	'm12': 0.01,
	'm13': 0.01,
	'm23': 0.01,
	'N': 10,
	'rep': 100
});

export const step = (params, _step, t) => { // allele frequencies in repeats
	const variables = Object.assign({}, params, _step);
	if (t === 0) {
		variables['p'] = Array.from(new Array(params['rep'])).map(v => variables['p0'].slice());
		const hs = range(0, 3).map(i => variables['p'].map(pop => 2*pop[i]*(1-pop[i])));
		const ht = variables['p'].map(pop => 2*mean(pop)*(1-mean(pop)));
		const fst = ht.map((ht,i) => ht === 0 ? 0 : 1 - mean([hs[0][i], hs[1][i], hs[2][i]]) / ht);
		const fst_pairwise = [[0, 1], [0, 2], [1, 2]].map(i => variables['p'].map((pop, rep) => {
			const pij = mean([pop[i[0]], pop[i[1]]]);
			return pij === 0 || pij === 1 ? 0 : 1 - mean([hs[i[0]][rep], hs[i[1]][rep]])/(2 * pij * (1 - pij))
		}));
		return {'p': variables['p'], 'hs1': hs[0], 'hs2': hs[1], 'hs3': hs[2], ht, fst, fst12: fst_pairwise[0], fst13: fst_pairwise[1], fst23: fst_pairwise[2]};
	}
	const migration = variables['p'].map(pop => [
		(1 - variables['m12'] - variables['m13']) * pop[0] + variables['m12'] * pop[1] + variables['m13'] * pop[2],
		(1 - variables['m12'] - variables['m23']) * pop[1] + variables['m12'] * pop[0] + variables['m23'] * pop[2],
		(1 - variables['m23'] - variables['m13']) * pop[2] + variables['m13'] * pop[0] + variables['m23'] * pop[1]
	]);
	const drift = migration.map(repeat => repeat.map(pop => binomial(2*variables['N'], pop) / (2 * variables['N'])));
	const hs = range(0, 3).map(i => drift.map(pop => 2*pop[i]*(1-pop[i])));
	const ht = drift.map(pop => 2*mean(pop)*(1-mean(pop)));
	const fst = ht.map((ht,i) => ht === 0 ? 0 : 1 - mean([hs[0][i], hs[1][i], hs[2][i]]) / ht);
	const fst_pairwise = [[0, 1], [0, 2], [1, 2]].map(i => drift.map((pop, rep) => {
		const pij = mean([pop[i[0]], pop[i[1]]]);
		return pij === 0 || pij === 1 ? 0 : 1 - mean([hs[i[0]][rep], hs[i[1]][rep]])/(2 * pij * (1 - pij))
	}));
	return {'p': drift, 'hs1': hs[0], 'hs2': hs[1], 'hs3': hs[2], ht, fst, fst12: fst_pairwise[0], fst13: fst_pairwise[1], fst23: fst_pairwise[2]};
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
