
const sum = arr => arr.reduce((a,v) => a+v, 0);
const mean = arr => arr.length === 0 ? 0 : arr.reduce((a,v) => a+v, 0) / arr.length;
const cumsum = arr => { let sum = 0; const out = []; for (const v of arr) { sum += v; out.push(sum) } return out; }
const binomial = (trials, p) => sum(Array.from(new Array(trials)).map(() => Math.random() < p));
const multi = (n, p) => {
	const s = cumsum(p);
	const res = [];
	while(res.length < n) {
		const r = Math.random();
		let i = 0;
		while(i < s.length) {
			if (r < s[i]) {
				res.push(i);
				break;
			}
			i++;
		}
	}
	return res;
};

export const defaults = () => ({
	'target_steps': 100,
	'p0': [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1],
	'rep': 100,
	'N0': 10,
	'mu': 0.01,
	'lambda': 1.05
});

export const step = (params, _step, t) => { // allele frequencies in repeats
	const variables = Object.assign({}, params, _step);
	const N = Math.round(Math.min(10000, variables['N0'] * variables['lambda']**t));
	if (N === 0)
		return false;
	if (t === 0) {
		variables['p'] = Array.from(new Array(params['rep'])).map(v => variables['p0'].slice());
		const heterozygosity = variables['p'].map(repeat => (1 - sum(repeat.map(freq => freq**2))));
		return {'p': variables['p'], 'N': N, 'alleles': variables['p'].map(repeat => repeat.filter(allele => allele !== 0).length), 'heterozygosity': heterozygosity, 'effective_alleles': heterozygosity.map(het => 1/(1-het))};
	}
	const drift = variables['p'].map(repeat => multi(2*N, repeat).reduce((a,v) => { a[v] += 1; return a }, new Array(repeat.length).fill(0)).map(v => v / (2*N)));
	const mutation = drift.map(repeat => {
		const mutations = binomial(2*N, variables['mu']);
		for(let i=0;i<mutations;i++) {
			const source = multi(1, repeat);
			repeat[source] -= 1/(2*N);
			repeat.push(1/(2*N));
		}
		return repeat;
	});
	const heterozygosity = mutation.map(repeat => (1 - sum(repeat.map(freq => freq**2))));
	return {'p': mutation, 'N': N, 'alleles': mutation.map(repeat => repeat.filter(allele => allele !== 0).length), 'heterozygosity': heterozygosity, 'effective_alleles': heterozygosity.map(het => 1/(1-het))};
};

export const result = (params, steps) => {
	const variables = Object.assign({}, params, steps[steps.length - 1]);
	return {'alleles': mean(variables['alleles'])};
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
