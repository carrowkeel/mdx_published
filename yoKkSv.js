const round = (n,p) => { var f = Math.pow(10, p); return Math.round(n * f) / f };
const sum = (arr) => arr.reduce((a,v)=>a+v,0);
const mean = (arr) => arr.length === 0 ? 0 : arr.reduce((a,v)=>a+v,0)/arr.length;
const range = (start,end,step=1) => Array.from(Array(Math.floor((end-start)/step))).map((v,i)=>start+i*step);
const cumsum = arr => { let sum = 0; const out = []; for (const v of arr) { sum += v; out.push(sum) } return out; }
const common_value = values => Object.entries(values.reduce((a,value) => Object.assign(a, {[value]: (a[value] || 0) + 1}), {})).sort((a,b)=>b[1]-a[1])[0][0];
const poisson = (l, L=Math.exp(-l), n=0, p=1) => {
	while(true) {
		if (p <= L)
			return n - 1;
		n += 1;
		p *= Math.random();
	}
};

const model_outcomes = ['Collapse', 'Spillover', 'Failure (resistance)', 'Rescue (resistance)', 'Failure', 'Gene swamping', 'Oscillations', 'Differential targeting', 'Suppression'];

const initial = params => {
	return {q1: params.q1, q2: params.q2, p1: 0, p2: 0, N1: params.K * (1 - params.eps), N2: params.K * (1 - params.eps)};
};

const r_s = (params, r0, s) => {
	if (params.r_s === 'Fixed')
		return 1 / r0;
	else if (params.N_model === 'Logistic difference equation') {
		if (params.r_s === 'Eq. S1')
			return -s * params.d;
		else
			return r0 - (1 - Math.E**(-s * params.d)) * (r0 + s);
	} else if (params.N_model === 'Beverton-Holt') {
		if (params.r_s === 'Eq. S1')
			return 1 - s * params.d;
		else
			return r0 - (1 - Math.E**(-s * params.d)) * (r0 + s - 1);
	}
};

const step = (params, _step, t) => {
	if (t === 0)
		return initial(params);
	const s = params.s;
	const c = params.c;
	const h = params.h;
	const m = params.m;
	const K = params.K;
	const s_b = params.s_b;
	const post_migration_freq = [{
		q: ((1 - m) * _step.q1 + m * (_step.N2/_step.N1) * _step.q2) * 1/(1 - m + m * (_step.N2/_step.N1)),
		p: ((1 - m) * _step.p1 + m * (_step.N2/_step.N1) * _step.p2) * 1/(1 - m + m * (_step.N2/_step.N1))
	},
	{
		q: ((1 - m) * _step.q2 + m * (_step.N1/_step.N2) * _step.q1) * 1/(1 - m + m * (_step.N1/_step.N2)),
		p: ((1 - m) * _step.p2 + m * (_step.N1/_step.N2) * _step.p1) * 1/(1 - m + m * (_step.N1/_step.N2))
	}];

	const N_m = [
		_step.N1 * (1 - m) + _step.N2 * m,
		_step.N1 * m + _step.N2 * (1 - m)
	];

	if (params['mu'] > 0 || params['nhej'] > 0) {
		for (const pop_i of [0, 1]) {
			const rate = params['mu'] * 2 * N_m[pop_i] * ((1 - post_migration_freq[pop_i].q)**2 + post_migration_freq[pop_i].q * (1 - post_migration_freq[pop_i].q) * (1 - c)) +
				     params['nhej'] * 2 * N_m[pop_i] * c * post_migration_freq[pop_i].q * (1 - post_migration_freq[pop_i].q);
			const resistance_alleles = rate === 0 ? 0 : poisson(rate);
			if (resistance_alleles > 0)
				post_migration_freq[pop_i].p += resistance_alleles / (2 * N_m[pop_i])
		}
	}

	const w_bar = post_migration_freq.map(({q, p}) => params.conversion === 'Zygotic' ?
		2 * q * (1 - q - p) * (1 - c) * (1 - h * s) + // Zygotic (Eq. S4)
		2 * p * (1 - q - p) * (1 - s_b / 2) +
		2 * p * q * (1 - (s * h + s_b * (1 - h))) +
		(q ** 2 + 2 * q * (1 - q - p) * c) * (1 - s) +
		p ** 2 * (1 - s_b) +
		(1 - q - p) ** 2 :

		2 * q * (1 - q - p) * (1 - h * s) + // Gametic (Eq. 3)
		2 * p * (1 - q - p) * (1 - s_b / 2) +
		2 * p * q * (1 - (s * h + s_b * (1 - h))) +
		q ** 2 * (1 - s) +
		p ** 2 * (1 - s_b) +
		(1 - q - p) ** 2
	);

	const post_selection_freq = post_migration_freq.map(({q, p}, pop_i) => ({
		q: params.conversion === 'Zygotic' ?
			(q * (1 - q - p) * (1 - c) * (1 - h * s) + // Zygotic (Eq. S4)
			p * q * (1 - (s * h + s_b * (1 - h))) +
			(q ** 2 + 2 * q * (1 - q - p) * c) * (1 - s)) / w_bar[pop_i] :
	
			(q * (1 - q - p) * (1 + c) * (1 - h * s) + // Gametic (Eq. 3)
			p * q * (1 - (s * h + s_b * (1 - h))) +
			q ** 2 * (1 - s)) / w_bar[pop_i],

		p: (p * (1 - q - p) * (1 - s_b / 2) +
		p * q * (1 - (s * h + s_b * (1 - h))) +
		p ** 2 * (1 - s_b)) / w_bar[pop_i]
	}));

	const r0 = params.N_model === 'Beverton-Holt' ? params.r0 : params.r0 / 2;

	const R = post_migration_freq.map(({q, p}, pop_i) => params.conversion === 'Zygotic' ?
		(2 * q * (1 - q - p) * (1 - c) * (1 - h * s) * r_s(params, r0, h * s) + // Zygotic
		2 * p * (1 - q - p) * (1 - s_b / 2) * r_s(params, r0, s_b / 2) +
		2 * p * q * (1 - (s * h + s_b * (1 - h))) * r_s(params, r0, (s * h + s_b * (1 - h))) +
		(q ** 2 + 2 * q * (1 - q - p) * c) * (1 - s) * r_s(params, r0, s) +
		p ** 2 * (1 - s_b) * r_s(params, r0, s_b) +
		(1 - q - p) ** 2 * r0) / w_bar[pop_i] :

		(2 * q * (1 - q - p) * (1 - h * s) * r_s(params, r0, h * s) + // Gametic
		2 * p * (1 - q - p) * (1 - s_b / 2) * r_s(params, r0, s_b / 2) +
		2 * p * q * (1 - (s * h + s_b * (1 - h))) * r_s(params, r0, (s * h + s_b * (1 - h))) +
		q ** 2 * (1 - s) * r_s(params, r0, s) +
		p ** 2 * (1 - s_b) * r_s(params, r0, s_b) +
		(1 - q - p) ** 2 * r0) / w_bar[pop_i]
	);

	const N_n = N_m.map((N, pop_i) => {
		if (params.N_model == 'Logistic difference equation')
			return N + R[pop_i] * N * (1 - (N / K)) ;
		else
			return (R[pop_i] * K * N) / (K + (R[pop_i] - 1) * N);
	});

	return {'q1': post_selection_freq[0].q, 'q2': post_selection_freq[1].q, 'p1': post_selection_freq[0].p, 'p2': post_selection_freq[1].p, 'N1': N_n[0], 'N2': N_n[1]}
};

const outcome = (params, steps) => {
	const freq = [steps[steps.length-1].q1, steps[steps.length-1].q2];
	const resistance_freq = [steps[steps.length-1].p1, steps[steps.length-1].p2];
	const population = steps.reduce((a,step) => range(0, 2).map(i => a[i].concat(step[`N${i+1}`] / params.K)), [[], []]);
	const max_values = population.map(p => p.filter((v,i) => i > 0 && i < p.length-1 && ((p[i-1] < v && v > p[i+1]) || (p[i-1] > v && v < p[i+1]))));
	const oscillating = max_values[0].length <= 1 ? false : Math.abs(max_values[0][max_values[0].length - 1] - max_values[0][max_values[0].length - 2]) > 0.5;
	const min = population.map(pop => Math.min.apply(null, pop));
	const pop_end = population.map(v => v[v.length-1]);
	const gen_threshold = 0.5;
	const pop_threshold = 0.9;
	const res_threshold = 0.1;
	switch(true) {
		case min[0] < pop_threshold && min[1] < pop_threshold: 
			return 'Collapse';
		case freq[0] >= gen_threshold && freq[1] >= gen_threshold:
			return 'Spillover';
		case freq[0] < gen_threshold && freq[1] < gen_threshold && (resistance_freq[0] >= res_threshold || resistance_freq[1] >= res_threshold) && min[0] >= pop_threshold && min[1] >= pop_threshold:
			return 'Failure (resistance)';
		case freq[0] < gen_threshold && freq[1] < gen_threshold && (resistance_freq[0] >= res_threshold || resistance_freq[1] >= res_threshold) && min[0] < pop_threshold && min[1] >= pop_threshold:
			return 'Rescue (resistance)';
		case freq[0] < gen_threshold && freq[1] < gen_threshold && min[0] < pop_threshold && pop_end[0] >= pop_threshold:
			return 'Gene swamping';
		case freq[0] < gen_threshold && freq[1] < gen_threshold && min[0] >= pop_threshold && min[1] >= pop_threshold:
			return 'Failure';
		case freq[0] > 0 && freq[1] < gen_threshold && oscillating:
			return 'Oscillations';
		case freq[1] < gen_threshold && min[0] < pop_threshold:
			return 'Suppression';
		case freq[0] >= gen_threshold && freq[1] < gen_threshold && min[0] >= pop_threshold && min[1] >= pop_threshold:
			return 'Differential targeting';
		default:
			return 'undefined';
	}
};

const sim = (params) => {
	const storage = [];
	let s = undefined;
	while (s !== false && storage.length <= params.target_steps) {
		s = step(params, s, storage.length);
		if (s)
			storage.push(s);
	}
	return result(params, storage);
};

const mcrit = (params) => {
	let m = 0;
	let incr = -2;
	const precision = -5;
	while (true) {
		const result = sim(Object.assign({}, params, {m, stat: 'outcome'}));
		if (!['failure', 'spillover', 'collapse'].includes(result.outcome)) {
			m += Math.pow(10, incr);
		} else if (m !== 0) {
			m -= Math.pow(10, incr);
			if (incr === precision)
				return {mcrit: Math.min(1, m)};
			incr--;
			m += Math.pow(10, incr);
		} else {
			return {mcrit: Math.min(1, m)};
		}
	}
};

const defaults = () => ({
	s: 0.5,
	c: 1,
	h: 1,
	d: 10,
	r0: 2,
	r_s: 'Eq. 6',
	m: 0.01,
	q1: 0.8,
	q2: 0,
	mu: 0,
	nhej: 0,
	s_b: 0,
	K: 1,
	conversion: 'Gametic',
	N_model: 'Beverton-Holt',
	eps: 0.01,
	target_steps: 100,
	repeats: 1,
	outcome: 'Suppression',
	output: 'Outcome'
});

const outcomeProp = (params) => {
	const s = range(0, 1.01, 0.01);
	const c = range(0, 1.01, 0.01);
	const outcomes = Object.fromEntries(model_outcomes.map(outcome => [outcome, 0]));
	for (const s_i of s) {
		for (const c_i of c) {
			outcomes[sim(Object.assign({}, params, {s: s_i, c: c_i}))['outcome']] += 1 / (s.length * c.length);
		}
	}
	return outcomes;
};

const run = (params) => {
	switch(params.output) {
		case 'm*':
			return mcrit(params);
		case 'Outcome proportion':
			return outcomeProp(params);
		default:
			if (params.repeats === 1)
				return sim(params);
			const runs = Array.from(new Array(params.repeats)).map(_ => sim(params)['outcome']);
			return {outcome_proportion: runs.filter(run => run === params.outcome).length / params.repeats};
	}
};

const result = (params, steps) => {
	const model_outcome = outcome(params, steps);
	const non_target = steps.sort((a,b) => b.q_2 - a.q_2)[0].q_2;
	const population = steps.reduce((a,step) => range(0, 2).map(i => a[i].concat(step[`N${i+1}`] / params.K)), [[], []]);
	return {
		outcome: model_outcome,
		max_suppression: model_outcome === 'Gene swamping' ? 1 - Math.min.apply(null, population[0]) : 0,
		sum_suppression: model_outcome === 'Gene swamping' ? sum(population[0].map(v => 1 - v)) / population[0].length : 0,
		dur_suppression: model_outcome === 'Gene swamping' ? Math.min(1, population[0].filter(v => v<0.9).length / population[0].length) : 0,
		non_target_freq: non_target
	};
};

export { step, run, defaults }