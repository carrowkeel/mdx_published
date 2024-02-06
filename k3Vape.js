'use strict';

const randint = (m,m1,g) => Math.floor(g() * (m1 - m)) + m;
const round = (n,p) => { var f = Math.pow(10, p); return Math.round(n * f) / f };
const sum = (arr) => arr.reduce((a,v)=>a+v,0);
const mean = (arr) => arr.length === 0 ? 0 : arr.reduce((a,v)=>a+v,0)/arr.length;
const range = (start,end) => Array.from(Array(end-start)).map((v,i)=>i+start);
const cumsum = arr => { let sum = 0; const out = []; for (const v of arr) { sum += v; out.push(sum) } return out; }
const vshuffle = (arr, g) => {
	var ridx = Array(arr.length).fill(0).map((v,i)=>[i,g()]).sort((a,b)=>a[1]-b[1]); // Annoying
	return ridx.map(v=>arr[v[0]]);
};
const rand_generator = (seed) => () => { let t = seed += 0x6D2B79F5; t = Math.imul(t ^ t >>> 15, t | 1); t ^= t + Math.imul(t ^ t >>> 7, t | 61); return ((t ^ t >>> 14) >>> 0) / 4294967296; };

const initial = params => {
	return {N_1: 2, N_2: 2, I_1: 0, I_2: 0, a_1: params.a1, a_2: params.a2};
};

const defaults = () => ({
	K: 1,
	l: 0.15,
	b11: 0.5,
	b22: 0.5,
	ba: 0.1,
	mu: 0.05,
	c: 0.1,
	a1: 1,
	a2: 1,
	d: 1,
	target_steps: 200
});

const step = (params, _step, t) => {
	if (t === 0)
		return stats(params, initial(params));
	const _N = [_step.N_1, _step.N_2];
	const _I = [_step.I_1, _step.I_2];
	const _a = [_step.a_1, _step.a_2];
	const N = _N.map((N,i) => N + (N*(params.l * params.K) / (params.K + N) - params.mu * N - _a[i] * _I[i]) / params.d);
	const I = _I.map((I,i) => I + ((_N[i] - I) * (params.ba*_N[i ? 0 : 1] + (i === 0 ? params.b11 : params.b22) * I) - (params.mu + _a[i]) * I) / params.d);
	const a = _a.map((a,i) => Math.max(0, a <= 0 ? 0 : a - (params.c * params.ba * _N[i ? 0 : 1]) / params.d));
	return stats(params, {N_1: N[0], N_2: N[1], I_1: I[0], I_2: I[1], a_1: a[0], a_2: a[1]});
};

const stats = (params, step) => {
	return Object.assign(step, {
		burden_1: step.I_1*step.a_1,
		burden_2: step.I_2*step.a_2
	});
};

const result = (params, steps) => {
	return {population_end: steps[steps.length-1].N_1};
};

const run = (params) => {
	let prev = undefined;
	for(let t = 0;t <= params.target_steps;t++) {
		const r = step(params, prev, t);
		if (!r)
			break;
		prev = r;
	}
	return result(params, [prev]);
};

export { defaults, step, result, run }