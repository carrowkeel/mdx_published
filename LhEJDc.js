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
	return {P_1: params.p1, P_2: params.p2, B_11: params.b11, B_12: params.b12, B_22: params.b22, B_21: params.b21, burden_1: 0, burden_2: 0};
};

const epidemic_step = (B, gamma, _I, _S) => {
	const S = _S.map((S,i) => S - B[i === 0 ? 0 : 2] * S * _I[i] - B[i === 0 ? 3 : 1] * S * _I[i ? 0 : 1]);
	const I = _I.map((I,i) => I + B[i === 0 ? 0 : 2] * _S[i] * I + B[i === 0 ? 3 : 1] * _S[i] * _I[i ? 0 : 1] - gamma * I);
	return [I, S];
};

const epidemic = (B, gamma) => {
	let I = [0, 0.01];
	let S = I.map(I => 1 - I);
	let R = 0;
	let step = 0;
	while (++step) {
		[I, S] = epidemic_step(B, gamma, I, S);
		const dR = gamma * I[0];
		R += dR;
		if (step > 1000)
			break;
	}
	return R;
};

const defaults = () => ({
	a: 0.00005,
	b: 0.02,
	b11: 0.05,
	b12: 0.01,
	b21: 0.01,
	b22: 0.05,
	c: 100,
	d: 1,
	gamma: 0.045,
	p1: 100,
	p2: 100,
	target_steps: 200
});

const step = (params, _step, t) => {
	if (t === 0)
		return initial(params);
	const initial_contact_rates = [params.b11, params.b12, params.b22, params.b21];
	const _P = [_step.P_1, _step.P_2];
	const _B = [_step.B_11, _step.B_12, _step.B_22, _step.B_21];
	const D = _P.map((P,i) => P * epidemic(i === 0 ? _B : [_B[2], _B[3], _B[0], _B[1]], params.gamma));
	const B = _B.map((B,i) => B - params.a * D[i === 0 || i === 3 ? 0 : 1] + params.b * (initial_contact_rates[i] - B));
	const P = _P.map((P,i) => Math.max(P - params.c * _B[i === 0 ? 3 : 1], 0));
	return {P_1: P[0], P_2: P[1], B_11: B[0], B_12: B[1], B_22: B[2], B_21: B[3], burden_1: D[0], burden_2: D[1]};
};

const result = (params, stats) => {
	return {};
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