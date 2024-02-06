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
const vcomp = (a,b,m) => a.map((v,i) => Array.isArray(v)?vcomp(v,b[i],m):b[i]!=null?(m==2?v*b[i]:m==1?v-b[i]:v+b[i]):v);
const stddev = arr => Math.pow(vcomp(Array(arr.length).fill(-mean(arr)),arr).map((v)=>Math.pow(v,2)).reduce((a,b)=>a+b) / (arr.length - 1), 0.5);
const rand_generator = (seed) => () => { let t = seed += 0x6D2B79F5; t = Math.imul(t ^ t >>> 15, t | 1); t ^= t + Math.imul(t ^ t >>> 7, t | 61); return ((t ^ t >>> 14) >>> 0) / 4294967296; };

const color_cycle = [
	[0, 255, 255],
	[255, 0, 255],
	[31,119,180],
	[255,127,14],
	[44,160,44],
	[214,39,40],
	[148,103,189],
	[140,86,75],
	[227,119,194],
	[127,127,127],
	[188,189,34],
	[23,190,207]
];

const vcol = (arr, col=0, f=(v)=>v) => { let o = [];for(let i of arr){let v=f(i[col]);if(v!==false){o.push(v)}} return o};
const vfilter = (arr, f) => { let o = []; for (let i of arr) { let v=f(i); if(v){o.push(i)} } return o };

const randchoice = (total, select, rand) => vshuffle(range(0, total), rand).slice(0, select);
const randof = (choices, rand) => { let o = []; let i = choices.length; for (let x in choices[0]) { o.push(choices[randint(0, i, rand)][x]); } return o; };

const initial = params => {
	const signalsn = 2;
	Object.assign(params, {
		signals: [params.signal_1, params.signal_2],
		signalsn,
		_qi: 1,
		_ri: 2,
		_ai: 3,
		_sl: 3 + signalsn,
		_mi: 3 + signalsn * 2,
		maxchoice: 15,
		muts: 20,
		recmuts: 5,
		qualityspread: 0.25,
		mutationstep: 0.25,
		offspring: 3,
		popsize: 500,
		initialadvertising: 0.5,
		g: rand_generator(Math.round(Math.random()*1e9))
	});
	return {agents: Array.from(new Array(params.popsize)).map(v => [
		randint(0, 2, params.g),
		(params.g() / 2) + 0.25,
		randint(0, params.signalsn, params.g),
		...Array(params.signalsn).fill(0).map((v,i)=>params.initialadvertising),
		...Array(params.signalsn).fill(0),
		0,
	])};
};

const siglevel = (params, agent, cue, signalcost) => {
	const sig = agent[params._qi] * agent[params._ai+cue];
	const sigcost = agent[params._qi] * (agent[params._ai+cue] - params.initialadvertising) * signalcost;
	const out = (v) => (v < params.signallimit ? v : params.signallimit);
	return [out(sig), sigcost];
};

const step = (params, _step, t) => {
	if (t === 0)
		return stats(params, initial(params));

	const _agents = _step.agents;

	if (_agents.length === 0)
		return false;

	let mn = vcol(_agents, 0, (v)=>v===0?true:false).length;
	let fn = vcol(_agents, 0, (v)=>v===1?true:false).length;

	if (mn === 0 || fn === 0)
		return false;

	const survived = [];
	for (const a of _agents) {
		let q_a = a[params._qi];
		if (a[0] === 0) {
			for (let s of range(0, params.signalsn)) {
				let [sig, cost] = siglevel(params, a, s, params.signals[s]);
				a[params._sl+s] = sig;
				q_a -= cost;
			}
		}
		if (params.g() > 1 - q_a)
			survived.push(a);
	}

	const agents = survived.slice(0);

	if (agents.length === 0)
		return false;

	const males = vfilter(agents, (v)=>v[0]===0);
	const females = vfilter(agents, (v)=>v[0]===1);
	if (males.length === 0 || females.length === 0)
		return false;

	const next = [];

	for (const a of females) {
		if (males.length === 0)
			break;
		const h = males.length < params.maxchoice ? males.length : params.maxchoice;
		const ind = randchoice(males.length, h, params.g);
		const signals = [];
		const index = {};
		for (const p of ind) {
			signals.push(males[p][params._sl+a[params._ri]]);
			index[males[p][params._sl+a[params._ri]]] = p;
		}

		const threshold = Math.max(...signals) - params.prefresolution;
		const potential = signals.filter(v=>v>=threshold);
		const chosen = index[potential[randint(0, potential.length, params.g)]];

		males[chosen][params._mi] += 1;

		for (let _x of range(0, params.offspring)) {
			const f1 = randof([a, males[chosen]], params.g);
			f1[params._mi] = 0;
			f1[params._qi] = (a[params._qi] + males[chosen][params._qi]) / 2;
			f1[params._qi] += params.qualityspread * (params.g()-0.5);
			f1[params._qi] = f1[params._qi] > 1 ? 1 : f1[params._qi];
			next.push(f1);
		}
		if (males[chosen][params._mi] >= params.maxmated)
			males.splice(chosen,1);
	}

	const shuffled = vshuffle(next, params.g).slice(0, params.popsize);

	for (const a of randchoice(shuffled.length, params.recmuts, params.g)) {
		shuffled[a][params._ri] = randint(0, params.signalsn, params.g);
	}
	for (const a of randchoice(shuffled.length, params.muts, params.g)) {
		const mutsig = randint(params._ai, params._ai+params.signalsn, params.g);
		shuffled[a][mutsig] += params.mutationstep * (params.g()-0.5);
		shuffled[a][mutsig] = shuffled[a][mutsig] < params.initialadvertising ? params.initialadvertising : shuffled[a][mutsig];
	}

	return stats(params, {agents: shuffled});
};

const defaults = () => ({
	maxmated: 10,
	prefresolution: 0.1,
	repeats: 1,
	signal_1: 0,
	signal_2: 2,
	signallimit: 1,
	target_steps: 1000
});

const stats = (params, step) => {
	const rechist = Array(params.signalsn).fill(0);
	for (const agent of step.agents)
		rechist[agent[params._ri]] += 1;
	const ahist = Array(params.signalsn).fill(0);
	for (const agent of step.agents)
		for (const s of range(0, params.signalsn))
			ahist[s] += agent[params._ai+s];
	return Object.assign(step, {
		quality_signal: range(0, params.signalsn).map(s => step.agents.map(agent => [agent[params._qi], Math.min(1, agent[params._ai+s] * agent[params._qi])])),
		a: range(0, params.signalsn).map(s => ahist[s]/step.agents.length),
		pref: range(0, params.signalsn).map(s => rechist[s]/step.agents.length)
	});
};

const result = (params, stats) => {
	const pref = stats.map(stat => stat.pref);
	const means = range(0, params.signalsn).map(i => mean(pref.map(v => v[i])));
	const sorted = means.slice().sort((a,b) => b-a);
	const signal = (sorted[0] - sorted[1] < 0.1 ? 0 : means.indexOf(sorted[0]) + 1);
	return {signal};
};

const run = (params) => {
	const stats = [];
	let _step = undefined;
	while (_step !== false && stats.length <= params.target_steps) {
		_step = step(params, _step, stats.length);
		if (_step)
			stats.push(_step);
	}
	return result(params, stats);
};

export { defaults, step, result, run }