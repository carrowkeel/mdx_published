const randint = (m,m1,g=Math.random) => Math.floor(g() * (m1 - m)) + m;
const choose = arr => arr[randint(0, arr.length)];
const round = (n,p) => { var f = Math.pow(10, p); return Math.round(n * f) / f };
const sum = (arr) => arr.reduce((a,v)=>a+v,0);
const mean = (arr) => arr.length === 0 ? 0 : arr.reduce((a,v)=>a+v,0)/arr.length;
const range = (start,end) => Array.from(Array(end-start)).map((v,i)=>i+start);
const cumsum = arr => { let sum = 0; const out = []; for (const v of arr) { sum += v; out.push(sum) } return out; }
const vshuffle = (arr, g) => {
	var ridx = Array(arr.length).fill(0).map((v,i)=>[i,g()]).sort((a,b)=>a[1]-b[1]); // Annoying
	return ridx.map(v=>arr[v[0]]);
};

const stats = (params, step) => {
	return Object.assign(step, {
		patches_population: step.patches.reduce((a, patch) => [a[0] + patch[0], a[1] + patch[1]], [0, 0]).map(v => v / params.patches),
		particles_population: step.particles.reduce((a, patch) => [a[0] + patch[0], a[1] + patch[1]], [0, 0]).map(v => v / (params.gridsize ** 2)),
		particles_spatial: step.particles.map((patch, i) => [i % params.gridsize, Math.floor(i / params.gridsize), 1, 1, patch[0] + patch[1] === 0 ? [0, 0, 0] : (patch[0] > patch[1] ? [31, 119, 180] : (patch[0] < patch[1] ? [255, 127, 14] : [255, 255, 255])), 1])
	});
};

const initial = params => {
	return {
		dynamical: [params.h0, params.d0],
		patches: new Array(params.patches).fill([Math.ceil(params.N * params.h0), Math.ceil(params.N * params.d0)]),
		particles: new Array(params.gridsize**2).fill([Math.ceil(params.N * params.h0), Math.ceil(params.N * params.d0)])
	};
};

const dynamical_system = (params, [_h, _d]) => {
	return [
		Math.max(0, _h + (_h * (params.a * (_h / (_h + _d)) + params.b * (_d / (_h + _d)) - params.K * (_h + _d))) / params.timestep),
		Math.max(0, _d + (_d * (params.c * (_h / (_h + _d)) + params.d * (_d / (_h + _d)) - params.K * (_h + _d))) / params.timestep)
	];
};

const patches_random = (params, _patches) => {
	range(0, params.patches).forEach(patch => {
		const [h, d] = _patches[patch];
		for (const i in range(0, h)) {
			if (Math.random() < params.mu / params.timestep) {
				_patches[patch][0] -= 1;
				_patches[randint(0, params.patches)][0] += 1;
			}
		}
		for (const i in range(0, d)) {
			if (Math.random() < params.mu / params.timestep) {
				_patches[patch][1] -= 1;
				_patches[randint(0, params.patches)][1] += 1;
			}
		}
	});
	const crowding = _patches.map(([_h, _d]) => {
		const p = _h / (_h + _d);
		let [h, d] = [_h, _d];
		for (const i in range(0, _h))
			if (Math.random() < params.K * (_h + _d) / params.timestep)
				h -= 1;
		for (const i in range(0, _d))
			if (Math.random() < params.K * (_h + _d) / params.timestep)
				d -= 1;
		return [h, d];
	});
	return crowding.map(([_h, _d]) => {
		const p = _h / (_h + _d);
		let [h, d] = [_h, _d];
		const birth = [params.a * p + params.b * (1 - p), params.c * p + params.d * (1 - p)];
		for (const i in range(0, _h))
			if (Math.random() < Math.abs(birth[0]) / params.timestep)
				h += birth[0] >= 0 ? 1 : -1;
		for (const i in range(0, _d))
			if (Math.random() < Math.abs(birth[1]) / params.timestep)
				d += birth[1] >= 0 ? 1 : -1;
		return [h, d];
	});
};

const neighborhood = (grid, i, neighborhood) => {
	switch(neighborhood) {
		case 'nearest':
			return [i, i % grid === 0 ? i + grid - 1 : i - 1, (i + 1) % grid === 0 ? i - grid + 1 : i + 1, i - grid < 0 ? grid**2 + i - grid : i - grid, i + grid >= grid * grid ? i % grid : i + grid];
	}
};

const particle_system = (params, _patches) => {
	_patches.forEach(([h, d], patch) => {
		const n = neighborhood(params.gridsize, patch, 'nearest').slice(1);
		for (const i in range(0, h)) {
			if (Math.random() < params.mu / params.timestep) {
				_patches[patch][0] -= 1;
				_patches[choose(n)][0] += 1;
			}
		}
		for (const i in range(0, d)) {
			if (Math.random() < params.mu / params.timestep) {
				_patches[patch][1] -= 1;
				_patches[choose(n)][1] += 1;
			}
		}
	});
	const crowding = _patches.map(([_h, _d], patch) => {
		let [h, d] = [_h, _d];
		for (const i in range(0, _h))
			if (Math.random() < params.K * (_h + _d) / params.timestep)
				h -= 1;
		for (const i in range(0, _d))
			if (Math.random() < params.K * (_h + _d) / params.timestep)
				d -= 1;
		return [h, d];
	});
	return crowding.map(([_h, _d], patch) => {
		const n = neighborhood(params.gridsize, patch, params.neighborhood);
		const [nh, nd] = n.reduce(([h, d], i) => [h + crowding[i][0], d + crowding[i][1]], [0, 0]);
		const p = nh / (nh + nd);
		let [h, d] = [_h, _d];
		const birth = [params.a * p + params.b * (1 - p), params.c * p + params.d * (1 - p)];
		for (const i in range(0, _h))
			if (Math.random() < Math.abs(birth[0]) / params.timestep)
				h += birth[0] >= 0 ? 1 : -1;
		for (const i in range(0, _d))
			if (Math.random() < Math.abs(birth[1]) / params.timestep)
				d += birth[1] >= 0 ? 1 : -1;
		return [h, d];
	});
}

export const defaults = () => {
	return {
		a: 0.4,
		b: 0.8,
		c: 0.6,
		d: 0.3,
		K: 0.1,
		h0: 2,
		d0: 2,
		patches: 1000,
		gridsize: 100,
		neighborhood: 'nearest',
		mu: 0.01,
		N: 1,
		timestep: 1,
		target_steps: 100
	};
};

export const step = (params, _step, t) => {
	if (t === 0)
		return stats(params, initial(params));
	let step = _step;
	for (let i=0;i<params.timestep;i++)
		step = {dynamical: dynamical_system(params, step.dynamical), patches: patches_random(params, step.patches), particles: particle_system(params, step.particles)};
	return stats(params, step);
};

const outcome = steps => {
	const zero = 0.05;
	const collapsed = ['dynamical', 'patches_population', 'particles_population'].map(stat => (steps[steps.length-1][stat][0] > zero ? 2 : 0) | (steps[steps.length-1][stat][1] > zero));
	switch(true) {
		case collapsed.reduce((a,v) => a === v ? v : false) !== false:
			return 'agreement';
		case collapsed[0] === collapsed[1]:
			return 'spatial';
		case collapsed[1] === collapsed[2]:
			return 'discrete';
		default:
			return 'spatial_discrete';
	}
};

export const result = (params, steps) => {
	return {
		outcome: outcome(steps)
	};
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