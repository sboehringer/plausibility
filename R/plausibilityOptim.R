#
#	plausibilityOptim.R
#Thu Jun 17 17:07:02 CEST 2021

#
#	<p> approximtion
#

evaluateCubesFor = function(g, v, Ngrid, cache, ..., f, conditions) {
	Ndim = ncol(g);
	# iterate sub-cubes
	Ncubes = (Ngrid - 1)^Ndim;
	# spacing between indeces, tb cached
	cubeI = sapply(0:(Ndim - 1), function(i)Ngrid^i)
	# where are the indeces of the other cube corners
	cubeIdxSpacing = t_(cubeI * sapply(0:(2^Ndim - 1), ord2bin, digits = Ndim)) %*% rep(1, Ndim);

	rCubes = lapply(1:Ncubes, function(i) {
		# lower left corner
		cubeLL = (ord2base(i - 1, digits = Ndim, base = Ngrid - 1) %*% cubeI)[1, 1] + 1;
		# nodes belonging to current cube
		cubeIdcs = cubeLL + cubeIdxSpacing;
		# extract coordinates
		coords = g[cubeIdcs, , drop = F];
		approximateFunctionFor(f,
			mn = apply(coords, 2, min), mx = apply(coords, 2, max), Ngrid = Ngrid,
			cache = list(args = coords, v = v[cubeIdcs]), ..., conditions = conditions
		)
	});
	return(unique(do.call(rbind, rCubes)));
}
approximateFunctionFor = function(f, mn, mx, ..., Ngrid = 4, cache = NULL, conditions) {
	Ndim = length(mn);
	# <p> stopping condition
	if (conditions$mnMx(mn, mx)) return(NULL);

	g = gridTensor(mn, mx, N = Ngrid);
	v = evaluateFunctionCached(g, f, cache, ...);
	# <p> stopping condition, do not stop on first run
	if (notE(cache) && conditions$value(v, g, cache, mn, mx)) return(cbind(g, v));
	return(evaluateCubesFor(g, v, Ngrid = Ngrid, cache = list(args = g, v = v), ...,
		f = f, conditions = conditions));
}

approximated = function(r, Nnn = 1, aggr = mean) {
	X = t(r[, -ncol(r), drop = F]);
	Y = r[, ncol(r)];
	return(function(xs) {
		if (is.vector(xs)) xs = t_(xs);
		apply(xs, 1, function(x) {
			D = apply((X - x)^2, 2, sum);
			o = order(D);
			return(aggr(Y[o[1:Nnn]]));
	})})
}

conditionsContour = function(contour = .05, epsX = 1e-2, epsY = contour/2 * 1.25, tolerance = 2.5) {
	return(list(
		mnMx = function(mn, mx)(any(abs(mx - mn) < epsX)),
		value = function(v, grid, cache, ...) (
			max(v) < contour - epsY || min(v) > contour + epsY || max(v) - min(v) < epsY / tolerance
		)
	));
}

conditionsMax = function(sign = 1, epsX = 1e-3, epsY = 1e-3, tolerance = 1/epsY * .8) {
	return(list(
		mnMx = function(mn, mx)(any(abs(mx - mn) < epsX)),
		value = function(v, grid, cache, mn, mx) {
			if (is.null(cache)) return(FALSE);
			Icache = DfSearch(Df_(cache$args), Df_(grid));	# old values
			v = v[-Icache];
			grid = grid[-Icache, , drop = F];	# can be commented out <d>
			delta = min(mx - mn);
			#print(list(v, grid, cache))
			#dp = list(delta, max(sign*v), max(sign * cache$v) - tolerance * delta * epsY, abs(max(v) - max(cache$v)))
			#print(list2df(dp))
			return(	max(sign * v) <= max(sign * cache$v) - tolerance * delta * epsY ||
					abs(max(sign * v) - max(sign * cache$v)) <= epsY);
		}
	));
}

#
#	<p> optimization
#

plausibilityOptimizeOptim = function(this, objectiveFunction, start,
	method = 'Nelder-Mead', control = list(fnscale = -1)) {

	if (is.null(start)) start = plausibilityStart(this);
	#print(list(start = start));
	#if (start == 0) browser();
	lower = start - plausibilityDelta(this);
	upper = start + plausibilityDelta(this);
	o = optimFn(start, objectiveFunction, method = method, control = control,
		lower = lower, upper = upper, this = this);
	P = o$value;
    pl = new('plausibilityResult', par = o$par, value = min(P, 1), valueRaw = P, optimizer = o);
	return(pl);
}

plausibilityOptimizeGrid = function(this, objectiveFunction, start,
	# pars searchOptimum
	Ngrid = 3, epsX = 1e-2, epsY = 1e-3) {

	if (is.null(start)) start = plausibilityStart(this);
	delta = plausibilityDelta(this);
	r = approximateFunctionFor(objectiveFunction,
		start - delta, start + delta, Ngrid = Ngrid,
		conditions = conditionsMax(epsX = epsX, epsY = epsY), this = this
	);
	#o  = searchOptimum(objectiveFunction, start, this = this,
	#	delta = plausibilityDelta(this), gridGen = gridGen, eps = eps);
	#value = objectiveFunction(o, this = this);
	#pl = new('plausibilityResult', par = o, value = value, optimizer = list(par = o));
	#return(pl);
	Imax = which.max(r[, 'v']);
	P = r[Imax, 'v'];

	pl = new('plausibilityResult',
		par = r[Imax, -ncol(r)],
		value = min(P, 1),
		valueRaw = P,
		optimizer = list(grid = r, approximation = approximated(r)),
		object = this
		#optimizer = list(approximation = approximated(r))
	);
	return(pl);
}

plausibilityOptimizers = list(
	grid = plausibilityOptimizeGrid,
	optim = plausibilityOptimizeOptim
);



#
#	<p> region weighting function
#


weightingRegion = function(f0, f1, data, parF) {
	# pre-processing
	#y = data[[formula.response(f1)]];	# assume standardized response
	mm0 = model_matrix_from_formula(f0, data)$mm
	mm1 = model_matrix_from_formula(f1, data)$mm
	NsFree = setdiff(dimnames(mm1)[[2]], dimnames(mm0)[[2]]);
	NsF0 = dimnames(mm0)[[2]];
	NpF1 = length(dimnames(mm1)[[2]]);
	NpF0 = length(dimnames(mm0)[[2]]);
	# reorder columns
	Mm1 = cbind(mm1[, NsF0], mm1[, NsFree]);
	dimnames(Mm1)[[2]] = c(NsF0, NsFree);

	# MLEs
	r1 = glmModelBinomialMM(list(y = data$y, X = Mm1), Nbinom = 1);
	if (missing(parF)) parF = r1$par[(NpF0 + 1):NpF1];
	r0 = glmModelBinomialMM(list(y = data$y, X = Mm1), Nbinom = 1, parF = parF);
	return(list(ll1 = r1$ll, ll0 = r0$ll, lr = r1$ll - r0$ll));
}

weightingFunctionRegion = function(parF) {
	return(function(this, data, ...)-weightingRegion(this@f0, this@f1, data, ..., parF = parF)$lr);
}

plausibilityRegionOF = function(this, par) {
	pl = new('plausibilityModelBinomialWeightedClone', this, weightingFunction = weightingFunctionRegion(par))
	P = plausibilityOptimizeGrid(pl, cumProbSIcompBinomial, this@par)
	#print(P);
	return(P@value);
}

freeParameters = function(f0, f1, data) {
	ns0 = dimnames(model_matrix_from_formula(f0, data)$mm)[[2]];
	ns1 = dimnames(model_matrix_from_formula(f1, data)$mm)[[2]];
	setdiff(ns1, ns0);
}
nullParameters = function(f0, f1, data) {
	ns0 = dimnames(model_matrix_from_formula(f0, data)$mm)[[2]];
	ns1 = dimnames(model_matrix_from_formula(f1, data)$mm)[[2]];
	intersect(ns1, ns0);
}


# return vertices in lexicographic order
gridTensor = function(min = 0, max = 1, N = 3) {
	gridMarginal = lapply(seq_along(min), function(i) {
		seq(min[i], max[i], length.out = N)
	});
	grid = as.matrix(Rbind(merge.multi.list(gridMarginal)));
	return(grid);
}


evaluateFunctionCached = function(args, f, cache = NULL, ...) {
	if (is.null(cache)) return(apply(args, 1, f, ...));
	iAll = 1:nrow(args);
	iDupl = DfSearch(Df_(args), Df_(cache$args), returnIdcs = T);	# vertices already computed
	iNew = setdiff(iAll, iDupl);
	v = rep(NA, nrow(args));
	v[iDupl[, 2]] = cache$v[iDupl[, 1]];
	v[iNew] = sapply(iNew, function(i)f(args[i, ], ...));
	return(v);
}

evaluateCubes = function(g, v, Ngrid, cache, ...,
	f, epsX, epsY, contour) {
	Ndim = ncol(g);
	# iterate sub-cubes
	Ncubes = (Ngrid - 1)^Ndim;
	# spacing between indeces, tb cached
	cubeI = sapply(0:(Ndim - 1), function(i)Ngrid^i)
	# where are the indeces of the other cube corners
	cubeIdxSpacing = t_(cubeI * sapply(0:(2^Ndim - 1), ord2bin, digits = Ndim)) %*% rep(1, Ndim);

	rCubes = lapply(1:Ncubes, function(i) {
		# lower left corner
		cubeLL = (ord2base(i - 1, digits = Ndim, base = Ngrid - 1) %*% cubeI)[1, 1] + 1;
		# nodes belonging to current cube
		cubeIdcs = cubeLL + cubeIdxSpacing;
		# extract coordinates
		coords = g[cubeIdcs, , drop = F];
#if (any(min(apply(coords, 2, min)) > -4)) browser();
		approximateFunctionAtContour(f,
			mn = apply(coords, 2, min), mx = apply(coords, 2, max), Ngrid = Ngrid,
			cache = list(args = coords, v = v[cubeIdcs]), ...,
			epsX = epsX, epsY = epsY, contour = contour
		)
	});
	return(unique(do.call(rbind, rCubes)));
}
approximateFunctionAtContour = function(f, mn, mx, ...,
	epsX = 1e-3, epsY = 1e-2, Ngrid = 4, cache = NULL, contour = 0.05) {
	Ndim = length(mn);
	# <p> stopping condition
	if (any(abs(mx - mn) < epsX)) return(NULL);

	g = gridTensor(mn, mx, N = Ngrid);
	v = evaluateFunctionCached(g, f, cache, ...);
	#print(cbind(g, v))
	#print(rbind(mn, mx, c(min(v), max(v), rep(NA, Ndim-2)) ));
	# <p> stopping condition, do not stop on first run
	#dprint(max(v), min(v), rng= max(v) - min(v));
	if (notE(cache) && (
		max(v) < contour - epsY || min(v) > contour + epsY || max(v) - min(v) < epsY
	))
		return(cbind(g, v));

	return(evaluateCubes(g, v, Ngrid = Ngrid, cache = list(args = g, v = v), ...,
		f = f, epsX = epsX, epsY = epsY, contour = contour));
}

findRegion = function(f, fixed, range = c(-5, 5), Ngrid = 1e3, sel = function(x, contour = .05)x < contour) {
	Ivar = which(is.na(fixed));
	rng = seq(range[1], range[2], length.out = Ngrid);
	m = matrix(rep(fixed, length(rng)), ncol = length(fixed), byrow = T);
	m[, Ivar] = rng;
	y = apply(m, 1, f);
	Iregion = which(sel(y));
	return(list(lower = rng[min(Iregion)], upper = rng[max(Iregion)]));
}
