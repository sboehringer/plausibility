#
#	plausibilitiyStochastic.R
#Fri May 31 14:59:59 CEST 2019
#
#	derived from plausibilityPenalized.R

packageDefinition = list(
	name = 'plausibility',
	files = c('Rmeta.R', 'Rdata.R', 'Rsystem.R', 'Rfunctions.R', 'RpropertyList.R'),
	description = list(
		title = 'Plausibility based estimation and inference',
		author = 'Stefan B\uf6hringer <r-packages@s-boehringer.org>',
		description = 'R-package to implement important methods of the plausibility framework, including goodness-of-fit testing, test model comparisons and computation of marginal and joint plausibility regions.',
		depends = c(),
		suggests = c(),
		license = 'LGPL',
		news = "0.2-0	Implementation Plausibility region\n0.1-0   Initial release"
	),
	git = list(
		readme = '## Installation\n```{r}\nlibrary(devtools);\ninstall_github("sboehringer/plausibility")\n```\n',
		push = F,
		pushOnNewVersion = T,
		remote = 'https://github.com/sboehringer/plausibility.git'
	)
);

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
#   <p> S4 classes
#

dummyObjective = function(par, this)-Inf;
regionObjective = function(objective, level = .95) {
	return(function(par, this)(regionObjective(par, this) - (1 - level))^2);
}

setClass("plausibilityFamily", representation = list(
	f0 = 'formula',
	f1 = 'formula',
	data = 'data.frame',
	objectiveFunction = 'function'
), prototype = list(
	f0 = as.formula(NULL),
	f1 = as.formula(NULL),
	data = data.frame(),
	objectiveFunction = dummyObjective
));

setMethod('initialize', 'plausibilityFamily', function(.Object, f0, f1, data,
	objectiveFunction = dummyObjective) {
	Data = DfNames2std(data, vars.as.rhs(formula.response(f0)), ~ y);
	f0 = formula.add.response(formula.expand(f0, data), y ~ 1); # expand ~ ., std outcome name
	f1 = if (notE(f1)) formula.add.response(formula.expand(f1, data), y ~ 1) else NULL;

	#.Object@f0 = formula.add.response(f0, y ~ 1);
	#.Object@f1 = if (length(f1) > 0) formula.add.response(f1, y ~ 1) else as.formula(NULL);
	#Data = data;
	.Object@f0 = f0;
	.Object@f1 = if (length(f1) > 0) f1 else as.formula(NULL);

	.Object@data = Data;
	.Object@objectiveFunction = objectiveFunction;
	return(.Object);
});

setClass("plausibilityResult", representation = list(
	par = 'numeric',
	value = 'numeric',
	optimizer = 'list',
	object = 'plausibilityFamily'
), prototype = list(
	par = as.numeric(NA),
	value = as.numeric(NA),
	optimizer = list()
));


setGeneric("plausibility", function(this, optMethod = 'grid', Niter = 1L, ...)this)
setGeneric("plausibilityStart", function(this)this)
setGeneric("plausibilityDelta", function(this)NULL)
setGeneric("region", function(this, level = .95, ...)this)

setClass("plausibilityFamilySI", contains = 'plausibilityFamily', representation = list(
	Nsi = 'integer'
), prototype = list(
	Nsi = 1e3L
));
setMethod('initialize', 'plausibilityFamilySI', function(.Object, f0, f1, data, objectiveFunction, Nsi) {
	.Object = callNextMethod(.Object, f0, f1, data, objectiveFunction);
	return(.Object);
});

#
#	<p> glm models
#


setClass("plausibilityGlm", contains = 'plausibilityFamilySI', representation = list(
	mm = 'matrix'		# model matrix null
), prototype = list());

setMethod('initialize', 'plausibilityGlm', function(.Object, f0, f1, data, objectiveFunction, Nsi) {
	.Object = callNextMethod(.Object, f0, f1, data, objectiveFunction, Nsi);

	.Object@mm = model_matrix_from_formula(.Object@f0, .Object@data)$mm
	return(.Object);
});

#	<p> debugging
# Iout = which(exp(eventWeights) > 40)
# Nevents = apply(this@sim[, Iout], 2, sum)

# parF: fixed part of parameter vector
glmLLmm = function(par, X, y, Nbinom = 1, parF = c()) {
	parFull = c(par, parF);
	ps = plogis(X %*% parFull)[, 1];
	lli = dbinom(y, Nbinom, ps, log = T);
	#ll = sum(lli) + lchoose(nrow(mm), sum(y));
	ll = sum(lli);
	return(ll);
}
glmLL = function(par, data, f, ...)glmLLmm(par, model_matrix_from_formula(f, data)$mm, ...);

glmModelBinomial = function(f, data, Nbinom = 1) {
	m = glm(f, data, family = 'binomial');
	par = coefficients(summary(m))[, 'Estimate'];
	mm = model_matrix_from_formula(f, data)$mm;
	y = data[[formula.response(f)]];
	ll = glmLLmm(par, mm, y = y, Nbinom = Nbinom);

	if (ll == -Inf) {
		o = try(optim(rep(0, length(par)), glmLLmm,
			method = 'BFGS', control = list(fnscale = -1),
			X = mm, y = data$y, Nbinom = Nbinom), silent = T);
		if (class(o) == 'try-error') browser();
		par = o$par;
		ll = o$value;
	}
	return(list(model = m, par = par, ll = ll, sds = coefficients(summary(m))[, 'Std. Error'], formula = f));
}

# equivalence classes for rows
uniqueEqCl = function(m) {
	mu = t(unique(m));
	return(apply(m, 1, function(r)which(apply(mu == r, 2, all))));
}
# binomial coefficient for identical rows
modelMatrixBinom = function(m) {
	cls = uniqueEqCl(m);
}
weightingFunctionLR = function(this, data, ...) {
	-glmCompBinomial(this@f0, this@f1, data, ...)$lr;
}


setClass('plausibilityBinomial', contains = 'plausibilityGlm', representation = list(
	Nbinom = 'integer',	# #{binomial outcomes} == 1 for standard logistic regression
	par = 'numeric',	# starting par (MLE under f0, plausibility estimate)
	mdl0 = 'list',		# glm model + corresponing ll
	sim = 'matrix',		# stochastic outcomes under f0, N x Nsi
	pSI = 'numeric'	# 
), prototype = list(Nbinom = 1L));

cumProbSIbinomial = function(par, this, epsPs = 1e-5) {
	N = length(this@data$y);
	Nsim = ncol(this@sim);
	ps = plogis(this@mm %*% par)[, 1];
	ps[ps == 1] = 1 - epsPs;	# help optimizer
	TsRaw = apply(this@sim, 2, function(r)sum(dbinom(r, this@Nbinom, ps, log = T)));

	# <p> Test statistic + IS correction: reweight events
	Tdata = sum(dbinom(this@data$y, this@Nbinom, ps, log = T));
	sel = TsRaw <= Tdata;
	Ts = TsRaw[sel] - this@pSI[sel];

	# <p> P-value
	P = if (any(Ts == -Inf)) .5 / (Nsim + 1) else (sumExp(Ts) + 1) / (Nsim + 1);
	#print(t(c(par, P, sum(sel))));
	return(P);
}

setMethod('plausibilityDelta', 'plausibilityBinomial', function(this)this@mdl0$sds * 6)
setGeneric("plausibilityStart", function(this)this@par)

setMethod('initialize', 'plausibilityBinomial', function(.Object, f0, data, Nsi = 1e3L, Nbinom = 1L,
	start = NULL, weightingFunction = weightingFunctionLR, objectiveFunction = cumProbSIbinomial, sim = NULL, sampleFromAlt = TRUE) {
	Data = dataBinomialSplit(formula.response(f0), data, Nbinom);
	.Object = callNextMethod(.Object, f0, NULL, Data, objectiveFunction, Nsi);
	.Object@Nbinom = Nbinom;
	.Object@mdl0 = suppressWarnings(glmModelBinomial(.Object@f0, .Object@data));
	.Object@par = .Object@mdl0$par;

	par = if (sampleFromAlt) .Object@mdl0$par else .Object@par;
	ps = plogis(.Object@mm %*% par)[, 1];
	# <p> stochastic integration sample
	if (is.null(sim)) {
		Msim = matrix(runif(length(ps) * Nsi), ncol = Nsi);
		sim = Msim <= ps;	# column-wise comparison
		mode(sim) = 'integer';
	}
	.Object@sim = sim;

	# <p> probabilities stochastic integration sample under starting pars, needed for importance correction
	.Object@pSI = apply(sim, 2, function(r)sum(dbinom(r, Nbinom, ps, log = T)));
	return(.Object);
});

setClass('plausibilityBinomialWeighted', contains = 'plausibilityFamilySI', representation = list(
	Nbinom = 'integer',	# #{binomial outcomes} == 1 for standard logistic regression
	par = 'numeric',	# starting par (MLE under f0, plausibility estimate)
	mdl = 'list',		# glm model + corresponing ll
	mdl0 = 'list',		# glm model + corresponing ll
	mm = 'matrix',		# model matrix null
	mm1 = 'matrix',		# model matrix alternative
	Nsi = 'integer',	# #{stochastic integration steps}
	sim = 'matrix',		# stochastic outcomes under f0, N x Nsi
	pSI = 'numeric',	# 
	lbcData = 'numeric',# log binomial coefficient (data)
	lbcSim = 'numeric',	# log binomial coefficient (stochastic sample)
	weights = 'numeric',# weighing function applied to integeration sample
	weight = 'numeric',	# weight of data
	dataE = 'data.frame'	# expanded data (bernoulli outcome)
), prototype = list(
	Nbinom = 1L,
	mdl = list(),
	mm = matrix(numeric(), ncol = 1L),
	mm1 = matrix(numeric(), ncol = 1L),
	Nsi = 1e3L,
	sim = matrix(integer(), ncol = 1e3L),
	pSI = numeric(),
	lbcData = numeric(),
	lbcSim = numeric()
));
setMethod('plausibilityDelta', 'plausibilityBinomialWeighted', function(this)this@mdl0$sds * 2)

cumProbSIcompBinomial = function(par, this) {
	#print(par);
	N = length(this@data$y);
	Nsim = ncol(this@sim);
	ps = plogis(this@mm %*% par)[, 1];
	TsRaw = apply(this@sim, 2, function(r)sum(dbinom(r, this@Nbinom, ps, log = T)));
	# IS correction: reweight events
	eventWeights = (TsRaw - this@pSI)[this@weight > this@weights];
	#eventWeights = rep(0, Nsim)[this@weight > this@weights];
	P = if (length(eventWeights) == 0) 1/(Nsim + 1) else (sumExp(eventWeights) / Nsim);
	#dprint(par, Nsel = sum(this@weight > this@weights), sumExp = P * Nsim, P);
	return(P);
}
#	<p> debugging
# Iout = which(exp(eventWeights) > 40)
# Nevents = apply(this@sim[, Iout], 2, sum)

# parF: fixed part of parameter vector
glmLLmm = function(par, X, y, Nbinom = 1, parF = c()) {
	parFull = c(par, parF);
	ps = plogis(X %*% parFull)[, 1];
	lli = dbinom(y, Nbinom, ps, log = T);
	#ll = sum(lli) + lchoose(nrow(mm), sum(y));
	ll = sum(lli);
	return(ll);
}
glmLL = function(par, data, f, ...)glmLLmm(par, model_matrix_from_formula(f, data)$mm, ...);

dataBinomialSplit = function(y, data, Nbinom = 2) {
	dataBlocks = by(data, factor(1:nrow(data)), function(d) {
		d1 = rep.each.row(d, Nbinom);
		d1[[y]] = c(rep(0, Nbinom - d[[y]]), rep(1, d[[y]]));
		return(d1);
	});
	return(do.call(rbind, dataBlocks));
}

glmModelBinomial = function(f, data, Nbinom = 1) {
	response = formula.response(f);
	dataE = if (Nbinom > 1L) dataBinomialSplit(response, data, Nbinom) else data;
	m = glm(f, dataE, family = 'binomial');
	par = coefficients(summary(m))[, 'Estimate'];
	mm = model_matrix_from_formula(f, data)$mm;
	y = data[[response]];
	ll = glmLLmm(par, mm, y = y, Nbinom = Nbinom);

	if (ll == -Inf) {
		o = try(optim(rep(0, length(par)), glmLLmm,
			method = 'BFGS', control = list(fnscale = -1),
			X = mm, y = data$y, Nbinom = Nbinom), silent = T);
		if (class(o) == 'try-error') browser();
		par = o$par;
		ll = o$value;
	}
	return(list(
		model = m, par = par, ll = ll, sds = coefficients(summary(m))[, 'Std. Error'],
		parsFree = Seq(1, ncol(mm))
	));
}

glmCompBinomial = cumProbSIcompModelBinomial = function(f0, f1, data, Nbinom = 1L) {
	m0 = suppressWarnings(glmModelBinomial(f0, data, Nbinom));
	m1 = suppressWarnings(glmModelBinomial(f1, data, Nbinom));
	return(list(m0 = m0, m1 = m1, lr = 2*(m1$ll - m0$ll), family = 'binomial'));
}

modelCharacteristics = function(f0, f1, data, Nbinom = 1L) {
	m0 = glmModelBinomial(f0, data, Nbinom);
	m1 = glmModelBinomial(f1, data, Nbinom);

	cfs = coefficients(summary(m1$m));
	sds = cfs[, 'Std. Error'];
	pars = cfs[, 'Estimate'];
	parsN = nullParameters(f0, f1, data);
	parsF = freeParameters(f0, f1, data);

	return(list(m0 = m0, m1 = m1, ll = m0$ll, family = 'binomial',
		# under full model
		sds = sds, pars = pars, parsNull = parsN, parsFree = parsF
	));
}

# equivalence classes for rows
uniqueEqCl = function(m) {
	mu = t(unique(m));
	return(apply(m, 1, function(r)which(apply(mu == r, 2, all))));
}
# binomial coefficient for identical rows
modelMatrixBinom = function(m) {
	cls = uniqueEqCl(m);
}
weightingFunctionLR = function(this, data, ...) {
	-glmCompBinomial(this@f0, this@f1, data, ..., Nbinom = this@Nbinom)$lr;
}

setMethod('initialize', 'plausibilityBinomialWeighted', function(.Object, f0, f1, data,
	Nsi = 1e3L, Nbinom = 1L,
	start = NULL, weightingFunction = weightingFunctionLR, objectiveFunction = cumProbSIcompBinomial,
	sim = NULL, sampleFromAlt = TRUE) {
	dataC = completeData(f1, data);
	.Object = callNextMethod(.Object, f0, f1, dataC, objectiveFunction, Nsi);

	.Object@mdl0 = suppressWarnings(glmModelBinomial(.Object@f0, .Object@data, Nbinom));
	.Object@mdl = mdl = modelCharacteristics(.Object@f0, .Object@f1, .Object@data, Nbinom);
	.Object@par = if (notE(start)) start else mdl$m0$par;
	#print(.Object@par);
	.Object@Nsi = Nsi;
	.Object@Nbinom = as.integer(Nbinom);
	.Object@mm = model_matrix_from_formula(.Object@f0, .Object@data)$mm;
	if (sampleFromAlt) {
		mm = model_matrix_from_formula(.Object@f1, .Object@data)$mm;
		ps = plogis(mm %*% mdl$m1$par)[, 1];
	} else {
		ps = plogis(.Object@mm %*% .Object@par)[, 1];
	}
	# <p> stochastic integration sample
	if (is.null(sim)) {
		Msim0 = matrix(runif(length(ps) * Nsi * Nbinom), ncol = Nbinom);
		Msim = Msim0 <= ps;	# column-wise comparison, row-recycling
		sim = matrix(apply(Msim, 1, sum), ncol = Nsi);
		mode(sim) = 'integer';
	}
	.Object@sim = sim;

	# <p> probabilities stochastic integration sample, we sample from alternative but want to 
	#	evaluate under the null hypothesis
	.Object@pSI = apply(sim, 2, function(r)sum(dbinom(r, Nbinom, ps, log = T)));
	#N = nrow(.Object@data);
	.Object@weights = apply(sim, 2, function(r) {
		weightingFunction(.Object, DfRepl(.Object@data, Df(y = r)))
	});
	.Object@weight = weightingFunction(.Object, .Object@data);
	return(.Object);
});
setMethod('plausibilityStart', 'plausibilityBinomialWeighted', function(this)return(this@par));

setClass('plausibilityBinomialWeightedClone', contains = 'plausibilityBinomialWeighted',
	representation = list(), prototype = list()
);

setMethod('initialize', 'plausibilityBinomialWeightedClone', function(.Object, other, weightingFunction) {
	.Object = callNextMethod(.Object,
		other@f0, other@f1, other@data,
		objectiveFunction = other@objectiveFunction, Nsi = other@Nsi, sim = other@sim
	);
	.Object@weights = apply(.Object@sim, 2, function(r)
		weightingFunction(.Object, DfRepl(.Object@data, Df(y = r)))
	);
	.Object@weight = weightingFunction(.Object, .Object@data);

	return(.Object);
});

#
#	<p> plausiblity methods
#

plausibilityOptimizeOptim = function(this, objectiveFunction, start,
	method = 'Nelder-Mead', control = list(startScale = 1.5)) {

	#if (start == 0) browser();
	o = Optim(start, objectiveFunction, method = method, control = control, this = this);
    pl = new('plausibilityResult', par = o$par, value = o$value, optimizer = o);
	return(pl);
}

plausibilityOptimizeGrid = function(this, objectiveFunction, start,
	# pars searchOptimum
	Ngrid = 3, epsX = 1e-2, epsY = 1e-3) {

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
	pl = new('plausibilityResult', par = r[Imax, -ncol(r)], value = r[Imax, 'v'],
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

#' Compute plausibility
#'
#' Take a plausibility-family object and compute the plausibility estimate and P-value
#'
#' @param this Plausibility family object
#' @param optMethod Optimization method to use. Either 'grid' (default) or 'optim'
#' @param Niter Number of iterations used when finding the plausibility estimate
#' @param ... Arguments passed to the optimizer
#'
#' @export plausibility
setMethod('plausibility', 'plausibilityFamily',
	function(this, start = NULL, optMethod = 'grid', ...) {
	optimizer = plausibilityOptimizers[[ optMethod ]];
	print(list(start = start));
	if (is.null(start)) start = plausibilityStart(this);
	rPl = optimizer(this, this@objectiveFunction, start, ...);
	print(list(par = rPl@par, value = rPl@value));
	return(rPl);
});

#' Compute weighted plausibility for a model comparison of binomial outcome
#'
#' A model comparison of binomial outcome for nested models is computed. The plausibility estimate and P-value are returned together with a glm-fit.
#'
#' @param f0 Null model given as R formula
#' @param f1 Alternative model given as R formula
#' @param data data frame
#' @param optMethod Optimization method to use. Either 'grid' (default) or 'optim'
#' @param Niter Number of iterations used when finding the plausibility estimate
#' @param Nsi Number of stochastic integration samples (defaults to 1e3L)
#' @param Nbinom Number of repetions in the binomial outcome
#' @param ... Arguments passed to the optimizer
#'
#' @export PlausibilityBinomialWeighted
PlausibilityBinomialWeighted = function(f0, f1, data, Nsi = 1e3L, Niter = 2L, ...,
	optMethod = 'grid', Nbinom = 1L) {
	start = NULL;
	for (i in Seq(1, Niter)) {
		pl = new('plausibilityBinomialWeighted',
			f0 = f0, f1 = f1, data = data,
			Nsi = Nsi, Nbinom = Nbinom);
		rPl = plausibility(pl,  start = start, ..., optMethod = optMethod);
		start = rPl@par;
	}
	return(rPl);
}

#
#	<p> region weighting function
#


glmModelBinomialMM = function(d, Nbinom = 1, parF = c()) {
	# only works for Nbinom == 1
	m = glm.fit(d$X[, 1:(ncol(d$X) - length(parF)), drop = F], d$y, family = binomial());
	par = m$coefficients;

	o = try(optim(
		par, glmLLmm, method = 'BFGS', control = list(fnscale = -1),
		X = d$X, y = d$y, Nbinom = Nbinom, parF = parF
	), silent = T);
	#if (class(o) == 'try-error') browser();
	par = o$par;
	ll = o$value;

	return(list(model = m, par = par, ll = ll));
}

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
	pl = new('plausibilityBinomialWeightedClone', this, weightingFunction = weightingFunctionRegion(par))
	P = plausibilityOptimizeGrid(pl, cumProbSIcompBinomial, this@par)
	print(P);
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

setMethod('region', 'plausibilityFamily',
	function(this, start, level = .95, ...) {

	# search aroudn 3 SDs for a level of .95
	# <i> calibrate to level
	mn = this@mdl0$par[this@mdl0$parsFree];
	delta = 3 * this@mdl0$sds[this@mdl0$parsFree];
	print(c(mn, delta));
	fa = approximateFunctionFor(this@objectiveFunction, mn - delta, mn + delta,
		Ngrid = 7, conditions = conditionsContour(contour = 1 - level),
		this = this
	);
	fx = approximated(fa, Nnn = 2, aggr = max);
	ci = findRegion(fx, fixed = rep(NA, ncol(fa) - 1), sel = function(x, contour = 1 - level)x < contour);

	return(list(region = ci, fa = fa));
});

RegionBinomial = function(f0, data, Nsi = 1e3L, Niter = 1L, ..., optMethod = 'grid', level = .95) {
	rPl = PlausibilityBinomial(f0, data, Nsi = Nsi, Niter = Niter, ..., optMethod = optMethod);
	plRegion = region(rPl@object, start = rPl@start, level = level);
	return(plRegion);
}



