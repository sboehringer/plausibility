#
#	plausibilitiyStochastic.R
#Fri May 31 14:59:59 CEST 2019
#
#	derived from plausibilityPenalized.R

packageDefinition = list(
	name = 'plausibility',
	files = c('~/src/Rprivate/Rmeta.R', '~/src/Rprivate/Rdata.R', '~/src/Rprivate/Rsystem.R', '~/src/Rprivate/Rfunctions.R', '~/src/Rprivate/RpropertyList.R'),
	description = list(
		title = 'Plausibility based estimation and inference',
		author = 'Stefan B\uf6hringer <r-packages@s-boehringer.org>',
		description = 'R-package to test model comparisons and to compute marginal and joint plausibility regions using the plausibility framework.',
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
packageDefinitionTemplate = list(
	name = 'package',
	files = c('Rmeta.R', 'Rdata.R', 'Rsystem.R', 'Rfunctions.R', 'RpropertyList.R'),
	instFiles = list(Rscripts = 'Dev/pkg-minimal.R'),
	testing = list(
		doInstall = TRUE,
		tests = c('RtestsPackages/package/package.R')
	),
	description = list(
		title = 'Create packages from R-code directly',
		# version to be documented in news section
		#version = '0.1-0',
		author = 'Stefan B\uf6hringer <r-packages@s-boehringer.org>',
		description = 'This package simplifies package generation by automating the use of `devtools` and `roxygen`. It also makes the development workflow more efficient by allowing ad-hoc development of packages. Use `?"package-package"` for a tutorial or visit the project wiki (belonging to the source repository).',
		depends = c('roxygen2', 'devtools', 'methods'),
		suggests = c('testme', 'jsonlite', 'yaml', 'knitr'),
		news = "0.9-0	Finished vignette. Clean test. RC1\n0.8-1	Bug fix automatic dependency addition.\n0.8-0	Vignette building finished. Project vignette in progress.\n0.7-0	Vignette building. Started vignette for package `package`\n0.6-0	Clean CRAN check\n0.5-1	Resolved documentation\n0.5-0	Error free CRAN check. Warnings left.\n0.4-4	bugfix NAMESPACE generation\n0.4-3	`createPackage` fully documented.\n0.4-2	More documentation\n0.4-1	Bug fix NEWS file\n0.4-0	Self-contained example\n0.3-1	bug fix for missing files\n0.3-0	Beta, self-contained\n0.2-0	Alpha version\n0.1-0	Initial release",
		license = 'LGPL-2',
		vignettes = "vignettes/vignette-package.Rmd"
	)
);

#
#   <p> S4 classes
#

dummyObjective = function(par, this)-Inf;
regionObjective = function(objective, level = .95) {
	return(function(par, this)(regionObjective(par, this) - (1 - level))^2);
}

setClass("plausibilityFamily", representation = list(
	objectiveFunction = 'function',
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
	f1 = formula.add.response(formula.expand(f1, data), y ~ 1);

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
	optimizer = 'list'
), prototype = list(
	par = as.numeric(NA),
	value = as.numeric(NA),
	optimizer = list()
));


setGeneric("plausibility", function(this, optMethod = 'grid', Niter = 1L, ...)this)
setGeneric("plausibilityStart", function(this)this)
setGeneric("region", function(this, level = .95, ...)this)

setClass("plausibilityFamilySI", contains = 'plausibilityFamily', representation = list(
	Nsi = 'integer'
), prototype = list(
	Nsi = 1e3L
));
setMethod('initialize', 'plausibilityFamilySI', function(.Object, f0, f1, data, objectiveFunction, Nsi) {
	.Object = callNextMethod(.Object, f0, f1, data, objectiveFunction);
	.Object@Nsi = as.integer(Nsi);
	return(.Object);
});

setClass('plausibilityBinomial', contains = 'plausibilityFamilySI', representation = list(
	Nbinom = 'integer',	# #{binomial outcomes} == 1 for standard logistic regression
	par = 'numeric',	# starting par (MLE under f0, plausibility estimate)
	mdl = 'list',		# glm model + corresponing ll
	mm = 'matrix',		# model matrix
	Nsi = 'integer',	# #{stochastic integration steps}
	sim = 'matrix',		# stochastic outcomes under f0, N x Nsi
	pSI = 'numeric',	# 
	lbcData = 'numeric',# log binomial coefficient (data)
	lbcSim = 'numeric',	# log binomial coefficient (stochastic sample)
	weights = 'numeric',# weighing function applied to integeration sample
	weight = 'numeric'	# weight of data
), prototype = list(
	Nbinom = 1L,
	mdl = list(),
	mm = matrix(numeric(), ncol = 1L),
	Nsi = 1e3L,
	sim = matrix(integer(), ncol = 1e3L),
	pSI = numeric(),
	lbcData = numeric(),
	lbcSim = numeric()
));

cumProbSIcompBinomial = function(par, this) {
	N = length(this@data$y);
	Nsim = ncol(this@sim);
	ps = plogis(this@mm %*% par)[, 1];
	TsRaw = apply(this@sim, 2, function(r)sum(dbinom(r, this@Nbinom, ps, log = T)));
	# IS correction: reweight events
	eventWeights = (TsRaw - this@pSI)[this@weight > this@weights];
	#eventWeights = rep(0, Nsim)[this@weight > this@weights];
	P = if (length(eventWeights) == 0) 1/(Nsim + 1) else (sumExp(eventWeights) / Nsim);
	#dprint(par, Nsel = sum(this@weight > this@weights), sumExp = P * Nsim, P);
	#if (P > 1) browser();
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
	return(list(model = m, par = par, ll = ll, sds = coefficients(summary(m))[, 'Std. Error']));
}

glmCompBinomial = cumProbSIcompModelBinomial = function(f0, f1, data) {
	m0 = suppressWarnings(glmModelBinomial(f0, data));
	m1 = suppressWarnings(glmModelBinomial(f1, data));
	return(list(m0 = m0, m1 = m1, lr = 2*(m1$ll - m0$ll), family = 'binomial'));
}

modelCharacteristics = function(f0, f1, data) {
	m0 = glmModelBinomial(f0, data);
	m1 = glmModelBinomial(f1, data);

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
weightingFunctionLR = function(this, data, ...)-glmCompBinomial(this@f0, this@f1, data, ...)$lr;

setMethod('initialize', 'plausibilityBinomial', function(.Object, f0, f1, data, Nsi = 1e3L, Nbinom = 1L,
	start = NULL, weightingFunction = weightingFunctionLR, objectiveFunction = cumProbSIcompBinomial, sim = NULL) {
	.Object = callNextMethod(.Object, f0, f1, data, objectiveFunction, Nsi);

	.Object@mdl = mdl = modelCharacteristics(.Object@f0, .Object@f1, .Object@data);
	.Object@par = par = if (notE(start)) start else mdl$m0$par;
	#print(.Object@par);
	.Object@Nsi = Nsi;
	.Object@mm = mm = model_matrix_from_formula(.Object@f0, .Object@data)$mm;
	ps = plogis(mm %*% par)[, 1];

	# <p> stochastic integration sample
	if (is.null(sim)) {
		Msim = matrix(runif(length(ps) * Nsi), ncol = Nsi);
		sim = Msim <= ps;	# column-wise comparison
		mode(sim) = 'integer';
	}
	.Object@sim = sim;

	# <p> probabilities stochastic integration sample
	.Object@pSI = apply(sim, 2, function(r)sum(dbinom(r, Nbinom, ps, log = T)));
	N = nrow(.Object@data);
	.Object@lbcData = lchoose(N, sum(.Object@data$y));
	.Object@lbcSim = apply(sim, 2, function(r)lchoose(N, sum(r)));
	.Object@weights = apply(sim, 2, function(r)
		weightingFunction(.Object, DfRepl(.Object@data, Df(y = r)))
	);
	.Object@weight = weightingFunction(.Object, .Object@data);

	return(.Object);
});
setMethod('plausibilityStart', 'plausibilityBinomial', function(this)return(this@par));

setClass('plausibilityBinomialClone', contains = 'plausibilityBinomial',
	representation = list(), prototype = list()
);

setMethod('initialize', 'plausibilityBinomialClone', function(.Object, other, weightingFunction) {
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
	searchDelta = .5, gridGen = gridBounding(Ngrid = 3), eps = 1e-4) {

	o  = searchOptimum(objectiveFunction, start, this = this,
		delta = this@mdl$m0$sds * 2, gridGen = gridGen, eps = eps);
	value = objectiveFunction(o, this = this);
	pl = new('plausibilityResult', par = o, value = value, optimizer = list(par = o));
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
	function(this, optMethod = 'grid', Niter = 1L, ...) {
	optimizer = plausibilityOptimizers[[ optMethod ]];

	rPl = NULL;
	for (i in Seq(1, Niter)) {
		LogS(4, 'Iteration %{i}d');
		start = if (i == 1) plausibilityStart(this) else rPl$par;
		rPl = optimizer(this, this@objectiveFunction, start, ...);
	}
	return(rPl);
});

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
	pl = new('plausibilityBinomialClone', this, weightingFunction = weightingFunctionRegion(par))
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


setMethod('region', 'plausibilityFamily',
	function(this, start, level = .95, ...) {

	# search aroudn 3 SDs for a level of .95
	# <i> calibrate to level
	delta = 4 * this@mdl$sds[this@mdl$parsFree];

	ct = searchContour(this@objectiveFunction, this@mdl$pars[this@mdl$parsFree], this = this,
		contour = 1 - level,
		delta = delta, gridGen = gridBounding(Ngrid = 3), eps = 1e-2, lower = TRUE
	);
	return(ct);
});
