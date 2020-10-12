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
		author = 'Stefan Böhringer <r-packages@s-boehringer.org>',
		description = 'R-package to test model comparisons and to compute marginal and joint plausibility regions using the plausibility framework.',
		depends = c(),
		suggests = c(),
		license = 'LGPL',
		news = "0.1-0   Initial release"
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
	f0 = formula.expand(f0, data); # expand ~ .
	f1 = formula.expand(f1, data);
	Data = DfNames2std(data, vars.as.rhs(formula.response(f0)), ~ y)

	.Object@f0 = formula.add.response(f0, y ~ 1);
	.Object@f1 = if (length(f1) > 0) formula.add.response(f1, y ~ 1) else as.formula(NULL);
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
	P = if (length(eventWeights) == 0) 1/(Nsim + 1) else sumExp(eventWeights) / Nsim;
	dprint(par, Nsel = sum(this@weight > this@weights), sumExp = sumExp(eventWeights), P);
	#if (P > 1) browser();
	return(P);
}
#	<p> debugging
# Iout = which(exp(eventWeights) > 40)
# Nevents = apply(this@sim[, Iout], 2, sum)

cumProbSIModelBinomial = function(f0, f1, data) {
	m0 = glmModelBinomial(f0, data);
	return(list(m0 = m0, ll = m0$ll, family = 'binomial'));
}
glmLLmm = function(par, mm, y, Nbinom = 1) {
	ps = plogis(mm %*% par)[, 1];
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
	ll = glmLLmm(par, mm, y = data$y, Nbinom = Nbinom);

	if (ll == -Inf) {
		o = optim(rep(0, length(par)), glmLLmm,
			method = 'BFGS', control = list(fnscale = -1),
			mm = mm, y = data$y, Nbinom = Nbinom);
		par = o$par;
		ll = o$value;
	}
	return(list(model = m, par = par, ll = ll));
}

glmCompBinomial = cumProbSIcompModelBinomial = function(f0, f1, data) {
	m0 = suppressWarnings(glmModelBinomial(f0, data));
	m1 = suppressWarnings(glmModelBinomial(f1, data));
	return(list(m0 = m0, m1 = m1, lr = 2*(m1$ll - m0$ll), family = 'binomial'));
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
weightingFunctionLR = function(f0, f1, data, ...)-glmCompBinomial(...);

setMethod('initialize', 'plausibilityBinomial', function(.Object, f0, f1, data, Nsi = 1e3L, Nbinom = 1L,
	start = NULL, weightingFunction = weightingFunctionLR) {
	.Object = callNextMethod(.Object, f0, f1, data, cumProbSIcompBinomial, Nsi);

	.Object@mdl = mdl = cumProbSIModelBinomial(.Object@f0, .Object@f1, .Object@data);
	.Object@par = par = if (notE(start)) start else mdl$m0$par;
	#print(.Object@par);
	.Object@Nsi = Nsi;
	.Object@mm = mm = model_matrix_from_formula(.Object@f0, .Object@data)$mm;

	# <p> stochastic integration sample
	ps = plogis(mm %*% par)[, 1];
	Msim = matrix(runif(length(ps) * Nsi), ncol = Nsi);
	sim = Msim <= ps;	# column-wise comparison
	mode(sim) = 'integer';
	.Object@sim = sim;

	# <p> probabilities stochastic integration sample
	.Object@pSI = apply(sim, 2, function(r)sum(dbinom(r, Nbinom, ps, log = T)));

	N = nrow(.Object@data);
	.Object@lbcData = lchoose(N, sum(.Object@data$y));
	.Object@lbcSim = apply(sim, 2, function(r)lchoose(N, sum(r)));
	.Object@weights = apply(sim, 2, function(r)
		-glmCompBinomial(.Object@f0, .Object@f1, DfRepl(.Object@data, Df(y = r)))$lr);
	.Object@weight = -glmCompBinomial(.Object@f0, .Object@f1, .Object@data)$lr;

	return(.Object);
});

setMethod('plausibilityStart', 'plausibilityBinomial', function(this)return(this@par));

#
#	<p> plausiblity methods
#

plausibilityOptimizeOptim = function(this, start,
	method = 'Nelder-Mead', control = list(startScale = 1.5)) {

	#if (start == 0) browser();
	o = Optim(start, this@objectiveFunction, method = method, control = control, this = this);
    pl = new('plausibilityResult', par = o$par, value = o$value, optimizer = o);
	return(pl);
}

plausibilityOptimizeGrid = function(this, start,
	# pars searchOptimum
	searchDelta = .5, gridGen = gridBounding(Ngrid = 3), eps = 1e-4) {

	o  = searchOptimum(this@objectiveFunction, start, this = this,
		delta = searchDelta, gridGen = gridGen, eps = eps);
	pl = new('plausibilityResult', par = o, value = as.numeric(NA), optimizer = list(par = o));
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
		rPl = optimizer(this, start, ...);
	}
	return(rPl);
});

setMethod('region', 'plausibilityFamily',
	function(this, level = .95, ...) {

	pl = plausibility(this, ...);


	browser();
});

# 
# #
# #	<p> generic plausibility functions
# #
# 
# #dmvnIndep = function(x, mu = 0, sd = 1, log = T)sum(dnorm(x, mu, sd, log = log));
# dmvnIndep = function(x, mu = 0, sd = 1, log = T, combine = if (log) sum else prod)
# 	combine(dnorm(x, mu, sd, log = log));
# 
# plausibilityNormLR = function(data) {
# 	m0 = lm(y ~ 1, data = data);
# 	m0S = summary(m0);
# 	m1 = lm(y ~ ., data = data);
# 	ll1 = sum(dnorm(data$y, predict(m1, data), summary(m1)$sigma, log = T));
# 	ll0 = sum(dnorm(data$y, coefficients(m0S)['(Intercept)', 'Estimate'], m0S$sigma, log = T));
# 	return(ll1 - ll0);
# }
# 
# plausibilityMLE = function(data, f0 = y ~ 1) {
# 	m0 = summary(lm(f0, data = data));
# 	mn0 = coefficients(m0)['(Intercept)', 'Estimate'];
# 	sd0 = m0$sigma;
# 	return(list(mean = mn0, sd = sd0));
# }
# 
# plausibilityNormLRsim = function(data, N = 1, Simplify = T, par = plausibilityMLE(data)) {
# 	y = sapply(1:N, function(i)rnorm(nrow(data), par$mean, par$sd));
# 	if (Simplify && dim(y)[2] == 1) y = avu(y);
# 	return(y);
# }
# 
# plausibilityNormLRint = function(data, par = plausibilityMLE(data)) {
# 	y = plausibilityNormLRsim(data, par = par);
# 	dataS = Df(Df_(data, min_ = 'y'), y = y);
# 	return(plausibilityNormLR(dataS));
# }
# 
# #sapplyWoI = function(v, f, ...)sapply(v, function(i, ...)f(...), ...)
# #applyWoI = function(v, f, ...)lapply(v, function(i, ...)f(...), ...)
# 
# plausibilityNorm = function(data, Nboot = 1e2, par = plausibilityMLE(data)) {
# 	Tdata = -plausibilityNormLR(data);
# 	Tint = -sapplyWoI(1:Nboot, plausibilityNormLRint, data = data, par = par);
# 	return(list(P = mean(Tdata < Tint), Tdata = Tdata, Tint = Tint));
# }
# 
# plausibilityNormLRintUnrolled = function(data, y) {
# 	dataS = Df(Df_(data, min_ = 'y'), y = y);
# 	return(plausibilityNormLR(dataS));
# }
# 
# plausibilityNormUnrolled = function(data, dataSim) {
# 	Tdata = -plausibilityNormLR(data);
# 	Tint = -apply(dataSim, 2, plausibilityNormLRintUnrolled, data = data);
# 	return(list(P = mean(Tdata < Tint), Tdata = Tdata, Tint = Tint));
# }
# 
# #
# #	<p> generic functions, high-dimensional (i.e. data y, X)
# #
# 
# dataLmNull = function(m0, X = NULL) {
# 	N = length(m0$residuals);
# 	if (is.null(X)) X = data.frame(dummy = rep(NA, N));
# 	y0 = rnorm(N);
# 	y1 = y0 / sd(y0) * summary(m0)$sigma + predict(m0, Df_(X));
# 	return(list(y = y1, X = X));
# }
# dataLmModel = function(d, f0 = y ~ 1) {
# 	m0 = lm(f0, data = Df(y = d$y));
# 	# optimization: detect intercept model
# 	d0 = if (length(f0[[3]]) == 1 && f0[[3]] == '1') Df(x = rep(NA, length(m0$residuals))) else Df_(d$X);
# 	mu0 = predict(m0, d0);
# 	return(list(model = m0, pred = mu0));
# }
# 
# 
# plausibilityNormLRintHD = function(data, ..., f0, plausibilityModel) {
# 	d0 = dataLmNull(dataLmModel(data, f0)$model, data$X);
# 	return(plausibilityModel(d0, ..., f0 = f0));
# }
# 
# plausibilityNormHD = function(data, ..., f0 = y ~ 1,
# 	Nboot = 1e2, plausibilityModel = plausibilityGlmnetRaw, key = 'T', digits = 6) {
# 
# 	Tdata = plausibilityModel(data, ..., f0 = f0);
# 	Tint = lapplyWoI(Seq(1, Nboot), plausibilityNormLRintHD,
# 		data = data, f0 = f0, plausibilityModel = plausibilityModel, ...);
# 	Tcomp = signif(-list.kpu(Tint, key), digits) <= signif(-Tdata[[key]], digits);
# 	return(list(P = mean(Tcomp), Tdata = Tdata, Tint = Tint));
# }
# 
# #
# #	<p> simulation
# #
# 
# CovMatrix = function(P, corr = 0) {
# 	Cov = diag(1, P);
# 	Cov[lower.tri(Cov)] = Cov[upper.tri(Cov)] = corr;
# 	return(Cov);
# }
# simulateHighDim = function(N = 5e3, p = 5e3, Nblocks = 10, corr = .9, beta = 0, sd = 1) {
# 	Np = floor(p / Nblocks);
# 	Cov = CovMatrix(Np, corr);
# 	Xs = lapply(1:Nblocks, function(i)mvrnorm(N, rep(0, Np), Cov));
# 	X = do.call(cbind, Xs);
# 	#Beta = recycle(beta, 1:p)[[1]];
# 	b0 = unlist(lapply(beta, function(b)c(b, rep(0, Np -1))));
# 	Beta = c(b0, rep(0, p - length(b0)));
# 	y = X %*% Beta + rnorm(N, sd = 1);
# 	data = list(y = y, X = X);
# 	return(data);
# }
# 
# #
# #	<p> penalized models
# #
# 
# plausibilityGlmnetRaw = function(d, glmnet_model = cv.glmnet, ..., f0 = y ~ 1) {
# 	# null
# 	m0 = dataLmModel(d, f0);
# 
# 	# alternative
# 	m1 = glmnet_model(d$X, d$y, ...);
# 	cfs = coefficients(m1);
# 	mu1 = predict(m1, d$X);
# 
# 	# model comparison
# 	res1 = d$y - mu1;
# 	ll1 = sum(dnorm(d$y, mu1, sd(res1), log = T));
# 	ll0 = sum(dnorm(d$y, m0$pred, summary(m0$model)$sigma, log = T));
# 	return(list(T = 2*(ll1 - ll0),
# 		Nsel = sum(cfs != 0) - 1,
# 		lambda.min = m1$lambda.min, lambda.1se = m1$lambda.1se)
# 	);
# }
# 
# plausibilityGlmnetLambda = function(d, NlambdaSel = 10, lambdaKey = 'lambda.min', ...) {
# 	m0 = lm(y ~ 1, data = Df(y = d$y));
# 	lambdas = list.kpu(lapply(1:NlambdaSel, function(i) {
# 		plausibilityGlmnetRaw(list(y = dataLmNull(m0, NULL)$y, X = d$X), ...);
# 	}), lambdaKey);	# %in% c('lambda.min', 'lambda.min')
# 	lambda = median(lambdas);
# 	return(lambda);
# }
# 
# plausibilityGlmnet = function(d, NlambdaSel = 10, Nboot = 50, f0 = y ~ 1, alpha = 1) {
# 	lambda = plausibilityGlmnetLambda(d, NlambdaSel, alpha = alpha);
# 
# 	r = plausibilityNormHD(d, glmnet_model = glmnet, lambda = lambda, f0 = f0,
# 		plausibilityModel = plausibilityGlmnetRaw, Nboot = Nboot, alpha = alpha);
# 	return(list(results = r, lambda = lambda));
# }
# 
# #
# #	<p> penalized models, general case
# #
# 
# dataGlmNullGaussian = function(m0, X = NULL) {
# 	N = length(m0$residuals);
# 	if (is.null(X)) X = data.frame(dummy = rep(NA, N));
# 	y0 = rnorm(N);
# 	y1 = y0 / sd(y0) * summary(m0)$sigma + predict(m0, Df_(X));
# 	return(list(y = y1, X = X));
# }
# # dataGlmNullBinomial = function(m0, X) {
# # 	N = length(m0$residuals);
# # 	if (is.null(X)) X = data.frame(dummy = rep(NA, N));
# # 	ps = predict(m0, Df_(X), type = 'response');
# # 	y = as.integer(runif(length(ps)) > ps);
# # 	return(Df(y = y));
# # }
# dataGlmNullBinomial = function(m0, X, Nbinom = 1) {
# 	N = length(m0$residuals);
# 	if (is.null(X)) X = data.frame(dummy = rep(NA, N));
# 	ps = predict(m0, Df_(X), type = 'response');
# 	# binomial case for Nbinom >= 1
# 	#Ps = t2r(sapply(ps, function(p)cumsum(dbinom(0:(Nbinom - 1), Nbinom, p))));
# 	#y = apply(apply(Ps, 1, function(P)R > P), 1, function(i)max(c(0, which(i))));
# 	# binomial case for Nbinom == 1
# 	y = as.integer(runif(length(ps)) > ps);
# 	return(Df(y = y));
# }
# 
# dataGlmNull = function(m0, X = NULL) {
# 	callDelegate('dataGlmNull', m0$family$family, list(m0 = m0, X = X));
# }
# 
# plGlmnetWeightGaussian = function(d, m0, m1) {
# 	mu1 = predict(m1, d$X);
# 	res1 = d$y - mu1;
# 	ll1 = sum(dnorm(d$y, mu1, sd(res1), log = T));
# 	ll0 = sum(dnorm(d$y, m0$pred, summary(m0$model)$sigma, log = T));
# 	return(2 * (ll1 - ll0));
# }
# plGlmnetWeightBinomial = function(d, m0, m1) {
# 	mu0 = predict(m0, Df_(d$X), type = 'response');
# 	mu1 = predict(m1, d$X, type = 'response');
# 	ll1 = sum(dbinom(d$y, 1, mu1, log = T));
# 	ll0 = sum(dbinom(d$y, 1, mu0, log = T));
# 	return(2 * (ll1 - ll0));
# }
# plGlmnetWeight = function(d, m0, m1) {
# 	callDelegate('plGlmnetWeight', m0$family$family, list(d = d, m0 = m0, m1 = m1));
# }
# 
# plausibilityGlmnetRawGlm = function(d, glmnet_model = cv.glmnet, ..., f0 = y ~ 1, family = 'gaussian') {
# 	# null
# 	d0 = if (length(f0[[3]]) == 1 && f0[[3]] == '1') Df(y = d$y) else Df(y = d$y, Df_(d$X));
# 	m0 = glm(f0, data = d0, family = family);
# 
# 	# alternative
# 	m1 = glmnet_model(d$X, d$y, ..., family = family);
# 
# 	# model comparison
# 	T = plGlmnetWeight(d, m0, m1);
# 	return(list(
# 		T = T, Nsel = sum(coefficients(m1) != 0) - 1,
# 		lambda.min = m1$lambda.min, lambda.1se = m1$lambda.1se)
# 	);
# }
# 
# plausibilityGlmnetLambdaGlm = function(d, NlambdaSel = 10, lambdaKey = 'lambda.min', ...,
# 	family = 'gaussian') {
# 	m0 = glm(y ~ 1, data = Df(y = d$y), family = family);
# 	lambdas = list.kpu(lapply(1:NlambdaSel, function(i) {
# 		plausibilityGlmnetRawGlm(list(y = dataGlmNull(m0, NULL)$y, X = d$X), ..., family = family);
# 	}), lambdaKey);	# %in% c('lambda.min', 'lambda.min')
# 	lambda = median(lambdas);
# 	return(lambda);
# }
# 
# plausibilityGlmnetGlmLambda = function(d, NlambdaSel = 10, f0 = y ~ 1,
# 	lambdaKey = 'lambda.min', ..., family = 'gaussian') {
# 	data = Df(y = d$y, Df_(d$X));
# 	m0 = glm(f0, data = data, family = family);
# 	lambdas = list.kpu(lapply(Seq(1, NlambdaSel), function(i) {
# 		dNull = list(y = dataGlmNull(m0, NULL)$y, X = d$X);
# 		plausibilityGlmnetRawGlm(dNull, ..., family = family);
# 	}), lambdaKey);	# %in% c('lambda.min', 'lambda.min')
# 	lambda = median(lambdas);
# 	return(lambda);
# }
# 
# plausibilityGlmLRintHD = function(data, ..., f0, plausibilityModel, family = 'gaussian') {
# 	# <N><o> detect intercept model and create smaller data frame
# 	d0 = if (length(f0[[3]]) == 1 && f0[[3]] == '1') Df(y = d$y) else Df(y = d$y, Df_(d$X));
# 	m0 = glm(f0, data = d0, family = family);
# 	d0 = dataGlmNull(m0, data$X);
# 	return(plausibilityModel(list(y = d0$y, X = data$X), ..., f0 = f0, family = family));
# }
# 
# plausibilityGlmHD = function(data, ..., f0 = y ~ 1,
# 	Nboot = 1e2, plausibilityModel = plausibilityGlmnetRawGlm, family = 'gaussian', key = 'T', digits = 6) {
# 
# 	Tdata = plausibilityModel(data, ..., f0 = f0, family = family);
# 	Tint = lapplyWoI(Seq(1, Nboot), plausibilityGlmLRintHD,
# 		data = data, f0 = f0, plausibilityModel = plausibilityModel, ..., family = family);
# 	Tcomp = signif(-list.kpu(Tint, key), digits) <= signif(-Tdata[[key]], digits);
# 	return(list(P = mean(Tcomp), Tdata = Tdata, Tint = Tint));
# }
# 
# # implicitely assume that f1 = f0 + all.vars(d$X)
# plausibilityGlmnetGlm = function(d, NlambdaSel = 10, Nboot = 50, f0 = y ~ 1,
# 	alpha = 1, family = 'gaussian', ...) {
# 	lambda = plausibilityGlmnetGlmLambda(d, NlambdaSel, f0 = f0, alpha = alpha, family = family);
# 
# 	r = plausibilityGlmHD(d, glmnet_model = glmnet, lambda = lambda, f0 = f0,
# 		plausibilityModel = plausibilityGlmnetRawGlm, Nboot = Nboot, alpha = alpha, family = family, ...);
# 	return(list(results = r, lambda = lambda));
# }
# 
# #
# #	<p> binomial model
# #
# 
# cumProbSIBinomial_old = function(par, dataSI, mm, d0) {V
# 	ps = plogis(mm %*% par);
# 	T = sum(dbinom(d0$y, 1, ps, log = T));
# 	Ts = apply(dataSI$sim, 2, function(r)dbinom(r, 1, ps, log = T));
# 	Ts = apply(2* Ts - dataSI$pSI, 2, sum);
# 	#Ts = apply(dataSI$sim, 2, function(r)sum(dbinom(r, 1, ps, log = T)));
# 	pl = mean(c(1, Ts <= T));
# 	#pl = mean(as.integer(Ts <= T) * exp(Ts));
# 	#dprint(par, pl);
# 	return(pl);
# }
# cumProbSIBinomialSI_old = function(par, f0, data, Nsi = 1e4, dataSI) {
# 	d0 = glmPlData2d0(f0, data);
# 	if (missing(dataSI)) dataSI = simSIBinomial(par, d0, Nsi);
# 	cumProbSIBinomial_old(par, dataSI, model_matrix_from_formula(f0, d0)$mm, d0)
# }
# 
# cumProbSIBinomial = function(par, dataSI, mm, d0, Nbinom = 1) {
# 	ps = plogis(mm %*% par);
# 	N = length(d0$y);
# 	#print(ps);
# 	T = sum(dbinom(d0$y, Nbinom, ps, log = T)) + dataSI$lbcData;
# 	#Tlc = lchoosemult(Table(d0$y, 0, Nbinom));
# 	TsRaw = apply(dataSI$sim, 2, function(r)sum(dbinom(r, Nbinom, ps, log = T)));
# 	# IS correction (M: matrix)
# 	Ts = 2*TsRaw - dataSI$pSI + dataSI$lbcSim;
# 	#TsLc = apply(dataSI$sim, 2, function(r)lchoosemult(Table(r, 0, Nbinom)));
# 	#print(sum(Ts0 - dataSI$pSI));
# 	#Ts = apply(dataSI$sim, 2, function(r)sum(dbinom(r, 1, ps, log = T)));
# 	#pl = mean(as.integer(Ts0 <= T) * exp(Ts));
# 	#pl = mean(exp(Ts[Ts0 <= T]));
# 	#pl = meanExp(Ts[Ts0 <= T]);
# 	#pl = sum(exp(Ts[Ts0 <= T])) / length(Ts);
# 	#pl = sumExp(Ts[Ts0 + TsLc <= T + Tlc]) / length(Ts);
# 	#pl = sumExp(Ts[Ts <= T]) / sumExp(Ts);
# 	pl = fraction(Ts <= T);
# 	#dprint(fraction = fraction(Ts <= T), pl = pl);
# 	#dprint(par, pl);
# 	#rbind(dataSI$sim[, which.max(Ts0)], exp(dbinom(dataSI$sim[, which.max(Ts0)], 1, ps, log = T)))
# 	#rbind(d0$y, dbinom(d0$y, 1, ps, log = T))
# dprint(par, pl);
# 	return(pl);
# }
# 
# 
# #
# #	<p> plausibility with stochastic integration
# #
# 
# simSIBinomial = function(par, d0, Nsi, ..., Nbinom = 1) {
# 	#ps = predict(m0, d0, type = 'response');
# 	ps = plogis(model_matrix_from_formula(y ~ 1, d0)$mm %*% par)[, 1];
# 	Msim = matrix(runif(length(ps) * Nsi), ncol = Nsi);
# 	sim = Msim <= ps;	# column-wise comparison
# 	mode(sim) = 'integer';
# 	pSI = apply(sim, 2, function(r)sum(dbinom(r, Nbinom, ps, log = T)));
# 	lbcData = lchoose(N, sum(d0$y));
# 	lbcSim = apply(sim, 2, function(r)lchoose(N, sum(r)));
# 	return(list(sim = sim, pSI = pSI, lbcData = lbcData, lbcSim = lbcSim));
# }
# 
# glmPlData2d0 = function(f0, data) {
# 	d0 = if (is.data.frame(data$y)) data$y else
# 		Df(y = if (is.factor(data$y)) as.integer(data$y) - 1 else data$y);
# }
# 
# cumProbSIBinomialSI = function(par, f0, data, Nsi = 1e4, dataSI) {
# 	d0 = glmPlData2d0(f0, data);
# 	if (missing(dataSI)) dataSI = simSIBinomial(par, d0, Nsi);
# 	cumProbSIBinomial(par, dataSI, model_matrix_from_formula(f0, d0)$mm, d0)
# }
# 
# 
# glmPlSI = function(f0, d0, par, family = 'binomial',
# 	# pars searchOptimum
# 	searchDelta = .5, gridGen = gridBounding(Ngrid = 3), eps = 1e-4,
# 	# stochastic integration
# 	Nsi = 1e4, useOptim = T, optimMethod = 'BFGS'
# 	) {
# 
# 	dataSI = callDelegate('simSI', family, list(par = par, d0 = d0, Nsi = Nsi))
# 	fn = get(Sprintf('cumProbSI%{family}u'));
# 
# 	plSI = if (!useOptim) searchOptimum(fn,
# 		start = par, delta = searchDelta, gridGen = gridGen, scale = 1, returnOpt = T, eps = eps,
# 		dataSI = dataSI, mm = model_matrix_from_formula(f0, d0)$mm, d0 = d0
# 	) else optim(par, fn, dataSI = dataSI, mm = mm, d0 = d0, method = optimMethod,
# 		control = list(fnscale = -1));
# 	return(plSI);
# }
# 
# # y: vector or data.frame (with nuisance covariates)
# glmPl_old = function(f0, f1 = NULL, data, family = 'binomial', Niter = 1,
# 	# pars searchOptimum
# 	searchDelta = .5, gridGen = gridBounding(Ngrid = 5), eps = 1e-4,
# 	# stochastic integration
# 	Nsi = 1e4, type = 'response'
# 	) {
# 
# 	d0 = glmPlData2d0(f0, data);
# 	m0 = glm(f0, d0, family = family);
# 	# <p> starting value
# 	par = coefficients(summary(m0))[, 'Estimate'];
# 
# 	for (i in Seq(1, Niter)) {
# 		LogS(3, 'Iteration %{i}d');
# 		plI = glmPlSI(f0, d0, par, family, searchDelta, gridGen, eps, Nsi);
# 		par = plI$par;
# 	}
# 
# 	return(list(o = plI, glm = summary(m0)));
# }
# 
# 
# #
# #	<p> start new unified approach
# #
# 
# #
# #	<p> helper functions
# #
# 
# glmLLmm = function(par, mm, y, Nbinom = 1) {
# 	ps = plogis(mm %*% par)[, 1];
# 	lli = dbinom(y, Nbinom, ps, log = T);
# 	#ll = sum(lli) + lchoose(nrow(mm), sum(y));
# 	ll = sum(lli);
# 	return(ll);
# }
# glmLL = function(par, data, f, ...)glmLLmm(par, model_matrix_from_formula(f, data)$mm, ...);
# 
# glmModelBinomial = function(f, data, Nbinom = 1) {
# 	m = glm(f, data, family = 'binomial');
# 	par = coefficients(summary(m))[, 'Estimate'];
# 	mm = model_matrix_from_formula(f, data)$mm;
# 	ll = glmLLmm(par, mm, y = data$y, Nbinom = Nbinom);
# 	if (ll == -Inf) {
# 		o = optim(rep(0, length(par)), glmLLmm,
# 			method = 'BFGS', control = list(fnscale = -1),
# 			mm = mm, y = data$y, Nbinom = Nbinom);
# 		par = o$par;
# 		ll = glmLLmm(par, mm, y = data$y, Nbinom = Nbinom);
# 	}
# 	return(list(model = m, par = par, ll = ll));
# }
# 
# #
# #	<p> model comparison
# #
# 
# glmCompBinomial = cumProbSIcompModelBinomial = function(f0, f1, data) {
# 	m0 = suppressWarnings(glmModelBinomial(f0, data));
# 	m1 = suppressWarnings(glmModelBinomial(f1, data));
# 	return(list(m0 = m0, m1 = m1, lr = 2*(m1$ll - m0$ll), family = 'binomial'));
# }
# glmCompBinomialAnova = function(f0, f1, data) {
# 	mdls = glmCompBinomial(f0, f1, data);
# 	return(c(mdls, list(anova = anova(mdls$m0$model, mdls$m1$model, test = 'Chisq'))));
# }
# 
# cumProbSIcompInitBinomial = function(f0, f1, d0, Nsi, ..., Nbinom = 1) {
# 	# refactored not to be called in 
# 	mdl = glmCompBinomial(f0, f1, d0);
# 	par = mdl$m0$par;
# 
# 	#ps = predict(m0, d0, type = 'response');
# 	ps = plogis(model_matrix_from_formula(f0, d0)$mm %*% par)[, 1];
# 	Msim = matrix(runif(length(ps) * Nsi), ncol = Nsi);
# 	sim = Msim <= ps;	# column-wise comparison
# 	mode(sim) = 'integer';
# 	pSI = apply(sim, 2, function(r)sum(dbinom(r, Nbinom, ps, log = T)));
# 	N = nrow(d0);
# 	lbcData = lchoose(N, sum(d0$y));
# 	lbcSim = apply(sim, 2, function(r)lchoose(N, sum(r)));
# 	weights = apply(sim, 2, function(r)glmCompBinomial(f0, f1, DfRepl(d0, Df(y = r)))$lr);
# 	return(list(sim = sim, pSI = pSI, lbcData = lbcData, lbcSim = lbcSim,
# 		weights = weights, weight = mdl$lr, mdl = mdl,
# 		PlStart = mean(mdl$lr < weights), start = par
# 	));
# }
# 
# cumProbSIcomp_Binomial = function(par, mm, data, sim, pSI, lbcData, lbcSim, 
# 	weights, weight, mdl, ..., Nbinom = 1) {
# 	ps = plogis(mm %*% par);
# 	N = length(data$y);
# 	#T = sum(dbinom(data$y, Nbinom, ps, log = T)) + lbcData;
# 	TsRaw = apply(sim, 2, function(r)sum(dbinom(r, Nbinom, ps, log = T)));
# 	# IS correction (M: matrix)
# 	Ts = 2*TsRaw - pSI + lbcSim;
# 	#pl = fraction(Ts <= T);
# 	#PlLog = cdfl(Ts, weights)(T);
# 	#return(exp(PlLog));
# 	pChange = (TsRaw - pSI)[weight < weights];
# 	#return(PlStart - mean(pChange));
# 	#return(PlStart - mean(pChange));
# 	#return(PlStart * exp(mean(pChange)));
# 	return(if (length(pChange) == 0) 0 else mean(pChange));
# }
# cumProbSIcompPostBinomial = function(pl, dataSI) {
# 	return(list(
# 	plausibilitiy = list(par = pl$par, value = exp(pl$value) * dataSI$PlStart, Nsi = ncol(dataSI$sim)),
# 	lh = list(
# 		glm = list(m1 = summary(dataSI$mdl$m1$model),
# 		anova = anova(dataSI$mdl$m0$model, dataSI$mdl$m1$model, test = 'Chisq')
# 	))));
# }
# 
# #
# #	<p> goodness-of-fit
# #
# 
# cumProbSIModelBinomial = function(f0, f1, data) {
# 	m0 = glmModelBinomial(f0, data);
# 	return(list(m0 = m0, ll = m0$ll, family = 'binomial'));
# }
# 
# cumProbSIInitBinomial = function(f0, f1, d0, Nsi, ..., Nbinom = 1) {
# 	# refactored not to be called in 
# 	mdl = cumProbSIModelBinomial(f0, f1, d0);
# 	par = mdl$m0$par;
# 
# 	#ps = predict(m0, d0, type = 'response');
# 	ps = plogis(model_matrix_from_formula(f0, d0)$mm %*% par)[, 1];
# 	Msim = matrix(runif(length(ps) * Nsi), ncol = Nsi);
# 	sim = Msim <= ps;	# column-wise comparison
# 	mode(sim) = 'integer';
# 	pSI = apply(sim, 2, function(r)sum(dbinom(r, Nbinom, ps, log = T)));
# 	N = nrow(d0);
# 	lbcData = lchoose(N, sum(d0$y));
# 	lbcSim = apply(sim, 2, function(r)lchoose(N, sum(r)));
# 	return(list(sim = sim, pSI = pSI, lbcData = lbcData, lbcSim = lbcSim, start = par, mdl = mdl,
# 		call = list(f0 = f0, f1 = f1, data = d0)));
# }
# 
# cumProbSIBinomial = function(par, mm, data, sim, pSI, lbcData, lbcSim, mdl, ..., Nbinom = 1) {
# 	ps = plogis(mm %*% par);
# 	N = length(data$y);
# 	T = sum(dbinom(data$y, Nbinom, ps, log = TRUE)) + lbcData;
# 	TsRaw = apply(sim, 2, function(r)sum(dbinom(r, Nbinom, ps, log = TRUE)));
# 	Ts = 2*TsRaw - pSI + lbcSim;
# # 	PlLog = cdfl(Ts)(T);
# # print(exp(logSumExp(TsRaw)));
# # print(PlLog);
# # 	return(exp(PlLog));
# 	Pl = mean(T <= Ts);
# 	return(Pl);
# }
# 
# #
# #	<p> generic functions
# #
# 
# glmPlPostStandard = function(pl, dataSI) {
# 	return(list(plausibilitiy = c(pl, list(Nsi = ncol(dataSI$sim))),
# 		lh = with(dataSI$call, glmCompBinomialAnova(f0, f1, data))));
# }
# 
# glmPlSIoptimize = function(f0, f1, d0, par,
# 	family = 'binomial', f = 'cumProbSIcomp_',
# 	# pars searchOptimum
# 	searchDelta = .5, gridGen = gridBounding(Ngrid = 3), eps = 1e-4, dataSI, ...,
# 	useOptim = T, optimMethod = 'Nelder-Mead'
# 	) {
# 
# 	of = get(Sprintf('%{f}s%{family}u'));
# 	mm = model_matrix_from_formula(f0, d0)$mm;
# 	# arguments searchOptimum; <A> start comes from dataSI
# 	a = list(of, delta = searchDelta, gridGen = gridGen, scale = 1, returnOpt = T, eps = eps);
# 	# arguments objective function
# 	aOf = c(dataSI, list(mm = mm, data = d0), list(...));
# 	plSI = if (!useOptim)
# 		do.call(searchOptimum, c(a, aOf)) else
# 		do.call(optim, c(
# 			list(par = dataSI$start, of),
# 			aOf,
# 			list(method = optimMethod, control = list(fnscale = -1))
# 		));
# 	plSI = rget(Sprintf('%{f}sPost%{family}u'), glmPlPostStandard)(plSI, dataSI);
# 	return(plSI);
# }
# 
# 
# 
# # y: vector or data.frame (with nuisance covariates)
# glmPlComp = function(f0, f1, data, family = 'binomial', Niter = 1,
# 	# pars searchOptimum
# 	searchDelta = .5, gridGen = gridBounding(Ngrid = 5), eps = 1e-4,
# 	# stochastic integration
# 	Nsi = 1e4, ...
# 	) {
# 	glmPlgeneric(f0, f1, data, family, 'cumProbSIcomp_', Niter, searchDelta, gridGen, eps, Nsi, ...)
# }
# 
# # y: vector or data.frame (with nuisance covariates)
# plausibilityGeneric = function(f0, f1, data, family = 'binomial', f = 'cumProbSIcomp_', Niter = 1,
# 	# pars searchOptimum
# 	searchDelta = .5, gridGen = gridBounding(Ngrid = 5), eps = 1e-4,
# 	# stochastic integration
# 	Nsi = 1e4, ...
# 	) {
# 
# 	for (i in Seq(1, Niter)) {
# 		LogS(4, 'Iteration %{i}d');
# 		dataSI = CallDelegate(c(f, 'Init'), family, f0 = f0, f1 = f1, d0 = data, Nsi = Nsi, ...);
# 		#plI = glmPlSIcomp(f0, f1, data, mdl$m0$par, family, searchDelta, gridGen, eps, Nsi);
# 		plI = glmPlSIoptimize(f0, f1, data, if (i == 1) dataSI$start else plI$par, family, f,
# 			searchDelta, gridGen, eps, dataSI = dataSI);
# 	}
# 	return(plI);
# }
# 
# plausibilityMethods = list(SI = 'cumProbSIcomp');
# 
# plausibility_ = function(f0, f1 = NULL, data, family = 'binomial', Niter = 1,
# 	# pars searchOptimum
# 	searchDelta = .5, gridGen = gridBounding(Ngrid = 5), eps = 1e-4,
# 	# stochastic integration
# 	Nsi = 1e4, method = 'SI', region = F, ...
# 	) {
# 
# 	#	<p> pre-processing
# 	# expand ~ . <i> -> move to more generic function
# 	if (all(as.character(formula.rhs(f1)) == c('~', '.'))) {
# 		covs = setdiff(names(data), as.character(formula.response(f1)));
# 		f1 = formula.set.rhs(f1, vars.as.rhs(covs));
# 	}
# 	Data = DfNames2std(data, vars.as.rhs(formula.response(f0)), ~ y)
# 	f0 = formula.add.response(f0, y ~ 1);
# 	f1 = formula.add.response(f1, y ~ 1);
# 
# 	plausibilityGeneric(f0, f1, Data, family,
# 		plausibilityMethods[method], Niter, searchDelta, gridGen, eps, Nsi, ...
# 	)
# }
# 
# 
