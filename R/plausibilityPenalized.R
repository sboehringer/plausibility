#
#	plausibilityPenalized.R
#Thu Jun 17 18:49:49 CEST 2021

require('glmnet');

#
#	<p> helper functions
#

foldIds = function(N, Nfolds, doSample = TRUE) {
	ids = if (doSample) sample(1:Nfolds) else 1:Nfolds;
	Nfold = splitSeatsForFractions(N, vn(rep(1, Nfolds)));
	foldIds = unlist(sapply(1:Nfolds, function(i)rep(ids[i], Nfold[i])));
	return(if (doSample) sample(foldIds) else foldIds);
}

#lambdaKey %in% c('lambda.min', 'lambda.1se')

lambdaAlphaSelect = function(X, y, lp0, model, NlambdaSel, nfolds, lambdaKey, alphaSel) {
	foldids = sapply(1:NlambdaSel, function(.) foldIds(nrow(X), nfolds));
	ns = c('lambda', 'mse');
	lambdaAlpha = lapply(alphaSel, function(alpha) {
		lambdas = lapply(1:NlambdaSel, \(i) {
			rGlmnet = cv.glmnet(X, y, offset = lp0, family = model@family,
				alpha = alpha, foldids = foldids[, i], standardize = FALSE);
			SetNames(c(rGlmnet[[lambdaKey]], mean(rGlmnet$cvm)), ns)
		});
	});
	# alpha-lambda-mse
	alm = aperm(array(unlist(lambdaAlpha), c(2, NlambdaSel, length(alphaSel))), c(2, 1, 3));
	dimnames(alm) = list(NULL, ns, alphaSel);
	# select lambda
	lambdas = apply(alm, 3, \(m)m[which.min((m[, 'lambda'] - median(m[, 'lambda']))^2), , drop = F]);
	dimnames(lambdas)[[1]] = ns;
	# alpha + lambda
	Ialpha = which.min(lambdas['mse', ]);
	rAlpha = c(alpha = alphaSel[Ialpha], lambda = lambdas['lambda', Ialpha]);
	return(rAlpha);
}

lambdaSelect = function(X, y, lp0, model, NlambdaSel, alpha, nfolds, lambdaKey) {
	if (is.null(lp0)) lp0 = rep(0, length(y));
	family = if (is.character(model)) model else model@family;
	lambdas = list.kpu(lapply(1:NlambdaSel, function(i) {
		cv.glmnet(X, y, offset = lp0, family = family, alpha = alpha, nfolds = nfolds, standardize = FALSE);
	}), lambdaKey);	# %in% c('lambda.min', 'lambda.min')
	lambda = median(lambdas);
	return(c(alpha = alpha, lambda = lambda));
}

plausibilityGlmnetTuning = function(y, mm0, X, model, alpha,
	NlambdaSel = 10, nfolds = 20, lambdaKey = 'lambda.1se', alphaSel = seq(0, 1, length.out = 11)) {
	# <p> create offset
	m0 = plausFit(model, y, mm0);
	lp0 =  (mm0 %*% m0$par)[, 1];

	# <p. select tuning paramters
	alphaLambda = if (is.na(alpha)) {
		lambdaAlphaSelect(X, y, lp0, model, NlambdaSel, nfolds, lambdaKey, alphaSel);
	} else {
		lambdaSelect(X, y, lp0, model, alpha, nfolds, lambdaKey)
	}
	return(alphaLambda);
}

plausibilityPenalizedFit = function(y, X, family = 'gaussian', offset = NULL,
	alpha = 1, lambda = 1e-3) {

	r = glmnet(X, y, family = family, offset = offset,
		alpha = alpha, lambda = lambda, standardize = FALSE)

	beta = coefficients(r)[, 1];
	lp = cbind(1, X) %*% beta + (if (notE(offset)) offset else 0);
	return(list(beta = beta, lp = lp, Nsel = sum(beta != 0) - 1));
}

# <!> cave: first arbument is plaubility object, not model object
weightingFunctionLRpenalized = function(p, y, mm0, X, alpha) {
	model = p@model;

	# <p> models
	m0 = plausFit(model, y, mm0);
	lp0 =  (mm0 %*% m0$par)[, 1];
	m1 = plausibilityPenalizedFit(y, X, family = model@family, offset = lp0, lambda = p@tuning['lambda'], alpha = p@tuning['alpha']);
	#print(m1$beta);
	ll1 = plausDensity(model, y, m1$lp, plausAncil(model, y, m1$lp))

	# ordering in favor of alternative has to be small values first
	return( ll1 - m0$ll );
}


#
#	<p> S4 classes
#

setClass("plausibilityPenalized", contains = 'plausibilityFamilyWeightedSI', representation = list(
	pSI = 'numeric',
	model = 'plausibilityModel',
	mdl0 = 'list',		# glm model + corresponing ll
	sim = 'matrix',
	weights = 'numeric',	# simulation weights
	weight = 'numeric',		# data weight
	tuning = 'numeric'		# mixing, penalty parameter
), prototype = list());

# alpha == NULL: perform model selection for alpha

setMethod('initialize', 'plausibilityPenalized', function(.Object, f0, f1, data, X, Nsi = 1e3, model,
	start = NULL, objectiveFunction = cumProbSIcomp, sim = NULL, NmaxIS = 5,
	weightingFunction = weightingFunctionLRpenalized, sampleFromAlt = TRUE, lp = NULL, fudge = NULL,
	NlambdaSel = 1e1, Nfolds = 10, standardize = TRUE, alpha = .5) {

	# <p> initialize
	#	f1 used to complete data, here: alternative given as matrix X
	.Object = callNextMethod(.Object, f0, f1 = f0, data, Nsi,
		objectiveFunction = objectiveFunction,
		weightingFunction = weightingFunction, NmaxIS = NmaxIS);

	if (standardize) X = standardize(X);

	#modelNm = Sprintf('plausibilityModel%{family}s');
	.Object@model = model;
	.Object@mdl0 = plausFit(model, .Object@y, .Object@mm);
	.Object@tuning = plausibilityGlmnetTuning(.Object@y, .Object@mm, X, model, alpha, NlambdaSel, Nfolds);
	print(.Object@tuning);

	lp0 =  (.Object@mm %*% .Object@mdl0$par)[, 1];
	# <p> sampling from alternative
	if (is.null(sim) && is.null(lp) && sampleFromAlt) {
		lp = plausibilityPenalizedFit(.Object@y, X, model@family,
			offset = lp0, lambda = .Object@tuning['lambda'], alpha = .Object@tuning['alpha'])$lp[, 1];
	}

	# <p> linear predictor
	if (is.null(lp)) lp = lp0;
	# <p> ancillary parameters
	parAncil = plausAncil(model, .Object@y, lp);
	# <p> stochastic integration sample	# plausSample(.Object, ., lp) %.% .
	if (is.null(sim))
		sim = apply(matrix(runif(length(lp) * Nsi), ncol = Nsi), 2, function(u)
			plausSample(model, u, lp, parAncil));

	# <p> data weight
	weight = weightingFunction(.Object, .Object@y, .Object@mm, X, .Object@tuning['alpha']);
	# <p> simulation weights
	weights = apply(sim, 2, function(y)weightingFunction(.Object, y, .Object@mm, X, .Object@tuning['alpha']));
	# filter simulation for required sub-sample
	.Object@sim = sim[, weights > weight, drop = FALSE];

	# <p> probabilities stochastic integration sample under starting pars, needed for importance correction
	.Object@pSI = apply(.Object@sim, 2, function(r)plausDensity(model, r, lp, parAncil));

	return(.Object);
});
setMethod('plausibilityStart', 'plausibilityPenalized', function(this)this@mdl0$par);
setMethod('plausibilityDelta', 'plausibilityPenalized', function(this)this@mdl0$sds * 2)

