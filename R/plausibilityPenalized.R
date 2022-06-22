#
#	plausibilityPenalized.R
#Thu Jun 17 18:49:49 CEST 2021

require('glmnet');

#
#	<p> helper functions
#

#lambdaKey %in% c('lambda.min', 'lambda.1se')

plausibilityGlmnetLambda = function(y, mm0, X, model, alpha,
	NlambdaSel = 10, nfolds = 20, lambdaKey = 'lambda.1se') {
	# <p> create offset
	m0 = plausFit(model, y, mm0);
	lp0 =  (mm0 %*% m0$par)[, 1];

	# <p. select lambda
	lambdas = list.kpu(lapply(1:NlambdaSel, function(i) {
		cv.glmnet(X, y, offset = lp0, family = model@family, alpha = alpha, nfolds = nfolds, standardize = FALSE);
	}), lambdaKey);	# %in% c('lambda.min', 'lambda.min')
	#print(stem(lambdas));
	lambda = median(lambdas);
	return(lambda);
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
	m1 = plausibilityPenalizedFit(y, X, family = model@family, offset = lp0, lambda = p@lambda, alpha = alpha);
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
	lambda = 'numeric'		# penalty parameter
), prototype = list());

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
	.Object@lambda = plausibilityGlmnetLambda(.Object@y, .Object@mm, X, model, alpha, NlambdaSel, Nfolds);
	#print(list(lambda = .Object@lambda));

	lp0 =  (.Object@mm %*% .Object@mdl0$par)[, 1];
	# <p> sampling from alternative
	if (is.null(sim) && is.null(lp) && sampleFromAlt) {
		lp = plausibilityPenalizedFit(.Object@y, X, model@family,
			offset = lp0, lambda = .Object@lambda, alpha = alpha)$lp[, 1];
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
	weight = weightingFunction(.Object, .Object@y, .Object@mm, X, alpha);
	# <p> simulation weights
	weights = apply(sim, 2, function(y)weightingFunction(.Object, y, .Object@mm, X, alpha));

	# filter simulation for required sub-sample
	.Object@sim = sim[, weights > weight ];

	# <p> probabilities stochastic integration sample under starting pars, needed for importance correction
	.Object@pSI = apply(.Object@sim, 2, function(r)plausDensity(model, r, lp, parAncil));

	return(.Object);
});
setMethod('plausibilityStart', 'plausibilityPenalized', function(this)this@mdl0$par);
setMethod('plausibilityDelta', 'plausibilityPenalized', function(this)this@mdl0$sds * 2)

