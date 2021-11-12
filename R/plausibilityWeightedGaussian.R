#
#	plausbilityWeightedGaussian.R
#Thu Jun 17 17:11:33 CEST 2021

#
#	<p> standard weighted plausbility
#

setClass("plausibilityGlmWeightedGaussian", contains = 'plausibilityGlmWeighted', representation = list(
));

setMethod('initialize', 'plausibilityGlmWeightedGaussian', function(.Object, f0, f1, data, Nsi = 1e3, model,
	start = NULL, objectiveFunction = cumProbSIcomp, sim = NULL, NmaxIS = 5,
	weightingFunction = weightingFunctionLR, sampleFromAlt = TRUE, lp = NULL, fudge = 4) {

	# <p> initialize
	.Object = callNextMethod(.Object, f0, f1, data, Nsi, model,
		start, objectiveFunction, sim, NmaxIS, weightingFunction, sampleFromAlt, lp, fudge);
	return(.Object);
});
setMethod('plausibilityStart', 'plausibilityGlmWeightedGaussian', function(this)this@mdl0$par);
setMethod('plausibilityDelta', 'plausibilityGlmWeightedGaussian', function(this)this@mdl0$sds * 2)

setMethod('plausibility', 'plausibilityGlmWeightedGaussian',
	function(this, start = NULL, optMethod = NULL, ...) {
	if (is.null(start)) start = plausibilityStart(this);
	P = this@objectiveFunction(start, this);
    rPl = new('plausibilityResult', par = start, value = min(P, 1), valueRaw = P, optimizer = list());
	return(rPl);
});



#
#	<p> penalized weighted plausbility
#

setClass("plausibilityPenalizedGaussian", contains = 'plausibilityPenalized', representation = list(
));

setMethod('initialize', 'plausibilityPenalizedGaussian', function(.Object, f0, f1, data, X, Nsi = 1e3, model,
	start = NULL, objectiveFunction = cumProbSIcomp, sim = NULL, NmaxIS = 5,
	weightingFunction = weightingFunctionLRpenalized, sampleFromAlt = TRUE, lp = NULL, fudge = NULL,
	NlambdaSel = 1e1, Nfolds = 10, standardize = TRUE, alpha = 0) {
	# <p> initialize
	.Object = callNextMethod(.Object, f0, f1, data, X, Nsi, model,
		start, objectiveFunction, sim, NmaxIS, weightingFunction, sampleFromAlt, lp, fudge,
		NlambdaSel, Nfolds, standardize, alpha
	);
	return(.Object);
});

setMethod('plausibility', 'plausibilityPenalizedGaussian',
	function(this, start = NULL, optMethod = NULL, ...) {
	if (is.null(start)) start = plausibilityStart(this);
	P = this@objectiveFunction(start, this);
    rPl = new('plausibilityResult', par = start, value = min(P, 1), valueRaw = P, optimizer = list());
	return(rPl);
});
