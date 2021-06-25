#
#	plausbilityWeighted.R
#Thu Jun 17 17:11:33 CEST 2021

require('glmnet');

#
#	<p> weighted models
#

cumProbSIcomp = function(par, this) {
	N = length(this@y);
	Nsim = ncol(this@sim);
	lp = (this@mm %*% par)[, 1];
	# help optimizer, if necessary (e.g. bound away form 0, 1), defaults to identity
	lp = plausBounder(this@model, lp);
	parAncil = plausAncil(this@model, this@y, lp)
	TsRaw = apply(this@sim, 2, function(y)plausDensity(this@model, y, lp, parAncil));
	# importance sampling, cancel out uniform factor (max)
	TsIS0 = TsRaw - this@pSI;
	#TsIS = TsIS0 - max(TsIS0);
	TsIS = TsIS0;

	#<A> comparison based on ll instead of -ll, keep on thies scale for IS correction
	sel = this@weight > this@weights;       # this is constant <N>
	Ts = TsIS[sel];
	#Nsim = sumExp(TsIS);    #Nsim = ncol(this@sim);
	#if (length(Nsim) == 0) return(1e50 * .Machine$double.xmin);

	# <p> P-value
	eps = .1 / ncol(this@sim);
	#eps = 0;
	P = if (any(Ts == -Inf) || length(Ts) == 0)
			(eps / 2) / (Nsim + eps) else
			((sumExp(Ts) + eps) / (Nsim + eps));
	#if (P > 1) browser();
	#dprint(par, P);
	return(P);
}

#
#	<p> weighted GLM models
#

weightingFunctionLR = function(model, y, mm0, mm1) {
	m0 = plausFit(model, y, mm0);
	m1 = plausFit(model, y, mm1);
	# ordering in favor of alternative has to be small values first
	return(-2*( m1$ll - m0$ll));
}

setClass("plausibilityGlmWeighted", contains = 'plausibilityFamilyWeightedSI', representation = list(
	pSI = 'numeric',
	model = 'plausibilityModel',
	mdl0 = 'list',		# glm model + corresponing ll
	sim = 'matrix',
	weights = 'numeric',	# simulation weights
	weight = 'numeric'		# data weight
), prototype = list());

setMethod('initialize', 'plausibilityGlmWeighted', function(.Object, f0, f1, data, Nsi = 1e3, model,
	start = NULL, objectiveFunction = cumProbSIcomp, sim = NULL, NmaxIS = 5,
	weightingFunction = weightingFunctionLR, sampleFromAlt = TRUE, lp = NULL) {

	mm1 = model_matrix_from_formula(f1, data)$mm;
	if (is.null(sim) && is.null(lp) && sampleFromAlt) {
		y = data[, formula.response(f1)];
		lp = (mm1 %*% plausFit(model, y, mm1)$par)[, 1];
	}
	# <p> initialize
	.Object = callNextMethod(.Object, f0, f1, data, Nsi,
		objectiveFunction = objectiveFunction,
		weightingFunction = weightingFunction, NmaxIS = NmaxIS);

	#modelNm = Sprintf('plausibilityModel%{family}s');
    .Object@model = model;
	.Object@mdl0 = plausFit(model, .Object@y, .Object@mm);

	# <p> linear predictor
	if (is.null(lp)) lp = (.Object@mm %*% .Object@mdl0$par)[, 1];
	# <p> ancillary parameters
	parAncil = plausAncil(model, .Object@y, lp);
	# <p> stochastic integration sample	# plausSample(.Object, ., lp) %.% .
	if (is.null(sim))
		sim = apply(matrix(runif(length(lp) * Nsi), ncol = Nsi), 2, function(u)
			plausSample(model, u, lp, parAncil));
	.Object@sim = sim;

	# <p> probabilities stochastic integration sample under starting pars, needed for importance correction
	.Object@pSI = apply(sim, 2, function(r)plausDensity(model, r, lp, parAncil));
	# <p> simulation weights
	.Object@weights = apply(sim, 2, function(r)weightingFunction(model, r, .Object@mm, mm1));
	# <p> data weights
	.Object@weight = weightingFunction(model, .Object@y, .Object@mm, mm1);

	return(.Object);
});
setMethod('plausibilityStart', 'plausibilityGlmWeighted', function(this)this@mdl0$par);
setMethod('plausibilityDelta', 'plausibilityGlmWeighted', function(this)this@mdl0$sds * 2)
