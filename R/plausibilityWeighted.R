#
#	plausbilityWeighted.R
#Thu Jun 17 17:11:33 CEST 2021

require('glmnet');

#
#	<p> weighted models
#

log1pExp = function(x, scale = 1) ifelse (x * scale < -5, x, log(log1p(exp(x * scale))) / scale) 

cumProbSIcomp = function(par, this) {
	N = length(this@y);
	Nsim = this@Nsi;
	lp = (this@mm %*% par)[, 1];
	# help optimizer, if necessary (e.g. bound away form 0, 1), defaults to identity
	lp = plausBounder(this@model, lp);
	parAncil = plausAncil(this@model, this@y, lp)
	TsRaw = apply(this@sim, 2, function(y)plausDensity(this@model, y, lp, parAncil));
	# importance sampling correction
	TsIS0 = TsRaw - this@pSI;
	Ts = log1pExp( TsIS0, .5 );	# dampen weighing
	#Ts = TsIS0;

	#<A> comparison based on ll instead of -ll, keep on thies scale for IS correction
	#sel = this@weight > this@weights;       # this is constant <N>
	#Ts = TsIS[sel];
	#Nsim = sumExp(TsIS);    #Nsim = ncol(this@sim);
	#if (length(Nsim) == 0) return(1e50 * .Machine$double.xmin);

	# <p> P-value
	eps = 1e-3 / Nsim;
	#eps = 0;
	if (any(Ts == -Inf) || length(Ts) == 0) return( eps / (Nsim + eps) );
	P = ((sumExp(Ts) + eps) / (Nsim + eps));
	#dprint(par, P);
	#print(list(par, P));
	#if (P > 1) browser();
	return(P);
}

#
#	<p> weighted GLM models
#

weightingFunctionLR = function(model, y, mm0, mm1) {
	m0 = plausFit(model, y, mm0);
	m1 = plausFit(model, y, mm1);
	# ordering in favor of alternative has to be small values first
	#return(-2*( m1$ll - m0$ll ));
	# <!> 28.6.2021 changed to integrate the tail to achieve more stable P-values,
	#	dispense with factor 2
	return( m1$ll - m0$ll );
}

fudgeLp = function(model, mm1, this, scale = 100, fudge = NULL) {
	lp0 = (this@mm %*% this@mdl0$par)[, 1];
	lp1 = (mm1 %*% plausFit(model, this@y, mm1)$par)[, 1];
	if (is.null(fudge)) fudge = dim(mm1)[2] / dim(mm1)[1] * scale;
	return( lp0 + (lp1 - lp0) / fudge);
}

setClass("plausibilityGlmWeighted", contains = 'plausibilityFamilyWeightedSI', representation = list(
	pSI = 'numeric',
	model = 'plausibilityModel',
	mdl0 = 'list',		# glm model + corresponing ll
	sim = 'matrix'
	#weights = 'numeric',	# simulation weights
	#weight = 'numeric'		# data weight
), prototype = list());

setMethod('initialize', 'plausibilityGlmWeighted', function(.Object, f0, f1, data, Nsi = 1e3, model,
	start = NULL, objectiveFunction = cumProbSIcomp, sim = NULL, NmaxIS = 5,
	weightingFunction = weightingFunctionLR, sampleFromAlt = TRUE, lp = NULL, fudge = 8) {

	# <p> initialize
	.Object = callNextMethod(.Object, f0, f1, data, Nsi,
		objectiveFunction = objectiveFunction,
		weightingFunction = weightingFunction, NmaxIS = NmaxIS);

	#modelNm = Sprintf('plausibilityModel%{family}s');
    .Object@model = model;
	.Object@mdl0 = plausFit(model, .Object@y, .Object@mm);

	# <p> sample from null/alternative
	mm1 = model_matrix_from_formula(.Object@f1, data)$mm;
	#lp = (mm1 %*% (plausFit(model, y, mm1)$par))[, 1];
	lp = if (is.null(sim) && is.null(lp) && sampleFromAlt)
		fudgeLp(model, mm1, .Object, fudge = fudge) else
		(.Object@mm %*% .Object@mdl0$par)[, 1];

	# <p> ancillary parameters
	parAncil = plausAncil(model, .Object@y, lp);
	# <p> stochastic integration sample	# plausSample(.Object, ., lp) %.% .
	if (is.null(sim))
		sim = apply(matrix(runif(length(lp) * Nsi), ncol = Nsi), 2, function(u)
			plausSample(model, u, lp, parAncil));

	# <p> simulation weights
	weights = apply(sim, 2, function(r)weightingFunction(model, r, .Object@mm, mm1));
	# <p> data weights
	weight = weightingFunction(model, .Object@y, .Object@mm, mm1);

	# filter simulation for required sub-sample
	.Object@sim = sim[, weights > weight ];

	# <p> probabilities stochastic integration sample under starting pars, needed for importance correction
	.Object@pSI = apply(.Object@sim, 2, function(r)plausDensity(model, r, lp, parAncil));
	#print(list(N = sum(weights > weight), ll = .Object@pSI));
	return(.Object);
});
setMethod('plausibilityStart', 'plausibilityGlmWeighted', function(this)this@mdl0$par);
setMethod('plausibilityDelta', 'plausibilityGlmWeighted', function(this)this@mdl0$sds * 2)
