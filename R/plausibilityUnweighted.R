#
#	plausibilityUnweighted.R
#Thu Jun 17 17:12:03 CEST 2021

#
#	<p> plausibility for regression models
#

#
#	<p> glm model helper classes/functions
#


# bounder@binomial:
#	function(lp, epsPs = 1e-5) { lp[lp == 1] = 1 - epsPs; return(lp); }

cumProbSI = function(par, this) {
	N = length(this@y);
	lp = (this@mm %*% par[1:ncol(this@mm)])[, 1];
	# help optimizer, if necessary (e.g. bound away form 0, 1), defaults to identity
	lp = plausBounder(this@model, lp);
	parAncil = plausAncil(this@model, this@y, lp)
	TsRaw = apply(this@sim, 2, function(y)plausDensity(this@model, y, lp, parAncil, par));
	# importance sampling, cancel out uniform factor (max)
	TsIS0 = TsRaw - this@pSI;
	TsIS = TsIS0 - max(TsIS0);

	# <p> Test statistic + IS correction: reweight events
	Tdata = plausDensity(this@model, this@y, lp, parAncil, par);
	#<A> comparison based on ll instead of -ll, keep on this scale for IS correction
	sel = TsRaw <= Tdata;
	Ts = TsIS[sel];
	#Ts = minimax(Ts0, -Inf, log(this@NmaxIS));	# limit influence of single observations
	Nsim = sumExp(TsIS);	#Nsim = ncol(this@sim);

	# <p> P-value
	eps = Nsim / ncol(this@sim);
	eps = 0;
	#if (length(Ts) == 0) browser(); #return(1);	#<A> should not happen
	P = if (any(Ts == -Inf) || length(Ts) == 0) eps / 2 / (Nsim + eps) else ((sumExp(Ts) + eps) / (Nsim + eps));
	#if (any(Ts0 > log(this@NmaxIS))) print(list(Ptg = P)) else print(list(P = P))
	#if (P > 1.001) { print(list(current = par)); }
	return(P);
}

setClass("plausibilityGlm", contains = 'plausibilityFamilySI', representation = list(
	pSI = 'numeric',
	model = 'plausibilityModel',
	mdl0 = 'list',		# glm model + corresponing ll
	NmaxIS = 'numeric',	# maximal upweighting of observations during IS
	sim = 'matrix'
), prototype = list());

setMethod('initialize', 'plausibilityGlm', function(.Object,
	f0, data, Nsi = 1e3L,
	model, start = NULL, objectiveFunction = cumProbSI, sim = NULL, NmaxIS = 5L, lp = NULL) {

	# <p> initialize
	.Object = callNextMethod(.Object, f0, data, objectiveFunction, Nsi);
	#modelNm = Sprintf('plausibilityModel%{family}s');
    .Object@model = model;
	.Object@NmaxIS = NmaxIS;
	.Object@mdl0 = plausFit(model, .Object@y, .Object@mm);

	# <p> linear predictor
	if (is.null(lp)) lp = (.Object@mm %*% .Object@mdl0$par[1:ncol(.Object@mm)])[, 1];
	# <p> ancillary parameters
	parAncil = plausAncil(model, .Object@y, lp);
	# <p> stochastic integration sample	# plausSample(.Object, ., lp) %.% .
	if (is.null(sim))
		sim = apply(matrix(runif(length(lp) * Nsi), ncol = Nsi), 2, function(u)
			plausSample(model, u, lp, parAncil, .Object@mdl0$par));
	.Object@sim = sim;

	# <p> probabilities stochastic integration sample under starting pars, needed for importance correction
	.Object@pSI = apply(sim, 2, function(r)plausDensity(model, r, lp, parAncil, .Object@mdl0$par));
	return(.Object);
});

setMethod('plausibilityStart', 'plausibilityGlm', function(this)this@mdl0$par);
setMethod('plausibilityDelta', 'plausibilityGlm', function(this)this@mdl0$sds * 6)
#setMethod("plausFit", 'plausibilityGlm', function(this, y, mm)glmModel(mm, y, this))


