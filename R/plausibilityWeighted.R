#
#	plausbilityWeighted.R
#Thu Jun 17 17:11:33 CEST 2021

#
#	<p> weighted models
#

log1pExp = function(x, scale = 1) ifelse (x * scale < -5, x, log(log1p(exp(x * scale))) / scale) 
logAtOff = function(x, off = 3) { ifelse (x <= off, x, log1p(ifelse(x < off, off, x) - off) + off)  }
dampen = log1pExp;

cumProbSIcomp = function(par, this, dampen = logAtOff) {
	N = length(this@y);
	Nsim = this@Nsi;
	#lp = (this@mm %*% par)[, 1];
	lp = plausLp(this@model, par, this@mm);
	# help optimizer, if necessary (e.g. bound away form 0, 1), defaults to identity
	lp = plausBounder(this@model, lp);
	parAncil = plausAncil(this@model, this@y, lp)
	TsRaw = apply(this@sim, 2, function(y)plausDensity(this@model, y, lp, parAncil, par));
	# importance sampling correction
	TsIS0 = TsRaw - this@pSI;
	Ts = dampen( TsIS0 );	# dampen weighing
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
	#if (P > 1.2) browser();
# 	if (P < 1e-6 && abs(par - this@mdl0$par) < .05) {
# 		lr = NILR(df2LRmatrix(Df(y = this@y, x = this@mm[, 2])), this@model@fnull);
# 		parLr = c(lr$o[[2]]$par[1], niH02h1(lr$o[[2]]$par, this@model@fnull));
# 		print(c(par = par, parLr = parLr, P = P, T = exp(Ts[1:5]), lp = c(lp[1:5], lp[201:205])));
# 		#browser();
# 	}
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
	lr =  m1$ll - m0$ll;
	return( max(lr, 0) );
}

fudgeLp = function(model, mm1, this, scale = 100, fudge = NULL) {
	#lp0 = (this@mm %*% this@mdl0$par)[, 1];
	lp0 = plausLp(model, this@mdl0$par, this@mm);
	#lp1 = (mm1 %*% plausFit(model, this@y, mm1)$par)[, 1];
	lp1 = plausLpAlt(model, plausFitAlt(model, this@y, mm1)$par, mm1);
	if (is.null(fudge)) fudge = dim(mm1)[2] / dim(mm1)[1] * scale;
	if (fudge < 1) warning(Sprintf('Fudge factor %{fudge}f < 1'));
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
		#(.Object@mm %*% .Object@mdl0$par)[, 1];
		plausLp(this, .Object@mdl0$par, .Object@mm);

	# <p> ancillary parameters
	parAncil = plausAncil(model, .Object@y, lp);
	# <p> stochastic integration sample	# plausSample(.Object, ., lp) %.% .
	if (is.null(sim))
		sim = apply(matrix(runif(length(lp) * Nsi), ncol = Nsi), 2, function(u)
			# <i> sampleFromAt && negBinomial
			plausSample(model, u, lp, parAncil, .Object@mdl0$par));

	# <p> simulation weights
	#weights = apply(sim, 2, function(r)weightingFunction(model, r, .Object@mm, mm1));
	weights = unlist(apply(sim, 2, function(r)weightingFunction(model, r, .Object@mm, mm1)));
	# <p> data weights
	weight = weightingFunction(model, .Object@y, .Object@mm, mm1);

	# filter simulation for required sub-sample
	#	assume weights tb LRs -> large values in rejection region
	#	">=" instead of ">" necessary for point masses
	.Object@sim = sim[, weights >= weight, drop = FALSE ];
	#dprint(sum(weights > weight));

	if (sum(weights >= weight) == 0) warning('Plausib: no samples generated more extreme than data');
	#if (sum(weights >= weight) == 0) browser();
	#lr = NILR(df2LRmatrix(Df(y = .Object@y, x = .Object@mm[, 2])), .Object@model@fnull);
	#parLr = c(lr$o[[2]]$par[1], niH02h1(lr$o[[2]]$par, .Object@model@fnull));
	#if (sum(weights >= weight) == Nsi && lr$o[[2]]$par[2] > 0) browser();

	# <p> probabilities stochastic integration sample under starting pars, needed for importance correction
	.Object@pSI = apply(.Object@sim, 2, function(r)plausDensity(model, r, lp, parAncil, .Object@mdl0$par));
	#print(df2LRmatrix(Df(y = .Object@y, x = .Object@mm[, -1])));
	#print(c(Nweights = sum(weights >= weight), weight = weight, par = .Object@mdl0$par, parAlt = plausFitAlt(model, .Object@y, mm1)$par), lp = c(lp[1:5], lp[201:205]));
	#print(list(N = sum(weights > weight), ll = .Object@pSI));
	return(.Object);
});
setMethod('plausibilityStart', 'plausibilityGlmWeighted', function(this)this@mdl0$par);
setMethod('plausibilityDelta', 'plausibilityGlmWeighted', function(this)this@mdl0$sds * 2)
