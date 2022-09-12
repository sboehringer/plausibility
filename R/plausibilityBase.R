#
#	plausibilityBase.R
#Fri 12 Nov 2021 03:41:45 PM CET
# moved from plausibility dt to avoid references

# Class Structure
#	plausibilityFamily [base class]
#		plausibilityFamilySI
#			plausibilityGlm [glm models: handles mm, but not specific LH]
#		plausibilityFamilyWeighted [weighted analsyes]
#			plausibilityFamilyWeightedSI
#				plausibilityGlmWeighted [weighted analsyes, glm models]
#
#
#	plaubilityModelGlm [base class for GLM models]
#		plaubilityModelBinomial


#
#   <p> S4 classes
#

#
#	<p> generic unweighted plausibility classes
#

dummyObjective = function(par, this)-Inf;
regionObjective = function(objective, level = .95) {
	return(function(par, this)(regionObjective(par, this) - (1 - level))^2);
}

setClass("plausibilityFamily", representation = list(
	f0 = 'formula',
	mm = 'matrix',
	y = 'numeric',
	objectiveFunction = 'function'
), prototype = list(
	f0 = as.formula(NULL),
	mm = matrix(),
	objectiveFunction = dummyObjective
));
setGeneric("plausibility", function(this, optMethod = 'grid', Niter = 1L, ...)this)
setGeneric("plausibilityStart", function(this)this)
setGeneric("plausibilityDelta", function(this)NULL)
setGeneric("region", function(this, level = .95, ...)this)

setMethod('initialize', 'plausibilityFamily', function(.Object, f0, data,
	objectiveFunction = dummyObjective) {
	Data = DfNames2std(data, vars.as.rhs(formula.response(f0)), ~ y);
	f0 = formula.add.response(formula.expand(f0, data), y ~ 1); # expand ~ ., std outcome name
	.Object@f0 = f0;
	.Object@mm = model_matrix_from_formula(f0, Data)$mm;
	.Object@y = Data[, 'y'];
	.Object@objectiveFunction = objectiveFunction;
	return(.Object);
});

#
#	<p> stochastic integration
#

setClass("plausibilityFamilySI", contains = 'plausibilityFamily', representation = list(
	Nsi = 'integer'
), prototype = list(
	Nsi = 1e3L
));
setMethod('initialize', 'plausibilityFamilySI', function(.Object, f0, data, objectiveFunction, Nsi) {
	.Object = callNextMethod(.Object, f0, data, objectiveFunction);
	.Object@Nsi = Nsi;
	return(.Object);
});


#
#	<p> result class
#

setClass("plausibilityResult", representation = list(
	par = 'numeric',
	value = 'numeric',
	valueRaw = 'numeric',
	optimizer = 'list',
	object = 'plausibilityFamily'
	#object = 'plausibilityFamily'	 # forces empty initialization
), prototype = list(
	par = as.numeric(NA),
	value = as.numeric(NA),
	valueRaw = as.numeric(NA),
	optimizer = list()
));

#
#	<p> generic weighted plausibility classes
#

setClass("plausibilityFamilyWeighted", contains = 'plausibilityFamily', representation = list(
	f1 = 'formula',
	weightingFunction = 'function'
), prototype = list(
	weightingFunction = identity
));

# objectiveFunction can be cumProbSIcomp in most cases, dummyObjective used here to avoid
#	forward reference (see plausibilityWeighted.R)
setMethod('initialize', 'plausibilityFamilyWeighted', function(.Object, f0, f1, data,
	objectiveFunction = dummyObjective, weightingFunction = identity) {

	.Object = callNextMethod(.Object, f0, data, objectiveFunction);
	.Object@f1 = formula.expand(f1, data);
	.Object@weightingFunction = weightingFunction;
	return(.Object);
});

setClass("plausibilityFamilyWeightedSI", contains = 'plausibilityFamilyWeighted', representation = list(
	Nsi = 'integer',
	NmaxIS = 'numeric'	# maximal upweighting of observations during IS
), prototype = list(
	weightingFunction = identity
));

setMethod('initialize', 'plausibilityFamilyWeightedSI', function(.Object, f0, f1, data, Nsi = 1e3L,
	objectiveFunction = dummyObjective, weightingFunction = identity, NmaxIS = 5) {

	.Object = callNextMethod(.Object, f0, f1, data, objectiveFunction, weightingFunction);
	.Object@Nsi = Nsi;
	.Object@NmaxIS = NmaxIS;
	return(.Object);
});


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
	if (is.null(start)) start = plausibilityStart(this);
	#print(list(start = start));
	rPl = optimizer(this, this@objectiveFunction, start, ...);
	return(rPl);
});

#
#	<p> internal interface
#

PlausibilityUnweighted = function(f0, data, family,
	Nsi = 1e3L, Niter = 1L, optMethod, ...) {
	modelNm = Sprintf('plausibilityModel%{family}u');
	model = new(modelNm, family = family, ...);
	Data = completeData(f0, data);
	start = NULL;
	for (i in Seq(1, Niter)) {
		pl = new('plausibilityGlm', f0 = f0, data = Data, Nsi = Nsi, model = model);
		rPl = plausibility(pl, start = start, optMethod = optMethod);
		start = rPl@par;
	}
	return(rPl);
}

plClasses = list(
	binomial = 'plausibilityGlmWeighted',
	gaussian = 'plausibilityGlmWeightedGaussian'
);
plFudgFactors = list(
	binomial = 6,
	gaussian = 4,
	NIbinomial = 1
);

PlausibilityWeighted = function(f0, f1 = NULL, data, family, 
	Nsi = 1e3L, Niter = 2L, optMethod, fudge = NULL, ..., plClass, initArgs = list()) {
	if (missing(plClass)) plClass = plClasses[[family]];
	if (is.null(fudge)) fudge = plFudgFactors[[family]];
	modelNm = Sprintf('plausibilityModel%{family}u');
	model = new(modelNm, family = family, ...);
	Data = completeData(if (class(f1) != 'formula') f0 else f1, data);

	start = NULL;
	for (i in Seq(1, Niter)) {
		args = merge.lists(
			list(f0 = f0, f1 = f1, data = Data, Nsi = Nsi, model = model, fudge = fudge),
			initArgs);
		pl = do.call('new', c(list(Class = plClass), args));
		rPl = plausibility(pl, start = start, optMethod = optMethod);
		start = rPl@par;
	}
	return(rPl);
}

#
#	<p> plausibility region
#

setMethod('region', 'plausibilityFamily',
	function(this, start = NULL, level = .95, Napprox = 20, ..., Nretries = 2, sigmaScale = 3, simplify = TRUE) {

	# search aroudn 3 SDs for a level of .95
	# <i> calibrate to level
	mn = if (is.null(start)) this@mdl0$par else start;
	delta = sigmaScale * this@mdl0$sds
	r = NA;
	for (i in 1:Nretries) {
		ranges = apply(cbind(mn - delta, mn + delta), 1, identity, simplify =- F);	# list fo ranges

		# 	cls = contourLinesStacked(ranges, Napprox, this@objectiveFunction);
		# 	contour = completeLevelSets(cls, contour = 1 - level);
		contour = try(FindRegion(this@objectiveFunction, ranges, level = level, Nout = Napprox, this = this));
		nr = lapply(contour, \(.)nrow(.@polygons[[1]]@Polygons[[1]]@coords));
		if (!any(unlist(sapply(contour, class)) == 'try-error') && all(nr > 1)) {
			r = contour;
			break;
		}
		delta = delta * sigmaScale;
	}
	if (simplify && length(level) == 1 && !is.na(r)) r = r[[1]];
	return(list(region = r));
});

PlausibilityRegionFull = function(f0, data, family = 'gaussian', level = .95, Nsi = 1e3L, Niter = 1L, Napprox = 20L, ..., optMethod = 'grid', plClass = 'plausibilityGlm', initArgs = list(), sigmaScale = 3) {
	modelNm = Sprintf('plausibilityModel%{family}u');
	model = new(modelNm, family = family, ...);
	Data = completeData(f0, data);

	args = merge.lists(
		list(f0 = f0, data = Data, Nsi = Nsi, model = model),
		initArgs);
	pl = do.call('new', c(list(Class = plClass), args));

	plRegion = region(pl, start = NULL, Napprox = Napprox, optMethod = optMethod, level = level, sigmaScale = sigmaScale);
	return(plRegion);
}

#
#	<p> external interface
#

#' Compute weighted or unweighted plausibility function
#'
#' A model comparison of binomial outcome for nested models is computed. The plausibility estimate and P-value are returned together with a glm-fit.
#'
#' @param f0 Null model given as R formula
#' @param f1 Alternative model given as R formula
#' @param data data frame
#' @param optMethod Optimization method to use. Either 'grid' (default) or 'optim'
#' @param Niter Number of iterations used when finding the plausibility estimate
#' @param Nsi Number of stochastic integration samples (defaults to 1e3L)
#' @param ... Arguments passed to the optimizer
#'
#' @export Plausibility
Plausibility = function(f0, f1 = NULL, data, family = 'binomial',
	Nsi = 1e3L, Niter = 2L, optMethod = 'optim', ..., plClass, initArgs = list()) {

	rPl = if (notE(f1))
		PlausibilityWeighted(f0, f1, data, family, Nsi, Niter, optMethod, ...,
			plClass = plClass, initArgs = initArgs) else
		PlausibilityUnweighted(f0, data, family, Nsi, Niter, optMethod, ...);
	return(rPl);
}


#' Compute full or marginal plausibility region
#'
#' glm-fit.
#'
#' @param f0 Null model given as R formula
#' @param f1 optional formula containing variables to be marginalized over
#' @param data data frame
#' @param optMethod Optimization method to use. Either 'grid' (default) or 'optim'
#' @param Niter Number of iterations used when finding the plausibility estimate
#' @param Nsi Number of stochastic integration samples (defaults to 1e3L)
#' @param ... Arguments passed to the optimizer
#'
#' @export Plausibility
PlausibilityRegion = function(f0, f1 = NULL, data, family = 'binomial', level = .95,
	Nsi = 1e3L, Niter = 1L, Napprox = 20L, optMethod = 'optim', ..., plClass = 'plausibilityGlm', initArgs = list(),
	sigmaScale = 3) {

	rPl = if (notE(f1))
		PlausibilityRegionMarginal(f0, f1, data, family, Nsi, Niter, Napprox, optMethod, ...,
			plClass = plClass, initArgs = initArgs, sigmaScale = sigmaScale) else
		PlausibilityRegionFull(f0, data, family, level, Nsi, Niter, Napprox, optMethod = optMethod, ...,
			plClass = plClass, initArgs = initArgs, sigmaScale = sigmaScale);
	return(rPl);
}
