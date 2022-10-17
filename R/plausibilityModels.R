#
#	plausibilityModels.R
#Thu Jun 17 17:03:40 CEST 2021

#
#	<p> helpers
#


glmLL = function(par, X, y, this, offset = 0) plausDensity(this, y,  (X %*% par)[, 1] + offset);
glmModelOptimize = function(X, y, this, offset) {
	start = rep(0, ncol(X));
	o = try(
		optim(start, glmLL, method = 'BFGS', control = list(fnscale = -1),
			X = X, y = y, this = this, offset = offset)
	, silent = T);
	#if (class(o) == 'try-error') browser();
	if (class(o) == 'try-error') return(list(par = rep(NA, ncol(X)), ll = -Inf));
	return(list(par = o$par, ll = o$value))
}

#
#	<p> generic model class
#

# <p> core methods
setGeneric("plausSample", function(this, u, lp, parAncil, par)NULL)
setGeneric("plausDensityS", function(this, x, lp, parAncil, par)-Inf)

plausFitModelGlm = function(this, X, y, offset) {
	r = glm(y ~ . + 0, as.data.frame(cbind(X, y)), family = this@family, offset = offset);
	return(list(par = r$coefficients, sds = sqrt(diag(vcov(r))), model = r));
}
setGeneric("plausFitModel", plausFitModelGlm);

plausFitGeneric = function(this, X, y, offset = 0) {
	m = plausFitModel(this, X, y, recycle(offset, y)[[1]]);
	parAncil = plausAncil(this, y, m$model$linear.predictors)
	ll = plausDensity(this, y, m$model$linear.predictors, parAncil, m$par);
	if (any(is.na(ll))) return(list(par = NA));
	if (ll == -Inf) return(glmModelOptimize(X, y, this, offset));
	return(list(par = m$par, ll = ll, sds = m$sds, model = m$model));
}
setGeneric("plausFit", function(this, y, mm, offset = 0)plausFitGeneric(this, mm, y, offset))

# <p> methods with "good" defaults
setGeneric("plausDensity", function(this, x, lp, parAncil, par)sum(plausDensityS(this, x, lp, parAncil, par)))
setGeneric("plausBounder", function(this, lp)lp)
setGeneric("plausLp", function(this, par, mm){ (mm %*% par[1:ncol(mm)])[, 1] })

# <p> estimate ancillory paramaters
setGeneric("plausAncil", function(this, y, lp, ...)NULL)

# methods alternative model, used for equivalence hyptotheses, default to non-Alt methods
setGeneric("plausFitAlt", function(this, y, mm, offset = 0)plausFit(this, y, mm, offset = 0))
setGeneric("plausLpAlt", function(this, par, mm)plausLp(this, par, mm))

setClass('plausibilityModel', representation = list(
	family = 'character'
), prototype = list(family = 'gaussian'));

setMethod('initialize', 'plausibilityModel', function(.Object, family = 'gaussian') {
	.Object@family = family;
	return(.Object);
});
setMethod("plausBounder", 'plausibilityModel', function(this, lp)lp)


#
#	<p> glm models
#


#
#	<p> normal
#

setClass('plausibilityModelGaussian', contains = 'plausibilityModel', representation = list());
setMethod("plausSample", 'plausibilityModelGaussian', function(this, u, lp, parAncil, par) {
	qnorm(u, mean = lp, sd = parAncil)
})
setMethod("plausDensityS", 'plausibilityModelGaussian', function(this, x, lp, parAncil, par)
	dnorm(x, lp, sd = parAncil, log = T))
setMethod('plausAncil', 'plausibilityModelGaussian', function(this, y, lp) {
	return(sd(y - lp))
})
setMethod('initialize', 'plausibilityModelGaussian', function(.Object, family = 'gaussian')
	return(callNextMethod(.Object, family))
);

#
#	<p> binomial
#

dataBinomialSplit = function(y, data, Nbinom = 2) {
	d1 = rep.each.row(data, Nbinom);
	yBernoulli = sapply(data[, y], function(y)c(rep(0, Nbinom - y), rep(1, y)));
	d1[, y] = as.vector(yBernoulli);
	return(d1);
}

setClass('plausibilityModelBinomial', contains = 'plausibilityModel', representation = list(
	Nbinom = 'integer',	# #{binomial outcomes} == 1 for standard logistic regression
	eps = 'numeric'
), prototype = list(Nbinom = 1L, eps = 1e-4));


setMethod('initialize', 'plausibilityModelBinomial', function(.Object, family, Nbinom = 1L) {
	.Object = callNextMethod(.Object, family);
	.Object@Nbinom = as.integer(Nbinom);
	return(.Object);
});
setMethod("plausBounder", 'plausibilityModelBinomial', function(this, lp) {
	ps = plogis(lp);
	lp[ps == 1] = qlogis(1 - this@eps);	# help optimizer
	return(lp);
})

setMethod("plausSample", 'plausibilityModelBinomial', function(this, u, lp, parAncil, par)
	{ qbinom(u, this@Nbinom, plogis(lp)) })
setMethod("plausDensityS", 'plausibilityModelBinomial', function(this, x, lp, parAncil, par) {
	dbinom(x, this@Nbinom, plogis(lp), log = T)
})
setMethod("plausFit", 'plausibilityModelBinomial', function(this, y, mm, offset) {
	if (this@Nbinom > 1L) {
#if (length(y) != dim(mm)[1]) browser();
		dSp = dataBinomialSplit('y', cbind(y, mm), this@Nbinom);
		y = dSp[, 1];
		mm = dSp[, -1, drop = FALSE];
		if (!missing(offset)) offset = rep.each(offset, this@Nbinom);
	}
	callNextMethod(this, y, mm, offset);
});

#
#	<p> negative binomial
#

setClass('plausibilityModelNegativeBinomial', contains = 'plausibilityModel', representation = list(), prototype = list());

setMethod('initialize', 'plausibilityModelNegativeBinomial', function(.Object, family = 'negativeBinomial') {
	.Object = callNextMethod(.Object, family);
	return(.Object);
});


# <!> offset needs to be implemented
glmFitNb = function(this, X, y, offset) {
	#glm.nb(y ~ . + 0, as.data.frame(cbind(X, y)), offset = offset);
	r = glm.nb(y ~ . + 0, as.data.frame(cbind(X, y)));
	par = c(r$coefficients, log(r$theta));
	sds = c(sqrt(diag(vcov(r))), r$SE.theta);
	return(list(par = par, sds = sds, model = r));
}
setMethod("plausFitModel", 'plausibilityModelNegativeBinomial', glmFitNb);

s2fromMuSize = function(mu, size)(mu/size + 1)*mu
probFromMuS2 = function(mu, s2)(1 - mu/s2)
probFromMuSize = function(mu, size)probFromMuS2(mu, s2fromMuSize(mu, size))

# theta/size coded as first parameter
setMethod("plausSample", 'plausibilityModelNegativeBinomial', function(this, u, lp, parAncil, par) {
	size = last(par);
	#qnbinom(u, size, probFromMuSize(exp(lp), size))
	qnbinom(u, exp(last(par)), mu = exp(lp))
})
setMethod("plausDensityS", 'plausibilityModelNegativeBinomial', function(this, x, lp, parAncil, par) {
	size = last(par);
	#dnbinom(x, size, probFromMuSize(exp(lp), size), log = TRUE)
	dnbinom(x, exp(last(par)), mu = exp(lp), log = TRUE)
})
