#
#	plausibilityModels.R
#Thu Jun 17 17:03:40 CEST 2021

#
#	<p> generic model class
#

#
#	<p> glm models
#

setGeneric("plausSample", function(this, u, lp, parAncil)NULL)
setGeneric("plausDensityS", function(this, x, lp, parAncil)-Inf)
setGeneric("plausDensity", function(this, x, lp, parAncil)sum(plausDensityS(this, x, lp, parAncil)))
setGeneric("plausFit", function(this, y, mm, offset = 0)glmModel(mm, y, this, offset))
setGeneric("plausBounder", function(this, lp)lp)
# estimate ancillory paramaters
setGeneric("plausAncil", function(this, y, lp, ...)NULL)

setClass('plausibilityModel', representation = list(
	family = 'character'
), prototype = list(family = 'gaussian'));

setMethod('initialize', 'plausibilityModel', function(.Object, family = 'gaussian') {
	.Object@family = family;
	return(.Object);
});
setMethod("plausBounder", 'plausibilityModel', function(this, lp)lp)


#
#	<p> normal
#

setClass('plausibilityModelGaussian', contains = 'plausibilityModel', representation = list());
setMethod("plausSample", 'plausibilityModelGaussian', function(this, u, lp, parAncil) {
	qnorm(u, mean = lp, sd = parAncil)
})
setMethod("plausDensityS", 'plausibilityModelGaussian', function(this, x, lp, parAncil)
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

setMethod("plausSample", 'plausibilityModelBinomial', function(this, u, lp, parAncil)
	{ qbinom(u, this@Nbinom, plogis(lp)) })
setMethod("plausDensityS", 'plausibilityModelBinomial', function(this, x, lp, parAncil) {
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
