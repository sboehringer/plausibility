#
#	Rstatistic.R
#Fri 19 Jan 2007 11:06:44 PM CET 

# contains simple statistics to evaluate consulting questions

sizesDesc = function(s) {
	col.frame(list(
		mean = mean(s),
		median = median(s),
		stddev = sqrt(var(s)),
		quantiles = quantile(s)
	), do.paste = " ", digits = 1)
}

compareSamples = function(l) {
	desc = data.frame(lapply(l, function(e)sizesDesc(e)));
	print(desc);

	tests = col.frame(list(
		test.t = t.test(l[[1]], l[[2]])$p.value,
		test.wilcoxon = wilcox.test(l[[1]], l[[2]])$p.value
	));
	print(tests);
}

df2numeric = function(df) apply(df, 2, function(col)as.numeric(as.vector(col)));
expandCounts = function(tab)  unlist(apply(tab, 1, function(r){rep(r[1], r[2])}));

chisq.test.robust = function(tab, bootstrapCellCount = 5, B = 5e3) {
	# reduce table by removing 0-marginals and check for degeneration
	tab = tab[, !apply(tab, 2, function(c)all(c == 0))];
	if (is.vector(tab)) return(list(p.value = NA));
	tab = tab[!apply(tab, 1, function(r)all(r == 0)), ];
	if (is.vector(tab)) return(list(p.value = NA));
	# determine whether to bootstrap
	r = if (any(tab < bootstrapCellCount))
		chisq.test(tab, simulate.p.value = T, B = B) else
		chisq.test(tab);
	r
}

# depends on coin package <!>, unfinished
armitage.test.robust = function(formula, df, scores) {
	tab = table(df);
	# only eliminate 0-rows of table from score vector
	zRows = sapply(1:dim(tab)[1], function(i){ all(tab[i,] == 0) });
	scores[[1]] = scores[[1]][!zRows];
	r =	independence_test(formula, df, teststat = "quad", scores = scores);
	r
}

# simulations in project 2014-02-Microsatellites
logSumExpRaw = function(v, pivot = median(v))(log(sum(exp(v - pivot))) + pivot)
logSumExpPivot = logSumExpMax = function(v)logSumExpRaw(v, pivot = max(v))
logSumExp = function(x) {
	Imx = which.max(x);
	log1p(sum(exp(x[-Imx] - x[Imx]))) + x[Imx]
}
logSumExp_add = function(v, bins) {
	diff = abs(bins - v);
	diffB = pop(bins) - shift(bins);
	if (max(diff) > max(diffB)) {	# bin-reshuffeling
		Ibin = which.min(diffB);
		bins[Ibin] = logSumExp(bins[c(Ibin, Ibin + 1)]);
		bins[Ibin + 1] = v;
	} else {
		Ibin = which.min(abs(diff));
		bins[Ibin] = logSumExp(c(v, bins[Ibin]));
	}
	return(sort(bins, na.last = T, decreasing = T));
}

logSumExpIter_raw = function(x, Nbins = 30) {
	if (length(x) < Nbins) return(x);
	lse = x[1:Nbins];
	for (e in x[-(1:Nbins)]) lse = logSumExp_add(e, lse);
	return(logSumExp(lse));
}
logSumExpIter = function(x, Nbins = 30) {
	return(logSumExp(logSumExpIter_raw(x)));
}
# sumExp = function(x) {
# 	Imx = which.max(x);
# 	sum(exp(x - x[Imx])) * exp(x[Imx])
# }
sumExp = function(x)exp(logSumExp(x))
meanExp = function(x)(sumExp(x)/length(x))


# rejFrac = function(x, alpha = 0.05) {
# 	x = na.omit(x);
# 	f = count(x <= alpha) / length(x);
# 	f
# }
rejFrac = function(x, alpha = 0.05)mean(x <= alpha, na.rm = T);
vector.std = function(v, C = 1)(C * v / sum(v));
vector.std.log = function(v, C = 0)(v - (logSumExp(v) - C));

#
#	<p> ml methods
#

lhWrapperFunctions = c("initialize",
	"parsScale", "parsMap", "parsMapInv", "parsStart", "parsNames", "lh", "null2alt", "alt2null"
);

# <!> transition to S4-objects
lhGetWrapper = function(prefix, self, ...) {
	createNullWrapper = F;
	f = list();
	if (substr(prefix, nchar(prefix) - 3, nchar(prefix)) == "null") {
		createNullWrapper = T;
		prefix = substr(prefix, 1, nchar(prefix) - 5);
	}
	for (n in lhWrapperFunctions) {
		f[[n]] = mget(sprintf("%s.%s", prefix, n), envir = globalenv(), ifnotfound=list(NULL))[[1]];
	}
	f$self = if (is.null(self)) { if (is.null(f$initialize)) list(...) else f$initialize(...) } else self;
	if (createNullWrapper) {
		f1 = f;
		self = f1$self = f$self;
		f1$parsStart = function(self){ f$alt2null(self, f$parsStart(self)) };
		f1$parsScale = function(self){ f$alt2null(self, f$parsScale(self)) };
		f1$parsMap = function(self, p){ f$alt2null(self, f$parsMap(self, f$null2alt(self, p))) };
		f1$parsMapInv = function(self, p){ f$alt2null(self, f$parsMapInv(self, f$null2alt(self, p))) };
		f1$lh = function(self){ lhRaw = f$lh(self); function(p)lhRaw(f$null2alt(self, p)) };
		return(f1);
	}
	f
}
lhCopyWrapper = function(name, template) {
	for (f in lhWrapperFunctions) {
		g = mget(sprintf("%s.%s", template, f), envir = globalenv(), ifnotfound=list(NULL))[[1]];
		if (!is.null(g)) eval.parent(parse(text = sprintf("%s.%s = %s.%s;", name, f, template, f)));
	}
}

lhInit = function(lhWrapper) {
	
}

mapU = function(y){ -log(1/y - 1) }
map2U = function(z){ 1/(1 + exp(-z)) }
# one-dimensional estimation
lhMlEstimatorOD = function(lhWrapper = NULL, start = NULL, c = NULL, ...) {
	if (is.null(c)) c = list(tol = .Machine$double.eps^0.25);
	f = lhGetWrapper(lhWrapper, c$self, ...);

	lhRaw = f$lh(f$self);
	lh = function(p) { lhRaw(mapU(f$parsMap(f$self, p))) }
	o = try(optimize(lh, lower = 0, upper = 1, tol = c$tol, maximum = T));
	r = list(par = mapU(f$parsMap(f$self, o$maximum)), par.os = o$maximum, value = o$objective);
	r
}

# multi-dimensional estimation
lhMlEstimatorMD = function(lhWrapper = NULL, start = NULL, c = NULL, ...) {
	if (is.null(c)) c = list(do.sann = F, sann.cycles = 1000);
	f = lhGetWrapper(lhWrapper, c$self, ...);
	eps = 1e-5;
	#if (!is.null(start)) { starts = matrix(start, nrow = 1); }
	if (is.null(start)) start = f$parsStart(f$self);

	starts = if (!is.matrix(start)) matrix(as.numeric(unlist(start)), nrow = 1) else start;
	parscale = f$parsScale(f$self);
	lhRaw = f$lh(f$self);
	lh = function(p) { lhRaw(f$parsMap(f$self, p)) }
	os = apply(starts, 1, function(s) {
		s = f$parsMapInv(f$self, s);
		o = try(optim(s, lh, method = "Nelder-Mead",
			control = list(fnscale = -1, parscale = parscale, maxit = 1000),
		));
		if (class(o) == "try-error") return(NA);
		if (0) { # if (o$convergence > 0 || c$do.sann) {	# Nelder-Mead failed to converged
		o1 = try(optim(s, lh, method = "SANN",
			control = list(fnscale = -1, parscale = parscale, maxit = c$sann.cycles),
		));
		#if (class(o1) == "try-error") return(NA);
		if (o$convergence > 0 || o1$value > o$value) o = o1;
		}
		o$par.os = o$par;	# parameter values in optimiztation space
		o$par = f$parsMap(f$self, o$par);
		o
	});

	if (all(is.na(os))) return(NA);
	vs = sapply(os, function(o){o$value});
	arg.max = which.max(vs);
	estimate = os[[arg.max[1]]];
	fisher = list();
	#if (!is.null(c$computeFisher) & c$computeFisher)
	if (!is.null(c$computeFisher)) fisher = estimate.fisher(d, estimate, fisher.eps = 1e-1);
	r = c(estimate, fisher);
	r
}

lhMlEstimator = function(lhWrapper = NULL, start = NULL, c = NULL, ...) {
	f = lhGetWrapper(lhWrapper, c$self, ...);
	r = if (length(f$parsStart(f$self)) > 1) {
		lhMlEstimatorMD(lhWrapper, start, c, ...);
	} else if (length(f$parsStart(f$self)) == 1) {
		lhMlEstimatorOD(lhWrapper, start, c, ...);
	} else { # null hypothesis w/o nuisance parameters
		r = f$lh(f$self)();
	}
	r
}


lhLRtest = function(lhWrapper = NULL, start = NULL, c = list(do.sann = F, sann.cycles = 1000), ...) {
	f = lhGetWrapper(lhWrapper, NULL, c$self, ...);	# f$self is likelihood object and absorbs ellipsis parameters
	self = f$self;
	if (is.null(start)) start = f$parsStart(self);

	startNull = if (is.matrix(start))
		t(apply(start, 1, function(r)f$alt2null(self, r))) else
		f$alt2null(self, start);

	e.null = lhMlEstimator(sprintf("%s.%s", lhWrapper, "null"), startNull, c(c, list(self = self)));

	start = rbind(start, f$null2alt(self, e.null$par));
	e.alt = lhMlEstimator(lhWrapper, start, c(c, list(self = self)));

	# <p> calcualte degrees of freedom
	st = f$parsStart(self);
	df = length(st) - length(f$alt2null(self, st));
	stat =  2 * (e.alt$value - e.null$value);
	list(ll.null = e.null$value, ll.alt = e.alt$value,
		test.stat = stat, p = 1 - pchisq(stat, df), df = df, par.null = e.null$par, par.alt = e.alt$par
	)
}

#
#	lh-functions based on likelihood specification
#

# Example: see dataAnalysis.R in hwe project
# Example: binomial distribution
# lhBin = function(p, k, N)dbinom(k, N, p)
# spec_lhBin = list(
# 	ll = "lhBin",
# 	alt = list(
# 		start = c(.5),	# also specifies number of parameters
# 		pars = list(list(name = "rho", type = "freq"))
# 	),
# 	null = list(
# 		start = c(.5),	# assume same likelihood and therefore #pars from alternative
# 		parsFree = 0	# alternative: list free parameters or specify tail from alt
# 	)
# );
# r = lhMl(spec_lhBin)

#define a function
toF = function(expr, args, env = parent.frame()) {
	as.function(c(args, expr), env)
}

logitI = expit = function(x, min = 0, max = 1) { (max - min)/(1 + exp(-x)) + min }
expitD = toF(D(expression((max - min)/(1 + exp(-x)) + min), 'x'), list(x = NULL, min = 0, max = 1));
logit = function(x, min = 0, max = 1) { log((x - min)/(max - x)) }
logitD = toF(D(expression(log((x - min)/(max - x))), 'x'), list(x = NULL, min = 0, max = 1));
# templates assuming X as argument, p as parameter description list
lhArgMappers = list(
	freq =		"expit(X)",
	int =		"expit(X, min, max)",
	real =		"X",
	positive =	"log1p(exp(X))"
);
lhArgMappersD = list(
	freq =		NULL, #D(expression(expit(x), 'x')),
	int =		"expit(X, min, max)",
	real =		"X",
	positive =	"log1p(exp(X))"
);
lhArgMappersI = list(
	freq =	"logit(X)",
	int =	"logit(X, min, max)",
	real =	"X",
	positive =	"log(expm1(X))"
);
lhSpecificationDefaults = list(
	# embed null-parameter into alt-parameter space: variables: npars, parsFree, s (specification),
	#	p: input parameters
	#	<i>: optimization: substitute literals from start
	default = list(mapper = 'c(c(ARGS_FREE), c(ARGS_BOUND))', mapperInv = 'c(ARGS_FREE)')
);
# richest: richest parametrization of the likelihood
# lhInterface: call the likelihood function with a vector (vector) or with separate arguments formula
#	the paramters (inline)
lhSpecificationDefault = list(richest = 'alt', lhInterface = 'vector');
lhSpecificationInterfaces = list(
	vector = 'function(p, ...) { pm = mapper(p); if (any(abs(pm) > 1e10)) return(-Inf); lf(pm, ...) }',
	inline = 'function(p, ...) { pm = mapper(p); if (any(abs(pm) > 1e10)) return(-Inf); lf(ARGS_INLINE, ...) }'
);

#
#	<p> logit derivatives
#simulations in 2014-07-Borstkanker/src/borstKankerExp.R
logExpit1 = function(x)log(expit(x))
logExpit = logExpit2 = function(x)-log1p(exp(-x))
logitExp1 = function(x)logit(exp(x))
logitExp = logitExp2 = function(x)-log(expm1(-x))

logExpit1m1 = function(x)log(1 - expit(x))
logExpit1m = logExpit1m2 = function(x)-log1p(exp(x))
logit1mExp1 = function(x)logit(1 - exp(x))
logit1mExp = logit1mExp2 = function(x)log(expm1(-x))


#
#	<p> helper functions
#

# mappers for individual parameters
# ps: list of parameters
# mappers: mapper templates to used
# target: name of variable on which to apply
# idcs: indeces to iterate
lhMapperPars = function(ps, mappers, target = 'p', idcs = 1:length(ps)) {
	maps = if (length(idcs) == 0) c() else sapply(idcs, function(i) {
		p = ps[[i]];
		a = gsub("X", sprintf("%s[%s]", target,
			deparse(if (length(p$entries)) p$entries else i)), mappers[[p$type]]);
		a = mergeDictToString(ps[[i]]$args, a);
		a
	});
	r = paste(maps, collapse = ", ");
	r
}

# <!> auto inverse mapping has to heed mapperPost time of application
# mappers map individual arguments, mapper sub-sequently maps the whole vector
lhMapperFunction = function(s, mappers, mapper) {
	free = 1:s$parsFree;	# idcs of free variables
	bound = if(s$parsFree < s$npars) (s$parsFree + 1):s$npars else c(); # idcs of bound variables
	mStr = sprintf('function(p){%s}',
		mergeDictToString(list(
			ARGS_FREE = lhMapperPars(s$pars, mappers, 'p', free),
			ARGS_BOUND = lhMapperPars(s$pars, mappers, 'start', bound)
	), mapper));
	mf = with(s, eval(parse(text = mStr)));
	mf
}

lhMapperFunctions = function(s) {
	r = list(
		mapper = lhMapperFunction(s, lhArgMappers, s$mapper),
		mapperInv = lhMapperFunction(s, lhArgMappersI, s$mapperInv)
	);
	r
}

#' Build wrapper function around likelihood
#'
#' @param template parameter specification used as template (usually richest parametrization tb reduced
#'	for other hypotheses)
lhPreparePars = function(pars, defaults = lhSpecificationDefaults$default, spec = lhSpecificationDefault,
	template = pars) {
	# <p> determine free parameters
	t = merge.lists(defaults, pars);
	npars = length(template$pars);
	if (!is.null(t$parsFree)) {
		t$pars = if(t$parsFree == 0) list() else template$pars[(npars - t$parsFree): npars];
	}
	if (is.null(t$start)) t$start = template$start;
	if (is.null(t$parsFree)) t$parsFree = length(t$pars);

	# <p> construct mapped likelihood function
	fs = mergeDictToString(
		list(ARGS_INLINE =
			paste(sapply(1:npars, function(i) { sprintf("pm[%s]",
				deparse(if (length(template$pars[[i]]$entries)) template$pars[[i]]$entries else i)) }
			), collapse = ', ')),
		lhSpecificationInterfaces[[spec$lhInterface]]
	);
	t = merge.lists(t, list(npars = npars));
	t = merge.lists(t, lhMapperFunctions(t), list(lf = get(spec$ll)));
	f = with(t, eval(parse(text = fs)));
	t = merge.lists(t, list(npars = npars, lh = f));
	t
}

# types: names of specifications for which to define wrapped functions
# richest: name of specification for model that includes a superset of parameters of all other types
lhPrepare = function(s, types = c('null', 'alt')) {
	# <p> preparation
	s = merge.lists(lhSpecificationDefault, s);
	ri = s[[s$richest]];
	# number of parameter groups
	npars = length(ri$pars);
	# number of parameters of the likelihood function
	#Npar = sum(list.kp(ri$pars, 'entries', template = 1));
	# <p> build wrappers
	m = nlapply(types, function(type) {
		defaults = merge.lists(lhSpecificationDefaults$default, lhSpecificationDefaults[[type]]);
		lhPreparePars(s[[type]], defaults, s, template = ri)
	});
	m = merge.lists(s, m);
	m
}
# <N> free parameters come first
lhFreePars = function(s, p)with(s, {
	r = if (parsFree > 0) {
		idcs = unlist(list.kp(s$pars[1:parsFree], 'entries'));
		if (length(idcs) == 0) idcs = 1:parsFree;
		p[idcs]
	} else c();
	r
})

# second numeric derivative of x
Dn2f = function(f, x, ..., eps = 1e-5) {
	(f(x + 2*eps, ...) + f(x - 2*eps, ...) - 2*f(x, ...))/(4*eps^2)
}

lhEstMLRaw = function(t, start = NULL, ..., optim_method = 'BFGS') {
	if (is.null(start)) start = t$start;
	for (method in optim_method) {
		o = try(OptimizeMultiStart(t$mapperInv(start), t$lh, method = method, ...));
		if (!('try-error' %in% class(o))) break();
	}
	o$par = t$mapper(o$par);
	o$ci$ciL = t$mapper(o$ci$ciL);
	o$ci$ciU = t$mapper(o$ci$ciU);
	o
}

lhEstML = lhMl = function(s, start = NULL, type = 'alt', ..., optim_method = 'BFGS') {
	# <p> mapping of parameters
	s = lhPrepare(s, types = type);
	lhEstMLRaw(s[[type]], start = start, ..., optim_method = optim_method)
}

lfPrepare = function(s, ...) {
	lhParsOrig = list(...);
	prepare = sprintf('%s%s', s$ll, c('prepare', '_prepare'));
	prepareExists = min(which(sapply(prepare, exists)));
	lhPars = if (prepareExists < Inf) get(prepare[prepareExists])(...) else lhParsOrig;
	lhPars
}

# specification based LR-test
lhTestLR = function(s, startNull = NULL, startAlt = NULL, types = c('null', 'alt'), ...,
	optim_method = 'BFGS', addTypeArg = F) {
	# <p> general preparation
	s = lhPrepare(s, types = types);
	null = s[[types[1]]];
	alt = s[[types[2]]];

	# <p> specific preparation (user defined)
	lhPars = lfPrepare(s, ...);
	
	# <p> null hypothesis
	if (is.null(startNull))
		startNull = if(null$parsFree == 0) NULL else matrix(lhFreePars(null, null$start), nrow = 1);
	lhEstMLRawArgs = c(list(t = null, start = startNull), lhPars, list(optim_method = optim_method));
	if (addTypeArg) lhEstMLRawArgs = c(lhEstMLRawArgs, list(lh_type__ = 'null'));
	o0 = do.call(lhEstMLRaw, lhEstMLRawArgs);

	# <p> alternative hypothesis
	if (is.null(startAlt)) {
		# build from fit under the null
		parNull = lhFreePars(null, o0$par);
		startAlt = matrix(c(parNull, alt$start[(length(parNull) + 1):length(alt$start)]), nrow = 1);
	}
	lhEstMLRawArgs = c(list(t = alt, start = startAlt), lhPars, list(optim_method = optim_method));
	if (addTypeArg) lhEstMLRawArgs = c(lhEstMLRawArgs, list(lh_type__ = 'alt'));
	o1 = do.call(lhEstMLRaw, lhEstMLRawArgs);

	# <p> calcualte degrees of freedom
	df = length(alt$start) - length(lhFreePars(null, o0$par));
	stat =  2 * (o1$value - o0$value);
	r = list(ll.null = o0$value, ll.alt = o1$value,
		test.stat = stat, p = 1 - pchisq(stat, df), df = df, par.null = o0$par, par.alt = o1$par,
		lh.pars = lhPars, lh.pars.orig = lhParsOrig
	);
	r
}

#
#	<p> latest iteration of LH wrapper
#

lhPrepareFormula = function(s, type, formula, data, ...) {
	# <o> compute on subset of data <N> cave: missingness
	X = model.matrix(model.frame(formula, data = data), data = data);

	# <p> expand paramters
	t = s[[type]];
	ps = t$pars;
	fparsI = which(list.key(ps, 'name') == 'formula');
	fpars = ps[[fparsI]];	# formula pars
	ps[[fparsI]] = merge.lists(ps[[fparsI]], list(name = 'beta', count = ncol(X)));
	# <p> determine slots
	counts = cumsum(list.key(ps, 'count'));
	countsStart = pop(c(1, counts + 1));
	ps = lapply(seq_along(ps), function(i)merge.lists(ps[[i]], list(entries = countsStart[i]:counts[i])));
	# <p> determine start
	start = avu(sapply(ps, function(p)rep(p$start, p$count)));
	# <p> map pars
	t$pars = ps;
	t = lhPreparePars(t, spec = merge.lists(lhSpecificationDefault, s));
	t$start = start;
	t
}

lhMlFormula = function(s, formula, data, type = 'formula', ..., optim_method = 'BFGS') {
	# <p> mapping of parameters
	t = lhPrepareFormula(s, type, formula, data, ...);
	# <p> extra args
	lhPars = lfPrepare(s, formula = formula, data = data, ...);
	# <p> call optimizer
	lhEstMLRawArgs = c(list(t = t, start = s$start), lhPars, list(optim_method = optim_method));
	r = try(do.call(lhEstMLRaw, lhEstMLRawArgs), silent = T);
print(r);
	if (class(r) == 'try-error') r = list(par = rep(NA, length(t$start)), value = NA, convergence = 1);
	r
}

#
# <p> model manipulation
#

response.is.binary = function(r) {
	vs = sort(unique(r));
	if (length(vs) != 2) F else all(vs == c(0, 1));
}



#
#	<p> clustered data
#

#
# <p> describe relationships (genetic) given a relational (database) model
#

# given relatedness in a data frame of ids and clusterIds, return a list of clusters containing ids
# clusterRelation2list_old = function(r, idName = "id", idClusterName = "idFam", byIndex = T) {
# 	r = r[, c(idName, idClusterName)];
# 	ns = sort(unique(r[, 2]));
# 	# <p> build clusters
# 	clusters = sapply(ns, function(e)list());		# holds members of clusters
# 	names(clusters) = ns;
# 	# <!> we can iterate the list, given it is ordered lexicographically
# 	for (i in 1:(dim(r)[1])) {
# 		clN = as.character(r[i, 2]);
# 		clusters[[clN]] = unlist(c(clusters[[clN]], ifelse(byIndex, i, as.character(r[i, 1]))));
# 	}
# 	clusters
# }
clusterRelation2list = function(r, idName = "id", idClusterName = "idFam", byIndex = T) {
	r = r[, c(idName, idClusterName)];
	clusters = nlapply(sort(unique(r[[idClusterName]])), function(n) {
		idcs = which(r[[idClusterName]] == n);
		c = if (byIndex) idcs else r[[idName]][idcs];
		c
	});
	clusters
}

# permute clusters of identical size and within clusters
# cluster specification as given by clusterRelation2list assuming byIndex = T
# returned permutation is relative to refIds
permuteClusters = function(cls, refIds = NULL, selectIds = NULL) {
	# allow to filter ids from cluster specification
	if (!is.null(selectIds)) {
		cls = lapply(cls, function(cl)intersect(cl, selectIds));
		cls = clusters[sapply(cls, length) > 0];
	}
	cSizes = sapply(cls, function(e)length(e));
	# which cluster sizes are present in the data set?
	sizes = unique(cSizes);
	# indexable list of ids
	refIds = if (is.null(refIds)) sort(unlist(cls));
	# final permutation of refIds, such that refIds[perm] gives new order
	perm = 1:length(refIds);

	for (s in sort(sizes, decreasing = T)) {	# permute cluster of same size, permute within cluster
		clsS = which(cSizes == s);
		p1 = sample(1:length(clsS));	# permute clusters
		for (i in 1:length(clsS)) {
			p2 = sample(1:s);
			# <p> indeces that are to be replaced
			indT = which.indeces(cls[[clsS[i]]], refIds);
			# <p> indeces where the replacement comes from
			indF = which.indeces(cls[[clsS[p1[i]]]][p2], refIds);
			# <p> save partial permutation
			perm[indT] = indF;
		}
	}
	perm
}

# clusters is a vector with cluster ids
clustersPermute = function(cls) {
	permuteClusters(clusterRelation2list(data.frame(id = 1:length(cls), idFam = cls)))
}

#
#	<p> wrap model fitting for lm/glm/gee fitters
#

#library("geepack");	# <i> move to init method
regressionMethods = list(
	# assume formula to contain random effect
	glmr = list(
		fit = function(formula, data, clusterCol = NULL, ...) {
			glmer(formula, data = data, ...)
		},
		compare = function(m1, m0){
			a = anova(m0$r, m1$r, test = "Chisq");
			list(anova = a, m0 = m0, m1 = m1,
				#p.value = a[["P(>|Chi|)"]][2],
				p.value = a[['Pr(>Chisq)']][2],	# as of R 2.15.1
				effects0 = coefficients(summary(m0$r))[, "Estimate"],
				sdevs0 = coefficients(summary(m0$r))[, "Std. Error"],
				effects1 = coefficients(summary(m1$r))[, "Estimate"],
				sdevs1 = coefficients(summary(m1$r))[, "Std. Error"]
			)
		}
	),
	# use cluster column <!> untested
	glmrcl = list(
		fit = function(formula, data, clusterCol = NULL, ...) {
			f = update(formula, as.formula(Sprintf('~ . + (1|%{clusterCol}s)')));
			glmer(f, data = data, ...)
		},
		compare = function(m1, m0){
			a = anova(m0$r, m1$r, test = "Chisq");
			list(anova = a, m0 = m0, m1 = m1,
				#p.value = a[["P(>|Chi|)"]][2],
				p.value = a[['Pr(>Chisq)']][2],	# as of R 2.15.1
				effects0 = coefficients(summary(m0$r))[, "Estimate"],
				sdevs0 = coefficients(summary(m0$r))[, "Std. Error"],
				effects1 = coefficients(summary(m1$r))[, "Estimate"],
				sdevs1 = coefficients(summary(m1$r))[, "Std. Error"]
			)
		}
	),
	glm = list(
		fit = function(formula, data, clusterCol = NULL, ...)glm(formula, data = data, ...),
		compare = function(m1, m0){
			a = anova(m0$r, m1$r, test = "Chisq");
			list(anova = a, m0 = m0, m1 = m1,
				#p.value = a[["P(>|Chi|)"]][2],
				p.value = a[['Pr(>Chi)']][2],	# as of R 2.15.1
				effects0 = coefficients(summary(m0$r))[, "Estimate"],
				sdevs0 = coefficients(summary(m0$r))[, "Std. Error"],
				effects1 = coefficients(summary(m1$r))[, "Estimate"],
				sdevs1 = coefficients(summary(m1$r))[, "Std. Error"]
			)
		}
	),
	lm = list(
		fit = function(formula, data, clusterCol = NULL, ...)lm(formula, data = data, ...),
		compare = function(m1, m0){
			a = anova(m0$r, m1$r);
			list(anova = a, m0 = m0, m1 = m1, p.value = a[["Pr(>F)"]][2],
				effects0 = coefficients(summary(m0$r))[, "Estimate"],
				sdevs0 = coefficients(summary(m0$r))[, "Std. Error"],
				effects1 = coefficients(summary(m1$r))[, "Estimate"],
				sdevs1 = coefficients(summary(m1$r))[, "Std. Error"]
			)
		}
	),
	gee = list(
		fit = function(formula, data, clusterCol, family = binomial(), ...) {
			if (!length(formula.covariates(formula))) return(NULL);
			# geeglm needs ordered clusterIds <!> <N> now pre-processed
			#data = data[order(data[[clusterCol]]), , drop = F];
			#data = Df_(data, headerMap = listKeyValue(clusterCol, 'regrCompClusters'));
			#names(data)[which(names(data) == clusterCol)] = "regrCompClusters";	# hack to make geeglm work
			#geeClusters = as.integer(data[[clusterCol]]);
			myFamily = get(family$family)();
			r = geeglm(formula, data = data, id = regrCompClusters, ...);
			r
		},
		compare = function(m1, m0){
			#a = if (is.null(m0)) anova(m1$r) else anova.geeglm(m0$r, m1$r);
			a = if (is.null(m0)) anova(m1$r) else { anova(m0$r, m1$r); }
			list(anova = a, m0 = m0, m1 = m1, p.value = a[["P(>|Chi|)"]][1],
				effects0 = coefficients(summary(m0$r))[, "Estimate"],
				sdevs0 = coefficients(summary(m0$r))[, "Std.err"],
				effects1 = coefficients(summary(m1$r))[, "Estimate"],
				sdevs1 = coefficients(summary(m1$r))[, "Std.err"]
			)
		}
	),
	geebin = list(
		fit = function(formula, data, clusterCol, ...) {
			if (!length(formula.covariates(formula))) return(NULL);
			# <!> workaround for generic approach that gives obscure error message
			#	gee with family = binomial()
			r = geeglm(formula, data = data, id = regrCompClusters, family = binomial(), ...);
			r
		},
		compare = function(m1, m0){
			#a = if (is.null(m0)) anova(m1$r) else anova.geeglm(m0$r, m1$r);
			a = if (is.null(m0)) anova(m1$r) else { anova(m0$r, m1$r); }
			list(anova = a, m0 = m0, m1 = m1, p.value = a[["P(>|Chi|)"]][1],
				effects0 = coefficients(summary(m0$r))[, "Estimate"],
				sdevs0 = coefficients(summary(m0$r))[, "Std.err"],
				effects1 = coefficients(summary(m1$r))[, "Estimate"],
				sdevs1 = coefficients(summary(m1$r))[, "Std.err"]
			)
		}
	),
	surv = list(
		fit = function(formula, data, clusterCol, ...)
			list(r = coxph(formula, data = data), formula = formula, data = data),
		compare = function(m1, m0){
			vars0 = formula.predictors(m0$r$formula, m0$r$data);
			Nvars0 = length(vars0);
			vars1 = formula.predictors(m1$r$formula, m1$r$data);
			Nvars1 = length(vars1);
			p.value = if (Nvars0 == 0 || (Nvars1 - Nvars0) == 1) {
				v = setdiff(vars1, vars0);
				summary(m1$r$r)$coefficients[v, 'Pr(>|z|)']
			} else {
				anova(m0$r$r, m1$r$r)[2, 'P(>|Chi|)']
			}
			r = list(
				p.value = p.value,
				effects0 = coefficients(summary(m0$r$r)),
				effects1 = coefficients(summary(m1$r$r))
			);
			r
		}
	)
);

CompleteRows = function(f1, data) {
	vars = if (is.character(f1)) f1 else all.vars(as.formula(f1));
	rows = apply(data[, vars, drop = F], 1, function(r)all(!is.na(r)));
	r = which(rows);
	r
}

# <!> clusterIds is needed as argument although just forwarded
regressionFit = function(f, data, type, ...) {
	r = regressionMethods[[type]]$fit(f, data, ...);
	list(type = type, r = r)
}

regressionCompare = function(m1, m0) {
	r = regressionMethods[[m1$type]]$compare(m1, m0);
	r
}

regressionCompareModelsRaw = function(f1, f0, data, type = "lm", clusterCol = NULL, ...) {
	f0 = as.formula(f0);
	f1 = as.formula(f1);
	# <p> jointly trim data according to missing data
	#rows = which(apply(data[, c(formula.vars(f1), clusterCol)], 1, function(r)all(!is.na(r))));
	# more robust version
	row.names(data) = NULL;
	#rows = as.integer(row.names(model.frame(f1, data = data)));
	# robust for random effects
	#rows = apply(data[, all.vars(as.formula(f1)), drop = F], 1, function(r)!any(is.na(r)));
	#d0 = data[rows, ];
	# <p> general cluster processing: order (<!> required for geeglm, rename <!> required for geeglm
	if (Nif(clusterCol)) {
		data = data[order(data[[clusterCol]]), , drop = F];
		data = Df_(data, headerMap = listKeyValue(clusterCol, 'regrCompClusters'));
	}
	d0 = model_matrix_from_formula(f1, data, returnComplete = T)$data;

	# <p> fit and compare models
	m1 = regressionFit(f1, data = d0, type = type, clusterCol = clusterCol, ...);
	m0 = regressionFit(f0, data = d0, type = type, clusterCol = clusterCol, ...);
	a = regressionCompare(m1, m0);
	a
}

permuteDefault = list(
	p.value = 0, sdev.rel = .3, Nchunk = 1e3,
	nuisanceCovariates = NULL, .clRunLocal = T
);
# idCol: used for permutation: column specifying identiy of individuals: could be filled automatically <i>
# permute:
#	sdev.rel: sdev relative to p.value to decide how often to permute
regressionCompareModels = function(f1, f0, data, type = "lm", clusterCol = NULL, ...,
	permute = permuteDefault) {
	permute = merge.lists(permuteDefault, permute);
	r = regressionCompareModelsRaw(f1, f0, data, type, clusterCol, ...);

	if (!is.null(r) && !is.null(r$p.value) && !is.na(r$p.value) && r$p.value < permute$p.value)
		r = regressionCompareModelsPermuted(f1, f0, data, type, clusterCol, ..., permute = permute);
	r
}

#
#	<p> permuted cluster regression
#

regressionCompareModelsPermuted = function(f1, f0, data, type = "lm", clusterCol = "cluster", ...,
	idCol = "id", permute = permuteDefault, Nmax = 1e5) {
	# <p> data p-value
	a.data = regressionCompareModelsRaw(f1, f0, data, type, clusterCol = clusterCol, ...);
	p.data = a.data$p.value;
	# <p> logging
	Log(sprintf("Permuting Regression: %s [p = %.2e]", paste(as.character(f1), collapse = " "), p.data), 4);
	# <p> permutation variables indeces
	pvs = setdiff(formula.covariates(f1), permute$nuisanceCovariates);

	# <p> precompute cluster data structure
	cls = clusterRelation2list(data.frame(id = 1:length(data[[clusterCol]]), idFam = data[[clusterCol]]))
	ps = NULL;
	d0 = data;
	# adaptive permutation
	repeat {
		ps0 = clapply(1:permute$Nchunk, function(i, f1, f0, data, type, clusterCol, cls, pvs){
			d0[, pvs] = if (is.null(clusterCol)) d0[sample(1:(dim(data)[1])), pvs] else
				d0[permuteClusters(cls), pvs];
			r = regressionCompareModelsRaw(f1, f0, d0, type, clusterCol, ...);
			r$p.value
		}, f1 = f1, f0 = f0, data = data, type = type, clusterCol = clusterCol, cls = cls, pvs = pvs,
		.clRunLocal = permute$.clRunLocal);
		ps0 = na.exclude(as.numeric(ps0));
		ps = c(ps, ps0);
		#print(ps[1:100]);
		p.emp = fraction(ps <= p.data);
		# <p> stopping criterion
		p.break = if (p.emp == 0) 1 / length(ps) else p.emp;
		sdev.rel = sqrt(p.break * (1 - p.break) / length(ps)) / p.break;
		#print(list(sd = sdev.rel * p.break, sd.rel = sdev.rel, p = p.emp));
		if (sdev.rel <= permute$sdev.rel) break;
		# <p> final stop
		if (length(ps) >= Nmax) break;
	};
	r = list(f1 = f1, f0 = f0, p.value = p.emp, p.data = p.data, anova = a.data$anova, ps = ps);
	r
}

# permute covariates in order to obtain empirical p-values
#	f1: model formula alternative
#	f0: model formula hypothesis
#	M: number of permutations
regressionCompareModelsEmp = function(f1, f0, data, nuisanceCovariates = c(), type = "lm", M = 1e3, ...,
	idName = "id", idClusterName = "cluster", .clRunLocal = T) {
	r = regressionCompareModelsPermuted(f1, f0, type, ..., clusterCol = idClusterName, idCol = idName,
		permute = list(Nchunk = M, nuisanceCovariates = nuisanceCovariates, .clRunLocal = .clRunLocal));
	r
}

# data: data.frame
# stat: function computing test statistic
# vars: formula for permuation
# Nperm: number of permutations
# Pvalue: c('upper', 'lower', 'two.tailed')
permute = function(data, stat, vars, ..., Nperm = 5e3, Pvalue = 'lower', na.rm = T, fracBadStatThres = .01,
	returnT = TRUE) {
	perm.vars = all.vars(as.formula(vars));
	f = function(i, ...) {
	};
	Ts = Sapply(0:Nperm, function(i, data, ...) {
		if (i > 0) data[, perm.vars] = data[sample(nrow(data)), perm.vars];
		stat(data, ...)
	}, data = data, ...);

	fracBadStatistics = mean(is.na(Ts[-1]));
	if (is.na(Ts[1]) || fracBadStatistics >= fracBadStatThres) return(list(p.value = NA));
	Ts = Ts[!is.na(Ts)];
	Tdata = Ts[1];
	Ts = Ts[-1];
	Plower = (1 + sum(Ts <= Tdata)) / Nperm;
	Pupper = (1 + sum(Ts >= Tdata)) / Nperm;
	p.value = switch(Pvalue,
		lower = Plower,
		upper = Pupper,
		two.tailed = 2 * min(Plower, Pupper)
	);
	r = if (returnT)
		list(p.value = p.value, t.data = Tdata, t.perm = Ts) else
		list(p.value = p.value, t.data = Tdata);
	r
}

calibrate = function(f1, data, regression = glm, family = binomial(), quantiles = seq(0, 1, 0.25), ...) {
	stat0 = regression(f1, data, family = family);
	scores = predict(stat0, data, type = 'link');
	ps = predict(stat0, data, type = 'response');
	y = data[[formula.response(f1)]];
	groups = quantile(scores, probs = quantiles, na.rm = FALSE);

	calibration = lapply(2:length(quantiles), function(i) {
		group = scores > groups[i - 1] & scores <= groups[i];
		c(mean(y[group]), mean(ps[group]), sum(group))
	});
	c0 = Df(do.call(rbind, calibration),
		row.names = paste0(sprintf('%.1f', quantiles[-1] * 100), '%'),
			names = c('truth', 'prediction', 'N'));
	r = list(calibration = c0, model = summary(stat0), quantiles = groups[-1], scores = scores);
	r
}

#
#	<p> P-values
#

alphaLevels = function(ps, Nalpha = 40, alphas = seq(1/Nalpha, 1 - 1/Nalpha, length.out = Nalpha)) {
	lvls = setNames(sapply(alphas, function(alpha)mean(ps < alpha)), round(alphas, 2));
	lvls
}

testSizes = function(ps, levels = seq(0, 1, length.out = 21), asDf = FALSE, scale = 1, round = NULL) {
	sizes = sapply(levels, function(l) fraction(ps <= l));
	r = if (asDf)
		Df(level = levels, size = sizes) else
		SetNames(sizes, Sprintf("%{l}.1f%%", l = levels*100));
	r = r * scale;
	if (notE(round)) r = round(r, round);
	return(r);
}
TestSizes = function(ps, levels = seq(0, 1, length.out = 21), asDf = FALSE, scale = 100, round = 2)
	testSizes(ps, levels, asDf, scale, round)
testSizePlot = function(ps, ..., diagonal = 'red') {
	psP = testSizes(ps, ..., asDf = TRUE);
	ggplot(psP, aes(x = levels, y = size)) + geom_line() +
		geom_abline(intercept = 0, slope = 1, col = diagonal)

}

#
#	<p> error propagation
#

# as derived from the RB project and tested therein

errProd = function(x, sdx, y, sdy, covxy = 0) {
	sdp = (x * y) * sqrt((sdx/x)^2 + (sdy/y)^2 + 2 * sdx * sdy * covxy);
	sdp
}

errFrac = function(x, sdx, y, sdy, covxy = 0) {
	sdp = (x / y) * sqrt((sdx/x)^2 + (sdy/y)^2 - 2 * sdx * sdy * covxy);
	sdp
}

errSum = function(sdx, cx = 1, sdy = 0, cy = 1, covxy = 0) {
	sds = sqrt((cx *sdx)^2 + (cy * sdy)^2 + 2 * cx * cy * covxy);
	sds
}

#
#	<§> some general statistical transformations
#

# convert confidence interval to standard dev based on a normality assumption
ciToSd_old = function(ci.lo, ci.up, level = .95) {
	# upper centered limit
	ciU = ci.up - mean(c(ci.lo, ci.up));
	span = ci.up - ci.lo;
	# corresponding sd
	sd = Vectorize(inverse(function(s)qnorm(1 - (1 - level)/2, 0, s), interval = c(0, span * 8)))(ciU);
	sd
}
ciToSd = function(ci.lo, ci.up, level = .95) {
	if (missing(ci.up)) {
		if (is.null(dim(ci.lo))) {
			ci.up = ci.lo[2];
			ci.lo = ci.lo[1];
		} else {
			ci.up = ci.lo[, 2];
			ci.lo = ci.lo[, 1];
		}
	}
	alpha = 1 - level;
	widthStd = qnorm(1 - alpha/2) * 2;
	(ci.up - ci.lo)/widthStd
}
ciToP = function(ci.lo, ci.up, level = .95, one.sided = F, against = 0) {
	sd = ciToSd(ci.lo, ci.up, level)
	P = peSdToP((ci.lo + ci.up)/2 - against, sd, one.sided);
	P
}
# convert point estimate and SD to p-value (assuming normality)
peSdToP = function(beta, sd, one.sided = F) {
	pnorm(-abs(beta), 0, sd, lower.tail = T) * ifelse(one.sided, 1, 2);
}

betaSdToCi = ciFromBetaSdev = function(beta, sdev, level = .95) {
	r = list(effect = beta,
		lower = qnorm((1 - level)/2, beta, sdev, lower.tail = T),
		upper = qnorm((1 - level)/2, beta, sdev, lower.tail = F)
	);
	r
}

ciFromSummary = function(s, var, level = .95) {
	cs = as.matrix(coefficients(s)[var, ]);
	ciFromBetaSdev(cs[, 'Estimate'], cs[, 'Std. Error'], level = level);
}

pFromBetaSd = function(beta, sd, null = 0)pnorm(null, abs(beta), sd)
sdFromBetaP = function(beta, p)Vectorize(inverse(function(s)peSdToP(beta, s), interval = c(0, 10)))(p);
betaPtoCi = function(beta, p) {
	sd = sdFromBetaP(beta, p);
	ciFromBetaSdev(beta, sd)
}
betaVarToCi = function(beta, var)ciFromBetaSdev(beta, sqrt(var))

coefficientsOR = function(s, level = .95) {
	cis = do.call(cbind, ciFromSummary(s, level = level));

	# <p> ORs
	cisOR = Df_(exp(cis));
	cisOR['(Intercept)', ] = NA;	# intercept does not have an OR interpretation
	names(cisOR) = paste(names(cisOR), 'OR', sep = '');

	# <p> result
	cfs = coefficients(s);
	r = cbind(cfs, cis[, -1], as.matrix(cisOR));
	r
}
coefficientsORmodel = function(model, level = .95)coefficientsOR(summary(model), level)

#
#	meta analysis
#

# meta analysis row-wise
metaPvalue = function(ps) {
	if (!is.matrix(ps)) ps = matrix(ps, nrow = 1);
	if (!all(is.numeric(ps))) ps = apply(ps, 1:2, as.numeric);
	cs = apply(ps, 1, function(r)sum(-2*log(r)));
	psM = pchisq(cs, 2*dim(ps)[2], lower.tail = F);
	psM
}

#
#	data imputation
#

Sample = function(x, ...)if (length(x) == 1)x else sample(x, ...);

mi.simple = function(data, n.imp = 20) {
	r = lapply(1:n.imp, function(i) {
		for (v in names(data)) {
			data[is.na(data[, v]), v] = Sample(na.omit(data[, v]), count(is.na(data[, v])));
		}
		data
	})
	r
}

cross.imputer = function(imputationData, imputationVars = NULL, doExpandFactors = T) {
	if (is.null(imputationVars)) imputationVars = names(imputationData);
	f = function(data) {
		d0 = data;
		for (v in imputationVars) { # cross impute from imputationData
			d0[is.na(d0[, v]), v] = Sample(na.omit(imputationData[[v]]), count(is.na(d0[, v])));
		}
		if (doExpandFactors) d0 = dataExpandFactors(d0)[, vars];
		d0
	};
	f
}

imputeMeanVar = function(col) {
	mn = mean(col, na.rm = T);
	col[is.na(col)] = mn;
	col
}
imputeMean = function(data) {
	d1 = apply(data, 2, imputeMeanVar);
	d1
}

# impute using mice
imputeVariables = function(data, formula = NULL, vars = NULL, Nimp = 10, returnOriginal = FALSE,
	imputationColumn = 'Imputation_', returnMice = FALSE) {
	require('mice');

	# <p> preparation
	if (notE(formula)) vars = all.vars(formula.re(formula, data = data));
	invalidVars = vars[!lapply(data[, vars], class) %in% c('factor', 'numeric')];
	if (length(invalidVars) > 0) {
		stop(Sprintf('some columns neither factor nor numeric: %{v}s', v = join(invalidVars, ', ')));
	}
	s = descriptives(data[, vars]);
	
	# <p> imputation
	dataImputed = mice(data[, vars, drop = F], m = Nimp);
	dataL = Df_(complete(dataImputed, "long", include = returnOriginal),
		min_ = c('.id'), headerMap = list(`.imp` = imputationColumn), as_numeric = imputationColumn);

	dataC = DfStack(data, (Nimp + returnOriginal));
	dataC[, vars] = Df_(dataL, min_ = imputationColumn);
	dataR = Df(dataL[, imputationColumn], dataC, names = imputationColumn);
	r = list(data = dataR, summary = s, vars = vars);
	if (returnMice) r = c(r, list(mice = dataImputed));
	#r = list(data = dataImputed, summary = s, vars = vars);
	r
}

missingness = function(d, f0) {
	ns = if (missing(f0)) names(d) else
		if (class(f0) == 'formula') all.vars(f0) else f0;
	r = lapply(ns, function(n) {
		list(N = sum(is.na(d[[n]])), freq = mean(is.na(d[[n]])))
	});
	return(SetNames(do.call(rbind, r), rnames = ns));
}

Missingness = function(d, f0) {
	dM = Df_(missingness(d, f0), deep_simplify_ = TRUE);
	dM$variable = row.names(dM);
	dM$perc = Sprintf('%{perc}.1f%%', perc = dM$freq * 100);
	dM$freq = Sprintf('%{freq}.3f', dM);
	row.names(dM) = NULL;
	
	return(dM[, c('variable', 'N', 'freq', 'perc')]);
}


#
#	<p> cross validation
#

# cross validation partitions for classification data
crossValidationPartitionsClassification = function(responses, K = 15, minEls = 3, maxTries = 15) {
	N = length(responses);
	cats = unique(responses);

	for (i in 1:maxTries) {
		# random permutation
		perm = sample(1:N, N);
		# compute partitions
		parts = splitListEls(perm, K, returnElements = T);
		counts = data.frame.types(
			lapply(parts, function(p)table.n(responses[-p], categories = cats)),
			names = cats, do.rbind = T
		);
		doReject = any(apply(counts, 1, function(r)any(r < minEls)));
		if (!doReject) break;
	}
	r = if (i < maxTries) parts else {
		Log("Error: failed to find suitable cross validation partition!");
		NULL
	}
	r
}

# cross validation parititions for clustered data
# return indeces into cluster vector (cluster identities assumed to be given by integers)
# so far do not heed cluster sizes
crossValidationPartitionsClusters = function(clusters, K = 20) {
	N = length(clusters);
	# unique cluster ids
	cls = unique(clusters);
	# random permutation
	perm = Sample(cls, length(cls));
	# compute partitions
	parts = splitListEls(perm, K, returnElements = T);
	r = lapply(parts, function(p)which.indeces(p, clusters, match.multi = T));
	r
}

#
#	<p> optimization
#

nested.search = function(f, ..., key = NULL, parameters = list(p1 = c(0, 10)),
	steps = 3, Ngrid = 4, rscale = 1, return.grid = F, par.as.vector = F, .clRunLocal = rget('.clRunLocal')) {
	ps = ps0 = parameters;

	for (i in 1:steps) {
		# <p> create serach grid
		pars = lapply(ps, function(p)seq(p[1], p[2], length.out = Ngrid));
		grid = merge.multi.list(pars);
		# <p> apply function
		r = clapply(1:dim(grid)[1], function(j, grid, ...) {
			args = if (par.as.vector) list(as.vector(grid[j, ]), ...) else c(as.list(grid[j, ]), list(...));
			do.call(f, args);
		}, grid = grid, ..., .clRunLocal = .clRunLocal);
		# <p> search optimum
		values = if (is.null(key)) r else list.kp(r, key, do.unlist = T);
		opt = which.min(values * rscale);
		pt = grid[opt, ];	# optimal point in the grid search

		ps = lapply(1:length(ps), function(j){
			from = max(pt[j] - (ps[[j]][2] - ps[[j]][1])/Ngrid, ps0[[j]][1]);
			to =   min(pt[j] + (ps[[j]][2] - ps[[j]][1])/Ngrid, ps0[[j]][2]);
			c(from, to)
		});
		names(ps) = names(ps0);
	}
	r = if (return.grid) list(value = values[[opt]], par = pt, grid = r) else
		list(value = values[[opt]], par = pt, r = r[[opt]]);
	r
}
optim.nested.defaults = list(steps = 5, Ngrid = 4, rscale = 1, return.grid = F);
optim.nested = function(par = NULL, f, ..., lower = -Inf, upper = Inf, control = list())
	with(merge.lists(optim.nested.defaults, control), {
	parameters = apply(cbind(lower, upper), 1, function(r)list(r));
	r = nested.search(f, ..., parameters,
		steps = steps, Ngrid = Ngrid, rscale = rscale, return.grid = return.grid, par.as.vector = T);
	r
})

#
#	<p> correlation in data
#

Df.corr = function(df, eps = 1e-2) {
	N = dim(df)[2];
	rc = rcorr(df);
	pairs = t(sapply(which(abs(rc$r) > (1 - eps)), function(e) {
		row = ((e - 1) %/% N) + 1;
		col = ((e - 1) %% N) + 1;
		r = c(row, col);
		r
	}));
	pairs = pairs[pairs[, 1] < pairs[, 2], ];
	clusters = sub.graph(pairs);
	remove = unlist(lapply(clusters, function(e)e[-1]));
	r = list(clusters = clusters, cols.remove = remove);
	r
}

identity = function(e)e
seq.transf = function(from = 0, to = 1, length.out = 1e1, ..., transf = log, transfI = exp, eps = 1e-5) {
	s = transfI(seq(from = transf(from + eps), to = transf(to - eps), length.out = length.out, ...));
	s
}

#
#	<p> bug fixes for packages
#

completeRows = function(data)!apply(data, 1, function(r)any(is.na(r)))
completeData = function(f = ~ ., data, collapse = F, varExclude = c()) {
	vsRaw = if (is.character(f)) f else all.vars(f);
	if (any(vsRaw == '.')) vsRaw = names(data);
	vs = setdiff(vsRaw, varExclude);
	data[completeRows(data[, vs, drop = F]), if (collapse) vs else 1:ncol(data), drop = F]
}

model_matrix_from_formula = function(f, data, offset = NULL, ignore.case = F, remove.intercept = F,
	returnComplete = F, formulaAsIs = F) {
	# <p> prepare data matrices
	f1 = if (formulaAsIs) f else formula.re(f, data = data, ignore.case = ignore.case);
	f1vars = all.vars(f1);
	response = formula.response(f1);
	responseVars = all.vars(as.formula(con(response, ' ~ 1')));

	# <p> complete data
	complete = completeRows(data[, f1vars, drop = F]);
	data = droplevels(data[complete, , drop = F]);
	
	# <p> response
	responseData = if (notE(responseVars)) with(data, Eval(response)) else NULL;
	# <p> offset
	offset = if (!is.null(offset)) offset[complete] else NULL;
	# <p> model matrix
	#mm = if (exists('lFormula') && length(unlist(Regex('[+]\\s*[(]', Deparse(f1)))) > 0)
	mm = if (exists('lFormula') && length(unlist(Regex('[+]\\s*[(].*?[|].*[)]', Deparse(f1)))) > 0)
		lFormula(f1, data)$X else	# lme4
		model.matrix(f1,  model.frame(f1, data = data));
	if (remove.intercept) mm = mm[, !(dimnames(mm)[[2]] == '(Intercept)'), drop = F];

	# <p> return
	r = list(mm = mm, response = responseData, offset = offset, indeces = which(complete),
		varsCov = setdiff(f1vars, responseVars), varsResponse = responseVars);
	if (returnComplete) r = c(r, list(data = data));
	r
}
complete_from_formula = function(f, data, offset = NULL, ignore.case = F, remove.intercept = F) {
	model_matrix_from_formula(f, data, offset, ignore.case, remove.intercept)$indeces
}
complete_from_vars = function(vars, data, offset = NULL, ignore.case = F, remove.intercept = F) {
	f = as.formula(Sprintf('~ %{vars}s', vars = join(vars, ' + ')));
	model_matrix_from_formula(f, data, offset, ignore.case, remove.intercept)$indeces
}

# <!><t> tb tested refactoring as of 19.2.2019

glmnet_re_mm = function(y, mm, ..., offset = NULL, ignore.case = F, remove.intercept = F,
	lambdas = NULL, cv = T, returnGlmnet = F) {
	# <p> fit model
	r = if (cv) {
		r0 = cv.glmnet(x = mm, y = y, lambda = lambdas, ..., offset = offset);
		args = c(List(..., min_ = c('foldid', 'nfolds', 'grouped')),
			list(x = mm, y = y, lambda = r0$lambda.min, offset = offset));
# 			list(x = mm, y = y, lambda = (3*r0$lambda.min + r0$lambda.1se)/4, offset = offset));
#			list(x = mm, y = y, lambda = (r0$lambda.min), offset = offset));
		do.call('glmnet', args);
	} else glmnet(x = mm, y = y, lambda = lambdas, ..., offset = offset);
	if (returnGlmnet) r = list(glmnet = r);
	return(r);
}

glmnet_re = function(f, data, ..., offset = NULL, ignore.case = F, remove.intercept = F,
	lambdas = NULL, cv = T, returnGlmnet = F) {
	d = model_matrix_from_formula(f, data, offset, ignore.case, remove.intercept);
	r = glmnet_re_mm(d$response, d$mm, ..., offset = d$offset,
		ignore.case = ignore.case, remove.intercept = remove.intercept, lambdas = lambdas,
		cv = cv, returnGlmnet = returnGlmnet);
	return(c(r, list(formula = f)));
}
glmnet_re_refit = function(model, data, ..., var_cutoff =  1e-6, intercept = '1', impute = NULL) {
	response = formula.response(model$formula);

	if (model$df <= 1) return(list());
	# <p> scrutinize model
	coefs = model$beta;
	varsSel = row.names(coefs)[abs(as.vector(coefs)) > var_cutoff];
	varsSel = setdiff(varsSel, '(Intercept)');
	
	if (!is.null(impute) && impute == 'mean') {
		# <!> use model matrix <i>
		d0 = sapply(varsSel, function(var) {
			data[[var]][is.na(data[[var]])] = mean(data[[var]], na.rm = T);
		});
		data[, varsSel] = d0;
	}
	# <p> refit
	f = as.formula(sprintf('%s ~ %s', response, paste(c(intercept, varsSel), collapse = ' + ')));
	glm1 = glm(f, data = data, ...);
	r0 = list(glm = glm1, score = as.vector(predict(glm1, data, type = 'link')))
	r0
}

#library('glmnet');
grid.glmnet.raw = function(..., glmnet.f = cv.glmnet, max.tries = 3) {
	for (i in 1:max.tries) {
		fit = try(glmnet.f(...), silent = T);
		if (all(class(fit) != 'try-error')) break();
	}
	if (any(class(fit) == 'try-error')) stop(fit[1]);
	fit
}

grid.glmnet.control = list(steps = 4, Ngrid = 50, from = .01,  to = .8, eps = 1e-5,
		transf = identity, transfI = identity);

grid.glmnet = function(..., control = grid.glmnet.control)
	with (merge.lists(grid.glmnet.control, control), {
	# initialize
	fit = NULL;
	fromO = from;
	toO = to;
	options(warn = -1);
	for (i in 1:steps) {
		lambda = seq.transf(from, to, length.out = Ngrid + 1, eps = eps,
			transf = transf, transfI = transfI);
		fit = grid.glmnet.raw(..., lambda = sort(lambda, decreasing = T));
		from = max(fit$lambda.min - (to - from)/Ngrid, 0);
		to = fit$lambda.min + (to - from)/Ngrid;
	}
	options(warn = 0);
	# choose lambdas to contain lambda.min also covering the range between from and to
	lambda = c(
		seq.transf(fromO, toO, length.out = Ngrid + 1, eps = eps,
			transf = transf, transfI = transfI),
		fit$lambda.min
	);
	fit0 = do.call('grid.glmnet.raw', c(list(...), list(lambda = sort(lambda, decreasing = T))));
	args = List(..., min_ = c('nfolds', 'grouped'));
	fit1 = do.call('grid.glmnet.raw', c(args, list(lambda = fit$lambda.min, glmnet.f = glmnet)));
	r = fit0;
	r$glmnet.fit = fit1;
	r
})

#	f: formula, passed through formula.re
#	data: data frame
grid.glmnet.re = function(f, data, ..., offset = NULL, control = grid.glmnet.control,
	ignore.case = F, remove.intercept = T)
	with (merge.lists(grid.glmnet.control, control), {

	# <p> prepare data matrices
	f1 = formula.re(f, data = data, ignore.case = ignore.case);
	f1vars = all.vars(f1);
	response = formula.response(f1);
	complete = !apply(data[, f1vars], 1, function(r)any(is.na(r)));
	d1 = data[complete, ];
	if (!is.null(offset)) offset = offset[complete];
	mm = model.matrix(f1, model.frame(f1, data = d1));
	if (remove.intercept) mm = mm[, !(dimnames(mm)[[2]] == '(Intercept)')];
	# <p> fit model
	r = grid.glmnet(x = mm, y = d1[[response]], ..., offset = offset, control = control);
	r = c(r, list(formula = f1));
	r
})

grid_glmnet_re_refit = function(model, data, ..., var_cutoff =  1e-6, intercept = '1', impute = NULL) {
	# <p> scrutinize model
	coefs = coefficients(model$glmnet.fit);
	varsSel = row.names(coefs)[abs(as.vector(coefs)) > var_cutoff];
	varsSel = setdiff(varsSel, '(Intercept)');
	response = formula.response(model$formula);

	if (!is.null(impute) && impute == 'mean') {
		# <!> use model matrix <i>
		d0 = sapply(varsSel, function(var) {
			data[[var]][is.na(data[[var]])] = mean(data[[var]], na.rm = T);
		});
		data[, varsSel] = d0;
	}
	# <p> refit
	f = as.formula(sprintf('%s ~ %s', response, paste(c(intercept, varsSel), collapse = ' + ')));
	glm1 = glm(f, data = data, ...);
	r0 = list(glm = glm1, score = as.vector(predict(glm1, data, type = 'link')))
	r0
}

refitModel = function(model, f1, f0, data, ..., var_cutoff =  1e-6, ignore.case = F, intercept = '0') {
	# <p> prepare formula and data set
	f1 = formula.re(f1, data = data, ignore.case = ignore.case);
	f0 = formula.re(f0, data = data, ignore.case = ignore.case);
	f0covs = formula.covariates(f0);
	f1vars = all.vars(f1);
	response = formula.response(f1);
	complete = complete.cases(data[, f1vars]);	#!apply(data[, f1vars], 1, function(r)(any(is.na(r))));
	d1 = data[complete, ];
	# <p> extract data set according to model
	coefs = coefficients(model);
	varsSel = row.names(coefs)[abs(as.vector(coefs)) > var_cutoff];
	varsSel = setdiff(varsSel, '(Intercept)');
	varsSel0 = intersect(varsSel, f0covs);
	if (!length(varsSel0)) return(
		list(coefficients = coefs, anova = NA, r2 = NA, r20 = NA, raw = NA, model1 = NA, model0 = NA)
	);
	# <p> re-fit glm
	f1 = as.formula(sprintf('%s ~ %s', response, paste(c(intercept, varsSel), collapse = ' + ')));
	glm1 = glm(f1, data = d1, ...);
	f0 = as.formula(sprintf('%s ~ %s', response, paste(c(intercept, varsSel0), collapse = ' + ')));
	glm0 = glm(f0, data = d1, ...);
	# <p> anova
	a = anova(glm0, glm1, test = 'Chisq');
	# <p> R^2
	mn = mean(d1[[response]]);
	#mm = model.matrix(f1, model.frame(f1, data = d1));
	pr = as.vector(predict(glm1, d1, type = 'response'));
	#r2 = cor((pr - mn), (d1[[response]] - mn));
	r2 = cor(pr, d1[[response]]);
	pr0 = as.vector(predict(glm0, d1, type = 'response'));
	#r20 = cor((pr0 - mn), (d1[[response]] - mn));
	r20 = cor(pr0, d1[[response]]);
	# <p> raw-model fit
	fScore = as.formula(sprintf('y ~ score + %s', paste(c(intercept, varsSel0), collapse = ' + ')));
	d2 = data.frame(
		d1[, varsSel0], y = d1[[response]], score = as.vector(predict(glm1, d1))
	);
	if (length(varsSel0)) names(d2)[1:length(varsSel0)] = varsSel0;
	raw = glm(fScore, data = d2, ...);
	r = list(coefficients = coefs, anova = a, r2 = r2, r20 = r20,
		raw = summary(raw), model1 = glm1, model0 = glm0);
	r
}

#
#	<p> crossvalidation
#

# <!> tb implemented
cv_summary_lm = function(model, pred, data, ...) {
	summary(r0)$fstatistic[1]
	r = mean( (pred - data)^2 );
	r
}

cv_test_glm = function(model, formula, data, ...) {
	response = formula.response(formula);
	responseP = predict(model, data, type = 'response');
	responseD = data[[response]];
	ll = sum(log(responseP));
	ll
}

cv_train_glmnet = function(f1, data, ...) {
	glmnet_re(f1, data, ..., returnGlmnet = T);
}
cv_test_glmnet = function(model, f1, data, ...) {
	preds = predict(model$glmnet, model_matrix_from_formula(f1, data)$mm, type = 'response');
	# detect logistic model, constrain probabilities
	if (notE(model$glmnet$classnames) && length(model$glmnet$classnames) == 2) preds = minimax(preds, 0, 1);
	return(preds);
}

cv_train_glmnet_mm = function(data, ...) {
	i = FirstDef(which(dimnames(data)[[2]] == 'y'), 1);
	glmnet_re_mm(data[, i], data[, -i], ..., returnGlmnet = T);
}
cv_test_glmnet_mm = function(model, data, ...) {
	i = FirstDef(which(dimnames(data)[[2]] == 'y'), 1);
	preds = predict(model$glmnet, data[, -i], type = 'response');
	# detect logistic model, constrain probabilities
	if (notE(model$glmnet$classnames) && length(model$glmnet$classnames) == 2) preds = minimax(preds, 0, 1);
	return(preds);
}


cv_train_glm = function(f1, data, ..., family = binomial()) {
	glm(f1, data = data, family = family)
}
cv_test_glm_predict = function(model, f1, data, ..., type = 'response') {
	predict(model, data, type = type)
}


# cv_prepare = function(data, argsFrom...)
# cv_train = function(data, argsFrom...)
# cv_test = function(model, data, argsFrom...)
# @arg cv_fold number of crossvalidation folds, denotes leave -cv_fold out if negative
# <!> 19.2.2019: removed cv_prepare for lack of functionality
# <i> argument routing to train/test

crossvalidate = function(cv_train, cv_test,
	data, cv_fold = 20, cv_repeats = 1, ..., parallel = F, align_order = TRUE) {
	if (cv_fold == 0) stop('crossvalidate: cv_fold must be an integer != 0');
	if (!parallel) Lapply = lapply;
	N = nrow(data);
	r = Lapply(1:cv_repeats, function(i, ...) {
			perm = Sample(1:N, N);
			# compute partitions
			fold = if (cv_fold > 0) cv_fold else as.integer(N/-cv_fold);
			parts = splitListEls(perm, fold, returnElements = T);
			o = order(unlist(parts));
			r = Lapply(parts, function(part, cv_train, cv_test, data, cv_repeats, ...) {
				d0 = data[-part, , drop = F];
				d1 = data[part, , drop = F];
				model = cv_train(..., data = d0);
				r = cv_test(model = model, ..., data = d1);
				gc();
				r
			}, cv_train = cv_train, cv_test = cv_test,
				data = data, cv_repeats = cv_repeats, ...);
			# re-establish order
			r = if (align_order
				&& all(sapply(r, class) %in% c('numeric', 'integer'))
				#&& all(sapply(r, length) == 1)) {
				&& sum(sapply(r, length)) == nrow(data)) {
				unlist(r)[o];
			} else if (align_order && all(sapply(r, class) %in% c('data.frame', 'matrix')) &&
				sum(sapply(r, nrow)) == nrow(data)) {
				#<!> untested
				#r = rbindDataFrames(r, colsFromFirstDf = T);
				r = do.call(rbind, r);
				r[o, ]
			} else if (align_order) stop("Crossvalidate: didn't know how to align order.") else r;
			gc();
			r
	}, ...);
	r
}

#
#	<p> data standardization
#

standardize = function (x, ...)UseMethod("standardize", x);
standardize.numeric = function(v, na.rm = TRUE, orig.mean = FALSE) {
	vM = mean(v, na.rm = na.rm);
	vC = v - vM;
	vS = vC / sd(vC, na.rm = na.rm);
	if (orig.mean) vS + vM else vS
}
standardize.data.frame = function(df, na.rm = TRUE)apply(df, 2, standardize.numeric, na.rm = na.rm);
standardize.matrix = function(m, na.rm = TRUE)apply(m, 2, standardize.numeric, na.rm = na.rm);

df2z = function(data, vars = names(as.data.frame(data))) {
	data = as.data.frame(data);
	df = data.frame.types(sapply(vars, function(v) {
		(data[[v]] - mean(data[[v]], na.rm = T)) / sd(data[[v]], na.rm = T)
	}), do.transpose = F);
	i = which.indeces(vars, names(data));
	d0 = data.frame(data[, -i], df);
	d0
}

lumpFactor = function(factor, minFreq = NULL, minN = 20, levelPrefix = 'l') {
	# <p> preparation
	f0 = as.factor(factor);
	t0 = table(f0);
	ls = levels(f0);
	N = length(f0);
	if (!is.null(minFreq)) minN = as.integer(minFreq * N + 0.5);
	
	# <p> lumping
	map = listKeyValue(ls, ls);
	for (i in 1:length(t0)) {
		t0 = table(factor);
		if (all(t0 >= minN) || length(t0) < 2)  break;
		# combine two smallest groups
		t1 = sort(t0);
		newLevel = sprintf('%s%d', levelPrefix, i);
		factor = as.character(factor);
		factor[factor == names(t1)[1] | factor == names(t1)[2]] = newLevel;
		map[[names(t1)[1]]] = map[[names(t1)[2]]] = newLevel;
		map[[newLevel]] = newLevel;
	}
	# <p> normalize map
	lsNew = as.character(ls);
	repeat {
		lsNew0 = lsNew;
		lsNew = as.character(map[lsNew]);
		if (all(lsNew == lsNew0)) break;
	}
	return(list(map = listKeyValue(ls, lsNew), factor = factor));
}

# lump a variable after checking other variables for non-missingness
lumpVariableOnVariables = function(data, var, vars, postfix = '_lump', minN = 20) {
	# prepare confounder afkomst
	lump = sapply(vars, function(v) {
		dvar = data[[var]][!is.na(data[[v]])];
		lump = lumpFactor(dvar, minN = minN);
		dvarNew = as.character(lump$map[as.factor(data[[var]])]);
		dvarNew[dvarNew == 'NULL'] = NA;
		as.factor(dvarNew)
	});
	d = data.frame(lump);
	names(d) = paste(var, paste(vars, postfix, sep = ''), sep = '_');
	d
}

#
#	<p> descriptive
#

compareVectors = function(l) {
	sets = names(l);
	# marginals
	r0 = nlapply(sets, function(n)c(n, length(l[[n]])));
	r1 = nlapply(sets, function(n)c(sprintf('%s-unique', n), length(unique(l[[n]]))));
	r2 = nlapply(sets, function(n)c(sprintf("%s-NA", n), sum(is.na(l[[n]]))));

	modelList = list(A = sets, B = sets);
	r3 = iterateModels(modelList, .constraint = function(A, B)(A < B), function(i, A, B) {
		r = list(
			c(sprintf("%s inter %s", A, B), length(intersect(l[[A]], l[[B]]))),
			c(sprintf("%s union %s", A, B), length(union(l[[A]], l[[B]]))),
			c(sprintf("%s min %s", A, B), length(setdiff(l[[A]], l[[B]]))),
			c(sprintf("%s min %s", B, A), length(setdiff(l[[B]], l[[A]])))
		);
		r
	}, lapply__ = lapply)$results;

	r = c(r0, r1, r2, unlist.n(r3, 1));
	r = data.frame.types(r, do.rbind = T, names = c('type', 'count'));
	r
}

pairs_std.panel.hist <- function(x, ...) {
	usr <- par("usr"); on.exit(par(usr));	# restore settings
	par(usr = c(usr[1:2], 0, 1.5) )
	h <- hist(x, plot = FALSE)
	breaks <- h$breaks; nB <- length(breaks)
	y <- h$counts; y <- y/max(y)
	rect(breaks[-nB], 0, breaks[-1], y, col = "grey", ...)
}
pairs_std.panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...) {
	usr <- par("usr"); on.exit(par(usr));
	# <p> hack to get column names
	nsAll = eval(expression(names(x)), parent.frame());
	#nsV = nsAll[c(eval(expression(i), parent.frame()), eval(expression(j), parent.frame()))];
	# <!> seems to be order (j, i) in terms of names 16.3.2021
	nsV = nsAll[c(eval(expression(j), parent.frame()), eval(expression(i), parent.frame()))];

	par(usr = c(0, 1, 0, 1))
	c0 = cor.test(x, y);
	Sx = summary(x);
	Sy = summary(y);
	Sums = RbindDfs(list(list2df(Sx), list2df(Sy)), embed = TRUE);	# summaries
	if (!('NA.s' %in% names(Sums))) Sums = Df(Sums, NA.s = NA);
	#sum = apply(cbind(Sy, Sx), 1:2, sprintf, fmt = '%.3f');
	sum = apply(Sums, 1:2, sprintf, fmt = '%.3f');
	ns = c('Min', '1st Q', 'Median', 'Mean', '3rd Q', 'Max', 'NAs'); 
	#ns = names(summary(x));
	txt <- paste0(prefix,
		sprintf("Cor: %.2f (%.2f, %.2f)\n", c0$estimate, c0$conf.int[1], c0$conf.int[2]),
		sprintf("P: %.2e\n", c0$p.value),
		sprintf("Vars: %s  %s\n", nsV[1], nsV[2]),
		#join(apply(cbind(ns, sum), 1, join, ' '), "\n")
		join(apply(cbind(ns, t(sum)), 1, join, ' '), "\n")
	)
	if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
	text(.9, 0.5, txt, adj = c(1, .5), cex = cex.cor * 1)	# used tb cex.cor * r
}

unique_fraction = function(x) (length(unique(x))/length(x))
guess_is_factor = function(x, na.rm = T, uniqueFraction = .1, Nlevels = 10) {
	all(as.integer(x) == x, na.rm = na.rm) &&
		unique_fraction(x) < uniqueFraction &&
		length(unique(x)) <= Nlevels
}


panel.scatter = function(x, y, col = par("col"), bg = NA, pch = par("pch"), 
	cex = 1, col.smooth = "red", span = 2/3, iter = 3, blotSize = cm2in(.3), ...) 
{
	#call_ = match.call(envir = parent.frame(1L));
	# it is difficult to get information on the plotted variables, for now a heuristics is used
	#	to determine factor status
	if (guess_is_factor(x) && guess_is_factor(y)) {
		pars = par(c('usr', 'xaxt', 'yaxt')); on.exit(par(pars));
		t0 = table(x, y);
		#plot(t0);
		mx = max(as.integer(t0));
		d0 = Df(merge.multi.list(dimnames(t0)), w = sqrt(as.integer(t(t0))), as_numeric = c('x', 'y'));
		#par(new = TRUE);
		#plot(d0[, c('x', 'y')], type = 'n', xaxt='n', yaxt = 'n');
		par(new = TRUE);
		with(d0,
			symbols(x = x, y = y, circles = w, inches = blotSize, ann = F, bg = "steelblue2", fg = NULL,
			 xaxt='n', yaxt = 'n')
		)
		#return();
	} else {
		# <i> call forwarding
		panel.smooth(x, y, col, bg, pch, cex, span, iter);
	}
}

pairs_std = function(...,
	lower.panel = panel.scatter, diag.panel = pairs_std.panel.hist, upper.panel = pairs_std.panel.cor) {
	pairs(...,
		lower.panel = lower.panel, diag.panel = diag.panel, upper.panel = upper.panel)
}

pairs_std_save = function(data, ..., path = NULL) {
	if (is.null(path)) path = tempfile(fileext = '.png');
		plot_save_eval(pairs_std(data), ..., plot_path = path);
	return(path);
}

pairs_wide_list = function(blocks = NULL, NperBlock = 8) {
	# <p> names
	ns = names(blocks);
	nsAll = unlist(blocks);
	blkNs = if (notE(ns)) ns else paste('Block', seq_along(blocks));
	# <p> block tiling
	Iblocks = seq_along(blocks);
	blocksPair = merge.multi(Iblocks, Iblocks);
	# <p> reorder
	Idiag = blocksPair[, 1] == blocksPair[, 2];
	blocksPair = rbind(blocksPair[Idiag, ], blocksPair[!Idiag, ]);
	# <p> tabnames
	pairNs = apply(blocksPair, 1, function(r)paste(unique(c(blkNs[r[1]], blkNs[r[2]])), collapse = ' - '));
	# <p> indeces
	horIdcs = lapply(blocksPair[, 1], function(b)which.indeces(blocks[[b]], nsAll));
	verIdcs = lapply(blocksPair[, 2], function(b)which.indeces(blocks[[b]], nsAll));
	# <p> by pair of blocks
	r = lapply(1:nrow(blocksPair), function(i)list(
		coord = blocksPair[i, ], name = pairNs[i],
		names = list(blocks[[ blocksPair[i, 1] ]], blocks[[ blocksPair[i, 2] ]]),
		horIdcs = horIdcs[[blocksPair[i, 1]]], verIdcs = verIdcs[[blocksPair[i, 2]]]
	));
	return(r);
}

pairs_wide_old = function(data, ..., blocks = NULL, NperBlock = 8, consumer = Identity1) {
	nsPairs = unlist(blocks);
	dPairs = data[, nsPairs];
	Iblocks = seq_along(blocks);
	blocksPair = merge.multi(Iblocks, Iblocks);
	plts = apply(blocksPair, 1, function(r) {
		horInd = which.indeces(blocks[[ r[1] ]], nsPairs);
		verInd = which.indeces(blocks[[ r[2] ]], nsPairs);
		pl = consumer(pairs_std(dPairs, horInd = horInd, verInd = verInd), r);
		return(list(pl))
	});
	return(unlist.n(plts, 1));
}

pairs_wide_single = function(data, b) { pairs_std(data, horInd = b$horIdcs, verInd = b$verIdcs); }
pairs_wide = function(data, ..., blocks = NULL, NperBlock = 8, consumer = Identity1) {
	prs = pairs_wide_list(blocks, NperBlock);
	dPairs = data[, unlist(blocks)];
	plts = lapply(prs, function(b)list(consumer(pairs_wide_single(dPairs, b), b)));
	return(unlist.n(plts, 1));
}


pairs_wide_save = function(data, ..., blocks = NULL, path = NULL, ext = 'jpeg') {
	prefix = if (is.null(path)) tempfile() else path;
	r = pairs_wide(data, ..., blocks = blocks, consumer = function(plot, b) {
		plot_save_eval(plot, plot_path = Sprintf('%{prefix}s-%{name}s.%{ext}', b));
	});
	return(r);
}

#
#	<p> omics data
#

quantileData = function(d, p) {
	dS = sort(d);
	q = dS[ceiling(p * length(dS))];
	q
}

quantileReference = function(reference, direction = 2, center = TRUE) {
	if (is.matrix(reference) && center) {
		refC =  matrixCenter(reference, direction);
		reference = matrixDeCenter(refC$matrix, mean(refC$center), direction);
	}
	ref = na.omit(as.vector(as.matrix(reference)));
	ref
}

Skewness = function(x, na.rm = T) {
	x0 = if(na.rm) na.omit(x) else x;
    N = length(x0);
    x1 = x0 - mean(x0)
    y = sqrt(N) * sum(x1^3) / (sum(x1^2)^(3/2))
    s = y * ((1 - 1/N))^(3/2);
    s
}
Noutliers = function(x, coef = 1.5)length(boxplot.stats(x, coef = coef)$out)

#' Quantile normalization of frame/matrix with respect to reference distribution
#'
#' Distribution to be normalized are represented as columns or rows of a matrix/data frame.
#' Each value is replaced by the quantile of the reference distribution as given by the value of the
#' empirical distribution function of the given value.
#'
#' @param reference numeric vector with realizations from the target distribution
#' @param data data frame or matrix with data to be normalized
#' @param direction is \code{data} organized per row or column?
#'
#' @examples
#' d = sapply(1:20, rnorm(1e4));
#' dNorm = quantileNormalization(as.vector(d), d)
quantileNormalization = function(reference, data, direction = 2,
	impute = TRUE, ties = 'random', center = TRUE, referenceDirection = direction) {
	ref = quantileReference(reference, referenceDirection, center);
	if (impute) mns = apply(data, 3 - direction, median, na.rm = T);
	dN = apply(data, direction, function(d) {
		d0 = d;
		if (impute) d[is.na(d0)] = mns[is.na(d0)];
		r = quantile(ref, probs = rank(d, na.last = 'keep', ties = ties)/length(na.omit(d)))
		if (impute) r[is.na(d0)] = NA;
		r
	});
	if (direction == 1) dN = t(dN);
	dimnames(dN) = dimnames(data);
	dN
}

# quantile normalization based on samples picked on the basis of their medians (around the medians)
# Nqn: number of reference samples
quantileNormalizationMedians = function(data, direction = 2, Nqn = 5, impute = TRUE) {
	# <p> determine median of medians, corresponding median, IQR
	medians = apply(data, direction, median);
	mediansO = order(medians);
	medianOI = as.integer(length(mediansO)/2 + .5);
	medianI = mediansO[medianOI];
	refMed = summary(data[, medianI]);
	refIQR = refMed[['3rd Qu.']] - refMed[['1st Qu.']];

	# <p> reference samples
	refL = as.integer(medianOI - Nqn/2 + .5);
	refU = refL + Nqn - 1;
	refSamples = mediansO[refL:refU];
	#print(refSamples)

	# <p> standardize reference samples wrt median, IQR
	refSampleValues = sapply(refSamples, function(i) {
		refI = summary(data[, i]);
		refIIQR = refI[['3rd Qu.']] - refI[['1st Qu.']];
		E = (data[, i] - refI[['Median']]) * refIQR/refIIQR + refMed[['Median']];
		#refIE = summary(E);
		#refIEIQR = refIE[['3rd Qu.']] - refIE[['1st Qu.']];
		#print(list(refI = refI, refIIQR = refIIQR, refIE = refIE, refIEIQR = refIEIQR));
		E
	});
	eQn = quantileNormalization(refSampleValues, data,
		direction = direction, impute = impute, center = FALSE);
	eQn
}

dataCentered = function(d, na.rm = T) {
	dC = apply(d, 2, function(col)col - mean(col, na.rm = na.rm));
	dC
}

#
#	<p> distributions
#

qtruncnorm = function(p, mean = 0, sd = 1, lower = -Inf, upper = Inf) {
	plow = pnorm(lower, mean, sd);
	pupp = pnorm(upper, mean, sd);
	qnorm(p * (pupp - plow) + plow, mean, sd = sd)
}
rtruncnorm = function(N, mean = 0, sd = 1, lower = -Inf, upper = Inf) {
	qtruncnorm(runif(N), mean, sd = sd, lower, upper)
}

qqDist = function(Nqts = 1e2, qdist, ...) {
	qtls = (1:Nqtls)/(Nqtls + 1);
	qtlsExp = qdist(qtls, ...);
	qtlsObs = quantile(rtn, qtls);
	qq = qplot(qtlsExp, qtlsObs) + theme_bw();
	qq
}

qqSim = function(Nsim, dist = 'truncnorm', Nqts = Nsim/10, ...) {
	require('ggplot2');
	rdist = get(Sprintf('r%{dist}s'));
	r = rdist(Nsim, ...);
	qdist = get(Sprintf('q%{dist}s'));
	qq = qqDist(Nqts, qdist, ...);
	qq
}

#
#	<p> entropy
#

table.entropy = function(d) {
	p = table.freq(d);
	p = p[p != 0];
	H = - sum(p * log(p));
}

#
#	<p> qvalue
#

Qvalue = function(P.value, ...) {
	Require('qvalue');
	P.valuesNotNA = na.omit(P.value);
	qv = qvalue(P.valuesNotNA, ...);
	r = qv;
	r$qvalue = vector.embed(rep(NA, sum(is.na(P.value))), which(!is.na(P.value)), qv$qvalue);
	r
}

#
#	<p> math functions
#

Ceiling = function(x, N = 0)(ceiling(x * 10^N) / 10^N)
Floor = function(x, N = 0)(floor(x * 10^N) / 10^N)

#
#	<p> descriptive statistics
#

SummaryColNum = function(col, qtls = c(0, 0.25, .5, .75, 1)) {
	qs = quantile(col, qtls, na.rm = T);
	mn = vector.named(mean(col, na.rm = T), 'mean');
	na = vector.named(fraction(is.na(col)), 'NA');
	Sd = vector.named(sd(col, na.rm = T), 'sd');
	medianI = which(qtls == .5);
	qs1 = if (length(medianI) > 0) vector.embed(qs, medianI + 1, mn) else c(qs, mean);
	c(qs1, Sd, na)
}
SummaryColFactor = function(col) {
	na = vector.named(fraction(is.na(col)), 'NA');
	na
}
SummaryCol = function(col, qtls = c(0, 0.25, .5, .75, 1)) {
	r = if (class(col) == 'numeric') SummaryColNum(col, qtls) else
		if (class(col) == 'factor') c(rep(NA, length(qtls) + 2), SummaryColFactor(col)) else
			stop("Could not summarize column");
	r
}

Summary = function(d) {
	r0 = lapply(d, SummaryCol)
	r = t(as.matrix(do.call(cbind, r0)));
	r
}


descriptivesNumeric = function(data) {
	bp = apply(data, 2, boxplot.stats);
	N = list.kp(bp, 'n', do.unlist = T);
	r = data.frame(
		N = N,
		Nna = nrow(data) - N,
		mean = apply(data, 2, mean, na.rm = T),
		median = apply(data, 2, median, na.rm = T),
		se = apply(data, 2, sd, na.rm = T) / sqrt(N),
		Nout = sapply(bp, function(v)length(v$out))
	);
	r
}
descriptivesFactor = function(data) {
	Nnas = apply(data, 2, function(c)sum(is.na(c)));
	r = data.frame(
		N = nrow(data) - Nnas,
		Nna = Nnas,
		Ncat = sapply(lapply(data, levels), length)
	);
	r
}

descriptives = function(data) {
	r1 = descriptivesFactor(data[, DfClasses(data) == 'factor']);
	r2 = descriptivesNumeric(data[, DfClasses(data) == 'numeric']);
	r0 = if (nrow(r2) == 0) list(r1) else if (nrow(r1) == 0) list(r2) else list(r1, r2);
	r = listOfDataFrames2data.frame(r0, idColumn = NULL, colsFromUnion = T, row.names = T);
	r$`NAperc` = with(r, Nna / nrow(data) * 100);
	r
}


#
#	<p> survival
#

# adapted from survival:::print.survfit
survfit2m = function (x, scale = 1, ..., rmean = 'none') {
	if (inherits(x, "survfitms")) {
		x$surv <- 1 - x$prev
		if (is.matrix(x$surv)) dimnames(x$surv) <- list(NULL, x$states)
		if (!is.null(x$lower)) {
			x$lower <- 1 - x$lower
			x$upper <- 1 - x$upper
		}
		rmean <- "none"
	}
	omit <- x$na.action
	temp <- survival:::survmean(x, scale = scale, rmean)
	mtemp <- if (is.matrix(temp$matrix)) 
		temp$matrix
	else matrix(temp$matrix, nrow = 1, dimnames = list(NULL, 
		names(temp$matrix)))
	if (all(mtemp[, 2] == mtemp[, 3])) {
		cname <- dimnames(mtemp)[[2]]
		mtemp <- mtemp[, -2, drop = FALSE]
		cname <- cname[-2]
		cname[2] <- "n"
		dimnames(mtemp)[[2]] <- cname
	}
	if (all(mtemp[, 1] == mtemp[, 2])) 
		mtemp <- mtemp[, -1, drop = FALSE]
	temp$matrix <- drop(mtemp)
	temp$matrix
}

#
#	<p> mixed models
#

lmer2ci = function(object, method = 'profile')confint.merMod(object, method = method);

lmer2p = function(object, method = 'profile') {
	cis = lmer2ci(object, method = method);
	apply(cis, 1, function(ci)ciToP(ci[1], ci[2]))
}

polr2p = function(rPoly) {
 	ci = t2r(confint(rPoly));
	cfs = coefficients(summary(rPoly));	
	cfs1 = cbind(cfs, P = c(ciToP(ci[, 1], ci[, 2]), rep(NA, nrow(cfs) - nrow(ci))));
	return(cfs1);
}

# confint plus P-value
confintP = function(m, level = .95) {
	ci = confint(m, level = level);
	P = apply(ci, 1, function(ci)ciToP(ci[1], ci[2], level = level));
	cbind(ci, P)
}

coefficientsLmer = function(object, method = 'Wald') {
	cfs = coefficients(summary(object));
	P = lmer2p(object, method);
	Iinter = which(names(P) == '(Intercept)');
	r = cbind(cfs, P = P[Iinter:length(P)]);
	return(r);
}


# alternative check for convergence of glmer models
#https://stackoverflow.com/questions/21344555/convergence-error-for-development-version-of-lme4
glmerFit = function(m) {
	relgrad = with(m@optinfo$derivs, solve(Hessian, gradient));
	return(max(abs(relgrad)));
}
GlmerFit = function(m, cutoff = 1e-3)(glmerFit(m) < cutoff)

#
#	<p> simulation studies
#

simulationRun = dataSimulation = function(Niteration, fSim, fStat = function(data, ...)summary(data),
	..., parSim = list(), parStat = list()) {
	parShared = list(...);
	r = Lapply(1:Niteration, function(i) {
		data = do.call(fSim, c(parSim, parShared));
		r = do.call(fStat, c(list(data = data), parStat, parShared));
	});
}

simulationRunBySpec = function(fSim, fStat = function(data, ...)summary(data),
	spec = list(Nsim = 1e1, sim = list(), stat = list(), pars = list())) {

	args = c(
		list(Niteration = spec$Nsim, fSim = fSim, fStat = fStat,
			parSim = spec$sim, parStat = spec$stat),
		spec$pars
	);
	r = do.call(simulationRun, args);
	return(r);
}

simulationDataBySpec = function(fSim,
	spec = list(Nsim = 1e1L, sim = list(), stat = list(), pars = list())) {
	data = do.call(fSim, c(spec$pars, spec$sim));
	return(data);
}
simulationStatBySpec = function(fTest, data,
	spec = list(Nsim = 1e1L, sim = list(), stat = list(), pars = list())) {
	r = do.call(fTest, c(list(data = data), spec$pars, spec$stat));
	return(r);
}

#
#	<p> delta rule
#

#
#	Numeric delta rule
#

gradient = function(f, x, ..., eps = 1e-5) {
	sapply(seq_along(x), function(i) {
		epsI = vector.assign(0, i, eps, N = length(x));
		(f(x + epsI, ...) - f(x - epsI, ...)) / (2 * eps)
	});
}

totalDeriv = function(f, x, ..., eps = 1e-5) {
	y = f(x, ...);
	t(sapply(seq_along(y), function(i)gradient(function(x)f(x, ...)[i], x, eps = eps)));
}

# delta rule for functions of binomial data
#N: list of pairs Nevents, Ntotal, or nx2 matrix of these numbers
#f: function to be computed on the frequencies
deltaFreqs = function(N, f, logScale = T, ..., eps = 1e-5, alpha = .95, pseudo = 0) {
	#if (is.list(N)) N = sapply(N, identity);
	if (is.list(N)) N = do.call(cbind, N);
	if (pseudo != 0) {
		N[1, ] = N[1, ] + pseudo;
		N[2, ] = N[2, ] + 2*pseudo;
	}
	freqs = N[1, ] / N[2, ];
	var = freqs * (1 - freqs) / N[2, ];

	fDeriv = if (logScale) function(x, ...)log(f(x, ...)) else f;
	fGradient = gradient(fDeriv, freqs, ..., eps = eps);

	varf = (t(fGradient) %*% diag(var) %*% fGradient)[1, 1];
	limit = qnorm(1 - (1 - alpha)/2, 0, sqrt(varf));

	y = fDeriv(freqs, ...);
	yCi = y + c(-limit, limit);

	r = if (logScale)
		list(y = exp(y), ci = exp(yCi), var = varf, se = sqrt(varf), logCi = yCi) else
		list(y = y, ci = yCi, var = varf, se = sqrt(varf));
	r
}

delta = function(x, Cov, f, logScale = F, ..., eps = 1e-5, alpha = .95) {
	fDeriv = if (logScale) function(x, ...)log(f(x, ...)) else f;
	fGradient = gradient(fDeriv, x, ..., eps = eps);

	varf = (t(fGradient) %*% Cov %*% fGradient)[1, 1];
	y = fDeriv(x, ...);

	limit = qnorm(1 - (1 - alpha)/2, 0, sqrt(varf));
	yCi = y + c(-limit, limit);

	r = if (logScale)
		list(y = exp(y), ci = exp(yCi), var = varf, se = sqrt(varf), logCi = yCi) else
		list(y = y, ci = yCi, var = varf, se = sqrt(varf));
	r
}

# return freq-vector, cov matrix
deltaMultinomFromCnts = function(counts) {
	N = sum(counts);
	freqs = counts / N;
	Var = freqs * (1 - freqs);
	Cov = freqs %*% t(freqs);
	diag(Cov) = Var;
	return(list(par = freqs, cov = Cov));
}

# compose lists of par/cov structures
deltaCompose = function(parCov) {
	return(list(par = list.kpu(parCov, 'par'), cov = .bdiag(list.kp(parCov, 'cov'))));
}

#
#	<p> sensitvity/specificity/ROC/AUC
#


# assume TRUE/FALSE
binaryMeasure = function(tab, strict = TRUE) {
	if (strict && !all(names(dimnames(tab)) == c('pred', 'truth'))) stop('dimnames != pred, truth');
	if (strict && !all(unique(unlist(dimnames(tab))) == c('FALSE', 'TRUE'))) stop('values != FALSE, TRUE');
	r = list(
		tp = sum(tab[, 'TRUE']) / sum(tab),
		tn = sum(tab[, 'FALSE']) / sum(tab),
		sensitivity = tab['TRUE', 'TRUE'] / sum(tab[, 'TRUE']),
		specificity = tab['FALSE', 'FALSE'] / sum(tab[, 'FALSE']),
		ppv = tab['TRUE', 'TRUE'] / sum(tab['TRUE', ]),
		npv = tab['FALSE', 'FALSE'] / sum(tab['FALSE', ])
	);
	return(r);
}

binaryMeasures = function(prob, labels) {
	cutoffs = sort(unique(prob));
	measures = lapply(cutoffs, function(co) {
		tab = Table(list(pred = prob >= co, truth = labels), cats = list(pred = c(F, T), truth = c(F, T)));
		bm = c(list(cutoff = co), binaryMeasure(tab));
		bm
	});
	r = apply(do.call(rbind, measures), 2, unlist);
	r
}

cvROC = function(f1, data,
	cv_fold = 10, cv_repeats = 1, cv_train = cv_train_glmnet, cv_test = cv_test_glmnet) {
	Library('pROC');
	Library('AUC');

	# <p> data preparation
	mm = model_matrix_from_formula(f1, data);
	d0 = data[mm$indeces, ];
	d0resp = d0[[formula.response(f1)]];

	# <p> crossvalidation
	pred = crossvalidate(cv_train, cv_test, data = d0, cv_fold = cv_fold, cv_repeats = cv_repeats, f1 = f1);

	# <p> results
	# <i><N> only first repeat headed at the moment
	roc = AUC::roc(pred[[1]], as.factor(d0resp));
	auc = AUC::auc(roc);
	rocProc = pROC::roc(as.factor(d0resp), pred[[1]]);
	aucProc = pROC::auc(rocProc);
	P = cor.test(d0resp, pred[[1]])

	r = list(roc = roc, auc = auc, rocProc = rocProc, aucProc = aucProc, aucCi = ci(rocProc), P = P,
		data = Df(response = d0resp, prediction = pred[[1]]));
	return(r);
}

plotRocs = function(data, interval = 0.2, breaks = seq(0, 1, interval)) {
	d0 = do.call(rbind, nlapply(data, function(n) {
		d = if (all( c('tpr', 'fpr') %in% names(data[[n]]))) data[[n]] else	# roc already computed
			with(data[[n]], AUC::roc(prediction, as.factor(response)));
		Df(fpr = 1 - d$fpr, tpr = d$tpr, type = n)
	}));
	p = ggplot(data = d0) +
		#geom_line(aes(x = fpr, y = tpr)) +
		geom_step(aes(x = fpr, y = tpr, col = type), lwd = 1) +
		geom_segment(aes(x = 0, y = 1, xend = 1, yend = 0)) +	#diagonal
		scale_x_reverse(name = "Specificity", limits = c(1,0), breaks = breaks,  expand = c(0.001,0.001)) + 
		scale_y_continuous(name = "Sensitivity", limits = c(0,1), breaks = breaks, expand = c(0.001, 0.001))+
		labs(color = 'ROC type') +
		coord_equal() + theme_bw();
	p
}


#
#	<p> multiple testing
#
tailStrength = function(ps, M = length(ps)) {
	N = length(ps);
	# standardize p-values to 0
	psS = 1 - sort(ps) * (N + 1)/(1:N);
	r = - sum(psS[1:M])/M;	# return negative as small values are deemed relevant
	r
}

#
#	<p> power calculations
#

findBinomN = function(pTrue, pLower, level = 0.95, Nmax = 1e3, lower = T) {
	Library('binom');
	fBinom = if (lower)
		function(N)binom.logit(N*pTrue, N, conf.level = level)$lower else
		function(N)binom.logit(N*pTrue, N, conf.level = level)$upper
	inverse(fBinom, c(1, Nmax))(pLower);
}

#
#	<p> pseudo counts
#

# implement pseudo-count, needed as observations c(0, 0, X) cause diverging estimates
mulitnomCoeff = function(x, log = TRUE) {
	N = sum(x);
	muco = lgamma(N + 1) - sum(sapply(x + 1, lgamma));
	if (!log) muco = exp(muco);
	return(muco);
}

dmultinomPseudo = function(x, prob, log = FALSE, pseudo = .1, factors = 1:1e2, F = NULL) {
	# find multiplication factor
	if (is.null(F)) {
		pseudoF = t(t_(pseudo) %*% factors);
		pseudoR = abs(pseudoF - round(pseudoF));
		F = factors[apply(pseudoR, 2, which.min)];
	}
	# standardize probabilities
	logP = log(v2freq(prob));
	# raw likelihood contributions
	pRaw = sum(round((x + pseudo) * F) * logP / F);
	# add multinomial coefficient
	pCooked = pRaw + mulitnomCoeff(x);
	if (!log) pCooked = exp(pCooked);
	return(pCooked);
}

#
#	<p> random numbers
#

generateSeed = function(preSeed = Sys.time()) {
	ints = hex2ints(md5sumString(as.integer(preSeed)));
	r = 0
	for (i in 1:length(ints)) r = bitwXor(r, ints[i]);
	return(r);
}

# run expression with fixed seed but leave the RNG stream unchanged
#	(bracket expression with save/restore of seed)
#	as.integer(Sys.time())
fixSeedForExpr = function(seed, expr, envir = parent.frame()) {
	if (!exists('.Random.seed')) runif(1);	# vivify RNG
	currentSeed = .Random.seed;
	on.exit(.Random.seed <<- currentSeed);
	set.seed(seed);
	r = eval(expr, envir = envir);
	#.Random.seed <<- currentSeed;
	return(r);
}

#
#	<p> general distributions
#

lchoosemult = function(x) {
	N = sum(x);
	K = length(x);
	lmult = dmultinom(x, prob = rep(1, K), log = T) + N*log(K);
	return(lmult);
}

#
#	<p> sampling
#

# individual probabilities, sample single event
# P: N x (K - 1) matrix with class probabilities
# size == 1 only supported case
rmultinomial = function(P, size = 1) {
	Ps = t2r(apply(P, 1, cumsum));
	R = runif(nrow(P));	# case of size == 1
	y = apply(apply(Ps, 1, function(P)R > P), 1, function(i)max(c(0, which(i))));
	return(y);
}

rbinomial = function(P, Nbinom = 1) {
	# binomial case for Nbinom >= 1
	P = t_(sapply(P, function(p)cumsum(dbinom(0:(Nbinom - 1), Nbinom, p))));
	return(rmultinomial(P));
}
