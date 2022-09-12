#
#	Rdata.R
#Mon 27 Jun 2005 10:49:06 AM CEST
#system("~/src/Rprivate/exportR.sh");
#system("~/src/Rprivate/exportR.sh"); source("RgenericAll.R"); source("Rgenetics.R"); loadLibraries();
#system('~/src/Rprivate/exportR.sh ; cp ~/src/Rprivate/RgenericAllRaw.R .');
# <!> copied to Rmeta.R as being the first file to be exported by now (26.3.2017)

#
#	<§> abstract data functions
#

defined = function(x) exists(as.character(substitute(x)));
defined.by.name = function(name) { class(try(get(name), silent = TRUE)) != 'try-error' }
# equivalent to i %in% v
is.in = function(i, v)(length((1:length(v))[v == i])>0)
rget = function(name, default = NULL, ..., pos = -1, envir = as.environment(pos)) {
	#obj = try(get(name, ...), silent = TRUE);
	#r = if(class(obj) == 'try-error') default else obj;
	#r = if (exists(name, where = pos, envir = envir)) get(name, ..., pos = pos, envir = envir) else default;
	r = if (exists(name, envir = envir)) get(name, ..., envir = envir) else default;
	r
}
# .fdE: use notE
firstDef = function(..., .fdInterpolate = FALSE, .fdIgnoreErrors = FALSE, .fdE = FALSE) {
	l = if (.fdInterpolate) c(...) else list(...);
	for (i in l) {
		if ((!is.null(i) && (!.fdE || notE(i))) && (!.fdIgnoreErrors || class(i) != 'try-error'))
			return(i)
	};
	NULL
}
FirstDef = function(..., .fdInterpolate = FALSE, .fdIgnoreErrors = FALSE, .fdE = TRUE)
	firstDef(..., .fdInterpolate = .fdInterpolate, .fdIgnoreErrors = .fdIgnoreErrors, .fdE = .fdE)
firstDefNA = function(..., .fdInterpolate = FALSE) {
	l = if (.fdInterpolate) c(...) else list(...);
	for (i in l) { if (!is.na(i)) return(i)};
	NULL
}
# <N> NULL behaviour
to.list = function(..., .remove.factors = TRUE){
	r = if(is.null(...)) NULL else if (is.list(...)) c(...) else list(...);
	if (.remove.factors) {
		r = sapply(r, function(e)ifelse(is.factor(e), levels(e)[e], e));
	}
	r
}
# clean list/vector
Avu = function(v)as.vector(unlist(v))
is.Null = nullOrLengthNull = function(e)(is.null(e) || length(e) == 0);
# pretty much force everything to be a vector
# toNA: length(0) or NULL converted to NA
avu = function(v, recursive = TRUE, toNA = TRUE) {
	transform = if (toNA)
		function(e, condition)(if (condition) NA else avu(e, toNA = TRUE, recursive = TRUE)) else
		function(e, ...)avu(e, toNA = FALSE, recursive = TRUE);

	r = if (is.list(v)) {
		nls = sapply(v, is.Null);	# detects nulls
		# unlist removes NULL values -> NA
		unlist(sapply(seq_along(v), function(i)transform(v[[i]], nls[i])));
	} else as.vector(v);
	if (!length(r)) return(NULL);
	r
}
#pop = function(v)rev(rev(v)[-1]);

assign.list = function(l, pos = -1, envir = as.environment(pos), inherits = FALSE, immediate = TRUE) {
	for (n in names(l)) {
		assign(n, l[[n]], pos, envir, inherits, immediate);
	}
}
eval.text = function(text, envir = parent.frame())eval(parse(text = text), envir= envir);

nullomit = function(r)r[!sapply(r, is.null)]

# replace elements base on list
# l may be a list of lists with elements f (from) and t (to), when f is replaced with t
# if both, f and t arguments are not NULL, l will be ignored and f is replaced with t
vector.replace = function(v, l, regex = FALSE, ..., f = NULL, t = NULL) {
# 	if (!is.null(f) & !is.null(t)) l = list(list(f = f, t = t));
# 	# replacments are given in f/t pairs
# 	if (all(sapply(l, length) == 2)) {
# 		from = list.key(l, "f");
# 		to = list.key(l, "t");
# 	} else {
# 		from = names(l);
# 		to = unlist(l);
# 	}
# 	for (i in 1:length(from)) {
# 		if (regex) {
# 			idcs = which(sapply(v, function(e)(length(fetchRegexpr(from[i], e, ...)) > 0)));
# 			v[idcs] = sapply(v[idcs], function(e)gsub(from[i], to[i], e));
# 		} else v[which(v == from[i])] = to[i];
# 	}
 	repl = if (!is.null(f) & !is.null(t)) listKeyValue(f, t) else l;
	# <!> tb tested
 	v = if (!regex) {
		raw = repl[v];
		unlist(ifelse(sapply(repl[v], is.null), v, raw))
	} else {
		sapply(v, function(e){
			# first match takes precedent
			j = which(sapply(names(repl), function(f)length(fetchRegexpr(f, e, ...)) > 0))[1];
			if (is.na(j)) e else gsub(names(repl)[j], repl[[j]], e)
		})
	}
	v
}

vector.with.names = function(v, all_names, default = 0) {
	r = rep(default, length(all_names));
	names(r) = all_names;
	is = which.indeces(names(v), all_names, ret.na = TRUE);
	r[is[!is.na(is)]] = v[!is.na(is)];
	r
}
vector.named = function(v, Names) {
	names(v) = Names;
	v
}

# dir: direction of selection: 1: select rows, 2: select columns
mat.sel = function(m, v, dir = 1) {
	r = if (dir == 1)
		sapply(1:length(v), function(i)m[v[i], i]) else
		sapply(1:length(v), function(i)m[i, v[i]]);
	r
}

# rbind on list
simplify = sapplyId = function(l)sapply(l, identity);
Simplify = function(l)unlist(simplify(l));

listFind = function(lsed, lsee) {
	values = sapply(names(lsee), function(n)list.key(lsed, n), simplify = FALSE, USE.NAMES = FALSE);
	values = sapply(values, identity);
	found = apply(values, 1, function(r) all(r == lsee));
	r = unlist.n(lsed[found], 1);
	r
}

#same.vector = function(v)(unique(v) == 1)
same.vector = function(v)all(v == v[1])

# in vector v, find index min j \in 1, ..., N so that v[1:j] contains at least U unique elements
uniqueIndex = function(v, U) {
	#Nu = sapply(seq_along(v), function(i)length(unique(data$chr[1:i])));
	# more efficient version
	u = c();
	for (i in seq_along(v)) {
		u = unique(c(u, v[i]));
		if (length(u) == U) return(i);
	}
	return(NA);
}


#
#	<§> string manipulation
#

#join = function(v, sep = " ")if (length(v) == 0) '' else paste(v, collapse = sep);
join = function(v, sep = " ")paste(v, collapse = sep);
con = function(..., Sep_ = '')paste(..., sep = Sep_);
Con = function(..., Sep_ = '')paste(unlist(list(...)), collapse = Sep_);
# pastem = function(a, b, ..., revsort = TRUE) {
# 	if (revsort)
# 		as.vector(apply(merge(data.frame(a = b), data.frame(b = a), sort = FALSE), 1,
# 			function(e)paste(e[2], e[1], ...))) else
# 		as.vector(apply(merge(data.frame(a = a), data.frame(b = b), sort = FALSE), 1,
# 			function(e)paste(e[1], e[2], ...)))
# }
pastem = function(a, b, ..., revsort = TRUE) {
	df = merge.multi.list(list(Df(a = a), Df(b = b)), .first.constant = revsort);
	paste(df[, 1], df[, 2], ...)
}

r.output.to.vector.int = function(s) {
	matches = gregexpr("(?<![\\[\\d])\\d+", s, perl=TRUE);
	starts = as.vector(matches[[1]]);
	lengthes = attr(matches[[1]], "match.length");
	v = sapply(1:length(starts), function(i){ substr(s, starts[i], starts[i] + lengthes[i] -1) });
	as.integer(v)
}
r.output.to.vector.numeric = function(s) {
	matches = gregexpr("\\d*\\.\\d+", s, perl=TRUE);
	starts = as.vector(matches[[1]]);
	lengthes = attr(matches[[1]], "match.length");
	v = sapply(1:length(starts), function(i){ substr(s, starts[i], starts[i] + lengthes[i] -1) });
	as.numeric(v)
}
readFile = function(path) { join(scan(path, what = "raw", sep = "\n", quiet = TRUE), sep = "\n") };
circumfix = function(s, post = NULL, pre = NULL) {
	if (is.null(s) || length(s) == 0) return('');
	sapply(s, function(s)if (s == '') s else con(pre, s, post))
}
abbr = function(s, Nchar = 20, ellipsis = '...') {
	ifelse(nchar(s) > Nchar, paste(substr(s, 1, Nchar - nchar(ellipsis)), ellipsis, sep = ''), s)
}
wrapStr = function(s, Nchar = 60, regex = '\\s+', indent = "\n") {
	r = '';
	while (nchar(s) > Nchar) {
		R = gregexpr('\\s+', s, perl = TRUE);
		Iws = R[[1]][R[[1]] <= Nchar];
		Ichr = max(Iws);
		# <i> handle Ichr = 1
		r = con(r, substr(s, 1, Ichr - 1), indent);
		s = substr(s, Ichr + attr(R[[1]], 'match.length')[length(Iws)], nchar(s));

	}
	r = con(r, s);
	return(r);
}

Which.max = function(l, last.max = TRUE, default = NA) {
	if (is.logical(l) && all(!l)) return(default);
	r = if (last.max) (length(l) - which.max(rev(l)) + 1) else which.max(l);
	r
}
Which.min = function(l, last.min = FALSE, default = NA) {
	if (is.logical(l) && all(!l)) return(default);
	r = if (last.min) (length(l) - which.min(rev(l)) + 1) else which.min(l);
	r
}
# capturesN: named captures; for each name in captureN put the captured value assuming names to be ordered
# captures: fetch only first capture per match <!> deprecated
# capturesAll: fetch all caputers for each match
fetchRegexpr = function(re, str, ..., ret.all = FALSE, globally = TRUE, captures = FALSE, captureN = c(),
	capturesAll = FALSE, maxCaptures = 9, returnMatchPositions = FALSE) {
	if (length(re) == 0) return(c());
	r = if (globally)
		gregexpr(re, str, perl = TRUE, ...)[[1]] else
		regexpr(re, str, perl = TRUE, ...);
	if (all(r < 0)) return(NULL);
	l = sapply(1:length(r), function(i)substr(str, r[i], r[i] + attr(r, "match.length")[i] - 1));
	if (captures) {
		l = sapply(l, function(e)gsub(re, '\\1', e, perl = TRUE, fixed = FALSE));
	} else if (length(captureN) > 0) {
		l = lapply(l, function(e) {
			r = sapply(1:length(captureN), function(i) {
				list(gsub(re, sprintf('\\%d', i), e, perl = TRUE, fixed = FALSE))
			});
			names(r) = captureN;
			r
		});
	} else if (capturesAll) {
		l = lapply(l, function(e) {
			cs = c();	# captures
			# <!> hack to remove zero-width assertions (no nested grouping!)
			#re = gsub('(\\(\\?<=.*?\\))|(\\(\\?=.*?\\))', '', re, perl = TRUE, fixed = FALSE);
			for (i in 1:maxCaptures) {
				n = gsub(re, sprintf('\\%d', i), e, perl = TRUE, fixed = FALSE);
				cs = c(cs, n);
			}
			cs
		});

		# trim list
		#maxEls = maxCaptures - min(c(maxCaptures + 1, sapply(l, function(e)Which.max(rev(e != ''))))
		#	, na.rm = TRUE) + 1;
		maxEls = max(c(sapply(l, function(e)Which.max(e != '', default = 1)), 1));
		l = lapply(l, function(e)(if (maxEls > 0) e[1:maxEls] else NULL));
	}
	if (!ret.all) l = l[l != ""];
	ret = if (returnMatchPositions) list(match = l, positions = r) else l;
	ret
}
# improved multistring version
FetchRegexpr = function(re, str, ..., ret.all = FALSE, globally = TRUE, captures = FALSE, captureN = c(),
	capturesAll = FALSE, maxCaptures = 9, returnMatchPositions = FALSE) {
	if (length(re) == 0) return(c());
	r = if (globally)
		gregexpr(re, str, perl = TRUE, ...) else
		list(regexpr(re, str, perl = TRUE, ...));
	if (all(unlist(r) < 0)) return(NULL);
	l = sapply(seq_along(r),
		function(j) {
			r0 = r[[j]];
			sapply(1:length(r0),
				function(i)substr(str[j], r0[i], r0[i] + attr(r0, "match.length")[i] - 1))
	});
	if (captures) {
		l = sapply(l, function(e)gsub(re, '\\1', e, perl = TRUE, fixed = FALSE));
		#print(l);
	} else if (length(captureN) > 0) {
		l = lapply(l, function(e) {
			r = sapply(1:length(captureN), function(i) {
				list(gsub(re, sprintf('\\%d', i), e, perl = TRUE, fixed = FALSE))
			});
			names(r) = captureN;
			r
		});
	} else if (capturesAll) {
		l = lapply(l, function(e) {
			cs = c();	# captures
			# <!> hack to remove zero-width assertions (no nested grouping!)
			#re = gsub('(\\(\\?<=.*?\\))|(\\(\\?=.*?\\))', '', re, perl = TRUE, fixed = FALSE);
			for (i in 1:maxCaptures) {
				n = gsub(re, sprintf('\\%d', i), e, perl = TRUE, fixed = FALSE);
				cs = c(cs, n);
			}
			cs
		});

		# trim list
		#maxEls = maxCaptures - min(c(maxCaptures + 1, sapply(l, function(e)Which.max(rev(e != ''))))
		#	, na.rm = TRUE) + 1;
		maxEls = max(c(sapply(l, function(e)Which.max(e != '', default = 1)), 1));
		l = lapply(l, function(e)(if (maxEls > 0) e[1:maxEls] else NULL));
	}
	if (!ret.all) l = l[l != ""];
	ret = if (returnMatchPositions) list(match = l, positions = r) else l;
	ret
}

regex = Vectorize(fetchRegexpr, 'str', SIMPLIFY = TRUE, USE.NAMES = TRUE);
Regex = Vectorize(FetchRegexpr, 're', SIMPLIFY = TRUE, USE.NAMES = TRUE);
RegexL = Vectorize(FetchRegexpr, 're', SIMPLIFY = FALSE, USE.NAMES = TRUE);
regexIdcs = function(re, s, ...)vectorIdcs(regex(re, s, ...), is.null, not = TRUE)

# unify capture extraction for gregexpr, regexpr
# pos == 0: grexepr, regexpr else by iterating pos as index into str
matchRegexCapture = function(reg, str, pos = NULL) {
	if (is.null(attr(reg, 'capture.start'))) return(NULL);
	if (!is.null(pos)) str = str[pos] else pos = seq_along(reg);
	captures = lapply(1:ncol(attr(reg, 'capture.start')), function(i) {
		vs = sapply(pos, function(j)Substr(str,
			attr(reg, 'capture.start')[j, i], attr(reg, 'capture.length')[j, i]))
		vs
	});
	names(captures) = attr(reg, 'capture.names');
	captures
}
MatchRegexExtract = function(m, s, pos = seq_along(m)) {
	matches = ifelse(m[pos] < 0, character(0),
		sapply(pos, function(i)Substr(s[i], m[i], attr(m, 'match.length')[i])));
	matches
}
matchRegexExtract = function(reg, str, pos = NULL) {
	if (!is.null(pos)) str = str[pos] else pos = seq_along(reg);
	matches = ifelse(reg[pos] < 0, character(0),
		sapply(pos, function(i)Substr(str, reg[i], attr(reg, 'match.length')[i])));
	matches
}
# <i> re nested list with sub-res for named captures
# <!> globally == FALSE, removeNonMatch == FALSE
matchRegex = function(re, str, ..., globally = TRUE, simplify = TRUE,
	positions = FALSE, removeNonMatch = FALSE) {
	if (length(re) == 0) return(NULL);
	reg = if (globally) gregexpr(re, str, perl = TRUE, ...) else regexpr(re, str, perl = TRUE, ...);
	ms = if (globally)
		lapply(seq_along(reg), function(i)matchRegexExtract(reg[[i]], str[i])) else
		lapply(seq_along(str), function(i)matchRegexExtract(reg, str, pos = i));
	#	regmatches(str, reg);
	captures = if (globally)
		lapply(seq_along(reg), function(i)matchRegexCapture(reg[[i]], str[i])) else
		lapply(seq_along(str), function(i)matchRegexCapture(reg, str, pos = i));
	if (removeNonMatch) {
		nonmatch = sapply(ms, length) == 0 | is.na(ms);
		ms = ms[!nonmatch];
		captures = captures[!nonmatch];
		reg = reg[!nonmatch];
	}
	if (simplify && length(str) == 1) {
		ms = ms[[1]];
		captures = captures[[1]];
		reg = reg[[1]];
	}
	r = if(positions) list(match = ms, capture = captures, positions = reg) else
		list(match = ms, capture = captures);
	r
}

#
#	<p> final interface as of 2016/04
#
MatchRegex = function(re, str, mode = 'return') {
	r = regexpr(re, str);
	if (mode == 'return') {
		r = str[which(r > 0)];
	}
	r
}
# handle attributes
# As.list assumes attributes and vector elements to be paired
#	corresponding values/attributes will be put into the list
As.list = function(v) {
	as = Recycle(attributes(v));
	l = lapply(seq_along(v), function(i) {
		# does not preserve matrices
		#attrs = list.kp(as, Sprintf('[[%{i}d]]'));
		# should become <i>
		#attrs = list.kp(as, Sprintf('[[%{i}d]]', accessor = function(e)accessIdx(e, i)));
		attrs = lapply(seq_along(as), function(j)accessIdx(as[[j]], i));
		Attr(v[i], SetNames(attrs, names(as)))
	});
	l
}

# transform results from Regexpr captures = TRUE
list.transpose = function(l)lapply(seq_along(l[[1]]), function(i)list.kp(l, Sprintf('[[%{i}d]]')));

# interface as of 2018/06
# if re is vector, iterate over
# by default, return matches
RegexprSingle = function(re, s, captures = FALSE, global = TRUE, simplify = TRUE, concatMatches = TRUE, drop = TRUE) {
	matches = if (global) gregexpr(re, s, perl = TRUE) else As.list(regexpr(re, s, perl = TRUE));
	#print(gregexpr(re, s, perl = TRUE));
	#print(regexpr(re, s, perl = TRUE));
	#print(matches);

	#matches = if (global) gregexpr(re, s, perl = TRUE) else list(regexpr(re, s, perl = TRUE));
	r = pairslapply(matches, s, function(m, s) {	# iterate strings
		if (captures) {
			r = matchRegexCapture(m, s);
			if (concatMatches) r = apply(do.call(cbind, r), 1, join, sep = '');
		} else {
			r = MatchRegexExtract(m, s);
			if (drop) r = r[!is.na(r)];
		}
		r
	});
	if (simplify && (
		(length(s) == 1 && captures && concatMatches)
	)) r = r[[1]];
	if (simplify && !global) r = Simplify(r);
	return(r);
}

Regexpr = function(re, s, ..., reSimplify = TRUE) {
	r = lapply(re, RegexprSingle, s = unlist(s), ...);
	if (length(re) == 1 && reSimplify) r = r[[1]];
	return(r);
}
RegexprM = function(re, s, ..., reSimplify = TRUE) {
	r = sapply(Regexpr(re, s, ..., reSimplify = reSimplify), function(e)length(e) > 0);
	return(r);
}

splitString = function(re, str, ..., simplify = TRUE) {
	l = lapply(str, function(str) {
		if (is.na(str)) return(NA);
		r = gregexpr(re, str, perl = TRUE, ...)[[1]];
		if (r[1] < 0) return(str);
		l = sapply(1:(length(r) + 1), function(i) {
			substr(str, ifelse(i == 1, 1, r[i - 1] + attr(r, "match.length")[i - 1]),
				ifelse(i > length(r), nchar(str), r[i] - 1))
		});
	});
	if (length(l) == 1 && simplify) l = l[[1]];
	l
}
# modeled after perl's qq
reString = '(?:([_\\/\\-a-zA-Z0-9.]+)|(?:\\"((?:\\\\\\\\.)*(?:[^"\\\\]+(?:\\\\\\\\.)*)*)\\"))';
# use reSep = '\\s+' to split based on a separator RE
qw = function(s, re = reString, reSep = NULL, names = NULL, byrow = TRUE) {
	r = if (notE(reSep)) unlist(splitString(reSep, s)) else {
	#r = if (TRUE) unlist(splitString('\\s+', s)) else
		unlist(Regexpr(re, unlist(s), captures = TRUE));
	}
	if (notE(names)) r = Df_(matrix(r, ncol = length(names), byrow = byrow), names = names);
	r
}
quoteString = function(s)sprintf('"%s"', s)
trimString = function(s) {
	sapply(s, function(e)
		if (is.na(e)) NA else FetchRegexpr('^\\s*(.*?)\\s*$', e, captures = TRUE)
	)
}
qwi = function(...)as.integer(qw(...))
qwn = function(...)as.numeric(qw(...))

valueMapperRaw = function(n, d)d[[n]]
valueMapperStandard = function(n, d) {
	if (is.na(d[[n]])) '{\\bf Value missing}' else (if (is.null(d[[n]])) n else d[[n]])
}

# <N> maxIterations needs to be large as a new iteration is entered after each successful substitution
#	this is necessary, as 
mergeDictToString = function(d, s,
	valueMapper = valueMapperStandard,
	#valueMapper = function(s)ifelse(is.na(d[[n]]), '{\\bf Value missing}', d[[n]]),
	iterative = FALSE, re = FALSE, maxIterations = 1e4, doApplyValueMap = TRUE, doOrderKeys = TRUE, maxLength = 1e7) {
	ns = names(d);
	# proceed in order of decreasing key lengthes
	if (doOrderKeys) ns = ns[rev(order(sapply(ns, nchar)))];
	for (i in 1:maxIterations) {
		s0 = s;
		for (n in ns) {
			# counteract undocumented string interpolation
			subst = if (doApplyValueMap)
				gsub("[\\\\]", "\\\\\\\\", valueMapper(n, d), perl = TRUE)
				else d[[n]];
			# <!> quoting
			if (!re) n = sprintf("\\Q%s\\E", n);
			s = gsub(n, firstDef(subst, ""), s, perl = TRUE, fixed = FALSE);
			# <A> if any substitution was made, it is nescessary to reiterate ns to preserver order
			#	of substitutions
			if (iterative && s != s0) break;
		}
		if (!iterative || s == s0 || nchar(s) > maxLength) break;
	}
	s
}
mergeDictToStringV = Vectorize(mergeDictToString, 's', SIMPLIFY = TRUE, USE.NAMES = TRUE);

mergeDictToVector = function(d, v) { unlist(ifelse(is.na(names(d[v])), v, d[v])) }

mergeDictToDict = function(dMap, dValues, ..., recursive = TRUE) {
	r = lapply(dValues, function(v) {
		r = if (class(v) == 'list') {
			if (recursive) mergeDictToDict(dMap, v, ...) else v
		} else if (class(v) == 'character') mergeDictToString(dMap, v, ...) else v;
		r
	});
	r
}

# double quote if needed
qsSingle = function(s, force = FALSE) {
	# <N> better implementation possible: detect unquoted white-space
	if (force || length(fetchRegexpr('[ \t"()\\[\\]:,]', s)) > 0) {
		s = gsub('([\\"])', '\\\\\\1', s);
		s = sprintf('"%s"', s);
	} else {
		s0 = gsub("([\\'])", '\\\\\\1', s);
		if (s0 != s) s = sprintf("$'%s'", s0);
	}
	s
}
qs = function(s, ...)sapply(s, qsSingle, ...)
# single quote if needed
qssSingle = function(s, force = FALSE) {
	# <N> better implementation possible: detect unquoted white-space
	if (force || nchar(s) == 0 || length(fetchRegexpr("[ \t'\"()\\[\\]:,]", s)) > 0) {
		s = gsub("(['])", "'\"'\"'", s);
		s = sprintf("'%s'", s);
	}
	s
}
qss = function(s, ...)sapply(s, qssSingle, ...)
# include special case for home folder expansion: do not quote initial '~'
qsSinglePath = function(s, ...) {
	if (s == '~')
		s else
	if (nchar(s) >= 2 && substring(s, 1, 2) == '~/')
		con('~/', qsSingle(substring(s, 3), ...)) else
		qsSingle(s, ...)
}
# include special case for home folder expansion"
qsPath = function(s, ...)sapply(s, qsSinglePath, ...)

#' Return sub-strings indicated by positions or produce a string by substituting those strings with
#'	replacements
#'
#' The function behaves similar to sprintf, except that character sequences to be substituted are
#' indicated by name.
#'
#' @param s template string
#' @param start vector of start positions of substrings to substitute
#' @param length vector of lengthes of substrings to substitute
#' @param replacement vector of strings to subsitute. If missing, \code{Substr} returns sub-strings indicated
#'	by start/length
#' @return character vector containing extracted sub-strings
#'
# #' @examples
# #' \dontrun{
# #' print(Substr("abc", c(2, 3), c(1, 1), c("def", 'jkl')));
# #' print(Substr("abcdef", c(2, 3, 5), c(1, 1, 1), c("123", '456', '789')));
# #' print(Substr("abcdef", c(1, 3, 5), c(1, 1, 1), c("123", '456', '789')));
# #' print(Substr("abcdef", c(1, 3, 5), c(0, 1, 0), c("123", '456', '789')));
# #' }
Substr = function(s, start, length, replacement) {
	if (missing(replacement)) return(substr(s, start, start + length - 1));
	start = c(start, nchar(s) + 1);
	l = sapply(seq_along(replacement), function(i)c(
		replacement[i],
		substr(s, start[i] + length[i], start[i + 1] - 1)
	));
	l = c(substr(s, 1, start[1] - 1), as.vector(l));
	r = join(as.vector(l), sep = '');
	r
}

sprintfIgnoreEscapes = function(r) {
	m = r$match;
	L = attr(r$positions, 'capture.length');
	if (!(any(L[, 1] == 0 & L[, 2] == 0))) return(r);
	Is = which(L[, 1] == 0 & L[, 2] == 0);
	r0 = r;
	r$match = r0$match[-Is];
	r$positions = r0$positions[-Is];
	attr(r$positions, 'match.length') = attr(r0$positions, 'match.length')[-Is];
	attr(r$positions, 'capture.start') = attr(r0$positions, 'capture.start')[-Is, , drop = FALSE];
	attr(r$positions, 'capture.length') = attr(r0$positions, 'capture.length')[-Is, , drop = FALSE];
	attr(r$positions, 'capture.names') = attr(r0$positions, 'capture.names')[-Is];
	return(r);
}

# <!> quoting
#'	Produce string by substituting placeholders
#'
#' The function behaves similar to sprintf, except that character sequences to be substituted are
#' indicated by name. To be implemented: *-specifications
#'
#' #@param s template string
#' #@param d values to substitute into \code{s}
#' #@param template template for substitution pattern. Within this pattern \code{__DICT_KEY__} is
#' # substituted for a key in \code{d}. This string \code{k} is substituted in \code{s} with \code{d[[k]]}.
#' @param .fmt formatting string into which values are interpolated (see details)
#' @param values list or vector of values to be used for interpolation
#' @param sprintf_cartesian boolean to indicate whether cartesian product of values should be used.
#'   Otherwise standard recyling rules apply.
#' @param envir environment in which values are to be evaluated
#' @return Interpolated character string
#'
# #' @examples
# #' \dontrun{
# #' Sprintf('These are N %{N} characters.', list(N = 10));
# #' Sprintf('These are N %{N}d characters.', list(N = 10));
# #' Sprintf('These are N %{N}02d characters.', list(N = 10));
# #' }
Sprintfl = function(.fmt, values, sprintf_cartesian = FALSE, envir = parent.frame()) {
	dict = extraValues = list();
	for (i in seq_along(values)) {
		if (is.list(values[[i]]))
			dict = merge.lists(dict, values[[i]]) else
		if (!is.null(names(values)[i]) && names(values)[i] != '')
			dict = merge.lists(dict, values[i]) else
			extraValues = c(extraValues, values[i]);
	}
# 	re = '(?x)(?:
# 		(?:^|[^%]|(?:%%)+)\\K
# 		[%]
# 			(?:[{]([^{}\\*\'"]*)[}])?
# 		((?:[-]?[*\\d]*[.]?[*\\d]*)?(?:[sdfegG]|))(?=[^%sdfegG]|$)
# 	)';
	# <!> new, untested regexpr as of 22.5.2014
	# un-interpolated formats do no longer work
# 	re = '(?xs)(?:
# 		(?:[^%]+|(?:%%)+)*\\K
# 		[%]
# 			(?:[{]([^{}\\*\'"]*)[}])?
# 		((?:[-]?[*\\d]*[.]?[*\\d]*)?(?:[sdfegGDQqu]|))(?=[^sdfegGDQqu]|$)
# 	)';

	re = '(?xs)(?:
		(?:[^%]+|(?:%%)+)*
		\\K[%]
			(?:[{]([^{}\\*\'"]*)[}])?
		((?:[-]?[*\\d]*[.]?[*\\d]*)?(?:[stdfegGDQqu]|))(?=[^stdfegGDQqu]|$)
	)';

# 	re = '(?xs)(?:
# 		(?:(?:[^%]+)(?:(?:%%)+(?:[^%]+))*)
# 		[%]
# 			(?:[{]([^{}\\*\'"]*)[}])?
# 		((?:[-]?[*\\d]*[.]?[*\\d]*)?(?:[sdfegGDQqu]|))(?=[^sdfegGDQqu]|$)
# 	)';

	r = fetchRegexpr(re, .fmt, capturesAll = TRUE, returnMatchPositions = TRUE);
	r = sprintfIgnoreEscapes(r);
	# <p> nothing to format
	if (length(r$match) == 0) return(.fmt);
	typesRaw = sapply(r$match, function(m)ifelse(m[2] == '', 's', m[2]));
	types = ifelse(typesRaw %in% c('D', 'Q'), 's', typesRaw);
	fmts = sapply(r$match, function(m)sprintf('%%%s',
		ifelse(m[2] %in% c('', 'D', 'Q', 'q', 't', 'u'), 's', m[2])));
	fmt1 = Substr(.fmt, r$positions, attr(r$positions, 'match.length'), fmts);

	keys = sapply(r$match, function(i)i[1]);
	nonKeysI = cumsum(keys == '');	# indeces of values not passed by name
	nonKeysIdcs = which(keys == '');

	# <p> collect all values
	allValues = c(extraValues, dict);
	# get interpolation variables
	interpolation = nlapply(keys[keys != ''], function(k)
		if (!is.null(allValues[[k]])) NULL else rget(k, default = NA, envir = envir)
	);
	# <p> handle %D: current day
	keys[typesRaw == 'D'] = '..Sprintf.date..';
	dateValue = if (sum(typesRaw == 'D'))
		list(`..Sprintf.date..` = format(Sys.time(), '%Y%m%d')) else
		list();
	allValues = c(allValues, dateValue, List_(interpolation, rm.null = TRUE));

	# 14.9.2015 -> convert to indeces
	# build value combinations
	listedValues = lapply(keys, function(k)allValues[[k]]);
	dictDf = if (!sprintf_cartesian) Df_(listedValues) else merge.multi.list(listedValues);
	# fill names of anonymous formats
	keys[keys == ''] = names(dictDf)[Seq(1, sum(nonKeysI != 0))];
	# due to repeat rules of R vectors might have been converted to factors
	#dictDf = Df_(dictDf, as_character = unique(keys[types == 's']));
	dictDf = Df_(dictDf, as_character = which(types == 's'));
	
	# <p> conversion <i>: new function
	#colsQ = keys[typesRaw == 'Q'];
	# <!> switch to index based transformation on account of duplicate keys
	colsQ = which(typesRaw == 'Q');
	dictDf[, colsQ] = apply(dictDf[, colsQ, drop = FALSE], 2, qsPath, force = TRUE);
	#colsq = keys[typesRaw == 'q'];
	colsq = which(typesRaw == 'q');;
	dictDf[, colsq] = apply(dictDf[, colsq, drop = FALSE], 2, qss);
	colst = which(typesRaw == 't');;
	dictDf[, colst] = apply(dictDf[, colst, drop = FALSE], 2, qss, force = TRUE);

	colsu = which(typesRaw == 'u');;
	dictDf[, colsu] = apply(dictDf[, colsu, drop = FALSE], 2, uc.first);

	colsd = which(typesRaw == 'd');;
	dictDf[, colsd] = apply(dictDf[, colsd, drop = FALSE], 2, as.integer);
	s = sapply(1:nrow(dictDf), function(i) {
		valueDict = as.list(dictDf[i, , drop = FALSE]);
# 		sprintfValues = lapply(seq_along(keys), function(i)
# 			ifelse(keys[i] == '', extraValues[[nonKeysI[i]]],
# 				firstDef(valueDict[[keys[i]]], rget(keys[i], default = '__no value__'), pos = -2)));
# 		sprintfValues = lapply(seq_along(keys), function(i)
# 			firstDef(valueDict[[keys[i]]], rget(keys[i], default = '__no value__', envir = envir)));
		#sprintfValues = lapply(seq_along(keys), function(i)valueDict[[keys[i]]]);
		#do.call(sprintf, c(list(fmt = fmt1), sprintfValues))
		# <!> simplify above two lines, now robust against duplicated entries -> <i> needs unit tests
		names(valueDict) = NULL;
		do.call(sprintf, c(list(fmt = fmt1), valueDict))
	});
	s
}

# 18.10.2019: fmt -> .fmt to avoid confusion with abbreviated named arguments (e.g. f = x substitutes fmt)
Sprintf = sprintd = function(.fmt, ..., sprintf_cartesian = FALSE, envir = parent.frame(),
	resetNames = TRUE, drop = TRUE) {
	r = sapply(.fmt, function(.fmt)
		Sprintfl(.fmt, list(...), sprintf_cartesian = sprintf_cartesian, envir = envir),
		USE.NAMES = !resetNames);
	# <!> special case when a single .fmt is provided -> do not return matrix for several values
	if (drop && length(.fmt) == 1) r = avu(r);
	r
}

#r = getPatternFromStrings(DOC, '(?:\\nDOCUMENTATION_BEGIN:)([^\\n]+)\\n(.*?)(?:\\nDOCUMENTATION_END\\n)');
getPatternFromStrings = function(strings, pattern, keyIndex = 1) {
	r = lapply(strings, function(s) {
		ps = fetchRegexpr(pattern, s, capturesAll = TRUE);
		listKeyValue(sapply(ps, function(e)e[[keyIndex]]), sapply(ps, function(e)e[-keyIndex]));
	});
	r
}

getPatternFromFiles = function(files, locations = NULL, ...) {
	strings = sapply(files, function(f)readFile(f, prefixes = locations));
	getPatternFromStrings(strings, ...);
}

#
#	hex strings
#

asc = function(x)strtoi(charToRaw(x), 16L);
character.as.characters = function(str) {
	sapply(str, function(s) sapply(1:nchar(s), function(i)substr(str, i, i)));
}

# bit_most_sig in bits
hex2int = function(str, bit_most_sig = 32) {
	cs = rev(sapply(character.as.characters(tolower(str)), asc));
	cms = bit_most_sig / 4;	# character containing most significant bit
	is = ifelse(cs >= asc('a'), cs - asc('a') + 10, cs - asc('0'));
	flipSign = (length(is) >= cms && is[cms] >= 8);
	if (flipSign) is[cms] = is[cms] - 8;
	r = sum(sapply(1:length(is), function(i)(is[i] * 16^(i-1))));
	if (flipSign) r = r - 2^(bit_most_sig - 1);
	r = if (r == - 2^(bit_most_sig - 1)) NA else as.integer(r);
	r
}

# chunk_size in bits
hex2ints = function(str, chunk_size = 32) {
	l = nchar(str);
	csc = chunk_size / 4;	# chunk_size in characters
	chunks = (l + csc - 1) %/% csc;
	r = sapply(1:chunks, function(i)hex2int(substr(str, (i - 1)*csc + 1, min(l, i*csc))));
	r
}

#
#	<§> binary numbers/n-adic numbers
#

ord2base = dec2base = function(o, digits = 5, base = 2) {
	sapply(1:digits, function(i){(o %/% base^(i-1)) %% base})
}
base2ord = base2dec = function(v, base = 2) {
	sum(sapply(1:length(v), function(i)v[i] * base^(i-1)))
}

ord2bin = dec.to.bin = function(number, digits = 5) ord2base(number, digits, base = 2);
bin2ord = bin.to.dec = function(bin) base2ord(bin, base = 2);

# mixed base calculations
#sapply(1:length(base), function(i)((n %/% div[i]) %% base[i]));
cumprod1 = function(v)c(1, cumprod(pop(v)))
# ord2adic = function(n, base = rep(2, 5)) {
# 	div = cumprod1(base);
# 	(n %/% div) %% base
# }
# adic2ord = function(v, base = rep(2, 5)) {
# 	mult = cumprod1(base);
# 	(v %*% mult)[1, 1]
# }
ord2adic = function(n, base = rep(2, 5))((n %/% cumprod1(base)) %% base)
adic2ord = function(v, base = rep(2, 5))((v %*% cumprod1(base))[1, 1])

#
#	<Par> sequences
#

#'	Produce constrained sequences
#'
#' This is a wrapper around seq that adds constraints. Setting ascending, descending to NA reverts to
#' standard \code{seq} behaviour.
#'
#' @param ascending restrict sequences to be ascending; return empty list if to < from
#' @param descending restrict sequences to be descending; return empty list if from < to
#' @param from starting value
#' @param to ending value
#' @param neg boolean to indicate wheter sequence should be negated before return
#' @param ... parameters passed on to \code{seq}
#' @return sequence from \code{from} to \code{to}
# #' @examples
# #' \dontrun{
# #' Seq(1, 10, ascending = TRUE)
# #' Seq(1, 10, descending = TRUE)
# #' Seq(10, 1, ascending = NA)
# #' }
Seq = function(from, to, ..., ascending = TRUE, descending = !ascending, neg = FALSE) {
	# <!> order matters: if called with only descending == TRUE
	if (nif(descending) && to > from) return(if (neg) TRUE else c()) else
	if (nif(ascending) && from > to) return(if (neg) TRUE else c());
	s = seq(from, to, ...);
	r = if (neg) -s else s;
	r
}
SeqRows = function(o)Seq(1, nrow(o))

#' Produce index pairs for vector of counts
#'
#' @param counts vector of integers specifying counts
#' @return vector of pairs of indeces indicating the first and last element in a vector for the blocks 
#'  specified by \code{counts}
#' @keywords internal
# #' @examples
# #' \dontrun{
# #' count2blocks(c(1, 5, 3))
# #' }
count2blocks = function(counts) {
	ccts = cumsum(counts);
	fidcs = c(1, ccts[-length(ccts)] + 1);
	blks = as.vector(rbind(fidcs, fidcs + counts - 1));
	blks
}

#
#	expand a block list - for example as from count2blocks - to a list of integers
#
expandBlocks = function(blks) {
	applyL(matrix(blks, ncol = 2, byrow = TRUE), 1, function(r) { r[1]:r[2] } )
}

# split 1:M into N partitions, return row-wise range
splitListIndcs = function(M, N = 1, .compact = FALSE, .truncate = TRUE) {
	if (.truncate & M < N) N = M;
	if (.compact) {
		n = rep(ceiling(M / N), N);	# size of parts
		idcs = c(0, cumsum(n));
		idcs = idcs[idcs < M];
		idcs = c(idcs, M);
	} else {
		n = rep(floor(M / N), N);		# size of parts
		R = M - n[1] * N;
		n = n + c(rep(1, R), rep(0, N - R));
		idcs = c(0, cumsum(n));
	}
	idcs = cbind(idcs + 1, c(idcs[-1], 0))[-length(idcs), ];	# from:to in a row
	# <!> usual R degeneracy
	if (!is.matrix(idcs)) idcs = matrix(idcs, nrow = 1);
	idcs
}
splitListEls = function(l, N, returnElements = FALSE) {
	idcs = splitListIndcs(length(l), N);
	li = apply(idcs, 1, function(r)(if (returnElements) l[r[1]:r[2]] else r[1]:r[2]));
	# <!> R ambiguity of apply return type
	if (is.matrix(li)) li = lapply(1:(dim(li)[2]), function(i)li[, i]);
	if (is.vector(li)) li = as.list(li);;
	li
}

# @param l list of index positions from another object
# @return return vector indicating to which list element an index was assigned
# Example: glmnet accepts fold numbers per index (as opposed to a partitioning of elements)
index2listPosition = function(l) {
	N = sum(sapply(l, length));
	na = rep(NA, N);
	m = sapply(1:length(l), function(i)vector.assign(na, l[[i]], i, na.rm = NA));
	r = apply(m, 1, na.omit);
	r
}

# idcs start positions in ragged list, converted to ranges
idcsStart2range = function(idcs, N = max(idcs)) {
	if (length(idcs) == 0) return(NULL);
	vector.intercalate(idcs, c(shift(idcs - 1), N))
}

# splitting based on fractions
# voting percentages to seats
#	simple algorithm based on size of residuals
# tiePreferHigh: for tied residuals add/subtract seats to high indeces (TRUE) or low ones (FALSE)
splitSeatsForFractions = function(Nseats, fractions = vn(rep(1, Nfractions)), Nfractions,
	tiePreferHigh = TRUE) {
	# number of parties
	Nparties = length(fractions);
	# fractional seats
	Nseats0 = fractions * Nseats;
	# garuantee one seat, otherwise round to nearest
	Nseats1 = ifelse (Nseats0 < 1, 1, round(Nseats0));
	# individual mismatch
	Nresid = Nseats0 - Nseats1;
	# mismatch total
	diff = sum(Nseats1) - Nseats;
	# redistribute deficit/overshoot
	if (diff != 0) {
		Nresid1 = ifelse(Nresid < 0, 1, Nresid);	# too few vs too many, too few -> maximal value of 1
		# take seats from whom? We need abs(diff) seats.
		#subtr = order(Nresid1, decreasing = diff < 0)[1:abs(diff)];
		prio = if (tiePreferHigh) 1:Nparties else rev(1:Nparties);
		subtr = Order(Df(Nresid1, prio))[1:abs(diff)];
		# assume one round of correction is always sufficient <!>
		Nseats1[subtr] = Nseats1[subtr] - sign(diff);
	}
	Nseats1
}

# tranform number of elements (as from splitSeatsForFractions) into from:to per row in a matrix
counts2idcs = function(counts) {
	idcs = c(0, cumsum(counts));
	idcs = cbind(idcs + 1, c(idcs[-1], 0))[-length(idcs), ];
	if (is.null(counts)) return(idcs);	# matrix w/ 0 rows
	t2r(idcs)	# fails on counts == NULL
}

# N is partitioned into fractions from p, where each element of p partitions the remaining part of N
# procedure makes sure to leave space for length(p) elements
cumpartition = function(N, p) {
	I = c();	# indeces within 1:N
	for (i in 1:length(p)) {
		# partition remaining space (ifelse), leave room for subsequent indeces
		Ii = floor(p[i] * (ifelse(i == 1, N, N - I[i - 1]) - (length(p) - i))) + 1;
		I = c(I, ifelse(i == 1, Ii, I[i - 1] + Ii));
	}
	as.integer(I)
}

#' Extract parts of a nested structure based on the range from..to
#'
#'
#' @param Ns Vector of integers that specify the size of the substructures
#' @param from absolute index where to start extraction
#' @param to absolute index where to stop extraction
#' @return Return list of lists, where each basic list contains key \code{segment}
#'  (which of the elements of Ns) and key \code{range}, a list with elements \code{from} and \code{to},
#'  specifying which elements to use from
#'  that segment.
# #' @examples
# #' \dontrun{
# #'    # TestMe: TRUE1
# #'    subListFromRaggedIdcs(c(2, 4, 10, 15), 1, 20)
# #' }
subListFromRaggedIdcs = function(Ns, from = 1, to) {
	NsCS = cumsum(Ns);
	NsCSs = c(0, pop(NsCS));	# shifted cumsum
	segments = which(from <= NsCS & to > NsCSs);
	if (missing(to)) to = sum(segments);
	r = lapply(segments, function(segment){
		N = Ns[segment];	# list-call
		from_ = 1;
		to_ = N;
		if (segment == segments[1]) from_ = from - NsCSs[segment];
		if (segment == rev(segments)[1]) to_ = to - NsCSs[segment];
		r = list(segment = segment, range = list(from = from_, to = to_));
		r
	});
	r
}

#' Extract parts of nested lists based on the range from..to
#'
#'
#' @param from absolute index where to start extraction
#' @param to absolute index where to stop extraction
#' @param ls nested list structure (currently only two levels supported)
#' @return Return list of list, where each basic list contains key \code{segment}
#'  (which of the elements of Ns) and key \code{range}, a list with elements \code{from} and \code{to},
#'  specifying which elements to use from
#'  that segment.
subListFromRaggedLists = function(ls, from = 1, to = sum(sapply(ls, length))) {
	sl = subListFromRaggedIdcs(sapply(ls, length), from = from, to = to);
	r = lapply(sl, function(s) with(s, {
		r = ls[[segment]][range$from: range$to];
		r
	}));
	r = unlist.n(r, 1);
	r
}


#
#	<§> vector functions
#

# does the position exists in vector v
exists.pos = function(v, i)(is.vector(v) && !is.na(v[i]))

# for a vector blocked by blockSize N, return indeces of elements of block i
rangeBlock = function(i, N)(((i - 1)*N + 1):(i * N))

#
#	<par> lists
#

merge.lists = function(..., ignore.nulls = TRUE, listOfLists = FALSE, concat = FALSE, useIndeces = FALSE) {
	lists = if (listOfLists) c(...) else list(...);
	l1 = lists[[1]];
	if (length(lists) > 1) for (i in 2:length(lists)) {
		l2 = lists[[i]];
		ns = if (useIndeces) 1L:length(l2) else names(l2);
		for(n in ns) {
			if (is.null(n)) print("Warning: tried to merge NULL key");
			if (!is.null(n) & (!ignore.nulls | !is.null(l2[[n]]))) {
				if (concat) l1[[n]] = c(l1[[n]], l2[[n]]) else l1[[n]] = l2[[n]];
			}
		}
	}
	l1
}

merge.lists.recursive = function(..., ignore.nulls = TRUE, listOfLists = FALSE, concat = FALSE) {
	lists = if (listOfLists) c(...) else list(...);
	l1 = lists[[1]];
	if (length(lists) > 1) for (i in 2:length(lists)) {
		l2 = lists[[i]];
		for(n in names(l2)) {
			if (is.null(n)) print("Warning: tried to merge NULL key");
			if (!is.null(n) & (!ignore.nulls | !is.null(l2[[n]])))
				l1[[n]] = if (is.list(l1[[n]]))
					merge.lists.recursive(l1[[n]], l2[[n]]) else
					(if (concat) c(l1[[n]], l2[[n]]) else l2[[n]])
		}
	}
	l1
}

unshift = function(l, listOfList = TRUE) {
	if (!listOfList) l = list(l);
	e1 = lapply(l, function(l0)if (is.list(l0)) l0[[1]] else l0[1]);
	r1 = lapply(l, function(l0)l0[-1]);
	r = list(elements = e1, remainder = r1);
	r
}

Merge.lists.raw = function(lists, ignore.nulls = TRUE, recursive = FALSE, keys = NULL) {
	if (!is.null(keys)) keys = unshift(keys);

	l1 = lists[[1]];
	if (length(lists) > 1) for (i in 2:length(lists)) {
		l2 = lists[[i]];
		for(n in names(l2)) {
			if (is.null(n)) print("Warning: tried to merge NULL key");
			if (!is.null(n) & (!ignore.nulls | !is.null(l2[[n]])))
				l1[[n]] = if (recursive && is.list(l1[[n]]) && (is.null(keys) || n %in% keys$elements))
					Merge.lists.raw(list(l1[[n]], l2[[n]]), ignore.nulls, recursive,
						if (is.null(keys)) NULL else keys$remainder) else
					l2[[n]]
		}
	}
	l1
}

Merge.lists = function(..., ignore.nulls = TRUE, listOfLists = FALSE, recursive = FALSE, keyPathes = NULL) {
	lists = if (listOfLists) c(...) else list(...);
	keys = if (!is.null(keyPathes)) splitString("[$]", keyPathes, simplify = FALSE) else NULL; 
	l = Merge.lists.raw(lists, ignore.nulls = ignore.nulls, recursive = recursive, keys = keys);
	l
}

# l: list of lists
# take parallel elements from l (1, ...) after recycling
list.combine = function(l, byRow = TRUE, names = NULL, doMerge = FALSE) {
	lR = Recycle(l, byRow = byRow);
	# <p> number of final elements
	N =	length(lR[[1]]);
	lC = lapply(1:N, function(i) {
		lol = list.kp(lR, Sprintf('[[%{i}d]]'));
		if (notE(names)) names(lol) = names;
		return(if (doMerge) merge.lists(lol, listOfLists = TRUE) else lol);
	});
	return(lC);
}
# inverse of unlist.n(, 1)
list.embed = function(l, key = 'key')lapply(l, function(e)SetNames(list(e), key));

# use.names preserves names and concatenates with lower level names
# reset sets names to top level names
unlist.n = function(l, n = 1, use.names = TRUE, reset = FALSE) {
	if (n > 0) for (i in 1:n) {
		ns = names(l);
		#names(l) = rep(NULL, length(l));	# <!> untested removal Tue Oct 19 17:11:53 2010
		l = unlist(l, recursive = FALSE, use.names = use.names);
		if (reset) names(l) = ns;
	}
	l
}

# <N> obsolete, better: with(l, { ...})
instantiate.list = function(l, n = 1) {
	for (nm in names(l)) {
 		eval.parent(parse(file = "", text = sprintf("%s = %s", nm, deparse(l[[nm]]))), n = n);
# 		if (is.integer(l[[nm]])) {
# 			eval.parent(parse(file = "", text = sprintf("%s = %d", nm, l[[nm]])), n = n);
# 		} else if (is.numeric(l[[nm]])) {
# 			eval.parent(parse(file = "", text = sprintf("%s = %f", nm, l[[nm]])), n = n);
# 		} else {
# 			eval.parent(parse(file = "", text = sprintf("%s = \"%s\"", nm, l[[nm]])), n = n);
# 		};
	}
}
# for use in testing code
instantiate = function(l, ..., envir = parent.frame()) {
	l0 = c(l, list(...));
	for (i in seq_along(l0)) assign(names(l0)[i], l0[[i]], envir = envir);
	invisible(l0)
}

# assume a list of lists (aka vector of dicts) and extract a certain key from each of the lists
list.key = function(v, key, unlist = TRUE, template = NULL, null2na = FALSE) {
	l = lapply(v, function(i){
		if (is.list(i)) {
			if (is.null(i[[key]])) { if (null2na) NA else NULL } else i[[key]]
		} else template});
	if (unlist) l = unlist(l);
	l
}

# iterative gsub
# substs: list with pairs re, sub
gsubI = function(substs, s, ...) {
	for (r in substs) s = gsub(r[[1]], r[[2]], s, ...);
	s
}

# concatenate lists, leave out lists only containing a single NULL
cList = function(...) {
	r = list(...);
	isNull = sapply(r, is.null);
	do.call(c, r[!isNull])
}

keyPathParse = function(kp) {
	kp = gsubI(
		list(c('\\*', '.'), c('\\[\\[(\\d+)\\]\\]', 'INDEX__\\1'), c('(\\$|^)\\s*[(]', '\\1PAR__('))
	, kp);
	Parse(kp)
}
keyPathExpression2key = function(e) {
	s = as.character(e);
	i = FetchRegexpr('INDEX__(\\d+)', s, captures = TRUE);
	r = if (s == '.') '*' else
		if (s == 'PAR__') NULL else
		if (!is.null(i)) return(as.integer(i)) else s;
	if (is.null(r)) NULL else list(r)
}

# level: level of nesting for parallel key pathes
# output: list or vector, a parallel pattern induces two more levels of list nesting
#	the parallel keyPathes each of which is forced to be a list
keyPathAstRaw = function(e, level = 0) {
	r = if (is.call(e)) {
		isPar = e[[1]] == '|';
		isPlain = e[[1]] == '$';
		levelN = ifelse(isPar, level + 1, 0);
		# simple walk through
		r = if (isPlain | isPar)
			cList(keyPathAstRaw(e[[2]], level = levelN), keyPathAstRaw(e[[3]], level = levelN))
		# this is a fake call generated by the PAR__ construct, the path leading to PAR__ is seen as an
		#	anonymous function
		else cList(keyPathAstRaw(e[[1]]), list(keyPathAstRaw(e[[2]], level = 1)));
		if (level & isPlain) list(r) else r
	} else if (is.name(e)) keyPathExpression2key(e) else {
		stop('malformed keyPath');
	}
	r
}

keyPathAst = function(kp) {
	unlist.n(keyPathAstRaw(keyPathParse(kp)[[1]]), n = 0);
}

list.kp.unquote = function(key) {
	# un-quote: remove single backslashes
	key = sub('(?<![\\\\])[\\\\](?![\\\\])', '', key, perl = TRUE);
	# de-quote: double backslashes become single backslashes
	key = sub('\\\\', '\\', key, fixed = TRUE);
	as.character(key)
}

# extract key path from list, general, recursive version
#	key path recursive worker
list.kprw = function(l, keys, unlist.pats, template, null2na, carryNames, test, keyAccess) {
	if (!length(keys)) return(l);
	key = keys[1];
	# <p> extract key
	r = if (key != "*") {
		index = fetchRegexpr("\\A\\[\\[(\\d+)\\]\\]\\Z", key, captures = TRUE);
		if (length(index) > 0) key = as.integer(index[[1]]);
		if (keyAccess[1] == '@') {
			r = slot(l, key);
			list.kprw(r, keys[-1], unlist.pats[-1], template, null2na, carryNames, test, keyAccess[-1]);
		} else if (is.list(l)) {
			# <N> logical(0) seen as NULL by second condition
			r = if (is.null(l[[key]]) || length(l[[key]]) == 0) {
					if (null2na) { NA } else firstDef(template, NULL)
				} else l[[key]];
			if (length(keys) > 1)
				list.kprw(r, keys[-1], unlist.pats[-1], template, null2na, carryNames, test, keyAccess[-1]) else
				if (test) !(is.null(r) || all(is.na(r))) else r;
		} else if (class(l) %in% c('character')) {
			if (notE(names(l))) l[names(l) %in% key] else l[key]
		} else if (class(l) %in% c('data.frame', 'matrix')) {
			l[, key]
		} else if (class(l) %in% c('numeric', 'integer')) {
			l[key]
		} else return(template);
# 		{
# 			r = template;
# 			attr(r, 'names') = keys[last(keys)];
# 			print(c(keys, r));
# 			return(r);
# 		}
	} else {
		if (length(keys) > 1)
			lapply(l, function(sl)
				list.kprw(sl, keys[-1], unlist.pats[-1], template, null2na, carryNames, test, keyAccess[-1])
			) else l;
	}
	# <p> unlisting
	if (notE(unlist.pats)) if (unlist.pats[1]) r = unlist.n(r, 1, reset = carryNames);
	r
}

# extract key path from list, general, recursive version
#	key path recursive worker: parallel keys
#	iterate over recursive keys
list.kprwPar = function(l, keys, ...) {
	key = keys[1];
	r = if (length(fetchRegexpr("\\|", key)) > 0) {
		parKeys = sapply(splitString('\\|', key), list.kp.unquote);
		r = lapply(parKeys, function(key)list.kprwkp(l, c(key, keys[-1]), ...));
		unlist.n(r, 1);
	} else list.kprw(l, keys, ...);
	r
}

# worker: keypath
list.kprwkp = function(l, keyPath, ..., keyAccess) {
	keysNew = fetchRegexpr("(?:[a-zA-Z0-9_.|\\[\\]*]+(?:\\\\[$@])?)+", keyPath[1]);
	keys = c(keysNew, keyPath[-1]);
	r = list.kprwPar(l, keys, ...);
	r
}

list.kp.keys = function(keyPath) fetchRegexpr("[^$@]+", keyPath);
list.kp.method = function(keyPath) fetchRegexpr("[$@]", keyPath);

# wrapper for list.kprw
# keyPath obeys EL1 $ EL2 $ ..., where ELn is '*' or a literal
# unlist.pat is pattern of truth values TR1 $ TR2 $..., where TRn is in 'T|F' and specifies unlist actions
# carryNames determines names to be carried over from the top level in case of unlist
list.kpr = function(l, keyPath, do.unlist = FALSE, template = NULL,
	null2na = FALSE, unlist.pat = NULL, carryNames = TRUE, as.matrix = FALSE, test = FALSE) {
	keys = list.kp.keys(keyPath);
	# list or slot?
	keyAccess = list.kp.method(keyPath);
	# if first element is '*', assume list
	if (length(keyAccess) < length(keys)) keyAccess = c('$', keyAccess);
	unlist.pats = if (notE(unlist.pat)) as.logical(fetchRegexpr("[^$]+", unlist.pat)) else NULL;

	# parallel keys
	#r = list.kprwkp(l, keyPath, unlist.pats, template, null2na, carryNames, test = test);
	r = list.kprw(l, keys, unlist.pats, template, null2na, carryNames, test = test, keyAccess = keyAccess);
	if (do.unlist) { r = unlist(r); }
	if (as.matrix) r = t(sapply(r, function(e)e));
	r
}
# extract key path from list
# <!> interface change: unlist -> do.unlist (Wed Sep 29 18:16:05 2010)
# test: test existance instead of returning value
list.kp = function(l, keyPath, do.unlist = FALSE, template = NULL, null2na = FALSE, test = FALSE, n,
	pathAsIs = FALSE) {
	fullPath = if (pathAsIs) keyPath else sprintf("*$%s", keyPath);
	r = list.kpr(l, fullPath, do.unlist = do.unlist,
		template = template, null2na = null2na, test = test);
	if (!missing(n)) r = unlist.n(r, n);
	r
}

list.kpu = function(..., do.unlist = TRUE)list.kp(..., do.unlist = do.unlist);
# allow for slot access
list.Kpu = function(..., do.unlist = TRUE)list.kp(..., do.unlist = do.unlist, pathAsIs = TRUE);

list.keys = function(l, keys, default = NA) {
	l = as.list(l);
	r = lapply(unlist(keys), function(key) if (is.null(l[[key]])) default else l[[key]]);
	r
}

# make A > B into B > A
listReverseHierarchy = function(l, unlist = FALSE) {
	ns = names(l[[1]]);
	r = lapply(ns, function(n)list.kp(l, n, do.unlist = unlist));
	names(r) = ns;
	return(r);
}

# return pair of lists: keys to access element, value
list.flatten_raw = function(lol, prefix = list(), ignoreIndeces = TRUE) {
	#ns = names(lol);
	#if (is.null(ns)) ns = 1:length(lol);
	# <A> duplicate names, handled by nelapply
	r = unlist.n(Nelapply(lol, function(n, e) {
		# generate list of lists, each entry being list(key path), value
		if (is.list(e)) {
			prefix = if (ignoreIndeces && is.integer(n)) prefix else c(prefix, list(n));
			list.flatten_raw(e, prefix = prefix)
		} else list(list(c(prefix, list(n)), e));
	}), 1);
	return(r);
}

# return flat list with keys representing keyPath
#	alternative: return pair with keys, values, keys joined or list
list.flatten = function(lol, joinby = NULL, aspairs = TRUE, ignoreIndeces = TRUE) {
	r = list.flatten_raw(lol, ignoreIndeces = ignoreIndeces);
	if (aspairs && is.null(joinby)) return(r);
	keys = list.kp(r, '[[1]]');
	values = list.kp(r, '[[2]]');
	if (is.null(joinby)) return(list(keys = keys, values = values));
	return(listKeyValue(sapply(keys, join, sep = joinby), values));
}


# analogous to list.merge; instead of merging, append entries
list.accrue = function(lol, sep = ':', ignoreIndeces = TRUE) {
	es = list.flatten(lol, joinby = sep, ignoreIndeces = ignoreIndeces);
	ns = names(es)
	r = list();
	for (i in seq_along(es)) {
		r[[ns[i]]] = union(r[[ns[i]]], es[[i]])
	}
	return(r);
}

null2na = function(l) {
	if (!length(l)) return(l);
	l[sapply(l, is.null)] = NA;
	return(l);
}


# return list without listed keys
list.min  = function(l, keys) {
	l[-which.indeces(keys, names(l))]
}
# list generation on steroids (wraps other functions)
.list = function(l, .min = NULL) {
	if (!is.null(.min)) l = list.min(l, .min);
	l
}
# get apply
gapply = function(l, key, unlist = FALSE)list.key(l, key, unlist)
# construct list as a dictionary for given keys and values
listKV = listKeyValue = function(keys, values, doRecycle = TRUE) {
	if (length(keys) != length(values) && doRecycle) {
		r = recycle(keys, values);
		keys = r[[1]];
		values = r[[2]];
	}
	if (length(keys) != length(values))
		stop("listKeyValue: number of provided keys does not match that of values");

	l = as.list(values);
	names(l) = keys;
	l
}
listNamed = function(l, names)setNames(l, names)
vectorNamed = function(v, names) {
	if (length(names) > length(v)) stop("vectorNamed: more names than vector elements");
	names(v) = names;
	v
}

vn = vectorNormed = function(v, type = 'O') {
	v0 = as.matrix(v);
	v0n = apply(v0, 2, function(v)norm(as.matrix(v), type = type));
	r = if (!is.matrix(v)) v/v0n else sapply(1:ncol(v), function(i) v[, i] / v0n[i]);
	r
}

#listInverse = function(l)listKeyValue(avu(l), names(l));
listInverse = function(l, toNA = FALSE) {
	n = sapply(l, length);
	# <p> values of inverse map
	vs = rep.each(names(l), n);
	# <p> construct list
	r = listKeyValue(avu(l, recursive = FALSE, toNA = toNA), vs);
	r
}

# name the list elements by the iterated vector elements ns (names)
nlapply = function(ns, f, ...) {
	if (is.list(ns)) ns = names(ns);
	r = lapply(ns, f, ...);
	names(r) = ns;
	r
}
nelapply = function(l, f, ..., name = '') {
	ns = names(l);
	if (is.null(ns)) ns = ( if (notE(name)) rep(name, length(l)) else seq_along(l) );
	r = lapply(seq_along(l), function(i, ...)f(ns[i], l[[i]], ...), ...);
	names(r) = ns;
	r
}
Nelapply = function(l, f, ..., name = NULL)nelapply(l, f, ..., name = name)

ilapply = function(l, f, ...) {
	r = lapply(1:length(l), function(i)f(l[[i]], i, ...));
	if (!is.null(names(l))) names(r) = names(l);
	r
}
einXapply = function(v, f, ..., einXapplyIterator = lapply) {
	l = as.list(v);
	ns = names(l);
	r = einXapplyIterator(seq_along(l), function(i)f(l[[i]], i, ns[i], ...));
	if (length(r) > 0) names(r) = ns;
	r
}

# pass element, index, name
einlapply = function(l, f = Identity, ...)einXapply(l, f, ..., einXapplyIterator = lapply);

# pass element, index
eilapply = function(l, f, ...) {
	r = lapply(seq_along(l), function(i)f(l[[i]], i, ...));
	names(r) = names(l);
	r
}
eisapply = function(v, f, ...) {
	l = as.list(v);
	r = sapply(seq_along(l), function(i)f(l[[i]], i, ...));
	names(r) = names(v);
	r
}
ensapply = function(l0, f, ...) {
	l = as.list(l0);
	ns = names(l);
	r = sapply(seq_along(l), function(i, ...)f(l[[i]], ns[i], ...), ...);
	names(r) = ns;
	r
}
einsapply = function(v, f = Identity, ...)einXapply(v, f, ..., einXapplyIterator = sapply)

kvlapply = function(l, f, ...) {
	ns = names(l);
	r = lapply(1:length(l), function(i)f(ns[i], l[[i]], ...));
	names(r) = ns;
	r
}
pairsapply = pairsapplyVL = function(l1, l2, f, ..., simplify = TRUE, USE.NAMES = TRUE) {
	if (length(l1) != length(l2)) stop('pairsapply: pair of collections of unequal length.');
	r = sapply(seq_along(l1), function(i)f(l1[i], l2[[i]], ...),
		simplify = simplify, USE.NAMES = USE.NAMES);
	r
}
pairsapplyLV = function(l1, l2, f, ..., simplify = TRUE, USE.NAMES = TRUE) {
	if (length(l1) != length(l2)) stop('pairsapply: pair of collections of unequal length.');
	r = sapply(seq_along(l1), function(i)f(l1[[i]], l2[i], ...),
		simplify = simplify, USE.NAMES = USE.NAMES);
	r
}
pairslapply = function(l1, l2, f, ...) {
	if (length(l1) != length(l2)) stop('pairslapply: pair of collections of unequal length.');
	r = lapply(seq_along(l1), function(i)f(l1[[i]], l2[[i]], ...));
	names(r) = names(l1);
	r
}

sapplyWoI = function(v, f, ...)sapply(v, function(i, ...)f(...), ...)
lapplyWoI = function(v, f, ...)lapply(v, function(i, ...)f(...), ...)

dfapply = function(Df__, f__) {
	r = lapply(1:nrow(Df__), function(i) {
		r = Df__[i, ];
		return(Df_(f__(as.list(r))));
	});
	Dfr = do.call(rbind, r);
	return(Dfr);
}

list.grep = sublist = filterList = function(o, f, ...) {
	l = if (!is.function(f)) f else sapply(o, f, ...);
	if (length(l) == 0) l = NULL;	#list corner case
	r = o[l];
	return(r);
}



# <i> copy MARGIN handling from apply (aperm)
lapplyDir = function(m, MARGIN, f_, ..., drop = FALSE) {
	selector = if (MARGIN == 1)
		function(m, i)m[i, , drop = drop] else
		function(m, i)m[, i, drop = drop];
	setNames(lapply(1:dim(m)[MARGIN], function(i)f_(selector(m, i), ...)), Dimnames(m, MARGIN))
}

# <!> as matrix to avoid warning
#lapplyRows = function(m, ...)lapply(split(as.matrix(m), row(m)), ...)
# lapplyRows = function(m, f_, ..., drop = FALSE)
# 	setNames(lapply(1:nrow(m), function(i)f_(m[i, , drop = drop], ...), ...), Row.names(m))
lapplyRows = function(m, f_, ..., drop = FALSE)lapplyDir(m, 1, f_ = f_, ..., drop = drop)
lapplyCols = function(m, f_, ..., drop = FALSE)lapplyDir(m, 2, f_ = f_, ..., drop = drop)

	
getElement = function(v, i)if (is.list(v)) v[[i]] else v[i];
# unify w/ list.takenFrom -> tests
List.takenFrom = function(listOfLists, v)
	lapply(1:length(listOfLists), function(j)getElement(listOfLists[[j]], v[j]));
# tuple-apply
tuapply = function(..., fct = Identity, args = list(), names = NULL) {
	tupels = list(...);
	M = length(tupels);
	Ns = sapply(tupels, length);
	N = Ns[1];
	if (any(Ns != N)) stop('Indexable elements not of same length');
	r = lapply(1:N, function(i)do.call(fct, c(List.takenFrom(tupels, rep(i, M)), args)));
	if (is.null(names) && !is.null(base::names(tupels[[1]]))) names = base::names(tupels[[1]]);
	if (!is.null(names)) base::names(r) = names;
	r
}

undrop2row = function(e)(if (is.vector(e)) matrix(e, ncol = length(e)) else e);
Lundrop2row = function(l)lapply(l, undrop2row);

undrop2col = function(e)(if (is.vector(e)) matrix(e, nrow = length(e)) else e);
Lundrop2col = function(l)lapply(l, undrop2col);

# return list from apply (undo simplify)
applyL = function(X, MARGIN, FUN, ...) {
	r = apply(X, MARGIN, FUN, ...);
	if (is.matrix(r)) return(lapply(1:ncol(r), function(i)r[, i]));
	if (!is.list(r) && is.vector(r)) return(lapply(1:length(r), function(i)r[i]));
	return(r);
}
# USE.NAMES logic reversed for sapply
sapplyn = function(l, f, ...)sapply(l, f, ..., USE.NAMES = FALSE);
list.with.names = function(..., .key = 'name') {
	l = list(...);
	ns = names(l);
	r = nlapply(l, function(n) c(l[[n]], listKeyValue(.key, n)));
	r
}

#
#	<p> names
#

Row.names = function(o, vivify = TRUE) {
	rn = row.names(o);
	if (is.null(rn) && vivify) 1:nrow(o) else rn
}
Col.names = function(o, vivify = TRUE) {
	rn = if (is.matrix(o)) dimnames(o)[[2]] else names(o);
	if (is.null(rn) && vivify) 1:ncol(o) else rn
}
# <i> implement general MARGINs
Dimnames = function(o, MARGIN, vivify = TRUE) {
	if (MARGIN == 1) Row.names(o, vivify) else Col.names(o, vivify)
}

SetNames = function(o, names, rnames, cnames, Dimnames, embed = FALSE) {
	if (!missing(Dimnames)) dimnames(o) = Dimnames;
	if (!missing(rnames)) row.names(o) = rnames;
	if (!missing(cnames)) dimnames(o)[[2]] = cnames;
	if (!missing(names)) {
		if (any(class(o) == 'matrix')) {
			if (embed) dimnames(o)[[2]][seq_along(names)] = names else {
				ns = if (is.list(names)) vector.replace(dimnames(o)[[2]], names) else names;
				if (is.null(dimnames(o))) dimnames(o) = list(NULL, ns) else dimnames(o)[[2]] = ns;
			}
		} else {
			if (embed) names(o)[seq_along(names)] = names else {
				names(o) = if (is.list(names)) vector.replace(names(o), names) else names;
			}
		}
	}
	o
}


#
#	<p> attributes
#

Attr = function(o, plus_, min_ = NULL) {
	if (!missing(plus_)) for (n in names(plus_)) { attr(o, n) = plus_[[n]]; }
	if (Nif(min_)) for (a in min_) { attr(o, a) = NULL; }
	o
}

#
#	<par> data type conversions
#

# assure m has at least 1 column
to.col = function(m) { if (is.null(dim(m))) t(t(m)) else m }
col.frame = function(l, col.name = 'value', minus = NULL, ignore.null = TRUE,
	do.paste = NULL, do.format = TRUE, digits = 3, plus = NULL) {
	if (ignore.null) { for (n in names(l)) { if (is.null(l[[n]])) l[[n]] = NULL; } }
	if (!is.null(minus)) { for (n in minus) { l[[n]] = NULL; } }
	my.names = if (!is.null(plus)) plus else names(l);
	digits = if (length(digits) > 1) digits else rep(digits, length(l));
	if (!is.null(do.paste)) {
		if (do.format) {
			i = 1;
			for (n in my.names) { if (is.vector(l[[n]])) {
				l[[n]] = paste(sapply(l[[n]],
						function(e){if (is.numeric(e)) sprintf("%.*f", digits[i], e) else e}
					), collapse = do.paste)
				i = i + 1;
			}}
		} else {
			for (n in my.names) { if (is.vector(l[[n]])) l[[n]] = paste(l[[n]], collapse = do.paste) }
		}
	}
	f = as.data.frame(l);
	if (dim(f)[2] > length(col.name) && length(col.name) == 1)
		row.names(f) = paste(col.name, 1:dim(f)[1], sep = "")
	else row.names(f) = c(col.name);
	t(f)
}

# <i> collect recursively until list or data.frame
# convert list of lists to data frame (assuming identical keys for each sub list)
#	also works on list of vectors
listOfLists2data.frame = function(l, idColumn = "id", .names = NULL) {
	# collect keys
	keys = if (is.list(l[[1]]))
		sort(unique(as.vector(unlist(sapply(l, function(e)names(e)))))) else 1:length(l[[1]]);
	if (is.null(.names)) .names = keys;
	# row names
	rows = names(l);
	if (is.null(rows)) rows = 1:length(l);
	# build df

	#df = t(sapply(rows, function(r) { unlist(l[[r]][keys]) }));
	df = t(sapply(rows, function(r)list2df(l[[r]], keys)));
	df = if (!is.null(idColumn)) {
		data.frame.types(data.frame(..idColumn.. = rows, df),
			row.names = 1:length(rows), names = c(idColumn, .names));
	} else {
		data.frame.types(df, row.names = rows, names = .names);
	}
	df
}

# resetColNames: reset column names to names of first data frame
# colsFromFirstDf: take columns from the first data frame
# <i> improved algorithm: unlist everything, bind together: cave: data types,
#	strictly valid only for matrices
# Use cases:
#	list with named vectors: get data frame that contains all vectors with all possible names represented
#		listOfDataFrames2data.frame(cfs, colsFromUnion = TRUE, do.transpose = TRUE, idColumn = NULL);
listOfDataFrames2data.frame = function(l, idColumn = "id", do.unlist = TRUE, direction = rbind,
	resetColNames = TRUE, colsFromFirstDf = FALSE, colsFromUnion = FALSE, do.transpose = FALSE, idAsFactor = FALSE,
	row.names = FALSE) {
	# row names
	# <!> 2009-11-20 changed from: rows = firstDef(names(l), list(1:length(l)));
	rows = firstDef(names(l), 1:length(l));
	# columns
	ns = NULL;
	if (colsFromUnion) {
		ns = unique(unlist(lapply(l, names)));
		# get data.frame names
		ns = names(do.call(data.frame, listKeyValue(ns, rep(NA, length(ns)))));
		resetColNames = FALSE;	# <!> mutually exclusive
	}
	# build df
	df = NULL;
	for (i in 1:length(rows)) {
		if (is.null(l[[i]])) next;	# ignore empty entries
		# <p> force to data frame
		df0 = if (do.transpose) as.data.frame(t(l[[i]])) else as.data.frame(l[[i]]);
		# <p> homogenize columns
		if (colsFromUnion) {
			# add missing columns
			ns0 = setdiff(ns, names(df0));
			df0 = do.call(data.frame, c(list(df0), listKeyValue(ns0, rep(NA, length(ns0)))));
			# correct order of columns
			df0 = df0[, ns];
		}
		if (!is.null(df)) {
			if (colsFromFirstDf) df0 = df0[, names(df)] else
			if (resetColNames) {
				names(df0) = if (is.null(idColumn)) names(df) else names(df)[-1];
			}
		}
		# <p> add id column
		df0 = if (is.null(idColumn)) df0 else cbind(rep(rows[i], dim(df0)[1]), df0);
		# <A> case differentiation should not me necessary
		df = if (i == 1) df0 else direction(df, df0);
	}
	if (!is.null(idColumn)) names(df)[1] = idColumn;
	if (do.unlist) for (n in names(df)) { df[[n]] = unlist(df[[n]]); }
	if (idAsFactor) df[[idColumn]] = as.factor(df[[idColumn]]);
	if (!row.names) row.names(df) = NULL;
	df
}
cbindDataFrames = function(l, do.unlist = FALSE, colsFromUnion = FALSE) {
	listOfDataFrames2data.frame(l, idColumn = NULL, do.unlist = do.unlist, direction = cbind,
		resetColNames = FALSE, colsFromUnion = colsFromUnion)
}
# @param embed corresponds to colsFromUnion in listOfDataFrames2data.frame
RbindDfs = function(dfl, namesFromFirst = TRUE, embed = FALSE) {
	if (namesFromFirst && !embed) dfl = lapply(dfl, setNames, nm = names(dfl[[1]]));
	if (embed) {
		ns = unique(unlist(sapply(dfl, names)));
		df0 = Df_(listKeyValue(ns, rep(NA, length(ns))));
		dfl = lapply(dfl, function(d)cbind(d, df0[, setdiff(ns, names(d)), drop = FALSE]));
	}
	do.call(rbind, dfl)
}

rbindDataFrames = function(l, do.unlist = FALSE, useDisk = FALSE, idColumn = NULL, transpose = FALSE,
	resetColNames = FALSE, colsFromFirstDf = FALSE, idAsFactor = FALSE) {
	r = if (useDisk) {
		tempTable = tempfile();
		for (i in 1:length(l)) {
			d0 = l[[i]];
			if (class(d0) != 'data.frame') d0 = as.data.frame(d0);
			if (transpose) d0 = t(d0);
			if (!is.null(idColumn)) {
				d0 = data.frame(idColumn = names(l)[i], d0);
				names(d0)[1] = idColumn;
			}
			write.table(d0, file = tempTable, col.names = i == 1, append = i != 1, row.names = FALSE);
		}
		read.table(tempTable, header = TRUE, as.is = TRUE);
	} else {
		listOfDataFrames2data.frame(l, idColumn = idColumn, do.unlist = do.unlist,
			direction = rbind, resetColNames = resetColNames, colsFromFirstDf = colsFromFirstDf,
			idAsFactor = idAsFactor)
	}
	r
}

# names2col assigns names of the list to a column of the data frame and values to the valueCol
list2df = function(l, cols = names(l), row.name = NULL, names2col = NULL, valueCol = 'value') {
	idcs = if (is.null(cols)) 1:length(l) else
		if (all(is.integer(cols))) cols else which.indeces(names(l), cols);
	if (is.null(cols) || all(is.integer(cols))) cols = paste('C', 1:length(l), sep = '');
	r = as.list(rep(NA, length(cols)));
	names(r) = cols;
	r[idcs] = l;
	r = as.data.frame(r, stringsAsFactors = FALSE);
	if (!is.null(row.name)) row.names(r)[1] = row.name;
	if (!is.null(names2col)) {
		r = data.frame(name = names(r), value = unlist(r[1, ]), row.names = NULL, stringsAsFactors = FALSE);
		names(r) = c(names2col, valueCol);
	}
	r
}

be.numeric = function(v)
	sapply(v, function(e)grepl('^-?\\d*(\\.\\d+)?(e-?\\d+)?$', e, ignore.case = TRUE, perl = TRUE));

list2df.print = function(l, valueCol = 'value', names2col = NULL, ..., digits = 3, scientific = 3) {
	l1 = list2df(l, valueCol = valueCol, names2col = names2col, ...);
	numericRows = be.numeric(l1[[valueCol]]);
	numbers = as.numeric(l1[[valueCol]][numericRows]);
	log10range = max(floor(log10(numbers))) - min(floor(log10(numbers)));
	#fmt = if (log10range > digits + 1) '%.*e' else '%.*f';
	numbers = sprintf(ifelse(abs(floor(log10(numbers))) > scientific, '%.*e', '%.*f'), digits, numbers);
	#numbers = sapply(numbers, function(n)sprintf(fmt, digits, n));
	separators = as.vector(names(l) == '' & is.na(l));
	l1[separators, names2col] = '-';
	l1[separators, valueCol] = '';
	l1[numericRows, valueCol] = numbers;
	print(l1);
}


rbind.list2df = function(d, l, row.name = NULL) {
	d = as.data.frame(d);
	r = list2df(l, names(d), row.name);
	r0 = rbind(d, r);
	r0
}

# take list of lists
#	names of list elements become column-names
listOfLists2df = function(l, columnNames = names(l[[1]])) {
	colV = lapply(columnNames, function(n)Df_(list.kp(l, n, do.unlist = TRUE)));
	r = Df_(do.call(cbind, colV), names = columnNames);
	r
}

ListOfLists2df_extract = function(l, kp, template) {
	l1 = list.kp(l, kp, null2na = TRUE, do.unlist = FALSE, template = template);
	do.call(rbind, l1);
}
# advanced version of the above
ListOfLists2df = function(l,
	keyPath = '*', columnNames = names(list.kp(l[1], keyPath)[[1]]),
	reverseKeys = FALSE, keySep = '-', template = NA) {
	colV = lapply(columnNames, function (n) {
		kp = Sprintf('%{keyPath}s$%{n}s');
		# <A> robustly choose name (assume first element is proper template)
		#name = if (collapse) names(ListOfLists2df_extract(l[1], kp, template, collapse)) else NULL;
		r = ListOfLists2df_extract(l, kp, template);
		# names
		kpk = list.kp.keys(Sprintf('%{n}s'));
		cns = Col.names(r, vivify = FALSE);
		if (is.null(cns)) keySep = '';
		ns = if (reverseKeys)
			paste(cns, join(rev(kpk), keySep), sep = keySep) else
			paste(join(kpk, keySep), cns, sep = keySep);
		r = SetNames(r, ns);
		r
	});
	r = do.call(cbind, colV);
	# <!> Df_ applies as.data.frame -> normalization of column names
	#r = if (collapse == 0) Df_(r0, names = columnNames) else r0;
	r
}

# select elements from iterable (<=> perl grep)
Select = function(l, f, ...)return(l[sapply(l, f, ...)]);

# # d: data frame, l: list with names corresponding to cols, values to be searched for in columns
searchDataFrame = function(d, l, .remove.factors = TRUE) {
	ns = names(l);
	d = d[, ns, drop = FALSE];
	if (.remove.factors) {
		l = sapply(l, function(e)ifelse(is.factor(e), levels(e)[e], e));
		#d = apply(d, 2, function(col)(if (is.factor(col)) levels(col)[col] else col));
	}
	rs = which(as.vector(apply(apply(d, 1, function(r)(r == l)), 2, all)));
	rs
}

Which.cols = function(d, cols, regex = FALSE) {
	which.indeces(cols[is.character(cols)], names(d), regex = regex)
}
.df.cols = which.cols = function(d, cols, regex = FALSE) {
	cols[is.numeric(cols)] = as.integer(cols[is.numeric(cols)]);
	cols[is.character(cols)] = which.indeces(cols[is.character(cols)], names(d), regex = regex);
	as.integer(cols)
}
DfColsAfter = function(d, col, regex = FALSE, inclusive = FALSE) {
	I = which.cols(d, col, regex = regex);
	if (is.null(I)) return(NULL);
	return( (I + ifelse(inclusive, 0, 1)):ncol(d));
}
DfColsBetween = function(d, col1, col2, regex = FALSE, inclusive = TRUE) {
	I1 = if (missing(col1)) 1 else which.cols(d, col1, regex = regex);
	I2 = if (missing(col2)) ncol(d) else which.cols(d, col2, regex = regex);
	if (is.null(I1) || is.null(I2)) return(NULL);
	return( (I1 + ifelse(inclusive, 0, 1)):(I2 + ifelse(inclusive, 0, -1)) );
}

# select columns by name
.df = function(d, names, regex = TRUE, as.matrix = FALSE) {
	cols = which.indeces(names, names(d), regex = regex);
	d0 = d[, cols, drop = FALSE];
	# <t> simpler version:
	# d0 = d[, .df.cols(d, names, regex)];
	if (as.matrix) d0 = as.matrix(d0);
	d0
}
.df.reorder = function(d, names, regex = TRUE) {
	cols = .df.cols(d, names, regex);
	d0 = d[, c(cols, setdiff(1:dim(d)[2], cols))];
	d0
}
# remove columns by name
.dfm = function(d, names, regex = FALSE, as.matrix = FALSE) {
	cols = if (all(is.numeric(names))) as.integer(names) else which.indeces(names, names(d), regex = regex);
	d0 = d[, -cols, drop = FALSE];
	if (as.matrix) d0 = as.matrix(d0);
	d0
}
# remove rows by name
.dfrmr = function(d, names, regex = FALSE, as.matrix = FALSE) {
	rows = if (all(is.numeric(names)))
		as.integer(names) else
		which.indeces(names, row.names(d), regex = regex);
	d0 = d[-rows, , drop = FALSE];
	if (as.matrix) d0 = as.matrix(d0);
	d0
}

# remove rows/columns by name
.dfrm = function(d, rows = NULL, cols = NULL, regex = FALSE, as.matrix = FALSE) {
	d = as.data.frame(d);	# enforce data frame
	rows = if (is.null(rows)) 1:dim(d)[1] else
		-(if (all(is.numeric(rows))) as.integer(rows) else which.indeces(rows, row.names(d), regex = regex));
	cols = if (is.null(cols)) 1:dim(d)[2] else 
		-(if (all(is.numeric(cols))) as.integer(cols) else which.indeces(cols, names(d), regex = regex));
	d0 = d[rows, cols, drop = FALSE];
	if (as.matrix) d0 = as.matrix(d0);
	d0
}

# alignByRowNames: logical: use row.names from first element, else use provided vector
Cbind = function(..., stringsAsFactors = FALSE, deparse.level = 0, alignByRowNames = NULL) {
	l = list(...);
	if (notE(alignByRowNames)) {
		if (is.null(row.names(l[[1]]))) stop('Cbind[alignByRowNames]: No row names @ 1');
		ref = if (is.logical(alignByRowNames) && alignByRowNames)
			row.names(l[[1]]) else
			alignByRowNames;
		l = pairslapply(l, seq_along(l), function(e, i) {
			if (is.null(row.names(e))) stop('Cbind[alignByRowNames]: No row names @ %{i}d');
			e[order_align(ref, row.names(e)), , drop = FALSE]
		});
	}
	if (length(l) == 1)
		# <p> special case vector
		t_(l[[1]]) else
		# <p> standard invocation
		do.call(cbind, c(l, list(deparse.level = deparse.level)))
}
Rbind = function(..., stringsAsFactors = FALSE) {
	l = list(...);
	r = if (length(l) == 1) t_(t_(l[[1]])) else {
		if (class(l[[1]]) == 'data.frame')
			rbind(..., stringsAsFactors = stringsAsFactors) else
			rbind(...);
	}
	r
}


subsetTop = function(obj, sel, N = 1) {
	d0 = subset(obj, sel);
	d1 = d0[1:min(nrow(d0), N), ];
	d1
}

# transpose to create column vector for vector
t_ = function(m)(if (is.vector(m)) t(t(m)) else t(m))
# double transpose aka transpose to row -> vector to 1 x N matrix, otherwise identity
t2r = function(m)t(t_(m))

# convert strings to data frame names
#	<i> create a data frame and extract names
.dfns = function(ns)gsub(':', '.', ns);

# manipulate list of vectors
# vectors i = 1,.., n with entries v_ij are represented as vector v_11, ..., v_n1, v_21, ...
vector.intercalate = meshVectors = function(...) {
	l = list(...);
	if (length(l) == 1) l = l[[1]];
	v = as.vector(t(sapply(l, function(v)unlist(v))));
	# <N> preferred implementation
	# No unlist -> should be part of input sanitization
	# v = as.vector(do.call(rbind, l));
	v
}

is.sorted = function(...)(!is.unsorted(...))
is.ascending = function(v) {
	if (length(v) < 2) return(TRUE);
	for (i in 2:length(v)) if (v[i] <= v[i - 1]) return(FALSE);
	return(TRUE);
}

# pad a vector to length N
pad = function(v, N, value = NA)c(v, rep(value, N - length(v)));

#
#	<par> number sequences
#

rep.each.vector = function(v, n)as.vector(matrix(rep(v, n), n, byrow = TRUE))
rep.each = function(l, n, simplify = unlist) {
	l = Avu(l);
	if (length(n) == 1) rep.each.vector(l, n) else simplify(pairsapply(l, n, rep))
}
factorWithLevels = function(f, levels_) {
	f = as.factor(f);
	levels(f) = levels_;
	f
}
Rep.each = function(v, n) {
	r = rep.each(v, n);
	return(if (is.factor(v)) factorWithLevels(r, levels(v)) else r)
}
copyFactorStructure = function(dS, dD) {
	factors = which(lapply(dS, class) == 'factor');
	for (f in factors) dD[[f]] = factorWithLevels(dD[[f]], levels(dS[[f]]));
	dD
}
rep.each.row = function(m, n) {
# 	r = matrix(rep.each(m, n), ncol = ncol(m));
# 	if (class(m) == 'data.frame') {
# 		r = Df_(r, names = names(m));
# 		r = copyFactorStructure(m, r);
# 	}
	r = if (is.data.frame(m))
		Df_(lapply(m, Rep.each, n = n)) else
		m[rep.each(Seq(1, nrow(m)), n), , drop = FALSE]
	r
}

#rep.list = function(l, n) lapply(1:length(l), function(e)l);
# <!> changed as of 23.8.2016; n not used before
rep.each.list = rep.list = function(l, n) lapply(1:n, function(e)l);

matrix.intercalate = function(..., direction = 1, listOfMatrices = FALSE) {
	l = list(...);
	if (listOfMatrices) l = l[[1]];
	# <!> assume same dimension
	d = dim(l[[1]]);
	N = prod(d);
	# <p> create new matrix
	v = c(if (direction == 1) sapply(l, as.vector) else sapply(sapply(l, t), as.vector), recursive = TRUE);
	vN = as.vector(matrix(v, ncol = N, byrow = TRUE));
	r = if (direction == 1)
		matrix(vN, nrow = d[1] * length(l)) else
		matrix(vN, ncol = d[2] * length(l), byrow = TRUE);
	# <p> return value
	if (class(l[[1]]) == 'data.frame') r = Df_(r, names = names(l[[1]]));
	r
}

matrixSearch = function(mSearch, mSearched, cols = 1:ncol(mSearch)) {
	df1 = Df_(mSearch, names = paste0('c', cols));
	df2 = Df_(mSearched[, cols, drop = FALSE], names = paste0('c', cols));
	return(DfSearch(df1, df2, returnIdcs = TRUE));
}

arrayFromRowPairs = function(m, halves = FALSE) {
	if (halves)
		aperm(array(t(m), dim = c(2, dim(m)[1]/2, dim(m)[2])), c(2, 1, 3)) else
		# adjecent pairs
		aperm(array(t(m), dim = c(2, dim(m)[2], dim(m)[1]/2)), c(3, 1, 2))
}

data.frame.expandWeigths = function(data, weights = 'weights') {
	w = data[[weights]];
	weightsCol = which(names(data) == weights);
	df0 = lapply(1:length(w), function(i) {
		if (w[i] > 0) rep.each.row(data[i, -weightsCol], w[i]) else list();
	});
	df1 = rbindDataFrames(df0);
	df1
}

# spread/fill vector to indeces
vector.spread = function(v, idcs, N, default = 0) {
	r = rep(default, N);
	r[idcs] = v;
	r
}

# search vector for value, fill value elements with elements prior to it
#	e.g. 1, NA, NA, 2, NA -> 1, 1, 1, 2, 2
vector.propagateValuesForward = function(v, value = NA, vs) {
	idcs = if (is.na(value)) which(!is.na(v)) else which(v != value);
	Idcs = c(idcs, length(v) + 1);	# padded version
	# assign positions
	iA = lapply(seq_along(idcs), function(i)Seq(idcs[i] + 1, Idcs[i + 1] - 1));
	# indeces of values tb assigned
	iV = lapply(seq_along(idcs), function(i)rep(idcs[i], length(iA[[i]])));
	# fill in values
	#r = vector.assign(v, unlist(iA), v[unlist(iV)]);
	v[unlist(iA)] = v[unlist(iV)];
	return(v);
}

# create new vector with length == length(v) + length(idcs)
# idcs are positions in the final vector
vector.embed = function(v, idcs, e, idcsResult = TRUE) {
	if (!idcsResult) idcs = idcs + (1:length(idcs)) - 1;
	N = length(v) + length(idcs);
	r = rep(NA, N);
	r[setdiff(1:N, idcs)] = v;
	r[idcs] = e;

	# <p> names
	if (!is.null(names(v)) || !is.null(names(e))) {
		ns = rep(NA, N);
		if (!is.null(names(v))) ns[setdiff(1:N, idcs)] = names(v);
		if (!is.null(names(e))) ns[idcs] = names(e);
		names(r) = ns;
	}
	r
}
# set values at idcs
vector.assign = function(v, idcs, e, na.rm = 0, N) {
	if (!missing(N)) v = rep(v, N);
	v[idcs] = e;
	if (!is.na(na.rm)) v[is.na(v)] = na.rm;
	v
}
# names based assignment
Vector.assign = function(v, e, na.rm = NA) {
	idcs = which.indeces(names(e), names(v));
	vector.assign(v, idcs, e, na.rm = na.rm)
}

matrix.assign = function(m, idcs, e, byrow = TRUE) {
	if (length(dim(idcs)) > 1) {
		m[as.matrix(idcs)] = e
	} else if (byrow)
		m[idcs, ] = e else
		m[, idcs] = e
	m
}

# extract elements from array/matrix indexed in a row-wise manner by ...
#	array.extract(m, c(1, 2), c(1, 2)) -> c(m[1, 1], m[2, 2])
array.extract = function(a, ...) {
	r = mapply(function(...)do.call('[', c(list(a), list(...))), ...);
	return(r);
}

# are columns/rows same values in matrix
matrix.same = function(m, direction = 1) {
	apply(m, direction, function(e)all(e[1] == e))
}

vectorIdcs = function(v, f, ..., not = FALSE) {
	r = sapply(v, f, ...);
	which(if (not) !r else r)
}

is.seq = function(v, offset = 1)all( (v - offset + 1) == seq_along(v))

# produce indeces for indeces positioned into blocks of blocksize of which count units exists
# example: expand.block(2, 10, 1:2) == c(1, 2, 11, 12)
expand.block = function(count, blocksize, indeces) {
	blks = Seq(1,count);
	if (is.null(blks)) return(NULL);
	as.vector(apply(to.col(blks), 1,
		function(i){ (i - 1) * blocksize + t(to.col(indeces)) }
	));
}

search.block = function(l, s) {
	b.sz = length(s);
	which(sapply(
		1:(length(l)/b.sz), function(i){all(l[((i - 1) * b.sz + 1):(i * b.sz)] == s)}
	));
}

#
#	<par> matrix functions
#

# <!> assumes same indeces for rows/columns
matrixFromIndexedDf = function(df, idx.r = 'idx.r', idx.c = 'idx.c', value = 'value', referenceOrder = NULL) {
	id = unique(c(df[[idx.r]], df[[idx.c]]));
	# matrix indeces
	# <A> canonical order is by repeating vector id for row index, constant for columns within repetition
	#	-> matrix filled by columns
	midcs = merge(data.frame(id = id), data.frame(id = id), by = NULL);
	midcs = data.frame(midcs, mfid.i = 1:nrow(midcs));
	map = merge(df[, c(idx.r, idx.c, value)], midcs,
		by.x = c(idx.r, idx.c), by.y = c('id.x', 'id.y'), all.y = TRUE);
	# return to midcs order
	map = map[order(map$mfid.i), ];
	# filled by rows
	m = matrix(map[[value]], nrow = length(id));
	# reorder matrix
	o = order_align(firstDef(referenceOrder, id), id);
	# reorder in two steps -> out of mem otherwise
	m1 = m[o, ];
	m2 = m1[, o];
	m2
}

symmetrizeMatrix = function(m) {
	m[is.na(m)] = t(m)[is.na(m)];
	m
}

which.row = function(m, row) {
	cols = names(as.list(row));
	if (is.null(cols)) cols = 1:length(row);
	rows = 1:(dim(m)[1]);
	rows.found = rows[sapply(rows, function(i){ all(m[i, cols] == row) })];
	rows.found
}

# lsee:	list with searchees
# lsed:	list with searched objects
# inverse: lsed are regexes matched against lsee; pre-condition: length(lsee) == 1
# ret.list: for match.multi return list by lsee
# <!><t> cave: semantics changed as of 17.8.2009: return NA entries for unfound lsee-entries
# <!> match multi only implemented for merge = TRUE
which.indeces = function(lsee, lsed, regex = FALSE, ret.na = FALSE, merge = TRUE, match.multi = FALSE, ...,
	inverse = FALSE, ret.list = FALSE) {
	if (!length(lsed) || !length(lsee)) return(c());
	v = if (is.list(lsed)) names(lsed) else lsed;
	idcs = if (regex) {
		which(sapply(lsed, function(e)(
			if (inverse) length(fetchRegexpr(e, lsee, ...)) > 0 else
				any(sapply(lsee, function(see)(length(fetchRegexpr(see, e, ...)) > 0)))
		)))
	} else if (merge) {
		d0 = merge(
			data.frame(d = lsed, ix = 1:length(lsed)),
			data.frame(d = lsee, iy = 1:length(lsee)), all.y = TRUE);
		d0 = d0[order(d0$iy), ];
		idcs = if (match.multi) {
				#d0$ix[unlist(sapply(lsee, function(e)which(d0$d == e)))]
				#na.omit(sort(d0$ix))
				r = if (ret.list)
					unlist.n(by(d0, d0$d, function(d)list(na.omit(d$ix)), simplify = FALSE)) else
					na.omit(d0$ix);
				r
			} else {
				d0$ix[pop(which(c(d0$iy, 0) - c(0, d0$iy) != 0))];
			}
		# less efficient version
#		} else d0$ix[unlist(sapply(lsee, function(e)which(d0$d == e)[1]))];
#		} else d0$ix[order(d0$iy)]
		if (!ret.na) idcs = idcs[!is.na(idcs)];
		idcs
	} else {
		unlist(as.vector(sapply(lsee, function(e){
			w = which(e == v);
			if (!ret.na) return(w);
			ifelse(length(w), w, NA)
		})))
	};
	r = if (ret.list) idcs else as.integer(idcs);
	r
}

grep.vector = function(lsee, lsed, regex = FALSE, ret.na = FALSE, merge = TRUE, match.multi = FALSE, ..., inverse = FALSE) {
	lsed[which.indeces(lsee, lsed, regex, ret.na, merge, match.multi, ..., inverse = inverse)]
}
grep.infixes = function(lsee, lsed, ...) {
	r = grep.vector(sapply(lsee, function(v)sprintf('^%s.*', v)), lsed, regex = TRUE, inverse = FALSE, ... );
	r
}

# force structure to be matrix (arrange vector into a row)
MR = function(m) {
	if (!is.matrix(m)) m = matrix(m, byrow = TRUE, ncol = length(m));
	m
}
# force structure to be matrix (arrange vector into a columns)
MC = function(m) {
	if (!is.matrix(m)) m = matrix(m, byrow = FALSE, nrow = length(m));
	m
}

#
#	<par> data processing
#

# like table but produce columns for all numbers 1..n (not only for counts > 0)
# cats are the expected categories
table.n = function(v, n, min = 1, categories = NULL, useNA = 'no') {
	if (is.null(categories)) categories = min:n;
	t = as.vector(table(c(categories, v), useNA = useNA) - rep(1, length(categories)));
	t
}

tableFreqMarg = function(tab, margin = 2)apply(as.matrix(tab), margin, vn)

table.freq = function(v, byCol = TRUE, useNA = 'no') {
	t0 = table(v, useNA = useNA);
	r = if (is.vector(v) || is.factor(v) || is.numeric(v) || ncol(v) == 1) { t0 / sum(t0) } else {
		if (byCol) tableFreqMarg(t0) else tableFreqMarg(t0, 1)
	}
	r
}
table.n.freq = function(...) {
	t0 = table.n(...);
	r = t0 / sum(t0);
	r
}
table2df = function(tab, perc = FALSE, total = FALSE) {
	df0 = Df_(tab);
	nms = setdiff(names(df0), 'Freq');
	f = do.call(formulaWith, as.list(nms));
	df1 = dcast(df0, f, value.var= 'Freq');
	tot = apply(df1[, -1], 1, sum, na.rm = T);
	df2 = Df_(df1);
	if (perc) df2 = cbind(df2, Df_(df1[, -1] / tot, names = paste0(names(df1)[-1], 'perc')));
	if (total) df2$Total = tot;
	return(df2);
}
Table = function(v, min, max, ..., cats, asDf = FALSE, perc = FALSE, total = FALSE) {
	if (missing(min) && missing(max) && missing(cats)) {
		t0 = table(v);
		return(if (asDf) table2df(t0, perc, total) else t0);
	}
	if (!missing(cats)) {
		d = Df_(lapply(v, Avu));
		catsV = SetNames(Df_(merge.multi.list(cats)), names(d));
		t0 = table(rbind(d, catsV)) - 1;
		return(if (asDf) table2df(t0) else t0);
	} else {
		if (missing(min)) min = min(v);
		if (missing(max)) max = max(v);
		t0 = table.n(v, n = max, min = min);
		return(if (asDf) table2df(t0, perc, total) else t0);
	}
}
TableDf = function(v, min, max, ..., cats, asDf = TRUE, perc = FALSE, total = FALSE)
	Table(v, min, max, ..., cats = cats, asDf = asDf, perc = perc, total = total);
v2freq = function(v)(v/sum(v))

TableReshaped = function(df, reshape = names(df), replace = list(Freq. = '')) {
	tab = table(df)
	tabReshaped = reshape.wide(Df_(tab), reshape[1], reshape[-1]);
	if (notE(replace)) names(tabReshaped) = vector.replace(names(tabReshaped), list(Freq. = ''), regex = T);
	return(tabReshaped);
}


#
#	<p> numeric function
#

to.numeric = function(x) { SetNames(suppressWarnings(as.numeric(x)), names(x)) }
minFloor = function(x)(x - floor(x))

#
#	<par> data types
#


# set types for columns: numeric: as.numeric
data.frame.types = function(df, numeric = c(), character = c(), factor = c(), integer = c(),
	do.unlist = TRUE, names = NULL, row.names = NULL, reset.row.names = FALSE, do.rbind = FALSE, do.transpose = FALSE,
	stringsAsFactors = FALSE) {
	if (do.rbind) {
		#old code: df = t(sapply(df, function(e)e));
		lengthes = sapply(df, length);
		maxL = max(lengthes);
		df = t(sapply(1:length(df), function(i)c(df[[i]], rep(NA, maxL - lengthes[i]))));
	}
	if (do.transpose) df = t(df);
	df = as.data.frame(df, stringsAsFactors = stringsAsFactors);
	# set or replace column names
	if (!is.null(names)) {
		if (class(names) == "character") names(df)[1:length(names)] = names;
		if (class(names) == "list") names(df) = vector.replace(names(df), names);
	}
	if (do.unlist) for (n in names(df)) { df[[n]] = unlist(df[[n]]); }
	for (n in numeric) { df[[n]] = as.numeric(df[[n]]); }
	for (n in integer) { df[[n]] = as.integer(df[[n]]); }
	for (n in character) { df[[n]] = as.character(df[[n]]); }
	for (n in factor) { df[[n]] = as.factor(df[[n]]); }
	if (reset.row.names) row.names(df) = NULL;
	if (length(row.names) > 0) row.names(df) = row.names;
	df
}

DfStack = function(df0, N)do.call(rbind, rep.list(df0, N));
DfClasses = function(dataFrame)nlapply(dataFrame, function(n)class(dataFrame[[n]]));
DfAsInteger = function(dataFrame, as_integer) {
	#dfn = apply(dataFrame[, as_integer, drop = FALSE], 2, function(col)as.integer(avu(col)));
	# <!> 6.6.2016 as.integer first needed to retain factor status on factors
	dfn = nlapply(as_integer, function(col)avu(as.integer(dataFrame[[col]])));
	dataFrame[, as_integer] = as.data.frame(do.call(cbind, dfn));
	dataFrame
}
DfAsLogical = function(dataFrame, as_logical) {
	dfn = nlapply(as_logical, function(n) {
		col = dataFrame[[n]];
		if (is.factor(col)) (col == levels(col)[1]) else avu(as.logical(col));
	});
	dataFrame[, as_logical] = as.data.frame(do.call(cbind, dfn));
	dataFrame
}
DfAsCharacter = function(dataFrame, as_character) {
	#dfn = apply(dataFrame[, as_character, drop = FALSE], 2, function(col)as.character(avu(col)));
	#dataFrame[, as_character] = as.data.frame(dfn, stringsAsFactors = FALSE);
	dfn = nlapply(as_character, function(col)avu(as.character(dataFrame[[col]])));
	dataFrame[, as_character] = as.data.frame(do.call(cbind, dfn), stringsAsFactors = FALSE);
	dataFrame
}
DfFac2num = function(dataFrame) {
	return(do.call(data.frame, lapply(dataFrame, function(e)if (is.factor(e)) as.numeric(e) else e)))
}
DfApplyValueMap = function(r, valueMap, Df_doTrimValues = FALSE,
	Df_mapping_value = '__df_mapping_value__',
	Df_mapping_empty = '__DF_EMPTY__', Do_Df_mapping_empty = TRUE) {

	for (n in names(valueMap)) {
		vs = if (Df_doTrimValues && class(r[[n]]) %in% c('character', 'factor'))
			nina(trimString(as.character(r[[n]])), Df_mapping_value) else
			as.character(r[[n]]);
		vm = valueMap[[n]];
		if (Do_Df_mapping_empty) {
			vs = ifelse(nit(vs == ''), Df_mapping_empty, vs);
			if (!(Df_mapping_empty %in% names(vm)))
				vm = c(vm, listKeyValue(Df_mapping_empty, NA));
		}
		vs = nina(valueMap[[n]][vs], Df_mapping_value);
		vs = ifelse(vs == Df_mapping_value, as.character(r[[n]]), vs);
		r[[n]] = vs;
	}
	return(r);
}
# copy over factor structure from other data frame (tentamen/bw for example)
DfEmbed = function(d, dSource) {
	cols = nlapply(d, function(n) {
		if (class(dSource[[n]]) == 'factor')factor(d[[n]], levels(dSource[[n]]))else d[[n]]
	})
	return(Df_(cols));
}

#	r = strptime(col, dateFormat[1], tz = firstDef(dateFormat[2], defaultTz));

DfDate = function(dataFrame, as_date, format = '%F', tz = 'UTC') {
	dfn = nlapply(as_date, function(n)Df_(strptime(dataFrame[[n]], format, tz), names = n));
	dataFrame[, as_date] = do.call(cbind, dfn);
	dataFrame
}

# as of 22.7.2013 <!>: min_ applied before names/headerMap
# as of 19.12.2013 <!>: as.numeric -> as_numeric
# as of 22.5.2014 <!>: t -> t_
# as of 13.11.2014 <!>: sapply -> simplify_
# Create data frames with more options than \code{data.frame}
Df_ = function(df0, headerMap = NULL, names = NULL, min_ = NULL,
	as_numeric = NULL, as_character = NULL, as_factor = NULL, as_integer = NULL, as_logical = NULL,
	row.names = NA, valueMap = NULL, Df_as_is = TRUE, simplify_ = FALSE,
	deep_simplify_ = FALSE, t_ = FALSE, unlist_cols = FALSE, transf_log = NULL, transf_m1 = NULL,
	Df_doTrimValues = FALSE, Df_mapping_value = '__df_mapping_value__',
	Df_mapping_empty = '__DF_EMPTY__', Do_Df_mapping_empty = TRUE, apply_ = FALSE,
	as_date = NULL, date_format = '%F', date_tz = 'UTC') {
	# <p> input sanitation
	#r = as.data.frame(df0);
	# for a vector with identical names for each entry, use this as a column name
	if (length(unique(names(df0))) == 1 && !Nif(names)) names = unique(names(df0));
	# sanitize row.names
	dn = dimnames(df0);
	if (Nif(dn) && any(duplicated(dn[[1]]))) dimnames(df0)[[1]] = NULL;
	# <!> commented out on 4.4.2019, test implemented to fix this behavior
	#if (length(row.names) == 0 || !all(is.na(row.names))) base::row.names(df0) = row.names;

	if (apply_) df0 = as.data.frame(apply(df0, 2, identity));
	#if (!Nif(Apply_)) df0 = as.data.frame(apply(df0, 2, Apply_));
	if (t_) df0 = t(df0);
	# reset_row_names breaks unit tests (27.9.2017)
	#r = data.frame(df0, stringsAsFactors = !Df_as_is, row.names = if (reset_row_names) NA else NULL);
	r = data.frame(df0, stringsAsFactors = !Df_as_is);
	if (notE(min_)) {
		is = which.indeces(min_, names(r));
		if (length(is) > 0) r = r[, -is, drop = FALSE];
	}
	if (simplify_) r = as.data.frame(sapply(r, identity));
	if (deep_simplify_) r = as.data.frame(
		nlapply(r, function(col)sapply(r[[col]], unlist)), stringsAsFactors = !Df_as_is
	);

	#
	#	<p> column names
	#
	if (notE(names)) {
		if (class(names) == 'character') names(r)[1:length(names)] = names;
		if (class(names) == 'list') names(r) = vector.replace(names(r), names);
	}
	if (notE(headerMap)) names(r) = vector.replace(names(r), headerMap);
	#
	#	<p> column types
	#
#if (class(df0) == 'data.frame' && ncol(df0) >= 3) browser();
	if (notE(as_numeric)) {
		#dfn = apply(r[, as_numeric, drop = FALSE], 2, function(col)as.numeric(avu(col)));
		dfn = lapply(r[, as_numeric, drop = FALSE], function(col)avu(as.numeric(col)));
		r[, as_numeric] = as.data.frame(do.call(cbind, dfn));
	}
	if (notE(as_logical)) r = DfAsLogical(r, as_logical);
	if (notE(as_integer)) r = DfAsInteger(r, as_integer);
	if (notE(as_character)) r = DfAsCharacter(r, as_character);
	if (notE(as_factor)) {
		# <N> does not work
		#dfn = apply(r[, as_factor, drop = FALSE], 2, function(col)as.factor(col));
		#r[, as_factor] = dfn;
		for (f in as_factor) r[, f] = as.factor(r[[f]]);
	}
	if (notE(as_date)) r = DfDate(r, as_date, date_format, date_tz);
	#
	#	<p> value map
	#
	if (notE(valueMap)) {
# 		for (n in names(valueMap)) {
# 			vs = if (Df_doTrimValues && class(r[[n]]) %in% c('character', 'factor'))
# 				nina(trimString(as.character(r[[n]])), Df_mapping_value) else
# 				as.character(r[[n]]);
# 			vm = valueMap[[n]];
# 			if (Do_Df_mapping_empty) {
# 				vs = ifelse(nit(vs == ''), Df_mapping_empty, vs);
# 				if (!(Df_mapping_empty %in% names(vm)))
# 					vm = c(vm, listKeyValue(Df_mapping_empty, NA));
# 			}
# 			vs = nina(valueMap[[n]][vs], Df_mapping_value);
# 			vs = ifelse(vs == Df_mapping_value, as.character(r[[n]]), vs);
# 			r[[n]] = vs;
# 		}
		r = DfApplyValueMap(r, valueMap,
			Df_doTrimValues, Df_mapping_value, Df_mapping_empty, Do_Df_mapping_empty);
	}
	#
	#	<p> transformations
	#
	if (notE(transf_log)) r[, transf_log] = log(r[, transf_log, drop = FALSE]);
	if (notE(transf_m1)) r[, transf_m1] = r[, transf_m1, drop = FALSE] - 1;
	if (length(row.names) == 0 || !all(is.na(row.names))) base::row.names(r) = row.names;
	if (unlist_cols) for (n in names(r)) r[[n]] = avu(r[[n]]);
	r
}

Df = function(..., headerMap = NULL, names = NULL, min_ = NULL, row.names = NA, Df_as_is = TRUE,
	as_numeric = NULL, as_character = NULL, as_factor = NULL, t_ = FALSE, unlist_cols = FALSE) {
	r = data.frame(...);
	Df_(r, headerMap = headerMap, names = names, min_ = min_, row.names = row.names,
		as_numeric = as_numeric,
		as_character = as_character,
		as_factor = as_factor,
		Df_as_is = Df_as_is,
		t_ = t_,
		unlist_cols = unlist_cols
	);
}
Df2list = function(df) {
	df = as.data.frame(df);
	nlapply(names(df), function(n)df[[n]]);
}
Dfselect = function(data, l, na.rm = nif) {
	sel = apply(sapply(nlapply(l, function(n)data[[n]] == l[[n]]), identity), 1, all);
	r = data[na.rm(sel), ];
	r
}
DfSearch = function(dfSearch, dfSearched,
	colNamesReset = 'col', colNameIdx = '.dfSearchIdx', returnIdcs = FALSE) {

	if (is.null(dfSearched)) return(NULL);
	nms = if (notE(colNamesReset)) {
		nms = paste(colNamesReset, 1:ncol(dfSearched), sep = '');
		names(dfSearch) = names(dfSearched) = nms;
	} else names(dfSearched);
	dfm = merge(
		Df(1:nrow(dfSearched), dfSearched, names = colNameIdx),
		Df(1:nrow(dfSearch), dfSearch, names = colNameIdx), by = nms);
	if (returnIdcs)
		return(dfm[, paste0(colNameIdx, c('.x', '.y')), drop = FALSE]) else
		return(dfm[[paste0(colNameIdx, '.x')]]);
}

DfDiff = function(d1, d2) {
	dC = rbind(d2, d1);
	row.names(dC) = NULL;
	dCu = unique(dC);
	# d2 comes first, non-unique rows left out from d1, sames as ones diffed out
	r = if (nrow(dCu) == nrow(d2)) dCu[c(), ] else dCu[(nrow(d2) + 1):nrow(dCu), , drop = FALSE];
	r
}
# replace columns in data.frame
DfRepl = function(d0, d1) {
	d0[, names(d1)] = d1;
	return(d0);
}

DfRound = function(df0, cols_round = names(df0), digits = 2, as_numeric = FALSE) {
	rounder = if (as_numeric)
		function(col)round(as.numeric(df0[[col]]), digits) else
		function(col)round(df0[[col]], digits)
	df0[, cols_round] = do.call(cbind, lapply(cols_round, rounder));
	df0
}


# standardize df names using formulas
# f: formula with names used in the Dataframe ~ x1 + x2 ... or x1 ~ x2 + ...
#	corresponding positionally to the standard names
# nmsStd: formula or vector with standard names
# d: data.frame with column names tb transformed
dfNmsStd = function(f, nmsStd, d) {
	nmsUsed = all.vars(f);
	#if (is.formula(nmsStd)) nmsStd = all.vars(nmsStd);
	# import from plyr (is.formula) leads to test failures <!>
	if (class(nmsStd) == 'formula') nmsStd = all.vars(nmsStd);
	if (length(nmsUsed) != length(nmsStd))
		stop(Sprintf('Formula names [%{f}s] do not match standard names [%{nm}s]',
			f = formula.to.character(f), nm = join(nmsStd, ', ')));
	d1 = Df_(d, headerMap = listKeyValue(nmsUsed, nmsStd));
	return(d1);
}
# DfNames2std = function(d, nmsFormula, nmsStandard) {
# 	d1 = Df_(d, headerMap = listKeyValue(all.vars(nmsFormula), nmsStandard));
# 	d1
# }
DfNames2std = function(d, nmsFormula, nmsStandard)dfNmsStd(nmsFormula, nmsStandard, d)

charRange = characterRange = function(ns, range, indeces = TRUE, invert = FALSE) {
	N = length(ns);
	r = if (class(range) == 'character') {
		(if (is.na(range)[1])1 else which(range[1] == ns)) :
			(if (is.na(range)[2])N else which(range[2] == ns))
	} else if (class(range) == 'integer') {
		(if (is.na(range)[1])1 else range[1]) :
			(if (is.na(range)[2])N else range[2])
	} else c();
	if (invert) r = setdiff(1:length(ns), r);
	if (!indeces) r = ns[r];
	return(r);
}

DfCol = function(d, range) {
	d = d[, characterRange(names(d), range), drop = F];
	return(d);
}

List_ = .List = function(l, min_ = NULL, sel_ = NULL,
	rm.null = FALSE, names_ = NULL, null2na = FALSE, simplify_ = FALSE, rm.na = FALSE) {
	if (!is.null(min_)) {
		i = which.indeces(min_, names(l));
		if (length(i) > 0) l = l[-i];
	}
	if (!is.null(sel_)) {
		i = which.indeces(sel_, names(l));
		if (length(i) > 0) l = l[i];
	}
	if (rm.null) {
		remove = -which(sapply(l, is.null));
		if (length(remove) > 0) l = l[remove];
	}
	if (null2na) {
		nullI = which(sapply(l, is.null));
		l[nullI] = NA;
	}
	if (rm.na) {
		l = l[!is.na(l)];
	}
	if (notE(names_)) {
		if (is.character(names_)) names(l)[Seq(1, length(names_))] = names_;
		if (is.list(names_)) names(l) = vector.replace(names(l), names_);
		if (is.na(names_)) names(l) = NULL;
	}
	if (simplify_) l = sapply(l, identity);
	l
}
List = function(..., min_ = NULL, envir = parent.frame(), names_ = NULL) {
	l = eval(list(...), envir = envir);
	.List(l, min_ = min_, names_ = names_);
}

Unlist = function(l, ..., null2na_ = FALSE) {
	if (null2na_) l[sapply(l, is.null)] = NA;
	unlist(l, ...)
}

#last = function(v)(rev(v)[1])
last = function(v)(v[length(v)])
pop = function(v)(v[-length(v)])
shift = function(v)(v[-1])
# differences between successive elements, first diff is first element with start
vectorLag = function(v, start = 0)pop(c(v, start) - c(start, v))
splitN = function(N, by = 4) vectorLag(round(cumsum(rep(N/by, by))));
splitToMax = function(N, max = 4) vectorLag(round(cumsum(rep(N/ceiling(N/max), ceiling(N/max)))));
# split into fixed block sizes + last incomplete block
splitBy = function(N, NperBlock = 4) {
	Nlast = N %% NperBlock;
	return(c(rep(NperBlock, N %/% NperBlock), if (Nlast == 0) c() else Nlast));
}

# cumsum returning indeces for numbers given in Ns
cumsumI = function(Ns, offset = 1, do.pop = FALSE) {
	cs = vectorNamed(c(0, cumsum(Ns)) + offset, c(names(Ns), 'N'));
	if (do.pop) cs = pop(cs);
	cs
}
# recursive cumsum (one level)
cumsumR = function(l, offset = 1) {
	cs0 = if (is.list(l)) lapply(l, cumsumR, offset = 0) else rev(cumsum(l))[1];
	cs = vectorNamed(c(0, pop(unlist(cs0))) + offset, names(cs0));
	cs
}

countsExtract = function(v, Ns, simplify = FALSE) {
	cnts = counts2idcs(Ns);
	r = apply(cnts, 1, function(r) {
		r = v[ r[1] : r[2] ];
		if (simplify) r else list(r)
	});
	return(if (!simplify) unlist.n(r, 1) else r);
}

#
#	<par> sets and permutations
#

#' @title wrapper for order to allow multivariate ordering
#'
#' @param v object (vector or data frame) for which order is to be calculated
#' @param ... additional arguemnts passed on to \code{order}
#' @return order of the object
#' @seealso {order{}} which this function wraps around
Order = function(v, ...) {
	if (is.data.frame(v)) do.call(order, lapply(v, identity), ...) else
	if (is.list(v)) do.call(order, v, ...) else
	order(v, ...)
}

#' @title Return all value combinations appearing in a data frame
#'
#' @param d data frame for which value combinations are to be caclulated
#' @return list with all value combinations present in \code{d}
# #' @examples
# #'
# #' combs = valueCombinations(iris);
# #'
valueCombinations = function(d) merge.multi.list(dimnames(table(d)));

#' @title Computes order so that inverseOrder after order is the identity
#'
#' Caculate ranks for arguemnt \code{p}. Works on vactors and data frames.
#'
#' @param p object for which ranks are to be comptued
#' @return vector of ranks of elements of \code{p}
#'
# #' @examples
# #' v = runif(1e2);
# #' print(all(sort(v)[inverseOrder(v)] == v))
Rank = inverseOrder = inversePermutation = function(p) {
	## <p> naive version
	# 	o = order(p);
	# 	i = rep(NA, length(o));
	# 	for (j in 1:length(o)) { i[o[j]] = j};
	# 	i
	## <p> build-in version (not working for multivariate case)
	#rank(v, ties.method = 'first')
	## <p> better version
	which.indeces(1:(if (class(p) == 'data.frame') nrow(p) else length(p)), Order(p))
}

#' @title Calculates inverseOrder, assuming that the argument is already an \code{order}-vector.
#'
#' @param p obect for which the inverse order is to be calculated
#' @return vector with integers representing the inverse order
inverseOrder_fromOrder = function(p)which.indeces(1:length(p), p)

#' @title Return vector that reorders v to equal reference.
#'
#' Assuming that two arguments are permutaions of each other, return a vector of indeces such that \code{all(reference == v[order_align(reference, v)]) == TRUE} for all vectors \code{reference, v}.
#'
#' @param reference vector with the reference ordering
#' @param v vector that is to be ordered the same way as \code{reference}
#' @return vector of indeces so that \code{v[return_value]} is the same as \code{reference}
#'
# #' @examples
# #' sapply(1:10, function(i){v = sample(1:5); v[order_align(5:1, v)]})
# #' sapply(1:10, function(i){
# #'    v = runif(1e2); v1 = sample(v, length(v));
# #'    all(v1[order_align(v, v1)] == v)
# #' })
order_align = function(reference, v)Order(v)[inverseOrder(reference)];

#' @title Calculates \code{order_align}, assuming that the both arguments are already orders.
#'
#' Analogous to \code{order_align} under the assumption that provided arguments are orders.
#'
#' @param reference order of a reference vector
#' @param v order of vector that is to be brought into the order of \code{reference}
#' @return order that can be applied to the orignal vector (from which \code{v} was calculated) to make it identical to the vector underlying \code{reference}
#'
# # ' @examples
# # ' \dontrun{
# # '   sapply(1:40, function(i){
# # '     v = runif(1e2);
# # '     v1 = sample(v, length(v));
# # '     all(v1[order_align_fromOrder(order(v), order(v1))] == v)
# # '   })
# # ' }
order_align_fromOrder = function(reference, v)v[inverseOrder_fromOrder(reference)];

# permutation is in terms of elements of l (not indeces)

applyPermutation = function(l, perm, from = 'from', to = 'to', returnIndeces = TRUE) {
	# 1. bring perm[[from]] in the same order as l
	# 2. apply this order to perm[[to]]
	r0 = perm[[to]][order(perm[[from]])[inverseOrder(l)]];
	# 3. determine permutation going from l to r0
	r = order(l)[inverseOrder(r0)]
	if (!returnIndeces) r = l[r];
	r
}

order.df = function(df, cols = NULL, decreasing = FALSE, na.last = FALSE) {
	if (is.null(cols)) cols = 1:ncol(df);
	if (!is.numeric(cols)) cols = which.indeces(cols, names(df));
	orderText = sprintf("order(%s, decreasing = %s, na.last = %s)",
		paste(sapply(cols, function(i) { sprintf("df[, %d]", i) }), collapse = ", "
		), as.character(decreasing), as.character(na.last)
#		paste(sapply(cols, function(i) {
#			if (is.numeric(i)) sprintf("df[, %d]", i) else sprintf("df$%s", i) }), collapse = ", "
#		), as.character(decreasing), as.character(na.last)
	);
	o = eval(parse(text = orderText));
	#print(list(text = orderText, order = o, df=df));
	o
}

order.df.maps = function(d, maps, ..., regex = FALSE) {
	cols = NULL;
	for (i in 1:length(maps)) {
		m = names(maps)[i];
		map = maps[[i]];
		keys = names(map);
		cols = c(cols, if (is.list(map)) {
			tempColName = sprintf("..order.df.maps.%04d", i);
			col = if (regex)
				sapply(d[[m]], function(e){ j = which.indeces(e, keys, regex = TRUE, inverse = TRUE)
					if (length(j) == 0) NA else map[[j]]
				}) else	as.character(map[d[[m]]]);
			col[col == "NULL"] = NA;
			d = data.frame(col, d, stringsAsFactors = FALSE);
			names(d)[1] = tempColName;
		} else { m });
	}
	o = order.df(d, cols, ...);
	o
}

data.frame.union = function(l) {
	dfu = NULL;
	for (n in names(l)) {
		df = l[[n]];
		factor = rep(n, dim(df)[1]);
		dfu = rbind(dfu, cbind(df, factor));
	}
	dfu
}

#
#	<p> factors
#

# levels: take levels in that order, unmentioned levels are appended
# setLevels: restrict to these levels, else set to NA
# setLevelsTo: set names of levels to argument, set excess levels to NA
# group: group levels, set names to concatenations
#	recodeLevels(as.factor(c('AA', 'AG', 'GG')), group = list(1:2, 3))
recodeLevels = function(f, map = NULL, others2na = TRUE, levels = NULL, setLevels = NULL,
	setLevelsTo = NULL, sortLevelsByMap = TRUE, group = NULL) {
	r = f;
	# <!> overwrites map
	# <!><i> does not implement grouping by level spec
	if (notE(group)) {
		lvls = levels(f);
		map = unlist.n(lapply(group, function(e) {
			newLevel = join(lvls[e], ' ');
			mapEl = recycle(lvls[e], newLevel);
			listKeyValue(mapEl[[1]], mapEl[[2]])
		}), 1);
	}
	if (!is.null(map)) {
		# map others to NA
		if (others2na) {
			nonmentioned = setdiff(if (is.factor(f)) levels(f) else unique(f), names(map));
			map = c(map, listKeyValue(nonmentioned, rep(NA, length(nonmentioned))));
		}
		v = vector.replace(as.character(f), map);
		# test for integer before and after
		# special case eliminated as of 14.9.2018
		#if (is.integer(f)) v = as.integer(v);
		#if (is.factor(f)) v = factor(v, levels = unique(as.character(map)));
		#v = factor(v, levels = unique(as.character(map)));
		r = if (sortLevelsByMap)
			factor(v, levels = union(unique(map), setdiff(unique(v), unique(map)))) else
			as.factor(v);
		# <!> r = v, r <- v do not work here, remain local
	}
	if (!is.null(levels) || !is.null(setLevels)) {
		# <p> preparation
		fact0 = as.factor(r);
		levls = levels(fact0);
		r = levls[fact0];

		# <p> new levels
		levlsN0 = firstDef(setLevels, levels, levls);
		levlsN = c(levlsN0, setdiff(levls, levlsN0));

		# <p> remove unwanted levels
		if (!is.null(setLevels)) r = ifelse(r %in% setLevels, r, NA);
		# <p> rename levels
		if (!is.null(setLevelsTo)) {
			#r = drop.levels(ifelse(as.integer(r) <= length(setLevels), r, NA));
			# 14.1.2020
			r = droplevels(ifelse(as.integer(r) <= length(setLevels), r, NA));
			levels(r) = setLevelsTo;
		}
		r = factor(r, levels = if (!is.null(setLevels)) levlsN0 else levlsN);
	}
	r
}
levelsSort = function(fac)recodeLevels(fac, listKeyValue(sort(levels(fac)), sort(levels(fac))))


factor2int = function(f)as.integer(as.character(f))
factor2numeric = function(f)as.numeric(as.character(f))

#
#	</p> factors
#

Union = function(..., .drop = TRUE, as.list = FALSE) {
	l = if (as.list) list(...)[[1]] else list(...);
	l = list(...);
	# auto-detect list of values
	if (.drop && length(l) == 1 && is.list(l[[1]])) l = l[[1]];
	r = NULL;
	for (e in l) { r = union(r, e); }
	r
}
Intersect = function(..., .drop = TRUE, as.list = FALSE) {
	l = if (as.list) list(...)[[1]] else list(...);
	# auto-detect list of values
	if (.drop && length(l) == 1 && is.list(l[[1]])) l = l[[1]];
	r = l[[1]];
	for (e in l[-1]) { r = intersect(r, e); }
	r
}

intersectSetsCount = function(sets) {
	i = iterateModels(list(s1 = names(sets), s2 = names(sets)), function(s1, s2) {
		length(intersect(sets[[s1]], sets[[s2]]))
	}, lapply__ = lapply);
	#r = reshape.wide(Df(i$models_symbolic, count = unlist(i$results)), 's1', 's2');
	rM = matrix(i$results, nrow = length(sets), byrow = TRUE);
	dimnames(rM) = list(names(sets), names(sets));
	rM
}
unionCum = function(..., .drop = TRUE) {
	l = list(...);
	# auto-detect list of values
	if (.drop && length(l) == 1 && is.list(l[[1]])) l = l[[1]];
	r = l[1];
	if (length(l) > 1)
		for (n in names(l)[-1]) { r = c(r, List(union(r[[length(r)]], l[[n]]), names_ = n)); }
	r
}

# row bind of data.frames/matrices with equal number of cols
lrbind = function(l, as.data.frame = FALSE, names = NULL) {
	d = dim(l[[1]])[2];
	v = unlist(sapply(l, function(m) unlist(t(m))));
	m = matrix(v, byrow = TRUE, ncol = d);
	dimnames(m) = list(NULL, names(l[[1]]));
	if (as.data.frame) {
		m = data.frame(m);
		if (!is.null(names)) names(m) = names;
	}
	m
}

#
#	logic arrays/function on list properties
#

# old versions:
#	if (na.rm) v = v[!is.na(v)];
#	sum(v)	# old version: length((1:length(v))[v])
# same as in Rlab
count = function(v, na.rm = TRUE)sum(v, na.rm = na.rm)
# old versions:
#	if (na.rm) v = v[!is.na(v)]; (sum(v)/length(v))
#	{ length(v[v]) / length(v) }
# v assumed to be logical
fraction = function(v, na.rm = TRUE)mean(v, na.rm = na.rm);
# treat v as set
set.card = function(v)count(unique(v))
# cardinality of a set
size = function(set)length(unique(set));

# null is false
#nif = function(b)(!(is.null(b) | is.na(b) | !b))
#nif = function(b)sapply(b, function(b)(!(is.null(b) || is.na(b) || !b)))
nif = function(b) {
	if (length(b) == 0) return(FALSE);
	if (class(b) %in% c('formula', 'function', 'list', 'data.frame')) return(TRUE);
	!(is.null(b) | is.na(b) | !b)
}
Nif = function(b, allnif = TRUE, nonLogicalIsTrue = TRUE) {
	if (is.null(b)) return(FALSE);
	if (class(b) %in% c('formula', 'function')) return(TRUE);
	bLog = sapply(b, as.logical);
	b = ifelse(is.na(b) | sapply(b, class) == 'logical', bLog, nonLogicalIsTrue);
	summ = (if (allnif) all else any);
	r = summ(sapply(b, nif));
	r
}
# null is true
#nit = function(b)(is.null(b) | is.na (b) | b)
#nit = function(b)sapply(b, function(b)(is.null(b) || is.na (b) || b))
nit = function(b) {
	if (length(b) == 0) return(TRUE);
	is.null(b) | is.na (b) | b
}
# null is zero
#niz = function(e)ifelse(is.null(e) | is.na(e), 0, e)
niz = function(e)ifelse(is.null(e) | is.na(e), 0, e)

# null is na (or other special value
#niz = function(e)ifelse(is.null(e) | is.na(e), 0, e)
nina = function(e, value = NA)sapply(e, function(e)ifelse(is.null(e), value, e))
Nina = function(e, value = NA)if (length(e) == 0) value else nina(e, value);

# not empty
notE = function(e)(length(e) > 0);

plus = function(x)ifelse(x > 0, x, 0)
minus = function(x)ifelse(x < 0, x, 0)

#
#	<p> complex structures
#

#
# Averaging a list of data frames per entry over list elements
#

# meanMatrices = function(d) {
# 	df = as.data.frame(d[[1]]);
# 	ns = names(df);
# 	# iterate columns
# 	dfMean = sapply(ns, function(n) {
# 		m = sapply(d, function(e)as.numeric(as.data.frame(e)[[n]]));
# 		mn = apply(as.matrix(m), 1, mean, na.rm = TRUE);
# 		mn
# 	});
# 	dfMean
# }
meanMatrices = function(d) {
	dm = dim(d[[1]]);
	good = sapply(d, function(m)(length(dim(m)) == 2 && all(dim(m) == dm)));
	if (any(!good)) warning('meanMatrices: malformed/incompatible matrices in list, ignored');
	d = d[good];
	m0 = sapply(d, function(e)avu(e));
	m1 = apply(m0, 1, mean, na.rm = TRUE);
	r = matrix(m1, ncol = dm[2], dimnames = dimnames(d[[1]]));
	r
}
meanVectors = function(d) {
	ns = names(d[[1]]);
	mn = apply(as.matrix(sapply(d, function(e)e)), 1, mean, na.rm = TRUE);
	mn
}
meanList = function(l)mean(as.numeric(l));

meanStructure = function(l) {
	r = nlapply(names(l[[1]]), function(n) {
		meanFct =
			if (is.matrix(l[[1]][[n]])) meanMatrices else
			if (length(l[[1]][[n]]) > 1) meanVectors else
				meanList;
		meanFct(list.key(l, n, unlist = FALSE));
	});
	r
}

matrixCenter = function(m, direction = 2, centerBy = median) {
	center = apply(m, direction, centerBy, na.rm = TRUE);
	m = if (direction == 1) (m - center) else t(t(m) - center);
	list(matrix = m, center = center)
}

matrixDeCenter = function(m, center, direction = 2) {
	m = if (direction == 1) t(t(m) + center) else (m + center);
	m
}


#
#	<p> combinatorial functions
#

# form all combinations of input arguments as after being constraint to lists
# .first.constant designates whether the first list changes slowest (TRUE) or fastest (FALSE)
#	in the resulting data frame,
#	i.e. all other factors are iterated for a fixed value of l[[1]] (TRUE) or not
# .constraint provides a function to filter the resulting data frame
merge.multi.list = function(l, .col.names = NULL, .col.names.prefix = "X",
	.return.lists = FALSE, .first.constant = TRUE, stringsAsFactors = FALSE, .cols.asAre = FALSE, .constraint = NULL, ...) {
	# <p> determine column names of final data frame
	.col.names.generic = paste(.col.names.prefix, 1:length(l), sep = "");
	if (is.null(.col.names)) .col.names = names(l);
	if (is.null(.col.names)) .col.names = .col.names.generic;
	.col.names[.col.names == ""] = .col.names.generic[.col.names == ""];
	names(l) = .col.names;		# overwrite names
	# <p> construct combinations
	if (.first.constant) l = rev(l);
	df0 = data.frame();
	if (length(l) >= 1) for (i in 1:length(l)) {
		newNames = if (.cols.asAre) names(l[[i]]) else names(l)[i];
		# <p> prepare data.frame: handle lists as well as data.frames
		# <!> changed 22.3.2016
		#dfi = if (is.list(l[[i]])) unlist(l[[i]]) else l[[i]];
		dfi = if (!is.data.frame(l[[i]])) unlist(l[[i]]) else l[[i]];
		df1 = data.frame.types(dfi, names = newNames, stringsAsFactors = stringsAsFactors);
		# <p> perform merge
		df0 = if (i > 1) merge(df0, df1, ...) else df1;
	}
	if (.first.constant) df0 = df0[, rev(names(df0)), drop = FALSE];
	if (.return.lists) df0 = apply(df0, 1, as.list);
	if (!is.null(.constraint)) {
		df0 = df0[apply(df0, 1, function(r).do.call(.constraint, as.list(r))), ];
	}
	df0
}

# list of list, vector contains index for each of these lists to select elements from
#	these elements are merged and return
#	if sub-element is not a list, take name of sub-element and contruct list therefrom
#	namesOfLists controls whether, if a selected element is a list, its name is used instead
#		can be used to produce printable summaries
list.takenFrom = function(listOfLists, v) {
	ns = names(listOfLists);
	if (any(ns != names(v))) v = v[order_align(ns, names(v))];
	l = lapply(1:length(v), function(i) {
		new = if (!is.list(listOfLists[[i]]))
			listKeyValue(ns[i], listOfLists[[i]][v[i]]) else {
				t = listOfLists[[i]][[v[i]]];
				# list of vectors
				t = (if (!is.list(t)) {
					# define name from higher level
					listKeyValue(firstDef(
						names(listOfLists[[i]])[v[i]], ns[i]
					), list(t))
					# <A> probably better and correct
					#listKeyValue(ns[i], list(t))
				} else if (is.null(names(t))) listKeyValue(ns[i], t) else t);
				t
			}
	});
	names(l) = names(v);
	l
}
# simplified version of list.takenFrom
list.extract = function(lol, idcs)pairsapplyLV(lol, idcs, function(l, i)l[i])
list.extractRows = function(lol, idcs)t(pairsapplyLV(lol, idcs, function(l, i)l[i, ]))


merge.lists.takenFrom = function(listOfLists, v) {
	merge.lists(list.takenFrom(listOfLists, v), listOfLists = TRUE);
}

merge.lists.takenFrom_old = function(listOfLists, v) {
	l = list();
	ns = names(listOfLists);
	if (any(ns != names(v))) v = v[order_align(ns, names(v))];
	for (i in 1:length(v)) {
		new = if (!is.list(listOfLists[[i]]))
			listKeyValue(ns[i], listOfLists[[i]][v[i]]) else {
				t = listOfLists[[i]][[v[i]]];
				# list of vectors
				t = (if (!is.list(t)) {
					# define name from higher level
					listKeyValue(firstDef(
						names(listOfLists[[i]])[v[i]], ns[i]
					), list(t))
					# <A> probably better and correct
					#listKeyValue(ns[i], list(t))
				} else if (is.null(names(t))) listKeyValue(ns[i], t) else t);
				t
			}
		l = merge.lists(l, new);
	}
	l
}

# take indeces given by v from a nested list
# namesOfLists: take the name of the list at the position in v
#	if null, take first element or leave aggregation to the function aggregator
# aggregator: called with the final result, should flatten existing lists into characters
lists.splice = function(listOfLists, v, namesOfLists = FALSE, aggregator = NULL, null2na = TRUE) {
	ns = names(listOfLists);
	l = lapply(1:length(ns), function(i) {
		name = ns[i];
		e = listOfLists[[i]][v[i]];
		r = if (!is.list(e)) e else {
			f = if (namesOfLists) {
				g = names(e)[1];
				# handle name == NULL
				if (is.null(g)) {
					# make an attempt later to print element
					#if (!is.null(aggregator)) e[[1]] else e[[1]][[1]]
					if (!is.null(aggregator))
						e[[1]] else
						join(as.character(e[[1]][[1]]), ", ")
				} else g
			} else e[[1]];
		}
		r
	});
	if (null2na) l = lapply(l, function(e)ifelse(is.null(e), NA, e));
	if (!is.null(aggregator)) l = aggregator(listKeyValue(ns, l), v, l);
	l
}

# dictionary produced by lists.splice, v: splice vector, l: aggregated list (w/o names)
merge.multi.symbolizer = function(d, v, l)unlist.n(d, 1);

merge.multi.list.symbolic = function(modelList, ..., symbolizer = NULL) {
	modelSize = lapply(modelList, function(m)1:length(m));
	models = merge.multi.list(modelSize, ...);
	namesDf = if (is.null(symbolizer)) names(modelList) else NULL;
	df0 = sapply(1:nrow(models), function(i, ...) {
		r = lists.splice(modelList, unlist(models[i, ]),
			namesOfLists = TRUE, aggregator = symbolizer);
		r
	});
	r = Df_(df0, t_ = TRUE, names = namesDf);
	r
}

inlist = function(l)lapply(l, function(e)list(e));
Inlist = function(...)inlist(list(...));

Do.callIm = function(im__f, args, ..., restrictArgs = TRUE, callMode = 'inline') {
	if (callMode == 'inlist') {
		.do.call(im__f, c(args, list(...)), restrictArgs = restrictArgs)
	} else if (callMode == 'list') {
		im__f(unlist.n(args, 1, reset = TRUE), ...)
	} else if (callMode == 'inline') {
		args = c(merge.lists(args, listOfLists = TRUE), list(...));
		.do.call(im__f, args, restrictArgs = restrictArgs)
	} else stop('Unknown call mode');
}

Kronecker = function(l, ...) {
	if (length(l) == 1) return(l[[1]]);
	kronecker(l[[1]], Kronecker(l[-1], ...), ...);
}


# <!> should be backwards compatible with iterateModels_old, not tested
# modelList: list of lists/vectors; encapuslate blocks of parameters in another level of lists
# Example:
#
#' Iterate combinations of parameters
#'
#' This function takes a list of parameters for which several values are to be evaluated. These values can be vectors of numbers or lists that contain blocks of parameters. All combinations are formed and passed to a user supplied function \code{f_iterate()}. This functions takes an index of the combination together with parameter values. Argument \code{callWithList} controls whether there is exactly one argument per parameter position or wether one more step of unlisting takes place. In case that a block of parameters is supplied, all values of the block are passed as individual arguments to \code{f_iterate()} in case \code{callWithList == FALSE}.
#'
#' #@param selectIdcs restrict models to the given indeces
#' @param modelList list specifying the models (see details)
#' @param models matrix containing indeces to sub-models (see details)
#' @param f_iterate function to be iterated across models
#' @param callWithList boolean to indicate whether model combination is to be supplied as a list.
#'   Otherwise model specification is inlined as arguments (see details)
#' @param callMode 'inline', 'list', 'inlist'
#' @param restrictArgs boolean to indicate whether over-supplied arguments (with respect to \code{f_iterate()})
#"   should be ignored. Otherwise, an error will be raised.
#' @param parallel boolean to inidcate whether iteration should be parallelized with
#'    \code{parallelize.dynamic}
#' @param lapply__ the iterator to be used (ignored at this moment)
#' @param ... extra arguments to be passed to \code{f_iterate()}
#' @return list containing the result of \code{f_iterate()} for all paramter combinations
#'
# #' @examples
# #' \dontrun{
# #' modelList = list(global = list(list(a=1, b=2)), N = c(1, 2, 3));
# #' print(iterateModels(modelList));
# #' modelList = list(N = c(1, 2, 3), parsAsBlock = list(list(list(c = 1, d = 2)),
# #'   list(list(c = 3, d = 4))));
# #' print(iterateModels(modelList));
# #' # ensure elements on A are given as a block (list)
# #' A = list(list(a = 1, b = 2), list(a = 3, b = 5));
# #' modelList = list(N = inlist(A), parsAsBlock = list(list(list(c = 1, d = 2)),
# #'   list(list(c = 3, d = 4))));
# #' print(iterateModels(modelList));
# #' # shorter version of the above
# #' modelList = list(N = Inlist(list(a = 1, b = 2), list(a = 3, b = 5)),
# #'   parsAsBlock = Inlist(list(c = 1, d = 2), list(c = 3, d = 4)));
# #' print(iterateModels(modelList));
# #' # inline calling
# #' modelList = list(N = list(list(a = 1, b = 2), list(a = 3, b = 5)),
# #'   parsAsBlock = list(list(c = 1, d = 2), list(c = 3, d = 4)));
# #' print(iterateModels(modelList));
# #' }
iterateModels_raw = function(modelList, models, f_iterate = function(...)list(...), ...,
	callWithList = FALSE, callMode = NULL, restrictArgs = TRUE, parallel = FALSE, lapply__) {
	if (!parallel) Lapply = lapply;
	if (is.null(callMode)) callMode = if (callWithList) 'list' else 'inline';
	# model indeces contains the original positions in models
	# this allows reordering of execution, eg with reverseEvaluationOrder
	r = Lapply(1:nrow(models), function(i, ..., im__f, im__model_idcs) {
		args = c(list(i = list(i = im__model_idcs[i])), list.takenFrom(modelList, unlist(models[i, ])));
#if (callMode == 'list') browser();
		Do.callIm(im__f, args, ..., restrictArgs = restrictArgs, callMode = callMode);
	}, ..., im__f = f_iterate, im__model_idcs = as.integer(row.names(models)));
	r
}

# <i> refactor iterateModels to use iterateModels_prepare
iterateModels_prepare = function(modelList, .constraint = NULL,
	callWithList = FALSE, callMode = NULL, restrictArgs = TRUE, selectIdcs = NULL, .first.constant = TRUE) {
	# <p> preparation
	if (is.null(callMode)) callMode = if (callWithList) 'list' else 'inline';

	modelSize = lapply(modelList, function(m)1:length(m));
	models = merge.multi.list(modelSize, .first.constant = .first.constant);

	# <p> handle constraints
	selC = if (is.null(.constraint)) TRUE else
		unlist(iterateModels_raw(modelList, models, f_iterate = .constraint,
			parallel = FALSE, callMode = callMode, restrictArgs = restrictArgs));
	selI = if (is.null(selectIdcs)) TRUE else 1:nrow(models) %in% selectIdcs;
	#	apply constraints
	models = models[selC & selI, , drop = FALSE];
	r = list(
		modelsRaw = models,
		selection = selC & selI,
		models = models
	);
	r
}

iterateModelsDefaultSymbolizer = function(i, ...) {
	l = list(...);
	r = lapply(l, function(e)unlist(as.character(unlist(e)[1])));
	r
}
iterateModelsJoinSymbolizer = function(i, ..., sep = ':') {
	l = list(...);
	r = lapply(l, function(e)join(unlist(as.character(unlist(e))), sep));
	r
}
iterateModelsSymbolizer = function(i, ..., im_symbolizer, im_symbolizerMode) {
	l = list(...);
	l0 = iterateModelsDefaultSymbolizer(i, ...);
	l1 = .do.call(im_symbolizer, c(list(i = i), list(...)), restrictArgs = TRUE);
	r = merge.lists(l0, l1);
	r
}

# <i> make name of supplied model index, currently 'i', configurable
iterateModels = function(modelList, f = function(...)list(...), ...,
	.constraint = NULL, .clRunLocal = TRUE, .resultsOnly = FALSE, .unlist = 0,
	callWithList = FALSE, callMode = NULL,
	symbolizer = iterateModelsDefaultSymbolizer, symbolizerMode = 'inlist',
	restrictArgs = TRUE, selectIdcs = NULL,
	.first.constant = TRUE, parallel = FALSE, lapply__, reverseEvaluationOrder = TRUE,
	modelTags = FALSE) {
	# <p> pre-conditions
	nsDupl = duplicated(names(modelList));
	if (any(nsDupl))
		stop(con('iterateModels: duplicated modelList entries: ', join(names(modelList)[nsDupl], ', ')));

	# <p> preparation
	if (is.null(callMode)) callMode = if (callWithList) 'list' else 'inline';

	# <p> produce raw combinations
	modelSize = lapply(modelList, function(m)1:length(m));
	models = merge.multi.list(modelSize, .first.constant = .first.constant);
# 	models_symbolic = merge.multi.list.symbolic(modelList,
# 		symbolizer = symbolizer, .first.constant = .first.constant);
	models_symbolic = do.call(rbind, iterateModels_raw(modelList, models, iterateModelsSymbolizer,
		callMode = 'inlist', parallel = FALSE,
		im_symbolizerMode = symbolizerMode, im_symbolizer = symbolizer));

	# <p> handle constraints
	selC = if (is.null(.constraint)) TRUE else
		unlist(iterateModels_raw(modelList, models, f_iterate = .constraint,
			callMode = callMode, restrictArgs = restrictArgs, ..., parallel = FALSE));
	selI = if (is.null(selectIdcs)) TRUE else 1:nrow(models) %in% selectIdcs;
	# <p> apply constraints
	models = models[selC & selI, , drop = FALSE];
	models_symbolic = models_symbolic[selC & selI, , drop = FALSE];

	# <p> models to be iterated
	modelsIt = if (reverseEvaluationOrder) models[rev(1:nrow(models)), , drop = FALSE] else models;
	metaStartTime = Sys.time();
	r = iterateModels_raw(modelList, modelsIt, f_iterate = f,
		callMode = callMode, restrictArgs = restrictArgs, ..., parallel = parallel);
	metaStopTime = Sys.time();
	if (reverseEvaluationOrder) r = rev(r);
	if (modelTags) {
		ns = dimnames(models_symbolic)[[2]]
		names(r) = apply(models_symbolic, 1, function(mr)
			join(apply(cbind(ns, mr), 1, join, sep = ':'), sep = '-'));
		#dimnames(models_symbolic)[[1]] = names(r);
	}
	r = if (.resultsOnly) r else list(
		models = models,
		results = r,
		models_symbolic = models_symbolic,
		meta = list(
			time = list(start = metaStartTime, stop = metaStopTime, duration = metaStopTime - metaStartTime)
		)
	);
	r = unlist.n(r, .unlist);
	r
}

iterateModelsExpand = function(modelList, .constraint = NULL) {
	modelSize = lapply(modelList, function(m)1:length(m));
	models = merge.multi.list(modelSize, .constraint = .constraint);
	r = list(
		models = models,
		models_symbolic = merge.multi.list.symbolic(modelList, .constraint = .constraint)
	);
	r
}

IterateModelsExpand = function(modelList, .constraint = NULL) {
	iterateModels(modelList, identity, .constraint = .constraint, callWithList = TRUE)$results
}

# reverse effect of .retern.lists = TRUE
#	list.to.df(merge.multi.list(..., .return.lists = TRUE)) === merge.multi.list(..., .return.lists = FALSE)
list.to.df = function(l)t(sapply(l, function(e)e))

merge.multi = function(..., .col.names = NULL, .col.names.prefix = "X",
	.return.lists = FALSE, stringsAsFactors = FALSE, .constraint = NULL, .first.constant = TRUE) {
	merge.multi.list(list(...), .col.names = .col.names, .return.lists = .return.lists,
		stringsAsFactors = stringsAsFactors, .constraint = .constraint, .first.constant = .first.constant)
}

merge.multi.dfs = function(l, .first.constant = TRUE, all = TRUE, stringsAsFactors = FALSE, ...) {
	if (.first.constant) l = rev(l);
	if (length(l) >= 1) for (i in 1:length(l)) {
		df1 = data.frame.types(l[[i]], stringsAsFactors = stringsAsFactors);
		df0 = if (i > 1) merge(df0, df1, all = all, ...) else df1;
	}
	if (.first.constant) df0 = df0[, rev(names(df0)), drop = FALSE];
	df0
}

Merge = function(x, y, by = intersect(names(x), names(y)), ..., safemerge = TRUE, stableByX = FALSE) {
	if (stableByX) x = data.frame(x, MergeStableByX = 1:nrow(x));
	if (safemerge && length(by) == 0) {
		stop(sprintf('Merge: safemerge triggered. No common columns between "%s" and "%s"',
			join(names(x), sep = ','), join(names(y), sep = ',')))
	}
	r = merge(x = x, y = y, by = by, ...);
	if (stableByX) {
		indexCol = which(names(r) == 'MergeStableByX');
		r = r[order(r$MergeStableByX), -indexCol, drop = FALSE];
	}
	r
}

MergeByRowNames = function(x, y, ...) {
	dMerge = Merge(Df(x, ROW_NAMES__ = row.names(x)), Df(y, ROW_NAMES__ = row.names(y)),
		by = 'ROW_NAMES__', ...)
	Df_(dMerge, min_ = 'ROW_NAMES__', row.names = dMerge$ROW_NAMES__)
}

# ids: variables identifying rows in final table
# vars: each combination of vars gets transformed to an own column
# <!> not tested for length(ids) > 1 || ength(rvars) > 1
# blockVars: should the repeated vars go in blocks or be meshed for vars
#
# Examples:
# intersection table
# i = intersectSetsCount(sets);
# reshape.wide(Df(i$models_symbolic, count = unlist(i$results)), 's1', 's2');
reshape.wide = function(d, ids, vars, blockVars = FALSE, reverseNames = FALSE, sort.by.ids = TRUE) {
	# remaining vars
	rvars = setdiff(names(d), union(ids, vars));
	# levels of variables used in the long expansion
	levls = lapply(vars, function(v)unique(as.character(d[[v]])));
	# combinations at the varying vars as passed to vars
	cbs = merge.multi.list(levls, .col.names = vars, .first.constant = !blockVars);
	# repvars: repeated variables
	repvars = merge.multi.list(c(list(rvars), levls),
		.first.constant = !blockVars, .col.names = c("..var", vars));
	varnames = apply(repvars, 1, function(r)join(if (reverseNames) rev(r) else r, "."));

	r0 = data.frame.types(unique(d[, ids], drop = FALSE), names = ids);
	r1 = data.frame.types(apply(r0, 1, function(r) {
		# <p> isolate rows which match to current id columns
		ids = which(apply(d[, ids, drop = FALSE], 1, function(id)all(id == r)));
		d1 = d[ids, ];
		# <p> construct vector of repeated values
		vs = sapply(1:dim(cbs)[1], function(i) {
			# <A> should be equal to one
			row = which(apply(d1[, vars, drop = FALSE], 1, function(r)all(r == cbs[i, ])));
			v = if (length(row) != 1) rep(NA, length(rvars)) else d1[row, rvars];
			v
		});
		# heed blockVars
		vs = as.vector(unlist(if (!blockVars) t(vs) else vs));
		vs
	}), do.transpose = TRUE, names = varnames);
	r = data.frame(r0, r1);
	if (sort.by.ids) r = r[order.df(r, ids), ];
	row.names(r) = NULL;
	r
}

#' Convert data in wide format to long format
#' 
#' Long format duplicates certain columns and adds rows for which one new column hold values coming
#' from a set of columns in wide format. Does not allow for parallel reshaping.
#'
#' @param d data frame with columns in wide format
#' @param vars columns in wide format by name or index
#' @param factors \code{vars} can be grouped. For each level of \code{factor} a new row is created. Implies
#'			that \code{length(vars)} is a multiple of \code{length(levels(factor))}
#' @param factorColumn name of the column to be created for the factor
#' @param valueColumn name of the new column of values that were in wide format
# factors: provide factor combinations explicitly for vars (otherwise split by '.', <i>)
#' @param rowNamesAs name of the column that should contain row names
#' @return data frame in long format
# #' @examples
# #' \dontrun{
# #'	#reshape variables 2:9 (forming two groups: case/ctr), value of which is named 'group'
# #'	# the shortened columns will get names valueColumn
# #'	d0 = reshape.long(d, vars = 2:9, factors = c('case', 'ctr'), factorColumn = 'group',
# #'		valueColumn = c('AA', 'AG', 'GG', 'tot'));
# #'
# #'	# reshape several grouped columns
# #' 	d2 = reshape.long(d1, vars = avu(vs),
# #'		factorColumn = 'time', valueColumn = valueNames, factors = as.factor(1:3));
# #'	}
reshape.long = function(d, vars = NULL, factorColumn = 'factor', valueColumn = 'value',
	factors = as.factor(vars), rowNamesAs = NULL) {
	if (is.null(vars)) vars = names(d);
	# make rownames an extra column
	if (!is.null(rowNamesAs)) {
		d = data.frame(reshape_row_names__ = rownames(d), d);
		names(d)[1] = rowNamesAs;
	}
	# indeces of columns vars
	Ivars = .df.cols(d, vars);
	# remaining vars
	rvars = setdiff(1:length(names(d)), Ivars);
	# names thereof
	Nrvars = names(d)[rvars];

	# how wide are the blocks?
	S = length(vars) / length(factors);
	# columns of intermediate data.frame
	N = length(rvars);
	# create list of data frames
	dfs = lapply(1:nrow(d), function(i) {
		st = d[i, rvars];	# start of the new row
		df0 = data.frame(factors, value =  matrix(d[i, vars], nrow = length(factors), byrow = TRUE));
		df1 = data.frame(st, df0, row.names = NULL);
		names(df1) = c(Nrvars, factorColumn, valueColumn);
		df1
	});
	#r = rbindDataFrames(dfs, do.unlist = TRUE, useDisk = useDisk);
	r = do.call(rbind, dfs);
	r
}

DfUniqueRowsByCols = function(d, cols) {
	row.names(d) = NULL;
	as.integer(row.names(unique(d[, cols, drop = FALSE])))
}

#' Reduce data frame to be unique on subset of columns
#'
#' Reduce data frame by picking the first row of blocks for which \code{cols} has the same values.
#'
#' @param d data frame to be made unique
#' @param cols columns for which the reduced data frame has to be unique
#' @param drop argument passed to subset selection \code{`[`}
#' @return the reduced data frame
DfUniqueByCols = uniqueByCols = function(d, cols, drop = FALSE) {
	d[DfUniqueRowsByCols(d, cols), , drop = drop]
}

# robustly access columns: if column name is NA, add column of NAs
# changed as of 28.11.2018 <!> rely on only use by Reshape.long.raw
DfSelectCols = function(d, vars) {
	# changed as of 28.11.2018
	#d0 = do.call(cbind, lapply(vars, function(v)if (is.na(v)) NA else d[, v]));
	d0 = do.call(cbind, lapply(vars, function(v)if (is.na(v)) NA else d[, v, drop = FALSE]));
	d0
}

#	vars:	columns for which are to be in long format
#	lvMap:	mapping from levels to columns in d (wide columns)
Reshape.long.raw = function(d, vars, lvMap, factorColumn = 'repeat',
	valuePostfix = '_long', varsLong = paste(vars, valuePostfix, sep = '')) {
	if (any(sapply(lvMap, length) != length(vars)))
		stop('selected variables per level of different length to result columns');
	# remaining vars
	rvars = setdiff(names(d), na.omit(Avu(lvMap)));
	# levels
	lvls = names(lvMap);
	# create list of data frames
	dfs = lapply(1:nrow(d), function(i) {
		dR = d[i, rvars, drop = FALSE];	# fixed, repeated part of the data set
		d0L = lapply(lvls, function(l)DfSelectCols(d[i, , drop = FALSE], lvMap[[l]]));
		d0 = do.call(rbind, lapply(d0L, setNames, varsLong));
		d1 = Df(index = lvls,
			Df_(d0, row.names = NULL), Df_(dR, row.names = NULL), names = c(factorColumn, varsLong));
		#d1 = Df(index = lvls, d0, dR, names = c(factorColumn, varsLong));
		d1
	});
	r = do.call(rbind, dfs);
	r
}

Reshape.levelMap_re = function(ns, vars, factorsRe) {
	# regular expressions for columns to be reshaped
	Res = sapply(vars, function(v)Sprintf(factorsRe, COLIDENT = v));
	# perform RE search
	lvlsRaw = Regex(Res, ns);
	lvlsRawL = sapply(lvlsRaw, length);

	# levels of index/reshape column
	#cols = Df_(lvlsRaw, names = vars);
	# level belonging to column (non-simplifying Regex)
	lvCol_old = RegexL(Res, ns, captures = TRUE);
	# allow to concat matches (several captures per Re)
	lvCap = lapply(lapply(Regexpr(Res, ns, captures = TRUE, reSimplify = FALSE), setNames, ns), unlist);
	lvCol = lapply(lvCap, filterList, f = function(e)e != '');
	# prepare level -> column mapping
	names(lvCol) = vars;
	lvls = unique(Avu(lvCol));
	if (any(sapply(lvCol, length) < sapply(lvlsRaw, length))) {
		print(list(lvlsRaw = lvlsRaw, lvls = lvCol));
		stop('Could not extract values for levels for all variables');
	}
	# 	# 28.11.2018: allow levels to be embedded
	# 	if (!all(lvlsRawL == lvlsRawL[1])) {
	# 		print(lvlsRaw);
	# 		stop('Different number of levels per factor');
	# 	}
	# map from level to columns
	lvMap = nlapply(lvls, function(l)Avu(nina(lapply(lvCol, function(c)names(c)[which(c == l)]))));
	return(lvMap);
}

Reshape.levelMap_list = function(ns, vars, factorsRe) {
	lels = sapply(vars, is.list);	# list elements
	vL = vars[lels];

	# <p> check input
	unmatched = unlist(vL)[which(!(unlist(vL) %in% ns))];
	if (length(unmatched) > 0)
		stop(Sprintf('reshape variables [%{r}s] do not exist in data', r = join(unmatched, ', ')));

	# <p> process re-matched variables
	lvMapRe = if (any(!lels)) Reshape.levelMap_re(ns, unlist(vars[!lels]), factorsRe) else list();
	# <p> take levels of repetition from Re variables if available, else enumerate
	levels = if (length(lvMapRe) > 0) names(lvMapRe) else as.character(1:length(vL[[1]]));

	# <p> construct level-map for explicit variable names
	lvMapL = lapply(Df(sapply(vL, identity), t_ = TRUE), unlist);
	names(lvMapL) = levels;
	lvMap = merge.lists(lvMapRe, lvMapL, concat = TRUE)
	lvMap
}

Reshape.levelMap = function(ns, vars, factorsRe) {
	if (is.character(vars)) return(Reshape.levelMap_re(ns, vars, factorsRe));
	if (!is.list(vars)) stop('invalid variable specification');
	lvMap = Reshape.levelMap_list(ns, vars, factorsRe);
	lvMap
}

# allow parallel re-shaping, i.e. take columns of form 'prefix.\d' and take \d as the value for the new
#	index column (reshape-column)
#	vars: prefix of columns to be reshaped
#	factorsRe: re to append to vars to identify wide columns
Reshape.long = function(d, vars, factorColumn = 'repeat', valuePostfix = '_long',
	factors = NULL, factorsRe = '^%{COLIDENT}s[._]?(\\d+)', useDisk = FALSE, rowNamesAs = NULL,
	varsLong = paste(vars, valuePostfix, sep = '')) {

	lvMap = Reshape.levelMap(names(d), vars, factorsRe);
	if (is.list(vars)) {
		nsVars = names(vars);
		varsLong[nsVars != ''] = nsVars[nsVars != ''];
	}
	Reshape.long.raw(d, vars, lvMap,
		factorColumn = factorColumn, valuePostfix = valuePostfix, varsLong = varsLong);
}

# reshape rows in blocks to avoid memory exhaustion
Reshape.long.byParts = function(d, ..., N = 1e4, path = tempfile(), filter = NULL) {
	Nrow = nrow(d);
	Nparts = ceiling(Nrow / N);

	#Nparts = 2;
	for (i in 1:Nparts) {
		dP = d[ (N*(i - 1) + 1):min((N*i), Nrow), ];
		dL = Reshape.long(dP, ...);
		if (notE(filter)) dL = filter(dL);
		write.table(dL, file = path, col.names = i == 1, append = i != 1, row.names = F);
	}
	gc();
	return(readTable(Sprintf('[SEP=S,HEADER=T]:%{path}s')));
}

byParts = function(d, fn, ..., N = 1e4, path = tempfile(), filter = NULL) {
	Nrow = nrow(d);
	Nparts = ceiling(Nrow / N);

	#Nparts = 2;
	for (i in 1:Nparts) {
		dP = d[ (N*(i - 1) + 1):min((N*i), Nrow), ];
		dL = fn(dP, ...);
		if (notE(filter)) dL = filter(dL);
		write.table(apply(dL, 2, as.character), file = path, col.names = i == 1, append = i != 1, row.names = F);
	}
	gc();
	return(readTable(Sprintf('[SEP=S,HEADER=T]:%{path}s')));
}

reshape.long.byParts = function(d, ..., N = 1e4, path = tempfile(), filter = NULL) {
	byParts(d, reshape.long, ..., N = N, path = path, filter = filter);
}

#
# <p> string functions
#

uc.first = firstUpper = function(s) {
	paste(toupper(substring(s, 1, 1)), substring(s, 2), sep = "", collapse = "");
}
substrM1 = function(s)substr(s, 1, nchar(s) - 1)
strAbbr = function(s, N = 15, ellipsis = '...') {
	r = if (nchar(s) <= N) s else join(c(substr(s, 1, N - nchar(ellipsis)), ellipsis));
	r
}

deduplicateLabels = function(v, labels = v[duplicated(v)], sep = '-', firstUntouched = TRUE) {
	for (label in labels) {
		idcs = which(v == label);
		if (firstUntouched) idcs = shift(idcs);
		v[idcs] = paste(label, 1:length(idcs), sep = sep);
	}
	v
}

Trimws = function(s)join(sub('^\t', '', splitString("\n", s)), '\n');

# non-empty string
nonEmpty = function(e)nif(e != '')
# non-empty, unique
vecCharNEU = function(v)unique(Filter(nonEmpty, v))

#
#	<p> factor transformations for data frames
#

dataExpandedNames = function(data) {
	dnames = unlist(lapply(names(data), function(v){
		if (is.factor(data[[v]])) paste(v, 1:(length(levels(data[[v]])) - 1), sep = "") else v;
	}));
	dnames
}
# model.matrix removes missing columns and could not be tweaked into working
dataExpandFactors = function(data, vars =  NULL) {
	if (is.null(vars)) vars = names(data);
	d0 = lapply(vars, function(v) {
		if (is.factor(data[[v]])) {
			ls = levels(data[[v]]);
			dcNa = rep(NA, length(ls) - 1);	# missing data coding
			dc = rep(0, length(ls) - 1);	# dummy coding
			sapply(data[[v]], function(e) {
				if (is.na(e)) return(dcNa);
				i = which(e == ls);
				if (i == 1) return(dc);
				dc[i - 1] = 1;
				return(dc);
			});
		} else data[[v]];
	});
	d0names = dataExpandedNames(data[, vars]);
	# re-transform data
	d1 = data.frame(matrix(unlist(lapply(d0, function(e)t(e))), ncol = length(d0names), byrow = FALSE));
	names(d1) = d0names;
	d1
}
coefficientNamesForData = function(vars, data) {
	lnames = dataExpandedNames(data);	# names of levels of factors
	cnames = lnames[unlist(sapply(vars, function(v)which.indeces(v, lnames, regex = TRUE)))];
	cnames
}

#
# <p> statistic oriented data frame manipulation
#

variableIndecesForData = function(d, vars, varsArePrefixes = TRUE, varRegex = '%s.*') {
	if (varsArePrefixes) vars = sapply(vars, function(e)sprintf(varRegex, e));
	which.indeces(vars, names(d), regex = TRUE, match.multi = TRUE)
}
variablesForData = function(d, vars, varsArePrefixes = TRUE, varRegex = '%s.*') {
	names(d)[variableIndecesForData(d, vars, varsArePrefixes, varRegex)]
}

subData = function(d, vars, varsArePrefixes = TRUE) {
	dfr = d[, variableIndecesForData(d, vars, varsArePrefixes), drop = FALSE];
	dfr
}

subDataFromFormula = function(d, formula, responseIsPrefix = TRUE, covariateIsPrefix = TRUE) {
	resp = formula.response(formula);
	cov = formula.covariates(formula);
	ns = names(d);
	r = list(
		response = subData(d, resp, responseIsPrefix),
		covariate = subData(d, cov, covariateIsPrefix)
	);
	r
}

#
#	<p> graph functions
#

sub.graph.merge = function(df, leader, follower) {
	# next transitive step
	r0 = merge(df, data.frame(leader = leader, follower = follower), by = 'follower');
	# add new connections
	r1 = rbind(df, data.frame(follower = r0$leader.y, leader = r0$leader.x, cluster = r0$cluster));
	# symmetric closure
	r1 = rbind(r1, data.frame(follower = r1$leader, leader = r1$follower, cluster = r1$cluster))
	# form clusters by selecting min cluster number per connection
	r1 = r1[order(r1$cluster), ];
	row.names(r1) = 1:dim(r1)[1];
	r2 = unique(r1[, c('leader', 'follower')]);
	# select unique rows (first occurunce selects cluster)
	r = r1[as.integer(row.names(r2)), ];
	# pretty sort data frame
	r = r[order(r$cluster), ];
	r
}
# form clusters from a relationally defined hierarchy
sub.graph = function(df) {
	df = as.data.frame(df);
	names(df)[1:2] = c('follower', 'leader');
	df = df[order(df$follower), ];
	# seed clusters
	ids = sort(unique(df$follower));
	idsC = as.character(ids);
	counts = lapply(ids, function(id)sum(df$follower == id));
	names(counts) = idsC;
	clusters = unlist(sapply(idsC, function(id){ rep(as.integer(id), counts[[id]]) }));

	df = cbind(df, data.frame(cluster = rep(clusters, 2)));
	df = unique(rbind(df, data.frame(follower = df$leader, leader = df$follower, cluster = df$cluster)));
	# receiving frame
	df0 = df;
	# results with clusters
	i = 1;
	repeat {
		Nrows = dim(df0)[1];
		cls = df0$clusters;
		# add transitive connections
		df0 = sub.graph.merge(df0, follower = df0$leader, leader = df0$follower);
		if (dim(df0)[1] == Nrows && all(cls == df0$clusters)) break();
	}
	df0 = df0[order(df0$cluster), ];
	cIds = unique(df0$cluster);
	cls = lapply(cIds, function(id)unique(avu(df0[df0$cluster == id, c('follower', 'leader')])));
	cls
}

#
#	<p> formulas
#

# formula: formula as a character string with wildcard character '%'
# 	<!>: assume whitespace separation in formula between terms
#	<!>: write interaction with spaces <!> such as in:
#		f = 'MTOTLOS_binair ~ ZRES% + sq(ZRes%) + ( ZRES% )^2';
formula.re = function(formula, data, ignore.case = FALSE, re.string = '.*') {
	vars = names(data);
	#regex = '(?:([A-Za-z_.]+[A-Za-z0-9_.]*)[(])?([A-Za-z.]+[%][A-Za-z0-9.%_]*)(?:[)])?';
	#			function names				(    regex						   )
	#regex = '(?:([A-Za-z_.]+[A-Za-z0-9_.]*)[(])?([A-Za-z%.]+[A-Za-z0-9.%_]*)(?:[)])?';
	# allow backslash quoting
	regex = '(?:([A-Za-z_.\\\\]+[A-Za-z0-9_.\\\\]*)[(])?([A-Za-z%.\\\\]+[A-Za-z0-9.%_\\\\]*)(?:[)])?';
	patterns = unique(fetchRegexpr(regex, formula, ignore.case = ignore.case));
	subst = nlapply(patterns, function(p) {
		comps = fetchRegexpr(regex, p, captureN = c('fct', 'var'), ignore.case = ignore.case)[[1]];
		p = sprintf("^%s$", gsub('%', re.string, comps$var));
		mvars = vars[sapply(vars, function(v)regexpr(p, v, perl = TRUE, ignore.case = ignore.case)>=0)];
		if (comps$fct != '') {
			varf = sprintf('%s', paste(sapply(mvars, function(v)sprintf('%s(%s)', comps$fct, v)),
				collapse = " + "));
		} else {
			varf = sprintf('%s', paste(mvars, collapse = " + "));
		}
		varf
	});
	formula1 = mergeDictToString(subst, formula);
	formulaExp = as.formula(formula1);
	formulaExp
}

formula.response = function(f) {
	#r = fetchRegexpr('[^\\s~][^~]*?(?=\\s*~)', if (is.formula(f)) deparse(f) else f);
	f = if (class(f) == 'formula') Deparse(f) else f;
	r = as.character(fetchRegexpr('^\\s*([^~]*?)(?:\\s*~)', f, captures = TRUE));
	# <p> version 2
	#fs = as.character(as.formula(as.character(f)));	# "~" "response" "covs"
	#r = fs[2];
	# <p> version 1
	#f = as.formula(f);
	#r = all.vars(f)[attr(terms(f), "response")];	# fails to work on 'response ~ .'
	r
}
formula.rhs = function(f, noTilde = FALSE, as_character = FALSE) {
	rhs = fetchRegexpr('[~](.*)', if (!is.character(f)) formula.to.character(f) else f, captures = TRUE);
	r = if (noTilde) rhs else con('~', rhs);
	r = if (as_character) r else as.formula(r);
	r
}
formula.covariates = function(f) {
	covs = all.vars(formula.rhs(f));
	#covs = setdiff(all.vars(as.formula(f)), formula.response(f));
	covs
}
formula.vars = function(f)union(formula.response(f), formula.covariates(f));
#formula.vars = function(f)all.vars(as.formula(f));

formula.nullModel = function(f) {
	r = formula.response(f);
	fn = as.formula(sprintf("%s ~ 1", r));
	fn
}
formula.to.character = function(f)join(deparse(as.formula(f)), '');
Formula.to.character = function(f)ifelse(is.character(f), f, formula.to.character(f));

formula.expand = function(f, data) {
	if (is.null(f)) return(NULL);
	if (any(all.vars(formula.rhs(f)) == '.')) {
		covs = setdiff(names(data), as.character(formula.response(f)));
		f = formula.set.rhs(f, vars.as.rhs(covs));
	}
	return(f);
}

formula2filename = function(f) {
	fs = join(f, sep = '');
	filename = mergeDictToString(list(
		`\\s+` = '',
		`_` = '-',
		`%` = ':',
		`Surv\\(.*\\)` = 'surv',
		MARKER = 'snp'
		# other components
	), fs, re = TRUE, doApplyValueMap = FALSE, doOrderKeys = FALSE);
	filename
}
data.vars = function(data, formula, re.string = '.*', ignore.case = FALSE) {
	all.vars(formula.re(formula = formula, data = data, re.string = re.string, ignore.case = ignore.case));
}
data.vars.after = function(data, col, skip = TRUE) {
	ns = names(data);
	ns[(which(ns == col) + skip):length(ns)]
}

dataColRange = function(data, from = NULL, to = NULL) {
	ns = names(data);
	start = if (is.integer(from)) from else (if (notE(from)) which(ns == from) else 1);
	stop = if (is.integer(to)) to else (if (notE(to)) which(ns == to) else ncol(data));
	data[, start:stop, drop = FALSE]
}


# select column names based on res, negation or literal names
# dataSelectVars(data, ~ cg + ab, ~ 0)
# dataSelectVars(data, list(~ cg, ~ !ab))
dataSelectVars = function(data, prefix = list(), fixed = list()) {
	# <p> input sanitation
	ns = names(data);
	if (class(prefix) == 'formula') prefix = list(prefix);
	if (class(fixed) == 'formula') fixed = list(fixed);

	# <p> fixed
	vsF = do.call(Union, lapply(fixed, function(e)setdiff(all.vars(e), '0')));
	# <p> prefix
	vsP = do.call(Union, lapply(prefix, function(e) {
		negation = class(as.list(e)[[2]]) == 'call' && as.list(e)[[2]][[1]] == '!';
		formula = join(c('~', join(sapply(all.vars(e), function(v)sprintf('%s%%', v)), ' + ')));
		vars = all.vars(formula.re(formula, data));
		if (negation) setdiff(ns, vars) else vars
	}));

	r = Union(vsF, vsP);
	r
}
dataSelectCols = function(data, prefix, fixed = ~ 0) {
	data[, dataSelectVars(data, prefix, fixed), drop = FALSE]
}


formula.add.rhs = function(f0, f1, envir = parent.frame()) {
	if (is.character(f1)) f1 = as.formula(paste0(c('~', paste(f1, collapse = ' + ')), collapse = ''));
	as.formula(join(c(
		formula.to.character(f0),
		formula.rhs(f1, noTilde = TRUE, as_character = TRUE)), '+'), env = envir)
}
vars.as.rhs = function(v)as.formula(Sprintf('~ %{vars}s', vars = join(v, '+')))
formula.set.rhs = function(f0, f1, envir = parent.frame()) {
	as.formula(join(c(formula.response(f0), formula.rhs(f1, as_character = TRUE))), env = envir)
}
formula.add.responseByName = function(f0, response, envir = parent.frame()) {
	formula = join(c(response, formula.rhs(f0, noTilde = FALSE)), ' ');
	as.formula(formula, env = envir)
}
formula.add.response = function(f0, f1, envir = parent.frame()) {
	formula.add.responseByName(f0, formula.response(f1), envir = envir)
}
# <!><t> w/ transformations in formula
formula.predictors = function(f, data, dataFrameNames = TRUE) {
	if (formula.rhs(f) == ~ 1) return('(Intercept)');
	#mm = model.matrix(model.frame(formula.rhs(f), data), data);
	mm = model.matrix(formula.rhs(f), data);
	ns = dimnames(mm)[[2]];

	# <p> create data frame to extract proper names
# 	if (dataFrameNames) {
# 		df0 = as.data.frame(t(rep(1, length(ns))));
# 		names(df0) = ns;
#		ns = names(df0);
# 	}
	ns
}

# <!> cave survival
formulaRemoveTransformation = function(model) {
	respVar = setdiff(all.vars(model), all.vars(formula.rhs(model)));
	formula.add.response(formula.rhs(model), as.formula(Sprintf('%{respVar}s ~ 1')))
}

formulas.free = function(f1, f0, data) {
	setdiff(formula.predictors(f1, data), formula.predictors(f0, data))
}


# <i> use terms.formula from a (a + ... + z)^2 formula
# <i> merge.multi.list(rep.list(covs, 2), .constraint = is.ascending)
covariatePairs = function(covs) {
	pairs = merge(data.frame(c1 = 1:length(covs)), data.frame(c2 = 1:length(covs)));
	pairs = pairs[pairs[, 1] > pairs[ ,2], ];
	df = data.frame(c1 = covs[pairs[, 1]], c2 = covs[pairs[, 2]]);
	df
}

formulaWith = function(response = "y", covariates = "x") {
	if (!Nif(response)) response = '';
	as.formula(sprintf("%s ~ %s", response,  paste(covariates, collapse = "+")))
}

#
#	<p> set operations
#

minimax = function(v, min = -Inf, max = Inf) {
	r = ifelse(v < min, min, ifelse(v > max, max, v));
	r
}

#
#	<p> recycling
#

accessIdx = function(e, i, byRow = TRUE) {
	if (class(e) != 'matrix' || is.na(byRow)) e[i] else
		(if (byRow) e[i, , drop = FALSE] else e[, i, drop = FALSE])
}

# fixed as of 10.8.2018: different types not correctly handles <f>
Recycle = function(l, byRow = TRUE) {
	# determine recyling pattern
	# old version would not preserve type
	# Recycle = function(l)lapply(apply(do.call(cbind, l), 2, as.list), unlist)
	lTmp = lapply(l, function(e)1:length(e));
	rTmp = 	lapply(apply(do.call(cbind, lTmp), 2, as.list), unlist)
	# extract values per component
	r = lapply(seq_along(l), function(i)accessIdx(l[[i]], rTmp[[i]], byRow = byRow));
	return(setNames(r, names(l)));
}
recycle = function(...)Recycle(list(...));
recycleTo = function(..., to, simplify = TRUE) {
	r = recycle(to, ...)[-1];
	if (simplify && length(r) == 1) r[[1]] else r
}
