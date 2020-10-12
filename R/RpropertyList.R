#
#	RpropertyList.R
#Fri Jan  7 17:40:12 2011

# wrap string for property list
ws = function(s) {
	s = if (length(grep('^([_/\\a-zA-Z0-9.]+)$', s)) > 0) { s } else {
		s = gsub('([\\"])', '\\\\\\1', s);
		sprintf('"%s"', s);
	}
	s
}

# can a string be condensed into a single line
condense = function(s, ident, o) {
	if (nchar(s) + ident * o$tabWidth - nchar(grep("\t", s)) < o$screenWidth) {
		s = gsub("\n", ' ', s);
		s = gsub("\t", '', s);
	}
	s
}

stringFromPropertyI = function(obj, ident, o) {
	str = '';
	inS = join(rep("\t", ident), '');
	in1S = join(rep("\t", ident + 1), '');

	if ( class(obj) == 'function' ) {
		str = sprintf('%s%s', str, ws(join(deparse(obj), "\n")))
	} else if ( class(obj) != 'list' & length(obj) == 1 & !(o$kp %in% o$forceVectors)) {
		# <i> data support
		str = sprintf('%s%s', str, ws(obj));
	} else if (class(obj) == 'list' && !is.null(names(obj))) {
		hash = sprintf("{\n%s%s;\n%s}", in1S, paste(sapply(names(obj), function(k) {
			o = merge.lists(o, list(kp = sprintf('%s.%s', o$kp, k)));
			r = sprintf('%s = %s', ws(k), stringFromPropertyI(obj[[k]], ident+1, o))
			r
		}), collapse = sprintf(";\n%s", in1S)), inS);
		if (!o$noFormatting) hash = condense(hash, ident, o);
		str = sprintf('%s%s', str, hash);
	} else { # vector or anonymous list
		obj = as.list(obj);
		array = sprintf("(\n%s%s\n%s)", in1S, if (length(obj) < 1) '' else paste(
			sapply(1:length(obj), function(i) {
			e = obj[[i]];
			o = merge.lists(o, list(kp = sprintf('%s.[%d]', o$kp, i)));
			stringFromPropertyI(e, ident+1, o)
		}), collapse = sprintf(",\n%s", in1S)), inS);
		if (!o$noFormatting) array = condense(array, ident, o);
		str = sprintf('%s%s', str, array);
	}
	str
}

stringFromPropertyDefaults = list(screenWidth = 80, tabWidth = 4, noFormatting = F, kp = '');
stringFromProperty = function(obj, o = list()) {
	o = merge.lists(stringFromPropertyDefaults, o);
	s = stringFromPropertyI(obj, 0, o);
	if (o$noFormatting) {
		s = gsub("[\n\t]", '', s);
	}
	s
}

# tokens: character vector of tokens
# ti: current token cursor (token index)
propertyFromStringRaw = function(tokens, ti = 1) {
	if (length(tokens) < 1) stop("propertyFromString: out of tokens");
	pl = if (tokens[ti] == '(') {	# we have an array here 	# ')' (bracket)
		a = NULL;
		repeat {
			ti = ti + 1;	# advance to next token
			if (ti > length(tokens) || tokens[ti] == ')') break;	# <A> empty list
			r = propertyFromStringRaw(tokens, ti);	# sub propertyList
			if (is.list(r$pl)) r$pl = list(r$pl);	# <A> concatanating of lists
			a = c(a, r$pl);
			ti = r$ti + 1;
			if (ti > length(tokens) || tokens[ti] == ')') break;	# <A> returning to list end
			if (tokens[ti] != ',') stop("propertyFromString: expected ',' or ')'");
		}
		if (ti > length(tokens) || tokens[ti] != ')') stop("propertyFromString: no array termination");
		a
	} else if (tokens[ti] == '{') {
		dict = list();
		repeat {
			ti = ti + 1;	# advance to next token
			if (ti > length(tokens) || tokens[ti] == '}') break;
			key = tokens[ti];
			if (tokens[ti + 1] != '=') stop("propertyFromString: expected '='");
			r = propertyFromStringRaw(tokens, ti + 2);
			dict[[key]] = r$pl;
			ti = r$ti + 1;
			if (tokens[ti] != ';') stop("propertyFromString: expected ';'");
		}
		if (ti > length(tokens) || tokens[ti] != '}') stop("propertyFromString: no dict termination");;
		dict
	#} elsif ($token =~ /^<(.*)>$/so) {		# we encountered data
	# <N> data not supported
	} else {	# string
		s = tokens[ti];
		if (substr(s, 1, 1) == '"') s = substr(s, 2, nchar(s) - 1);
		s
	}
	r = list(pl = pl, ti = ti);
	r
}

plStringRE = '(?:(?:[_\\/\\-a-zA-Z0-9.]+)|(?:\"(?:(?:\\\\.)*(?:[^"\\\\]+(?:\\\\.)*)*)\"))';
plCommentRE = '(?:/\\*(?:.*?)\\*/)';

propertyFromString = function(plistString, o = list()) {
	plistString = gsub(plCommentRE, '', plistString, perl = T);
	tokens = fetchRegexpr(sprintf('%s|[(]|[)]|[{]|[}]|[=]|[,]|[;]|<.*?>', plStringRE), plistString);
	pl = propertyFromStringRaw(tokens);
	pl$pl
}

#	pattern for extended plist properties
# #------------------------------------------------------------------------------------------------------------
# #PROPERTY_TAG
plExtDefaults = list( REstart = '^#[-]{80,}', REname = '^#(.*)' );
splitExtendedPlist = function(s, c = list()) {
	c = merge.lists(plExtDefaults, c);
	lines = splitString("\\n", s);
	m = Regexpr(c$REstart, lines);

	# <p> separate items
	start = which(sapply(m, length) > 0);
	startA = c(start, length(lines));
	plist = join(lines[Seq(1, start[1] - 1)], "\n");

	plistE = unlist.n(eilapply(start, function(st, i) {
		n = Regexpr(c$REname, lines[st + 1], captures  = T);
		SetNames(list(join(lines[Seq(startA[i] + 2, startA[i + 1] - 1)], "\n")), n);
	}), 1);
	return(list(propertyList = plist, interpolation = plistE));
}


propertyListTraverseRaw = function(plist, fString = identity, fData = NULL, ..., logLevel = 7) {
	Log("Plist traversal enter...", logLevel);
	plistI = if ( is.character(plist) && length(plist) == 1 ) {
		LogS(logLevel, "Plist traversal reached String [%{plist}s]");
		if (RegexprM('[\\x00-\\x08\\x80-\\x9f\\x7f-\\xff]', plist) && notE(fData)) fData(plist, ...) else
		if (notE(fString)) fString(plist, ...) else plist
	} else if (is.character(plist)) {
		LogS(logLevel, "Plist traversal reached array [1..%{N}d]", length(plist));
		sapply(plist, propertyListTraverseRaw, fString = fString, fData = fData, ...)
	} else if (is.list(plist)) {
		LogS(logLevel, "Plist traversal reached dictionary {%{ks}s}", ks = join(names(plist), ','));
		ns = sapply(names(plist), propertyListTraverseRaw, fString = fString, fData = fData, ...);
		values = lapply(plist, propertyListTraverseRaw, fString = fString, fData = fData, ...);
		SetNames(values, ns)
	} else stopS("Unknown Plist type: %{t}s", t = class(plist));
	return(plistI);
}

# c: fString => function handling strings, fArray, fDict, c: context
propertyListTraverse = function(plist, fString = identity, ...)propertyListTraverseRaw(plist, fString, ...);

plistInterpolateString = function(s, i, logLevel = 7) {
	LogS(logLevel, "Plist interpolation: %{s}s --> '%{d}s'", d = if (notE(i[[s]])) i[[s]] else s);
 	return(if (notE(i[[s]])) i[[s]] else s);
}
plistInterpolate = function(plist, i) {
	return(propertyListTraverse(plist, fString = plistInterpolateString, i = i));
}

propertyFromStringExt = function(s, c) {
	plRaw = splitExtendedPlist(s);
	return(plistInterpolate(propertyFromString(plRaw$propertyList), plRaw$interpolation));
}

