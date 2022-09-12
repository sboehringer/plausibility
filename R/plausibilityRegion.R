#
#	plausibilityRegion.R
#Wed 08 Jun 2022 12:17:48 PM CEST

extendNodeBy = function(n, nBy, by = 'x')merge(n, nBy, by = by)
extendNodeTry = function(n, ns, by = 'x')
	do.call(rbind, lapply(ns, function(nBy)extendNodeBy(n, nBy, by)))
first = function(v)v[1]
reverse = function(m)return(m[Seq(nrow(m), 1, ascending = FALSE), , drop = F])

# make first and last coordinate same
closePathAppend = function(coords) {
	if (all(coords[1, ] == coords[nrow(coords), ])) return(coords);	# already closed
	return(rbind(unique(coords), coords[1, , drop = F]));
}

pathIsClosed = function(nodes) return(all(nodes[1, ] == nodes[nrow(nodes), ]));

# <i> concatenate contours
filterContours = function(ct, level, fromBelow = TRUE) {
	lvls = list.kpu(ct, 'level');
	# <!><N> assume levels to be ordered (tb checked)
	I = if (fromBelow) which.max(lvls <= level) else which.min(lvls >= level);
	return(ct[[I]]);
}

countourStack = function(cls, level, fromBelow = TRUE) {
	r = lapply(cls, function(ct) {
		filterContours(ct, level)
	});
}
contour2shape = function(ct) Polygon(closePathAppend(ct[, c('x', 'y'), drop = F]));

polygonFromCoords = function(coords, id)Polygons(list(Polygon(coords)), as.character(id))
polygonsFromCoords = function(coordsList)
	SpatialPolygons(ilapply(coordsList, \(coords, i)polygonFromCoords(coords, i)));

PolyUnion = function(plgs) {
	poly = plgs[[1]];
	for (p in  plgs[-1]) poly = gUnion(poly, p);
	return(poly);
}

# returns list
#	rev: did the order of nodes reverse
nodesConnect_old = function(n1, n2, by = 'x') {
	if (last(n1[[ by ]]) == first(n2[[ by ]]))
		return(list(rev = FALSE, coords = rbind(n1, n2)));
	if (last(n1[[ by ]]) == last(n2[[ by ]]))
		return(list(rev = FALSE, coords = rbind(n1, reverse(n2))));
	if (first(n1[[ by ]]) == last(n2[[ by ]]))
		return(list(rev = TRUE, coords = rbind(n2, n1)));
	if (first(n1[[ by ]]) == first(n2[[ by ]]))
		return(list(rev = TRUE, coords = rbind(reverse(n1), reverse(n2))));
	stop('Could not connnect nodes');
}

nodesConnect = function(n1, n2, by = 'x') {
	if (last(n1[[ by ]]) == first(n2[[ by ]])) return(rbind(n1, n2));
	if (last(n1[[ by ]]) == last(n2[[ by ]])) return(rbind(n1, reverse(n2)));
	stop('Could not connnect nodes');
}
nodesExtendBy = function(nodes, by = 'x') {
	if (length(nodes) <= 1) return(nodes);	# nothing to extend
	nodeL = lapply(seq_along(nodes), function(i)Df(id = i, nodes[[i]]));
	nodesF = nodeL[[1]];	# first node
	nodesR = nodeL[-1];		# all other nodes

	extendsBy = extendNodeTry(nodesF[nrow(nodesF), , drop = F], nodesR, by = by);
	if (nrow(extendsBy) > 0) {
		# next node, +1 tp trams;ate frp, mpdesR tp mpdes
		nodeN = which(extendsBy$id.y[1] == list.kpu(nodesR, 'id$[[1]]')) + 1;
		nodesC = nodesConnect(nodes[[1]], nodes[[nodeN]], by);	# use original list w/o id
		nodesU = nodes[ - c(1, nodeN) ];
		#print(list(nodesExtendBy_nodesC = nodesC, nodesExtendBy_nodesU = nodesU, nodeN = nodeN ));
		return( closePathRaw(c(list(nodesC), nodesU)) )
	}
	return(NULL);
}

# assume closed loop exists
closePathRaw = function(nodes) {
	if (length(nodes) <= 1) return(nodes);	# nothing to extend

	# try to extend by x, assume path exists from last coordinate
	extension = nodesExtendBy(nodes, 'x');
	if (notE(extension)) return(extension);

	# try to extend by y, assume path exists from last coordinate
	extension = nodesExtendBy(nodes, 'y');
	if (notE(extension)) return(extension);
	stop('Could not extend path');
}

# nodes list of segments (df each)
# assume x or y coordinate same for segments
#	connect these segments
closePath = function(nodes) {
	if (length(nodes) == 0) return(nodes);
	isClosed = sapply(nodes, pathIsClosed);
	if (!all(isClosed)) warning('Contour leaves boundingbox. Region might undercover.');
	pathClosed = if (all(!isClosed)) list() else nodes[ isClosed ];
	pathClosedNew = if (!all(isClosed))
		list(closePathAppend(closePathRaw(nodes[ !isClosed ])[[1]])) else list();
	path = c(pathClosed, pathClosedNew);
	#print(list(closePath_path = path));
	return(path);
}

completedLevelSet = function(cl, contour, Nround = 4) {
	lvlsI = which.max(round(list.kpu(cl, 'level') - contour, Nround) >= 0);
	print(list(levels = round(list.kpu(cl, 'level'), Nround), I = lvlsI));
	nodeList = lapply(lvlsI, function(i)Df(x = cl[[i]]$x, y = cl[[i]]$y));
	print(list(completedLevelSet_nodeList = nodeList));
	return(closePath(nodeList));
}

length0.omit = function(l)filterList(l, function(e)length(e) > 0)

completeLevelSets = function(cts, contour = 40) {
	cts = length0.omit(unlist.n(lapply(cts, completedLevelSet, contour = contour), 1));
	contour = polygonsFromCoords(cts);
# 	splgs = SpatialPolygons(list(Polygons(plgs, ID = 'contour')));
# 	contourP = unionSpatialPolygons(splgs, 'contour');
# 	#contour = contourP@polygons[[1]]@Polygons[[1]]@coords;
# 	contour = contourP@polygons[[1]]@Polygons[[1]];
	return(contour);
}

array2vectorWrapper = function(x, y, myFun, coordsFixed = NULL, coordsIdcs = NULL, ...) {
	Npar = (length(coordsFixed) + 2);
	idcsFree = setdiff(1:Npar, coordsIdcs);
	pTempl = vector.assign(rep(NA, Npar), coordsIdcs, coordsFixed);

	apply(cbind(x, y), 1, function(r) {
		p = vector.assign(pTempl, idcsFree, r);
		myFun(p, ...)
	})
}
contour2DforFixed = function(fixedCoords, fn, seqs, free, fixedIdcs, nlevels, ..., Nround = 4) {
	x = seqs[[free[1]]];
	y = seqs[[free[2]]];
	z = outer(x, y, array2vectorWrapper, myFun = fn, coordsFixed = fixedCoords, coordsIdcs = fixedIdcs, ...);
	cl = contourLines(x, y, z, nlevels = nlevels);
	return(cl);
}

contourLinesStacked = function(ranges, N, fn, free, nlevels = 40, ...) {
	seqs = lapply(ranges, function(r)seq(r[1], r[2], length.out = N));
	if (missing(free)) free = c(length(ranges) - 1, length(ranges));
	fixed = setdiff(1:length(ranges), free);
	rs = ranges[free];

	fixedCoords = merge.multi.list(seqs[fixed]);
	cls = if (length(fixed) > 0) apply(fixedCoords, 1, contour2DforFixed, simplify = FALSE,
		fn = fn, seqs = seqs, free = free, fixcedCoords = fixed, nlevels = nlevels, ...) else
		list(contour2DforFixed(NULL, fn, seqs, free, fixed, nlevels, ...))
	return(cls);
}

FindRegion = function(f, range, level = .95, Nout = 20, ...) {
	Nargs = length(range);
	cls = contourLinesStacked(range, Nout, f, ...);
	contour = SetNames(lapply(level, \(level)completeLevelSets(cls, contour = 1 - level)), level);
	return(contour);
}

#isPointInRegion = function(p, region)
#	return(point.in.polygon(p[1], p[2], region$region@coords[, 'x'], region$region@coords[, 'y']))
isPointInRegion = function(p, region)
	if (class(region) == 'SpatialPolygons')gContains(region, SpatialPoints(matrix(p, nrow = 1))) else NA
