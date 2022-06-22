#
#	plausibilityRegion.R
#Wed 08 Jun 2022 12:17:48 PM CEST

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
contour2shape = function(ct) {
	coords0 = cbind(ct$x, ct$y);
	coords = rbind(coords0, coords0[1, , drop = F]);
	Polygon(coords)
}


extendNodeBy = function(n, nBy, by = 'x')merge(n, nBy, by = by)
extendNodeTry = function(n, ns, by = 'x')
	do.call(rbind, lapply(ns, function(nBy)extendNodeBy(n, nBy, by)))
first = function(v)v[1]
reverse = function(m)return(m[Seq(nrow(m), 1, ascending = FALSE), , drop = F])

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
print(list(nodesExtendBy_nodesC = nodesC, nodesExtendBy_nodesU = nodesU, nodeN = nodeN ));
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

# make first and last coordinate same
closePathAppend = function(coords) {
	if (all(coords[1, ] == coords[nrow(coords), ])) return(coords);	# already closed
	return(rbind(unique(coords), coords[1, , drop = F]));
}

# nodes list of segments (df each)
# assume x or y coordinate same for segments
#	connect these segments
closePath = function(nodes) {
	if (length(nodes) == 0) return(nodes);
	if (length(nodes) > 1) warning('Contour leaves boundingbox. Region might undercover.');
	path = closePathRaw(nodes)[[ 1 ]];
	print(list(closePath_path = path));
	return(closePathAppend(path));
}

completeeLevelSet = function(cl, contour) {
	lvlsI = which(list.kpu(cl, 'level') == contour);
	nodeList = lapply(lvlsI, function(i)Df(x = cl[[i]]$x, y = cl[[i]]$y));
	print(list(completeeLevelSet_nodeList = nodeList));
	return(closePath(nodeList));
}

lengthnullomit = function(l)filterList(l, function(e)length(e) > 0)

completeLevelSets = function(cts, contour = 40) {
	cts = lengthnullomit(lapply(cts, completeeLevelSet, contour = contour));
	if (length(cts) == 0) return(NULL);
	plgs = lapply(cts, contour2shape);
	splgs = SpatialPolygons(list(Polygons(plgs, ID = 'contour')));
	contourP = unionSpatialPolygons(splgs, 'contour');
	contour = contourP@polygons[[1]]@Polygons[[1]]@coords;
	return(contour);
}
