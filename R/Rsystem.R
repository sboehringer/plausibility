#
#	Rsystem.R
#Mon 27 Jun 2005 10:51:30 AM CEST 

#
#	<par> file handling
#

# <!><N> works only on atomic path
# <!> 5.1.2016: trailing slash leads to basename of ""
splitPath = function(path, removeQualifier = TRUE, ssh = FALSE, skipExists = FALSE) {
	if (is.null(path)) return(NULL);
	if (removeQualifier) {
		q = fetchRegexpr('(?<=^\\[).*?(?=\\]:)', path);
		if (length(q) > 0) path = substr(path, nchar(q) + 4, nchar(path));
	}
	sshm = list(user = '', host = '', userhost = '');
	if (ssh) {
		sshm = fetchRegexpr('^(?:(?:([a-z]\\w*)(?:@))?([a-z][\\w.]*):)?(.*)', path,
			ignore.case = TRUE, captureN = c('user', 'host', 'path'))[[1]];
		sshm$userhost = if (sshm$user != '') sprintf('%s@%s', sshm$user, sshm$host) else sshm$host;
		path = sshm$path;
	}

	#path = "abc/def.ext";
	#r.base = basename(path);
	#re = "([^.]*$)";
	#r = gregexpr(re, r.base)[[1]];
	#ext = substr(r.base, r[1], r[1] + attr(r, "match.length")[1] - 1);
	#ext = firstDef(fetchRegexpr('(?<=\\.)[^/.]+\\Z', path), '');
	ext = fetchRegexpr('(?<=\\.)[^/.]+\\Z', path);
	# take everything before ext and handle possible absence of '.'
	#base = substr(r.base, 1, r[1] - 1 - (ifelse(substr(r.base, r[1] - 1, r[1] - 1) == '.', 1, 0)));
	# reduce to file.ext
	Nchar = nchar(path);
	# replace leading '~'
	if (path == '~') {
		path = Sys.getenv('HOME');
	} else if (Nchar > 1 && substr(path, 1, 2) == '~/') {
		path = join(c(Sys.getenv('HOME'), substring(path, 2)), '');
	}
	Nchar = nchar(path);
	if (Nchar != 0 && substr(path, Nchar, Nchar) == '/') {
		base = '';
		dir = substr(path, 1, Nchar - 1);
	} else {
		base = basename(path);
		dir = dirname(path);
	}
	# base as yet still contains the file extension
	file = base;
	# chop off extension if present
	if (length(fetchRegexpr('\\.', base)) > 0) base = fetchRegexpr('\\A.*(?=\\.)', base);
	
	#pieces = regexpr(re, path, perl = TRUE);
	pieces = fetchRegexpr('([^.]+)', path);
	isAbsolute = Nchar != 0 && substr(path, 1, 1) == '/';
	# <N> disk is accessed
	exists = if (!skipExists) File.exists(path, host = sshm$userhost, ssh = FALSE) else NA;
	nonempty = exists && (file.info(path)$size > 0);
	ret = c(list(
		dir = dir,
		base = base,
		path = path,
		fullbase = sprintf("%s/%s", dir, base),
		ext = ext,
		file = file,
		isAbsolute = isAbsolute,
		absolute = if (isAbsolute) path else sprintf('%s/%s', getwd(), path),
		# fs properties
		exists = exists, nonempty = nonempty,
		# remote
		is.remote = !(sshm$user == '' && sshm$host == ''),
			user = sshm$user, host = sshm$host, userhost = sshm$userhost
	), if (removeQualifier && length(q) > 0)
		list(qualifier = q, qualifierFull = Sprintf('[%{q}s]:')) else
		list(qualifier = NA, qualifierFull = ''));
	ret
}
path.absolute = absolutePath = function(path, home.dir = TRUE, ssh = TRUE) {
	path = splitPath(path, ssh = ssh)$path;
	if (home.dir && nchar(path) >= 2 && substr(path, 1, 2) == "~/")
		path = sprintf("%s/%s", Sys.getenv('HOME'), substr(path, 3, nchar(path)));
	if (nchar(path) > 0 && substr(path, 1, 1) == "/") path else sprintf("%s/%s", getwd(), path)
}
pathSimplify = function(path) {
	path = gsub('(^./|/.(?=/|$))', '', path, perl = TRUE);
	return(path);
}

tempFileName = function(prefix, extension = NULL, digits = 6, retries = 5, inRtmp = FALSE,
	createDir = FALSE, home.dir = TRUE, doNotTouch = FALSE) {
	ext = if (is.null(extension)) '' else sprintf('.%s', extension);
	path = NULL;
	if (inRtmp) prefix = sprintf('%s/%s', tempdir(), prefix);
	if (home.dir) prefix = path.absolute(prefix, home.dir = home.dir);
	for (i in 1:retries) {
		path = sprintf('%s%0*d%s', prefix, digits, floor(runif(1) * 10^digits), ext);
		LogS(5, 'tempFileName trying path: %{path}s');
		if (!File.exists(path)) break;
	}
	if (File.exists(path))
		stop(sprintf('Could not create tempfile with prefix "%s" after %d retries', prefix, retries));
	# potential race condition <N>
	if (createDir)
		Dir.create(path, recursive = TRUE) else
		if (!doNotTouch) writeFile(path, '', mkpath = TRUE, ssh = TRUE);
	# # old implementation
	#path = tempfile(prefix);
	#cat('', file = path);	# touch path to lock name
	#path = sprintf("%s%s%s", path, ifelse(is.null(extension), "", "."),
	#	ifelse(is.null(extension), "", extension));
	Log(sprintf('Tempfilename:%s', path), 5);
	path
}
dirList = function(dir, regex = TRUE, case = TRUE, path = FALSE, absolute = FALSE) {
	base = splitPath(dir)$dir;
	files = list.files(base);
	if (regex) {
		re = splitPath(dir)$file;
		files = files[grep(re, files, perl = TRUE, ignore.case = !case)];
	}
	prefix = if (absolute) splitPath(base)$absolute else base;
	if (absolute || path) files = paste(prefix, files, sep = '/');
	files
}
list_files_with_exts = function(path, exts, full.names = TRUE)
	list.files(path, pattern = Sprintf('.(%{Exts}s)$', Exts = join(exts, '|')), full.names = full.names);

list_files_with_base = function(path, exts, full.names = TRUE) {
	sp = splitPath(path);
	list.files(sp$dir,
		pattern = Sprintf('^%{base}s.(%{Exts}s)$', base = sp$base, Exts = join(exts, '|')),
		full.names = full.names
	);
}

#
#	<p> file manipulation
#

File.exists = function(path, host = '', agent = 'ssh', ssh = TRUE) {
	if (ssh) {
		sp = splitPath(path, skipExists = TRUE, ssh = TRUE);
		host = sp$userhost;
		path = sp$path;
	}
	r = if (!is.null(host) && host != '') {
		ret = system(sprintf('%s %s stat %s >/dev/null 2>&1', agent, host, qs(path)));
		ret == 0
	} else file.exists(path);
	r
}

File.copy_raw = function(from, to, ...,
	overwrite = FALSE, recursive = FALSE, agent = 'scp', logLevel = 6, ignore.shell = TRUE,
	symbolicLinkIfLocal = TRUE) {
	spF = splitPath(from, ssh = TRUE);
	spT = splitPath(to, ssh = TRUE);
	is.remote.f = spF$is.remote || spF$host == 'localhost';
	is.remote.t = spT$is.remote || spT$host == 'localhost';

	r = if (!is.remote.f && !is.remote.t) {
		if (symbolicLinkIfLocal) {
			LogS(4, 'Symlinking "%{from}s --> %{to}s', from = spF$path, to = spT$path);
			file.symlink(spF$path, spT$path, ...);
		} else {
			LogS(4, 'Copy "%{from}s --> %{to}s', from = spF$path, to = spT$path);
			file.copy(spF$path, spT$path, recursive = recursive, ..., overwrite = overwrite);
		}
	} else {
		# <A> assume 'to' to be atomic
		cmd = sprintf('%s %s %s %s %s',
			agent,
			ifelse(recursive, '-r', ''),
			paste(sapply(from, qs), collapse = ' '),
			qs(to),
			ifelse(ignore.shell, '>/dev/null', '')
		);
		System(cmd, logLevel);
	}
	r
}

File.copy = function(from, to, ..., recursive = FALSE, agent = 'scp', logLevel = 6, ignore.shell = TRUE,
	symbolicLinkIfLocal = TRUE) {
	if (is.null(from)) return(NULL);
	pairs = cbind(from, to);
	r = apply(pairs, 1, function(r) {
		File.copy_raw(r[1], r[2], ...,
			recursive = recursive, agent = agent, logLevel = logLevel,
			ignore.shell = ignore.shell, symbolicLinkIfLocal = symbolicLinkIfLocal)
	})
	r
}

File.remove = function(path, ..., agent = 'ssh', ssh = TRUE, logLevel = 6) {
	r = if (ssh) {
		sp = splitPath(path, skipExists = TRUE, ssh = TRUE);
		host = sp$userhost;
		rpath = sp$path;
		if (File.exists(path, ssh = TRUE))
			System(sprintf('rm %s', join(sapply(rpath, qs))), pattern = agent,
				ssh_host = host, logLevel = logLevel);
	} else if (file.exists(path)) file.remove(path, ...);
	r
}

# <i> remote operations
File.symlink = function(from, to, replace = TRUE, agent = 'ssh', ssh = FALSE, logLevel = 6,
	warnings = FALSE) {
	r = if (ssh) {
		sp = splitPath(from, skipExists = TRUE, ssh = TRUE);
		host = sp$userhost;
		rpath = sp$path;
		# <!><i>
		stop('not implmenented');
	} else {
		Log(sprintf('symlink %s -> %s', qs(from), qs(to)), logLevel);
		if (replace && file.exists(to)) file.remove(to);
		if (warnings)
			file.symlink(from, to) else
			suppressWarnings(file.symlink(from, to))
	}
	r
}


# <!> only atomic path
#	treatAsFile: causes Dir.create to split off last path-component
Dir.create = function(path, ..., recursive = FALSE, agent = 'ssh', logLevel = 6,
	ignore.shell = TRUE, allow.exists = TRUE, treatPathAsFile = FALSE) {
	sp = splitPath(path, ssh = TRUE);
	# ignore last path-component
	if (treatPathAsFile) {
		sp$path = sp$dir;
		Log(sprintf('creating path %s', sp$path), 4);
	}
	if (sp$is.remote) {
		System(sprintf('ssh %s mkdir %s %s %s',
			sp$userhost,
			if (recursive) '--parents' else '',
			paste(sapply(sp$path, qs), collapse = ' '),
			if (ignore.shell) '2>/dev/null' else ''
		), logLevel);
	} else {
		if (allow.exists && !file.exists(sp$path)) dir.create(sp$path, ..., recursive = recursive);
	}
}

Save = function(..., file = NULL, symbolsAsVectors = FALSE, mkpath = TRUE, envir = parent.frame(1)) {
	sp = splitPath(file, ssh = TRUE);
	localPath = if (sp$is.remote) tempfile() else file;
	if (mkpath) { Dir.create(file, recursive = TRUE, treatPathAsFile = TRUE); }
	r = if (symbolsAsVectors) {
		do.call('save', c(as.list(c(...)), list(file = localPath)), envir = envir);
	} else save(..., file = localPath, envir = envir);
	if (sp$is.remote) File.copy(localPath, file);
	r
}
Load = function(..., file = NULL, Load_sleep = 0, Load_retries = 3, envir = parent.frame(1), logLevel = 6) {
	sp = splitPath(file, ssh = TRUE);
	localPath = if (sp$is.remote) tempfile() else file;
	r = NULL;
	for (i in 1:Load_retries) {
		if (sp$is.remote) {
			if (!File.exists(file)) {
				Sys.sleep(Load_sleep);
				next;
			}
			File.copy(file, localPath, logLevel = logLevel);
		}
		r = try(load(..., file = localPath, envir = envir));
		if (class(r) == 'try-error' && Load_sleep > 0) Sys.sleep(Load_sleep) else break;
	}
	if (is.null(r)) stop(sprintf('could not Load %s', file));
	if (class(r) == 'try-error') stop(r[1]);
	r
}

#
#	create output file names
# output = list(prefix = "results/pch", extension = "pdf", tag = "20100727");
fileName = function(output, extension = NULL, subtype = NULL) {
	if (is.null(output)) return(NULL);
	if (is.null(output$prefix)) return(NULL);
	subtype = firstDef(subtype, output$subtype, "");
	if (subtype != "") subtype =  sprintf("%s-", subtype);
	r = sprintf("%s-%s%s.%s", output$prefix, subtype, output$tag,
		firstDef(extension, output$extension, ""));
	Log(r, 4);
	r
}
#.globalOutput = list(prefix = 'results/20120126-');
#save(r, file = .fn('simulation', 'RData'))
.globalOutputDefault = .globalOutput = list(prefix = '', tag = NULL, tagFirst = FALSE);
GlobalOutput_env__ = new.env();
# .fn.set(prefix = 'results/predictionTesting-')
# @par prefix character, start path name with this character string
# @par tag character, add dashed string to all files (defaults to appending to filename)
# @par tagFirst boolean, put tag as a prefix to the file name instead
.fn.set = function(...) {
	.globalOutput = merge.lists(.globalOutputDefault, list(...));
	assign('.globalOutput', .globalOutput, envir = GlobalOutput_env__);
}
# create output file name on globalOptions
.fn = function(name, extension = '', options = NULL) {
	o = merge.lists(.globalOutputDefault, .globalOutput,
		get('.globalOutput', envir = GlobalOutput_env__), options);
	# construct plain filename
	pathes = sprintf('%s%s%s%s', o$prefix, name, ifelse(extension == '', '', '.'), extension);
	fn = sapply(pathes, function(path) {
		sp = splitPath(path);
		# <p> dir
		if (!file.exists(sp$dir)) dir.create(sp$dir);
		# <p> tag
		ext = firstDef(sp$ext, '');
		fn = if (!is.null(o$tag)) {
			if (o$tagFirst) {
				sprintf('%s/%s-%s%s%s', sp$dir, o$tag, sp$base, ifelse(ext == '', '', '.'), ext)
			} else { sprintf('%s/%s-%s%s%s', sp$dir, sp$base, o$tag, ifelse(ext == '', '', '.'), ext) };
		} else sprintf('%s/%s%s%s', sp$dir, sp$base, ifelse(ext == '', '', '.'), ext);
		fn
	});
	avu(fn)
}
.fn.pushPrefix = function(prefix) {
	output = merge.lists(.globalOutput, list(prefix = sprintf('%s%s', .globalOutput$prefix, prefix)));
	assign('.globalOutput', output, envir = GlobalOutput_env__);
	.globalOutput
}
.fn.popPrefix = function(prefix) {
	output = merge.lists(.globalOutput, list(prefix = sprintf('%s/', splitPath(.globalOutput$prefix)$dir)));
	assign('.globalOutput', output, envir = GlobalOutput_env__);
	.globalOutput
}

exprInDir = function(expr, dir = '.', envir = parent.frame()) {
	prev = setwd(dir);
	on.exit(setwd(prev));
	return(eval(expr, envir = envir));
}

#
#	create consecutive files
#
# findNextFile = function(path, N = 1e2) {
# 	sp = splitPath(path);
# 	for (i in 0:N) {
# 		path = if (i > 0) with(sp, Sprintf('%{fullbase}s-%{i}d.%{ext}')) else path;
# 		if (!file.exists(path)) return(path);
# 	}
# 	stop(Sprintf('No path could be crated from base path: %{path}s'));
# }

findLastVersion = function(path, retAll = FALSE) {
	sp = splitPath(path);
	E = if (is.null(sp$ext)) '' else with(sp, Sprintf('.%{ext}s'));
	re = with(sp, Sprintf('^%{base}s-(\\d+)%{E}s$'));
	files = list.files(sp$dir, pattern = re);
	i = max(c(0, as.integer(Regexpr(re, files, captures = TRUE))));
	highest = if (i == 0) path else with(sp, Sprintf('%{fullbase}s-%{i}d%{E}s'));
	return(if (retAll) list(path = highest, version = i, ext = E) else highest);
}

findNextFile = function(path, Nmax = 1e2) {
# 	sp = splitPath(path);
# 	re = with(sp, Sprintf('^%{base}s-(\\d+).%{ext}s$'));
# 	files = list.files(sp$dir, pattern = re);
	lv = findLastVersion(path, retAll = TRUE);
	v = lv$version;
	if (v == 0 && !file.exists(path)) return(path);
	if (v >= Nmax)
		stop(Sprintf('No path could be crated from base path: %{path}s [Maximum versions exhausted: %{Nmax}d]'));
	r = with(splitPath(path), Sprintf('%{fullbase}s-%{i}d%{ext}s', i = v + 1, ext = lv$ext));
	LogS(6, 'findNextFile: %{r}s');
	return(r);
}


#
#	command argument handling
#

# default args: command line call minus command
evaluateArgs = function(c = commandArgs()[-1]) {
	is.no.option = is.na(as.integer(sapply(c, function(a)grep("^--", a))));
	#c = c[!(c == "--vanilla")];	# eliminate '--vanilla' arguments
	c = c[is.no.option];
	if (length(c) > 0) {
		eval.parent(parse(text = c[1]));
		argListString = gsub(";", ",", gsub(";$", "", c[1]));
		print(argListString);
		return(eval(parse(text = sprintf("list(%s)", argListString))));
	}
	return(NULL);
}

# default args: command line call minus command
getCommandOptions = function(c = commandArgs()[-1]) {
	is.no.option = is.na(as.integer(sapply(c, function(a)grep("^--", a))));
	#c = c[!(c == "--vanilla")];	# eliminate '--vanilla' arguments
	c = c[is.no.option];
	o = lapply(c, function(e) {
		eval(parse(text = e));
		nlapply(setdiff(ls(), 'e'), function(n)get(n))
	});
	o = unlist.n(o, 1);
	o
}

# R.pl interface

handleTriggers = function(o, triggerDefinition = NULL) {
	if (is.null(triggerDefinition)) triggerDefinition = rget('.globalTriggers');
	if (!is.list(o) || is.null(triggerDefinition)) return(NULL);
	for (n in names(triggerDefinition)) {
		if (!is.null(o[[n]])) triggerDefinition[[n]](o$args, o);
	}

}

#
#	<p> extended system call
#

# Example of patterns:
# 	System(cmd, 5, patterns = c('cwd', 'qsub', 'ssh'),
# 		cwd = sp$path, ssh_host = sp$userhost,
# 		qsubPath = sprintf('%s/qsub', sp$path), qsubMemory = self@config$qsubRampUpMemory);


.System.fileSystem = list(
	#tempfile = function(prefix, ...)tempfile(splitPath(prefix)$base, tmpdir = splitPath(prefix)$dir, ...),
	tempfile = function(prefix, ...)tempFileName(prefix, ...),
	readFile = function(...)readFile(...)
);
.System.patterns = list(
	default = list(pre = function(cmd, ...)cmd, post = function(spec, ret, ...)list()	),
	qsub = list(pre = function(cmd, spec,
		jidFile = spec$fs$tempfile(sprintf('/tmp/R_%s/qsub_pattern', Sys.getenv('USER'))),
		qsubOptions = '',
		waitForJids = NULL, ...) {
		Dir.create(jidFile, treatPathAsFile = TRUE);
		waitOption = if (is.null(waitForJids)) '' else
			sprintf('--waitForJids %s', join(waitForJids, sep = ','));
		message(cmd);
		ncmd = sprintf('qsub.pl --type ogs --jidReplace %s %s --unquote %s -- %s',
			jidFile, waitOption, qsubOptions, qs(cmd));
		message(ncmd);
		spec = list(cmd = ncmd, jidFile = jidFile);
		spec
	},
	post = function(spec, ret, ...) { list(jid = as.integer(spec$fs$readFile(spec$jidFile))) }
	),

	qsub_slurm = list(pre = function(cmd, spec,
		jidFile = spec$fs$tempfile(sprintf('/tmp/R_%s/qsub_pattern', Sys.getenv('USER'))),
		qsubOptions = '',
		waitForJids = NULL, ...) {
		Dir.create(jidFile, treatPathAsFile = TRUE);
		waitOption = if (is.null(waitForJids)) '' else
			sprintf('--waitForJids %s', join(waitForJids, sep = ','));
		message(cmd);
		ncmd = sprintf('qsub.pl --type slurm --jidReplace %s %s --unquote %s -- %s',
			jidFile, waitOption, qsubOptions, qs(cmd));
		message(ncmd);
		spec = list(cmd = ncmd, jidFile = jidFile);
		spec
	},
	post = function(spec, ret, ...) { list(jid = as.integer(spec$fs$readFile(spec$jidFile))) }
	),

	cwd = list(pre = function(cmd, spec, cwd = '.', ...) {
		ncmd = sprintf('cd %s ; %s', qs(cwd), cmd);
		spec = list(cmd = ncmd);
		spec
	},
	post = function(spec, ret, ...) { list() }
	),
	# <i> stdout/stderr handling
	ssh = list(pre = function(cmd, spec, ssh_host = 'localhost', ssh_source_file = NULL, ...,
		ssh_single_quote = TRUE) {
		if (!is.null(ssh_source_file)) {
			cmd = sprintf('%s ; %s',
				join(paste('source', qs(ssh_source_file), sep = ' '), ' ; '), cmd);
		}
		fmt = if (ssh_single_quote) 'ssh %{ssh_host}s %{cmd}q' else 'ssh %{ssh_host}s %{cmd}Q';
		spec = list(cmd = Sprintf(fmt));
		spec
	},
	fs = function(fs, ..., ssh_host) {
		list(
			tempfile = function(prefix, ...) {
				Log(sprintf('tempfile ssh:%s', prefix), 1);
				r = splitPath(tempFileName(sprintf('%s:%s', ssh_host, prefix), ...), ssh = TRUE)$path;
				Log(sprintf('tempfile ssh-remote:%s', r), 1);
				r
			},
			readFile = function(path, ...)readFile(sprintf('%s:%s', ssh_host, path), ..., ssh = TRUE)
		);
	},
	post = function(spec, ret, ...) { list() }
	)
);
#
#	a system call (c.f. privatePerl/TempFilenames::System)
#
System_env__ <- new.env();
assign(".system.doLogOnly", FALSE, envir = System_env__);

System = function(cmd, logLevel = get('DefaultLogLevel', envir = Log_env__),
	doLog = TRUE, printOnly = NULL, return.output = FALSE,
	pattern = NULL, patterns = NULL, ..., return.cmd = FALSE, return.error = FALSE) {
	# prepare
	if (!exists(".system.doLogOnly", envir = System_env__))
		assign(".system.doLogOnly", FALSE, envir = System_env__);
	doLogOnly = ifelse (!is.null(printOnly), printOnly, get('.system.doLogOnly', envir = System_env__));

	# pattern mapping
	fs = .System.fileSystem;
	if (!is.null(patterns)) {
		spec = list();
		# map file accesses
		for (pattern in rev(patterns)) {
			fsMapper = .System.patterns[[pattern]]$fs;
			if (!is.null(fsMapper)) fs = fsMapper(fs, ...);
			spec[[length(spec) + 1]] = list(fs = fs);
		}
		# wrap commands into each other
		for (i in 1:length(patterns)) {
			spec[[i]] = merge.lists(spec[[i]], .System.patterns[[patterns[[i]]]]$pre(cmd, spec[[i]], ...));
			cmd = spec[[i]]$cmd;
		}
	} else if (!is.null(pattern)) {
		spec = .System.patterns[[pattern]]$pre(cmd, list(fs = fs), ...);
		spec$fs = fs;	# manually install fs
		cmd = spec$cmd;
	}
	# redirection (after patterns) <A>
	if (return.output & !doLogOnly) {
		tmpOutput = tempfile();
		cmd = sprintf("%s > %s", cmd, tmpOutput);
	}
	if (return.error & !doLogOnly) {
		tmpError = tempfile();
		cmd = sprintf("%s 2> %s", cmd, tmpError);
	}
	# logging
	if (doLog){ Log(sprintf("system: %s", cmd), logLevel); }
	# system call
	ret = NULL;
	if (!doLogOnly) ret = system(cmd);
	# return value
	r = list(error = ret);
	if (return.output & !doLogOnly) {
		r = merge.lists(r, list(output = readFile(tmpOutput)));
	}
	if (return.error & !doLogOnly) {
		r = merge.lists(r, list(output.err = readFile(tmpError)));
	}
	# postprocess
	if (!doLogOnly) if (!is.null(patterns)) {
		for (i in rev(1:length(patterns))) {
			r = merge.lists(r, .System.patterns[[patterns[[i]]]]$post(spec[[i]], ret, ...));
		}
	} else if (!is.null(pattern)) {
		r = merge.lists(r, .System.patterns[[pattern]]$post(spec, ret, ...));
	}
	if (return.cmd) r$command = cmd;
	# simplified output
	if (!return.output && !return.cmd && !return.error && is.null(pattern)) r = r$error;
	r
}
SystemS = function(cmd, logLevel = get('DefaultLogLevel', envir = Log_env__),
	doLog = TRUE, printOnly = NULL, return.output = FALSE, return.cmd = FALSE, ..., envir = parent.frame()) {

	cmd = Sprintf(cmd, ..., envir = envir);
	System(cmd, logLevel, doLog, printOnly, return.output, return.cmd = return.cmd);
}

qsub_wait_function = function(r, ...) {
	ids = if (is.list(r[[1]]) & !is.null(r[[1]]$jid)) list.kp(r, 'jid', do.unlist = TRUE) else r$jid;
	idsS = if (length(ids) == 0) '' else paste(ids, collapse = ' ');
	System(sprintf('qwait.pl %s', idsS), ...);
}
	
# wait on job submitted by system
.System.wait.patterns = list(
	default = function(r, ...)(NULL),
	qsub = qsub_wait_function,
	qsub_slurm = qsub_wait_function
);
System.wait = function(rsystem, pattern = NULL, ...) {
	r = if (!is.null(pattern)) .System.wait.patterns[[pattern]](rsystem, ...) else NULL;
	r
}

System.SetDoLogOnly = function(doLogOnly = FALSE) {
	assign(".system.doLogOnly", doLogOnly, envir = System_env__);
}

#
#	<p> io
#

# Capture.ouput(..., type = c('input', 'output', 'merged', 'all', 'none', 'discard'), split = c('input', 'output', 'all', 'none'), append = c('input', 'output', 'all', 'none'), return = TRUE)
silence = function(expr, verbose = FALSE) {
	if (verbose || Sys.info()['sysname'] == 'Windows') eval(expr) else {
		sink('/dev/null', type = 'output');
		sink(stdout(), type = 'message');
		on.exit({ sink(type = 'message'); sink(type = 'output'); });
		r = eval(expr);
		r
	}
}

#
#	<p> calls
#

evalCall = function(call) {
	call = callEvalArgs(call);
	do.call(call$f, call$args, envir = call$envir)
}

# envirArgs: non-functional, depracated
Do.call = function(what, args, quote = FALSE, envir = parent.frame(),
	defaultEnvir = .GlobalEnv, envirArgs = NULL, do_evaluate_args = FALSE) {
	if (is.null(envir)) envir = defaultEnvir;
	if (do_evaluate_args) args = nlapply(args, function(e)eval(args[[e]], envir = envir));
	do.call(what = what, args = args, quote = quote, envir = envir)
}

#
#	<p> file operations
#

# <A> overlap with Source; avoid dependecy with RCurl
SourceLocal = function(file, ...,
	locations = c('', '.', sprintf('%s/src/Rscripts', Sys.getenv('HOME'))),
	envir = NULL) {
	sapply(file, function(file) {
		file0 = file.locate(file, prefixes = locations);
			if (notE(envir)) sys.source(file = file0, envir = envir, ...) else source(file = file0, ...)
	})
}

# on the fly activation of package w/o installation
#	SourcePackage('~/src/Rprivate/Packages/plausibility/plausibility.R');
SourcePackage = function(defFile, ...,
	locations = c('', '.', sprintf('%s/src/Rscripts', Sys.getenv('HOME'))),
	envir = NULL) {

	dir = splitPath(defFile)$dir;
	tmpEnv = new.env();
	SourceLocal(defFile, locations = locations, envir = tmpEnv);
	files = c(defFile, get('packageDefinition', tmpEnv)$files);

	SourceLocal(files, locations = c(dir, locations), envir = envir);
}


#' Return absolute path for name searched in search-pathes
#'
#' Search for pathes.
#'
#' @param path path (segment) to be located in standard locations
#' @param prefixes prefixes to be prepended to path to check existance
#' @param normalize boolean to inidcate whether a normalized path should be returned (absolute)
#' @param home boolean to indicate whether starting prefix '~' should be interpolated to the home folder
#' @param as.dirs assume that prefixes are pathes, i.e. a slash will be put between path and prefix
#' @param force enforces that path and prefix are always joined, otherwise if path is absolute no prefixing is performed
#' @return character vector with path to file or NULL if no file could be located
file.locate = function(path, prefixes = NULL, normalize = TRUE, as.dirs = TRUE, force = FALSE, home = TRUE) {
	if (!force && substr(path, 1, 1) == '/') return(path);
	if (substr(path, 1, 1) == '~' && home) {
		path = path.absolute(path, home.dir = TRUE);
		if (!force) return(path);
	}
	if (is.null(prefixes)) prefixes = if (as.dirs) '.' else '';
	sep = ifelse(as.dirs, '/', '');
	for (prefix in prefixes) {
		npath = sprintf('%s%s%s', prefix, sep, path);
		if (normalize) npath = path.absolute(npath);
		if (file.exists(npath)) return(npath);
	}
	NULL
}

#' Read content of file and return as character object.
#' 
#' @param path Path to the file to be read.
#' @param prefixes Search for file by prepending character strings from
#' prefixes.
#' @param normalize Standardize pathes.
#' @param ssh Allow pathes to remote files in \code{scp} notation.
#' @author Stefan Böhringer <r-packages@@s-boehringer.org>
#' @return character vector containing the file content
#' @keywords io input
# #' @examples
# #' \dontrun{
# #'   parallel8 = function(e) log(1:e) %*% log(1:e);
# #'   cat(readFile(tempcodefile(parallel8)));
# #' }
# prefixes only supported locally <!>
readFile = function(path, prefixes = NULL, normalize = TRUE, ssh = FALSE) {
	s = splitPath(path, ssh = ssh);
	r = if (s$is.remote) {
		tf = tempfile();
		File.copy(path, tf);
		readChar(tf, nchars = as.list(file.info(tf)[1,])$size);
	} else {
		if (!is.null(prefixes)) path = file.locate(path, prefixes, normalize);
		readChar(path, nchars = as.list(file.info(path)[1,])$size);
	}
	r
}

writeFile = function(path, str, mkpath = FALSE, ssh = FALSE) {
	s = splitPath(path, ssh = ssh);
	if (s$is.remote) {
		Dir.create(sprintf('%s:%s', s$userhost, s$dir), recursive = mkpath);
		tf = tempfile();
		out = file(description = tf, open = 'w', encoding='UTF-8');
			cat(str, file = out, sep = "");
		close(con = out);
		File.copy(tf, path);
	} else {
		if (mkpath) {
			if (!file.exists(s$dir)) dir.create(s$dir, recursive = TRUE);
		}
		out = file(description = path, open = 'w', encoding='UTF-8');
			cat(str, file = out, sep = "");
		close(con = out);
	}
	path
}

isURL = function(path)(length(grep("^(ftp|http|https|file)://", path)) > 0L)

#
#	<p> helper functions readTable/writeTable
#

compressPathBz2 = function(pathRaw, path, doRemoveOrig = TRUE) {
	cmd = Sprintf("cat %{pathRaw}q | bzip2 -9 > %{path}q");
	r = System(cmd, 5);
	if (doRemoveOrig && !get('.system.doLogOnly', envir = System_env__)) file.remove(pathRaw);
	r
}
compressPath = function(pathRaw, path, extension = NULL, doRemoveOrig = TRUE) {
	if (is.null(extension)) return(path);
	compressor = get(Sprintf('compressPath%{extension}u'));
	r = compressor(pathRaw, path, doRemoveOrig = doRemoveOrig);
	r
}
decompressPathBz2 = function(path, pathTmp, doRemoveOrig = FALSE) {
	cmd = Sprintf("cat %{path}q | bunzip2 > %{pathTmp}q");
	r = System(cmd, 5);
	if (doRemoveOrig && !get('.system.doLogOnly', envir = System_env__)) file.remove(path);
	r
}
decompressPath = function(path, pathTmp, extension = NULL, doRemoveOrig = FALSE) {
	if (is.null(extension)) return(path);
	decompressor = get(Sprintf('decompressPath%{extension}u'));
	r0 = decompressor(path, pathTmp, doRemoveOrig = doRemoveOrig);
	r = list(destination = pathTmp, pathOrig = path, return = r0);
	r
}

compressedConnectionBz2 = function(path, mode = '') {
	#r = Sprintf('%{path}s.bz2');
	bzfile(path, open = mode)
}
compressedConnectionGz = function(path, mode = '') {
	gzfile(path, open = mode)
}
compressedConnection = function(path, extension = NULL, mode = '') {
	if (is.null(extension)) return(path);
	compressor = get(Sprintf('compressedConnection%{extension}u'));
	compressor(path, mode = mode)
}
compressedConnectionPath = function(conn) {
	if ('connection' %in% class(conn)) summary(conn)$description else conn
}

#
#	<p> readTable
#


#
#	<p> print
#

stdOutFromCall = function(call_) {
	tf = tempfile();
	sink(tf);
		eval.parent(call_, n = 2);
	sink();
	readFile(tf)
}

#
#	crypotgraphy/checksumming
#

# md5sumString = function(s, prefix = 'md5generator') {
# 	Require('tools');
# 	path = tempfile('md5generator');
# 	writeFile(path, s);
# 	md5 = avu(md5sum(path));
# 
# 	md5
# }
# same as above, less dpendencies
md5sumString = function(s, length = 32, ..., logLevel = 5)
	substr(
		SystemS('echo -n %{s}q | md5sum', return.output = TRUE, logLevel = logLevel)$output
	, 1, min(length, 32))
sha256sumString = function(s, length = 32, ..., logLevel = 5)
	substr(
		SystemS('echo -n %{s}q | sha256sum', return.output = TRUE, logLevel = logLevel)$output
	, 1, min(length, 64))
sha256sumPath = function(path, length = 64, ..., logLevel = 5)
	substr(
		SystemS('sha256sum %{path}q', return.output = TRUE, logLevel = logLevel)$output
	, 1, min(length, 64))

hashPathContent = function(path, type = 'sha256', length = 64,..., logLevel = 5) {
	f = get(paste0(type, 'sumPath'));
	f(path, length, ..., logLevel = logLevel)
}
hashStringContent = function(path, type = 'sha256', length = 64,..., logLevel = 5) {
	f = get(paste0(type, 'sumString'));
	f(path, length, ..., logLevel = logLevel)
}

#
#	<p> package documentation
#

#
#	<p> Rcpp helpers
#

#
#	<p> sqlite
#

#
#	<p> publishing
#

#
#	<p> quick pdf generation
#

#
#	<p> workarounds
#

clearWarnings = function()assign('last.warning', NULL, envir = baseenv())

# fix broken install from dir: create tarball -> install_local
Install_local = function(path, ..., tarDir = tempdir()) {
	sp = splitPath(path);
	pkgPath = Sprintf('%{tarDir}s/%{base}s.tar.gz', sp);
	# dir component is containing folder
	System(Sprintf('cd %{dir}Q ; tar czf %{pkgPath}Q %{file}Q', sp), 2);
	#lib = list(...)$lib;
	#libLocation = if (is.null(lib)) 'default location' else lib;
	#LogS(4, 'Installing to lib:%{libLocation}s');
	#print(Sprintf('Installing to lib:%{libLocation}s'));
	install_local(pkgPath, ...);
}

#
#	<p> packages
#

#
#	<p> misc linux system stuff
#

#
#	<p> random numbers
#

getRandomSeed = function(tag = date()) {
	md5 = md5sumString(join(c(getwd(), tag)));
	is = hex2ints(md5);
	seed = is[1];
	for (i in 2:length(is)) { seed = bitwXor(seed, is[i]); }
	seed
}

#
#	<p> Reporting
#

#
#	<p> stop
#

stopS = function(str, ...)stop(Sprintf(str, ...));

#
#	<p> debugging
#

# r__: return printed values as list
dprint = function(..., r__ = TRUE) {
	vs = as.character(as.list(substitute(list(...)))[-1]);
	ns = names(list(...));
	Ns = if (is.null(ns)) vs else ifelse(ns == '', vs, ns);
	l = listKeyValue(Ns, c(...));
	print(list2df(l));
	if (r__) return(l);
}

debugOn = function()options(error = recover);

#
#	<p> file system
#

normalizePath = function(p) {
	p = gsub('^~', Sys.getenv('HOME'), p);
	p = gsub('(?:g)//', '/', p, perl = TRUE);
	return(p);
}

# recursive version of gsbu
gsubR = function(pattern, replacement, x, ..., Nmax = 1e3) {
	for (i in 1:Nmax) {
		xNew = gsub(pattern, replacement, x, ...);
		if (all(xNew == x)) return(x);
		x = xNew;
	}
	return(NA);	# trigger special case (consider stop) <N>
}

# <i><!> unify with normalizePath after testing 8/2020
NormalizePath = function(p) {
	p = gsub('^~', Sys.getenv('HOME'), p);
	p = gsub('//+', '/', p, perl = TRUE);
	p = gsubR('(^|/)[^/]+/[.][.]/', '/', p, perl = TRUE);
	return(p);
}
pathToHome = function(path)
	gsub(Sprintf('^%{home}s((?=/)|$)', home = Sys.getenv('HOME')), '~', path, perl = TRUE)

# how to refer to to from within from
relativePathSingle = function(from, to) {
	from = normalizePath(from);
	to = normalizePath(to);
	spF = splitPath(from);
	spT = splitPath(to);
	if (spT$isAbsolute) return(to);
	join(c(rep('..', length(splitString('/', spF$dir)) + 0), to), '/');
}
relativePath = Vectorize(relativePathSingle, c('from', 'to'));
SplitPath = function(path, ...)lapply(path, splitPath, ...);
absolutePathSingle = function(path)splitPath(path)$absolute
absolutePath = Vectorize(absolutePathSingle, c('path'));
pathSimplify = function(p)gsub('[:]', '_', p)
pathInsertPostfix = function(path, postfix, sep = '-')
	Sprintf('%{fullbase}s%{sep}s%{postfix}s.%{ext}s', splitPath(path))

# keys of input-list are folder names, use folderstring in key to create subfolders
# 	createZip(list(results = c('r/ref1.html', 'r/ref2.html')), 'r/myZip.zip', doCopy = TRUE);
# 	createZip(list(`results::sub` = c('r/ref1.html', 'r/ref2.html')), 'r/myZip.zip', doCopy = TRUE);
#	values are slash-dependend dest = 'source' copies into dest/source; dest = source/, copies into dest

createZip = function(input, output, pword, doCopy = FALSE, readmeText, readme, logOnly = FALSE,
	absoluteSymlink = FALSE, simplifyFileNames = FALSE, folderString = '::') {
	destDir = splitPath(output)$fullbase;
	Dir.create(destDir);
	if (!missing(readmeText)) writeFile(Sprintf('%{destDir}s/README'), readmeText);
	nelapply(input, function(n, e) {
		if (notE(folderString)) n = gsub(folderString, '/', n);
		subdir = join(c(destDir, n, ''), '/');
		Dir.create(subdir);
		toFiles = list.kpu(SplitPath(e), 'file');
		if (simplifyFileNames) toFiles = sapply(toFiles, pathSimplify);
		to = paste(subdir, toFiles, sep = '/');
		LogS(4, 'Copy: %{e}s --> %{to}s');
		if (doCopy) file.copy(e, to, recursive = TRUE) else {
			#from = relativePath(subdir, e);
			from = absolutePath(e);
			if (absoluteSymlink) from = NormalizePath(paste(splitPath(subdir)$absolute, from, sep = '/'));
			print(list(from = from, to = to));
			file.symlink(from, to);
		}
	});
	dir = splitPath(output)$dir;
	zip = splitPath(output)$base;
	options = '';
	if (!missing(pword)) options = Sprintf('%{options}s -P %{pword}q');
	SystemS('cd %{dir}q ; zip %{options}s -r %{zip}q.zip %{zip}q', logLevel = 1, printOnly = logOnly);
}
