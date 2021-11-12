#
#	plausibilitiy.R
#Fri May 31 14:59:59 CEST 2019
#
#	derived from plausibilityPenalized.R

packageDefinition = list(
	name = 'plausibility',
	files = c(
		c(	'plausibilityBase.R', 'plausibilityModels.R', 'plausibilityOptim.R',
			'plausibilityUnweighted.R', 'plausibilityWeighted.R', 'plausibilityPenalized.R',
			'plausibilityWeightedGaussian.R'
		),
		c('Rmeta.R', 'Rdata.R', 'Rsystem.R', 'Rfunctions.R', 'RpropertyList.R', 'Rstatistic.R')
	),
	description = list(
		title = 'Plausibility based estimation and inference',
		author = 'Stefan B\uf6hringer <r-packages@s-boehringer.org>',
		description = 'R-package to implement important methods of the plausibility framework, including goodness-of-fit testing, test model comparisons and computation of marginal and joint plausibility regions.',
		depends = c(),
		suggests = c(),
		license = 'LGPL',
		news = "0.5-4	Cleanup\n0.5-3	Code refactoring\n0.5-2	Collate fix\n0.5-1	Package dependencies\n0.5-0	Gaussian models, Penalized models\n0.4-1	Installation fixes\n0.4-0	Weighted plausibility, penalized models.\n0.3-0	Unweighted plausibility with stochastic integeration\n0.2-0	Implementation Plausibility region\n0.1-0   Initial release"
	),
	git = list(
		readme = '## Installation\n```{r}\nlibrary(devtools);\ninstall_github("sboehringer/plausibility")\n```\n',
		push = F,
		pushOnNewVersion = T,
		#remote = 'https://github.com/sboehringer/plausibility.git'
		remote = 'git@github.com:sboehringer/plausibility.git'
	)
);
