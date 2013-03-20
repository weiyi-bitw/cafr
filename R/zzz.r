.onAttach <- function(libname, pkgname){
	DBVer <- read.dcf(file=system.file("DESCRIPTION", package=pkgname), fields="Version")
	packageStartupMessage(
	"\n",
	"===============================\n",
	"\n",
	paste(pkgname, DBVer), "\n",
	"\n",
	"===============================\n",
	"\n",
	"Welcome to cafr package!\n\n",
	" -- Brought to you by Wei-Yi Cheng and Team Attractor Metagenes.\n\n",
	"===============================\n\n"
	)
}


