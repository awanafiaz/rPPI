## create_package("~/My Drive (aafiaz@uw.edu)/Github_Awan/rPPI")

library(devtools)


#use_git()

use_r("simdat")

load_all()

simdat()

check()

#use_mit_license()

use_package(package = "stats")
use_package(package = "gam")

document()

install()

library(rPPI)

df <- rPPI::simdat()

use_testthat()

use_test("simdat")

test()

## Adding to github
usethis::create_from_github(
  repo_spec = "https://github.com/awanafiaz/rPPI.git",
  destdir = "~/My Drive (aafiaz@uw.edu)/Github_Awan/rPPI"
)



