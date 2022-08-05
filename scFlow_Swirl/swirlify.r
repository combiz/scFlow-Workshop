setwd("~/Documents/scFlow-Workshop")

library(devtools)
install_github("swirldev/swirlify", ref = "dev")

library(swirl)
library(swirlify)
new_lesson("SingleCellExperiment", "scFlow Swirl")

add_to_manifest()
pack_course()

install_course(
    swc_path = file.path(getwd(), "scFlow_Swirl.swc")
)

swirl()


