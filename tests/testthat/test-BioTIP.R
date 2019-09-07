context("protein coding regions")
library("BioTIP")

data("gencode")
data("intron")
data("ILEF")
data("cod")

gencode_gr = GRanges(gencode)
ILEF_gr = GRanges(ILEF)
cod_gr = GRanges(cod)
intron_gr =GRanges(intron)

test_that("protein coding regions", {
    expect_equal(getBiotypes(ILEF_gr,gencode_gr,intron_gr),
                 getBiotypes(ILEF_gr,gencode_gr,intron_gr))

})
