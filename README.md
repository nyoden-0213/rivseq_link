# rivseq_link
generate river sequence considering water transfers
# calc_rivseq_link.F90
main program of FTCS algorithm
# testdata
canal_info.csv: list of information of IBWTs in Indus basin (-999 is the missing value)

diminfo_test-1deg.txt: the number of grids of the model and lontitude & latitude

nextxy.bin: MERIT Hydro (Yamazaki et al., 2019) river network data specifying the destination of each grid

# How to run the program
1. download all the files
2. compaile by running Makefile
3. run calc_rivseq_link
