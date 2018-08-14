library(ISBF)

data(CGHDisease1)
cgh = cghISBF(CGHDisease1$data,CGHDisease1$chromosome,CGHDisease1$nucposi,s=1,K=100)
