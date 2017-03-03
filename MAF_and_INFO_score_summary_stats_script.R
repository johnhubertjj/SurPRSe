
args=(commandArgs(trailingOnly = T))

if (length(args) == 0) {
  print ("No arguments supplied.")
}else{
  for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
  }
}
print(a)