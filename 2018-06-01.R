


get_range <- function(chr,lan,width){
  m = lan%/%width
  if(!lan%%width==0){
    start = 0:m*width+1
    end   = c(1:m*width,lan)
  }else{
    start = 0:(m-1)*width+1
    end   = 1:m*width
  }
  return(GRanges(chr,IRanges(start = start,end =end)))
}
get_range(chr = "1",lan = 249250621, width = 5000)

?options()


"
warn如果是负数，则所有warning message都被忽略
warn = 0 (0为默认值)，则所有warning messages会被储存起来直到上级函数(此例中则是repeat()函数)运行结束
warn = 1，则一旦产生warning message，这条信息会被立即显示出来
warn = 2 或更大的数值， 则warning message会被立即显示并转换成error message。此例中，如果warn = 2，整个自定义函数都会被中断，提示warning message的内容，
但是会以error message的形式弹出。
"