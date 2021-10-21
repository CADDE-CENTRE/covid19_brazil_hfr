# Write function to put commas for numbers
# This could be improved: check if as.integer of input makes sense
# Then determine the length, and divide it by 3. Use rep() for the gsub pattern
write_commas <- function(vec)
{
  # Doesn t work for decimal, we can get there if need be
  # Just split about '.' and perform write_commas on the left
  
  out <- rep('', length(vec))
  tmp <- which(is.na(as.numeric(vec)))
  out[tmp] <- vec[tmp]
  idx <- which(!is.na(as.numeric(vec)))
  for (i in idx)
  {
    
    tmp <- (nchar(vec[i]) - 1)  %/% 3 
    if (tmp > 0){
      tmp1 <- paste0(rep('([0-9][0-9][0-9])', tmp), collapse = '')
      tmp2 <- paste0(',\\',seq(1:tmp), collapse = '')
      out[i] <- gsub(paste0(tmp1,'$'), tmp2, vec[i])
    }else{out[i] <- vec[i]}
  }
  return(out)
}