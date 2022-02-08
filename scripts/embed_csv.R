embed_csv=function(x= mtcars, filename= "file.csv", label= "Download as CSV") {
  
  # Create encoded Base64 datastream 
  encode_data= function(x){
    write.csv2(x, "./file.csv")
    enc= sprintf('data:text/csv;base64,%s', openssl::base64_encode(paste0(readLines("./file.csv"), collapse="\n")) )
    unlink("./file.csv")
    return(enc)
  }
  
  # String result ready to be placed in rmarkdown
  paste0("<a download='", filename, "' href=", encode_data(x), ">", label, "</a>")
  
}
