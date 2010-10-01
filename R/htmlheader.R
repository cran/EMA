`htmlheader` <-
function(towrite, filename) {
  write(paste("<b", towrite, "</b>"), file = filename)
}

