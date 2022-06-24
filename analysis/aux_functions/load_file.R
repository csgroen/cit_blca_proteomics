load_file <- function(file) {
    data <- eval(parse(text = load(file)))
    return(data)
}