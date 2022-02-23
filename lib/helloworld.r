suppressWarnings(library(optparse))

### Functions ###
sayHello <- function () {
    print('hello!')
}

sayHelloTo <- function (name) {
    print(paste0('hello', ', ', name, '!'))
}

### OptParse ###
option_list <- list(make_option(c("-n", "--name"), type="character", default="world",
		                help="Enter a name to say hello! [default: world]"),
		    make_option(c("-p", "--plot"), type="logical", default=FALSE,
				help="T/F create plot of name popularity.") )
opts <- parse_args(OptionParser(option_list=option_list), args=commandArgs(trailing=T))

### Output ###
sayHello()
sayHelloTo(opts$name)

### Create Plot ###
if (opts$plot) {
    x <- seq(1980, 2015, by=5)
    y <- runif(length(x), 0, 1)*100
    png('helloworld_plot.png')
    plot(range(x), range(y), type='n',  main=paste0('Popularity of ', opts$name, ' over time'), xlab='Year', ylab='Percent Approval')
    lines(x, y, type='b', lwd=2)
    suppressMessages(dev.off())
    browseURL('helloworld_plot.png')
    print('I created a plot for you.')
}
