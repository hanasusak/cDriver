# cDriver
cDriver R package for finding candidate driver genes in cancers. To be able to use cDriver, user must have installation of R programing language.

Current version is a *Beta* version. We will regularly improve this package based on a users comments.

Manuscript for cDriver is in submition process, and we have preprent at bioRxiv: http://biorxiv.org/content/early/2017/04/13/059360.


# Installation from R

```Rscript
install.packages("devtools")

library(devtools)

install_github("hanasusak/cDriver")

library(cDriver)

# And now is read to use!
```

# Installation from Command line (terminal)
```Shell
curl -L https://api.github.com/repos/hanasusak/cDriver/tarball > cDriver.tar.gz

R CMD INSTALL cDriver.tar.gz

# cDriver package should be installed now.
```

or

Click at this [link](https://api.github.com/repos/hanasusak/cDriver/tarball) and download cDriver package.
Then you rename it as cDriver.tar.gz and you can install it from Shell as mentioned before:
```Shell
R CMD INSTALL cDriver.tar.gz

```

or

Go to this [link](https://github.com/hanasusak/cDriver_tools), what is another github repository.
There you can download wrapper for cDriver to run analysis from command line wihout any programing knowledge.

# License
This project is licensed under the terms of the MIT license.

