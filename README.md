# GAMBLR - an R package with convenience functions for working with GAMBL results

## Getting Started


This package relies heavily on a MariaDB/MySQL database hosted on the GSC network. Because of this, you can't use many of the functions outside of the network. You will also need a config file that contains the credentials that give you access to one of the databases. If you have the appropriate approvals, you can make a copy of the credentials file in the home directory of another GAMBL user. 

```
cp /home/rmorin/.my.cnf ~/.my.cnf && chgrp morinlab ~/.my.cnf && chmod 750 ~/.my.cnf
```
