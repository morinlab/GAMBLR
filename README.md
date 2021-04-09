# GAMBLR - an R package with convenience functions for working with GAMBL results

## Getting Started


This package relies heavily on a MariaDB/MySQL database hosted on the GSC network. Because of this, you can't use many of the functions outside of the network. You will also need a config file that contains the credentials that give you access to one of the databases. If you have the appropriate approvals, you can make a copy of the credentials file in the home directory of another GAMBL user. 

```
cp /home/rmorin/.my.cnf ~/.my.cnf && chgrp morinlab ~/.my.cnf && chmod 750 ~/.my.cnf
```

It is extremely important that you ensure the file has the right unix group (`morinlab`) and is not readable by "all". The output of `ls -l` should look exactly as shown below (with your username shown instead of rmorin). 

```
-rwxr-x--- 1 rmorin morinlab 54 Mar 31 08:41 /home/rmorin/.my.cnf
```
